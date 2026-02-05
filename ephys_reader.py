import os
from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# ----------------------------
# Minimal "tree node" objects
# ----------------------------

class _Node:
    """Mimics the parts of heka_reader TreeNode API your app uses."""
    def __init__(self, **fields: Any):
        self.children: List[_Node] = []
        self.fields: Dict[str, Any] = dict(fields)
        for k, v in self.fields.items():
            setattr(self, k, v)

    def __getitem__(self, i: int) -> "_Node":
        return self.children[i]

    def __len__(self) -> int:
        return len(self.children)

    def __iter__(self):
        return iter(self.children)

    def get_fields(self) -> Dict[str, Any]:
        out = dict(self.fields)
        out["children"] = [ch.get_fields() for ch in self.children]
        return out


@dataclass
class _AbfSeriesDef:
    column_name: str          # e.g. "ap_rheo_series"
    code: str                 # e.g. "003"
    abf_path: str             # full path to "24o10003.abf"


# ----------------------------
# Unit conversion helpers
# ----------------------------

def _scale_to_si(values: np.ndarray, unit: str, target: str) -> np.ndarray:
    """
    Convert common electrophysiology units to SI base:
      - voltage to V
      - current to A
    """
    if unit is None:
        unit = ""
    u = unit.strip()

    # Normalize micro symbol variants
    u = u.replace("µ", "u")

    if target == "V":
        if u == "V":
            return values
        if u == "mV":
            return values * 1e-3
        if u == "uV":
            return values * 1e-6
        # Unknown: assume already V (best effort)
        return values

    if target == "A":
        if u == "A":
            return values
        if u == "mA":
            return values * 1e-3
        if u == "uA":
            return values * 1e-6
        if u == "nA":
            return values * 1e-9
        if u == "pA":
            return values * 1e-12
        # Unknown: assume already A (best effort)
        return values

    return values


# ----------------------------
# ABF-backed Bundle
# ----------------------------

class _AbfComposedBundle:
    """
    Presents multiple ABF files (one per "series") as a single HEKA-like bundle:
      pul[group][series][sweep][trace]
      data[group, series, sweep, trace] -> np.ndarray
    """
    def __init__(self, root_name: str, data_folder: str, metadata_file: str):
        self.file_name = root_name  # keep consistent with your analysis_points keys
        self._data_folder = data_folder
        self._metadata_file = metadata_file

        row = self._get_metadata_row(root_name)
        self._series_defs = self._extract_series_defs(row)

        self._pul = self._build_pul_tree(self._series_defs)
        self._data = _AbfData(self, self._series_defs)

        # Optional mapping: metadata column -> synthetic series index
        self.series_map: Dict[str, int] = {sd.column_name: i for i, sd in enumerate(self._series_defs)}

        # Eagerly populate sweep/trace children so the browser tree can show them
        # without needing to access bundle.data first.
        if len(self._series_defs) == 0:
            raise ValueError(
                f"No ABF series columns found for file_name={root_name!r}. "
                f"Expected at least one of: rmp_series, rin_series, ap_rheo_series, ap_max_series, ap_broadening_series."
            )

        for series_id in range(len(self._series_defs)):
            # This loads the ABF (cached) and creates sweep/trace nodes + NumberSweeps.
            self._data._get_abf(series_id)

    @property
    def pul(self) -> _Node:
        return self._pul

    @property
    def data(self) -> "_AbfData":
        return self._data

    def _get_metadata_row(self, root_name: str) -> pd.Series:
        df = pd.read_excel(self._metadata_file)
        if "file_name" not in df.columns:
            raise ValueError("metadata.xlsx must contain a 'file_name' column.")
        matches = df[df["file_name"].astype(str) == str(root_name)]
        if len(matches) == 0:
            raise ValueError(f"Could not find file_name={root_name!r} in metadata.")
        if len(matches) > 1:
            raise ValueError(f"Multiple rows found for file_name={root_name!r}; please make it unique.")
        return matches.iloc[0]

    def _extract_series_defs(self, row: pd.Series) -> List[_AbfSeriesDef]:
        root = str(row["file_name"])

        # Decide which columns count as “series columns”
        # Keep it strict to your known names to avoid accidentally picking numeric parameter columns.
        candidate_cols = [
            "rmp_series",
            "rin_series",
            "ap_rheo_series",
            "ap_max_series",
            "ap_broadening_series",
        ]

        series_defs: List[_AbfSeriesDef] = []

        for col in candidate_cols:
            if col not in row.index:
                continue
            val = row[col]
            if pd.isna(val):
                continue
            code = str(val).strip()

            # Accept values like 3, 3.0, "003"
            try:
                if code.replace(".", "", 1).isdigit():
                    code_int = int(float(code))
                    code = f"{code_int:03d}"
            except Exception:
                # if it's non-numeric text, keep as-is
                pass

            if code == "" or code.lower() == "nan":
                continue

            abf_name = f"{root}{code}.abf"
            abf_path = os.path.join(self._data_folder, abf_name)
            series_defs.append(_AbfSeriesDef(column_name=col, code=code, abf_path=abf_path))

        return series_defs

    def _build_pul_tree(self, series_defs: List[_AbfSeriesDef]) -> _Node:
        # Root pul node
        pul = _Node(Label="ABF_Composed", VersionName="ABF", RootText="", StartTime=0.0)
        group = _Node(Label="Group1", GroupCount=1)
        pul.children.append(group)

        # Create series nodes; each series nodes creates sweep/trace children lazily via Data object.
        for i, sd in enumerate(series_defs):
            series_label = f"Series{i + 1}"
            series = _Node(Label=series_label, NumberSweeps=0, _abf_path=sd.abf_path, _series_def=sd)
            group.children.append(series)

        return pul


class _AbfData:
    def __init__(self, bundle: _AbfComposedBundle, series_defs: List[_AbfSeriesDef]):
        self.bundle = bundle
        self._series_defs = series_defs
        self._abf_cache: Dict[str, Any] = {}  # path -> pyabf.ABF

        # Cache sweep-level arrays to avoid repeated disk access:
        # key: (series_id, sweep_id, trace_id) -> np.ndarray
        self._trace_cache: Dict[Tuple[int, int, int], np.ndarray] = {}

        # Cache sampling info per series (XStart/XInterval/DataPoints)
        self._series_info: Dict[int, Dict[str, Any]] = {}

    def __getitem__(self, *args):
        index = args[0]
        assert len(index) == 4
        group_id, series_id, sweep_id, trace_id = index
        if group_id != 0:
            raise IndexError("Only group_id=0 is supported for ABF composed bundles.")

        key = (series_id, sweep_id, trace_id)
        if key in self._trace_cache:
            return self._trace_cache[key]

        abf = self._get_abf(series_id)
        sweep_numbers = list(abf.sweepList)
        if sweep_id < 0 or sweep_id >= len(sweep_numbers):
            raise IndexError(f"sweep_id {sweep_id} out of range for series_id {series_id}.")
        abf.setSweep(sweep_numbers[sweep_id])

        if trace_id == 0:
            y = np.array(abf.sweepY, dtype=float)
            y_unit = getattr(abf, "sweepUnitsY", "")
            out = _scale_to_si(y, y_unit, target="V")
        elif trace_id == 1:
            c = np.array(abf.sweepC, dtype=float)
            c_unit = getattr(abf, "sweepUnitsC", "")
            out = _scale_to_si(c, c_unit, target="A")
        else:
            raise IndexError("ABF composed bundles expose exactly 2 traces: trace_id 0 (V) and 1 (A).")

        self._trace_cache[key] = out
        return out

    def _get_abf(self, series_id: int):
        try:
            series_node = self.bundle.pul[0][series_id]
            abf_path = getattr(series_node, "_abf_path")
        except Exception as e:
            raise RuntimeError(f"Could not resolve ABF path for series_id {series_id}") from e

        if abf_path in self._abf_cache:
            return self._abf_cache[abf_path]

        if not os.path.exists(abf_path):
            raise FileNotFoundError(f"ABF file not found: {abf_path}")

        import pyabf  # local import so HEKA users don't need it installed
        abf = pyabf.ABF(abf_path)
        self._abf_cache[abf_path] = abf

        # Build sweep/trace nodes now that we know sweep count and sampling
        self._populate_series_tree(series_id, abf)

        return abf

    def _populate_series_tree(self, series_id: int, abf: Any) -> None:
        if series_id in self._series_info:
            return

        series_node = self.bundle.pul[0][series_id]

        sweep_numbers = list(abf.sweepList)
        series_node.NumberSweeps = len(sweep_numbers)
        series_node.fields["NumberSweeps"] = len(sweep_numbers)

        dt = float(getattr(abf, "dataSecPerPoint", np.nan))
        if not np.isfinite(dt) or dt <= 0:
            # fallback: infer from sweepX
            abf.setSweep(sweep_numbers[0])
            x = np.array(abf.sweepX, dtype=float)
            if len(x) < 2:
                raise RuntimeError("Cannot infer XInterval from ABF sweepX.")
            dt = float(x[1] - x[0])

        # Two trace nodes per sweep
        for sweep_id, sweep_num in enumerate(sweep_numbers):
            sweep_node = _Node(Label=f"Sweep{sweep_id + 1}", SweepCount=sweep_id, Time=0.0)
            series_node.children.append(sweep_node)

            # Trace0: voltage
            tr0 = _Node(
                Label="Trace1",
                YUnit="V",
                XUnit="s",
                XStart=0.0,
                XInterval=dt,
            )
            # Trace1: command
            tr1 = _Node(
                Label="Trace2",
                YUnit="A",
                XUnit="s",
                XStart=0.0,
                XInterval=dt,
            )
            sweep_node.children.extend([tr0, tr1])

        self._series_info[series_id] = {"dt": dt, "sweeps": len(sweep_numbers)}


# ----------------------------
# Public Bundle factory
# ----------------------------

class Bundle:
    """
    Drop-in-ish replacement for heka_reader.Bundle:
      - If given a .dat path -> uses heka_reader.Bundle
      - Otherwise treats input as ABF root name and composes an ABF-backed bundle
    """
    def __init__(self, file_name: str, *, data_folder: Optional[str] = None, metadata_file: Optional[str] = None):
        self._impl = None

        if str(file_name).lower().endswith(".dat"):
            import heka_reader
            self._impl = heka_reader.Bundle(file_name)
            self.file_name = self._impl.file_name
        else:
            if data_folder is None or metadata_file is None:
                raise ValueError("ABF mode requires data_folder and metadata_file.")
            self._impl = _AbfComposedBundle(root_name=str(file_name), data_folder=data_folder, metadata_file=metadata_file)
            self.file_name = self._impl.file_name

        # forward extra attributes (like series_map) if present
        self.series_map = getattr(self._impl, "series_map", {})

    @property
    def pul(self):
        return self._impl.pul

    @property
    def data(self):
        return self._impl.data