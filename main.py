try:
    import config
except ImportError:
    print(
        "\nERROR: 'config.py' not found.\n"
        "Please create a local 'config.py' by copying 'config_template.py' and "
        "adjusting the paths for your system.\n"
    )
    raise SystemExit(1)
import os
import sys
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import find_peaks
import heka_reader
import git_save as myGit
from collections import defaultdict
import json
from statistics import median

from calculate_filtered_trace_and_derivatives import calculate_filtered_trace_and_derivatives


# --- parameters ---
A_to_pA = 1e12
V_to_mV = 1e3

def get_first_phase_peak(time, voltage_filt, d1_filt, peak_time, phase_plot_t1, phase_plot_t2,
                         phase_plot_range_start, phase_plot_range_end):
    """Find the first local peak in the phase plot within the requested voltage range.

    The search is restricted to phase_plot_t1 ms before and phase_plot_t2 ms after the AP peak.
    The phase plot is d1_filt vs. voltage_filt.
    """
    if peak_time is None:
        return None

    plot_start = peak_time - phase_plot_t1 / 1000.0
    plot_end = peak_time + phase_plot_t2 / 1000.0

    time_mask = (time >= plot_start) & (time <= plot_end)
    voltage_range_mask = (
        (voltage_filt >= phase_plot_range_start) &
        (voltage_filt <= phase_plot_range_end)
    )
    search_mask = time_mask & voltage_range_mask

    search_indices = np.where(search_mask)[0]
    if len(search_indices) < 3:
        return None

    search_d1 = d1_filt[search_indices]
    peaks, _ = find_peaks(search_d1)

    if len(peaks) == 0:
        return None

    phase_peak_idx = search_indices[peaks[0]]
    return {
        "time": float(time[phase_peak_idx]),
        "voltage_filt": float(voltage_filt[phase_peak_idx]),
        "d1_filt": float(d1_filt[phase_peak_idx]),
    }

def _points_inside_xlim(t_points, y_points, x_min, x_max):
    """Return only points whose x-coordinate lies inside the intended x-axis range."""
    t_points = np.asarray(t_points)
    y_points = np.asarray(y_points)
    inside = (t_points >= x_min) & (t_points <= x_max)
    return t_points[inside], y_points[inside]

def valid_point_pairs(t_values, y_values):
    """Return only point pairs where both time and y value are present."""
    return [
        (t, y)
        for t, y in zip(t_values, y_values)
        if t is not None and y is not None
    ]

def plot_analysis_points_on_phase_page(ax_voltage, ax_d1, ax_d2, ax_phase, points, phase_peak=None):
    """Overlay AP analysis points using the same visual language as the browser."""
    if not points:
        return

    voltage_symbols = {
        'threshold': ('o', 'red', 'Threshold'),
        'threshold_2nd': ('o', 'brown', 'Threshold 2nd'),
        'half_duration_start': ('s', 'blue', 'Half duration start'),
        'half_duration_end': ('s', 'green', 'Half duration end'),
        'peak': ('^', 'gold', 'Peak'),
        'ahp': ('D', 'purple', 'AHP'),
        'dvdt_max': ('P', 'cyan', 'Max dV/dt'),
    }

    d1_symbols = {
        'threshold_v1': ('o', 'red', 'Threshold'),
    }

    d2_symbols = {
        'd2_threshold': ('o', 'brown', 'Threshold 2nd'),
        'd2_peak1': ('^', 'purple', 'Peak1'),
        'd2_peak2': ('D', 'gold', 'Peak2'),
    }

    voltage_x_min, voltage_x_max = ax_voltage.get_xlim()
    d1_x_min, d1_x_max = ax_d1.get_xlim()
    d2_x_min, d2_x_max = ax_d2.get_xlim()
    phase_x_min, phase_x_max = ax_phase.get_xlim()

    for point_type, (marker, color, label) in voltage_symbols.items():
        if points.get(point_type):
            t_points, v_points = zip(*points[point_type])
            t_points, v_points = _points_inside_xlim(t_points, V_to_mV * np.asarray(v_points),
                                                     voltage_x_min, voltage_x_max)
            if len(t_points) > 0:
                ax_voltage.scatter(t_points, v_points,
                                   marker=marker, color=color, s=25, label=label, zorder=5)

    for point_type, (marker, color, label) in d1_symbols.items():
        if points.get(point_type):
            t_points, y_points = zip(*points[point_type])
            t_points, y_points = _points_inside_xlim(t_points, y_points, d1_x_min, d1_x_max)
            if len(t_points) > 0:
                ax_d1.scatter(t_points, y_points,
                              marker=marker, color=color, s=25, label=label, zorder=5)

    for point_type, (marker, color, label) in d2_symbols.items():
        if points.get(point_type):
            t_points, y_points = zip(*points[point_type])
            t_points, y_points = _points_inside_xlim(t_points, y_points, d2_x_min, d2_x_max)
            if len(t_points) > 0:
                ax_d2.scatter(t_points, y_points,
                              marker=marker, color=color, s=25, label=label, zorder=5)

    if points.get('threshold'):
        _, threshold_v = points['threshold'][0]
        threshold_mv = V_to_mV * threshold_v
        if phase_x_min <= threshold_mv <= phase_x_max:
            ax_phase.axvline(threshold_mv, color='red', linestyle='--', alpha=0.4)

    if points.get('threshold_2nd'):
        _, threshold_2nd_v = points['threshold_2nd'][0]
        threshold_2nd_mv = V_to_mV * threshold_2nd_v
        if phase_x_min <= threshold_2nd_mv <= phase_x_max:
            ax_phase.axvline(threshold_2nd_mv, color='brown', linestyle='--', alpha=0.4)

    if phase_peak is not None:
        phase_peak_voltage_mv = V_to_mV * phase_peak["voltage_filt"]
        if phase_x_min <= phase_peak_voltage_mv <= phase_x_max:
            ax_phase.scatter(
                phase_peak_voltage_mv,
                phase_peak["d1_filt"],
                marker='*',
                color='magenta',
                s=80,
                label='Phase peak AIS',
                zorder=6,
            )

def plot_phase_analysis_column(axs_column, phase_data, title_prefix):
    """Plot voltage, d1, d2 and phase plot for one selected first AP."""
    ax_voltage, ax_d1, ax_d2, ax_phase = axs_column

    if phase_data is None:
        for ax in axs_column:
            ax.set_title(f"{title_prefix}: SKIPPED")
            ax.grid(True)
        return

    time = phase_data["time"]
    voltage = phase_data["voltage"]
    voltage_filt = phase_data["voltage_filt"]
    d1 = phase_data["d1"]
    d1_filt = phase_data["d1_filt"]
    d2 = phase_data["d2"]
    d2_filt = phase_data["d2_filt"]
    plot_mask = phase_data["plot_mask"]
    points = phase_data["points"]
    phase_peak = phase_data["phase_peak"]

    ax_voltage.plot(time[plot_mask], V_to_mV * voltage[plot_mask], color='0.65', label='voltage')
    ax_voltage.plot(time[plot_mask], V_to_mV * voltage_filt[plot_mask], color='black', label='voltage_filt')
    ax_voltage.set_title(f"{title_prefix}: AP")
    ax_voltage.set_ylabel("Voltage (mV)")
    ax_voltage.set_xlabel("Time (s)")
    ax_voltage.grid(True)

    ax_d1.plot(time[plot_mask], d1[plot_mask], color='tab:blue', label='d1')
    ax_d1.plot(time[plot_mask], d1_filt[plot_mask], color='black', label='d1_filt')
    ax_d1.set_title(f"{title_prefix}: 1st derivative")
    ax_d1.set_ylabel("dV/dt (V/s)")
    ax_d1.set_xlabel("Time (s)")
    ax_d1.grid(True)

    ax_d2.plot(time[plot_mask], d2[plot_mask], color='tab:blue', label='d2')
    ax_d2.plot(time[plot_mask], d2_filt[plot_mask], color='black', label='d2_filt')
    ax_d2.set_title(f"{title_prefix}: 2nd derivative")
    ax_d2.set_ylabel("d²V/dt² (V/s²)")
    ax_d2.set_xlabel("Time (s)")
    ax_d2.grid(True)

    ax_phase.plot(V_to_mV * voltage_filt[plot_mask], d1_filt[plot_mask], color='black')
    ax_phase.set_title(f"{title_prefix}: phase plot")
    ax_phase.set_ylabel("dV/dt (V/s)")
    ax_phase.set_xlabel("voltage_filt (mV)")
    ax_phase.grid(True)


    # Set the intended x-axis limits explicitly.
    # Time-based plots: use the AP zoom window.
    # Phase plot: use the voltage range of the displayed phase trace.
    if np.any(plot_mask):
        time_x_min = float(time[plot_mask][0])
        time_x_max = float(time[plot_mask][-1])
        ax_voltage.set_xlim(time_x_min, time_x_max)
        ax_d1.set_xlim(time_x_min, time_x_max)
        ax_d2.set_xlim(time_x_min, time_x_max)

        #phase_x_values = V_to_mV * voltage_filt[plot_mask]
        #ax_phase.set_xlim(float(np.nanmin(phase_x_values)), float(np.nanmax(phase_x_values)))

    plot_analysis_points_on_phase_page(ax_voltage, ax_d1, ax_d2, ax_phase, points, phase_peak)

    for ax in axs_column:
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(fontsize=6, loc='best')

def prepare_phase_analysis_data(time, voltage, smooth_window, points, phase_plot_t1, phase_plot_t2,
                                phase_plot_range_start, phase_plot_range_end):
    """Prepare all arrays and detected phase peak for one first AP."""
    if not points or not points.get("peak"):
        return None

    peak_time = points["peak"][0][0]

    voltage_filt, d1, d1_filt, d2, d2_filt = calculate_filtered_trace_and_derivatives(
        time, voltage, smooth_window
    )

    plot_start = peak_time - phase_plot_t1 / 1000.0
    plot_end = peak_time + phase_plot_t2 / 1000.0
    plot_mask = (time >= plot_start) & (time <= plot_end)

    phase_peak = get_first_phase_peak(
        time,
        voltage_filt,
        d1_filt,
        peak_time,
        phase_plot_t1,
        phase_plot_t2,
        phase_plot_range_start,
        phase_plot_range_end,
    )

    return {
        "time": time,
        "voltage": voltage,
        "voltage_filt": voltage_filt,
        "d1": d1,
        "d1_filt": d1_filt,
        "d2": d2,
        "d2_filt": d2_filt,
        "plot_mask": plot_mask,
        "points": points,
        "phase_peak": phase_peak,
    }



def ap_analysis(time, voltage, v_threshold, dvdt_threshold, smooth_window, fraction_of_max_of_2nd_derivative,
                window_for_searching_threshold, window_for_searching_ahp,
                minimal_ap_interval, minimal_ap_duration, maximal_ap_duration,
                maximal_relative_amplitude_decline):

    sg_polyorder = 3
    dt = time[1] - time[0]

    voltage_filt, d1, d1_filt, d2, d2_filt = \
        calculate_filtered_trace_and_derivatives(time, voltage, smooth_window)

    # Find threshold crossings where voltage crosses v_threshold
    threshold_idx = np.where((voltage[:-1] < v_threshold) & (voltage[1:] >= v_threshold))[0]

    # Initialize results
    th_v = []
    th_v1 = []
    th_t = []
    th_v_2nd = []
    th_t_2nd = []
    th_d2_2nd = []
    p_v = []
    p_t = []
    ahp_v = []
    ahp_t = []
    hd_start_v = []
    hd_start_t = []
    hd_end_v = []
    hd_end_t = []
    maxdvdt_v = []
    maxdvdt_d1 = []
    maxdvdt_t = []
    d2_peak1_v = []
    d2_peak1_t = []
    d2_peak2_v = []
    d2_peak2_t = []

    for idx in threshold_idx:
        # Get window indices
        start_idx = max(0, idx - int(window_for_searching_threshold / dt))
        end_idx = min(len(voltage) - 1, idx + int(window_for_searching_threshold / dt))

        # Find precise threshold where dV/dt crosses threshold
        window_d1 = d1_filt[start_idx:end_idx]
        th_idx_candidates = np.where(window_d1 >= dvdt_threshold)[0]
        if len(th_idx_candidates) == 0:
            # Handle case where no threshold crossing is found
            continue  # Skip this iteration of the loop
        # Find threshold based on dV/dt
        th_idx = start_idx + th_idx_candidates[0]

        # Find AP peak first.
        # The 2nd-derivative peak search should only use the AP rising phase,
        # i.e. it must not extend beyond the voltage peak.
        peak_search_end_idx = min(len(voltage), th_idx + int(maximal_ap_duration / dt))
        if peak_search_end_idx <= th_idx:
            print("Warning: Invalid AP peak search window. Skipped this AP.")
            continue

        peak_idx = th_idx + np.argmax(voltage[th_idx:peak_search_end_idx])

        # Restrict 2nd-derivative analysis to the rising phase only:
        # from the threshold-search start up to and including the AP peak.
        d2_end_idx = peak_idx + 1
        if d2_end_idx - start_idx < 3:
            print("Warning: 2nd-derivative rising-phase window too short. Skipped this AP.")
            continue

        window_d2 = d2_filt[start_idx:d2_end_idx]

        # Find threshold based on 2nd derivative crossing
        # fraction_of_max_of_2nd_derivative of the rising-phase d2 maximum.
        d2_max_idx = np.argmax(window_d2)
        d2_max = window_d2[d2_max_idx]
        d2_threshold = fraction_of_max_of_2nd_derivative * d2_max

        # Find first crossing of the 2nd-derivative threshold.
        crossing_indices = []
        for i in range(1, len(window_d2)):
            if window_d2[i - 1] < d2_threshold and window_d2[i] >= d2_threshold:
                crossing_indices.append(i)

        if len(crossing_indices) > 0:
            th_idx_2nd_window = crossing_indices[0]
            th_idx_2nd = start_idx + th_idx_2nd_window
        else:
            print("Warning: Could not find 2nd derivative threshold crossing. Using d2 maximum instead.")
            th_idx_2nd_window = d2_max_idx
            th_idx_2nd = start_idx + th_idx_2nd_window

        # Find candidate peaks in the 2nd derivative.
        # Because window_d2 ends at peak_idx + 1, peak2 cannot be later than the AP peak.
        min_peak_height = np.max(window_d2) * 0.05  # 5% of max value
        peaks, _ = find_peaks(window_d2, height=min_peak_height)

        # Only peaks at/after the 2nd-derivative threshold are valid.
        # Peaks before the 2nd-derivative threshold are discarded.
        valid_peaks = np.array([peak for peak in peaks if peak >= th_idx_2nd_window])

        d2_peak1_idx = None
        d2_peak2_idx = None
        d2_peak1 = None
        d2_peak2 = None

        if len(valid_peaks) >= 2:
            # Use the first two valid peaks in temporal order.
            valid_peaks = np.sort(valid_peaks)
            d2_peak1_idx = valid_peaks[0]
            d2_peak2_idx = valid_peaks[1]
            d2_peak1 = window_d2[d2_peak1_idx]
            d2_peak2 = window_d2[d2_peak2_idx]

        elif len(valid_peaks) == 1:
            # If only one valid peak exists, it is peak2.
            # peak1 is intentionally left unset.
            d2_peak2_idx = valid_peaks[0]
            d2_peak2 = window_d2[d2_peak2_idx]

        else:
            # No valid peak after the 2nd-derivative threshold.
            # Leave both peaks unset.
            print("Warning: No valid 2nd derivative peak found between 2nd-derivative threshold and AP peak.")

        # Find AP peak
        peak_idx = th_idx + np.argmax(voltage[th_idx:th_idx + int(maximal_ap_duration / dt)])

        # Find AHP
        ahp_window = voltage[peak_idx:peak_idx + int(window_for_searching_ahp / dt)]
        ahp_idx = peak_idx + np.argmin(ahp_window)

        # Max dV/dt
        dvdt_idx = th_idx + np.argmax(d1_filt[th_idx:peak_idx])

        # Calculate half-duration points
        half_amplitude = (voltage[peak_idx] - voltage[th_idx]) / 2 + voltage[th_idx]

        hd_start_time = None
        hd_end_time = None

        # Find first crossing before peak
        for i in range(peak_idx, th_idx, -1):
            if voltage[i] <= half_amplitude:
                # --- Linear interpolation between i and i+1 ---
                v1, v2 = voltage[i], voltage[i + 1]
                t1, t2 = time[i], time[i + 1]
                frac = (half_amplitude - v1) / (v2 - v1) if v2 != v1 else 0.0
                hd_start_time = t1 + frac * (t2 - t1)
                break

        # Find first crossing after peak
        for i in range(peak_idx, ahp_idx):
            if voltage[i] <= half_amplitude:
                v1, v2 = voltage[i - 1], voltage[i]
                t1, t2 = time[i - 1], time[i]
                frac = (half_amplitude - v1) / (v2 - v1) if v2 != v1 else 0.0
                hd_end_time = t1 + frac * (t2 - t1)
                break

        # Skip this AP if we couldn't find valid half-duration points
        if hd_start_time is None or hd_end_time is None:
            print(f"Warning: Could not find valid half-duration points. Skipped this AP.")
            continue

        # For validation of the AP before appending
        half_duration = hd_end_time - hd_start_time
        ap_amplitude = voltage[peak_idx] - voltage[th_idx]

        # Calculate relative amplitude (relative to first AP)
        relative_amplitude = 1.0  # Default for first AP
        if len(p_v) > 0:  # If we already have APs
            first_ap_amplitude = p_v[0] - th_v[0]
            relative_amplitude = ap_amplitude / first_ap_amplitude

        # Check if this is first AP or if interval from previous AP is sufficient
        is_valid_interval = len(p_t) == 0 or (time[peak_idx] - p_t[-1] > minimal_ap_interval)

        if (is_valid_interval and
                minimal_ap_duration <= half_duration <= maximal_ap_duration and
                relative_amplitude >= maximal_relative_amplitude_decline):
            th_v.append(voltage[th_idx])
            th_v1.append(d1_filt[th_idx])
            th_t.append(time[th_idx])
            th_v_2nd.append(voltage[th_idx_2nd])
            th_t_2nd.append(time[th_idx_2nd])
            th_d2_2nd.append(d2_filt[th_idx_2nd])
            p_v.append(voltage[peak_idx])
            p_t.append(time[peak_idx])
            ahp_v.append(voltage[ahp_idx])
            ahp_t.append(time[ahp_idx])
            maxdvdt_v.append(voltage[dvdt_idx])
            maxdvdt_d1.append(d1_filt[dvdt_idx])
            maxdvdt_t.append(time[dvdt_idx])
            hd_start_v.append(half_amplitude)
            hd_start_t.append(hd_start_time)
            hd_end_v.append(half_amplitude)
            hd_end_t.append(hd_end_time)

            if d2_peak1_idx is not None:
                d2_peak1_v.append(d2_peak1)
                d2_peak1_t.append(time[start_idx + d2_peak1_idx])
            else:
                d2_peak1_v.append(None)
                d2_peak1_t.append(None)

            if d2_peak2_idx is not None:
                d2_peak2_v.append(d2_peak2)
                d2_peak2_t.append(time[start_idx + d2_peak2_idx])
            else:
                d2_peak2_v.append(None)
                d2_peak2_t.append(None)

    ap_number = len(p_v)

    return ap_number, th_v, th_v1, th_t, th_v_2nd, th_t_2nd, th_d2_2nd, hd_start_t, hd_start_v, hd_end_t, hd_end_v, p_v, p_t, ahp_v, ahp_t, maxdvdt_v, maxdvdt_d1, maxdvdt_t, d2_peak1_v, d2_peak1_t, d2_peak2_v, d2_peak2_t


def CC_eval():

    # --- Create Output Folders ---
    timestamp = datetime.now().strftime("%Y-%m-%d__%H-%M-%S")
    output_folder = os.path.join(config.ROOT_FOLDER, f"output_{config.MY_INITIAL}_{timestamp}")
    os.makedirs(output_folder, exist_ok=True)

    output_folder_results = os.path.join(output_folder, "results")
    os.makedirs(output_folder_results, exist_ok=True)

    output_folder_traces = os.path.join(output_folder, "traces")
    os.makedirs(output_folder_traces, exist_ok=True)

    output_folder_used_data_and_code = os.path.join(output_folder, "used_data_and_code")
    os.makedirs(output_folder_used_data_and_code, exist_ok=True)

    # --- Load Metadata ---
    metadata_df = pd.read_excel(config.METADATA_FILE)
    # save used data
    metadata_df.to_excel(os.path.join(output_folder_used_data_and_code, "my_data.xlsx"), index=False)

    # initialize analysis points structure
    analysis_points = defaultdict(
        lambda: defaultdict(
            lambda: defaultdict(
                lambda: defaultdict(
                    lambda: defaultdict(list)
                )
            )
        )
    )

    # === GIT SAVE ===
    # Provide the current script path (only works in .py, not notebooks)
    script_path = __file__ if '__file__' in globals() else None
    myGit.save_git_info(output_folder_used_data_and_code, script_path)

    # --- Results Storage ---
    results = []


    # --- Process Each Cell ---
    for cell_count, row in metadata_df.iterrows():
        print(f"Processing cell {cell_count + 1}: {row['file_name']}")

        file_name = row['file_name']

        # Check if series values are provided and valid (not NaN, not empty, and numeric)
        def is_valid_series(value):
            if pd.isna(value):
                return False
            try:
                # Try to convert to int - if it fails, it's not a valid series number
                int_val = int(float(value))  # float() first in case it's a string like "1.0"
                return int_val > 0  # Series numbers should be positive
            except (ValueError, TypeError):
                return False

        rmp_series_available = is_valid_series(row['rmp_series'])
        rin_series_available = is_valid_series(row['rin_series'])
        ap_rheo_series_available = is_valid_series(row['ap_rheo_series'])
        ap_max_series_available = is_valid_series(row['ap_max_series'])
        ap_broadening_series_available = is_valid_series(row['ap_broadening_series'])

        # Only convert to int and subtract 1 if value is available
        rmp_series = int(float(row['rmp_series'])) - 1 if rmp_series_available else None
        rin_series = int(float(row['rin_series'])) - 1 if rin_series_available else None
        ap_rheo_series = int(float(row['ap_rheo_series'])) - 1 if ap_rheo_series_available else None
        ap_max_series = int(float(row['ap_max_series'])) - 1 if ap_max_series_available else None
        ap_broadening_series = int(float(row['ap_broadening_series'])) - 1 if ap_broadening_series_available else None

        # Debug print to see what's happening
        #print(f"Series availability: RMP={rmp_series_available}, Rin={rin_series_available}, AP_rheo={ap_rheo_series_available}, AP_max={ap_max_series_available}")
        #print(f"Raw values: RMP={row['rmp_series']}, Rin={row['rin_series']}, AP_rheo={row['ap_rheo_series']}, AP_max={row['ap_max_series']}")

        # ... existing parameter loading code ...
        rin_series_from_current = row['rin_series_from_current']
        rin_series_to_current = row['rin_series_to_current']

        window1_rin_start = row['window1_rin_start']
        window1_rin_end = row['window1_rin_end']
        window2_rin_start = row['window2_rin_start']
        window2_rin_end = row['window2_rin_end']

        window1_ap_rheo_start = row['window1_ap_rheo_start']
        window1_ap_rheo_end = row['window1_ap_rheo_end']
        window2_ap_rheo_start = row['window2_ap_rheo_start']
        window2_ap_rheo_end = row['window2_ap_rheo_end']

        window1_ap_max_start = row['window1_ap_max_start']
        window1_ap_max_end = row['window1_ap_max_end']
        window2_ap_max_start = row['window2_ap_max_start']
        window2_ap_max_end = row['window2_ap_max_end']

        window1_ap_broadening_start = row['window1_ap_max_start']
        window1_ap_broadening_end = row['window1_ap_max_end']

        v_threshold = row['v_threshold'] / V_to_mV  #is provided in mV
        dvdt_threshold = row['dvdt_threshold']
        smooth_window = row['smooth_window']
        fraction_of_max_of_2nd_derivative = row['fraction_of_max_of_2nd_derivative']
        window_for_searching_threshold = row['window_for_searching_threshold']
        window_for_searching_ahp = row['window_for_searching_ahp']
        minimal_ap_interval = row['minimal_ap_interval']
        minimal_ap_duration = row['minimal_ap_duration']
        maximal_ap_duration = row['maximal_ap_duration']
        maximal_relative_amplitude_decline = row['maximal_relative_amplitude_decline']

        phase_plot_t1 = row['phase_plot_t1']
        phase_plot_t2 = row['phase_plot_t2']
        phase_plot_range_start = row['phase_plot_range_start'] / V_to_mV
        phase_plot_range_end = row['phase_plot_range_end'] / V_to_mV

        dat_path = os.path.join(config.EXTERNAL_DATA_FOLDER, file_name)
        try:
            bundle = heka_reader.Bundle(dat_path)
        except Exception as e:
            print(f"Error reading {file_name}: {e}")
            continue

        group_id = 0
        fig, axs = plt.subplots(5, 2, figsize=(8, 12))
        axs = axs.flatten()

        # Initialize all result variables to None (will appear as blank in Excel)
        rmp_mean = None
        rmp_min = None
        rin_fit = None
        rheobase = None
        max_ap_number = None

        ap_rheo_baseline_voltage = None
        ap_rheo_first_ap_delay = None

        ap_rheo_half_duration_1st = None
        ap_rheo_threshold_1st = None
        ap_rheo_threshold_2nd_1st = None
        ap_rheo_amplitude_1st = None

        ap_rheo_half_duration_av = None
        ap_rheo_threshold_av = None
        ap_rheo_threshold_2nd_av = None
        ap_rheo_amplitude_av = None

        ap_rheo_peak1_peak2_interval_1st = None
        ap_rheo_maxdvdt_1st = None

        ap_rheo_phase_peak_d1_filt = None
        ap_rheo_phase_peak_voltage_filt = None
        ap_rheo_phase_data = None

        ap_max_baseline_voltage = None

        ap_max_half_duration_1st = None
        ap_max_threshold_1st = None
        ap_max_threshold_2nd_1st = None
        ap_max_amplitude_1st = None

        ap_max_half_duration_av = None
        ap_max_threshold_av = None
        ap_max_threshold_2nd_av = None
        ap_max_amplitude_av = None

        ap_max_instantaneous_freq_1_2 = None
        ap_max_instantaneous_freq_last = None
        ap_max_freq_adaptation = None
        ap_max_average_freq_all = None
        ap_max_peak1_peak2_interval_1st = None
        ap_max_peak1_peak2_interval_av = None
        ap_max_peak1_peak2_interval_median = None

        ap_max_phase_peak_d1_filt = None
        ap_max_phase_peak_voltage_filt = None
        ap_max_phase_data = None

        ap_max_list_current_steps = []
        ap_max_list_ap_numbers = []
        ap_max_list_half_duration_1st = []
        ap_max_list_threshold_1st = []
        ap_max_list_amplitude_1st = []
        ap_max_list_voltage_baseline = []
        ap_max_list_average_frequency = []
        ap_max_list_instantaneous_freq_1_2 = []
        ap_max_list_instantaneous_last = []
        ap_max_list_freq_adaptation = []

        ap_broadening_list_voltage_baseline = []
        ap_broadening_list_half_duration_1st = []
        ap_broadening_list_threshold_1st = []
        ap_broadening_list_threshold_2nd_1st = []
        ap_broadening_list_amplitude_1st = []


        # ==========================================================================================
        # --- RMP ---
        # ==========================================================================================

        if rmp_series_available:
            series_id = rmp_series

            n_sweeps = bundle.pul[group_id][series_id].NumberSweeps
            voltage_trace = bundle.data[group_id, series_id, 0, 0]

            n_points = len(voltage_trace)
            sampling_interval = bundle.pul[group_id][series_id][0][0].XInterval
            time = np.arange(n_points) * sampling_interval

            rmp_mean = voltage_trace.mean()
            rmp_min = voltage_trace.min()

            # Plot RMP trace with custom axis
            axs[0].plot(time, V_to_mV * voltage_trace)
            axs[0].set_title("RMP trace")
            axs[0].set_ylabel("Voltage (mV)")
            axs[0].set_xlabel("time (s)")
            axs[0].grid(True)

            if n_sweeps != 1:
                print(f"Warning: RMP series in {file_name} has {n_sweeps} sweeps. Only number 1 was used.")
        else:
            # Skip RMP analysis - create empty plot
            axs[0].set_title("SKIPPED")
            axs[0].set_ylabel("Voltage (mV)")
            axs[0].set_xlabel("time (s)")
            axs[0].grid(True)
            print(f"        Skipping RMP analysis")

        # ==========================================================================================
        # --- Rin ---
        # ==========================================================================================

        if rin_series_available:
            series_id = rin_series
            n_sweeps = bundle.pul[group_id][series_id].NumberSweeps

            delta_vs = []
            delta_is = []

            # ... existing Rin analysis code ...
            # get time base
            current = bundle.data[group_id, series_id, 0, 1]
            n_points = len(current)
            sampling_interval = bundle.pul[group_id][series_id][0][0].XInterval
            time = np.arange(n_points) * sampling_interval

            # Convert windows to indices
            idx1 = (time >= window1_rin_start) & (time <= window1_rin_end)
            idx2 = (time >= window2_rin_start) & (time <= window2_rin_end)

            # Plot superposition of all current traces
            axs[2].set_title("Rin Current Traces Superposition")
            axs[2].set_ylabel("Current (pA)")
            axs[2].set_xlabel("Time (s)")

            # Plot superposition of all voltage traces
            axs[3].set_title("Rin Voltage Traces Superposition")
            axs[3].set_ylabel("Voltage (mV)")
            axs[3].set_xlabel("Time (s)")

            # Determine sweep range
            for sweep_id in range(n_sweeps):
                voltage = bundle.data[group_id, series_id, sweep_id, 0]
                current = A_to_pA * bundle.data[group_id, series_id, sweep_id, 1]

                dv = voltage[idx2].mean() - voltage[idx1].mean()
                di = current[idx2].mean() - current[idx1].mean()

                if not (rin_series_from_current <= di <= rin_series_to_current):
                    continue

                # Plot current trace
                axs[2].plot(time, current, alpha=0.5)

                # Plot voltage trace
                axs[3].plot(time, V_to_mV * voltage, alpha=0.5)

                delta_vs.append(dv)
                delta_is.append(di)

            # Add vertical lines for analysis windows in both current and voltage plots
            for ax in [axs[2], axs[3]]:
                ax.axvline(x=window1_rin_start, color='r', linestyle='--', alpha=0.3)
                ax.axvline(x=window1_rin_end, color='r', linestyle='--', alpha=0.3)
                ax.axvline(x=window2_rin_start, color='r', linestyle='--', alpha=0.3)
                ax.axvline(x=window2_rin_end, color='r', linestyle='--', alpha=0.3)
                ax.grid(True)

            delta_vs = V_to_mV * np.array(delta_vs)
            delta_is = np.array(delta_is)

            # Linear fit: ΔV = R * ΔI (forced through zero)
            slope = np.sum(delta_is * delta_vs) / np.sum(delta_is * delta_is)  # Least squares through origin
            rin_fit = 1000 * slope  # input resistance in MOhm (1e6 = 1000 * milli / pico)

            # Plot I-V relationship
            axs[1].scatter(delta_is, delta_vs, label="data")
            axs[1].plot(delta_is, slope * delta_is, color="red",
                        label=f"fit (R={rin_fit:.1f} MΩ)")
            axs[1].set_xlabel("ΔI (pA)")
            axs[1].set_ylabel("ΔV (mV)")
            axs[1].set_title("Input Resistance")
            axs[1].legend()
            axs[1].grid(True)
        else:
            # Skip Rin analysis - create empty plots
            axs[1].set_title("SKIPPED")
            axs[1].set_xlabel("ΔI (pA)")
            axs[1].set_ylabel("ΔV (mV)")
            axs[1].grid(True)

            axs[2].set_title("SKIPPED")
            axs[2].set_ylabel("Current (pA)")
            axs[2].set_xlabel("Time (s)")
            axs[2].grid(True)

            axs[3].set_title("SKIPPED")
            axs[3].set_ylabel("Voltage (mV)")
            axs[3].set_xlabel("Time (s)")
            axs[3].grid(True)
            print(f"        Skipping Rin analysis")

        # ==========================================================================================
        # --- AP rheo ---
        # ==========================================================================================

        if ap_rheo_series_available:
            series_id = ap_rheo_series
            n_sweeps = bundle.pul[group_id][series_id].NumberSweeps

            # get time base
            current = bundle.data[group_id, series_id, 0, 1]
            n_points = len(current)
            sampling_interval = bundle.pul[group_id][series_id][0][0].XInterval
            time = np.arange(n_points) * sampling_interval

            # Convert windows to indices
            idx1 = (time >= window1_ap_rheo_start) & (time <= window1_ap_rheo_end)
            idx2 = (time >= window2_ap_rheo_start) & (time <= window2_ap_rheo_end)

            # Plot superposition of all voltage traces
            axs[4].set_title("AP rheo Voltage Traces Superposition")
            axs[4].set_ylabel("Voltage (mV)")
            axs[4].set_xlabel("Time (s)")

            sweep_points = None  # Default value when no APs are detected

            for sweep_id in range(n_sweeps):
                voltage = bundle.data[group_id, series_id, sweep_id, 0]
                current = bundle.data[group_id, series_id, sweep_id, 1]

                # Plot voltage trace in superposition plot
                axs[4].plot(time, V_to_mV * voltage, alpha=0.5, label=f'Sweep {sweep_id + 1}')

                di = current[idx2].mean() - current[idx1].mean()

                ap_number, th_v, th_v1, th_t, th_v_2nd, th_t_2nd, th_d2_2nd, hd_start_t, hd_start_v, hd_end_t, hd_end_v, p_v, p_t, ahp_v, ahp_t, maxdvdt_v, maxdvdt_d1, maxdvdt_t, d2_peak1_v, d2_peak1_t, d2_peak2_v, d2_peak2_t = ap_analysis(
                    time, voltage, v_threshold, dvdt_threshold, smooth_window, fraction_of_max_of_2nd_derivative,
                    window_for_searching_threshold, window_for_searching_ahp,
                    minimal_ap_interval, minimal_ap_duration, maximal_ap_duration,
                    maximal_relative_amplitude_decline)

                #print(f"Sweep {sweep_id + 1}: ap_number = {ap_number}")
                if ap_number > 0:
                    sweep_points = {
                        'threshold': list(zip(th_t, th_v)),
                        'threshold_v1': list(zip(th_t, th_v1)),
                        'threshold_2nd': list(zip(th_t_2nd, th_v_2nd)),
                        'half_duration_start': list(zip(hd_start_t, hd_start_v)),
                        'half_duration_end': list(zip(hd_end_t, hd_end_v)),
                        'peak': list(zip(p_t, p_v)),
                        'ahp': list(zip(ahp_t, ahp_v)),
                        'dvdt_max': list(zip(maxdvdt_t, maxdvdt_v)),
                        'd2_peak1': valid_point_pairs(d2_peak1_t, d2_peak1_v),
                        'd2_peak2': valid_point_pairs(d2_peak2_t, d2_peak2_v),
                        'd2_threshold': list(zip(th_t_2nd, th_d2_2nd)),
                        'smooth_window': smooth_window
                    }
                    # Store the analysis points in the nested dictionary
                    analysis_points[file_name][group_id][series_id][sweep_id][0] = sweep_points

                if rheobase is None and ap_number > 0:
                    rheobase_voltage_trace = voltage
                    rheobase = di
                    # Calculate baseline voltage from window1 (before stimulus)
                    ap_rheo_baseline_voltage = voltage[idx1].mean()
                    # Calculate delay of first AP (time from stimulus start to first AP threshold)
                    # Assuming stimulus starts at window2_ap_rheo_start
                    ap_rheo_first_ap_delay = th_t[0]
                    
                    ap_rheo_half_duration_1st = hd_end_t[0] - hd_start_t[0]
                    ap_rheo_threshold_1st = V_to_mV * th_v[0]
                    ap_rheo_threshold_2nd_1st = V_to_mV * th_v_2nd[0]
                    ap_rheo_amplitude_1st = V_to_mV * p_v[0] - th_v[0]
                    if d2_peak1_t[0] is not None and d2_peak2_t[0] is not None:
                        ap_rheo_peak1_peak2_interval_1st = d2_peak2_t[0] - d2_peak1_t[0]
                    else:
                        ap_rheo_peak1_peak2_interval_1st = None
                    ap_rheo_maxdvdt_1st = maxdvdt_d1[0]

                    ap_rheo_phase_data = prepare_phase_analysis_data(
                        time,
                        voltage,
                        smooth_window,
                        sweep_points,
                        phase_plot_t1,
                        phase_plot_t2,
                        phase_plot_range_start,
                        phase_plot_range_end,
                    )
                    if ap_rheo_phase_data is not None and ap_rheo_phase_data["phase_peak"] is not None:
                        ap_rheo_phase_peak_d1_filt = ap_rheo_phase_data["phase_peak"]["d1_filt"]
                        ap_rheo_phase_peak_voltage_filt = ap_rheo_phase_data["phase_peak"]["voltage_filt"]

                    # Calculate average AP parameters
                    ap_rheo_half_duration_av = sum(hd_end_t[i] - hd_start_t[i] for i in range(ap_number)) / ap_number
                    ap_rheo_threshold_av = V_to_mV * sum(th_v) / len(th_v)
                    ap_rheo_threshold_2nd_av = V_to_mV * sum(th_v_2nd) / len(th_v_2nd)
                    ap_rheo_amplitude_av = V_to_mV * sum(p_v[i] - th_v[i] for i in range(len(th_v))) / len(th_v)

            # Formatting for AP rheo superposition plot
            axs[4].grid(True)

            # Setup rheobase trace plot
            axs[5].set_title("Rheobase Voltage Trace")
            axs[5].set_ylabel("Voltage (mV)")
            axs[5].set_xlabel("Time (s)")
            if rheobase is not None:
                axs[5].plot(time, V_to_mV * rheobase_voltage_trace, label=f'Rheobase: {rheobase:.1f} pA')
                axs[5].legend()
            axs[5].grid(True)
        else:
            # Skip AP rheo analysis - create empty plots
            axs[4].set_title("SKIPPED")
            axs[4].set_ylabel("Voltage (mV)")
            axs[4].set_xlabel("Time (s)")
            axs[4].grid(True)

            axs[5].set_title("SKIPPED")
            axs[5].set_ylabel("Voltage (mV)")
            axs[5].set_xlabel("Time (s)")
            axs[5].grid(True)
            print(f"        Skipping AP rheo analysis")

        # ==========================================================================================
        # --- AP max ---
        # ==========================================================================================

        if ap_max_series_available:
            series_id = ap_max_series
            n_sweeps = bundle.pul[group_id][series_id].NumberSweeps

            # get time base
            current = bundle.data[group_id, series_id, 0, 1]
            n_points = len(current)
            sampling_interval = bundle.pul[group_id][series_id][0][0].XInterval
            time = np.arange(n_points) * sampling_interval

            # Convert windows to indices
            idx1 = (time >= window1_ap_max_start) & (time <= window1_ap_max_end)
            idx2 = (time >= window2_ap_max_start) & (time <= window2_ap_max_end)

            # Plot superposition of all voltage traces
            axs[6].set_title("AP max Voltage Traces Superposition")
            axs[6].set_ylabel("Voltage (mV)")
            axs[6].set_xlabel("Time (s)")

            max_ap_number = 0
            max_ap_di = None
            max_ap_voltage_trace = None

            for sweep_id in range(n_sweeps):
                voltage = bundle.data[group_id, series_id, sweep_id, 0]
                current = bundle.data[group_id, series_id, sweep_id, 1]

                # Plot voltage trace in superposition plot
                axs[6].plot(time, V_to_mV * voltage, alpha=0.5, label=f'Sweep {sweep_id + 1}')

                di = current[idx2].mean() - current[idx1].mean()
                ap_max_list_current_steps.append(float(di))  # Convert numpy.float64 to Python float
                
                # Calculate baseline voltage for this sweep
                baseline_voltage = voltage[idx1].mean()
                ap_max_list_voltage_baseline.append(float(baseline_voltage))

                ap_number, th_v, th_v1, th_t, th_v_2nd, th_t_2nd, th_d2_2nd, hd_start_t, hd_start_v, hd_end_t, hd_end_v, p_v, p_t, ahp_v, ahp_t, maxdvdt_v, maxdvdt_d1, maxdvdt_t, d2_peak1_v, d2_peak1_t, d2_peak2_v, d2_peak2_t = ap_analysis(
                    time, voltage, v_threshold, dvdt_threshold, smooth_window, fraction_of_max_of_2nd_derivative,
                    window_for_searching_threshold, window_for_searching_ahp,
                    minimal_ap_interval, minimal_ap_duration, maximal_ap_duration,
                    maximal_relative_amplitude_decline)

                ap_max_list_ap_numbers.append(int(ap_number))  # Convert to Python int
                #print(f"Sweep {sweep_id + 1}: ap_number = {ap_number}")

                if ap_number > 0:
                    sweep_points = {
                        'threshold': list(zip(th_t, th_v)),
                        'threshold_v1': list(zip(th_t, th_v1)),
                        'threshold_2nd': list(zip(th_t_2nd, th_v_2nd)),
                        'half_duration_start': list(zip(hd_start_t, hd_start_v)),
                        'half_duration_end': list(zip(hd_end_t, hd_end_v)),
                        'peak': list(zip(p_t, p_v)),
                        'ahp': list(zip(ahp_t, ahp_v)),
                        'dvdt_max': list(zip(maxdvdt_t, maxdvdt_v)),
                        'd2_peak1': valid_point_pairs(d2_peak1_t, d2_peak1_v),
                        'd2_peak2': valid_point_pairs(d2_peak2_t, d2_peak2_v),
                        'd2_threshold': list(zip(th_t_2nd, th_d2_2nd)),
                        'smooth_window': smooth_window
                    }
                    # Store the analysis points in the nested dictionary
                    analysis_points[file_name][group_id][series_id][sweep_id][0] = sweep_points
                    # Log the structure of analysis points:
                    # print(f"Storing points: file_name={file_name}, group_id={group_id}, series_id={series_id}, sweep_id={sweep_id}, trace_id={0}")
                    # print(f"Points to store: {json.dumps(sweep_points, indent=2)}")

                    ap_max_list_half_duration_1st.append(float(hd_end_t[0] - hd_start_t[0]))
                    ap_max_list_threshold_1st.append(V_to_mV * float(th_v[0]))
                    ap_max_list_amplitude_1st.append(V_to_mV * float(p_v[0] - th_v[0]))
                else:
                    ap_max_list_half_duration_1st.append(None)
                    ap_max_list_threshold_1st.append(None)
                    ap_max_list_amplitude_1st.append(None)

                # Calculate various average frequencies for this sweep
                if ap_number >= 2:
                    instantaneous_freq_1_2 = 1.0 / (p_t[1] - p_t[0]) # Instantaneous frequency of first two APs (Hz)
                    instantaneous_freq_last = 1.0 / (p_t[-1] - p_t[-2]) # Instantaneous frequency of last two APs (Hz)
                    total_time = p_t[-1] - p_t[0]  # Time from first to last AP
                    average_frequency = (ap_number - 1) / total_time  # Number of intervals / total time
                    ap_max_list_average_frequency.append(float(average_frequency))
                    ap_max_list_instantaneous_freq_1_2.append(float(instantaneous_freq_1_2))
                    ap_max_list_instantaneous_last.append(float(instantaneous_freq_last))
                    ap_max_list_freq_adaptation.append(float(instantaneous_freq_last/instantaneous_freq_1_2))
                else:
                    ap_max_list_average_frequency.append(None)
                    ap_max_list_instantaneous_freq_1_2.append(None)
                    ap_max_list_instantaneous_last.append(None)
                    ap_max_list_freq_adaptation.append(None)

                # find trace with maximal number of APs
                # if the next trace again has this amount of APs, then use these values. I.e. the last trace with the most number of APs is used for saving the parameters.
                if ap_number >= max_ap_number and ap_number > 0:
                    max_ap_number = ap_number
                    max_ap_voltage_trace = voltage
                    max_ap_di = di
                    ap_max_baseline_voltage = V_to_mV * voltage[idx1].mean() # Calculate baseline voltage from window1 (before stimulus) for max AP trace
                    ap_max_half_duration_1st = hd_end_t[0] - hd_start_t[0]
                    ap_max_threshold_1st = V_to_mV * th_v[0]
                    ap_max_threshold_2nd_1st = V_to_mV * th_v_2nd[0]
                    ap_max_amplitude_1st = V_to_mV * p_v[0] - th_v[0]
                    if d2_peak1_t[0] is not None and d2_peak2_t[0] is not None:
                        ap_max_peak1_peak2_interval_1st = d2_peak2_t[0] - d2_peak1_t[0]
                    else:
                        ap_max_peak1_peak2_interval_1st = None

                    ap_max_phase_data = prepare_phase_analysis_data(
                        time,
                        voltage,
                        smooth_window,
                        sweep_points,
                        phase_plot_t1,
                        phase_plot_t2,
                        phase_plot_range_start,
                        phase_plot_range_end,
                    )
                    if ap_max_phase_data is not None and ap_max_phase_data["phase_peak"] is not None:
                        ap_max_phase_peak_d1_filt = ap_max_phase_data["phase_peak"]["d1_filt"]
                        ap_max_phase_peak_voltage_filt = ap_max_phase_data["phase_peak"]["voltage_filt"]
                    else:
                        ap_max_phase_peak_d1_filt = None
                        ap_max_phase_peak_voltage_filt = None

                    # Calculate average AP parameters
                    ap_max_half_duration_av = sum(hd_end_t[i] - hd_start_t[i] for i in range(ap_number)) / ap_number
                    ap_max_threshold_av = V_to_mV * sum(th_v) / ap_number
                    ap_max_threshold_2nd_av = V_to_mV * sum(th_v_2nd) / ap_number
                    ap_max_amplitude_av = V_to_mV * sum(p_v[i] - th_v[i] for i in range(ap_number)) / ap_number

                    valid_peak1_peak2_intervals = [
                        d2_peak2_t[i] - d2_peak1_t[i]
                        for i in range(ap_number)
                        if d2_peak1_t[i] is not None and d2_peak2_t[i] is not None
                    ]

                    if valid_peak1_peak2_intervals:
                        ap_max_peak1_peak2_interval_av = sum(valid_peak1_peak2_intervals) / len(valid_peak1_peak2_intervals)
                        ap_max_peak1_peak2_interval_median = median(valid_peak1_peak2_intervals)
                    else:
                        ap_max_peak1_peak2_interval_av = None
                        ap_max_peak1_peak2_interval_median = None

                    if ap_number >= 2:      #within the if for finding the trace with most APs we check if there are more then 2 APs
                        ap_max_average_freq_all = average_frequency # calculated above
                        ap_max_instantaneous_freq_1_2 = instantaneous_freq_1_2
                        ap_max_instantaneous_freq_last = instantaneous_freq_last
                        ap_max_freq_adaptation = instantaneous_freq_last/instantaneous_freq_1_2
                    else:
                        ap_max_average_freq_all = None
                        ap_max_instantaneous_freq_1_2 = None
                        ap_max_instantaneous_freq_last = None
                        ap_max_freq_adaptation = None

            axs[6].grid(True)

            # Setup max AP trace plot
            axs[7].set_title("Max AP Voltage Trace")
            axs[7].set_ylabel("Voltage (mV)")
            axs[7].set_xlabel("Time (s)")
            if max_ap_voltage_trace is not None:
                axs[7].plot(time, V_to_mV * max_ap_voltage_trace, label=f'Max APs: {max_ap_number} at {max_ap_di:.1f} pA')
                axs[7].legend()
            axs[7].grid(True)

            # Plot number of APs vs current
            axs[8].set_title("AP Frequency vs Current")
            axs[8].set_ylabel("Number of APs")
            axs[8].set_xlabel("Current (pA)")
            axs[8].plot(ap_max_list_current_steps, ap_max_list_ap_numbers, 'o-')
            axs[8].grid(True)

        else:
            # Skip AP max analysis - create empty plots
            axs[6].set_title("SKIPPED")
            axs[6].set_ylabel("Voltage (mV)")
            axs[6].set_xlabel("Time (s)")
            axs[6].grid(True)

            axs[7].set_title("SKIPPED")
            axs[7].set_ylabel("Voltage (mV)")
            axs[7].set_xlabel("Time (s)")
            axs[7].grid(True)

            axs[8].set_title("SKIPPED")
            axs[8].set_ylabel("Number of APs")
            axs[8].set_xlabel("Current (pA)")
            axs[8].grid(True)
            print(f"        Skipping AP max analysis")

        # ==========================================================================================
        # --- AP broadening ---
        # ==========================================================================================

        if ap_broadening_series_available:
            series_id = ap_broadening_series
            n_sweeps = bundle.pul[group_id][series_id].NumberSweeps

            # get time base
            current = bundle.data[group_id, series_id, 0, 1]
            n_points = len(current)
            sampling_interval = bundle.pul[group_id][series_id][0][0].XInterval
            time = np.arange(n_points) * sampling_interval

            # Convert windows to indices
            idx1 = (time >= window1_ap_broadening_start) & (time <= window1_ap_broadening_end)

            for sweep_id in range(n_sweeps):
                voltage = bundle.data[group_id, series_id, sweep_id, 0]
                #current = bundle.data[group_id, series_id, sweep_id, 1]

                baseline_voltage = voltage[idx1].mean()
                ap_broadening_list_voltage_baseline.append(float(baseline_voltage))  # Convert numpy.float64 to Python float

                ap_number, th_v, th_v1, th_t, th_v_2nd, th_t_2nd, th_d2_2nd, hd_start_t, hd_start_v, hd_end_t, hd_end_v, p_v, p_t, ahp_v, ahp_t, maxdvdt_v, maxdvdt_d1, maxdvdt_t, d2_peak1_v, d2_peak1_t, d2_peak2_v, d2_peak2_t = ap_analysis(
                    time, voltage, v_threshold, dvdt_threshold, smooth_window, fraction_of_max_of_2nd_derivative,
                    window_for_searching_threshold, window_for_searching_ahp,
                    minimal_ap_interval, minimal_ap_duration, maximal_ap_duration,
                    maximal_relative_amplitude_decline)

                # Store half duration of first AP (or None if no AP detected)
                if ap_number > 0:
                    # Store analysis points
                    sweep_points = {
                        'threshold': list(zip(th_t, th_v)),
                        'threshold_v1': list(zip(th_t, th_v1)),
                        'threshold_2nd': list(zip(th_t_2nd, th_v_2nd)),
                        'half_duration_start': list(zip(hd_start_t, hd_start_v)),
                        'half_duration_end': list(zip(hd_end_t, hd_end_v)),
                        'peak': list(zip(p_t, p_v)),
                        'ahp': list(zip(ahp_t, ahp_v)),
                        'dvdt_max': list(zip(maxdvdt_t, maxdvdt_v)),
                        'd2_peak1': valid_point_pairs(d2_peak1_t, d2_peak1_v),
                        'd2_peak2': valid_point_pairs(d2_peak2_t, d2_peak2_v),
                        'd2_threshold': list(zip(th_t_2nd, th_d2_2nd)),
                        'smooth_window': smooth_window
                    }
                    analysis_points[file_name][group_id][series_id][sweep_id][0] = sweep_points

                    ap_broadening_list_half_duration_1st.append(float(hd_end_t[0] - hd_start_t[0]))
                    ap_broadening_list_threshold_1st.append(V_to_mV * float(th_v[0]))
                    ap_broadening_list_threshold_2nd_1st.append(V_to_mV * float(th_v_2nd[0]))
                    ap_broadening_list_amplitude_1st.append(V_to_mV * float(p_v[0] - th_v[0]))
                else:
                    ap_broadening_list_half_duration_1st.append(None)
                    ap_broadening_list_threshold_1st.append(None)
                    ap_broadening_list_threshold_2nd_1st.append(None)
                    ap_broadening_list_amplitude_1st.append(None)

            # Plot half duration of first AP vs sweep number in axs[9]
            sweep_numbers = list(range(1, len(ap_broadening_list_half_duration_1st) + 1))
            valid_indices = [i for i, hd in enumerate(ap_broadening_list_half_duration_1st) if hd is not None]
            valid_sweep_numbers = [sweep_numbers[i] for i in valid_indices]
            valid_half_durations = [ap_broadening_list_half_duration_1st[i] for i in valid_indices]

            axs[9].set_title("AP Broadening (1st AP Half Duration)")
            axs[9].set_ylabel("Half Duration (s)")
            axs[9].set_xlabel("Sweep Number")
            if valid_half_durations:
                axs[9].plot(valid_sweep_numbers, valid_half_durations, 'o-')
            axs[9].grid(True)

        else:
            # Skip AP broadening analysis - create empty plot
            axs[9].set_title("SKIPPED")
            axs[9].set_ylabel("Half Duration (s)")
            axs[9].set_xlabel("Sweep Number")
            axs[9].grid(True)
            print(f"        Skipping AP broadening analysis")

        # ==========================================================================================
        # --- Store results within the loop ---
        # ==========================================================================================
        results.append({
            "cell_count": cell_count + 1,
            "file_name": file_name,
            "RMP_mean": rmp_mean,
            "RMP_min": rmp_min,
            "Rin_fit": rin_fit,
            "Rheobase": rheobase,
            "max_ap_number": max_ap_number,

            "ap_rheo_baseline_voltage": ap_rheo_baseline_voltage,
            "ap_rheo_first_ap_delay": ap_rheo_first_ap_delay,

            "ap_rheo_half_duration_1st": ap_rheo_half_duration_1st,
            "ap_rheo_threshold_1st": ap_rheo_threshold_1st,
            "ap_rheo_threshold_2nd_1st": ap_rheo_threshold_2nd_1st,
            "ap_rheo_amplitude_1st": ap_rheo_amplitude_1st,
            
            "ap_rheo_half_duration_av": ap_rheo_half_duration_av,
            "ap_rheo_threshold_av": ap_rheo_threshold_av,
            "ap_rheo_threshold_2nd_av": ap_rheo_threshold_2nd_av,
            "ap_rheo_amplitude_av": ap_rheo_amplitude_av,
            "ap_rheo_maxdvdt_1st": ap_rheo_maxdvdt_1st,
            "ap_rheo_peak1_peak2_interval_1st": ap_rheo_peak1_peak2_interval_1st,
            "ap_rheo_phase_peak_d1_filt": ap_rheo_phase_peak_d1_filt,
            "ap_rheo_phase_peak_voltage_filt": ap_rheo_phase_peak_voltage_filt,

            "ap_max_half_duration_1st": ap_max_half_duration_1st,
            "ap_max_threshold_1st": ap_max_threshold_1st,
            "ap_max_threshold_2nd_1st": ap_max_threshold_2nd_1st,
            "ap_max_amplitude_1st": ap_max_amplitude_1st,
            "ap_max_baseline_voltage": ap_max_baseline_voltage,

            "ap_max_half_duration_av": ap_max_half_duration_av,
            "ap_max_threshold_av": ap_max_threshold_av,
            "ap_max_threshold_2nd_av": ap_max_threshold_2nd_av,
            "ap_max_amplitude_av": ap_max_amplitude_av,
            "ap_max_peak1_peak2_interval_1st": ap_max_peak1_peak2_interval_1st,
            "ap_max_peak1_peak2_interval_av": ap_max_peak1_peak2_interval_av,
            "ap_max_peak1_peak2_interval_median": ap_max_peak1_peak2_interval_median,
            "ap_max_average_freq_all": ap_max_average_freq_all,
            "ap_max_instantaneous_freq_1_2": ap_max_instantaneous_freq_1_2,
            "ap_max_instantaneous_freq_last": ap_max_instantaneous_freq_last,
            "ap_max_freq_adaptation": ap_max_freq_adaptation,
            "ap_max_phase_peak_d1_filt": ap_max_phase_peak_d1_filt,
            "ap_max_phase_peak_voltage_filt": ap_max_phase_peak_voltage_filt,

            "ap_max_list_current_steps": ap_max_list_current_steps,
            "ap_max_list_ap_numbers": ap_max_list_ap_numbers,
            "ap_max_list_half_duration_1st": ap_max_list_half_duration_1st,
            "ap_max_list_threshold_1st": ap_max_list_threshold_1st,
            "ap_max_list_amplitude_1st": ap_max_list_amplitude_1st,
            "ap_max_list_voltage_baseline": ap_max_list_voltage_baseline,
            "ap_max_list_instantaneous_freq_1_2": ap_max_list_instantaneous_freq_1_2,
            "ap_max_list_average_frequency": ap_max_list_average_frequency,
            "ap_max_list_instantaneous_last": ap_max_list_instantaneous_last,
            "ap_max_list_freq_adaptation": ap_max_list_freq_adaptation,

            "ap_broadening_list_voltage_baseline": ap_broadening_list_voltage_baseline,
            "ap_broadening_list_half_duration_1st": ap_broadening_list_half_duration_1st,
            "ap_broadening_list_threshold_1st": ap_broadening_list_threshold_1st,
            "ap_broadening_list_threshold_2nd_1st": ap_broadening_list_threshold_2nd_1st,
            "ap_broadening_list_amplitude_1st": ap_broadening_list_amplitude_1st
        })

        # --- Save PDF ---
        plt.tight_layout()
        pdf_filename = f"{cell_count + 1:03d}_{os.path.splitext(file_name)[0]}.pdf"
        pdf_path = os.path.join(output_folder_traces, pdf_filename)

        with PdfPages(pdf_path) as pdf:
            pdf.savefig(fig)
            plt.close(fig)

            phase_fig, phase_axs = plt.subplots(4, 2, figsize=(8, 12))

            plot_phase_analysis_column(
                phase_axs[:, 0],
                ap_rheo_phase_data,
                "Rheobase first AP",
            )
            plot_phase_analysis_column(
                phase_axs[:, 1],
                ap_max_phase_data,
                "Max-AP trace first AP",
            )

            phase_fig.suptitle("First-AP derivative and phase-plot analysis", fontsize=12)
            phase_fig.tight_layout(rect=[0, 0, 1, 0.97])
            pdf.savefig(phase_fig)
            plt.close(phase_fig)

    # ==========================================================================================
    # --- Save all results to Excel ---
    # ==========================================================================================
    # results EXCEL file
    results_df = pd.DataFrame(results)
    excel_output_path = os.path.join(output_folder_results, "results.xlsx")

    list_columns = [
        "ap_max_list_current_steps",
        "ap_max_list_ap_numbers",
        "ap_max_list_half_duration_1st",
        "ap_max_list_threshold_1st",
        "ap_max_list_amplitude_1st",
        "ap_max_list_voltage_baseline",
        "ap_max_list_instantaneous_freq_1_2",
        "ap_max_list_average_frequency",
        "ap_max_list_instantaneous_last",
        "ap_max_list_freq_adaptation",
        "ap_broadening_list_voltage_baseline",
        "ap_broadening_list_half_duration_1st",
        "ap_broadening_list_threshold_1st",
        "ap_broadening_list_threshold_2nd_1st",
        "ap_broadening_list_amplitude_1st",
    ]
    export_df = results_df.drop(columns=list_columns, errors="ignore")
    export_df.to_excel(excel_output_path, index=False)

    # lists of ap max =========================================
    # Create DataFrames for lists
    ap_max_list_current_steps_data = []
    ap_max_list_ap_numbers_data = []
    ap_max_list_half_duration_1st_data = []
    ap_max_list_threshold_1st_data = []
    ap_max_list_amplitude_1st_data = []
    ap_max_list_voltage_baseline_data = []
    ap_max_list_average_frequency_data = []
    ap_max_list_instantaneous_freq_1_2_data = []
    ap_max_list_instantaneous_last_data = []
    ap_max_list_freq_adaptation_data = []

    # Process each cell's data
    for cell_count, row in metadata_df.iterrows():
        file_name = row['file_name']

        # Check if ap_max_series is valid before converting
        def is_valid_series(value):
            if pd.isna(value):
                return False
            try:
                int_val = int(float(value))
                return int_val > 0
            except (ValueError, TypeError):
                return False

        # Only process if ap_max_series is available and valid
        if not is_valid_series(row['ap_max_series']):
            # Skip this cell if ap_max_series is not valid
            # print(f"Skipping ap_max_list_current_steps/ap_max_list_ap_numbers export for {file_name} - invalid ap_max_series")
            continue

        series_id = int(float(row['ap_max_series'])) - 1

        # Get the stored ap_max_list_current_steps and ap_max_list_ap_numbers for this cell
        ap_max_list_current_steps_row = {'cell_count': cell_count + 1, 'file_name': file_name}
        ap_max_list_ap_numbers_row = {'cell_count': cell_count + 1, 'file_name': file_name}
        ap_max_list_instantaneous_freq_1_2_row = {'cell_count': cell_count + 1, 'file_name': file_name}
        ap_max_list_half_duration_1st_row = {'cell_count': cell_count + 1, 'file_name': file_name}
        ap_max_list_threshold_1st_row = {'cell_count': cell_count + 1, 'file_name': file_name}
        ap_max_list_amplitude_1st_row = {'cell_count': cell_count + 1, 'file_name': file_name}
        ap_max_list_voltage_baseline_row = {'cell_count': cell_count + 1, 'file_name': file_name}
        ap_max_list_average_frequency_row = {'cell_count': cell_count + 1, 'file_name': file_name}
        ap_max_list_instantaneous_last_row = {'cell_count': cell_count + 1, 'file_name': file_name}
        ap_max_list_freq_adaptation_row = {'cell_count': cell_count + 1, 'file_name': file_name}

        for i, (tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10) in enumerate(
                zip(results[cell_count]['ap_max_list_current_steps'],
                    results[cell_count]['ap_max_list_ap_numbers'],
                    results[cell_count]['ap_max_list_instantaneous_freq_1_2'],
                    results[cell_count]['ap_max_list_half_duration_1st'],
                    results[cell_count]['ap_max_list_threshold_1st'],
                    results[cell_count]['ap_max_list_amplitude_1st'],
                    results[cell_count]['ap_max_list_voltage_baseline'],
                    results[cell_count]['ap_max_list_average_frequency'],
                    results[cell_count]['ap_max_list_instantaneous_last'],
                    results[cell_count]['ap_max_list_freq_adaptation']), 1):
            ap_max_list_current_steps_row[f'sweep_{i}'] = tmp1
            ap_max_list_ap_numbers_row[f'sweep_{i}'] = tmp2
            ap_max_list_instantaneous_freq_1_2_row[f'sweep_{i}'] = tmp3
            ap_max_list_half_duration_1st_row[f'sweep_{i}'] = tmp4
            ap_max_list_threshold_1st_row[f'sweep_{i}'] = tmp5
            ap_max_list_amplitude_1st_row[f'sweep_{i}'] = tmp6
            ap_max_list_voltage_baseline_row[f'sweep_{i}'] = tmp7
            ap_max_list_average_frequency_row[f'sweep_{i}'] = tmp8
            ap_max_list_instantaneous_last_row[f'sweep_{i}'] = tmp9
            ap_max_list_freq_adaptation_row[f'sweep_{i}'] = tmp10

        ap_max_list_current_steps_data.append(ap_max_list_current_steps_row)
        ap_max_list_ap_numbers_data.append(ap_max_list_ap_numbers_row)
        ap_max_list_instantaneous_freq_1_2_data.append(ap_max_list_instantaneous_freq_1_2_row)
        ap_max_list_half_duration_1st_data.append(ap_max_list_half_duration_1st_row)
        ap_max_list_threshold_1st_data.append(ap_max_list_threshold_1st_row)
        ap_max_list_amplitude_1st_data.append(ap_max_list_amplitude_1st_row)
        ap_max_list_voltage_baseline_data.append(ap_max_list_voltage_baseline_row)
        ap_max_list_average_frequency_data.append(ap_max_list_average_frequency_row)
        ap_max_list_instantaneous_last_data.append(ap_max_list_instantaneous_last_row)
        ap_max_list_freq_adaptation_data.append(ap_max_list_freq_adaptation_row)

    # ap_broadening series
    ap_broadening_list_voltage_baseline_data = []
    ap_broadening_list_half_duration_1st_data = []
    ap_broadening_list_threshold_1st_data = []
    ap_broadening_list_threshold_2nd_1st_data = []
    ap_broadening_list_amplitude_1st_data = []

    # Process each cell's data

    for cell_count, row in metadata_df.iterrows():
        file_name = row['file_name']

        # Only process if ap_broadening_series is available and valid
        if not is_valid_series(row['ap_broadening_series']):
            # Skip this cell if ap_broadening_series is not valid
            continue

        series_id = int(float(row['ap_broadening_series'])) - 1

        # Get the stored ap_broadening_list_current_steps and ap_broadening_list_ap_numbers for this cell
        ap_broadening_list_voltage_baseline_row = {'cell_count': cell_count + 1, 'file_name': file_name}
        ap_broadening_list_half_duration_1st_row = {'cell_count': cell_count + 1, 'file_name': file_name}
        ap_broadening_list_threshold_1st_row = {'cell_count': cell_count + 1, 'file_name': file_name}
        ap_broadening_list_threshold_2nd_1st_row = {'cell_count': cell_count + 1, 'file_name': file_name}
        ap_broadening_list_amplitude_1st_row = {'cell_count': cell_count + 1, 'file_name': file_name}

        for i, (tmp1, tmp2 ,tmp3, tmp4, tmp5) in enumerate(zip(results[cell_count]['ap_broadening_list_voltage_baseline'],
                                                                    results[cell_count]['ap_broadening_list_half_duration_1st'],
                                                                    results[cell_count]['ap_broadening_list_threshold_1st'],
                                                                    results[cell_count]['ap_broadening_list_threshold_2nd_1st'],
                                                                    results[cell_count]['ap_broadening_list_amplitude_1st']), 1):
            ap_broadening_list_voltage_baseline_row[f'sweep_{i}'] = tmp1
            ap_broadening_list_half_duration_1st_row[f'sweep_{i}'] = tmp2
            ap_broadening_list_threshold_1st_row[f'sweep_{i}'] = tmp3
            ap_broadening_list_threshold_2nd_1st_row[f'sweep_{i}'] = tmp4
            ap_broadening_list_amplitude_1st_row[f'sweep_{i}'] = tmp5

        ap_broadening_list_voltage_baseline_data.append(ap_broadening_list_voltage_baseline_row)
        ap_broadening_list_half_duration_1st_data.append(ap_broadening_list_half_duration_1st_row)
        ap_broadening_list_threshold_1st_data.append(ap_broadening_list_threshold_1st_row)
        ap_broadening_list_threshold_2nd_1st_data.append(ap_broadening_list_threshold_2nd_1st_row)
        ap_broadening_list_amplitude_1st_data.append(ap_broadening_list_amplitude_1st_row)

    # exporting all lists of ap max and ap broadening =========================================
    # ap_max series
    ap_max_list_current_steps_df = pd.DataFrame(ap_max_list_current_steps_data)
    ap_max_list_current_steps_excel_path = os.path.join(output_folder_results, "ap_max_list_current_steps.xlsx")
    ap_max_list_current_steps_df.to_excel(ap_max_list_current_steps_excel_path, index=False)

    ap_max_list_ap_numbers_df = pd.DataFrame(ap_max_list_ap_numbers_data)
    ap_max_list_ap_numbers_excel_path = os.path.join(output_folder_results, "ap_max_list_ap_numbers.xlsx")
    ap_max_list_ap_numbers_df.to_excel(ap_max_list_ap_numbers_excel_path, index=False)

    ap_max_list_instantaneous_freq_1_2_df = pd.DataFrame(ap_max_list_instantaneous_freq_1_2_data)
    ap_max_list_instantaneous_freq_1_2_excel_path = os.path.join(output_folder_results,
                                                                 "ap_max_list_instantaneous_freq_1_2.xlsx")
    ap_max_list_instantaneous_freq_1_2_df.to_excel(ap_max_list_instantaneous_freq_1_2_excel_path, index=False)

    ap_max_list_half_duration_1st_df = pd.DataFrame(ap_max_list_half_duration_1st_data)
    ap_max_list_half_duration_1st_excel_path = os.path.join(output_folder_results, "ap_max_list_half_duration_1st.xlsx")
    ap_max_list_half_duration_1st_df.to_excel(ap_max_list_half_duration_1st_excel_path, index=False)

    ap_max_list_threshold_1st_df = pd.DataFrame(ap_max_list_threshold_1st_data)
    ap_max_list_threshold_1st_excel_path = os.path.join(output_folder_results, "ap_max_list_threshold_1st.xlsx")
    ap_max_list_threshold_1st_df.to_excel(ap_max_list_threshold_1st_excel_path, index=False)

    ap_max_list_amplitude_1st_df = pd.DataFrame(ap_max_list_amplitude_1st_data)
    ap_max_list_amplitude_1st_excel_path = os.path.join(output_folder_results, "ap_max_list_amplitude_1st.xlsx")
    ap_max_list_amplitude_1st_df.to_excel(ap_max_list_amplitude_1st_excel_path, index=False)

    # ap_broadening series
    ap_broadening_list_voltage_baseline_df = pd.DataFrame(ap_broadening_list_voltage_baseline_data)
    ap_broadening_list_voltage_baseline_excel_path = os.path.join(output_folder_results,
                                                                       "ap_broadening_list_voltage_baseline.xlsx")
    ap_broadening_list_voltage_baseline_df.to_excel(ap_broadening_list_voltage_baseline_excel_path,
                                                         index=False)

    ap_broadening_list_half_duration_1st_df = pd.DataFrame(ap_broadening_list_half_duration_1st_data)
    ap_broadening_list_half_duration_1st_excel_path = os.path.join(output_folder_results,
                                                                   "ap_broadening_list_half_duration_1st.xlsx")
    ap_broadening_list_half_duration_1st_df.to_excel(ap_broadening_list_half_duration_1st_excel_path, index=False)

    ap_broadening_list_threshold_1st_df = pd.DataFrame(ap_broadening_list_threshold_1st_data)
    ap_broadening_list_threshold_1st_excel_path = os.path.join(output_folder_results, "ap_broadening_list_threshold_1st.xlsx")
    ap_broadening_list_threshold_1st_df.to_excel(ap_broadening_list_threshold_1st_excel_path, index=False)

    ap_broadening_list_threshold_2nd_1st_df = pd.DataFrame(ap_broadening_list_threshold_2nd_1st_data)
    ap_broadening_list_threshold_2nd_1st_excel_path = os.path.join(output_folder_results, "ap_broadening_list_threshold_2nd_1st.xlsx")
    ap_broadening_list_threshold_2nd_1st_df.to_excel(ap_broadening_list_threshold_2nd_1st_excel_path, index=False)

    ap_broadening_list_amplitude_1st_df = pd.DataFrame(ap_broadening_list_amplitude_1st_data)
    ap_broadening_list_amplitude_1st_excel_path = os.path.join(output_folder_results, "ap_broadening_list_amplitude_1st.xlsx")
    ap_broadening_list_amplitude_1st_df.to_excel(ap_broadening_list_amplitude_1st_excel_path, index=False)

    ap_max_list_average_frequency_df = pd.DataFrame(ap_max_list_average_frequency_data)
    ap_max_list_average_frequency_excel_path = os.path.join(output_folder_results, "ap_max_list_average_frequency.xlsx")
    ap_max_list_average_frequency_df.to_excel(ap_max_list_average_frequency_excel_path, index=False)

    ap_max_list_voltage_baseline_df = pd.DataFrame(ap_max_list_voltage_baseline_data)
    ap_max_list_voltage_baseline_excel_path = os.path.join(output_folder_results, "ap_max_list_voltage_baseline.xlsx")
    ap_max_list_voltage_baseline_df.to_excel(ap_max_list_voltage_baseline_excel_path, index=False)

    ap_max_list_instantaneous_last_df = pd.DataFrame(ap_max_list_instantaneous_last_data)
    ap_max_list_instantaneous_last_excel_path = os.path.join(output_folder_results,
                                                             "ap_max_list_instantaneous_last.xlsx")
    ap_max_list_instantaneous_last_df.to_excel(ap_max_list_instantaneous_last_excel_path, index=False)

    ap_max_list_freq_adaptation_df = pd.DataFrame(ap_max_list_freq_adaptation_data)
    ap_max_list_freq_adaptation_excel_path = os.path.join(output_folder_results, "ap_max_list_freq_adaptation.xlsx")
    ap_max_list_freq_adaptation_df.to_excel(ap_max_list_freq_adaptation_excel_path, index=False)

    # save points =========================================
    # Save analysis points to JSON file
    analysis_points_path = os.path.join(config.IMPORT_FOLDER, "analysis_points.json")
    analysis_points_output_path = os.path.join(output_folder_used_data_and_code, "analysis_points.json")

    # Convert defaultdict to regular dict for JSON serialization
    analysis_points_dict = json.loads(json.dumps(analysis_points, default=lambda o: {str(k): v for k, v in o.items()}))

    # Save to both locations
    with open(analysis_points_path, 'w') as f:
        json.dump(analysis_points_dict, f, indent=2)
    with open(analysis_points_output_path, 'w') as f:
        json.dump(analysis_points_dict, f, indent=2)

    print(f"Analysis complete. Results saved to {excel_output_path}")
    print(f"Analysis points saved to {analysis_points_path}")


def start_browser():
    # Import and start the browser
    from browser import app, win
    win.show()
    if sys.flags.interactive == 0:
        app.exec_()


if __name__ == '__main__':
    CC_eval()
    start_browser()  # This will start the browser after CC_eval completes