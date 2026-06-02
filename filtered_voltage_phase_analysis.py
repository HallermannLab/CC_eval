# ... existing code ...
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import savgol_filter, find_peaks
import heka_reader
import git_save as myGit
from collections import defaultdict
import json
from statistics import median


# --- parameters ---
A_to_pA = 1e12
V_to_mV = 1e3


def calculate_filtered_trace_and_derivatives(time, voltage, smooth_window):
    """Calculate filtered voltage plus first and second derivatives for plotting/phase analysis.

    voltage is expected in mV.
    d1 is therefore initially mV/s, and d1_in_V_per_s is V/s.
    d2 is initially V/s² after gradient of d1_in_V_per_s; d2_in_V_per_s_s is V/s².
    """
    sg_polyorder = 3
    dt = time[1] - time[0]

    sg_window_s = smooth_window / 1000.0
    win_samples = int(round(sg_window_s / dt))
    if win_samples <= sg_polyorder + 1:
        win_samples = sg_polyorder + 3
    if win_samples % 2 == 0:
        win_samples += 1
    if win_samples >= len(voltage):
        win_samples = len(voltage) - 1 if len(voltage) % 2 == 0 else len(voltage)
    if win_samples <= sg_polyorder:
        win_samples = sg_polyorder + 2
        if win_samples % 2 == 0:
            win_samples += 1

    voltage_filt = savgol_filter(voltage, window_length=win_samples, polyorder=sg_polyorder)

    d1 = np.gradient(voltage_filt, dt)
    d1_in_V_per_s = savgol_filter(d1, window_length=win_samples, polyorder=sg_polyorder) / V_to_mV

    d2 = np.gradient(d1_in_V_per_s, dt)
    d2_in_V_per_s_s = savgol_filter(d2, window_length=win_samples, polyorder=sg_polyorder)

    return voltage_filt, d1, d1_in_V_per_s, d2, d2_in_V_per_s_s


def get_first_phase_peak(time, voltage_filt, d1_in_V_per_s, peak_time, phase_plot_t1, phase_plot_t2,
                         phase_plot_range_start, phase_plot_range_end):
    """Find the first local peak in the phase plot within the requested voltage range.

    The search is restricted to phase_plot_t1 ms before and phase_plot_t2 ms after the AP peak.
    The phase plot is d1_in_V_per_s vs. voltage_filt.
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

    search_d1 = d1_in_V_per_s[search_indices]
    peaks, _ = find_peaks(search_d1)

    if len(peaks) == 0:
        return None

    phase_peak_idx = search_indices[peaks[0]]
    return {
        "time": float(time[phase_peak_idx]),
        "voltage_filt": float(voltage_filt[phase_peak_idx]),
        "d1_in_V_per_s": float(d1_in_V_per_s[phase_peak_idx]),
    }


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

    for point_type, (marker, color, label) in voltage_symbols.items():
        if points.get(point_type):
            t_points, v_points = zip(*points[point_type])
            ax_voltage.scatter(t_points, [v * V_to_mV for v in v_points],
                               marker=marker, color=color, s=25, label=label, zorder=5)

    for point_type, (marker, color, label) in d1_symbols.items():
        if points.get(point_type):
            t_points, y_points = zip(*points[point_type])
            ax_d1.scatter(t_points, y_points, marker=marker, color=color, s=25, label=label, zorder=5)

    for point_type, (marker, color, label) in d2_symbols.items():
        if points.get(point_type):
            t_points, y_points = zip(*points[point_type])
            ax_d2.scatter(t_points, y_points, marker=marker, color=color, s=25, label=label, zorder=5)

    if points.get('threshold'):
        _, threshold_v = points['threshold'][0]
        ax_phase.axvline(threshold_v * V_to_mV, color='red', linestyle='--', alpha=0.4)

    if points.get('threshold_2nd'):
        _, threshold_2nd_v = points['threshold_2nd'][0]
        ax_phase.axvline(threshold_2nd_v * V_to_mV, color='brown', linestyle='--', alpha=0.4)

    if phase_peak is not None:
        ax_phase.scatter(
            phase_peak["voltage_filt"],
            phase_peak["d1_in_V_per_s"],
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
    d1_in_V_per_s = phase_data["d1_in_V_per_s"]
    d2 = phase_data["d2"]
    d2_in_V_per_s_s = phase_data["d2_in_V_per_s_s"]
    plot_mask = phase_data["plot_mask"]
    points = phase_data["points"]
    phase_peak = phase_data["phase_peak"]

    ax_voltage.plot(time[plot_mask], voltage[plot_mask], color='0.65', label='voltage')
    ax_voltage.plot(time[plot_mask], voltage_filt[plot_mask], color='black', label='voltage_filt')
    ax_voltage.set_title(f"{title_prefix}: AP")
    ax_voltage.set_ylabel("Voltage (mV)")
    ax_voltage.set_xlabel("Time (s)")
    ax_voltage.grid(True)

    ax_d1.plot(time[plot_mask], (d1 / V_to_mV)[plot_mask], color='tab:blue', label='d1 / V_to_mV')
    ax_d1.plot(time[plot_mask], d1_in_V_per_s[plot_mask], color='black', label='d1_in_V_per_s')
    ax_d1.set_title(f"{title_prefix}: 1st derivative")
    ax_d1.set_ylabel("dV/dt (V/s)")
    ax_d1.set_xlabel("Time (s)")
    ax_d1.grid(True)

    ax_d2.plot(time[plot_mask], (d2 / V_to_mV)[plot_mask], color='tab:blue', label='d2 / V_to_mV')
    ax_d2.plot(time[plot_mask], d2_in_V_per_s_s[plot_mask], color='black', label='d2_in_V_per_s_s')
    ax_d2.set_title(f"{title_prefix}: 2nd derivative")
    ax_d2.set_ylabel("d²V/dt²")
    ax_d2.set_xlabel("Time (s)")
    ax_d2.grid(True)

    ax_phase.plot(voltage_filt[plot_mask], d1_in_V_per_s[plot_mask], color='black')
    ax_phase.set_title(f"{title_prefix}: phase plot")
    ax_phase.set_ylabel("dV/dt (V/s)")
    ax_phase.set_xlabel("voltage_filt (mV)")
    ax_phase.grid(True)

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

    voltage_filt, d1, d1_in_V_per_s, d2, d2_in_V_per_s_s = calculate_filtered_trace_and_derivatives(
        time, voltage, smooth_window
    )

    plot_start = peak_time - phase_plot_t1 / 1000.0
    plot_end = peak_time + phase_plot_t2 / 1000.0
    plot_mask = (time >= plot_start) & (time <= plot_end)

    phase_peak = get_first_phase_peak(
        time,
        voltage_filt,
        d1_in_V_per_s,
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
        "d1_in_V_per_s": d1_in_V_per_s,
        "d2": d2,
        "d2_in_V_per_s_s": d2_in_V_per_s_s,
        "plot_mask": plot_mask,
        "points": points,
        "phase_peak": phase_peak,
    }


def ap_analysis(time, voltage, v_threshold, dvdt_threshold, smooth_window, fraction_of_max_of_2nd_derivative,
                window_for_searching_threshold, window_for_searching_ahp,
                minimal_ap_interval, minimal_ap_duration, maximal_ap_duration,
                maximal_relative_amplitude_decline):

    sg_polyorder = 3
# ... existing code ...
