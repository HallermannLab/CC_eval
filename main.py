import os
import sys
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#from scipy.misc import derivative
from scipy.stats import linregress
#from scipy import signal
from scipy.signal import savgol_filter
import heka_reader
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtCore
import git_save as myGit
from config import ROOT_FOLDER, IMPORT_FOLDER, METADATA_FILE, EXTERNAL_DATA_FOLDER
#from config import analysis_points
from collections import defaultdict
import json
from trace_viewer import show_three_traces



# --- parameters ---
A_to_pA = 1e12
V_to_mV = 1e3


"""
These and some more values are now imported from the metadata file (metadata.xlsx)
# windows for analysis estimating the current jump
window1_rin_start = 0.01
window1_rin_end = 0.09
window2_rin_start = 0.30
window2_rin_end = 0.39

window1_ap_rheo_start = 0.01
window1_ap_rheo_end = 0.09
window2_ap_rheo_start = 0.30
window2_ap_rheo_end = 0.39

window1_ap_max_start = 0.01
window1_ap_max_end = 0.09
window2_ap_max_start = 0.30
window2_ap_max_end = 0.39

# parameters for AP detection
v_threshold = -20
dvdt_threshold = 10
smooth_window = 1500 #Hz
fraction_of_max_of_2nd_derivative = 0.3
window_for_searching_threshold = 0.002
window_for_searching_ahp = 0.01
minimal_ap_interval = 0.001
minimal_ap_duration = 0.0005
maximal_ap_duration = 0.005
"""


def ap_analysis(time, voltage, v_threshold, dvdt_threshold, smooth_window, fraction_of_max_of_2nd_derivative,
                window_for_searching_threshold, window_for_searching_ahp,
                minimal_ap_interval, minimal_ap_duration, maximal_ap_duration,
                maximal_relative_amplitude_decline):

    sg_polyorder = 3

    """
    Analyze action potential features from voltage trace

    Args:
        time: time points array (seconds)
        voltage: voltage trace array (mV)
        v_threshold: voltage threshold for AP detection (mV)
        dvdt_threshold: dV/dt threshold for AP detection (V/s)
        smooth_window: smoothing window for voltage traces (ms)
        window_for_searching_threshold: time window for finding threshold (s)
        window_for_searching_ahp: time window for finding AHP (s)
        minimal_ap_interval: minimum time between APs (s)
        minimal_ap_duration: minimum AP duration (s)
        maximal_ap_duration: maximum AP duration (s)
        maximal_relative_amplitude_decline: maximum allowed amplitude decline

    Returns:
        ap_number: number of APs detected
        th_v: threshold voltages
        th_t: threshold times
        th_v_2nd: threshold voltages from 2nd derivative
        th_t_2nd: threshold times from 2nd derivative
        hd_start_t: half-duration start times
        hd_start_v: half-duration start voltages
        hd_end_t: half-duration end times
        hd_end_v: half-duration end voltages
        p_v: peak voltages
        p_t: peak times
        ahp_v: after-hyperpolarization voltages
        ahp_t: after-hyperpolarization times
        dvdt_v: max dV/dt values
        dvdt_t: max dV/dt times
    """

    """
    # use Bessel filter for low-pass filtering
    # Apply low-pass filter to voltage trace
    dt = time[1] - time[0]
    nyquist = 0.5 / dt
    cutoff_normalized = smooth_window / nyquist
    b, a = signal.butter(2, cutoff_normalized, 'low')
    voltage_filt = signal.filtfilt(b, a, voltage)

    # Calculate 1st and 2nd derivatives using filtered trace
    d1 = np.diff(voltage_filt) / dt
    d1_filt = signal.filtfilt(b, a, d1) # CAVE overwrites the trace with the filtered one
    d1_in_V_per_s = d1_filt / V_to_mV
    d2 = np.diff(d1_in_V_per_s) / dt
    d2_filt = signal.filtfilt(b, a, d2)  # CAVE overwrites the trace with the filtered one
    d2_in_V_per_s_s = d2_filt / V_to_mV
    """
    # use Savitzky-Golay filter for smoothing
    dt = time[1] - time[0]
    sg_window_s = smooth_window / 1000.0
    win_samples = int(round(sg_window_s / dt))
    # make odd and at least polyorder+2
    if win_samples <= sg_polyorder + 1:
        win_samples = sg_polyorder + 3
    if win_samples % 2 == 0:
        win_samples += 1

    # smooth voltages and use numerical derivatives (central differences via np.gradient)
    voltage_filt = savgol_filter(voltage, window_length=win_samples, polyorder=sg_polyorder)
    d1 = np.gradient(voltage_filt, dt)
    d1_filt = savgol_filter(d1, window_length=win_samples, polyorder=sg_polyorder)
    d1_in_V_per_s = d1_filt / V_to_mV
    d2 = np.gradient(d1_in_V_per_s, dt)
    d2_filt = savgol_filter(d2, window_length=win_samples, polyorder=sg_polyorder)  # CAVE overwrites the trace with the filtered one
    d2_in_V_per_s_s = d2_filt / V_to_mV

    """
    to check difference with browser.py
    voltage_filt = savgol_filter(data, window_length=win_samples, polyorder=sg_polyorder)
    d1 = np.gradient(voltage_filt, dt)
    d1_in_V_per_s = savgol_filter(d1, window_length=win_samples, polyorder=sg_polyorder)
    d2 = np.gradient(d1_in_V_per_s, dt)
    d2 = d2 / V_to_mV
    d2 = d2 / V_to_mV
    # the 1st is V/s and 2nd is V/ms^2
    d2_in_V_per_s_s = savgol_filter(d2, window_length=win_samples,
                                    polyorder=sg_polyorder)  # CAVE overwrites the trace with the filtered one
    """

    # Find threshold crossings where voltage crosses v_threshold
    threshold_idx = np.where((voltage[:-1] < v_threshold) & (voltage[1:] >= v_threshold))[0]

    # Initialize results
    th_v = []
    th_t = []
    th_v_2nd = []
    th_t_2nd = []
    p_v = []
    p_t = []
    ahp_v = []
    ahp_t = []
    hd_start_v = []
    hd_start_t = []
    hd_end_v = []
    hd_end_t = []
    dvdt_v = []
    dvdt_t = []

    for idx in threshold_idx:
        # Get window indices
        start_idx = max(0, idx - int(window_for_searching_threshold / dt))
        end_idx = min(len(voltage) - 1, idx + int(window_for_searching_threshold / dt))

        # Find precise threshold where dV/dt crosses threshold
        window_d1 = d1_in_V_per_s[start_idx:end_idx]
        th_idx_candidates = np.where(window_d1 >= dvdt_threshold)[0]
        if len(th_idx_candidates) == 0:
            # Handle case where no threshold crossing is found
            continue  # Skip this iteration of the loop
        # Find threshold based on dV/dt
        th_idx = start_idx + th_idx_candidates[0]

        # Find threshold based on 2nd derivative crossing 30% of max
        window_d2 = d2_in_V_per_s_s[start_idx:end_idx]
        d2_max = np.max(window_d2)
        d2_threshold = fraction_of_max_of_2nd_derivative * d2_max

        # Find first crossing of 30% threshold - but ensure it's actually a crossing
        # Look for points that cross the threshold from below
        crossing_indices = []
        for i in range(1, len(window_d2)):
            if window_d2[i-1] < d2_threshold and window_d2[i] >= d2_threshold:
                crossing_indices.append(i)
        
        if len(crossing_indices) > 0:
            th_idx_2nd = start_idx + crossing_indices[0]
        else:
            print(f"Warning: Could not find 2nd derivative threshold crossing in. Using maximum instead.")
            th_idx_2nd = start_idx + np.argmax(window_d2)  # Fallback to max if no crossing found

        peak_idx = th_idx + np.argmax(voltage[th_idx:th_idx + int(maximal_ap_duration / dt)])

        # Find AHP
        ahp_window = voltage[peak_idx:peak_idx + int(window_for_searching_ahp / dt)]
        ahp_idx = peak_idx + np.argmin(ahp_window)

        # Max dV/dt
        dvdt_idx = th_idx + np.argmax(d1_in_V_per_s[th_idx:peak_idx])

        # Calculate half-duration points
        half_amplitude = (voltage[peak_idx] - voltage[th_idx]) / 2 + voltage[th_idx]

        hd_start_idx = None
        hd_end_idx = None
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
                hd_start_idx = i
                break

        # Find first crossing after peak
        for i in range(peak_idx, ahp_idx):
            if voltage[i] <= half_amplitude:
                v1, v2 = voltage[i - 1], voltage[i]
                t1, t2 = time[i - 1], time[i]
                frac = (half_amplitude - v1) / (v2 - v1) if v2 != v1 else 0.0
                hd_end_time = t1 + frac * (t2 - t1)
                hd_end_idx = i
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
            th_t.append(time[th_idx])
            th_v_2nd.append(voltage[th_idx_2nd])
            th_t_2nd.append(time[th_idx_2nd])
            p_v.append(voltage[peak_idx])
            p_t.append(time[peak_idx])
            ahp_v.append(voltage[ahp_idx])
            ahp_t.append(time[ahp_idx])
            dvdt_v.append(voltage[dvdt_idx])
            dvdt_t.append(time[dvdt_idx])
            hd_start_v.append(half_amplitude)  
            hd_start_t.append(hd_start_time)
            hd_end_v.append(half_amplitude)
            hd_end_t.append(hd_end_time)

    ap_number = len(p_v)

    return ap_number, th_v, th_t, th_v_2nd, th_t_2nd, hd_start_t, hd_start_v, hd_end_t, hd_end_v, p_v, p_t, ahp_v, ahp_t, dvdt_v, dvdt_t


def CC_eval():

    # --- Create Output Folders ---
    timestamp = datetime.now().strftime("%Y-%m-%d__%H-%M-%S")
    output_folder = os.path.join(ROOT_FOLDER, f"output_SH_{timestamp}")
    os.makedirs(output_folder, exist_ok=True)

    output_folder_results = os.path.join(output_folder, "results")
    os.makedirs(output_folder_results, exist_ok=True)

    output_folder_traces = os.path.join(output_folder, "traces")
    os.makedirs(output_folder_traces, exist_ok=True)

    output_folder_used_data_and_code = os.path.join(output_folder, "used_data_and_code")
    os.makedirs(output_folder_used_data_and_code, exist_ok=True)

    # --- Load Metadata ---
    metadata_df = pd.read_excel(METADATA_FILE)
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

        v_threshold = row['v_threshold']
        dvdt_threshold = row['dvdt_threshold']
        smooth_window = row['smooth_window']
        fraction_of_max_of_2nd_derivative = row['fraction_of_max_of_2nd_derivative']
        window_for_searching_threshold = row['window_for_searching_threshold']
        window_for_searching_ahp = row['window_for_searching_ahp']
        minimal_ap_interval = row['minimal_ap_interval']
        minimal_ap_duration = row['minimal_ap_duration']
        maximal_ap_duration = row['maximal_ap_duration']
        maximal_relative_amplitude_decline = row['maximal_relative_amplitude_decline']

        dat_path = os.path.join(EXTERNAL_DATA_FOLDER, file_name)
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

        ap_rheo_half_duration_1st = None
        ap_rheo_threshold_1st = None
        ap_rheo_threshold_2nd_1st = None
        ap_rheo_amplitude_1st = None

        ap_rheo_half_duration_av = None
        ap_rheo_threshold_av = None
        ap_rheo_threshold_2nd_av = None
        ap_rheo_amplitude_av = None

        ap_max_half_duration_1st = None
        ap_max_threshold_1st = None
        ap_max_threshold_2nd_1st = None
        ap_max_amplitude_1st = None

        ap_max_half_duration_av = None
        ap_max_threshold_av = None
        ap_max_threshold_2nd_av = None
        ap_max_amplitude_av = None

        ap_max_instantaneous_freq_1_2 = None
        ap_max_average_freq_all = None

        ap_max_list_current_steps = []
        ap_max_list_ap_numbers = []
        ap_max_list_instantaneous_freq_1_2 = []
        ap_max_list_half_duration_1st = []
        ap_max_list_threshold_1st = []
        ap_max_list_amplitude_1st = []

        ap_broadening_list_voltage_baseline = []
        ap_broadening_list_half_duration_1st = []
        ap_broadening_list_threshold_1st = []
        ap_broadening_list_amplitude_1st = []


        # ==========================================================================================
        # --- RMP ---
        # ==========================================================================================

        if rmp_series_available:
            series_id = rmp_series
            
            n_sweeps = bundle.pul[group_id][series_id].NumberSweeps
            voltage_trace = V_to_mV * bundle.data[group_id, series_id, 0, 0]

            n_points = len(voltage_trace)
            sampling_interval = bundle.pul[group_id][series_id][0][0].XInterval
            time = np.arange(n_points) * sampling_interval

            rmp_mean = voltage_trace.mean()
            rmp_min = voltage_trace.min()

            # Plot RMP trace with custom axis
            axs[0].plot(time, voltage_trace)
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
                voltage = V_to_mV * bundle.data[group_id, series_id, sweep_id, 0]
                current = A_to_pA * bundle.data[group_id, series_id, sweep_id, 1]

                dv = voltage[idx2].mean() - voltage[idx1].mean()
                di = current[idx2].mean() - current[idx1].mean()

                if not (rin_series_from_current <= di <= rin_series_to_current):
                    continue

                # Plot current trace
                axs[2].plot(time, current, alpha=0.5)

                # Plot voltage trace
                axs[3].plot(time, voltage, alpha=0.5)

                delta_vs.append(dv)
                delta_is.append(di)

            # Add vertical lines for analysis windows in both current and voltage plots
            for ax in [axs[2], axs[3]]:
                ax.axvline(x=window1_rin_start, color='r', linestyle='--', alpha=0.3)
                ax.axvline(x=window1_rin_end, color='r', linestyle='--', alpha=0.3)
                ax.axvline(x=window2_rin_start, color='r', linestyle='--', alpha=0.3)
                ax.axvline(x=window2_rin_end, color='r', linestyle='--', alpha=0.3)
                ax.grid(True)

            delta_vs = np.array(delta_vs)
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
                voltage = V_to_mV * bundle.data[group_id, series_id, sweep_id, 0]
                current = A_to_pA * bundle.data[group_id, series_id, sweep_id, 1]

                # Plot voltage trace in superposition plot
                axs[4].plot(time, voltage, alpha=0.5, label=f'Sweep {sweep_id + 1}')

                di = current[idx2].mean() - current[idx1].mean()

                ap_number, th_v, th_t, th_v_2nd, th_t_2nd, hd_start_t, hd_start_v, hd_end_t, hd_end_v, p_v, p_t, ahp_v, ahp_t, dvdt_v, dvdt_t = ap_analysis(
                    time, voltage, v_threshold, dvdt_threshold, smooth_window, fraction_of_max_of_2nd_derivative,
                    window_for_searching_threshold, window_for_searching_ahp,
                    minimal_ap_interval, minimal_ap_duration, maximal_ap_duration,
                    maximal_relative_amplitude_decline)

                #print(f"Sweep {sweep_id + 1}: ap_number = {ap_number}")
                if ap_number > 0:
                    sweep_points = {
                        'threshold': list(zip(th_t, [v / V_to_mV for v in th_v])),
                        'threshold_2nd': list(zip(th_t_2nd, [v / V_to_mV for v in th_v_2nd])),
                        'half_duration_start': list(zip(hd_start_t, [v / V_to_mV for v in hd_start_v])),
                        'half_duration_end': list(zip(hd_end_t, [v / V_to_mV for v in hd_end_v])),
                        'peak': list(zip(p_t, [v / V_to_mV for v in p_v])),
                        'ahp': list(zip(ahp_t, [v / V_to_mV for v in ahp_v])),
                        'dvdt_max': list(zip(dvdt_t, [v / V_to_mV for v in dvdt_v])),
                        'smooth_window': smooth_window
                    }
                    # Store the analysis points in the nested dictionary
                    analysis_points[file_name][group_id][series_id][sweep_id][0] = sweep_points

                if rheobase is None and ap_number > 0:
                    rheobase_voltage_trace = voltage
                    rheobase = di
                    ap_rheo_half_duration_1st = hd_end_t[0] - hd_start_t[0]
                    ap_rheo_threshold_1st = th_v[0]
                    ap_rheo_threshold_2nd_1st = th_v_2nd[0]
                    ap_rheo_amplitude_1st = p_v[0] - th_v[0]
                    # Calculate average AP parameters
                    ap_rheo_half_duration_av = sum(hd_end_t[i] - hd_start_t[i] for i in range(ap_number)) / ap_number
                    ap_rheo_threshold_av = sum(th_v) / len(th_v)
                    ap_rheo_threshold_2nd_av = sum(th_v_2nd) / len(th_v_2nd)
                    ap_rheo_amplitude_av = sum(p_v[i] - th_v[i] for i in range(len(th_v))) / len(th_v)

            # Formatting for AP rheo superposition plot
            axs[4].grid(True)

            # Setup rheobase trace plot
            axs[5].set_title("Rheobase Voltage Trace")
            axs[5].set_ylabel("Voltage (mV)")
            axs[5].set_xlabel("Time (s)")
            if rheobase is not None:
                axs[5].plot(time, rheobase_voltage_trace, label=f'Rheobase: {rheobase:.1f} pA')
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
                voltage = V_to_mV * bundle.data[group_id, series_id, sweep_id, 0]
                current = A_to_pA * bundle.data[group_id, series_id, sweep_id, 1]

                # Plot voltage trace in superposition plot
                axs[6].plot(time, voltage, alpha=0.5, label=f'Sweep {sweep_id + 1}')

                di = current[idx2].mean() - current[idx1].mean()
                ap_max_list_current_steps.append(float(di))  # Convert numpy.float64 to Python float

                ap_number, th_v, th_t, th_v_2nd, th_t_2nd, hd_start_t, hd_start_v, hd_end_t, hd_end_v, p_v, p_t, ahp_v, ahp_t, dvdt_v, dvdt_t = ap_analysis(
                    time, voltage, v_threshold, dvdt_threshold, smooth_window, fraction_of_max_of_2nd_derivative,
                    window_for_searching_threshold, window_for_searching_ahp,
                    minimal_ap_interval, minimal_ap_duration, maximal_ap_duration,
                    maximal_relative_amplitude_decline)

                ap_max_list_ap_numbers.append(int(ap_number))  # Convert to Python int
                #print(f"Sweep {sweep_id + 1}: ap_number = {ap_number}")

                if ap_number > 0:
                    sweep_points = {
                        'threshold': list(zip(th_t, [v / V_to_mV for v in th_v])),
                        'threshold_2nd': list(zip(th_t_2nd, [v / V_to_mV for v in th_v_2nd])),
                        'half_duration_start': list(zip(hd_start_t, [v / V_to_mV for v in hd_start_v])),
                        'half_duration_end': list(zip(hd_end_t, [v / V_to_mV for v in hd_end_v])),
                        'peak': list(zip(p_t, [v / V_to_mV for v in p_v])),
                        'ahp': list(zip(ahp_t, [v / V_to_mV for v in ahp_v])),
                        'dvdt_max': list(zip(dvdt_t, [v / V_to_mV for v in dvdt_v])),
                        'smooth_window': smooth_window
                    }
                    # Store the analysis points in the nested dictionary
                    analysis_points[file_name][group_id][series_id][sweep_id][0] = sweep_points
                    # Log the structure of analysis points:
                    # print(f"Storing points: file_name={file_name}, group_id={group_id}, series_id={series_id}, sweep_id={sweep_id}, trace_id={0}")
                    # print(f"Points to store: {json.dumps(sweep_points, indent=2)}")

                    ap_max_list_half_duration_1st.append(float(hd_end_t[0] - hd_start_t[0]))
                    ap_max_list_threshold_1st.append(float(th_v[0]))
                    ap_max_list_amplitude_1st.append(float(p_v[0] - th_v[0]))
                else:
                    ap_max_list_half_duration_1st.append(None)
                    ap_max_list_threshold_1st.append(None)
                    ap_max_list_amplitude_1st.append(None)

                # Calculate AP frequencies
                if ap_number >= 2:
                    # Instantaneous frequency of first two APs (Hz)
                    ap_max_instantaneous_freq_1_2 = 1.0 / (p_t[1] - p_t[0])
                    # Average frequency of all APs (Hz)
                    total_time = p_t[-1] - p_t[0]  # Time from first to last AP
                    ap_max_average_freq_all = (ap_number - 1) / total_time  # Number of intervals / total time
                    ap_max_list_instantaneous_freq_1_2.append(float(ap_max_instantaneous_freq_1_2))
                else:
                    ap_max_instantaneous_freq_1_2 = None
                    ap_max_average_freq_all = None
                    ap_max_list_instantaneous_freq_1_2.append(None)

                # find trace with maximal number of APs
                if ap_number > max_ap_number:
                    max_ap_number = ap_number
                    max_ap_voltage_trace = voltage
                    max_ap_di = di
                    ap_max_half_duration_1st = hd_end_t[0] - hd_start_t[0]
                    ap_max_threshold_1st = th_v[0]
                    ap_max_threshold_2nd_1st = th_v_2nd[0]
                    ap_max_amplitude_1st = p_v[0] - th_v[0]
                    # Calculate average AP parameters
                    ap_max_half_duration_av = sum(hd_end_t[i] - hd_start_t[i] for i in range(ap_number)) / ap_number
                    ap_max_threshold_av = sum(th_v) / ap_number
                    ap_max_threshold_2nd_av = sum(th_v_2nd) / ap_number
                    ap_max_amplitude_av = sum(p_v[i] - th_v[i] for i in range(ap_number)) / ap_number

            axs[6].grid(True)

            # Setup max AP trace plot
            axs[7].set_title("Max AP Voltage Trace")
            axs[7].set_ylabel("Voltage (mV)")
            axs[7].set_xlabel("Time (s)")
            if max_ap_voltage_trace is not None:
                axs[7].plot(time, max_ap_voltage_trace, label=f'Max APs: {max_ap_number} at {max_ap_di:.1f} pA')
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
                voltage = V_to_mV * bundle.data[group_id, series_id, sweep_id, 0]
                #current = A_to_pA * bundle.data[group_id, series_id, sweep_id, 1]

                baseline_voltage = voltage[idx1].mean()
                ap_broadening_list_voltage_baseline.append(float(baseline_voltage))  # Convert numpy.float64 to Python float

                ap_number, th_v, th_t, th_v_2nd, th_t_2nd, hd_start_t, hd_start_v, hd_end_t, hd_end_v, p_v, p_t, ahp_v, ahp_t, dvdt_v, dvdt_t = ap_analysis(
                    time, voltage, v_threshold, dvdt_threshold, smooth_window, fraction_of_max_of_2nd_derivative,
                    window_for_searching_threshold, window_for_searching_ahp,
                    minimal_ap_interval, minimal_ap_duration, maximal_ap_duration,
                    maximal_relative_amplitude_decline)

                # Store half duration of first AP (or None if no AP detected)
                if ap_number > 0:
                    # Store analysis points
                    sweep_points = {
                        'threshold': list(zip(th_t, [v / V_to_mV for v in th_v])),
                        'threshold_2nd': list(zip(th_t_2nd, [v / V_to_mV for v in th_v_2nd])),
                        'half_duration_start': list(zip(hd_start_t, [v / V_to_mV for v in hd_start_v])),
                        'half_duration_end': list(zip(hd_end_t, [v / V_to_mV for v in hd_end_v])),
                        'peak': list(zip(p_t, [v / V_to_mV for v in p_v])),
                        'ahp': list(zip(ahp_t, [v / V_to_mV for v in ahp_v])),
                        'dvdt_max': list(zip(dvdt_t, [v / V_to_mV for v in dvdt_v])),
                        'smooth_window': smooth_window
                    }
                    analysis_points[file_name][group_id][series_id][sweep_id][0] = sweep_points

                    ap_broadening_list_half_duration_1st.append(float(hd_end_t[0] - hd_start_t[0]))
                    ap_broadening_list_threshold_1st.append(float(th_v[0]))
                    ap_broadening_list_amplitude_1st.append(float(p_v[0] - th_v[0]))
                else:
                    ap_broadening_list_half_duration_1st.append(None)
                    ap_broadening_list_threshold_1st.append(None)
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

            "ap_rheo_half_duration_1st": ap_rheo_half_duration_1st,
            "ap_rheo_threshold_1st": ap_rheo_threshold_1st,
            "ap_rheo_threshold_2nd_1st": ap_rheo_threshold_2nd_1st,
            "ap_rheo_amplitude_1st": ap_rheo_amplitude_1st,
            
            "ap_rheo_half_duration_av": ap_rheo_half_duration_av,
            "ap_rheo_threshold_av": ap_rheo_threshold_av,
            "ap_rheo_threshold_2nd_av": ap_rheo_threshold_2nd_av,
            "ap_rheo_amplitude_av": ap_rheo_amplitude_av,
            
            "ap_max_half_duration_1st": ap_max_half_duration_1st,
            "ap_max_threshold_1st": ap_max_threshold_1st,
            "ap_max_threshold_2nd_1st": ap_max_threshold_2nd_1st,
            "ap_max_amplitude_1st": ap_max_amplitude_1st,
            
            "ap_max_half_duration_av": ap_max_half_duration_av,
            "ap_max_threshold_av": ap_max_threshold_av,
            "ap_max_threshold_2nd_av": ap_max_threshold_2nd_av,
            "ap_max_amplitude_av": ap_max_amplitude_av,
            
            "ap_max_instantaneous_freq_1_2": ap_max_instantaneous_freq_1_2,
            "ap_max_average_freq_all": ap_max_average_freq_all,

            "ap_max_list_current_steps": ap_max_list_current_steps,
            "ap_max_list_ap_numbers": ap_max_list_ap_numbers,
            "ap_max_list_instantaneous_freq_1_2": ap_max_list_instantaneous_freq_1_2,
            "ap_max_list_half_duration_1st": ap_max_list_half_duration_1st,
            "ap_max_list_threshold_1st": ap_max_list_threshold_1st,
            "ap_max_list_amplitude_1st": ap_max_list_amplitude_1st,

            "ap_broadening_list_voltage_baseline": ap_broadening_list_voltage_baseline,
            "ap_broadening_list_half_duration_1st": ap_broadening_list_half_duration_1st,
            "ap_broadening_list_threshold_1st": ap_broadening_list_threshold_1st,
            "ap_broadening_list_amplitude_1st": ap_broadening_list_amplitude_1st
        })

        # --- Save PDF ---
        plt.tight_layout()
        pdf_filename = f"{cell_count + 1:03d}_{os.path.splitext(file_name)[0]}.pdf"
        plt.savefig(os.path.join(output_folder_traces, pdf_filename))
        plt.close()

    # ==========================================================================================
    # --- Save all results to Excel ---
    # ==========================================================================================
    # results EXCEL file
    results_df = pd.DataFrame(results)
    excel_output_path = os.path.join(output_folder_results, "results.xlsx")
    # Create copy without the last 10 columns for export
    export_df = results_df.iloc[:, :-10]
    export_df.to_excel(excel_output_path, index=False)

    # lists of ap max =========================================
    # Create DataFrames for lists
    ap_max_list_current_steps_data = []
    ap_max_list_ap_numbers_data = []
    ap_max_list_instantaneous_freq_1_2_data = []
    ap_max_list_half_duration_1st_data = []
    ap_max_list_threshold_1st_data = []
    ap_max_list_amplitude_1st_data = []

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

        for i, (tmp1, tmp2 ,tmp3, tmp4, tmp5, tmp6) in enumerate(zip(results[cell_count]['ap_max_list_current_steps'],
                                                                    results[cell_count]['ap_max_list_ap_numbers'],
                                                                    results[cell_count]['ap_max_list_instantaneous_freq_1_2'],
                                                                    results[cell_count]['ap_max_list_half_duration_1st'],
                                                                    results[cell_count]['ap_max_list_threshold_1st'],
                                                                    results[cell_count]['ap_max_list_amplitude_1st']), 1):
            ap_max_list_current_steps_row[f'sweep_{i}'] = tmp1
            ap_max_list_ap_numbers_row[f'sweep_{i}'] = tmp2
            ap_max_list_instantaneous_freq_1_2_row[f'sweep_{i}'] = tmp3
            ap_max_list_half_duration_1st_row[f'sweep_{i}'] = tmp4
            ap_max_list_threshold_1st_row[f'sweep_{i}'] = tmp5
            ap_max_list_amplitude_1st_row[f'sweep_{i}'] = tmp6

        ap_max_list_current_steps_data.append(ap_max_list_current_steps_row)
        ap_max_list_ap_numbers_data.append(ap_max_list_ap_numbers_row)
        ap_max_list_instantaneous_freq_1_2_data.append(ap_max_list_instantaneous_freq_1_2_row)
        ap_max_list_half_duration_1st_data.append(ap_max_list_half_duration_1st_row)
        ap_max_list_threshold_1st_data.append(ap_max_list_threshold_1st_row)
        ap_max_list_amplitude_1st_data.append(ap_max_list_amplitude_1st_row)

    # lists of ap broadening =========================================
    # Create DataFrames for lists
    ap_broadening_list_voltage_baseline_data = []
    ap_broadening_list_half_duration_1st_data = []
    ap_broadening_list_threshold_1st_data = []
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
        ap_broadening_list_amplitude_1st_row = {'cell_count': cell_count + 1, 'file_name': file_name}

        for i, (tmp1, tmp2 ,tmp3, tmp4) in enumerate(zip(results[cell_count]['ap_broadening_list_voltage_baseline'],
                                                                    results[cell_count]['ap_broadening_list_half_duration_1st'],
                                                                    results[cell_count]['ap_broadening_list_threshold_1st'],
                                                                    results[cell_count]['ap_broadening_list_amplitude_1st']), 1):
            ap_broadening_list_voltage_baseline_row[f'sweep_{i}'] = tmp1
            ap_broadening_list_half_duration_1st_row[f'sweep_{i}'] = tmp2
            ap_broadening_list_threshold_1st_row[f'sweep_{i}'] = tmp3
            ap_broadening_list_amplitude_1st_row[f'sweep_{i}'] = tmp4

        ap_broadening_list_voltage_baseline_data.append(ap_broadening_list_voltage_baseline_row)
        ap_broadening_list_half_duration_1st_data.append(ap_broadening_list_half_duration_1st_row)
        ap_broadening_list_threshold_1st_data.append(ap_broadening_list_threshold_1st_row)
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

    ap_broadening_list_amplitude_1st_df = pd.DataFrame(ap_broadening_list_amplitude_1st_data)
    ap_broadening_list_amplitude_1st_excel_path = os.path.join(output_folder_results, "ap_broadening_list_amplitude_1st.xlsx")
    ap_broadening_list_amplitude_1st_df.to_excel(ap_broadening_list_amplitude_1st_excel_path, index=False)

    # save points =========================================
    # Save analysis points to JSON file
    analysis_points_path = os.path.join(IMPORT_FOLDER, "analysis_points.json")
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