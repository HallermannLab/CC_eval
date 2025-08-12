import os
import sys
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy import signal
import heka_reader
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtCore
import git_save as myGit
from config import ROOT_FOLDER, IMPORT_FOLDER, METADATA_FILE
from config import analysis_points
from collections import defaultdict
import json


# --- parameters ---
A_to_pA = 1e12
V_to_mV = 1e3

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
filter_cut_off = 3000 #Hz
window_for_searching_threshold = 0.002
window_for_searching_ahp = 0.01
minimal_ap_interval = 0.001
minimal_ap_duration = 0.0005
maximal_ap_duration = 0.005
maximal_relative_amplitude_decline = 0.3


def ap_analysis(time, voltage):
    """
    Analyze action potential features from voltage trace

    Args:
        time: time points array (seconds)
        voltage: voltage trace array (mV)

    Returns:
        ap_number: number of APs detected
        th_v: threshold voltages
        th_t: threshold times 
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
    # Apply low-pass filter to voltage trace
    dt = time[1] - time[0]
    nyquist = 0.5 / dt
    cutoff_normalized = filter_cut_off / nyquist
    b, a = signal.butter(2, cutoff_normalized, 'low')
    voltage = signal.filtfilt(b, a, voltage) # CAVE overwrites the trace with the filtered one

    # Calculate derivative using filtered trace
    d1 = np.diff(voltage) / dt
    d1_in_V_per_s = d1 / V_to_mV
    


    # Find threshold crossings where voltage crosses v_threshold
    threshold_idx = np.where((voltage[:-1] < v_threshold) & (voltage[1:] >= v_threshold))[0]

    # Initialize results
    th_v = []
    th_t = []
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
        th_idx = start_idx + th_idx_candidates[0]

        # Find peak
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

        # Find first crossing before peak
        for i in range(peak_idx, th_idx, -1):
            if voltage[i] <= half_amplitude:
                hd_start_idx = i
                break

        # Find first crossing after peak
        for i in range(peak_idx, ahp_idx):
            if voltage[i] <= half_amplitude:
                hd_end_idx = i
                break

        # Validate AP before appending
        half_duration = time[hd_end_idx] - time[hd_start_idx]
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
            p_v.append(voltage[peak_idx])
            p_t.append(time[peak_idx])
            ahp_v.append(voltage[ahp_idx])
            ahp_t.append(time[ahp_idx])
            dvdt_v.append(voltage[dvdt_idx])
            dvdt_t.append(time[dvdt_idx])
            hd_start_v.append(voltage[hd_start_idx])
            hd_start_t.append(time[hd_start_idx])
            hd_end_v.append(voltage[hd_end_idx])
            hd_end_t.append(time[hd_end_idx])

    ap_number = len(p_v)

    return ap_number, th_v, th_t, hd_start_t, hd_start_v, hd_end_t, hd_end_v, p_v, p_t, ahp_v, ahp_t, dvdt_v, dvdt_t


def CC_eval():

    # --- Create Output Folders ---
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output_folder = os.path.join(ROOT_FOLDER, f"output_SH_{timestamp}")
    os.makedirs(output_folder, exist_ok=True)

    output_folder_traces = os.path.join(output_folder, "traces")
    os.makedirs(output_folder_traces, exist_ok=True)
    output_folder_used_data_and_code = os.path.join(output_folder, "used_data_and_code")
    os.makedirs(output_folder_used_data_and_code, exist_ok=True)

    # --- Load Metadata ---
    metadata_df = pd.read_excel(METADATA_FILE)
    # save used data
    metadata_df.to_excel(os.path.join(output_folder_used_data_and_code, "my_data.xlsx"), index=False)

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
        rmp_series = row['rmp_series']
        rin_series = row['rin_series']
        ap_rheo_series = row['ap_rheo_series']
        ap_max_series = row['ap_max_series']

        dat_path = os.path.join(IMPORT_FOLDER, file_name)
        try:
            bundle = heka_reader.Bundle(dat_path)
        except Exception as e:
            print(f"Error reading {file_name}: {e}")
            continue

        group_id = 0
        fig, axs = plt.subplots(4, 2, figsize=(8, 10))
        axs = axs.flatten()

        rmp_mean = None
        rmp_min = None
        rin_fit = None

        # ==========================================================================================
        # --- RMP ---
        # ==========================================================================================
        try:
            series_id = int(rmp_series)

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
        except Exception as e:
            print(f"RMP error in {file_name}: {e}")

        # ==========================================================================================
        # --- Rin ---
        # ==========================================================================================
        try:
            series_id = int(rin_series)
            n_sweeps = bundle.pul[group_id][series_id].NumberSweeps

            delta_vs = []
            delta_is = []

            #get time base
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

            for sweep_id in range(n_sweeps):
                voltage = V_to_mV * bundle.data[group_id, series_id, sweep_id, 0]
                current = A_to_pA * bundle.data[group_id, series_id, sweep_id, 1]

                # Plot current trace
                axs[2].plot(time, current, alpha=0.5)

                # Plot voltage trace
                axs[3].plot(time, voltage, alpha=0.5)

                dv = voltage[idx2].mean() - voltage[idx1].mean()
                di = current[idx2].mean() - current[idx1].mean()

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

            # Linear fit: ΔV = R * ΔI + offset
            slope, intercept, r_value, p_value, std_err = linregress(delta_is, delta_vs)
            rin_fit = 1000 * slope  # input resistance in MOhm (1e6 = 1000 * milli / pico)

            # Plot I-V relationship
            axs[1].scatter(delta_is, delta_vs, label="data")
            axs[1].plot(delta_is, slope * delta_is + intercept, color="red",
                        label=f"fit (R={rin_fit:.1f} MΩ)")
            axs[1].set_xlabel("ΔI (pA)")
            axs[1].set_ylabel("ΔV (mV)")
            axs[1].set_title("Input Resistance")
            axs[1].legend()
            axs[1].grid(True)

        except Exception as e:
            print(f"Rin error in {file_name}: {e}")

        # ==========================================================================================
        # --- AP rheo ---
        # ==========================================================================================
        series_id = int(ap_rheo_series)
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

        rheobase = None
        sweep_points = None  # Default value when no APs are detected

        for sweep_id in range(n_sweeps):
            voltage = V_to_mV * bundle.data[group_id, series_id, sweep_id, 0]
            current = A_to_pA * bundle.data[group_id, series_id, sweep_id, 1]

            # Plot voltage trace in superposition plot
            axs[4].plot(time, voltage, alpha=0.5, label=f'Sweep {sweep_id + 1}')

            di = current[idx2].mean() - current[idx1].mean()

            ap_number, th_v, th_t, hd_start_t, hd_start_v, hd_end_t, hd_end_v, p_v, p_t, ahp_v, ahp_t, dvdt_v, dvdt_t = ap_analysis(
                time, voltage)
            #print(f"Sweep {sweep_id + 1}: ap_number = {ap_number}")
            if ap_number > 0:
                sweep_points = {
                    'threshold': list(zip(th_t, [v / V_to_mV for v in th_v])),
                    'half_duration_start': list(zip(hd_start_t, [v / V_to_mV for v in hd_start_v])),
                    'half_duration_end': list(zip(hd_end_t, [v / V_to_mV for v in hd_end_v])),
                    'peak': list(zip(p_t, [v / V_to_mV for v in p_v])),
                    'ahp': list(zip(ahp_t, [v / V_to_mV for v in ahp_v])),
                    'dvdt_max': list(zip(dvdt_t, [v / V_to_mV for v in dvdt_v]))
                }
                # Store the analysis points in the nested dictionary
                analysis_points[file_name][group_id][series_id][sweep_id][0] = sweep_points
                # Log the structure of analysis points:
                #print(f"Storing points: file_name={file_name}, group_id={group_id}, series_id={series_id}, sweep_id={sweep_id}, trace_id={0}")
                #print(f"Points to store: {json.dumps(sweep_points, indent=2)}")

            if rheobase is None and ap_number > 0:
                rheobase = di
                rheobase_voltage = voltage
                ap_rheo_half_duration = hd_end_t[0] - hd_start_t[0]
                ap_rheo_threshold = th_v[0]
                ap_rheo_amplitude = p_v[0] - th_v[0]

        # Formatting for AP rheo superposition plot
        axs[4].grid(True)

        # Setup rheobase trace plot
        axs[5].set_title("Rheobase Voltage Trace")
        axs[5].set_ylabel("Voltage (mV)")
        axs[5].set_xlabel("Time (s)")
        axs[5].plot(time, rheobase_voltage, label=f'Rheobase: {rheobase:.1f} pA')
        axs[5].legend()
        axs[5].grid(True)

        # ==========================================================================================
        # --- AP max ---
        # ==========================================================================================
        series_id = int(ap_max_series)
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
        max_ap_voltage = None

        for sweep_id in range(n_sweeps):
            voltage = V_to_mV * bundle.data[group_id, series_id, sweep_id, 0]
            current = A_to_pA * bundle.data[group_id, series_id, sweep_id, 1]

            # Plot voltage trace in superposition plot
            axs[6].plot(time, voltage, alpha=0.5, label=f'Sweep {sweep_id + 1}')

            di = current[idx2].mean() - current[idx1].mean()

            ap_number, th_v, th_t, hd_start_t, hd_start_v, hd_end_t, hd_end_v, p_v, p_t, ahp_v, ahp_t, dvdt_v, dvdt_t = ap_analysis(
                time, voltage)
            #print(f"Sweep {sweep_id + 1}: ap_number = {ap_number}")
            if ap_number > 0:
                sweep_points = {
                    'threshold': list(zip(th_t, [v / V_to_mV for v in th_v])),
                    'half_duration_start': list(zip(hd_start_t, [v / V_to_mV for v in hd_start_v])),
                    'half_duration_end': list(zip(hd_end_t, [v / V_to_mV for v in hd_end_v])),
                    'peak': list(zip(p_t, [v / V_to_mV for v in p_v])),
                    'ahp': list(zip(ahp_t, [v / V_to_mV for v in ahp_v])),
                    'dvdt_max': list(zip(dvdt_t, [v / V_to_mV for v in dvdt_v]))
                }
                # Store the analysis points in the nested dictionary
                analysis_points[file_name][group_id][series_id][sweep_id][0] = sweep_points
                # Log the structure of analysis points:
                #print(f"Storing points: file_name={file_name}, group_id={group_id}, series_id={series_id}, sweep_id={sweep_id}, trace_id={0}")
                #print(f"Points to store: {json.dumps(sweep_points, indent=2)}")

            if ap_number > max_ap_number:
                max_ap_number = ap_number
                max_ap_di = di
                max_ap_voltage = voltage
                # Calculate average AP parameters
                ap_max_half_duration = sum(hd_end_t[i] - hd_start_t[i] for i in range(ap_number)) / ap_number
                ap_max_threshold = sum(th_v) / ap_number
                ap_max_amplitude = sum(p_v[i] - th_v[i] for i in range(ap_number)) / ap_number

        # Add vertical lines and formatting for AP max superposition plot
        axs[6].grid(True)

        # Setup max AP trace plot
        axs[7].set_title("Max AP Voltage Trace")
        axs[7].set_ylabel("Voltage (mV)")
        axs[7].set_xlabel("Time (s)")
        if max_ap_voltage is not None:
            axs[7].plot(time, max_ap_voltage, label=f'Max APs: {max_ap_number} at {max_ap_di:.1f} pA')
            axs[7].legend()
        axs[7].grid(True)


        # Turn off unused subplots
        # for i in range(4, 8):
        #    axs[i].set_visible(False)

        # --- Store results ---
        results.append({
            "cell_count": cell_count + 1,
            "file_name": file_name,
            "RMP_mean": rmp_mean,
            "RMP_min": rmp_min,
            "Rin_fit": rin_fit,
            "Rheobase": rheobase,
            "max_ap_number": max_ap_number,
            "ap_rheo_half_duration": ap_rheo_half_duration,
            "ap_rheo_threshold": ap_rheo_threshold,
            "ap_rheo_amplitude": ap_rheo_amplitude,
            "ap_max_half_duration": ap_max_half_duration if max_ap_number > 0 else None,
            "ap_max_threshold": ap_max_threshold if max_ap_number > 0 else None,
            "ap_max_amplitude": ap_max_amplitude if max_ap_number > 0 else None,
        })



        # --- Save PDF ---
        plt.tight_layout()
        pdf_filename = f"{cell_count + 1:03d}_{os.path.splitext(file_name)[0]}.pdf"
        plt.savefig(os.path.join(output_folder_traces, pdf_filename))
        plt.close()

    # --- Save Results to Excel ---
    results_df = pd.DataFrame(results)
    excel_output_path = os.path.join(output_folder, "results.xlsx")
    results_df.to_excel(excel_output_path, index=False)

    print(f"Analysis complete. Results saved to {excel_output_path}")

def start_browser():
    # Import and start the browser
    from browser import app, win
    win.show()
    if sys.flags.interactive == 0:
        app.exec_()

if __name__ == '__main__':
    CC_eval()
    start_browser()  # This will start the browser after CC_eval completes


