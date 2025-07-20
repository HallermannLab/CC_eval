import os
import sys
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import heka_reader
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets, QtCore

# --- parameters ---
A_to_pA = 1e12
V_to_mV = 1e3

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
dvdt_threshold = 30
window_for_searching_threshold = 0.001
window_for_searching_ahp = 0.001
minimal_ap_interval = 0.001
minimal_ap_duration = 0.0005
maximal_ap_duration = 0.005
maximal_relative_amplitude_decline = 0.7


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
    # Calculate derivative for threshold detection
    dt = time[1] - time[0]
    d1 = np.diff(voltage) / dt

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
        window_d1 = d1[start_idx:end_idx]
        th_idx = start_idx + np.where(window_d1 >= dvdt_threshold)[0][0]

        th_v.append(voltage[th_idx])
        th_t.append(time[th_idx])

        # Find peak
        peak_idx = th_idx + np.argmax(voltage[th_idx:th_idx + int(maximal_ap_duration / dt)])
        p_v.append(voltage[peak_idx])
        p_t.append(time[peak_idx])

        # Find AHP
        ahp_window = voltage[peak_idx:peak_idx + int(window_for_searching_ahp / dt)]
        ahp_idx = peak_idx + np.argmin(ahp_window)
        ahp_v.append(voltage[ahp_idx])
        ahp_t.append(time[ahp_idx])

        # Max dV/dt
        dvdt_idx = th_idx + np.argmax(d1[th_idx:peak_idx])
        dvdt_v.append(d1[dvdt_idx])
        dvdt_t.append(time[dvdt_idx])

        # Calculate half-duration points
        half_amplitude = (voltage[peak_idx] - voltage[th_idx]) / 2 + voltage[th_idx]

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

        hd_start_v.append(voltage[hd_start_idx])
        hd_start_t.append(time[hd_start_idx])
        hd_end_v.append(voltage[hd_end_idx])
        hd_end_t.append(time[hd_end_idx])

    ap_number = len(th_v)

    return ap_number, th_v, th_t, hd_start_t, hd_start_v, hd_end_t, hd_end_v, p_v, p_t, ahp_v, ahp_t, dvdt_v, dvdt_t


def CC_eval():
    # --- Configuration ---
    ROOT_FOLDER = "/Users/stefanhallermann/Library/CloudStorage/Dropbox/tmp/Rinako"
    import_folder = os.path.join(ROOT_FOLDER, "in")
    metadata_file = os.path.join(import_folder, "metadata.xlsx")

    # --- Create Output Folders ---
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output_folder = os.path.join(ROOT_FOLDER, f"output_SH_{timestamp}")
    os.makedirs(output_folder, exist_ok=True)

    output_folder_traces = os.path.join(output_folder, "traces")
    os.makedirs(output_folder_traces, exist_ok=True)

    # --- Load Metadata ---
    metadata_df = pd.read_excel(metadata_file)

    # --- Results Storage ---
    results = []


    # --- Process Each Cell ---
    for cell_count, row in metadata_df.iterrows():
        file_name = row['file_name']
        rmp_series = row['rmp_series']
        rin_series = row['rin_series']
        ap_rheo_series = row['ap_rheo_series']
        ap_max_series = row['ap_max_series']

        dat_path = os.path.join(import_folder, file_name)
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
        for sweep_id in range(n_sweeps):
            voltage = V_to_mV * bundle.data[group_id, series_id, sweep_id, 0]
            current = A_to_pA * bundle.data[group_id, series_id, sweep_id, 1]

            # Plot voltage trace in superposition plot
            axs[4].plot(time, voltage, alpha=0.5, label=f'Sweep {sweep_id + 1}')

            di = current[idx2].mean() - current[idx1].mean()

            ap_number, th_v, th_t, hd_start_t, hd_start_v, hd_end_t, hd_end_v, p_v, p_t, ahp_v, ahp_t, dvdt_v, dvdt_t = ap_analysis(
                time, voltage)
            #print(f"Sweep {sweep_id + 1}: ap_number = {ap_number}")

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

if __name__ == '__main__':
    CC_eval()