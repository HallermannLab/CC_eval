>         d1 = np.gradient(voltage_filt, dt)
> ```
>
> It should be aligned with `voltage_filt = ...`:
>
<llm-snippet-file>main.py</llm-snippet-file>
```python
# ... existing code ...
    # smooth voltages and use numerical derivatives (central differences via np.gradient)
    # NOTE difference with browser.py, where voltage is in V (to show unmodified traces) and dV/dt is in V/s (also used as thershold value) but d2V/dt2 is in V/ms^2 (to remove 1e6 in the number on the axes)
    voltage_filt = savgol_filter(voltage, window_length=win_samples, polyorder=sg_polyorder)
    d1 = np.gradient(voltage_filt, dt)
    d1_filt = savgol_filter(d1, window_length=win_samples, polyorder=sg_polyorder)
    d1_in_V_per_s = d1_filt / V_to_mV
    d2 = np.gradient(d1_in_V_per_s, dt)
# ... existing code ...
```

Now add the new metadata values inside `CC_eval()`, after your existing AP-analysis parameter loading.

<llm-snippet-file>main.py</llm-snippet-file>
```python
# ... existing code ...
        minimal_ap_interval = row['minimal_ap_interval']
        minimal_ap_duration = row['minimal_ap_duration']
        maximal_ap_duration = row['maximal_ap_duration']
        maximal_relative_amplitude_decline = row['maximal_relative_amplitude_decline']

        phase_plot_t1 = row['phase_plot_t1']
        phase_plot_t2 = row['phase_plot_t2']
        phase_plot_range_start = row['phase_plot_range_start']
        phase_plot_range_end = row['phase_plot_range_end']

        dat_path = os.path.join(config.EXTERNAL_DATA_FOLDER, file_name)
        try:
            bundle = heka_reader.Bundle(dat_path)
# ... existing code ...
```

Initialize the new phase-analysis variables and storage objects near the other result variables.

<llm-snippet-file>main.py</llm-snippet-file>
```python
# ... existing code ...
        ap_rheo_peak1_peak2_interval_1st = None
        ap_rheo_maxdvdt_1st = None

        ap_rheo_phase_peak_d1_in_V_per_s = None
        ap_rheo_phase_peak_voltage_filt = None
        ap_rheo_phase_data = None

        ap_max_baseline_voltage = None

        ap_max_half_duration_1st = None
# ... existing code ...
        ap_max_peak1_peak2_interval_av = None
        ap_max_peak1_peak2_interval_median = None

        ap_max_phase_peak_d1_in_V_per_s = None
        ap_max_phase_peak_voltage_filt = None
        ap_max_phase_data = None

        ap_max_list_current_steps = []
        ap_max_list_ap_numbers = []
# ... existing code ...
```

Store the rheobase first-AP phase-analysis data at the moment the rheobase trace is identified.

<llm-snippet-file>main.py</llm-snippet-file>
```python
# ... existing code ...
                if rheobase is None and ap_number > 0:
                    rheobase_voltage_trace = voltage
                    rheobase = di
                    # Calculate baseline voltage from window1 (before stimulus)
                    ap_rheo_baseline_voltage = voltage[idx1].mean()
                    # Calculate delay of first AP (time from stimulus start to first AP threshold)
                    # Assuming stimulus starts at window2_ap_rheo_start
                    ap_rheo_first_ap_delay = th_t[0]
                    
                    ap_rheo_half_duration_1st = hd_end_t[0] - hd_start_t[0]
                    ap_rheo_threshold_1st = th_v[0]
                    ap_rheo_threshold_2nd_1st = th_v_2nd[0]
                    ap_rheo_amplitude_1st = p_v[0] - th_v[0]
                    ap_rheo_peak1_peak2_interval_1st = d2_peak2_t[0] - d2_peak1_t[0]
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
                        ap_rheo_phase_peak_d1_in_V_per_s = ap_rheo_phase_data["phase_peak"]["d1_in_V_per_s"]
                        ap_rheo_phase_peak_voltage_filt = ap_rheo_phase_data["phase_peak"]["voltage_filt"]

                    # Calculate average AP parameters
                    ap_rheo_half_duration_av = sum(hd_end_t[i] - hd_start_t[i] for i in range(ap_number)) / ap_number
                    ap_rheo_threshold_av = sum(th_v) / len(th_v)
                    ap_rheo_threshold_2nd_av = sum(th_v_2nd) / len(th_v_2nd)
                    ap_rheo_amplitude_av = sum(p_v[i] - th_v[i] for i in range(len(th_v))) / len(th_v)
# ... existing code ...
```

Store the max-AP trace first-AP phase-analysis data when the trace with the maximal AP number is updated.

<llm-snippet-file>main.py</llm-snippet-file>
```python
# ... existing code ...
                if ap_number >= max_ap_number and ap_number > 0:
                    max_ap_number = ap_number
                    max_ap_voltage_trace = voltage
                    max_ap_di = di
                    ap_max_baseline_voltage = voltage[idx1].mean() # Calculate baseline voltage from window1 (before stimulus) for max AP trace
                    ap_max_half_duration_1st = hd_end_t[0] - hd_start_t[0]
                    ap_max_threshold_1st = th_v[0]
                    ap_max_threshold_2nd_1st = th_v_2nd[0]
                    ap_max_amplitude_1st = p_v[0] - th_v[0]
                    ap_max_peak1_peak2_interval_1st = d2_peak2_t[0] - d2_peak1_t[0]

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
                        ap_max_phase_peak_d1_in_V_per_s = ap_max_phase_data["phase_peak"]["d1_in_V_per_s"]
                        ap_max_phase_peak_voltage_filt = ap_max_phase_data["phase_peak"]["voltage_filt"]
                    else:
                        ap_max_phase_peak_d1_in_V_per_s = None
                        ap_max_phase_peak_voltage_filt = None

                    # Calculate average AP parameters
                    ap_max_half_duration_av = sum(hd_end_t[i] - hd_start_t[i] for i in range(ap_number)) / ap_number
                    ap_max_threshold_av = sum(th_v) / ap_number
                    ap_max_threshold_2nd_av = sum(th_v_2nd) / ap_number
# ... existing code ...
```

Add the new Excel-output fields to the results dictionary.

<llm-snippet-file>main.py</llm-snippet-file>
```python
# ... existing code ...
            "ap_rheo_amplitude_av": ap_rheo_amplitude_av,
            "ap_rheo_maxdvdt_1st": ap_rheo_maxdvdt_1st,
            "ap_rheo_peak1_peak2_interval_1st": ap_rheo_peak1_peak2_interval_1st,
            "ap_rheo_phase_peak_d1_in_V_per_s": ap_rheo_phase_peak_d1_in_V_per_s,
            "ap_rheo_phase_peak_voltage_filt": ap_rheo_phase_peak_voltage_filt,

            "ap_max_half_duration_1st": ap_max_half_duration_1st,
            "ap_max_threshold_1st": ap_max_threshold_1st,
# ... existing code ...
            "ap_max_instantaneous_freq_1_2": ap_max_instantaneous_freq_1_2,
            "ap_max_instantaneous_freq_last": ap_max_instantaneous_freq_last,
            "ap_max_freq_adaptation": ap_max_freq_adaptation,
            "ap_max_phase_peak_d1_in_V_per_s": ap_max_phase_peak_d1_in_V_per_s,
            "ap_max_phase_peak_voltage_filt": ap_max_phase_peak_voltage_filt,

            "ap_max_list_current_steps": ap_max_list_current_steps,
            "ap_max_list_ap_numbers": ap_max_list_ap_numbers,
# ... existing code ...
```

Finally, replace the single-page `plt.savefig(...)` block with a multipage PDF writer.

<llm-snippet-file>main.py</llm-snippet-file>
```python
# ... existing code ...
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
# ... existing code ...
```

One more small but important change: because you are adding two new scalar result columns, the current `export_df = results_df.iloc[:, :-15]` becomes fragile. It depends on the number of list-like columns at the end. Replace it with an explicit exclusion of list columns.

<llm-snippet-file>main.py</llm-snippet-file>
```python
# ... existing code ...
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
# ... existing code ...
```

This keeps the main `results.xlsx` stable and ensures the new phase-peak values are exported instead of accidentally being chopped off by column position.
