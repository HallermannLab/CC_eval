import numpy as np
from scipy.signal import savgol_filter, find_peaks

def calculate_filtered_trace_and_derivatives(time, voltage, smooth_window):
    """Calculate filtered voltage plus first and second derivatives for plotting/phase analysis.

    voltage is expected in V.
    voltage_filt is in V.
    d1 and d1_filt is in V/s = mV/ms
    d2 and d2_filt is V/s².
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
    d1_filt = savgol_filter(d1, window_length=win_samples, polyorder=sg_polyorder)

    d2 = np.gradient(d1_filt, dt)
    d2_filt = savgol_filter(d2, window_length=win_samples, polyorder=sg_polyorder)

    return voltage_filt, d1, d1_filt, d2, d2_filt

