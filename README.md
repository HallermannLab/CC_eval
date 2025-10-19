# Electrophysiology Analysis Tool

A Python-based tool for analyzing electrophysiological recordings, specifically designed for processing and analyzing current-clamp recordings from neurons. 

## Features

- Resting Membrane Potential (RMP) analysis
- Input Resistance (Rin) calculation
- Rheobase determination
- Action Potential (AP) detection and analysis including low-pass filtering of voltage trace:
  - AP threshold detection (based on both a fixed dV/dt value of the 1st derivative and dual peak detection of the 2nd derivative)
  - AP amplitude calculation
  - AP half-width measurement
  - Maximum firing rate analysis
  - Peak1-Peak2 interval analysis in 2nd derivative for AP initiation kinetics
- Action Potential Broadening Analysis: Tracks changes in AP half-duration across multiple sweeps to detect activity-dependent broadening
- Baseline voltage tracking across different current injection levels
- Individual cell parameter customization for optimized analysis
- Automated batch processing of multiple recordings
- Comprehensive Excel export with detailed parameter lists
- Interactive browser interface for visualizing action potential parameters and 1st and 2nd derivatives

## Installation

1. Clone this repository
2. Install required packages:
   - See `requirements.txt`

## Configuration

### Data Folder Setup
The tool supports external data folders for flexible file organization. Configure the paths in `config.py`:

### Metadata File Requirements
Create a metadata Excel file (`metadata.xlsx`) in your import folder with all of the following columns. Each cell requires individual parameter values for optimal analysis:

#### Example Metadata File Structure (see provided `example metadata.xlsx`):
| file_name | rmp_series | rin_series | ap_rheo_series | ap_max_series | ap_broadening_series | rin_series_from_current | rin_series_to_current | window1_rin_start | window1_rin_end | window2_rin_start | window2_rin_end | window1_ap_rheo_start | window1_ap_rheo_end | window2_ap_rheo_start | window2_ap_rheo_end | window1_ap_max_start | window1_ap_max_end | window2_ap_max_start | window2_ap_max_end | v_threshold | dvdt_threshold | smooth_window | fraction_of_max_of_2nd_derivative | window_for_searching_threshold | window_for_searching_ahp | minimal_ap_interval | minimal_ap_duration | maximal_ap_duration | maximal_relative_amplitude_decline |
|-----------|------------|------------|----------------|---------------|---------------------|------------------------|----------------------|-------------------|-----------------|-------------------|------------------|----------------------|---------------------|----------------------|---------------------|----------------------|--------------------|--------------------|--------------------|-----------|--------------|--------------|---------------------------------|------------------------------|-------------------------|--------------------|--------------------|--------------------|---------------------------------|
| cell_001.dat | 1 | 2 | 3 | 4 | 5 | -50 | 50 | 0.01 | 0.09 | 0.30 | 0.39 | 0.01 | 0.09 | 0.30 | 0.39 | 0.01 | 0.09 | 0.30 | 0.39 | -25 | 15 | 1500 | 0.3 | 0.002 | 0.01 | 0.001 | 0.0005 | 0.005 | 0.8 |
| cell_002.dat | 1 | 2 | 3 | 4 | 5 | -50 | 50 | 0.01 | 0.09 | 0.30 | 0.39 | 0.01 | 0.09 | 0.30 | 0.39 | 0.01 | 0.09 | 0.30 | 0.39 | -30 | 12 | 2000 | 0.3 | 0.002 | 0.01 | 0.001 | 0.0007 | 0.005 | 0.8 |
| cell_003.dat | 1 | 2 | 3 | 4 | 5 | -50 | 50 | 0.01 | 0.09 | 0.30 | 0.39 | 0.01 | 0.09 | 0.30 | 0.39 | 0.01 | 0.09 | 0.30 | 0.39 | -20 | 18 | 1200 | 0.3 | 0.002 | 0.01 | 0.001 | 0.0004 | 0.005 | 0.8 |


#### Columns Explained:

**Basic Recording Information:**
- `file_name`: Name of the HEKA data file
- `rmp_series`: Series number for RMP analysis
- `rin_series`: Series number for input resistance analysis
- `ap_rheo_series`: Series number for rheobase analysis
- `ap_max_series`: Series number for maximum AP analysis
- `ap_broadening_series`: Series number for AP broadening analysis (tracks half-duration of first AP in each sweep)

You can skip analysis for a given series by leaving the cell empty in the metadata excel file. 

**Input Resistance Analysis:**
- `rin_series_from_current`: Lower current bound for Rin analysis (pA)
- `rin_series_to_current`: Upper current bound for Rin analysis (pA)
- `window1_rin_start/end`: First time window for Rin analysis (s)
- `window2_rin_start/end`: Second time window for Rin analysis (s)

**Action Potential Analysis Windows:**
- `window1_ap_rheo_start/end`: Time windows for rheobase analysis (s)
- `window2_ap_rheo_start/end`
- `window1_ap_max_start/end`: Time windows for maximum AP analysis (s)
- `window2_ap_max_start/end`

**Action Potential Detection Parameters:**
- `v_threshold`: Voltage threshold for AP detection (mV, e.g., -25)
- `dvdt_threshold`: dV/dt threshold for AP detection (V/s, e.g., 15)
- `smooth_window`: Smoothing window for filtering using Savitzky-Golay filter (s, e.g. 0.1)
- `fraction_of_max_of_2nd_derivative`: Fraction of max 2nd derivative for threshold detection (e.g., 0.3)
- `window_for_searching_threshold`: Time window for finding threshold (s, e.g., 0.002)
- `window_for_searching_ahp`: Time window for finding AHP (s, e.g., 0.01)
- `minimal_ap_interval`: Minimum time between APs (s, e.g., 0.001)
- `minimal_ap_duration`: Minimum AP duration (s, e.g., 0.0005)
- `maximal_ap_duration`: Maximum AP duration (s, e.g., 0.005)
- `maximal_relative_amplitude_decline`: Maximum allowed amplitude decline (e.g., 0.8)

# Usage

1. **Prepare your data:**
   - Place your HEKA data files in the external data folder (as configured in `config.py`)
   - Create the metadata Excel file with all required columns and individual parameters for each cell (see `example metadata.xlsx`)

3. **Run the analysis:**
   - `python main.py`
   - Performs full analysis and automatically launches the browser interface

4. **Browser interface:**
   - Interactive visualization of action potential parameters
   - Zoom and navigation controls
   - Select files from a dropdown menu based on metadata
   - Visualize detected analysis points overlaid on traces
   - Can be launched separately: `python browser.py`

## Output

The analysis generates:

- **Timestamped output folder** containing:
  - Individual PDF files with analysis plots for each recording ([see example here PDF](003_Rinako_20250709_003.pdf))
  - Main Excel file (`results.xlsx`) with compiled analysis results including:
    - Basic parameters (RMP, Rin, Rheobase)
    - AP properties for rheobase and maximum current traces including baseline voltages, AP delays, and peak1-Peak2 intervals from 2nd derivative analysis
  - Detailed Excel files with sweep-by-sweep data:
    - `ap_max_list_current_steps.xlsx` - Current injection levels
    - `ap_max_list_ap_numbers.xlsx` - Number of APs per sweep
    - `ap_max_list_instantaneous_freq_1_2.xlsx` - Instantaneous frequency between first two APs
    - `ap_max_list_half_duration_1st.xlsx` - Half-duration of first AP
    - `ap_max_list_threshold_1st.xlsx` - Threshold voltage of first AP
    - `ap_max_list_amplitude_1st.xlsx` - Amplitude of first AP
    - `ap_max_list_voltage_baseline.xlsx` - Baseline voltage for each sweep
    - `ap_max_list_average_frequency.xlsx` - Average firing frequency for each sweep
    - `ap_broadening_list_*` files for AP broadening analysis
  - Copy of `analysis_points.json` for documentation

- **`analysis_points.json`** in the input folder:
  - Contains all detected action potential parameters for each trace
  - Used by browser interface for visualization
  - Overwritten with each analysis run

## Analysis Results Details

### Main Results Excel File
The `results.xlsx` file contains the following key parameters:

**Basic Parameters:**
- `RMP_mean`, `RMP_min`: Resting membrane potential statistics
- `Rin_fit`: Input resistance from linear fit
- `Rheobase`: Minimum current to elicit an AP
- `max_ap_number`: Maximum number of APs fired in any sweep

**Rheobase Analysis:**
- `ap_rheo_*_1st`: Properties of the first AP at rheobase
- `ap_rheo_*_av`: Average properties of all APs at rheobase
- `ap_rheo_baseline_voltage`: Baseline voltage at rheobase trace
- `ap_rheo_first_ap_delay`: Delay from stimulus start to first AP

**Maximum AP Analysis:**
- `ap_max_*_1st`: Properties of the first AP in trace with most APs
- `ap_max_*_av`: Average properties of all APs in trace with most APs
- `ap_max_baseline_voltage`: Baseline voltage in trace with most APs

## Browser Interface Features

The interactive browser provides:

- **File selection**: Choose from files listed in metadata via dropdown menu
- **Interactive navigation through the HEKA traces tree**: Zoom, pan, and explore traces with overlaid analysis points
- **Analysis point visualization**: View detected thresholds, peaks, AHP, and other parameters with color-coded markers:
  - **Red circles**: Threshold (1st derivative)
  - **Brown circles**: Threshold (2nd derivative)
  - **Blue squares**: Half-duration start
  - **Green squares**: Half-duration end
  - **Yellow triangles**: Peak
  - **Purple diamonds**: AHP
  - **Cyan pentagons**: Maximum dV/dt
- **Derivative plotting**: Automatic display of 1st and 2nd derivatives when analysis points are available
  - **Red circles**: Threshold shown on 1st derivative
  - **Brown circles**: Threshold of crossing 'fraction_of_max_of_2nd_derivative' shown on 2nd derivative
  - **Purple triangles**: Peak1 shown on 2nd derivative
  - **Yellow diamonds**: Peak2 shown on 2nd derivative
- **Interactive navigation**: Zoom, pan, and explore traces with overlaid analysis points

## Advanced Features

### Dual Peak Detection in 2nd Derivative
The tool now implements peak detection in the 2nd derivative of voltage during AP initiation:
- Uses `scipy.signal.find_peaks` for peak identification
- Identifies Peak1 (left) and Peak2 (right) based on temporal position
- Calculates intervals between peaks for kinetic analysis
- Provides fallback mechanisms for edge cases

## Acknowledgements

I am grateful for using the HEKA reader and the trace browser provided by: [https://github.com/campagnola/heka_reader](https://github.com/campagnola/heka_reader)

The trace browser has been modified to:
- Allow file selection from metadata
- Superimpose analysis points on action potential traces
- Display 1st and 2nd derivative plots with detected peaks
- Provide interactive visualization of analysis results with color-coded parameter markers