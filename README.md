# Electrophysiology Analysis Tool

A Python-based tool for analyzing electrophysiological recordings, specifically designed for processing and analyzing current-clamp recordings from neurons. 

## Features

- Resting Membrane Potential (RMP) analysis
- Input Resistance (Rin) calculation
- Rheobase determination
- Action Potential (AP) detection and analysis including low-pass filtering of voltage trace:
  - AP threshold detection (based on both a fixed dV/dt value of the 1st derivative and the maximum of the 2nd derivative)
  - AP amplitude calculation
  - AP half-width measurement
  - Maximum firing rate analysis
- Action Potential Broadening Analysis: Tracks changes in AP half-duration across multiple sweeps to detect activity-dependent broadening
- Individual cell parameter customization for optimized analysis
- Automated batch processing of multiple recordings
- Generation of detailed visualization plots with optional derivative plotting
- Export of results to Excel
- Interactive browser interface for visualizing action potential parameters

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
| file_name | rmp_series | rin_series | ap_rheo_series | ap_max_series | ap_broadening_series | rin_series_from_current | rin_series_to_current | window1_rin_start | window1_rin_end | window2_rin_start | window2_rin_end | window1_ap_rheo_start | window1_ap_rheo_end | window2_ap_rheo_start | window2_ap_rheo_end | window1_ap_max_start | window1_ap_max_end | window2_ap_max_start | window2_ap_max_end | v_threshold | dvdt_threshold | filter_cut_off | fraction_of_max_of_2nd_derivative | window_for_searching_threshold | window_for_searching_ahp | minimal_ap_interval | minimal_ap_duration | maximal_ap_duration | maximal_relative_amplitude_decline |
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
- `filter_cut_off`: Cutoff frequency for low-pass filter (Hz, e.g., 1500)
- `fraction_of_max_of_2nd_derivative`: Fraction of max 2nd derivative for threshold detection (e.g., 0.3)
- `window_for_searching_threshold`: Time window for finding threshold (s, e.g., 0.002)
- `window_for_searching_ahp`: Time window for finding AHP (s, e.g., 0.01)
- `minimal_ap_interval`: Minimum time between APs (s, e.g., 0.001)
- `minimal_ap_duration`: Minimum AP duration (s, e.g., 0.0005)
- `maximal_ap_duration`: Maximum AP duration (s, e.g., 0.005)
- `maximal_relative_amplitude_decline`: Maximum allowed amplitude decline (e.g., 0.8)

### Plotting Options for Debugging

The tool includes visualization options for derivatives to help with parameter optimization:

- **Derivative Plotting**: Set `which_sweep_to_plot_derivatives` in `main.py` to visualize the 1st and 2nd derivatives
  - Set to `-1` (default): No derivative plotting
  - Set to sweep number (e.g., `0`, `1`, `2`): Plot derivatives for that specific sweep
  - **Purpose**: Shows the 1st and 2nd derivatives overlaid with the voltage trace to help understand and optimize threshold detection parameters
  - **Use case**: Essential for debugging AP detection and fine-tuning `v_threshold`, `dvdt_threshold`, and `fraction_of_max_of_2nd_derivative` parameters

## Usage

1. **Prepare your data:**
   - Place your HEKA data files in the external data folder (as configured in `config.py`)
   - Create the metadata Excel file with all required columns and individual parameters for each cell (see `example metadata.xlsx`)

2. **Configure analysis parameters:**
   - Set individual cell parameters in the metadata file for cell-specific optimization
   - Set `which_sweep_to_plot_derivatives` in `main.py` if you want to visualize derivatives for debugging

3. **Run the analysis:**
   ```bash
   python main.py
   ```
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
  - Individual PDF files with analysis plots for each recording ([see example here PDF](003_Rinako_20250709_003.pdf)
)
  - Excel file (`results.xlsx`) with compiled analysis results
  - Additional Excel files (`ap_max_currents.xlsx`, `ap_max_ap_numbers.xlsx`, `ap_broadening.xlsx`) with current, AP number and AP duration data
  - Copy of `analysis_points.json` for documentation

- **`analysis_points.json`** in the input folder:
  - Contains all detected action potential parameters for each trace
  - Used by browser interface for visualization
  - Overwritten with each analysis run

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
- **Interactive navigation**: Zoom, pan, and explore traces with overlaid analysis points

## Acknowledgements

I am grateful for using the HEKA reader and the trace browser provided by: [https://github.com/campagnola/heka_reader](https://github.com/campagnola/heka_reader)

The trace browser has been modified to:
- Allow file selection from metadata
- Superimpose analysis points on action potential traces
- Provide interactive visualization of analysis results with color-coded parameter markers