# Electrophysiology Analysis Tool

A Python-based tool for analyzing electrophysiological recordings, specifically designed for processing and analyzing current-clamp recordings from neurons. This tool provides comprehensive analysis of membrane properties and action potential characteristics.

## Features

- Resting Membrane Potential (RMP) analysis
- Input Resistance (Rin) calculation
- Rheobase determination
- Action Potential (AP) detection and analysis including low-pass filtering of voltage trace:
  - AP threshold detection (based on both a fixed dV/dt value of the 1st derivative and the maximum of the 2nd derivative)
  - AP amplitude calculation
  - AP half-width measurement
  - Maximum firing rate analysis
- Automated batch processing of multiple recordings
- Generation of detailed visualization plots
- Export of results to Excel
- Visualization of action potential parameters for individual action potentials in a browser

## Installation

1. Clone this repository
2. Install required packages:
   - See `requirements.txt`.

## Usage

1. Prepare your data:
   - Place your HEKA data files in the input folder.
   - Create a metadata Excel file (`metadata.xlsx`) with the following columns:
     - file_name
     - rmp_series
     - rin_series
     - ap_rheo_series
     - ap_max_series

2. Configure parameters in the script:
   - Adjust analysis windows if needed.
   - Modify AP detection parameters if required.

3. Run the analysis:
   - Execute `main.py` to perform the full analysis and automatically launch the browser interface.

4. Browser interface:
   - After running `main.py`, the browser interface will automatically start, in which various options for zoom in and out exist.
   - The etermined action potential parameters can be visualized for each action potential in an interactive browser interface.
   - The browser can also be launched separately by running `browser.py` directly (useful for reviewing results without re-running the analysis).

## Output

The script generates:
- A timestamped output folder containing:
  - Individual PDF files with analysis plots for each recording.
  - An Excel file (`results.xlsx`) with compiled analysis results.
  - A copy of `analysis_points.json` for documentation purposes.
- An `analysis_points.json` file in the input folder containing all detected action potential parameters (this file is overwritten with each analysis run).

## Acknowledgement

I am greatful for using the HEKA reader and the trace browser provided by: [https://github.com/campagnola/heka_reader](https://github.com/campagnola/heka_reader)

I have modified the trace browser for allowing the selection of the files from the metadata and for superposition of the analysis point of the action potentials.

## Individual Cell Parameters (individual-parameters branch)

The software supports two parameter configuration modes:

### Default Mode (main branch)
Uses uniform parameters defined globally for all cells.

### Individual Parameters Mode (individual-parameters branch)

The `individual-parameters` branch allows **cell-specific parameter customization**, which is particularly useful when:

- Different cell types require different analysis parameters
- Individual recordings have varying signal characteristics  
- You need to optimize parameters for specific experimental conditions
- Recording quality varies between cells

#### How to Use Individual Parameters

1. **Switch to the individual-parameters branch:**
   ```bash
   git checkout individual-parameters
   ```

2. **Update your metadata file:**
   The metadata Excel file should include individual parameter columns. Each row can have its own parameter values:

   | file_name | v_threshold | dvdt_threshold | filter_cut_off | minimal_ap_duration | window1_rin_start | ... |
   |-----------|-------------|----------------|----------------|---------------------|-------------------|-----|
   | cell_001.dat | -25 | 15 | 1500 | 0.0005 | 0.01 | ... |
   | cell_002.dat | -30 | 12 | 2000 | 0.0007 | 0.01 | ... |
   | cell_003.dat | -20 | 18 | 1200 | 0.0004 | 0.01 | ... |

#### Supported Individual Parameters

**Action Potential Detection:**
- `v_threshold`: Voltage threshold for AP detection (mV)
- `dvdt_threshold`: dV/dt threshold for AP detection (V/s)  
- `filter_cut_off`: Cutoff frequency for low-pass filter (Hz)
- `window_for_searching_threshold`: Time window for finding threshold (s)
- `window_for_searching_ahp`: Time window for finding AHP (s)
- `minimal_ap_interval`: Minimum time between APs (s)
- `minimal_ap_duration`: Minimum AP duration (s)
- `maximal_ap_duration`: Maximum AP duration (s)
- `maximal_relative_amplitude_decline`: Maximum allowed amplitude decline

**Analysis Windows:**
- `window1_rin_start/end`: Time windows for input resistance analysis
- `window2_rin_start/end`
- `window1_ap_rheo_start/end`: Time windows for rheobase analysis
- `window2_ap_rheo_start/end` 
- `window1_ap_max_start/end`: Time windows for maximum AP analysis
- `window2_ap_max_start/end`

