# Electrophysiology Analysis Tool

A Python-based tool for analyzing electrophysiological recordings, specifically designed for processing and analyzing current-clamp recordings from neurons. This tool provides comprehensive analysis of membrane properties and action potential characteristics.

## Features

- Resting Membrane Potential (RMP) analysis
- Input Resistance (Rin) calculation
- Rheobase determination
- Action Potential (AP) detection and analysis including low-pass filtering of voltage trace:
  - AP threshold detection
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