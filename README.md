# ASHA: Automatic Single-molecule High content screening Analysis

![Python](https://img.shields.io/badge/Python-3.8%2B-blue)
![R](https://img.shields.io/badge/R-4.0%2B-blue)
![Jupyter](https://img.shields.io/badge/Jupyter-Interactive-orange)

A comprehensive, interactive pipeline for processing, analyzing, and visualizing High-Content Screening Single-Molecule Localization Microscopy (HCS-SMLM) data. 

This toolkit bridges automated data extraction from PoCA and PALMTracer outputs with powerful visualization modules, allowing researchers to evaluate fluorescent protein performance across 96-well plates.

## Key Features

* **Automated Data Parsing:** Seamlessly recursively scan and read `.txt` and `.csv` outputs from PoCA and PALMTracer (`locPALMTracer_merged.txt`, `statMIA.txt`).
* **Photophysical Parameter Calculation:** Automatically compute critical metrics such as `intensity`, `photon_loc`, `blinks`, `total ON`, `avg_on`, and `avg_off`, handling outliers and infinities robustly.
* **Interactive UI (Jupyter Widgets):** A user-friendly, widget-based graphical interface inside Jupyter Notebooks for zero-code data exploration.
* **Multi-Level Visualizations:**
    * **Plate Heatmaps:** Map median/mean photophysical parameters across a standard 96-well format.
    * **Intra-well Boxplots:** Assess Field of View (FOV) variability within a single well.
    * **Inter-well Boxplots:** Compare technical replicates (3 wells) for a specific target protein.
    * **Comparison Matrices (R):** Generate normalized Log2 ratio matrices to benchmark multiple proteins against a reference (e.g., mEos3.2) across multiple plates.
* **Batch Excel Export:** Pool thousands of localizations across multiple experimental plates and export aligned statistics into multi-sheet Excel files.

## Project Structure

```text
ASHA/
│
├── notebooks/                  # Interactive exploration and pre-processing
│   ├── hcs_plate_exploration.ipynb
│   └── preprocessing.ipynb
│
├── scripts/                    # Ready-to-run Python execution scripts
│   ├── gradient_meos_*.py      # Gradients and linear/log fits for SpyCatcher
│   ├── spycatcher_*.py
│   ├── *_heatmap.py            # Plate and control heatmap generators
│   └── *_histograms.py         # Advanced ON/OFF duration analysis (frame.csv)
│
├── src/                        # Python Backend Core
│   ├── io_utils.py             # File parsing and data cleaning (PoCA, PALMTracer)
│   ├── preprocessing/          # Data filtering pipelines
│   ├── visualization/          # Heatmaps and boxplots rendering mechanics
│   └── save/                   # Data pooling and Excel export modules
│
└── R/                          # R Statistical & Multi-plot Suite
    ├── utils.R                 # I/O utilities for PoCA file extractions
    ├── table_*.R               # Specialized comparison matrices (Log2 ratios)
    └── all_protein_wells_*.R   # 2x4 and mega-grid visualization with stat tests
```

## Dependencies

### Python Requirements
* `pandas`
* `numpy`
* `matplotlib`
* `seaborn`
* `ipywidgets`
* `ipyfilechooser`
* `openpyxl` (for Excel export)

*Install via:* 
```bash
python -m venv venv

venv/source/activate

pip install -r requirements.txt
```


### R Requirements
* `ggplot2`
* `dplyr`
* `reshape2`
* `patchwork`

## Usage Guide

### 1. Interactive Data Exploration (Python)
Launch your Jupyter Notebook located in the `notebooks/` directory. The notebook utilizes `ipywidgets` to render a graphical interface.
1. Select your target **Plate Directory**.
2. Choose an **Action** (Experiment Statistics, Boxplots, or Plate Heatmap).
3. Select the **Photophysical Parameter** (e.g., `photon_loc`).
4. Set the **Scale** (Compare FOVs vs Compare 3 wells).
5. Click **Run** to generate and automatically save the PDF plots to the `results/` folder.

### 2. Batch Excel Export
Using the *Excel Exporter UI* in the notebook:
1. Click **+ Add Plate** to dynamically add multiple experiment directories using the interactive File Chooser.
2. Select the statistical method (**Median** or **Mean**).
3. Click **Run Extraction**. The script will pool all localizations per protein, compute the statistics, separate distinct plates with empty rows, and save everything into a clean, multi-sheet `results_combined.xlsx` file.

### 3. Protein Comparison Matrix (R)
To benchmark proteins against a reference across multiple experiments:
1. Open `R/comparison_matrix.R`.
2. Update the `base_dirs` vector with your experiment paths.
3. Run the script. It will extract the data, normalize it against the reference (default: mEos3.2), calculate the Log2 ratio, and output a high-resolution heatmap matrix (`comparison_matrix_[param].png`).