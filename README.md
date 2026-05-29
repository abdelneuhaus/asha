# HCS-SMLM Data Analysis Pipeline

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

## 📂 Project Structure

```text
HCS-SMLM-Pipeline/
│
├── notebooks/              # Jupyter Notebooks containing the interactive UIs
│   └── interactive_dashboard.ipynb
│
├── src/                    # Python Backend
│   ├── __init__.py
│   ├── io_utils.py         # File parsing and data cleaning (get_poca_files, etc.)
│   ├── visualization.py    # Matplotlib/Seaborn functions for heatmaps and boxplots
│   └── save_as_excel.py    # Logic for pooling and exporting data to Excel
│
├── R/                      # R Backend for statistical matrices
│   ├── io_utils.R          # R equivalent for file parsing
│   └── comparison_matrix.R # Log2 ratio matrix generation with ggplot2
│
├── data/                   # (Optional) Place raw experimental data here
└── results/                # Output directory for exported PDFs and Excel files