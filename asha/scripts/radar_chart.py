import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import pi

def plot_protein_radar_chart(excel_file):
    """
    Generates a grid of individual radar charts comparing protein photophysical properties.

    Args:
        excel_file (str): Absolute path to the Excel file containing the analysis data. Each sheet name must correspond to a key in the internal mapping

    Returns:
        None: Displays the generated matplotlib figure directly.    
    """
    
    name_mapping = {
        "photon_loc": "Photons/Loc",
        "intensity": "Photons/Mol",
        "total ON": "Total ON time",
        "avg_on": "ON duration",
        "avg_off": "OFF duration",
        "blinks": "NBlinks"
    }

    xls = pd.read_excel(excel_file, sheet_name=None)
    if 'global_metadata' in xls: del xls['global_metadata']
    
    available_params = [p for p in name_mapping.keys() if p in xls]
    proteins = xls[available_params[0]].columns.tolist()
    labels = [name_mapping[p] for p in available_params]
    
    normalized_means = {prot: [] for prot in proteins}
    for param in available_params:
        param_means = xls[param].mean(numeric_only=True)
        max_val = param_means.max()
        for prot in proteins:
            val = param_means[prot] / max_val if max_val > 0 else 0
            if param in ["avg_off", "blinks"]: val = 1.0 - val
            normalized_means[prot].append(val)

    N_params = len(labels)
    angles = [n / float(N_params) * 2 * pi for n in range(N_params)] + [0]
    
    fig, axes = plt.subplots(2, 4, figsize=(16, 8), subplot_kw=dict(polar=True))
    axes = axes.flatten() 

    for i, prot in enumerate(proteins):
        ax = axes[i]
        values = normalized_means[prot] + normalized_means[prot][:1]
        ax.plot(angles, values, linewidth=2, color='tab:blue')
        ax.fill(angles, values, alpha=0.2, color='tab:blue')
        
        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(labels, size=8)
        ax.set_ylim(0, 1)
        ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
        ax.set_yticklabels([])
        ax.set_title(prot, size=12, fontweight='bold', pad=15)
        ax.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)

    axes[7].set_visible(False)
    
    plt.tight_layout()
    plt.show()

raw_data = plot_protein_radar_chart("/Users/aneuhaus/Desktop/asha/DATA/figure3d_&_3e.xlsx")