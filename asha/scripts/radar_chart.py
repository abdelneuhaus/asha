import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import pi

def plot_protein_radar_chart(excel_file):
    """
    Lit un fichier Excel multiparamétrique, calcule la moyenne par protéine,
    et génère un graphique en radar normalisé.
    """
    # Mapping pour renommer les paramètres
    name_mapping = {
        "photon_loc": "Photons/Loc",
        "intensity": "Photons/Mol",
        "total ON": "Total ON time",
        "avg_on": "ON duration",
        "avg_off": "OFF duration",
        "blinks": "NBlinks"
    }

    # 1. Lecture
    xls = pd.read_excel(excel_file, sheet_name=None)
    if 'global_metadata' in xls:
        del xls['global_metadata']
    
    # On filtre les feuilles selon le mapping
    available_params = [p for p in name_mapping.keys() if p in xls]
    
    # 2. Calcul et Normalisation
    proteins = xls[available_params[0]].columns.tolist()
    normalized_means = {prot: [] for prot in proteins}

    for param in available_params:
        df = xls[param]
        param_means = df.mean(numeric_only=True)
        max_val = param_means.max()
        for prot in proteins:
            val = param_means[prot] / max_val if max_val > 0 else 0
            normalized_means[prot].append(val)

    # 3. Préparation Radar Chart
    labels = [name_mapping[p] for p in available_params]
    N = len(labels)
    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1] 

    fig, ax = plt.subplots(figsize=(5, 5), subplot_kw=dict(polar=True))
    
    # Configuration des axes (Labels)
    plt.xticks(angles[:-1], labels, size=11, fontweight='bold')
    ax.tick_params(axis='x', pad=40)
    # Ajout des valeurs 0 à 1 sur l'axe radial
    ax.set_rlabel_position(0)
    plt.yticks([0.2, 0.4, 0.6, 0.8, 1.0], ["0.2", "0.4", "0.6", "0.8", "1.0"], color="grey", size=8)
    plt.ylim(0, 1)

    # Style
    ax.spines['polar'].set_visible(False)
    ax.grid(color='grey', linestyle='--', linewidth=0.5, alpha=0.5)

    # Tracé
    for prot in proteins:
        values = normalized_means[prot] + normalized_means[prot][:1]
        ax.plot(angles, values, linewidth=2, linestyle='solid', label=prot)
        ax.fill(angles, values, alpha=0.1)

    # Légende et Titre (English)
    plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), title="Proteins", frameon=False)
    
    plt.tight_layout()
    plt.show()

# Utilisation :
raw_data = plot_protein_radar_chart("/Users/aneuhaus/Desktop/asha/DATA/figure3d_&_3e.xlsx")