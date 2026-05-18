import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PAPER.utils import read_poca_files, get_poca_files, get_PALMTracer_files, read_locPALMTracer_file

# 240417_W1_FPs
# 240924_W1_FPs
# 241113_W1_FPs

# pathway = "//filer3/TEAM_M/everyone/_Transferts/Transfert Neuhaus/STABLE-FP/230627_stable_fp/A3/561-405.PT"
pathway = "D:/ANALYSIS_PAPER/241113_W1_FPs/D4/P 04"
# pathway = "D:/DPP/mAple3"
poca_files = sorted(get_poca_files(pathway))
pt_files = sorted(get_PALMTracer_files(pathway))
parameters = ['blinks', 'intensity', 'photon_loc', 'lifetime', 'avg_on', 'avg_off']

# --- 1. Chargement et Séparation par fichier ---
first_20_list = []
last_20_list = []

for poca_file, pt_file in zip(poca_files, pt_files):
    # 1. Lecture
    df_poca = read_poca_files(poca_file)
    
    # 2. Calculs AVANT filtrage (pour ne pas perdre de données)
    df_poca["avg_on"] = df_poca['total ON'] / df_poca['# seq ON']
    df_poca["avg_off"] = df_poca['total OFF'] / df_poca['# seq OFF']
    
    # Ton nouveau calcul de photon_loc (basé sur POCA uniquement)
    # On le fait AVANT le filtrage pour que l'alignement soit parfait
    df_poca["photon_loc"] = df_poca['intensity'] / df_poca['total ON']

    # 3. FILTRAGE (On garde les molécules < 500)
    df_poca = df_poca[df_poca['total ON'] <= 500].copy()
    
    # Conversion d'intensité (après filtrage c'est bon aussi)
    df_poca['intensity'] *= 0.012
    # On convertit aussi photon_loc si l'intensité n'était pas encore convertie
    df_poca['photon_loc'] *= 0.012 
    
    # 4. Seuils temporels 20%
    f_min, f_max = df_poca['frame'].min(), df_poca['frame'].max()
    threshold = (f_max - f_min) * 0.20
    
    # 5. Extraction First et Last (Tout est dans df_poca désormais)
    poca_f20 = df_poca[df_poca['frame'] <= (f_min + threshold)].copy()
    poca_l20 = df_poca[df_poca['frame'] >= (f_max - threshold)].copy()
    
    # On ajoute à la liste globale
    first_20_list.append(poca_f20)
    last_20_list.append(poca_l20)

# Fusion finale (Pas besoin de s'occuper de df_pt)
first_20 = pd.concat(first_20_list, ignore_index=True)
last_20 = pd.concat(last_20_list, ignore_index=True)

n_files = len(poca_files)

parameters_map = {
    'intensity': 'Photon budget',
    'photon_loc': 'Photons per loc',
    'total ON': 'Total ON Time (frames)',
    'avg_on': 'ON times',
    'avg_off': 'OFF times',
    'blinks': 'Number of blinks'
}

# --- 2. Plotting ---
fig, axes = plt.subplots(2, 3, figsize=(15, 11))
fig.suptitle(f"Données Concaténées (n = {n_files} fichiers)\nSource : POCA & PALMTracer (.PT)", fontsize=16)
axes_flat = axes.flatten()

for i, (col_name, title) in enumerate(parameters_map.items()):
    if col_name in first_20.columns:
        data_first = first_20[col_name].dropna()
        data_last = last_20[col_name].dropna()

        means = [data_first.mean(), data_last.mean()]
        labels = ['First 20%', 'Last 20%']
        
        # Barres (Vert et Bleu comme demandé)
        axes_flat[i].bar(labels, means, color=["#23C941", "#006AFF"], edgecolor='black', alpha=0.7, width=0.6)
        
        # Erreur (SEM)
        sem = [data_first.sem(), data_last.sem()]
        axes_flat[i].errorbar(labels, means, yerr=sem, fmt='none', c='black', capsize=5, lw=1.5)

        # Cosmétique
        axes_flat[i].set_title(title, fontweight='bold')
        axes_flat[i].set_ylabel("Mean Value")
        axes_flat[i].grid(axis='y', linestyle=':', alpha=0.5)
        
        # Ajustement de l'axe Y
        current_max = axes_flat[i].get_ylim()[1]
        axes_flat[i].set_ylim(0, current_max * 1.1)

plt.tight_layout(rect=[0, 0.03, 1, 0.93])
plt.show()