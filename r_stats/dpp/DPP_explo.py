import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PAPER.utils import read_poca_files, get_poca_files, get_PALMTracer_files, read_locPALMTracer_file


# 240417_W1_FPs
# 240924_W1_FPs
# 241113_W1_FPs


# --- 1. Liste modulable des chemins d'accès ---
pathways_list = [
    "D:/ANALYSIS_PAPER/240924_W1_FPs/D5",
    # "D:/ANALYSIS_PAPER/241113_W1_FPs/D6",
    # "D:/ANALYSIS_PAPER/241113_W1_FPs/E6"
]

parameters_map = {
    'intensity': 'Photon budget',
    'photon_loc': 'Photons per loc',
    'total ON': 'Total ON Time (frames)',
    'avg_on': 'ON times',
    'avg_off': 'OFF times',
    'blinks': 'Number of blinks'
}

# --- 2. Chargement et calculs par fichier ---
file_level_data = []
total_files = 0

for pathway in pathways_list:
    if not os.path.exists(pathway):
        print(f"Attention : Le chemin {pathway} n'existe pas et sera ignoré.")
        continue
        
    poca_files = sorted(get_poca_files(pathway))
    pt_files = sorted(get_PALMTracer_files(pathway))
    
    # On itère sur les deux types de fichiers en même temps
    for poca_file, pt_file in zip(poca_files, pt_files):
        # --- TRAITEMENT POCA ---
        df_poca = read_poca_files(poca_file)
        
        if df_poca.empty:
            continue
            
        df_poca["avg_on"] = df_poca['total ON'] / df_poca['# seq ON']
        df_poca["avg_off"] = df_poca['total OFF'] / df_poca['# seq OFF']
        
        # Filtrage POCA
        df_poca = df_poca[df_poca['total ON'] <= 500].copy()
        if df_poca.empty:
            continue
            
        df_poca['intensity'] *= 0.012
        
        # Découpage temporel POCA
        f_min_p, f_max_p = df_poca['frame'].min(), df_poca['frame'].max()
        threshold_p = (f_max_p - f_min_p) * 0.20
        poca_f20 = df_poca[df_poca['frame'] <= (f_min_p + threshold_p)]
        poca_l20 = df_poca[df_poca['frame'] >= (f_max_p - threshold_p)]
        
        # --- TRAITEMENT PALMTracer ---
        df_pt = read_locPALMTracer_file(pt_file)
        
        if not df_pt.empty:
            # 1. Filtrage sur SigmaX(px) pour retirer les extremums (< 20% et > 80%)
            q20 = df_pt['SigmaX(px)'].quantile(0.20)
            q80 = df_pt['SigmaX(px)'].quantile(0.80)
            df_pt = df_pt[(df_pt['SigmaX(px)'] >= q20) & (df_pt['SigmaX(px)'] <= q80)].copy()
            
            # 2. Conversion (Optionnel : à retirer si Integrated_Intensity est déjà en photons)
            df_pt['Integrated_Intensity'] *= 0.012
            
            # 3. Découpage temporel PALMTracer
            # (Assure-toi que la colonne s'appelle bien 'frame' et non 'Frame' dans tes fichiers PT)
            f_min_pt, f_max_pt = df_pt['Plane'].min(), df_pt['Plane'].max()
            threshold_pt = (f_max_pt - f_min_pt) * 0.20
            pt_f20 = df_pt[df_pt['Plane'] <= (f_min_pt + threshold_pt)]
            pt_l20 = df_pt[df_pt['Plane'] >= (f_max_pt - threshold_pt)]
        else:
            # Sécurité si le fichier pt est vide
            pt_f20 = pd.DataFrame()
            pt_l20 = pd.DataFrame()

        # --- COMPILATION DES MOYENNES DU FICHIER ---
        file_record = {}
        for col_name in parameters_map.keys():
            if col_name == 'photon_loc':
                # On utilise PALMTracer pour photon_loc
                if not pt_f20.empty and 'Integrated_Intensity' in pt_f20.columns:
                    file_record[f'{col_name}_First20'] = pt_f20['Integrated_Intensity'].median()
                    file_record[f'{col_name}_Last20'] = pt_l20['Integrated_Intensity'].median()
                else:
                    file_record[f'{col_name}_First20'] = np.nan
                    file_record[f'{col_name}_Last20'] = np.nan
            else:
                # On utilise POCA pour le reste
                if col_name in df_poca.columns:
                    file_record[f'{col_name}_First20'] = poca_f20[col_name].median()
                    file_record[f'{col_name}_Last20'] = poca_l20[col_name].median()
                
        file_level_data.append(file_record)
        total_files += 1

# Création du DataFrame global
df_results = pd.DataFrame(file_level_data)

# --- 3. Plotting ---
fig, axes = plt.subplots(2, 3, figsize=(15, 11))
fig.suptitle(f"Données moyennées par fichier (n = {total_files} fichiers)\nSource : POCA & PALMTracer", fontsize=16)
axes_flat = axes.flatten()

for i, (col_name, title) in enumerate(parameters_map.items()):
    col_f20 = f'{col_name}_First20'
    col_l20 = f'{col_name}_Last20'
    
    if col_f20 in df_results.columns and col_l20 in df_results.columns:
        # Moyenne et SD calculées sur les fichiers (en omettant les NaNs éventuels)
        mean_f20 = df_results[col_f20].median()
        mean_l20 = df_results[col_l20].median()
        
        sd_f20 = df_results[col_f20].std()
        sd_l20 = df_results[col_l20].std()

        means = [mean_f20, mean_l20]
        sds = [sd_f20, sd_l20]
        labels = ['First 20%', 'Last 20%']
        
        axes_flat[i].bar(labels, means, color=["#23C941", "#006AFF"], edgecolor='black', alpha=0.7, width=0.6)
        axes_flat[i].errorbar(labels, means, yerr=sds, fmt='none', c='black', capsize=5, lw=1.5)

        axes_flat[i].set_title(title, fontweight='bold')
        axes_flat[i].set_ylabel("Mean Value per file")
        axes_flat[i].grid(axis='y', linestyle=':', alpha=0.5)
        
        # Sécurité pour ajuster l'axe Y (évite les bugs si toutes les valeurs sont NaN)
        if not np.isnan(mean_f20) and not np.isnan(sd_f20):
            max_val_with_error = max(mean_f20 + sd_f20, mean_l20 + sd_l20)
            axes_flat[i].set_ylim(0, max_val_with_error * 1.15)

plt.tight_layout(rect=[0, 0.03, 1, 0.93])
plt.show()