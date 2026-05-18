import os
import pandas as pd
import numpy as np
from PAPER.utils import read_poca_files, get_poca_files, get_PALMTracer_files, read_locPALMTracer_file

# --- Configuration des chemins ---
# Rendu modulable sous forme de liste pour traiter un ou plusieurs dossiers facilement
pathways_list = [
    "D:/ANALYSIS_PAPER/240417_W1_FPs/C5",
    "D:/ANALYSIS_PAPER/240417_W1_FPs/D7",
    "D:/ANALYSIS_PAPER/240924_W1_FPs/C5"
]

# Le fichier Excel sera sauvegardé dans le premier dossier de la liste
output_excel = os.path.join(pathways_list[0], "D:/DPP/mEOSEM.xlsx")

all_data = {
    'Photon budget': {'first': [], 'last': []},
    'Photons per loc': {'first': [], 'last': []}, # Calculé via PALMTracer
    'Bleach time': {'first': [], 'last': []},
    'On times': {'first': [], 'last': []},
    'Off times': {'first': [], 'last': []},
    'Number of blinks': {'first': [], 'last': []}
}

total_files = 0

for pathway in pathways_list:
    if not os.path.exists(pathway):
        print(f"Attention : Le chemin {pathway} n'existe pas et sera ignoré.")
        continue
        
    poca_files = sorted(get_poca_files(pathway))
    pt_files = sorted(get_PALMTracer_files(pathway))

    print(f"Traitement de {len(poca_files)} fichiers dans {pathway}...")

    for poca_file, pt_file in zip(poca_files, pt_files):
        # --- 1. TRAITEMENT POCA ---
        df_poca = read_poca_files(poca_file)
        
        if df_poca.empty:
            continue
            
        df_poca["avg_on"] = df_poca['total ON'] / df_poca['# seq ON']
        df_poca["avg_off"] = df_poca['total OFF'] / df_poca['# seq OFF']

        # Application du filtre (<= 500)
        df_poca = df_poca[df_poca['total ON'] <= 100].copy()
        if df_poca.empty:
            continue

        df_poca['intensity'] *= 0.012

        # Seuils temporels POCA
        f_min_p, f_max_p = df_poca['frame'].min(), df_poca['frame'].max()
        threshold_p = (f_max_p - f_min_p) * 0.20
        p_f20 = df_poca[df_poca['frame'] <= (f_min_p + threshold_p)]
        p_l20 = df_poca[df_poca['frame'] >= (f_max_p - threshold_p)]

        # --- 2. TRAITEMENT PALMTracer (pour photon_loc) ---
        df_pt = read_locPALMTracer_file(pt_file)
        if not df_pt.empty:
            # Filtre sur SigmaX(px) (crop des bords <20% et >80%)
            q20 = df_pt['SigmaX(px)'].quantile(0.20)
            q80 = df_pt['SigmaX(px)'].quantile(0.80)
            df_pt = df_pt[(df_pt['SigmaX(px)'] >= q20) & (df_pt['SigmaX(px)'] <= q80)].copy()
            
            df_pt['Integrated_Intensity'] *= 0.012
            
            # Seuils temporels PALMTracer
            f_min_pt, f_max_pt = df_pt['Plane'].min(), df_pt['Plane'].max()
            threshold_pt = (f_max_pt - f_min_pt) * 0.20
            pt_f20 = df_pt[df_pt['Plane'] <= (f_min_pt + threshold_pt)]
            pt_l20 = df_pt[df_pt['Plane'] >= (f_max_pt - threshold_pt)]
        else:
            pt_f20 = pd.DataFrame()
            pt_l20 = pd.DataFrame()

        # --- 3. ACCUMULATION DES MOYENNES (1 valeur par fichier) ---
        
        # Photon budget (POCA : intensity)
        all_data['Photon budget']['first'].append(p_f20['intensity'].mean() if not p_f20.empty else np.nan)
        all_data['Photon budget']['last'].append(p_l20['intensity'].mean() if not p_l20.empty else np.nan)
        
        # Photons per loc (PALMTracer : Integrated_Intensity)
        if not pt_f20.empty and 'Integrated_Intensity' in pt_f20.columns:
            all_data['Photons per loc']['first'].append(pt_f20['Integrated_Intensity'].mean())
        else:
            all_data['Photons per loc']['first'].append(np.nan)
            
        if not pt_l20.empty and 'Integrated_Intensity' in pt_l20.columns:
            all_data['Photons per loc']['last'].append(pt_l20['Integrated_Intensity'].mean())
        else:
            all_data['Photons per loc']['last'].append(np.nan)
        
        # Bleach time (POCA : total ON)
        all_data['Bleach time']['first'].append(p_f20['total ON'].mean() if not p_f20.empty else np.nan)
        all_data['Bleach time']['last'].append(p_l20['total ON'].mean() if not p_l20.empty else np.nan)
        
        # On times (POCA : avg_on)
        all_data['On times']['first'].append(p_f20['avg_on'].mean() if not p_f20.empty else np.nan)
        all_data['On times']['last'].append(p_l20['avg_on'].mean() if not p_l20.empty else np.nan)
        
        # Off times (POCA : avg_off)
        all_data['Off times']['first'].append(p_f20['avg_off'].mean() if not p_f20.empty else np.nan)
        all_data['Off times']['last'].append(p_l20['avg_off'].mean() if not p_l20.empty else np.nan)
        
        # Number of blinks (POCA : blinks)
        all_data['Number of blinks']['first'].append(p_f20['blinks'].mean() if not p_f20.empty else np.nan)
        all_data['Number of blinks']['last'].append(p_l20['blinks'].mean() if not p_l20.empty else np.nan)

        total_files += 1

# --- 4. EXPORT VERS EXCEL ---
print(f"Exportation vers : {output_excel}")
with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
    for sheet_name, values in all_data.items():
        # Création du DataFrame pour l'onglet courant
        df_export = pd.DataFrame({
            'First 20%': pd.Series(values['first']),
            'Last 20%': pd.Series(values['last'])
        })
        # sheet_name[:31] assure que le nom de l'onglet ne dépasse pas la limite d'Excel
        df_export.to_excel(writer, sheet_name=sheet_name[:31], index=False)

print(f"Processus terminé. Les moyennes de {total_files} fichiers ont été exportées.")