import os
import pandas as pd
import numpy as np
from PAPER.utils import read_poca_files, get_poca_files, get_PALMTracer_files, read_locPALMTracer_file

# --- Configuration des chemins ---
pathway = "D:/DPP/mAple3"
output_excel = os.path.join(pathway, "Export_Stats_First_vs_Last_20.xlsx")

poca_files = sorted(get_poca_files(pathway))
pt_files = sorted(get_PALMTracer_files(pathway)) # On les garde si tu en as besoin ailleurs

all_data = {
    'Photon budget': {'first': [], 'last': []},
    'Photons per loc': {'first': [], 'last': []}, # Sera calculé via POCA
    'Bleach time': {'first': [], 'last': []},
    'On times': {'first': [], 'last': []},
    'Off times': {'first': [], 'last': []},
    'Number of blinks': {'first': [], 'last': []}
}

print(f"Traitement de {len(poca_files)} fichiers...")

for poca_file, pt_file in zip(poca_files, pt_files):
    # 1. Lecture du fichier POCA uniquement (puisque tout vient de là maintenant)
    df_poca = read_poca_files(poca_file)
    
    # 2. SEUILS TEMPORELS (Sur l'acquisition totale avant filtre)
    f_min, f_max = df_poca['frame'].min(), df_poca['frame'].max()
    threshold = (f_max - f_min) * 0.20
    limit_first = f_min + threshold
    limit_last = f_max - threshold

    # 3. CALCULS DES PARAMÈTRES
    # On convertit l'intensité d'abord
    df_poca['intensity_conv'] = df_poca['intensity'] * 0.012
    
    # Calcul de photon_loc : Intégralité des photons / Nombre de détections
    df_poca["photon_loc_calc"] = df_poca['intensity_conv'] / df_poca['total ON']
    
    df_poca["avg_on"] = df_poca['total ON'] / df_poca['# seq ON']
    df_poca["avg_off"] = df_poca['total OFF'] / df_poca['# seq OFF']

    # 4. APPLICATION DU FILTRE (<= 500)
    # Ce filtre s'applique maintenant à TOUTES les colonnes calculées au-dessus
    df_poca_f = df_poca[df_poca['total ON'] <= 500].copy()

    # 5. RÉPARTITION FIRST / LAST
    p_f20 = df_poca_f[df_poca_f['frame'] <= limit_first]
    p_l20 = df_poca_f[df_poca_f['frame'] >= limit_last]

    # 6. ACCUMULATION DANS LE DICTIONNAIRE
    all_data['Photon budget']['first'].extend(p_f20['intensity_conv'].tolist())
    all_data['Photon budget']['last'].extend(p_l20['intensity_conv'].tolist())
    
    all_data['Photons per loc']['first'].extend(p_f20['photon_loc_calc'].tolist())
    all_data['Photons per loc']['last'].extend(p_l20['photon_loc_calc'].tolist())
    
    all_data['Bleach time']['first'].extend(p_f20['total ON'].tolist())
    all_data['Bleach time']['last'].extend(p_l20['total ON'].tolist())
    
    all_data['On times']['first'].extend(p_f20['avg_on'].tolist())
    all_data['On times']['last'].extend(p_l20['avg_on'].tolist())
    
    all_data['Off times']['first'].extend(p_f20['avg_off'].tolist())
    all_data['Off times']['last'].extend(p_l20['avg_off'].tolist())
    
    all_data['Number of blinks']['first'].extend(p_f20['blinks'].tolist())
    all_data['Number of blinks']['last'].extend(p_l20['blinks'].tolist())

# --- 7. EXPORT VERS EXCEL ---
print(f"Exportation vers : {output_excel}")
with pd.ExcelWriter(output_excel, engine='openpyxl') as writer:
    for sheet_name, values in all_data.items():
        df_export = pd.DataFrame({
            'First 20%': pd.Series(values['first']),
            'Last 20%': pd.Series(values['last'])
        })
        df_export.to_excel(writer, sheet_name=sheet_name[:31], index=False)

print("Processus terminé. Toutes les données sont synchronisées via POCA.")