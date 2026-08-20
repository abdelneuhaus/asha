import math
import os
import pandas as pd

from src.io_utils import get_statMIA_files, read_statMIA, get_poca_files, read_poca_files


def count_filtering(plate_path:str, nlocs:int=60, len_window:int=300, 
                              sigma_min:float=0.5, sigma_max:float=1.5,
                              intensity_min:int=0, intensity_max:int=10e6):
    """
    Analyse les fichiers statMIA pour le filtrage des localisations et les fichiers 
    PoCA pour le décompte post-clustering.
    """
    
    # ==========================================
    # ÉTAPE 1 : Filtrage des fichiers statMIA
    # ==========================================
    list_of_mia_files = get_statMIA_files(plate_path)

    total_initial = 0
    total_filtered_params = 0
    total_filtered_window = 0
    total_final = 0

    for mia_file in list_of_mia_files:
        raw_data = read_statMIA(mia_file)
        
        initial_locs = len(raw_data)
    
        raw_data['Intensity Gauss'] *= 2 * math.pi * raw_data['SigmaX'] * raw_data['SigmaY']
        raw_data = raw_data[(raw_data['SigmaX'] > sigma_min) & (raw_data['SigmaX'] < sigma_max)]
        raw_data = raw_data[(raw_data['Intensity Gauss'] >= intensity_min) & (raw_data['Intensity Gauss'] <= intensity_max)]
        
        locs_after_params = len(raw_data)
        filtered_by_params = initial_locs - locs_after_params
        
        loc_per_frame = raw_data.groupby('Plane').size()

        try:
            frame_with_maximum = 0  # loc_per_frame.idxmax()
        except ValueError:
            pass

        raw_data = raw_data[(raw_data['Plane'] >= frame_with_maximum)]
        loc_per_frame = raw_data.groupby('Plane').size()
        rolling_avg = loc_per_frame.rolling(window=len_window).mean()
        valid_frames = rolling_avg[rolling_avg < nlocs].index
        
        locs_after_window = 0
        filtered_by_window = locs_after_params 

        if not valid_frames.empty:
            start_frame = valid_frames[0]
            raw_data = raw_data[(raw_data['Plane'] >= start_frame)]
            locs_after_window = len(raw_data)
            filtered_by_window = locs_after_params - locs_after_window
            loc_per_frame = loc_per_frame[loc_per_frame.index >= start_frame]
            rolling_avg = rolling_avg[rolling_avg.index >= start_frame]
            
        total_initial += initial_locs
        total_filtered_params += filtered_by_params
        total_filtered_window += filtered_by_window
        total_final += locs_after_window


    # ==========================================
    # ÉTAPE 2 : Décompte dans les fichiers PoCA
    # ==========================================
    list_of_poca_files = get_poca_files(plate_path)
    
    total_poca_molecules = 0
    total_poca_locs = 0
    
    for poca_file in list_of_poca_files:
        df_poca = read_poca_files(poca_file)
        
        if 'total ON' in df_poca.columns:
            total_poca_molecules += len(df_poca)
            total_poca_locs += df_poca['total ON'].sum()
        else:
            print(f"Attention : Colonne 'total ON' manquante dans {os.path.basename(poca_file)}")


    # ==========================================
    # AFFICHAGE DU RÉSUMÉ
    # ==========================================
    print("\n" + "="*50)
    print("RÉSUMÉ GLOBAL DU FILTRAGE ET CLUSTERING")
    print("="*50)
    print("ÉTAPE 1 : Localisations brutes (statMIA)")
    print(f"Total initialement détecté           : {total_initial}")
    print(f"Retirées par les filtres (Sigma/Int) : {total_filtered_params}")
    print(f"Retirées par le Sliding Window       : {total_filtered_window}")
    print(f"TOTAL LOCALISATIONS CONSERVÉES       : {total_final}")
    print("-" * 50)
    print("ÉTAPE 2 : Clusters validés (PoCA)")
    print(f"Nombre de molécules (clusters)       : {total_poca_molecules}")
    print(f"Localisations dans les clusters (ON) : {int(total_poca_locs)}")
    print("-" * 50)
    
    # Calcul des localisations rejetées par le DBSCAN (bruit)
    if total_final > 0:
        locs_perdues = total_final - int(total_poca_locs)
        pourcentage_perte = (locs_perdues / total_final) * 100
        print(f"Localisations non-clusterisées       : {locs_perdues} ({pourcentage_perte:.1f}%)")
    
    print("="*50 + "\n")
    
    print("Done")