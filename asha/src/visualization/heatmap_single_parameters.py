import os

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from src.io_utils import get_poca_files, read_poca_files, creer_matrice_et_medianes, fusion


def plot_plate_heatmap(pathway, param):
    """
    Prépare les données et lance la heatmap pour la plaque SMLM-HCS.
    """
    exp_name = os.path.basename(os.path.normpath(pathway))
    list_of_poca = get_poca_files(pathway) 
    if not list_of_poca:
        print(f"Error : No PoCA files found in {pathway}")
        return
    heatmap_single_parameters(exp_name=exp_name, param=param, list_of_poca_files=list_of_poca, stats=np.median)


def heatmap_single_parameters(exp_name, param, list_of_poca_files, stats=np.median):
    """
    Génère une heatmap à partir des fichiers PoCA pour un paramètre photophysique donné.
    """
    heatmap_data = []
    
    for f in list_of_poca_files:
        df = read_poca_files(f)
        
        if param == 'avg_on':
            values = np.divide(df['total ON'], df['# seq ON'], out=np.zeros_like(df['total ON'], dtype=float), where=df['# seq ON']!=0)
        
        elif param == 'avg_off':
            values = np.divide(df['total OFF'], df['# seq OFF'], out=np.zeros_like(df['total OFF'], dtype=float), where=df['# seq OFF']!=0)
            
        elif param == 'photon_loc':
            values = np.divide(df['intensity'], df['total ON'], out=np.zeros_like(df['intensity'], dtype=float), where=df['total ON']!=0)
            
        else:
            values = df[param] 
            
        clean_values = pd.Series(values).replace([np.inf, -np.inf], np.nan).dropna()
        if len(clean_values) > 0:
            heatmap_data.append(float(stats(clean_values)))
        else:
            heatmap_data.append(np.nan) # empty wells
            
    idx = list(map(str, range(1, 13)))
    cols = list(map(str, [chr(i) for i in range(ord('A'), ord('H')+1)]))
    all_wells = fusion(cols, idx)
    matrice_resultante = creer_matrice_et_medianes(list_of_poca_files, heatmap_data, all_wells)
    df_heat = pd.DataFrame(np.array(matrice_resultante).reshape(8,12), index=cols, columns=idx)

    plt.figure(figsize=(12, 5))
    sns.heatmap(df_heat, annot=True, fmt='.1f', cmap="coolwarm", linewidths=1, linecolor='black')
    plt.yticks(rotation=0)
    plt.title(f"Heatmap: {param.upper()} (Plate: {exp_name})")    
    plt.show()
    plt.close('all')
