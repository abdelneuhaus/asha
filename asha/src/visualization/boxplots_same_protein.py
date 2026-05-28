import re
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from src.io_utils import get_poca_files, read_poca_files, create_database



def plot_well_boxplots(pathway, param, protein):
    """
    Génère un boxplot comparant les réplicats (puits individuels) d'une protéine spécifique.
    """

    db = create_database()
    if protein not in db:
        print(f"Error : Protein '{protein}' is not in the experiment.")
        return
        
    target_wells = db[protein][0]
    
    all_poca = get_poca_files(pathway)
    if not all_poca:
        print("No PoCA files found in this directory.")
        return
        
    data_frames = []
    
    for f in all_poca:
        match = re.search(r'([A-H]\d+)', f)
        if not match:
            continue
        well_name = match.group(1)
        
        if well_name in target_wells:
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
            
            temp_df = pd.DataFrame({'Value': clean_values, 'Puits': well_name})
            data_frames.append(temp_df)
            
    if not data_frames:
        print(f"No files for wells {target_wells}.")
        return

    final_df = pd.concat(data_frames, ignore_index=True)
    

    plt.figure(figsize=(8, 6))
    sns.boxplot(data=final_df, x='Puits', y='Value', order=target_wells, hue='Puits', legend=False, palette="Set2", showfliers=False)
    plt.title(f"Comparaison des réplicats : {protein}\n(Paramètre : {param.upper()})", fontsize=14)
    plt.ylabel(param)
    plt.xlabel("Puits")
    plt.grid(axis='y', linestyle='--', alpha=0.7)    
    plt.show()
    plt.close()