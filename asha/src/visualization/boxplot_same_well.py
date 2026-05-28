import re

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from src.io_utils import get_poca_files, read_poca_files


def plot_fov_boxplots(pathway, param, protein, well):
    """
    Generates a boxplot comparing individual Fields of View (FOV) within a single well.

    Args:
        pathway (str): Path to the experiment directory.
        param (str): Photophysical parameter to visualize (e.g., 'intensity').
        protein (str): Protein name for the plot title.
        well (str): The identifier of the well to analyze (e.g., 'C3').

    Returns:
        None: Displays the plot and saves it as a PDF in the results folder.    
    """

    all_poca = get_poca_files(pathway)
    if not all_poca:
        print("Error: No PoCa files found.")
        return
        
    data_frames = []
    fov_count = 1
    
    for f in all_poca:
        match = re.search(r'([A-H]\d+)', f)
        if not match:
            continue
            
        current_well = match.group(1)
        
        if current_well == well:
            df = read_poca_files(f)
            df['intensity'] *= (3.6 / 300)
            if param == 'avg_on':
                values = np.divide(df['total ON'], df['# seq ON'], out=np.zeros_like(df['total ON'], dtype=float), where=df['# seq ON']!=0)
            elif param == 'avg_off':
                values = np.divide(df['total OFF'], df['# seq OFF'], out=np.zeros_like(df['total OFF'], dtype=float), where=df['# seq OFF']!=0)
            elif param == 'photon_loc':
                values = np.divide(df['intensity'], df['total ON'], out=np.zeros_like(df['intensity'], dtype=float), where=df['total ON']!=0)
            else:
                values = df[param]
                
            clean_values = pd.Series(values).replace([np.inf, -np.inf], np.nan).dropna()
            
            # Pour identifier le FOV, on peut utiliser le nom du dossier parent (ex: 'SR_001.MIA' ou juste un compteur)
            fov_name = f"FOV {fov_count}"
            
            temp_df = pd.DataFrame({
                'Value': clean_values,
                'FOV': fov_name
            })
            data_frames.append(temp_df)
            fov_count += 1
            
    if not data_frames:
        print(f"No FOV for well {well}.")
        return

    final_df = pd.concat(data_frames, ignore_index=True)
    
    plt.figure(figsize=(8, 6))
    sns.boxplot(data=final_df, x='FOV', y='Value', hue='FOV', legend=False, palette="Pastel1", showfliers=False)
    plt.title(f"Intra-well variability : {protein} (Well {well})\n(Parameter : {param.upper()})", fontsize=14)
    plt.ylabel(param)
    plt.xlabel("Replicates (FOVs)")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.show()
    plt.close()