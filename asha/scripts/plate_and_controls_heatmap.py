import os

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from asha.src.io_utils import get_poca_files, read_poca_files, get_PALMTracer_files, read_locPALMTracer_file, creer_matrice_et_medianes, fusion


def plot_controlled_plate_heatmap(pathway_plate, pathway_control, param):
    """
    Combine protein and control data on one heatmap
    
    Args:
        pathway_plate (str): pathway to protein repertory.
        pathway_control (str): pathway to controls repertory.
        param (str): 'clusters' (total number) ou 'locs_per_plane' (mean per plane).
    """
    exp_name = f"{os.path.basename(os.path.normpath(pathway_plate))} + Controls"
    
    list_of_files = []
    
    if param == 'clusters':
        list_of_files = get_poca_files(pathway_plate) + get_poca_files(pathway_control)
    elif param in ['locs_per_plane', 'density']:
        list_of_files = get_PALMTracer_files(pathway_plate) + get_PALMTracer_files(pathway_control)
    else:
        print("Error: param must be 'clusters' or 'locs_per_plane'.")
        return
        
    if not list_of_files:
        print(f"Error : No files found in provided pathways for param='{param}'.")
        return
        
    heatmap_count_per_plane(exp_name=exp_name, param=param, list_of_files=list_of_files)


def heatmap_count_per_plane(exp_name, param, list_of_files):
    """
    Generates and displays a heatmap of clusters or molecular localizations per plane.

    Args:
        exp_name (str): identifier for the experiment, used in the plot title.
        param (str): metric to visualize.
        list_of_files (list[str]): A list of absolute file paths to the PoCA or PALMTracer data files to be analyzed.

    Returns:
        None: The function generates and displays the matplotlib figure directly and closes the plot instance upon completion.    
    """
    heatmap_data = []
    
    for f in list_of_files:
        if param == 'clusters':
            df = read_poca_files(f)
            heatmap_data.append(len(df))
            
        elif param in ['locs_per_plane', 'density']:
            df = read_locPALMTracer_file(f)
            
            if not df.empty and 'Plane' in df.columns:
                total_locs = len(df)
                unique_planes = df['Plane'].nunique()
                
                if unique_planes > 0:
                    heatmap_data.append(total_locs / unique_planes)
                else:
                    heatmap_data.append(np.nan)
            else:
                heatmap_data.append(np.nan)
                
    idx = list(map(str, range(1, 13)))
    cols = list(map(str, [chr(i) for i in range(ord('A'), ord('H')+1)]))
    all_wells = fusion(cols, idx)
    
    matrice_resultante = creer_matrice_et_medianes(list_of_files, heatmap_data, all_wells)
    df_heat = pd.DataFrame(np.array(matrice_resultante).reshape(8,12), index=cols, columns=idx)

    plt.figure(figsize=(12, 5))
    
    cmap = sns.color_palette("coolwarm", as_cmap=True)
    cmap.set_bad(color='#f0f0f0') 
    annotation_format = '.0f' if param == 'clusters' else '.1f'
    sns.heatmap(df_heat, annot=True, fmt=annotation_format, cmap=cmap, linewidths=1, linecolor='black')
    plt.yticks(rotation=0)
    
    title_suffix = "Average Locs / Plane" if param in ['locs_per_plane', 'density'] else "Total Clusters"
    plt.title(f"Heatmap: {title_suffix}\n(Plate: {exp_name})")    
    plt.show()
    plt.close('all')
    


if __name__ == '__main__':
    path_prot = "D:/ANALYSIS_PAPER/new_threshold/241113_W1_FPs"
    path_ctrl = "D:/ANALYSIS_PAPER/new_threshold/241113_W1_control"

    print("Generating Locs/Plane Heatmap...")
    plot_controlled_plate_heatmap(path_prot, path_ctrl, param='locs_per_plane')