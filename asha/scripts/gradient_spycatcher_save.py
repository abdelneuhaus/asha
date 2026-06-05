import os
import numpy as np
import pandas as pd
from src.io_utils import read_poca_files, get_poca_files, get_PALMTracer_files, read_locPALMTracer_file



list_of_pt_files = get_PALMTracer_files('D:/ANALYSIS_PAPER/gamme/SPYCATCHER')
list_of_poca_files = get_poca_files('D:/ANALYSIS_PAPER/gamme/SPYCATCHER')

spycatcher_wells = ["C9", "C10", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "E3", "E4", "E5", "E6"]
concentration = [10, 7.5, 5, 2, 1, 0.75, 0.5, 0.2, 0.1, 0.075, 0.05, 0.02, 0.01, 0]

surface = (256 * 0.16) * (256 * 0.16)
export_data = []

print("Extracting data...")

for position, conc in zip(spycatcher_wells, concentration):
    position_files_pt = [f for f in list_of_pt_files if f"/{position}/" in f]
    position_files_poca = [f for f in list_of_poca_files if f"/{position}/" in f]
    
    density_per_file = []
    cluster_per_file = []
    
    for file in position_files_pt:
        raw_file_pt = read_locPALMTracer_file(file)
        if not raw_file_pt.empty:
            total_locs = len(raw_file_pt)
            unique_planes = raw_file_pt['Plane'].nunique()
            
            if unique_planes > 0:
                mean_density = (total_locs / unique_planes) / surface
                if mean_density > 0:
                    density_per_file.append(mean_density)
            
    for file in position_files_poca:
        raw_file_poca = read_poca_files(file)
        cluster_per_file.append(len(raw_file_poca))

    val_density_finale = sum(density_per_file) / len(density_per_file) if density_per_file else np.nan
    val_cluster_final = sum(cluster_per_file) / len(cluster_per_file) if cluster_per_file else np.nan
    
    export_data.append({
        "Concentration %": conc,
        "Mean loc density (loc/um2)": val_density_finale,
        "Mean #Clusters": val_cluster_final
    })

output_excel = "spycatcher_gradient.xlsx"
df_results = pd.DataFrame(export_data)
df_results.to_excel(output_excel, index=False)

print(f"Finished! File '{output_excel}' created.")