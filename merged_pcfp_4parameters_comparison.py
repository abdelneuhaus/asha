from PAPER.utils import read_poca_files, get_poca_files, get_csv_poca_frame_files, get_csv_poca_intensity_files, get_csv_poca_sigma_files
from preprocessing import photon_calculation, pre_process_off_frame_csv, pre_process_on_frame_csv, pre_process_sigma, pre_process_single_intensity, get_and_save_well_and_FOV

import os
import seaborn as sns
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt

pathway = 'D:/240925_W1_FPs'
list_of_poca_files = get_poca_files(pathway)
csv_intensity_files = get_csv_poca_intensity_files(pathway)
csv_frame_files = get_csv_poca_frame_files(pathway)
csv_sigma_files = get_csv_poca_sigma_files(pathway)

heatmap_data = []

# Calcul des statistiques pour chaque fichier
for f in range(len(list_of_poca_files)):
    raw_file_poca = read_poca_files(list_of_poca_files[f])
    raw_file_poca = raw_file_poca[raw_file_poca['blinks'] <  np.quantile(raw_file_poca['blinks'], 1)]
    raw_file_poca = raw_file_poca[raw_file_poca['total ON'] < np.quantile(raw_file_poca['total ON'], 1)]
    
    well_data = dict()
    on_times = pre_process_on_frame_csv(csv_frame_files[f])
    off_times = pre_process_off_frame_csv(csv_frame_files[f])
    well_data['_num_blinks'] = raw_file_poca['blinks'].values
    well_data['_num_on_times'] = raw_file_poca['# seq ON'].values.tolist()
    well_data['_num_off_times'] = raw_file_poca['# seq OFF'].values.tolist()
    well_data['_on_times'] = on_times
    well_data['_off_times'] = off_times
    well_data['_phot_per_cluster'] = photon_calculation(raw_file_poca['intensity'].values.tolist())
    well_data['_phot_per_loc'] = (np.array(raw_file_poca['intensity']) / np.array(raw_file_poca['total ON'])) * (3.6/300) 
    well_data['_total_on'] = raw_file_poca['total ON'].values
    heatmap_data.append(well_data)

# Définir les groupes de puits et leurs noms associés
# well_groups = [
#     ["C3", "D9", "E3"], ["C4", "D8", "E4"], ["C5", "D7", "E5"],
#     ["C6", "D6", "E6"], ["C7", "D5", "E7"], ["C8", "D4", "E8"],
#     ["C9", "D3", "E9"]
# ]
# group_names = ["mEos4b-L93M", "mEos3.2", "mEosEM", "pcSTAR", "mAple", "Dendra", "mEos4b"]


"""
mEos Family only
"""
# Définir les groupes de puits et leurs noms associés
well_groups = [
    ["C3", "D9", "E3"], ["C4", "D8", "E4"], ["C5", "D7", "E5"],
    ["C6", "D6", "E6"], ["C7", "D5", "E7"], ["C8", "D4", "E8"], ["C9", "D3", "E9"]
]
group_names = ["mEos4b-L93M", "mEos3.2", "mEosEM", "pcSTAR", "Dendra2", "mMaple3", "mEos4b"]

# Structure pour stocker les données pour chaque métrique
boxplot_data_phot_per_cluster = []
boxplot_data_phot_per_loc = []
boxplot_data_total_on = []
boxplot_data_num_blinks = []

for wells in well_groups:
    group_phot_per_cluster = []
    group_phot_per_loc = []
    group_total_on = []
    group_num_blinks = []
    
    for well in wells:
        well_files = [file for file in list_of_poca_files if well in file]
        
        for f in well_files:
            raw_file_poca = read_poca_files(f)
            
            total_on_filtered = raw_file_poca[raw_file_poca['total ON'] < np.quantile(raw_file_poca['total ON'], 1)]['total ON']
            group_total_on.extend(total_on_filtered.values.tolist())
            
            phot_per_cluster = photon_calculation(raw_file_poca['intensity'].values.tolist())
            group_phot_per_cluster.extend(phot_per_cluster)
            
            phot_per_loc = (np.array(raw_file_poca['intensity']) / np.array(raw_file_poca['total ON'])) * (3.6/300)
            group_phot_per_loc.extend(phot_per_loc.tolist())
            
            num_blinks = raw_file_poca[raw_file_poca['blinks'] < np.quantile(raw_file_poca['blinks'], 1)]['blinks']
            group_num_blinks.extend(num_blinks.values.tolist())
    
    boxplot_data_phot_per_cluster.append(group_phot_per_cluster)
    boxplot_data_phot_per_loc.append(group_phot_per_loc)
    boxplot_data_total_on.append(group_total_on)
    boxplot_data_num_blinks.append(group_num_blinks)

# Création des sous-figures
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
custom_palette = ['dodgerblue', 'salmon', 'palegreen', "tomato", "black", "red", "orchid"]

# Phot per cluster
sns.boxplot(data=boxplot_data_phot_per_cluster, ax=axes[0, 0], showfliers=False, palette=custom_palette)
axes[0, 0].set_xticklabels(group_names, rotation=0)
axes[0, 0].set_ylim(-3000, 22000)

# Phot per loc
sns.boxplot(data=boxplot_data_phot_per_loc, ax=axes[0, 1], showfliers=False, palette=custom_palette)
axes[0, 1].set_xticklabels(group_names, rotation=0)
axes[0, 1].set_ylim(50, 175)

# Total ON
sns.boxplot(data=boxplot_data_total_on, ax=axes[1, 0], showfliers=False, palette=custom_palette)
axes[1, 0].set_xticklabels(group_names, rotation=0)
axes[1, 0].set_ylim(-100, 180)

# Num blinks
sns.boxplot(data=boxplot_data_num_blinks, ax=axes[1, 1], showfliers=False, palette=custom_palette)
axes[1, 1].set_xticklabels(group_names, rotation=0)
axes[1, 1].set_ylim(-50, 70)

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()