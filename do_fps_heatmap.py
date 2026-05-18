from PAPER.utils import read_poca_files
from preprocessing import photon_calculation, loc_prec_calculation, pre_process_off_frame_csv, pre_process_on_frame_csv, pre_process_sigma, pre_process_single_intensity, get_and_save_well_and_FOV

import os
import seaborn as sns
import numpy as np
import statistics
import pandas as pd 
import matplotlib.pyplot as plt


def do_fps_heatmap(list_of_poca_files, list_of_frame_csv, stats=statistics.mean):
    heatmap_data = []

    # Calcul des statistiques pour chaque fichier
    for f in range(len(list_of_poca_files)):
        raw_file_poca = read_poca_files(list_of_poca_files[f])
        raw_file_poca = raw_file_poca[raw_file_poca['blinks'] <  np.quantile(raw_file_poca['blinks'], 1)]
        raw_file_poca = raw_file_poca[raw_file_poca['total ON'] < np.quantile(raw_file_poca['total ON'], 1)]
        
        well_data = dict()
        on_times = pre_process_on_frame_csv(list_of_frame_csv[f])
        off_times = pre_process_off_frame_csv(list_of_frame_csv[f])
        well_data['_num_blinks'] = float(stats(raw_file_poca['blinks'].values))
        well_data['_num_on_times'] = float(stats(raw_file_poca['# seq ON'].values.tolist()))
        well_data['_num_off_times'] = float(stats(raw_file_poca['# seq OFF'].values.tolist()))
        well_data['_on_times'] = float(stats(on_times))
        well_data['_off_times'] = float(stats(off_times))
        well_data['_phot_per_cluster'] = float(stats(photon_calculation(raw_file_poca['intensity'].values.tolist())))
        well_data['_phot_per_loc'] = float(stats(np.array(raw_file_poca['intensity']) / np.array(raw_file_poca['total ON'])) * (3.6/300)) 
        well_data['_total_on'] = float(statistics.mean(raw_file_poca['total ON'].values))
        # Calcul du duty cycle
        # total_on = raw_file_poca['total ON'].values
        # total_off = raw_file_poca['total OFF'].values
        # duty_cycle = total_on / (total_on + total_off)
        heatmap_data.append(well_data)
    
    protein_data = pd.DataFrame(heatmap_data).mean()
    return protein_data  


def generate_fps_heatmap(all_protein_data, exp_name):
    # Convertir les données en DataFrame
    df = pd.DataFrame(all_protein_data)#.T  # Transposer pour avoir les protéines en colonnes

    # Normalisation des données pour la heatmap
    tab_n = df.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=1)
    tab_orig = df.copy()

    # Génération de la heatmap
    fig, ax = plt.subplots(figsize=(14, 7))
    sns.heatmap(tab_n, annot=tab_orig, cmap="coolwarm", fmt=".1f", vmin=0, vmax=1)
    cbar = ax.collections[0].colorbar
    cbar.set_ticks([0, 0.5, 1])
    cbar.ax.set_yticklabels(['Low', 'Median', 'High'])
    ax.xaxis.tick_top()

    # Labels des lignes (statistiques) et colonnes (protéines)
    row_labels = ["Number of blinks", "Number ON times", "Number OFF times", 
                  "Length ON times (frames)", "Length OFF times (frames)", 
                  "Photons/Molecule", "Photons/Localization", "Total ON time (frames)"]#, "Duty Cycle"]
    
    ax.set_yticklabels(row_labels, minor=False)
    ax.set_xticklabels(all_protein_data.keys())#, rotation=90)

    # Sauvegarde de la heatmap
    results_dir = os.path.join('results', exp_name)
    sample_file = 'experiment_global_heatmap.pdf'
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    plt.savefig(os.path.join(results_dir, sample_file))
    plt.show()
    plt.close('all')
