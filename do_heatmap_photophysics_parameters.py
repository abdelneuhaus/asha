from PAPER.utils import read_poca_files
from preprocessing import photon_calculation, loc_prec_calculation, pre_process_off_frame_csv, pre_process_on_frame_csv, pre_process_sigma, pre_process_single_intensity, get_and_save_well_and_FOV

import os
import seaborn as sns
import numpy as np
import statistics
import pandas as pd 
import matplotlib.pyplot as plt



def do_heatmap_photophysics_parameters(exp, list_of_poca_files, list_of_frame_csv, list_of_int_csv, list_of_sigma_csv, 
                                       isPT=True, stats=statistics.median):
    heatmap_data = []
    cpt = 0
    tmp_pho_loc = list()
    for f in range(len(list_of_poca_files)):
        raw_file_poca = read_poca_files(list_of_poca_files[f])
        raw_file_poca = raw_file_poca[raw_file_poca['blinks'] < np.quantile(raw_file_poca["blinks"], 1)]
        raw_file_poca = raw_file_poca[raw_file_poca['total ON'] < np.quantile(raw_file_poca["total ON"], 1)]
        well_data = dict()
        well_data['_on_times'] = int(stats(pre_process_on_frame_csv(list_of_frame_csv[f])))
        well_data['_off_times'] = int(stats(pre_process_off_frame_csv(list_of_frame_csv[f])))
        tmp = photon_calculation(pre_process_single_intensity(list_of_int_csv[f]))
        tmp_pho_loc.append(tmp)
        well_data['_phot_per_loc'] = int(stats(photon_calculation(pre_process_single_intensity(list_of_int_csv[f]))))
        well_data['_total_on'] = int(stats(raw_file_poca.loc[:, 'total ON'].values.tolist()))
        well_data['_num_blinks'] = int(stats(raw_file_poca.loc[:, 'blinks'].values.tolist()))
        well_data['_phot_per_cluster'] = int(stats(photon_calculation(raw_file_poca.loc[:, 'intensity'].values.tolist())))
        well_data['_num_on_times'] = int(stats(raw_file_poca.loc[:, '# seq ON'].values.tolist()))
        well_data['_num_off_times'] = int(stats(raw_file_poca.loc[:, '# seq OFF'].values.tolist()))
        # well_data['_sigma'] = float(stats(loc_prec_calculation(pre_process_sigma(list_of_sigma_csv[f], tmp_pho_loc[cpt]))))
        cpt += 1
        heatmap_data.append(well_data)

    # Convert dict to pd.dataframe
    data = []
    if isPT == True:
        for d in range(len(heatmap_data)):
            name = get_and_save_well_and_FOV(list_of_poca_files[d],'/SR_001.MIA/locPALMTracer_merged.txt')
            b=pd.DataFrame.from_dict(heatmap_data[d], orient='index').rename({0:name}, axis='columns')
            data.append(b)
    else:
        name = list()
        for d in range(len(heatmap_data)):
            name = os.path.basename(os.path.normpath(list_of_poca_files[d].replace('/locPALMTracer_merged.txt', '')))
            b=pd.DataFrame.from_dict(heatmap_data[d], orient='index').rename({0:name}, axis='columns')
            data.append(b)
            
    # Create figure and convert dataframe to heatmap data
    fig, ax = plt.subplots(figsize=(14, 7))
    heatmap_data = pd.concat(data)
    heatmap_data = heatmap_data.groupby(heatmap_data.index).sum()
    heatmap_data = heatmap_data.replace(np.nan, 0)
    
    # Normalise les MOYENNES (ON PEUT FAIRE LES MEDIANES) car sinon on a un trop gros écart sur une même plaque
    tab_n = heatmap_data.div(heatmap_data.max(axis=1), axis=0)
    tab_orig = heatmap_data.copy()
    sns.heatmap(tab_n, annot=tab_orig, cmap="mako",  fmt=".1f", vmin=0, vmax=1)
    cbar = ax.collections[0].colorbar
    cbar.set_ticks([0, 0.5, 1])
    cbar.ax.set_yticklabels(['Low', 'Median', 'High'])    
    ax.xaxis.tick_top()

    # Columns & Rows labels
    column_labels = heatmap_data.columns
    row_labels = ["Number of blinks", "Number OFF times", "Number ON times", 
                "Length OFF times (frames)", "Length ON times (frames)", "Photons/Molecule", 
                "Photons/Localization", "Total ON time (frames)"]
    ax.set_yticklabels(row_labels, minor=False)
    plt.show()
    # Save figure
    results_dir = os.path.join('results/'+exp+'/')
    sample_file = 'experiment_heatmap.pdf'
    sample_file = sample_file.replace('.PT', '')
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    plt.savefig(results_dir+sample_file)
    plt.close('all')
    
