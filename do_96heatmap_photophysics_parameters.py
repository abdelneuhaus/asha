from PAPER.utils import read_poca_files
from preprocessing import photon_calculation, loc_prec_calculation, pre_process_off_frame_csv, pre_process_on_frame_csv, pre_process_sigma, pre_process_single_intensity, get_and_save_well_and_FOV

import os
import seaborn as sns
import numpy as np
import statistics
import pandas as pd 
import matplotlib.pyplot as plt


def do_96heatmap_photophysics_parameters(exp, list_of_poca_files, list_of_frame_csv, list_of_int_csv, list_of_sigma_csv, 
                                       col, isPT=True, stats=statistics.median):
    heatmap_data = []
    for f in range(len(list_of_poca_files)):
        raw_file_poca = read_poca_files(list_of_poca_files[f])
        raw_file_poca = raw_file_poca[raw_file_poca['blinks'] < np.quantile(raw_file_poca["blinks"], 1)]
        raw_file_poca = raw_file_poca[raw_file_poca['total ON'] < np.quantile(raw_file_poca["total ON"], 1)]

        well_data = dict()
        on_times = pre_process_on_frame_csv(list_of_frame_csv[f])
        off_times = pre_process_off_frame_csv(list_of_frame_csv[f])
        well_data['_on_times'] = float(stats(on_times))
        well_data['_off_times'] = float(stats(off_times))
        well_data['_phot_per_loc'] = float(stats(np.array(raw_file_poca['intensity']) / np.array(raw_file_poca['total ON']))* (3.6/300)) 
        well_data['_total_on'] = float(stats(raw_file_poca['total ON'].values))
        well_data['_num_blinks'] = float(stats(raw_file_poca['blinks'].values))
        well_data['_phot_per_cluster'] = float(stats(photon_calculation(raw_file_poca['intensity'].values.tolist())))
        well_data['_num_on_times'] = float(stats(raw_file_poca['# seq ON'].values.tolist()))
        well_data['_num_off_times'] = float(stats(raw_file_poca['# seq OFF'].values.tolist()))

        heatmap_data.append(well_data)

    # Convert to DataFrame and group by columns
    data = []
    if isPT:
        for d in range(len(heatmap_data)):
            name = get_and_save_well_and_FOV(list_of_poca_files[d], '/SR_001.MIA/locPALMTracer_merged.txt')
            b = pd.DataFrame.from_dict(heatmap_data[d], orient='index').rename({0: name}, axis='columns')
            data.append(b)
    else:
        for d in range(len(heatmap_data)):
            name = os.path.basename(os.path.normpath(list_of_poca_files[d].replace('/locPALMTracer_merged.txt', '')))
            b = pd.DataFrame.from_dict(heatmap_data[d], orient='index').rename({0: name}, axis='columns')
            data.append(b)

    heatmap_data = pd.concat(data)
    heatmap_data = heatmap_data.groupby(heatmap_data.index).mean()
    # heatmap_data = heatmap_data.replace(np.nan, 0)
    heatmap_data.columns = heatmap_data.columns.str.replace(r'\\', '/', regex=True)
    
    # Group data by column prefixes
    grouped_data = pd.DataFrame()
    for prefix in sorted(set(col.split('/')[0] for col in heatmap_data.columns)):
        cols_to_group = [col for col in heatmap_data.columns if col.startswith(prefix)]
        grouped_data[prefix] = heatmap_data[cols_to_group].mean(axis=1)
    
    # Normalize data
    tab_n = grouped_data.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=1)
    tab_orig = grouped_data.copy()

    fig, ax = plt.subplots(figsize=(14, 7))
    # Change the colormap to something with more nuance
    sns.heatmap(tab_n, annot=tab_orig, cmap="coolwarm", fmt=".1f", vmin=0, vmax=1)
    cbar = ax.collections[0].colorbar
    cbar.set_ticks([0, 0.5, 1])
    cbar.ax.set_yticklabels(['Low', 'Median', 'High'])
    ax.xaxis.tick_top()

    column_labels = grouped_data.columns
    row_labels = ["Number of blinks", "Number OFF times", "Number ON times", 
                "Length OFF times (frames)", "Length ON times (frames)", "Photons/Molecule", 
                "Photons/Localization", "Total ON time (frames)"]

    ax.set_xticklabels(column_labels, minor=False)
    ax.set_yticklabels(row_labels, minor=False)

    results_dir = os.path.join('results/'+exp+'/')
    sample_file = 'experiment_heatmap_row_'+ col +'.pdf'
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    plt.savefig(results_dir+sample_file)
    plt.show()
    plt.close('all')