from PAPER.utils import read_poca_files
from preprocessing import fusion, photon_calculation, loc_prec_calculation, pre_process_off_frame_csv, pre_process_on_frame_csv, pre_process_sigma, pre_process_single_intensity, get_and_save_well_and_FOV
from statistical_test_same_well import statistical_test_same_well

import os

import numpy as np
import matplotlib.pyplot as plt
import numpy as np


def do_stat_tests_fovs(exp, well, index, list_of_poca_files, list_of_frame_csv, list_of_int_csv, list_of_sigma_csv):
    csv_frame_label  = ['ON times', "OFF times"]
    csv_int_label =  "Intensity_loc"
    csv_sigma_label = "Loc_Precision"
    
    idx = list(map(str, range(1, 13)))
    cols = list(map(str, [chr(i) for i in range(ord('A'), ord('H')+1)]))
    all_wells = fusion(cols, idx)

    for i in index:
        boxplot_data = []
        # Case where index is 'ON times' or 'OFF times'
        if i in csv_frame_label:
            if i == 'ON times':
                for f in range(len(list_of_frame_csv)):
                    boxplot_data.append(pre_process_on_frame_csv(list_of_frame_csv[f]))
            else:
                for f in range(len(list_of_frame_csv)):
                    boxplot_data.append(pre_process_off_frame_csv(list_of_frame_csv[f]))
        # Case where index is 'intensity per loc'
        elif i == csv_int_label:
            for f in range(len(list_of_int_csv)):
                raw_file_poca = read_poca_files(list_of_poca_files[f])
                filtered = raw_file_poca[raw_file_poca['blinks'] < np.quantile(raw_file_poca["blinks"], 1)]
                filtered = filtered[filtered['total ON'] < np.quantile(filtered["total ON"], 1)]
                boxplot_data.append(np.array(filtered['intensity'])/np.array(filtered['total ON'])*(3.6/300))
        # Case where we compute localization precision       
        elif i == csv_sigma_label:
            for f in range(len(list_of_sigma_csv)):
                boxplot_data.append(pre_process_sigma(list_of_sigma_csv[f]))
        # Case where we read from locPALMTracer_merged file
        else:
            for f in range(len(list_of_poca_files)):
                raw_file_poca = read_poca_files(list_of_poca_files[f])
                raw_file_poca = raw_file_poca[raw_file_poca['blinks'] < np.quantile(raw_file_poca["blinks"], 1)]
                raw_file_poca = raw_file_poca[raw_file_poca['total ON'] < np.quantile(raw_file_poca["total ON"], 1)]
                if i == 'intensity':
                    boxplot_data.append(photon_calculation((raw_file_poca.loc[:, i].values.tolist())))
                else:
                    boxplot_data.append(raw_file_poca.loc[:, i].values.tolist())
                    
        fig, ax = plt.subplots()
        boxplot = ax.boxplot(boxplot_data, showfliers=False)
        statistical_test_same_well(ax, boxplot_data)
        results_dir = os.path.join('results/'+exp+'/'+well+'/')
        sample_file = i+'_stats.pdf'
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)
        fig.suptitle(well+'_'+i+'_stats')
        plt.savefig(results_dir+sample_file)
        plt.show()
        plt.close('all')