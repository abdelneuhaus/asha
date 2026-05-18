from PAPER.utils import read_poca_files
from preprocessing import creer_matrice_et_moyennes, creer_data_boxplot, convertir_coordonnees, pad_list, fusion, fusion_position, photon_calculation, loc_prec_calculation, pre_process_off_frame_csv, pre_process_on_frame_csv, pre_process_sigma, pre_process_single_intensity, get_and_save_well_and_FOV
from statistical_test_same_well import statistical_test_same_well

import os
import seaborn as sns
import pandas as pd
import numpy as np
import statistics
import matplotlib.pyplot as plt
import numpy as np


def do_96heatmap_one_photophysics_parameter(exp, index, list_of_poca_files, list_of_frame_csv, list_of_int_csv, list_of_sigma_csv,
                                          isPT=True, stats=statistics.mean, get_boxplot=True):
    csv_frame_label  = ['ON times', "OFF times"]
    csv_int_label =  "Intensity_loc"
    csv_sigma_label = "Loc_Precision"
    
    idx = list(map(str, range(1, 13)))
    cols = list(map(str, [chr(i) for i in range(ord('A'), ord('H')+1)]))
    all_wells = fusion(cols, idx)

    for i in index:
        heatmap_data = []
        # Case where index is 'ON times' or 'OFF times'
        if i in csv_frame_label:
            if i == 'ON times':
                for f in range(len(list_of_frame_csv)):
                    heatmap_data.append(float(stats(pre_process_on_frame_csv(list_of_frame_csv[f]))))
            else:
                for f in range(len(list_of_frame_csv)):
                    heatmap_data.append(float(stats(pre_process_off_frame_csv(list_of_frame_csv[f]))))
        # Case where index is 'intensity per loc'
        elif i == csv_int_label:
            for f in range(len(list_of_int_csv)):
                raw_file_poca = read_poca_files(list_of_poca_files[f])
                filtered = raw_file_poca[raw_file_poca['blinks'] < np.quantile(raw_file_poca["blinks"], 1)]
                filtered = filtered[filtered['total ON'] < np.quantile(raw_file_poca["total ON"], 1)]
                heatmap_data.append(int(stats(np.array(filtered['intensity'])/np.array(filtered['total ON'])))*(3.6/300))

        # Case where we compute localization precision       
        elif i == csv_sigma_label:
            for f in range(len(list_of_sigma_csv)):
                heatmap_data.append(float(160*stats(loc_prec_calculation(pre_process_sigma(list_of_sigma_csv[f]), photon_calculation(pre_process_single_intensity(list_of_int_csv[f]))))))
        # Case where we read from locPALMTracer_merged file
        else:
            for f in range(len(list_of_poca_files)):
                raw_file_poca = read_poca_files(list_of_poca_files[f])
                raw_file_poca = raw_file_poca[raw_file_poca['blinks'] < np.quantile(raw_file_poca["blinks"], 1)]
                raw_file_poca = raw_file_poca[raw_file_poca['total ON'] < np.quantile(raw_file_poca["total ON"], 1)]
                if i == 'intensity':
                    heatmap_data.append(int(stats(photon_calculation((raw_file_poca.loc[:, i].values.tolist())))))
                else:
                    heatmap_data.append(float(stats(raw_file_poca.loc[:, i].values.tolist())))
                    
        if get_boxplot == True:
            for f in range(len(list_of_poca_files)):
                raw_file_poca = read_poca_files(list_of_poca_files[f])

            # We initialize well names
            well_name = []
            if isPT == True:
                for d in range(len(list_of_poca_files)):
                    name = get_and_save_well_and_FOV(list_of_poca_files[d],'/SR_001.MIA/locPALMTracer_merged.txt')
                    well_name.append(name)
            else:
                for d in list_of_poca_files:
                    well_name.append(os.path.basename(os.path.normpath(d.replace('/locPALMTracer_merged.txt', ''))))
            
            matrice_resultante = creer_matrice_et_moyennes(list_of_poca_files, heatmap_data, all_wells)
            df = pd.DataFrame(np.array(matrice_resultante).reshape(8,12), index=cols, columns=idx)

            sns.heatmap(df, annot=True, fmt='g', cmap="coolwarm", linewidths=1, linecolor='black')
            plt.yticks(rotation=0)
            plt.gcf().set_size_inches((12, 5))
            plt.title(i)    
            # Sauvegarder la figure
            results_dir = os.path.join('results/'+exp+'/')
            sample_file = i+'.pdf'
            sample_file = sample_file.replace('.PT', '')
            if not os.path.isdir(results_dir):
                os.makedirs(results_dir)
            plt.savefig(results_dir+sample_file)
            plt.show()
            plt.close('all')
