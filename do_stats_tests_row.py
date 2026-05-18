from PAPER.utils import read_poca_files
from preprocessing import fusion, photon_calculation, loc_prec_calculation, pre_process_off_frame_csv, pre_process_on_frame_csv, pre_process_sigma, pre_process_single_intensity, get_and_save_well_and_FOV
from statistical_test_same_well import statistical_test_same_well

import os
import numpy as np
import matplotlib.pyplot as plt
import numpy as np


import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

def do_stats_tests_row(exp, col, index, wells_data):
    csv_frame_label = ['ON times', "OFF times"]
    csv_int_label = "Intensity_loc"
    csv_sigma_label = "Loc_Precision"
    
    cmap = sns.color_palette("coolwarm", as_cmap=True)  # Utiliser la palette de couleurs coolwarm
    
    for i in index:
        boxplot_data = []
        well_names = []
        mean_values = []  # Stocker les valeurs moyennes pour la coloration
        
        for well, data in wells_data.items():
            well_poca_data = data["poca_data"]
            well_frame_data = data["frame_data"]
            well_intensity_data = data["intensity_data"]
            well_sigma_data = data["sigma_data"]
            
            # Cas où l'index est 'ON times' ou 'OFF times'
            if i in csv_frame_label:
                for f in well_frame_data:
                    if i == 'ON times':
                        processed_data = pre_process_on_frame_csv(f)
                    else:
                        processed_data = pre_process_off_frame_csv(f)
                    boxplot_data.append(processed_data)
                    mean_values.append(np.mean(processed_data))  # Calcul de la valeur moyenne
                    
            # Cas où l'index est 'intensity per loc'
            elif i == csv_int_label:
                filtered = well_poca_data[well_poca_data['blinks'] < np.quantile(well_poca_data["blinks"], 1)]
                filtered = filtered[filtered['total ON'] < np.quantile(filtered["total ON"], 1)]
                processed_data = np.array(filtered['intensity']) / np.array(filtered['total ON']) * (3.6 / 300)
                boxplot_data.append(processed_data)
                mean_values.append(np.mean(processed_data))  # Calcul de la valeur moyenne
                
            # Cas où on calcule la précision de localisation       
            elif i == csv_sigma_label:
                for f, int_f in zip(well_sigma_data, well_intensity_data):
                    processed_data = 160 * loc_prec_calculation(pre_process_sigma(f), photon_calculation(pre_process_single_intensity(int_f)))
                    boxplot_data.append(processed_data)
                    mean_values.append(np.mean(processed_data))  # Calcul de la valeur moyenne
                    
            # Cas où on lit depuis le fichier locPALMTracer_merged
            else:
                filtered_poca_data = well_poca_data[well_poca_data['blinks'] < np.quantile(well_poca_data["blinks"], 1)]
                filtered_poca_data = filtered_poca_data[filtered_poca_data['total ON'] < np.quantile(filtered_poca_data["total ON"], 1)]
                if i == 'intensity':
                    processed_data = photon_calculation(filtered_poca_data.loc[:, i].values.tolist())
                else:
                    processed_data = filtered_poca_data.loc[:, i].values.tolist()
                boxplot_data.append(processed_data)
                mean_values.append(np.mean(processed_data))  # Calcul de la valeur moyenne
            
            well_names.append(well)

        # Normaliser les valeurs moyennes pour qu'elles soient comprises entre 0 et 1
        norm_mean_values = (np.array(mean_values) - np.min(mean_values)) / (np.max(mean_values) - np.min(mean_values))
        
        fig, ax = plt.subplots()
        boxprops = dict(linestyle='-', linewidth=1.5, color='black')  # Propriétés de contour des boîtes
        bp = ax.boxplot(boxplot_data, showfliers=False, labels=well_names, patch_artist=True, boxprops=boxprops)
        
        # Appliquer les couleurs aux boxplots
        for patch, norm_value in zip(bp['boxes'], norm_mean_values):
            color = cmap(norm_value)  # Obtenir la couleur à partir de la palette
            patch.set_facecolor(color)
        ax.set_ylim(0, 200)
        # statistical_test_same_well(ax, boxplot_data)
        results_dir = os.path.join('results/' + exp + '/' + col + '/')
        sample_file = i + '_stats.pdf'
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)
        fig.suptitle(col + '_' + i + '_stats')
        plt.savefig(results_dir + sample_file)
        plt.show()
        plt.close('all')

