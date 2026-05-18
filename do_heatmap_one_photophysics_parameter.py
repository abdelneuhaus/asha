from PAPER.utils import read_poca_files
from preprocessing import pad_list, fusion, fusion_position, photon_calculation, loc_prec_calculation, pre_process_off_frame_csv, pre_process_on_frame_csv, pre_process_sigma, pre_process_single_intensity, get_and_save_well_and_FOV

import os
import itertools
import seaborn as sns
import pandas as pd
import numpy as np
import statistics
import matplotlib.pyplot as plt
import numpy as np

from scipy.stats import f_oneway, mannwhitneyu


def add_statistical_annotations(ax, data, significance_level=0.05):
    pairs = list(itertools.combinations(range(len(data)), 2))
    y_max = max(max(sublist) for sublist in data)
    height = y_max * 0.05  # Adjust height of the significance marker

    for pair in pairs:
        group1 = data[pair[0]]
        group2 = data[pair[1]]
        
        # Perform Mann-Whitney U test
        stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
        
        if p < significance_level:
            x1, x2 = pair
            y, h = y_max + height, height

            # Determine the significance level
            if p < 0.001:
                significance = '***'
            elif p < 0.01:
                significance = '**'
            elif p < 0.05:
                significance = '*'
            else:
                significance = ''

            # Plot significance line and asterisk
            if significance:
                ax.plot([x1 + 1, x1 + 1, x2 + 1, x2 + 1], [y, y + h, y + h, y], lw=1.5, c='black')
                ax.text((x1 + x2 + 2) * .5, y + h, significance, ha='center', va='bottom', color='black')

            # Increase y_max for the next line
            y_max += height * 2

def do_heatmap_one_photophysics_parameter(exp, index, list_of_poca_files, list_of_frame_csv, list_of_int_csv, list_of_sigma_csv,
                                          isPT=True, stats=statistics.median, get_boxplot=False):
    csv_frame_label  = ['ON times', "OFF times"]
    csv_int_label =  "Intensity_loc"
    csv_sigma_label = "Loc_Precision"
    idx = ['1', '2', '3', '4']
    cols = ['A', 'B']
    all_wells = fusion(cols, idx)
    
    for i in index:
        heatmap_data = []
        boxplot_data = []
        # Case where index is 'ON times' or 'OFF times'
        if i in csv_frame_label:
            if i == 'ON times':
                for f in range(len(list_of_frame_csv)):
                    heatmap_data.append(int(stats(pre_process_on_frame_csv(list_of_frame_csv[f]))))
                    boxplot_data.append(pre_process_on_frame_csv(list_of_frame_csv[f]))
            else:
                for f in range(len(list_of_frame_csv)):
                    heatmap_data.append(int(stats(pre_process_off_frame_csv(list_of_frame_csv[f]))))
                    boxplot_data.append(pre_process_off_frame_csv(list_of_frame_csv[f]))
        # Case where index is 'intensity per loc'
        elif i == csv_int_label:
            for f in range(len(list_of_int_csv)):
                heatmap_data.append(int(stats(photon_calculation(pre_process_single_intensity(list_of_int_csv[f])))))
                boxplot_data.append(photon_calculation(pre_process_single_intensity(list_of_int_csv[f])))
        # Case where we compute localization precision       
        elif i == csv_sigma_label:
            for f in range(len(list_of_sigma_csv)):
                heatmap_data.append(float(stats(loc_prec_calculation(pre_process_sigma(list_of_sigma_csv[f]), photon_calculation(pre_process_single_intensity(list_of_int_csv[f]))))))
                boxplot_data.append(loc_prec_calculation(pre_process_sigma(list_of_sigma_csv[f]), photon_calculation(pre_process_single_intensity(list_of_int_csv[f]))))
        # Case where we read from locPALMTracer_cleaned file
        else:
            for f in range(len(list_of_poca_files)):
                raw_file_poca = read_poca_files(list_of_poca_files[f])
                
                if i == 'intensity':
                    heatmap_data.append(int(stats(photon_calculation((raw_file_poca.loc[:, i].values.tolist())))))
                    boxplot_data.append(photon_calculation((raw_file_poca.loc[:, i].values.tolist())))
                else:
                    heatmap_data.append(int(stats(raw_file_poca.loc[:, i].values.tolist())))
                    boxplot_data.append(stats(raw_file_poca.loc[:, i].values.tolist()))
        
        duty_cycle = []
        if get_boxplot == True:
            for f in range(len(list_of_poca_files)):
                raw_file_poca = read_poca_files(list_of_poca_files[f])
                duty_cycle.append(np.array(raw_file_poca.loc[:, 'total ON'])/np.array(raw_file_poca.loc[:, 'lifetime']))
        # We initialize well names
        well_name = []
        if isPT == True:
            for d in range(len(list_of_poca_files)):
                name = get_and_save_well_and_FOV(list_of_poca_files[d],'/561.PT/locPALMTracer_merged.txt', '/561-405.PT/locPALMTracer_merged.txt')
                well_name.append(name)
            legend = fusion_position(all_wells, well_name)
        else:
            legend = list()
            for d in list_of_poca_files:
                well_name.append(os.path.basename(os.path.normpath(d.replace('.PT/locPALMTracer_merged.txt', ''))))
                legend.append(os.path.basename(os.path.normpath(d.replace('.PT/locPALMTracer_merged.txt', ''))))
            # heatmap_data = pad_list(heatmap_data)

        # Rotation of the 8x1 data to 2x4 and plot it on the heatmap
        heatmap_data = pad_list(heatmap_data)
        df = pd.DataFrame(np.array(heatmap_data).reshape(2,4), index=cols, columns=idx)
        sns.heatmap(df, annot=True, fmt='g', cmap="YlGnBu", linewidths=1, linecolor='black')
        plt.yticks(rotation=0)

        # Récupération de la position de l'échelle de couleur
        colorbar = plt.gcf().axes[-1]
        colorbar_pos = colorbar.get_position()

        # Ajout de la boîte de texte à droite de l'échelle de couleur
        plt.text(colorbar_pos.x1 + 0.05, colorbar_pos.y1, "\n".join(legend),
                transform=plt.gcf().transFigure, fontsize=10,
                verticalalignment='top', horizontalalignment='left',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        plt.gcf().set_size_inches((12, 5))
        plt.title(i + ' median')        

        # Save figure
        results_dir = os.path.join('results/'+exp+'/')

        sample_file = i+'.pdf'
        sample_file = sample_file.replace('.PT', '')
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)
        plt.savefig(results_dir+sample_file)
        
        fig, ax = plt.subplots()
        boxplot = ax.boxplot(boxplot_data, showfliers=False)
        # ax.set_xticklabels(well_name, rotation=45, ha='right')
        
        # Add statistical annotations
        # add_statistical_annotations(ax, boxplot_data)
        
        plt.show()

        fig, ax = plt.subplots()
        boxplot = ax.boxplot(duty_cycle, showfliers=False)
        ax.set_xticklabels(well_name, rotation=45, ha='right')
        plt.title('duty cycle boxplots (no outliers)')        
        plt.show()
        
        plt.close('all')
