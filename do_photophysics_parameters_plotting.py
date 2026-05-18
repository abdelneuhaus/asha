from PAPER.utils import read_poca_files
from preprocessing import photon_calculation, loc_prec_calculation, pre_process_off_frame_csv, pre_process_on_frame_csv, pre_process_sigma, pre_process_single_intensity

import matplotlib.pyplot as plt
import os
import statistics
import numpy as np


def do_photophysics_parameters_plotting(list_of_poca_files, list_of_frame_csv, list_of_int_csv, list_of_sigma_csv, exp=None, isPT=True,
                                        on_times=True, off_times=True, total_on=True, num_blinks=True,
                                        phot_per_loc=True, phot_per_cluster=True, num_on_times=True, 
                                        num_off_times=True, sigma=True, boxplot=True):
    tmp_pho_loc = list()
    cpt = 0
    for j in range(len(list_of_frame_csv)):
        # length of each ON time
        _on_times = pre_process_on_frame_csv(list_of_frame_csv[j]) if on_times else None
        # length of each OFF time
        _off_times = pre_process_off_frame_csv(list_of_frame_csv[j]) if off_times else None
        # photon per localization
        _phot_per_loc = photon_calculation(pre_process_single_intensity(list_of_int_csv[j])) if phot_per_loc else None
        tmp_pho_loc.append(photon_calculation(pre_process_single_intensity(list_of_int_csv[j])))
        
        raw_file_poca = read_poca_files(list_of_poca_files[j])
        raw_file_poca = raw_file_poca[raw_file_poca['blinks'] < np.quantile(raw_file_poca["blinks"], 1)]
        raw_file_poca = raw_file_poca[raw_file_poca['total ON'] < np.quantile(raw_file_poca["total ON"], 1)]

        # = bleachtime or total ON in frame number
        _total_on = raw_file_poca.loc[:, 'total ON'].values.tolist() if total_on else None
        # num blinks
        _num_blinks = raw_file_poca.loc[:, 'blinks'].values.tolist() if num_blinks else None
        # photon per cluster
        _phot_per_cluster = photon_calculation(raw_file_poca.loc[:, 'intensity'].values.tolist()) if phot_per_cluster else None
        # number of ON times per cluster
        _num_on_times = raw_file_poca.loc[:, '# seq ON'].values.tolist() if num_on_times else None
        # number of OFF times per cluster
        _num_off_times = raw_file_poca.loc[:, '# seq OFF'].values.tolist() if num_off_times else None
        _sigma = loc_prec_calculation(pre_process_sigma(list_of_sigma_csv[j]), tmp_pho_loc[cpt]) if sigma else None
        cpt += 1
        
        non_none_elements = [_on_times, _off_times, _total_on, _num_blinks, _phot_per_loc, _phot_per_cluster, _num_on_times, _num_off_times, _sigma]
        non_none_elements = [elem for elem in non_none_elements if elem is not None]
        label = ['ON times', "OFF times", "Total ON (Bleachtime)", "#Blinks", "Photon/loc", "Photon/cluster", "#ON times", "#OFF times", "Loc. Precision"]
        
        fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(15,8))
        plt.subplots_adjust(wspace=0.3, hspace=0.5)
        for m, df in enumerate(non_none_elements):
            row = m // 4
            col = m % 4
            ax = axes[row][col]
            ax.set_title(label[m]+', median:'+ str(np.round(statistics.median(df), 2)))
            
            if boxplot ==  True:
            # Boxplot version
                ax.boxplot(df, showfliers=False)

            # Histogram version
            elif boxplot == False:
                median = np.median(df)
                left = np.percentile(df, 5)
                right = np.percentile(df, 95)
                cropped_data = [x for x in df if left <= x <= right]
                ax.hist(cropped_data, bins=50, density=True)
                ax.axvline(median, color='red', linestyle='dashed', linewidth=1)
                                
        title_plot = os.sep.join(list_of_poca_files[j].replace('/SR_001.MIA/locPALMTracer_merged.txt', '').split('/')[-2:])
        results_dir = os.path.join('results/'+exp+'/'+title_plot+'/')
        plt.suptitle(title_plot, fontsize=14)

        sample_file = '/photophysics_plots.pdf'
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)
        plt.savefig(results_dir+sample_file)
        plt.close('all')