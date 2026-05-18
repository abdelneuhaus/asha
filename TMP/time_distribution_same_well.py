import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import statistics as statistics

from PAPER.utils import read_poca_files, get_poca_files
from preprocessing import photon_calculation


# Function to process a list of POCA files and create plots
def time_distribution_same_well(pathway, well, expname):
    list_of_poca_files = get_poca_files(pathway)
    
    # Initialize global trends
    all_trends_photloc = []
    all_trends_total_on = []
    all_trends_blinks = []
    filenames = []

    # Iterate over each POCA file
    for poca_file in list_of_poca_files:
        raw_file_poca = read_poca_files(poca_file)
        filtered_data = raw_file_poca[raw_file_poca['total ON'] < np.quantile(raw_file_poca["blinks"], 1)]
        filtered_data = filtered_data[filtered_data['blinks'] < np.quantile(filtered_data["blinks"], 1)]
        filtered_data["avg_on"] = np.array(filtered_data['total ON'] / filtered_data['# seq ON'])
        filtered_data["photon_loc"] = np.array(filtered_data['intensity'] / filtered_data['total ON'])

        grouped_data = filtered_data.groupby(pd.cut(filtered_data['frame'], list(range(0, 5001, 500)))).mean()

        all_trends_photloc.append(photon_calculation(grouped_data['photon_loc']))
        all_trends_total_on.append(grouped_data['total ON'])
        all_trends_blinks.append(grouped_data['blinks'])
        filenames.append(os.sep.join(poca_file.replace('/SR_001.MIA/locPALMTracer_merged.txt', '').split('/')[-2:]))

    frames = np.array(list(range(0, 5000, 500)))

    # Create the figure with 1x3 subplots
    fig = plt.figure(layout='constrained', figsize=(18, 6))
    fig.patch.set_facecolor('white')
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    # Plot all photon_loc trends
    for trend_photloc, filename in zip(all_trends_photloc, filenames):
        ax1.plot(frames, trend_photloc[:len(frames)], 'o-', label=filename)
    ax1.set_title('Photon Loc Over Time')
    ax1.set_xlabel('Frame Intervals')
    ax1.set_ylabel('Mean Photon per Loc')

    # Plot all total ON time trends
    for trend_total_on, filename in zip(all_trends_total_on, filenames):
        ax2.plot(frames, trend_total_on[:len(frames)], 'o-', label=filename)
    ax2.set_title('Total ON Time Over Time')
    ax2.set_xlabel('Frame Intervals')
    ax2.set_ylabel('Mean Total ON Time')

    # Plot all intensity trends
    for trend_blinks, filename in zip(all_trends_blinks, filenames):
        ax3.plot(frames, trend_blinks[:len(frames)], 'o-', label=filename)
    ax3.set_title('nBlinks Over Time')
    ax3.set_xlabel('Frame Intervals')
    ax3.set_ylabel('Mean nBlinks')
    ax3.legend()
    
    # Sauvegarder la figure
    results_dir = os.path.join('results/'+expname+'/'+well+'/')
    sample_file = 'time_distribution_FOV_of_well.pdf'
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    plt.savefig(results_dir+sample_file)
    plt.close(fig)