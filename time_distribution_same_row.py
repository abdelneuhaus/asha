import os
import matplotlib.pyplot as plt
import numpy as np
import statistics as statistics

from preprocessing import photon_calculation

# Function to process a list of POCA files and create plots
def time_distribution_same_row(super_table, expname, column, well):
    all_trends_photloc = []
    all_trends_total_on = []
    all_trends_blinks = []
    filenames = []

    for position in super_table:
        all_trends_photloc.append(photon_calculation(position['photon_loc']))
        all_trends_total_on.append(position['total ON'])
        all_trends_blinks.append(position['blinks'])
        filenames.append(well)

    frames = np.array(list(range(0, 5000, 500)))

    # Create the figure with 1x3 subplots
    fig = plt.figure(layout='constrained', figsize=(18, 6))
    fig.patch.set_facecolor('white')
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)

    # Plot all photon_loc trends
    for trend_photloc, filename in zip(all_trends_photloc, well):
        ax1.plot(frames, trend_photloc[:len(frames)], 'o-', label=filename)
    ax1.set_title('Photon Loc Over Time')
    ax1.set_xlabel('Frame Intervals')
    ax1.set_ylabel('Mean Photon per Loc')

    # Plot all total ON time trends
    for trend_total_on, filename in zip(all_trends_total_on, well):
        ax2.plot(frames, trend_total_on[:len(frames)], 'o-', label=filename)
    ax2.set_title('Total ON Time Over Time')
    ax2.set_xlabel('Frame Intervals')
    ax2.set_ylabel('Mean Total ON Time')

    # Plot all intensity trends
    for trend_blinks, filename in zip(all_trends_blinks, well):
        ax3.plot(frames, trend_blinks[:len(frames)], 'o-', label=filename)
    ax3.set_title('nBlinks Over Time')
    ax3.set_xlabel('Frame Intervals')
    ax3.set_ylabel('Mean nBlinks')
    ax3.legend()
    
    # Sauvegarder la figure
    results_dir = os.path.join('results/'+expname+'/')
    sample_file = 'time_distribution_row_'+column+'.pdf'
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    plt.savefig(results_dir+sample_file)
    plt.close(fig)