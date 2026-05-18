from utils import read_poca_files
from preprocessing import get_and_save_well_and_FOV

import matplotlib.pyplot as plt
import os


def do_cumulative_number_clusters(list_of_poca_files, exp, isPT=True):
    tmp = list()
    for i in list_of_poca_files:
        raw_file_poca = read_poca_files(i)
        loc_per_frame = raw_file_poca.groupby(['frame']).size()
        cum_loc_per_frame = loc_per_frame.cumsum()
        tmp.append(cum_loc_per_frame)
 
    fig = plt.figure()
    fig.tight_layout()
    ax = plt.subplot(111)
    for i in tmp:
        i.plot(ax=ax)
    plt.xlabel('Time (in frames)')
    plt.ylabel('Cumulative Number of Clusters')
    
    legend = list()
    if isPT == True:
        for d in range(len(list_of_poca_files)):
            name = get_and_save_well_and_FOV(list_of_poca_files[d],'/SR_001.MIA/locPALMTracer_merged.txt')
            legend.append(name)
    else:
        for d in list_of_poca_files:
            legend.append(os.path.basename(os.path.normpath(d.replace('/locPALMTracer_merged.txt', ''))))
    
    
    plt.legend(legend, loc='upper left')
    plt.grid(linestyle='-', linewidth=1)
    results_dir = os.path.join('results/'+exp+'/')
    sample_file = 'cumulative_clusters.pdf'
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    plt.savefig(results_dir+sample_file)
    plt.close('all')


