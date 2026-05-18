import numpy as np
import matplotlib.pyplot as plt

from setup_plate import setup_plate
from asha.utils import read_poca_files, read_locPALMTracer_file, create_database


def generate_data_plate_boxplot(plate_path, protein, parameter, prefix):
    database = create_database()
    group_poca_files, group_palm_tracer_files = setup_plate(plate_path, prefix)
    protein_index = database[protein][1]
    group_poca_files = group_poca_files[protein_index]
    group_palm_tracer_files = group_palm_tracer_files[protein_index]
    boxplot_data = []
    
    for i in range(len(group_palm_tracer_files)):
        raw_file_poca = read_poca_files(group_poca_files[i])
        raw_pt_file = read_locPALMTracer_file(group_palm_tracer_files[i])
        raw_pt_file['Integrated_Intensity'] *= 3.6/300
        raw_file_poca['intensity'] *= 3.6/300
        raw_file_poca = raw_file_poca[raw_file_poca['intensity'] <= np.quantile(raw_file_poca["intensity"], 1)]
        raw_file_poca = raw_file_poca[raw_file_poca['total ON'] <= np.quantile(raw_file_poca["total ON"], 1)]
        raw_file_poca["avg_on"] = np.array(raw_file_poca['total ON'] / raw_file_poca['# seq ON'])
        raw_file_poca["photon_loc"] = np.array(raw_file_poca['intensity'] / raw_file_poca['total ON'])

        if parameter == "Integrated_Intensity":
            boxplot_data.append(raw_pt_file["Integrated_Intensity"])
        else:
            boxplot_data.append(raw_file_poca[parameter])
    boxplot_data = np.concatenate(boxplot_data)
    return boxplot_data



def boxplots_between_plate(boxplots):
    fig, ax = plt.subplots()
    bp = ax.boxplot(boxplots, patch_artist=True, showfliers=False)
    base_colors = ['orangered', 'dodgerblue', 'palegreen', 'gold', 'mediumorchid']
    num_groups = len(boxplots)
    colors = base_colors[:num_groups]

    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)

    for median in bp['medians']:
        median.set_color('white')
        median.set_linewidth(1)
        
    ax.set_xticklabels([f'Plate {i+1}' for i in range(num_groups)])

    plt.tight_layout()
    plt.show()

