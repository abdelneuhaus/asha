import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

from boxplots_between_plate import generate_data_plate_boxplot
from utils import read_poca_files, read_locPALMTracer_file, create_database


def well_scale(target_protein, target_property, group_poca_files, group_palm_tracer_files):
    database = create_database()
    target_well = random.choice(database[target_protein][0])
    protein_index = database[target_protein][1]
    target_well_poca_files = [group for group in group_poca_files[protein_index] if target_well in group]
    target_well_palm_tracer_files = [group for group in group_palm_tracer_files[protein_index] if target_well in group]
    poca_data = [read_poca_files(f) for f in target_well_poca_files]
    palm_tracer_data = [read_locPALMTracer_file(f) for f in target_well_palm_tracer_files]
    
    if target_property == 'intensity':
        intensity_values = [df["intensity"]* (3.6/300) for df in poca_data]
        return intensity_values
    elif target_property == 'total ON':
        total_on_values = [df["total ON"] for df in poca_data]
        return total_on_values
    elif target_property == 'blinks':
        blinks_values = [df["blinks"] for df in poca_data]
        return blinks_values
    elif target_property == 'Integrated_Intensity':
        integrated_intensity_values = [df["Integrated_Intensity"]* (3.6/300) for df in palm_tracer_data]
        return integrated_intensity_values
    else:
        raise ValueError("Unknown photophysical parameter: {}".format(target_property))


def replicate_scale(target_protein, target_property, group_poca_files, group_palm_tracer_files):
    database = create_database()
    wells = database[target_protein][0]
    protein_index = database[target_protein][1]
    all_well_data = []
    all_replicate_medians = []
    for well in wells:
        target_well_poca_files = [f for f in group_poca_files[protein_index] if well in f]
        target_well_palm_tracer_files = [f for f in group_palm_tracer_files[protein_index] if well in f]
        poca_data = [read_poca_files(f) for f in target_well_poca_files]
        palm_tracer_data = [read_locPALMTracer_file(f) for f in target_well_palm_tracer_files]

        if target_property == 'intensity':
            values = [df["intensity"] * (3.6 / 300) for df in poca_data]
        elif target_property == 'total ON':
            values = [df["total ON"] for df in poca_data]
        elif target_property == 'blinks':
            values = [df["blinks"] for df in poca_data]
        elif target_property == 'Integrated_Intensity':
            values = [df["Integrated_Intensity"] * (3.6 / 300) for df in palm_tracer_data]
        else:
            raise ValueError("Unknown photophysical parameter: {}".format(target_property))
        
        well_flat_values = pd.concat(values).tolist() if values else []
        all_well_data.append(well_flat_values)
        replicate_medians = [v.median() for v in values if not v.empty]
        all_replicate_medians.append(replicate_medians)
    
    return all_well_data, all_replicate_medians


def experiment_scale(plate_lists, protein, parameter, prefix):
    stockage = []
    for i in plate_lists:
        tmp = generate_data_plate_boxplot(i, protein, parameter, prefix)
        stockage.append(tmp)
    return stockage


def figure2(well_boxplots, replicate_boxplots, replicate_medians, experiment_boxplots):
    fig, axs = plt.subplots(1, 3, figsize=(6, 2), sharex=False, sharey=True, dpi=150)
    # Well scale boxplot
    axs[0].boxplot(well_boxplots, showfliers=False)
    axs[0].set_xticks(range(1, len(well_boxplots) + 1))
    axs[0].set_xticklabels([f'FOV {i+1}' for i in range(len(well_boxplots))], fontsize=8)

    # Replicate scale boxplot
    axs[1].boxplot(replicate_boxplots, showfliers=False)
    axs[1].set_xticks(range(1, len(replicate_boxplots) + 1))
    axs[1].set_xticklabels([f'Well {i+1}' for i in range(len(replicate_boxplots))], fontsize=8)
    for i, medians in enumerate(replicate_medians):
        x = [i + 1] * len(medians)
        axs[1].scatter(x, medians, facecolors='none', edgecolors='black', linewidths=1.2, zorder=3)
        
    # Experiment scale boxplot
    boxprops = axs[2].boxplot(experiment_boxplots, showfliers=False, patch_artist=True)
    colors = plt.cm.tab10.colors
    num_boxes = len(boxprops['boxes'])
    for i in range(num_boxes):
        boxprops['boxes'][i].set_facecolor(colors[i % len(colors)])
        boxprops['medians'][i].set_color('white')
        boxprops['medians'][i].set_linewidth(1)
    axs[2].set_xticks(range(1, len(experiment_boxplots) + 1))
    axs[2].set_xticklabels([f'Plate {i+1}' for i in range(len(experiment_boxplots))], fontsize=8)
    
    plt.tight_layout()
    plt.show()
    
    
    
def save_individual_boxplots(well_boxplots, replicate_boxplots, replicate_medians, experiment_boxplots, save_dir="."):
    # Définir un bleu pour le premier boxplot de la figure 2
    base_blue = '#1f77b4'  # plt default blue

     # === FIGURE 1: WELL BOXPLOTS ===
    fig1, ax1 = plt.subplots(figsize=(2, 2), dpi=150)
    box1 = ax1.boxplot(well_boxplots, showfliers=False, patch_artist=True)

    # Boxplots : fond blanc, bord noir
    for patch in box1['boxes']:
        patch.set_facecolor('white')
        patch.set_edgecolor('black')
    for element in ['whiskers', 'caps', 'medians']:
        for item in box1[element]:
            item.set_color('black')

    # Cadre (axes) en base_blue
    for spine in ax1.spines.values():
        spine.set_color(base_blue)

    ax1.set_xticks(range(1, len(well_boxplots) + 1))
    ax1.set_xticklabels([f'FOV {i+1}' for i in range(len(well_boxplots))], fontsize=8)
    plt.tight_layout()
    fig1.savefig(f"{save_dir}/well_boxplot.png")
    plt.close(fig1)


    # === REPLICATE BOXPLOT ===
    fig2, ax2 = plt.subplots(figsize=(2, 2), dpi=150)

    num_boxes = len(replicate_boxplots)
    blue_cmap = cm.get_cmap('BuPu', num_boxes + 2)
    # Le premier est base_blue, les suivants sont des variantes
    blue_colors = [base_blue] + [mcolors.to_hex(blue_cmap(i)) for i in range(2, num_boxes + 1)]

    box2 = ax2.boxplot(replicate_boxplots, showfliers=False, patch_artist=True)

    for i, patch in enumerate(box2['boxes']):
        patch.set_facecolor('white')
        patch.set_edgecolor(blue_colors[i])
    for i, (whisker1, whisker2) in enumerate(zip(box2['whiskers'][::2], box2['whiskers'][1::2])):
        whisker1.set_color(blue_colors[i])
        whisker2.set_color(blue_colors[i])
    for i, (cap1, cap2) in enumerate(zip(box2['caps'][::2], box2['caps'][1::2])):
        cap1.set_color(blue_colors[i])
        cap2.set_color(blue_colors[i])
    for i, median in enumerate(box2['medians']):
        median.set_color(blue_colors[i])

    ax2.set_xticks(range(1, num_boxes + 1))
    ax2.set_xticklabels([f'Well {i+1}' for i in range(num_boxes)], fontsize=8)
    for i, medians in enumerate(replicate_medians):
        x = [i + 1] * len(medians)
        ax2.scatter(x, medians, facecolors='none', edgecolors='black', linewidths=1.2, zorder=3)

    plt.tight_layout()
    fig2.savefig(f"{save_dir}/replicate_boxplot.png")
    plt.close(fig2)
    
    
    # EXPERIMENT boxplot
    fig3, ax3 = plt.subplots(figsize=(2, 2), dpi=150)
    boxprops = ax3.boxplot(experiment_boxplots, showfliers=False, patch_artist=True)
    colors = plt.cm.tab10.colors
    num_boxes = len(boxprops['boxes'])
    for i in range(num_boxes):
        boxprops['boxes'][i].set_facecolor(colors[i % len(colors)])
        boxprops['medians'][i].set_color('white')
        boxprops['medians'][i].set_linewidth(1)
    ax3.set_xticks(range(1, len(experiment_boxplots) + 1))
    ax3.set_xticklabels([f'Plate {i+1}' for i in range(len(experiment_boxplots))], fontsize=8)
    plt.tight_layout()
    fig3.savefig(f"{save_dir}/experiment_boxplot.png")
    plt.close(fig3)



from setup_plate import setup_plate

# plate0 = "D:/ANALYSIS_PAPER/240313_W1_FPs"
# plate1 = "D:/ANALYSIS_PAPER/240326_W1_FPs"
# plate2 = "D:/ANALYSIS_PAPER/240417_W1_FPs"
# plate3 = "D:/ANALYSIS_PAPER/240924_W1_FPs"
# plate4 = "D:/ANALYSIS_PAPER/241113_W1_FPs"
# plate5 = "D:/ANALYSIS_PAPER/250326_W1_FPs"
plate2bis = "D:/ANALYSIS_PAPER/new_threshold/240417_W1_FPs"
plate3bis = "D:/ANALYSIS_PAPER/new_threshold/240924_W1_FPs"
plate4bis = "D:/ANALYSIS_PAPER/new_threshold/241113_W1_FPs"
prefix = 'SR_001.MIA'

group_poca_files, group_palm_tracer_files = setup_plate(plate4bis, prefix)

protein = "mEos4b"
parameter = "intensity"  # or "intensity", "avg_on", "photon_loc", "total ON", "blinks", "Integrated_Intensity"

# Generate boxplots for well scale
well_boxplots = well_scale(protein, parameter, group_poca_files, group_palm_tracer_files)
replicate_boxplots, replicate_medians = replicate_scale(protein, parameter, group_poca_files, group_palm_tracer_files)
experiment_boxplots = experiment_scale([plate2bis, plate3bis, plate4bis], protein, parameter, prefix)
# figure2(well_boxplots, replicate_boxplots, replicate_medians, experiment_boxplots)
save_individual_boxplots(well_boxplots, replicate_boxplots, replicate_medians, experiment_boxplots)
