import matplotlib.pyplot as plt
import numpy as np
import random
from asha.utils import read_poca_files, read_locPALMTracer_file, create_database


def boxplots_inside_a_well(target_protein, group_poca_files, group_palm_tracer_files):
    
    database = create_database()
    target_well = random.choice(database[target_protein][0])
    protein_index = database[target_protein][1]

    target_well_poca_files = [group for group in group_poca_files[protein_index] if target_well in group]
    target_well_palm_tracer_files = [group for group in group_palm_tracer_files[protein_index] if target_well in group]

    poca_data = [read_poca_files(f) for f in target_well_poca_files]
    palm_tracer_data = [read_locPALMTracer_file(f) for f in target_well_palm_tracer_files]

    intensity_values = [df["intensity"]* (3.6/300) for df in poca_data]
    total_on_values = [df["total ON"] for df in poca_data]
    blinks_values = [df["blinks"] for df in poca_data]
    integrated_intensity_values = [df["Integrated_Intensity"]* (3.6/300) for df in palm_tracer_data]

    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    axs = axs.flatten()

    axs[0].boxplot(intensity_values, labels=[f"FOV {i+1}" for i in range(len(intensity_values))], showfliers=False)
    axs[1].boxplot(total_on_values, labels=[f"FOV {i+1}" for i in range(len(total_on_values))], showfliers=False)
    axs[2].boxplot(blinks_values, labels=[f"FOV {i+1}" for i in range(len(blinks_values))], showfliers=False)
    axs[3].boxplot(integrated_intensity_values, labels=[f"FOV {i+1}" for i in range(len(integrated_intensity_values))], showfliers=False)
    axs[0].set_ylabel("Photons/Molecule")
    axs[1].set_ylabel("Total ON (frames)")
    axs[2].set_ylabel("NBlinks")
    axs[3].set_ylabel("Photons/Localization")

    plt.tight_layout()
    plt.show()