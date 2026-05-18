import matplotlib.pyplot as plt
import numpy as np 
import random
import os

from asha.utils import read_poca_files, read_locPALMTracer_file

SURFACE = (256*0.16) * (256*0.16)

def save_plot(conc, dens, clus, gradient_name):
    fig, ax1 = plt.subplots(figsize=(5,4))
    ax2 = ax1.twinx()

    color_densities = {
        "gradient0": "steelblue",
        "gradient1": "lightblue"
    }
    color_clusters = {
        "gradient0": "indianred",
        "gradient1": "lightcoral"
    }
    markers = {
        "gradient0": "o",
        "gradient1": "s"
    }

    ax1.scatter(conc, dens, color=color_densities[gradient_name],
                label=f'Densité - {gradient_name}', marker=markers[gradient_name])
    ax2.scatter(conc, clus, color=color_clusters[gradient_name],
                label=f'Clusters - {gradient_name}', marker=markers[gradient_name])

    # ax1.set_xlabel('Concentration mEos3.2 (log scale)', fontsize=12)
    # ax1.set_ylabel('Densité de locs (loc/µm²)', color='steelblue', fontsize=12)
    ax1.set_xscale('log')
    ax1.tick_params(axis='y', labelcolor='steelblue')
    ax1.grid(True, which="major", linestyle="--", linewidth=0.7, alpha=0.7)
    ax1.tick_params(axis='both', which='minor', length=4, color='gray')

    # ax2.set_ylabel('Nombre de clusters', color='indianred', fontsize=12)
    ax2.tick_params(axis='y', labelcolor='indianred')

    fig.tight_layout()

    output_dir = "gradient_plots/meos"
    os.makedirs(output_dir, exist_ok=True)
    filepath = os.path.join(output_dir, f"{gradient_name}.png")
    plt.savefig(filepath, dpi=300)
    plt.show()


def plot_gradient_meos(list_of_pt_files, list_of_poca_files, display_mode="random"):
    gradients = {
        "gradient0": {
            "positions": ["F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10"],
            "concentration": [1, 0.5, 0.2, 0.15, 0.1, 0.075, 0.05, 0.01]
        },
        "gradient1": {
            "positions": ["G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10"],
            "concentration": [1, 0.5, 0.2, 0.15, 0.1, 0.075, 0.05, 0.01]
        }
    }

    if display_mode == "random":
        display_mode = random.choice(["gradient0", "gradient1"])

    selected_gradients = gradients.keys() if display_mode == "both" else [display_mode]

    for gradient_name in selected_gradients:
        data = gradients[gradient_name]
        average_density_per_position = {}
        number_of_clusters_per_position = {}

        for position in data["positions"]:
            position_files_pt = [f for f in list_of_pt_files if f"/{position}/" in f]
            position_files_poca = [f for f in list_of_poca_files if f"/{position}/" in f]
            density_per_file = []
            cluster_per_file = []

            for file in position_files_pt:
                raw_file_pt = read_locPALMTracer_file(file)
                loc_per_frame = raw_file_pt.groupby(['Plane']).size()
                density_per_frame = loc_per_frame / SURFACE
                mean_density = density_per_frame.mean()
                if mean_density > 0:
                    density_per_file.append(mean_density)

            for file in position_files_poca:
                raw_file_poca = read_poca_files(file)
                raw_file_poca = raw_file_poca[raw_file_poca['total ON'] <= np.quantile(raw_file_poca["total ON"], 1)]
                raw_file_poca = raw_file_poca[raw_file_poca['intensity'] <= np.quantile(raw_file_poca["intensity"], 1)]
                raw_file_poca["intensity"] = raw_file_poca['intensity'] * (3.6 / 300)
                clusters_num = len(raw_file_poca)
                cluster_per_file.append(clusters_num)

            if density_per_file:
                average_density_per_position[position] = sum(density_per_file) / len(density_per_file)
            if cluster_per_file:
                number_of_clusters_per_position[position] = sum(cluster_per_file) / len(cluster_per_file)

        densities = list(average_density_per_position.values())
        clusters = list(number_of_clusters_per_position.values())
        concentrations = data["concentration"]
        save_plot(concentrations, densities, clusters, gradient_name)
