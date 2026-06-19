import os

import numpy as np
import matplotlib.pyplot as plt
import concurrent.futures

from asha.src.io_utils import read_poca_files, get_poca_files, get_PALMTracer_files, read_locPALMTracer_file


def plot_linear_fits(conc_densities, densities, conc_clusters, clusters, fit_mode="full"):
    """
    Generates a dual-axis plot with linear regression fits for density and cluster data.

    Args:
        conc_densities (list[float] or np.ndarray): concentration values corresponding to density measurements.
        densities (list[float] or np.ndarray): Mean molecular density values.
        conc_clusters (list[float] or np.ndarray): concentration values corresponding to cluster measurements.
        clusters (list[float] or np.ndarray): Mean cluster count values.
        fit_mode (str, optional): The fitting range mode.

    Returns:
        None: Displays the matplotlib figure directly.
    """
    fig, ax1 = plt.subplots(figsize=(7, 5))

    x_den, y_den = np.array(conc_densities), np.array(densities)
    sort_den = np.argsort(x_den)
    x_den, y_den = x_den[sort_den], y_den[sort_den]

    x_clu, y_clu = np.array(conc_clusters), np.array(clusters)
    sort_clu = np.argsort(x_clu)
    x_clu, y_clu = x_clu[sort_clu], y_clu[sort_clu]

    if fit_mode == "0-2":
        mask_den = x_den <= 2.0
        mask_clu = x_clu <= 2.0
        fit_label_suffix = "(0 à 2)"
    else:
        mask_den = np.ones_like(x_den, dtype=bool)
        mask_clu = np.ones_like(x_clu, dtype=bool)
        fit_label_suffix = "(Complet)"

    a_den, b_den = np.polyfit(x_den[mask_den], y_den[mask_den], 1)
    fit_den_line = a_den * x_den[mask_den] + b_den
    r2_den = np.corrcoef(x_den[mask_den], y_den[mask_den])[0, 1]**2

    ax1.scatter(x_den, y_den, color='steelblue', label='Densité data')
    ax1.plot(x_den[mask_den], fit_den_line, color='steelblue', linestyle='--', 
             label=f'Fit Densité {fit_label_suffix} : y={a_den:.2f}x + {b_den:.2f} ($R^2$={r2_den:.3f})')
    
    ax1.set_xlabel('SpyCatcher Concentration', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Loc Density (loc/µm²)', color='steelblue', fontsize=12, fontweight='bold')
    ax1.tick_params(axis='y', labelcolor='steelblue')

    ax2 = ax1.twinx()
    a_clu, b_clu = np.polyfit(x_clu[mask_clu], y_clu[mask_clu], 1)
    fit_clu_line = a_clu * x_clu[mask_clu] + b_clu
    r2_clu = np.corrcoef(x_clu[mask_clu], y_clu[mask_clu])[0, 1]**2

    ax2.scatter(x_clu, y_clu, color='salmon', label='Clusters data')
    ax2.plot(x_clu[mask_clu], fit_clu_line, color='salmon', linestyle='--', 
             label=f'Fit Clusters {fit_label_suffix} : y={a_clu:.1f}x + {b_clu:.1f} ($R^2$={r2_clu:.3f})')
    
    ax2.set_ylabel('#Clusters', color='indianred', fontsize=12, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor='indianred')

    ax1.grid(True, which="major", linestyle="--", linewidth=0.7, alpha=0.7)
    
    lines_1, labels_1 = ax1.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='upper left', framealpha=0.9)

    fig.tight_layout()
    plt.show()



def process_single_well(args):
    """
    Processes a single well to extract mean molecular density and cluster count.

    Args:
        args (tuple): A tuple containing the following elements in order:
            - position (str): The well identifier (e.g., "C9").
            - conc (float): The corresponding concentration value.
            - list_pt (list[str]): List of all available PALMTracer file paths.
            - list_poca (list[str]): List of all available PoCA file paths.
            - surface (float): The conversion factor for surface area (µm²).

    Returns:
        tuple: A tuple containing (conc, avg_density, avg_clusters). 
    """
    position, conc, list_pt, list_poca, surface = args
    
    pos_files_pt = [f for f in list_pt if f"/{position}/" in f.replace('\\', '/')]
    pos_files_poca = [f for f in list_poca if f"/{position}/" in f.replace('\\', '/')]
    
    density_per_file = []
    cluster_per_file = []
    
    for file in pos_files_pt:
        raw_file_pt = read_locPALMTracer_file(file)
        if not raw_file_pt.empty:
            total_locs = len(raw_file_pt)
            unique_planes = raw_file_pt['Plane'].nunique()
            if unique_planes > 0:
                mean_density = (total_locs / unique_planes) / surface
                if mean_density > 0:
                    density_per_file.append(mean_density)
            
    for file in pos_files_poca:
        raw_file_poca = read_poca_files(file)
        cluster_per_file.append(len(raw_file_poca))

    avg_density = sum(density_per_file) / len(density_per_file) if density_per_file else None
    avg_clusters = sum(cluster_per_file) / len(cluster_per_file) if cluster_per_file else None

    return conc, avg_density, avg_clusters


if __name__ == '__main__':
    print("Extracting data (Multiprocessing enabled)...")

    list_of_pt_files = get_PALMTracer_files('D:/ANALYSIS_PAPER/gammes/SPYCATCHER')
    list_of_poca_files = get_poca_files('D:/ANALYSIS_PAPER/gammes/SPYCATCHER')

    spycatcher_wells = ["C9", "C10", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "E3", "E4", "E5", "E6"]
    concentration = [10, 7.5, 5, 2, 1, 0.75, 0.5, 0.2, 0.1, 0.075, 0.05, 0.02, 0.01, 0]
    surface = (256 * 0.16) * (256 * 0.16)

    plot_conc_densities = []
    plot_densities = []
    plot_conc_clusters = []
    plot_clusters = []

    tasks = [(pos, conc, list_of_pt_files, list_of_poca_files, surface) for pos, conc in zip(spycatcher_wells, concentration)]

    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = list(executor.map(process_single_well, tasks))

    for conc, avg_density, avg_clusters in results:
        if avg_density is not None:
            plot_conc_densities.append(conc)
            plot_densities.append(avg_density)
        if avg_clusters is not None:
            plot_conc_clusters.append(conc)
            plot_clusters.append(avg_clusters)


    print("Generating linear fit plot...")
    plot_linear_fits(plot_conc_densities, plot_densities, plot_conc_clusters, plot_clusters, fit_mode="0-2")