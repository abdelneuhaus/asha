import matplotlib.pyplot as plt

from src.io_utils import read_poca_files, get_poca_files, get_PALMTracer_files, read_locPALMTracer_file


list_of_pt_files = get_PALMTracer_files('D:/ANALYSIS_PAPER/gamme/MEOS')
list_of_poca_files = get_poca_files('D:/ANALYSIS_PAPER/gamme/MEOS')

gradient_spycatcher1 = ["C9", "C10", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "E3", "E4", "E5", "E6"]
concentration = [1, 0.5, 0.2, 0.15, 0.1, 0.075, 0.05, 0.01]
surface = (256 * 0.16) * (256 * 0.16)

plot_conc_densities = []
plot_densities = []
plot_conc_clusters = []
plot_clusters = []

print("Extracting data...")

for position, conc in zip(gradient_spycatcher1, concentration):
    position_files_pt = [f for f in list_of_pt_files if f"/{position}/" in f]
    position_files_poca = [f for f in list_of_poca_files if f"/{position}/" in f]
    
    density_per_file = []
    cluster_per_file = []
    
    for file in position_files_pt:
        raw_file_pt = read_locPALMTracer_file(file)
        if not raw_file_pt.empty:
            total_locs = len(raw_file_pt)
            unique_planes = raw_file_pt['Plane'].nunique()
            
            if unique_planes > 0:
                mean_density = (total_locs / unique_planes) / surface
                if mean_density > 0:
                    density_per_file.append(mean_density)
            
    for file in position_files_poca:
        raw_file_poca = read_poca_files(file)
        cluster_per_file.append(len(raw_file_poca))

    if density_per_file:
        plot_conc_densities.append(conc)
        plot_densities.append(sum(density_per_file) / len(density_per_file))
        
    if cluster_per_file:
        plot_conc_clusters.append(conc)
        plot_clusters.append(sum(cluster_per_file) / len(cluster_per_file))

print("Generating plot...")

fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.scatter(plot_conc_densities, plot_densities, color='steelblue', label='Densité de locs')
ax1.set_xlabel('mEOS3.2 Concentration (log scale)', fontsize=12)
ax1.set_ylabel('Loc Density (loc/µm²)', color='steelblue', fontsize=12)
ax1.tick_params(axis='y', labelcolor='steelblue')
ax1.set_xscale('log')

ax2 = ax1.twinx()
ax2.scatter(plot_conc_clusters, plot_clusters, color='salmon', label='Nombre de clusters')
ax2.set_ylabel('#Clusters', color='indianred', fontsize=12)
ax2.tick_params(axis='y', labelcolor='indianred')

ax1.grid(True, which="major", linestyle="--", linewidth=0.7, alpha=0.7) 
ax1.tick_params(axis='both', which='minor', length=4, color='gray')

fig.tight_layout()
plt.show()