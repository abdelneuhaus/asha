from PAPER.utils import read_poca_files, get_poca_files, get_PALMTracer_files, read_locPALMTracer_file
import matplotlib.pyplot as plt
import numpy as np 

# Charger les fichiers
list_of_pt_files = get_PALMTracer_files('D:/250205_HCS_gamme/W1/')
list_of_poca_files = get_poca_files('D:/250205_HCS_gamme/W1/')
gradient_spycatcher0 = ["B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "C3", "C4", "C5", "C6", "C7", "C8"]
gradient_spycatcher1 = ["C9", "C10", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "E3", "E4", "E5", "E6"]
gradient_meos0 = ["F3", "F4", "F5", "F6", "F7", "F8", "F9", "F10"]
gradient_meos1 = ["G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10"]
# meos1 = ["B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10"]
# meos2 = ["C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"]
# meos3 = ["D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10"]

surface = (256 * 0.16) * (256 * 0.16)


average_density_per_position = {}
number_of_clusters_per_position = {}
density_per_file = []
cluster_per_file = []

# Calculer les densités et le nombre de clusters par position
for position in gradient_meos1:
    position_files_pt = [f for f in list_of_pt_files if f"/{position}/" in f]
    position_files_poca = [f for f in list_of_poca_files if f"/{position}/" in f]
    
#     for file in position_files_pt:
#         raw_file_pt = read_locPALMTracer_file(file)
#         loc_per_frame = raw_file_pt.groupby(['Plane']).size()
#         density_per_frame = loc_per_frame / surface
#         mean_density = density_per_frame.median()
#         print(file,"density:", mean_density)
#         if (mean_density) > 0:
#             density_per_file.append(mean_density)
# print(density_per_file)
    for file in position_files_poca:
        raw_file_poca = read_poca_files(file)
        clusters_num = len(raw_file_poca)
        cluster_per_file.append(clusters_num)
        print(file,"clusters:", clusters_num)
print(cluster_per_file)

        
        
#     if density_per_file:
#         average_density_per_position[position] = sum(density_per_file) / len(density_per_file)
#     if cluster_per_file:
#         number_of_clusters_per_position[position] = sum(cluster_per_file) / len(cluster_per_file)

# # # Exemple de données
# concentration = [10, 7.5, 5, 2, 1, 0.75, 0.5, 0.2, 0.1, 0.075, 0.05, 0.02, 0.01, 0]
# concentration0 = [1, 0.5, 0.2, 0.15, 0.1, 0.075, 0.05, 0.01]
# densities = list(average_density_per_position.values())
# print(densities)
# clusters = list(number_of_clusters_per_position.values())
# print(clusters)