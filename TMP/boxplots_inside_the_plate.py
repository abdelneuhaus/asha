import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np
import random
from asha.utils import read_poca_files, read_locPALMTracer_file, create_database


def boxplots_inside_the_plate(group_poca_files, group_palm_tracer_files):
    """
    This function generates boxplots for a specific parameter inside a plate.
    Args:
        target_protein (str): The target protein to analyze.
        group_poca_files (list): List of POCA files grouped by protein.
        group_palm_tracer_files (list): List of PALM Tracer files grouped by protein.
    """
    
    database = create_database()
    protein_names = list(database.keys())

    # Ces listes contiendront les données et noms valides uniquement
    valid_protein_names = []
    intensity_values_all = []
    total_on_values_all = []
    blinks_values_all = []
    integrated_intensity_values_all = []

    for i, protein in enumerate(protein_names):
        poca_files = group_poca_files[i]
        palm_tracer_files = group_palm_tracer_files[i]

        poca_data = [read_poca_files(f) for f in poca_files]
        palm_tracer_data = [read_locPALMTracer_file(f) for f in palm_tracer_files]

        try:
            intensity = np.concatenate([df["intensity"] * (3.6/300) for df in poca_data if not df.empty])
            total_on = np.concatenate([df["total ON"] for df in poca_data if not df.empty])
            blinks = np.concatenate([df["blinks"] for df in poca_data if not df.empty])
            integrated_intensity = np.concatenate([df["Integrated_Intensity"] * (3.6/300) for df in palm_tracer_data if not df.empty])

            # Ne garder que si les données sont non vides
            if all(len(arr) > 0 for arr in [intensity, total_on, blinks, integrated_intensity]):
                intensity_values_all.append(intensity)
                total_on_values_all.append(total_on)
                blinks_values_all.append(blinks)
                integrated_intensity_values_all.append(integrated_intensity)
                valid_protein_names.append(protein)

        except Exception as e:
            print(f"Données invalides pour {protein} : {e}")
            continue

    # Plot
    # Génère une palette de couleurs distinctes
    num_proteins = len(valid_protein_names)
    colors = cm.get_cmap('Paired', num_proteins).colors  # 'tab20' donne des couleurs bien distinctes

    fig, axs = plt.subplots(2, 2, figsize=(16, 12))
    axs = axs.flatten()

    # Liste des valeurs à tracer
    all_values = [
        intensity_values_all,
        total_on_values_all,
        blinks_values_all,
        integrated_intensity_values_all
    ]

    titles = [
        "Photons/Molecule",
        "Total ON (frames)",
        "NBlinks",
        "Photons/Localization"
    ]

    # Tracer les boxplots un par un avec couleurs personnalisées
    for i in range(4):
        bp = axs[i].boxplot(all_values[i], patch_artist=True, labels=valid_protein_names, showfliers=False)
        axs[i].set_title(titles[i])

        # Appliquer une couleur différente à chaque boîte
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)

        axs[i].tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.show()
