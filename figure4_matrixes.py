import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from utils import get_PALMTracer_files, get_poca_files, read_locPALMTracer_file, read_poca_files



def process_files(poca_files, pt_files, param):
    """
    Process a list of poca and pt files, computing photon_loc values.
    """
    data_list = []
    for poca_file, pt_file in zip(poca_files, pt_files):
        raw_file_poca = read_poca_files(poca_file)
        raw_pt_file = read_locPALMTracer_file(pt_file)
        raw_file_poca = raw_file_poca[raw_file_poca['total ON'] <= np.quantile(raw_file_poca["total ON"], 1)]
        raw_file_poca = raw_file_poca[raw_file_poca['intensity'] <= np.quantile(raw_file_poca["intensity"], 1)]
        raw_file_poca["photon_loc"] = raw_file_poca['intensity'] / raw_file_poca['total ON']*(3.6/300)
        raw_file_poca["intensity"] = raw_file_poca['intensity']*(3.6/300)
        data_list.append(raw_file_poca[param].values)
    return np.concatenate(data_list)


def normalize_data(data_dict, reference_key):
    """
    Normalize the data in data_dict using the values from reference_key as 100%.
    """
    reference_data = data_dict[reference_key]
    normalized_data = {key: (values / np.mean(reference_data)) * 100 for key, values in data_dict.items()}
    return normalized_data



def create_comparison_matrix(normalized_data, parameter):
    """
    Create a symmetric double-entry table comparing proteins based on the chosen parameter.
    If protein1 has 10% more than protein2 in the upper triangle, the lower triangle will show -10%.
    """
    proteins = list(normalized_data.keys())
    comparison_matrix = pd.DataFrame(index=proteins, columns=proteins, dtype=float)

    for i, prot1 in enumerate(proteins):
        for j, prot2 in enumerate(proteins):
            if i < j:  # Fill upper triangle
                mean1 = np.mean(normalized_data[prot1])
                mean2 = np.mean(normalized_data[prot2])
                if mean2 != 0:
                    ratio = (mean1 / mean2 - 1)  # Percentage difference
                    comparison_matrix.loc[prot1, prot2] = ratio
                    comparison_matrix.loc[prot2, prot1] = -ratio  # Symmetric value
                else:
                    comparison_matrix.loc[prot1, prot2] = np.nan
                    comparison_matrix.loc[prot2, prot1] = np.nan
            elif i == j:
                comparison_matrix.loc[prot1, prot2] = 0  # Identity values
    return comparison_matrix




def plot_boxplots(normalized_data, parameters):
    """
    Plot boxplots for all parameters for each protein.
    """
    for param in parameters:
        plt.figure(figsize=(12, 6))
        plt.title(f"Boxplot of Proteins for Parameter: {param}", fontsize=16)
        data_to_plot = [normalized_data[protein] for protein in normalized_data.keys()]
        plt.boxplot(data_to_plot, labels=normalized_data.keys(), patch_artist=True, showfliers=False)
        plt.xlabel("Proteins", fontsize=14)
        plt.ylabel(f"{param} (Normalized)", fontsize=14)
        plt.xticks(rotation=45, fontsize=12)
        plt.tight_layout()
        plt.show()



# Define the directories for each protein
proteins = {
    "mEos3.2": [],
    "mEos4b": [],
    "mEos4b-L93M": [],
    "pcSTAR": [],
    "mEosEM": [],
    "Dendra2": [],
    "Maple3": []
}

plates = [
    # "D:/ANALYSIS_PAPER/240313_W1_FPs",
    # "D:/ANALYSIS_PAPER/240326_W1_FPs",
    "D:/ANALYSIS_PAPER/new_threshold/240417_W1_FPs",
    "D:/ANALYSIS_PAPER/new_threshold/240924_W1_FPs",
    "D:/ANALYSIS_PAPER/new_threshold/241113_W1_FPs"
    # "D:/ANALYSIS_PAPER/250326_W1_FPs"
]

def create_database():
    return {
        "mEos4b-L93M": (['C3', 'D9', 'E3'], 0),
        "mEos3.2": (['C4', 'D8', 'E4'], 1),
        "mEosEM": (['C5', 'D7', 'E5'], 2),
        "pcSTAR": (['C6', 'D6', 'E6'], 3),
        "Dendra2": (['C7', 'D5', 'E7'], 4),
        "Maple3": (['C8', 'D4', 'E8'], 5),
        "mEos4b": (['C9', 'D3', 'E9'], 6)
    }

# Mise à jour des chemins pour chaque protéine
db = create_database()

for protein_name in proteins:
    if protein_name in db:
        wells, _ = db[protein_name]
        proteins[protein_name] = [
            f"{plate}/{well}" for plate in plates for well in wells
        ]


# Process all proteins
# Process all proteins with multiple directories
poca_data = {}
pt_data = {}
parameters = ['total ON','intensity', 'blinks', 'photon_loc']
# parameters = ['intensity']#,'blinks']
for a in parameters:
    for protein, paths in proteins.items():
        protein_poca_data = []
        protein_pt_data = []

        for path in paths:
            poca_files = get_poca_files(path)
            pt_files = get_PALMTracer_files(path)
            
            if not poca_files:
                print(f"Warning: No POCA files found for {protein} in {path}")
            if not pt_files:
                print(f"Warning: No PT files found for {protein} in {path}")

            try:
                processed_data = process_files(poca_files, pt_files, a)
                if processed_data.size > 0:  # Check if any data was processed
                    protein_poca_data.append(processed_data)
            except Exception as e:
                print(f"Error processing files for {protein} in {path}: {e}")

        # Combine data from all directories for this protein
        if protein_poca_data:
            poca_data[protein] = np.concatenate(protein_poca_data)
        else:
            print(f"Warning: No valid POCA data for {protein}")

    # Normalize data using mEos3.2 as the reference
    if "mEos3.2" not in poca_data:
        raise ValueError("mEos3.2 data is missing and cannot be used as a reference for normalization.")

    normalized_poca_data = normalize_data(poca_data, reference_key="mEos3.2")
    # plot_boxplots(normalized_poca_data, [a])
    # Create comparison matrix for photon_loc
    comparison_matrix = create_comparison_matrix(normalized_poca_data, parameter=a)

    # Display the matrix
    print(comparison_matrix)
    print(a)

    # Optionally save the matrix to a CSV file
    # comparison_matrix.to_csv("protein_comparison_matrix.csv")

    # Visualize the comparison matrix as a heatmap
    plt.figure(figsize=(10, 8))
    plt.title("Protein Comparison Matrix " + a)
    plt.imshow(comparison_matrix.astype(float), cmap="RdBu")#, interpolation='none')
    plt.colorbar(label="Ratio")
    plt.xticks(range(len(comparison_matrix.columns)), comparison_matrix.columns, rotation=45)
    plt.yticks(range(len(comparison_matrix.index)), comparison_matrix.index)
    plt.tight_layout()
    plt.show()
