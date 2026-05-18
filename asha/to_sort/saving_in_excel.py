import os
import pandas as pd
import numpy as np

def get_poca_files(repertory):
    a = [os.path.join(dirpath, filename) for dirpath, _, filenames in os.walk(repertory) for filename in filenames if filename.endswith('merged.txt')]
    return [x.replace("\\", "/") for x in a]

def read_poca_files(file):
    df = pd.read_csv(file)
    df.loc[df['# seq OFF'] > 5000, '# seq OFF'] = df['blinks']
    return df.iloc[:, :-1]

def process_files(file_list, names):
    data_dict = {name: {"photon_loc": [], "total ON": [], "intensity": [], "blinks": []} for name in names}
    for file in file_list:
        raw_file_poca = read_poca_files(file)
        raw_file_poca = raw_file_poca[raw_file_poca['blinks'] < np.quantile(raw_file_poca["blinks"], 1)]
        raw_file_poca = raw_file_poca[raw_file_poca['total ON'] < np.quantile(raw_file_poca["total ON"], 1)]
        raw_file_poca["photon_loc"] = raw_file_poca['intensity'] / raw_file_poca['total ON']
        for name in names:
            if name in file:
                print(name)
                data_dict[name]["photon_loc"].append(np.median(raw_file_poca["photon_loc"].values*(3.6/300)))
                data_dict[name]["total ON"].append(np.median(raw_file_poca["total ON"].values))
                data_dict[name]["intensity"].append(np.median(raw_file_poca["intensity"].values*(3.6/300)))
                data_dict[name]["blinks"].append(np.median(raw_file_poca["blinks"].values))
    return data_dict

def save_results(data_dict, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for param in ["photon_loc", "total ON", "intensity", "blinks"]:
            max_length = max(len(data_dict[name][param]) for name in data_dict)
            aligned_data = {
                name: data_dict[name][param] + [np.nan] * (max_length - len(data_dict[name][param]))
                for name in data_dict
            }
            df = pd.DataFrame(aligned_data)
            df.to_excel(writer, index=False, sheet_name=param)

names = ["mEos3.2", "mEos4b_", "mEos4b-L93M", "pcSTAR", "mEosEM", "mMaple3", "Dendra2"]

base_dirs = [
    "D:/ANALYSIS_PAPER/new_threshold/240417_W1_FPs",
    "D:/ANALYSIS_PAPER/new_threshold/240924_W1_FPs",
    "D:/ANALYSIS_PAPER/new_threshold/241113_W1_FPs",
]

fixed_subfolders = {
    "mEos4b_": ["C9", "D3", "E9"],
    "mEos4b-L93M": ["C3", "D9", "E3"],
    "mEos3.2": ["C4", "D8", "E4"],
    "mEosEM": ["C5", "D7", "E5"],
    "pcSTAR": ["C6", "D6", "E6"],
    "Dendra2": ["C7", "D5", "E7"],
    "mMaple3": ["C8", "D4", "E8"],
}

proteins = {
    protein: [
        f"{base}/{sub}"
        for base in base_dirs
        for sub in subs
    ]
    for protein, subs in fixed_subfolders.items()
}

# Traiter tous les fichiers et fusionner les résultats
global_data = {name: {"photon_loc": [], "total ON": [], "intensity": [], "blinks": []} for name in proteins}

for name in proteins:
    file_list = []
    for folder_path in proteins[name]:
        file_list.extend(get_poca_files(folder_path))
    processed_data = process_files(file_list, names=[name])
    for param in ["photon_loc", "total ON", "intensity", "blinks"]:
        global_data[name][param].extend(processed_data[name][param])
        global_data[name][param].append(np.nan)  # séparation entre expériences

# Sauvegarder le fichier final
save_results(global_data, "./results/combined/results_combined.xlsx")
