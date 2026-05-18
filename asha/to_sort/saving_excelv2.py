import os
import pandas as pd
import numpy as np

# ---------- I/O ----------

def get_poca_files(repertory):
    files = [
        os.path.join(dirpath, filename)
        for dirpath, _, filenames in os.walk(repertory)
        for filename in filenames
        if filename.endswith("merged.txt")
    ]
    return [f.replace("\\", "/") for f in files]


def read_poca_files(file):
    df = pd.read_csv(file)
    df.loc[df["# seq OFF"] > 5000, "# seq OFF"] = df["blinks"]
    return df.iloc[:, :-1]


# ---------- PROCESS ----------

def process_files(file_list):
    data = {
        "photon_loc": [],
        "total ON": [],
        "intensity": [],
        "blinks": []
    }

    for file in file_list:
        raw = read_poca_files(file)

        raw = raw[raw["blinks"] < np.quantile(raw["blinks"], 1)]
        raw = raw[raw["total ON"] < np.quantile(raw["total ON"], 1)]
        raw["photon_loc"] = raw["intensity"] / raw["total ON"]

        data["photon_loc"].append(np.median(raw["photon_loc"].values * (3.6 / 300)))
        data["total ON"].append(np.median(raw["total ON"].values))
        data["intensity"].append(np.median(raw["intensity"].values * (3.6 / 300)))
        data["blinks"].append(np.median(raw["blinks"].values))

    return data


# ---------- SAVE ----------

def save_results(data_dict, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with pd.ExcelWriter(output_file, engine="openpyxl") as writer:
        for param in ["photon_loc", "total ON", "intensity", "blinks"]:
            max_len = max(len(data_dict[name][param]) for name in data_dict)

            aligned = {
                name: data_dict[name][param] + [np.nan] * (max_len - len(data_dict[name][param]))
                for name in data_dict
            }

            pd.DataFrame(aligned).to_excel(
                writer, index=False, sheet_name=param
            )


# ---------- CONFIG ----------

names = [
    "mEos3.2",
    "mEos4b_",
    "mEos4b-L93M",
    "pcSTAR",
    "mEosEM",
    "mMaple3",
    "Dendra2"
]

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


# ---------- MAIN ----------

global_data = {
    name: {"photon_loc": [], "total ON": [], "intensity": [], "blinks": []}
    for name in names
}

for base in base_dirs:

    # stock temporaire du bloc courant
    block_data = {
        name: {"photon_loc": [], "total ON": [], "intensity": [], "blinks": []}
        for name in names
    }

    # traitement du dossier
    for name in names:
        file_list = []
        for sub in fixed_subfolders[name]:
            file_list.extend(get_poca_files(f"{base}/{sub}"))

        processed = process_files(file_list)

        for param in processed:
            block_data[name][param].extend(processed[param])

    # hauteur max du bloc (toutes protéines confondues)
    block_height = max(
        len(block_data[name][param])
        for name in names
        for param in block_data[name]
    )

    # padding + concaténation
    for name in names:
        for param in block_data[name]:
            current_len = len(block_data[name][param])
            block_data[name][param].extend([np.nan] * (block_height - current_len))
            global_data[name][param].extend(block_data[name][param])

    # ligne vide de séparation ENTRE dossiers
    for name in names:
        for param in global_data[name]:
            global_data[name][param].append(np.nan)


save_results(global_data, "./results/combined/results_combined.xlsx")
