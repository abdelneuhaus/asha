import os
import pandas as pd
import numpy as np

from src.io_utils import get_poca_files, read_poca_files


def process_and_export_plates(plate_paths, stat_choice, output_file="./results/hcs_excel.xlsx"):
    """Processes multiple HCS-SMLM plates and exports the chosen statistics to Excel.

    Args:
        plate_paths (list of str): List of absolute paths to the plate directories.
        stat_choice (str): The statistical metric to compute ('Median' or 'Mean').
        output_file (str, optional): The destination path for the Excel file.

    Returns:
        None: Saves the compiled data directly to the specified Excel file.
    """
    fixed_subfolders = {
        "mEos4b_": ["C9", "D3", "E9"],
        "mEos4b-L93M": ["C3", "D9", "E3"],
        "mEos3.2": ["C4", "D8", "E4"],
        "mEosEM": ["C5", "D7", "E5"],
        "pcSTAR": ["C6", "D6", "E6"],
        "Dendra2": ["C7", "D5", "E7"],
        "mMaple3": ["C8", "D4", "E8"],
    }
    
    parameters = ["photon_loc", "total ON", "intensity", "blinks", "avg_on", "avg_off"]
    global_data = {param: {prot: [] for prot in fixed_subfolders} for param in parameters}
    stat_func = np.median if stat_choice == 'Median' else np.mean

    print("Extracting data...")
    
    for plate_path in plate_paths:
        if not plate_path.strip():
            continue
            
        print(f"Plate processing : {os.path.basename(plate_path)}")
        
        for protein, wells in fixed_subfolders.items():
            plate_prot_data = {param: [] for param in parameters}
            
            for well in wells:
                well_path = os.path.join(plate_path, well)
                if not os.path.exists(well_path):
                    continue
                    
                poca_files = get_poca_files(well_path) 
                
                for file in poca_files:
                    try:
                        df = read_poca_files(file) 
                        photon_loc = (df['intensity'] / df['total ON']) * (3.6/300)
                        intensity = df['intensity'] * (3.6/300)
                        avg_on = df['total ON'] / df['# seq ON']
                        avg_off = df['total OFF'] / df['# seq OFF']

                        plate_prot_data["photon_loc"].append(stat_func(photon_loc.dropna()))
                        plate_prot_data["total ON"].append(stat_func(df["total ON"].dropna()))
                        plate_prot_data["intensity"].append(stat_func(intensity.dropna()))
                        plate_prot_data["blinks"].append(stat_func(df["blinks"].dropna()))
                        plate_prot_data["avg_on"].append(stat_func(avg_on.dropna()))
                        plate_prot_data["avg_off"].append(stat_func(avg_off.dropna()))

                    except Exception as e:
                        print(f"Error on file {file}: {e}")
            
            for param in parameters:
                global_data[param][protein].extend(plate_prot_data[param])
                global_data[param][protein].append(np.nan)

    # excel export
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    print(f"Saving in {output_file}...")
    
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for param in parameters:
            max_length = max(len(global_data[param][prot]) for prot in fixed_subfolders)
            aligned_data = {
                prot: global_data[param][prot] + [np.nan] * (max_length - len(global_data[param][prot]))
                for prot in fixed_subfolders
            }
            df_export = pd.DataFrame(aligned_data)
            
            if not df_export.empty and df_export.iloc[-1].isna().all():
                df_export = df_export.iloc[:-1]
                
            df_export.to_excel(writer, index=False, sheet_name=param)
            
    print("Export worked properly!") 