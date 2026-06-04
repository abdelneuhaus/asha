import os
import numpy as np
import pandas as pd

from src.io_utils import get_poca_files, read_poca_files, fusion, creer_matrice_et_medianes



def save_heatmap_as_excel(pathway, output_file="./results/heatmap_matrices.xlsx", stats=np.median):
    """
    Extracts all photophysical parameters from a plate and exports a 8x12 matrix for each parameter into an Excel file.

    Args:
        pathway (str): Absolute path to the plate directory containing the PoCA files.
        output_file (str): Destination path for the generated Excel file.
        stats (callable, optional): The statistical function to apply (default is np.median).

    Returns:
        None: The function saves an Excel file, with each parameter in a separate sheet.
    """
    exp_name = os.path.basename(os.path.normpath(pathway))
    list_of_poca = get_poca_files(pathway) 
    
    if not list_of_poca:
        print(f"Error : No PoCA files found in {pathway}")
        return

    parameters = ["photon_loc", "total ON", "intensity", "blinks", "avg_on", "avg_off"]
    heatmap_data_dict = {param: [] for param in parameters}
    
    print(f"Extraction for plate : {exp_name}...")
    
    for f in list_of_poca:
        df = read_poca_files(f)
        intensity_corr = df['intensity'] * (3.6 / 300)
        
        for param in parameters:
            if param == 'avg_on':
                values = np.divide(df['total ON'], df['# seq ON'], out=np.zeros_like(df['total ON'], dtype=float), where=df['# seq ON']!=0)
            elif param == 'avg_off':
                values = np.divide(df['total OFF'], df['# seq OFF'], out=np.zeros_like(df['total OFF'], dtype=float), where=df['# seq OFF']!=0)
            elif param == 'photon_loc':
                values = np.divide(intensity_corr, df['total ON'], out=np.zeros_like(intensity_corr, dtype=float), where=df['total ON']!=0)
            elif param == 'intensity':
                values = intensity_corr
            else:
                values = df[param] 
                
            clean_values = pd.Series(values).replace([np.inf, -np.inf], np.nan).dropna()
            
            if len(clean_values) > 0:
                heatmap_data_dict[param].append(float(stats(clean_values)))
            else:
                heatmap_data_dict[param].append(np.nan)
                
    idx = list(map(str, range(1, 13)))
    cols = list(map(str, [chr(i) for i in range(ord('A'), ord('H')+1)]))
    all_wells = fusion(cols, idx)
    
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    print(f"Save in {output_file}...")
    
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for param in parameters:
            matrice_resultante = creer_matrice_et_medianes(list_of_poca, heatmap_data_dict[param], all_wells)
            df_heat = pd.DataFrame(np.array(matrice_resultante).reshape(8, 12), index=cols, columns=idx)
            df_heat.to_excel(writer, sheet_name=param)
            
    print(f"Excel export done!")