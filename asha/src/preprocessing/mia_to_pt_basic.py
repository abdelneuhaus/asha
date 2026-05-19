import pandas as pd
import math

from src.io_utils import get_statMIA_files, read_statMIA


def mia_to_pt_basic(plate_path:str):
    """
    Convert all statMIA files from a directory to PALMTracer format in the same directory.

    Args:
        plate_path (str): Path to the plate directory.

    Returns:
        None. Save the converted PALMTracer file in the same place as the statMIA file.
    """
    list_of_pt_files = get_statMIA_files(plate_path)

    for mia_file in list_of_pt_files:
        raw_data = read_statMIA(mia_file)
        raw_data['Intensity Gauss'] *= 2 * math.pi * raw_data['SigmaX'] * raw_data['SigmaY']
        filtered_data = filtered_data.rename(columns={'Chi2(Gauss)': 'MSE(Gauss)',
                                                    'Intensity Gauss': 'Integrated_Intensity',
                                                    'CentroidX(pix)': 'CentroidX(px)', 
                                                    'CentroidY(pix)': 'CentroidY(px)',
                                                    'SigmaX': 'SigmaX(px)',
                                                    'SigmaY': 'SigmaY(px)',
                                                    'AngleRad': 'Angle(rad)',
                                                    'MSE(Z)': 'MSE_Z(um)',
                                                    'CentroidZ(µm)': 'CentroidZ(um)'})
        
        filtered_data = filtered_data.drop(columns=['Mean', 'Surface', 'Intensity', 'Perimeter', 'Morpho', 'AngleDeg', 'maxX'])
        filtered_data['id'] = range(1, len(filtered_data) + 1)
        filtered_data['Channel'] = 0.0
        filtered_data = filtered_data[['id', 'Plane', 'Index', 'Channel', 'Integrated_Intensity', 
                                    'CentroidX(px)', 'CentroidY(px)', 'SigmaX(px)', 'SigmaY(px)', 
                                    'Angle(rad)', 'MSE(Gauss)', 'CentroidZ(um)', 'MSE_Z(um)']]
        
        meta_data = pd.read_csv(mia_file, nrows=1, sep='\t')
        meta_data = meta_data.rename(columns={'Planes count': 'nb_Planes', 'Points count': 'nb_Points'})
        meta_data = meta_data[['Width', 'Height', 'nb_Planes', 'nb_Points']]
        meta_data['Pixel_Size(um)'] = 0.160
        meta_data['Frame_Duration(s)'] = 0.05
        meta_data['Gaussian_Fit'] = str(None)
        meta_data['Spectral'] = False
        meta_data['nb_Points'] = len(filtered_data)
        
        output_file = mia_file.replace('statMIA', 'locPALMTracer')
        
        with open(output_file, 'w', newline='') as file:
            file.write('\t'.join(meta_data.columns) + '\n')
            file.write(meta_data.to_csv(sep='\t', index=False, header=False))
            file.write('\t'.join(filtered_data.columns) + '\n')
            file.write(filtered_data.to_csv(sep='\t', index=False, header=False))