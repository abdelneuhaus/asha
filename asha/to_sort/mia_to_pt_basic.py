import pandas as pd
import os
import math

def get_PALMTracer_files(repertory):
    return [os.path.join(dirpath, filename).replace("\\", "/") for dirpath, _, filenames in os.walk(repertory) for filename in filenames if filename.endswith('statMIA.txt')]

def read_statMIA(file):
    data = pd.read_csv(file, sep='\t', skiprows=2)
    data['Plane'] = data['Plane'].astype(int)
    data['Index'] = data['Index'].astype(int)
    return data

list_of_pt_files = get_PALMTracer_files('D:/DPP/')

for mia_file in list_of_pt_files:
    raw_data = read_statMIA(mia_file)
    raw_data['Intensity Gauss'] *= 2 * math.pi * raw_data['SigmaX'] * raw_data['SigmaY']
    filtered_data = raw_data[(raw_data['SigmaX'] < 1.15) & (raw_data['SigmaX'] > 0.85) & 
                             (raw_data['Intensity Gauss'] <= 13500) & (raw_data['Intensity Gauss'] > 0.0)]

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