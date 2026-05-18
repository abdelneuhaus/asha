import pandas as pd
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats 

def get_PALMTracer_files(repertory):
    return [os.path.join(dirpath, filename).replace("\\", "/") for dirpath, _, filenames in os.walk(repertory) for filename in filenames if filename.endswith('locPALMTracer.txt')]

def read_locPALMTracer_file(file):
    data=pd.read_csv(file, sep='\t', skiprows=2)
    data['id'] = [int(i) for i in data['id']]
    data['Plane'] = [int(i) for i in data['Plane']]
    data['Index'] = [int(i) for i in data['Index']]
    return data

list_of_pt_files = get_PALMTracer_files('D:/DPP/230531/meos4b_f')

for mia_file in list_of_pt_files:
    raw_data = read_locPALMTracer_file(mia_file)
    lower_quantile = 0.6
    upper_quantile = 1.1
    # raw_data = raw_data[(raw_data['SigmaX(px)'] > 0.75) & (raw_data['SigmaX(px)'] < 1.25)]
    raw_data = raw_data[(raw_data['MSE(Gauss)'] > 0.6)]
    # plt.hist(raw_data['SigmaX(px)'], bins=100, color='steelblue', edgecolor='black')
    # plt.show()
    filtered_data = raw_data[(raw_data['SigmaX(px)'] > lower_quantile) & (raw_data['SigmaX(px)'] < upper_quantile)]
    # filtered_data = filtered_data[(filtered_data['Plane'] > 1499)]
    # filtered_data = filtered_data[(filtered_data['Integrated_Intensity'] >= 3375)]
    # filtered_data = raw_data[(raw_data['Integrated_Intensity'] <= 20000)]
    
    filtered_data = raw_data.copy()
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
    output_file = mia_file.replace('locPALMTracer', 'locPALMTracer')
    
    with open(output_file, 'w', newline='') as file:
        file.write('\t'.join(meta_data.columns) + '\n')
        file.write(meta_data.to_csv(sep='\t', index=False, header=False))
        file.write('\t'.join(filtered_data.columns) + '\n')
        file.write(filtered_data.to_csv(sep='\t', index=False, header=False))
print("Done")