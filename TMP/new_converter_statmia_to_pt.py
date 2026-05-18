import pandas as pd
import os
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats 

def get_statMIA_files(repertory):
    return [os.path.join(dirpath, filename).replace("\\", "/") for dirpath, _, filenames in os.walk(repertory) for filename in filenames if filename.endswith('statMIA.txt')]

def read_statMIA(file):
    data = pd.read_csv(file, sep='\t', skiprows=2)
    data['Plane'] = data['Plane'].astype(int)
    data['Index'] = data['Index'].astype(int)
    return data

# 240417_W1_FPs, 240924_W1_FPs, 241113_W1_FPs
list_of_mia_files = get_statMIA_files('D:/ANALYSIS_PAPER/test_gammes/MEOS_BIS')
nlocs = 60
len_window = 00

for mia_file in list_of_mia_files:
    raw_data = read_statMIA(mia_file)
    raw_data['Intensity Gauss'] *= 2 * math.pi * raw_data['SigmaX'] * raw_data['SigmaY']
    # raw_data = raw_data[(raw_data['SigmaX'] > 0.5) & (raw_data['SigmaX'] < 1.5)]
    # raw_data = raw_data[(raw_data['Intensity Gauss'] >= 4200) & (raw_data['Intensity Gauss'] <= 16666)]
    loc_per_frame = raw_data.groupby('Plane').size()
    # plt.plot(loc_per_frame)
    # plt.title(mia_file)
    # plt.show()
    try:
        frame_with_maximum = 0#loc_per_frame.idxmax()
    except ValueError:
        pass
    raw_data = raw_data[(raw_data['Plane'] >= frame_with_maximum)]
    loc_per_frame = raw_data.groupby('Plane').size()
    rolling_avg = loc_per_frame.rolling(window=len_window).mean()
    valid_frames = rolling_avg[rolling_avg < nlocs].index
    
    if not valid_frames.empty:
        start_frame = valid_frames[0]
        raw_data = raw_data[(raw_data['Plane'] >= start_frame)]
        loc_per_frame = loc_per_frame[loc_per_frame.index >= start_frame]
        rolling_avg = rolling_avg[rolling_avg.index >= start_frame]
        # rd = raw_data[(raw_data['Plane'] == start_frame)]
        # print(rd['Intensity Gauss']*(3.6/300))
        # # afficher les points sur une figure 256x256
        # plt.figure(figsize=(8, 8))
        # plt.scatter(rd['Index'], rd['Intensity Gauss']*(3.6/300), color='blue', label='Intensity Gauss')
        # plt.xlabel('Index')
        # plt.ylabel('Intensity Gauss (scaled)')
        # plt.title(f'Intensity Gauss for {os.path.basename(mia_file)} at Frame {start_frame}')
        # plt.legend()
        # plt.show()

        # Plot the density and rolling average
        # plt.figure(figsize=(10, 6))
        # plt.plot(loc_per_frame.index, loc_per_frame, label='Density (Localizations per Frame)', alpha=0.7)
        # plt.plot(rolling_avg.index, rolling_avg, label='Rolling Average', color='red', linewidth=2)
        # plt.axhline(nlocs, color='gray', linestyle='--', label='Threshold')
        # plt.xlabel('Frame')
        # plt.ylabel('Density (Localizations per Frame)')
        # plt.title(f'Density Over Time for {os.path.basename(mia_file)}')
        # plt.legend()
        # plt.show()
    
        filtered_data = raw_data.rename(columns={'Chi2(Gauss)': 'MSE(Gauss)',
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
            
    else:
        print(f"No valid frames found for file: {mia_file}")
print("Done")
