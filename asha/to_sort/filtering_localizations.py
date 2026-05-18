import os
import math
import pandas as pd


def get_PALMTracer_files(repertory):
    return [os.path.join(dirpath, filename).replace("\\", "/") for dirpath, _, filenames in os.walk(repertory) for filename in filenames if filename.endswith('locPALMTracer.txt')]

def read_locPALMTracer_file(file):
    data = pd.read_csv(file, sep='\t', skiprows=2)
    data['id'] = [int(i) for i in data['id']]
    data['Plane'] = [int(i) for i in data['Plane']]
    data['Index'] = [int(i) for i in data['Index']]
    return data

def get_statMIA_files(repertory):
    return [os.path.join(dirpath, filename).replace("\\", "/") for dirpath, _, filenames in os.walk(repertory) for filename in filenames if filename.endswith('statMIA.txt')]

def read_statMIA(file):
    data = pd.read_csv(file, sep='\t', skiprows=2)
    data['Plane'] = data['Plane'].astype(int)
    data['Index'] = data['Index'].astype(int)
    return data



def true_filter_with_density(repertory):
    list_of_pt_files = get_statMIA_files(repertory)
    nlocs = 50 # good density (between 30 and 65 usually)
    len_window = 300 # = 15 secondes @ 50 ms/frame

    for mia_file in list_of_pt_files:
        raw_data = read_statMIA(mia_file)
        raw_data['Intensity Gauss'] *= 2 * math.pi * raw_data['SigmaX'] * raw_data['SigmaY']
        raw_data = raw_data[(raw_data['Intensity Gauss'] <= 25000)]
        raw_data = raw_data[(raw_data['SigmaX'] > 0) & (raw_data['SigmaX'] < 2)]

        # Count the number of localizations per frame
        loc_per_frame = raw_data.groupby('Plane').size()
        frame_with_maximum = loc_per_frame.idxmax()
        raw_data = raw_data[(raw_data['Plane'] >= frame_with_maximum)]
        loc_per_frame = raw_data.groupby('Plane').size()

        # Calculate the rolling average over X frames
        rolling_avg = loc_per_frame.rolling(window=len_window).mean()
        # Find the first valid frame where the rolling average is less than density loc
        valid_frames = rolling_avg[rolling_avg < nlocs].index
        # Start plotting from the first valid frame
        if not valid_frames.empty:
            start_frame = valid_frames[0]
            raw_data = raw_data[(raw_data['Plane'] >= start_frame)]
            loc_per_frame = loc_per_frame[loc_per_frame.index >= start_frame]
            rolling_avg = rolling_avg[rolling_avg.index >= start_frame]

            filtered_data = raw_data.rename(columns={'Chi2(Gauss)': 'MSE(Gauss)',
                                                    'Intensity Gauss': 'Integrated_Intensity',
                                                    'CentroidX(pix)': 'CentroidX(px)', 
                                                    'CentroidY(pix)': 'CentroidY(px)',
                                                    'SigmaX': 'SigmaX(px)',
                                                    'SigmaY': 'SigmaY(px)',
                                                    'AngleRad': 'Angle(rad)',
                                                    'MSE(Z)': 'MSE_Z(um)',
                                                    'CentroidZ(µm)': 'CentroidZ(um)'})
        try:
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
        except KeyError as e:
            print(f"KeyError: {e} in file: {mia_file}")
            continue
    else:
        print(f"No valid frames found for file: {mia_file}")



repertory = 'D:/241114_W1_FPs'
true_filter_with_density(repertory)
print("FINISHED")