import math

import pandas as pd

from src.io_utils import get_statMIA_files, read_statMIA


def mia_to_pt_sliding_windows(plate_path:str, nlocs:int=60, len_window:int=300, 
                              sigma_min:float=0.0, sigma_max:float=2.0,
                              intensity_min:int=0, intensity_max:int=10e6):
    """
    Converts statMIA files to PALMTracer format using a sliding window density filter.

    Args:
        plate_path (str): The root directory containing the statMIA files.
        nlocs (int, optional): Density threshold for localizations per frame. 
            Defaults to 60.
        len_window (int, optional): Window size for rolling average calculation. 
            Defaults to 300.
        sigma_min (float, optional): Minimum threshold for SigmaX. Defaults to 0.0.
        sigma_max (float, optional): Maximum threshold for SigmaX. Defaults to 2.0.
        intensity_min (int, optional): Minimum integrated intensity threshold. Defaults to 0.
        intensity_max (int, optional): Maximum integrated intensity threshold. 
            Defaults to 10e6.

    Returns:
        None: Saves processed data to disk in the same directory as the input files.
    """
    
    list_of_mia_files = get_statMIA_files(plate_path)

    for mia_file in list_of_mia_files:
        raw_data = read_statMIA(mia_file)
        raw_data['Intensity Gauss'] *= 2 * math.pi * raw_data['SigmaX'] * raw_data['SigmaY']
        raw_data = raw_data[(raw_data['SigmaX'] > sigma_min) & (raw_data['SigmaX'] < sigma_max)]
        raw_data = raw_data[(raw_data['Intensity Gauss'] >= intensity_min) & (raw_data['Intensity Gauss'] <= intensity_max)]
        loc_per_frame = raw_data.groupby('Plane').size()

        try:
            frame_with_maximum = 0  # loc_per_frame.idxmax()
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
            
            ## Plot the density and rolling average
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
