import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_plate_on_off_histograms(plate_path):
    """
    Extracts and plots outlier-free distributions for molecular ON and OFF times.

    Args:
        plate_path (str): The absolute path to the plate or well directory containing the 'frame.csv' files.

    Returns:
        None: The function generates and displays the three-panel histogram directly using matplotlib.
    """
    all_total_on = []
    all_avg_on = []
    all_avg_off = []
    
    print(f"Recherche des fichiers frame.csv dans : {plate_path}")
    
    frame_files = []
    for root, _, files in os.walk(plate_path):
        if 'frame.csv' in files:
            frame_files.append(os.path.join(root, 'frame.csv').replace('\\', '/'))

    if not frame_files:
        print("error: no frame.csv files")
        return

    print(f"{len(frame_files)} files found.")

    for frame_file in frame_files:
        try:
            with open(frame_file, 'r') as f:
                for line in f:
                    parts = line.strip().split(',')
                    
                    if len(parts) < 2:
                        continue
                        
                    try:
                        frames = sorted([int(p) for p in parts[1:] if p.strip()])
                    except ValueError:
                        continue 
                        
                    if not frames:
                        continue
                        
                    on_lengths = []
                    off_lengths = []
                    current_on_length = 1
                    
                    for i in range(1, len(frames)):
                        diff = frames[i] - frames[i-1]
                        
                        if diff == 1:
                            current_on_length += 1
                        elif diff > 1:
                            on_lengths.append(current_on_length)
                            off_lengths.append(diff - 1) 
                            current_on_length = 1 
                            
                    on_lengths.append(current_on_length)

                    all_total_on.append(sum(on_lengths))
                    all_avg_on.append(np.mean(on_lengths))
                    
                    if off_lengths:
                        all_avg_off.append(np.mean(off_lengths))
                        
        except Exception as e:
            print(f"error reading file {frame_file}: {e}")

    if not all_total_on or not all_avg_on:
        print("no valid data.")
        return

    def filter_outliers(data_list):
        if not data_list:
            return []
        series = pd.Series(data_list)
        q1 = series.quantile(0.25)
        q3 = series.quantile(0.75)
        iqr = q3 - q1
        lower_bound = q1 - 1.5 * iqr
        upper_bound = q3 + 1.5 * iqr
        return series[(series >= lower_bound) & (series <= upper_bound)].tolist()

    clean_total_on = filter_outliers(all_total_on)
    clean_avg_on = filter_outliers(all_avg_on)
    clean_avg_off = filter_outliers(all_avg_off)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    
    if clean_total_on:
        min_val = int(np.floor(min(clean_total_on)))
        max_val = int(np.ceil(max(clean_total_on)))
        bins_tot = np.arange(min_val, max_val + 2) - 0.5
        axes[0].hist(clean_total_on, bins=bins_tot, color='lightgreen', edgecolor='white', linewidth=0.5, alpha=0.8)
        
    axes[0].set_title('Distribution of Total ON Time', fontweight='bold', size=11)
    axes[0].set_xlabel('Total ON time (frames)')
    axes[0].set_ylabel('Count')
    
    if clean_avg_on:
        min_val = int(np.floor(min(clean_avg_on)))
        max_val = int(np.ceil(max(clean_avg_on)))
        bins_avg_on = np.arange(min_val, max_val + 2) - 0.5
        axes[1].hist(clean_avg_on, bins=bins_avg_on, color='skyblue', edgecolor='white', linewidth=0.5, alpha=0.8)
        
    axes[1].set_title('Distribution of ON Duration', fontweight='bold', size=11)
    axes[1].set_xlabel('Avg. ON duration (frames)')
    axes[1].set_ylabel('Count')
    
    if clean_avg_off:
        min_val = int(np.floor(min(clean_avg_off)))
        max_val = int(np.ceil(max(clean_avg_off)))
        bins_avg_off = np.arange(min_val, max_val + 2) - 0.5
        axes[2].hist(clean_avg_off, bins=bins_avg_off, color='lightcoral', edgecolor='white', linewidth=0.5, alpha=0.8)
        
    axes[2].set_title('Distribution of OFF Duration', fontweight='bold', size=11)
    axes[2].set_xlabel('Avg. OFF duration (frames)')
    axes[2].set_ylabel('Count')
    
    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.3)
        ax.xaxis.get_major_locator().set_params(integer=True)
        
    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    plot_plate_on_off_histograms("D:/ANALYSIS_PAPER/new_threshold/240417_W1_FPs/C4")