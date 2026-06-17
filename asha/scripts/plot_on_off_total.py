import os
import pandas as pd
import matplotlib.pyplot as plt

from asha.src.io_utils import get_poca_files, read_poca_files


def plot_plate_on_off_histograms(plate_path, bins=100):
    """
    Extracts and plots outlier-free distributions for total_on, avg_on, and avg_off.

    Traverses all POCA files in the plate, calculates durations per molecule,
    removes outliers using the Interquartile Range (IQR) method (the core of a boxplot),
    and displays the cleaned data across three histograms.

    Args:
        plate_path (str): Absolute path to the plate directory.
        bins (int, optional): Number of bins for the histograms. Defaults to 100.

    Returns:
        None: Displays the matplotlib figure directly.
    """
    
    all_total_on = []
    all_avg_on = []
    all_avg_off = []
    
    print(f"Reading plate data from: {os.path.basename(plate_path)}")
    poca_files = get_poca_files(plate_path)
    
    if not poca_files:
        print("Warning: No POCA files found in this plate directory.")
        return

    for file in poca_files:
        try:
            df = read_poca_files(file)
            if df.empty:
                continue
                
            avg_on = df['total ON'] / df['# seq ON']
            avg_off = df['total OFF'] / df['# seq OFF']
            
            all_total_on.extend(df['total ON'].dropna().tolist())
            all_avg_on.extend(avg_on.dropna().tolist())
            all_avg_off.extend(avg_off.dropna().tolist())
            
        except Exception as e:
            print(f"Error reading file {file}: {e}")

    if not all_total_on or not all_avg_on or not all_avg_off:
        print("No valid data could be extracted.")
        return

    def filter_outliers(data_list):
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
    
    # total ON time
    axes[0].hist(clean_total_on, bins=bins, color='lightgreen', edgecolor='none', alpha=0.8)
    axes[0].set_title('Distribution of Total ON Time', fontweight='bold', size=11)
    axes[0].set_xlabel('Total ON time (frames)')
    axes[0].set_ylabel('Count')
    
    # avg ON duration
    axes[1].hist(clean_avg_on, bins=bins, color='skyblue', edgecolor='none', alpha=0.8)
    axes[1].set_title('Distribution of ON Duration', fontweight='bold', size=11)
    axes[1].set_xlabel('Avg. ON duration (frames)')
    axes[1].set_ylabel('Count')
    
    # avg OFF duration
    axes[2].hist(clean_avg_off, bins=bins, color='lightcoral', edgecolor='none', alpha=0.8)
    axes[2].set_title('Distribution of OFF Duration', fontweight='bold', size=11)
    axes[2].set_xlabel('Avg. OFF duration (frames)')
    axes[2].set_ylabel('Count')
    
    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(axis='y', linestyle='--', alpha=0.3)
        
    plt.tight_layout()
    plt.show()

plot_plate_on_off_histograms("D:/ANALYSIS_PAPER/new_threshold/240417_W1_FPs/C4")