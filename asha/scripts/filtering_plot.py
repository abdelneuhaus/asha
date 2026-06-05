import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from asha.src.io_utils import read_statMIA


def format_ticks(x, pos):
    if x >= 1000:
        return f'{x*1e-3:g}k'
    return f'{x:g}'


def plot_filtering_steps_minimal(mia_file: str, nlocs: int=60, len_window: int=300, 
                                 sigma_min: float=0.65, sigma_max: float=1.3,
                                 intensity_min: int=30, intensity_max: int=300):
    """
    Generates a diagnostic plot to visualize SMLM data filtering steps.

    Args:
        mia_file (str): Absolute path to the 'statMIA.txt' file.
        nlocs (int, optional): The threshold for localizations per frame to identify the signal decay onset. Defaults to 60.
        len_window (int, optional): The size of the rolling window for the temporal average. Defaults to 300.
        sigma_min (float, optional): Minimum threshold for SigmaX values. 
        sigma_max (float, optional): Maximum threshold for SigmaX values. 
        intensity_min (int, optional): Minimum intensity threshold in photons. 
        intensity_max (int, optional): Maximum intensity threshold in photons. 

    Returns:
        None: Displays the matplotlib figure directly.
    """
    
    raw_data = read_statMIA(mia_file)
    raw_data['Integrated_Intensity'] = raw_data['Intensity Gauss'] * 2 * math.pi * raw_data['SigmaX'] * raw_data['SigmaY']
    raw_data['Integrated_Intensity'] *= (3.6/300)
    
    fig, axes = plt.subplots(3, 1, figsize=(5, 5))
    
    axes[0].hist(raw_data['SigmaX'], bins=100, color='skyblue', edgecolor='none')
    axes[0].axvline(sigma_min, color='red', linestyle='dashed', linewidth=1.5)
    axes[0].axvline(sigma_max, color='red', linestyle='dashed', linewidth=1.5)
    axes[0].set_ylabel('Count')
    axes[0].set_xlabel('Sigma (pixel)')
    
    upper_lim = 1500
    axes[1].hist(raw_data['Integrated_Intensity'], bins=100, range=(0, upper_lim), color='lightgreen', edgecolor='none')
    axes[1].axvline(intensity_min, color='red', linestyle='dashed', linewidth=1.5)
    if intensity_max <= upper_lim:
        axes[1].axvline(intensity_max, color='red', linestyle='dashed', linewidth=1.5)
    axes[1].set_ylabel('Count')
    axes[1].set_xlabel('Photons per localisation')
    
    filtered_data = raw_data[(raw_data['SigmaX'] > sigma_min) & (raw_data['SigmaX'] < sigma_max)]
    filtered_data = filtered_data[(filtered_data['Integrated_Intensity'] >= intensity_min) & (filtered_data['Integrated_Intensity'] <= intensity_max)]
    
    loc_per_frame = filtered_data.groupby('Plane').size()
    rolling_avg = loc_per_frame.rolling(window=len_window).mean()
    
    peak_frame = rolling_avg.idxmax()
    valid_frames = rolling_avg[(rolling_avg.index > peak_frame) & (rolling_avg < nlocs)].index
    
    axes[2].plot(loc_per_frame.index, loc_per_frame, alpha=0.3, color='gray', linewidth=1)
    axes[2].plot(rolling_avg.index, rolling_avg, color='orange', linewidth=1.5)
    axes[2].axhline(nlocs, color='red', linestyle='--', linewidth=1.5)
    
    if not valid_frames.empty:
        start_frame = valid_frames[0]
        axes[2].axvline(start_frame, color='purple', linestyle='-', linewidth=1.5)
        axes[2].fill_between(loc_per_frame.index, 0, loc_per_frame.max(), 
                             where=(loc_per_frame.index >= start_frame), 
                             color='purple', alpha=0.1)
        
    axes[2].set_ylabel('Locs/frame')
    axes[2].set_xlabel('Frame')
    
    for ax in axes:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.tick_params(axis='both', which='both', length=0)
        
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(format_ticks))

    fig.align_ylabels(axes)

    plt.tight_layout(pad=1.0)
    plt.show()


plot_filtering_steps_minimal("/Users/aneuhaus/Desktop/asha/DATA/W1_bis/C9/P 02/SR_001.MIA/statMIA.txt", nlocs=50, len_window=300)