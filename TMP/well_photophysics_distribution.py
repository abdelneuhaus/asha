import os
import matplotlib.pyplot as plt
import numpy as np
import statistics as statistics

from preprocessing import photon_calculation


def well_photophysics_distribution(grouped_data, well, expname):

    global_trend_total_on = []
    global_trend_intensity = []
    global_trend_photloc = []
    global_trend_on_lifetime = []
    
    # Create subplots in a 2x2 grid
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 7))
    fig.patch.set_facecolor('white')
    for ax in [ax1, ax2, ax3, ax4]:
        ax.set_facecolor('white')

    # Plot total ON time in the first subplot (ax1)
    ax1.plot(list(grouped_data['avg_on']))
    ax1.set_xlabel('Frame Intervals')
    ax1.set_ylabel('Mean ON time duration')
    
    # Plot photon per loc in the second subplot (ax2)
    ax2.plot(photon_calculation(list(grouped_data['photon_loc'])))
    ax2.set_xlabel('Frame Intervals')
    ax2.set_ylabel('Mean Photon per Loc')
    
    # Plot total ON time in the third subplot (ax3)
    ax3.plot(list(grouped_data['total ON']))
    ax3.set_xlabel('Frame Intervals')
    ax3.set_ylabel('Mean total ON time')
    
    # Plot intensity in the fourth subplot (ax4)
    ax4.plot(photon_calculation(list(grouped_data['intensity'])))
    ax4.set_xlabel('Frame Intervals')
    ax4.set_ylabel('Mean photon budget')
    
    # Append global trends
    global_trend_total_on.extend(grouped_data['avg_on'])
    global_trend_photloc.extend(photon_calculation(grouped_data['photon_loc']))
    global_trend_on_lifetime.extend(grouped_data['total ON'])
    global_trend_intensity.extend(photon_calculation(grouped_data['intensity']))
    
    # Calculate the global trend lines using linear regression
    frames = np.array(list(range(0, 5001, 500)))
    coeff_total_on = np.polyfit(frames, global_trend_total_on[:len(frames)], 1)
    coeff_photon_loc = np.polyfit(frames, global_trend_photloc[:len(frames)], 1)
    coeff_on_lifetime = np.polyfit(frames, global_trend_on_lifetime[:len(frames)], 1)
    coeff_intensity = np.polyfit(frames, global_trend_intensity[:len(frames)], 1)
    
    trendline_total_on = np.polyval(coeff_total_on, frames)
    trendline_photon_loc = np.polyval(coeff_photon_loc, frames)
    trendline_on_lifetime = np.polyval(coeff_on_lifetime, frames)
    trendline_intensity = np.polyval(coeff_intensity, frames)
    
    # Plot the trend lines
    ax1.plot(trendline_total_on, '--', color='gray', label='Trend Line')
    ax2.plot(trendline_photon_loc, '--', color='gray', label='Trend Line')
    ax3.plot(trendline_on_lifetime, '--', color='gray', label='Trend Line')
    ax4.plot(trendline_intensity, '--', color='gray', label='Trend Line')
    
    results_dir = os.path.join('results/'+expname+'/'+well+'/')
    sample_file = 'well_fitted_distribution_photophysic_parameters.pdf'
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    plt.savefig(results_dir+sample_file)
    plt.close(fig)