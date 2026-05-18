import matplotlib.pyplot as plt
import os

from PAPER.utils import read_locPALMTracer_file


def save_loc_as_pdf(locpalmtracer_file, idx, density_per_frame, cum_loc_per_frame, avg_density, exp):
    raw_file = read_locPALMTracer_file(locpalmtracer_file)
    fig = plt.figure(layout='constrained')
    fig.suptitle('Results of '+idx)
    axis1 = fig.add_subplot(211)
    axis2 = fig.add_subplot(212)
    axis1.plot(range(max(raw_file['Plane'])), density_per_frame)
    axis2.plot(range(max(raw_file['Plane'])), cum_loc_per_frame)
    axis1.set_title('Density through time')
    axis2.set_title('Cumulative number of localizations through time')
    axis1.set_ylabel("Density (mol/µm²)")
    axis2.set_ylabel("Cumul. number of localisations")
    axis2.set_xlabel('Time (frames)\n\nAverage density per frame: '+str(round(avg_density, 5))+' mol/µm²')
    
    # save in results/name_of_exp
    results_dir = os.path.join('results/'+exp+'/'+idx+'/')
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    plt.savefig(results_dir)
    plt.close('all')