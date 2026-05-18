from PAPER.utils import read_locPALMTracer_file
import pandas as pd

def do_analysis_for_one_acquisition(locpalmtracer_file, x_size_image=256, y_size_image=256):
    raw_file = read_locPALMTracer_file(locpalmtracer_file)
    loc_per_frame = raw_file.groupby(['Plane']).count()['id'].to_frame() # num of loc per frame
    if max(raw_file['Plane']) != len(loc_per_frame):
        for j in range(1,max(raw_file['Plane'])):
            if j not in loc_per_frame.index.values:
                df = pd.DataFrame(0, columns=['id'], index=[j])
                loc_per_frame = pd.concat([loc_per_frame, df]).sort_index()
    
    # Calculate cumulative number of frame, density per frame and average density over the acquisition
    cum_loc_per_frame = loc_per_frame.cumsum()
    surface = (x_size_image*0.16)*(y_size_image*0.16)
    density_per_frame = loc_per_frame/surface
    avg_density = density_per_frame['id'].mean()
    return cum_loc_per_frame, density_per_frame, avg_density