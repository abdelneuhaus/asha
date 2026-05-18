import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_statMIA(file):
    data = pd.read_csv(file, sep='\t', skiprows=2)
    data['Plane'] = data['Plane'].astype(int)
    data['Index'] = data['Index'].astype(int)
    return data

def get_poca_files(repertory):
    a = [os.path.join(dirpath,filename) for dirpath, _, filenames in os.walk(repertory) for filename in filenames if filename.endswith('merged.txt')]
    return [x.replace("\\", "/") for x in a]

def read_poca_files(file):
    df = pd.read_csv(file)
    df.loc[df['# seq OFF'] > 5000, '# seq OFF'] = df['blinks']
    return df.iloc[:,:-1]

def read_locPALMTracer_file(file):
    data=pd.read_csv(file, sep='\t', skiprows=2)
    data['id'] = [int(i) for i in data['id']]
    data['Plane'] = [int(i) for i in data['Plane']]
    data['Index'] = [int(i) for i in data['Index']]
    return data

def get_PALMTracer_files(repertory):
    a = [os.path.join(dirpath,filename) for dirpath, _, filenames in os.walk(repertory) for filename in filenames if filename.endswith('locPALMTracer.txt')]
    return [x.replace("\\", "/") for x in a]


def get_statMIA_files(repertory):
    return [os.path.join(dirpath, filename).replace("\\", "/") for dirpath, _, filenames in os.walk(repertory) for filename in filenames if filename.endswith('statMIA.txt')]


def calculate_median_values(repertory):
    median_values = {}
    clusters = []
    for root, dirs, files in os.walk(repertory):
        for dir in dirs:
            if dir.startswith('P 0'):
                poca_files = get_poca_files(os.path.join(root, dir))
                pt_files = get_PALMTracer_files(os.path.join(root, dir))
                medians = []
                for f in range(len(poca_files)):
                    raw_file_poca = read_poca_files(poca_files[f])
                    raw_pt_file = read_locPALMTracer_file(pt_files[f])
                    # raw_file_poca = raw_file_poca[raw_file_poca['blinks'] < np.quantile(raw_file_poca["blinks"], 1)]
                    # raw_file_poca = raw_file_poca[raw_file_poca['total ON'] < np.quantile(raw_file_poca["total ON"], 1)]
                    raw_file_poca["avg_on"] = np.array(raw_file_poca['total ON'] / raw_file_poca['# seq ON'])
                    raw_file_poca["photon_loc"] = np.array(raw_file_poca['intensity'] / raw_file_poca['total ON'])
                    medians.append(np.median(raw_file_poca["total ON"]))#*(3.6/300))
                    # medians.append(len(raw_file_poca["total ON"]))
                top_level_folder = os.path.basename(((os.path.dirname(root))))
                if top_level_folder not in median_values:
                    median_values[top_level_folder] = []
                median_values[top_level_folder].append(np.median(medians))
    return median_values

def plot_median_values(median_values):
    fig, ax = plt.subplots()
    for key, values in median_values.items():
        ax.scatter([key]*len(values), values, label=key)
    ax.set_xlabel('Parent Folder')
    ax.set_ylabel('Median Values')
    ax.legend()
    plt.show()

def calculate_cumulative_localizations(repertory):
    cumulative_localizations = []
    lenn = 0
    for root, dirs, files in os.walk(repertory):
        for dir in dirs:
            if dir.startswith('P 0'):
                    pt_files = get_PALMTracer_files(os.path.join(root, dir))
                    for f in range(len(pt_files)):
                        raw_pt_file = read_locPALMTracer_file(pt_files[f])
                        lenn += len(raw_pt_file)
                        if pt_files[f] != 'D:/2503_comparison/15sec/2409/C4/P 03/SR.PT/locPALMTracer.txt':
                        # raw_pt_file = raw_pt_file[(raw_pt_file['Plane'] >= 1500) & (raw_pt_file['Plane'] <= 5000)]
                            cumulative_sum = raw_pt_file.groupby('Plane').size()#.cumsum()
                            cumulative_localizations.append(cumulative_sum)
    return cumulative_localizations

def plot_cumulative_localizations(cumulative_localizations):
    fig, ax = plt.subplots()
    for  cum_sum in cumulative_localizations:
        ax.plot(cum_sum)#, label=label)
    ax.set_xlabel('Frame')
    ax.set_ylabel('Localizations')
    ax.legend()
    plt.show()
    
    
    
def show_locs_density_before_filtering(repertory):
    mia_files = get_statMIA_files(repertory)
    plt.figure(figsize=(10, 6))
    for f in mia_files:
        raw_mia_file = read_statMIA(f)
        locs_per_plane = raw_mia_file.groupby('Plane').size()
        plt.plot(locs_per_plane.index, locs_per_plane.values)  # Nom du fichier en légende
    plt.xlabel("Frame (time)")
    plt.ylabel("NLocs")
    plt.title('Before filtering')
    plt.grid(True)
    plt.show()

    
def show_locs_density_after_filtering(repertory):
    pt_files = get_PALMTracer_files(repertory)
    plt.figure(figsize=(10, 6))
    for f in pt_files:
        raw_pt_file = read_locPALMTracer_file(f)
        locs_per_plane = raw_pt_file.groupby('Plane').size()
        plt.plot(locs_per_plane.index, locs_per_plane.values)  # Nom du fichier en légende
    plt.xlabel("Frame (time)")
    plt.ylabel("NLocs")
    plt.title('After filtering')
    plt.grid(True)
    plt.show()


def show_clusters(repertory):
    poca_files = get_poca_files(repertory)
    plt.figure(figsize=(10, 6))
    for f in poca_files:
        raw_file_poca = read_poca_files(f)
        locs_per_plane = raw_file_poca['total ON']
        plt.hist(locs_per_plane, bins=100)
        plt.xlabel("Frame (time)")
        plt.ylabel("NLocs")
        plt.title('')
        plt.grid(True)
        plt.show()
    
    
# Main execution
# repertory = 'D:/241114_W1_FPs/C4'
# show_locs_density_before_filtering(repertory)
# show_locs_density_after_filtering(repertory)

repertory = 'D:/ANALYSIS_PAPER/241113_W1_FPs'
# show_locs_density_before_filtering(repertory)
# show_locs_density_after_filtering(repertory)
median_values = calculate_median_values(repertory)
# plot_median_values(median_values)

# cumulative_localizations = calculate_cumulative_localizations(repertory)
# plot_cumulative_localizations(cumulative_localizations)
