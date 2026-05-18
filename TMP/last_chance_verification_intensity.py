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

def plot_total_ON(meos32, base_path):
    colors = ['red', 'green', 'blue', "yellow", "orange"]  # Couleurs pour chaque dossier de meos32
    plt.figure(figsize=(10, 6))
    for idx, folder in enumerate(meos32):
        folder_path = os.path.join(base_path, folder)
        poca_files = get_poca_files(folder_path)
        pt_files = get_PALMTracer_files(folder_path)       
        total_on_values = []
        for f in range(len(poca_files)):
            poca_data = read_poca_files(poca_files[f])
            raw_pt_file = read_locPALMTracer_file(pt_files[f])
            poca_data["avg_on"] = np.array(poca_data['total ON'] / poca_data['# seq ON'])
            poca_data["photon_loc"] = (np.array(poca_data['intensity'] / poca_data['total ON']))*(3.6/300)
            poca_data['intensity'] = poca_data['intensity']*(3.6/300)
            raw_pt_file["Integrated_Intensity"] = raw_pt_file["Integrated_Intensity"]*(3.6/300)
            plt.hist(raw_pt_file["Integrated_Intensity"], bins=100, alpha=0.6, color=colors[idx], label=f"{folder}")
            plt.show()
    #         total_on_values.extend(raw_pt_file["Integrated_Intensity"].values)
    #     # Tracer les valeurs de total ON
    #     # plt.hist(total_on_values, bins=100, alpha=0.6, color=colors[idx], label=f"{folder}")
    #     # plt.xlim(0, 1000)  # Limiter l'axe des x à 1000
    #     plt.boxplot(total_on_values, positions=[idx], widths=0.6, patch_artist=True, boxprops=dict(facecolor=colors[idx], color=colors[idx]), medianprops=dict(color='black'), showfliers=False)
    # plt.xlabel('Total ON')
    # plt.ylabel('Frequency')
    # plt.title('Distribution of Total ON Values')
    # plt.show()
    
    
def plot_median_points(meos32, base_path):    
    colors = ['red', 'green', 'blue', 'yellow', 'orange']  # Couleurs pour chaque groupe
    for idx, folder in enumerate(meos32):
        folder_path = os.path.join(base_path, folder)
        poca_files = get_poca_files(folder_path)
        
        median_values = []  # Liste pour stocker les valeurs médianes de chaque fichier poca_file
        for poca_file in poca_files:
            poca_data = read_poca_files(poca_file)
            poca_data["avg_on"] = np.array(poca_data['total ON'] / poca_data['# seq ON'])
            poca_data["photon_loc"] = (np.array(poca_data['intensity'] / poca_data['total ON'])) * (3.6 / 300)
            poca_data['intensity'] = poca_data['intensity'] * (3.6 / 300)
            
            # Calculer la médiane pour ce fichier poca_file
            median_value = np.median(poca_data['total ON'].values)
            median_values.append(median_value)
        
        # Ajouter les points pour ce groupe
        plt.scatter([idx] * len(median_values), median_values, color=colors[idx], label=folder, alpha=0.7)
    
    # Ajouter des labels et personnaliser le graphique
    meos32 = ["10 mW", "20 mW", "30 mW", "50 mW", "70 mW"]
    plt.xticks(range(len(meos32)), meos32)
    # plt.legend(loc='upper right')
    plt.show()

# Main execution
base_path = 'D:/250408_HCS_High power/'
meos32 = ["10mW/C4", "405-10microW_561-20mW/W1/D8", "405-10microW_561-30mW/W1/E4", "405-10microW_561-50mW/W1/F8", "405-10microW_561-70mW/W1/G4"]
meos4b = ["10mW/C9", "405-10microW_561-20mW/W1/D3", "405-10microW_561-30mW/W1/E9", "405-10microW_561-50mW/W1/F3", "405-10microW_561-70mW/W1/G9"]
# plot_total_ON(meos4b, base_path)
plot_median_points(meos4b, base_path)