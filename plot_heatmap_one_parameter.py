import os
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numpy as np

from utils import read_poca_files, fusion, creer_matrice_et_medianes, get_poca_files, get_PALMTracer_files, read_locPALMTracer_file


def plot_heatmap_one_parameter(plate_pathway, parameter, prefix):
    
    idx = list(map(str, range(1, 13)))
    cols = list(map(str, [chr(i) for i in range(ord('A'), ord('H')+1)]))
    all_wells = fusion(cols, idx)
    
    list_of_poca_files = get_poca_files(plate_pathway, prefix)
    list_of_palm_tracer_files = get_PALMTracer_files(plate_pathway, prefix)
    heatmap_data = []
    
    if parameter == "Integrated_Intensity":
        for file in list_of_palm_tracer_files:
            df = read_locPALMTracer_file(file)
            df['Integrated_Intensity'] *= 3.6/300
            if parameter in df.columns:
                heatmap_data.append(int(np.median(((df.loc[:, parameter].values.tolist())))))
        matrice_resultante = creer_matrice_et_medianes(list_of_palm_tracer_files, heatmap_data, all_wells)
        df = pd.DataFrame(np.array(matrice_resultante).reshape(8,12), index=cols, columns=idx)
        sns.heatmap(df, annot=True, fmt='g', cmap="coolwarm", linewidths=1, linecolor='black')
        plt.yticks(rotation=0)
        plt.gcf().set_size_inches((12, 5))
        plt.title("Heatmap")
        plt.show()
        
    else:
        for file in list_of_poca_files:
            df = read_poca_files(file)
            df['intensity'] *= 3.6/300
            df["avg_on"] = np.array(df['total ON'] / df['# seq ON'])
            df = df[df['intensity'] < np.quantile(df["intensity"], 1)]
            df = df[df['total ON'] < np.quantile(df["total ON"], 1)]
            if parameter in df.columns:
                heatmap_data.append(int(np.median(((df.loc[:, parameter].values.tolist())))))
                
        matrice_resultante = creer_matrice_et_medianes(list_of_poca_files, heatmap_data, all_wells)
        df = pd.DataFrame(np.array(matrice_resultante).reshape(8,12), index=cols, columns=idx)
        sns.heatmap(df, annot=True, fmt='g', cmap="coolwarm", linewidths=1, linecolor='black')
        plt.yticks(rotation=0)
        plt.gcf().set_size_inches((12, 5))
        plt.title("Heatmap")    
        plt.show()
