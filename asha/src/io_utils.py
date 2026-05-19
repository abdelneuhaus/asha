import os
import re

import numpy as np
import pandas as pd

from tkinter import Tk, filedialog
from collections import defaultdict


def get_poca_files(repertory, prefix='SR_001.MIA'):
        return [os.path.join(dirpath, f).replace("\\", "/") for dirpath, _, files in os.walk(repertory) if os.path.basename(dirpath) == prefix for f in files if f == 'locPALMTracer_merged.txt']


def get_PALMTracer_files(repertory, prefix='SR_001.MIA'):
    return [os.path.join(dirpath, f).replace("\\", "/") for dirpath, _, files in os.walk(repertory) if os.path.basename(dirpath) == prefix for f in files if f == 'locPALMTracer.txt']


def get_statMIA_files(repertory):
    return [os.path.join(dirpath, filename).replace("\\", "/") for dirpath, _, filenames in os.walk(repertory) for filename in filenames if filename.endswith('statMIA.txt')]


def read_poca_files(file):
    df = pd.read_csv(file)
    df.loc[df['# seq OFF'] > 5000, '# seq OFF'] = df['blinks']
    return df.iloc[:,:-1]


def read_locPALMTracer_file(file):
    data = pd.read_csv(file, sep='\t', skiprows=2)
    data['id'] = [int(i) for i in data['id']]
    data['Plane'] = [int(i) for i in data['Plane']]
    data['Index'] = [int(i) for i in data['Index']]
    return data


def read_statMIA(file):
    data = pd.read_csv(file, sep='\t', skiprows=2)
    data['Plane'] = data['Plane'].astype(int)
    data['Index'] = data['Index'].astype(int)
    return data



def create_database():
    """
    Create a database mapping target proteins to their corresponding wells and indices.
    Returns:
        dict: A dictionary mapping target proteins to a tuple containing a list of wells and their index.
    """
    database = dict()
    database["mEos4b-L93M"] = ['C3', 'D9', 'E3'], 0
    database["mEos3.2"] = ['C4', 'D8', 'E4'], 1
    database["mEosEM"] = ['C5', 'D7', 'E5'], 2
    database["pcSTAR"] = ['C6', 'D6', 'E6'], 3
    database["Dendra2"] = ['C7', 'D5', 'E7'], 4
    database["Maple3"] = ['C8', 'D4', 'E8'], 5
    database["mEos4b"] = ['C9', 'D3', 'E9'], 6
    return database



def creer_matrice_et_medianes(noms_fichiers, valeurs, positions):
    matrice_puits = np.full((len(positions),), float('nan'))
    valeurs_par_position = defaultdict(list)
    for nom_fichier, valeur in zip(noms_fichiers, valeurs):
        match = re.search(r'([A-H]\d+)', nom_fichier)
        if match:
            position = match.group(1)
            valeurs_par_position[position].append(valeur)
    for i, position in enumerate(positions):
        if position in valeurs_par_position:
            mediane = np.median(valeurs_par_position[position])
            matrice_puits[i] = round(float(mediane), 1)
    return matrice_puits



