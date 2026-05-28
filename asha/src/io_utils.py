import os
import re

import numpy as np
import pandas as pd

from collections import defaultdict


def get_poca_files(repertory, prefix='SR_001.MIA'):
    """
    Finds all 'locPALMTracer_merged.txt' files within specific directories.

    Args:
        repertory (str): The root directory to start the recursive search.
        prefix (str, optional): The name of the sub-directory containing the files. 
            Defaults to 'SR_001.MIA'.

    Returns:
        list of str: A list of absolute file paths to all found PoCA merge files.
    """
    return [os.path.join(dirpath, f).replace("\\", "/") for dirpath, _, files in os.walk(repertory) if os.path.basename(dirpath) == prefix for f in files if f == 'locPALMTracer_merged.txt']


def get_PALMTracer_files(repertory, prefix='SR_001.MIA'):
    """
    Finds all 'locPALMTracer.txt' files within specific directories.

    Args:
        repertory (str): The root directory to start the recursive search.
        prefix (str, optional): The name of the sub-directory containing the files. 
            Defaults to 'SR_001.MIA'.

    Returns:
        list of str: A list of absolute file paths to all found PALMTracer files.
    """
    return [os.path.join(dirpath, f).replace("\\", "/") for dirpath, _, files in os.walk(repertory) if os.path.basename(dirpath) == prefix for f in files if f == 'locPALMTracer.txt']


def get_statMIA_files(repertory):
    """
    Finds all files ending with 'statMIA.txt' recursively.

    Args:
        repertory (str): The root directory to start the recursive search.

    Returns:
        list of str: A list of absolute file paths to all found statMIA files.
    """
    return [os.path.join(dirpath, filename).replace("\\", "/") for dirpath, _, filenames in os.walk(repertory) for filename in filenames if filename.endswith('statMIA.txt')]


def read_poca_files(file):
    """
    Reads a PoCA CSV file and applies data cleaning to the OFF sequence column.

    This function reads the CSV, corrects outliers in '# seq OFF' by capping 
    them with 'blinks' values, and drops the last column (metadata/filler).

    Args:
        file (str): Path to the PoCA CSV file.

    Returns:
        pd.DataFrame: The cleaned and processed dataframe.
    """
    df = pd.read_csv(file)
    df.loc[df['# seq OFF'] > 5000, '# seq OFF'] = df['blinks']
    return df.iloc[:,:-1]


def read_locPALMTracer_file(file):
    """
    Reads a tab-separated PALMTracer file with custom data type conversion.

    Skips the first two header rows and enforces integer types for 'id', 
    'Plane', and 'Index' columns.

    Args:
        file (str): Path to the tab-separated PALMTracer text file.

    Returns:
        pd.DataFrame: The formatted dataframe with correctly typed columns.
    """
    data = pd.read_csv(file, sep='\t', skiprows=2)
    data['id'] = [int(i) for i in data['id']]
    data['Plane'] = [int(i) for i in data['Plane']]
    data['Index'] = [int(i) for i in data['Index']]
    return data


def read_statMIA(file):
    """
    Reads a statMIA tab-separated file and prepares it for density analysis.

    Args:
        file (str): Path to the statMIA tab-separated text file.

    Returns:
        pd.DataFrame: The processed dataframe with correctly typed columns.
    """
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
    """
    Computes the median value for each well and maps them to a standard plate layout.

    Args:
        noms_fichiers (list of str): List of file paths used to extract well IDs.
        valeurs (list of float): List of values corresponding to each file.
        positions (list of str): List of all expected well positions (e.g., ['A1', 'A2', ...]).

    Returns:
        np.ndarray: A 1D numpy array of length `len(positions)` containing 
            the median values (rounded to 1 decimal place), with NaNs for empty wells.
    """
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


def fusion(cols, idx):
    """
    Generates all combinations of rows and columns to represent a 96-well plate.

    Args:
        cols (list of str): A list of row identifiers (e.g., ['A', 'B', ...]).
        idx (list of str): A list of column identifiers (e.g., ['1', '2', ...]).

    Returns:
        list of str: A sorted list of all well identifiers in the order 
            defined by the input lists (e.g., ['A1', 'A2', ..., 'H12']).
    """
    return [f"{c}{i}" for c in cols for i in idx]

