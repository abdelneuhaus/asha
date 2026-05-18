import csv 
import os
import statistics
import numpy as np
import re
import pandas as pd

from collections import defaultdict
from get_length_on_off import get_length_off, get_length_on
from localization_precision import localization_precision


def lire_csv(nom_fichier):
    lignes = []
    with open(nom_fichier, 'r') as f:
        lecteur = csv.reader(f)
        for ligne in lecteur:
            # Convertir chaque élément de la ligne en entier (integer)
            ligne = [float(element) for element in ligne]
            if len(ligne) < 51:
                lignes.append(ligne[1:])
    return lignes


def pre_process_single_intensity(file):
    tmp = list()
    read_file = lire_csv(file)
    for line in read_file:
        tmp.append(line)
    return [j for i in tmp for j in i]


def pre_process_on_frame_csv(file):
    tmp = list()
    read_file = lire_csv(file)
    for line in read_file:
        tmp.append(get_length_on(line))
    return [j for i in tmp for j in i]


def pre_process_off_frame_csv(file):
    tmp = list()
    read_file = lire_csv(file)
    for line in read_file:
        tmp.append(get_length_off(line)) 
    return [j for i in tmp for j in i]


def pre_process_sigma(file):
    tmp = list()
    read_file = lire_csv(file)
    for line in read_file:
        tmp.append(line)
    return [j for i in tmp for j in i]


def get_and_save_well_and_FOV(i, PT):
    if 'SR' in i and PT in i:
        title_plot = os.sep.join(i.replace(PT, '').split('/')[-2:])
    else:
        idx = os.path.basename(os.path.normpath(i.replace(PT, '')))
        title_plot = os.path.join(idx)
    return title_plot


def pad_list(lst):
    while len(lst) < 8:
        lst.append(float('nan'))
    return lst

def fusion(liste1, liste2):
    fusions = []
    for mot1 in liste1:
        for mot2 in liste2:
            fusions.append(mot1 + mot2)
    return fusions

def fusion_position(liste1, liste2):
    resultat = []
    for i in range(len(liste2)):
        resultat.append(liste1[i] + ': ' + liste2[i])
    return resultat

def photon_calculation(liste, gain=3.6, emgain=300, qe=0.95):
    exp_liste = []
    otp = gain/emgain
    for valeur in liste:
        exp_liste.append(valeur*otp/qe)
    return exp_liste


def loc_prec_calculation(sigma, photon_loc):
    otp = []
    median = statistics.median(sigma)
    for i in range(len(sigma)):
        otp.append(localization_precision(photon_loc[i], sigma[i], median=median))
    return otp


def creer_matrice_et_moyennes(noms_fichiers, valeurs, positions):
    matrice_puits = np.full((len(positions),), float('nan'))
    valeurs_par_position = defaultdict(list)
    for nom_fichier, valeur in zip(noms_fichiers, valeurs):
        match = re.search(r'([A-H]\d+)', nom_fichier)
        if match:
            position = match.group(1)
            valeurs_par_position[position].append(valeur)

    for i, position in enumerate(positions):
        if position in valeurs_par_position:
            moyenne = sum(valeurs_par_position[position]) / len(valeurs_par_position[position])
            matrice_puits[i] = round(float(moyenne),2)

    return matrice_puits



def creer_data_boxplot(noms_fichiers, valeurs, positions):
    valeurs_par_position = defaultdict(list)
    
    # Ajouter les valeurs dans le defaultdict par position
    for nom_fichier, valeur in zip(noms_fichiers, valeurs):
        match = re.search(r'([A-H]\d+)', nom_fichier)
        if match:
            position = match.group(1)
            valeurs_par_position[position].append(valeur)

    matrice_puits = []

    # Concaténer les listes de valeurs pour chaque position
    for position in positions:
        if position in valeurs_par_position:
            concat_values = []
            for data in valeurs_par_position[position]:
                concat_values.extend(data)
            matrice_puits.append(concat_values)

    return matrice_puits


# Fonction pour convertir les coordonnees de la plaque (ex: A1) en indices de matrice (ex: (0, 0))
def convertir_coordonnees(coordonnees):
    colonne = ord(coordonnees[0]) - ord('A')
    ligne = int(coordonnees[1:]) - 1
    return ligne, colonne