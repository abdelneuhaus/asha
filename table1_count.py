from utils import get_PALMTracer_files, get_poca_files, read_poca_files, read_locPALMTracer_file
import numpy as np

def table1_count(plate_path):
    """ Count the number of localization and clusters per plate"""
    list_of_pt_files = get_PALMTracer_files(plate_path)
    list_of_poca_files = get_poca_files(plate_path)
    total_localizations = 0
    total_clusters = 0
    for pt_file in list_of_pt_files:
        raw_file_pt = read_locPALMTracer_file(pt_file)
        total_localizations += len(raw_file_pt)
    for poca_file in list_of_poca_files:
        raw_file_poca = read_poca_files(poca_file)
        total_clusters += len(raw_file_poca)
    return total_localizations, total_clusters


import pandas as pd
import os

def average_density(plate_path):
    """Calculate average localization density per frame and mean length from poca files for wells C4, D8, and E4."""

    wells_of_interest = ['C4', 'D8', 'E4']
    surface = (256 * 0.16) ** 2
    list_of_pt_files = [f for f in get_PALMTracer_files(plate_path) if any(well in f for well in wells_of_interest)]
    list_of_poca_files = [f for f in get_poca_files(plate_path) if any(well in f for well in wells_of_interest)]
    total_density = 0
    total_frames = 0

    for pt_file in list_of_pt_files:
        df = read_locPALMTracer_file(pt_file)  # ou pd.read_table selon le format
        grouped = df.groupby("Plane")
        for _, group in grouped:
            loc_count = len(group)
            density = loc_count / surface
            total_density += density
            total_frames += 1
    average_density_loc = total_density / total_frames if total_frames > 0 else 0
    
    # Pour les poca files
    individual_means = []
    for poca_file in list_of_poca_files:
        df = read_poca_files(poca_file)
        df = df[df['total ON'] <= np.quantile(df["total ON"], 1)]
        df = df[df['intensity'] <= np.quantile(df["intensity"], 1)]
        length = len(df)# ou pd.read_csv / pd.read_table selon le format
        individual_means.append(length)
    average_density_cluster = sum(individual_means) / len(individual_means) if individual_means else 0

    return average_density_loc, average_density_cluster


plate2bis = "D:/ANALYSIS_PAPER/new_threshold/240417_W1_FPs"
plate3bis = "D:/ANALYSIS_PAPER/new_threshold/240924_W1_FPs"
plate4bis = "D:/ANALYSIS_PAPER/new_threshold/241113_W1_FPs"

# print(table1_count(plate2bis))
# print(table1_count(plate3bis))
# print(table1_count(plate4bis))
# print(average_density(plate2bis))
# print(average_density(plate3bis))
# print(average_density(plate4bis))