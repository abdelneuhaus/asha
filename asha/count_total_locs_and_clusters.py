from utils import get_PALMTracer_files, get_poca_files, read_poca_files, read_locPALMTracer_file
import numpy as np

def count_total_locs_and_clusters(plate_path):
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


plate = "DATA/W1/"
print(count_total_locs_and_clusters(plate))