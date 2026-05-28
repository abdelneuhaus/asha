import os
import math

from src.io_utils import get_PALMTracer_files, get_poca_files, read_poca_files, read_locPALMTracer_file



def convert_size(size_bytes):
    """
    Converts a size in bytes into a human-readable string format. Inspired from https://stackoverflow.com/a/1392549

    Args:
        size_bytes (int or float): The total number of bytes.

    Returns:
        str: A human-readable string representation (e.g., "1.25 MB").
    """
    if size_bytes <= 0:
        return "0 B"
    size_name = ("B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB")
    
    i = int(math.floor(math.log(size_bytes, 1024)))
    i = min(i, len(size_name) - 1)
    p = math.pow(1024, i)
    s = round(size_bytes / p, 2)
    return f"{s} {size_name[i]}"


def get_size(start_path = '.'):
    """
    Calculates the total size of a directory in bytes.

    Args:
        start_path (str, optional): The directory path to calculate. 
            Defaults to the current directory ('.').

    Returns:
        int: The total size of all files in bytes.
    """
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            if not os.path.islink(fp):
                total_size += os.path.getsize(fp)
    return total_size


def count_total_locs_and_clusters(plate_path):
    """
    Count the number of localizations and clusters found after processing the HCS plate.

    Args:
        plate_path (str): Path to the plate directory.

    Returns:
        tuple: A tuple containing:
            - total_localizations (int): Number of localizations in one plate.
            - total_clusters (int): Number of clusters in one plate.
    """
    
    list_of_pt_files = get_PALMTracer_files(plate_path)
    list_of_poca_files = get_poca_files(plate_path)
    total_localizations = sum(len(read_locPALMTracer_file(f)) for f in list_of_pt_files)
    total_clusters = sum(len(read_poca_files(f)) for f in list_of_poca_files)
    total_size = convert_size(get_size(plate_path))
    
    return total_localizations, total_clusters, total_size