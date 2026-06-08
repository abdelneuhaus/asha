import os
import math

from src.io_utils import get_PALMTracer_files, get_poca_files, get_statMIA_files, read_statMIA, read_poca_files, read_locPALMTracer_file 


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


def get_size(start_path='.'):
    """
    Calculates the total size of a directory in bytes, separated by specific file types.

    Args:
        start_path (str, optional): The directory path to calculate. 
            Defaults to the current directory ('.').

    Returns:
        tuple: A tuple containing:
            - total_size (int): Total size of all files.
            - smf_size (int): Size of image files (.smf).
            - processed_size (int): Size of processed data files (.csv, .txt).
    """
    total_size = 0
    smf_size = 0
    processed_size = 0
    
    for dirpath, dirnames, filenames in os.walk(start_path):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            if not os.path.islink(fp):
                file_size = os.path.getsize(fp)
                total_size += file_size
                
                # Check extension to categorize
                ext = os.path.splitext(f)[1].lower()
                if ext == '.smf':
                    smf_size += file_size
                elif ext in ['.csv', '.txt']:
                    processed_size += file_size
                    
    return total_size, smf_size, processed_size


def count_total_locs_and_clusters(plate_path):
    """
    Count the number of localizations and clusters found after processing the HCS plate,
    along with disk usage detailed by file type.

    Args:
        plate_path (str): Path to the plate directory.

    Returns:
        tuple: A tuple containing:
            - total_localizations (int): Number of localizations in one plate.
            - total_clusters (int): Number of clusters in one plate.
            - total_size (str): Total disk size.
            - smf_size (str): Disk size occupied by .smf images.
            - processed_size (str): Disk size occupied by .csv and .txt files.
    """
    list_of_mia_files = get_statMIA_files(plate_path)
    list_of_pt_files = get_PALMTracer_files(plate_path)
    list_of_poca_files = get_poca_files(plate_path)
    
    total_localizations_pre = sum(len(read_statMIA(f)) for f in list_of_mia_files)
    total_localizations = sum(len(read_locPALMTracer_file(f)) for f in list_of_pt_files)
    total_clusters = sum(len(read_poca_files(f)) for f in list_of_poca_files)
    total_bytes, smf_bytes, processed_bytes = get_size(plate_path)
    total_size = convert_size(total_bytes)
    smf_size = convert_size(smf_bytes)
    processed_size = convert_size(processed_bytes)
    
    return total_localizations_pre, total_localizations, total_clusters, total_size, smf_size, processed_size