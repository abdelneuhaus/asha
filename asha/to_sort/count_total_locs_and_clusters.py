from utils import get_PALMTracer_files, get_poca_files, read_poca_files, read_locPALMTracer_file


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
    
    return total_localizations, total_clusters