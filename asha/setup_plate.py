from asha.utils import get_PALMTracer_files, get_poca_files


def setup_plate(repertory, prefix='SR_001.MIA'):
    """
    Set up the plate for analysis by organizing file paths.

    Args:
        repertory (str): The directory containing the files.

    Returns:
        tuple: A tuple containing two lists - the first list contains the POCA files,
               and the second list contains the PALMTracer files.
    """
    well_groups = [["C3", "D9", "E3"], 
                   ["C4", "D8", "E4"], 
                   ["C5", "D7", "E5"], 
                   ["C6", "D6", "E6"], 
                   ["C7", "D5", "E7"], 
                   ["C8", "D4", "E8"], 
                   ["C9", "D3", "E9"]]
    
    poca_files = get_poca_files(repertory, prefix)
    palm_tracer_files = get_PALMTracer_files(repertory, prefix)
    group_poca_files = []
    group_palm_tracer_files = []
    
    for wells in well_groups:
        group_poca_files.append([f for f in poca_files if any(well in f for well in wells)])
        group_palm_tracer_files.append([f for f in palm_tracer_files if any(well in f for well in wells)])

    return group_poca_files, group_palm_tracer_files