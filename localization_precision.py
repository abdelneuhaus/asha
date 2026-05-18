import math

def localization_precision(photon, sigma, a=160, median=0, bckg=0):
    """Implementation of Localization Precision (nm) for EMCCD Camera (Shot noise limited)

    Args:
        photon: number of photon (count)
        sigma: standard deviation (in nm)
        a: pixel size (nm)
    """
    if photon == 0:
        photon = median
    return sigma/math.sqrt(photon)