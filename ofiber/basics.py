"""
Useful basic routines for optical fibers

A handful of functions to calculate commonly needed optical fiber parameters.  Many are so simple that one might reasonably question if making them a function has any utility.  

Todo:
    * add proper documentation
    * move MFD here?

Scott Prahl
Mar 2018
"""

import numpy as np

__all__ = ['acceptance_angle',
           'critical_angle',
           'cutoff_wavelength',
           'numerical_aperture',
           'numerical_aperture_graded_index',
           'relative_refractive_index',
           'v_parameter']


def acceptance_angle(NA):
    """
    Calculates the acceptance angle of an optical fiber
    
    This is the half-angle in air.
    
    Args:
        NA : numerical aperture of the fiber  [--]
    Returns:
        maximum entrance/exit half-angle of the fiber [radians]
    """
    return np.arcsin(NA)


def critical_angle(n_core, n_clad):
    """
    Calculates the angle (from the normal) for total internal reflection
        
    Args:
        n_core : the index of refraction of the fiber core  [--]
        n_core : the index of refraction of the fiber cladding  [--]
    Returns:
        angle of total internal reflection [radians]
    """
    return np.arcsin(n_clad / n_core)


def cutoff_wavelength(a, NA, ell=0, q=np.inf):
    Vc = jn_zeros(ell, 1)
    if np.isfinite(q):       # graded index fiber
        Vc *= np.sqrt(1 + 2 / q)
    return 2 * np.pi * a * NA / Vc


def numerical_aperture(n_core, n_clad):
    return np.sqrt(n_core**2 - n_clad**2)


def numerical_aperture_graded_index(n_core, n_clad, q, r, a):
    return np.sqrt(n_core**2 - n_clad**2) * np.sqrt(1 - (r / a)**q)


def relative_refractive_index(n_core, n_clad):
    return (n_core**2 - n_clad**2) / (2 * n_core)


def v_parameter(a, NA, lambda0):
    return 2 * np.pi / lambda0 * a * NA
