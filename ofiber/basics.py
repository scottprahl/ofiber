"""
Useful basic routines for optical fibers

A handful of functions to calculate commonly needed optical fiber parameters.
Many are so simple that one might reasonably question if making them a
function has any utility.

Todo:
    * move MFD here?

Scott Prahl
Mar 2018
"""

import numpy as np
from scipy.special import jn_zeros


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
    """
    Calculates the cutoff wavelength for an optical fiber

    The default operation is for this function to calculate the cutoff
    wavelength for the fundamental mode of a step-index fiber.  The cutoff
    wavelength for higher order modes may be found by specifying a different
    value of ell.

    If the cutoff wavelength for a graded index fiber is desired, then specify
    a different value for q.

    Args:
        a :   radius of the fiber                               [m]
        NA :  numerical aperture of the fiber                   [-]
        ell : (optional) mode number                            [-]
        q :   (optional) parameter for graded index fiberr      [-]
    Returns:
        shortest wavelength for operation in the specified mode [m]
    """
    Vc, = jn_zeros(int(ell), 1)
    if np.isfinite(q):       # graded index fiber
        Vc *= np.sqrt(1 + 2 / q)
    return 2 * np.pi * a * NA / Vc


def numerical_aperture(n_core, n_clad):
    """
    Calculates the numerical aperture of an optical fiber

    Args:
        n_core : the index of refraction of the fiber core      [-]
        n_core : the index of refraction of the fiber cladding  [-]
    Returns:
        numerical aperture                                      [-]
    """
    return np.sqrt(n_core**2 - n_clad**2)


def numerical_aperture_graded_index(n_core, n_clad, q, r_over_a):
    """
    Calculates the numerical aperture of a graded-index optical fiber

    The numerical aperture varies across the face of a graded-index fiber.
    This give the result at the fractional distance across the fiber core.

    Args:
        n_core :  the index of refraction of the fiber core      [-]
        n_core :  the index of refraction of the fiber cladding  [-]
        q :       (optional) parameter for graded index fiber    [-]
        r_over_a : ratio of radius to the core radius            [-]
    Returns:
        numerical aperture at r_over_a                           [-]
    """
    return np.sqrt(n_core**2 - n_clad**2) * np.sqrt(1 - r_over_a**q)


def relative_refractive_index(n_core, n_clad):
    """
    Calculates the relative refractive index (Delta) for an optical fiber

    Args:
        n_core :  the index of refraction of the fiber core      [-]
        n_core :  the index of refraction of the fiber cladding  [-]
    Returns:
        the relative refractive index (Delta)                    [-]
    """
    return (n_core**2 - n_clad**2) / (2 * n_core**2)


def v_parameter(a, NA, lambda0):
    """
    Calculates the V-parameter for an optical fiber

    The default operation is for this function to calculate the cutoff
    wavelength for the fundamental mode of a step-index fiber.  The cutoff
    wavelength for higher order modes may be found by specifying a different
    value of ell.

    If the cutoff wavelength for a graded index fiber is desired, then specify
    a different value for q.

    Args:
        a :       radius of the fiber                 [m]
        NA :      numerical aperture of the fiber     [-]
        lambda0 : wavelength in vacuum                [m]
    Returns:
        V-parameter                                   [-]
    """
    return 2 * np.pi / lambda0 * a * NA
