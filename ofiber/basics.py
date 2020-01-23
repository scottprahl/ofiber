# pylint: disable=invalid-name


"""
Useful basic routines for optical fibers.

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


__all__ = ('acceptance_angle',
           'critical_angle',
           'cutoff_wavelength',
           'esi_delta',
           'esi_radius',
           'esi_v_parameter',
           'numerical_aperture',
           'numerical_aperture_graded_index',
           'relative_refractive_index',
           'numerical_aperture_from_Delta',
           'r_par',
           'r_per',
           'r_unpolarized',
           'v_parameter')


def acceptance_angle(NA):
    """
    Calculate the acceptance angle of an optical fiber.

    This is the half-angle in air.

    Args:
        NA : numerical aperture of the fiber  [--]

    Returns:
        maximum entrance/exit half-angle of the fiber [radians]
    """
    return np.arcsin(NA)


def critical_angle(n_core, n_clad):
    """
    Calculate the angle (from the normal) for total internal reflection.

    Args:
        n_core : the index of refraction of the fiber core  [--]
        n_core : the index of refraction of the fiber cladding  [--]

    Returns:
        angle of total internal reflection [radians]
    """
    return np.arcsin(n_clad / n_core)


def cutoff_wavelength(a, NA, ell=0, q=np.inf):
    """
    Calculate the cutoff wavelength for an optical fiber.

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
        q :   (optional) parameter for graded index fiber       [-]

    Returns:
        shortest wavelength for operation in the specified mode [m]
    """
    Vc, = jn_zeros(int(ell), 1)
    if np.isfinite(q):       # graded index fiber
        Vc *= np.sqrt(1 + 2 / q)
    return 2 * np.pi * a * NA / Vc


def esi_delta(Delta, q):
    """
    Calculate equivalent step index (esi) Delta for a graded-index fiber.

    Args:
        Delta :  relative refractive index     [-]

    Returns:
        equivalent relative refractive index   [-]
    """
    return q * (2 + q) / (1 + q)**2 * Delta


def esi_radius(a, q):
    """
    Calculate equivalent step index (esi) radius for a graded-index fiber.

    Args:
        a :   radius of the fiber                  [m]
        q :   parameter for graded index fiber     [-]

    Returns:
        equivalent step index radius               [m]
    """
    return a * (1 + q) / (2 + q)


def esi_v_parameter(V, q):
    """
    Calculate equivalent step index (esi) V for a graded-index fiber.

    Args:
        V :       V parameter                       [-]
        q :       parameter for graded index fiber  [-]

    Returns:
        equivalent step index V-parameter           [-]
    """
    return V * np.sqrt(q / (q + 2))


def numerical_aperture(n_core, n_clad):
    """
    Calculate the numerical aperture of an optical fiber.

    Args:
        n_core : the index of refraction of the fiber core      [-]
        n_clad : the index of refraction of the fiber cladding  [-]

    Returns:
        numerical aperture                                      [-]
    """
    return np.sqrt(n_core**2 - n_clad**2)


def numerical_aperture_from_Delta(n_core, Delta):
    """
    Calculate the numerical aperture of an optical fiber.

    Args:
        n_core : the index of refraction of the fiber core      [-]
        Delta : relative index of refraction                    [-]

    Returns:
        numerical aperture                                      [-]
    """
    return n_core * np.sqrt(2*Delta)


def numerical_aperture_graded_index(n_core, n_clad, q, r_over_a):
    """
    Calculate the numerical aperture of a graded-index optical fiber.

    The numerical aperture varies across the face of a graded-index fiber.
    This give the result at the fractional distance across the fiber core.

    Args:
        n_core :  the index of refraction of the fiber core      [-]
        n_clad :  the index of refraction of the fiber cladding  [-]
        q :       parameter for graded index fiber               [-]
        r_over_a : ratio of radius to the core radius            [-]

    Returns:
        numerical aperture at r_over_a                           [-]
    """
    return np.sqrt(n_core**2 - n_clad**2) * np.sqrt(1 - r_over_a**q)


def relative_refractive_index(n_core, n_clad):
    """
    Calculate the relative refractive index (Delta) for an optical fiber.

    Args:
        n_core :  the index of refraction of the fiber core      [-]
        n_clad:  the index of refraction of the fiber cladding   [-]

    Returns:
        the relative refractive index (Delta)                    [-]
    """
    return (n_core**2 - n_clad**2) / (2 * n_core**2)


def r_par(m, theta):
    """
    Calculate the Fresnel reflection for parallel polarized light.

    Args:
        m :     complex index of refraction   [-]
        theta : angle from normal to surface  [radians]

    Returns:
        reflected power                       [-]
    """
    m2 = m * m
    c = np.cos(theta)
    s = np.sin(theta)
    d = np.sqrt(m2 - s * s)
    return abs((m2 * c - d) / (m2 * c + d))**2


def r_per(m, theta):
    """
    Calculate the Fresnel reflection for perpendicular polarized light.

    Args:
        m :     complex index of refraction   [-]
        theta : angle from normal to surface  [radians]

    Returns:
        reflected power                       [-]
    """
    m2 = m * m
    c = np.cos(theta)
    s = np.sin(theta)
    d = np.sqrt(m2 - s * s)
    return abs((c - d) / (c + d))**2


def r_unpolarized(m, theta):
    """
    Calculate the Fresnel reflection for unpolarized incident light.

    Args:
        m :     complex index of refraction   [-]
        theta : angle from normal to surface  [radians]

    Returns:
        reflected power                       [-]
    """
    return (r_par(m, theta) + r_per(m, theta)) / 2


def v_parameter(a, NA, lambda0):
    """
    Calculate the V-parameter for an optical fiber.

    Args:
        a :       radius of the fiber              [m]
        NA :      numerical aperture of the fiber  [-]
        lambda0 : wavelength in vacuum             [m]

    Returns:
        V-parameter                                [-]
    """
    V = 2 * np.pi / lambda0 * a * NA
    return V
