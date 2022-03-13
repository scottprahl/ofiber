# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=non-ascii-name
"""
Useful routines for dispersion in cylindrical waveguides.

See <https://ofiber.readthedocs.io> for usage examples.

Based on chapter 10 of A. Ghatak, K. Thyagarajan, An Introduction to Fiber
Optics, Cambridge University Press, 1998
"""

import scipy.constants
import numpy as np
import ofiber.refraction as ofr
import ofiber.cylinder_step as ofc
import ofiber.basics as ofb


__all__ = ('Material_Dispersion',
           'Waveguide_Dispersion',
           'Dispersion',
           )


def Material_Dispersion(core, λ):
    """
    Calculate the material dispersion using Sellmeier coefficients.

    Args:
        core: array of Sellmeier coefficients for core
        λ:    wavelength [m]

    Returns:
        material dispersion [s/m**2]   (multiply by 1e6 to get ps/(km*nm))
    """
    c = scipy.constants.speed_of_light
    return -λ * ofr.d2n(core, λ) / c


def Waveguide_Dispersion(n_core, n_clad, r_core, λ, q=1e20, approx=False):
    """
    Calculate the waveguide dispersion of a fiber.

    The default value of q represents a step index fiber.  Other values
    allow parabolic (q=2) or triangular (q=1) profiles.

    The waveguide dispersion is for the fundamental mode of the fiber.

    The approximation is reasonably good for values from 1.6<V<2.6

    Args:
        n_core: core index of refraction          [-]
        n_clad: cladding index of refraction      [-]
        r_core: radius of the fiber core          [m]
        λ:      wavelength in vacuum              [m]
        q:      power in graded index fiber       [-]
        approx: approximate when True             [True/False]

    Returns:
        waveguide dispersion [s/m**2]   (multiply by 1e6 to get ps/(km*nm))
    """
    c = scipy.constants.speed_of_light
    Δ = (n_core**2 - n_clad**2) / 2 / n_core**2
    V = 2 * np.pi / λ * r_core * np.sqrt(n_core**2 - n_clad**2)

    # Find the equivalent step index fiber parameters
    esi_Delta = ofb.esi_Delta(Δ, q)
    esi_V = ofb.esi_V_parameter(V, q)

    if approx:
        vtemp = ofc.V_d2bV_by_V_Approx(esi_V)
    else:
        vtemp = ofc.V_d2bV_by_V(esi_V, 0)
    return -n_clad * esi_Delta / c / λ * vtemp


def Dispersion(core, n_clad, r_core, λ, q=1e20, approx=False):
    """
    Calculate the material and waveguide dispersion.

    This is a convenience routine that finds the total dispersion for
    a specific type of core glass and refractive index difference.

    The returned dispersion is in units of [s/m**2].  To convert to
    [ps/km/nm], multiply by 1e6.

    Args:
        core:   Sellmeier coefficients for core   [various]
        n_clad: index of cladding                 [-]
        r_core: radius of fiber core              [m]
        λ:      wavelength in vacuum              [m]
        q:      power in graded index fiber       [-]
        approx: approximate when True             [True/False]

    Returns:
        material, waveguide, and total dispersion [s/m**2]
    """
    n_core = ofr.n(core, λ)
    Dm = Material_Dispersion(core, λ)
    Dw = Waveguide_Dispersion(n_core, n_clad, r_core, λ, q, approx)
    return Dm, Dw
