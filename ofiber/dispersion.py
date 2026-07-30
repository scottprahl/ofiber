# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=non-ascii-name
"""
Useful routines for dispersion in cylindrical waveguides.

See <https://ofiber.readthedocs.io> for usage examples.

The three routines are::

    Material_Dispersion(core, λ)
    Waveguide_Dispersion(n_core, n_clad, r_core, λ, q=np.inf, approx=False)
    Dispersion(core, n_clad, r_core, λ, q=np.inf, approx=False)

All return [s/m**2]; multiply by 1e6 for the more familiar [ps/(km*nm)].
The total dispersion is the sum of the material and waveguide terms, and
`Dispersion` returns the two separately so they can be plotted apart.

A step-index fiber is signalled by q=np.inf, the same sentinel
`ofiber.cutoff_wavelength` uses.

Based on chapter 10 of A. Ghatak, K. Thyagarajan, An Introduction to Fiber
Optics, Cambridge University Press, 1998
"""

import numpy as np
import scipy.constants
import ofiber.refraction as ofr
import ofiber.cylinder_step as ofc
import ofiber.basics as ofb

__all__ = (
    "Material_Dispersion",
    "Waveguide_Dispersion",
    "Dispersion",
)


def Material_Dispersion(core, λ):
    """
    Calculate the material dispersion using Sellmeier coefficients.

    This is -λ/c times the second derivative of the index, so it vanishes
    where the group index is stationary.  For pure silica that is near
    1272 nm, and adding GeO2 moves it to longer wavelengths.

    Args:
        core: Sellmeier coefficients, from `ofiber.glass` or `ofiber.doped_glass`
        λ:    wavelength [m]

    Returns:
        material dispersion [s/m**2]   (multiply by 1e6 to get ps/(km*nm))
    """
    c = scipy.constants.speed_of_light
    return -λ * ofr.d2n(core, λ) / c


def Waveguide_Dispersion(n_core, n_clad, r_core, λ, q=np.inf, approx=False):
    """
    Calculate the waveguide dispersion of a fiber.

    The default value of q represents a step index fiber, matching the
    sentinel `ofiber.cutoff_wavelength` uses.  Other values allow parabolic
    (q=2) or triangular (q=1) profiles.

    The waveguide dispersion is for the fundamental mode of the fiber, and is
    negative for an ordinary fiber, so it pulls the zero of the total
    dispersion to a longer wavelength than the material term alone.

    Setting approx=True swaps the exact V*d2(bV)/dV2 for Marcuse's quadratic.
    That is within about 1% for 1.35<V<2.08 but drifts to roughly 6% near
    V=2.44 before crossing back, so the endpoints of 1.6<V<2.6 look better
    than the middle.

    The two indices are used as given, so any wavelength dependence of the
    core and cladding belongs to the caller; `Dispersion` evaluates them from
    the Sellmeier coefficients at each wavelength.

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
    Δ = ofb.relative_refractive_index(n_core, n_clad)
    V = ofb.V_parameter(r_core, ofb.numerical_aperture(n_core, n_clad), λ)

    # Find the equivalent step index fiber parameters
    esi_Delta = ofb.esi_Delta(Δ, q)
    esi_V = ofb.esi_V_parameter(V, q)

    if approx:
        vtemp = ofc.V_d2bV_by_V_Approx(esi_V)
    else:
        vtemp = ofc.V_d2bV_by_V(esi_V, 0)
    return -n_clad * esi_Delta / c / λ * vtemp


def Dispersion(core, n_clad, r_core, λ, q=np.inf, approx=False):
    """
    Calculate the material and waveguide dispersion.

    This is a convenience routine that finds the total dispersion for
    a specific type of core glass and refractive index difference.  The core
    index is evaluated from the Sellmeier coefficients at every wavelength,
    unlike `Waveguide_Dispersion`, which takes the indices as given.

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
        material and waveguide dispersion, as a pair; their sum is the
        total dispersion                          [s/m**2]
    """
    n_core = ofr.n(core, λ)
    Dm = Material_Dispersion(core, λ)
    Dw = Waveguide_Dispersion(n_core, n_clad, r_core, λ, q, approx)
    return Dm, Dw
