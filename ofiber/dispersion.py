"""
Useful routines for cylindrical waveguides

   Based on chapter 10 of A. Ghatak, K. Thyagarajan, An Introduction to Fiber
   Optics, Cambridge University Press, 1998

Todo:
    * adhere to google function docstrings
    * lowercase function names?

"""

import numpy as np
import ofiber.refraction as ofr
import ofiber.cylinder_step as ofc


__all__ = ['Material_Dispersion',
           'Waveguide_Dispersion',
           'V_d2bV_by_V_Approx',
           'Waveguide_Dispersion_Approx',
           'Waveguide_Dispersion_Delta',
           'Total_Dispersion']


def Material_Dispersion(glass, lambda0):
    """
    Return the material dispersion using Sellmeier coefficients

    Args:
        glass  : an array of Sellmeier coefficients.
        lambda0: wavelength [m]
    Returns:
        material dispersion [s/m**2]   (multiply by 1e6 to get ps/(km*nm))
    """
    c = 2.997e8                                    # m/s
    return -lambda0 * ofr.d2n(glass, lambda0) / c  # s/m**2


def Waveguide_Dispersion(n1, n2, a, lambda0):
    """
    Return the waveguide dispersion of a step index fiber.

    The waveguide dispersion is for the fundamental mode of the fiber.

    Args:
        n1:      core index of refraction     [--]
        n2:      cladding index of refraction [--]
        a:       radius of the fiber          [m]
        lambda0: wavelength                   [m]
    Returns:
        waveguide dispersion [s/m**2]   (multiply by 1e6 to get ps/(km*nm))
    """
    Delta = (n1**2 - n2**2) / 2 / n1**2
    V = 2 * np.pi / lambda0 * a * np.sqrt(n1**2 - n2**2)
    c = 2.997e8
    dw = -n2 * Delta / c / lambda0 * ofc.V_d2bV_by_V(V, 0)
    return dw


def V_d2bV_by_V_Approx(V):
    """
    Returns approximate V*d^2(bV)/dV^2 when 1.4<V<2.4
    for the fundamental mode of a step index fiber
    Approximation by Marcuse (1979)
    """
    return 0.080 + 0.549 * (2.834 - V)**2


def Waveguide_Dispersion_Approx(n1, n2, a, lambda0):
    """
    Return the waveguide dispersion of the fundamental
    mode of a step index fiber (core n1, cladding n2) with
    radius a [m] and wavelength lambda0 [m]
    returned dispersion is in [s/m**2]   (multiply by 1e6 to get ps/(km*nm))
    """
    Delta = (n1**2 - n2**2) / 2 / n1**2
    V = 2 * np.pi / lambda0 * a * np.sqrt(n1**2 - n2**2)
    c = 2.997e8
    dw = -n2 * Delta / c / lambda0 * V_d2bV_by_V_Approx(V)
    return dw


def Waveguide_Dispersion_Delta(glass, Delta, a, lambda0):
    """
    Returns waveguide dispersion for glass with refractive index

    glass[0:3] is in [um**2] and glass[3:6] is in [1/um**2]
    Delta is the refractive index difference between core and cladding
    a is the fiber radius [m]
    lambda0 is the wavelength in [m]
    returned dispersion is in [s/m**2]   (multiply by 1e6 to get ps/(km*nm))
    """
    n1 = ofr.n(glass, lambda0)
    n2 = n1 * (1 - Delta)
    return Waveguide_Dispersion(n1, n2, a, lambda0)


def Total_Dispersion(glass, Delta, a, lambda0):
    """
    Returns total dispersion for specified glass

    glass[0:3] is in [um**2] and glass[3:6] is in [1/um**2]
    Delta is the refractive index difference between core and cladding
    a is the fiber radius [m]
    lambda0 is the wavelength in [m]
    returned dispersion is in [s/m**2]   (multiply by 1e6 to get ps/(km*nm))
    """
    Dm = Material_Dispersion(glass, lambda0)                   # [s/m**2]
    Dw = Waveguide_Dispersion_Delta(glass, Delta, a, lambda0)  # [s/m**2]
    return Dw + Dm
