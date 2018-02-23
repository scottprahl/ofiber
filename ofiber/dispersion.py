# Useful routines for cylindrical waveguides based on chapter 10 of
#
#    A. Ghatak, K. Thyagarajan, An Introduction to Fiber Optics,
#    Cambridge University Press, 1998
#
# Scott Prahl
# Feb 2018

from __future__ import division
import numpy as np
import ofiber.refraction

def approx_V_d2bV_by_V(V):
    """
    Returns approximate V*d^2(bV)/dV^2 when 1.4<V<2.4
    for the fundamental mode of a step index fiber
    Approximation by Marcuse (1979)
    """
    return 0.080 + 0.549*(2.834-V)**2
    
def waveguide_dispersion(n1,n2,a,lambda0):
    """
    Return the waveguide dispersion of the fundamental
    mode of a step index fiber (core n1, cladding n2) with
    radius a [m] and wavelength lambda0 [m]
    result has units of s/m**2  (multiply by 1e6 to get ps/(km*nm))
    """
    Delta = (n1**2-n2**2)/2/n1**2
    V = 2*np.pi/lambda0 * a * np.sqrt(n1**2-n2**2)
    c = 3e8
    disp = -n2*Delta/c/lambda0 * approx_V_d2bV_by_V(V)
    return disp


def waveguide_dispersion_delta(glass,Delta,a,lambda0):
    """
    Returns waveguide dispersion for glass with refractive index 

    glass is an array of 6 Sellmeier coefficients [um**2 or 1/um**2]
    Delta is the refractive index difference between core and cladding
    a is the fiber radius [m]
    lambda0 is the wavelength in [m]
    returned dispersion is in [s/m**2]
    """
    n1 = sellmeier_refractive_index(glass,lambda0)
    n2 = n1*(1-Delta)
    if np.isscalar(lambda0):
        yy = waveguide_dispersion(n1,n2,a,lambda0)
    else :
        yy = np.empty_like(lambda0)
        for i in range(len(yy)):
            yy[i] = waveguide_dispersion(n1[i],n2[i],a,lambda0[i])
    return yy

def total_dispersion(glass,Delta,a,lambda0):
    """
    Returns total dispersion for (x)GeO2:(1-x)SiO2 glass 

    glass is an array of 6 Sellmeier coefficients [um**2 or 1/um**2]
    Delta is the refractive index difference between core and cladding
    a is the fiber radius [m]
    lambda0 is the wavelength in [m]
    returned dispersion is in [s/m**2]
    """
    Dm = sellmeier_material_dispersion(glass,lambda0)      # [s/m**2]
    Dw = waveguide_dispersion_delta(glass,Delta,a,lambda0) # [s/m**2]
    return Dw+Dm

def total_dispersion_geo2_sio2_mixture(x, *args):
    """
    Returns total dispersion for (x)GeO2:(1-x)SiO2 glass 

    wrapper function for use with brentq()
    """
    a = args[0]
    Delta = args[1]
    lambda0 = args[2]
    glass = doped_glass(x)
    return total_dispersion(glass,Delta,a,lambda0)

