# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals

"""
Planar parabolic graded-index optical waveguides.

See <https://ofiber.readthedocs.io> for usage examples.

A graded-index planar waveguide is an infinite waveguide in which confinement
takes place in the x-direction because of a parabolic index of refraction
gradient.

Let z be the direction of light propagation through the waveguide.  Let x be the
direction in which the index of refraction varies, and then the y-direction is
perpendicular to both x and z.  The index of refraction n(x) is given by::

    n(x)**2 = n_1**2 * ( 1 - 2*Δ*(x/a)**2 )

where n_1 is the index of refraction in the center of the waveguide (x=0)
and Δ is the usual relative refractive index.  In this case there are not
cladding layers, and instead, the light is guided by the refractive index
variation.

Based on section 7.9 of Ghatak and Thyagarajan, *An Introduction to
Fiber Optics*, Cambridge University Press, 1998.
"""

import numpy as np
import ofiber.basics

__all__ = ('parabolic_propagation_constant',
           'parabolic_propagation_constants',
           'TE_planar_parabolic_field',
           )


def _herm(n, x):
    """
    Calculate the Hermite polynomial H_n(x).

    Args:
        n: order of the Hermite polynomial    [-]
        x: argument of the Hermite polynomial [-]

    Returns:
        H_n(x)
    """
    c = np.zeros(n + 1, dtype="i4")
    c[n] = 1
    return np.polynomial.hermite.hermval(x, c)


def parabolic_propagation_constant(m, lambda0, n1, a, V):
    """
    Calculate the waveguide propagation constant beta.

    Args:
        m: mode number                                  [-]
        lambda0: wavelength in vacuum                   [m]
        n1: centerline index of refraction of waveguide [-]
        Delta: relative refractive index                [-]
        a: half thickness of the waveguide              [m]

    Returns:
        beta for the mth mode
    """
    gamma = np.sqrt(V) / a
    k = 2 * np.pi * n1 / lambda0
    return np.sqrt(k**2 - gamma**2 * (2 * m + 1))


def parabolic_propagation_constants(lambda0, n1, a, V):
    """
    Return all the betas for a parabolic planar waveguide.

    Args:
        lambda0: wavelength in vacuum                   [m]
        n1: centerline index of refraction of waveguide [-]
        a: half thickness of the waveguide              [m]
        V: the V-parameter for the waveguide
    Returns:
        an array eigenvalues.
    """
    Delta = np.sqrt(V / n1) / 2
    modes = np.floor(V / Delta - 1 / 2)
    betas = np.empty(modes)

    for mode in range(modes):
        betas[mode] = parabolic_propagation_constant(mode, lambda0, n1, a, V)

    return betas


def TE_planar_parabolic_field(m, lambda0, n1, Delta, a, x):
    """
    Calculate the transverse electric field E_y(x) for mode m.

    This is normalized such that the integral of E_y(x)**2 over
    all x is one.

    Args:
        m: mode number                                  [-]
        lambda0: wavelength in vacuum                   [m]
        n1: centerline index of refraction of waveguide [-]
        Delta: relative refractive index                [-]
        a: half thickness of the waveguide                   [m]
        x: position in the waveguide                    [m]

    Returns:
         E_y(x) the transverse electric field at x      [m**-0.5]
    """
    NA = ofiber.basics.numerical_aperture_from_Delta(n1, Delta)
    V = ofiber.basics.V_parameter(a, NA, lambda0)
    gamma = np.sqrt(V) / a
    Nm = np.sqrt(gamma / (2**m * np.math.factorial(m) * np.sqrt(np.pi)))
    xi = gamma * x
    E_y = Nm * _herm(m, xi) * np.exp(-0.5 * xi**2)
    return E_y
