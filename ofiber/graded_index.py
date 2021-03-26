# pylint: disable=invalid-name
# pylint: disable=too-many-arguments

"""
Useful basic routines for graded-index optical fibers.

See <https://ofiber.readthedocs.io> for usage examples.

This needs more testing.
"""
import scipy.constants
import numpy as np


__all__ = ('power_law_profile',
           'first_derivative',
           'curvature',
           'transverse_location')


def power_law_profile(ncore, nclad, q, a, x):
    """Calculate the index of refraction at a particular radius."""
    delta = (ncore**2 - nclad**2) / 2 / ncore**2
    nsqr = ncore**2 * (1 - 2 * delta * (np.abs(x) / a)**q)
    if not np.isscalar(x):
        np.place(nsqr, abs(x) >= a, nclad**2)
    return np.sqrt(nsqr)


def first_derivative(x, f):
    """Return the first derivative of the array f w.r.t. x but the same length as f."""
    deriv = np.diff(f, n=1)
    deriv = np.insert(deriv, -1, 2 * deriv[-1] - deriv[-2])
    dx = x[1] - x[2]
    return deriv / dx


def curvature(ncore, nclad, q, a, x, theta):
    """Return the curvature at a position x on the fiber profile."""
    beta = ncore * np.cos(theta)
    nsqr = power_law_profile(ncore, nclad, q, a, x)**2
    curve = first_derivative(x, nsqr) / 2 / beta**2
    return curve


def transverse_location(n1, theta1, Delta, a, z):
    """Equation 4.13 from Ghatak."""
    beta = n1 * np.cos(theta1)
    Gamma = n1 * np.sqrt(2 * Delta) / beta / a
    A = a * np.sin(theta1) / np.sqrt(2 * Delta)
    return A * np.sin(Gamma * z)

def velocity(ncore, q, beta_invariant, material_dispersion=None):
    """
    The velocity for a power-law circular fiber.
    
    Equation 5.4 in Ghatak.
    """
    c = scipy.constants.speed_of_light
    if glass is None:
        A = 2/c/(2+q)
        B = q * ncore**2/c/(2+q)
    else:
        N1 = ncore + material_dispersion
        y = 2* ncore/N1
        A = 2 * N1/ncore  * (1+0.25*y)/c/(q+2)
        B = q *ncore**2*A-1/4/c*N1*ncore*y
    
    return A*beta_invariant + B/beta_invariant
    return (A * beta_invariant + B/beta_invariant)