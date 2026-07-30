# pylint: disable=invalid-name
# pylint: disable=too-many-arguments

"""
Useful basic routines for graded-index optical fibers.

See <https://ofiber.readthedocs.io> for usage examples.

This needs more testing.
"""
import scipy.constants
import numpy as np


__all__ = ("power_law_profile", "first_derivative", "curvature", "transverse_location", "ray_delay")


def power_law_profile(n_core, nclad, q, a, x):
    """
    Calculate the index of refraction at a particular radius.

    Outside the core the index is the cladding index.  Clamping |x|/a at 1
    delivers that for scalars and arrays alike, because n_core**2*(1-2*delta)
    is exactly nclad**2, and it avoids overflow when q is large.

    Args:
        n_core: index on the fiber axis           [-]
        nclad:  index of the cladding             [-]
        q:      power in the graded index profile [-]
        a:      radius of the fiber core          [m]
        x:      radial position                   [m]

    Returns:
        index of refraction at x                  [-]
    """
    delta = (n_core**2 - nclad**2) / 2 / n_core**2
    r_over_a = np.minimum(np.abs(x) / a, 1.0)
    nsqr = n_core**2 * (1 - 2 * delta * r_over_a**q)
    return np.sqrt(nsqr)


def first_derivative(x, f):
    """
    Return the first derivative of the array f w.r.t. x but the same length as f.

    This is a forward difference, so element i approximates df/dx midway between
    x[i] and x[i+1] and carries an O(dx) offset.  The final element is a linear
    extrapolation of the last two differences, which keeps the result the same
    length as f.  The spacing of x is assumed to be uniform.

    Args:
        x: uniformly spaced positions          [m]
        f: values of the function at each x

    Returns:
        df/dx, the same length as f
    """
    deriv = np.diff(f, n=1)
    deriv = np.append(deriv, 2 * deriv[-1] - deriv[-2])
    dx = x[1] - x[0]
    return deriv / dx


def curvature(n_core, nclad, q, a, x, theta):
    """Return the curvature at a position x on the fiber profile."""
    beta = n_core * np.cos(theta)
    nsqr = power_law_profile(n_core, nclad, q, a, x) ** 2
    curve = first_derivative(x, nsqr) / 2 / beta**2
    return curve


def transverse_location(n1, theta1, Delta, a, z):
    """Equation 4.13 from Ghatak."""
    beta = n1 * np.cos(theta1)
    Gamma = n1 * np.sqrt(2 * Delta) / beta / a
    A = a * np.sin(theta1) / np.sqrt(2 * Delta)
    return A * np.sin(Gamma * z)


def ray_delay(n_core, q, beta_invariant, material_dispersion=None, profile_dispersion=0):
    """
    Find the delay per unit length for a ray in a power-law circular fiber.

    Equation 5.4 in Ghatak.  This is a time per unit length in [s/m], the
    reciprocal of a velocity: an axial ray comes back as N1/c whatever the
    profile.

    A ray of invariant beta accumulates delay at a rate N(r)*n(r)/(beta*c) per
    unit length, and for a power-law profile the path average of n**2 is
    (2*beta**2 + q*n_core**2)/(2+q).  Writing the group index profile as

        N(r) = (N1/n_core) * [(1 - P) * n(r) + n_core * P]

    and using <n> = (n_core**2 + <n**2>)/(2*n_core), which holds to first order
    in Delta, collects the delay into A*beta + B/beta with the coefficients
    below.  The exponent that equalises the axial and edge-ray delays then
    falls out as q = 2 - 2P.

    `material_dispersion` is N1 - n_core, the amount by which the group index
    on axis exceeds the phase index, so passing 0 is the same as omitting it.
    `profile_dispersion` is P, which needs dDelta/dlambda and so cannot be
    recovered from the material dispersion alone.  Beware that other texts use
    a parameter twice this size, and of opposite sign, for the same quantity.

    Args:
        n_core: index on the fiber axis                        [-]
        q:      power in the graded index profile              [-]
        beta_invariant: ray invariant n(r)*cos(theta(r))       [-]
        material_dispersion: N1 - n_core, defaults to none     [-]
        profile_dispersion: P = (n_core/N1)*(lambda/Delta)*(dDelta/dlambda) [-]

    Returns:
        delay per unit length                                  [s/m]
    """
    c = scipy.constants.speed_of_light
    N1 = n_core if material_dispersion is None else n_core + material_dispersion
    P = profile_dispersion

    A = N1 / n_core * (2 - P) / c / (2 + q)
    B = N1 * n_core * (q + P) / c / (2 + q)

    return A * beta_invariant + B / beta_invariant
