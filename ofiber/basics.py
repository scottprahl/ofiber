# pylint: disable=invalid-name


"""
Calculate simple optical fiber parameters.

See <https://ofiber.readthedocs.io> for usage examples.

Basic parameters that can be found are::

    acceptance_angle(NA, n_outside=1)
    critical_angle(n_core, n_clad)
    cutoff_wavelength(a, NA, ell=0, q=np.inf)
    numerical_aperture(n_core, n_clad)
    numerical_aperture_from_Delta(n_core, Delta)
    relative_refractive_index(n_core, n_clad)
    V_parameter(a, NA, lambda0)

If you want Δ (Delta), then use `relative_refractive_index`

Some are just generic Fresnel equations::

    R_par(m, theta)
    R_per(m, theta)
    R_unpolarized(m, theta)

And finally, some apply to graded index fibers where 'esi' is short for
'equivalent step index'::

    esi_Delta(Delta, q)
    esi_radius(a, q)
    esi_V_parameter(V, q)
    numerical_aperture_graded_index(n_core, n_clad, q, r_over_a)

The three esi routines agree with one another: building a V-parameter from
`esi_radius` and `esi_Delta` reproduces `esi_V_parameter` exactly.
"""

import numpy as np
from scipy.special import jn_zeros


__all__ = (
    "acceptance_angle",
    "critical_angle",
    "cutoff_wavelength",
    "esi_Delta",
    "esi_radius",
    "esi_V_parameter",
    "numerical_aperture",
    "numerical_aperture_graded_index",
    "relative_refractive_index",
    "numerical_aperture_from_Delta",
    "R_par",
    "R_per",
    "R_unpolarized",
    "V_parameter",
)


def acceptance_angle(NA, n_outside=1):
    """
    Find the acceptance angle for a cone of light in/out of an optical fiber.

    This is the half-angle measured from the normal to the fiber face
    to the edge of the entering (or exiting) cone of light.

    The face of the optical fiber is in a medium that defaults to
    air, but whose index can be specified.

    If NA exceeds n_outside then no such cone exists and nan is returned.

    Args:
        NA : numerical aperture of the fiber  [--]
        n_outside : (optional) refractive index of medium outside fiber [--]

    Returns:
        maximum entrance/exit half-angle of the fiber [radians]
    """
    return np.arcsin(NA / n_outside)


def critical_angle(n_core, n_clad):
    """
    Calculate the angle (from the normal) for total internal reflection.

    Args:
        n_core : the index of refraction of the fiber core  [--]
        n_clad : the index of refraction of the fiber cladding  [--]

    Returns:
        angle of total internal reflection [radians]
    """
    return np.arcsin(n_clad / n_core)


def cutoff_wavelength(a, NA, ell=0, q=np.inf):
    """
    Calculate the cutoff wavelength for an optical fiber.

    By default this is the wavelength above which a step-index fiber carries
    only the fundamental mode.  The fundamental mode LP_01 itself has no
    cutoff; what the default returns is the cutoff of LP_11, found from the
    first zero of J_0 at V=2.405, and below that wavelength LP_11 also
    propagates.

    More generally, ell selects the Bessel function whose first zero sets the
    cutoff, so ell gives the cutoff wavelength of mode LP_(ell+1),1.  Using
    ell=1 and the first zero of J_1 at V=3.832 gives the cutoff of LP_21.

    If the cutoff wavelength for a graded index fiber is desired, then specify
    a different value for q.  The cutoff V is then larger by sqrt(1+2/q), so a
    parabolic fiber (q=2) has a cutoff wavelength sqrt(2) times shorter than a
    step-index fiber of the same radius and numerical aperture.

    Args:
        a :   radius of the fiber                               [m]
        NA :  numerical aperture of the fiber                   [-]
        ell : (optional) mode number, truncated to an integer   [-]
        q :   (optional) parameter for graded index fiber       [-]

    Returns:
        shortest wavelength for operation in the specified mode [m]
    """
    (Vc,) = jn_zeros(int(ell), 1)
    if np.isfinite(q):  # graded index fiber
        Vc *= np.sqrt(1 + 2 / q)
    return 2 * np.pi * a * NA / Vc


def esi_Delta(Delta, q):
    """
    Calculate equivalent step index (esi) Delta for a graded-index fiber.

    Args:
        Delta :  relative refractive index         [-]
        q :      parameter for graded index fiber  [-]

    Returns:
        equivalent relative refractive index   [-]
    """
    return q * (2 + q) / (1 + q) ** 2 * Delta


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


def esi_V_parameter(V, q):
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

    Just a convenience function.

    Args:
        n_core : the index of refraction of the fiber core      [-]
        Delta : relative index of refraction                    [-]

    Returns:
        numerical aperture                                      [-]
    """
    return n_core * np.sqrt(2 * Delta)


def numerical_aperture_graded_index(n_core, n_clad, q, r_over_a):
    """
    Calculate the numerical aperture of a graded-index optical fiber.

    The numerical aperture varies across the face of a graded-index fiber.
    This gives the result at the fractional distance across the fiber core,
    falling from the full numerical aperture on axis to zero at r_over_a=1.

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
        n_clad :  the index of refraction of the fiber cladding   [-]

    Returns:
        the relative refractive index (Delta)                    [-]
    """
    return (n_core**2 - n_clad**2) / (2 * n_core**2)


def R_par(m, theta):
    """
    Calculate the Fresnel reflection for parallel polarized light.

    This is the fraction of reflected intensity (not field) for light with an
    electric field parallel to the plane of incidence.

    When m is less than 1 and theta exceeds the critical angle the square root
    turns imaginary, np.emath.sqrt follows it onto the complex branch, and the
    reflectance comes back as exactly 1 for total internal reflection.

    Args:
        m :     index of the far medium relative to the near one,
                which may be complex                  [-]
        theta : angle from normal to surface  [radians]

    Returns:
        reflected power                       [-]
    """
    m2 = m * m
    c = np.cos(theta)
    s = np.sin(theta)
    d = np.emath.sqrt(m2 - s * s)  # complex past the critical angle
    return abs((m2 * c - d) / (m2 * c + d)) ** 2


def R_per(m, theta):
    """
    Calculate the Fresnel reflection for perpendicular polarized light.

    This is the fraction of reflected intensity (not field) for light with an
    electric field perpendicular to the plane of incidence.

    When m is less than 1 and theta exceeds the critical angle the square root
    turns imaginary, np.emath.sqrt follows it onto the complex branch, and the
    reflectance comes back as exactly 1 for total internal reflection.

    Args:
        m :     index of the far medium relative to the near one,
                which may be complex                  [-]
        theta : angle from normal to surface  [radians]

    Returns:
        reflected power                       [-]
    """
    m2 = m * m
    c = np.cos(theta)
    s = np.sin(theta)
    d = np.emath.sqrt(m2 - s * s)  # complex past the critical angle
    return abs((c - d) / (c + d)) ** 2


def R_unpolarized(m, theta):
    """
    Calculate the Fresnel reflection for unpolarized incident light.

    This is the fraction of reflected intensity (not field) for unpolarized
    incident light.

    Args:
        m :     index of the far medium relative to the near one,
                which may be complex                  [-]
        theta : angle from normal to surface  [radians]

    Returns:
        reflected power                       [-]
    """
    return (R_par(m, theta) + R_per(m, theta)) / 2


def V_parameter(a, NA, lambda0):
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
