# pylint: disable=invalid-name
# pylint: disable=no-name-in-module
# pylint: disable=consider-using-f-string
# pylint: disable=too-many-arguments
# pylint: disable=no-member

"""
Useful routines for step-index cylindrical waveguides.

See <https://ofiber.readthedocs.io> for usage examples.

Based on chapter 8 of A. Ghatak, K. Thyagarajan, An Introduction to Fiber
Optics, Cambridge University Press, 1998

Functions to calculate and plot modes for step index fibers.  Specifically::

    LP_mode_value(V, ell, em)
    LP_mode_values(V, ell)
    LP_core_irradiance(V, b, ell)
    LP_clad_irradiance(V, b, ell)
    LP_total_irradiance(V, b, ell)
    LP_radial_field(V, b, ell, r_over_a)
    LP_radial_irradiance(V, b, ell, r_over_a)
    gaussian_envelope_Omega(V)
    gaussian_radial_irradiance(V, r_over_a)
    plot_LP_modes(V, ell)

Functions to estimate losses::

    angular_misalignment_loss_db(n, w, theta, lambda0)
    bending_loss_db(n1, Delta, a, Rc, lambda0)
    longitudinal_misalignment_loss_db(n1, w, D, lambda0)
    transverse_misalignment_loss_db(w1, w2, u)

Functions to find equivalent core diameters::

    MFR(V)
    MFD(V)
    PetermannW(V)
    PetermannW_Approx(V)

And finally, a couple of routines to help with waveguide dispersion
calculations::

    V_d2bV_by_V(V, ell)
    V_d2bV_by_V_Approx(V, ell)
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from scipy import special

_all_ = ('LP_mode_value',
         'LP_mode_values',
         'plot_LP_modes',
         'LP_core_irradiance',
         'LP_clad_irradiance',
         'LP_total_irradiance',
         'LP_radial_field',
         'LP_radial_irradiance',
         'gaussian_envelope_Omega',
         'gaussian_radial_irradiance',
         'transverse_misalignment_loss_db',
         'angular_misalignment_loss_db',
         'longitudinal_misalignment_loss_db',
         'bending_loss_db',
         'MFR',
         'MFD',
         'PetermannW',
         'PetermannW_Approx',
         'V_d2bV_by_V',
         'V_d2bV_by_V_Approx',
         'FF_polar_irradiance_x',
         'FF_irradiance_x',
          )


def _LHS_eqn_8_40(b, V, ell):
    """
    Calculate the left hand side of the eigenvalue eqn 8.40 in Ghatak.

    Also works for ell=0 (but is multiplied by -1 relative to eqn 8.41).
    This is private method that should not be needed outside this module.

    Args:
        b:      normalized propagation constant  [-]
        V:      V-parameter for fiber            [-]
        ell:    desired fiber mode               [-]
    Returns:
        LHS of equation 8.40                     [-]
    """
    U = V * np.sqrt(1 - b)
    return U * special.jv(ell - 1, U) / special.jv(ell, U)


def _RHS_eqn_8_40(b, V, ell):
    """
    Calculate the right hand side of the eigenvalue eqn 8.40 in Ghatak.

    Also works for ell=0 (but is multiplied by -1 relative to eqn 8.41).
    This is private method that should not be needed outside this module.

    Args:
        b:      normalized propagation constant  [-]
        V:      V-parameter for fiber            [-]
        ell:    desired fiber mode               [-]
    Returns:
        RHS of equation 8.40                     [-]
    """
    W = V * np.sqrt(b)
    return -W * special.kn(ell - 1, W) / special.kn(ell, W)


def _cyl_mode_eqn(b, *args):
    """
    Return the difference of RHS and LHS of 8.40 in Ghatak.

    This function is zero when a guided mode exists in the step index fiber.
    This is a private function and should not be needed outside this module.

    Args:
        b:      normalized propagation constant  [-]
        arg[0]: V-parameter for optical fiber    [-]
        arg[1]: desired fiber mode               [-]
    Returns:
        LHS-RHS of equation 8.40                 [-]
    """
    V = args[0]
    ell = args[1]
    g1 = _LHS_eqn_8_40(b, V, ell)
    g2 = _RHS_eqn_8_40(b, V, ell)
    return g1 - g2


def _LP_mode_value(V, ell, em):
    """
    Calculate guided b for mode (ell,em) in a circular step-index fiber.

    b is the normalized propagation constant.  Each guided mode in an optical
    fiber has a specific value of b that depends on the fiber parameter V
    and the mode number.

    If no mode exists, a value of None is returned

    The LP_lm is specified by the (ell,em) to avoid confusion between the
    number 1 and the letter l.

    For cylindrical fibers, em is a positive integer: thus there are modes
    LP_01, LP_02, but not LP_10.

    Args:
        V:   V-parameter for optical fiber    [-]
        ell: primary fiber mode   (integer)   [-]
        em:  secondary fiber mode (integer>0) [-]
    Returns:
        guided normalized propagation constant for mode (ell,em)  [-]
    """
    if ell < 0:
        ell *= -1   # negative ells are same as positive ones

    if em <= 0:
        return None    # modes start with 1, e.g., LP_01

    if V <= 0:
        return None    # V must be positive

    abit = 1e-5

    # set up bounds for this mode
    jnz = special.jn_zeros(ell, em)
    lo = max(0, 1 - (jnz[em - 1] / V)**2) + abit

    if em == 1:
        hi = 1 - abit
    else:
        hi = 1 - (jnz[em - 2] / V)**2 - abit

    if hi < lo:
        return None  # no such mode

    try:
        b = scipy.optimize.brentq(_cyl_mode_eqn, lo, hi, args=(V, ell))
    except ValueError:  # happens when both hi and lo values have same sign
        return None     # therefore no such mode exists

    return b


def LP_mode_value(V, ell, em):
    """
    Calculate guided b for mode (ell,em) in a circular step-index fiber.

    b is the normalized propagation constant.  Each guided mode in an optical
    fiber has a specific value of b that depends on the fiber parameter V
    and the mode numbers ell and em.
    
    This is a wrapper function that handles V, ell, or em being possible
    arrays.  See `LP_mode_value` for details.

    Args:
        V:   V-parameter for optical fiber    [-]
        ell: primary fiber mode   (integer)   [-]
        em:  secondary fiber mode (integer>0) [-]
    Returns:
        guided normalized propagation constant for mode (ell,em)  [-]
    """
    if em < 1:
        raise ValueError("The mode number 'em' must be one or greater.")

    if np.isscalar(V) + np.isscalar(ell) + np.isscalar(em) < 2:
        raise ValueError("only one of V, ell, and em can be an array")

    if not np.isscalar(V):
        b = np.zeros_like(V)
        for i, VV in enumerate(V):
            b[i] = LP_mode_value(VV, ell, em)
        return b

    if not np.isscalar(ell):
        b = np.zeros_like(ell)
        for i, ll in enumerate(ell):
            b[i] = LP_mode_value(V, ll, em)
        return b

    if not np.isscalar(em):
        b = np.zeros_like(em)
        for i, mm in enumerate(em):
            b[i] = LP_mode_value(V, ell, mm)
        return b

    return _LP_mode_value(V, ell, em)


def LP_mode_values(V, ell):
    """
    Calculate all guided b for mode ell in a circular step-index fiber.

    If there is no such mode, returns an empty array

    Note that in the returned array b[0] will correspond to LP_ell,1

    Args:
        V:   V-parameter for optical fiber    [-]
        ell: primary fiber mode   (integer)   [-]
    Returns:
        array of normalized propagation constant for mode ell  [-]
    """
    all_b = np.array([])
    for em in range(1, 10):
        b = LP_mode_value(V, ell, em)
        if b is None:
            break
        all_b = np.append(all_b, b)

    return all_b


def plot_LP_modes(V, ell):
    """
    Produce a plot show possible eigenvalue solutions for step index fiber.

    The solutions correspond to places where the curves cross one another.  No
    crossing means that there is no guided mode for that mode value.

    Args:
        V:   V-parameter for optical fiber    [-]
        ell: primary fiber mode   (integer)   [-]
    Returns:
        graph for mode ell   [matplotlib.pyplot object]
    """
    abit = 1e-5
    pltmin = -2 * V
    pltmax = 2 * V

    b = np.linspace(abit, 1 - abit, 251)

    g1 = _LHS_eqn_8_40(b, V, ell)
    g2 = _RHS_eqn_8_40(b, V, ell)

    # remove points so confusing vertical retrace lines are not shown
    np.place(g1, g1 < pltmin, np.nan)
    np.place(g2, g2 < pltmin, np.nan)

    plt.plot([0, 1], [0, 0], ':k')
    plt.plot(b, g1)
    plt.plot(b, g2)

    # plot and label all the crossings
    all_b = LP_mode_values(V, ell)
    for i, bb in enumerate(all_b):
        y = _LHS_eqn_8_40(bb, V, ell)
        plt.scatter([bb], [y], s=30)
        plt.annotate(r'   LP$_{%d%d}$' % (ell, i + 1), xy=(bb, y), va='top')

    plt.title(r'Modes for $\ell$=%d when V=%.3f' % (ell, V))
    plt.xlabel('b')
    plt.ylim(pltmin, pltmax)
    plt.xlim(0, 1)

    return plt


def LP_core_irradiance(V, b, ell):
    """
    Calculate the core irradiance for a step-index fiber.

    See Ghatak equation 8.56.  The returned value is the total
    core power divided by the area of the core.

    Args:
        V:      V-parameter for fiber            [-]
        b:      normalized propagation constant  [-]
        ell:    desired fiber mode               [-]
    Returns:
        total core power over core area          [-]
    """
    U = V * np.sqrt(1 - b)
    return 1 - special.jv(ell + 1, U) * special.jv(ell - 1, U) / special.jv(ell, U)**2


def LP_clad_irradiance(V, b, ell):
    """
    Calculate the cladding irradiance for a step-index fiber.

    See Ghatak equation 8.57.  The returned value is the total
    cladding power divided by the area of the core.

    Args:
        V:      V-parameter for fiber            [-]
        b:      normalized propagation constant  [-]
        ell:    desired fiber mode               [-]
    Returns:
        total cladding power over core area      [-]
    """
    W = V * np.sqrt(b)
    return special.kn(ell + 1, W) * special.kn(ell - 1, W) / special.kn(ell, W)**2 - 1


def LP_total_irradiance(V, b, ell):
    """
    Calculate the total irradiance for a step-index fiber.

    See Ghatak equation 8.58.  The returned value is the total
    power (cladding + core) divided by the area of the core.

    Args:
        V:      V-parameter for fiber            [-]
        b:      normalized propagation constant  [-]
        ell:    desired fiber mode               [-]
    Returns:
        total power over core area               [-]
    """
    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)
    val = V**2 / U**2 * special.kn(ell + 1, W)
    val *= special.kn(ell - 1, W) / special.kn(ell, W)**2
    return val


def LP_radial_field(V, b, ell, r_over_a):
    """
    Calculate the normalized field in a step-index fiber.

    Args:
        V:        V-parameter for fiber            [-]
        b:        normalized propagation constant  [-]
        ell:      desired fiber mode               [-]
        r_over_a: (radial position)/(core radius)  [-]
    Returns:
        normalized field at point r_over_a         [-]
    """
    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)
    r = abs(r_over_a)  # same value for negative radii

    A = special.jv(ell, U * r) / special.jv(ell, U)
    B = special.kn(ell, W * r) / special.kn(ell, W)
    values = np.where(r < 1, A, B)
    return values / np.sqrt(LP_total_irradiance(V, b, ell))


def LP_radial_irradiance(V, b, ell, r_over_a):
    """
    Calculate the normalized irradiance in a step-index fiber.

    The normalization is done such that
    integral_over_space/(area of core) = 1
    or
    2*np.trapz(LP(r_over_a)*r_over_a, r_over_a) =1

    Args:
        V:        V-parameter for fiber            [-]
        b:        normalized propagation constant  [-]
        ell:      desired fiber mode               [-]
        r_over_a: (radial position)/(core radius)  [-]
    Returns:
        normalized irradiance at points r_over_a   [-]
    """
    field = LP_radial_field(V, b, ell, r_over_a)
    return field**2


def gaussian_envelope_Omega(V):
    """
    Calculate the normalized irradiance in a step-index fiber.

    The normalization is done assuming
    the Gaussian envelope approximation for the LP_01 mode.

    Args:
        V:        V-parameter for fiber            [-]
    Returns:
        Omega_over_core_radius                     [-]
    """
    b = LP_mode_value(V, 0, 1)
    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)
    Omega_over_a = special.jv(0, U) * V / U * special.kn(1, W) / special.kn(0, W)
    return Omega_over_a


def gaussian_radial_irradiance(V, r_over_a):
    """
    Calculate the normalized irradiance in a step-index fiber.

    The normalization is done assuming
    the Gaussian envelope approximation for the LP_01 mode. The result
    is normalized such that
    np.trapz(Gaussian(r_over_a)*r_over_a, r_over_a) = 1/2

    Args:
        V:        V-parameter for fiber            [-]
        r_over_a: (radial position)/(core radius)  [-]
    Returns:
        normalized irradiance at points r_over_a   [-]
    """
    Omega_over_a = gaussian_envelope_Omega(V)
    return 1 / Omega_over_a**2 * np.exp(-r_over_a**2 / Omega_over_a**2)


def transverse_misalignment_loss_db(w1, w2, u):
    """
    Calculate the loss due to transverse fiber misalignment.

    See Ghatak eqn 8.69

    Args:
        w1:      mode field radius of first fiber  [m]
        w2:      mode field radius of second fiber [m]
        u:       transverse misalignment           [m]
    Returns:
        transverse misalignment loss in dB         [-]
    """
    sq = w1**2 + w2**2
    loss = (2 * w1 * w2 / sq)**2 * np.exp(-2 * u**2 / sq)
    return -10 * np.log10(loss)


def angular_misalignment_loss_db(n, w, theta, lambda0):
    """
    Calculate the loss due to angular fiber misalignment.

    See Ghatak eqn 8.75

    Args:
        n:        index between fiber ends [-]
        w:        mode field radius        [m]
        theta:    angular misalignment     [radians]
        lambda0:  wavelength in vacuum     [m]
    Returns:
        angular misalignment loss in dB    [-]
    """
    return 4.34 * (np.pi * w * theta * n / lambda0)**2


def longitudinal_misalignment_loss_db(n1, w, D, lambda0):
    """
    Calculate the loss due to longitudinal fiber misalignment.

    See Ghatak eqn 8.81

    Args:
        n:        index between fiber ends      [-]
        w:        mode field radius             [m]
        D:        longitudinal fiber separation [m]
        lambda0:  wavelength in vacuum          [m]
    Returns:
        longitudinal misalignment loss dB       [-]
    """
    dhat = D * lambda0 / (2 * np.pi * n1 * w**2)
    return 10 * np.log10(1 + dhat**2)


def _bending_loss_db_scalar(n1, Delta, a, Rc, lambda0):
    """
    Calculate the bending loss in dB/m.

    The bending loss is given by eqn 10.29 in Ghatak.  This private method
    only works for scalar values.

    Args:
        a:        core radius                 [m]
        n1:       core index                  [-]
        Delta:    refractive index difference [-]
        Rc:       radius of curvature in      [m]
        lambda0:  wavelength in vacuum in     [m]
    Returns:
        bending loss in dB/m                  [1/m]
    """
    k0 = 2 * np.pi / lambda0
    V = k0 * a * n1 * np.sqrt(2 * Delta)
    b = LP_mode_value(V, 0, 1)
    if b is None:
        return np.nan
    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)
    val = 4.343 * np.sqrt(np.pi / 4 / a / Rc)
    val *= (U / V / special.kn(1, W))**2
    val *= W**-1.5
    val *= np.exp(-2 * W**3 * Rc / 3 / k0**2 / a**3 / n1**2)
    return val


def bending_loss_db(n1, Delta, a, Rc, lambda0):
    """
    Calculate the bending loss in dB/m.

    This is a convenience method that works when a is an array.

    Args:
        a:        core radius                 [m]
        n1:       core index                  [-]
        Delta:    refractive index difference [-]
        Rc:       radius of curvature in      [m]
        lambda0:  wavelength in vacuum in     [m]
    Returns:
        bending loss in dB/m                  [1/m]
    """
    if np.isscalar(a):
        alpha = _bending_loss_db_scalar(n1, Delta, a, Rc, lambda0)
    else:
        alpha = np.empty_like(a)
        for i, aa in enumerate(a):
            alpha[i] = _bending_loss_db_scalar(n1, Delta, aa, Rc, lambda0)
    return alpha


def MFR(V):
    """
    Approximate the mode field radius for a step-index single mode fiber.

    The approximation is fairly accurate for V>1. In the multimode range
    (V > 2.405), it applies to the fundamental mode.

    D. Marcuse, "Loss analysis of single-mode fiber splices", Bell Syst.
    Tech. J., 56, 703 (1977)

    Args:
        V:      V-parameter of the fiber                            [--]
    Returns:
        approximate mode field radius normalized by the core radius [--]
    """
    return 0.65 + 1.619 * V**-1.5 + 2.879 * V**-6


def MFD(V):
    """
    Approximate the mode field diameter for a step-index single mode fiber.

    See MFR() for details.

    Args:
        V:      V-parameter of the fiber                              [--]
    Returns:
        approximate mode field diameter normalized by the core radius [--]
    """
    return 2 * MFR(V)


def _PetermannW_scalar(V):
    """
    Calculate the Petermann-2 radius for a step-index fiber.

    This private method only works when V is a scalar.

    Args:
        V:      V-parameter of the fiber                          [--]
    Returns:
        approximate Petermann-2 radius normalized by core radius  [--]
    """
    b = LP_mode_value(V, 0, 1)
    if b is None:
        return np.nan
    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)
    denom = W * special.jv(0, U)
    return np.sqrt(2) * special.jv(1, U) / denom


def PetermannW(V):
    """
    Calculate the Petermann-2 radius for a step-index fiber.

    This is a convenience function that works when V is an array.

    Args:
        V:      V-parameter of the fiber                          [--]
    Returns:
        approximate Petermann-2 radius normalized by core radius  [--]
    """
    if np.isscalar(V):
        wp = _PetermannW_scalar(V)
    else:
        wp = np.empty_like(V)
        for i, VV in enumerate(V):
            wp[i] = _PetermannW_scalar(VV)
    return wp


def PetermannW_Approx(V):
    """
    Approximate the Petermann-2 radius for a step-index fiber.

    The approximation is valid for single mode fibers (1.5<V<2.5).  The result
    is the ratio of the Petermann-2 radius to the core radius.

    C. D. Hussey and F. Martinez, “Approximate analytical forms for
       the propagation characteristics of single-mode optical fibres”,
       Electron. Lett. 21, 1103 (1985).

    Args:
        V:      V-parameter of the fiber                          [--]
    Returns:
        approximate Petermann-2 radius normalized by core radius  [--]
    """
    return MFR(V) - 0.016 - 1.567 * V**-7


def _V_d2bV_by_V_scalar(V, ell):
    """
    Calculate V*d^2(bV)/dV^2 for mode ell of a step-index fiber.

    This private function only works for scalar values of V and ell. It
    finds V*d^2(bV)/dV^2 for mode ell of a step-index fiber using eqn 10.14

    Args:
        V:      V-parameter of the fiber     [--]
    Returns:
        V*d^2(bV)/dV^2                       [--]
    """
    b = LP_mode_value(V, ell, 1)
    if b is None:
        return 0

    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)

    kappa_ell = special.kn(ell, W)**2 / special.kn(ell - 1, W)
    kappa_ell /= special.kn(ell + 1, W)
    summ = 3 * W**2 - 2 * kappa_ell * (W**2 - U**2)
    val = W * (W**2 + U**2 * kappa_ell) * (kappa_ell - 1)
    val *= (special.kn(ell - 1, W) + special.kn(ell + 1, W))
    val /= special.kn(ell, W)
    summ += val
    return 2 * U**2 * kappa_ell / V**2 / W**2 * summ


def V_d2bV_by_V(V, ell):
    """
    Calculate V*d^2(bV)/dV^2 for mode ell of a step-index fiber.

    This value is needed to determine the waveguide dispersion.  This
    routine is a convenience function that works when V is an array.

    Args:
        V:      V-parameter of the fiber     [--]
    Returns:
        V*d^2(bV)/dV^2                       [--]
    """
    if np.isscalar(V):
        return _V_d2bV_by_V_scalar(V, ell)

    v_by_v = np.empty_like(V)
    for i, VV in enumerate(V):
        v_by_v[i] = _V_d2bV_by_V_scalar(VV, ell)

    return v_by_v


def V_d2bV_by_V_Approx(V):
    """
    Approximate V*d^2(bV)/dV^2 for single mode fiber.

    This value is needed to determine the waveguide dispersion.  This
    approximation is for the fundamental mode in the fiber and is good
    to 1% when 1.4<V<2.4.  Approximation by Marcuse (1979)

    Args:
        V:      V-parameter of the fiber     [--]
    Returns:
        V*d^2(bV)/dV^2                       [--]
    """
    return 0.080 + 0.549 * (2.834 - V)**2


def _FF_polar_x(ell, ka, theta, V, b):
    """
    Polar component of the far-field irradiance polarized in the x-direction.
    
    This implements equation 10.13 from Chen, Foundations for Guided-Wave Optics,
    Wiley-Interscience, 2007.  Use FF_irradiance_x() below to get the far-field
    irradiance.
    
    The polar angle is measured from the optical axis of the fiber.
    
    The product ka is dimensionless and is the product k * a = 2𝜋a/λ

    Args:
        ell: azimuthal mode number.
        ka: wavenumber * fiber radius
        theta: polar angle to calculate the field (radians).
        V: normalized frequency parameter of the fiber.
        b: normalized propagation constant.
    Returns:
        The calculated azimuthal field distribution
    """
    kasin = ka * np.sin(theta)
    Vb = V * np.sqrt(1 - b)
    ell1 = ell + 1

    Flnumer = kasin * special.jv(ell1, kasin) - (special.jv(ell, kasin) \
              / special.jv(ell, Vb)) * Vb * special.jv(ell1, Vb)
    Fldenom = (Vb**2 - kasin**2) * (V**2 * b + kasin**2)

    return Flnumer / Fldenom


def FF_irradiance_x(r, theta, phi, ell, lambda0, a, V, b):
    """
    Normalized far-field irradiance polarized in the x-direction.

    This calculation is based eqn 10.12 for the far-field from Chen, 
    Foundations for Guided-Wave Optics, Wiley-Interscience, 2007.  

    The magnitude of the field is squared and normalized by the square
    of the electric field magnitude.

    Args:
        r: radial distance from the fiber axis (microns).
        theta: polar angle in radians.
        phi: azimuthal angle in radians.
        ell: azimuthal mode number.
        lambda0: wavelength of light in the medium (microns).
        a: Radius of the fiber core (microns).
        V: normalized frequency parameter of the fiber.
        b: normalized propagation constant.
    Returns:
        The irradiance pattern as a function of r, theta, and phi.
    """
    k = 2 * np.pi / lambda0
    FF_ell = _FF_polar_x(ell, k*a, theta, V, b)
    FF = FF_ell * (k * a * V)**2 * np.cos(ell * phi) / (k * r)
    return FF**2


def FF_polar_irradiance_x(r, theta, ell, lambda0, a, V, b):
    """
    Integral of x-polarized far-field irradiance over all azimuthal angles.

    The magnitude of the field is squared and normalized by the square
    of the electric field magnitude. The integral of cos^2(ell*phi) over
    phi from 0 to 2𝜋 is 𝜋 and we get a slightly simpler result than
    the function above.

    Args:
        r: radial distance from the fiber axis (microns).
        theta: polar angle in radians.
        ell: azimuthal mode number.
        lambda0: Wavelength of light in the medium (microns).
        a: Radius of the fiber core (microns).
        V: Normalized frequency parameter of the fiber.
        b: Normalized propagation constant.
    Returns:
        Azimuthally integrated irradiance pattern.
    """
    k = 2 * np.pi / lambda0
    FF_ell = _FF_polar_x(ell, k*a, theta, V, b)
    return np.pi * (FF_ell * (k * a * V)**2 / (k * r))**2
