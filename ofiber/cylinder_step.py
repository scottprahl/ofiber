# pylint: disable=invalid-name
# pylint: disable=no-name-in-module

"""
Useful routines for step-index cylindrical waveguides.

Based on chapter 8 of A. Ghatak, K. Thyagarajan, An Introduction to Fiber
Optics, Cambridge University Press, 1998

Todo:
    * lowercase functions
    * rename Plot_LP_modes
    * sort out normalization for fields and irradiances
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.special import jn
from scipy.special import jn_zeros
from scipy.special import kn

__all__ = ('LP_mode_value',
           'LP_mode_values',
           'Plot_LP_modes',
           'LP_core_irradiance',
           'LP_clad_irradiance',
           'LP_total_irradiance',
           'LP_radial_field',
           'LP_radial_irradiance',
           'Gaussian_envelope_Omega',
           'Gaussian_radial_irradiance',
           'Tranverse_misalignment_loss_db',
           'Angular_misalignment_loss_db',
           'Longitudinal_misalignment_loss_db',
           'Bending_loss_db',
           'MFR',
           'MFD',
           'PetermannW',
           'PetermannW_Approx',
           'V_d2bV_by_V',
           'V_d2bV_by_V_Approx')


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
    return U * jn(ell - 1, U) / jn(ell, U)


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
    return -W * kn(ell - 1, W) / kn(ell, W)


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


def LP_mode_value(V, ell, em):
    """
    Calculate guided b for mode (ell,em) in a circular step-index fiber.

    b is the normalized propagation constant.  Each guided mode in an optical
    fiber has a specific value of b that depends on the fiber parameter V
    and the mode number.

    If no mode exists, a value of 0 is returned

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
        return 0    # modes start with 1, e.g., LP_01

    if V <= 0:
        return 0    # V must be positive

    abit = 1e-3

    # set up bounds for this mode
    jnz = jn_zeros(ell, em)
    lo = max(0, 1 - (jnz[em - 1] / V)**2) + abit

    if em == 1:
        hi = 1 - abit
    else:
        hi = 1 - (jnz[em - 2] / V)**2 - abit

    if hi < lo:
        return 0  # no such mode

    try:
        b = brentq(_cyl_mode_eqn, lo, hi, args=(V, ell))
    except ValueError:  # happens when both hi and lo values have same sign
        return 0        # therefore no such mode exists

    return b


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
        if b == 0:
            break
        all_b = np.append(all_b, b)

    return all_b


def Plot_LP_modes(V, ell):
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
    return 1 - jn(ell + 1, U) * jn(ell - 1, U) / jn(ell, U)**2


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
    return kn(ell + 1, W) * kn(ell - 1, W) / kn(ell, W)**2 - 1


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
    val = V**2 / U**2 * kn(ell + 1, W)
    val *= kn(ell - 1, W) / kn(ell, W)**2
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

    A = jn(ell, U * r) / jn(ell, U)
    B = kn(ell, W * r) / kn(ell, W)
    values = np.where(r < 1, A, B)
    return values / np.sqrt(LP_total_irradiance(V, b, ell))


def LP_radial_irradiance(V, b, ell, r_over_a):
    """
    Calculate the normalized irradiance in a step-index fiber.

    The normalization is done such that
           2*np.trapz(LP(r_over_a)*r_over_a, r_over_a) =1
    or
          integral_over_space/(area of core) = 1
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


def Gaussian_envelope_Omega(V):
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
    Omega_over_a = jn(0, U) * V/U * kn(1, W)/kn(0, W)
    return Omega_over_a


def Gaussian_radial_irradiance(V, r_over_a):
    """
    Calculate the normalized irradiance in a step-index fiber.

    The normalization is done assuming
    the Gaussian envelope approximation for the LP_01 mode. The result
    is normalized such that
           np.trapz(Gaussian(r_over_a)*r_over_a, r_over_a) =1/2

    Args:
        V:        V-parameter for fiber            [-]
        r_over_a: (radial position)/(core radius)  [-]
    Returns:
        normalized irradiance at points r_over_a   [-]
    """
    Omega_over_a = Gaussian_envelope_Omega(V)
    return 1/Omega_over_a**2 * np.exp(-r_over_a**2/Omega_over_a**2)


def Tranverse_misalignment_loss_db(w1, w2, u):
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


def Angular_misalignment_loss_db(n, w, theta, lambda0):
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


def Longitudinal_misalignment_loss_db(n1, w, D, lambda0):
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


def _Bending_loss_db_scalar(n1, Delta, a, Rc, lambda0):
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
    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)
    if W == 0:
        return np.nan
    val = 4.343 * np.sqrt(np.pi / 4 / a / Rc)
    val *= (U / V / kn(1, W))**2
    val *= W**-1.5
    val *= np.exp(-2 * W**3 * Rc / 3 / k0**2 / a**3 / n1**2)
    return val


def Bending_loss_db(n1, Delta, a, Rc, lambda0):
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
        alpha = _Bending_loss_db_scalar(n1, Delta, a, Rc, lambda0)
    else:
        alpha = np.empty_like(a)
        for i, aa in enumerate(a):
            alpha[i] = _Bending_loss_db_scalar(n1, Delta, aa, Rc, lambda0)
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
    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)
    denom = W * jn(0, U)
    if denom == 0:
        return np.nan
    return np.sqrt(2) * jn(1, U) / denom


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
    if b == 0:
        return 0

    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)

    kappa_ell = kn(ell, W)**2 / kn(ell - 1, W)
    kappa_ell /= kn(ell + 1, W)
    summ = 3 * W**2 - 2 * kappa_ell * (W**2 - U**2)
    val = W * (W**2 + U**2 * kappa_ell) * (kappa_ell - 1)
    val *= (kn(ell - 1, W) + kn(ell + 1, W))
    val /= kn(ell, W)
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
