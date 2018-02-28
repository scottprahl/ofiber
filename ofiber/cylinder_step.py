# Useful routines for step-index cylindrical waveguides based on chapter 8 of
#
#    A. Ghatak, K. Thyagarajan, An Introduction to Fiber Optics,
#    Cambridge University Press, 1998
#
# Scott Prahl
# Feb 2018

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.special import jn
from scipy.special import jn_zeros
from scipy.special import kn

__all__ = [ 'LP_mode_value',
			'LP_mode_values',
			'Plot_LP_modes',
			'LP_core_irradiance',
			'LP_clad_irradiance',
			'LP_total_irradiance',
			'LP_radial_field',
			'LP_radial_irradiance',
			'Tranverse_misalignment_loss_db',
			'Angular_misalignment_loss_db',
			'Longitudinal_misalignment_loss_db',
			'Bending_loss_db',
			'MFD',
			'PetermannW',
			'PetermannW_Approx',
			'V_d2bV_by_V']

def _LHS_eqn_8_40(b, V, ell):
    """
    Returns the left hand side of equation 8.40 in Ghatak
    This also works for eqn 8.40 (but is multiplied by -1)
    """
    U = V * np.sqrt(1 - b)
    return U * jn(ell - 1, U) / jn(ell, U)


def _RHS_eqn_8_40(b, V, ell):
    """
    Returns the right hand side of equation 8.40 in Ghatak
    This also works for eqn 8.40 (but is multiplied by -1)
    """
    W = V * np.sqrt(b)
    return -W * kn(ell - 1, W) / kn(ell, W)


def _cyl_mode_eqn(b, *args):
    """
    Returns the difference of RHS and LHS of 8.40 in Ghatak
    This function is zero when a mode exists in a step index fiber
    """
    V = args[0]
    ell = args[1]
    g1 = _LHS_eqn_8_40(b, V, ell)
    g2 = _RHS_eqn_8_40(b, V, ell)
    return g1 - g2


def LP_mode_value(V, ell, em):
    """
    Returns b for the specified mode ell,em in a circular step-index fiber
    If there is no such mode, returns zero
    """

    if em <= 0:  # modes start with 1, e.g., LP_01
        return 0

    if V <= 0:    # V must be positive
        return 0

    abit = 1e-3

    # set up bounds for this mode
    jnz = jn_zeros(ell, em)
    lo = max(0, 1 - (jnz[em - 1] / V)**2) + abit

    if em == 1:
        hi = 1 - abit
    else:
        hi = 1 - (jnz[em - 2] / V)**2 - abit
        if hi < lo:   # no such mode
            return 0

    try:
        b = brentq(_cyl_mode_eqn, lo, hi, args=(V, ell))
    except ValueError:  # happens when both hi and lo values have same sign
        b = 0           # therefore no such mode exists

    return b


def LP_mode_values(V, ell):
    """
    Returns array of b for the mode ell in a circular step-index fiber
    If there is no such mode, returns zero
    b[0] will correspond to LP_ell,1
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
    Returns a plot showing possible modes for the specified mode ell
    and V-parameter V
    """
    abit = 1e-5
    pltmin = -2 * V
    pltmax = 2 * V

    b = np.linspace(abit, 1 - abit, 251)

    g1 = _LHS_eqn_8_40(b, V, ell)
    g2 = _RHS_eqn_8_40(b, V, ell)

    # remove points so vertical lines are not drawn
    np.place(g1, g1 < pltmin, np.nan)
    np.place(g2, g2 < pltmin, np.nan)

    plt.plot([0, 1], [0, 0], ':k')
    plt.plot(b, g1)
    plt.plot(b, g2)

    # plot and label all the crossings
    all_b = LP_mode_values(V, ell)
    for i in range(len(all_b)):
        bb = all_b[i]
        y = _LHS_eqn_8_40(bb, V, ell)
        plt.scatter([bb], [y], s=30)
        plt.annotate('   LP%d%d' % (ell, i + 1), xy=(bb, y),va='top')

    plt.title('Modes for $\ell$=%d when V=%.3f' % (ell, V))
    plt.xlabel('b')
    plt.ylim(pltmin, pltmax)
    plt.xlim(0, 1)

    return plt


def LP_core_irradiance(V, b, ell):
    """
    core irradiance = core_power/fiber_core_area
    """
    U = V * np.sqrt(1 - b)
    return 1 - jn(ell + 1, U) * jn(ell - 1, U) / jn(ell, U)**2


def LP_clad_irradiance(V, b, ell):
    """
    clad irradiance = clad_power/fiber_core_area
    """
    W = V * np.sqrt(b)
    return kn(ell + 1, W) * kn(ell - 1, W) / kn(ell, W)**2 - 1


def LP_total_irradiance(V, b, ell):
    """
    total irradiance = total_power/fiber_core_area
    """
    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)
    val = V**2 / U**2 * kn(ell + 1, W)
    val *= kn(ell - 1, W) / kn(ell, W)**2
    return val 


def LP_radial_field(V, b, ell, r_over_a):
    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)

    values = np.empty_like(r_over_a)
    for i in range(len(r_over_a)):
        r = abs(r_over_a[i])
        if r < 1:
            values[i] = jn(ell, U * r) / jn(ell, U)
        else:
            values[i] = kn(ell, W * r) / kn(ell, W)

    return values / np.sqrt(LP_total_irradiance(V, b, ell))


def LP_radial_irradiance(V, b, ell, r_over_a):

    field = LP_radial_field(V, b, ell, r_over_a)
    return field**2


def Tranverse_misalignment_loss_db(w1, w2, u):
    sq = w1**2 + w2**2
    loss = (2 * w1 * w2 / sq)**2 * np.exp(-2 * u**2 / sq)
    return -10 * np.log10(loss)


def Angular_misalignment_loss_db(n, w, theta, lambda0):
    return 4.34 * (np.pi * w * theta * n / lambda0)**2


def Longitudinal_misalignment_loss_db(n, w, D, lambda0):
    dhat = D * lambda0 / (2 * np.pi * n * w**2)
    return 10 * np.log10(1 + dhat**2)


def _Bending_loss_db_scalar(n1, Delta, a, Rc, lambda0):
    """
    returns the bending loss in dB/m
    using eqn below eqn 10.29 in Ghatak 
    a is core radius in [m]
    n1 is core index
    Delta is core-cladding difference
    Rc is radius of curvature in [m]
    lambda0 is wavelength in vacuum in [m]
    """
    k0 = 2*np.pi/lambda0
    V = k0 * a * n1 * np.sqrt(2*Delta)
    b = LP_mode_value(V, 0, 1)
    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)
    if W == 0:
    	return np.nan
    val = 4.343*np.sqrt(np.pi/4/a/Rc)
    val *= (U/V/kn(1, W))**2
    val *= W**-1.5
    val *= np.exp(-2*W**3*Rc/3/k0**2/a**3/n1**2)
    return val


def Bending_loss_db(n1, Delta, a, Rc, lambda0):
    """
    returns the bending loss in dB/m
    using eqn below eqn 10.29 in Ghatak 
    a is core radius in [m]
    n1 is core index
    Delta is core-cladding difference
    Rc is radius of curvature in [m]
    lambda0 is wavelength in vacuum in [m]
    """
    if np.isscalar(a):
        alpha = _Bending_loss_db_scalar(n1, Delta, a, Rc, lambda0)
    else:
        alpha = np.empty_like(a)
        for i in range(len(alpha)):
            alpha[i] = _Bending_loss_db_scalar(n1, Delta, a[i], Rc, lambda0)
    return alpha


def MFD(V):
    return 0.65 + 1.619 * V**-1.5 + 2.879 * V**-6


def _PetermannW_scalar(V):
    b = LP_mode_value(V, 0, 1)
    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)
    denom =  W * jn(0, U)
    if denom == 0 :
        return np.nan
    return np.sqrt(2) * jn(1, U) / denom


def PetermannW(V):
    """
    returns the Petermann-2 radius (divided by core radius) for a step index fiber
    using eqn 8.86 in Ghatak
    result has no units
    """
    if np.isscalar(V):
        wp = _PetermannW_scalar(V)
    else:
        wp = np.empty_like(V)
        for i in range(len(wp)):
            wp[i] = _PetermannW_scalar(V[i])
    return wp


def PetermannW_Approx(V):
    """
    returns approximate Petermann-2 radius for a step index fiber
    using eqn below eqn 8.86 in Ghatak (good for 1.5<V<2.5)
    result has no units because it is the spot size divided by core radius
    """
    return MFD(V) - 0.016 - 1.567 * V**-7


def _V_d2bV_by_V_scalar(V, ell):
    """
    returns V*d^2(bV)/dV^2 for mode ell of a step index fiber
    using eqn 10.14 in Ghatak
    result has no units
    """
    b = LP_mode_value(V, ell, 1)
    if b == 0:
        return 0

    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)

    kappa_ell = kn(ell, W)**2 / kn(ell - 1, W) 
    kappa_ell /= kn(ell + 1, W)
    sum = 3 * W**2 - 2 * kappa_ell * (W**2 - U**2)
    val = W * (W**2 + U**2 * kappa_ell) * (kappa_ell - 1)
    val *= (kn(ell - 1, W) + kn(ell + 1, W))
    val /= kn(ell, W)
    sum += val
    return 2 * U**2 * kappa_ell / V**2 / W**2 * sum


def V_d2bV_by_V(V, ell):
    """
    returns V*d^2(bV)/dV^2 for mode ell of a step index fiber
    using eqn 10.14 in Ghatak
    result has no units
    """
    if np.isscalar(V):
        return _V_d2bV_by_V_scalar(V,ell)
    else:
        v_by_v = np.empty_like(V)
        for i in range(len(v_by_v)):
            v_by_v[i] = _V_d2bV_by_V_scalar(V[i],ell)

    return v_by_v

