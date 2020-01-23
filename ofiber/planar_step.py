# pylint: disable=invalid-name

"""
Useful routines for step-index planar waveguides

Based on chapter 7 of A. Ghatak, K. Thyagarajan, An Introduction to 
Fiber Optics, Cambridge University Press, 1998
	
Todo:
    * Normalize field functions
    * Plot and label points at eigenvalue crossings
    * sort out language about eigenvalues
    * check TM solutions

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

__all__ = ('TE_crossing',
           'TE_crossings',
           'TE_field',
           'TE_mode_plot',
           'TE_propagation_constant',
           'TM_crossing',
           'TM_crossings',
           'TM_field',
           'TM_mode_plot',
           'TM_propagation_constant')


def _base_mode_plot(V):
    """
    Plot the symmetric and asymmetric modes in a planar waveguide.

    Args:
        V: the V-parameter for the waveguide

    Returns:
        a matplotlib.pyplot object
    """
    abit = 1e-5

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_aspect('equal')

    xi = np.linspace(abit, np.pi / 2 - abit, 50)
    while xi[0] < V / 2:
        plt.plot(xi, xi * np.tan(xi), color='red')
        xi += np.pi / 2

    xi = np.linspace(np.pi / 2, np.pi - abit, 50)
    while xi[0] < V / 2:
        plt.plot(xi, -xi / np.tan(xi), '--b')
        xi += np.pi

    plt.xlabel(r'$\xi=(d/2)\sqrt{k_0^2n_1^2-\beta^2}=(d/2)\kappa$')

    return plt


def TE_mode_plot(V):
    """
    Plot the symmetric and asymmetric TE modes in a planar waveguide.

    Args:
        V: the V-parameter for the waveguide

    Returns:
        a matplotlib.pyplot object
    """
    abit = 1e-5
    xi = np.linspace(0, V / 2 - abit, 100)
    circle = np.sqrt((V / 2)**2 - xi**2)

    aplt = _base_mode_plot(V)
    aplt.plot(xi, circle, ':k')
    ystr = r'$\xi\,\tan\xi$   or   $-\xi\,\cot\xi$   or   $\sqrt{V^2/4-\xi^2}$'
    aplt.ylabel(ystr)
    aplt.title('TE Modes in Planar Film Waveguide (V=%.2f)' % V)
    aplt.ylim(0, V / 2 + 1)
    aplt.xlim(0, V / 2 + 1)

    return aplt


def _TE_mode(xi, *args):
    """
    This is the eigenvalue equation for TE modes.  The zeros of this function
    can be used to determine the propagation factor beta for a particular mode.

    Args:
        xi    : non-dimensional propagation value in waveguide
        arg[0]: the V-parameter for the waveguide
        arg[1]: the desired mode (even values are symmetric)

    Returns:
        the distance from a solution
    """
    V = args[0]
    mode = args[1]
    if mode % 2 == 0:
        return xi * np.tan(xi) - np.sqrt(V**2 / 4 - xi * xi)
    else:
        return xi / np.tan(xi) + np.sqrt(V**2 / 4 - xi * xi)


def TE_crossing(V, mode):
    """
    Finds the xi value for a planar waveguide that solves the eigenvalue
    equation for the specified TE mode.

    Args:
        V    : the V-parameter for the waveguide
        mode : the mode

    Returns:
        the eigenvalue.  If no eigenvalue, 0 is returned.
    """
    abit = 1e-5
    lo = abit + mode * np.pi/2
    hi = min(np.pi/2 - abit + mode * np.pi/2, V/2)
    if lo > V / 2:
        return 0   # mode does not exist

    try:
        b = brentq(_TE_mode, lo, hi, args=(V, mode))
    except ValueError:  # happens when both hi and lo values have same sign
        return 0        # therefore no such mode exists

    return b


def TE_crossings(V):
    """
    Finds all xi values that for a planar waveguide that are solutions
    to the eigenvalue equation.

    Args:
        V    : the V-parameter for the waveguide

    Returns:
        an array eigenvalues.
    """
    ncross = int(V / np.pi) + 1
    crossings = np.empty(ncross)

    for mode in range(ncross):
        crossings[mode] = TE_crossing(V, mode)

    return crossings


def TM_mode_plot(V, n1, n2):
    """
    Plot the symmetric and asymmetric TM modes in a planar waveguide.

    Args:
        V:  the V-parameter for the waveguide
        n1: index of refraction inside the waveguide
        n2: index of refraction of the cladding

    Returns:
        a matplotlib.pyplot object
    """
    abit = 1e-5
    xi = np.linspace(0, V / 2 - abit, 100)
    ellipse = (n1 / n2)**2 * np.sqrt((V / 2)**2 - xi**2)

    aplt = _base_mode_plot(V)
    aplt.plot(xi, ellipse, ':k')
    ystr = r'$\xi\,\tan\xi$   or   $-\xi\,\cot\xi$   or'
    ystr = ystr + r'$(n_1/n_2)^2\sqrt{V^2/4-\xi^2}$'
    aplt.ylabel(ystr)
    aplt.title('TM Modes in Planar Film Waveguide (V=%.2f)' % V)
    ymax = (n1 / n2)**2 * V / 2
    aplt.ylim(0, ymax + 1)
    aplt.xlim(0, ymax + 1)

    return aplt


def _TM_mode(xi, *args):
    """
    This is the eigenvalue equation for TM modes.  The zeros of this function
    can be used to determine the propagation factor beta for a particular mode.

    Args:
        xi    : non-dimensional propagation value in waveguide
        arg[0]: the V-parameter for the waveguide
        arg[1]: index of the planar waveguide
        arg[2]: index of the cladding
        arg[3]: the desired mode (even values are symmetric)

    Returns:
        the distance from a solution
    """
    V = args[0]
    n1 = args[1]
    n2 = args[2]
    mode = args[3]
    if mode % 2 == 0:
        return xi * np.tan(xi) - (n1 / n2)**2 * np.sqrt(V**2 / 4 - xi**2)
    else:
        return xi / np.tan(xi) + (n1 / n2)**2 * np.sqrt(V**2 / 4 - xi**2)


def TM_crossing(V, n1, n2, mode):
    """
    Finds the xi value for a planar waveguide that solves the eigenvalue
    tequation for the specified TM mode.

    Args:
        V    : the V-parameter for the waveguide
        n1   : index of waveguide
        n2   : index of material next to waveguide
        mode : desired propagation mode

    Returns:
        the eigenvalue.  Returns 0 if no eigenvalue.
    """
    abit = 1e-5
    lo = abit + mode * np.pi / 2
    hi = min(np.pi / 2 - abit + mode * np.pi / 2, V / 2)
    if lo > V / 2:
        return 0  # mode does not exist

    try:
        b = brentq(_TM_mode, lo, hi, args=(V, n1, n2, mode))
    except ValueError:  # happens when both hi and lo values have same sign
        return 0        # therefore no such mode exists

    return b


def TM_crossings(V, n1, n2):
    """
    Finds all TM eigenvalues for a planar waveguide

    Args:
        V    : the V-parameter for the waveguide
        n1   : index of waveguide
        n2   : index of material next to waveguide

    Returns:
        an array of eigenvalues
    """
    ncross = int(V / 2 / (np.pi / 2)) + 1
    crossings = np.empty(ncross)

    for mode in range(ncross):
        crossings[mode] = TM_crossing(V, n1, n2, mode)

    return crossings


def _basic_field(V, d, x, mode, xi):
    """
    Calculates a field at location(s) x in a planar waveguide

    Args:
        V    : the V-parameter for the waveguide
        d    : thickness of the waveguide
        x    : desired field positions
        mode : specific mode
        xi   : kappa*d/2 for this mode

    Returns:
        the field at each position x
    """
    gdby2 = np.sqrt((V / 2)**2 - xi**2)   # gamma*d/2
    xgamma = 2 / d * gdby2 * abs(x)       # gamma*x
    kappa = 2 / d * xi

    if mode % 2 == 0:
        A = np.cos(kappa*x)
        B = np.cos(xi) * np.exp(gdby2-xgamma)
    else:
        A = np.sin(kappa*x)
        B = np.sin(xi) * np.sign(x) * np.exp(gdby2-xgamma)

    return np.where(abs(x) < d/2, A, B)


def TE_field(V, d, x, mode):
    """
    Calculates the TE field at location(s) x in a planar waveguide

    Args:
        V    : the V-parameter for the waveguide
        d    : thickness of the waveguide
        x    : desired field positions
        mode : specific mode

    Returns:
        the field at each position x
    """
    xi = TE_crossing(V, mode)
    E_y = _basic_field(V, d, x, mode, xi)
    return E_y


def TM_field(V, n1, n2, d, x, mode):
    """
    Calculates the TM field at location(s) x in a planar waveguide

    Args:
        V    : the V-parameter for the waveguide
        n1   : index of waveguide
        n2   : index of material next to waveguide
        d    : thickness of the waveguide
        x    : desired field positions
        mode : specific mode

    Returns:
        the field at each position x
    """
    xi = TM_crossing(V, n1, n2, mode)
    H_y = _basic_field(V, d, x, mode, xi)
    return H_y


def TE_propagation_constant(V, mode):
    """
    Calculates the dimensionless propagation constants for
    TE modes in a planar waveguide

    Args:
        V    : the V-parameter for the waveguide
        mode : specific mode

    Returns:
        an array of propagation constants for each value of V
    """
    b = np.empty_like(V)
    for i, VV in enumerate(V):
        xi = TE_crossing(VV, mode)
        if xi == 0:
            b[i] = 0
        else:
            b[i] = 1 - (2 * xi / VV)**2
    return b


def TM_propagation_constant(V, n1, n2, mode):
    """
    Calculates the dimensionless propagation constants for
    TM modes in a planar waveguide

    Args:
        V    : the V-parameter for the waveguide
        mode : specific mode

    Returns:
        an array of propagation constants for each value of V
    """
    b = np.empty_like(V)
    for i, VV in enumerate(V):
        xi = TM_crossing(VV, n1, n2, mode)
        if xi == 0:
            b[i] = 0
        else:
            b[i] = 1 - (2 * xi / VV)**2
    return b
