# Useful routines for step-index planar waveguides based on chapter 7 of
#
#    A. Ghatak, K. Thyagarajan, An Introduction to Fiber Optics,
#    Cambridge University Press, 1998
#
# Scott Prahl
# Feb 2018

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

__all__ = [ 'TE_mode_plot',
			'TE_crossing',
			'TE_crossings',
			'TM_mode_plot',
			'TM_crossing',
			'TM_crossings',
			'TE_field',
			'TE_propagation_constant',
			'TM_propagation_constant']

def _base_mode_plot(V):
    """ Create a plot with symmetric and asymmetric modes
        This is used for both TE and TM plots
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
    V = args[0]
    mode = args[1]
    if mode % 2 == 0:
        return xi * np.tan(xi) - np.sqrt((V / 2)**2 - xi * xi)
    else:
        return xi / np.tan(xi) + np.sqrt((V / 2)**2 - xi * xi)


def TE_crossing(V, mode):
    abit = 1e-5
    lo = abit + mode * np.pi / 2
    hi = min(np.pi / 2 - abit + mode * np.pi / 2, V / 2)
    if lo > V / 2:
        return 0   # mode does not exist

    try:
        b = brentq(_TE_mode, lo, hi, args=(V, mode))
    except ValueError:  # happens when both hi and lo values have same sign
        b = 0           # therefore no such mode exists

    return b


def TE_crossings(V):
    ncross = int(V / 2 / (np.pi / 2)) + 1
    crossings = np.empty(ncross)

    for mode in range(ncross):
        crossings[mode] = TE_crossing(V, mode)

    return crossings


def TM_mode_plot(V, n1, n2):
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
    V = args[0]
    n1 = args[1]
    n2 = args[2]
    mode = args[3]
    if mode % 2 == 0:
        return xi * np.tan(xi) - (n1/n2)**2 * np.sqrt((V/2)**2 - xi**2)
    else:
        return xi / np.tan(xi) + (n1/n2)**2 * np.sqrt((V/2)**2 - xi**2)


def TM_crossing(V, n1, n2, mode):
    abit = 1e-5
    lo = abit + mode * np.pi / 2
    hi = min(np.pi / 2 - abit + mode * np.pi / 2, V / 2)
    if lo > V / 2:
        return 0  # mode does not exist

    try:
        b = brentq(_TM_mode, lo, hi, args=(V, n1, n2, mode))
    except ValueError:  # happens when both hi and lo values have same sign
        b = 0           # therefore no such mode exists

    return b


def TM_crossings(V, n1, n2):
    ncross = int(V / 2 / (np.pi / 2)) + 1
    crossings = np.empty(ncross)

    for mode in range(ncross):
        crossings[mode] = TM_crossing(V, n1, n2, mode)

    return crossings


def TE_field(V, d, x, mode):
    xi = TE_crossing(V, mode)
    gdby2 = np.sqrt((V / 2)**2 - xi**2)
    gamma = 2 / d * gdby2
    kappa = 2 / d * xi
    Ey = np.empty(len(x))
    if mode % 2 == 0:
        A = 1
        C = A * np.cos(xi) / np.exp(-gdby2)
        for j in range(len(x)):
            if abs(x[j]) < d / 2:
                Ey[j] = A * np.cos(x[j] * kappa)
            else:
                Ey[j] = C * np.exp(-gamma * abs(x[j]))
    else:
        B = 1
        D = B * np.sin(xi) / np.exp(-gdby2)
        for j in range(len(x)):
            if abs(x[j]) < d / 2:
                Ey[j] = B * np.sin(kappa * x[j])
            else:
                Ey[j] = D * np.sign(x[j]) * np.exp(-gamma * abs(x[j]))

    return Ey


def TE_propagation_constant(V, mode):
    b = np.empty(len(V))
    for i in range(len(V)):
        xi = TE_crossing(V[i], mode)
        if xi == 0:
            b[i] = 0
        else:
            b[i] = 1 - (2 * xi / V[i])**2
    return b


def TM_propagation_constant(V, n1, n2, mode):
    b = np.empty(len(V))
    for i in range(len(V)):
        xi = TM_crossing(V[i], n1, n2, mode)
        if xi == 0 :
            b[i] = 0
        else :
            b[i] = 1 - (2 * xi / V[i])**2
    return b


