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

# a constant that to avoid infinity in tan(x)
abit = 1e-5

# make everything larger
plt.rcParams.update({'font.size': 14})

def base_mode_plot(V):
    """ Create a plot with symmetric and asymmetric modes
        This is used for both TE and TM plots
    """

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

    xi = np.linspace(0, V/2 - abit, 100)
    circle = np.sqrt((V/2)**2 - xi**2)

    aplt = base_mode_plot(V)
    aplt.plot(xi, circle, ':k')
    aplt.ylabel(r'$\xi\,\tan\xi$   or   $-\xi\,\cot\xi$   or   $\sqrt{V^2/4-\xi^2}$')
    aplt.title('TE Modes in Planar Film Waveguide (V=%.2f)' % V)
    aplt.ylim(0, V / 2 + 1)
    aplt.xlim(0, V / 2 + 1)

    return aplt


def TE_mode(xi, *args):
    V = args[0]
    mode = args[1]
    if mode % 2 == 0:
        return xi * np.tan(xi) - np.sqrt((V / 2)**2 - xi * xi)
    else:
        return xi / np.tan(xi) + np.sqrt((V / 2)**2 - xi * xi)


def TE_crossing(V, i):
    lower = abit + i * np.pi / 2
    upper = min(np.pi / 2 - abit + i * np.pi / 2, V / 2)
    if lower > V / 2:
        return 0   # mode does not exist

    return brentq(TE_mode, lower, upper, args=(V, i))


def TE_crossings(V):
    ncross = int(V / 2 / (np.pi / 2)) + 1
    crossings = np.empty(ncross)

    for i in range(ncross):
        crossings[i] = TE_crossing(V, i)

    return crossings


def TM_mode_plot(V, n1, n2):

    xi = np.linspace(0, V/2 - abit, 100)
    ellipse = (n1 / n2)**2 * np.sqrt((V / 2)**2 - xi**2)

    aplt = base_mode_plot(V)
    aplt.plot(xi, ellipse, ':k')
    aplt.ylabel(r'$\xi\,\tan\xi$   or   $-\xi\,\cot\xi$   or   $(n_1/n_2)^2\sqrt{V^2/4-\xi^2}$')
    aplt.title('TM Modes in Planar Film Waveguide (V=%.2f)' % V)
    ymax = (n1 / n2)**2 * V / 2
    aplt.ylim(0, ymax + 1)
    aplt.xlim(0, ymax + 1)

    return aplt


def TM_mode(xi, *args):
    V = args[0]
    n1 = args[1]
    n2 = args[2]
    mode = args[3]
    if mode % 2 == 0:
        return xi * np.tan(xi) - (n1 / n2)**2 * np.sqrt((V / 2)**2 - xi**2)
    else:
        return xi / np.tan(xi) + (n1 / n2)**2 * np.sqrt((V / 2)**2 - xi**2)


def TM_crossing(V, n1, n2, i):
    lower = abit + i * np.pi / 2
    upper = min(np.pi / 2 - abit + i * np.pi / 2, V / 2)
    if lower > V / 2:
        return 0  # mode does not exist

    return brentq(TM_mode, lower, upper, args=(V, n1, n2, i))


def TM_crossings(V, n1, n2):
    ncross = int(V / 2 / (np.pi / 2)) + 1
    crossings = np.empty(ncross)

    for i in range(ncross):
        crossings[i] = TM_crossing(V, n1, n2, i)

    return crossings


def TE_field(V, d, x, i):
    xi = TE_crossing(V, i)
    gdby2 = np.sqrt((V / 2)**2 - xi**2)
    gamma = 2 / d * gdby2
    kappa = 2 / d * xi
    Ey = np.empty(len(x))
    if i % 2 == 0:
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