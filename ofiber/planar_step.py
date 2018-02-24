# Useful routines for step-index planar waveguides based on chapter 7 of
#
#    A. Ghatak, K. Thyagarajan, An Introduction to Fiber Optics,
#    Cambridge University Press, 1998
#
# Scott Prahl
# Feb 2018

import numpy
import matplotlib.pyplot
import scipy.optimize


def _base_mode_plot(V):
    """ Create a plot with symmetric and asymmetric modes
        This is used for both TE and TM plots
    """
    abit = 1e-5
    
    fig, ax = matplotlib.pyplot.subplots(figsize=(8, 8))
    ax.set_aspect('equal')

    xi = numpy.linspace(abit, numpy.pi / 2 - abit, 50)
    while xi[0] < V / 2:
        matplotlib.pyplot.plot(xi, xi * numpy.tan(xi), color='red')
        xi += numpy.pi / 2

    xi = numpy.linspace(numpy.pi / 2, numpy.pi - abit, 50)
    while xi[0] < V / 2:
        matplotlib.pyplot.plot(xi, -xi / numpy.tan(xi), '--b')
        xi += numpy.pi

    matplotlib.pyplot.xlabel(r'$\xi=(d/2)\sqrt{k_0^2n_1^2-\beta^2}=(d/2)\kappa$')

    return matplotlib.pyplot


def TE_mode_plot(V):
    abit = 1e-5
    xi = numpy.linspace(0, V / 2 - abit, 100)
    circle = numpy.sqrt((V / 2)**2 - xi**2)

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
        return xi * numpy.tan(xi) - numpy.sqrt((V / 2)**2 - xi * xi)
    else:
        return xi / numpy.tan(xi) + numpy.sqrt((V / 2)**2 - xi * xi)


def TE_crossing(V, mode):
    abit = 1e-5
    lo = abit + mode * numpy.pi / 2
    hi = min(numpy.pi / 2 - abit + mode * numpy.pi / 2, V / 2)
    if lo > V / 2:
        return 0   # mode does not exist

    try:
        b = scipy.optimize.brentq(_TE_mode, lo, hi, args=(V, mode))
    except ValueError:  # happens when both hi and lo values have same sign
        b = 0           # therefore no such mode exists

    return b


def TE_crossings(V):
    ncross = int(V / 2 / (numpy.pi / 2)) + 1
    crossings = numpy.empty(ncross)

    for mode in range(ncross):
        crossings[mode] = TE_crossing(V, mode)

    return crossings


def TM_mode_plot(V, n1, n2):
    abit = 1e-5
    xi = numpy.linspace(0, V / 2 - abit, 100)
    ellipse = (n1 / n2)**2 * numpy.sqrt((V / 2)**2 - xi**2)

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
        return xi * numpy.tan(xi) - (n1/n2)**2 * numpy.sqrt((V/2)**2 - xi**2)
    else:
        return xi / numpy.tan(xi) + (n1/n2)**2 * numpy.sqrt((V/2)**2 - xi**2)


def TM_crossing(V, n1, n2, mode):
    abit = 1e-5
    lo = abit + mode * numpy.pi / 2
    hi = min(numpy.pi / 2 - abit + mode * numpy.pi / 2, V / 2)
    if lo > V / 2:
        return 0  # mode does not exist

    try:
        b = scipy.optimize.brentq(_TM_mode, lo, hi, args=(V, n1, n2, mode))
    except ValueError:  # happens when both hi and lo values have same sign
        b = 0           # therefore no such mode exists

    return b


def TM_crossings(V, n1, n2):
    ncross = int(V / 2 / (numpy.pi / 2)) + 1
    crossings = numpy.empty(ncross)

    for mode in range(ncross):
        crossings[mode] = TM_crossing(V, n1, n2, mode)

    return crossings


def TE_field(V, d, x, mode):
    xi = TE_crossing(V, mode)
    gdby2 = numpy.sqrt((V / 2)**2 - xi**2)
    gamma = 2 / d * gdby2
    kappa = 2 / d * xi
    Ey = numpy.empty(len(x))
    if mode % 2 == 0:
        A = 1
        C = A * numpy.cos(xi) / numpy.exp(-gdby2)
        for j in range(len(x)):
            if abs(x[j]) < d / 2:
                Ey[j] = A * numpy.cos(x[j] * kappa)
            else:
                Ey[j] = C * numpy.exp(-gamma * abs(x[j]))
    else:
        B = 1
        D = B * numpy.sin(xi) / numpy.exp(-gdby2)
        for j in range(len(x)):
            if abs(x[j]) < d / 2:
                Ey[j] = B * numpy.sin(kappa * x[j])
            else:
                Ey[j] = D * numpy.sign(x[j]) * numpy.exp(-gamma * abs(x[j]))

    return Ey


def TE_propagation_constant(V, mode):
    b = numpy.empty(len(V))
    for i in range(len(V)):
        xi = TE_crossing(V[i], mode)
        if xi == 0:
            b[i] = 0
        else:
            b[i] = 1 - (2 * xi / V[i])**2
    return b


def TM_propagation_constant(V, n1, n2, mode):
    b = numpy.empty(len(V))
    for i in range(len(V)):
        xi = TM_crossing(V[i], n1, n2, mode)
        if xi == 0 :
            b[i] = 0
        else :
            b[i] = 1 - (2 * xi / V[i])**2
    return b


