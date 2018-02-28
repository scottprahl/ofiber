# Useful routines for index of refraction based on Sellmeier coefficients
#
# Based on the equation
#
#    $$ n^2 = 1 + \sum_{i=0}^2 b_i \lambda_0^2/(\lambda_0^2-c_i) $$
#
#   See Fleming 1978
#   See Book on Optical Properties
#
# Scott Prahl
# Feb 2018

__all__ = [ 'all_glass_names',
			'd2n',
			'dn',
			'doped_glass',
			'doped_glass_name',
			'glass',
			'glass_name',
			'n']

import numpy as np

_glass = [
    # format is [c1, c2, c3, b1, b2, b3]
    # where b1, b2, b3 are unitless and c1, c2, c3 have units of [microns**2]
    [4.67914826e-3, 1.35120631e-2, 9.79340025e1,
     6.96166300e-1, 4.07942600e-1, 8.97479400e-1],  # [0] SiO2
    [4.75722038e-3, 2.37055446e-2, 1.40231330e2,
     8.06866420e-1, 7.18158480e-1, 8.54168310e-1],  # [1] GeO2
    [3.79061862e-3, 1.43810462e-2, 7.49374334e1,
     6.95790000e-1, 4.52497000e-1, 7.12513000e-1],  # [2] P2O2
    [3.83161000e-3, 1.52922902e-2, 8.27910731e1,
     6.90618000e-1, 4.01996000e-1, 8.98817000e-1],  # [3] 13.3% B2O3
    [4.65492353e-3, 1.35629316e-2, 9.98741796e1,
     6.91116000e-1, 3.99166000e-1, 8.90423000e-1],  # [4] F
    [8.90362088e-3, 8.72094500e-3, 3.59958241e1,
     7.96468000e-1, 4.97614000e-1, 3.58924000e-1],  # [5] NaO:B2O3
    [4.54079433e-2, -1.35241091e-5, 3.09549568e2,
     1.30355147e-1, 9.13764925e-1, 1.14207828],     # [6] ABCY
    [-1.55905386e-4, 7.32455962e-3, 5.96822762e2,
     2.03974072e-5, 1.25885153, 2.11857374],        # [7] HBL
    [3.31382978e-4, 9.53013988e-3, 3.85595295e2,
     3.50883275e-1, 9.36323861e-1, 1.45963548],     # [8] ZBG
    [1.49169281e-8, 8.95628044e-3, 2.39968296e2,
     3.28391032e-2, 1.25579928, 8.97176663e-1],     # [9] ZBLA
    [-2.40488039e-2, 1.73740457e-2, 4.02611805e2,
     3.05900633e-1, 9.18318740e-1, 1.50695421],     # [10] ZBLAN
    [0.004981838, 0.01375664, 97.93353,
     0.6910021, 0.4022430, 0.9439644],              # [11] 5.2% B2O3
    [0.005202431, 0.01287730, 97.93401,             
     0.7058489, 0.4176021, 0.8952753],              # [12] 10.5% P2O2
    [6.00069867e-3, 2.00179144e-2, 103.560653, 
     1.03961212, 0.231792344, 1.01046945],          # [13] BK7
    [4.67914826e-3, 1.35120631e-2, 97.9340025, 
     0.696166300, 0.407942600, 0.897479400],        # [14] fused silica
    [5.2799261e-3, 1.42382647e-2, 325.017834, 
     1.43134930, 0.65054713, 5.3414021],            # [15] sapphire (ord. wave)
    [5.48041129e-3, 1.47994281e-2, 402.89514, 
     1.5039759, 0.55069141, 6.5927379]              # [16] sapphire (eo wave)
]


all_glass_names = np.array([
    "SiO$_2$",             "GeO$_2$", "9.1% P$_2$O$_2$", 
    "13.3% B$_2$O$_3$",    "1.0% F",  "16.9% Na$_2$O : 32.5% B$_2$O$_3$",
    "ABCY",                "HBL",     "ZBG",  
    "ZBLA",                "ZBLAN",   "5.2% B$_2$O$_3$",
    "10.5% P$_2$O$_2$",    "BK7",     "fused silica",
    "sapphire (ordinary)", "sapphire (exordinary)"
    ])


def glass(i):
    """
    return an array of Sellmeier coefficients for glass with index i
    Use like this
        lambda0 = np.linspace(1000,1700,50)*1e-9 # [m]
        glass = ofiber.refraction.glass(0)       # SiO2
        n = ofiber.refraction.n(glass,lambda0)    
        plt.plot(lambda0*1e9, n)
    """
    return _glass[i]


def glass_name(i):
    """
    return the name of the glass with index i
    (A list of all possible names is in the array 
        ofiber.refraction.all_glass_names
    )
    """
    return all_glass_names[i]


def doped_glass(x):
    """
    return Sellmeier coefficients for mixed system 
             x GeO_2 : (1 - x)SiO_2 
    where x is the molar fraction of GeO_2 in the mixture
    """
    SA = np.array([0.6961663, 0.4079426, 0.8974794])
    SL = np.array([0.0684043, 0.1162414, 9.896161])
    GA = np.array([0.80686642, 0.71815848, 0.85416831])
    GL = np.array([0.068972606, 0.15396605, 11.841931])
    a = (SL + x * (GL - SL))**2
    b = abs(SA + x * (GA - SA))
    return np.array([a, b]).reshape(-1)


def doped_glass_name(x):
    """
    return name for the mixed system 
             x GeO_2 : (1 - x)SiO_2 
    where x is the molar fraction of GeO_2 in the mixture
    """
    if x == 0:
        return r'SiO$_2$'
    if x == 1:
        return r'GeO$_2$'

    return r'%.2f GeO$_2$ : %.2f SiO$_2$' % (x, 1 - x)


def _sellmeier(b, c, lambda0):
    """
    returns the index of refraction using the Sellmeier equation
    $$ n^2 = 1 + \sum_{i=0}^2 b_i \lambda_0^2/(\lambda_0^2-c_i) $$

    b is unitless
    c is in [um**2]
    lambda0 is in [m]
    """
    lam2 = lambda0**2 * 1e12                   # um**2
    nsq = 1
    for i in range(3):
        nsq += b[i] * lam2 / (lam2 - c[i])

    return np.sqrt(nsq)


def _d_sellmeier(b, c, lambda0):
    """
    returns the first derivative (wrt to wavelength) of the Sellmeier equation

    b is unitless
    c is in [um**2]
    lambda0 is in [m]
    returned value dn/dlambda is in [1/m]
    """
    n1 = _sellmeier(b, c, lambda0)
    lam = lambda0 * 1e6  # microns
    lam2 = lam**2

    dy = 0
    for i in range(3):
        dy -= b[i] * c[i] / (lam2 - c[i])**2  # 1/um**2

    dy *= lam / n1                            # 1/um
    dy *= 1e6                                 # 1/m

    return dy


def _d2_sellmeier(b, c, lambda0):
    """
    returns the second derivative (wrt to wavelength) of the Sellmeier equation

    b is unitless
    c is in [um**2]
    lambda0 is in [m]
    returned value d^2 n/dlambda^2 is in [1/m**2]
    """
    n = _sellmeier(b, c, lambda0)
    lam = lambda0 * 1e6  # to match Sellmeier Coefficients
    lam2 = lam**2

    dy = 0
    d2y = 0
    for i in range(3):
        dy = b[i] * c[i] / (lam2 - c[i])**2                        # 1/um
        d2y += b[i] * c[i] * (3 * lam2 + c[i]) / (lam2 - c[i])**3  # 1/um**2

    d2n = d2y / n - lam2 * dy**2 / n**3                            # 1/um**2
    d2n *= 1e12                                                    # 1/m**2

    return d2n


def n(glass, lambda0):
    """
    returns the index of refraction for Sellmeier array, glass, at wavelength lambda0

    glass is an array obtained from ofiber.refraction.glass(i)
    lambda0 is in [m]
    returned value is in [s/m**]  ... just multiply be 1e6 to get ps/(km nm)
    """
    return _sellmeier(glass[3:6], glass[0:3], lambda0)


def dn(glass, lambda0):
    """
    returns the first derivative (wrt to wavelength) of the Sellmeier equation

    glass is an array obtained from ofiber.refraction.glass(i)
    lambda0 is in [m]
    returned value is in [1/m]
    """
    return _d_sellmeier(glass[3:6], glass[0:3], lambda0)


def d2n(glass, lambda0):
    """
    returns the second derivative (wrt to wavelength) of the Sellmeier equation

    glass is an array obtained from ofiber.refraction.glass(i)
    lambda0 is in [m]
    returned value is in [1/m**2]
    """
    return _d2_sellmeier(glass[3:6], glass[0:3], lambda0)


# code used to generate the Sellmeier coefficients for other_glass below
# data straight from fleming 1978
# flemingzz=np.array([
#     [0.696750, 0.069066, 0.408218, 0.115662, 0.890815, 9.900559],
#     [0.711040, 0.064270, 0.451885, 0.129408, 0.704048, 9.425478],
#     [0.695790, 0.061568, 0.452497, 0.119921, 0.712513, 8.656641],
#     [0.690618, 0.061900, 0.401996, 0.123662, 0.898817, 9.098960],
#     [0.691116, 0.068227, 0.399166, 0.116460, 0.890423, 9.993707],
#     [0.796468, 0.094359, 0.497614, 0.093386, 0.358924, 5.999652]]).T
#
# rearrange so each row is a type of glass with Sellmeier coefficients needed
# b1,l1,b2,l2,b3,l3=flemingzz
# a1=l1**2
# a2=l2**2
# a3=l3**2
# fleming_all=np.array([a1,a2,a3,b1,b2,b3]).T
# fleming_names = np.array(["Quenched SiO$_2$","13.5% Ge$O_2$",\
#     "9.1% P$_2$O$_2$","13.3% B$_2$O$_3$","1.0% F",\
#     "16.9% Na$_2$O : 32.5% B$_2$O$_3$"])
#
# extract glasses that are not SiO2:GeO2 mixtures
# glass = fleming_all[2:6]
# glass_names = fleming_names[2:6]


# data straight from Ghatak
# pure_sio2 = [0.004679148,0.01351206,97.93400,0.6961663,0.4079426,0.8974794]
# geo2_63 = [0.007290464,0.01050294,97.93428,0.7083952,0.4203993,0.8663412]
# geo2_193 = [0.005847345,0.01552717,97.93484,0.7347008,0.4461191,0.8081698]
# b203 = [0.004981838,0.01375664,97.93353,0.6910021,0.4022430,0.9439644]
# p2o3 = [0.005202431,0.01287730,97.93401,0.7058489,0.4176021,0.8952753]


# def paek_refraction(lambda0):
#     """
#     Return the index of refraction of silica fiber at
#     vacuum wavelength lambda0 [m]
#     result is dimensionless
#     """
#     ell = 0.035       # um**2
#     c0 = 1.4508554
#     c1 = -0.0031268
#     c2 = -0.0000381
#     c3 = 0.0030270
#     c4 = -0.0000779
#     c5 = 0.0000018
#     lam = lambda0 * 1e6   # now in microns
#     den = lam**2 - ell
#     n = c0 + c1 * lam**2 + c2 * lam**4 + c3 / den + c4 / den**2 + c5 / den**3
#     return n
#
#
# def d_paek_refraction(lambda0):
#     """
#     Return the first derivative of index of refraction of silica fiber at
#     vacuum wavelength lambda0 [m] with respect to wavelength
#     result has units of [1/m]
#     """
#     ell = 0.035
#     c1 = -0.0031268
#     c2 = -0.0000381
#     c3 = 0.0030270
#     c4 = -0.0000779
#     c5 = 0.0000018
#     lam = lambda0 * 1e6   # microns
#     den = lam**2 - ell
#     n = 2 * c1 * lam + 4 * c2 * lam**3 - 2 * c3 * lam / den**2
#     n += - 4 * c4 * lam / den**3 - 6 * c5 * lam / den**4      # um**-1
#     n *= 1e6                                                  # m**-1
#     return n
#
#
# def d2_paek_refraction(lambda0):
#     """
#     Return the second derivative of index of refraction of silica fiber at
#     vacuum wavelength lambda0 [m] with respect to wavelength
#     result has units of [1/m**2]
#     """
#     ell = 0.035
#     c1 = -0.0031268
#     c2 = -0.0000381
#     c3 = 0.0030270
#     c4 = -0.0000779
#     c5 = 0.0000018
#     lam = lambda0 * 1e6   # microns
#     den = lam**2 - ell
#     n = 2 * c1 + 12 * c2 * lam**2 - 2 * c3 / \
#         den**2 - 4 * c4 / den**3 - 6 * c5 / den**4
#     n += 48 * c5 * lam**2 / den**5 + 24 * c4 * \
#         lam**2 / den**4 + 8 * c3 * lam**2 / den**3
#     n *= 1e12           # m**-2
#     return n
#
# ### method for converting from mendez coefficients to sellmeier coefficients
# ### just a simple curve fit
#
# from scipy.optimize import curve_fit
#
# def sell(x,a1,a2,a3,b1,b2,b3):
#     x2 = x**2
#     return np.sqrt(1 + b1*x2/(x2 - a1) + b2*x2/(x2 - a2) + b3*x2/(x2 - a3))
#
#
# lambda0 = np.linspace(1500,3500,50)*1e-9 # [m]
# for i in range(len(ofr.mendez_glass)):
#     glass = ofr.mendez_glass[i]
#     n = ofr.mendez_refraction(glass,lambda0)
#     popt, pcov = curve_fit(sell, lambda0*1e6, n, p0=ofr.sell_glass[0])
#     plt.plot(lambda0*1e9, n-sell(lambda0*1e6, *popt), 'r-', label='fit')
#     print(popt)

# mendez_glass =[
# [7.67742e-6, 2.16195e-3, 1.42969, -1.28304e-3, -5.35487e-6],
# [-28.61020e-6, 3.11470e-3, 1.50294, -1.17821e-3, -2.64123e-6],
# [93.67070e-6, 2.94329e-3, 1.51236, -1.25045e-3, -4.01026e-6],
# [-300.80370e-6, 4.03214e-3, 1.51272, -1.21921e-3, -6.77630e-6],
# [93.67070e-6, 2.94329e-3, 1.49136, -1.25045e-3, -4.01026e-6]
# ]
# mendezall_glass_namess=["ABCY", "HBL",  "ZBG",  "ZBLA", "ZBLAN"]
#
# def mendez_refraction(glass, lambda0):
#     """
#     returns the index of refraction using the Mendez equation
#     lambda0 is in [m]
#     """
#     lam = lambda0 * 1e6   # microns
#     sum = glass[0] * lam**-4
#     sum += glass[1] * lam**-2
#     sum += glass[2]
#     sum += glass[3] * lam**2
#     sum += glass[4] * lam**4
#     return sum
#
#
# def d_mendez_refraction(glass, lambda0):
#     """
#     returns first derivative of the  index of refraction using the Mendez eqn
#     lambda0 is in [m]
#     """
#     sum = 0
#     lam = lambda0 * 1e6   # microns
#     sum = -4 * glass[0] * lam**-5
#     sum += -2 * glass[1] * lam**-3
#     sum += 2 * glass[3] * lam
#     sum += 4 * glass[4] * lam**3
#     sum *= 1e6  # [1/m]
#     return sum
#
#
# def d2_mendez_refraction(glass, lambda0):
#     """
#     returns second derivative of the  index of refraction using the Mendez eqn
#     lambda0 is in [m]
#     """
#     sum = 0
#     lam = lambda0 * 1e6   # microns
#     sum = 20 * glass[0] * lam**-6
#     sum += 6 * glass[1] * lam**-4
#     sum += 2 * glass[3]
#     sum += 12 * glass[4] * lam**2
#     sum *= 1e12  # [1/m]
#     return sum
