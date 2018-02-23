# Useful routines for index of refraction based on Sellmeier coefficients
#
#   See Fleming 1978
#   See Book on Optical Properties
#
# Scott Prahl
# Feb 2018


import numpy as np

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
# fleming_names = np.array(["Quenched SiO$_2$","13.5% Ge$O_2$","9.1% P$_2$O$_2$","13.3% B$_2$O$_3$","1.0% F","16.9% Na$_2$O : 32.5% B$_2$O$_3$"])
# 
# extract glasses that are not SiO2:GeO2 mixtures
# glass = fleming_all[2:6]
# glass_names = fleming_names[2:6]

other_glass=[
[3.79061862e-03,1.43810462e-02,7.49374334e+01,6.95790000e-01,4.52497000e-01,7.12513000e-01],
[3.83161000e-03,1.52922902e-02,8.27910731e+01,6.90618000e-01,4.01996000e-01,8.98817000e-01],
[4.65492353e-03,1.35629316e-02,9.98741796e+01,6.91116000e-01,3.99166000e-01,8.90423000e-01],
[8.90362088e-03,8.72094500e-03,3.59958241e+01,7.96468000e-01,4.97614000e-01,3.58924000e-01]
]
other_glass_names = np.array(["9.1% P$_2$O$_2$","13.3% B$_2$O$_3$",
                              "1.0% F","16.9% Na$_2$O : 32.5% B$_2$O$_3$"])

def geo2_sio2_glass(x):
    """
    return Sellmeier coefficients for mixed system 
             x GeO_2 : (1 - x)SiO_2 
    where x is the molar fraction of GeO_2 in the mixture
    """
    SA = np.array([0.6961663,0.4079426,0.8974794])
    SL = np.array([0.0684043,0.1162414,9.896161])
    GA = np.array([0.80686642,0.71815848,0.85416831])
    GL = np.array([0.068972606,0.15396605,11.841931])
    a = (SL + x*(GL-SL))**2
    b = abs(SA + x*(GA-SA))
    return np.array([a,b]).reshape(-1)


def sellmeier(a,b,lambda0):
    """
    returns the index of refraction using the Sellmeier equation
    $$ n^2 = 1 + \sum_{i=0}^2 b_i \lambda_0^2/(\lambda_0^2-a_i) $$
    
    a is in [um**2]
    b is in [1/um**2]
    lambda0 is in [um]
    """
    lam2 = lambda0**2 *1e12                   # um**2
    nsq  = 1
    for i in range(3):
        nsq += b[i]*lam2/(lam2-a[i])
        
    return np.sqrt(nsq)


def d_sellmeier(a,b,lambda0):
    """
    returns the first derivative (wrt to wavelength) of the Sellmeier equation
    
    a is in [um**2]
    b is in [1/um**2]
    lambda0 is in [um]
    returned value dn/dlambda is in [1/m]
    """
    n1   = sellmeier_indexx(a,b,lambda0)
    lam  = lambda0*1e6  # microns
    lam2 = lam**2
    
    dy  = 0
    for i in range(3):
        dy -= a[i]*b[i]/(lam2-a[i])**2  # 1/um**2

    dy *= lam/n1                        # 1/um
    dy *= 1e6                           # 1/m
    
    return dy


def d2_sellmeier(a,b,lambda0):
    """
    returns the second derivative (wrt to wavelength) of the Sellmeier equation
    
    a is in [um**2]
    b is in [1/um**2]
    lambda0 is in [um]
    returned value dn/dlambda is in [1/m]
    """
    n   = sellmeier(a,b,lambda0)
    lam  = lambda0*1e6  # to match Sellmeier Coefficients
    lam2 = lam**2
    
    dy  = 0
    d2y = 0
    for i in range(3):
        dy = a[i]*b[i]/(lam2-a[i])**2                 # 1/um
        d2y += a[i]*b[i]*(3*lam2+a[i])/(lam2-a[i])**3 # 1/um**2

    d2n  = d2y/n - lam2*dy**2/n**3                    # 1/um**2
    d2n *= 1e12                                       # 1/m**2
    
    return d2n
    

def sellmeier_refractive_index(sell,lambda0):
    """
    returns the index of refraction for Sellmeier array, sell, at wavelength lambda0
    
    sell[0:3] is in [um**2]
    sell[3:6] is in [1/um**2]
    lambda0 is in [um]
    returned value is in [s/m**]  ... just multiply be 1e6 to get ps/(km nm)
    """
    return sellmeier(sell[0:3],sell[3:6],lambda0)

