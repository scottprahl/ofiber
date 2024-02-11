'''
Far field radiation of cylindrical step-index fiber LP modes.
Adapted from Chin-Lin Chen, 'Foundations for Guided-Wave Optics', (2006), section 10.2.1.
'''
import numpy as np
import matplotlib.pyplot as plt
import scipy.special

def FFazim(l,lmbd,a,Theta,V,b):
    '''
    Calculate $F_{l}(\Theta)$ from Chen Eq. (10.13)
    l: azimuthal mode number $l$
    Theta: polar angle $\Theta$ in radians'''
    k = 2*np.pi/lmbd
    # V = ofiber.basics.V_parameter(a, NA, lmbd)
    ka = np.multiply(k,a)
    kasin = np.multiply(ka,np.sin(Theta))
    Vsqrt1mb = V*np.sqrt(1-b)
    Jl = scipy.special.jv(l, kasin)
    Jlp1 = scipy.special.jv(l+1,kasin)
    kassq = np.square(kasin)
    Flden1 = np.subtract(np.square(Vsqrt1mb),kassq)
    Flden2 = np.add(V**2*b, kassq)
    Flden = np.multiply(Flden1,Flden2)
    Flnum1 = np.multiply(kasin,Jlp1)
    VsqrtJlp1 = np.multiply(Vsqrt1mb,scipy.special.jv(l+1,Vsqrt1mb))
    JlkovJl = np.divide(Jl,scipy.special.jv(l, Vsqrt1mb))
    Flnum2 = np.multiply(JlkovJl,VsqrtJlp1)
    Flnum = np.subtract(Flnum1,Flnum2)
    return np.divide(Flnum,Flden)

def IrradFFx(R,Theta,Phi,l,lmbd,a,V,b):
    '''
    Chen Eq. (10.12), absolute value squared.
    We divide by $E_{l}^2$.
    FFl : output of FFazim
    '''
    FFl = np.square(FFazim(l,lmbd,a,Theta,V,b))
    k = 2*np.pi/lmbd
    kaVcub = np.power(np.multiply(k*a,V),4)
    cossq = np.square(np.cos(np.multiply(l,Phi)))
    return np.divide(np.multiply(FFl,np.multiply(kaVcub,cossq)),np.square(k*R))


def IrradFFxazint(R,Theta,l,lmbd,a,V,b):
    '''
    Chen Eq. (10.12), absolute value squared.
    We divide by $E_{l}^2$.
    FFl : output of FFazim
    The integral from 0 to $2\pi$ of $\cos^2 l\Phi$ is simply $\pi$.
    '''
    FFl = np.square(FFazim(l,lmbd,a,Theta,V,b))
    k = 2*np.pi/lmbd
    kaVcub = np.power(np.multiply(k*a,V),4)
    return np.divide(np.multiply(FFl,np.multiply(kaVcub,np.pi)),np.square(k*R))
