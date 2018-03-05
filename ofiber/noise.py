import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc
from scipy.special import erfcinv

__all__ = [ 'shot_noise',
            'thermal_noise',
            'NEP',
            'best_APD_gain',
            'BER_at_SNR',
            'SNR_at_BER',
            'thermal_min_power',
            'quantum_min_power',
            ]


def shot_noise(I,Idark,bandwidth, M=1, x=0):
    """
    Return the noise current associated with shot/poisson noise
    Inputs:
        I     = current         (A)
        Idark = dark current    (A)
        M     = APD Gain factor (--)
        x     = excess noise (0.3 for Si, 0.7 for InGaAs, 1.0 for Ge APDs)
    Returns:
        shot_noise       (A)
    """
    q = 1.602e-19 #coulomb
    return np.sqrt(2*q*(I/M+Idark)*bandwidth*M**(2+x))


def thermal_noise(T,Rload,bandwidth):
    """
    Return the noise current associated with thermal processes
    Inputs:
        T     = temperature (Kelvin)
        Rload = resistance  (Ohms)
        BW    = bandwidth   (Hz)
    Returns:
        thermal_noise       (A)
    """
    k = 1.38e-23   #J/K
    return np.sqrt(4*k*T*bandwidth/Rload)


def NEP(Responsivity,Rload,Idark,T):
    """
    Return noise equivalent power
    Inputs:
        responsivity = photodetector response (A/W)
        Rload        = resistance             (Ohms)
        Idark        = dark current           (A)
        T            = temperature            (Kelvin)
    Returns:
        power       (W/sqrt(Hz))
    """
    q = 1.602e-19 #coulomb
    k = 1.38e-23   #J/K
    return 1/Responsivity*np.sqrt(4*k*T/Rload + 2*q*Idark)


def best_APD_gain(I,Rload,Idark,x,T):
    q = 1.602e-19  # coulomb
    k = 1.38e-23   # J/K
    return (4*k*T/(x*q*Rload*(I+Idark)))**(1/(x+2))


def BER_at_SNR(snr):
    """
    Return the bit error rate for a particular signal-to-noise ratio
    Inputs:
        snr = signal to noise ratio (--)
    Returns:
        bit error rate              (--)
    """
    return 0.5*erfc(np.sqrt(snr/8))


def SNR_at_BER(ber):
    """
    Return the necessary signal-to-noise ratio to achieve specified BER
    Inputs:
        ber = bit error rate (--)
    Returns:
        signal to noise ratio (--)
    """
    return 8*erfcinv(2*ber)**2


def thermal_min_power(bitrate,responsivity,capacitance,T,snr):
    """
    Return the minimum optical power needed to achieve a signal-to-noise
    ratio assuming that thermal noise dominates
    Inputs:
        bitrate      = bits per second           (Hz)
        responsivity = photodetector response    (A/W)
        capacitance  = photodetector capacitance (Farads)
        T            = temperature               (Kelvin)
        snr          = signal to noise ratio     (--)
    Returns:
        optical power                            (W)
    """
    k = 1.38e-23   #J/K
    val = 2*np.pi*k*T*C*snr
    return bitrate/responsivity*np.sqrt(val)


def quantum_min_power(bitrate, ber, lambda0):
    """
    Return the minimum optical power needed to achieve a bit error rate
    Inputs:
        bitrate  = bits per second      (Hz)
        ber      = bit error rate       (--)
        lambda0  = wavelength in vacuum (m)
    Returns:
        optical power                   (W)
    """
    h = 6.626e-34  # J*s
    c = 2.998e8    # m/s
    nu = c/lambda0 # Hz
    Np = -np.log(2*ber)
    return h*nu*Np*bitrate
