# pylint: disable=invalid-name
# pylint: disable=no-member

"""
Useful routines for noise in optical communications.

See <https://ofiber.readthedocs.io> for usage examples.

Noise currents in a photodetector::

    shot_noise(I0, Idark, bandwidth, M=1, x=0)
    thermal_noise(T, Rload, bandwidth)
    NEP(Responsivity, Rload, Idark, T)
    best_APD_gain(I0, Rload, Idark, x, T)

Bit error rates and the powers needed to reach them::

    BER_at_SNR(snr)
    SNR_at_BER(ber)
    thermal_min_power(bitrate, responsivity, capacitance, T, snr)
    quantum_min_power(bitrate, ber, lambda0)

`BER_at_SNR` and `SNR_at_BER` are exact inverses of one another.

Beware that I0 means different things in two of these.  `shot_noise` wants
the current after avalanche gain, since it divides by M to recover the
primary current, while `best_APD_gain` wants the primary photocurrent
before any gain is applied.

Based on chapter 13 of A. Ghatak, K. Thyagarajan, An Introduction to
Fiber Optics, Cambridge University Press, 1998
"""

import scipy.constants
import scipy.special
import numpy as np

__all__ = (
    "shot_noise",
    "thermal_noise",
    "NEP",
    "best_APD_gain",
    "BER_at_SNR",
    "SNR_at_BER",
    "thermal_min_power",
    "quantum_min_power",
)


def shot_noise(I0, Idark, bandwidth, M=1, x=0):
    """
    Return the noise current associated with shot/poisson noise.

    The avalanche gain M multiplies the primary photocurrent, so I0 is the
    current measured after gain and I0/M recovers the primary current.  Idark
    is the primary dark current and is multiplied along with it.  The factor
    M**(2+x) is the usual M**2 gain together with the excess noise factor
    F(M) = M**x.  With the defaults M=1 and x=0 this is the shot noise of a
    PIN detector, sqrt(2 q (I0 + Idark) B).

    Args:
        I0:    current after avalanche gain  (A)
        Idark:  primary dark current    (A)
        bandwidth:  bandwidth   (Hz)
        M:      APD Gain factor (--)
        x:      excess noise (0.3 for Si, 0.7 for InGaAs, 1.0 for Ge APDs)

    Returns:
        shot_noise       (A)
    """
    q = scipy.constants.e  # Coulomb
    return np.sqrt(2 * q * (I0 / M + Idark) * bandwidth * M ** (2 + x))


def thermal_noise(T, Rload, bandwidth):
    """
    Return the noise current associated with thermal processes.

    Args:
        T:      temperature (Kelvin)
        Rload:  resistance  (Ohms)
        bandwidth:     bandwidth   (Hz)

    Returns:
        thermal_noise       (A)
    """
    k = scipy.constants.k  # J/K
    return np.sqrt(4 * k * T * bandwidth / Rload)


def NEP(Responsivity, Rload, Idark, T):
    """
    Return noise equivalent power.

    Args:
        Responsivity: photodetector response (A/W)
        Rload:        resistance             (Ohms)
        Idark:        dark current           (A)
        T:            temperature            (Kelvin)

    Returns:
        power       (W/sqrt(Hz))
    """
    q = scipy.constants.e  # Coulomb
    k = scipy.constants.k  # J/K
    return 1 / Responsivity * np.sqrt(4 * k * T / Rload + 2 * q * Idark)


def best_APD_gain(I0, Rload, Idark, x, T):
    """
    Return best gain for an avalanche photodiode.

    This is the gain that maximises the signal-to-noise ratio, where more gain
    buys signal until the excess noise F(M) = M**x overtakes the fixed thermal
    noise of the load.  Unlike `shot_noise`, I0 here is the primary
    photocurrent before any gain.

    The excess noise x must be greater than zero; with x=0 more gain always
    helps and there is no optimum.

    Args:
        I0:    primary photocurrent, before gain  (A)
        Rload: load resistance (Ohms)
        Idark: dark current    (A)
        x:     excess noise (0.3 for Si, 0.7 for InGaAs, 1.0 for Ge APDs)
        T:     temperature     (Kelvin)

    Returns:
        optimal gain (--)
    """
    q = scipy.constants.e  # Coulomb
    k = scipy.constants.k  # J/K
    return (4 * k * T / (x * q * Rload * (I0 + Idark))) ** (1 / (x + 2))


def BER_at_SNR(snr):
    """
    Return the bit error rate for a particular signal-to-noise ratio.

    Args:
        snr: signal to noise ratio (--)

    Returns:
        bit error rate              (--)
    """
    return 0.5 * scipy.special.erfc(np.sqrt(snr / 8))


def SNR_at_BER(ber):
    """
    Return the necessary signal-to-noise ratio to achieve specified BER.

    Args:
        ber: bit error rate (--)

    Returns:
        signal to noise ratio (--)
    """
    return 8 * scipy.special.erfcinv(2 * ber) ** 2


def thermal_min_power(bitrate, responsivity, capacitance, T, snr):
    """
    Return the minimum optical power needed to reach a signal-to-noise ratio.

    The assumption is that thermal noise dominates, with the load resistance
    set so the detector bandwidth matches the bit rate.

    Args:
        bitrate:       bits per second           (Hz)
        responsivity:  photodetector response    (A/W)
        capacitance:   photodetector capacitance (Farads)
        T:             temperature               (Kelvin)
        snr:           signal to noise ratio     (--)

    Returns:
        optical power                            (W)
    """
    k = scipy.constants.k  # J/K
    val = 2 * np.pi * k * T * capacitance * snr
    return bitrate / responsivity * np.sqrt(val)


def quantum_min_power(bitrate, ber, lambda0):
    """
    Return the minimum optical power needed to achieve a bit error rate.

    In the quantum limit a receiver errs only when no photon arrives, so a bit
    carrying Np photons has an error rate of exp(-Np)/2.  Inverting that gives
    Np = -ln(2*ber) photons per bit, and the power follows as h*nu*Np*bitrate.

    Args:
        bitrate:   bits per second      (Hz)
        ber:       bit error rate       (--)
        lambda0:   wavelength in vacuum (m)

    Returns:
        optical power                   (W)
    """
    h = scipy.constants.h  # J*s
    c = scipy.constants.c  # m/s
    nu = c / lambda0  # Hz
    Np = -np.log(2 * ber)
    return h * nu * Np * bitrate
