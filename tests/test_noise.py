# pylint: disable=invalid-name
# pylint: disable=no-member
"""
Tests for the ofiber.noise module.

No notebook prints numbers from this module, so the references are textbook
identities from chapter 13 of Ghatak & Thyagarajan, *An Introduction to Fiber
Optics*, plus a few relations the routines have to satisfy jointly:

* `BER_at_SNR` and `SNR_at_BER` are exact inverses
* `shot_noise` collapses to the PIN result sqrt(2 q (I0+Idark) B) at M=1, x=0
* `best_APD_gain` really is the gain that maximises the signal-to-noise ratio,
  checked by sweeping M through `shot_noise` and `thermal_noise` and finding
  the maximum numerically
* `quantum_min_power` delivers exactly the photons per bit that give the
  requested error rate under ber = exp(-Np)/2

The last two are the strongest, since neither can pass unless several
routines agree with one another.
"""

import numpy as np
import pytest
import scipy.constants
import scipy.special
import ofiber

BANDWIDTH = 1e9  # Hz
T_ROOM = 300.0  # K
R_LOAD = 1e4  # ohm
I_PRIMARY = 1e-7  # A
I_DARK = 1e-9  # A
RESPONSIVITY = 0.8  # A/W


# --------------------------------------------------------------------------
# shot noise
# --------------------------------------------------------------------------


def test_shot_noise_reduces_to_the_pin_result():
    """Verify the M=1, x=0 defaults give sqrt(2 q (I0+Idark) B)."""
    expected = np.sqrt(2 * scipy.constants.e * (I_PRIMARY + I_DARK) * BANDWIDTH)
    assert ofiber.shot_noise(I_PRIMARY, I_DARK, BANDWIDTH) == pytest.approx(expected, rel=1e-12, abs=0)


def test_shot_noise_grows_with_current_and_bandwidth():
    """Verify shot noise rises as the square root of current and of bandwidth."""
    assert ofiber.shot_noise(4e-7, 0, BANDWIDTH) == pytest.approx(2 * ofiber.shot_noise(1e-7, 0, BANDWIDTH))
    assert ofiber.shot_noise(1e-7, 0, 4 * BANDWIDTH) == pytest.approx(2 * ofiber.shot_noise(1e-7, 0, BANDWIDTH))


def test_dark_current_adds_shot_noise():
    """Verify a dark current can only increase the noise."""
    assert ofiber.shot_noise(I_PRIMARY, I_DARK, BANDWIDTH) > ofiber.shot_noise(I_PRIMARY, 0, BANDWIDTH)


def test_shot_noise_divides_the_output_current_by_the_gain():
    """Verify I0 is the current after gain, so I0/M is the primary current."""
    M = 10.0
    with_gain = ofiber.shot_noise(M * I_PRIMARY, I_DARK, BANDWIDTH, M, 0)
    expected = np.sqrt(2 * scipy.constants.e * (I_PRIMARY + I_DARK) * BANDWIDTH * M**2)
    assert with_gain == pytest.approx(expected, rel=1e-12, abs=0)


@pytest.mark.parametrize("x", [0.0, 0.3, 0.7, 1.0])
def test_excess_noise_raises_the_gain_exponent(x):
    """Verify the M**(2+x) scaling, with x the excess noise factor exponent."""
    M = 8.0
    ratio = ofiber.shot_noise(M * I_PRIMARY, I_DARK, BANDWIDTH, M, x) / ofiber.shot_noise(
        M * I_PRIMARY, I_DARK, BANDWIDTH, M, 0
    )
    assert ratio == pytest.approx(M ** (x / 2))


def test_noisier_detectors_have_larger_excess_noise():
    """Verify a Ge APD is noisier than silicon at the same gain."""
    M = 20.0
    silicon = ofiber.shot_noise(M * I_PRIMARY, I_DARK, BANDWIDTH, M, 0.3)
    germanium = ofiber.shot_noise(M * I_PRIMARY, I_DARK, BANDWIDTH, M, 1.0)
    assert germanium > silicon


# --------------------------------------------------------------------------
# thermal noise
# --------------------------------------------------------------------------


def test_thermal_noise_matches_its_closed_form():
    """Verify the Johnson noise expression sqrt(4 k T B / R)."""
    expected = np.sqrt(4 * scipy.constants.k * T_ROOM * BANDWIDTH / R_LOAD)
    assert ofiber.thermal_noise(T_ROOM, R_LOAD, BANDWIDTH) == pytest.approx(expected, rel=1e-12, abs=0)


def test_thermal_noise_scales_as_the_square_root_of_temperature():
    """Verify quadrupling the temperature doubles the noise current."""
    assert ofiber.thermal_noise(4 * T_ROOM, R_LOAD, BANDWIDTH) == pytest.approx(
        2 * ofiber.thermal_noise(T_ROOM, R_LOAD, BANDWIDTH)
    )


def test_a_larger_load_resistor_is_quieter():
    """Verify thermal noise falls as the load resistance rises."""
    R = np.array([1e3, 1e4, 1e5, 1e6])
    assert np.all(np.diff(ofiber.thermal_noise(T_ROOM, R, BANDWIDTH)) < 0)


def test_thermal_noise_vanishes_at_absolute_zero():
    """Verify no thermal agitation means no thermal noise."""
    assert ofiber.thermal_noise(0.0, R_LOAD, BANDWIDTH) == 0


# --------------------------------------------------------------------------
# noise equivalent power
# --------------------------------------------------------------------------


def test_nep_matches_its_closed_form():
    """Verify NEP is the total noise current density divided by the responsivity."""
    expected = np.sqrt(4 * scipy.constants.k * T_ROOM / R_LOAD + 2 * scipy.constants.e * I_DARK) / RESPONSIVITY
    assert ofiber.NEP(RESPONSIVITY, R_LOAD, I_DARK, T_ROOM) == pytest.approx(expected, rel=1e-12, abs=0)


def test_nep_improves_with_responsivity():
    """Verify a more responsive detector needs less power to match its noise."""
    assert ofiber.NEP(2 * RESPONSIVITY, R_LOAD, I_DARK, T_ROOM) == pytest.approx(
        ofiber.NEP(RESPONSIVITY, R_LOAD, I_DARK, T_ROOM) / 2
    )


def test_nep_worsens_with_dark_current_and_temperature():
    """Verify both noise sources push the noise equivalent power up."""
    base = ofiber.NEP(RESPONSIVITY, R_LOAD, I_DARK, T_ROOM)
    assert ofiber.NEP(RESPONSIVITY, R_LOAD, 10 * I_DARK, T_ROOM) > base
    assert ofiber.NEP(RESPONSIVITY, R_LOAD, I_DARK, 2 * T_ROOM) > base


# --------------------------------------------------------------------------
# optimum avalanche gain
# --------------------------------------------------------------------------


@pytest.mark.parametrize("x", [0.3, 0.7, 1.0])
def test_best_apd_gain_maximises_the_signal_to_noise_ratio(x):
    """
    Verify the reported gain is where a numerical sweep finds the SNR peak.

    This ties best_APD_gain, shot_noise and thermal_noise together: the
    optimum is where the growing excess noise M**x overtakes the fixed
    thermal noise of the load.
    """
    M = np.linspace(1, 400, 400001)
    signal = M * I_PRIMARY
    shot = ofiber.shot_noise(signal, I_DARK, BANDWIDTH, M, x)
    thermal = ofiber.thermal_noise(T_ROOM, R_LOAD, BANDWIDTH)
    snr = signal**2 / (shot**2 + thermal**2)
    assert M[np.argmax(snr)] == pytest.approx(ofiber.best_APD_gain(I_PRIMARY, R_LOAD, I_DARK, x, T_ROOM), rel=1e-3)


def test_best_apd_gain_matches_its_closed_form():
    """Verify the documented (4kT/(x q R (I0+Id)))**(1/(x+2))."""
    x = 0.7
    thermal = 4 * scipy.constants.k * T_ROOM
    shot = x * scipy.constants.e * R_LOAD * (I_PRIMARY + I_DARK)
    assert ofiber.best_APD_gain(I_PRIMARY, R_LOAD, I_DARK, x, T_ROOM) == pytest.approx(
        (thermal / shot) ** (1 / (x + 2))
    )


def test_noisier_detectors_want_less_gain():
    """Verify a larger excess noise exponent lowers the optimum gain."""
    gains = [ofiber.best_APD_gain(I_PRIMARY, R_LOAD, I_DARK, x, T_ROOM) for x in (0.3, 0.7, 1.0)]
    assert np.all(np.diff(gains) < 0)


def test_a_brighter_signal_wants_less_gain():
    """Verify more photocurrent means thermal noise matters less, so less gain helps."""
    gains = [ofiber.best_APD_gain(I0, R_LOAD, I_DARK, 0.7, T_ROOM) for I0 in (1e-8, 1e-7, 1e-6)]
    assert np.all(np.diff(gains) < 0)


def test_a_noisier_load_wants_more_gain():
    """Verify a small load resistor, which is noisy, calls for more gain."""
    assert ofiber.best_APD_gain(I_PRIMARY, 1e3, I_DARK, 0.7, T_ROOM) > ofiber.best_APD_gain(
        I_PRIMARY, 1e5, I_DARK, 0.7, T_ROOM
    )


# --------------------------------------------------------------------------
# bit error rate
# --------------------------------------------------------------------------


@pytest.mark.parametrize("snr", [10.0, 50.0, 100.0, 144.0, 200.0])
def test_ber_and_snr_are_exact_inverses(snr):
    """Verify SNR_at_BER undoes BER_at_SNR."""
    assert ofiber.SNR_at_BER(ofiber.BER_at_SNR(snr)) == pytest.approx(snr, rel=1e-9)


@pytest.mark.parametrize("ber", [1e-3, 1e-6, 1e-9, 1e-12])
def test_snr_and_ber_are_exact_inverses_the_other_way(ber):
    """Verify BER_at_SNR undoes SNR_at_BER."""
    assert ofiber.BER_at_SNR(ofiber.SNR_at_BER(ber)) == pytest.approx(ber, rel=1e-9, abs=0)


def test_ber_falls_as_the_signal_to_noise_ratio_rises():
    """Verify a cleaner signal makes fewer errors."""
    snr = np.linspace(1, 200, 100)
    assert np.all(np.diff(ofiber.BER_at_SNR(snr)) < 0)


def test_ber_is_one_half_without_any_signal():
    """Verify a receiver with no signal guesses, erring half the time."""
    assert ofiber.BER_at_SNR(0.0) == pytest.approx(0.5)


def test_a_standard_error_rate_needs_the_textbook_signal_to_noise_ratio():
    """Check the SNR of 144, or 21.6 dB, that a BER of 1e-9 famously requires."""
    assert ofiber.SNR_at_BER(1e-9) == pytest.approx(144, rel=0.02)
    assert 10 * np.log10(ofiber.SNR_at_BER(1e-9)) == pytest.approx(21.6, abs=0.1)


def test_ber_matches_the_complementary_error_function():
    """Verify the documented 0.5*erfc(sqrt(snr/8))."""
    assert ofiber.BER_at_SNR(64.0) == pytest.approx(0.5 * scipy.special.erfc(np.sqrt(8.0)))


# --------------------------------------------------------------------------
# minimum detectable power
# --------------------------------------------------------------------------


def test_thermal_min_power_matches_its_closed_form():
    """Verify the documented bitrate/responsivity * sqrt(2 pi k T C snr)."""
    C, snr = 1e-12, 144.0
    expected = 1e9 / RESPONSIVITY * np.sqrt(2 * np.pi * scipy.constants.k * T_ROOM * C * snr)
    assert ofiber.thermal_min_power(1e9, RESPONSIVITY, C, T_ROOM, snr) == pytest.approx(expected)


def test_thermal_min_power_rises_with_bitrate_and_snr():
    """Verify a faster link, or a stricter error rate, needs more power."""
    base = ofiber.thermal_min_power(1e9, RESPONSIVITY, 1e-12, T_ROOM, 144.0)
    assert ofiber.thermal_min_power(4e9, RESPONSIVITY, 1e-12, T_ROOM, 144.0) == pytest.approx(4 * base)
    assert ofiber.thermal_min_power(1e9, RESPONSIVITY, 1e-12, T_ROOM, 4 * 144.0) == pytest.approx(2 * base)


def test_thermal_min_power_falls_with_responsivity():
    """Verify a better detector needs less light."""
    base = ofiber.thermal_min_power(1e9, RESPONSIVITY, 1e-12, T_ROOM, 144.0)
    assert ofiber.thermal_min_power(1e9, 2 * RESPONSIVITY, 1e-12, T_ROOM, 144.0) == pytest.approx(base / 2)


@pytest.mark.parametrize("ber", [1e-6, 1e-9, 1e-12])
def test_quantum_min_power_delivers_the_photons_that_error_rate_needs(ber):
    """
    Verify the power carries exactly Np photons per bit, with ber = exp(-Np)/2.

    This inverts the physics the routine is built on rather than repeating its
    arithmetic.
    """
    bitrate, lambda0 = 1e9, 1.55e-6
    photon_energy = scipy.constants.h * scipy.constants.c / lambda0
    Np = ofiber.quantum_min_power(bitrate, ber, lambda0) / photon_energy / bitrate
    assert 0.5 * np.exp(-Np) == pytest.approx(ber, rel=1e-9, abs=0)


def test_quantum_min_power_is_proportional_to_bitrate():
    """Verify twice the bits need twice the power at a fixed error rate."""
    base = ofiber.quantum_min_power(1e9, 1e-9, 1.55e-6)
    assert ofiber.quantum_min_power(2e9, 1e-9, 1.55e-6) == pytest.approx(2 * base)


def test_shorter_wavelengths_carry_more_energy_per_photon():
    """Verify the quantum limit rises as the wavelength shortens."""
    lam = np.array([1.55e-6, 1.31e-6, 0.85e-6])
    powers = [ofiber.quantum_min_power(1e9, 1e-9, wavelength) for wavelength in lam]
    assert np.all(np.diff(powers) > 0)


def test_a_stricter_error_rate_needs_more_photons():
    """Verify demanding fewer errors costs power."""
    powers = [ofiber.quantum_min_power(1e9, ber, 1.55e-6) for ber in (1e-6, 1e-9, 1e-12)]
    assert np.all(np.diff(powers) > 0)


def test_the_quantum_limit_is_below_the_thermal_one():
    """Verify a real thermally limited receiver needs more power than the quantum limit."""
    quantum = ofiber.quantum_min_power(1e9, 1e-9, 1.55e-6)
    thermal = ofiber.thermal_min_power(1e9, RESPONSIVITY, 1e-12, T_ROOM, ofiber.SNR_at_BER(1e-9))
    assert quantum < thermal


# --------------------------------------------------------------------------
# module surface
# --------------------------------------------------------------------------


def test_everything_in_all_is_reachable_from_the_package():
    """Verify each exported name survives the star import in __init__.py."""
    for name in ofiber.noise.__all__:
        assert hasattr(ofiber, name), "%s is in __all__ but missing from ofiber" % name


# --------------------------------------------------------------------------
# physical constants
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "routine, args, exact",
    [
        (ofiber.thermal_noise, (T_ROOM, R_LOAD, BANDWIDTH), scipy.constants.k),
        (ofiber.shot_noise, (I_PRIMARY, I_DARK, BANDWIDTH), scipy.constants.e),
    ],
)
def test_noise_currents_use_codata_constants(routine, args, exact):
    """Verify the constants match scipy.constants, which the rest of ofiber uses."""
    if exact == scipy.constants.k:
        expected = np.sqrt(4 * exact * T_ROOM * BANDWIDTH / R_LOAD)
    else:
        expected = np.sqrt(2 * exact * (I_PRIMARY + I_DARK) * BANDWIDTH)
    # abs=0 matters: these currents are ~1e-9, so approx's default abs=1e-12
    # floor would otherwise swallow the discrepancy entirely
    assert routine(*args) == pytest.approx(expected, rel=1e-9, abs=0)


def test_quantum_min_power_uses_codata_constants():
    """Verify the photon energy uses scipy.constants rather than four-digit values."""
    lambda0, bitrate, ber = 1.55e-6, 1e9, 1e-9
    expected = scipy.constants.h * scipy.constants.c / lambda0 * -np.log(2 * ber) * bitrate
    assert ofiber.quantum_min_power(bitrate, ber, lambda0) == pytest.approx(expected, rel=1e-9, abs=0)
