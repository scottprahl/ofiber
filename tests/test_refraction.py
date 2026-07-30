# pylint: disable=invalid-name
"""
Tests for the ofiber.refraction module.

Reference values come from three independent places:

* Malitson, J. Opt. Soc. Am. 55, 1205 (1965) for fused silica.
* The Schott optical glass catalogue for N-BK7.
* Ghatak & Thyagarajan figs. 6.1-6.6, which the notebooks
  `1-Refractive-Index.ipynb` and `2-Materials.ipynb` reproduce.

The derivative tests do not re-use ofiber's own derivative code.  They
difference an independent Sellmeier evaluation written below, so an error
in `_d_sellmeier` or `_d2_sellmeier` cannot cancel out of the comparison.
"""

import numpy as np
import pytest
import ofiber
from ofiber.refraction import _GLASS

# glass table indices used throughout
SIO2 = 0
GEO2 = 1
BK7 = 13

# wavelengths where the derivatives are far from zero, so a relative
# tolerance is meaningful
DERIV_CASES = [
    (SIO2, 500e-9),
    (SIO2, 850e-9),
    (SIO2, 1550e-9),
    (GEO2, 850e-9),
    (GEO2, 1310e-9),
    (BK7, 587.6e-9),
    (BK7, 1550e-9),
]


def _n_reference(coef, lam):
    """
    Evaluate the Sellmeier equation without using ofiber.

    Args:
        coef: six Sellmeier coefficients, c[0:3] in microns**2 then b[3:6]
        lam: wavelength in vacuum [m]

    Returns:
        index of refraction [-]
    """
    c = np.asarray(coef[0:3], dtype=float)
    b = np.asarray(coef[3:6], dtype=float)
    lam2 = (np.asarray(lam, dtype=float) * 1e6) ** 2
    nsq = 1.0
    for i in range(3):
        nsq = nsq + b[i] * lam2 / (lam2 - c[i])
    return np.sqrt(nsq)


def _fd1(coef, lam, rel_h=1e-4):
    """
    Differentiate `_n_reference` once with a five-point stencil.

    Args:
        coef: six Sellmeier coefficients
        lam: wavelength in vacuum [m]
        rel_h: step size as a fraction of lam

    Returns:
        dn/dlambda [1/m]
    """
    h = lam * rel_h

    def f(k):
        return _n_reference(coef, lam + k * h)

    return (-f(2) + 8 * f(1) - 8 * f(-1) + f(-2)) / (12 * h)


def _fd2(coef, lam, rel_h=1e-3):
    """
    Differentiate `_n_reference` twice with a five-point stencil.

    A relatively coarse step is deliberate: the second difference loses
    roughly 2/3 of the available digits to cancellation, so shrinking h
    makes the estimate worse rather than better.

    Args:
        coef: six Sellmeier coefficients
        lam: wavelength in vacuum [m]
        rel_h: step size as a fraction of lam

    Returns:
        d2n/dlambda2 [1/m**2]
    """
    h = lam * rel_h

    def f(k):
        return _n_reference(coef, lam + k * h)

    return (-f(2) + 16 * f(1) - 30 * f(0) + 16 * f(-1) - f(-2)) / (12 * h**2)


# --------------------------------------------------------------------------
# index of refraction
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "lam, expected",
    [
        (587.6e-9, 1.458464),  # Malitson, sodium d line
        (632.8e-9, 1.457018),  # Malitson, HeNe
        (1550e-9, 1.444022),  # Malitson, telecom C band
    ],
)
def test_n_fused_silica_matches_malitson(lam, expected):
    """Check fused silica against the published Malitson indices."""
    assert ofiber.n(ofiber.glass(SIO2), lam) == pytest.approx(expected, abs=5e-6)


def test_n_bk7_matches_schott_catalogue():
    """Check N-BK7 against the Schott catalogue value of n_d."""
    assert ofiber.n(ofiber.glass(BK7), 587.6e-9) == pytest.approx(1.51680, abs=5e-6)


def test_n_matches_independent_sellmeier():
    """Verify n() against a Sellmeier evaluation written in this file."""
    lam = np.linspace(500e-9, 1700e-9, 41)
    for i in range(len(_GLASS)):
        coef = ofiber.glass(i)
        assert np.allclose(ofiber.n(coef, lam), _n_reference(coef, lam), rtol=1e-12)


def test_n_accepts_arrays_and_scalars_alike():
    """Verify that array input gives the same answer as looping over scalars."""
    coef = ofiber.glass(SIO2)
    lam = np.linspace(600e-9, 1600e-9, 17)
    expected = np.array([ofiber.n(coef, single) for single in lam])
    assert np.allclose(ofiber.n(coef, lam), expected, rtol=1e-14)


def test_n_decreases_with_wavelength():
    """Verify normal dispersion for silica across the near infrared."""
    lam = np.linspace(600e-9, 1600e-9, 101)
    index = ofiber.n(ofiber.glass(SIO2), lam)
    assert np.all(np.diff(index) < 0)


def test_every_glass_has_a_physical_index():
    """Verify all tabulated glasses give a sane index in the visible and at 1550 nm."""
    for i in range(len(_GLASS)):
        for lam in (587.6e-9, 1550e-9):
            index = ofiber.n(ofiber.glass(i), lam)
            assert np.isfinite(index), "%s is not finite at %.1f nm" % (ofiber.glass_name(i), lam * 1e9)
            assert 1.0 < index < 2.5, "%s has index %f at %.1f nm" % (ofiber.glass_name(i), index, lam * 1e9)


# --------------------------------------------------------------------------
# first derivative
# --------------------------------------------------------------------------


@pytest.mark.parametrize("i, lam", DERIV_CASES)
def test_dn_matches_numerical_derivative(i, lam):
    """Verify dn() against a finite difference of the independent Sellmeier."""
    coef = ofiber.glass(i)
    assert ofiber.dn(coef, lam) == pytest.approx(_fd1(coef, lam), rel=1e-6)


def test_dn_is_negative_for_normal_dispersion():
    """Verify dn/dlambda stays negative for silica in the near infrared."""
    lam = np.linspace(600e-9, 1600e-9, 51)
    assert np.all(ofiber.dn(ofiber.glass(SIO2), lam) < 0)


def test_dn_accepts_arrays():
    """Verify array input to dn() matches scalar evaluation."""
    coef = ofiber.glass(SIO2)
    lam = np.linspace(1000e-9, 1700e-9, 13)
    expected = np.array([ofiber.dn(coef, single) for single in lam])
    assert np.allclose(ofiber.dn(coef, lam), expected, rtol=1e-14)


# --------------------------------------------------------------------------
# second derivative
# --------------------------------------------------------------------------


@pytest.mark.parametrize("i, lam", DERIV_CASES)
def test_d2n_matches_numerical_derivative(i, lam):
    """Verify d2n() against a finite difference of the independent Sellmeier."""
    coef = ofiber.glass(i)
    assert ofiber.d2n(coef, lam) == pytest.approx(_fd2(coef, lam), rel=1e-3)


def test_d2n_sums_all_three_sellmeier_terms():
    """
    Verify d2n() accumulates every Sellmeier term.

    A single dropped term in the first-derivative sum inside `_d2_sellmeier`
    shifts d2n by about half a percent, and by several percent near the zero
    crossing, so this compares against a deliberately term-by-term reference.
    """
    for i in (SIO2, GEO2, BK7):
        coef = ofiber.glass(i)
        c = np.asarray(coef[0:3], dtype=float)
        b = np.asarray(coef[3:6], dtype=float)
        for lam in (850e-9, 1310e-9, 1550e-9):
            lam2 = (lam * 1e6) ** 2
            index = _n_reference(coef, lam)
            first = sum(b[k] * c[k] / (lam2 - c[k]) ** 2 for k in range(3))
            second = sum(b[k] * c[k] * (3 * lam2 + c[k]) / (lam2 - c[k]) ** 3 for k in range(3))
            expected = (second / index - lam2 * first**2 / index**3) * 1e12
            assert ofiber.d2n(coef, lam) == pytest.approx(expected, rel=1e-10)


def test_d2n_accepts_arrays():
    """Verify array input to d2n() matches scalar evaluation."""
    coef = ofiber.glass(SIO2)
    lam = np.linspace(1000e-9, 1700e-9, 13)
    expected = np.array([ofiber.d2n(coef, single) for single in lam])
    assert np.allclose(ofiber.d2n(coef, lam), expected, rtol=1e-14)


def test_silica_material_zero_dispersion_wavelength():
    """
    Verify the wavelength where silica has zero material dispersion.

    Ghatak quotes roughly 1270 nm.  This pins the crossing tightly because it
    is the most sensitive consumer of d2n() in the package: dropping a term
    from the Sellmeier sum moved this crossing to 1274.5 nm.
    """
    lam = np.linspace(1.20e-6, 1.35e-6, 15001)
    d2 = ofiber.d2n(ofiber.glass(SIO2), lam)
    crossings = np.where(np.diff(np.sign(d2)) != 0)[0]
    assert len(crossings) == 1
    assert lam[crossings[0]] * 1e9 == pytest.approx(1272.75, abs=0.5)


def test_group_index_minimum_coincides_with_zero_d2n():
    """
    Verify d(n_group)/dlambda vanishes exactly where d2n does.

    Since n_group = n - lambda*dn/dlambda, its derivative is
    -lambda*d2n/dlambda2, so the two features must sit at the same
    wavelength.  This ties n(), dn() and d2n() together.
    """
    lam = np.linspace(1.20e-6, 1.35e-6, 15001)
    coef = ofiber.glass(SIO2)
    ng_min = lam[np.argmin(ofiber.n_group(coef, lam))]
    d2 = ofiber.d2n(coef, lam)
    zero = lam[np.where(np.diff(np.sign(d2)) != 0)[0][0]]
    assert ng_min == pytest.approx(zero, abs=2 * (lam[1] - lam[0]))


# --------------------------------------------------------------------------
# group index
# --------------------------------------------------------------------------


def test_n_group_follows_its_definition():
    """Verify n_group equals n - lambda*dn/dlambda."""
    coef = ofiber.glass(SIO2)
    lam = np.linspace(600e-9, 1600e-9, 21)
    expected = ofiber.n(coef, lam) - lam * ofiber.dn(coef, lam)
    assert np.allclose(ofiber.n_group(coef, lam), expected, rtol=1e-14)


def test_n_group_exceeds_phase_index():
    """Verify the group index is larger than the phase index in normal dispersion."""
    coef = ofiber.glass(SIO2)
    lam = np.linspace(600e-9, 1600e-9, 21)
    assert np.all(ofiber.n_group(coef, lam) > ofiber.n(coef, lam))


def test_n_group_of_silica_at_1550nm():
    """Check the silica group index at 1550 nm against the accepted 1.4626."""
    assert ofiber.n_group(ofiber.glass(SIO2), 1550e-9) == pytest.approx(1.4626, abs=1e-3)


# --------------------------------------------------------------------------
# glass table
# --------------------------------------------------------------------------


def test_glass_table_and_name_list_are_the_same_length():
    """Verify every row of Sellmeier coefficients has exactly one name."""
    assert len(_GLASS) == len(ofiber.ALL_GLASS_NAMES)


def test_every_glass_row_has_six_finite_coefficients():
    """Verify the shape and finiteness of every row in the glass table."""
    for i in range(len(_GLASS)):
        coef = ofiber.glass(i)
        assert len(coef) == 6, "%s has %d coefficients" % (ofiber.glass_name(i), len(coef))
        assert np.all(np.isfinite(coef)), "%s has a non-finite coefficient" % ofiber.glass_name(i)


@pytest.mark.parametrize("i, expected", [(SIO2, "SiO2"), (GEO2, "GeO2"), (BK7, "N-BK7")])
def test_glass_name_of_known_indices(i, expected):
    """Check the names of the glass indices the notebooks rely on."""
    assert ofiber.glass_name(i) == expected


def test_glass_returns_the_matching_row():
    """Verify glass() hands back the coefficient row for that index."""
    assert np.allclose(ofiber.glass(BK7), _GLASS[BK7])


# --------------------------------------------------------------------------
# find_glass
# --------------------------------------------------------------------------


@pytest.mark.parametrize("name, expected", [("SiO2", SIO2), ("GeO2", GEO2), ("BK7", BK7)])
def test_find_glass_locates_known_names(name, expected):
    """Check lookups used by 1-Refractive-Index.ipynb."""
    assert ofiber.find_glass(name) == expected


@pytest.mark.parametrize("name", ["sio2", "SIO2", "SiO2"])
def test_find_glass_ignores_case(name):
    """Verify the documented case-insensitive matching."""
    assert ofiber.find_glass(name) == SIO2


@pytest.mark.parametrize("name", ["SF1", "SF2", "SF4", "SF5", "SF10", "SF11", "F2", "K7", "N-BK7"])
def test_find_glass_prefers_an_exact_name_over_a_substring(name):
    """Verify an exact name wins over a longer name that contains it."""
    assert ofiber.glass_name(ofiber.find_glass(name)) == name


def test_find_glass_still_matches_on_a_substring():
    """Verify the documented fallback, which is how the notebooks find N-BK7."""
    assert ofiber.find_glass("BK7") == BK7


def test_find_glass_rejects_an_unknown_name():
    """Verify an unknown glass name raises rather than returning a valid index."""
    with pytest.raises(ValueError, match="not a known glass"):
        ofiber.find_glass("not-a-real-glass")


def test_glass_names_are_unique():
    """Verify every glass is reachable by an unambiguous name."""
    names = [str(name) for name in ofiber.ALL_GLASS_NAMES]
    duplicates = sorted({name for name in names if names.count(name) > 1})
    assert duplicates == []


def test_every_glass_name_round_trips_through_find_glass():
    """Verify find_glass recovers the index of every name in the table."""
    for i, name in enumerate(ofiber.ALL_GLASS_NAMES):
        found = ofiber.find_glass(str(name))
        assert found == i, "%s resolved to %s" % (name, ofiber.glass_name(found))


# --------------------------------------------------------------------------
# doped glass
# --------------------------------------------------------------------------


@pytest.mark.parametrize("x, i", [(0.0, SIO2), (1.0, GEO2)])
def test_doped_glass_endpoints_match_the_pure_glasses(x, i):
    """Verify the interpolation endpoints reproduce pure SiO2 and pure GeO2."""
    assert np.allclose(ofiber.doped_glass(x), ofiber.glass(i), rtol=1e-6)


def test_doped_glass_returns_six_coefficients():
    """Verify doped_glass produces a usable coefficient row."""
    assert len(ofiber.doped_glass(0.1)) == 6


def test_doped_glass_index_rises_with_germania():
    """Verify adding GeO2 raises the index, which is why it is the core dopant."""
    lam = 1550e-9
    fractions = np.array([0.0, 0.063, 0.193, 0.5, 1.0])  # Ghatak fig 6.4 uses the first three
    index = np.array([ofiber.n(ofiber.doped_glass(x), lam) for x in fractions])
    assert np.all(np.diff(index) > 0)


@pytest.mark.parametrize(
    "x, expected",
    [
        (0, r"SiO$_2$"),
        (1, r"GeO$_2$"),
        (0.1, r"0.10 GeO$_2$ : 0.90 SiO$_2$"),
    ],
)
def test_doped_glass_name(x, expected):
    """Check the labels the notebooks put in plot legends."""
    assert ofiber.doped_glass_name(x) == expected


# --------------------------------------------------------------------------
# air
# --------------------------------------------------------------------------


def test_n_air_in_the_visible():
    """Check air against the accepted 1.000277 at 550 nm and 15 C."""
    assert ofiber.n_air(550e-9) == pytest.approx(1.0002771, abs=2e-6)


def test_n_air_decreases_with_wavelength():
    """Verify air disperses normally."""
    lam = np.linspace(400e-9, 1600e-9, 25)
    assert np.all(np.diff([ofiber.n_air(single) for single in lam]) < 0)


def test_n_air_decreases_with_temperature():
    """Verify warmer air is less dense and so less refractive."""
    assert ofiber.n_air(550e-9, 30) < ofiber.n_air(550e-9, 15)


def test_n_air_is_close_to_unity():
    """Verify air stays just above vacuum across the tabulated range."""
    lam = np.linspace(400e-9, 1600e-9, 25)
    for single in lam:
        assert 1.0 < ofiber.n_air(single) < 1.001
