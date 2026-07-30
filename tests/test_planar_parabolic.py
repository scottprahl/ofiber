# pylint: disable=invalid-name
"""
Tests for the ofiber.planar_parabolic module.

No notebook exercises this module, so the reference values are closed forms
from section 7.9 of Ghatak & Thyagarajan, *An Introduction to Fiber Optics*.

A parabolic index turns the wave equation into a harmonic oscillator, so the
results are the textbook Hermite-Gaussian ones:

* gamma = sqrt(V)/a, which also equals sqrt(k0*n1*sqrt(2*Delta)/a)
* beta**2 = k**2 - (2m+1)*gamma**2, real while 2m+1 < (k*a)**2/V
* the fields are orthonormal, even or odd with m, and have exactly m zeros

The orthonormality check is the strongest of these: it exercises the Hermite
polynomial and the normalisation constant together, and neither can be wrong
on its own without the Gram matrix ceasing to be the identity.
"""

import numpy as np
import pytest
import ofiber
from ofiber.planar_parabolic import _herm

LAMBDA0 = 1.55e-6
N1 = 1.46
HALF_THICKNESS = 5e-6  # m
DELTA = 0.01

NA = N1 * np.sqrt(2 * DELTA)
V_PARAM = 2 * np.pi / LAMBDA0 * HALF_THICKNESS * NA
K_AXIS = 2 * np.pi * N1 / LAMBDA0
GAMMA = np.sqrt(V_PARAM) / HALF_THICKNESS


def _field(m, x):
    """
    Evaluate the mode field with the fixed waveguide used throughout.

    Args:
        m: mode number      [-]
        x: position         [m]

    Returns:
        E_y(x) [m**-0.5]
    """
    return ofiber.TE_planar_parabolic_field(m, LAMBDA0, N1, DELTA, HALF_THICKNESS, x)


# --------------------------------------------------------------------------
# Hermite polynomials
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "n, closed_form",
    [
        (0, np.ones_like),
        (1, lambda x: 2 * x),
        (2, lambda x: 4 * x**2 - 2),
        (3, lambda x: 8 * x**3 - 12 * x),
        (4, lambda x: 16 * x**4 - 48 * x**2 + 12),
    ],
)
def test_hermite_polynomials_match_their_closed_forms(n, closed_form):
    """Check the physicists' Hermite polynomials against their explicit forms."""
    x = np.linspace(-3, 3, 61)
    assert np.allclose(_herm(n, x), closed_form(x))


def test_hermite_polynomials_have_the_parity_of_their_order():
    """Verify H_n(-x) = (-1)**n H_n(x)."""
    x = np.linspace(0.1, 3, 30)
    for n in range(6):
        assert np.allclose(_herm(n, -x), (-1) ** n * _herm(n, x))


# --------------------------------------------------------------------------
# propagation constants
# --------------------------------------------------------------------------


def test_gamma_follows_from_the_waveguide_parameters():
    """Verify sqrt(V)/a is the harmonic oscillator constant of the index profile."""
    from_first_principles = np.sqrt(2 * np.pi / LAMBDA0 * N1 * np.sqrt(2 * DELTA) / HALF_THICKNESS)
    assert GAMMA == pytest.approx(from_first_principles, rel=1e-12)


@pytest.mark.parametrize("m", [0, 1, 2, 5, 20])
def test_propagation_constant_satisfies_the_eigenvalue_relation(m):
    """Verify beta**2 = k**2 - (2m+1) gamma**2."""
    beta = ofiber.parabolic_propagation_constant(m, LAMBDA0, N1, HALF_THICKNESS, V_PARAM)
    assert beta**2 == pytest.approx(K_AXIS**2 - (2 * m + 1) * GAMMA**2, rel=1e-12)


def test_propagation_constant_falls_with_mode_number():
    """Verify higher order modes are less tightly bound."""
    beta = np.array([ofiber.parabolic_propagation_constant(m, LAMBDA0, N1, HALF_THICKNESS, V_PARAM) for m in range(30)])
    assert np.all(np.diff(beta) < 0)


def test_fundamental_propagation_constant_is_just_below_the_axial_wavenumber():
    """Verify no mode travels faster than a plane wave on the waveguide axis."""
    beta = ofiber.parabolic_propagation_constant(0, LAMBDA0, N1, HALF_THICKNESS, V_PARAM)
    assert 0 < beta < K_AXIS
    assert beta == pytest.approx(K_AXIS, rel=1e-2)


def test_propagation_constant_is_nan_beyond_cutoff():
    """Verify an unguided mode number gives nan rather than a complex result."""
    unguided = int(((K_AXIS / GAMMA) ** 2 - 1) / 2) + 5
    with np.errstate(invalid="ignore"):
        beta = ofiber.parabolic_propagation_constant(unguided, LAMBDA0, N1, HALF_THICKNESS, V_PARAM)
    assert np.isnan(beta)


def test_mode_count_matches_the_guided_condition():
    """Verify one propagation constant is returned per mode with 2m+1 < (k a)**2/V."""
    betas = ofiber.parabolic_propagation_constants(LAMBDA0, N1, HALF_THICKNESS, V_PARAM)
    expected = int(np.floor(((K_AXIS / GAMMA) ** 2 - 1) / 2)) + 1
    assert len(betas) == expected


def test_every_returned_propagation_constant_is_guided():
    """Verify all the betas are real, positive and below the axial wavenumber."""
    betas = ofiber.parabolic_propagation_constants(LAMBDA0, N1, HALF_THICKNESS, V_PARAM)
    assert np.all(np.isfinite(betas))
    assert np.all(betas > 0)
    assert np.all(betas < K_AXIS)
    assert np.all(np.diff(betas) < 0)


def test_one_mode_past_the_last_is_not_guided():
    """Verify the mode count stops exactly where beta stops being real."""
    betas = ofiber.parabolic_propagation_constants(LAMBDA0, N1, HALF_THICKNESS, V_PARAM)
    assert K_AXIS**2 - (2 * len(betas) + 1) * GAMMA**2 < 0


def test_propagation_constants_agree_with_the_single_mode_routine():
    """Verify the array routine reproduces the scalar one term by term."""
    betas = ofiber.parabolic_propagation_constants(LAMBDA0, N1, HALF_THICKNESS, V_PARAM)
    one_by_one = [ofiber.parabolic_propagation_constant(m, LAMBDA0, N1, HALF_THICKNESS, V_PARAM) for m in range(6)]
    assert np.allclose(betas[:6], one_by_one)


def test_a_waveguide_guiding_nothing_returns_an_empty_array():
    """
    Verify the degenerate case gives an empty array rather than an error.

    This needs V above (k a)**2, which no real waveguide has, but the guard
    keeps np.empty from being handed a negative length.
    """
    huge_V = 10 * (K_AXIS * HALF_THICKNESS) ** 2
    assert len(ofiber.parabolic_propagation_constants(LAMBDA0, N1, HALF_THICKNESS, huge_V)) == 0


# --------------------------------------------------------------------------
# mode fields
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", [0, 1, 2, 3, 4, 5])
def test_mode_fields_are_normalized(m):
    """Verify the documented normalization, that the integral of E**2 is one."""
    x = np.linspace(-40e-6, 40e-6, 200001)
    assert np.trapezoid(_field(m, x) ** 2, x) == pytest.approx(1.0, rel=1e-6)


def test_mode_fields_are_orthogonal():
    """
    Verify distinct modes are orthogonal, so the Gram matrix is the identity.

    This pins the Hermite polynomial and the normalisation constant at once.
    """
    x = np.linspace(-40e-6, 40e-6, 200001)
    fields = [_field(m, x) for m in range(5)]
    for i, first in enumerate(fields):
        for j, second in enumerate(fields):
            overlap = np.trapezoid(first * second, x)
            assert overlap == pytest.approx(1.0 if i == j else 0.0, abs=1e-6)


def test_fundamental_mode_is_a_gaussian():
    """Verify the m=0 field is the Gaussian exp(-(gamma x)**2/2)."""
    x = np.linspace(-20e-6, 20e-6, 501)
    field = _field(0, x)
    assert np.allclose(field / field[len(x) // 2], np.exp(-0.5 * (GAMMA * x) ** 2), rtol=1e-9)


@pytest.mark.parametrize("m", [0, 1, 2, 3, 4, 5])
def test_mode_fields_have_the_parity_of_their_order(m):
    """Verify even modes are symmetric and odd modes antisymmetric."""
    x = np.linspace(1e-7, 20e-6, 200)
    assert np.allclose(_field(m, -x), (-1) ** m * _field(m, x))


@pytest.mark.parametrize("m", [0, 1, 2, 3, 4, 5])
def test_mode_m_has_m_zero_crossings(m):
    """
    Verify the mth mode has exactly m nodes, as a Hermite-Gaussian must.

    Samples that land exactly on a node are dropped first: an odd mode is
    exactly zero at x=0, and np.sign would report that single node as two
    changes of sign rather than one.
    """
    x = np.linspace(-25e-6, 25e-6, 20001)
    field = _field(m, x)
    assert np.count_nonzero(np.diff(np.sign(field[field != 0]))) == m


def test_mode_fields_decay_far_from_the_axis():
    """Verify the Gaussian envelope confines every mode."""
    x = np.linspace(-25e-6, 25e-6, 2001)
    for m in range(5):
        peak = np.max(np.abs(_field(m, x)))
        assert abs(_field(m, 60e-6)) < 1e-9 * peak


def test_field_accepts_a_scalar_position():
    """Verify a scalar position gives a scalar field."""
    assert np.isscalar(_field(0, 0.0)) or _field(0, 0.0).ndim == 0


def test_high_order_fields_overflow_above_about_m_155():
    """
    Document the largest mode number the normalization can express.

    Nm divides by 2**m * m!, which exceeds the range of a float once m passes
    roughly 155, and the Hermite polynomial itself overflows near m of 300.
    Waveguides needing modes that high are outside what this routine handles.
    """
    assert np.isfinite(_field(150, 0.0))
    with pytest.raises(OverflowError):
        _field(160, 0.0)


# --------------------------------------------------------------------------
# module surface
# --------------------------------------------------------------------------


def test_everything_in_all_is_reachable_from_the_package():
    """Verify each exported name survives the star import in __init__.py."""
    for name in ofiber.planar_parabolic.__all__:
        assert hasattr(ofiber, name), "%s is in __all__ but missing from ofiber" % name
