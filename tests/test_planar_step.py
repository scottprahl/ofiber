# pylint: disable=invalid-name
"""
Tests for the ofiber.planar_step module.

No notebook exercises this module directly, so the references are the closed
forms of chapter 7 of Ghatak & Thyagarajan, *An Introduction to Fiber Optics*.

For a symmetric three layer slab with xi = kappa*d/2 the eigenvalues satisfy

    xi tan(xi)  =              sqrt(V**2/4 - xi**2)   even TE modes
   -xi cot(xi)  =              sqrt(V**2/4 - xi**2)   odd TE modes
    xi tan(xi)  = (n1/n2)**2 * sqrt(V**2/4 - xi**2)   even TM modes
   -xi cot(xi)  = (n1/n2)**2 * sqrt(V**2/4 - xi**2)   odd TM modes

so mode m lives in (m*pi/2, (m+1)*pi/2) and a slab guides floor(V/pi)+1 modes.

Setting n1 = n2 makes the TM factor unity, so the TM eigenvalues must collapse
onto the TE ones exactly.  That is the strongest check here: it ties the two
solvers together without appealing to any tabulated number.
"""

import matplotlib
import numpy as np
import pytest
import ofiber
from ofiber.planar_step import _TE_mode, _TM_mode

matplotlib.use("Agg")

N1 = 1.5
N2 = 1.0
THICKNESS = 2e-6  # m


def _mode_count(V):
    """
    Return the number of guided modes of a symmetric slab.

    Args:
        V: the V-parameter of the waveguide [-]

    Returns:
        number of guided modes [-]
    """
    return int(V / np.pi) + 1


# --------------------------------------------------------------------------
# eigenvalue equations
# --------------------------------------------------------------------------


@pytest.mark.parametrize("V", [2.0, 5.0, 8.0, 12.0, 20.0])
def test_te_eigenvalues_solve_the_eigenvalue_equation(V):
    """Verify every TE eigenvalue is a root of its own eigenvalue equation."""
    for mode, xi in enumerate(ofiber.TE_crossings(V)):
        assert _TE_mode(xi, V, mode) == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize("V", [2.0, 5.0, 8.0, 12.0, 20.0])
def test_tm_eigenvalues_solve_the_eigenvalue_equation(V):
    """Verify every TM eigenvalue is a root of its own eigenvalue equation."""
    for mode, xi in enumerate(ofiber.TM_crossings(V, N1, N2)):
        assert _TM_mode(xi, V, N1, N2, mode) == pytest.approx(0.0, abs=1e-9)


@pytest.mark.parametrize("V", [2.0, 8.0, 20.0])
def test_tm_collapses_onto_te_when_the_indices_are_equal(V):
    """Verify the (n1/n2)**2 factor is the only thing separating TM from TE."""
    assert np.allclose(ofiber.TM_crossings(V, N1, N1), ofiber.TE_crossings(V), rtol=1e-12)


@pytest.mark.parametrize("V", [4.0, 8.0, 16.0])
def test_tm_modes_are_less_confined_than_te_modes(V):
    """Verify each TM eigenvalue exceeds its TE counterpart when n1 > n2."""
    assert np.all(ofiber.TM_crossings(V, N1, N2) > ofiber.TE_crossings(V))


@pytest.mark.parametrize("V", [1.0, 5.0, 12.0, 20.0])
def test_each_eigenvalue_lies_in_its_own_half_period(V):
    """Verify mode m has an eigenvalue between m*pi/2 and (m+1)*pi/2."""
    for mode, xi in enumerate(ofiber.TE_crossings(V)):
        assert mode * np.pi / 2 < xi < (mode + 1) * np.pi / 2


@pytest.mark.parametrize("V", [1.0, 5.0, 12.0, 20.0])
def test_eigenvalues_increase_with_mode_number(V):
    """Verify higher order modes have larger eigenvalues."""
    assert np.all(np.diff(ofiber.TE_crossings(V)) > 0)
    assert np.all(np.diff(ofiber.TM_crossings(V, N1, N2)) > 0)


@pytest.mark.parametrize("V", [0.5, 3.0, 3.15, 6.28, 6.3, 12.0, 20.0])
def test_mode_count_is_one_more_than_V_over_pi(V):
    """Verify the slab guides floor(V/pi)+1 modes, TE and TM alike."""
    assert len(ofiber.TE_crossings(V)) == _mode_count(V)
    assert len(ofiber.TM_crossings(V, N1, N2)) == _mode_count(V)


def test_every_eigenvalue_is_below_half_V():
    """Verify no eigenvalue exceeds V/2, where the cladding decay would vanish."""
    for V in (2.0, 7.0, 15.0):
        assert np.all(ofiber.TE_crossings(V) < V / 2)


def test_a_mode_above_cutoff_reports_no_eigenvalue():
    """Verify a mode number the slab cannot guide gives the zero sentinel."""
    assert ofiber.TE_crossing(5.0, 2) == 0
    assert ofiber.TM_crossing(5.0, N1, N2, 2) == 0


def test_a_mode_just_under_cutoff_reports_no_eigenvalue():
    """Verify a mode whose bracket starts past V/2 gives the zero sentinel."""
    assert ofiber.TE_crossing(np.pi + 1e-6, 1) == 0
    assert ofiber.TM_crossing(np.pi + 1e-6, N1, N1, 1) == 0


def test_a_degenerate_bracket_reports_no_eigenvalue():
    """
    Verify a bracket too narrow to hold a root is handled rather than raising.

    The search for mode 1 starts at pi/2 + 1e-5, so V = pi + 2e-5 makes the
    bracket collapse to a single point and brentq raises; the routine falls
    back to the zero sentinel.
    """
    V = np.pi + 2e-5
    assert ofiber.TE_crossing(V, 1) == 0
    assert ofiber.TM_crossing(V, N1, N1, 1) == 0


# --------------------------------------------------------------------------
# propagation constants
# --------------------------------------------------------------------------


def test_normalized_propagation_constant_is_between_zero_and_one():
    """Verify b stays between the cladding and core indices."""
    V = np.linspace(1.0, 20.0, 40)
    te = ofiber.TE_propagation_constant(V, 0)
    tm = ofiber.TM_propagation_constant(V, N1, N2, 0)
    assert np.all((0 < te) & (te < 1))
    assert np.all((0 < tm) & (tm < 1))


def test_propagation_constant_matches_its_eigenvalue():
    """Verify b = 1 - (2 xi / V)**2 for each solved eigenvalue."""
    V = np.linspace(2.0, 20.0, 25)
    b = ofiber.TE_propagation_constant(V, 0)
    expected = [1 - (2 * ofiber.TE_crossing(v, 0) / v) ** 2 for v in V]
    assert np.allclose(b, expected)


def test_propagation_constant_rises_with_V():
    """Verify a larger V confines the fundamental mode more tightly."""
    V = np.linspace(1.0, 20.0, 60)
    assert np.all(np.diff(ofiber.TE_propagation_constant(V, 0)) > 0)


def test_higher_modes_have_smaller_propagation_constants():
    """Verify b falls with mode number at fixed V."""
    V = np.array([20.0])
    b = [ofiber.TE_propagation_constant(V, mode)[0] for mode in range(_mode_count(20.0))]
    assert np.all(np.diff(b) < 0)


def test_propagation_constant_is_zero_for_an_unguided_mode():
    """Verify a mode the slab cannot guide reports zero rather than a number."""
    assert ofiber.TE_propagation_constant(np.array([5.0]), 2)[0] == 0
    assert ofiber.TM_propagation_constant(np.array([5.0]), N1, N2, 2)[0] == 0


# --------------------------------------------------------------------------
# fields
# --------------------------------------------------------------------------


@pytest.mark.parametrize("mode", [0, 1, 2, 3])
def test_te_field_is_continuous_at_the_core_boundary(mode):
    """Verify the sinusoidal and evanescent branches meet at x = d/2."""
    V = 12.0
    edge = THICKNESS / 2
    inside = ofiber.TE_field(V, THICKNESS, np.array([edge * (1 - 1e-9)]), mode)[0]
    outside = ofiber.TE_field(V, THICKNESS, np.array([edge * (1 + 1e-9)]), mode)[0]
    assert inside == pytest.approx(outside, abs=1e-8)


@pytest.mark.parametrize("mode", [0, 1, 2, 3])
def test_tm_field_is_continuous_at_the_core_boundary(mode):
    """Verify the TM field matches across the boundary as well."""
    V = 12.0
    edge = THICKNESS / 2
    inside = ofiber.TM_field(V, N1, N2, THICKNESS, np.array([edge * (1 - 1e-9)]), mode)[0]
    outside = ofiber.TM_field(V, N1, N2, THICKNESS, np.array([edge * (1 + 1e-9)]), mode)[0]
    assert inside == pytest.approx(outside, abs=1e-8)


@pytest.mark.parametrize("mode", [0, 1, 2, 3])
def test_field_has_the_parity_of_its_mode_number(mode):
    """Verify even modes are symmetric about the axis and odd modes antisymmetric."""
    x = np.linspace(1e-9, 3e-6, 200)
    field = ofiber.TE_field(12.0, THICKNESS, x, mode)
    mirrored = ofiber.TE_field(12.0, THICKNESS, -x, mode)
    assert np.allclose(mirrored, (-1) ** mode * field)


@pytest.mark.parametrize("mode", [0, 1, 2, 3])
def test_field_of_mode_m_has_m_nodes(mode):
    """
    Verify mode m crosses zero exactly m times inside the slab.

    Samples landing exactly on a node are dropped, since an odd mode is
    exactly zero at x=0 and np.sign would count that as two changes.
    """
    x = np.linspace(-THICKNESS / 2 * 0.999, THICKNESS / 2 * 0.999, 4001)
    field = ofiber.TE_field(12.0, THICKNESS, x, mode)
    assert np.count_nonzero(np.diff(np.sign(field[field != 0]))) == mode


@pytest.mark.parametrize("mode", [0, 1])
def test_field_decays_outside_the_slab(mode):
    """Verify the cladding field is evanescent."""
    x = np.linspace(THICKNESS / 2, 4e-6, 100)
    field = np.abs(ofiber.TE_field(12.0, THICKNESS, x, mode))
    assert np.all(np.diff(field) < 0)
    assert field[-1] < 1e-3 * field[0]


def test_te_and_tm_fields_agree_when_the_indices_are_equal():
    """Verify the two field routines differ only through their eigenvalues."""
    x = np.linspace(-3e-6, 3e-6, 101)
    te = ofiber.TE_field(8.0, THICKNESS, x, 0)
    tm = ofiber.TM_field(8.0, N1, N1, THICKNESS, x, 0)
    assert np.allclose(te, tm)


# --------------------------------------------------------------------------
# plots
# --------------------------------------------------------------------------


@pytest.mark.parametrize("V", [2.0, 12.0])
def test_te_mode_plot_builds(V):
    """Verify the TE eigenvalue plot runs for a slab with one and with several modes."""
    plot = ofiber.TE_mode_plot(V)
    assert plot is not None
    plot.close("all")


@pytest.mark.parametrize("V", [2.0, 12.0])
def test_tm_mode_plot_builds(V):
    """Verify the TM eigenvalue plot runs."""
    plot = ofiber.TM_mode_plot(V, N1, N2)
    assert plot is not None
    plot.close("all")


# --------------------------------------------------------------------------
# module surface
# --------------------------------------------------------------------------


def test_everything_in_all_is_reachable_from_the_package():
    """Verify each exported name survives the star import in __init__.py."""
    for name in ofiber.planar_step.__all__:
        assert hasattr(ofiber, name), "%s is in __all__ but missing from ofiber" % name


# --------------------------------------------------------------------------
# scalar and integer arguments
# --------------------------------------------------------------------------


def test_integer_array_of_V_gives_float_results():
    """Verify an integer array of V is not truncated to zeros."""
    assert np.allclose(
        ofiber.TE_propagation_constant(np.array([5, 8]), 0),
        ofiber.TE_propagation_constant(np.array([5.0, 8.0]), 0),
    )
    assert np.allclose(
        ofiber.TM_propagation_constant(np.array([5, 8]), N1, N2, 0),
        ofiber.TM_propagation_constant(np.array([5.0, 8.0]), N1, N2, 0),
    )


def test_propagation_constant_accepts_a_scalar_V():
    """Verify a scalar V works, as it does everywhere else in the package."""
    assert ofiber.TE_propagation_constant(3.0, 0) == pytest.approx(0.62801672, rel=1e-6)
    assert ofiber.TM_propagation_constant(3.0, N1, N2, 0) == pytest.approx(0.44908264, rel=1e-6)
    assert ofiber.TE_propagation_constant(5.0, 0) == pytest.approx(0.80268263, rel=1e-6)
    assert ofiber.TM_propagation_constant(5.0, N1, N2, 0) == pytest.approx(0.72744261, rel=1e-6)


def test_scalar_and_array_propagation_constants_agree():
    """Verify the scalar branch reproduces the array branch term by term."""
    V = np.linspace(2.0, 18.0, 17)
    assert np.allclose(ofiber.TE_propagation_constant(V, 0), [ofiber.TE_propagation_constant(v, 0) for v in V])
    assert np.allclose(
        ofiber.TM_propagation_constant(V, N1, N2, 1),
        [ofiber.TM_propagation_constant(v, N1, N2, 1) for v in V],
    )


def test_scalar_propagation_constant_is_zero_for_an_unguided_mode():
    """Verify the scalar branch reports an unguided mode the same way as the array one."""
    assert ofiber.TE_propagation_constant(5.0, 2) == 0
    assert ofiber.TM_propagation_constant(5.0, N1, N2, 2) == 0


# --------------------------------------------------------------------------
# unguided modes
# --------------------------------------------------------------------------


def test_an_unguided_mode_does_not_produce_a_field():
    """
    Verify a mode the slab cannot guide gives nan rather than a plausible field.

    A zero eigenvalue would otherwise give kappa=0, so cos(kappa x) is 1 right
    across the slab and the result is indistinguishable from a real mode.
    """
    x = np.linspace(-2e-6, 2e-6, 9)
    assert np.all(np.isnan(ofiber.TE_field(5.0, THICKNESS, x, 2)))
    assert np.all(np.isnan(ofiber.TM_field(5.0, N1, N2, THICKNESS, x, 2)))


def test_an_unguided_field_keeps_the_shape_of_its_input():
    """Verify the nan result matches the positions asked for, scalar or array."""
    x = np.linspace(-2e-6, 2e-6, 7)
    assert ofiber.TE_field(5.0, THICKNESS, x, 2).shape == x.shape
    assert np.isnan(ofiber.TE_field(5.0, THICKNESS, 0.0, 2))


def test_the_last_guided_mode_still_produces_a_field():
    """Verify the nan guard does not swallow the highest mode the slab does guide."""
    x = np.linspace(-2e-6, 2e-6, 9)
    field = ofiber.TE_field(5.0, THICKNESS, x, 1)
    assert np.all(np.isfinite(field))
    assert np.any(field != 0)
