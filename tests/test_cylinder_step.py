# pylint: disable=invalid-name
# pylint: disable=no-member
"""
Tests for the ofiber.cylinder_step module.

Two independent sources of truth are used.

The mode tables are the committed outputs of `4-Circular-Step-Index-Fiber.ipynb`,
so a change that moves an eigenvalue breaks these tests and the published
documentation together rather than one silently drifting from the other.

The eigenvalue check itself does not reuse the solver's own equation.  The
module solves Ghatak eq. 8.40, written with J_(l-1) and K_(l-1); the tests
below verify the returned b against the equivalent form

    U J_(l+1)(U) / J_l(U) = W K_(l+1)(W) / K_l(W)

which follows from different Bessel recurrences, so an error in the module's
LHS or RHS cannot cancel out of the comparison.
"""

import matplotlib
import numpy as np
import pytest
from scipy import special
import ofiber

matplotlib.use("Agg")

# committed output of 4-Circular-Step-Index-Fiber.ipynb for V = 6.5
MODES_AT_6p5 = {
    (0, 1): 0.8977,
    (0, 2): 0.4752,
    (1, 1): 0.7422,
    (1, 2): 0.1792,
}

# committed output of the same notebook for V = 12.5
MODES_AT_12p5 = {
    (0, 1): 0.9683, (0, 2): 0.8336, (0, 3): 0.5942, (0, 4): 0.2597,
    (1, 1): 0.9196, (1, 2): 0.7319, (1, 3): 0.4427, (1, 4): 0.0735,
    (2, 1): 0.8558, (2, 2): 0.6155, (2, 3): 0.2791,
    (3, 1): 0.7777, (3, 2): 0.4850, (3, 3): 0.1070,
    (4, 1): 0.6860, (4, 2): 0.3415,
    (5, 1): 0.5812, (5, 2): 0.1863,
}

# the first two zeros of J_0 and J_1 set the LP mode cutoffs
V_CUTOFF_LP11 = 2.404825557695773
V_CUTOFF_LP02 = 3.831705970207512


def _eigen_residual(V, b, ell):
    """
    Evaluate the LP eigenvalue equation in its J_(l+1)/K_(l+1) form.

    Args:
        V:   V-parameter of the fiber        [-]
        b:   normalized propagation constant [-]
        ell: azimuthal mode number           [-]

    Returns:
        difference of the two sides, zero for a guided mode
    """
    U = V * np.sqrt(1 - b)
    W = V * np.sqrt(b)
    return U * special.jv(ell + 1, U) / special.jv(ell, U) - W * special.kn(ell + 1, W) / special.kn(ell, W)


# --------------------------------------------------------------------------
# eigenvalues
# --------------------------------------------------------------------------


@pytest.mark.parametrize("mode, b", sorted(MODES_AT_6p5.items()))
def test_mode_values_match_the_notebook_at_V_6p5(mode, b):
    """Check each eigenvalue against the value published in the notebook."""
    ell, em = mode
    assert ofiber.LP_mode_value(6.5, ell, em) == pytest.approx(b, abs=5e-5)


@pytest.mark.parametrize("mode, b", sorted(MODES_AT_12p5.items()))
def test_mode_values_match_the_notebook_at_V_12p5(mode, b):
    """Check all eighteen guided modes published for a V of 12.5."""
    ell, em = mode
    assert ofiber.LP_mode_value(12.5, ell, em) == pytest.approx(b, abs=5e-5)


@pytest.mark.parametrize("V", [1.0, 2.5, 4.0, 6.5, 12.5, 25.0])
def test_every_mode_satisfies_the_independent_eigenvalue_equation(V):
    """Verify each solved b is a root of the J_(l+1)/K_(l+1) form of eq. 8.40."""
    found = 0
    for ell in range(8):
        for b in ofiber.LP_mode_values(V, ell):
            assert _eigen_residual(V, b, ell) == pytest.approx(0.0, abs=1e-7)
            found += 1
    assert found > 0


@pytest.mark.parametrize("V", [0.5, 1.0, 2.5, 6.5, 12.5])
def test_normalized_propagation_constant_is_between_zero_and_one(V):
    """Verify every guided mode has 0 < b < 1, so beta lies between the two indices."""
    for ell in range(6):
        for b in ofiber.LP_mode_values(V, ell):
            assert 0.0 < b < 1.0


def test_fundamental_mode_exists_at_any_practical_V():
    """Verify LP01 has no cutoff, unlike every other mode."""
    for V in (0.6, 1.0, 5.0, 50.0):
        assert ofiber.LP_mode_value(V, 0, 1) is not None


def test_modes_below_the_solver_floor_are_reported_as_absent():
    """
    Document that a mode with b under 1e-5 is not found.

    `_LP_mode_value` brackets the root in [abit, 1-abit] with abit = 1e-5, so a
    mode too close to cutoff falls outside the bracket and brentq reports no
    sign change.  LP01 has no true cutoff, yet it goes missing below V of about
    0.55, purely because b has dropped under that floor.
    """
    assert ofiber.LP_mode_value(0.5, 0, 1) is None
    assert ofiber.LP_mode_value(0.6, 0, 1) > 1e-5


def test_fundamental_mode_is_more_confined_as_V_grows():
    """Verify b rises monotonically with V and approaches one."""
    V = np.linspace(0.6, 40, 200)
    b = np.array([ofiber.LP_mode_value(v, 0, 1) for v in V])
    assert np.all(np.diff(b) > 0)
    assert b[-1] > 0.97


def test_lp11_appears_at_the_single_mode_cutoff():
    """Verify LP11 is absent below V=2.405 and present above it, the single-mode limit."""
    assert ofiber.LP_mode_value(V_CUTOFF_LP11 - 1e-3, 1, 1) is None
    assert ofiber.LP_mode_value(V_CUTOFF_LP11 + 1e-3, 1, 1) is not None


def test_second_mode_group_appears_at_the_next_bessel_zero():
    """
    Verify LP02 and LP21 both switch on at the first zero of J_1.

    The margin above cutoff differs between the two only because b climbs away
    from zero at different rates, and the solver cannot see b under 1e-5.
    """
    for mode, margin in (((0, 2), 0.1), ((2, 1), 1e-3)):
        assert ofiber.LP_mode_value(V_CUTOFF_LP02 - 1e-3, *mode) is None
        assert ofiber.LP_mode_value(V_CUTOFF_LP02 + margin, *mode) is not None


def test_a_fiber_below_cutoff_is_single_mode():
    """Verify only LP01 is guided just under the LP11 cutoff."""
    V = V_CUTOFF_LP11 - 1e-3
    assert len(ofiber.LP_mode_values(V, 0)) == 1
    assert len(ofiber.LP_mode_values(V, 1)) == 0


def test_mode_count_grows_with_V():
    """Verify a larger V guides at least as many modes of each order."""
    for ell in range(4):
        counts = [len(ofiber.LP_mode_values(V, ell)) for V in (2.0, 5.0, 10.0, 20.0)]
        assert counts == sorted(counts)


def test_mode_values_are_ordered_by_radial_number():
    """Verify b falls as the radial mode number rises, LP01 being the most confined."""
    for V in (6.5, 12.5):
        for ell in range(4):
            all_b = ofiber.LP_mode_values(V, ell)
            assert np.all(np.diff(all_b) < 0)


def test_negative_ell_matches_positive_ell():
    """Verify the documented equivalence of -ell and +ell."""
    assert ofiber.LP_mode_value(6.5, -1, 1) == ofiber.LP_mode_value(6.5, 1, 1)


def test_mode_number_em_must_be_positive():
    """Verify a nonsensical radial mode number is rejected rather than guessed at."""
    with pytest.raises(ValueError):
        ofiber.LP_mode_value(6.5, 0, 0)


def test_only_one_argument_may_be_an_array():
    """Verify the wrapper refuses ambiguous combinations of array arguments."""
    with pytest.raises(ValueError):
        ofiber.LP_mode_value(np.array([2.0, 3.0]), np.array([0, 1]), 1)


def test_mode_values_is_empty_when_no_mode_is_guided():
    """Verify a high order with a small V gives an empty array, not an error."""
    assert len(ofiber.LP_mode_values(2.0, 5)) == 0


def test_array_of_V_matches_scalar_calls():
    """Verify a float array of V gives the same answers as looping."""
    V = np.linspace(3.0, 9.0, 13)
    got = ofiber.LP_mode_value(V, 0, 1)
    assert np.allclose(got, [ofiber.LP_mode_value(v, 0, 1) for v in V])


# --------------------------------------------------------------------------
# irradiance and fields
# --------------------------------------------------------------------------


@pytest.mark.parametrize("V", [2.0, 4.0, 8.0])
def test_core_and_cladding_irradiance_sum_to_the_total(V):
    """Verify Ghatak eqs. 8.56 and 8.57 add up to eq. 8.58."""
    b = ofiber.LP_mode_value(V, 0, 1)
    core = ofiber.LP_core_irradiance(V, b, 0)
    clad = ofiber.LP_clad_irradiance(V, b, 0)
    assert core + clad == pytest.approx(ofiber.LP_total_irradiance(V, b, 0), rel=1e-12)


def test_confinement_improves_with_V():
    """Verify the fraction of power in the core rises toward one as V grows."""
    fractions = []
    for V in (1.5, 2.5, 5.0, 10.0, 20.0):
        b = ofiber.LP_mode_value(V, 0, 1)
        fractions.append(ofiber.LP_core_irradiance(V, b, 0) / ofiber.LP_total_irradiance(V, b, 0))
    assert np.all(np.diff(fractions) > 0)
    assert 0.0 < fractions[0] < fractions[-1] < 1.0


@pytest.mark.parametrize("V", [2.0, 4.0])
def test_radial_irradiance_is_normalized(V):
    """Verify the documented normalization, 2*trapezoid(I*r, r) = 1."""
    b = ofiber.LP_mode_value(V, 0, 1)
    r = np.linspace(0, 60, 200001)
    total = 2 * np.trapezoid(ofiber.LP_radial_irradiance(V, b, 0, r) * r, r)
    assert total == pytest.approx(1.0, rel=1e-4)


def test_radial_field_is_continuous_across_the_core_boundary():
    """Verify the Bessel and modified-Bessel branches meet at r = a."""
    V = 3.0
    b = ofiber.LP_mode_value(V, 0, 1)
    inside = ofiber.LP_radial_field(V, b, 0, 1 - 1e-7)
    outside = ofiber.LP_radial_field(V, b, 0, 1 + 1e-7)
    assert inside == pytest.approx(outside, rel=1e-5)


def test_radial_field_is_symmetric_in_radius():
    """Verify negative radii mirror positive ones."""
    V = 3.0
    b = ofiber.LP_mode_value(V, 0, 1)
    r = np.linspace(0.1, 3, 30)
    assert np.allclose(ofiber.LP_radial_field(V, b, 0, r), ofiber.LP_radial_field(V, b, 0, -r))


def test_radial_field_decays_outside_the_core():
    """Verify the cladding field is evanescent."""
    V = 3.0
    b = ofiber.LP_mode_value(V, 0, 1)
    r = np.linspace(1.0, 8.0, 60)
    field = ofiber.LP_radial_field(V, b, 0, r)
    assert np.all(np.diff(field) < 0)
    assert field[-1] < 1e-3 * field[0]


def test_radial_irradiance_is_the_square_of_the_field():
    """Verify the two radial routines stay consistent."""
    V = 3.0
    b = ofiber.LP_mode_value(V, 0, 1)
    r = np.linspace(0, 4, 50)
    assert np.allclose(ofiber.LP_radial_irradiance(V, b, 0, r), ofiber.LP_radial_field(V, b, 0, r) ** 2)


def test_gaussian_envelope_approximates_the_true_mode():
    """Verify the Gaussian envelope tracks the real LP01 irradiance near the axis."""
    V = 2.4
    b = ofiber.LP_mode_value(V, 0, 1)
    r = np.linspace(0, 1.0, 40)
    exact = ofiber.LP_radial_irradiance(V, b, 0, r)
    approx = ofiber.gaussian_radial_irradiance(V, r)
    assert np.allclose(exact, approx, rtol=0.15)


def test_gaussian_irradiance_is_normalized():
    """Verify the documented normalization of the Gaussian envelope."""
    V = 2.4
    r = np.linspace(0, 40, 200001)
    assert np.trapezoid(ofiber.gaussian_radial_irradiance(V, r) * r, r) == pytest.approx(0.5, rel=1e-5)


# --------------------------------------------------------------------------
# mode field radii
# --------------------------------------------------------------------------


def test_mode_field_diameter_is_twice_the_radius():
    """Verify MFD and MFR stay in step."""
    for V in (1.5, 2.0, 2.5):
        assert ofiber.MFD(V) == pytest.approx(2 * ofiber.MFR(V), rel=1e-15, abs=0)


def test_mode_field_radius_shrinks_as_V_grows():
    """Verify tighter confinement gives a smaller mode field."""
    V = np.linspace(1.2, 4.0, 40)
    assert np.all(np.diff(ofiber.MFR(V)) < 0)


def test_petermann_radius_matches_the_notebook():
    """Check the Petermann-2 radius and MFD published for V = 1.5."""
    assert ofiber.PetermannW(1.5) == pytest.approx(1.69342, abs=5e-5)
    assert ofiber.MFD(1.5) == pytest.approx(3.56805, abs=5e-5)


def test_petermann_radius_matches_its_approximation():
    """Verify the exact and approximate Petermann-2 radii agree over 1.5 < V < 2.5."""
    V = np.linspace(1.5, 2.5, 21)
    assert np.allclose(ofiber.PetermannW(V), ofiber.PetermannW_Approx(V), rtol=0.02)


def test_petermann_radius_is_nan_without_a_guided_mode():
    """Verify an unguided fiber gives nan rather than a plausible number."""
    assert np.isnan(ofiber.PetermannW(-1.0))


# --------------------------------------------------------------------------
# waveguide dispersion term
# --------------------------------------------------------------------------


def test_waveguide_dispersion_term_matches_its_approximation():
    """
    Verify the exact and approximate V d2(bV)/dV2 agree near the single mode range.

    The docstring's 1% claim holds over 1.35 < V < 2.08; this pins both that
    the claim is met inside the range and that it fails outside it, so the
    stated bounds cannot quietly drift.
    """
    inside = np.linspace(1.35, 2.08, 41)
    assert np.allclose(ofiber.V_d2bV_by_V(inside, 0), ofiber.V_d2bV_by_V_Approx(inside), rtol=0.01)
    assert abs(ofiber.V_d2bV_by_V_Approx(2.4) / ofiber.V_d2bV_by_V(2.4, 0) - 1) == pytest.approx(0.06, abs=0.005)


@pytest.mark.parametrize("V", [1.0, 1.5, 2.0, 2.405, 3.0, 4.0, 5.0, 8.0])
def test_waveguide_dispersion_term_matches_a_numerical_second_derivative(V):
    """
    Verify eqn 10.14 against differencing b(V)*V directly.

    This does not reuse the closed form at all, so it independently confirms
    the sign change near V=3.01: V*d2(bV)/dV2 really does go negative once
    the fiber is well past the LP11 cutoff, which is outside the single mode
    range this term is normally used in.
    """
    h = 1e-4

    def bV(x):
        return ofiber.LP_mode_value(x, 0, 1) * x

    numerical = V * (bV(V + h) - 2 * bV(V) + bV(V - h)) / h**2
    assert ofiber.V_d2bV_by_V(V, 0) == pytest.approx(numerical, abs=1e-5)


def test_waveguide_dispersion_term_changes_sign_past_the_single_mode_range():
    """Verify the term is positive where single mode fibers operate and negative beyond."""
    assert ofiber.V_d2bV_by_V(2.405, 0) > 0
    assert ofiber.V_d2bV_by_V(4.0, 0) < 0
    V = np.linspace(1.0, 8.0, 4001)
    crossings = V[np.where(np.diff(np.sign(ofiber.V_d2bV_by_V(V, 0))) != 0)[0]]
    assert len(crossings) == 1
    assert crossings[0] == pytest.approx(3.007, abs=0.01)


def test_waveguide_dispersion_term_falls_through_the_single_mode_range():
    """Verify the term decreases with V where the fiber is single mode."""
    V = np.linspace(1.4, 2.4, 30)
    assert np.all(np.diff(ofiber.V_d2bV_by_V(V, 0)) < 0)


def test_waveguide_dispersion_term_is_nan_without_a_mode():
    """
    Verify an unguided mode gives nan rather than zero.

    Zero is a legal value of V d2(bV)/dV2, so it could not be told apart from
    a fiber that simply has no waveguide dispersion at that V.
    """
    assert np.isnan(ofiber.V_d2bV_by_V(2.0, 5))
    assert np.all(np.isnan(ofiber.V_d2bV_by_V(np.array([2.0, 2.2]), 5)))


def test_waveguide_dispersion_term_is_finite_where_a_mode_exists():
    """Verify the nan guard does not swallow guided modes."""
    assert np.all(np.isfinite(ofiber.V_d2bV_by_V(np.linspace(1.4, 2.4, 11), 0)))


# --------------------------------------------------------------------------
# coupling losses
# --------------------------------------------------------------------------


def test_perfectly_aligned_fibers_lose_nothing():
    """Verify all three misalignment losses vanish when the fibers are aligned."""
    assert ofiber.transverse_misalignment_loss_db(5e-6, 5e-6, 0.0) == pytest.approx(0.0)
    assert ofiber.angular_misalignment_loss_db(1.0, 5e-6, 0.0, 1.55e-6) == pytest.approx(0.0)
    assert ofiber.longitudinal_misalignment_loss_db(1.0, 5e-6, 0.0, 1.55e-6) == pytest.approx(0.0)


@pytest.mark.parametrize(
    "loss, args",
    [
        (ofiber.transverse_misalignment_loss_db, (5e-6, 5e-6)),
        (ofiber.angular_misalignment_loss_db, None),
        (ofiber.longitudinal_misalignment_loss_db, None),
    ],
)
def test_misalignment_losses_grow_with_the_offset(loss, args):
    """Verify every misalignment loss increases monotonically away from alignment."""
    offsets = np.linspace(0, 2e-6, 20)
    if args is None:
        values = np.array([loss(1.0, 5e-6, d, 1.55e-6) for d in offsets])
    else:
        values = np.array([loss(*args, d) for d in offsets])
    assert np.all(np.diff(values) > 0)


def test_mismatched_mode_fields_lose_power_even_when_centred():
    """Verify unequal mode field radii cost power with no transverse offset."""
    assert ofiber.transverse_misalignment_loss_db(4e-6, 6e-6, 0.0) > 0


def test_bending_loss_falls_as_the_bend_radius_grows():
    """Verify a gentler bend leaks less power."""
    Rc = np.array([0.005, 0.01, 0.02, 0.05, 0.1])
    loss = ofiber.bending_loss_db(1.46, 0.005, 4e-6, Rc, 1.55e-6)
    assert np.all(np.diff(loss) < 0)
    assert np.all(loss > 0)


def test_bending_loss_accepts_an_array_of_core_radii():
    """Verify the scalar and array paths of bending_loss_db agree."""
    a = np.linspace(3e-6, 5e-6, 7)
    got = ofiber.bending_loss_db(1.46, 0.005, a, 0.01, 1.55e-6)
    assert np.allclose(got, [ofiber.bending_loss_db(1.46, 0.005, aa, 0.01, 1.55e-6) for aa in a])


# --------------------------------------------------------------------------
# far field
# --------------------------------------------------------------------------


def test_far_field_node_is_a_zero_of_the_polar_factor():
    """
    Verify the reported node really is a zero of the far-field polar factor.

    Note the returned quantity is k*a*sin(Theta), not the polar angle itself.
    """
    V = 3.0
    b = ofiber.LP_mode_value(V, 0, 1)
    kasin = ofiber.FF_node_polar_angle(V, 0, 1)
    # pylint: disable=protected-access
    assert ofiber.cylinder_step._FF_polar_x(kasin, V, 0, b) == pytest.approx(0.0, abs=1e-12)


def test_far_field_has_no_node_below_the_reported_one():
    """Verify the node found is the first one, as documented."""
    V = 3.0
    b = ofiber.LP_mode_value(V, 0, 1)
    kasin = ofiber.FF_node_polar_angle(V, 0, 1)
    # pylint: disable=protected-access
    sampled = ofiber.cylinder_step._FF_polar_x(np.linspace(1e-5, 0.999 * kasin, 500), V, 0, b)
    assert np.all(sampled > 0) or np.all(sampled < 0)


def test_far_field_irradiance_is_positive():
    """Verify the far-field irradiance is a squared magnitude."""
    V = 3.0
    b = ofiber.LP_mode_value(V, 0, 1)
    theta = np.linspace(0.01, 0.5, 40)
    assert np.all(ofiber.FF_irradiance_x(1.0, theta, 0.0, 0, 1.55e-6, 4e-6, V, b) >= 0)
    assert np.all(ofiber.FF_polar_irradiance_x(1.0, theta, 0, 1.55e-6, 4e-6, V, b) >= 0)


def test_far_field_irradiance_falls_off_with_distance():
    """Verify the inverse square dependence on range."""
    V = 3.0
    b = ofiber.LP_mode_value(V, 0, 1)
    near = ofiber.FF_polar_irradiance_x(1.0, 0.1, 0, 1.55e-6, 4e-6, V, b)
    far = ofiber.FF_polar_irradiance_x(2.0, 0.1, 0, 1.55e-6, 4e-6, V, b)
    assert far == pytest.approx(near / 4, rel=1e-12, abs=0)


def test_far_field_node_requires_a_guided_mode():
    """Verify an unguided combination gives nan rather than a plausible angle."""
    assert np.isnan(ofiber.FF_node_polar_angle(2.0, 5, 1))


# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------


def test_mode_plot_builds_without_error():
    """Verify the eigenvalue plot runs, since the notebooks depend on it."""
    plot = ofiber.plot_LP_modes(6.5, 0)
    assert plot is not None
    plot.close("all")


# --------------------------------------------------------------------------
# module surface
# --------------------------------------------------------------------------


def test_everything_in_all_is_reachable_from_the_package():
    """Verify each exported name survives the star import in __init__.py."""
    for name in ofiber.cylinder_step.__all__:
        assert hasattr(ofiber, name), "%s is in __all__ but missing from ofiber" % name


@pytest.mark.parametrize("name", ["np", "plt", "scipy", "special"])
def test_imported_modules_do_not_leak_into_the_package(name):
    """
    Verify star importing this module exports no third party names.

    The guard is __all__; while it was misspelled _all_ the module had no
    __all__ at all and ofiber.np, ofiber.plt, ofiber.scipy and ofiber.special
    all silently existed.
    """
    assert not hasattr(ofiber, name)


# --------------------------------------------------------------------------
# integer arrays
#
# Every one of these routines sized its output with np.zeros_like or
# np.empty_like on an argument that is naturally an integer array, so the
# float result was truncated on assignment and every entry came back 0.
# --------------------------------------------------------------------------


def test_integer_array_of_V_gives_float_results():
    """Verify an integer array of V is not truncated to zeros."""
    assert ofiber.LP_mode_value(np.array([3, 4]), 0, 1) == pytest.approx([0.651471, 0.772734], abs=1e-5)


def test_integer_array_of_ell_gives_float_results():
    """Verify an integer array of ell is not truncated to zeros."""
    got = ofiber.LP_mode_value(6.5, np.array([0, 1]), 1)
    assert got == pytest.approx([MODES_AT_6p5[(0, 1)], MODES_AT_6p5[(1, 1)]], abs=5e-5)


def test_integer_array_of_em_gives_float_results():
    """Verify an integer array of em is not truncated to zeros."""
    got = ofiber.LP_mode_value(6.5, 0, np.array([1, 2]))
    assert got == pytest.approx([MODES_AT_6p5[(0, 1)], MODES_AT_6p5[(0, 2)]], abs=5e-5)


@pytest.mark.parametrize("routine", [ofiber.PetermannW, lambda V: ofiber.V_d2bV_by_V(V, 0)])
def test_integer_array_of_V_is_not_truncated_by_the_single_argument_routines(routine):
    """Verify PetermannW and V_d2bV_by_V agree on integer and float arrays of V."""
    assert np.allclose(routine(np.array([3, 4, 5])), routine(np.array([3.0, 4.0, 5.0])))


def test_integer_array_of_core_radii_is_not_truncated():
    """Verify bending_loss_db handles an integer-valued array of core radii."""
    got = ofiber.bending_loss_db(1.46, 0.005, np.array([4, 5]), 1, 1)
    assert np.allclose(got, ofiber.bending_loss_db(1.46, 0.005, np.array([4.0, 5.0]), 1, 1))


def test_integer_array_of_ell_for_the_far_field_node():
    """Verify FF_node_polar_angle survives an integer array of ell."""
    got = ofiber.FF_node_polar_angle(6.5, np.array([0, 1]), 1)
    assert np.allclose(got, [ofiber.FF_node_polar_angle(6.5, 0, 1), ofiber.FF_node_polar_angle(6.5, 1, 1)])
    assert np.all(got > 0)


def test_integer_array_of_V_for_the_far_field_node():
    """Verify FF_node_polar_angle survives an integer array of V."""
    got = ofiber.FF_node_polar_angle(np.array([6, 7]), 0, 1)
    assert np.allclose(got, [ofiber.FF_node_polar_angle(6.0, 0, 1), ofiber.FF_node_polar_angle(7.0, 0, 1)])


def test_integer_array_of_em_for_the_far_field_node():
    """Verify FF_node_polar_angle survives an integer array of em."""
    got = ofiber.FF_node_polar_angle(6.5, 0, np.array([1, 2]))
    assert np.allclose(got, [ofiber.FF_node_polar_angle(6.5, 0, 1), ofiber.FF_node_polar_angle(6.5, 0, 2)])


# --------------------------------------------------------------------------
# argument checking and edge cases
# --------------------------------------------------------------------------


def test_far_field_node_rejects_a_bad_radial_mode_number():
    """Verify FF_node_polar_angle guards em the way LP_mode_value does."""
    with pytest.raises(ValueError):
        ofiber.FF_node_polar_angle(6.5, 0, 0)


def test_far_field_node_allows_only_one_array_argument():
    """Verify ambiguous combinations of array arguments are refused."""
    with pytest.raises(ValueError):
        ofiber.FF_node_polar_angle(np.array([6.5, 7.0]), np.array([0, 1]), 1)


def test_far_field_node_treats_negative_ell_as_positive():
    """Verify the documented equivalence of -ell and +ell."""
    assert ofiber.FF_node_polar_angle(6.5, -1, 1) == ofiber.FF_node_polar_angle(6.5, 1, 1)


def test_far_field_node_needs_a_positive_V():
    """Verify a nonphysical V gives None rather than a number."""
    assert ofiber.FF_node_polar_angle(-1.0, 0, 1) is None


def test_private_far_field_node_rejects_a_bad_radial_mode_number():
    """Verify the private routine guards em, since it bypasses the wrapper."""
    # pylint: disable=protected-access
    assert ofiber.cylinder_step._FF_node_polar_angle(6.5, 0, 0) is None


def test_private_mode_value_rejects_a_bad_radial_mode_number():
    """Verify the private solver guards em, since the wrapper raises before it can."""
    # pylint: disable=protected-access
    assert ofiber.cylinder_step._LP_mode_value(6.5, 0, 0) is None


def test_far_field_node_reports_a_failed_search(monkeypatch):
    """
    Verify a polar factor that never changes sign is reported, not silently wrong.

    No real mode does this, so the branch is unreachable with physical inputs
    and the polar factor is replaced with a constant to exercise it.
    """
    # pylint: disable=protected-access
    monkeypatch.setattr(ofiber.cylinder_step, "_FF_polar_x", lambda *args: 1.0)
    with pytest.raises(StopIteration):
        ofiber.FF_node_polar_angle(6.5, 0, 1)


def test_all_nine_searched_radial_modes_can_be_found():
    """
    Verify LP_mode_values fills its search range when every mode is guided.

    The loop only looks for nine radial modes, so a V above the LP09 cutoff
    exhausts it without ever hitting the break.
    """
    all_b = ofiber.LP_mode_values(26.0, 0)
    assert len(all_b) == 9
    assert np.all(np.diff(all_b) < 0)


def test_bending_loss_is_nan_without_a_guided_mode():
    """Verify a core too small to guide gives nan rather than a plausible loss."""
    assert np.isnan(ofiber.bending_loss_db(1.46, 0.005, 1e-7, 0.01, 1.55e-6))
