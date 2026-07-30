# pylint: disable=invalid-name
"""
Tests for the ofiber.graded_index module.

No notebook exercises this module, so unlike test_refraction.py there are no
published outputs to check against.  The reference values here are closed forms
from chapters 4 and 5 of Ghatak & Thyagarajan, *An Introduction to Fiber
Optics*, which the module docstrings cite:

* a power-law profile equals n_core on axis and n_clad at the core boundary
* a q=2 profile has d(n^2)/dx = -4 n_core^2 Delta x / a^2 inside the core
* eq. 4.13 launches a ray with slope tan(theta) and period 2 pi a cos(theta)/sqrt(2 Delta)
* eq. 5.4 gives an axial-ray delay of exactly n_core/c

`first_derivative` is a forward difference, so its error is a predictable
O(dx) offset rather than noise, and the tests below pin that offset exactly
instead of hiding it behind a loose tolerance.
"""

import numpy as np
import pytest
import scipy.constants
import ofiber
from ofiber.graded_index import ray_delay

N_CORE = 1.46
N_CLAD = 1.45
CORE_RADIUS = 25e-6  # m, a typical graded-index multimode core
DELTA = (N_CORE**2 - N_CLAD**2) / 2 / N_CORE**2


def _parabolic_index(x):
    """
    Evaluate a q=2 power-law profile from its closed form.

    Args:
        x: radial position [m]

    Returns:
        index of refraction [-]
    """
    return np.sqrt(N_CORE**2 * (1 - 2 * DELTA * (x / CORE_RADIUS) ** 2))


def _parabolic_dnsqr_dx(x):
    """
    Differentiate a q=2 profile analytically inside the core.

    Args:
        x: radial position [m]

    Returns:
        d(n**2)/dx [1/m]
    """
    return -4 * N_CORE**2 * DELTA * x / CORE_RADIUS**2


# --------------------------------------------------------------------------
# power_law_profile
# --------------------------------------------------------------------------


@pytest.mark.parametrize("q", [1.0, 2.0, 4.0, 10.0])
def test_profile_on_axis_is_the_core_index(q):
    """Verify every power-law profile starts at the core index."""
    assert ofiber.power_law_profile(N_CORE, N_CLAD, q, CORE_RADIUS, 0.0) == pytest.approx(N_CORE)


@pytest.mark.parametrize("q", [1.0, 2.0, 4.0, 10.0])
def test_profile_at_the_core_boundary_is_the_cladding_index(q):
    """Verify the profile meets the cladding index at r = a, for either sign of x."""
    for x in (CORE_RADIUS, -CORE_RADIUS):
        assert ofiber.power_law_profile(N_CORE, N_CLAD, q, CORE_RADIUS, x) == pytest.approx(N_CLAD)


def test_profile_is_symmetric_about_the_axis():
    """Verify n(x) = n(-x), since the profile depends only on |x|."""
    x = np.linspace(0, CORE_RADIUS, 25)
    plus = ofiber.power_law_profile(N_CORE, N_CLAD, 2.0, CORE_RADIUS, x)
    minus = ofiber.power_law_profile(N_CORE, N_CLAD, 2.0, CORE_RADIUS, -x)
    assert np.allclose(plus, minus)


def test_parabolic_profile_matches_its_closed_form():
    """Check the q=2 profile against the analytic parabolic expression."""
    x = np.linspace(-CORE_RADIUS, CORE_RADIUS, 101)
    got = ofiber.power_law_profile(N_CORE, N_CLAD, 2.0, CORE_RADIUS, x)
    assert np.allclose(got, _parabolic_index(x), rtol=1e-12)


def test_profile_decreases_across_the_core():
    """Verify the index falls monotonically from the axis to the cladding."""
    x = np.linspace(0, CORE_RADIUS, 101)
    index = ofiber.power_law_profile(N_CORE, N_CLAD, 2.0, CORE_RADIUS, x)
    assert np.all(np.diff(index) < 0)


def test_profile_outside_the_core_is_the_cladding_index():
    """Verify array input clamps to the cladding beyond the core."""
    x = np.array([1.5, 2.0, 10.0]) * CORE_RADIUS
    index = ofiber.power_law_profile(N_CORE, N_CLAD, 2.0, CORE_RADIUS, x)
    assert np.allclose(index, N_CLAD)


def test_large_q_approaches_a_step_index_profile():
    """Verify the profile tends to a step index as q grows."""
    x = np.linspace(0, 0.9 * CORE_RADIUS, 50)
    index = ofiber.power_law_profile(N_CORE, N_CLAD, 500.0, CORE_RADIUS, x)
    assert np.allclose(index, N_CORE, atol=1e-6)


@pytest.mark.parametrize("x", [1.0001, 1.5, 2.0, 100.0])
def test_scalar_profile_outside_the_core_is_the_cladding_index(x):
    """Verify a scalar position beyond the core clamps like an array does."""
    index = ofiber.power_law_profile(N_CORE, N_CLAD, 2.0, CORE_RADIUS, x * CORE_RADIUS)
    assert index == pytest.approx(N_CLAD)


def test_profile_outside_the_core_does_not_overflow_for_large_q():
    """Verify a large q outside the core stays finite rather than overflowing."""
    x = np.array([2.0, 50.0]) * CORE_RADIUS
    assert np.allclose(ofiber.power_law_profile(N_CORE, N_CLAD, 500.0, CORE_RADIUS, x), N_CLAD)


# --------------------------------------------------------------------------
# first_derivative
# --------------------------------------------------------------------------


def test_first_derivative_preserves_length():
    """Verify the result is as long as the input, which is the point of the padding."""
    x = np.linspace(0, 1, 17)
    assert len(ofiber.first_derivative(x, x**2)) == len(x)


def test_first_derivative_is_exact_for_a_straight_line():
    """Verify a forward difference reproduces a constant slope exactly."""
    x = np.linspace(0, 1, 9)
    assert np.allclose(ofiber.first_derivative(x, 3 * x + 1), 3.0)


def test_first_derivative_of_a_quadratic_carries_the_expected_offset():
    """
    Verify the forward difference of x**2 is 2x + dx.

    Pinning the offset exactly is what distinguishes a correct forward
    difference from one differenced against the wrong spacing.
    """
    x = np.linspace(0, 1, 21)
    dx = x[1] - x[0]
    assert np.allclose(ofiber.first_derivative(x, x**2), 2 * x + dx)


def test_first_derivative_has_the_sign_of_the_slope():
    """Verify a falling function gives a negative derivative."""
    x = np.linspace(0, 1, 11)
    assert np.all(ofiber.first_derivative(x, -2 * x) < 0)
    assert np.all(ofiber.first_derivative(x, 2 * x) > 0)


def test_first_derivative_is_ordered_like_the_true_derivative():
    """
    Verify the derivative of x**2 increases across the whole array.

    A padding value inserted anywhere but the end would leave the last two
    entries out of order while every individual value still looked plausible.
    """
    x = np.linspace(0, 1, 9)
    deriv = ofiber.first_derivative(x, x**2)
    assert np.all(np.diff(deriv) > 0)
    assert deriv[-1] == max(deriv)


def test_first_derivative_converges_as_the_grid_refines():
    """Verify the O(dx) error shrinks with the step size."""
    errors = []
    for n in (11, 101, 1001):
        x = np.linspace(0, 1, n)
        errors.append(np.max(np.abs(ofiber.first_derivative(x, x**2) - 2 * x)))
    assert errors[0] > errors[1] > errors[2]
    assert errors[-1] < 1e-2


def test_first_derivative_handles_a_negative_step():
    """Verify a descending x still gives the correct sign."""
    x = np.linspace(1, 0, 11)
    assert np.allclose(ofiber.first_derivative(x, 3 * x), 3.0)


# --------------------------------------------------------------------------
# curvature
# --------------------------------------------------------------------------


def test_parabolic_curvature_matches_the_analytic_derivative():
    """Check curvature against d(n^2)/dx / (2 beta^2) for a q=2 profile."""
    theta = 0.05
    x = np.linspace(-0.8 * CORE_RADIUS, 0.8 * CORE_RADIUS, 401)
    dx = x[1] - x[0]
    beta = N_CORE * np.cos(theta)
    got = ofiber.curvature(N_CORE, N_CLAD, 2.0, CORE_RADIUS, x, theta)
    # the forward difference lands midway between samples
    expected = _parabolic_dnsqr_dx(x + dx / 2) / 2 / beta**2
    assert np.allclose(got[:-1], expected[:-1], rtol=1e-9)


def test_curvature_bends_rays_back_toward_the_axis():
    """Verify the curvature is restoring, negative above the axis and positive below."""
    x = np.linspace(-0.8 * CORE_RADIUS, 0.8 * CORE_RADIUS, 201)
    curve = ofiber.curvature(N_CORE, N_CLAD, 2.0, CORE_RADIUS, x, 0.05)
    assert np.all(curve[x < -0.05 * CORE_RADIUS] > 0)
    assert np.all(curve[x > 0.05 * CORE_RADIUS] < 0)


def test_curvature_grows_with_launch_angle():
    """Verify curvature scales as 1/cos(theta)**2 through beta."""
    x = np.linspace(-0.5 * CORE_RADIUS, 0.5 * CORE_RADIUS, 101)
    shallow = ofiber.curvature(N_CORE, N_CLAD, 2.0, CORE_RADIUS, x, 0.02)
    steep = ofiber.curvature(N_CORE, N_CLAD, 2.0, CORE_RADIUS, x, 0.20)
    ratio = np.cos(0.02) ** 2 / np.cos(0.20) ** 2
    assert np.allclose(steep, shallow * ratio, rtol=1e-12)


# --------------------------------------------------------------------------
# transverse_location, Ghatak eq. 4.13
# --------------------------------------------------------------------------


def test_ray_starts_on_the_axis():
    """Verify the ray launches from x = 0."""
    assert ofiber.transverse_location(N_CORE, 0.1, DELTA, CORE_RADIUS, 0.0) == pytest.approx(0.0)


@pytest.mark.parametrize("theta", [0.01, 0.05, 0.1, 0.2])
def test_initial_ray_slope_is_tan_theta(theta):
    """
    Verify dx/dz at the launch point equals tan(theta).

    This ties the amplitude and the period together: either one alone being
    wrong would break the product that sets the initial slope.
    """
    dz = 1e-9
    plus = ofiber.transverse_location(N_CORE, theta, DELTA, CORE_RADIUS, dz)
    minus = ofiber.transverse_location(N_CORE, theta, DELTA, CORE_RADIUS, -dz)
    assert (plus - minus) / (2 * dz) == pytest.approx(np.tan(theta), rel=1e-6)


def test_ray_returns_to_the_axis_after_one_period():
    """Verify the sinusoidal path has period 2 pi a cos(theta)/sqrt(2 Delta)."""
    theta = 0.1
    period = 2 * np.pi * CORE_RADIUS * np.cos(theta) / np.sqrt(2 * DELTA)
    for z in (period, 2 * period, 3 * period):
        assert ofiber.transverse_location(N_CORE, theta, DELTA, CORE_RADIUS, z) == pytest.approx(0.0, abs=1e-12)


def test_ray_excursion_peaks_at_a_quarter_period():
    """Verify the maximum excursion is a sin(theta)/sqrt(2 Delta)."""
    theta = 0.1
    period = 2 * np.pi * CORE_RADIUS * np.cos(theta) / np.sqrt(2 * DELTA)
    peak = ofiber.transverse_location(N_CORE, theta, DELTA, CORE_RADIUS, period / 4)
    assert peak == pytest.approx(CORE_RADIUS * np.sin(theta) / np.sqrt(2 * DELTA))


def test_shallower_rays_stay_closer_to_the_axis():
    """Verify a smaller launch angle gives a smaller excursion."""
    z = np.linspace(0, 1e-3, 500)
    shallow = np.max(np.abs(ofiber.transverse_location(N_CORE, 0.02, DELTA, CORE_RADIUS, z)))
    steep = np.max(np.abs(ofiber.transverse_location(N_CORE, 0.20, DELTA, CORE_RADIUS, z)))
    assert shallow < steep


# --------------------------------------------------------------------------
# ray_delay, Ghatak eq. 5.4
# --------------------------------------------------------------------------
# This returns a delay per unit length, in s/m: an axial ray comes back as
# n_core/c, not c/n_core.


@pytest.mark.parametrize("q", [1.0, 2.0, 4.0, 10.0])
def test_axial_ray_delay_is_the_core_index_over_c(q):
    """
    Check the on-axis ray against n_core/c for every profile.

    An axial ray has beta = n_core and never leaves the axis, so its delay
    cannot depend on q.  The two q-dependent terms of eq. 5.4 must cancel.
    """
    assert ray_delay(N_CORE, q, N_CORE) == pytest.approx(N_CORE / scipy.constants.speed_of_light, rel=1e-12, abs=0)


@pytest.mark.parametrize("q", [1.0, 2.0, 4.0])
def test_axial_ray_delay_with_dispersion_is_the_group_index_over_c(q):
    """Verify a dispersive axial ray travels at the group velocity c/N1."""
    md = 0.021
    got = ray_delay(N_CORE, q, N_CORE, material_dispersion=md)
    assert got == pytest.approx((N_CORE + md) / scipy.constants.speed_of_light, rel=1e-12, abs=0)


def test_material_dispersion_scales_the_whole_delay():
    """Verify dispersion scales the delay by N1/n_core at every ray invariant."""
    md = 0.021
    beta = np.linspace(N_CLAD, N_CORE, 50)
    plain = ray_delay(N_CORE, 2.0, beta)
    dispersive = ray_delay(N_CORE, 2.0, beta, material_dispersion=md)
    assert np.allclose(dispersive, plain * (N_CORE + md) / N_CORE, rtol=1e-12)


def test_zero_material_dispersion_reduces_to_the_simple_form():
    """Verify passing a zero dispersion matches omitting it, since then N1 = n_core."""
    beta = np.linspace(N_CLAD, N_CORE, 50)
    assert np.allclose(ray_delay(N_CORE, 2.0, beta, material_dispersion=0.0), ray_delay(N_CORE, 2.0, beta), rtol=1e-15)


@pytest.mark.parametrize("q", [1.0, 2.0, 4.0])
def test_zero_profile_dispersion_is_the_default(q):
    """Verify P=0 matches omitting the argument."""
    beta = np.linspace(N_CLAD, N_CORE, 50)
    got = ray_delay(N_CORE, q, beta, material_dispersion=0.02, profile_dispersion=0.0)
    assert np.allclose(got, ray_delay(N_CORE, q, beta, material_dispersion=0.02), rtol=1e-15)


@pytest.mark.parametrize("q", [1.0, 2.0, 4.0])
@pytest.mark.parametrize("P", [-0.3, -0.1, 0.0, 0.15, 0.4])
def test_axial_ray_is_unaffected_by_profile_dispersion(q, P):
    """
    Verify the axial ray still takes N1/c however the profile disperses.

    A ray on the axis never samples r > 0, so no property of the profile shape
    can reach it.  This is the tightest constraint on the two coefficients.
    """
    md = 0.021
    got = ray_delay(N_CORE, q, N_CORE, material_dispersion=md, profile_dispersion=P)
    assert got == pytest.approx((N_CORE + md) / scipy.constants.speed_of_light, rel=1e-12, abs=0)


@pytest.mark.parametrize("P", [-0.2, -0.1, 0.0, 0.1, 0.2, 0.3])
def test_optimal_exponent_is_two_minus_twice_the_profile_dispersion(P):
    """
    Verify the delay-equalising exponent is q = 2 - 2P.

    Nothing in the coefficients was fitted to this: it is an independent
    consequence of them, and it reproduces the g_opt of the standard
    graded-index treatment.  A tiny Delta keeps the O(Delta) correction to
    the optimum below the search resolution.
    """
    delta = 1e-6
    n_edge = N_CORE * (1 - delta)
    q = np.linspace(0.5, 3.5, 60001)
    edge = ray_delay(N_CORE, q, n_edge, profile_dispersion=P)
    axial = ray_delay(N_CORE, q, N_CORE, profile_dispersion=P)
    assert q[np.argmin(np.abs(edge - axial))] == pytest.approx(2 - 2 * P, abs=1e-3)


def test_profile_dispersion_shifts_delay_in_opposite_directions_across_the_ray_fan():
    """Verify P tilts the delay curve, slowing edge rays while the axial ray is fixed."""
    beta = N_CORE * (1 - 0.01)
    slower = ray_delay(N_CORE, 2.0, beta, profile_dispersion=0.2)
    faster = ray_delay(N_CORE, 2.0, beta, profile_dispersion=-0.2)
    assert slower > ray_delay(N_CORE, 2.0, beta) > faster


def test_delay_is_minimised_at_the_predicted_invariant():
    """Verify the fastest ray has beta = n_core sqrt(q/2), where dtau/dbeta vanishes."""
    q = 2.0
    beta = np.linspace(1.20, 1.46, 20001)
    tau = ray_delay(N_CORE, q, beta)
    assert beta[np.argmin(tau)] == pytest.approx(N_CORE * np.sqrt(q / 2), abs=1e-4)


def test_delay_is_positive_for_every_guided_ray():
    """Verify the delay stays physical across the range of ray invariants."""
    beta = np.linspace(N_CLAD, N_CORE, 200)
    assert np.all(ray_delay(N_CORE, 2.0, beta) > 0)


# --------------------------------------------------------------------------
# module surface
# --------------------------------------------------------------------------


def test_everything_in_all_is_reachable_from_the_package():
    """Verify each exported name survives the star import in __init__.py."""
    for name in ofiber.graded_index.__all__:
        assert hasattr(ofiber, name), "%s is in __all__ but missing from ofiber" % name


def test_ray_delay_is_exported_from_the_package():
    """Verify ray_delay is reachable as ofiber.ray_delay like its neighbours."""
    assert hasattr(ofiber, "ray_delay")
    assert ofiber.ray_delay is ray_delay
