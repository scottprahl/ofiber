# pylint: disable=invalid-name
"""
Tests for the ofiber.basics module.

The cutoff wavelengths and index values checked here are the committed outputs
of `0-Basics.ipynb`, so a change that moves them breaks the tests and the
published documentation together.

Everything else is a closed form or an identity that ties two routines
together, for instance:

* NA from a pair of indices equals NA from the Delta those indices imply
* a V-parameter rebuilt from esi_radius and esi_Delta equals esi_V_parameter
* the esi routines all tend to their step-index values as q grows
* Fresnel reflection at normal incidence is ((m-1)/(m+1))**2 for both
  polarizations, and the parallel one vanishes at Brewster's angle
"""

import numpy as np
import pytest
from scipy.special import jn_zeros
import ofiber

N_CORE = 1.46
N_CLAD = 1.45


# --------------------------------------------------------------------------
# acceptance and critical angles
# --------------------------------------------------------------------------


def test_acceptance_angle_inverts_the_numerical_aperture():
    """Verify the acceptance angle is arcsin(NA) in air."""
    for NA in (0.1, 0.22, 0.5):
        assert np.sin(ofiber.acceptance_angle(NA)) == pytest.approx(NA)


def test_acceptance_angle_narrows_in_a_denser_medium():
    """Verify immersing the fiber face shrinks the acceptance cone."""
    assert ofiber.acceptance_angle(0.3, 1.33) < ofiber.acceptance_angle(0.3)


def test_acceptance_angle_is_nan_when_no_cone_exists():
    """Verify an NA above the outside index gives nan, as documented."""
    with np.errstate(invalid="ignore"):
        assert np.isnan(ofiber.acceptance_angle(1.6))


def test_critical_angle_of_a_typical_fiber():
    """Check the angle for total internal reflection at a core-cladding boundary."""
    assert ofiber.critical_angle(N_CORE, N_CLAD) == pytest.approx(np.arcsin(N_CLAD / N_CORE))


def test_critical_angle_is_ninety_degrees_without_an_index_step():
    """Verify equal indices guide nothing, so the critical angle is grazing."""
    assert ofiber.critical_angle(1.5, 1.5) == pytest.approx(np.pi / 2)


def test_a_larger_index_step_lowers_the_critical_angle():
    """Verify a stronger contrast accepts a wider range of internal rays."""
    angles = [ofiber.critical_angle(N_CORE, nc) for nc in (1.30, 1.40, 1.45)]
    assert np.all(np.diff(angles) > 0)


# --------------------------------------------------------------------------
# cutoff wavelength
# --------------------------------------------------------------------------


def test_cutoff_wavelength_matches_the_notebook():
    """Check the step-index and parabolic cutoffs published in 0-Basics.ipynb."""
    r_core = 4e-6
    NA = 2.1 / 2 / np.pi * 1300e-9 / r_core
    assert ofiber.cutoff_wavelength(r_core, NA) * 1e9 == pytest.approx(1135, abs=1)
    assert ofiber.cutoff_wavelength(r_core, NA, q=2) * 1e9 == pytest.approx(803, abs=1)


def test_a_second_notebook_cutoff():
    """Check the 863 nm cutoff published for a 2 micron core."""
    Delta = 0.0064
    n_core = 1.45 / (1 - Delta)
    NA = ofiber.numerical_aperture_from_Delta(n_core, Delta)
    assert ofiber.cutoff_wavelength(2e-6, NA) * 1e9 == pytest.approx(863, abs=1)


def test_cutoff_wavelength_uses_the_first_zero_of_J0():
    """Verify the default cutoff is set by V=2.405, the LP11 cutoff."""
    a, NA = 4e-6, 0.12
    assert ofiber.cutoff_wavelength(a, NA) == pytest.approx(2 * np.pi * a * NA / jn_zeros(0, 1)[0])


def test_higher_ell_uses_the_next_bessel_zero():
    """Verify ell selects the Bessel function whose first zero sets the cutoff."""
    a, NA = 4e-6, 0.12
    for ell in (0, 1, 2):
        expected = 2 * np.pi * a * NA / jn_zeros(ell, 1)[0]
        assert ofiber.cutoff_wavelength(a, NA, ell=ell) == pytest.approx(expected)


def test_higher_order_modes_cut_off_at_shorter_wavelengths():
    """Verify a larger ell gives a shorter cutoff wavelength."""
    a, NA = 4e-6, 0.12
    cutoffs = [ofiber.cutoff_wavelength(a, NA, ell=ell) for ell in range(4)]
    assert np.all(np.diff(cutoffs) < 0)


def test_a_fractional_ell_is_truncated():
    """Verify ell is truncated to an integer, as documented."""
    a, NA = 4e-6, 0.12
    assert ofiber.cutoff_wavelength(a, NA, ell=1.9) == ofiber.cutoff_wavelength(a, NA, ell=1)


@pytest.mark.parametrize("q", [1.0, 2.0, 4.0, 10.0])
def test_graded_index_cutoff_scales_by_sqrt_one_plus_two_over_q(q):
    """Verify the graded-index cutoff V is larger by sqrt(1+2/q)."""
    a, NA = 4e-6, 0.12
    ratio = ofiber.cutoff_wavelength(a, NA) / ofiber.cutoff_wavelength(a, NA, q=q)
    assert ratio == pytest.approx(np.sqrt(1 + 2 / q))


def test_a_huge_q_approaches_the_step_index_cutoff():
    """Verify a very large q recovers the step-index result."""
    a, NA = 4e-6, 0.12
    assert ofiber.cutoff_wavelength(a, NA, q=1e12) == pytest.approx(ofiber.cutoff_wavelength(a, NA), rel=1e-9)


def test_cutoff_wavelength_grows_with_core_and_aperture():
    """Verify a bigger or more strongly guiding core pushes the cutoff to longer waves."""
    assert ofiber.cutoff_wavelength(8e-6, 0.12) > ofiber.cutoff_wavelength(4e-6, 0.12)
    assert ofiber.cutoff_wavelength(4e-6, 0.24) > ofiber.cutoff_wavelength(4e-6, 0.12)


# --------------------------------------------------------------------------
# numerical aperture and Delta
# --------------------------------------------------------------------------


def test_numerical_aperture_matches_the_notebook():
    """Check the NA published in 0-Basics.ipynb for a 1.5/1.48 fiber."""
    assert ofiber.numerical_aperture(1.5, 1.48) == pytest.approx(0.2441, abs=5e-5)


def test_relative_refractive_index_matches_the_notebook():
    """Check the exact Delta published for the same fiber."""
    assert ofiber.relative_refractive_index(1.5, 1.48) == pytest.approx(0.01324, abs=5e-6)


def test_numerical_aperture_agrees_whichever_way_it_is_computed():
    """Verify NA from two indices equals NA from the Delta they imply."""
    Delta = ofiber.relative_refractive_index(N_CORE, N_CLAD)
    assert ofiber.numerical_aperture_from_Delta(N_CORE, Delta) == pytest.approx(
        ofiber.numerical_aperture(N_CORE, N_CLAD), rel=1e-12
    )


def test_numerical_aperture_is_zero_without_an_index_step():
    """Verify a fiber with no contrast has no acceptance cone."""
    assert ofiber.numerical_aperture(1.5, 1.5) == pytest.approx(0.0)
    assert ofiber.relative_refractive_index(1.5, 1.5) == pytest.approx(0.0)


def test_relative_refractive_index_is_close_to_the_usual_approximation():
    """Verify Delta agrees with (n_core-n_clad)/n_core for a small index step."""
    approx = (N_CORE - N_CLAD) / N_CORE
    assert ofiber.relative_refractive_index(N_CORE, N_CLAD) == pytest.approx(approx, rel=0.01)


def test_graded_index_numerical_aperture_at_the_axis_and_edge():
    """Verify the graded NA is the full NA on axis and zero at the core edge."""
    full = ofiber.numerical_aperture(N_CORE, N_CLAD)
    assert ofiber.numerical_aperture_graded_index(N_CORE, N_CLAD, 2, 0) == pytest.approx(full)
    assert ofiber.numerical_aperture_graded_index(N_CORE, N_CLAD, 2, 1) == pytest.approx(0.0)


def test_graded_index_numerical_aperture_falls_across_the_core():
    """Verify the acceptance cone narrows toward the cladding."""
    r = np.linspace(0, 1, 40)
    assert np.all(np.diff(ofiber.numerical_aperture_graded_index(N_CORE, N_CLAD, 2, r)) < 0)


def test_v_parameter_matches_its_definition():
    """Verify V = 2 pi a NA / lambda."""
    a, NA, lam = 4e-6, 0.12, 1.55e-6
    assert ofiber.V_parameter(a, NA, lam) == pytest.approx(2 * np.pi * a * NA / lam)


def test_v_parameter_is_one_at_the_cutoff_wavelength():
    """Verify the cutoff wavelength is where V reaches 2.405."""
    a, NA = 4e-6, 0.12
    lam_c = ofiber.cutoff_wavelength(a, NA)
    assert ofiber.V_parameter(a, NA, lam_c) == pytest.approx(jn_zeros(0, 1)[0])


# --------------------------------------------------------------------------
# equivalent step index
# --------------------------------------------------------------------------


@pytest.mark.parametrize("q", [1.0, 2.0, 4.0, 10.0, 100.0])
def test_the_three_esi_routines_are_mutually_consistent(q):
    """
    Verify a V rebuilt from esi_radius and esi_Delta equals esi_V_parameter.

    V is proportional to a*sqrt(2*Delta), so the radius and Delta corrections
    have to combine into exactly the V correction for the three to describe
    the same equivalent fiber.
    """
    a, Delta = 5e-6, 0.01
    scale = ofiber.esi_radius(a, q) / a * np.sqrt(ofiber.esi_Delta(Delta, q) / Delta)
    assert scale == pytest.approx(ofiber.esi_V_parameter(1.0, q), rel=1e-12)


@pytest.mark.parametrize("q", [1e6, 1e9])
def test_esi_routines_tend_to_the_step_index_values(q):
    """Verify a very large q leaves radius, Delta and V essentially unchanged."""
    assert ofiber.esi_radius(5e-6, q) == pytest.approx(5e-6, rel=1e-5)
    assert ofiber.esi_Delta(0.01, q) == pytest.approx(0.01, rel=1e-5)
    assert ofiber.esi_V_parameter(3.0, q) == pytest.approx(3.0, rel=1e-5)


@pytest.mark.parametrize(
    "q, radius, delta, v",
    [
        (1.0, 2 / 3, 3 / 4, 1 / np.sqrt(3)),
        (2.0, 3 / 4, 8 / 9, 1 / np.sqrt(2)),
    ],
)
def test_esi_closed_forms_for_triangular_and_parabolic_profiles(q, radius, delta, v):
    """Check the equivalent step index values for the two common profiles."""
    assert ofiber.esi_radius(1.0, q) == pytest.approx(radius)
    assert ofiber.esi_Delta(1.0, q) == pytest.approx(delta)
    assert ofiber.esi_V_parameter(1.0, q) == pytest.approx(v)


def test_esi_values_never_exceed_the_step_index_ones():
    """Verify grading always reduces the equivalent radius, Delta and V."""
    q = np.linspace(0.5, 50, 60)
    assert np.all(ofiber.esi_radius(1.0, q) < 1.0)
    assert np.all(ofiber.esi_Delta(1.0, q) < 1.0)
    assert np.all(ofiber.esi_V_parameter(1.0, q) < 1.0)


# --------------------------------------------------------------------------
# Fresnel reflection
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", [1.33, 1.5, 3.5])
def test_normal_incidence_reflection_is_the_same_for_both_polarizations(m):
    """Verify both polarizations give ((m-1)/(m+1))**2 at normal incidence."""
    expected = ((m - 1) / (m + 1)) ** 2
    assert ofiber.R_par(m, 0) == pytest.approx(expected)
    assert ofiber.R_per(m, 0) == pytest.approx(expected)
    assert ofiber.R_unpolarized(m, 0) == pytest.approx(expected)


@pytest.mark.parametrize("m", [1.33, 1.5, 3.5])
def test_parallel_reflection_vanishes_at_brewsters_angle(m):
    """Verify R_par is zero where tan(theta) = m."""
    assert ofiber.R_par(m, np.arctan(m)) == pytest.approx(0.0, abs=1e-15)


@pytest.mark.parametrize("m", [1.33, 1.5])
def test_perpendicular_reflection_never_vanishes(m):
    """Verify only the parallel polarization has a Brewster angle."""
    theta = np.linspace(0, np.pi / 2 - 1e-6, 200)
    assert np.all(ofiber.R_per(m, theta) > 0)


def test_reflection_approaches_unity_at_grazing_incidence():
    """Verify everything reflects at grazing incidence."""
    grazing = np.pi / 2 - 1e-7
    assert ofiber.R_par(1.5, grazing) == pytest.approx(1.0, abs=1e-4)
    assert ofiber.R_per(1.5, grazing) == pytest.approx(1.0, abs=1e-4)


def test_unpolarized_reflection_is_the_mean_of_the_two_polarizations():
    """Verify R_unpolarized averages the parallel and perpendicular results."""
    theta = np.linspace(0, 1.5, 40)
    expected = (ofiber.R_par(1.5, theta) + ofiber.R_per(1.5, theta)) / 2
    assert np.allclose(ofiber.R_unpolarized(1.5, theta), expected)


def test_reflection_stays_between_zero_and_one():
    """Verify reflectance is a fraction for a range of angles and indices."""
    theta = np.linspace(0, np.pi / 2 - 1e-6, 100)
    for m in (1.33, 1.5, 3.5):
        for R in (ofiber.R_par(m, theta), ofiber.R_per(m, theta), ofiber.R_unpolarized(m, theta)):
            assert np.all((0 <= R) & (R <= 1))


def test_an_absorbing_medium_reflects_more_than_a_transparent_one():
    """Verify a complex index raises the reflectance, as metals do."""
    assert ofiber.R_per(1.5 - 0.5j, 0.4) > ofiber.R_per(1.5, 0.4)


# --------------------------------------------------------------------------
# module surface
# --------------------------------------------------------------------------


def test_everything_in_all_is_reachable_from_the_package():
    """Verify each exported name survives the star import in __init__.py."""
    for name in ofiber.basics.__all__:
        assert hasattr(ofiber, name), "%s is in __all__ but missing from ofiber" % name


# --------------------------------------------------------------------------
# total internal reflection
# --------------------------------------------------------------------------


@pytest.mark.parametrize("m", [1 / 1.5, 1.45 / 1.46, 0.5])
def test_reflection_is_total_past_the_critical_angle(m):
    """
    Verify light beyond the critical angle is entirely reflected.

    Going from a dense medium to a rare one the relative index is below 1, so
    m**2 - sin(theta)**2 turns negative past the critical angle.  The
    discriminant is taken on the complex branch so the reflectance comes back
    as exactly 1 rather than nan.
    """
    critical = np.arcsin(m)
    theta = critical + np.linspace(1e-6, np.pi / 2 - critical - 1e-6, 25)
    assert np.allclose(ofiber.R_par(m, theta), 1.0)
    assert np.allclose(ofiber.R_per(m, theta), 1.0)
    assert np.allclose(ofiber.R_unpolarized(m, theta), 1.0)


def test_reflectance_is_finite_and_bounded_across_every_angle():
    """Verify sweeping from normal to grazing gives no nan and nothing outside [0,1]."""
    m = 1 / 1.5
    critical = np.arcsin(m)
    theta = np.linspace(0, np.pi / 2 - 1e-9, 4001)
    for R in (ofiber.R_par(m, theta), ofiber.R_per(m, theta)):
        assert np.all(np.isfinite(R))
        assert np.all((0 <= R) & (R <= 1 + 1e-12))
        assert R[theta > critical] == pytest.approx(1.0)


def test_reflectance_is_continuous_through_the_critical_angle():
    """
    Verify the reflectance has no jump where the discriminant changes branch.

    R approaches 1 with a diverging slope, so any fixed grid shows a visible
    step near the critical angle.  Continuity is checked by refining the grid
    instead: the largest step between neighbours has to shrink.
    """
    m = 1 / 1.5
    critical = np.arcsin(m)
    jumps = []
    for n in (2001, 8001, 32001):
        theta = np.linspace(critical - 0.05, critical + 0.05, n)
        jumps.append(np.max(np.abs(np.diff(ofiber.R_par(m, theta)))))
    assert jumps[0] > jumps[1] > jumps[2]


def test_reflectance_rises_monotonically_toward_the_critical_angle():
    """Verify R climbs steadily to 1 rather than wandering near the branch point."""
    m = 1 / 1.5
    theta = np.linspace(np.arctan(m), np.arcsin(m) - 1e-6, 500)
    assert np.all(np.diff(ofiber.R_par(m, theta)) > 0)
    assert np.all(np.diff(ofiber.R_per(m, theta)) > 0)


def test_a_fiber_core_totally_reflects_beyond_its_critical_angle():
    """Verify the core-cladding interface of a real fiber reflects everything."""
    m = N_CLAD / N_CORE
    just_past = ofiber.critical_angle(N_CORE, N_CLAD) + 1e-4
    assert ofiber.R_unpolarized(m, just_past) == pytest.approx(1.0)


def test_reflection_below_the_critical_angle_is_unchanged():
    """Verify the complex branch leaves ordinary partial reflection alone."""
    m = 1 / 1.5
    theta = np.linspace(0, np.arcsin(m) - 1e-3, 50)
    assert np.all(ofiber.R_par(m, theta) < 1)
    assert np.all(ofiber.R_per(m, theta) < 1)
    assert ofiber.R_par(m, 0) == pytest.approx(((m - 1) / (m + 1)) ** 2)


def test_reflection_results_are_real_valued():
    """Verify the complex discriminant does not leak into the returned power."""
    for m in (1 / 1.5, 1.5, 1.5 - 0.2j):
        for theta in (0.0, 0.8, 1.4):
            assert not isinstance(ofiber.R_par(m, theta), complex)
            assert np.isreal(ofiber.R_per(m, theta))
