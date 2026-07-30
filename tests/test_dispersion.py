# pylint: disable=invalid-name
"""
Tests for the ofiber.dispersion module.

The anchors are Ghatak & Thyagarajan tables 10.2 and 10.3, which
`5-Dispersion.ipynb` reproduces and prints.

One caution about those printed tables: the D_m row was computed before the
missing accumulation in `_d2_sellmeier` was fixed, so only the rows that do
not depend on the second derivative are used as published values here, namely
n_core, n_clad, V, V*d2(bV)/dV2 and D_w.  The material dispersion is instead
checked against its own definition, -λ/c times d2n, which test_refraction.py
already pins against finite differences of an independent Sellmeier.

This module is thin glue, so most of what is worth testing is that it wires
`refraction`, `cylinder_step` and `basics` together correctly.
"""

import numpy as np
import pytest
import scipy.constants
import ofiber

LAMBDAS = np.array([1100, 1300, 1560]) * 1e-9

# Ghatak table 10.2: 2.6% GeO2 core, silica cladding, 4.1 micron core radius
CORE_10p2 = ofiber.doped_glass(0.0260)
R_CORE_10p2 = 4.1e-6
N_CORE_10p2 = [1.45308, 1.45080, 1.44781]
N_CLAD_10p2 = [1.44920, 1.44692, 1.44390]
DW_10p2 = [-1.83, -3.76, -5.96]

# Ghatak table 10.3: 7.35% GeO2 core, 2.3 micron core radius
CORE_10p3 = ofiber.doped_glass(0.0735)
R_CORE_10p3 = 2.3e-6
V_10p3 = [2.344, 1.983, 1.656]
VBPP_10p3 = [0.224, 0.477, 0.843]
DW_10p3 = [-7.34, -13.27, -19.64]

CLAD = ofiber.doped_glass(0)
PS_PER_KM_NM = 1e6  # multiply [s/m**2] by this


def _dw(core, r_core, **kwargs):
    """
    Evaluate the waveguide dispersion at LAMBDAS in ps/(km nm).

    Args:
        core: Sellmeier coefficients of the core glass
        r_core: core radius [m]
        **kwargs: passed through to Waveguide_Dispersion

    Returns:
        waveguide dispersion at each wavelength [ps/(km nm)]
    """
    n_core = ofiber.n(core, LAMBDAS)
    n_clad = ofiber.n(CLAD, LAMBDAS)
    return np.array(
        [
            ofiber.Waveguide_Dispersion(n_core[i], n_clad[i], r_core, LAMBDAS[i], **kwargs) * PS_PER_KM_NM
            for i in range(len(LAMBDAS))
        ]
    )


# --------------------------------------------------------------------------
# material dispersion
# --------------------------------------------------------------------------


def test_material_dispersion_follows_its_definition():
    """Verify D_m is -lambda/c times the second derivative of the index."""
    lam = np.linspace(1.0e-6, 1.7e-6, 25)
    expected = -lam * ofiber.d2n(CORE_10p2, lam) / scipy.constants.speed_of_light
    assert np.allclose(ofiber.Material_Dispersion(CORE_10p2, lam), expected, rtol=1e-14)


def test_silica_material_dispersion_vanishes_near_1272nm():
    """Verify pure silica has zero material dispersion at the documented wavelength."""
    lam = np.linspace(1.2e-6, 1.35e-6, 15001)
    Dm = ofiber.Material_Dispersion(CLAD, lam)
    crossings = np.where(np.diff(np.sign(Dm)) != 0)[0]
    assert len(crossings) == 1
    assert lam[crossings[0]] * 1e9 == pytest.approx(1272.75, abs=0.5)


def test_germania_moves_the_material_zero_to_longer_wavelengths():
    """Verify doping the core with GeO2 pushes the zero of D_m upward."""
    lam = np.linspace(1.2e-6, 1.45e-6, 25001)
    zeros = []
    for x in (0.0, 0.05, 0.10, 0.20):
        Dm = ofiber.Material_Dispersion(ofiber.doped_glass(x), lam)
        zeros.append(lam[np.where(np.diff(np.sign(Dm)) != 0)[0][0]])
    assert np.all(np.diff(zeros) > 0)


def test_material_dispersion_changes_sign_across_the_zero():
    """Verify D_m is negative below the zero and positive above it."""
    assert ofiber.Material_Dispersion(CLAD, 1.1e-6) < 0
    assert ofiber.Material_Dispersion(CLAD, 1.5e-6) > 0


# --------------------------------------------------------------------------
# waveguide dispersion
# --------------------------------------------------------------------------


def test_indices_match_ghatak_table_10p2():
    """Check the core and cladding indices the notebook prints for table 10.2."""
    assert np.allclose(ofiber.n(CORE_10p2, LAMBDAS), N_CORE_10p2, atol=5e-6)
    assert np.allclose(ofiber.n(CLAD, LAMBDAS), N_CLAD_10p2, atol=5e-6)


def test_waveguide_dispersion_matches_ghatak_table_10p2():
    """Check D_w against the values published for table 10.2."""
    assert np.allclose(_dw(CORE_10p2, R_CORE_10p2), DW_10p2, atol=5e-3)


def test_waveguide_dispersion_matches_ghatak_table_10p3():
    """Check D_w against the values published for the smaller, more doped fiber."""
    assert np.allclose(_dw(CORE_10p3, R_CORE_10p3), DW_10p3, atol=5e-3)


def test_v_parameter_and_dispersion_term_match_table_10p3():
    """Check the V and V*d2(bV)/dV2 rows that feed D_w."""
    n_core, n_clad = ofiber.n(CORE_10p3, LAMBDAS), ofiber.n(CLAD, LAMBDAS)
    V = ofiber.V_parameter(R_CORE_10p3, ofiber.numerical_aperture(n_core, n_clad), LAMBDAS)
    assert np.allclose(V, V_10p3, atol=5e-4)
    assert np.allclose(ofiber.V_d2bV_by_V(V, 0), VBPP_10p3, atol=5e-4)


def test_waveguide_dispersion_matches_its_closed_form():
    """Verify D_w is -n_clad*esi_Delta/(c*lambda) times V*d2(bV)/dV2."""
    n_core, n_clad, lam, r = 1.4531, 1.4492, 1.3e-6, 4.1e-6
    Delta = ofiber.relative_refractive_index(n_core, n_clad)
    V = ofiber.V_parameter(r, ofiber.numerical_aperture(n_core, n_clad), lam)
    q = 1e20
    expected = (
        -n_clad
        * ofiber.esi_Delta(Delta, q)
        / scipy.constants.speed_of_light
        / lam
        * ofiber.V_d2bV_by_V(ofiber.esi_V_parameter(V, q), 0)
    )
    assert ofiber.Waveguide_Dispersion(n_core, n_clad, r, lam) == pytest.approx(expected, rel=1e-12, abs=0)


def test_waveguide_dispersion_is_negative_for_an_ordinary_fiber():
    """Verify D_w pulls the total dispersion zero to longer wavelengths."""
    assert np.all(_dw(CORE_10p2, R_CORE_10p2) < 0)
    assert np.all(_dw(CORE_10p3, R_CORE_10p3) < 0)


def test_waveguide_dispersion_grows_in_magnitude_with_wavelength():
    """Verify a longer wavelength is less well confined and disperses more."""
    assert np.all(np.diff(_dw(CORE_10p2, R_CORE_10p2)) < 0)


def test_a_smaller_core_disperses_more():
    """Verify tighter confinement strengthens the waveguide term."""
    n_core, n_clad = ofiber.n(CORE_10p3, 1.3e-6), ofiber.n(CLAD, 1.3e-6)
    small = ofiber.Waveguide_Dispersion(n_core, n_clad, 2.0e-6, 1.3e-6)
    large = ofiber.Waveguide_Dispersion(n_core, n_clad, 5.0e-6, 1.3e-6)
    assert abs(small) > abs(large)


def test_waveguide_dispersion_peaks_partway_through_the_single_mode_range():
    """
    Verify D_w is not monotonic in the index step but has a clear extremum.

    D_w goes as Delta times V*d2(bV)/dV2, and raising the index step raises
    both Delta and V while the dispersion term falls with V.  The product
    peaks near V=1.4 rather than growing without bound.
    """
    n_clad, r, lam = 1.4492, 4.1e-6, 1.3e-6
    n_core = np.linspace(1.4500, 1.4540, 60)
    V = ofiber.V_parameter(r, ofiber.numerical_aperture(n_core, n_clad), lam)
    assert np.all(V < 2.405), "stay inside the single mode range"
    Dw = np.array([ofiber.Waveguide_Dispersion(nc, n_clad, r, lam) for nc in n_core])
    peak = np.argmax(np.abs(Dw))
    assert 0 < peak < len(n_core) - 1
    assert V[peak] == pytest.approx(1.4, abs=0.2)


# --------------------------------------------------------------------------
# the graded index parameter
# --------------------------------------------------------------------------


def test_the_default_q_behaves_as_a_step_index_fiber():
    """Verify the huge default q is indistinguishable from an even larger one."""
    args = (1.4531, 1.4492, 4.1e-6, 1.3e-6)
    assert ofiber.Waveguide_Dispersion(*args) == pytest.approx(
        ofiber.Waveguide_Dispersion(*args, q=1e30), rel=1e-9, abs=0
    )


def test_infinite_q_is_not_a_usable_step_index_sentinel():
    """
    Verify np.inf gives nan, which is why the default is a large finite q.

    esi_Delta multiplies q*(2+q)/(1+q)**2, and that is nan at infinity.
    """
    assert np.isnan(ofiber.esi_Delta(0.0027, np.inf))
    assert np.isnan(ofiber.Waveguide_Dispersion(1.4531, 1.4492, 4.1e-6, 1.3e-6, q=np.inf))


def test_grading_the_profile_strengthens_the_waveguide_dispersion():
    """
    Verify a graded profile disperses more than the step index equivalent.

    Lowering q shrinks both the equivalent V and the equivalent Delta, but
    V*d2(bV)/dV2 climbs steeply as V falls, and that wins.  A triangular
    profile disperses roughly three times as strongly as a step index one.
    """
    args = (1.4531, 1.4492, 4.1e-6, 1.3e-6)
    magnitudes = [abs(ofiber.Waveguide_Dispersion(*args, q=q)) for q in (1e20, 100.0, 10.0, 4.0, 2.0, 1.0)]
    assert np.all(np.diff(magnitudes) > 0)
    assert magnitudes[-1] / magnitudes[0] == pytest.approx(2.9, abs=0.2)


def test_larger_q_approaches_the_step_index_result():
    """Verify the graded result tends to the step index one as q grows."""
    args = (1.4531, 1.4492, 4.1e-6, 1.3e-6)
    step = ofiber.Waveguide_Dispersion(*args)
    gaps = [abs(ofiber.Waveguide_Dispersion(*args, q=q) - step) for q in (2.0, 10.0, 100.0, 1000.0)]
    assert np.all(np.diff(gaps) < 0)


# --------------------------------------------------------------------------
# the Marcuse approximation
# --------------------------------------------------------------------------


def test_the_approximation_tracks_the_exact_waveguide_dispersion():
    """
    Verify approx=True stays within the documented few percent of the exact result.

    The gap peaks near V=2.44, where the quadratic is about 6% low, so the
    tolerance here reflects the worst case rather than the typical one.
    """
    exact = _dw(CORE_10p2, R_CORE_10p2)
    approx = _dw(CORE_10p2, R_CORE_10p2, approx=True)
    assert np.allclose(approx, exact, rtol=0.07)


def test_the_approximation_is_close_where_it_is_documented_to_be():
    """Verify the two agree to about 1% where V sits in the accurate range."""
    n_core, n_clad = ofiber.n(CORE_10p2, 1.4e-6), ofiber.n(CLAD, 1.4e-6)
    args = (n_core, n_clad, 4.1e-6, 1.4e-6)
    V = ofiber.V_parameter(4.1e-6, ofiber.numerical_aperture(n_core, n_clad), 1.4e-6)
    assert 1.35 < V < 2.08
    assert ofiber.Waveguide_Dispersion(*args, approx=True) == pytest.approx(
        ofiber.Waveguide_Dispersion(*args), rel=0.01
    )


# --------------------------------------------------------------------------
# total dispersion
# --------------------------------------------------------------------------


def test_dispersion_returns_the_two_terms_it_advertises():
    """Verify Dispersion is exactly its two components, computed separately."""
    Dm, Dw = ofiber.Dispersion(CORE_10p2, ofiber.n(CLAD, 1.3e-6), R_CORE_10p2, 1.3e-6)
    assert Dm == pytest.approx(ofiber.Material_Dispersion(CORE_10p2, 1.3e-6), rel=1e-14, abs=0)
    n_core = ofiber.n(CORE_10p2, 1.3e-6)
    expected = ofiber.Waveguide_Dispersion(n_core, ofiber.n(CLAD, 1.3e-6), R_CORE_10p2, 1.3e-6)
    assert Dw == pytest.approx(expected, rel=1e-14, abs=0)


def test_dispersion_passes_q_and_approx_through():
    """Verify the optional arguments reach Waveguide_Dispersion."""
    n_clad = ofiber.n(CLAD, 1.3e-6)
    _, plain = ofiber.Dispersion(CORE_10p2, n_clad, R_CORE_10p2, 1.3e-6)
    _, graded = ofiber.Dispersion(CORE_10p2, n_clad, R_CORE_10p2, 1.3e-6, q=2.0)
    _, approximated = ofiber.Dispersion(CORE_10p2, n_clad, R_CORE_10p2, 1.3e-6, approx=True)
    assert graded != plain
    assert approximated != plain


def test_the_waveguide_term_moves_the_zero_of_the_total_dispersion():
    """
    Verify adding the waveguide term pushes the zero past the material one.

    This is the whole point of the module: a standard fiber has its material
    zero near 1272 nm but disperses zero closer to 1310 nm.
    """
    lam = np.linspace(1.2e-6, 1.45e-6, 4001)
    n_clad = ofiber.n(CLAD, lam)
    Dm, Dw = ofiber.Dispersion(CORE_10p2, n_clad, R_CORE_10p2, lam)
    material_zero = lam[np.where(np.diff(np.sign(Dm)) != 0)[0][0]]
    total_zero = lam[np.where(np.diff(np.sign(Dm + Dw)) != 0)[0][0]]
    assert total_zero > material_zero
    assert 1290 < total_zero * 1e9 < 1340


def test_total_dispersion_accepts_an_array_of_wavelengths():
    """Verify the whole chain vectorizes over wavelength."""
    Dm, Dw = ofiber.Dispersion(CORE_10p2, ofiber.n(CLAD, LAMBDAS), R_CORE_10p2, LAMBDAS)
    assert Dm.shape == LAMBDAS.shape
    assert Dw.shape == LAMBDAS.shape
    assert np.all(np.isfinite(Dm)) and np.all(np.isfinite(Dw))


# --------------------------------------------------------------------------
# module surface
# --------------------------------------------------------------------------


def test_everything_in_all_is_reachable_from_the_package():
    """Verify each exported name survives the star import in __init__.py."""
    for name in ofiber.dispersion.__all__:
        assert hasattr(ofiber, name), "%s is in __all__ but missing from ofiber" % name
