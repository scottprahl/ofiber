#!/usr/bin/env python3
"""Generate README/docs SVG images from notebook plotting code.

Running this script writes all SVG images into the directory where the script is
executed (the current working directory).
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

# Ensure Matplotlib cache/config is writable in sandboxed/headless runs.
MPLCONFIGDIR = Path(tempfile.gettempdir()) / "ofiber-mplconfig"
MPLCONFIGDIR.mkdir(parents=True, exist_ok=True)
os.environ["MPLCONFIGDIR"] = str(MPLCONFIGDIR)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy

# Allow importing local package without installation.
REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

import ofiber


def dt_mixture(x: float, *args) -> float:
    """
    Total dispersion for (x)GeO2 : (1-x)SiO2 glass fiber.

    The cladding is assumed to be pure SiO2 surrounding a GeO2 doped core
    with specified radius.

    The returned dispersion is in ps/km/nm.
    """
    r_core = args[0]  # core radius [m]
    wavelength = args[1]  # core wavelength [m]

    cladding = ofiber.doped_glass(0)
    core = ofiber.doped_glass(x)

    n_clad = ofiber.n(cladding, wavelength)

    dm, dw = ofiber.Dispersion(core, n_clad, r_core, wavelength)
    return (dm + dw) * 1e6


def dB(x):
    return 10 * np.log10(x)


def _save(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(path)
    plt.close()
    return path


def make_dispersion_svg(output_dir: Path) -> Path:
    # Ghatak Section 3.A.1, page 85.
    c = 3e8  # [m/s]

    a = 4.1e-6  # [m] Fiber radius
    delta = 0.0027  # Fractional change in the index of refraction

    start = 1100
    finish = 1700
    resolution = 10

    wavelength = np.arange(start, finish, resolution) * 1e-9

    glass = ofiber.doped_glass(0)
    n1 = ofiber.n(glass, wavelength)
    d2n = ofiber.d2n(glass, wavelength)

    V = (2 * np.pi / wavelength) * a * n1 * np.sqrt(2 * delta)

    # Material Dispersion
    M = -(wavelength / c) * d2n

    # Waveguide dispersion
    Dw = -n1 * (1 + delta) * delta / c / wavelength * (0.080 + 0.549 * (2.834 - V) ** 2)

    ps_nm_km = 1e-12 / 1e-9 / 1e3
    plt.subplots(1, 1, figsize=(8, 4.5))
    plt.plot(wavelength * 1e9, M / ps_nm_km, color="blue")

    plt.plot(wavelength * 1e9, Dw / ps_nm_km, color="orange")
    plt.plot(wavelength * 1e9, (M + Dw) / ps_nm_km, color="red")
    plt.axhline(0, color="black", linewidth=0.5)
    plt.xlabel("Wavelength [microns]")
    plt.ylabel("Fiber Dispersion [ps/km/nm]")
    plt.title(r"SiO$_2$ core, Diameter=%.1f$\mu$m, $\Delta$=%.4f" % (a * 2e6, delta))
    plt.xlim(1100, 1700)
    plt.text(1510, 10, "Total Dispersion", color="red")
    plt.text(1600, 26, "Material Dispersion", color="blue", ha="right")
    plt.text(1400, -10, "Waveguide Dispersion", color="orange")
    plt.axvline(1315, color="black", linewidth=0.5)
    plt.text(1300, 10, r" zero dispersion $\lambda$", rotation=90)
    return _save(output_dir / "dispersion.svg")


def make_planarwaveguide_svg(output_dir: Path) -> Path:
    V = 15
    d = 1
    x = np.linspace(-1, 1, 100)

    plt.figure(figsize=(8, 4.5))
    m = 0
    plt.plot(x, ofiber.TE_field(V, d, x, m), color="blue")
    plt.text(-0.42, 0.7, "m=%d" % m, color="blue")

    m = 1
    plt.plot(x, ofiber.TE_field(V, d, x, m), color="red")
    plt.text(-0.10, -0.7, "m=%d" % m, color="red")

    m = 2
    plt.plot(x, ofiber.TE_field(V, d, x, m), color="orange")
    plt.text(0.25, -0.5, "m=%d" % m, color="orange")

    plt.plot(x, np.exp(-(x**2) / 0.01), ":b")
    plt.plot([-1, 1], [0, 0], "k")

    plt.axvspan(-0.5, 0.5, color="cyan", alpha=0.5)
    plt.text(-0.55, 0.25, "waveguide bottom", rotation=90, fontsize=8)
    plt.text(0.5, -1, "waveguide top", rotation=90, fontsize=8)

    plt.xlabel("Position (x/d)")
    plt.ylabel("$|E_y(x)|^2$ [Normalized]")
    plt.title("Modal Fields in symmetric planar waveguide V=%.2f" % V)
    return _save(output_dir / "planarwaveguide.svg")


def make_stepindexmodes_svg(output_dir: Path) -> Path:
    V = np.linspace(0.01, 10, 50)
    b = np.empty_like(V)

    plt.figure(figsize=(8, 4.5))
    for ell in range(3):
        for em in range(1, 4):
            for i in range(len(V)):
                b[i] = ofiber.LP_mode_value(V[i], ell, em)
            plt.plot(V, b)
            plt.annotate(r" LP$_{%d%d}$" % (ell, em), xy=(10, b[-1]), va="center")

    plt.xlim(0, 12)
    plt.xlabel("V-parameter")
    plt.ylabel("Propagation Constant b")
    plt.title("Modes in a step-index Fiber")
    return _save(output_dir / "stepindexmodes.svg")


def make_modeirradiance_svg(output_dir: Path) -> Path:
    plt.subplots(figsize=(8, 4.5))
    r_over_a = np.linspace(-1.2, 1.2, 100)
    V = 8
    for ell in range(3):
        for em in range(2, 4):
            b = ofiber.LP_mode_value(V, ell, em)
            if b is None:
                continue
            irradiance = ofiber.LP_radial_irradiance(V, b, ell, r_over_a)
            i = np.argmax(irradiance)

            plt.plot(r_over_a, irradiance)
            plt.annotate(
                r" LP$_{%d%d}$" % (ell, em),
                xy=(r_over_a[i], irradiance[i]),
                va="bottom",
                ha="center",
            )

    plt.xlabel("r / r$_{core}$")
    plt.ylabel("Normalized Irradiance")
    plt.plot([-1, -1], [0, 10], ":k")
    plt.plot([1, 1], [0, 10], ":k")
    plt.axvspan(-1, 1, color="cyan", alpha=0.5)
    plt.title("Step-index fiber, V=%.2f" % V)
    plt.ylim(0, 10)
    plt.xlim(-1.2, 1.2)
    return _save(output_dir / "modeirradiance.svg")


def make_internalmodes_svg(output_dir: Path) -> Path:
    V = 8

    r_over_a = np.linspace(0, 1.5, 50)
    phi = np.linspace(0, 2 * np.pi, 40)
    np.meshgrid(r_over_a, phi)  # keep parity with notebook setup

    clevs = np.linspace(-4, 4, 9)
    dlevs = np.linspace(0, 8, 9)

    fig = plt.figure(figsize=(8, 3.5))

    ell = 4
    em = 1

    phi_field = np.cos(ell * phi)

    b = ofiber.LP_mode_value(V, ell, em)
    r_field = ofiber.LP_radial_field(V, b, ell, r_over_a)
    R_FIELD, PHI_FIELD = np.meshgrid(phi_field, r_field)
    Z = R_FIELD * PHI_FIELD

    ax = plt.subplot(1, 2, 1, polar=True)
    cax = ax.contourf(phi, r_over_a, Z, levels=clevs)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    fig.colorbar(cax, ax=ax)
    ax.set_title(r"LP$_{%d%d}$ Field" % (ell, em))

    ax = plt.subplot(1, 2, 2, polar=True)
    cax = ax.contourf(phi, r_over_a, Z * Z, levels=dlevs)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    fig.colorbar(cax, ax=ax)
    ax.set_title(r"LP$_{%d%d}$ Irradiance" % (ell, em))

    return _save(output_dir / "internalmodes.svg")


def make_fiberdesign_svg(output_dir: Path) -> Path:
    r_core = 2.1e-6  # meters
    design_wavelength = 1550e-9  # meters

    cladding = ofiber.doped_glass(0)
    n_clad = ofiber.n(cladding, design_wavelength)

    x = np.linspace(0.05, 0.9)
    Dm = np.empty_like(x)
    Dw = np.empty_like(x)

    for i, xx in enumerate(x):
        core_glass = ofiber.doped_glass(xx)
        n_core = ofiber.n(core_glass, design_wavelength)
        Dm[i] = ofiber.Material_Dispersion(core_glass, design_wavelength) * 1e6  # ps/(km nm)
        Dw[i] = ofiber.Waveguide_Dispersion(n_core, n_clad, r_core, design_wavelength) * 1e6  # ps/(km nm)

    # plot waveguide, material, and total dispersion
    plt.figure(figsize=(8, 4.5))
    plt.plot(x * 100, Dw, "red")
    plt.text(0, Dw[0], " $D_w$", va="center", ha="left", color="red")

    plt.plot(x * 100, Dm, "blue")
    plt.text(0, Dm[0], " $D_m$", va="center", ha="left", color="blue")

    plt.plot(x * 100, Dm + Dw, "orange")
    plt.text(0, Dm[0] + Dw[0], " $D_t$", va="center", ha="left", color="orange")

    plt.axhline(0, color="black", ls=":")

    # find the zero dispersion concentration
    xx = scipy.optimize.brentq(dt_mixture, 0.1, 0.9, args=(r_core, design_wavelength))

    core = ofiber.doped_glass(xx)
    n_core = ofiber.n(core, design_wavelength)
    V = 2 * np.pi * r_core / design_wavelength * np.sqrt(n_core**2 - n_clad**2)
    print("V=%.3f" % V)

    # label the zero dispersion concentration
    Dtxx = dt_mixture(xx, r_core, design_wavelength)
    plt.plot([xx * 100], [Dtxx], "ok")
    plt.text(xx * 100, Dtxx, " %.2f%% GeO$_2$" % (xx * 100), va="bottom")

    plt.xlim(-5, 100)
    plt.xlabel("Molar Percentage of GeO$_2$")
    plt.ylabel("Dispersion [ps/(km nm)]")
    plt.title("GeO$_2$ added to SiO$_2$")
    return _save(output_dir / "fiberdesign.svg")


def make_farfieldirradiance_svg(output_dir: Path) -> Path:
    k = 1
    ell = 0  # azimuthal mode ell
    a = 20.3  # k=1 so ka=20.3
    V = 2.4
    b = 0.53303  # this value corresponds to the LP_01 mode
    lambda0 = 2 * np.pi  # so that k=1 and ka=a=20.3
    r = 100 * lambda0  # 100 wavelengths from the fiber face

    theta = np.radians(np.linspace(-90, 90, 501))

    # calculate the normalized irradiance
    IFF = ofiber.FF_polar_irradiance_x(r, theta, ell, lambda0, a, V, b)
    F = IFF / np.max(IFF)

    title = "$LP_{01}$ x-polarized far-field irradiance (normalized to max)\n"
    title += r"b=%.5f, V=%.2f, ka=%.2f" % (b, V, a)
    rlabel = r"$\log_{10}|E_x(\theta)/E_x(0)|^2 (dB)$"

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(7, 7))
    ax.plot(theta, dB(F), "b")
    ax.set_theta_zero_location("N")
    ax.set_thetalim(np.pi / 2, -np.pi / 2)
    ax.set_ylim(-100, 0)
    ax.grid(True)
    ax.set_title(title, va="bottom", pad=-40)
    plt.figtext(0.23, 0.23, rlabel, ha="left", va="bottom")
    return _save(output_dir / "farfieldirradiance.svg")


def make_theta01_svg(output_dir: Path) -> Path:
    ell = 0
    em = 1
    n_core = 1.589
    n_cladding = 1.48
    lambda0 = 0.550  # center wavelength over visible spectrum 0.400-0.700 um
    a = np.linspace(2, 15, num=131)  # core radii from 2-15 microns

    Delta = ofiber.relative_refractive_index(n_core, n_cladding)
    NA = ofiber.numerical_aperture_from_Delta(n_core, Delta)
    V = ofiber.V_parameter(a, NA, lambda0)

    b = np.zeros_like(V)
    b[0] = ofiber.LP_mode_value(V[0], ell, em)

    for i in range(1, len(V)):
        b[i] = scipy.optimize.root_scalar(
            f=ofiber.cylinder_step._cyl_mode_eqn,
            args=(V[i], 0),
            bracket=(b[i - 1], 1),
        ).root

    theta = np.radians(np.linspace(0, 90, 1801))
    I_a_theta = np.zeros((len(theta), len(a)))

    r = 100
    ell = 0
    lambda0 = 0.55
    for i in range(len(V)):
        I_a_theta[:, i] = ofiber.FF_polar_irradiance_x(r, theta, ell, lambda0, a[i], V[i], b[i])

    firstMin = np.zeros_like(a)

    for i in range(a.size):
        peaks, _ = scipy.signal.find_peaks(-I_a_theta[:, i])
        firstMin[i] = theta[peaks[0]]

    FF_node = ofiber.FF_node_polar_angle(V, ell, em)

    k = 2 * np.pi / lambda0  # free-space wavenumber
    firstMinDeg = np.degrees(firstMin)  # convert the first min angle peak to degrees
    FF_node_deg = np.degrees(np.arcsin(FF_node / (k * a)))  # convert k*a*sin theta to theta in degrees

    plt.figure(figsize=(8, 9))
    plt.subplot(2, 1, 1)
    # plt.xlabel(r'Fiber core radius $a$ ($\mu$m)')
    plt.gca().set_xticklabels([])
    plt.ylabel(r"Fundamental mode emission angle $\theta_{01}$ ($^{\circ}$)")
    plt.title(r"Emission angle $\theta_{01}$ vs fiber core radius")
    plt.plot(a, firstMinDeg, "ob", label=r"Angle $\theta_{01}$ using peak finding", markersize=1)
    plt.plot(a, FF_node_deg, label=r"Angle $\theta_{01}$ using zero bracketing", lw=0.5, color="red")
    plt.grid()
    plt.legend()

    # plt.figure(figsize=(8,4.5))
    plt.subplot(2, 1, 2)
    plt.xlabel(r"Fiber core radius ($\mu$m)")
    plt.ylabel("Error in fitted emission angle (degrees)")
    # plt.title(r'Emission angle $\theta_{01}$ vs core radius $a$')
    plt.plot(a, firstMinDeg - FF_node_deg, "ob", lw=0.5, markersize=2)
    plt.axhline(0)
    plt.ylim(-0.05, 0.05)
    plt.grid(True)
    return _save(output_dir / "theta01.svg")


def generate_all(output_dir: Path | None = None) -> list[Path]:
    output_dir = Path.cwd() if output_dir is None else Path(output_dir)
    generated = []

    generators = [
        make_dispersion_svg,
        make_planarwaveguide_svg,
        make_stepindexmodes_svg,
        make_modeirradiance_svg,
        make_internalmodes_svg,
        make_fiberdesign_svg,
        make_farfieldirradiance_svg,
        make_theta01_svg,
    ]

    for make_image in generators:
        image_path = make_image(output_dir)
        print(f"wrote {image_path}")
        generated.append(image_path)

    return generated


def main() -> int:
    generate_all()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
