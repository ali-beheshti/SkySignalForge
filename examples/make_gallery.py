#!/usr/bin/env python
"""Generate the example catalog, README figures, and animation."""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from flat_sky_signal_painter import (
    MapGeometry,
    MovingLensCatalog,
    paint_circular_masks,
    paint_moving_lens,
    paint_point_sources,
    paint_tsz_gaussian,
    single_moving_lens_template,
    truncated_nfw_profile,
    tsz_spectral_factor,
)

FIG = ROOT / "figures"


def make_catalog():
    df = pd.DataFrame({
        "Z": [0.30, 0.42, 0.55, 0.36, 0.48],
        "cNFW": [4.2, 5.0, 4.6, 3.9, 4.4],
        "Rs": [0.34, 0.28, 0.31, 0.25, 0.29],
        "rhoS": [1.00e15, 1.35e15, 1.15e15, 1.25e15, 1.10e15],
        "comovDist": [900.0, 1100.0, 1250.0, 1000.0, 1180.0],
        "vTh": [430.0, -280.0, 120.0, -350.0, 260.0],
        "vPh": [180.0, 360.0, -420.0, -160.0, -300.0],
        "RArad": np.deg2rad([-0.30, 0.28, 0.24, -0.20, 0.02]),
        "DECrad": np.deg2rad([0.20, 0.27, -0.25, -0.23, 0.02]),
        "tszGaussAmp": [1.4e-6, 1.0e-6, 1.8e-6, 0.8e-6, 1.2e-6],
        "Thetavir": np.deg2rad(np.array([5.0, 4.0, 5.5, 3.5, 4.5]) / 60.0),
        "pointAmp": [1.0, 0.7, 1.3, 0.9, 1.1],
        "dustAmp": [4.25, 3.2, 5.0, 2.8, 4.6],
    })
    df.to_csv(ROOT / "examples" / "demo_catalog.csv", index=False)
    return df


def save_map(arr, geometry, filename, title, label):
    fig, ax = plt.subplots(figsize=(6.0, 5.1))
    imshow_kwargs = {"origin": "lower", "extent": geometry.extent_deg, "cmap": "RdYlBu_r"}
    if np.nanmin(arr) < 0.0 < np.nanmax(arr):
        vmax = np.nanmax(np.abs(arr))
        imshow_kwargs.update(vmin=-vmax, vmax=vmax)
    im = ax.imshow(arr, **imshow_kwargs)
    ax.set_xlabel("RA offset [deg]")
    ax.set_ylabel("Dec offset [deg]")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=label)
    fig.tight_layout()
    fig.savefig(FIG / filename, dpi=180)
    plt.close(fig)


def plot_nfw():
    x = np.geomspace(0.011, 25.0, 700)
    y = truncated_nfw_profile(x, 4.5)
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(x, y)
    ax.axvline(4.5, linestyle="--", linewidth=1, label=r"$R_{\rm vir}/R_s$")
    ax.set_xscale("log")
    ax.set_xlabel(r"Projected radius $R/R_s$")
    ax.set_ylabel("Dimensionless deflection profile")
    ax.set_title("Truncated NFW deflection profile")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "nfw_profile.png", dpi=180)
    plt.close(fig)


def plot_tsz_spectrum():
    nu = np.linspace(20.0, 400.0, 600)
    g = tsz_spectral_factor(nu)
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(nu, g)
    ax.axhline(0.0, linewidth=1)
    ax.axvline(150.0, linestyle="--", linewidth=1, label="150 GHz")
    ax.set_xlabel("Frequency [GHz]")
    ax.set_ylabel(r"$g(x)$")
    ax.set_title("Non-relativistic tSZ spectral factor")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "tsz_spectral_factor.png", dpi=180)
    plt.close(fig)


def animate_velocity():
    geo = MapGeometry(npix=220, size_deg=0.95)
    n_frames = 32
    phase = np.linspace(0, 2 * np.pi, n_frames, endpoint=False)

    halos = [
        {
            "ra0": np.deg2rad(-0.18), "dec0": np.deg2rad(0.12),
            "dra": np.deg2rad(0.025), "ddec": np.deg2rad(0.018),
            "speed": 520.0, "phase0": 0.0,
            "redshift": 0.32, "concentration": 4.6,
            "rs_mpc": 0.32, "rho_s_msun_mpc3": 1.15e15,
            "comoving_distance_mpc": 1050.0,
        },
        {
            "ra0": np.deg2rad(0.16), "dec0": np.deg2rad(0.18),
            "dra": np.deg2rad(0.020), "ddec": np.deg2rad(-0.020),
            "speed": 430.0, "phase0": 1.4,
            "redshift": 0.45, "concentration": 5.0,
            "rs_mpc": 0.27, "rho_s_msun_mpc3": 1.30e15,
            "comoving_distance_mpc": 1300.0,
        },
        {
            "ra0": np.deg2rad(0.10), "dec0": np.deg2rad(-0.16),
            "dra": np.deg2rad(-0.022), "ddec": np.deg2rad(0.016),
            "speed": 480.0, "phase0": 3.0,
            "redshift": 0.52, "concentration": 4.2,
            "rs_mpc": 0.30, "rho_s_msun_mpc3": 1.10e15,
            "comoving_distance_mpc": 1500.0,
        },
        {
            "ra0": np.deg2rad(-0.06), "dec0": np.deg2rad(-0.06),
            "dra": np.deg2rad(0.015), "ddec": np.deg2rad(0.022),
            "speed": 370.0, "phase0": 4.2,
            "redshift": 0.38, "concentration": 3.9,
            "rs_mpc": 0.25, "rho_s_msun_mpc3": 1.20e15,
            "comoving_distance_mpc": 1180.0,
        },
    ]

    maps = []
    marker_positions = []
    for t in phase:
        total = np.zeros((geo.npix, geo.npix), dtype=float)
        frame_pos = []
        for h in halos:
            ra = h["ra0"] + h["dra"] * np.cos(t + h["phase0"])
            dec = h["dec0"] + h["ddec"] * np.sin(t + 0.8 * h["phase0"])
            angle = t + h["phase0"]
            speed = h["speed"] * (0.92 + 0.08 * np.cos(t + h["phase0"]))
            halo = {
                "ra_rad": ra,
                "dec_rad": dec,
                "redshift": h["redshift"],
                "concentration": h["concentration"],
                "rs_mpc": h["rs_mpc"],
                "rho_s_msun_mpc3": h["rho_s_msun_mpc3"],
                "comoving_distance_mpc": h["comoving_distance_mpc"],
                "v_theta_kms": speed * np.sin(angle),
                "v_phi_kms": speed * np.cos(angle),
            }
            total += single_moving_lens_template(halo, geo)
            frame_pos.append((np.rad2deg(ra), np.rad2deg(dec)))
        maps.append(total)
        marker_positions.append(np.array(frame_pos))

    vmax = max(np.nanmax(np.abs(m)) for m in maps)
    fig, ax = plt.subplots(figsize=(6.0, 5.2))
    im = ax.imshow(
        maps[0],
        origin="lower",
        extent=geo.extent_deg,
        cmap="RdYlBu_r",
        vmin=-vmax,
        vmax=vmax,
    )
    ax.set_xlabel("RA offset [deg]")
    ax.set_ylabel("Dec offset [deg]")
    title = ax.set_title("Multi-halo moving-lens evolution")
    scat = ax.scatter(
        marker_positions[0][:, 0],
        marker_positions[0][:, 1],
        s=22,
        facecolors="none",
        edgecolors="k",
        linewidths=0.8,
    )
    fig.colorbar(im, ax=ax, label=r"$\Delta T$ [$\mu$K]")
    fig.tight_layout()

    def update(i):
        im.set_data(maps[i])
        scat.set_offsets(marker_positions[i])
        title.set_text(f"Multi-halo moving-lens evolution — frame {i + 1}/{n_frames}")
        return im, scat, title

    ani = FuncAnimation(fig, update, frames=n_frames, interval=120, blit=False)
    ani.save(FIG / "velocity_rotation.gif", writer=PillowWriter(fps=8))
    plt.close(fig)


def main():
    FIG.mkdir(exist_ok=True)
    df = make_catalog()

    ml_geo = MapGeometry(npix=256, size_deg=0.45)
    one = MovingLensCatalog(
        ra_rad=[0.0], dec_rad=[0.0],
        redshift=[0.35], concentration=[4.5], rs_mpc=[0.32],
        rho_s_msun_mpc3=[1.1e15], comoving_distance_mpc=[1200.0],
        v_theta_kms=[450.0], v_phi_kms=[250.0],
    )
    ml = paint_moving_lens(one, ml_geo, show_progress=False)
    save_map(ml, ml_geo, "moving_lens_dipole.png",
             "Single-halo moving-lens temperature dipole", r"$\Delta T$ [$\mu$K]")

    tsz_geo = MapGeometry(npix=256, size_deg=0.45)
    tsz = paint_tsz_gaussian(
        [0.0], [0.0], [1.4e-6], tsz_geo,
        frequency_ghz=150.0, show_progress=False,
    )
    save_map(tsz, tsz_geo, "tsz_gaussian.png",
             "Gaussian tSZ profile at 150 GHz", r"$\Delta T_{\rm tSZ}$ [$\mu$K]")

    wide = MapGeometry(npix=260, size_deg=1.2)
    pts = paint_point_sources(
        df["RArad"], df["DECrad"], wide,
        amplitudes=df["pointAmp"], show_progress=False,
    )
    save_map(pts, wide, "point_sources.png",
             "Example weighted point-source map", "source amplitude")

    mask = paint_circular_masks(
        df["RArad"], df["DECrad"], df["Thetavir"], wide,
        combine="union", show_progress=False,
    )
    save_map(mask, wide, "halo_mask.png",
             "Circular halo mask", "mask")

    plot_nfw()
    plot_tsz_spectrum()
    animate_velocity()
    print("Wrote example catalog and gallery.")


if __name__ == "__main__":
    main()
