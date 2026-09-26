"""Command-line interface for the signal painters."""

from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .catalog import MovingLensCatalog, PositionCatalog
from .geometry import MapGeometry
from .masks import paint_circular_masks
from .moving_lens import paint_moving_lens
from .sources import paint_dust_test_sources, paint_point_sources
from .tsz import paint_tsz_gaussian


def _save(array, output, png, geometry, label):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.save(output, array)
    print(f"Saved {output}")
    if png:
        png = Path(png)
        png.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(6, 5))
        imshow_kwargs = {"origin": "lower", "extent": geometry.extent_deg, "cmap": "RdYlBu_r"}
        if np.nanmin(array) < 0.0 < np.nanmax(array):
            vmax = np.nanmax(np.abs(array))
            imshow_kwargs.update(vmin=-vmax, vmax=vmax)
        im = ax.imshow(array, **imshow_kwargs)
        ax.set_xlabel("RA offset [deg]")
        ax.set_ylabel("Dec offset [deg]")
        fig.colorbar(im, ax=ax, label=label)
        fig.tight_layout()
        fig.savefig(png, dpi=180)
        plt.close(fig)
        print(f"Saved {png}")


def main():
    p = argparse.ArgumentParser(description="Paint flat-sky astrophysical signal maps.")
    sub = p.add_subparsers(dest="command", required=True)

    for name in ["moving-lens", "tsz", "points", "dust", "mask"]:
        sp = sub.add_parser(name)
        sp.add_argument("catalog")
        sp.add_argument("--npix", type=int, default=256)
        sp.add_argument("--size-deg", type=float, default=4.0)
        sp.add_argument("--output", default=f"{name}.npy")
        sp.add_argument("--png", default=None)
        sp.add_argument("--no-progress", action="store_true")

    sub.choices["tsz"].add_argument("--frequency-ghz", type=float, default=150.0)
    sub.choices["mask"].add_argument("--union", action="store_true")

    args = p.parse_args()
    df = pd.read_csv(args.catalog)
    geometry = MapGeometry(npix=args.npix, size_deg=args.size_deg)
    show = not args.no_progress

    if args.command == "moving-lens":
        arr = paint_moving_lens(MovingLensCatalog.from_dataframe(df), geometry, show_progress=show)
        label = r"$\Delta T$ [$\mu$K]"
    else:
        pos = PositionCatalog.from_dataframe(df)
        if args.command == "tsz":
            arr = paint_tsz_gaussian(
                pos.ra_rad, pos.dec_rad, df["tszGaussAmp"].to_numpy(), geometry,
                frequency_ghz=args.frequency_ghz, show_progress=show,
            )
            label = r"$\Delta T_{\rm tSZ}$ [$\mu$K]"
        elif args.command == "points":
            amp = df["pointAmp"].to_numpy() if "pointAmp" in df else None
            arr = paint_point_sources(pos.ra_rad, pos.dec_rad, geometry, amplitudes=amp, show_progress=show)
            label = "source amplitude"
        elif args.command == "dust":
            amp = df["dustAmp"].to_numpy() if "dustAmp" in df else None
            arr = paint_dust_test_sources(pos.ra_rad, pos.dec_rad, geometry, amplitudes_uk=amp, show_progress=show)
            label = r"$\Delta T$ [$\mu$K]"
        else:
            arr = paint_circular_masks(
                pos.ra_rad, pos.dec_rad, df["Thetavir"].to_numpy(), geometry,
                combine="union" if args.union else "sum", show_progress=show,
            )
            label = "mask"

    _save(arr, args.output, args.png, geometry, label)


if __name__ == "__main__":
    main()
