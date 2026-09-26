"""Catalog containers and adapters for flat-patch catalog columns."""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np
import pandas as pd


def _arr(x):
    return np.atleast_1d(np.asarray(x, dtype=float))


@dataclass
class PositionCatalog:
    ra_rad: np.ndarray
    dec_rad: np.ndarray

    def __post_init__(self):
        self.ra_rad = _arr(self.ra_rad)
        self.dec_rad = _arr(self.dec_rad)
        if len(self.ra_rad) != len(self.dec_rad):
            raise ValueError("ra_rad and dec_rad must have equal length")

    def __len__(self):
        return len(self.ra_rad)

    @classmethod
    def from_dataframe(cls, df: pd.DataFrame):
        # Accept both standard and shifted flat-patch coordinate names.
        ra_name = "RArad" if "RArad" in df else "shiftedRArad"
        dec_name = "DECrad" if "DECrad" in df else "shiftedDECrad"
        return cls(df[ra_name].to_numpy(), df[dec_name].to_numpy())


@dataclass
class MovingLensCatalog(PositionCatalog):
    redshift: np.ndarray
    concentration: np.ndarray
    rs_mpc: np.ndarray
    rho_s_msun_mpc3: np.ndarray
    comoving_distance_mpc: np.ndarray
    v_theta_kms: np.ndarray
    v_phi_kms: np.ndarray

    def __post_init__(self):
        super().__post_init__()
        names = [
            "redshift", "concentration", "rs_mpc", "rho_s_msun_mpc3",
            "comoving_distance_mpc", "v_theta_kms", "v_phi_kms",
        ]
        for name in names:
            setattr(self, name, _arr(getattr(self, name)))
        lengths = {len(self.ra_rad)} | {len(getattr(self, n)) for n in names}
        if len(lengths) != 1:
            raise ValueError("All moving-lens catalog columns must have equal length")
        if np.any(self.redshift < 0):
            raise ValueError("redshift must be non-negative")
        if np.any(self.concentration <= 0):
            raise ValueError("concentration must be positive")
        if np.any(self.rs_mpc <= 0) or np.any(self.rho_s_msun_mpc3 <= 0):
            raise ValueError("NFW scale radius and density must be positive")
        if np.any(self.comoving_distance_mpc <= 0):
            raise ValueError("comoving distance must be positive")

    @property
    def scale_factor(self):
        return 1.0 / (1.0 + self.redshift)

    @property
    def transverse_speed_kms(self):
        return np.hypot(self.v_theta_kms, self.v_phi_kms)

    @property
    def transverse_angle_rad(self):
        return np.arctan2(self.v_theta_kms, self.v_phi_kms)

    def row(self, i):
        return {
            "ra_rad": float(self.ra_rad[i]),
            "dec_rad": float(self.dec_rad[i]),
            "redshift": float(self.redshift[i]),
            "concentration": float(self.concentration[i]),
            "rs_mpc": float(self.rs_mpc[i]),
            "rho_s_msun_mpc3": float(self.rho_s_msun_mpc3[i]),
            "comoving_distance_mpc": float(self.comoving_distance_mpc[i]),
            "v_theta_kms": float(self.v_theta_kms[i]),
            "v_phi_kms": float(self.v_phi_kms[i]),
        }

    @classmethod
    def from_dataframe(cls, df):
        pos = PositionCatalog.from_dataframe(df)
        return cls(
            ra_rad=pos.ra_rad,
            dec_rad=pos.dec_rad,
            redshift=df["Z"].to_numpy(),
            concentration=df["cNFW"].to_numpy(),
            rs_mpc=df["Rs"].to_numpy(),
            rho_s_msun_mpc3=df["rhoS"].to_numpy(),
            comoving_distance_mpc=df["comovDist"].to_numpy(),
            v_theta_kms=df["vTh"].to_numpy(),
            v_phi_kms=df["vPh"].to_numpy(),
        )


def read_csv(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path)
