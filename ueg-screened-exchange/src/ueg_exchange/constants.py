"""Units, Feynman-rule constant, derived scales.

Reduced units: momentum k_F, internal energy E0=hbar^2*k_F^2/m=2*E_F,
angular frequency omega0=E0/hbar, interaction v0=4*pi*e^2/k_F^2, and
density response 1/v0. Thus eps_k=k^2/2 and v(q)=1/q^2. Reported final
energies are Rydberg per electron; the internal and final units differ.
"""
import numpy as np
from math import gamma as _gamma

ALPHA_L = (4.0 / (9.0 * np.pi)) ** (1.0 / 3.0)
Q_TF_SQ_COEFF = 4.0 * ALPHA_L / np.pi          # q_TF^2 = Q_TF_SQ_COEFF * rs
A0 = 3.0 / (32.0 * np.pi ** 5)         # Feynman-rule constant, Eq. (A0)
ZETA3 = 1.2020569031595942854
E2X = np.log(2.0) / 3.0 - 3.0 * ZETA3 / (2.0 * np.pi ** 2)  # Onsager, Ry
VBALL = 4.0 * np.pi / 3.0

C_GAMMA = np.sqrt(2 * np.pi) / np.pi ** 2 * _gamma(0.75) / _gamma(0.25)
A1 = -4.0 / (7.0 * np.sqrt(2.0) * np.pi ** 1.5) * _gamma(0.25) / _gamma(0.75)
A2 = 8.0 / 15.0
A3 = -16408.0 / 3375.0 * C_GAMMA
C_LN2 = -ALPHA_L * (1.0 - np.log(2.0)) / (6.0 * np.pi ** 3)   # rs ln^2 rs, Ry
C_L_RING = 2.0 * (1.0 - np.log(2.0)) / np.pi ** 2             # GMB ring log, Ry


def q_tf_sq(rs):
    return Q_TF_SQ_COEFF * np.asarray(rs, float)


def q_star(rs):
    return (4.0 * q_tf_sq(rs) / 3.0) ** 0.25


def spin_radii(zeta):
    if not -1.0 <= zeta <= 1.0:
        raise ValueError("zeta must lie in [-1,1]")
    return (1.0 + zeta) ** (1.0 / 3.0), (1.0 - zeta) ** (1.0 / 3.0)


def spin_moment(n, zeta):
    xu, xd = spin_radii(zeta)
    return 0.5 * (xu ** n + xd ** n)


def asymptotic_kite(rs, zeta=0.0):
    """Three-term low-density expansion, Ry per electron (companion paper)."""
    qs = q_star(rs)
    M4, M5, M6 = (spin_moment(n, zeta) for n in (4, 5, 6))
    S3 = (1160 * M6 + 891 * M4 * M5) / 2051.0
    return A1 * M4 / qs + A2 * M5 / qs ** 2 + A3 * S3 / qs ** 3
