"""Spectral dispersion of the screening pole (JPA screening-pole analysis) and the leading low-density coefficient.

At large q the RPA line is lambda times a bare line minus a weighted Lorentzian kernel, with
    theta_lam(q) = lam q_*^4 / (q^4 + lam q_*^4),
    Lambda_lam(q) = (q/2) sqrt(1 + lam (q_*/q)^4).
Integrating radial_energy_large_q along that dispersion curve reproduces A1 in closed form.
"""
import numpy as np
from scipy.integrate import quad
from .constants import A1, q_star
from .lorentzian import p_loc


def theta(x, lam):
    """theta_lam at x = q/q_*."""
    return lam / (x ** 4 + lam)


def scaled_pole_frequency(x, lam):
    """scaled_pole_frequency = Lambda_lam(q)/q along the pole dispersion, as a function of x = q/q_*."""
    return 0.5 * np.sqrt(1.0 + lam * x ** -4)


def A1_dispersion(epsabs=1e-12):
    """A1 = -(12/pi^2) int_0^1 lam dlam int_0^inf dx x^-2 theta(x,lam) s_lambda/(1+2*s_lambda)^2."""
    def inner(lam):
        f = lambda x: x ** -2 * theta(x, lam) * scaled_pole_frequency(x, lam) / (1 + 2 * scaled_pole_frequency(x, lam)) ** 2
        s = 0.0
        for a, b in [(0, 0.3), (0.3, 1), (1, 3), (3, np.inf)]:
            s += quad(f, a, b, limit=200, epsabs=epsabs)[0]
        return s
    val = quad(lambda l: l * inner(l), 0.0, 1.0, limit=200, epsabs=epsabs)[0]
    return -12.0 / np.pi ** 2 * val


def Lambda_pole(q, rs, lam=1.0):
    """Width Lambda_lambda(q) in the scaled frequency nu/q; q is in units of k_F."""
    qs = q_star(rs)
    return 0.5 * q * np.sqrt(1.0 + lam * (qs / q) ** 4)
