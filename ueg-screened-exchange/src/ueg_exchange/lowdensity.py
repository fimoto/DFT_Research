"""Low-density (strong-screening) coefficients A1, A2, A3 from their profile integrals.

Companion paper, Sec. IV and App. A.  Reduced units; energies in Ry per electron;
e_kite = A1/q* + A2/q*^2 + A3/q*^3 + ..., q*^4 = (4/3) q_tf_sq, q_tf_sq = a rs.

The A3 pieces are computed here from the explicit profiles, not copied:
  cross  : 8 pi A0 V_B^2 * Fp int_0^inf x^-4 J_C(x^-4) dx            (closed form via Beta functions)
  recoil : 4 pi A0 (104 pi^2/27) int_0^inf x^-4 Phi_4(x^-4) dx
  Lindhard: 24 pi^3 A0 int_0^inf x^-8 Psi(x^-4) dx
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma, beta
from .constants import A0

C_GAMMA = np.sqrt(2 * np.pi) / np.pi ** 2 * gamma(0.75) / gamma(0.25)
V_B = 4 * np.pi / 3


def R(c):
    c = float(c)
    return 2 * (c - np.log1p(c)) / c ** 2 if c > 1e-6 else 1 - 2 * c / 3 + c * c / 2 - 2 * c ** 3 / 5


def dR(c):
    c = float(c)
    if c > 1e-4:
        return (2 * c / (1 + c) - 4 * (c - np.log1p(c)) / c) / c ** 2
    return -2.0 / 3 + c - 1.2 * c * c


def J_C(g):
    """Cross-ball frequency integral J_C(g) = 2 int lam dlam / sqrt(1+lam g)."""
    g = float(g)
    if g < 1e-2:   # below this the closed form loses digits to cancellation
        return 1 - g / 3 + 3 * g * g / 16 - g ** 3 / 8 + 35 * g ** 4 / 384
    return 4 / (3 * g * g) * ((1 + g) ** 1.5 - 3 * np.sqrt(1 + g) + 2)


def lindhard_coeff_4(u):
    """q^-4 coefficient of the large-q Lindhard function at nu = q^2 u."""
    return -16 * (12 * u * u - 1) / (15 * (1 + 4 * u * u) ** 3)


def _reA(u, n):
    return ((1j * u - 0.5) ** (-n)).real


def A3_cross_closed():
    """Fp int_0^inf x^-4 J_C dx = (2/3)[2^(1/2) B(3/4,-1/2) - 2^(-3/2) B(3/4,3/2)] = -(8/5) sqrt(2pi) G(3/4)/G(1/4)."""
    Fp = 2.0 / 3 * (2 ** 0.5 * beta(0.75, -0.5) - 2 ** -1.5 * beta(0.75, 1.5))
    return 8 * np.pi * A0 * V_B ** 2 * Fp


def A3_cross_numeric():
    f = lambda x: J_C(x ** -4) / x ** 4
    Fp = quad(lambda x: f(x) - 4 / (3 * x * x), 0, 1, limit=200)[0] - 4.0 / 3 + quad(f, 1, np.inf, limit=200)[0]
    return 8 * np.pi * A0 * V_B ** 2 * Fp


def Phi4(g):
    return quad(lambda u: R(g / (1 + 4 * u * u)) * _reA(u, 4), 0, np.inf, limit=400)[0] / np.pi


def A3_recoil_numeric():
    I = quad(lambda x: Phi4(x ** -4) / x ** 4, 0, np.inf, limit=200)[0]
    return 4 * np.pi * A0 * (104 * np.pi ** 2 / 27) * I


def Psi(g):
    return quad(lambda u: dR(g / (1 + 4 * u * u)) * lindhard_coeff_4(u) * _reA(u, 2), 0, np.inf, limit=400)[0] / np.pi


def A3_lindhard_numeric():
    I = quad(lambda x: Psi(x ** -4) / x ** 8, 0, np.inf, limit=200)[0]
    return 24 * np.pi ** 3 * A0 * I


# --- closed-form reductions (companion App. A) -------------------------------

def I_recoil(m):
    """int (du/pi) (1+4u^2) Re[(iu-1/2)^-4]/(1+4u^2+m) = -4m/[sqrt(1+m)(sqrt(1+m)+1)^4]."""
    m = np.asarray(m, float)
    s = np.sqrt(1 + m)
    return -4 * m / (s * (s + 1) ** 4)


def U_lindhard(m):
    """-(16/15) int (du/pi) (12u^2-1) Re[(iu-1/2)^-2]/[(1+4u^2)(1+4u^2+m)^2]
       = 8(2s^3-s^2+4s+1)/[15 s^3 (s+1)^4], s = sqrt(1+m)."""
    m = np.asarray(m, float)
    s = np.sqrt(1 + m)
    return 8 * (2 * s ** 3 - s * s + 4 * s + 1) / (15 * s ** 3 * (s + 1) ** 4)


def A3_recoil_closed():
    """4 pi A0 (104 pi^2/27) * (2/5) * (-2 sqrt2 B(7/4,3/2)) = -(416/675) C_Gamma."""
    K = -2 * np.sqrt(2) * beta(1.75, 1.5)
    return 4 * np.pi * A0 * (104 * np.pi ** 2 / 27) * (2.0 / 5) * K


def M_euler():
    """M = (16/15) int_0^1 w^-1/2 (1-w)^3/4 (1+w)^-13/4 (2-w+4w^2+w^3) dw  (companion Eq. MEuler)."""
    f = lambda w: w ** -0.5 * (1 - w) ** 0.75 * (1 + w) ** -3.25 * (2 - w + 4 * w * w + w ** 3)
    return 16.0 / 15 * quad(f, 0, 1, limit=400)[0]


def A3_lindhard_closed():
    """-(48/5) pi^3 A0 M, with M = (176/75) sqrt(2 pi) Gamma(3/4)/Gamma(1/4) = -(264/125) C_Gamma."""
    return -(48.0 / 5) * np.pi ** 3 * A0 * M_euler()


def A3_pieces_claimed():
    return {"cross": -32 / 15 * C_GAMMA, "recoil": -416 / 675 * C_GAMMA, "lindhard": -264 / 125 * C_GAMMA}


def A3_from_profiles(closed=True):
    """A3 from the three profiles.  closed=True uses the reduced 1-D forms of App. A,
    closed=False the raw multi-dimensional profile quadratures."""
    if closed:
        return A3_cross_closed() + A3_recoil_closed() + A3_lindhard_closed()
    return A3_cross_numeric() + A3_recoil_numeric() + A3_lindhard_numeric()
