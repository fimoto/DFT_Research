"""Static (z = 0) finite-q first-order exchange polarization of the 2D electron gas,
compared with Table I of S. Das Sarma, arXiv:2610.08309.

This is the crescent-coordinate code that produced the row Delta R/R(0) printed in
MPP SI Sec. S6 ("At z=0 we compare ..."): nphi = 400, nt = 10 reproduces
4.454e-4, 2.020e-3, 5.096e-3, 9.625e-3 (the SI prints these truncated to four digits).
It was written for z > 0 and is not adapted to the coincident-momentum singularity of the
static integrand, so it is NOT converged at z = 0 (see the nphi = 600 column).

Units: Hartree, k_F = 1, paramagnetic.
  chi1 = int d^2p d^2k/(2 pi)^4 v(p-k) Delta_p Delta_k Re[(g_p - g_k)^2],  v = 2 pi/|p-k|,
  g_p = 1/(i w - E_p),  E_p = p.q + q^2/2,  w = q z,
  R(q, z) = -chi1/N0,  N0 = 1/pi,  R_2(0) = 2/pi.
Crescent coordinates: rho = rho*(phi) + (1 - rho*(phi)) t; the signed measure covers Delta_p = +-1.
The phi integration is split at phi = +-pi/2 (Gauss-Legendre on each half).

usage: python3 finiteq2d_static.py            (about 10-20 min on one core)
"""
import numpy as np

# Das Sarma, arXiv:2610.08309, Table I, Delta P/P(0) (as quoted in MPP SI Sec. S6)
DAS_SARMA = {0.1: 4.524e-4, 0.2: 2.104e-3, 0.3: 5.141e-3, 0.4: 9.699e-3}
BETA = 6.6441  # small-x form x^2 (ln x + beta)/96, SI Sec. S6


def gl(n, a, b):
    x, w = np.polynomial.legendre.leggauss(n)
    return (b - a) / 2 * x + (a + b) / 2, (b - a) / 2 * w


def nodes(q, nphi, nt):
    ph, wp = [], []
    for a, b in [(-np.pi / 2, np.pi / 2), (np.pi / 2, 3 * np.pi / 2)]:
        x, w = gl(nphi, a, b)
        ph.append(x); wp.append(w)
    ph = np.concatenate(ph); wp = np.concatenate(wp)
    t, wt = gl(nt, 0, 1)
    PH, T = np.meshgrid(ph, t, indexing='ij')
    rs = -q * np.cos(PH) + np.sqrt(1 - q * q * np.sin(PH) ** 2)
    rho = rs + (1 - rs) * T
    W = (wp[:, None] * wt[None, :]) * (1 - rs) * rho
    P = np.stack([rho * np.cos(PH), rho * np.sin(PH)], -1).reshape(-1, 2)
    return P, W.ravel()


def R(q, z, nphi=400, nt=10, chunk=3000):
    P, W = nodes(q, nphi, nt)
    w = q * z
    g = 1 / (1j * w - (P[:, 0] * q + q * q / 2))
    tot = 0.0
    for i in range(0, len(W), chunk):
        dp = P[i:i + chunk, None, :] - P[None, :, :]
        r = np.sqrt(np.sum(dp * dp, -1)); r[r == 0] = np.inf   # drop the diagonal
        f = (2 * np.pi / r) * np.real((g[i:i + chunk, None] - g[None, :]) ** 2)
        tot += np.sum(W[i:i + chunk, None] * W[None, :] * f)
    return -tot / (2 * np.pi) ** 4 * np.pi


if __name__ == '__main__':
    R0 = 2 / np.pi
    runs = [(300, 6), (400, 6), (400, 10), (600, 10)]
    print('x     ' + '  '.join(f'nphi={a},nt={b}' for a, b in runs)
          + '   asympt x^2(ln x+b)/96   DasSarma TableI', flush=True)
    for x in (0.1, 0.2, 0.3, 0.4):
        vals = [(R(x, 0.0, nphi=a, nt=b) - R0) / R0 for a, b in runs]
        asym = x * x * (np.log(x) + BETA) / 96
        print(f'{x:.1f}  ' + '  '.join(f'{v:.5e}'.rjust(15) for v in vals)
              + f'   {asym:.4e}               {DAS_SARMA[x]:.3e}', flush=True)
