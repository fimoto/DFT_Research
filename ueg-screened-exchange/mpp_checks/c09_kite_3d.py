"""SI Sec. S19 ("Frequency integrals and normalization"): the checks of the
three-dimensional kite that are not already in kernels_checks/.

Reproduces: the closed form of the Lindhard function against direct
three-dimensional quadrature of its representation (1/4pi) Int d^3p Delta_p
delta_p/(delta_p^2 + nu^2) at (q, nu) = (0.7, 0.4) and (3, 2) with 30-digit working
precision (agreement to 1e-25); the symbolic kite prefactor A_0 = 3/(32 pi^5),
SI Eq. (S18); and a Monte Carlo evaluation of K[v] = e_2x, main Eq. (43), with the
sampling scheme described in SI Sec. S19 (p uniform in C_q, k from an equal mixture
of the uniform density on C_q and a density proportional to |l|^-2).

Run with --mc-samples N --mc-runs M to set the Monte Carlo size; the SI quotes
four runs of 1e8 samples.  The default here is four runs of 2.5e7 samples.
"""
import sys, argparse
import numpy as np, mpmath as mp, sympy as sp
from common import Recorder

ap = argparse.ArgumentParser(); ap.add_argument('--mc-samples', type=float, default=2.5e7); ap.add_argument('--mc-runs', type=int, default=4)
ap.add_argument('--seed', type=int, default=20261019)
args, _ = ap.parse_known_args()
R = Recorder("c09_kite_3d", "3D kite: Lindhard closed form, prefactor A_0 and Monte Carlo bare limit (SI Sec. S19)")

# ---------------------------------------------------------------- Lindhard function, 30 digits
mp.mp.dps = 30
def _F(p, q, b, n):
    """antiderivative of p log((q p + b)^2 + n^2) in p, written with x = q p + b"""
    x = q * p + b; Lg = mp.log(x * x + n * n)
    return (((x * x + n * n) * Lg - x * x) / 2 - b * (x * Lg - 2 * x + 2 * n * mp.atan(x / n))) / (q * q)
def L_closed(q, nu):
    """(1/(2q)) Int_0^1 p log[((pq+q^2/2)^2+nu^2)/((pq-q^2/2)^2+nu^2)] dp with the elementary antiderivative _F."""
    q, nu = mp.mpf(q), mp.mpf(nu)
    I1 = _F(1, q, q * q / 2, nu) - _F(0, q, q * q / 2, nu)
    I2 = _F(1, q, -q * q / 2, nu) - _F(0, q, -q * q / 2, nu)
    return (I1 - I2) / (2 * q)
def L_direct(q, nu):
    """(1/4pi) Int d^3p Delta_p delta_p/(delta_p^2+nu^2) over the two transition regions, in (rho, mu)."""
    q, nu = mp.mpf(q), mp.mpf(nu)
    f = lambda rho, mu: (rho * q * mu + q * q / 2) / ((rho * q * mu + q * q / 2) ** 2 + nu * nu) * 2 * mp.pi * rho * rho
    mus = lambda rho: (1 - rho * rho - q * q) / (2 * rho * q)      # |p+q| = 1
    def cplus(rho):   # |p| < 1 < |p+q|: mu > mu*
        lo = max(mp.mpf(-1), min(mp.mpf(1), mus(rho)))
        return mp.quad(lambda mu: f(rho, mu), [lo, 1])
    def cminus(rho):  # |p+q| < 1 < |p|: mu < mu*, rho in (1, 1+q)
        hi = max(mp.mpf(-1), min(mp.mpf(1), mus(rho)))
        return mp.quad(lambda mu: f(rho, mu), [-1, hi])
    rb = sorted(set([0, max(0, 1 - q), 1]))
    I = mp.quad(cplus, rb) - mp.quad(cminus, sorted(set([1, min(1 + q, max(1, q - 1)) if q > 2 else 1 + q, 1 + q])) if q <= 2 else [q - 1, q + 1])
    return I / (4 * mp.pi)
for q, nu in [('0.7', '0.4'), ('3', '2')]:
    R.close("SI Sec. S19", f"Lindhard closed form vs direct 3D quadrature, (q, nu) = ({q}, {nu}), 30 digits", L_direct(q, nu), L_closed(q, nu), mp.mpf('1e-25'))

# ---------------------------------------------------------------- A_0, SI Eq. (S18)
e, hb, m, kF = sp.symbols('e hbar m k_F', positive=True)
ne = kF ** 3 / (3 * sp.pi ** 2); v0 = 4 * sp.pi * e ** 2 / kF ** 2; E0 = hb ** 2 * kF ** 2 / m
A0 = sp.simplify(2 / (4 * ne) * kF ** 9 / (2 * sp.pi) ** 9 * v0 ** 2 / E0 / (m * e ** 4 / (2 * hb ** 2)))
R.exact("SI Eq. (S18), SI Sec. S19", "A_0 [Ry] = 3/(32 pi^5)", sp.simplify(A0 - sp.Rational(3, 32) / sp.pi ** 5), 0)

# ---------------------------------------------------------------- Monte Carlo of K[v]
rng = np.random.default_rng(args.seed)
A0n = 3 / (32 * np.pi ** 5)
def vol_C(q):
    return np.where(q >= 2, 4 * np.pi / 3, 4 * np.pi / 3 - np.pi * (4 + q) * (2 - np.minimum(q, 2)) ** 2 / 12)
def sample_C(q):
    """uniform points in C_q = {|p| < 1 < |p+q z^|}: propose in the shell max(0,1-q) < |p| < 1, accept |p+q| > 1."""
    n = q.size; out = np.empty((n, 3)); todo = np.arange(n)
    while todo.size:
        qq = q[todo]; r0 = np.maximum(0.0, 1 - qq)
        u = rng.random(todo.size); rad = (r0 ** 3 + u * (1 - r0 ** 3)) ** (1 / 3)
        ct = 2 * rng.random(todo.size) - 1; ph = 2 * np.pi * rng.random(todo.size); st = np.sqrt(1 - ct * ct)
        p = np.stack([rad * st * np.cos(ph), rad * st * np.sin(ph), rad * ct], axis=1)
        ok = (p[:, 0] ** 2 + p[:, 1] ** 2 + (p[:, 2] + qq) ** 2) > 1
        out[todo[ok]] = p[ok]; todo = todo[~ok]
    return out
def mc_run(n, batch=2_000_000):
    tot = 0.0; tot2 = 0.0; cnt = 0
    while cnt < n:
        b = int(min(batch, n - cnt))
        uq = rng.random(b); q = uq / (1 - uq); gq = (1 - uq) ** 2          # density g(q) = 1/(1+q)^2
        Vq = vol_C(q)
        p = sample_C(q)
        # k: mixture of uniform-in-C_q and |l|^-2 on the ball |l| < 2 + q
        Lmax = 2 + q
        pick = rng.random(b) < 0.5
        k = sample_C(q)
        lr = Lmax * rng.random(b); ct = 2 * rng.random(b) - 1; ph = 2 * np.pi * rng.random(b); st = np.sqrt(1 - ct * ct)
        lvec = np.stack([lr * st * np.cos(ph), lr * st * np.sin(ph), lr * ct], axis=1)
        k2 = lvec - p - np.stack([np.zeros(b), np.zeros(b), q], axis=1)
        k = np.where(pick[:, None], k, k2)
        l = p + k; l[:, 2] += q
        l2 = np.sum(l * l, axis=1)
        inC = (np.sum(k * k, axis=1) < 1) & ((k[:, 0] ** 2 + k[:, 1] ** 2 + (k[:, 2] + q) ** 2) > 1)
        mix = 0.5 / Vq + 0.5 / (4 * np.pi * Lmax * l2) * (np.sqrt(l2) < Lmax)
        w = np.where(inC, 2 * A0n * 4 * np.pi / (l2 * q * l[:, 2]) * Vq / mix / gq, 0.0)
        tot += w.sum(); tot2 += (w * w).sum(); cnt += b
    return tot / cnt
runs = [mc_run(int(args.mc_samples)) for _ in range(args.mc_runs)]
mean = float(np.mean(runs)); se = float(np.std(runs, ddof=1) / np.sqrt(len(runs)))
e2x = float(mp.log(2) / 3 - 3 * mp.zeta(3) / (2 * mp.pi ** 2))
R.note("SI Sec. S19", f"Monte Carlo K[v], {args.mc_runs} runs of {int(args.mc_samples)} samples (seed {args.seed})", f"{mean:.6f} +- {se:.6f} (runs: {', '.join(f'{x:.6f}' for x in runs)})")
R.close("SI Sec. S19, main Sec. 3.2", "Monte Carlo K[v] agrees with e_2x within 3 standard errors", mean, e2x, 3 * se, printed="0.04836(2) vs e_2x = 0.0483583")
R.digits("main Eq. (43)", "e_2x = log(2)/3 - 3 zeta(3)/(2 pi^2)", mp.log(2) / 3 - 3 * mp.zeta(3) / (2 * mp.pi ** 2), "0.0483583")
raise SystemExit(R.write(dict(mc=dict(runs=runs, mean=mean, se=se, samples_per_run=int(args.mc_samples), seed=args.seed))))
