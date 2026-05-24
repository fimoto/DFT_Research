"""
Monte Carlo evaluation of the momentum-resolved kite kernel X_c(q)
(static frequency kernel, units k_F = 1):

  X_c(q) = (1/16pi) * Int d3p d3k  Delta(p;q) Delta(k;q) W(dp,dk) / |p-k|^2,
  Delta(p;q) = Theta(1-|p|) - Theta(1-|p+q zhat|),
  dp = q*p_z + q^2/2,   W(a,b) = -1{ab<0} / (|a|+|b|),

normalized so that  Int_0^inf X_c(q) dq = D0 = pi^3/18*log2 - pi/4*zeta(3)
for the bare interaction (checked by run_scan.py) and X_c(q->0) = c0*q with
c0 = pi(2log2-1)/8 (checked by validate() below).

Sampling design (see the PRB Supplemental Material):
 * union-of-balls importance sampling for p AND k: uniform over
   ball(0,1) u ball(-q zhat,1); acceptance ~ 1 for every q, no variance
   blow-up at q ~ 2 where the balls separate;
 * additional Coulomb mixture component for k (k = p + r*omega, r ~ U(0,Rc))
   with weight beta, which cancels the leading Coulomb singularity in the
   importance weight; finiteness of the complete second moment is not proved;
 * batch means with a median-of-means cross-check (the residual
   1/(|dp|+|dk|) near-zero-plane tails are log-heavy).

run_to_precision(q, abs_err=..., rel_err=...) grows the sample size until
a target on the replicate scatter (not a confidence interval) is met.
"""
import sys
from functools import lru_cache
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"analytic"))

VBALL = 4*np.pi/3
ZETA3 = 1.2020569031595942854
D0_EXACT = np.pi**3/18*np.log(2) - np.pi/4*ZETA3       # 0.2499018971
C0_EXACT = np.pi*(2*np.log(2)-1)/8                     # 0.151697440877...


def _ball_points(rng, n):
    u = rng.random((n, 3))
    r = u[:, 0]**(1/3.)
    ct = 2*u[:, 1] - 1
    ph = 2*np.pi*u[:, 2]
    st = np.sqrt(1 - ct**2)
    return np.column_stack([r*st*np.cos(ph), r*st*np.sin(ph), r*ct])


def _f_union(x, q):
    """Density of the uniform-union sampler evaluated at points x."""
    inA = np.einsum('ij,ij->i', x, x) < 1.0
    xB = x.copy(); xB[:, 2] += q
    inB = np.einsum('ij,ij->i', xB, xB) < 1.0
    return (inA.astype(float) + inB.astype(float))/(2*VBALL)


def _sample_union(rng, q, n):
    pts = _ball_points(rng, n)
    which = rng.integers(0, 2, n)
    pts[:, 2] -= q*which
    return pts, _f_union(pts, q)


def _delta_occ(x, q):
    inA = np.einsum('ij,ij->i', x, x) < 1.0
    xB = x.copy(); xB[:, 2] += q
    inB = np.einsum('ij,ij->i', xB, xB) < 1.0
    return inA.astype(float) - inB.astype(float)


def make_dyn_kernel(q, rs, nxi=24, nseg=10, zeta=0.0, xspin=1.0,
                    corner_refined=True):
    """Precompute (xi_j, cw_j) for the dynamical (RPA-screened,
    coupling-constant-integrated) frequency kernel

      K(a,b) = W(a,b) + sum_j cw_j (a b - xi_j^2)/((a^2+xi_j^2)(b^2+xi_j^2)),
      cw_j = w_j [R(q,xi_j;rs) - 1]/pi,

    which equals (1/pi) Int_0^inf dxi u(q,i xi)(ab-xi^2)/(...) with the
    imaginary-axis screened-to-bare line ratio u = R.  This is exactly the
    paper's spectral (Lehmann) superposition of single-pole kernels
    J(omega;a,b), rewritten on the imaginary axis where it is smooth
    (R = 1, i.e. rs -> 0, gives back the static bare kernel W)."""
    from lindhard import R_ac
    xg, wg = np.polynomial.legendre.leggauss(nxi)
    s0 = max(q*(1 + q/2), 0.3)
    if corner_refined:
        # Extend log panels down to 1e-7*s0 so the
        # narrow xi ~ |a|,|b| structure of near-Fermi-surface samples is
        # resolved (the standard grid's smallest node is ~7e-5*s0)
        edges = np.concatenate([[0.0], s0*np.geomspace(1e-7, 1e-3, 9)[:-1],
                                s0*np.geomspace(1e-3, 1.0, nseg)])
    else:
        edges = np.concatenate([[0.0], s0*np.geomspace(1e-3, 1.0, nseg)])
    xi = (0.5*(edges[1:, None] - edges[:-1, None])*xg[None, :]
          + 0.5*(edges[1:, None] + edges[:-1, None])).ravel()
    w = (0.5*(edges[1:, None] - edges[:-1, None])*wg[None, :]).ravel()
    t = 0.5*(xg + 1)*(1 - 1e-9)
    xi_t = s0*(1 + t/(1 - t))
    w_t = 0.5*(1 - 1e-9)*wg*s0/(1 - t)**2
    xi = np.concatenate([xi, xi_t]); w = np.concatenate([w, w_t])
    # spin-sigma loop in scaled variables (Fermi radius xspin):
    # physical (q, xi) = (xspin*q~, xspin^2*xi~); zeta enters the
    # two-component dielectric.
    Rm1 = R_ac(xspin*q, xspin**2*xi, rs, zeta=zeta) - 1.0
    return xi, w*Rm1/np.pi


def _dyn_correction(dp, dk, kern):
    """sum_j cw_j (dp dk - xi_j^2)/((dp^2+xi_j^2)(dk^2+xi_j^2)),
    chunked over samples."""
    xi, cw = kern
    out = np.empty_like(dp)
    step = max(1, 2_000_000//max(len(xi), 1))
    for i in range(0, len(dp), step):
        a = dp[i:i + step, None]; b = dk[i:i + step, None]
        out[i:i + step] = ((a*b - xi[None, :]**2)
                           / ((a*a + xi[None, :]**2)
                              * (b*b + xi[None, :]**2)) @ cw)
    return out


@lru_cache(maxsize=256)
def _pole1_params_cached(q, rs, zeta, xspin):
    """Cached implementation with every physical input in the cache key."""
    if abs(zeta) > 1e-12 or abs(xspin - 1.0) > 1e-12:
        raise NotImplementedError(
            "kernel='pole1' is currently derived only for zeta=0 and "
            "xspin=1; polarized use requires polarized spectral moments")
    from lindhard import c_screen
    c0 = float(np.asarray(c_screen(q, np.atleast_1d(1e-9), rs)).ravel()[0])
    R0 = 2.0*(c0 - np.log1p(c0))/c0**2 if c0 > 1e-8 else 1.0 - 2*c0/3
    eta = 1.0 - R0
    xL = 60.0*max(1.0, q*q)
    m2a = float(np.asarray(c_screen(q, np.atleast_1d(xL), rs)).ravel()[0])*xL*xL
    m2b = float(np.asarray(c_screen(q, np.atleast_1d(2*xL), rs)).ravel()[0])*4*xL*xL
    m2 = (4*m2b - m2a)/3.0            # Richardson in 1/xi^2
    M2 = (2.0/3.0)*m2
    return eta, M2/max(eta, 1e-300)


def _pole1_params(q, rs, zeta=0.0, xspin=1.0):
    """Moment-conserving one-pole embedding of the coupling-averaged RPA
    line.  The rounded immutable key includes all physical inputs, and
    unsupported polarized calls fail explicitly rather than reusing an
    unpolarized cache entry silently."""
    return _pole1_params_cached(round(float(q), 10), round(float(rs), 10),
                                round(float(zeta), 10),
                                round(float(xspin), 10))


def _J_pole(w, a, b):
    """Exact closed form of J(w;a,b) = (1/pi) Int_0^inf (ab-xi^2) dxi /
    ((a^2+xi^2)(b^2+xi^2)(xi^2+w^2)), by sign branch:
      ab > 0:  1/(2w(al+w)(be+w))
      ab < 0: -(al+be+2w)/(2w(al+be)(al+w)(be+w))
    with al=|a|, be=|b|. All coincidences (al=be, al=w, ...) are removable;
    no regularization is required.
    """
    al, be = np.abs(a), np.abs(b)
    same = (a*b) > 0
    out = np.empty_like(al)
    out[same] = 1.0/(2*w*(al[same] + w)*(be[same] + w))
    o = ~same
    s_ab = al[o] + be[o]
    out[o] = -(s_ab + 2*w)/(2*w*np.maximum(s_ab, 1e-300)
                            * (al[o] + w)*(be[o] + w))
    return out


_pole_cache = {}


def xc_samples(q, n, rng, beta=0.5, Rc=None, kernel="static", rs=None,
               zeta=0.0, xspin=1.0, lam=None, _kern_cache={}, pole_N=1):
    """Return n one-sample estimates of the chosen momentum-resolved kernel.

    kernel = 'static' : bare static frequency kernel W (default);
             'dyn'    : full dynamical AC kernel K (requires rs);
             'diff'   : K - W only (the screening-induced correction;
                        use with D0 exact for maximal variance reduction).
             'lorentzian', 'lorentzian_same', 'lorentzian_opp': finite
                        positive scaled-frequency width lam; the last two
                        select equal/opposite signs of the energy denominators.
    Additional approximate screening rules are 'pole1', 'poleN' and 'poleNerr'.
    The static choice has the X_c normalization stated at the top of this file.
    """
    allowed_kernels = ('static', 'dyn', 'diff', 'lorentzian',
                       'lorentzian_same', 'lorentzian_opp', 'pole1', 'poleN', 'poleNerr')
    if kernel not in allowed_kernels:
        raise ValueError(f'Unknown kernel {kernel!r}; choose one of {allowed_kernels}')
    if kernel.startswith('lorentzian') and (lam is None or not np.isfinite(lam) or lam <= 0):
        raise ValueError('A Lorentzian kernel requires a finite positive width lam')
    if Rc is None:
        Rc = 2.0 + q
    p, fp = _sample_union(rng, q, n)

    # k: mixture of union sampler and Coulomb sampler around p
    use_c = rng.random(n) < beta
    ku, _ = _sample_union(rng, q, n)
    r = rng.random(n)*Rc
    ct = 2*rng.random(n) - 1
    ph = 2*np.pi*rng.random(n)
    st = np.sqrt(1 - ct**2)
    kc = p + np.column_stack([r*st*np.cos(ph), r*st*np.sin(ph), r*ct])
    k = np.where(use_c[:, None], kc, ku)

    r2 = np.einsum('ij,ij->i', k - p, k - p)
    fc = np.where(r2 < Rc*Rc, 1.0/(4*np.pi*np.maximum(r2, 1e-300)*Rc), 0.0)
    fk = (1 - beta)*_f_union(k, q) + beta*fc

    Dp = _delta_occ(p, q)
    Dk = _delta_occ(k, q)
    dp = q*p[:, 2] + q*q/2
    dk = q*k[:, 2] + q*q/2
    occ = (fp > 0) & (fk > 0) & (Dp != 0) & (Dk != 0)
    opp = ((dp*dk) < 0) & occ

    out = np.zeros(n)
    ker = np.zeros(n)
    if kernel in ("static", "dyn"):
        ker[opp] = -1.0/(np.abs(dp[opp]) + np.abs(dk[opp]))
    if kernel == "pole1":
        # one-pole correction kernel: -eta_q Omega_q^2 J(Omega_q; a, b),
        # the single-pole counterpart of the exact (R-1) correction 'diff'
        eta, Om2 = _pole1_params(q, rs, zeta=zeta, xspin=xspin)
        Om = np.sqrt(Om2)
        ker[occ] = -eta*Om2*_J_pole(Om, dp[occ], dk[occ])
    if kernel in ("lorentzian_same", "lorentzian_opp"):
        # Sector-resolved Lorentzian kernel: the same-sign sector
        # carries the pi*Lambda/q^3 tail, opposite-sign sector does not
        w0 = lam*q
        jv = w0*w0*_J_pole(w0, dp[occ], dk[occ])
        mask = (dp[occ]*dk[occ] > 0) if kernel == "lorentzian_same" else \
               (dp[occ]*dk[occ] < 0)
        jv = np.where(mask, jv, 0.0)
        ker[occ] = jv
    if kernel == "lorentzian":
        # finite-Lambda frequency kernel: (Lam q)^2 J(Lam q; a, b);
        # Lam -> inf recovers the static W exactly.
        w0 = lam*q
        ker[occ] = (w0*w0)*_J_pole(w0, dp[occ], dk[occ])
    if kernel in ("poleN", "poleNerr"):
        # positive N-pole spectral quadrature (mc/pole_rule.py):
        #   poleN    : correction  -sum_j eta_j t_j J(sqrt(t_j); a, b)
        #   poleNerr : paired error Delta X_N = corr_N - corr_exact per
        #              sample; static-mass matching cancels the leading
        #              1/(|a|+|b|) Fermi-surface singularity, so this
        #              difference estimator has a much reduced (not proven finite) variance.
        if zeta != 0.0 or xspin != 1.0:
            raise NotImplementedError("poleN/poleNerr: pole_rule builds moments from the unpolarized "
                                      "Lindhard function; polarized N-pole rules are not implemented "
                                      "(use kernel='diff' or 'exact' for zeta != 0).")
        keyp = (q, rs, pole_N)
        if keyp not in _pole_cache:
            from pole_rule import pole_rule
            _pole_cache.clear()
            _pole_cache[keyp] = pole_rule(q, rs, pole_N)[:2]
        tn, wn = _pole_cache[keyp]
        pc = np.zeros(occ.sum())
        for tj, wj in zip(tn, wn):
            pc += -wj*tj*_J_pole(np.sqrt(tj), dp[occ], dk[occ])
        if kernel == "poleN":
            ker[occ] = pc
        else:
            key = (q, rs, zeta, xspin)
            if key not in _kern_cache:
                _kern_cache.clear()
                _kern_cache[key] = make_dyn_kernel(q, rs, zeta=zeta,
                                                   xspin=xspin)
            ker[occ] = pc - _dyn_correction(dp[occ], dk[occ],
                                            _kern_cache[key])
    if kernel in ("dyn", "diff"):
        key = (q, rs, zeta, xspin)
        if key not in _kern_cache:
            _kern_cache.clear()
            _kern_cache[key] = make_dyn_kernel(q, rs, zeta=zeta,
                                               xspin=xspin)
        ker[occ] = ker[occ] + _dyn_correction(dp[occ], dk[occ],
                                              _kern_cache[key])
    out[occ] = Dp[occ]*Dk[occ]*ker[occ]/np.maximum(r2[occ], 1e-300) \
        / (fp[occ]*fk[occ])
    return out/(16*np.pi)


def estimate(q, n_total, rng, n_batch=32, return_batch_means=False, robust_guard=False, **kw):
    """(mean, sigma) from n_batch batch means.

    By default the plain batch mean is returned, so the estimand does not
    depend on the data: with ``robust_guard=False`` (the default) this is the
    ordinary unbiased sample mean and ``sigma`` its batch standard error,
    which for heavy-tailed kernels is a scatter measure, not a confidence
    interval.  Passing ``robust_guard=True`` restores the earlier
    data-dependent switch to the batch median when mean and median disagree;
    that estimand is a median, not a mean, and should not be quoted as one.

    kernel: 'static', 'dyn', 'diff', 'lorentzian' (finite-Lambda Lorentzian kernel,
    pass lam=...), 'lorentzian_same'/'lorentzian_opp' (sector-resolved lorentzian),
    'pole1' (moment-conserving one-pole correction; zeta=0 and xspin=1
    only; unsupported polarized calls raise explicitly).
    Set ``return_batch_means=True`` to append a copy of the raw batch-mean
    vector to the return tuple.  The default two-value return is unchanged;
    production drivers can opt in when they need covariance and heavy-tail
    diagnostics without changing the sampled estimator.
    """
    m = max(1, n_total//n_batch)
    means = np.array([xc_samples(q, m, rng, **kw).mean()
                      for _ in range(n_batch)])
    # sig is the batch standard error; for the heavy-tailed kernels here it is a
    # scatter measure, not a confidence half-width.
    mu, sig = means.mean(), means.std(ddof=1)/np.sqrt(n_batch)
    mom = np.median(means)
    mom_sig = 1.4826*np.median(np.abs(means - mom))/np.sqrt(n_batch)
    # heavy-tail guard: if plain and robust estimates disagree strongly,
    # report the (more conservative) robust pair
    if robust_guard and mom_sig > 0 and abs(mu - mom) > 4*mom_sig:
        # reporting the median: its SE exceeds the mean's by sqrt(pi/2)
        # under near-Gaussian batches.
        result = (mom, max(sig, mom_sig*np.sqrt(np.pi/2)))
    else:
        # the reported sigma is max(batch SE, 0.7 x MAD-scale), a deliberately
        # conservative scatter measure; it is not a confidence half-width
        result = (mu, max(sig, 0.7*mom_sig))
    if return_batch_means:
        return result + (means.copy(),)
    return result


def run_to_precision(q, abs_err=None, rel_err=None, n0=200_000,
                     n_max=50_000_000, seed=0, verbose=False, **kw):
    """Grow samples until a fresh confirmation estimate meets the target.

    Returns ``(mean, sigma, n)`` only after the independently drawn
    confirmation stage also satisfies the requested precision.  If the target
    is still missed at ``n_max``, raises ``RuntimeError`` instead of silently
    returning a nonconverged estimate.
    """
    if rel_err is not None and abs_err is None:
        raise ValueError("run_to_precision: a relative target alone is undefined where the "
                         "integrand crosses zero; supply abs_err as a floor.")

    def target(mu):
        if rel_err is not None and abs_err is not None:
            return max(rel_err*abs(mu), abs_err)
        if abs_err is not None:
            return abs_err
        return np.inf

    rng = np.random.default_rng([seed, int(q*10**6)])
    n = int(n0)
    while True:
        mu, sig = estimate(q, n, rng, **kw)
        tgt = target(mu)
        if verbose:
            print(f"  q={q:.3f} n={n:.2e} X_c={mu:.6f}({sig:.6f})",
                  flush=True)
        if sig <= tgt:
            # Discard the trigger-stage estimate and confirm at the same fixed
            # N with fresh draws.  This removes the trigger sample but not the
            # selection: only runs whose confirmation passes are returned, so
            # published numbers should come from fixed-budget replication
            # (see README) rather than from this loop.
            mu2, sig2 = estimate(q, n, rng, **kw)
            tgt2 = target(mu2)
            if verbose:
                print(f"    confirmation X_c={mu2:.6f}({sig2:.6f}); "
                      f"target={tgt2:.6f}", flush=True)
            if sig2 <= tgt2:
                return mu2, sig2, n
            mu, sig, tgt = mu2, sig2, tgt2
        if n >= n_max:
            raise RuntimeError(
                f"precision target not met at n_max={n_max}: "
                f"q={q}, sigma={sig:.6g}, target={tgt:.6g}")
        ratio = sig/tgt if tgt > 0 else np.inf
        n = min(int(n*min(4.0, max(2.0, ratio**2))), int(n_max))


def validate(n=2_000_000, seed=1):
    """Quick bare-limit checks against closed forms (see also run_scan.py)."""
    rng = np.random.default_rng(seed)
    print(f"c0 (exact) = {C0_EXACT:.5f}")
    ok = True
    for q in (0.1, 0.2, 0.4):
        mu, sig = estimate(q, n, rng)
        pull = (mu - C0_EXACT*q)/sig
        print(f"  q={q}: X_c = {mu:.5f} +- {sig:.5f} | c0*q = "
              f"{C0_EXACT*q:.5f}  ({pull:+.1f} sigma)")
        ok &= abs(pull) < 4
    # large-q tail C_inf/q^4, C_inf = 2 pi/9
    Cinf = 2*np.pi/9
    for q in (4.0, 6.0):
        mu, sig = estimate(q, n, rng)
        print(f"  q={q}: X_c*q^4 = {mu*q**4:.4f} +- {sig*q**4:.4f} | "
              f"C_inf = {Cinf:.4f}  (finite-q corrections expected)")
    print("PASS" if ok else "FAIL (small-q slope)")


if __name__ == "__main__":
    validate()
