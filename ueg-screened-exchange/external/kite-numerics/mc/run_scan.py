"""
Production driver for the Monte Carlo kernel X_c(q).

Modes
-----
1) q-scan at fixed one-sigma precision:
     python mc/run_scan.py scan --rel-err 0.005 --out data/xc_scan.csv
   Writes (q, X_c, err, n_samples) and prints the bare check
     D0 = Int X_c dq  vs  the closed form 0.2499018971,
   with the quadrature error propagated from the MC errors.

2) screened static factor vs density (physical Lindhard screening):
     python mc/run_scan.py screened --rs 0.5 1 2 5 10 20 50 --rel-err 0.01 \
            --out data/d0scr.csv
   Computes D0_scr(rs) = Int X_c(q) q^2/(q^2 + qTF^2 fL(q/2)) dq with
   qTF^2 = 4 alpha_L rs / pi, alpha_L = (4/9pi)^(1/3), reusing one q-scan
   (the screening enters only through the weight, so a single X_c table
   serves all rs -> fully correlated errors are propagated exactly).

3) convergence audit (half-sample + seed splitting):
     python mc/run_scan.py audit --q 1.0 --n 4000000
   Reports full vs half-sample estimates and the pull distribution over
   independent seeds: a standard check that the quoted sigma is honest.
"""
import argparse
import sys
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from xc_mc import run_to_precision, estimate, xc_samples, D0_EXACT

ALPHA_L = (4/(9*np.pi))**(1/3.)


def f_lindhard(x):
    x = np.asarray(x, float)
    out = np.full_like(x, 0.5)
    m = np.abs(1 - x) > 1e-12
    out[m] += (1 - x[m]**2)/(4*x[m])*np.log(np.abs((1 + x[m])/(1 - x[m])))
    return out


def default_qgrid():
    """Graded grid: dense at small q (slope/curvature), through the 2kF
    feature, log-spaced tail."""
    return np.unique(np.concatenate([
        np.linspace(0.05, 0.5, 10),
        np.linspace(0.55, 1.8, 13),
        np.linspace(1.85, 2.3, 10),          # Kohn/Friedel region
        np.geomspace(2.4, 12.0, 12),
    ]))


def scan(qs, rel_err, abs_floor, seed, out, verbose=True):
    rows = []
    for q in qs:
        mu, sig, n = run_to_precision(q, rel_err=rel_err, abs_err=abs_floor,
                                      seed=seed, verbose=verbose)
        rows.append((q, mu, sig, n))
    arr = np.array(rows)
    hdr = "q,Xc,err,n_samples"
    np.savetxt(out, arr, delimiter=",", header=hdr, comments="")
    print(f"wrote {out}")
    return arr


def integrate_with_errors(q, y, e, weight=None, n_resample=300, seed=11,
                          kappa2=None, _noise=None):
    """Int weight(q) X_c(q) dq on [q0,qmax] by cubic-spline quadrature,
    plus analytic endpoint tails (exact small-q slope c0; C/q^4 large q).

    Error budget (returned as total, with parts printed by callers):
      * statistical: std of the spline integral over n_resample synthetic
        data sets y + N(0,e)  (exact propagation of independent MC errors
        through the spline weights);
      * grid systematic: |spline - trapezoid|, an estimate of the
        discretization error of the q-grid (this term is what a plain
        trapezoid-with-error-bars analysis silently omits).
    """
    from scipy.interpolate import CubicSpline
    from xc_mc import C0_EXACT
    w = np.ones_like(q) if weight is None else weight(q)
    f = w*y

    def spline_int(vals):
        return CubicSpline(q, w*vals).integrate(q[0], q[-1])

    val = spline_int(y)
    trap = np.trapezoid(f, q)
    sys_grid = abs(val - trap)/2.0
    rng = np.random.default_rng(seed)
    if _noise is None:
        reals = np.array([spline_int(y + rng.normal(0, e))
                          for _ in range(n_resample)])
    else:                       # shared noise realizations -> cross-rs cov
        reals = np.array([spline_int(y + dz*e) for dz in _noise])
    stat = reals.std(ddof=1)
    # endpoint tails.
    # small q: exact closed form of Int_0^{q0} c0 q w(q) dq. For the bare
    # weight this is c0 q0^2/2; for the screened weight the caller passes
    # kappa2 = q_TF^2 and we use (c0/2)[q0^2 - k2 log(1+q0^2/k2)], which
    # retains the screening logarithm (the midpoint rule w(q0/2) used
    # avoids a factor-of-two endpoint underestimate from midpoint sampling.
    if kappa2 is None:
        val += 0.5*C0_EXACT*q[0]**2
    else:
        val += 0.5*C0_EXACT*(q[0]**2 - kappa2*np.log1p(q[0]**2/kappa2))
    # large q: exact C_inf = 2 pi/9 (weight -> 1 there for all forms used);
    # the last-point estimate is printed by callers as a consistency check
    # only, so the tail no longer correlates with the interior quadrature.
    C_INF = 2*np.pi/9
    val += C_INF/(3*q[-1]**3)
    return (val, np.hypot(stat, sys_grid)) if _noise is None else \
        (val, np.hypot(stat, sys_grid), reals + (val - reals.mean()), sys_grid)


def screened_from_scan(arr, rs_list, out):
    q, y, e = arr[:, 0], arr[:, 1], arr[:, 2]
    rows, reals_all, sys_all = [], [], []
    noise = np.random.default_rng(11).normal(size=(300, len(q)))
    for rs in rs_list:
        qtf2 = 4*ALPHA_L*rs/np.pi

        def wgt(qq):
            return qq**2/(qq**2 + qtf2*f_lindhard(qq/2))
        val, err, reals, sg = integrate_with_errors(
            q, y, e, weight=wgt, kappa2=qtf2, _noise=noise)
        rows.append((rs, val, err))
        reals_all.append(reals)
        sys_all.append(sg)
        print(f"rs={rs:7.3f}: D0_scr = {val:.6f} +- {err:.6f}")
    np.savetxt(out, np.array(rows), delimiter=",",
               header="rs,D0scr,err", comments="")
    # full cross-r_s covariance: statistical part from shared-noise
    # realizations; grid systematic treated as one fully correlated mode
    # (conservative: same-sign discretization error across r_s).
    R = np.array(reals_all)
    sg = np.array(sys_all)
    cov = np.cov(R) + np.outer(sg, sg)
    np.savez(out + ".cov.npz", rs=np.array(rs_list), cov=cov)
    print(f"wrote {out} and {out}.cov.npz (full cross-r_s covariance; "
          f"pass to fit_functional for a GLS fit)")


def audit(q, n, seed):
    rng = np.random.default_rng(seed)
    mu_full, s_full = estimate(q, n, rng)
    rng = np.random.default_rng(seed)
    mu_half, s_half = estimate(q, n//2, rng)
    print(f"q={q}: full N  : {mu_full:.6f} +- {s_full:.6f}")
    print(f"q={q}: half N  : {mu_half:.6f} +- {s_half:.6f}   "
          f"pull={(mu_full-mu_half)/np.hypot(s_full,s_half):+.2f} sigma")
    pulls = []
    for s in range(16):
        r = np.random.default_rng([seed + 1, s])
        m, sg = estimate(q, n//8, r)
        pulls.append((m - mu_full)/np.hypot(sg, s_full))
    pulls = np.array(pulls)
    print(f"seed-split pulls: mean={pulls.mean():+.2f}, "
          f"std={pulls.std(ddof=1):.2f}  (should be ~0 and ~1)")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("scan")
    s.add_argument("--rel-err", type=float, default=0.01)
    s.add_argument("--abs-floor", type=float, default=2e-4)
    s.add_argument("--seed", type=int, default=0)
    s.add_argument("--out", default="data/xc_scan.csv")

    c = sub.add_parser("screened")
    c.add_argument("--rs", type=float, nargs="+",
                   default=[0.2, 0.5, 1, 2, 5, 10, 20, 50, 100])
    c.add_argument("--scan-csv", default=None,
                   help="reuse an existing q-scan CSV")
    c.add_argument("--rel-err", type=float, default=0.01)
    c.add_argument("--abs-floor", type=float, default=2e-4)
    c.add_argument("--seed", type=int, default=0)
    c.add_argument("--out", default="data/d0scr.csv")

    dyn = sub.add_parser("dyn", help="dynamical (AC, RPA-screened) kite "
                         "via the correlated difference estimator: "
                         "D0dyn(rs) = D0_exact + Int [Xc_dyn - Xc_bare] dq")
    dyn.add_argument("--rs", type=float, nargs="+",
                     default=[0.2, 0.5, 1, 2, 5, 8, 12, 20])
    dyn.add_argument("--rel-err", type=float, default=0.02)
    dyn.add_argument("--abs-floor", type=float, default=3e-4)
    dyn.add_argument("--zeta", type=float, default=0.0,
                     help="spin polarization (0 or 1 validated); the "
                          "bare limit is zeta-independent by construction")
    dyn.add_argument("--seed", type=int, default=0)
    dyn.add_argument("--out", default="data/d0dyn.csv")

    a = sub.add_parser("audit")
    a.add_argument("--q", type=float, default=1.0)
    a.add_argument("--n", type=int, default=4_000_000)
    a.add_argument("--seed", type=int, default=7)

    args = ap.parse_args()
    if args.cmd == "scan":
        arr = scan(default_qgrid(), args.rel_err, args.abs_floor,
                   args.seed, args.out)
        d0, d0e = integrate_with_errors(arr[:, 0], arr[:, 1], arr[:, 2])
        pull = (d0 - D0_EXACT)/d0e
        print(f"\nD0 (MC)    = {d0:.6f} +- {d0e:.6f}")
        print(f"D0 (exact) = {D0_EXACT:.6f}   pull = {pull:+.2f} sigma")
    elif args.cmd == "screened":
        if args.scan_csv:
            arr = np.loadtxt(args.scan_csv, delimiter=",", skiprows=1)
        else:
            arr = scan(default_qgrid(), args.rel_err, args.abs_floor,
                       args.seed, "data/xc_scan.csv")
        screened_from_scan(arr, args.rs, args.out)
    elif args.cmd == "dyn":
        from xc_mc import run_to_precision as rtp
        qs = np.geomspace(0.02, 6.0, 20)
        rows = []
        spins = ([(1.0, 1.0)] if args.zeta == 0.0 else
                 [(s, (1 + s*args.zeta)**(1/3.)) for s in (+1.0, -1.0)
                  if (1 + s*args.zeta) > 1e-12])
        for rs in args.rs:
            dv, de2 = 0.0, 0.0
            for isp, (_, x) in enumerate(spins):
                wsp = 0.5*x**3 if args.zeta != 0.0 else 1.0
                ys, es = [], []
                for q in qs:
                    mu, sig, _ = rtp(q, rel_err=args.rel_err,
                                     abs_err=args.abs_floor,
                                     seed=args.seed + isp*10_000_019,
                                     kernel="diff", rs=rs,
                                     zeta=args.zeta, xspin=x)
                    ys.append(mu); es.append(sig)
                ys = np.array(ys); es = np.array(es)
                dv += wsp*(np.trapezoid(ys, qs) + 0.5*ys[0]*qs[0])
                tw = np.zeros_like(qs)
                tw[1:] += 0.5*np.diff(qs); tw[:-1] += 0.5*np.diff(qs)
                tw[0] += 0.5*qs[0]   # small-q endpoint enters the value,
                                     # so it must enter the error weight
                de2 += (wsp**2)*np.sum((tw*es)**2)
            de = np.sqrt(de2)
            d0 = D0_EXACT + dv
            rows.append((rs, d0, de))
            print(f"rs={rs:7.3f}: D0dyn = {d0:+.6f} +- {de:.6f}   "
                  f"e_kite_dyn = {6e3/np.pi**3*d0:+7.2f} "
                  f"+- {6e3/np.pi**3*de:.2f} mRy", flush=True)
        np.savetxt(args.out, np.array(rows), delimiter=",",
                   header="rs,D0dyn,err", comments="")
        print(f"wrote {args.out}")
    elif args.cmd == "audit":
        audit(args.q, args.n, args.seed)


if __name__ == "__main__":
    main()
