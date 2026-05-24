"""
RPA (ring) correlation energy per electron of the unpolarized UEG, from
the adiabatic-connection fluctuation--dissipation formula evaluated
entirely on the imaginary frequency axis (no plasmon-pole hunting; the
integrand is smooth there):

  e_c^RPA(rs) [Ry] = (3 k_F^2 / 2 pi) Int_0^inf dx x^2 Int_0^inf du
                     [ ln(1 + c(x,u)) - c(x,u) ],

  c(x,u) = (q_TF^2/x^2) G(x,u),   q_TF^2 = 4 alpha_L rs/pi,
  k_F = 1/(alpha_L rs),           G = reduced imaginary-axis Lindhard
                                      (analytic/lindhard.py, closed form),

with x = q/k_F and u the frequency in units where eps_k = k^2/2, k_F = 1.
The coupling-constant integral is already closed form in the
ln(1+c) - c combination.

Validation (run this file):
  * the exact Gell-Mann--Brueckner high-density limit
      e_c^RPA -> (2(1-ln2)/pi^2) ln rs - 0.1423  Ry
    (log coefficient checked to <0.3%, constant to ~2 mRy);
  * quadrature-refinement error estimate at each rs.

CLI:  python analytic/rpa_energy.py --rs 0.5 1 2 5 10 --out data/rpa.csv
writes (rs, e_c_RPA [Ry], quadrature error estimate).

Scope: zeta = 0 (unpolarized). The kite suite's polarization statements
(bare diagram polarization-independent) do not require the polarized RPA.
"""
import argparse
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lindhard import G, ALPHA_L


def _panels(edges, nodes, weights):
    e = np.asarray(edges)
    x = (0.5*(e[1:, None] - e[:-1, None])*nodes[None, :]
         + 0.5*(e[1:, None] + e[:-1, None])).ravel()
    w = (0.5*(e[1:, None] - e[:-1, None])*weights[None, :]).ravel()
    return x, w


def _c_two_spin(x, u, qtf2, zeta):
    """dimensionless screening c = v|chi0| with the two-spin polarization:
    chi0 for Fermi radius x_s scales as x_s * chi0^(1)(q/x_s, u/x_s^2)."""
    if zeta == 0.0:
        return qtf2/x**2*G(x, u)
    c = np.zeros_like(np.atleast_1d(u), dtype=float)
    for sgn in (+1.0, -1.0):
        xs = (1.0 + sgn*zeta)**(1/3.)
        if xs > 1e-12:
            c = c + qtf2/x**2*0.5*xs*G(x/xs, np.asarray(u, float)/xs**2)
    return c


def _u_integral(x, qtf2, ng, nseg, zeta=0.0):
    """Int_0^inf du [ln(1+c)-c] at fixed x; u-structure scale ~ x(x_max+x/2)."""
    xg, wg = np.polynomial.legendre.leggauss(ng)
    x_up = (1.0 + zeta)**(1/3.)
    ue = x*(x_up + x/2)
    edges = np.concatenate([[0.0], ue*np.geomspace(1e-4, 1.0, nseg)])
    u, w = _panels(edges, xg, wg)
    # mapped tail u = ue(1 + t/(1-t)):  du = ue dt/(1-t)^2
    t, wt = _panels(np.linspace(0, 1, nseg//2 + 1)[:-1].tolist()
                    + [1 - 1e-9], xg, wg)
    u_t = ue*(1 + t/(1 - t))
    w_t = wt*ue/(1 - t)**2
    u = np.concatenate([u, u_t]); w = np.concatenate([w, w_t])
    c = _c_two_spin(x, u, qtf2, zeta)
    return float(np.sum(w*(np.log1p(c) - c)))


def e_rpa(rs, ng=16, nx=28, nseg=12, zeta=0.0):
    """RPA correlation energy per electron in Ry; zeta via the two-spin
    dielectric (Kohn breakpoints at 2*x_sigma for both spins)."""
    kf = 1.0/(ALPHA_L*rs)
    qtf2 = 4*ALPHA_L*rs/np.pi
    xg, wg = np.polynomial.legendre.leggauss(ng)
    x_lo = min(1.0, np.sqrt(qtf2))*1e-3
    xkohn = sorted({2.0*(1.0 + sgn*zeta)**(1/3.) for sgn in (+1.0, -1.0)
                    if (1.0 + sgn*zeta) > 1e-12})
    edges = np.concatenate([x_lo*np.geomspace(1, xkohn[0]/x_lo, nx),
                            *[np.linspace(a, b, 5)[1:] for a, b in
                              zip(xkohn, xkohn[1:])],
                            np.linspace(xkohn[-1], 6.0, 7)[1:]])
    xs, wx = _panels(edges, xg, wg)
    val = sum(wxi*xi**2*_u_integral(xi, qtf2, ng, nseg, zeta)
              for xi, wxi in zip(xs, wx))
    # mapped x tail beyond 6
    t, wt = _panels([0, 1 - 1e-9], xg, wg)
    xt = 6.0/(1 - t)
    val += float(np.sum(wt*6.0/(1 - t)**2
                        * np.array([x*x*_u_integral(x, qtf2, ng, nseg, zeta)
                                    for x in xt])))
    return 3*kf**2/(2*np.pi)*val


def e_rpa_with_error(rs, zeta=0.0):
    lo = e_rpa(rs, zeta=zeta)
    hi = e_rpa(rs, ng=24, nx=40, nseg=18, zeta=zeta)
    return hi, abs(hi - lo)


def validate():
    A = 2*(1 - np.log(2))/np.pi**2          # 0.0621814 (exact GMB, Ry)
    e1, _ = e_rpa_with_error(0.01)
    e2, _ = e_rpa_with_error(0.02)
    slope = (e2 - e1)/np.log(2.0)
    const = e1 - A*np.log(0.01)
    ok1 = abs(slope - A)/A < 0.02
    ok2 = abs(const + 0.1423) < 0.004
    print(f"[{'PASS' if ok1 else 'FAIL'}] GMB log coefficient: "
          f"{slope:.5f} vs {A:.5f}")
    print(f"[{'PASS' if ok2 else 'FAIL'}] GMB RPA constant: "
          f"{const:.4f} vs -0.1423 Ry")
    for rs in (1.0, 4.0):
        e, de = e_rpa_with_error(rs)
        print(f"  e_c_RPA(rs={rs}) = {e:.5f} +- {de:.1e} Ry "
              f"({1e3*e/2:.2f} mHa)")
    return ok1 and ok2


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--rs", type=float, nargs="+", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    if a.rs is None:
        raise SystemExit(0 if validate() else 1)
    rows = []
    for rs in a.rs:
        e, de = e_rpa_with_error(rs)
        rows.append((rs, e, de))
        print(f"rs={rs:8.3f}: e_c_RPA = {e:.6f} +- {de:.1e} Ry")
    if a.out:
        np.savetxt(a.out, np.array(rows), delimiter=",",
                   header="rs,ecRPA_Ry,quad_err", comments="")
        print(f"wrote {a.out}")
