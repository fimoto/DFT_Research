"""Static (z = 0) finite-q first-order exchange polarization of the 2D gas, R(q,0), with the
same integrand and product quadrature as finiteq2d_static.py, but with the phi panels ending at
the crescent tips cos(phi0) = -q/2 (where E_p = p.q + q^2/2 vanishes and the crescent closes),
nodes graded towards the tips, and E_p and 1 - rho* evaluated without cancellation near the tips."""
import numpy as np, sys, time
def gl01(n):
    x, w = np.polynomial.legendre.leggauss(n); return (x + 1) / 2, w / 2
def nodes(q, nphi, nt, m):
    p0 = np.arccos(-q / 2); s0 = np.sin(p0)
    u, wu = gl01(nphi)
    dl, wd = [], []
    # each panel: delta = L u^m measured from the tip, sigma = -1 (C+, phi = p0 - delta) or +1 (C-, phi = p0 + delta)
    rows = []
    for sigma, L in ((-1, p0), (+1, np.pi - p0)):
        d = L * u ** m; w = L * m * u ** (m - 1) * wu
        rows.append((sigma, d, w))
    t, wt = gl01(nt)
    P, W, E = [], [], []
    for sigma, d, w in rows:
        phi = p0 + sigma * d
        c = q * np.sin(d / 2) ** 2 - sigma * s0 * np.sin(d)          # cos(phi) + q/2, accurate
        cphi = c - q / 2
        rs = -q * cphi + np.sqrt(1 - q * q * np.sin(phi) ** 2)
        one_m_rs = 2 * q * rs * c / (rs + 1 - q * q)                 # 1 - rho*, accurate
        for ysgn in (1, -1):
            PH = phi[:, None]; T = t[None, :]
            rho_m1 = -one_m_rs[:, None] * (1 - T)                    # rho - 1
            rho = 1 + rho_m1
            Ep = q * (rho_m1 * cphi[:, None] + c[:, None])           # = q (rho cos phi + q/2)
            Wt = (w[:, None] * wt[None, :]) * one_m_rs[:, None] * rho
            P.append(np.stack([rho * np.cos(PH), ysgn * rho * np.sin(PH)], -1).reshape(-1, 2))
            W.append(Wt.ravel()); E.append(Ep.ravel())
    return np.concatenate(P), np.concatenate(W), np.concatenate(E)
def R(q, nphi, nt, m, chunk=300):
    P, W, E = nodes(q, nphi, nt, m)
    tot = 0.0
    for i in range(0, len(W), chunk):
        dp = P[i:i + chunk, None, :] - P[None, :, :]
        r = np.sqrt(np.sum(dp * dp, -1)); r[r == 0] = np.inf
        Ei = E[i:i + chunk, None]; Ej = E[None, :]
        f = (2 * np.pi / r) * ((Ej - Ei) / (Ei * Ej)) ** 2          # Re[(g_p - g_k)^2] at nu = 0
        tot += np.sum(W[i:i + chunk, None] * W[None, :] * f)
    return -tot / (2 * np.pi) ** 4 * np.pi
if __name__ == '__main__':
    x = float(sys.argv[1]); m = float(sys.argv[2])
    for s in sys.argv[3:]:
        nphi, nt = map(int, s.split(','))
        t0 = time.time(); v = (R(x, nphi, nt, m) - 2 / np.pi) / (2 / np.pi)
        print(x, m, nphi, nt, f'{v:.6e}', f'{time.time()-t0:.0f}s', flush=True)
