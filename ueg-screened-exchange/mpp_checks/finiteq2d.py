"""Finite-q first-order exchange polarization of the two-dimensional gas (SI Sec. S6).

R(q,z) = -Pi^(1)(q, i q z)/[q_TF (r_s/sqrt2)] = -(1/(8 pi^2)) Int d^2p d^2k Delta_p Delta_k Re[(g_p - g_k)^2]/|p - k|,
g_p = (i nu - delta_p)^-1, delta_p = p.q + q^2/2, nu = q z, k_F = 1, from SI Eq. (S6) with d = 2.

Method (independent of the crescent quadrature of the SI):
  * (g_p - g_k)^2 = q^2 (p_x - k_x)^2 g_p^2 g_k^2; the joint map (p, k) -> (-p-q, -k-q) leaves the integrand
    invariant, so the p integral runs over the crescent C+ = {|p|<1<|p+q|} only, times 2.
  * For fixed p the k integral is done in polar coordinates about p, k = p + s e_alpha.  The factor
    (p_x-k_x)^2/|p-k| s ds = s^2 cos^2(alpha) ds removes the Coulomb singularity, and the s integral
    of s^2 (A - c s)^-2, A = i nu - delta_p, c = q cos(alpha), is elementary.  Delta_k = 1(|k|<1) - 1(|k+q|<1)
    gives the exit chord of the unit disk (all alpha, periodic trapezoid rule) minus the chord through the
    shifted disk (an arc of directions, mapped so that the chord length is cos(theta), Gauss--Legendre).
  * The p integral uses the crescent coordinates of SI Sec. S6 with Gauss--Legendre nodes.
"""
import numpy as np

def Phi(s, A, c):
    """Int_0^s t^2 (A - c t)^-2 dt, complex A with Im A > 0, real c."""
    s, A, c = np.broadcast_arrays(np.asarray(s, float), np.asarray(A, complex), np.asarray(c, float))
    x = c * s / A
    out = np.empty(s.shape, complex)
    small = np.abs(x) < 0.25
    # series: sum_n (n+1) c^n s^(n+3)/((n+3) A^(n+2)) = s^3/A^2 sum_n (n+1) x^n/(n+3)
    xs = x[small]; acc = np.zeros(xs.shape, complex); p = np.ones(xs.shape, complex)
    for n in range(60):
        acc += (n + 1) / (n + 3) * p
        p = p * xs
    out[small] = s[small] ** 3 / A[small] ** 2 * acc
    b = ~small
    Ab, cb, sb = A[b], c[b], s[b]
    u1 = Ab - cb * sb
    out[b] = -(1 / cb ** 3) * (-Ab ** 2 / u1 - 2 * Ab * (np.log(u1) - np.log(Ab)) + u1)
    return out

def R_finite(q, z, n_phi=96, n_t=48, n_alpha=2048, n_theta=256):
    nu = q * z
    # outer nodes on C+, phi in [0, phi0] doubled by the reflection y -> -y
    phi0 = np.arccos(-q / 2)
    xg, wg = np.polynomial.legendre.leggauss(n_phi)
    # two panels in phi: [0, phi0/2], [phi0/2, phi0]
    phis = np.concatenate([(xg + 1) / 2 * phi0 / 2, phi0 / 2 + (xg + 1) / 2 * phi0 / 2])
    wphi = np.concatenate([wg * phi0 / 4, wg * phi0 / 4])
    xt, wt = np.polynomial.legendre.leggauss(n_t)
    tt = (xt + 1) / 2; wtt = wt / 2
    rho_star = -q * np.cos(phis) + np.sqrt(1 - q * q * np.sin(phis) ** 2)
    rho = rho_star[:, None] + (1 - rho_star[:, None]) * tt[None, :]
    wp = (wphi[:, None] * wtt[None, :]) * rho * (1 - rho_star[:, None]) * 2      # x2 for phi < 0
    px = (rho * np.cos(phis)[:, None]).ravel(); py = (rho * np.sin(phis)[:, None]).ravel(); wp = wp.ravel()
    dp = px * q + q * q / 2
    A = 1j * nu - dp                                   # (Np,)
    gp2 = 1 / A ** 2
    # unit-disk chord, trapezoid in alpha
    al = 2 * np.pi * np.arange(n_alpha) / n_alpha; wal = 2 * np.pi / n_alpha
    ca, sa = np.cos(al), np.sin(al)
    I = np.zeros(px.shape, complex)
    chunk = max(1, 4_000_000 // n_alpha)
    for i0 in range(0, px.size, chunk):
        sl = slice(i0, i0 + chunk)
        pe = px[sl, None] * ca[None, :] + py[sl, None] * sa[None, :]
        s1 = -pe + np.sqrt(pe * pe + 1 - (px[sl, None] ** 2 + py[sl, None] ** 2))
        I[sl] += wal * np.sum(ca ** 2 * Phi(s1, A[sl, None], q * ca[None, :]), axis=1)
    # shifted-disk chord: P = p + q x^, |P| > 1; directions alpha = alpha_c + psi, sin(psi) = sin(theta)/|P|
    xth, wth = np.polynomial.legendre.leggauss(n_theta)
    th = xth * np.pi / 2; wthe = wth * np.pi / 2
    chunk = max(1, 4_000_000 // n_theta)
    for i0 in range(0, px.size, chunk):
        sl = slice(i0, i0 + chunk)
        Px, Py = px[sl] + q, py[sl]
        Pn = np.hypot(Px, Py)
        alc = np.arctan2(-Py, -Px)
        sinpsi = np.sin(th)[None, :] / Pn[:, None]
        psi = np.arcsin(sinpsi); cospsi = np.sqrt(1 - sinpsi ** 2)
        alpha = alc[:, None] + psi
        dpsi = (np.cos(th)[None, :] / Pn[:, None]) / cospsi
        sp_, sm_ = Pn[:, None] * cospsi + np.cos(th)[None, :], Pn[:, None] * cospsi - np.cos(th)[None, :]
        c = q * np.cos(alpha)
        val = np.cos(alpha) ** 2 * (Phi(sp_, A[sl, None], c) - Phi(sm_, A[sl, None], c)) * dpsi
        I[sl] -= np.sum(wthe[None, :] * val, axis=1)
    T = 2 * np.sum(wp * np.real(gp2 * q * q * I))
    return -T / (8 * np.pi ** 2)

def R2(z):
    s = np.sqrt(1 + z * z); r = s - z; a = z / s
    return 2 / np.pi * ((1 - z / s) - z / s ** 3 + 2 * a * a * r * np.arctan(r) / (1 + r * r))

