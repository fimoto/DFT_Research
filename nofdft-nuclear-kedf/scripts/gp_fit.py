#!/usr/bin/env python3
"""
Gaussian-process (reduced-rank) fit of the nonlocal-kernel correction delta-g.

Self-contained and deterministic: all inputs are read from ./gp_data/ (staged
alongside this script), no absolute paths, no hidden state. Reproduces the
GP results of the preprint:

  - evidence (log marginal likelihood) selection over (ell, sigma_f, kappa),
  - MAP kernel correction eg_GP.dat (+ local-correction shift),
  - fixed-active-set inverse-Hessian Gaussian approximation
    -> uncertainty band eg_GP_band.dat,
  - A=140 (not used in fitting or model selection) 2-sigma coverage diagnostic.

The Gaussian data blocks (chi values, exact EL gradients, and surface-tail
values) are linear operator observations of delta-g. Energy matching is imposed
as exact linear equalities. Violations of the tail and selected-path energy-rise
target margins enter only at the MAP stage as soft quadratic inequality
penalties; they are not hard inequalities. Their active rows are updated by an
active-set repeated penalized least-squares iteration.

Run:
    cd scripts && python3 gp_fit.py
Requires: numpy. Outputs are written to the working directory.
"""
import os
import numpy as np
import pickle

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gp_data")
def dpath(name): return os.path.join(DATA, name)

# ---- reduced-rank GP basis (Solin-Sarkka): sine eigenfunctions of -d^2/dy^2 on [0,L]
M = 20
L = 12.0
h2m = 20.73553
om = np.pi * np.arange(1, M + 1) / L           # angular frequencies omega_m = sqrt(lambda_m)

with open(dpath("dg_resp.pkl"), "rb") as f:   RESP = pickle.load(f)   # geometry (r, rho, ...)
with open(dpath("sine_resp.pkl"), "rb") as f: SR = pickle.load(f)     # sine-basis chi responses
with open(dpath("stab_rows_sine.pkl"), "rb") as f: _STAB_ALL = pickle.load(f)  # energy-rise-margin rows

def opath(name):
    """Write outputs next to this script (scripts/), independent of CWD."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), name)
# Strict data separation: keep only soft energy-rise-margin rows from the three
# fitting nuclei. The A=140 rows are inactive throughout in any case (removing
# them leaves the MAP and covariance bitwise-identical), but filtering here
# makes it explicit that no A=140-derived information enters fitting or model
# selection.
# The released filtered data contain exactly three such rows, all from A=80.
_TRAIN_NUC = ("O16", "Ca40", "A80")
STAB = [row for row in _STAB_ALL if row[0] in _TRAIN_NUC]

TEACH = {n: dpath("teacher_ddkf_%s.dat" % n) for n in ("O16", "Ca40", "A80", "A140")}
KSF   = {n: dpath("ks_nlf_%s.dat" % n)       for n in ("O16", "Ca40", "A80", "A140")}
PREF  = {"O16": "vdd", "Ca40": "vddCa", "A80": "vddZ"}
BANDS = {"O16": ((4.5, 8.0), (8.0, 11.5)), "Ca40": ((4.5, 8.0), (8.0, 11.5)),
         "A80": ((5.5, 9.0), (9.0, 13.0)), "A140": ((5.5, 9.0), (9.0, 13.0))}
BETA  = {"O16": 2.0, "Ca40": 2.0, "A80": 12.0}
ETAT  = {"O16": 1.2, "Ca40": 0.3, "A80": 0.3}


def load(nuc, need_grad):
    r5, rho5, k5, _ = RESP[nuc]
    t5 = np.loadtxt(TEACH[nuc]); k5f = np.loadtxt(KSF[nuc])
    base = np.loadtxt(dpath("%s_base.dat" % PREF[nuc])) if need_grad else np.loadtxt(dpath("vddY_base.dat"))
    r2 = base[:, 0]
    I = lambda x: np.interp(r2, r5, x)
    rho = I(rho5); tauTF = np.maximum(I(k5f[:, 2]), 1e-30)
    tauvW = I(k5f[:, 3]); tauKS = I(k5f[:, 4]); taudd = I(t5[:, 6])
    Bt = [I(bb) for bb in SR[nuc]]
    chiKS = (tauKS - tauvW) / tauTF; chidd = (taudd - tauvW) / tauTF
    rho53 = np.maximum(rho, 1e-30) ** (5 / 3)
    Xc = np.column_stack([bb / tauTF for bb in Bt] + [-rho53 / tauTF])
    dr = r2[1] - r2[0]; vol = 4 * np.pi * r2 ** 2 * dr
    g = np.array([np.sum(vol * bb) for bb in Bt] + [float(-np.sum(vol * rho53))])
    tgt = float(np.sum(vol * (tauKS - taudd)))
    boostw = 3.2 if nuc in ("O16", "Ca40") else 3.9
    J = 2 * np.abs(chiKS) / (1 + chiKS ** 2) ** 2
    wch = (r2 ** 2 + 0.7) * rho * (1 + 8 * np.exp(-(r2 / boostw) ** 2)) * J ** 2
    (t1, t2), (h1, h2) = BANDS[nuc]
    out = dict(r2=r2, rho=rho, Xc=Xc, ych=chiKS - chidd, chidd=chidd, chiKS=chiKS, wch=wch,
               g=g, tgt=tgt, mB=rho > 1e-6, mE=rho > 1e-5, mT=(r2 > t1) & (r2 < t2),
               mH=(r2 >= h1) & (r2 < h2),
               wT=r2 ** 2 * np.maximum(rho, 1e-30) ** 0.4 / (1 + (chiKS / 8.0) ** 2),
               wEL=(r2 ** 2 + 0.7) * rho)
    if need_grad:
        V0 = base[:, 15]; Vs = np.loadtxt(dpath("%s_shift.dat" % PREF[nuc]))[:, 15] - V0
        Vb = [np.loadtxt(dpath("vddS_%s_b%02d.dat" % (nuc, m)))[:, 15] - V0 for m in range(M)]
        out["V"] = np.column_stack(Vb + [Vs]); out["yv"] = -V0
    return out


train = ["O16", "Ca40", "A80"]
D = {n: load(n, True) for n in train}; D["A140"] = load("A140", False)
nb = M + 1; ntr = 3; nun = nb + ntr
# Soft-penalty weights and target margins: chi is dimensionless and eps_stab is
# in MeV. Finite weights discourage violations but do not guarantee either
# margin; exact energy-matching equalities below are treated separately by KKT.
eta_h = 12.0; margin = 5.0; eps_stab = 2.0

G = np.zeros((ntr, nun)); T = np.zeros(ntr)
for i, n in enumerate(train):
    G[i, :nb] = D[n]["g"]; T[i] = D[n]["tgt"]


def prior_prec(l, sf):
    """Diagonal prior precision from the SE spectral density S(omega_m)."""
    S = sf ** 2 * np.sqrt(2 * np.pi) * l * np.exp(-0.5 * l ** 2 * om ** 2)
    P = np.zeros(nun); P[:M] = 1.0 / np.maximum(S, 1e-300); P[M] = 1e-4; P[M + 1:] = 1e-10
    return P


# ---- Gaussian operator-observation blocks, assembled once
def gauss_blocks():
    RA, RY, RIV = [], [], []
    for i, n in enumerate(train):
        d = D[n]
        for (mask, cols, yv, w, eta, mu) in [
            (d["mB"], "X", d["ych"], d["wch"], 1.0, False),
            (d["mE"], "V", d["yv"], d["wEL"], BETA[n], True),
            (d["mT"], "X", d["ych"], d["wT"], ETAT[n], False)]:
            A = np.zeros((int(mask.sum()), nun))
            A[:, :nb] = d["Xc"][mask] if cols == "X" else d["V"][mask]
            if mu: A[:, nb + i] = -1.0
            y = yv[mask]; wm = w[mask]; s = np.sum(wm * y ** 2)
            RA.append(A); RY.append(y); RIV.append(eta * wm / s)
    return RA, RY, RIV


RA, RY, RIV = gauss_blocks()


def evidence(P, kap):
    """Log marginal likelihood of the LINEAR (Gaussian) blocks under prior P,
    with a global likelihood (noise) scale kappa. Soft quadratic inequality
    penalties are MAP-level terms and are excluded from the evidence, as stated
    in the paper."""
    H0 = np.diag(P).astype(float); b = np.zeros(nun); quad = 0.0; logw = 0.0
    for A, y, iv in zip(RA, RY, RIV):
        H0 += kap * A.T @ (A * iv[:, None]); b += kap * A.T @ (iv * y)
        quad += kap * np.sum(iv * y * y); logw += np.sum(np.log(kap * iv))
    _, ld = np.linalg.slogdet(H0)
    m = np.linalg.solve(H0, b)
    lp = np.sum(np.log(P[:M]))
    return 0.5 * (b @ m - quad) - 0.5 * ld + 0.5 * lp + 0.5 * logw


def map_fit(P, kap, nit=16):
    """MAP fit with exact energy equalities and soft quadratic penalties.

    Each update selects rows that violate either the tail target chi>=margin or
    the selected-path energy-rise target Delta E>=eps_stab at the current
    iterate. With those rows fixed, a penalized least-squares subproblem with
    exact energy equalities is solved by a bordered KKT system. This is an
    active-set repeated penalized least-squares iteration, not a hard-inequality
    optimizer: a target margin may remain violated at the finite penalty weight.

    The released implementation deliberately performs nit=16 fixed updates.
    For the released data (three energy-rise rows, all from A=80), both the
    active set and coefficients cease to change at update 14; updates 15 and 16
    therefore repeat the same quadratic subproblem and leave the result intact.
    """
    RIV2 = [kap * iv for iv in RIV]
    c = np.zeros(nun)
    for _ in range(nit):
        H = np.diag(P).astype(float); b = np.zeros(nun)
        for A, y, iv in zip(RA, RY, RIV2):
            H += A.T @ (A * iv[:, None]); b += A.T @ (iv * y)
        for i, n in enumerate(train):
            d = D[n]
            # Tail active set from the previous coefficient iterate. The
            # normalization sH is recomputed from the currently active rows, so
            # each update is a convex quadratic subproblem rather than a claim
            # that a hard lower bound on chi is enforced.
            chiH = d["chidd"][d["mH"]] + d["Xc"][d["mH"]][:, :nb] @ c[:nb]
            viol = chiH < margin
            if viol.any():
                AH = np.zeros((int(viol.sum()), nun)); AH[:, :nb] = d["Xc"][d["mH"]][viol]
                yH = (margin - d["chidd"][d["mH"]])[viol]; sH = max(np.sum(yH ** 2), 1e-12)
                H += kap * eta_h * (AH.T @ AH) / sH; b += kap * eta_h * (AH.T @ yH) / sH
        # Selected-path energy-rise penalty: weight 10 and normalization
        # 100 MeV^2=(10 MeV)^2. Only rows below the 2 MeV soft target are active;
        # finite weighting does not require the final energy rise to reach it.
        for (nuc, i, q, dEb) in STAB:
            if dEb + float(q @ c[:nb]) < eps_stab:
                row = np.zeros(nun); row[:nb] = q
                H += kap * 10.0 * np.outer(row, row) / 100.0
                b += kap * 10.0 * row * (eps_stab - dEb) / 100.0
        KKT = np.zeros((nun + ntr, nun + ntr))
        KKT[:nun, :nun] = H; KKT[:nun, nun:] = G.T; KKT[nun:, :nun] = G
        sol = np.linalg.solve(KKT, np.concatenate([b, T]))
        c = sol[:nun]
    return c, H


# ---- hyperparameter selection by evidence.
# Grid brackets the reported optimum (ell=1.5, sigma_f=0.3, kappa=300).
ELL   = [0.6, 1.0, 1.5, 2.2]
SIGF  = [0.3, 1.0, 3.0]
KAPPA = [30.0, 100.0, 300.0, 1000.0]
best = None
for l in ELL:
    for sf in SIGF:
        for kap in KAPPA:
            ev = evidence(prior_prec(l, sf), kap)
            if best is None or ev > best[0]:
                best = (ev, l, sf, kap)
ev, l, sf, kap = best
print("SELECTED by evidence: ell=%.2f sigma_f=%.2f kappa=%.0f (logEv=%.1f)" % (l, sf, kap, ev))

P = prior_prec(l, sf)
c, H = map_fit(P, kap)
np.save(opath("gp_map.npy"), c)
Hinv = np.linalg.inv(H)
# Fixed-active-set inverse-Hessian Gaussian approximation, projected onto the
# exact energy-equality subspace. It is not a probability distribution
# conditioned on satisfaction of the soft target margins.
Sig = Hinv - Hinv @ G.T @ np.linalg.inv(G @ Hinv @ G.T) @ G @ Hinv
np.save(opath("gp_cov.npy"), Sig)
cc = c[:nb]
print("local-correction shift = %.5f" % cc[M])

for n in train + ["A140"]:
    d = D[n]; chi_m = d["chidd"] + d["Xc"] @ cc
    nlfm = 1 / (1 + chi_m ** 2); nlfK = 1 / (1 + d["chiKS"] ** 2)
    wd = d["r2"] ** 2 * d["rho"]; md = d["r2"] < (7 if n in ("O16", "Ca40") else 9)
    dist = np.sqrt(np.sum(wd[md] * (nlfm[md] - nlfK[md]) ** 2) / np.sum(wd[md]))
    role = "not used" if n == "A140" else "used"
    print("  %-5s %-9s d_NLF = %.4f" % (n, role, dist))

# A=140 2-sigma coverage diagnostic (parameter uncertainty propagated to chi).
# A=140 was not used in fitting or model selection.
d = D["A140"]
Achi = np.zeros((len(d["r2"]), nun)); Achi[:, :nb] = d["Xc"]
sd = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", Achi, Sig, Achi), 0))
chi_m = d["chidd"] + d["Xc"] @ cc
mm = (d["r2"] > 1) & (d["r2"] < 8) & (d["rho"] > 1e-5)
cov = float(np.mean(np.abs(chi_m - d["chiKS"])[mm] <= 2 * sd[mm]))
print("  A=140 (not used in fit): 2-sigma coverage = %.2f" % cov)
print("    (small coverage indicates that model discrepancy is not represented by this band)")

# Export the fixed-density A=140 NLF diagnostic and its 2-sigma band (Fig. 6,
# right panel). The legacy output filename is retained for compatibility.
# Propagate the chi uncertainty sd through NLF = 1/(1+chi^2).
_nlf_m = 1.0 / (1.0 + chi_m ** 2)
_nlf_K = 1.0 / (1.0 + d["chiKS"] ** 2)
_dnlf_dchi = -2.0 * chi_m / (1.0 + chi_m ** 2) ** 2
_nlf_sd = np.abs(_dnlf_dchi) * sd
np.savetxt(opath("gp_holdout_A140_nlf.dat"),
           np.column_stack([d["r2"], _nlf_m, _nlf_K, 2.0 * _nlf_sd]),
           header="r  NLF_GP_model  NLF_KS  2sigma_band")

# Export the MAP kernel correction and fixed-active-set Gaussian band.
yy = np.linspace(0, 25, 2501)
PHI = np.where(yy[:, None] <= L, np.sqrt(2 / L) * np.sin(np.pi * np.arange(1, M + 1)[None, :] * yy[:, None] / L), 0.0)
dg = PHI @ cc[:M]
sd_dg = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", PHI, Sig[:M, :M], PHI), 0))
np.savetxt(opath("eg_GP.dat"), np.column_stack([yy, dg]))
np.savetxt(opath("eg_GP_band.dat"), np.column_stack([yy, dg, sd_dg]))
with open(opath("gp_shift.txt"), "w") as fh:
    fh.write("%.5f" % cc[M])
print("exported eg_GP.dat (max|delta_g|=%.3f), eg_GP_band.dat, gp_map.npy, gp_cov.npy" % abs(dg).max())
