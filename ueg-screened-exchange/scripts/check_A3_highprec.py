#!/usr/bin/env python3
"""High-precision check of the A3 reductions (companion App. A).

Requires mpmath (optional dependency; see requirements.txt).  Verifies, at 30
significant digits:
  * the two contour identities for the recoil and Lindhard frequency integrals,
  * the Euler integral M of Eq. (MEuler) against (176/75) sqrt(2 pi) G(3/4)/G(1/4),
  * the three A3 pieces and their sum against -16408/3375 C_Gamma.
"""
import mpmath as mp

FAILED = []
mp.mp.dps = 30
A0 = mp.mpf(3) / (32 * mp.pi ** 5)
CG = mp.sqrt(2 * mp.pi) / mp.pi ** 2 * mp.gamma(mp.mpf(3) / 4) / mp.gamma(mp.mpf(1) / 4)
B = lambda a, b: mp.gamma(a) * mp.gamma(b) / mp.gamma(a + b)
V_B = 4 * mp.pi / 3


def R(c):
    c = mp.mpf(c)
    return 2 * (c - mp.log(1 + c)) / c ** 2 if c > mp.mpf("1e-8") else 1 - 2 * c / 3 + c * c / 2


def line(label, a, b, tol=mp.mpf("1e-20")):
    passed = mp.isfinite(a) and mp.isfinite(b) and abs(a - b) < tol * max(1, abs(b))
    if not passed:
        FAILED.append(label)
    ok = "OK " if passed else "!! "
    print(f"{ok}{label:52s} {mp.nstr(a, 22)}  ref {mp.nstr(b, 22)}")


print("== frequency identities (residues) ==")
for m in (mp.mpf("0.3"), mp.mpf(1), mp.mpf(7)):
    lhs = mp.quad(lambda u: (1 + 4 * u * u) * mp.re((1j * u - mp.mpf(1) / 2) ** -4) / (1 + 4 * u * u + m),
                  [0, mp.mpf("0.5"), 2, 10, mp.inf]) / mp.pi
    s = mp.sqrt(1 + m)
    line(f"recoil, m={m}", lhs, -4 * m / (s * (s + 1) ** 4), mp.mpf("1e-18"))
    lindhard_coeff_4 = lambda u: -16 * (12 * u * u - 1) / (15 * (1 + 4 * u * u) ** 3)
    lhs2 = mp.quad(lambda u: lindhard_coeff_4(u) * mp.re((1j * u - mp.mpf(1) / 2) ** -2) * (1 + 4 * u * u) ** 2 / (1 + 4 * u * u + m) ** 2,
                   [0, mp.mpf("0.5"), 2, 10, mp.inf]) / mp.pi
    line(f"Lindhard, m={m}", lhs2, 8 * (2 * s ** 3 - s * s + 4 * s + 1) / (15 * s ** 3 * (s + 1) ** 4), mp.mpf("1e-18"))

print("\n== the three A3 pieces ==")
Fp = mp.mpf(2) / 3 * (2 ** mp.mpf("0.5") * B(mp.mpf(3) / 4, -mp.mpf(1) / 2) - 2 ** mp.mpf("-1.5") * B(mp.mpf(3) / 4, mp.mpf(3) / 2))
line("cross  = 8 pi A0 V_B^2 Fp", 8 * mp.pi * A0 * V_B ** 2 * Fp, -mp.mpf(32) / 15 * CG)
K = -2 * mp.sqrt(2) * B(mp.mpf(7) / 4, mp.mpf(3) / 2)
line("recoil = 4 pi A0 (104 pi^2/27)(2/5) K", 4 * mp.pi * A0 * (104 * mp.pi ** 2 / 27) * (mp.mpf(2) / 5) * K, -mp.mpf(416) / 675 * CG)
# w = t^2 removes the inverse-square-root endpoint; tanh-sinh handles (1-t^2)^{3/4}
M = (mp.mpf(32) / 15) * mp.quad(lambda t: (1 - t * t) ** mp.mpf("0.75") * (1 + t * t) ** mp.mpf("-3.25")
                                * (2 - t ** 2 + 4 * t ** 4 + t ** 6),
                                [0, 1], method="tanh-sinh")
line("M (Euler integral)", M, mp.mpf(176) / 75 * mp.sqrt(2 * mp.pi) * mp.gamma(mp.mpf(3) / 4) / mp.gamma(mp.mpf(1) / 4))
line("Lindhard = -(48/5) pi^3 A0 M", -(mp.mpf(48) / 5) * mp.pi ** 3 * A0 * M, -mp.mpf(264) / 125 * CG)

print("\n== total ==")
tot = 8 * mp.pi * A0 * V_B ** 2 * Fp + 4 * mp.pi * A0 * (104 * mp.pi ** 2 / 27) * (mp.mpf(2) / 5) * K - (mp.mpf(48) / 5) * mp.pi ** 3 * A0 * M
line("A3", tot, -mp.mpf(16408) / 3375 * CG)
print("C_Gamma =", mp.nstr(CG, 18))

if FAILED:
    raise SystemExit("A3 high-precision checks failed: " + ", ".join(FAILED))
print("All A3 high-precision checks passed.")
