"""Main-text Sec. 3.3 (low density) and Table 3: the leading low-density term of the
kite and its comparison with the direct values and with BRM.

Reproduces: A_1 from the (lambda, x) double integral of main Eq. (53) and from the
gamma-function formula; q_*^-1 = 1.0311348 r_s^-1/4; the Table 3 coefficients of
e_kite and F = (7/8) e_kite; A_1/q_* at r_s = 1e8; the momentum scale q_* k_F at
r_s = 2e5 and 1e6; and, from the stored data of the repository, the direct values
at r_s = 1e8, 2e5, 1e6 (production/raw/tableI) and the BRM values quoted in the main
text (production/reference_brm.json, an external reference, not computed here).
"""
import json, os
import mpmath as mp
from common import Recorder, HERE

R = Recorder("c11_lowdensity_kite", "low-density kite (main Sec. 3.3, Table 3)")
root = os.path.dirname(HERE)
mp.mp.dps = 30
aL = (4 / (9 * mp.pi)) ** (mp.mpf(1) / 3)
A1c = -4 * mp.gamma(mp.mpf(1) / 4) / (7 * mp.sqrt(2) * mp.pi ** mp.mpf(1.5) * mp.gamma(mp.mpf(3) / 4))
def integrand_x(lam):
    w = lambda x: mp.sqrt(1 + lam / x ** 4) / 2
    return mp.quad(lambda x: lam / (x ** 4 + lam) * w(x) / (1 + 2 * w(x)) ** 2 / x ** 2, [0, lam ** 0.25, 1, mp.inf])
A1q = -12 / mp.pi ** 2 * mp.quad(lambda lam: lam * integrand_x(lam), [0, 1])
R.close("main Eq. (53)", "A_1 from the (lambda, x) double integral equals the gamma-function formula", A1q, A1c, mp.mpf('1e-20'))
R.close("SI Eq. (S91)", "-(12/(7 sqrt2 pi^2)) B(1/2, 5/4) equals the gamma-function formula", -12 / (7 * mp.sqrt(2) * mp.pi ** 2) * mp.beta(mp.mpf(1) / 2, mp.mpf(5) / 4), A1c, mp.mpf('1e-28'))
R.digits("main Eq. (53)", "A_1", A1c, "-0.2146937609")
qs_coef = (3 * mp.pi / (16 * aL)) ** (mp.mpf(1) / 4)
R.digits("main Sec. 3.3", "q_*^-1 r_s^(1/4) = [3 pi/(16 alpha_L)]^(1/4)", qs_coef, "1.0311348")
R.digits("main Table 3", "e_kite: A_1/q_* coefficient of r_s^-1/4", A1c * qs_coef, "-0.2213782293")
R.digits("main Table 3", "F: (7/8) A_1/q_* coefficient", mp.mpf(7) / 8 * A1c * qs_coef, "-0.1937059506")
R.close("main Sec. 3.3, SI Sec. S12", "Mellin multiplier 2/(p+2) at p = -1/4 is 8/7", 2 / (mp.mpf(-1) / 4 + 2), mp.mpf(8) / 7, mp.mpf('1e-28'))
R.digits("main Sec. 3.3", "A_1/q_* at r_s = 1e8, in mRy", 1000 * A1c * qs_coef * mp.mpf(10) ** -2, "-2.214")
qstar = lambda rs: (16 * aL * rs / (3 * mp.pi)) ** (mp.mpf(1) / 4)
R.digits("main Sec. 3.3", "q_* at r_s = 2e5 (units of k_F)", qstar(mp.mpf(2e5)), "20")
R.digits("main Sec. 3.3", "q_* at r_s = 1e6 (units of k_F)", qstar(mp.mpf(1e6)), "31")
# stored direct values (production/raw/tableI) and BRM reference values
rows = {}
for rs, pm, pse in [(100000000, '-2.22', '0.21'), (200000, '-9.13', None), (1000000, '-6.35', None)]:
    f = os.path.join(root, 'production/raw/tableI', f"rs_{'1e+08' if rs == 100000000 else '1e+06' if rs == 1000000 else rs}_z_0.json")
    d = json.load(open(f)); rows[rs] = d
    R.digits("main Sec. 3.3", f"direct value at r_s = {rs:g} (mRy), {os.path.relpath(f, root)}", d['mean_mry'], pm)
    if pse: R.digits("main Sec. 3.3", f"its standard error at r_s = {rs:g}", d['se_mry'], pse)
brm = json.load(open(os.path.join(root, 'production/reference_brm.json')))
t1 = {r: (v, e) for r, v, e in brm['tableI']}; t2 = {r: (v, e) for r, v, e in brm['tableII']}
R.exact("main Sec. 3.3", "BRM Monte Carlo values at r_s = 2e5 and 1e6 quoted in the text (mRy)", (t1[200000.0][0], t1[1000000.0][0]), (-7.13, -3.89))
R.exact("main Sec. 3.3", "BRM value at r_s = 0.1 quoted in the text, 43.53(18) mRy", t2[0.1], (43.53, 0.18))
raise SystemExit(R.write())
