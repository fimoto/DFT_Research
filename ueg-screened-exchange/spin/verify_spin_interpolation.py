"""Algebraic checks of the added spin interpolation and its plotting data.

This checks consequences of the existing three energy coefficients, not
their underlying many-body derivation and not the accuracy of FRG-DFT.
"""
from pathlib import Path
import importlib.util
import json
import numpy as np
import sympy as s
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("spin_plot", ROOT / "figs/make_spin_interpolation.py")
plot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plot)


def run_checks():
    labels = json.loads((ROOT / "spin/equation_map.json").read_text())["labels"]
    expansion_number = labels["eq:spin-interpolation-expansion"]
    correction_number = labels["eq:spin-interpolation-correction"]
    # Derive the correction by differentiating the quotient with respect to
    # t=1/q_* at t=0, keeping all energy coefficients independent.
    t = s.Symbol("t", real=True)
    a1, a2, d4, d5 = s.symbols("a1 a2 d4 d5", nonzero=True)
    n4, n5, n3, d3 = s.symbols("n4 n5 n3 d3")
    ratio = (a1*n4 + a2*n5*t + n3*t*t)/(a1*d4 + a2*d5*t + d3*t*t)
    leading = n4/d4
    corr = a2/(a1*d4)*(n5-d5*leading)
    assert s.simplify(ratio.subs(t, 0)-leading) == 0
    assert s.simplify(s.diff(ratio, t).subs(t, 0)-corr) == 0
    h2 = s.simplify(s.diff(ratio, t, 2).subs(t, 0)/2)
    # Substituting the expansion back into numerator-denominator identity
    # provides a second assembly check through the next order.
    residual = s.expand((leading+corr*t+h2*t*t)*(a1*d4+a2*d5*t+d3*t*t)
                        -(a1*n4+a2*n5*t+n3*t*t))
    assert all(s.simplify(residual.coeff(t, k)) == 0 for k in range(3))
    print(f"EXACT PASS Eqs.({expansion_number}),({correction_number}) symbolic quotient gives f_x and H; the next coefficient also cancels.")

    x = s.real_root(2, 3)
    for endpoint, m4, m5 in [(0, s.Integer(1), s.Integer(1)),
                              (1, x, x*x), (-1, x, x*x)]:
        value = s.simplify((m5-1)-(x*x-1)*(m4-1)/(x-1))
        assert value == 0
    print(f"EXACT PASS Eq.({correction_number}) H(0)=H(+1)=H(-1)=0 exactly.")

    # Independent high precision evaluation of the original three-term
    # energy quotient, rather than the rearranged plotting expression.
    mp.mp.dps = 70
    a1m = -4/(7*mp.sqrt(2)*mp.pi**mp.mpf("1.5"))*mp.gamma(mp.mpf(".25"))/mp.gamma(mp.mpf(".75"))
    a2m = mp.mpf(8)/15
    pref = -mp.sqrt(2*mp.pi)/mp.pi**2*mp.gamma(mp.mpf(".75"))/mp.gamma(mp.mpf(".25"))
    def mom(n, z):
        return ((1+z)**(mp.mpf(n)/3)+(1-z)**(mp.mpf(n)/3))/2
    def energy(q, z):
        return a1m*mom(4,z)/q+a2m*mom(5,z)/q**2+pref*(mp.mpf(1856)/675*mom(6,z)+mp.mpf(264)/125*mom(4,z)*mom(5,z))/q**3
    for rs in (mp.mpf("1e4"), mp.mpf("1e6"), mp.mpf("1e8")):
        q = (mp.mpf(16)/(3*mp.pi)*(mp.mpf(4)/(9*mp.pi))**(mp.mpf(1)/3)*rs)**mp.mpf(".25")
        for z in (mp.mpf(".1"), mp.mpf(".5"), mp.mpf(".9"), mp.mpf("1")):
            fx = (mom(4,z)-1)/(2**(mp.mpf(1)/3)-1)
            direct = q*((energy(q,z)-energy(q,0))/(energy(q,1)-energy(q,0))-fx)
            numerical = plot.scaled_deviation(float(rs), float(z))
            assert abs(float(direct)-float(numerical)) < 2e-12
    print(f"NUM PASS Eq.({expansion_number}) all plotted formulas match the direct 70-digit energy quotient.")

    for z in (.1, .5, .9):
        # For the three-term ratio, q_*(f^(3)-f_x)-H is O(1/q_*).
        errors = [abs(float(plot.scaled_deviation(rs,z)-plot.correction(z)))
                  for rs in (1e12, 16e12, 256e12)]
        assert 1.98 < errors[0]/errors[1] < 2.02
        assert 1.98 < errors[1]/errors[2] < 2.02
    print(f"NUM PASS Eq.({expansion_number}) subtracted ratio converges at the predicted first-correction order.")

    data = np.genfromtxt(ROOT / "figs/data/fig_spin_interpolation.csv", delimiter=",", names=True)
    assert len(data) == 501
    assert np.max(np.abs(data["fx"]-plot.fx(data["zeta"]))) < 1e-14
    for rs, field in [(1e4,"qstar_f3_minus_fx_rs1e4"),(1e6,"qstar_f3_minus_fx_rs1e6"),(1e8,"qstar_f3_minus_fx_rs1e8")]:
        assert np.max(np.abs(data[field]-plot.scaled_deviation(rs,data["zeta"]))) < 1e-14
    assert np.max(np.abs(data["H"]-plot.correction(data["zeta"]))) < 1e-14
    print("NUM PASS spin-interpolation CSV: all 501 rows reproduce the analytic plotting functions.")


if __name__ == "__main__":
    run_checks()
