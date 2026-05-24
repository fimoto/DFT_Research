"""Verify the signed spin-interpolation bias implied by the kite expansion.

The analytic argument checks exact algebra and the ingredients of an
elementary sign proof on the whole interval 0 < |zeta| < 1.  The separate
70-digit numerical checks are supplemental consistency checks, not a proof
of that sign or of the underlying many-body energy expansion.

All expansions below are at fixed polarization.  No energy remainder is
differentiated with respect to polarization, and no spin susceptibility is
inferred.  Running this file writes spin_bias_algebra_validation.txt beside
it; importing and calling run_checks() only prints the validation results.
"""

from contextlib import redirect_stdout
from io import StringIO
import json
from pathlib import Path

import mpmath as mp
import sympy as sp


HERE = Path(__file__).resolve().parent
LABEL_NAMES = {
    "eq:spin-endpoint-interpolation": "endpoint-calibrated interpolation",
    "eq:spin-interpolation-bias": "signed spin-interpolation bias",
    "eq:spin-bias-derivative": "derivative factorization and sign proof",
}


def run_checks():
    """Print exact and supplemental numerical checks; return their counts."""
    map_path = HERE / "equation_map.json"
    labels = json.loads(map_path.read_text())["labels"] if map_path.exists() else {}
    counts = {"EXACT": 0, "NUM": 0}

    def passed(kind, label, description):
        number = labels.get(label)
        reference = f"Eq.({number})" if number else LABEL_NAMES[label]
        print(f"{kind} PASS {reference} [{label}] {description}")
        counts[kind] += 1

    endpoint_label = "eq:spin-endpoint-interpolation"
    bias_label = "eq:spin-interpolation-bias"
    derivative_label = "eq:spin-bias-derivative"

    # Keep all energy coefficients independent to check the interpolation
    # construction before imposing the particular spin moments.
    t = sp.Symbol("t", real=True)  # t = 1 / q_*; never a spin variable.
    a1, a2, d4, d5 = sp.symbols("a1 a2 d4 d5", nonzero=True)
    m4, m5, a3z, a30, a31 = sp.symbols("m4 m5 a3z a30 a31")
    fx = (m4 - 1) / d4
    energy = a1*m4*t + a2*m5*t**2 + a3z*t**3
    energy0 = a1*t + a2*t**2 + a30*t**3
    energy1 = a1*(1+d4)*t + a2*(1+d5)*t**2 + a31*t**3
    interpolated = energy0 + fx*(energy1-energy0)
    assert sp.simplify(interpolated.subs(m4, 1)-energy0) == 0
    assert sp.simplify(interpolated.subs(m4, 1+d4)-energy1) == 0
    passed("EXACT", endpoint_label,
           "epsilon_int = epsilon_0 + f_x*(epsilon_1-epsilon_0) matches both endpoints.")

    bracket = m5 - 1 - d5*fx
    residual = sp.expand(energy-interpolated)
    assert sp.simplify(residual.coeff(t, 1)) == 0
    assert sp.simplify(residual.coeff(t, 2)-a2*bracket) == 0
    assert sp.simplify(residual.coeff(t, 3)
                       -(a3z-a30-fx*(a31-a30))) == 0
    passed("EXACT", bias_label,
           "the 1/q_* term cancels and the 1/q_*^2 residual is A_2*B; the next coefficient is retained independently.")

    # A single density-independent shape f would have to match the first
    # two independent powers of t.  Their required shapes disagree by
    # B/d5, strictly positive by the interval proof below.
    f = sp.Symbol("f")
    trial_residual = sp.expand(energy-energy0-f*(energy1-energy0))
    f_first = sp.solve(trial_residual.coeff(t, 1), f)[0]
    f_second = sp.solve(trial_residual.coeff(t, 2), f)[0]
    assert sp.simplify(f_first-fx) == 0
    assert sp.simplify(f_second-f_first-bracket/d5) == 0
    passed("EXACT", bias_label,
           "matching the first two powers requires shapes differing by B/d_5; a density-independent shape cannot match both at interior polarization.")

    c = sp.root(2, 3)
    z = sp.Symbol("z", real=True)

    def moment(n):
        return ((1+z)**sp.Rational(n, 3)
                +(1-z)**sp.Rational(n, 3))/2

    b_z = moment(5) - (1+c)*moment(4) + c
    original_b = moment(5)-1-(c**2-1)*(moment(4)-1)/(c-1)
    assert sp.simplify(b_z-original_b) == 0
    assert sp.simplify((c**2-1)/(c-1)-(1+c)) == 0
    assert sp.simplify(b_z.subs(z, -z)-b_z) == 0
    for endpoint in (-1, 0, 1):
        assert sp.simplify(b_z.subs(z, endpoint)) == 0
    passed("EXACT", bias_label,
           "B is even, B(0)=B(+1)=B(-1)=0, and d_5/d_4=1+2^(1/3).")

    # On 0 < z < 1, x=(1+z)^(1/3), y=(1-z)^(1/3) are positive.
    # Use dx/dz=1/(3x^2), dy/dz=-1/(3y^2) to avoid fractional-power
    # simplifications that would require complex branch conventions.
    x, y, cc = sp.symbols("x y c", positive=True)
    b_xy = (x**5+y**5)/2 - (1+cc)*(x**4+y**4)/2 + cc
    derivative = sp.diff(b_xy, x)/(3*x**2) - sp.diff(b_xy, y)/(3*y**2)
    derivative_factored = (x-y)*(5*(x+y)-4*(1+cc))/6
    assert sp.factor(derivative-derivative_factored) == 0
    passed("EXACT", derivative_label,
           "B'=(x-y)*[5*(x+y)-4*(1+c)]/6 by the exact chain rule.")

    # The positive real cube root c satisfies 1<c<3/2, proved by cubing
    # these rational bounds.  This does not use floating-point estimates.
    assert sp.Integer(1)**3 < 2 < sp.Rational(3, 2)**3
    assert sp.simplify(c**3-2) == 0
    assert bool(c > 1) and bool(c < sp.Rational(3, 2))
    # x^3-y^3=2z>0 and x^2+xy+y^2>0 imply x-y>0.
    assert sp.factor(x**3-y**3) == (x-y)*(x**2+x*y+y**2)
    assert (x**2+x*y+y**2).is_positive
    # Once x>y, set x=y+h with h>0 to establish strict monotonicity.
    h = sp.Symbol("h", positive=True)
    sum_derivative = 1/(3*x**2)-1/(3*y**2)
    decreasing = sp.factor(sum_derivative.subs(x, y+h))
    assert decreasing == -h*(h+2*y)/(3*y**2*(h+y)**2)
    assert decreasing.is_negative
    # x+y has endpoint values 2 and c; the linear bracket is positive
    # initially and negative finally.  Strict monotonicity gives one zero.
    start_bracket = 10-4*(1+c)
    end_bracket = 5*c-4*(1+c)
    assert bool(start_bracket > 0)
    assert bool(end_bracket < 0)
    passed("EXACT", derivative_label,
           "analytic interval proof: 1<c<3/2, x>y, and x+y strictly decreases from 2 to c; B' changes sign exactly once (+ to -). With zero endpoint values, B>0 for 0<|zeta|<1.")

    aa, bb, dd, positive_b = sp.symbols("aa bb dd positive_b", positive=True)
    signed_h = bb/(-aa*dd)*positive_b
    signed_bias = bb*positive_b
    assert signed_h.is_negative and signed_bias.is_positive
    passed("EXACT", bias_label,
           "A_1<0<A_2 and d_4>0 imply H<0 and positive leading energy bias at each fixed interior polarization; no remainder is differentiated.")

    with mp.workdps(70):
        cm = mp.root(2, 3)
        a1m = (-4/(7*mp.sqrt(2)*mp.pi**mp.mpf("1.5"))
               * mp.gamma(mp.mpf(".25"))/mp.gamma(mp.mpf(".75")))
        a2m = mp.mpf(8)/15
        pref = (-mp.sqrt(2*mp.pi)/mp.pi**2
                * mp.gamma(mp.mpf(".75"))/mp.gamma(mp.mpf(".25")))

        def mom(n, zz):
            return ((1+zz)**(mp.mpf(n)/3)+(1-zz)**(mp.mpf(n)/3))/2

        def third(zz):
            return pref*(mp.mpf(1856)/675*mom(6, zz)
                         + mp.mpf(264)/125*mom(4, zz)*mom(5, zz))

        def energy3(qq, zz):
            return a1m*mom(4, zz)/qq+a2m*mom(5, zz)/qq**2+third(zz)/qq**3

        for index in range(1, 200):
            zz = mp.mpf(index)/200
            xm, ym = mp.root(1+zz, 3), mp.root(1-zz, 3)
            bm = mom(5, zz)-(1+cm)*mom(4, zz)+cm
            fxm = (mom(4, zz)-1)/(cm-1)
            assert bm > 0
            assert a2m/(a1m*(cm-1))*bm < 0
            expanded_derivative = (mp.mpf(5)/6*((1+zz)**(mp.mpf(2)/3)
                                   -(1-zz)**(mp.mpf(2)/3))
                                   - mp.mpf(2)/3*(1+cm)*(xm-ym))
            factored_derivative = (xm-ym)*(5*(xm+ym)-4*(1+cm))/6
            assert abs(expanded_derivative-factored_derivative) < mp.mpf("1e-65")
            assert abs(bm-(mom(5, zz)-1-(cm**2-1)*fxm)) < mp.mpf("1e-65")
        passed("NUM", derivative_label,
               "supplemental 70-digit evaluations confirm the identities and signs at 199 interior polarizations; the sign proof is analytic above.")

        ntests = 0
        for zz in map(mp.mpf, ("-.99", "-.5", "-.01", ".01", ".25", ".5", ".9", ".99")):
            fxm = (mom(4, zz)-1)/(cm-1)
            bm = mom(5, zz)-1-(cm**2-1)*fxm
            third_residual = third(zz)-third(0)-fxm*(third(1)-third(0))
            for qq in map(mp.mpf, ("100", "1000", "10000")):
                interpolation = energy3(qq, 0)+fxm*(energy3(qq, 1)-energy3(qq, 0))
                direct_bias = energy3(qq, zz)-interpolation
                assembled_bias = a2m*bm/qq**2+third_residual/qq**3
                assert abs(direct_bias-assembled_bias) < mp.mpf("1e-65")
                assert direct_bias > 0
                assert abs((direct_bias-a2m*bm/qq**2)*qq**3
                           - third_residual) < mp.mpf("1e-55")
                ntests += 1
        passed("NUM", bias_label,
               f"supplemental 70-digit direct assembly confirms {ntests} three-term energy residuals and their 1/q_*^3 corrections; this is not a full-integral energy test.")

    print(f"SUMMARY: {counts['EXACT']} EXACT PASS; {counts['NUM']} NUM PASS. "
          "Scope: algebraic consequences of the fixed-polarization kite expansion.")
    return counts


if __name__ == "__main__":
    output = StringIO()
    with redirect_stdout(output):
        run_checks()
    report = output.getvalue()
    (HERE / "spin_bias_algebra_validation.txt").write_text(report)
    print(report, end="")
