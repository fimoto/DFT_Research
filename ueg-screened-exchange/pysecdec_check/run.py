"""Evaluate the two pySecDec packages and compare with the closed forms.

Both the finite part and every negative order are checked: for these integrands
the analytic regulator must produce a vanishing pole, so a nonzero one is a
failure of the subtraction, not a detail. The script exits nonzero if any check
fails, so it can be used in a pipeline.
"""
import platform, re, sys, time
from math import gamma, sqrt, pi, isfinite

CG = sqrt(2 * pi) * gamma(0.75) / gamma(0.25)
TARGETS = {
    'fp_cross': ('Fp int x^-4 J_C dx = -(8/5) sqrt(2pi) G(3/4)/G(1/4)', -8 / 5 * CG),
    'me_euler': ('M_E = (176/75) sqrt(2pi) G(3/4)/G(1/4)', 176 / 75 * CG),
}
POLE_TOL = 1e-8          # absolute, on the coefficient of every negative order
MAX_ERR = 1e-6           # absolute; a larger reported error makes the comparison vacuous

TERM = re.compile(
    r"\(\(\s*(%(n)s)\s*,\s*(%(n)s)\s*\)\s*\+/-\s*"
    r"\(\s*(%(n)s)\s*,\s*(%(n)s)\s*\)\)(?:\*delta\^(-?\d+))?" % dict(n=r"[-+]?(?:nan|inf|[0-9.eE+-]+)"),
    re.IGNORECASE)

def parse(series):
    """Return {order: (re, im, err_re, err_im)} for every printed order."""
    out = {}
    matches = list(TERM.finditer(series))
    for m in matches:
        order = int(m.group(5)) if m.group(5) else 0
        if order in out:
            raise SystemExit(f'order delta^{order} printed twice in the pySecDec result:\n{series}')
        out[order] = tuple(float(m.group(i)) for i in (1, 2, 3, 4))
    # The whole string must have the form  + term + term ... + O(delta^k):
    # the text between parsed terms may only be '+' and whitespace, and the
    # string must end with the order symbol. Anything else (an extra term in
    # another format, leading or trailing text) is refused rather than ignored.
    rest, pos = [], 0
    for m in matches:
        rest.append(series[pos:m.start()]); pos = m.end()
    tail = series[pos:]
    gaps_ok = all(re.fullmatch(r'\s*\+\s*', g) for g in rest)
    tail_ok = re.fullmatch(r'\s*\+\s*O\(delta(?:\^-?\d+)?\)\s*', tail) is not None
    if not out or not gaps_ok or not tail_ok:
        raise SystemExit(f'could not parse the pySecDec result completely:\n{series}')
    return out


def check(orders, exact):
    """Return the list of reasons the parsed result fails; empty means it passes."""
    if 0 not in orders:
        return ['no finite (delta^0) part in the result']
    why = []
    for order, vals in sorted(orders.items()):
        if not all(isfinite(v) for v in vals):
            why.append(f'delta^{order} term is not finite: {vals}')
        elif vals[2] < 0 or vals[3] < 0:
            why.append(f'delta^{order} term has a negative error: {vals[2:]}')
    if why:
        return why
    val, imag, err, err_im = orders[0]
    if err > MAX_ERR or err_im > MAX_ERR:
        why.append(f'reported error {max(err, err_im):.1e} exceeds {MAX_ERR:.0e}')
    if abs(val - exact) > max(5 * err, 1e-8):
        why.append(f'finite part off by {abs(val - exact):.2e} against an error of {err:.1e}')
    if abs(imag) > POLE_TOL:
        why.append(f'finite part has an imaginary part {imag:.2e}')
    for order, (re_, im_, e_, ei_) in sorted(orders.items()):
        if order < 0 and (abs(re_) > POLE_TOL or abs(im_) > POLE_TOL):
            why.append(f'delta^{order} coefficient is {re_:+.2e}{im_:+.2e}i, not zero')
    return why

def main():
    from pySecDec.integral_interface import IntegralLibrary
    import pySecDec
    print('pySecDec', pySecDec.__version__, '| Python', platform.python_version())
    failures = []
    for name, (label, exact) in TARGETS.items():
        t0 = time.time()
        lib = IntegralLibrary(f'{name}/{name}_pylink.so')
        lib.use_Qmc(verbosity=0, minn=2 ** 16, maxeval=2 ** 24,
                    epsrel=1e-7, epsabs=1e-10, transform='korobov3')
        orders = parse(lib()[2])
        why = check(orders, exact)
        val, imag, err, err_im = orders.get(0, (float('nan'),) * 4)
        ok = not why
        if not ok:
            failures.append((name, why))
        print(('OK   ' if ok else 'FAIL ') + label)
        poles = ', '.join(f'delta^{o}: {orders[o][0]:+.1e}' for o in sorted(orders) if o < 0)
        print(f'      pySecDec {val:.11f} +- {err:.1e}   exact {exact:.11f}   '
              f'|diff| {abs(val - exact):.1e}   {poles or "no negative orders"}   [{time.time() - t0:.0f}s]')
        for w in why:
            print(f'      -> {w}')
    if failures:
        print('FAILURES PRESENT')
        raise SystemExit(1)
    print('ALL OK')

if __name__ == '__main__':
    main()
