#!/usr/bin/env python3
"""
Check that the raw-data pipeline (make_tables.py) reproduces the values recorded
in results/reference/paper_example_values.json, within the printed rounding.

paper_example_values.json is a frozen transcription of the printed per-nucleus and
shared multi-nucleus tables (paper Tables II and III). It must be kept in sync with the
manuscript .tex by hand: when a table number in the .tex changes, update this
file to match, and this validator then confirms the released raw data still
reproduces it. (The .tex itself is distributed separately from this archive, so
the validator does not parse it directly.)

Exit status is non-zero if any value disagrees beyond half a ULP of its printed
precision, so this doubles as a CI gate.
"""
import json, os, sys, argparse

# printed precision per field (decimal places shown in the manuscript)
PREC1 = dict(E_KS=3, E64=3, dE_abs=3, R_EL=4, rho_c=4, rho_c_KS=4, rrms=3, rrms_KS=3, dNLF=4)
PREC2 = dict(E64=2, dE_ref=1, E_grad=2, rho_c=4, rho_c_KS=4, rrms=3, rrms_KS=3, dNLF=4)


def check(gen, exp, prec, label, errors):
    for tag in exp:
        if tag not in gen:
            errors.append('%s: %s missing from generated table' % (label, tag))
            continue
        for key, p in prec.items():
            if key not in exp[tag]:
                continue
            tol = 0.5 * 10**(-p) + 1e-12
            gv, ev = gen[tag].get(key), exp[tag][key]
            if gv is None:
                errors.append('%s %s.%s missing from generated' % (label, tag, key))
                continue
            if abs(gv - ev) > tol:
                errors.append('%s %s.%s: generated %.*f vs manuscript %.*f (|d|=%.2e > %.2e)'
                              % (label, tag, key, p, gv, p, ev, abs(gv - ev), tol))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--tables', required=True, help='directory with generated table1.json/table2.json')
    ap.add_argument('--expected', required=True, help='paper_example_values.json')
    args = ap.parse_args()

    t1 = json.load(open(os.path.join(args.tables, 'table1.json')))
    t2 = json.load(open(os.path.join(args.tables, 'table2.json')))
    exp = json.load(open(args.expected))

    errors = []
    check(t1, exp['table1'], PREC1, 'Table II', errors)
    check(t2, exp['table2'], PREC2, 'Table III', errors)

    if errors:
        print('VALIDATION FAILED (%d discrepancies):' % len(errors))
        for e in errors:
            print('  -', e)
        sys.exit(1)
    print('VALIDATION PASSED: raw-data pipeline reproduces the recorded manuscript values (Tables II & III) within printed rounding.')


if __name__ == '__main__':
    main()
