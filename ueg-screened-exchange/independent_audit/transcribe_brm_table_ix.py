"""Record the zeta = 0 column of BRM Table IX from the text layer of the paper.

    python transcribe_brm_table_ix.py <text layer of page 26 of the BRM scan>

The page has two blocks. In the upper block each row reads
`r_s  E(zeta=0) (sd)  r_s'  E(zeta=0.2) ...`; the lower block repeats large r_s
with the zeta = 0.2 ... 1.0 columns only. A zeta = 0 row is therefore one whose
`value (sd)` is followed by an unsigned r_s' of the right-hand table. The result is
written to brm_table_ix_zeta0.json and compared with production/reference_brm.json.
The paper itself is not distributed with this archive.
"""
import json, os, re, sys
HERE = os.path.dirname(os.path.abspath(__file__))
NUM = r'[0-9]+(?:\.[0-9]+)?'
ROW = re.compile(rf'^\s*({NUM})\s+([−-]?{NUM})\s+\((\d+)\)\s+({NUM})\s+[−-]?{NUM}\s+\(')

def main(path):
    lines = open(path, encoding='utf-8', errors='replace').read().splitlines()
    caption = ' '.join(l.strip() for l in lines[1:4])
    rows, raw = [], []
    for l in lines:
        m = ROW.match(l)
        if m:
            val = m.group(2).replace('−', '-')
            dec = len(val.split('.')[1]) if '.' in val else 0
            rows.append([float(m.group(1)), float(val), int(m.group(3)) / 10 ** dec]); raw.append(l.strip())
    rec = dict(source='S. Benites, J. Rosado and E. Manousakis, Phys. Rev. B 110, 195151 (2024), Table IX, zeta = 0 column',
               method='Text layer of the publisher page image (page 26 of the supplied scan); rows of the upper block only. '
                      'The parenthesized standard deviation applies to the last printed digits. Extracted by this script, not typed.',
               caption=caption, rows=rows, raw_lines=raw)
    with open(os.path.join(HERE, 'brm_table_ix_zeta0.json'), 'w') as f:
        json.dump(rec, f, indent=1, ensure_ascii=False); f.write('\n')
    ref = json.load(open(os.path.join(HERE, '..', 'production', 'reference_brm.json')))
    tab = {r[0]: (r[1], r[2]) for r in rows}
    used = [r for key in ('tableI', 'tableII') for r in ref[key] if r[1] is not None]
    bad = [r for r in used if r[0] not in tab or abs(tab[r[0]][0] - r[1]) > 1e-9 or abs(tab[r[0]][1] - r[2]) > 1e-9]
    print(f'{len(rows)} zeta=0 rows transcribed; {len(used)} archive values checked; mismatches: {bad}')
    return 1 if bad else 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1]))
