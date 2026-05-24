"""Write the plotted series as plain ASCII, for gnuplot or any other tool.

Every figure script calls `write` once per data block.  The files land in
`figs/data/` as whitespace-separated columns with a commented header, which is
what gnuplot's `plot "file" using 1:2` expects; `#` lines are ignored by
gnuplot, awk and numpy.loadtxt alike.  Missing entries are written as `NaN`,
which gnuplot skips.

Analytic curves that the figures draw as smooth lines are written on the same
grid as they are plotted, so a gnuplot user can either reuse the column or
type the closed form in directly; the constants are repeated in each header.
"""
import os

import numpy as np

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def write(name, columns, header, notes=()):
    """Write `columns` (a list of (label, 1-D array)) to figs/data/<name>.

    `header` is a one-line description, `notes` extra comment lines.
    """
    os.makedirs(DATA, exist_ok=True)
    labels = [lab for lab, _ in columns]
    arrays = [np.asarray(col, dtype=float) for _, col in columns]
    n = len(arrays[0])
    if any(len(a) != n for a in arrays):
        raise ValueError("columns of unequal length in %s" % name)
    path = os.path.join(DATA, name)
    with open(path, "w") as fh:
        fh.write("# %s\n" % header)
        for line in notes:
            fh.write("# %s\n" % line)
        fh.write("# columns: " + "  ".join(
            "%d:%s" % (i + 1, lab) for i, lab in enumerate(labels)) + "\n")
        for row in zip(*arrays):
            fh.write("  ".join("% .10E" % v for v in row) + "\n")
    print("wrote figs/data/%s" % name)
    return path
