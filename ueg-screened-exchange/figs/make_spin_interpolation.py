"""Plot spin interpolation from the three analytic kite coefficients.

Run from any directory: python figs/make_spin_interpolation.py
Outputs are a deterministic PDF and CSV. No integral or FRG calculation
is performed. The coefficient derivation is in the PRB manuscript.
"""
from pathlib import Path
from datetime import datetime, timezone
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def moment(n, zeta):
    z = np.asarray(zeta, dtype=float)
    return ((1 + z)**(n / 3) + (1 - z)**(n / 3)) / 2


def coefficients():
    ratio = math.gamma(0.25) / math.gamma(0.75)
    a1 = -4 * ratio / (7 * math.sqrt(2) * math.pi**1.5)
    a2 = 8 / 15
    pref3 = -math.sqrt(2 * math.pi) / (math.pi**2 * ratio)
    return a1, a2, pref3


def qstar(rs):
    alpha = (4 / (9 * math.pi))**(1 / 3)
    a = 4 * alpha / math.pi
    return (4 * a * np.asarray(rs) / 3)**0.25


def fx(zeta):
    return (moment(4, zeta) - 1) / (2**(1 / 3) - 1)


def correction(zeta):
    a1, a2, _ = coefficients()
    d4, d5 = 2**(1 / 3) - 1, 2**(2 / 3) - 1
    return a2 / (a1 * d4) * (moment(5, zeta) - 1 - d5 * fx(zeta))


def scaled_deviation(rs, zeta):
    """q_* (f_kite^(3) - f_x), evaluated without subtracting nearby ratios."""
    q = qstar(rs)
    a1, a2, pref3 = coefficients()
    d4, d5 = 2**(1 / 3) - 1, 2**(2 / 3) - 1
    m4, m5, m6 = (moment(n, zeta) for n in (4, 5, 6))
    a3z = pref3 * (1856 * m6 / 675 + 264 * m4 * m5 / 125)
    a30 = pref3 * 16408 / 3375
    # At full polarization a3(1)=2*a3(0), hence its difference is a3(0).
    num = a2 * (m5 - 1 - d5 * fx(zeta))
    num += ((a3z - a30) - a30 * fx(zeta)) / q
    den = a1 * d4 + a2 * d5 / q + a30 / q**2
    return num / den


def main():
    zeta = np.linspace(0, 1, 501)
    densities = (1e4, 1e6, 1e8)
    curves = [scaled_deviation(rs, zeta) for rs in densities]
    data_dir = HERE / "data"
    data_dir.mkdir(exist_ok=True)
    np.savetxt(data_dir / "fig_spin_interpolation.csv",
               np.column_stack([zeta, fx(zeta), correction(zeta), *curves]),
               delimiter=",", fmt="%.17g",
               header="zeta,fx,H,qstar_f3_minus_fx_rs1e4,qstar_f3_minus_fx_rs1e6,qstar_f3_minus_fx_rs1e8",
               comments="")
    plt.rcParams.update({"font.family": "serif", "font.size": 8.5,
                         "mathtext.fontset": "dejavuserif", "pdf.fonttype": 42,
                         "axes.linewidth": 0.7, "xtick.direction": "in",
                         "ytick.direction": "in"})
    fig, ax = plt.subplots(figsize=(3.4, 2.8), layout="constrained")
    for vals, color, exponent in zip(curves, ("#0072B2", "#D55E00", "#009E73"), (4, 6, 8)):
        ax.plot(zeta, vals, color=color, lw=1.5, label=rf"$r_s=10^{exponent}$")
    ax.plot(zeta, correction(zeta), color="black", ls="--", lw=1.4,
            label=r"$H(\zeta)$, limiting curve")
    ax.set(xlim=(0, 1), xlabel=r"Spin polarization $\zeta$",
           ylabel=r"$q_*[f_{\rm kite}^{(3)}-f_x]$")
    ax.axhline(0, color="0.75", lw=0.5, zorder=0)
    ax.legend(loc="lower left", frameon=False, fontsize=7.6, handlelength=2.3)
    ax.tick_params(top=True, right=True)
    stamp = datetime(2026, 9, 27, tzinfo=timezone.utc)
    fig.savefig(HERE / "fig_spin_interpolation.pdf",
                metadata={"Creator": "Analytic kite spin interpolation",
                          "CreationDate": stamp, "ModDate": stamp})
    plt.close(fig)
    print("Wrote analytic figure and CSV; no new many-body integral evaluations.")


if __name__ == "__main__":
    main()
