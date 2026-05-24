"""Main-text Fig. 4: momentum regions of the kite loop for q > 2.

(a) Cross section through the axis of q for q = 3.2.  The color is the
excitation energy delta_p = p.q + q^2/2 in units of q^2/2, evaluated from this
closed form (the same quantity is tabulated in
../mpp_checks/figure_data/fig4_delta.csv).  (b) The four regions of the two
loop momenta.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle
import figstyle as st

NAME = 'fig4_ball_regions'
Q = 3.2                      # |q| in units of k_F
P = r'\boldsymbol{p}'
K = r'\boldsymbol{k}'
B = r'\mathbb{B}'
BQ = r'\mathbb{B}_{-\boldsymbol{q}}'


def delta_map(n=801):
    """2 delta_p / q^2 on the plane through the q axis (q along +x), inside the
    Fermi ball B (center 0) and its translate B_{-q} (center -q)."""
    x = np.linspace(-Q - 1, 1, int(n * (Q + 2) / 2))
    y = np.linspace(-1, 1, n)
    X, Y = np.meshgrid(x, y)
    val = 2 * X / Q + 1                      # (p.q + q^2/2)/(q^2/2)
    inside = (X ** 2 + Y ** 2 <= 1) | ((X + Q) ** 2 + Y ** 2 <= 1)
    return np.ma.masked_where(~inside, val), (x[0], x[-1], y[0], y[-1])


def panel_a(fig, rect, cbar_rect):
    ax = fig.add_axes(rect)
    val, ext = delta_map()
    vmax = 1 + 2 / Q
    im = ax.imshow(val, extent=ext, origin='lower', cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                   interpolation='bilinear', zorder=1)
    for c in (0.0, -Q):
        ax.add_patch(Circle((c, 0), 1, fill=False, ec='black', lw=st.LW['thin'], zorder=3))
        ax.plot(c, 0, 'o', ms=2.2, color='black', zorder=4)
    ax.plot([-Q - 1.25, 1.25], [0, 0], color='0.55', lw=0.5, ls=(0, (1.5, 1)), zorder=2)
    ax.annotate('', xy=(0, 0), xytext=(-Q, 0), zorder=5,
                arrowprops=dict(arrowstyle='-|>', color='black', lw=st.LW['thin'],
                                mutation_scale=7, shrinkA=0, shrinkB=1.5))
    ax.text(-Q / 2, 0.14, r'$\boldsymbol{q}$', ha='center', va='bottom', fontsize=st.FONT['label'], zorder=5)
    fa = st.FONT['annotation']
    ax.text(-Q, 1.12, rf'$\delta_{P}\simeq{{-q^2/2}}$', ha='center', va='bottom', fontsize=fa)
    ax.text(0, 1.12, rf'$\delta_{P}\simeq{{+q^2/2}}$', ha='center', va='bottom', fontsize=fa)
    ax.text(-Q, -1.12, rf'${BQ}$: $\Delta_{P}={{-1}}$', ha='center', va='top', fontsize=fa)
    ax.text(0, -1.12, rf'${B}$: $\Delta_{P}={{+1}}$', ha='center', va='top', fontsize=fa)
    ax.set_xlim(-Q - 1.3, 1.3)
    ax.set_ylim(-1.15, 1.15)
    ax.set_aspect('equal')
    ax.axis('off')
    cax = fig.add_axes(cbar_rect)
    cb = fig.colorbar(im, cax=cax, ticks=[-1.5, 0, 1.5])
    cb.ax.tick_params(labelsize=st.FONT['tick'] - 1, length=2, width=0.6)
    cb.outline.set_linewidth(0.5)
    cax.set_title(rf'$2\delta_{P}/q^2$', fontsize=fa, pad=4)


def panel_b(fig, rect):
    ax = fig.add_axes(rect)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    x0, x1, x2 = 0.24, 0.62, 1.0          # column edges
    y0, y1, y2 = 0.0, 0.43, 0.86          # row edges (bottom, middle, top)
    cells = {
        (0, 0): ('+1', '>', r'\leq2'), (0, 1): ('-1', '<', r'\simeq q'),
        (1, 0): ('-1', '<', r'\simeq q'), (1, 1): ('+1', '>', r'\leq2'),
    }
    cols = [(x0, x1), (x1, x2)]
    rows = [(y1, y2), (y0, y1)]            # row 0 (p in B) on top
    ft = st.FONT['table']
    for (i, j), (dd, sg, pk) in cells.items():
        (xa, xb), (ya, yb) = cols[j], rows[i]
        ax.add_patch(Rectangle((xa, ya), xb - xa, yb - ya, fc='0.88' if i == j else 'white',
                               ec='black', lw=st.LW['thin']))
        xc, yc, dy = (xa + xb) / 2, (ya + yb) / 2, 0.125
        ax.text(xc, yc + dy, rf'$\Delta_{P}\Delta_{K}={{{dd}}}$', ha='center', va='center', fontsize=ft)
        ax.text(xc, yc, rf'$\delta_{P}\delta_{K}{sg}0$', ha='center', va='center', fontsize=ft)
        ax.text(xc, yc - dy, rf'$|{P}-{K}|{pk}$', ha='center', va='center', fontsize=ft)
    for j, lab in enumerate((B, BQ)):
        ax.text(sum(cols[j]) / 2, y2 + 0.07, rf'${K}\in{lab}$', ha='center', va='center', fontsize=ft + 0.5)
    for i, lab in enumerate((B, BQ)):
        ax.text(x0 - 0.03, sum(rows[i]) / 2, rf'${P}\in{lab}$', ha='right', va='center', fontsize=ft + 0.5)


def main():
    st.apply()
    fig = plt.figure(figsize=st.figsize(NAME, 0.97))
    panel_a(fig, [0.03, 0.555, 0.76, 0.40], [0.865, 0.62, 0.028, 0.25])
    panel_b(fig, [0.0, 0.015, 0.985, 0.47])
    fig.text(0.015, 0.985, '(a)', ha='left', va='top', fontsize=st.FONT['panel'])
    fig.text(0.015, 0.50, '(b)', ha='left', va='top', fontsize=st.FONT['panel'])
    return st.save(fig, NAME)


if __name__ == '__main__':
    print(main())
