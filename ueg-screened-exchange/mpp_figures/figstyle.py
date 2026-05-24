"""Common style for the figures of the Mellin-pole paper (MPP).

All sizes are the sizes on the printed page.  Each figure is drawn at the
width at which the manuscript includes it, so that LaTeX does not rescale it:

    \\includegraphics[width=WIDTH\\textwidth]{figN_....pdf}

with WIDTH equal to the entry of WIDTH_FRACTION below.  To change a font size
or a figure width, edit this file (and, for a width, the \\includegraphics line
of the manuscript), then run make_figures.py.
"""
import os
import matplotlib as mpl

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'mpp_checks', 'figure_data')

# Text width of iopjournal.cls (2024/01/31): \textwidth = 153 mm.
TEXTWIDTH_IN = 153 / 25.4

# Fraction of \textwidth used by \includegraphics in MPP_main.tex.
WIDTH_FRACTION = {
    'fig1_poles': 0.56,
    'fig2_R2': 0.95,
    'fig3_exchange_response': 0.90,
    'fig4_ball_regions': 0.50,
}

# Font sizes in pt on the printed page.  The captions of iopjournal.cls are 8 pt
# and the body text is 12 pt.
FONT = {
    'label': 9.0,       # axis labels
    'tick': 8.0,        # tick labels
    'panel': 9.0,       # panel labels (a), (b)
    'legend': 7.5,      # legend entries
    'annotation': 7.5,  # labels inside the axes
    'table': 8.0,       # text of the table in Fig. 4(b)
}

# Line widths in pt.
LW = {'axis': 0.8, 'main': 1.6, 'secondary': 1.2, 'thin': 0.8}

# Colors (matplotlib tab10 blue and orange, grays).
BLUE = '#1f77b4'
ORANGE = '#ff7f0e'
GRAY_LIGHT = '0.75'
GRAY_DARK = '0.4'


def figsize(name, aspect):
    """Width from WIDTH_FRACTION, height = aspect * width (inches)."""
    w = WIDTH_FRACTION[name] * TEXTWIDTH_IN
    return (w, aspect * w)


def apply():
    mpl.rcParams.update({
        'font.family': 'serif',
        'font.serif': ['cmr10'],
        'mathtext.fontset': 'cm',
        'axes.formatter.use_mathtext': True,
        'axes.unicode_minus': True,
        'font.size': FONT['label'],
        'axes.labelsize': FONT['label'],
        'xtick.labelsize': FONT['tick'],
        'ytick.labelsize': FONT['tick'],
        'legend.fontsize': FONT['legend'],
        'axes.linewidth': LW['axis'],
        'xtick.major.width': LW['axis'],
        'ytick.major.width': LW['axis'],
        'xtick.minor.width': 0.6,
        'ytick.minor.width': 0.6,
        'xtick.major.size': 3.5,
        'ytick.major.size': 3.5,
        'xtick.minor.size': 2.0,
        'ytick.minor.size': 2.0,
        'legend.frameon': False,
        'legend.handlelength': 2.0,
        'pdf.fonttype': 42,
        'savefig.dpi': 600,
    })


def save(fig, name):
    """Save at the exact figure size (no tight bounding box), so that the
    printed scale is 1 when included at WIDTH_FRACTION[name] * \\textwidth."""
    path = os.path.join(HERE, name + '.pdf')
    fig.savefig(path, metadata={'CreationDate': None, 'ModDate': None})
    return path
