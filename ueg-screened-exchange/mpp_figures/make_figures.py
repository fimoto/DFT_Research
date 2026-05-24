"""Draw the four figures of the main text of the Mellin-pole paper (MPP).

    python mpp_figures/make_figures.py            # all figures
    python mpp_figures/make_figures.py 1 4        # Figs. 1 and 4 only

Figs. 1, 2 and 4 use matplotlib.  Fig. 3 is TikZ and needs pdflatex.
The PDF files are written next to this script under the names used by
\\includegraphics in MPP_main.tex.  Sizes and fonts are set in figstyle.py.
"""
import os, subprocess, sys, glob

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def fig3():
    tex = 'fig3_exchange_response.tex'
    # fixed dates so that the PDF file is the same on every run
    env = dict(os.environ, SOURCE_DATE_EPOCH='1767225600', FORCE_SOURCE_DATE='1')
    subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error', tex],
                   cwd=HERE, check=True, stdout=subprocess.DEVNULL, env=env)
    for ext in ('aux', 'log'):
        for f in glob.glob(os.path.join(HERE, '*.' + ext)):
            os.remove(f)
    return os.path.join(HERE, 'fig3_exchange_response.pdf')


def main(which):
    import fig1_poles, fig2_R2, fig4_ball_regions
    jobs = {'1': fig1_poles.main, '2': fig2_R2.main, '3': fig3, '4': fig4_ball_regions.main}
    for k in which or sorted(jobs):
        print(jobs[k]())


if __name__ == '__main__':
    main(sys.argv[1:])
