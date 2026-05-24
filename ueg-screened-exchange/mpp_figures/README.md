# Figures of the Mellin-pole paper

Scripts for the four figures of the main text of *Revisiting the $r_s$ expansion of the uniform electron gas in $d$ dimensions: $\log r_s$ terms from coincident Mellin poles*, and the PDF files they produce. The file names are those used by `\includegraphics` in the manuscript.

```sh
python -m pip install -r mpp_figures/requirements.txt
python mpp_figures/make_figures.py        # all four figures
python mpp_figures/make_figures.py 2      # Fig. 2 only
```

| Figure | Script | Input |
| --- | --- | --- |
| Fig. 1, positions of the Mellin poles | `fig1_poles.py` | `../mpp_checks/figure_data/fig1_poles.csv`, `fig1_coincidences.csv` |
| Fig. 2, $R_2(z)$ and $h_2(z)$ | `fig2_R2.py` | `../mpp_checks/figure_data/fig2_R2.csv`, `fig2_R2_loglog.csv` |
| Fig. 3, diagrams | `fig3_exchange_response.tex` (TikZ, needs `pdflatex`) | none |
| Fig. 4, momentum regions of the kite | `fig4_ball_regions.py` | closed form $\delta_{\bm p}=\bm p\cdot\bm q+q^2/2$, tabulated also in `../mpp_checks/figure_data/fig4_delta.csv` |

The CSV files are written by `../mpp_checks/fig_data.py`.

## Sizes and fonts

Each figure is drawn at the width at which the manuscript includes it, so that LaTeX does not rescale it and the font sizes set here are the printed sizes. The text width of `iopjournal.cls` is 153 mm.

| Figure | `\includegraphics` width | Printed width |
| --- | --- | --- |
| Fig. 1 | `0.56\textwidth` | 85.7 mm |
| Fig. 2 | `0.95\textwidth` | 145.4 mm |
| Fig. 3 | `0.9\textwidth` | 137.7 mm |
| Fig. 4 | `0.5\textwidth` | 76.5 mm |

To change the layout or a font size:

* **Font sizes and line widths** of Figs. 1, 2 and 4: edit `FONT` and `LW` in `figstyle.py`. The defaults are 9 pt for axis labels and panel labels, 8 pt for tick labels and for the table of Fig. 4(b), and 7.5 pt for legends and labels inside the axes. The captions of the journal class are 8 pt.
* **Width of a figure**: change `WIDTH_FRACTION` in `figstyle.py` and the `\includegraphics` line of the manuscript together. The height is set by the `aspect` argument of `figstyle.figsize` in each script.
* **Placement of panels, legends and labels**: the axes rectangles (`fig.add_axes([left, bottom, width, height])`, in fractions of the figure) and the text positions (data coordinates) are set in each script.
* **Fig. 3**: `\titlefont` and `\labelfont` at the top of `fig3_exchange_response.tex` set the font sizes, `\rx` and `\ry` the size of the bubble, and `\dx` the spacing of the columns.

The figures are saved without a tight bounding box, so their size is exactly the size set in the script.

With the pinned versions the scripts reproduce the stored PDF files byte for byte (the hashes are in `../release_manifest.json`). The PDF of Fig. 3 is the same byte for byte only with the same TeX Live version (2023 here). Other versions give the same drawing.
