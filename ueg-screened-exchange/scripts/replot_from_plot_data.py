"""Redraw every exported figure from plot_data/ alone, as a completeness test of the export.

Reads only spec.json and the CSV files of each figure directory (no physics code, no raw
records) and writes replot.png next to them; the published PDF is rasterized to
original.png for comparison. A reader re-plotting in another program can use this script
as a worked example of how the files map onto a figure.

    python scripts/replot_from_plot_data.py [plot_data]
"""
import ast, csv, json, subprocess, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
import matplotlib.colors as mcolors
import matplotlib.cm as cm

ROOT = Path(__file__).resolve().parents[1]

def read(path):
    with open(path) as f:
        r = csv.reader(f); head = next(r)
        rows = [[float(v) if v != '' else np.nan for v in row] for row in r]
    return head, np.array(rows, float).reshape(-1, len(head))

def transform(ax, coords):
    return {'data': ax.transData, 'axes fraction': ax.transAxes,
            'x data, y axes fraction': ax.get_xaxis_transform(), 'x axes fraction, y data': ax.get_yaxis_transform()}.get(coords, ax.transData)

def replot(d):
    spec = json.loads((d / 'spec.json').read_text())
    plt.rcParams.update(spec.get('rc', {}))
    fig = plt.figure(figsize=spec['size_inches'])
    meshes = []
    for a in spec['axes']:
        ax = fig.add_axes(a['position'])
        if a['is_colorbar']:
            c = a.get('colorbar')
            if c:
                fig.colorbar(cm.ScalarMappable(norm=mcolors.Normalize(*c['clim']), cmap=c['cmap']), cax=ax, orientation=c['orientation'])
                ax.set_ylabel(c['label']) if c['orientation'] == 'vertical' else ax.set_xlabel(c['label'])
            continue
        ax.set_xscale(a['xscale'], **a.get('xscale_kwargs', {}))
        ax.set_yscale(a['yscale'], **a.get('yscale_kwargs', {}))
        for e in a['elements']:
            _, v = read(d / e['file'])
            if e['kind'] == 'line':
                if not e.get('visible', True) or len(v) == 0: continue
                ls = e['linestyle'] if e['linestyle'] not in ('None', '') else 'none'
                mk = None if e['marker'] in ('None', 'none', '') else e['marker']
                ax.plot(v[:, 0], v[:, 1], transform=transform(ax, e['coordinates']), color=e['color'], linestyle=ls, label=e['label'],
                        linewidth=e['linewidth'], marker=mk, markersize=e['markersize'], markerfacecolor=e['markerfacecolor'],
                        drawstyle=e['drawstyle'], alpha=e['alpha'], zorder=e['zorder'], clip_on=e.get('clip_on', True),
                        markevery=None if e['markevery'] == 'None' else ast.literal_eval(e['markevery']))
            elif e['kind'] == 'scatter':
                ax.scatter(v[:, 0], v[:, 1], s=v[:, 2] if not np.isnan(v[:, 2]).all() else None,
                           c=(e['facecolors'][0] if e['facecolors'] else None), edgecolors=(e['edgecolors'][0] if e['edgecolors'] else None),
                           transform=transform(ax, e['coordinates']), zorder=e['zorder'], label=e['label'], clip_on=e.get('clip_on', True))
            elif e['kind'] == 'segments':
                segs = [v[v[:, 0] == i][:, 1:] for i in np.unique(v[:, 0])] if len(v) else []
                ax.add_collection(LineCollection(segs, colors=e['colors'] or None, linewidths=e['linewidths'] or None,
                                                 transform=transform(ax, e['coordinates'])))
            elif e['kind'] == 'polygons':
                vals = e.get('values'); cmap = plt.get_cmap(e['cmap']) if e.get('cmap') else None
                polys = [v[v[:, 0] == i][:, 1:] for i in np.unique(v[:, 0])] if len(v) else []
                for i, P in enumerate(polys):
                    fc = cmap(i / max(len(polys) - 1, 1)) if (cmap is not None and vals) else (e['facecolors'][min(i, len(e['facecolors']) - 1)] if e['facecolors'] else 'none')
                    ax.add_patch(Polygon(P, closed=True, facecolor=fc, edgecolor='none', alpha=e.get('alpha'), transform=transform(ax, e['coordinates'])))
            elif e['kind'] == 'mesh_nodes':
                nodes = v.reshape(e['shape'][0], e['shape'][1], 2); meshes.append(dict(nodes=nodes))
            elif e['kind'] == 'mesh_values':
                vals = np.full(e['cells'], np.nan)
                for r in v: vals[int(r[0]), int(r[1])] = r[2]
                nodes = meshes[-1]['nodes']
                ax.pcolormesh(nodes[..., 0], nodes[..., 1], vals, cmap=e['cmap'], vmin=e['clim'][0], vmax=e['clim'][1], shading='flat' if vals.shape[0] == nodes.shape[0] - 1 else 'nearest')
                meshes[-1].update(cmap=e['cmap'], clim=e['clim'])
            elif e['kind'] in ('contour_lines', 'contour_filled'):
                for il in np.unique(v[:, 0]).astype(int):
                    sel = v[v[:, 0] == il]
                    for sg in np.unique(sel[:, 2]):
                        P = sel[sel[:, 2] == sg][:, 3:]
                        col = e['colors'][min(il, len(e['colors']) - 1)] if e['colors'] else 'k'
                        if e['kind'] == 'contour_lines':
                            ls = e['linestyles'][min(il, len(e['linestyles']) - 1)] if e['linestyles'] else '-'
                            ax.plot(P[:, 0], P[:, 1], color=col, lw=e['linewidths'][min(il, len(e['linewidths']) - 1)] if e['linewidths'] else 1,
                                    linestyle='--' if 'dash' in ls or ls.startswith('(0, (') else '-')
                        else:
                            ax.add_patch(Polygon(P, closed=True, facecolor=col, edgecolor='none'))
            elif e['kind'] == 'image':
                n0, n1 = int(v[:, 0].max()) + 1, int(v[:, 1].max()) + 1
                img = v[:, 2].reshape(n0, n1)
                ax.imshow(img, extent=e['extent'], cmap=e['cmap'], vmin=e['clim'][0], vmax=e['clim'][1], origin='lower', aspect='auto')
                meshes.append(dict(cmap=e['cmap'], clim=e['clim']))
        for p in a.get('patches', []):
            if 'file' in p:
                _, v = read(d / p['file'])
                ax.add_patch(Polygon(v, closed=False, facecolor=p['facecolor'] or 'none', edgecolor=p['edgecolor'], linewidth=p['linewidth'], alpha=p['alpha']))
        for t in a.get('contour_labels', []):
            ax.text(t['x'], t['y'], t['text'], fontsize=t['fontsize'], rotation=t['rotation'], ha='center', va='center')
        for t in a.get('texts', []):
            if 'arrow_to' in t:
                ax.annotate(t['text'], xy=t['arrow_to'], xytext=(t['x'], t['y']), fontsize=t['fontsize'], ha=t['ha'], va=t['va'], color=t['color'],
                            arrowprops=dict(arrowstyle='->', lw=0.7))
            else:
                ax.text(t['x'], t['y'], t['text'], transform=transform(ax, t['coordinates']), fontsize=t['fontsize'], ha=t['ha'], va=t['va'], color=t['color'], rotation=t['rotation'])
        ax.set_xlim(a['xlim']); ax.set_ylim(a['ylim'])
        ax.set_xlabel(a['xlabel'], fontsize=a.get('xlabel_fontsize', plt.rcParams['axes.labelsize']))
        ax.set_ylabel(a['ylabel'], fontsize=a.get('ylabel_fontsize', plt.rcParams['axes.labelsize']))
        ax.set_title(a['title'], fontsize=a.get('title_fontsize', plt.rcParams['axes.titlesize']))
        if not a['axis_visible']: ax.set_axis_off()
        lines_drawn = [l for l in ax.get_lines()]
        for c in a.get('containers', []):
            i = c.get('data_line_index')
            if i is not None and i < len(lines_drawn) and c['label'] and not c['label'].startswith('_'):
                lines_drawn[i].set_label(c['label'])
        if a.get('legend') and a['legend']['texts']:
            hs, ls_ = ax.get_legend_handles_labels()
            keep = [(h, l) for h, l in zip(hs, ls_) if l in a['legend']['texts']]
            if keep: ax.legend(*zip(*keep), fontsize=7, frameon=False)
    for lg in spec.get('figure_legends', []):
        hs, ls_ = [], []
        for ax in fig.axes:
            for h, l in zip(*ax.get_legend_handles_labels()):
                if l in lg['texts'] and l not in ls_: hs.append(h); ls_.append(l)
        if hs: fig.legend(hs, ls_, loc='upper center', ncols=len(hs), frameon=False, fontsize=7)
    fig.savefig(d / 'replot.png', dpi=110, bbox_inches='tight', pad_inches=0.05)
    plt.close(fig)
    pdf = ROOT / json.loads((d / 'spec.json').read_text())['published_pdf']
    subprocess.run(['pdftoppm', '-r', '110', '-png', '-singlefile', str(pdf), str(d / 'original')], check=True)

def main():
    base = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / 'plot_data'
    for d in sorted(p.parent for p in base.glob('*/*/spec.json')):
        replot(d); print('redrawn', d.relative_to(base))

if __name__ == '__main__':
    main()
