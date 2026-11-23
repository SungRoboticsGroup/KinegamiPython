"""
Interactive viewer for every result on one chain of a parameter sweep.

Usage:
    python viewSweepChain.py "Trials Before Experiments/<sweep folder>" --seed 2

Shows the starting tree and the final tree of each settings label run on that chain, side by side in 3D,
each titled with its total length (and a warning if its final tree collides). Drag a panel to rotate it;
the panels rotate together. All trees are drawn in the neutral configuration.
"""
import argparse
import csv
import math
import os
import matplotlib.pyplot as plt
from KinematicTree import *

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sweep', help="sweep folder containing summary.csv and seed<k>/")
    parser.add_argument('--seed', type=int, required=True, help="which chain to show")
    args = parser.parse_args()

    rows = [r for r in csv.DictReader(open(os.path.join(args.sweep, 'summary.csv')))
            if int(r['seed']) == args.seed and not r['error']]
    if not rows:
        raise SystemExit(f"no finished runs for seed {args.seed} in {args.sweep}")
    rows.sort(key=lambda r: float(r['final_length']))
    seedDir = os.path.join(args.sweep, f"seed{args.seed}")
    lower = float(rows[0]['lower_bound'])

    panels = [(os.path.join(seedDir, 'start.tree'), f"Start: {float(rows[0]['start_length']):.1f}")]
    for r in rows:
        title = f"{r['label']}\n{float(r['final_length']):.1f} ({float(r['final_length']) / lower:.1f}x lower bound)"
        if r.get('final_colliding_configs') not in (None, '', '0'):
            title += "\nCOLLIDES in final tree"
        panels.append((os.path.join(seedDir, r['label'], 'final.tree'), title))

    cols = min(3, len(panels))
    nrows = math.ceil(len(panels) / cols)
    fig = plt.figure(figsize=(5.5 * cols, 5 * nrows))
    axes = []
    for k, (path, title) in enumerate(panels):
        ax = fig.add_subplot(nrows, cols, k + 1, projection='3d')
        loadKinematicTree(path).addToPlot(ax, showJointPoses=False, showLinkPath=False)
        ax.set_aspect('equal')
        ax.set_title(title, fontsize=9, color='red' if 'COLLIDES' in title else 'black')
        ax.set_axis_off()
        axes.append(ax)
    fig.suptitle(f"{os.path.basename(os.path.normpath(args.sweep))}, seed {args.seed}: lower bound {lower:.1f}",
                 fontsize=11)
    fig.tight_layout()

    def onMove(event):
        # keep all panels at the view of the one being dragged
        if event.inaxes in axes:
            for ax in axes:
                if ax is not event.inaxes:
                    ax.view_init(elev=event.inaxes.elev, azim=event.inaxes.azim)
            fig.canvas.draw_idle()
    fig.canvas.mpl_connect('motion_notify_event', onMove)
    plt.show()

if __name__ == "__main__":
    main()
