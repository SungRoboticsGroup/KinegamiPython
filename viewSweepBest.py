"""
Interactive viewer for the best result on each chain of a parameter sweep.

Usage:
    python viewSweepBest.py "Trials Before Experiments/<sweep folder>" [--power 2] [--baseline LABEL]

For the chain picked with the radio buttons on the left (or the left/right arrow keys), shows side by side:
the starting tree, the final tree with the baseline settings (default: the default settings), and the
final tree of the settings that finished shortest on that chain. Each panel's title gives its total length.
All three are drawn in the neutral configuration; to step through the other configurations, run
    python viewConfigurations.py <sweep folder>/seed<k>/<settings>/final.tree --configs <sweep folder>/seed<k>/configs.json
"""
import argparse
import csv
import os
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons
from KinematicTree import *

def loadBest(sweepDir, power, baseline):
    rows = [r for r in csv.DictReader(open(os.path.join(sweepDir, 'summary.csv')))
            if not r['error'] and (power is None or r['power'] == str(power))]
    seeds = sorted({int(r['seed']) for r in rows})
    chains = []
    for s in seeds:
        runs = [r for r in rows if int(r['seed']) == s]
        best = min(runs, key=lambda r: float(r['final_length']))
        base = next((r for r in runs if r['label'] == baseline), None)
        chains.append(dict(seed=s, best=best, base=base, start=float(best['start_length']),
                           lower=float(best['lower_bound'])))
    return chains

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sweep', help="sweep folder containing summary.csv and seed<k>/")
    parser.add_argument('--power', type=int, default=2, help="only consider runs with this power (default 2; 0 for all)")
    parser.add_argument('--baseline', default='box_nostraight_align10_power2', help="settings label for the middle panel")
    args = parser.parse_args()

    chains = loadBest(args.sweep, args.power or None, args.baseline)
    cache = {}
    def tree(path):
        if path not in cache:
            cache[path] = loadKinematicTree(path)
        return cache[path]

    fig = plt.figure(figsize=(16, 6))
    axes = [fig.add_subplot(1, 3, k + 1, projection='3d') for k in range(3)]
    fig.subplots_adjust(left=0.1, right=0.99, wspace=0.02)
    radioAx = fig.add_axes([0.005, 0.3, 0.07, 0.04 * len(chains)])
    radio = RadioButtons(radioAx, [f"seed {c['seed']}" for c in chains])
    state = {'index': 0}

    def draw(index):
        c = chains[index]
        seedDir = os.path.join(args.sweep, f"seed{c['seed']}")
        panels = [(os.path.join(seedDir, 'start.tree'), f"Start: {c['start']:.1f}")]
        if c['base'] is not None:
            panels.append((os.path.join(seedDir, args.baseline, 'final.tree'),
                           f"{args.baseline}: {float(c['base']['final_length']):.1f}"))
        else:
            panels.append((None, f"{args.baseline}: not run"))
        b = c['best']
        ratio = f" ({float(b['final_length']) / float(c['base']['final_length']):.2f} of baseline)" if c['base'] else ""
        panels.append((os.path.join(seedDir, b['label'], 'final.tree'),
                       f"Best: {b['label']}\n{float(b['final_length']):.1f}{ratio}"))
        for ax, (path, title) in zip(axes, panels):
            ax.clear()
            if path is not None:
                tree(path).addToPlot(ax, showJointPoses=False, showLinkPath=False)
                ax.set_aspect('equal')
            ax.set_title(title, fontsize=9)
            ax.set_axis_off()
        fig.suptitle(f"Seed {c['seed']}: lower bound {c['lower']:.1f}", fontsize=11)
        fig.canvas.draw_idle()

    def onRadio(label):
        state['index'] = [f"seed {c['seed']}" for c in chains].index(label)
        draw(state['index'])

    def onKey(event):
        step = {'right': 1, 'left': -1}.get(event.key)
        if step:
            radio.set_active((state['index'] + step) % len(chains))  # calls onRadio

    radio.on_clicked(onRadio)
    fig.canvas.mpl_connect('key_press_event', onKey)
    draw(0)
    plt.show()

if __name__ == "__main__":
    main()
