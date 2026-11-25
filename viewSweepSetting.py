"""
Interactive viewer for one setting's final trees across the chains of one or more parameter sweeps, in each of
the chain's configurations, for inspecting them (e.g. for collisions by eye).

Usage:
    python viewSweepSetting.py "Trials Before Experiments/<sweep>" ["Trials Before Experiments/<sweep2>" ...]
                               --label box_straight_align10_power1_parentlock

Pick the tree with the radio buttons on the left (or the left/right arrow keys) and the configuration with the
ones on the right (or the up/down arrow keys). The title gives the tree's total length and how many collisions
the checker finds in that configuration.
"""
import argparse
import copy
import csv
import json
import os
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons
from KinematicTree import *

def posed(tree, configuration):
    tree = copy.deepcopy(tree)
    for i in range(len(tree.Joints)):
        if not isWaypoint(tree.Joints[i]):
            tree.setJointState(i, configuration[i])
        tree.Joints[i].recomputeCollisionCapsules()
    return tree

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sweeps', nargs='+', help="sweep folders containing summary.csv and seed<k>/")
    parser.add_argument('--label', required=True, help="settings label whose final trees to show")
    args = parser.parse_args()

    entries = []
    for sweep in args.sweeps:
        joints = json.load(open(os.path.join(sweep, 'run_info.json')))['settings'].get('joints', '?')
        for r in csv.DictReader(open(os.path.join(sweep, 'summary.csv'))):
            if r['label'] == args.label and not r['error']:
                seedDir = os.path.join(sweep, f"seed{r['seed']}")
                entries.append(dict(name=f"{joints}J seed {r['seed']}", tree=os.path.join(seedDir, r['label'], 'final.tree'),
                                    configs=json.load(open(os.path.join(seedDir, 'configs.json'))),
                                    length=float(r['final_length']), lower=float(r['lower_bound'])))
    if not entries:
        raise SystemExit(f"no finished runs with label {args.label}")

    trees = {}
    state = {'tree': 0, 'config': 0}
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_axes([0.12, 0.02, 0.76, 0.9], projection='3d')
    treeRadio = RadioButtons(fig.add_axes([0.005, 0.3, 0.1, 0.035 * len(entries)]), [e['name'] for e in entries])
    maxConfigs = max(len(e['configs']) for e in entries)
    configRadio = RadioButtons(fig.add_axes([0.89, 0.4, 0.1, 0.035 * maxConfigs]),
                               ["neutral" if k == 0 else f"config {k}" for k in range(maxConfigs)])

    def draw():
        e = entries[state['tree']]
        k = min(state['config'], len(e['configs']) - 1)
        if e['tree'] not in trees:
            trees[e['tree']] = loadKinematicTree(e['tree'])
        tree = posed(trees[e['tree']], e['configs'][k])
        collisions = tree.detectCollisions()
        elev, azim = ax.elev, ax.azim
        ax.clear()
        tree.addToPlot(ax, showJointPoses=False, showLinkPath=False)
        ax.set_aspect('equal')
        ax.set_axis_off()
        ax.view_init(elev=elev, azim=azim)
        ax.set_title(f"{args.label}: {e['name']}, {'neutral' if k == 0 else f'config {k}'}\n"
                     f"length {e['length']:.1f} ({e['length'] / e['lower']:.1f}x lower bound); "
                     f"checker: {collisions} collision{'s' if collisions != 1 else ''}",
                     fontsize=10, color='red' if collisions else 'black')
        fig.canvas.draw_idle()

    def onTree(label):
        state['tree'] = [e['name'] for e in entries].index(label)
        draw()

    def onConfig(label):
        state['config'] = 0 if label == "neutral" else int(label.split()[-1])
        draw()

    def onKey(event):
        if event.key in ('left', 'right'):
            treeRadio.set_active((state['tree'] + (1 if event.key == 'right' else -1)) % len(entries))
        elif event.key in ('up', 'down'):
            configRadio.set_active((state['config'] + (1 if event.key == 'down' else -1)) % maxConfigs)

    treeRadio.on_clicked(onTree)
    configRadio.on_clicked(onConfig)
    fig.canvas.mpl_connect('key_press_event', onKey)
    draw()
    plt.show()

if __name__ == "__main__":
    main()
