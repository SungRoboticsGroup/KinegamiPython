"""
Debugging harness for bunched waypoints: runs single waypoint visits on saved examples and compares
the result with known better placements, without running a whole tree optimization.

Usage:
    python bunchingExamples.py [names...] [--sampling box|multiscale] [--straight] [--alignment W]
                               [--repeats N] [--seed S] [--serial]
    python bunchingExamples.py W8 --show [--sampling ...]   # plot one example (current / result / references)

The examples live in debugExamples/bunching/examples.json (see its "description"). Each records a tree, the
waypoint index, the configurations and failurePenalty from the run it came from, and two reference placements
found by a heavy search: the best under the waypoint objective with alignmentWeight 10 ("ref obj"), and the
best under plain linkLoss ("ref len"). Both are collision-free in every configuration.

For each example and repeat, the table reports linkLoss (power 2, the length term alone) and the waypoint
objective with the default alignmentWeight 10, before and after the visit, next to the references. "gap
closed" is how much of the distance from the current linkLoss to the "ref len" linkLoss the visit closed.
"""
import argparse
import json
import os
import io
import contextlib
import numpy as np
import matplotlib.pyplot as plt
from optimizationFunctions import *

EXAMPLES_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'debugExamples', 'bunching')

def loadExamples(names=None):
    with open(os.path.join(EXAMPLES_DIR, 'examples.json')) as f:
        data = json.load(f)
    trees = {}
    examples = []
    for e in data['examples']:
        if names and e['name'] not in names:
            continue
        if e['tree'] not in trees:
            trees[e['tree']] = loadKinematicTree(os.path.join(EXAMPLES_DIR, e['tree']))
        examples.append(dict(e, tree=trees[e['tree']], configurations=data['configurations'],
                             failurePenalty=data['failurePenalty']))
    if names:
        missing = set(names) - {e['name'] for e in examples}
        if missing:
            raise ValueError(f"Unknown examples {sorted(missing)}; available: {[e['name'] for e in data['examples']]}")
    return examples

def alignedOrientation(tree, index):
    waypoint, parent = tree.Joints[index], tree.Joints[tree.Parents[index]]
    return (parent.DistalDubinsFrame() * waypoint.ProximalDubinsFrame().inv()).R @ waypoint.Pose.R

def measure(tree, index, aligned, alignmentWeight=10):
    # linkLoss (power 2) and the waypoint objective with the given alignment weight, without collision penalty
    L = linkLoss(tree, index)
    misalignment = rotationAngle(tree.Joints[index].Pose.R.T @ aligned)
    return L, L + alignmentWeight * misalignment

def applyMove(tree, index, move):
    moved = tree.copyAbbreviatedSelf()
    if not moved.transformJoint(index, SE3(np.array(move), check=False), propogate=False, safe=True,
                                relative=True, recomputeBoundingBall=False):
        raise ValueError("reference move could not be applied")
    return moved

def runVisit(example, seed, parallelize=True, **settings):
    tree, index = example['tree'], example['index']
    np.random.seed(seed)
    record = {}
    with contextlib.redirect_stdout(io.StringIO()):
        result, _ = optimizeWaypointPlacement(tree, index, maxiter=50, tol=tree.r / 10,
                                              failurePenalty=example['failurePenalty'],
                                              configurations=example['configurations'], parallelize=parallelize,
                                              verbose=False, visitRecord=record, **settings)
    return result, record

def runAll(examples, repeats, seed, parallelize, settings):
    print(f"Settings: {settings or 'defaults'}; {repeats} repeat(s) per example, seeds {seed}..{seed + repeats - 1}\n")
    header = f"{'example':8} {'seed':>4}  {'linkLoss: before':>16} {'after':>8} {'ref len':>8} {'gap closed':>10}   " \
             f"{'objective(w=10): before':>23} {'after':>8} {'ref obj':>8}  {'retried/reverted':>16} {'time':>6}"
    print(header)
    print('-' * len(header))
    summary = []
    for e in examples:
        tree, index = e['tree'], e['index']
        aligned = alignedOrientation(tree, index)
        L0, obj0 = measure(tree, index, aligned)
        refLen, refObj = e['reference']['linkLoss']['linkLoss'], e['reference']['objective']['objective']
        closed = []
        for r in range(repeats):
            start = time.time()
            result, record = runVisit(e, seed + r, parallelize, **settings)
            L, obj = measure(result, index, aligned)
            gap = (L0 - L) / (L0 - refLen) if L0 - refLen > 1e-9 else float('nan')
            closed.append(gap)
            flags = f"{'R' if record.get('retried') else '-'}/{'V' if record.get('reverted') else '-'}"
            print(f"{e['name']:8} {seed + r:>4}  {L0:16.2f} {L:8.2f} {refLen:8.2f} {gap:10.0%}   "
                  f"{obj0:23.2f} {obj:8.2f} {refObj:8.2f}  {flags:>16} {time.time() - start:5.1f}s")
        summary.append((e['name'], np.nanmean(closed)))
    print("\nMean share of the linkLoss gap closed: " + ", ".join(f"{n} {g:.0%}" for n, g in summary)
          + f"; overall {np.nanmean([g for _, g in summary]):.0%}")

def show(example, seed, parallelize, settings, savePath=None):
    tree, index = example['tree'], example['index']
    result, _ = runVisit(example, seed, parallelize, **settings)
    versions = [('current', tree, 'gray'), (f"optimizer {settings or '(defaults)'}", result, 'tab:blue'),
                ('ref obj', applyMove(tree, index, example['reference']['objective']['move']), 'tab:green'),
                ('ref len', applyMove(tree, index, example['reference']['linkLoss']['move']), 'tab:purple')]
    # the parts whose placement depends on this waypoint: its incoming and outgoing links, and its neighbors
    parent = tree.Parents[index]
    links = [index] + list(tree.Children[index])
    joints = [parent, index] + list(tree.Children[index])
    aligned = alignedOrientation(tree, index)
    ax = plt.figure(figsize=(10, 9)).add_subplot(projection='3d')
    allPoints = []
    for label, t, color in versions:
        L, obj = measure(t, index, aligned)
        first = True
        for i in links:
            if t.Links[i].path.length > 1e-6:
                p = t.Links[i].interpolate(density=10)
                allPoints.append(p)
                ax.plot(p[:, 0], p[:, 1], p[:, 2], color=color, linewidth=3 if label == 'current' else 2,
                        label=f"{label}: linkLoss {L:.1f}, obj {obj:.1f}" if first else None)
                first = False
        pose = t.Joints[index].Pose
        ax.quiver(*pose.t, *(pose.R[:, 2] * 2 * tree.r), color=color)  # waypoint's z-axis
    for j in joints:
        if j == index:
            continue
        p = tree.Joints[j].Pose.t
        ax.scatter(*p, color='black', s=30)
        ax.text(*p, f" {'W' if isWaypoint(tree.Joints[j]) else 'J'}{j}", fontsize=9)
    points = np.vstack(allPoints)
    center = (points.max(axis=0) + points.min(axis=0)) / 2
    half = (points.max(axis=0) - points.min(axis=0)).max() / 2 + tree.r
    ax.set_xlim(center[0] - half, center[0] + half)
    ax.set_ylim(center[1] - half, center[1] + half)
    ax.set_zlim(center[2] - half, center[2] + half)
    ax.set_box_aspect((1, 1, 1))
    ax.set_title(f"{example['name']}: {example['note']}", fontsize=9, wrap=True)
    ax.legend(loc='upper left', fontsize=8)
    if savePath:
        plt.savefig(savePath, dpi=120, bbox_inches='tight')
    else:
        plt.show()

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('names', nargs='*', help="examples to run (default: all)")
    parser.add_argument('--sampling', choices=WAYPOINT_SAMPLINGS, help="waypointSampling (default: box)")
    parser.add_argument('--straight', action='store_true', help="straightThroughGuess=True")
    parser.add_argument('--alignment', type=float, help="alignmentWeight (default: 10)")
    parser.add_argument('--repeats', type=int, default=3, help="runs per example with different seeds (default 3)")
    parser.add_argument('--seed', type=int, default=0, help="first seed (default 0)")
    parser.add_argument('--serial', action='store_true', help="don't parallelize particle evaluation")
    parser.add_argument('--show', action='store_true', help="plot the (single) named example instead of the table")
    parser.add_argument('--save', help="with --show: save the plot to this file instead of opening a window")
    args = parser.parse_args()

    settings = {}
    if args.sampling:
        settings['waypointSampling'] = args.sampling
    if args.straight:
        settings['straightThroughGuess'] = True
    if args.alignment is not None:
        settings['alignmentWeight'] = args.alignment
    examples = loadExamples(args.names or None)
    if args.show:
        if len(examples) != 1:
            raise SystemExit("--show needs exactly one example name")
        show(examples[0], args.seed, not args.serial, settings, args.save)
    else:
        runAll(examples, args.repeats, args.seed, not args.serial, settings)

if __name__ == "__main__":
    main()
