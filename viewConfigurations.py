"""
Interactive viewer for a kinematic tree in each of the configurations it was optimized for.

Usage:
    python viewConfigurations.py <tree file> [--configs <configs json>] [--threshold <clearance>]

The configs file defaults to the <n>_configs.json in the run folder that contains the tree file
(e.g. for "Trials Before Experiments/<run>/0/<traversal>/final.tree" it finds "<run>/0_configs.json").

In the window:
- Pick a configuration with the radio buttons on the left (or the left/right arrow keys).
  "All overlaid" draws every configuration at once in different colors.
- Pairs of parts closer than the clearance threshold in the shown configuration are joined by a
  dashed red line between their closest points, labeled with the gap.
- "Tight in other configs" highlights in orange, in the current pose, the parts that come within the
  threshold of something in some other configuration, labeled with which one. Use this in the neutral
  configuration to see whether a bend that looks unnecessary there is avoiding a collision elsewhere.

Clearances are surface-to-surface distances between the pairs the optimizer's collision checking
compares (KinematicTree.collisionCandidates, with joints split into halves and the path window),
sampled along each part's centerline, so they're
approximate to within the sampling spacing.
"""
import argparse
import copy
import json
import os
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RadioButtons, CheckButtons
from KinematicTree import *

SAMPLE_DENSITY = 4  # centerline samples per unit length for clearance estimates
OVERLAY_COLORS = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple',
                  'tab:brown', 'tab:pink', 'tab:olive', 'tab:cyan', 'tab:gray']

def findConfigsFile(treePath):
    # Walk up from the tree file looking for <n>_configs.json, preferring the one whose n matches
    # the folder we came up through (runs with several chains have 0_configs.json, 1_configs.json, ...)
    directory = os.path.dirname(os.path.abspath(treePath))
    cameFrom = None
    for _ in range(4):
        candidates = [f for f in os.listdir(directory) if re.match(r'^\d+_configs\.json$', f)]
        if cameFrom is not None and f"{cameFrom}_configs.json" in candidates:
            return os.path.join(directory, f"{cameFrom}_configs.json")
        if len(candidates) == 1:
            return os.path.join(directory, candidates[0])
        cameFrom = os.path.basename(directory)
        directory = os.path.dirname(directory)
    raise FileNotFoundError(f"No <n>_configs.json found above {treePath}; pass --configs")

def configuredTree(tree, configuration):
    configured = copy.deepcopy(tree)
    for i in range(len(configured.Joints)):
        if not isWaypoint(configured.Joints[i]):
            configured.setJointState(i, configuration[i])
            configured.Joints[i].recomputeCollisionCapsules()
    return configured

def partName(tree, kind, i):
    if kind == 'Link':
        return f"L{i} ({tree.Parents[i]}->{i})"
    return f"{'W' if isWaypoint(tree.Joints[i]) else 'J'}{i}"

def isAttachment(obj1, points1, i1, obj2, points2, i2):
    # Same end test KinematicTree.collision uses: when both closest points are within r of an end, the
    # contact is real only if the end discs cross or a cross-section of either part near its closest point
    # reaches inside the other (KinematicTree.crossSectionsReach). Otherwise it's only the round caps of the
    # center-line distance meeting, e.g. where two connected parts attach.
    def nearestEnd(obj, points, i):
        startDist = np.linalg.norm(points[i] - points[0])
        endDist = np.linalg.norm(points[i] - points[-1])
        if min(startDist, endDist) >= obj.r:
            return None
        return obj.startCircle(forward=False) if startDist <= endDist else obj.endCircle(forward=True)
    circle1 = nearestEnd(obj1, points1, i1)
    circle2 = nearestEnd(obj2, points2, i2)
    if circle1 is None or circle2 is None:
        return False
    isTrueContact = discs_cross(circle1, circle2) \
                    or KinematicTree.crossSectionsReach(obj1, points1, points1[i1], obj2, 0) \
                    or KinematicTree.crossSectionsReach(obj2, points2, points2[i2], obj1, 0)
    return not isTrueContact

def clearances(tree, collisionMatrices, threshold):
    # Returns [(clearance, (kind1, i1), (kind2, i2), point1, point2)] for every collidable pair closer than threshold,
    # using the checker's pairs and path window (KinematicTree.collisionCandidates); a joint's halves are
    # reported as the joint, keeping the closest of their clearances
    candidates, structure = tree.collisionCandidates(None, collisionMatrices)
    best = {}
    for x, y, obj1, obj2, keep in candidates:
        if obj1.length() <= 1e-6 or obj2.length() <= 1e-6:
            continue
        points1 = obj1.interpolate(density=SAMPLE_DENSITY)
        allPoints2 = obj2.interpolate(density=SAMPLE_DENSITY)
        points2 = allPoints2
        if keep is not None:
            # only the stretch of obj2 outside the window (samples are evenly spaced by arc length)
            arc = np.linspace(0, obj2.length(), len(allPoints2))
            points2 = allPoints2[(arc >= keep[0]) & (arc <= keep[1])]
            if len(points2) == 0:
                continue
        # skip pairs whose sample bounding boxes are clearly farther apart than the threshold
        reach = obj1.r + obj2.r + threshold
        if np.any(points1.min(axis=0) - points2.max(axis=0) > reach) or np.any(points2.min(axis=0) - points1.max(axis=0) > reach):
            continue
        dists2 = obj1.sdf(points2) - obj2.r
        i2 = np.argmin(dists2)
        if keep is None:
            dists1 = obj2.sdf(points1) - obj1.r
            i1 = np.argmin(dists1)
            # like KinematicTree.collision, require both directions to agree, which avoids end-cap artifacts
            clearance = max(dists1[i1], dists2[i2])
            if isAttachment(obj1, points1, i1, obj2, allPoints2, int(np.argmin(np.linalg.norm(allPoints2 - points2[i2], axis=1)))):
                continue
        else:
            # as in KinematicTree.collision: from the checked stretch, paired with the nearest point of obj1
            clearance = dists2[i2]
            i1 = int(np.argmin(np.linalg.norm(points1 - points2[i2], axis=1)))
        if clearance < threshold:
            (n1, k1), (n2, k2) = structure.key(x), structure.key(y)
            key = ((k1, n1), (k2, n2))
            if key not in best or clearance < best[key][0]:
                best[key] = (clearance, key[0], key[1], points1[i1], points2[i2])
    close = sorted(best.values(), key=lambda c: c[0])
    return close

def centerline(tree, kind, i):
    obj = tree.Joints[i] if kind == 'Joint' else tree.Links[i]
    return obj.interpolate(density=SAMPLE_DENSITY)

class ConfigurationViewer:
    def __init__(self, tree, configurations, threshold, title):
        self.tree = tree
        self.configurations = configurations
        self.threshold = threshold
        self.trees = [configuredTree(tree, c) for c in configurations]
        collisionMatrices = tree.buildCollisionMatrices()
        print(f"Computing clearances in {len(configurations)} configurations...")
        # Compute out to a wider range so the closest pair can be reported even when nothing is within threshold
        wide = [clearances(t, collisionMatrices, max(threshold, 10 * tree.r)) for t in self.trees]
        self.closest = [close[0] if close else None for close in wide]
        self.close = [[c for c in close if c[0] < threshold] for close in wide]
        for k, close in enumerate(self.close):
            print(f"  {self.configLabel(k)}: {len(close)} pairs within {threshold}; closest {self.closestText(k)}")

        self.labels = [self.configLabel(k) for k in range(len(configurations))] + ["All overlaid"]
        self.current = 0
        self.showOther = True

        self.fig = plt.figure(figsize=(15, 9))
        self.fig.suptitle(title)
        self.ax = self.fig.add_axes([0.24, 0.02, 0.75, 0.92], projection='3d')
        self.setLimits()
        self.ax.view_init(elev=25, azim=-60)

        radioAx = self.fig.add_axes([0.01, 0.62, 0.21, 0.03 * len(self.labels) + 0.02])
        self.radio = RadioButtons(radioAx, self.labels, active=0)
        self.radio.on_clicked(self.onSelect)
        checkAx = self.fig.add_axes([0.01, 0.55, 0.21, 0.05])
        self.check = CheckButtons(checkAx, [f"Tight in other configs"], [self.showOther])
        self.check.on_clicked(self.onToggle)
        self.infoText = self.fig.text(0.01, 0.53, "", va='top', ha='left', family='monospace', fontsize=8)
        self.fig.canvas.mpl_connect('key_press_event', self.onKey)
        self.draw()

    def closestText(self, k):
        if self.closest[k] is None:
            return f"none within {max(self.threshold, 10 * self.tree.r):g}"
        clearance, key1, key2, _, _ = self.closest[k]
        return f"{clearance:.2f} ({partName(self.trees[k], *key1)} - {partName(self.trees[k], *key2)})"

    def configLabel(self, k):
        c = self.configurations[k]
        if all(abs(s) < 1e-12 for s in c):
            return f"{k}: neutral"
        moved = [f"{i}:{s:.2f}" for i, s in enumerate(c) if abs(s) > 1e-12]
        return f"{k}: " + ", ".join(moved[:3]) + (" ..." if len(moved) > 3 else "")

    def setLimits(self):
        # Fixed cube-shaped limits covering every configuration, so switching doesn't rescale the view
        points = np.vstack([centerline(t, kind, i)
                            for t in self.trees
                            for kind, objs in (('Joint', t.Joints), ('Link', t.Links))
                            for i, obj in enumerate(objs) if obj is not None and obj.length() > 1e-6])
        center = (points.max(axis=0) + points.min(axis=0)) / 2
        half = (points.max(axis=0) - points.min(axis=0)).max() / 2 + 2 * self.tree.r
        self.limits = [(c - half, c + half) for c in center]

    def draw(self):
        elev, azim = self.ax.elev, self.ax.azim
        self.ax.cla()
        lines = []
        if self.current == len(self.configurations):
            for k, t in enumerate(self.trees):
                color = OVERLAY_COLORS[k % len(OVERLAY_COLORS)]
                label = self.labels[k]  # attached to the first part drawn, so each config appears once in the legend
                for kind, objs in (('Joint', t.Joints), ('Link', t.Links)):
                    for i, obj in enumerate(objs):
                        if obj is not None and obj.length() > 1e-6:
                            p = centerline(t, kind, i)
                            self.ax.plot(p[:, 0], p[:, 1], p[:, 2], color=color, linewidth=2, label=label)
                            label = None
            self.ax.legend(loc='upper left', fontsize=8)
            lines.append("Centerlines of every configuration.")
            lines.append("Closest pair in each:")
            lines += [f"  {self.labels[k][:24]}: {self.closestText(k)}" for k in range(len(self.trees))]
        else:
            k = self.current
            t = self.trees[k]
            t.addToPlot(self.ax, showJointPoses=False)
            here = self.close[k]
            for clearance, key1, key2, p1, p2 in here:
                self.ax.plot(*zip(p1, p2), color='red', linestyle='--', linewidth=2)
                mid = (p1 + p2) / 2
                self.ax.text(*mid, f"{clearance:.2f}", color='red', fontsize=8)
            lines.append(f"Closest pair here: {self.closestText(k)}")
            lines.append(f"Within {self.threshold} here ({len(here)}):")
            lines += [f"  {clearance:6.2f}  {partName(t, *key1)} - {partName(t, *key2)}"
                      for clearance, key1, key2, _, _ in here[:12]]
            if len(here) > 12:
                lines.append(f"  ... {len(here) - 12} more")

            if self.showOther:
                # For each part, the tightest clearance it has in any other configuration
                tightest = {}
                for other, close in enumerate(self.close):
                    if other == k:
                        continue
                    for clearance, key1, key2, _, _ in close:
                        for key, partner in ((key1, key2), (key2, key1)):
                            if key not in tightest or clearance < tightest[key][0]:
                                tightest[key] = (clearance, partner, other)
                for (kind, i), (clearance, partner, other) in tightest.items():
                    p = centerline(t, kind, i)
                    self.ax.plot(p[:, 0], p[:, 1], p[:, 2], color='darkorange', linewidth=5, alpha=0.8)
                    mid = p[len(p) // 2]
                    self.ax.text(*mid, f" {partName(t, kind, i)}: cfg {other}", color='darkorange', fontsize=8)
                lines.append("")
                lines.append(f"Within {self.threshold} in other configs (orange):")
                shown = set()
                for (kind, i), (clearance, partner, other) in sorted(tightest.items(), key=lambda item: item[1][0]):
                    pairKey = frozenset([(kind, i), partner])
                    if pairKey in shown:
                        continue
                    shown.add(pairKey)
                    lines.append(f"  {clearance:6.2f}  {partName(t, kind, i)} - {partName(t, *partner)}  (cfg {other})")
                    if len(shown) >= 12:
                        lines.append("  ...")
                        break

        self.ax.set_xlim(*self.limits[0])
        self.ax.set_ylim(*self.limits[1])
        self.ax.set_zlim(*self.limits[2])
        self.ax.set_box_aspect((1, 1, 1))
        self.ax.view_init(elev=elev, azim=azim)
        self.ax.set_title(self.labels[self.current])
        self.infoText.set_text("\n".join(lines))
        self.fig.canvas.draw_idle()

    def onSelect(self, label):
        self.current = self.labels.index(label)
        self.draw()

    def onToggle(self, label):
        self.showOther = not self.showOther
        self.draw()

    def onKey(self, event):
        if event.key in ('right', 'left'):
            step = 1 if event.key == 'right' else -1
            self.radio.set_active((self.current + step) % len(self.labels))  # triggers onSelect

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('tree', help="path to a .tree file")
    parser.add_argument('--configs', help="path to the configurations json (default: found from the run folder)")
    parser.add_argument('--threshold', type=float, default=None,
                        help="clearance below which pairs are highlighted (default: the tree's tube radius)")
    args = parser.parse_args()

    tree = loadKinematicTree(args.tree)
    configsPath = args.configs or findConfigsFile(args.tree)
    with open(configsPath) as f:
        configurations = json.load(f)
    threshold = args.threshold if args.threshold is not None else tree.r
    print(f"Tree: {args.tree}\nConfigurations: {configsPath} ({len(configurations)})")
    title = f"{os.path.basename(args.tree)}: total length {tree.totalLength():.1f}"
    viewer = ConfigurationViewer(tree, configurations, threshold, title)
    plt.show()
    return viewer

if __name__ == "__main__":
    main()
