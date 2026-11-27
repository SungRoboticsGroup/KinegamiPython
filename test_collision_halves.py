"""
Tests for collision checking with joints split into halves and the path window (KinematicTree.collisionCandidates),
against ground truth from sampling the parts' solids.

Run: python test_collision_halves.py
"""
import json
import os
import numpy as np
from spatialmath import SE3
from KinematicTree import *

r = 1.0
RNG = np.random.default_rng(0)

def place(tree, parent, joint):
    return tree.addJoint(parent, joint, relative=False, safe=False, fixedPosition=True, fixedOrientation=True)

def transverse(h, pose):
    tmax = 2 * np.arctan(h / r)  # the range a transverse revolute of half-length h allows
    return TransverseRevolute(r=r, Pose=pose, minAngle=-tmax, maxAngle=tmax)

def adjacentRevolutes(hs, twists=None):
    # Transverse revolutes of half-lengths hs in a row along x, each starting where the previous one ends
    # (empty links between them); twists: rotation of each one's axis about the path (default 0: all bend about z)
    twists = twists if twists is not None else [0] * len(hs)
    tree = KinematicTree(transverse(hs[0], SE3()))
    x = hs[0]
    for i, h in enumerate(hs[1:]):
        x += h
        place(tree, i, transverse(h, SE3.Trans(x, 0, 0) @ SE3.Rx(twists[i + 1])))
        x += h
    return tree

def bent(tree, fractions):
    # each joint at the given fraction of its range (positive = the same side when the axes aren't twisted)
    for i, f in enumerate(fractions):
        tree.setJointState(i, f * tree.Joints[i].maxAngle)
    return tree

def solidPoints(half, ns=120, nr=8, nphi=36):
    # points filling a joint half's solid (its cross-section discs along its center line)
    a, b = half.endpoints()
    t = half.direction()
    u = np.cross(t, [0, 0, 1.0] if abs(t[2]) < 0.9 else [1.0, 0, 0])
    u /= np.linalg.norm(u)
    v = np.cross(t, u)
    s = np.linspace(0, 1, ns)[:, None, None, None]
    rho = (np.arange(1, nr + 1) / nr)[None, :, None, None] * (1 - 1e-9) * half.r
    phi = np.linspace(0, 2 * np.pi, nphi, endpoint=False)[None, None, :, None]
    points = (a + s * (b - a)) + rho * (np.cos(phi) * u + np.sin(phi) * v)
    return np.vstack([points.reshape(-1, 3), a + s.reshape(-1, 1) * (b - a)])

def solidOverlap(half1, half2):
    # ground truth: how deep the sampled solid of either half reaches into the other (positive = overlap)
    return max(-np.min(half2.sdf(solidPoints(half1))), -np.min(half1.sdf(solidPoints(half2))))

def truthPairs(tree, threshold):
    # pairs of joints (index order) with halves whose solids overlap by more than threshold
    halves = [(i, JointHalf(tree.Joints[i], w)) for i in range(len(tree.Joints)) for w in ('proximal', 'distal')]
    return sorted({(i, j) for i, a in halves for j, b in halves if i < j and solidOverlap(a, b) > threshold})

def checkerPairs(tree):
    return sorted({(a, b) for (a, _), (b, _) in tree.getCollidingPairs()})

def check(name, condition):
    print(("✓ " if condition else "✗ FAILED: ") + name)
    assert condition, name

def testFixedCases():
    # (half-lengths, fractions of range, expected colliding joint pairs)
    cases = [
        ([2, 2], [0.98, 0.98], []),      # two long joints folded toward each other: their outer arms stay apart
        ([2, 2], [1, 1], []),            # and only touch at full range
        ([2, 2], [0.98, -0.98], []),     # folded to opposite sides: the second's end comes back beside the first's
                                         # distal half (a false positive of the old end-plane test)
        ([1, 1], [1, 1], []),            # +-90 degree joints (as in the random test chains): only touch
        ([3, 1, 3], [0.98] * 3, [(0, 2)]),    # three in a row: the outer joints' far halves meet across the middle one
        ([3, 0.8, 3], [0.98] * 3, [(0, 2)]),
        ([3, 0.3, 3], [1] * 3, [(0, 1), (1, 2)]),  # a fully folded long joint against the far half of a very short one
    ]
    # Each expected pair must also overlap in the sampled solids; the solids can show more with a short middle
    # joint (h < r), whose elbow reaches into its neighbors (not counted, see testShortJointElbow)
    for hs, fractions, expected in cases:
        tree = bent(adjacentRevolutes(hs), fractions)
        found, truth = checkerPairs(tree), truthPairs(tree, 1e-3)
        check(f"joints {hs} at {fractions} of range: checker {found}, solids {truth}, expected {expected}",
              found == expected and set(expected) <= set(truth))

def testShortJointElbow():
    # A joint with range under +-90 degrees (half-length h < r) bent far: in the two-cylinder joint model, the
    # inner edge of its distal half swings back by r sin(angle) > h, past its proximal end, into its own
    # proximal half and the end of the part before it. Not counted, like a joint's two halves against each
    # other (the window covers it: it always lies within r of the joint's proximal end).
    tree = bent(adjacentRevolutes([0.6, 0.6]), [0, 1])
    elbow = solidOverlap(JointHalf(tree.Joints[0], 'distal'), JointHalf(tree.Joints[1], 'distal'))
    check(f"short joint (h = 0.6r) at full range: its elbow reaches {elbow:.3f}r into the part before it, "
          f"not counted ({tree.getCollidingPairs()})", elbow > 0.1 and tree.getCollidingPairs() == [])

def testRandomAgainstSolids(trials=60):
    # Random adjacent pairs and triples of joints with ranges of at least +-90 degrees (h >= r, where a bent
    # joint's elbow stays within the joint), random twists between their axes and random states: the checker
    # must agree with the sampled solids, except within a band around touching that sampling can't resolve.
    # Half the trials bend every joint far toward the same side with small twists, where overlaps happen.
    disagreements, overlapping = [], 0
    for trial in range(trials):
        n = 2 if trial % 2 == 0 else 3
        hs = list(RNG.uniform(1, 3, n))
        if n == 3:
            hs[1] = RNG.uniform(1, 1.5)
        if trial % 4 < 2:
            twists, fractions = RNG.uniform(-np.pi, np.pi, n), RNG.uniform(-1, 1, n)
        else:
            twists, fractions = RNG.uniform(-0.3, 0.3, n), RNG.uniform(0.7, 1, n)
        tree = bent(adjacentRevolutes(hs, [0] + list(twists[1:])), fractions)
        found = checkerPairs(tree)
        definite, possible = truthPairs(tree, 0.01), truthPairs(tree, -0.01)
        overlapping += bool(definite)
        if not (set(definite) <= set(found) <= set(possible)):
            disagreements.append((np.round(hs, 2), found, definite, possible))
    check(f"{trials} random adjacent pairs and triples ({overlapping} with solids overlapping by over 0.01r): "
          f"checker agrees with the solids ({disagreements})", not disagreements)

def testUTurnTangency():
    # A link that goes straight for 3, U-turns and ends at a joint heading back alongside it, exactly 2r from
    # the link's straight part (beyond the window): they only touch, no collision. The CSC path solver places
    # the straight part about 7e-6 r off, which CONTACT_TOLERANCE (1e-4 r) absorbs. (A tighter U-turn can't be
    # built: links bend with radius r.)
    tree = KinematicTree(TransverseRevolute(r=r, Pose=SE3(), minAngle=-np.pi/2, maxAngle=np.pi/2))
    place(tree, 0, TransverseRevolute(r=r, Pose=SE3.Trans(4 - 1, 2 * r, 0) @ SE3.Rz(np.pi),
                                      minAngle=-np.pi/2, maxAngle=np.pi/2))
    pairs = tree.getCollidingPairs()
    check(f"U-turn: joint 1 exactly 2r from the link's straight part: touching, no collision ({pairs})", pairs == [])

def testBranchSiblings():
    # A joint with two children: one directly at its distal end (empty link), the other far away. The link
    # to the far one leaves in the same direction, through the near child, which is a collision (the window
    # doesn't pass between branches). Moving the near child away clears it.
    for gap, expect in [(0.0, True), (6.0, False)]:
        tree = KinematicTree(TransverseRevolute(r=r, Pose=SE3(), minAngle=-np.pi/2, maxAngle=np.pi/2))
        place(tree, 0, TransverseRevolute(r=r, Pose=SE3.Trans(2 + gap, 0, 0), minAngle=-np.pi/2, maxAngle=np.pi/2))
        place(tree, 0, TransverseRevolute(r=r, Pose=SE3.Trans(3, 8, 0) @ SE3.Rz(np.pi/2), minAngle=-np.pi/2, maxAngle=np.pi/2))
        pairs = tree.getCollidingPairs()
        check(f"branch, near child {gap} from the branch point: "
              f"{'collides with the sibling link' if expect else 'clear'} ({pairs})",
              (((1, 'Joint'), (2, 'Link')) in pairs) if expect else pairs == [])

def testSavedTrees():
    # Hidden overlaps found by the closest-pair masking probe (2026-11-26): a link swinging back past the
    # proximal end of the joint it leaves, missed before the halves and the path window
    base = "Trials Before Experiments"
    for sweep, seed, label, c in [("2026.11.23_19.01.14_Sweep_Joints4_Seeds0-7", 5, "box_straight_align10_power2", 1),
                                  ("2026.11.25_04.23.23_Sweep_Joints6_Seeds0-2", 1, "box_straight_align10_power2", 3)]:
        path = os.path.join(base, sweep, f"seed{seed}", label, "final.tree")
        if not os.path.exists(path):
            print(f"- skipped (not here): {path}")
            continue
        tree = loadKinematicTree(path)
        config = json.load(open(os.path.join(base, sweep, f"seed{seed}", "configs.json")))[c]
        for i in range(len(tree.Joints)):
            if not isWaypoint(tree.Joints[i]):
                tree.setJointState(i, config[i])
        pairs = tree.getCollidingPairs()
        check(f"saved tree {sweep[:19]} seed {seed} {label} config {c}: joint 0 vs link 2 ({pairs})",
              pairs == [((0, 'Joint'), (2, 'Link'))])

def testCylinderSdf():
    # sdf_cylinder (the joint halves' sdf, and links' straight ends) against the distance to densely sampled
    # surface points, including beside the rims, where the former capsule-and-planes construction underestimated it
    from scipy.spatial import cKDTree
    a, b, rad = np.array([0.3, -0.2, 0.5]), np.array([1.3, 1.8, 0.0]), 0.8
    L = np.linalg.norm(b - a); n = (b - a) / L
    u = np.cross(n, [0, 0, 1.0]); u /= np.linalg.norm(u); v = np.cross(n, u)
    phi = np.linspace(0, 2 * np.pi, 360, endpoint=False)[None, :, None]
    ring = np.cos(phi) * u + np.sin(phi) * v
    side = (a + np.linspace(0, 1, 400)[:, None, None] * (b - a) + rad * ring).reshape(-1, 3)
    discs = [(c + np.linspace(0, rad, 120)[:, None, None] * ring).reshape(-1, 3) for c in (a, b)]
    surface = cKDTree(np.vstack([side] + discs))
    P = a + RNG.uniform(-2, 2, (3000, 3)) + (b - a) * RNG.uniform(-0.5, 1.5, (3000, 1))
    axial = (P - a) @ n; radial = np.linalg.norm((P - a) - axial[:, None] * n, axis=1)
    inside = (axial >= 0) & (axial <= L) & (radial <= rad)
    truth = np.where(inside, -1, 1) * surface.query(P)[0]
    err = np.abs(sdf_cylinder(np, P, a, b, rad) - truth).max()
    rim = sdf_cylinder(np, (a - 0.765 * n + (rad + 0.78) * u)[None], a, b, rad)[0]
    check(f"cylinder sdf: within {err:.4f} of sampled surface distances (spacing 0.014); beside a rim "
          f"{rim:.4f} (true {np.hypot(0.765, 0.78):.4f})", err < 0.01 and abs(rim - np.hypot(0.765, 0.78)) < 1e-9)

def testRimFalsePositive():
    # A link passing beside the rim of a joint's flat end, clear of it by 0.035r in the sampled solids, which the
    # former joint-half sdf counted as a collision (scan of saved trees, 2026-11-27)
    path = os.path.join("Trials Before Experiments", "2026.11.22_02.00.54_Sweep_Joints4_Seeds0-7", "seed2",
                        "multiscale_straight_align10_power2", "224.2343807220459_2.tree")
    if not os.path.exists(path):
        print(f"- skipped (not here): {path}")
        return
    tree = loadKinematicTree(path)
    config = json.load(open(os.path.join(os.path.dirname(os.path.dirname(path)), "configs.json")))[4]
    for i in range(len(tree.Joints)):
        if not isWaypoint(tree.Joints[i]):
            tree.setJointState(i, config[i])
    pairs = tree.getCollidingPairs()
    check(f"link 2 beside the rim of joint 0's proximal end, 0.035r clear: no collision ({pairs})",
          ((0, 'Joint'), (2, 'Link')) not in pairs)

if __name__ == "__main__":
    testFixedCases()
    testShortJointElbow()
    testRandomAgainstSolids()
    testUTurnTangency()
    testBranchSiblings()
    testSavedTrees()
    testCylinderSdf()
    testRimFalsePositive()
    print("\n✓ All tests passed!")
