"""
test_webGeometry.py — visual smoke test for webGeometry.py.

Builds a small PrintedKinematicTree and renders the extracted geometry
using matplotlib's Poly3DCollection (no pyqtgraph required at render time).

Run from the KinegamiPython directory:
    python test_webGeometry.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from spatialmath import SE3

# ---- Build tree ----
from PrintedTube import TransverseRDS3225, PrintedEndHemisphere, PrintedKinematicTree

root_pose = SE3.Ry(-np.pi / 2)
root = TransverseRDS3225(root_pose, version=270)

tree = PrintedKinematicTree(root)

# Add a second servo — addJoint returns the actual index of the new joint
# (may be >1 because intermediate Waypoints are inserted automatically)
second = TransverseRDS3225(SE3(), version=270)
second_idx = tree.addJoint(0, second)

# Add an end hemisphere after the second servo
tip = PrintedEndHemisphere(TransverseRDS3225.R, SE3(), pathIndex=0)
tree.addJoint(second_idx, tip)

print(f"Tree: {len(tree.Joints)} joints, {len(tree.Links)} links")
print(f"  second servo at index {second_idx}")

# ---- Extract geometry ----
from webGeometry import getTreeGeometry

geo = getTreeGeometry(tree,
                      showJointSurface=True,
                      showJointAxis=True,
                      showJointPoses=True,
                      showLinkSurface=True,
                      numSides=12)

print(f"Meshes: {len(geo['meshes'])}")
print(f"Lines:  {len(geo['lines'])}")
for m in geo["meshes"]:
    print(f"  [{m['label']:25s}]  verts={m['vertices'].shape}  faces={m['faces'].shape}  color={m['color']}")

# ---- Plot ----
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel('X (mm)')
ax.set_ylabel('Y (mm)')
ax.set_zlabel('Z (mm)')
ax.set_title('webGeometry smoke test')

# Draw meshes
for m in geo["meshes"]:
    verts = m["vertices"]
    faces = m["faces"]
    rgba  = m["color"]
    tri_verts = verts[faces]      # shape (F, 3, 3)
    poly = Poly3DCollection(tri_verts, alpha=rgba[3], linewidths=0)
    poly.set_facecolor(rgba[:3])
    ax.add_collection3d(poly)

# Draw lines
for ln in geo["lines"]:
    pts = ln["points"]
    col = ln["color"]
    ax.plot(pts[:, 0], pts[:, 1], pts[:, 2],
            color=col[:3], alpha=col[3], linewidth=1.2)


def set_equal_axes(ax):
    """Force equal aspect ratio in 3D."""
    pts = []
    for m in geo["meshes"]:
        pts.append(m["vertices"])
    for ln in geo["lines"]:
        pts.append(ln["points"])
    if not pts:
        return
    all_pts = np.vstack(pts)
    mins = all_pts.min(axis=0)
    maxs = all_pts.max(axis=0)
    center = (mins + maxs) / 2
    half_range = (maxs - mins).max() / 2
    ax.set_xlim(center[0] - half_range, center[0] + half_range)
    ax.set_ylim(center[1] - half_range, center[1] + half_range)
    ax.set_zlim(center[2] - half_range, center[2] + half_range)


set_equal_axes(ax)
plt.tight_layout()
plt.show()
