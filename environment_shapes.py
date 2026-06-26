"""
Reusable shape primitives for environment files.
Each function returns a single pyqtgraph.opengl item ready to add to the scene.

Built-in pyqtgraph factories used where available:
  gl.MeshData.sphere(rows, cols, radius)
  gl.MeshData.cylinder(rows, cols, radius, length)
Box geometry is computed with numpy (no built-in exists).
"""

import pyqtgraph.opengl as gl
import numpy as np


def circle_outline(cx, cy, radius, z=0.0,
                   color=(1.0, 0.8, 0.2, 1.0), width=2.0, n=128):
    """Closed circle outline on a horizontal plane at height z."""
    theta = np.linspace(0, 2 * np.pi, n, endpoint=True)
    pts = np.zeros((n, 3))
    pts[:, 0] = cx + radius * np.cos(theta)
    pts[:, 1] = cy + radius * np.sin(theta)
    pts[:, 2] = z
    return gl.GLLinePlotItem(pos=pts, color=color, width=width,
                             antialias=True, mode='line_strip')


def solid_box(cx, cy, cz, size,
              color=(0.2, 0.6, 1.0, 0.8)):
    """
    Solid axis-aligned box centered at (cx, cy, cz).
    size: scalar for a cube, or (sx, sy, sz) for a rectangular box.
    """
    if np.isscalar(size):
        hx = hy = hz = size / 2.0
    else:
        hx, hy, hz = size[0] / 2.0, size[1] / 2.0, size[2] / 2.0

    verts = np.array([
        [cx - hx, cy - hy, cz - hz],
        [cx + hx, cy - hy, cz - hz],
        [cx + hx, cy + hy, cz - hz],
        [cx - hx, cy + hy, cz - hz],
        [cx - hx, cy - hy, cz + hz],
        [cx + hx, cy - hy, cz + hz],
        [cx + hx, cy + hy, cz + hz],
        [cx - hx, cy + hy, cz + hz],
    ], dtype=float)

    faces = np.array([
        [0, 1, 2], [0, 2, 3],  # bottom
        [4, 5, 6], [4, 6, 7],  # top
        [0, 1, 5], [0, 5, 4],  # front
        [2, 3, 7], [2, 7, 6],  # back
        [1, 2, 6], [1, 6, 5],  # right
        [0, 3, 7], [0, 7, 4],  # left
    ], dtype=int)

    meshdata = gl.MeshData(vertexes=verts, faces=faces)
    return gl.GLMeshItem(meshdata=meshdata, smooth=False, shader='shaded',
                         color=color)


def solid_sphere(cx, cy, cz, radius,
                 color=(0.2, 0.8, 0.4, 0.8), rows=20, cols=20):
    """Solid sphere centered at (cx, cy, cz)."""
    meshdata = gl.MeshData.sphere(rows=rows, cols=cols, radius=radius)
    item = gl.GLMeshItem(meshdata=meshdata, smooth=True, shader='shaded',
                         color=color)
    item.translate(cx, cy, cz)
    return item


def solid_cylinder(cx, cy, cz, radius, length,
                   color=(0.8, 0.4, 0.2, 0.8), rows=2, cols=20):
    """
    Solid cylinder with its base centered at (cx, cy, cz), extending along +Z.
    Use item.rotate(...) after calling this if you need a different orientation.
    """
    meshdata = gl.MeshData.cylinder(rows=rows, cols=cols,
                                    radius=[radius, radius], length=length)
    item = gl.GLMeshItem(meshdata=meshdata, smooth=True, shader='shaded',
                         color=color)
    item.translate(cx, cy, cz)
    return item
