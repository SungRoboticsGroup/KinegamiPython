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


def flat_rectangle(cx, cy, sx, sy, z=0.0,
                   color=(0.8, 0.8, 0.8, 1.0), glOptions='opaque'):
    """
    Flat horizontal rectangle centered at (cx, cy) on the plane at height z,
    with side length sx along x and sy along y.
    """
    hx, hy = sx / 2.0, sy / 2.0
    verts = np.array([
        [cx - hx, cy - hy, z],
        [cx + hx, cy - hy, z],
        [cx + hx, cy + hy, z],
        [cx - hx, cy + hy, z],
    ], dtype=float)
    faces = np.array([[0, 1, 2], [0, 2, 3]], dtype=int)
    meshdata = gl.MeshData(vertexes=verts, faces=faces)
    # shader=None: flat, unlit color. The 'shaded' shader leaves an upward-
    # facing plane at ~20% brightness, so it would render very dark from above.
    return gl.GLMeshItem(meshdata=meshdata, smooth=False, shader=None,
                         color=color, glOptions=glOptions)


def solid_sphere(cx, cy, cz, radius,
                 color=(0.2, 0.8, 0.4, 0.8), rows=20, cols=20):
    """Solid sphere centered at (cx, cy, cz)."""
    meshdata = gl.MeshData.sphere(rows=rows, cols=cols, radius=radius)
    item = gl.GLMeshItem(meshdata=meshdata, smooth=True, shader='shaded',
                         color=color)
    item.translate(cx, cy, cz)
    return item


def solid_cylinder(cx, cy, cz, radius, length,
                   color=(0.8, 0.4, 0.2, 0.8), cols=20, caps=True,
                   glOptions='opaque'):
    """
    Solid cylinder with its base centered at (cx, cy, cz), extending along +Z.
    caps: if True, close the top and bottom with flat disks.
    glOptions: pass 'translucent' for the color's alpha to take effect.
    Use item.rotate(...) after calling this if you need a different orientation.
    """
    theta = np.linspace(0, 2 * np.pi, cols, endpoint=False)
    ring = np.column_stack([radius * np.cos(theta),
                            radius * np.sin(theta),
                            np.zeros(cols)])
    top = ring + [0, 0, length]
    idx = np.arange(cols)
    nxt = (idx + 1) % cols

    # Side wall: bottom ring = 0..cols-1, top ring = cols..2cols-1
    verts = [ring, top]
    faces = [np.column_stack([idx, nxt, cols + nxt]),
             np.column_stack([idx, cols + nxt, cols + idx])]

    if caps:
        # Caps get their own vertices so smooth shading doesn't blend
        # the flat cap normals with the curved side normals.
        b0 = 2 * cols                 # bottom cap ring start
        bc = 3 * cols                 # bottom cap center
        t0 = 3 * cols + 1             # top cap ring start
        tc = 4 * cols + 1             # top cap center
        verts += [ring, [[0, 0, 0]], top, [[0, 0, length]]]
        faces += [np.column_stack([np.full(cols, bc), b0 + nxt, b0 + idx]),
                  np.column_stack([np.full(cols, tc), t0 + idx, t0 + nxt])]

    meshdata = gl.MeshData(vertexes=np.vstack(verts).astype(float),
                           faces=np.vstack(faces).astype(int))
    item = gl.GLMeshItem(meshdata=meshdata, smooth=True, shader='shaded',
                         color=color, glOptions=glOptions)
    item.translate(cx, cy, cz)
    return item
