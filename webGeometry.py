"""
webGeometry.py — geometry extraction layer for the web port.

Returns plain numpy arrays (no pyqtgraph) in the format:
    {
        "meshes": [{"vertices": ndarray(N,3,f32), "faces": ndarray(M,3,i32),
                    "color": (r,g,b,a), "label": str}, ...],
        "lines":  [{"points": ndarray(N,3,f32), "color": (r,g,b,a),
                    "label": str}, ...]
    }

All units match the source data (mm for PrintedKinematicTree).
"""
from __future__ import annotations
import numpy as np
from numpy.linalg import norm

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _empty():
    return {"meshes": [], "lines": []}


def _merge(*geos):
    out = _empty()
    for g in geos:
        out["meshes"].extend(g["meshes"])
        out["lines"].extend(g["lines"])
    return out


def _rgba(color):
    """Normalise to a plain (r,g,b,a) tuple of Python floats."""
    return tuple(float(c) for c in color)


def _mesh(vertices, faces, color, label=""):
    return {
        "vertices": np.asarray(vertices, dtype=np.float32),
        "faces":    np.asarray(faces,    dtype=np.int32),
        "color":    _rgba(color),
        "label":    label,
    }


def _line(points, color, label=""):
    return {
        "points": np.asarray(points, dtype=np.float32),
        "color":  _rgba(color),
        "label":  label,
    }


# ---------------------------------------------------------------------------
# Low-level primitives (no external GL dependencies)
# ---------------------------------------------------------------------------

def _sphere_verts_faces(center, radius, rows=10, cols=10):
    """UV sphere matching pyqtgraph MeshData.sphere(rows, cols)."""
    center = np.asarray(center, dtype=np.float64)
    verts = []
    faces = []

    for r in range(rows + 1):
        theta = np.pi * r / rows          # 0 … pi
        sin_t = np.sin(theta)
        cos_t = np.cos(theta)
        for c in range(cols):
            phi = 2 * np.pi * c / cols    # 0 … 2pi (endpoint=False)
            x = sin_t * np.cos(phi)
            y = sin_t * np.sin(phi)
            z = cos_t
            verts.append(center + radius * np.array([x, y, z]))

    # quads → 2 triangles each
    for r in range(rows):
        for c in range(cols):
            c_next = (c + 1) % cols
            v00 = r * cols + c
            v01 = r * cols + c_next
            v10 = (r + 1) * cols + c
            v11 = (r + 1) * cols + c_next
            faces.append([v00, v01, v11])
            faces.append([v00, v11, v10])

    return np.array(verts, dtype=np.float32), np.array(faces, dtype=np.int32)


def _circle_disc_verts_faces(center, normal, radius, n=32):
    """Filled disc (for cylinder end caps)."""
    normal = np.asarray(normal, dtype=np.float64)
    normal = normal / norm(normal)
    if abs(normal[0]) < 0.9:
        u = np.cross(normal, [1.0, 0.0, 0.0])
    else:
        u = np.cross(normal, [0.0, 1.0, 0.0])
    u = u / norm(u)
    v = np.cross(normal, u)

    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ring = (np.asarray(center) +
            radius * (np.cos(angles)[:, None] * u + np.sin(angles)[:, None] * v))
    verts = np.vstack([np.asarray(center).reshape(1, 3), ring])
    tris = np.array([[0, i + 1, (i % n) + 2] for i in range(n - 1)] + [[0, n, 1]])
    return verts.astype(np.float32), tris.astype(np.int32)


def _box_verts_faces(box_center, xhat, yhat, zhat, bx, by, bz):
    """Axis-aligned box (arbitrary orientation via xhat/yhat/zhat)."""
    hx, hy, hz = bx / 2, by / 2, bz / 2
    local = np.array([
        [-hx, -hy, -hz], [+hx, -hy, -hz], [+hx, +hy, -hz], [-hx, +hy, -hz],
        [-hx, -hy, +hz], [+hx, -hy, +hz], [+hx, +hy, +hz], [-hx, +hy, +hz],
    ])
    R = np.column_stack([xhat, yhat, zhat])
    corners = (R @ local.T).T + np.asarray(box_center)
    box_tris = np.array([
        [0,1,2], [0,2,3],   # -z face
        [4,5,6], [4,6,7],   # +z face
        [0,1,5], [0,5,4],   # -y face
        [2,3,7], [2,7,6],   # +y face
        [0,3,7], [0,7,4],   # -x face
        [1,2,6], [1,6,5],   # +x face
    ], dtype=np.int32)
    return corners.astype(np.float32), box_tris


# ---------------------------------------------------------------------------
# Per-primitive geometry extractors
# ---------------------------------------------------------------------------

def getBallGeometry(ball, color):
    """Ball → one mesh entry."""
    verts, faces = _sphere_verts_faces(ball.c, ball.r, rows=10, cols=10)
    return _merge({"meshes": [_mesh(verts, faces, color, "ball")], "lines": []})


def getCylinderGeometry(cylinder, color, numSides=16, numCircles=2):
    """Cylinder → one mesh entry (no end caps)."""
    verts, faces = cylinder.interpolateQtCircles(numSides, numCircles)
    return {"meshes": [_mesh(verts, faces, color, "cylinder")], "lines": []}


def getElbowGeometry(elbow, color, numSides=16):
    """Single Elbow → one mesh entry."""
    verts, faces = elbow.circleEllipseCircleQT(numSides)
    return {"meshes": [_mesh(verts, faces, color, "elbow")], "lines": []}


def getCompoundElbowGeometry(compound_elbow, color, numSides=16):
    """CompoundElbow (list of elbows with vertex offsets) → one mesh entry."""
    verts_list, faces_list = compound_elbow.circleEllipseCircleQT(numSides)
    if len(verts_list) == 0:
        return _empty()
    verts = np.array(verts_list, dtype=np.float32)
    faces = np.array(faces_list, dtype=np.int32)
    return {"meshes": [_mesh(verts, faces, color, "compound_elbow")], "lines": []}


def getCircle3DGeometry(circle3d, color, count=64):
    """Circle3D → one line entry."""
    pts = circle3d.interpolate(count=count)
    return {"meshes": [], "lines": [_line(pts, color, "circle3d")]}


# ---------------------------------------------------------------------------
# Link geometry
# ---------------------------------------------------------------------------

def getLinkGeometry(link, color, numSides=8):
    """
    LinkCSC/PrintedLinkCSC → combined mesh of elbow1 + cylinder + elbow2.
    Mirrors LinkCSC.addToWidget logic exactly (including vertex offset handling).
    """
    DISTANCE_EPSILON = getattr(link, 'DISTANCE_EPSILON', 1e-4)

    all_verts = []
    all_faces = []
    v_offset = 0

    if link.elbow1 is not None:
        v, f = link.elbow1.circleEllipseCircleQT(numSides)
        v = np.asarray(v, dtype=np.float32)
        f = np.asarray(f, dtype=np.int32) + v_offset
        all_verts.append(v)
        all_faces.append(f)
        v_offset += len(v)

    if link.path.tMag > DISTANCE_EPSILON:
        v, f = link.cylinder.interpolateQtCircles(numSides, 2)
        v = np.asarray(v, dtype=np.float32)
        f = np.asarray(f, dtype=np.int32) + v_offset
        all_verts.append(v)
        all_faces.append(f)
        v_offset += len(v)

    if link.elbow2 is not None:
        v, f = link.elbow2.circleEllipseCircleQT(numSides)
        v = np.asarray(v, dtype=np.float32)
        f = np.asarray(f, dtype=np.int32) + v_offset
        all_verts.append(v)
        all_faces.append(f)
        v_offset += len(v)

    if not all_verts:
        return _empty()

    verts = np.vstack(all_verts)
    faces = np.vstack(all_faces)
    return {"meshes": [_mesh(verts, faces, color, "link")], "lines": []}


# ---------------------------------------------------------------------------
# Joint geometry helpers
# ---------------------------------------------------------------------------

def _pose_frame_lines(pose, axis_colors, axis_scale):
    """
    Returns 3 line entries for x/y/z axes of a pose frame.
    axis_colors: (xColor, yColor, zColor)
    """
    lines = []
    for i, col in enumerate(axis_colors):
        start = pose.t
        end = start + axis_scale * pose.R[:, i]
        lines.append(_line(np.array([start, end]), col, f"axis_{i}"))
    return lines


def _joint_axis_line(joint, color):
    """Long axis line through joint center (±10r along z)."""
    zhat = joint.Pose.R[:, 2]
    pts = np.array([joint.Pose.t - 10 * joint.r * zhat,
                    joint.Pose.t + 10 * joint.r * zhat])
    return _line(pts, color, "joint_axis")


def _tip_hemisphere_mesh(tip, color, n_lat=6, n_lon=8):
    """
    Stretched hemisphere mesh for a Tip joint.
    Mirrors Tip.addToWidget logic exactly.
    """
    if tip.forward:
        theta = np.linspace(0, np.pi / 2, n_lat)
    else:
        theta = np.linspace(np.pi / 2, np.pi, n_lat)
    phi = np.linspace(0, 2 * np.pi, n_lon, endpoint=False)

    vertices = []
    for i in range(n_lat):
        for j in range(n_lon):
            x_s = np.sin(theta[i]) * np.cos(phi[j])
            y_s = np.sin(theta[i]) * np.sin(phi[j])
            z_s = np.cos(theta[i])
            vertices.append([tip.r * x_s, tip.r * y_s, tip.neutralLength * z_s])

    # Pole vertex
    if tip.forward:
        vertices.append([0.0, 0.0, tip.neutralLength])
    else:
        vertices.append([0.0, 0.0, -tip.neutralLength])

    pole_idx = n_lat * n_lon
    vertices = np.array(vertices, dtype=np.float32)

    faces = []
    for i in range(n_lat - 1):
        for j in range(n_lon):
            j_next = (j + 1) % n_lon
            v00 = i * n_lon + j
            v01 = i * n_lon + j_next
            v10 = (i + 1) * n_lon + j
            v11 = (i + 1) * n_lon + j_next
            faces.append([v00, v01, v11])
            faces.append([v00, v11, v10])

    if tip.forward:
        for j in range(n_lon):
            j_next = (j + 1) % n_lon
            faces.append([pole_idx, j_next, j])
    else:
        last_row = (n_lat - 1) * n_lon
        for j in range(n_lon):
            j_next = (j + 1) % n_lon
            faces.append([pole_idx, last_row + j, last_row + j_next])

    faces = np.array(faces, dtype=np.int32)

    # Transform to world (mirrors Tip.addToWidget)
    if tip.forward:
        base_frame = tip.ProximalFrame()
    else:
        base_frame = tip.DistalFrame()
    R = base_frame.R
    x_sc, y_sc, z_sc = vertices[:, 0], vertices[:, 1], vertices[:, 2]
    pidx = tip.pidx
    if pidx == 0:
        xw = base_frame.t[0] + z_sc*R[0,0] + x_sc*R[0,1] + y_sc*R[0,2]
        yw = base_frame.t[1] + z_sc*R[1,0] + x_sc*R[1,1] + y_sc*R[1,2]
        zw = base_frame.t[2] + z_sc*R[2,0] + x_sc*R[2,1] + y_sc*R[2,2]
    elif pidx == 1:
        xw = base_frame.t[0] + x_sc*R[0,0] + z_sc*R[0,1] + y_sc*R[0,2]
        yw = base_frame.t[1] + x_sc*R[1,0] + z_sc*R[1,1] + y_sc*R[1,2]
        zw = base_frame.t[2] + x_sc*R[2,0] + z_sc*R[2,1] + y_sc*R[2,2]
    else:
        xw = base_frame.t[0] + x_sc*R[0,0] + y_sc*R[0,1] + z_sc*R[0,2]
        yw = base_frame.t[1] + x_sc*R[1,0] + y_sc*R[1,1] + z_sc*R[1,2]
        zw = base_frame.t[2] + x_sc*R[2,0] + y_sc*R[2,1] + z_sc*R[2,2]

    world_verts = np.column_stack([xw, yw, zw]).astype(np.float32)
    return _mesh(world_verts, faces, color, "tip")


# ---------------------------------------------------------------------------
# Per-joint-type geometry
# ---------------------------------------------------------------------------

def _getBaseJointLines(joint, xColor, yColor, zColor,
                       proximalColor, centerColor, distalColor,
                       showAxis, showPoses, axisScale,
                       showAxisColor=(0.75, 0.75, 0.75, 1)):
    lines = []
    if showAxis:
        lines.append(_joint_axis_line(joint, showAxisColor))
    if showPoses:
        for pose, col in zip(
                [joint.ProximalFrame(), joint.DistalFrame(), joint.Pose],
                [proximalColor, distalColor, centerColor]):
            lines.extend(_pose_frame_lines(pose, (xColor, yColor, zColor),
                                           axisScale * joint.r))
    return lines


def _getRevoluteGeometry(joint, color, xColor, yColor, zColor,
                         proximalColor, centerColor, distalColor,
                         showSurface, showAxis, showPoses, axisScale,
                         numSides=16):
    """Generic Revolute (not TransverseRDS3225) — proximal + center + distal cylinders + sphere."""
    geo = _empty()
    geo["lines"].extend(_getBaseJointLines(
        joint, xColor, yColor, zColor, proximalColor, centerColor, distalColor,
        showAxis, showPoses, axisScale))

    if showSurface:
        geo = _merge(geo,
                     getCylinderGeometry(joint.proximalCylinder(), color, numSides),
                     getBallGeometry(joint.centerSphere(), color),
                     getCylinderGeometry(joint.distalCylinder(), color, numSides))
    return geo


def _getPrismaticGeometry(joint, color, xColor, yColor, zColor,
                          proximalColor, centerColor, distalColor,
                          showSurface, showAxis, showPoses, axisScale,
                          numSides=16):
    geo = _empty()
    geo["lines"].extend(_getBaseJointLines(
        joint, xColor, yColor, zColor, proximalColor, centerColor, distalColor,
        showAxis, showPoses, axisScale))
    if showSurface:
        geo = _merge(geo,
                     getCylinderGeometry(joint.proximalCylinder(), color, numSides),
                     getCylinderGeometry(joint.distalCylinder(), color, numSides))
    return geo


def _getWaypointGeometry(joint, color, xColor, yColor, zColor,
                         proximalColor, centerColor, distalColor,
                         showSurface, showAxis, showPoses, axisScale):
    from geometryHelpers import Circle3D
    geo = _empty()
    geo["lines"].extend(_getBaseJointLines(
        joint, xColor, yColor, zColor, proximalColor, centerColor, distalColor,
        showAxis, showPoses, axisScale))
    if showSurface:
        circle = Circle3D(joint.r, joint.Pose.t, joint.Pose.R[:, joint.pidx])
        geo = _merge(geo, getCircle3DGeometry(circle, color))
        geo = _merge(geo, getBallGeometry(
            type('_B', (), {'c': joint.Pose.t, 'r': 0.05 * joint.r})(), color))
    # Full-size transparent sphere so the waypoint is raycasting-selectable
    verts, faces = _sphere_verts_faces(joint.Pose.t, joint.r, rows=8, cols=8)
    pick_mesh = _mesh(verts, faces, (0, 0, 0, 0), "pick_sphere")
    pick_mesh["pick_only"] = True
    geo["meshes"].append(pick_mesh)
    return geo


def _getTipGeometry(joint, color, xColor, yColor, zColor,
                    proximalColor, centerColor, distalColor,
                    showSurface, showAxis, showPoses, axisScale):
    geo = _empty()
    geo["lines"].extend(_getBaseJointLines(
        joint, xColor, yColor, zColor, proximalColor, centerColor, distalColor,
        showAxis, showPoses, axisScale))
    if showSurface:
        geo["meshes"].append(_tip_hemisphere_mesh(joint, color))
    return geo


def _getTransverseRDS3225Geometry(joint, color, xColor, yColor, zColor,
                                   proximalColor, centerColor, distalColor,
                                   showSurface, showAxis, showPoses, axisScale,
                                   numCylPoints=32):
    """
    Mirrors TransverseRDS3225.addToWidget: box + 2 cylinders + brackets.
    """
    geo = _empty()
    geo["lines"].extend(_getBaseJointLines(
        joint, xColor, yColor, zColor, proximalColor, centerColor, distalColor,
        showAxis, showPoses, axisScale))

    if not showSurface:
        return geo

    xhat = joint.Pose.R[:, 0]
    yhat = joint.Pose.R[:, 1]
    zhat = joint.Pose.R[:, 2]

    # --- Box ---
    bx, by, bz = 40.0, 20.0, 40.0
    box_center = joint.Pose.t + (-11.0) * xhat
    hx, hy, hz = bx / 2, by / 2, bz / 2

    box_v, box_f = _box_verts_faces(box_center, xhat, yhat, zhat, bx, by, bz)
    geo["meshes"].append(_mesh(box_v, box_f, (0.1, 0.1, 0.1, 0.6), "servo_box"))

    # --- Proximal cylinder + end cap ---
    cyl_length = 9.0
    from geometryHelpers import Cylinder
    proximal_pos = joint.ProximalFrame().t
    proximal_cyl = Cylinder(joint.r, proximal_pos, xhat, cyl_length)
    geo = _merge(geo, getCylinderGeometry(proximal_cyl, color, numCylPoints))

    prox_cap_center = proximal_pos + cyl_length * xhat
    cv, cf = _circle_disc_verts_faces(prox_cap_center, xhat, joint.r, numCylPoints)
    geo["meshes"].append(_mesh(cv, cf, color, "prox_cap"))

    # --- Distal cylinder + end cap ---
    distal_frame = joint.DistalFrame()
    distal_pos = distal_frame.t
    distal_xhat = distal_frame.R[:, 0]
    distal_yhat = distal_frame.R[:, 1]
    distal_zhat = distal_frame.R[:, 2]
    distal_cyl = Cylinder(joint.r, distal_pos, -distal_xhat, cyl_length)
    geo = _merge(geo, getCylinderGeometry(distal_cyl, color, numCylPoints))

    distal_inner = distal_pos - cyl_length * distal_xhat
    dv, df = _circle_disc_verts_faces(distal_inner, -distal_xhat, joint.r, numCylPoints)
    geo["meshes"].append(_mesh(dv, df, color, "dist_cap"))

    # --- Bracket lines (proximal [-bracket) ---
    line_color = (0.3, 0.3, 0.3, 1.0)
    prox_inner = prox_cap_center  # already computed above
    box_prox_face = box_center - hx * xhat

    for sign_y in [+1, -1]:
        y_off = sign_y * 10.0 * yhat
        pts = np.array([
            box_prox_face + hz * zhat + y_off,
            prox_inner    + hz * zhat + y_off,
            prox_inner    - hz * zhat + y_off,
            box_prox_face - hz * zhat + y_off,
        ])
        geo["lines"].append(_line(pts, line_color, "prox_bracket_outline"))

    for z_sign in [+1, -1]:
        for x_pos in [box_prox_face, prox_inner]:
            corner = x_pos + z_sign * hz * zhat
            edge = np.array([corner + 10.0 * yhat, corner - 10.0 * yhat])
            geo["lines"].append(_line(edge, line_color, "prox_bracket_edge"))

    # Proximal bracket face quads
    bracket_color = (0.5, 0.5, 0.5, 0.8)
    prox_quads = [
        [box_prox_face + hz*zhat + 10*yhat, prox_inner + hz*zhat + 10*yhat,
         prox_inner + hz*zhat - 10*yhat,   box_prox_face + hz*zhat - 10*yhat],
        [prox_inner + hz*zhat + 10*yhat,   prox_inner - hz*zhat + 10*yhat,
         prox_inner - hz*zhat - 10*yhat,   prox_inner + hz*zhat - 10*yhat],
        [prox_inner - hz*zhat + 10*yhat,   box_prox_face - hz*zhat + 10*yhat,
         box_prox_face - hz*zhat - 10*yhat, prox_inner - hz*zhat - 10*yhat],
    ]
    for quad in prox_quads:
        v = np.array(quad, dtype=np.float32)
        f = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int32)
        geo["meshes"].append(_mesh(v, f, bracket_color, "prox_bracket_face"))

    # --- Distal [-bracket ---
    for sign_y in [+1, -1]:
        y_off = sign_y * 10.0 * distal_yhat
        pts = np.array([
            joint.Pose.t + hz * zhat + y_off,
            distal_inner  + hz * zhat + y_off,
            distal_inner  - hz * zhat + y_off,
            joint.Pose.t  - hz * zhat + y_off,
        ])
        geo["lines"].append(_line(pts, line_color, "dist_bracket_outline"))

    for z_sign in [+1, -1]:
        for x_pos in [joint.Pose.t, distal_inner]:
            corner = x_pos + z_sign * hz * zhat
            edge = np.array([corner + 10.0 * distal_yhat, corner - 10.0 * distal_yhat])
            geo["lines"].append(_line(edge, line_color, "dist_bracket_edge"))

    dist_quads = [
        [joint.Pose.t + hz*zhat + 10*distal_yhat, distal_inner + hz*zhat + 10*distal_yhat,
         distal_inner + hz*zhat - 10*distal_yhat, joint.Pose.t + hz*zhat - 10*distal_yhat],
        [distal_inner + hz*zhat + 10*distal_yhat, distal_inner - hz*zhat + 10*distal_yhat,
         distal_inner - hz*zhat - 10*distal_yhat, distal_inner + hz*zhat - 10*distal_yhat],
        [distal_inner - hz*zhat + 10*distal_yhat, joint.Pose.t - hz*zhat + 10*distal_yhat,
         joint.Pose.t - hz*zhat - 10*distal_yhat, distal_inner - hz*zhat - 10*distal_yhat],
    ]
    for quad in dist_quads:
        v = np.array(quad, dtype=np.float32)
        f = np.array([[0, 1, 2], [0, 2, 3]], dtype=np.int32)
        geo["meshes"].append(_mesh(v, f, bracket_color, "dist_bracket_face"))

    return geo


def _getCoaxialRDS3225Geometry(joint, color, xColor, yColor, zColor,
                                proximalColor, centerColor, distalColor,
                                showSurface, showAxis, showPoses, axisScale,
                                numCylPoints=32):
    """
    Mirrors CoaxialRDS3225.addToWidget — falls back to generic revolute display
    since CoaxialRDS3225 doesn't override addToWidget in PrintedTube.py.
    """
    return _getRevoluteGeometry(
        joint, color, xColor, yColor, zColor, proximalColor, centerColor, distalColor,
        showSurface, showAxis, showPoses, axisScale, numCylPoints)


# ---------------------------------------------------------------------------
# Public dispatch
# ---------------------------------------------------------------------------

# Default colors (matching style.py)
_REVOLUTE_COLOR  = (0.0, 0.0, 0.545, 1.0)
_LINK_COLOR      = (0.5, 0.5, 0.5,   1.0)
_X_COLOR         = (1.0, 0.27, 0.0,  1.0)
_Y_COLOR         = (0.56, 0.93, 0.56, 1.0)
_Z_COLOR         = (0.0, 0.0, 0.55,  1.0)
_PROXIMAL_COLOR  = (0.0, 1.0, 1.0,   1.0)
_CENTER_COLOR    = (1.0, 0.0, 1.0,   1.0)
_DISTAL_COLOR    = (1.0, 1.0, 0.0,   1.0)
_AXIS_COLOR      = (0.75, 0.75, 0.75, 1.0)
_SELECTED_COLOR  = (1.0, 0.5, 0.0,   1.0)


def getJointGeometry(joint,
                     color=None,
                     xColor=_X_COLOR, yColor=_Y_COLOR, zColor=_Z_COLOR,
                     proximalColor=_PROXIMAL_COLOR, centerColor=_CENTER_COLOR,
                     distalColor=_DISTAL_COLOR,
                     showSurface=True, showAxis=True, showPoses=False,
                     axisScale=1):
    """Dispatch to the correct per-type geometry function."""
    # Import joint types here to avoid circular imports at module load
    from Joint import Revolute, Prismatic, Waypoint, Tip

    if color is None:
        color = _REVOLUTE_COLOR

    # Try most-specific PrintedTube types first
    try:
        from PrintedTube import TransverseRDS3225, CoaxialRDS3225
        if isinstance(joint, TransverseRDS3225):
            return _getTransverseRDS3225Geometry(
                joint, color, xColor, yColor, zColor, proximalColor, centerColor,
                distalColor, showSurface, showAxis, showPoses, axisScale)
        if isinstance(joint, CoaxialRDS3225):
            return _getCoaxialRDS3225Geometry(
                joint, color, xColor, yColor, zColor, proximalColor, centerColor,
                distalColor, showSurface, showAxis, showPoses, axisScale)
    except ImportError:
        pass

    if isinstance(joint, Tip):
        return _getTipGeometry(joint, color, xColor, yColor, zColor,
                               proximalColor, centerColor, distalColor,
                               showSurface, showAxis, showPoses, axisScale)
    if isinstance(joint, Waypoint):
        return _getWaypointGeometry(joint, color, xColor, yColor, zColor,
                                    proximalColor, centerColor, distalColor,
                                    showSurface, showAxis, showPoses, axisScale)
    if isinstance(joint, Revolute):
        return _getRevoluteGeometry(joint, color, xColor, yColor, zColor,
                                    proximalColor, centerColor, distalColor,
                                    showSurface, showAxis, showPoses, axisScale)
    if isinstance(joint, Prismatic):
        return _getPrismaticGeometry(joint, color, xColor, yColor, zColor,
                                     proximalColor, centerColor, distalColor,
                                     showSurface, showAxis, showPoses, axisScale)

    # Unknown joint type: fall back to pose frames only
    geo = _empty()
    geo["lines"].extend(_getBaseJointLines(
        joint, xColor, yColor, zColor, proximalColor, centerColor, distalColor,
        showAxis, showPoses, axisScale))
    return geo


def getTreeGeometry(tree,
                    selectedJoint=-1, selectedLink=-1,
                    jointColor=_REVOLUTE_COLOR, linkColor=_LINK_COLOR,
                    showJointSurface=True, showJointAxis=True, showJointPoses=False,
                    showLinkSurface=True,
                    xColor=_X_COLOR, yColor=_Y_COLOR, zColor=_Z_COLOR,
                    proximalColor=_PROXIMAL_COLOR, centerColor=_CENTER_COLOR,
                    distalColor=_DISTAL_COLOR,
                    numSides=8):
    """
    Top-level entry point. Build the full geometry for a KinematicTree.

    Returns {"meshes": [...], "lines": [...]} with all joints and links.
    """
    geo = _empty()

    for i, joint in enumerate(tree.Joints):
        jcolor = _SELECTED_COLOR if i == selectedJoint else jointColor
        joint_geo = getJointGeometry(
            joint, color=jcolor,
            xColor=xColor, yColor=yColor, zColor=zColor,
            proximalColor=proximalColor, centerColor=centerColor,
            distalColor=distalColor,
            showSurface=showJointSurface, showAxis=showJointAxis,
            showPoses=showJointPoses)
        for item in joint_geo["meshes"] + joint_geo["lines"]:
            item["joint_index"] = i
        geo = _merge(geo, joint_geo)

    for i, link in enumerate(tree.Links):
        lcolor = _SELECTED_COLOR if i == selectedLink else linkColor
        # Skip the degenerate "root" link (index 0 connects root to itself)
        if i == 0 and tree.Parents[0] == -1:
            continue
        link_geo = getLinkGeometry(link, lcolor, numSides)
        for item in link_geo["meshes"] + link_geo["lines"]:
            item["link_index"] = i
        geo = _merge(geo, link_geo)

    # Per-joint metadata so the client can position gizmos
    from scipy.spatial.transform import Rotation as _Rot
    geo["joints"] = [
        {
            "index": i,
            "type":  type(joint).__name__,
            "position":   joint.Pose.t.tolist(),
            "quaternion": _Rot.from_matrix(joint.Pose.R).as_quat().tolist(),
        }
        for i, joint in enumerate(tree.Joints)
    ]

    return geo


# ---------------------------------------------------------------------------
# JSON serialisation helper
# ---------------------------------------------------------------------------

def geometryToJson(geo):
    """
    Convert numpy arrays in a geometry dict to plain Python lists,
    ready for json.dumps().
    """
    def arr_to_list(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        return x

    return {
        "meshes": [
            {k: arr_to_list(v) for k, v in m.items()}
            for m in geo["meshes"]
        ],
        "lines": [
            {k: arr_to_list(v) for k, v in ln.items()}
            for ln in geo["lines"]
        ],
        "joints": geo.get("joints", []),
    }
