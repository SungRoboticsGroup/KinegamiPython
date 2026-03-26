"""
web_server.py — Kinegami web backend (FastAPI)

Run with:
    python web_server.py

Then open http://localhost:8000 in a browser.
"""
import os
import numpy as np
from fastapi import FastAPI
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from typing import List

app = FastAPI()
STATIC_DIR = os.path.join(os.path.dirname(__file__), "web_static")


# ── Persistent tree state ─────────────────────────────────────────────────────

_tree = None

def get_tree():
    """Build the demo tree once and reuse it across requests."""
    global _tree
    if _tree is None:
        from spatialmath import SE3
        from PrintedTube import (TransverseRDS3225, PrintedEndHemisphere,
                                  PrintedKinematicTree)
        root = TransverseRDS3225(SE3.Ry(-np.pi / 2), version=270)
        _tree = PrintedKinematicTree(root)
        second_idx = _tree.addJoint(0, TransverseRDS3225(SE3(), version=270))
        _tree.addJoint(second_idx,
                       PrintedEndHemisphere(TransverseRDS3225.R, SE3(), pathIndex=0))
    return _tree


# ── Geometry endpoint ─────────────────────────────────────────────────────────

@app.get("/geometry")
def get_geometry():
    from webGeometry import getTreeGeometry, geometryToJson
    return geometryToJson(getTreeGeometry(get_tree()))


# ── Warm-start link rebuilder ─────────────────────────────────────────────────

def _rebuild_link_warm(tree, link_idx, start_pose, end_pose):
    """
    Rebuild tree.Links[link_idx] using the existing path as a warm-start guess.
    Falls back to a cold solve if the warm-start produces an invalid path.
    """
    import warnings
    from scipy.optimize import fsolve as sp_fsolve
    from numpy.linalg import norm as np_norm

    mk_link = tree._get_link_constructor()
    r       = tree.r
    max_ang = tree.maxAnglePerElbow

    prev = getattr(tree.Links[link_idx], 'path', None)
    if prev is not None:
        try:
            from PathCSC import pathErrorCSC, PathCSC
            x0 = np.append(prev.tUnit, prev.tMag)
            sp = start_pose.t;  sd = start_pose.R[:, 0]
            ep = end_pose.t;    ed = end_pose.R[:, 0]
            with warnings.catch_warnings(action="ignore"):
                sol = sp_fsolve(pathErrorCSC, x0=x0,
                                args=(r, sp, sd, ep, ed,
                                      prev.circle1sign, prev.circle2sign))
            candidate = PathCSC(sol, r, sp, sd, ep, ed,
                                prev.circle1sign, prev.circle2sign)
            eps = 0.01
            if (np_norm(candidate.error) <= 0.005 * r and
                    candidate.theta1 >= -eps and candidate.theta1 < np.pi and
                    candidate.theta2 >= -eps and candidate.theta2 < np.pi):
                return mk_link(r, start_pose, end_pose, max_ang, path=candidate)
        except Exception:
            pass

    return mk_link(r, start_pose, end_pose, max_ang)


# ── Move-joint endpoint ───────────────────────────────────────────────────────

class MoveJointRequest(BaseModel):
    joint_index: int
    position:   List[float]   # world-space [x, y, z] in mm
    quaternion: List[float]   # [x, y, z, w]  (scipy / Three.js convention)


@app.post("/move_joint")
def move_joint(data: MoveJointRequest):
    import copy
    from scipy.spatial.transform import Rotation
    from spatialmath import SE3
    from webGeometry import getTreeGeometry, geometryToJson

    tree  = get_tree()
    joint = tree.Joints[data.joint_index]

    # Save previous pose so we can revert if any link rebuild fails
    old_pose           = copy.deepcopy(joint.Pose)
    old_proximal_dubins = getattr(joint, 'proximalDubins', None)
    old_distal_dubins   = getattr(joint, 'distalDubins',   None)

    # Apply new world pose
    R = Rotation.from_quat(data.quaternion).as_matrix()
    T = np.eye(4)
    T[:3, :3] = R
    T[:3,  3] = data.position
    joint.Pose = SE3(T)
    joint.proximalDubins = joint.ProximalDubinsFrame()
    joint.distalDubins   = joint.DistalDubinsFrame()

    # Attempt to rebuild all affected links (transactional: revert on failure)
    new_links = {}
    failure   = None

    parent_idx = tree.Parents[data.joint_index]
    if parent_idx >= 0:
        parent = tree.Joints[parent_idx]
        try:
            new_links[data.joint_index] = _rebuild_link_warm(
                tree, data.joint_index,
                parent.DistalDubinsFrame(), joint.ProximalDubinsFrame())
        except Exception as e:
            failure = str(e)

    if failure is None:
        for child_idx, p in enumerate(tree.Parents):
            if p != data.joint_index:
                continue
            child = tree.Joints[child_idx]
            try:
                new_links[child_idx] = _rebuild_link_warm(
                    tree, child_idx,
                    joint.DistalDubinsFrame(), child.ProximalDubinsFrame())
            except Exception as e:
                failure = str(e)
                break

    if failure is not None:
        # Revert joint pose — server state is unchanged
        joint.Pose = old_pose
        if old_proximal_dubins is not None:
            joint.proximalDubins = old_proximal_dubins
        if old_distal_dubins is not None:
            joint.distalDubins = old_distal_dubins
        print(f"[move_joint] rejected (link rebuild failed): {failure}")
        return JSONResponse(status_code=422, content={"error": failure})

    # All links rebuilt — commit
    for idx, link in new_links.items():
        tree.Links[idx] = link

    return geometryToJson(getTreeGeometry(tree))


# ── Live link-update endpoint (called every frame during drag) ────────────────

@app.post("/update_links")
def update_links(data: MoveJointRequest):
    from scipy.spatial.transform import Rotation
    from spatialmath import SE3
    from webGeometry import getAffectedLinksGeometry, geometryToJson

    tree  = get_tree()
    joint = tree.Joints[data.joint_index]

    R = Rotation.from_quat(data.quaternion).as_matrix()
    T = np.eye(4)
    T[:3, :3] = R
    T[:3,  3] = data.position
    joint.Pose = SE3(T)
    joint.proximalDubins = joint.ProximalDubinsFrame()
    joint.distalDubins   = joint.DistalDubinsFrame()

    affected = []

    parent_idx = tree.Parents[data.joint_index]
    if parent_idx >= 0:
        parent = tree.Joints[parent_idx]
        try:
            tree.Links[data.joint_index] = _rebuild_link_warm(
                tree, data.joint_index,
                parent.DistalDubinsFrame(), joint.ProximalDubinsFrame())
            affected.append(data.joint_index)
        except Exception as e:
            print(f"[update_links] parent link error: {e}")

    for child_idx, p in enumerate(tree.Parents):
        if p == data.joint_index:
            child = tree.Joints[child_idx]
            try:
                tree.Links[child_idx] = _rebuild_link_warm(
                    tree, child_idx,
                    joint.DistalDubinsFrame(), child.ProximalDubinsFrame())
                affected.append(child_idx)
            except Exception as e:
                print(f"[update_links] child link error: {e}")

    return geometryToJson(getAffectedLinksGeometry(tree, affected))


# ── Static files (index.html, gizmo_test.html, …) ────────────────────────────

@app.get("/")
def root():
    return FileResponse(os.path.join(STATIC_DIR, "index.html"))

app.mount("/", StaticFiles(directory=STATIC_DIR), name="static")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
