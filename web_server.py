"""
web_server.py — Kinegami web backend (FastAPI)

Run with:
    python web_server.py

Then open http://localhost:8000 in a browser.
"""
import os
import copy
import math
import numpy as np
from fastapi import FastAPI, UploadFile, File
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel
from typing import List

app = FastAPI()
STATIC_DIR = os.path.join(os.path.dirname(__file__), "web_static")


# ── Persistent tree state + undo/redo history ─────────────────────────────────

_tree = None
_history: list = []       # list of dill-serialised tree snapshots
_history_index: int = -1  # points to current position in _history


def _push_history(tree) -> None:
    """Snapshot tree into history, discarding any redo future."""
    global _history, _history_index
    import dill
    snap = dill.dumps(tree)
    _history = _history[:_history_index + 1]
    _history.append(snap)
    _history_index = len(_history) - 1


_EMPTY_GEO = {"meshes": [], "lines": [], "joints": [], "parents": [],
              "maxAnglePerElbow": math.pi / 12}


def _geometry_response(tree):
    """Return geometry JSON for tree, or empty geometry if tree is None."""
    if tree is None:
        return _EMPTY_GEO
    from webGeometry import getTreeGeometry, geometryToJson
    return geometryToJson(getTreeGeometry(tree))


# ── Geometry endpoint ─────────────────────────────────────────────────────────

@app.get("/geometry")
def get_geometry():
    return _geometry_response(_tree)


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

    tree  = _tree
    joint = tree.Joints[data.joint_index]

    old_pose           = copy.deepcopy(joint.Pose)
    old_proximal_dubins = getattr(joint, 'proximalDubins', None)
    old_distal_dubins   = getattr(joint, 'distalDubins',   None)

    R = Rotation.from_quat(data.quaternion).as_matrix()
    T = np.eye(4)
    T[:3, :3] = R
    T[:3,  3] = data.position
    joint.Pose = SE3(T)
    joint.proximalDubins = joint.ProximalDubinsFrame()
    joint.distalDubins   = joint.DistalDubinsFrame()

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
        joint.Pose = old_pose
        if old_proximal_dubins is not None:
            joint.proximalDubins = old_proximal_dubins
        if old_distal_dubins is not None:
            joint.distalDubins = old_distal_dubins
        print(f"[move_joint] rejected (link rebuild failed): {failure}")
        return JSONResponse(status_code=422, content={"error": failure})

    for idx, link in new_links.items():
        tree.Links[idx] = link

    _push_history(tree)
    return geometryToJson(getTreeGeometry(tree))


# ── Live link-update endpoint (called every frame during drag) ────────────────

@app.post("/update_links")
def update_links(data: MoveJointRequest):
    from scipy.spatial.transform import Rotation
    from spatialmath import SE3
    from webGeometry import getAffectedLinksGeometry, geometryToJson

    tree  = _tree
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


# ── Add joint endpoint ────────────────────────────────────────────────────────

class AddJointRequest(BaseModel):
    joint_type:   str   # "TransverseRevolute" | "CoaxialRevolute" | "Tip"
    parent_index: int = 0


@app.post("/add_joint")
def add_joint_ep(data: AddJointRequest):
    global _tree
    import math as _math
    from spatialmath import SE3
    from PrintedTube import (TransverseRDS3225, CoaxialRDS3225,
                              PrintedEndHemisphere, PrintedKinematicTree)

    try:
        tree_empty = (_tree is None)

        if data.joint_type == "TransverseRevolute":
            if tree_empty:
                # Root: point along +X, matching the demo tree convention
                pose = SE3.Ry(-_math.pi / 2)
            else:
                # Child: 4r + half neutral-length in front along distal dubins X axis
                distance = 4 * TransverseRDS3225.R + TransverseRDS3225.NEUTRAL_LENGTH / 2
                pose = SE3.Rt(np.eye(3), np.array([distance, 0, 0]))
            new_joint = TransverseRDS3225(pose, version=270)

        elif data.joint_type == "CoaxialRevolute":
            if tree_empty:
                pose = SE3()
            else:
                distance = 4 * CoaxialRDS3225.R + CoaxialRDS3225.NEUTRAL_LENGTH / 2
                pose = SE3.Rt(SE3.Ry(_math.pi / 2).R, np.array([distance, 0, 0]))
            new_joint = CoaxialRDS3225(pose, version=270)

        elif data.joint_type == "Tip":
            if tree_empty:
                pose = SE3()
            else:
                parent_r = _tree.Joints[data.parent_index].r
                distance = parent_r * 4 + TransverseRDS3225.R / 2
                pose = SE3.Rt(SE3.Ry(_math.pi / 2).R, np.array([distance, 0, 0]))
            new_joint = PrintedEndHemisphere(r=TransverseRDS3225.R, Pose=pose)

        else:
            return JSONResponse(status_code=400,
                                content={"error": f"Unknown joint type: {data.joint_type}"})

        if tree_empty:
            _tree = PrintedKinematicTree(new_joint)
        else:
            if data.parent_index < 0 or data.parent_index >= len(_tree.Joints):
                return JSONResponse(status_code=400, content={"error": "Invalid parent index"})
            _tree.addJoint(data.parent_index, new_joint,
                           relativeToDistalDubins=True,
                           fixedPosition=True, fixedOrientation=True,
                           safe=False)

        _push_history(_tree)
        return _geometry_response(_tree)
    except Exception as e:
        return JSONResponse(status_code=422, content={"error": str(e)})


# ── Delete joint endpoint ─────────────────────────────────────────────────────

class DeleteJointRequest(BaseModel):
    joint_index: int


@app.post("/delete_joint")
def delete_joint_ep(data: DeleteJointRequest):
    global _tree
    import dill

    if _tree is None:
        return JSONResponse(status_code=400, content={"error": "Tree is empty"})

    # Deleting the only remaining joint empties the tree
    if len(_tree.Joints) == 1 and data.joint_index == 0:
        _tree = None
        _push_history(_tree)
        return _EMPTY_GEO

    backup = dill.dumps(_tree)
    try:
        _tree.deleteJoint(data.joint_index)
        _push_history(_tree)
        return _geometry_response(_tree)
    except Exception as e:
        _tree = dill.loads(backup)
        return JSONResponse(status_code=422, content={"error": str(e)})


# ── Set joint state endpoint ──────────────────────────────────────────────────

class SetJointStateRequest(BaseModel):
    joint_index:  int
    state:        float   # degrees for Revolute, raw units for others
    push_history: bool = True


@app.post("/set_joint_state")
def set_joint_state_ep(data: SetJointStateRequest):
    global _tree
    import dill
    from Joint import Revolute
    from webGeometry import getTreeGeometry, geometryToJson

    tree   = _tree
    backup = dill.dumps(tree)
    try:
        joint = tree.Joints[data.joint_index]
        actual_state = (math.radians(data.state)
                        if isinstance(joint, Revolute) else data.state)

        success = tree.setJointState(data.joint_index, actual_state)
        if not success:
            return JSONResponse(status_code=422, content={"error": "State out of range"})

        tree.resyncFromLightweight()
        if data.push_history:
            _push_history(tree)
        return geometryToJson(getTreeGeometry(tree))
    except Exception as e:
        _tree = dill.loads(backup)
        return JSONResponse(status_code=422, content={"error": str(e)})


# ── Set all joint states endpoint (batch) ────────────────────────────────────

class SetAllJointStatesRequest(BaseModel):
    states:       dict   # {str_idx: float_degrees}
    push_history: bool = False


@app.post("/set_all_joint_states")
def set_all_joint_states_ep(data: SetAllJointStatesRequest):
    global _tree
    import dill
    from Joint import Revolute
    from webGeometry import getTreeGeometry, geometryToJson

    tree   = _tree
    backup = dill.dumps(tree)
    try:
        for idx_str, state_deg in data.states.items():
            idx   = int(idx_str)
            joint = tree.Joints[idx]
            actual_state = (math.radians(state_deg)
                            if isinstance(joint, Revolute) else state_deg)
            success = tree.setJointState(idx, actual_state)
            if not success:
                raise ValueError(f"Joint {idx} state {state_deg:.1f}° out of range")
        tree.resyncFromLightweight()
        if data.push_history:
            _push_history(tree)
        return geometryToJson(getTreeGeometry(tree))
    except Exception as e:
        _tree = dill.loads(backup)
        return JSONResponse(status_code=422, content={"error": str(e)})


# ── Clear tree endpoint ───────────────────────────────────────────────────────

@app.post("/clear_tree")
def clear_tree_ep():
    global _tree
    _tree = None
    _push_history(_tree)
    return _EMPTY_GEO


# ── Undo / redo endpoints ─────────────────────────────────────────────────────

@app.post("/undo")
def undo_ep():
    global _history_index, _tree
    import dill

    if _history_index <= 0:
        return JSONResponse(status_code=400, content={"error": "Nothing to undo"})
    _history_index -= 1
    _tree = dill.loads(_history[_history_index])
    return _geometry_response(_tree)


@app.post("/redo")
def redo_ep():
    global _history_index, _tree
    import dill

    if _history_index >= len(_history) - 1:
        return JSONResponse(status_code=400, content={"error": "Nothing to redo"})
    _history_index += 1
    _tree = dill.loads(_history[_history_index])
    return _geometry_response(_tree)


# ── Save / load tree endpoints ────────────────────────────────────────────────

@app.get("/save_tree")
def save_tree_ep():
    import dill
    data = dill.dumps(_tree)
    return Response(
        content=data,
        media_type="application/octet-stream",
        headers={"Content-Disposition": "attachment; filename=kinegami_tree.pkl"},
    )


@app.post("/load_tree")
async def load_tree_ep(file: UploadFile = File(...)):
    global _tree
    import dill
    from webGeometry import getTreeGeometry, geometryToJson

    try:
        data  = await file.read()
        _tree = dill.loads(data)
        _push_history(_tree)
        return geometryToJson(getTreeGeometry(_tree))
    except Exception as e:
        return JSONResponse(status_code=422, content={"error": str(e)})


# ── Export link modules endpoint (stub) ───────────────────────────────────────

@app.get("/export_link_modules")
def export_link_modules_ep():
    import io, zipfile
    from webGeometry import getTreeGeometry

    tree = _tree
    buf  = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as zf:
        for i, link in enumerate(tree.Links):
            if link is None:
                continue
            try:
                stl_bytes = link.toSTL()
                zf.writestr(f"link_{i:02d}.stl", stl_bytes)
            except Exception as e:
                zf.writestr(f"link_{i:02d}_error.txt", str(e))
    buf.seek(0)
    return Response(
        content=buf.read(),
        media_type="application/zip",
        headers={"Content-Disposition": "attachment; filename=link_modules.zip"},
    )


# ── Static files (index.html, …) ──────────────────────────────────────────────

@app.get("/")
def root():
    return FileResponse(os.path.join(STATIC_DIR, "index.html"))

app.mount("/", StaticFiles(directory=STATIC_DIR), name="static")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
