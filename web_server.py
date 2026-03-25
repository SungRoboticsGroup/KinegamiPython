"""
web_server.py — Kinegami web backend (FastAPI)

Run with:
    python web_server.py

Then open http://localhost:8000 in a browser.
"""
import os
import numpy as np
from fastapi import FastAPI
from fastapi.responses import FileResponse
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


# ── Move-joint endpoint ───────────────────────────────────────────────────────

class MoveJointRequest(BaseModel):
    joint_index: int
    position:   List[float]   # world-space [x, y, z] in mm
    quaternion: List[float]   # [x, y, z, w]  (scipy / Three.js convention)


@app.post("/move_joint")
def move_joint(data: MoveJointRequest):
    from scipy.spatial.transform import Rotation
    from spatialmath import SE3
    from webGeometry import getTreeGeometry, geometryToJson

    tree  = get_tree()
    joint = tree.Joints[data.joint_index]

    # Apply new world pose
    R = Rotation.from_quat(data.quaternion).as_matrix()
    T = np.eye(4)
    T[:3, :3] = R
    T[:3,  3] = data.position
    joint.Pose = SE3(T)
    # Refresh cached Dubins frames
    joint.proximalDubins = joint.ProximalDubinsFrame()
    joint.distalDubins   = joint.DistalDubinsFrame()

    mk_link = tree._get_link_constructor()

    # Rebuild link from parent → this joint
    parent_idx = tree.Parents[data.joint_index]
    if parent_idx >= 0:
        parent = tree.Joints[parent_idx]
        try:
            tree.Links[data.joint_index] = mk_link(
                tree.r,
                parent.DistalDubinsFrame(),
                joint.ProximalDubinsFrame(),
                tree.maxAnglePerElbow)
        except Exception as e:
            print(f"[move_joint] could not rebuild link to joint "
                  f"{data.joint_index}: {e}")

    # Rebuild links from this joint → each child
    for child_idx, p in enumerate(tree.Parents):
        if p == data.joint_index:
            child = tree.Joints[child_idx]
            try:
                tree.Links[child_idx] = mk_link(
                    tree.r,
                    joint.DistalDubinsFrame(),
                    child.ProximalDubinsFrame(),
                    tree.maxAnglePerElbow)
            except Exception as e:
                print(f"[move_joint] could not rebuild link to child "
                      f"{child_idx}: {e}")

    return geometryToJson(getTreeGeometry(tree))


# ── Static files (index.html, gizmo_test.html, …) ────────────────────────────

@app.get("/")
def root():
    return FileResponse(os.path.join(STATIC_DIR, "index.html"))

app.mount("/", StaticFiles(directory=STATIC_DIR), name="static")


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
