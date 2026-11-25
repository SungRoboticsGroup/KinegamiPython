# Known Issues

Known bugs and limitations that haven't been fixed yet. Each entry describes what goes wrong, what it affects, and how to fix it. Remove an entry once it's fixed.

- [End effector collision checking under-approximates the geometry](#end-effector-collision-checking-under-approximates-the-geometry)
- [`realJointIndices()` raises `NameError`, and related `realJointsOnly` bugs](#realjointindices-raises-nameerror-and-related-realjointsonly-bugs)

---

## End effector collision checking under-approximates the geometry

**Type:** deliberate limitation, to be fixed. Added with end effector support, 2026-11-11.

### What happens
Collision detection (`KinematicTree.collision`) treats every object as a centerline plus a radius `.r`. An `EndEffector`'s collision shape (`EndEffector.sdf`) is a flat-ended cylinder of radius `r` running from the attachment face to the tool tip. For `PrintedParallelGripper`, `r` = 33 mm, the radius of the incoming tube.

The real gripper is wider than that:
- the carriage reaches up to √(54² + 31.65²) ≈ 62.6 mm from the axis;
- the jaws reach up to √(33.75² + 13.5²) ≈ 36.4 mm when fully open.

So **contacts with the carriage or jaws outside the 33 mm tube are not detected**. Collision results also don't depend on the gripper's state (how far the jaws are open).

### Scope
- **Affected:**
  - collision highlighting in the GUI;
  - `getCollidingPairs`, `detectCollisions` and `collisionsCountAndError`;
  - anything built on those, such as an optimization that penalizes collisions.
- **Not affected:**
  - joint placement, which uses `boundingBall()`, a ball enclosing the real gripper geometry at every jaw opening;
  - rendering, which draws the real parts.

### Why `r` can't just be made larger
`r` is also the link radius and the Dubins turning radius:
- `shortestCSC(newJoint.r, …)` uses it for path planning;
- joint placement uses `4*r` spacing;
- `placeJointAndWayPoints` asserts `neighbor.r == jointToPlace.r`.

So the gripper's `r` has to stay 33. A wider collision model needs its own radius or shape, separate from `r`.

### Where it's marked in the code
- `TODO(collision)` comments:
  - `EndEffector` (class docstring and `sdf`), `ParallelJawGripper` and `KinematicTree.collision`, all with pointers to this section;
  - `PrintedParallelGripper`'s docstring (no pointer).
- **Quick way to check it:** put a `PrintedParallelGripper` at the origin and a straight link (radius 33) running along y at x = 70 mm, through the middle of the carriage (z = −14.25 mm). This link physically overlaps the carriage, since its surface reaches x = 37 mm, inside the carriage's 54 mm half-width. But `KinematicTree.collision(gripper, link)` currently returns `None`. Once this issue is fixed, that check should report a collision.

### How collision checking works today (relevant to any fix)
- `collision(a, b)` samples `a.interpolate()` (the centerline) and tests `b.sdf(points) − a.r`, then does the same the other way round. **Both** directions must come within the threshold for a collision to count.
- If the deepest points on both sides are near their objects' endpoints (within `.r`), a disc test on `startCircle` / `endCircle` decides whether it's a real collision or just the place where a joint attaches to a link.
- The GUI checks every real joint against every link, **including its own incoming link**. So that disc test is what keeps each joint's attachment from being flagged.

### Fix options
| Model | Infrastructure needed | Risk |
|---|---|---|
| **A. Conservative bounding tube** (recommended first step) | 1. Add a `collisionRadius()` method to `Joint` and `LinkCSC`, defaulting to `self.r`; `EndEffector` overrides it (≈ 62.6 mm for the gripper).<br>2. In `collision()`, replace `tube.r` with `tube.collisionRadius()` at its 6 uses. Results for every existing joint and link type stay identical.<br>3. Make `EndEffector.sdf` a flat-ended cylinder of radius `collisionRadius()`.<br>4. Keep `startCircle` / `endCircle` at radius `r`, since they describe the physical attachment face.<br>5. Size it for the jaws fully open, so collision results stay independent of state. | Low. Before adding any gripper, check that `getCollidingPairs()` gives identical results on existing trees. Expect false positives: the tube is 62.6 mm wide along its full length, including the 33 mm housing, so a tight bend right before the gripper may be flagged. |
| B. Bounding capsule (rounded ends) | Same as A, but the SDF has no flat cuts. | The rounded back end reaches 62.6 mm back over the incoming link and the parent joint. Avoiding false collisions there would rely on the attachment disc test over a much larger region. Not recommended over A. |
| C. Stepped tubes (33 mm housing, then 62.6 mm carriage and jaws) | Either a separate radius at each sample point (`interpolate` returns radii; `−r` becomes `−r_i`; endpoint windows use the local radius), or one joint becomes several collision objects. The second means changing the collision matrices, `getAllCollisionPairs`, `collisionPairsFromMovingJoint`, and how the GUI maps (index, type) pairs to highlights. | Medium. It touches code shared by every pair type. |
| D. Exact part-based shape | A second, one-sided collision path for non-tube objects, with its own handling of the incoming link (the existing disc test doesn't apply). The collision shape would also depend on the gripper's state, so the state slider would have to mark collisions dirty and recompute them; fast slider drags could show stale highlights. | Highest. |

---

## `realJointIndices()` raises `NameError`, and related `realJointsOnly` bugs

**Type:** bug. It predates the end effector work.

### What happens
`KinematicTree.realJointIndices()` ([KinematicTree.py](KinematicTree.py)) filters joints with:

```python
not isinstance(self.Joints[i], Waypoint)
and not isinstance(self.Joints[i], PrintedWaypoint)   # PrintedWaypoint is undefined
and not isinstance(self.Joints[i], EndTip)
```

`PrintedWaypoint` was a class in the old `PrintedJoint.py`. It was deleted in the printed-pipeline rework (commits `6f48995f` and `5f814cee`, early Feb 2026), but this method, added in `321f578a`, still refers to it.

The file still imports fine, because Python looks up the name only when the line runs. Because `and` short-circuits, the lookup is skipped for `Waypoint` joints. The method therefore raises `NameError` at the first non-`Waypoint` joint, which means on essentially every real tree.

### Scope / relevance
- **Not reached by anything in the repo right now.** It is only called from the `realJointsOnly=True` option of `configuration`, `setConfiguration`, `stateRanges` and `randomConfiguration`, and no code calls those with that option.
- The GUIs, `dancePartyViewer.py` and `generateArduino.py` each build their own list of joints to skip (`Waypoint` and `Tip`).
- **It would affect scripts**, for example random-configuration sampling or the multi-configuration optimization these methods appear to have been written for.

### Related problems
1. **`setConfiguration(newConfig, realJointsOnly=True)` sets the wrong joints.** Its loop calls `self.setJointState(i, newConfig[i])`, where `i` is the position in the list rather than `jointIndex`. Once the `NameError` is fixed, it would quietly write states onto the wrong joints whenever waypoints come earlier in the tree.
2. **The definition of "real joint" doesn't match the rest of the code.** This method skips `EndTip` but not `StartTip`, while the GUI's configuration panel (`_should_skip_joint_in_config`) and `generateArduino.py` skip every `Tip`. Under either rule, `EndEffector`s (grippers) count as real, which is correct.
3. **The same undefined name appears in three other places:**
   - `KinematicTree.export3DKinematicTree`: its only call is commented out in an example.
   - `origamiGUI.generate_stl`: its button connection is commented out.
   - `examples/trees/manualHandPrinted.py` builds a `PrintedWaypoint` at the root, so that example fails immediately.

   The first two look like leftovers from the old STL-export pipeline and need a decision on whether to keep that code.

### Suggested fix
```python
def realJointIndices(self) -> np.ndarray:
    """Indices of actuated joints: everything except Waypoints and Tips
    (same rule as the GUI's configuration panel and generateArduino.py)."""
    return np.array([i for i, j in enumerate(self.Joints)
                     if not isinstance(j, (Waypoint, Tip))], dtype=int)
```
- In `setConfiguration`, change `self.setJointState(i, newConfig[i])` to `self.setJointState(jointIndex, newConfig[i])` in the `realJointsOnly` branch.
- Add a test with a tree containing waypoints, a start tip, an end tip, and a gripper. It should check that `configuration(realJointsOnly=True)` → `setConfiguration(..., realJointsOnly=True)` round-trips and touches only the intended joints.
- Separately, decide whether to delete or update `export3DKinematicTree`, `origamiGUI.generate_stl` and `examples/trees/manualHandPrinted.py`.
