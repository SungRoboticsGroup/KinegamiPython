# PSO-based kinematic tree optimization

This document covers the code that shortens a kinematic tree by moving its joints and waypoints, using particle swarm optimization (PSO). It describes how the code is structured, which traversals and parameters are available, how to run and debug it, and where development stands.

*Status as of 2026-11-20, branch `pso-fixes`.*

## What it does

The goal is to make the tree as compact as possible without changing what it can do:
- Each joint is free to slide along its own z-axis and rotate about it. Neither move changes the robot's kinematics.
- Each waypoint is free to move anywhere (a full 6-DOF pose).
- The links between them are recomputed as CSC Dubins paths.

The tree is optimized **one joint or waypoint at a time**, following a traversal order. Each optimization of one joint or waypoint is a **visit**. A visit uses PSO to search the possible moves, minimizing a local objective that depends only on the links next to the moved joint or waypoint. A full traversal is a **pass**, and runs usually do several passes.

The result must be collision-free in every configuration in a given list of joint configurations (sets of joint states). The tree's current state doesn't matter, and neither does any configuration that isn't on the list.

## Code map

| File | What's in it |
|---|---|
| [optimizationFunctions.py](optimizationFunctions.py) | Everything below except the objective and traversals: `optimizeTree`, the joint and waypoint optimizers, the retry/revert logic, checkpointing |
| [treeTraversals.py](treeTraversals.py) | The local objective `linkLoss`, and the traversal generators `dfs`, `bfs`, `randomized` |
| [KinematicTree.py](KinematicTree.py) | Collision detection: `collision`, `collisionsCountAndError`, `detectCollisions`, `buildCollisionMatrices`. Also `transformJoint`, which every visit uses to apply its move |
| [randomChain.py](randomChain.py) | Experiment script for random chains built with the generalized gimbal construction |
| [randomTree.py](randomTree.py) | Experiment script for random branching trees |
| [manageCheckpoints.py](manageCheckpoints.py) | Examples of loading, inspecting and resuming checkpoints |
| [view.py](view.py), [viewBrokenExample.py](viewBrokenExample.py) | Scripts for viewing saved results and stepping through a failing example |

### Call structure

```
optimizeTree(subject, configurations, traversal, ...)
 ├─ failurePenalty = max over leaves of the sum of squared link lengths from root to leaf
 └─ for each pass, for each index in the traversal order:
      ├─ optimizeJointPlacement     (real joints)  ─┐
      └─ optimizeWaypointPlacement  (waypoints)    ─┴─ optimizeWithCollisionRetry
             1. optimize...Once(includeCollisionPenalty=False)
             2. if any configuration collides: optimize...Once(includeCollisionPenalty=True)
             3. if the result still collides, or loss > initial loss: revert to the previous
                placement and log a warning
      ├─ log the tree (.tree file) and its total length
      └─ save a checkpoint (.dill)
```

The no-penalty attempt runs first because the penalty is expensive: it runs collision checks in every configuration for every particle evaluation. Most visits don't need it.

### One joint visit (`optimizeJointPlacementOnce`)

- **Search space:** 2-D, `[translation along z, rotation about z]`, relative to the joint's current pose.
  - Each candidate is evaluated two ways, as is and with the joint's z-axis reversed, and the cheaper one counts.
- **Initial swarm (16 particles):**
  - Particle 0 is a guess that slides the joint up against its parent's distal frame.
  - Particle 1 is `[0, 0]`, i.e. the current placement. This guarantees the result is no worse than where the visit started.
  - The other 14 are the guess plus uniform noise, then clipped (translation) or wrapped (angle) to fit the bounds.
- **Bounds:**
  - Translation covers ±2·`localLinkScale` around both 0 and the guess. `localLinkScale` is the longest link touching the joint (incoming or outgoing), and never less than the tube radius r.
  - Rotation covers ±2π.
- **Optimizer:** `GlobalBestPSO` (c1 = 0.6, c2 = 0.7, w = 0.5) for 25 iterations, then 25 iterations of Nelder-Mead starting from the swarm's best result.

### One waypoint visit (`optimizeWaypointPlacementOnce`)

- **Search space:** 6-D, `[tx, ty, tz, ZYZ Euler angles]`. The move is expressed in the **waypoint's own frame** (`relative=True`), so rotations turn the waypoint in place.
- **Objective:** `linkLoss` plus 10 × the angle (in radians) between the waypoint's orientation and the orientation that would line it up with its parent's distal frame.
- **Initial swarm (24 particles):** particle 0 lines the waypoint up with its parent, particle 1 is the current placement, and the rest are random noise.
  - Translation bounds: ±2·`localLinkScale` per axis, widened to include the line-up guess.
  - Angle bounds: ±2π.
- **Optimizer:** `GlobalBestPSO` (c1 = 0.7, c2 = 0.5, w = 0.5) for 50 iterations. Unlike joint visits, there's no Nelder-Mead step afterward.

### Objective

`linkLoss`, in [treeTraversals.py](treeTraversals.py):

```
loss = (incoming link length)^power
     + childFraction * sum over outgoing links of (length)^power
     + [if includeCollisionPenalty] sum over configurations of collisionErrorWeight * collisionError
```

`collisionError` is a logistic function of the signed surface-to-surface distance, and it reaches about 1 within `fineDistanceThreshold` = 0.001. In practice any collision adds the full `collisionErrorWeight` (set to `failurePenalty`), however deep it is.

Two consequences:
- The penalty gives PSO no gradient to follow out of a collision. It only makes colliding placements look bad.
- The visits minimize **squared** link lengths, while the experiments report **total length**. A visit can lower its own objective and still slightly raise total length.

## Options

### `optimizeTree` arguments

| Argument | Values | Notes |
|---|---|---|
| `traversal` | `"dfs"`, `"bfs"`, `"randomized"` | See [Traversals](#traversals) |
| `direction` | `"outward"`, `"inward"` | dfs/bfs only. Outward is root → leaves; inward is the reverse of the outward order |
| `orderBy` | `"default"`, `"longest"`, `"shortest"` | dfs/bfs only. The order in which to visit siblings at each branch. **Has no effect on chains** |
| `repeatTraversal` | int, or `"n"` | Number of passes. `"n"` means one pass per entry in `tree.Joints` (joints and waypoints together) |
| `configurations` | list of state lists | Each state list has an entry for every joint and waypoint (0 for waypoints). The result must be collision-free in all of them. `None` means a single all-zeros configuration |
| `childFraction` | float ≥ 0 | Weight on outgoing links. **Currently only used by the randomized traversal**; see [Known issues](#known-issues) |
| `power` | float ≥ 0 | **Only used by the randomized traversal's weights.** Visits always use power 2 |
| `guarantee` | bool | **Always pass `True`.** `False` makes visits optimize a cut-down subtree; that path is untested and slated for removal |
| `parallelize` | bool | Evaluates particles across processes. Particle sampling is seeded by `np.random`, but results aren't guaranteed to be identical from run to run |
| `directory` | path or `None` | Where per-visit `.tree` files and `checkpoints/` go |
| `save_checkpoints` | bool | Save a `.dill` checkpoint after every visit (default `True`) |
| `evaluate` | bool | Return `(tree, times, lengths)` instead of just the tree |
| `showSteps`, `verbose` | bool | Display and printing only |
| `resetOnFail` | | Unused |

### Traversals

- **`dfs`:** depth-first from the root. With `longest`/`shortest`, it visits children in order of each child's longest/shortest path length down to a leaf.
- **`bfs`:** level by level. With `longest`/`shortest`, it sorts each level by incoming link length.
- **`randomized`:** each pass makes `len(tree.Joints)` random visits, excluding the root. A joint or waypoint is chosen with probability proportional to its `linkLoss` (using `power` and `childFraction`, no collision penalty).

The traversal order is computed from the tree passed to `optimizeTree`, i.e. the starting tree, and isn't updated as the tree changes. ([Known issues](#known-issues) covers what happens on resume.)

The experiment scripts use 9 variants: {DFS, BFS} × {outward, inward} × {longest, shortest}, plus randomized with power 3.

### PSO and sampling parameters

These are hard-coded in [optimizationFunctions.py](optimizationFunctions.py), not exposed as arguments:

| Parameter | Joint | Waypoint |
|---|---|---|
| Algorithm | `ps.single.GlobalBestPSO` | `ps.single.GlobalBestPSO` |
| Particles | 16 | 24 |
| c1, c2, w | 0.6, 0.7, 0.5 | 0.7, 0.5, 0.5 |
| PSO iterations | 25 (`(maxiter+1)/2`) | 50 (`maxiter`) |
| Local refinement | Nelder-Mead, 25 iterations | none |
| Tolerance | `r/10` (PSO `ftol` and Nelder-Mead `fatol`) | `r/10` |
| Translation bounds | ±2·`localLinkScale` around 0 and the guess | same, per axis |
| Angle bounds | ±2π, wrapped | ±2π per angle, wrapped |

`maxiter = 50` is set in `optimizeTree`. Spring 2026 meetings tried `LocalBestPSO` and other values of its neighbor count k and inertia w (see [Development status](#development-status)). That code is not in this branch.

## Running

### Experiment scripts

Settings are module-level variables at the top of each script. Edit them, then run the script:

```
python randomChain.py     # jointCount, seed, cubeSize, chainCount; the traversals to run are the `versions` list in test()
python randomTree.py      # jointCount, probabilityOfBranching, seed, treeCount; the `optimizations` list in testRandomTrees()
```

**Warning:** both scripts create their results folder as soon as they're imported, not just when run. So don't import from them; copy the function you need instead.

Each script:
1. Generates a random tree (seeded).
2. Picks collision-free configurations at random (`find_collision_free_configs`).
3. Runs each traversal variant on that tree, resetting the random seed before each one.
4. Plots time vs. length for every variant against `totalLengthLowerBound()`.

`randomChain.py` uses `repeatTraversal="n"`. `randomTree.py` uses the default of 1 pass.

**How long it takes:** an 8-joint gimbal chain has 23 entries in `tree.Joints` (joints plus waypoints). At `"n"` passes that's 529 visits, about 65–70 minutes per traversal. A visit that needs the penalized retry takes roughly 20 times longer than one that doesn't.

### Output layout

```
Trials Before Experiments/<timestamp>_Joints<N>_Chains<C>_Seed<S>/
  run_info.json                   git commit, branch, uncommitted changes, package versions, script settings
  0.tree                          starting tree (repr)
  0_configs.json                  configurations used
  collision_penalty_log.txt       penalized retries, reverts and their reasons, and any traversal that raised
  random_results_chkpt0.json      (times, lengths) per traversal, updated after each traversal (null = it failed)
  random_results.json, plot_0.png once all trees are done
  0/<traversal label>/
     visits.csv                       one row per visit (see below)
     <elapsed seconds>_<index>.tree   tree after each visit (index -1 = start)
     final.tree
     checkpoints/visit00012_iter_<pass>_<index>.dill
```

`visits.csv` has one row per visit, plus a `visit = 0` row for the starting tree:

| Column | Meaning |
|---|---|
| `visit`, `repeat`, `index`, `kind` | Visit number (matches the checkpoint name), pass, index in `tree.Joints`, and `joint` or `waypoint` |
| `initial_loss`, `final_loss` | The visit's objective at the current placement, and after the visit. If the penalized retry ran, both include the penalty |
| `retried`, `reverted`, `revert_reason` | Whether the penalized retry ran, and whether (and why) the visit was reverted |
| `nopenalty_pso_s`, `nopenalty_nm_s` | Seconds spent on PSO and Nelder-Mead in the no-penalty attempt (Nelder-Mead is 0 for waypoints) |
| `penalty_pso_s`, `penalty_nm_s` | The same for the penalized retry; blank if it didn't run |
| `collision_check_s` | Seconds spent on the collision checks after each attempt. Checks inside the penalized objective run in PSO's worker processes and are counted in `penalty_pso_s` instead |
| `visit_s`, `elapsed_s` | The visit's total time, and the time since the run started |
| `total_length`, `sum_squared_length` | The tree's total link length, and its sum of squared link lengths (the tree-wide version of the visits' objective) |

For example, with pandas (`pip install pandas`; it isn't in `requirements.txt`):

```python
import pandas as pd
visits = pd.read_csv(".../0/<traversal label>/visits.csv")
visits.groupby("kind")[["nopenalty_pso_s", "penalty_pso_s", "visit_s"]].sum()   # where the time goes
visits[visits.reverted == True]                                                # every revert and why (the visit-0 row is blank)
```

### Checkpoints: inspecting and resuming

Each checkpoint saves:
- the tree;
- the `np.random` state;
- the traversal's progress (pass, completed indices);
- the times, lengths and losses so far;
- the configurations, the `failurePenalty`, and the run's settings.

```python
from optimizationFunctions import listCheckpoints, loadCheckpoint, inspectCheckpoint, resumeFromCheckpoint

cps = listCheckpoints(".../checkpoints")       # sorted in visit order
state = loadCheckpoint(cps[-1]["path"])        # state.tree, state.configurations, state.failure_penalty, ...
resumeFromCheckpoint(cps[-1]["path"], directory="some/new/dir")
```

- **Older names:** checkpoints saved before visit numbering are named `iter_<pass>_<index>.dill`. Within a pass, their names sort by joint index, **not** by visit order. For inward traversals the last name alphabetically is the *first* visit of that pass.
- **Where resumed runs write:** pass a new `directory` to `resumeFromCheckpoint`, so the resumed run doesn't overwrite the original run's checkpoints.

### Debugging a single visit

To reproduce one visit exactly, load the checkpoint just before it with `restore_random_state=True`. Then call `optimizeJointPlacement` or `optimizeWaypointPlacement` directly with `state.tree`, `state.failure_penalty` and `state.configurations`.

To check collisions, call `collidingConfigurations(tree, index, configurations)`. Avoid `tree.detectCollisions()` for this, because it only checks the tree's current state, which usually isn't one of the configurations.

## Tests

There's no automated test suite for the optimizer, and pytest isn't installed. The `test_*.py` files in the repo root are standalone scripts (`python test_sdf_aabb.py`, etc.). The ones relevant here cover collision geometry:
- [test_sdf_aabb.py](test_sdf_aabb.py)
- [test_disc_intersection.py](test_disc_intersection.py)
- [test_endcap_filtering.py](test_endcap_filtering.py)
- [test_collision_free_configs.py](test_collision_free_configs.py)

Checks to run after changing optimizer code:
1. **Small run from start to finish:** a 4-joint random chain (`np.random.seed(42)`, 11 entries in `tree.Joints`, 4 configurations), 2 passes of inward DFS longest. On `pso-fixes` it reaches length 77.26 → 66.53 with 0 colliding configurations, and resuming from the middle checkpoint completes.
2. **Replaying a known-hard visit:** joint 9 from `Trials Before Experiments/2026.04.06_21.38.31_Joints8_Chains1_Seed44/0/DFS Inward Longest n repetitions/checkpoints/iter_18_10.dill`. It should take the penalized retry and finish collision-free, with loss 323.58 → 244.34.

## Development status

### History

- **Dec 2025 – Jan 2026:**
  - the SDF collision checking from the `collisions` branch was merged in;
  - collision checks were restricted to joint/link pairs that can actually collide;
  - multi-configuration optimization was added;
  - the gimbal chain generator was added.
- **Feb – Apr 2026 (Samhitha):**
  - collisions are measured as signed surface-to-surface distances;
  - `collisionErrorWeight` is tied to `failurePenalty`;
  - checkpoints are saved after every visit.
  - The 8-joint runs crashed on inward traversals.
- **May 2026 (not pushed):** Samhitha experimented with `LocalBestPSO` and with k and w, and added a revert when a retry still collided. Her notes suggest a lower k and higher w did better on one example; the 5/7 meeting notes have the details. This code isn't in the repo.
- **Nov 2026 (`pso-fixes`):**
  - **False collisions fixed.** The final collision checks used the tree's current state, where the gimbal construction leaves prismatic joints fully extended, instead of the configuration list. That flagged results that were collision-free in every configuration, and it caused the April crashes.
  - **Reverts.** A visit that still collides, or that raises its own objective, is now reverted with a logged warning instead of raising.
  - **Joint sampling.** Particles were clipped to ±2× the incoming link length, which is often about 0 right after the gimbal construction. That collapsed the whole swarm onto one point.
  - **Waypoint moves.** These used to be in the world frame, so rotations swung waypoints around the world origin.
  - **Checkpoints.** Names now sort in visit order.

### Known issues

- **`childFraction` has no effect on joint or waypoint placement.** Neither optimizer passes it to `linkLoss`, which defaults to 1. It needs wiring up before it can be an experiment parameter.
- **Visits always use power 2.** `optimizeTree`'s `power` only changes the randomized traversal's weights.
- **The root joint has a broken initial guess.** DFS and BFS visit the root (index 0), but the root has no parent: `Parents[0] = -1`, so the joint optimizer uses `Joints[-1]` (the last joint) to compute its guess. The current-placement particle keeps this from making things worse, but the guess is meaningless. The randomized traversal never visits the root.
- **Resuming can change the traversal order.** `optimizeTree` computes the order from the starting tree. `resumeFromCheckpoint` recomputes it from the checkpoint's tree, so `longest`/`shortest` orderings can differ after a resume.
- **The collision penalty is all-or-nothing** (see [Objective](#objective)), so PSO has no gradient to follow out of a collision.
- **Reverts haven't been seen in a real run yet.** They're implemented and checked by reading the code, but the small test run didn't trigger one.
- **Runs aren't exactly reproducible.** Earlier 3-joint runs with the same seed ended between 18.1 and 50.6, partly because the code changed between runs. Comparing variants needs several trials per setting.
- **Tree branches have fixes in the GUI branch.** `chainsGUIdevelopment` has collision fixes that `pso-fixes` doesn't: the sibling-waypoint clustering fix `af6b0359`, which matters for branching trees, and the short-tube end-disc fix `392f2061`. Port them before running trees.
- **Hard-coded paths.** [manageCheckpoints.py](manageCheckpoints.py) and [view.py](view.py) contain `/home/samhitha/...` paths.

### Before systematic experiments

1. Rerun the 8-joint seed-44 chain on all traversals, and confirm none crash and how often visits revert.
2. Wire up `childFraction`, and decide whether to also expose `power` for visits.
3. Decide which metric to report, total length or the visits' sum of squares.
4. Port the clustering fix from the GUI branch, then smoke-test `randomTree.py`.
5. Profile:
   - the split between PSO and Nelder-Mead time;
   - joint visits vs. waypoint visits;
   - how collision checking scales with the number of configurations.
6. Fix the experiment list. The parameter space from the 5/7 meeting notes:
   - **Inputs:** number of joints, number of configurations.
   - **Our parameters:** traversal, number of passes (or visits), `childFraction`.
   - **PSO parameters:** global vs. local best; number of particles; sampling bounds; initial particles (uniform random vs. evenly spaced, from [TODO_tests.md](TODO_tests.md)); k; w.

   The planned plots compare time per joint against the objective divided by the best objective found by a reference method.
