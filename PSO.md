# PSO-based kinematic tree optimization

This document covers the code that shortens a kinematic tree by moving its joints and waypoints, using particle swarm optimization (PSO). It describes how the code is structured, which traversals and parameters are available, how to run and debug it, and where development stands.

*Status as of 2026-11-22, branch `pso-fixes`. For why each experiment was run and what it showed, see [EXPERIMENTS.md](EXPERIMENTS.md).*

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
| [parameterSweep.py](parameterSweep.py) | Overnight sweep over the visit settings on small random chains; see [Parameter sweeps](#parameter-sweeps) |
| [bunchingExamples.py](bunchingExamples.py), [debugExamples/bunching/](debugExamples/bunching/) | Saved bunched waypoints with known better placements, a harness that runs single waypoint visits on them, and `surveyBunching.py`, which found the placements |
| [EXPERIMENTS.md](EXPERIMENTS.md) | Experiment log: the question, command, result and conclusion for each experiment |
| [viewConfigurations.py](viewConfigurations.py) | Interactive viewer that switches between a tree's configurations and highlights near-collisions |
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
- **Objective:** `linkLoss` plus `alignmentWeight` (default 10) × the angle (in radians) between the waypoint's orientation and the orientation that would line it up with its parent's distal frame.
- **Initial swarm (24 particles):** particle 0 lines the waypoint up with its parent, particle 1 is the current placement, and the rest are random.
  - With `straightThroughGuess=True`, two more particles point the waypoint along the line from the end of the previous part to the start of the next one: one turned in place, one also moved onto that line.
  - With `waypointSampling="box"` (default), each random particle is the line-up guess plus uniform noise: ±2·`localLinkScale` per translation axis and ±2π per ZYZ angle.
  - With `waypointSampling="multiscale"`, the random particles are split into thirds that rotate in place, move up to ±0.5·`localLinkScale` per axis, or move anywhere within the bounds. All of them get rotations sampled uniformly over orientations. Box sampling almost never pairs a large rotation with a small move, which is what fixing a bunched waypoint usually needs.
  - Translation bounds: ±2·`localLinkScale` per axis, widened to include the line-up guess.
  - Angle bounds: ±2π.
- **Optimizer:** `GlobalBestPSO` (c1 = 0.7, c2 = 0.5, w = 0.5) for 50 iterations. Unlike joint visits, there's no Nelder-Mead step afterward.

### Parent-locking

With `parentLocking=True` (the default whenever `guarantee=True`), a waypoint can be **parent-locked**: its proximal Dubins frame coincides with its parent's distal Dubins frame (within 1e-6), so its incoming link is empty. The tree then behaves as if the waypoint were deleted, but it stays in the tree in case it's needed again. Whether a waypoint is parent-locked is read from the geometry (`isParentLocked`), so it needs no extra saved state, and existing trees with empty links count as parent-locked.

- **Locked waypoints move with what they're locked to.** Every move a joint or waypoint visit considers is also applied to its parent-locked descendants (recursively), which are snapped back onto their parents' distal frames (`moveWithParentLocked`). The links leaving them count as the visited joint's outgoing links in `linkLoss`. Without this, a joint with a waypoint sitting on its distal frame can barely move: shifting the joint forward makes the path to the waypoint loop back.
- **Waypoint visits consider locking.** After PSO, the visit also evaluates the waypoint locked onto its parent, and locks it if that's better on `linkLoss` by at least 0.1% (`PARENT_LOCK_MARGIN`).
- **Visits to a locked waypoint consider unlocking it.** The swarm gets one more particle, the **unlocking guess**: the middle of the path leaving the waypoint, pointing along that path. Both halves of a feasible path are feasible, so this is a safe way out. The waypoint is unlocked if the swarm's best is unlocked and better on `linkLoss` by the margin, or else if the unlocking guess is.
- **Lock decisions use `linkLoss` alone,** without the alignment penalty, which is always 0 for a locked waypoint and would favor locking for reasons unrelated to length. (Inside PSO the alignment penalty still applies, which is why the unlocking guess is also judged separately: with `alignmentWeight` 10, the swarm almost never picks it.) When a visit unlocks via the guess, its `final_loss` is reported as that `linkLoss`.
- **What to expect at each power.** At power 2, merging links of lengths a and b costs (a+b)², more than a² + b², so visits rarely lock, and a locked waypoint is unlocked as soon as it's visited. At power 1, locking and splitting cost about the same, and the margin keeps waypoints from flipping back and forth.
- **Branching:** the unlocking guess uses the path to the waypoint's first child. Waypoints in the current tree construction have one child each.

### Objective

`linkLoss`, in [treeTraversals.py](treeTraversals.py):

```
loss = (incoming link length)^power
     + childFraction * sum over outgoing links of (length)^power
     + [if includeCollisionPenalty] sum over configurations of collisionErrorWeight * collisionError
```

A pair counts as a collision only when the tubes really overlap: when the surface-to-surface distance, measured both ways (from each tube's center line to the other's surface), is below 0 (more precisely, below −1e-9·r, to absorb floating-point error when tubes exactly touch). Tubes that just touch don't collide; for example, after a U-turn of radius r right before a joint, the straight part before the turn runs alongside the joint exactly touching it. The distances are found by sampling each center line every `fineDistanceThreshold` = 0.001 and then refining the sampled minima along the true center line (`refinedMinimum`), so the sampling can't hide a shallow overlap. For each counted pair, `collisionError` adds a logistic function of the distance centered on that threshold: 0.5 there, and about 1 at an overlap of 0.001 or more. So every counted collision adds at least half of `collisionErrorWeight` (set to `failurePenalty`), and in practice the full weight, however deep it is.

Which pairs are checked: every real joint against every other real joint and every link, and links against each other unless they're in the same **rigid part**. A rigid part is everything from a real joint's distal end to the next real joints, through any waypoints (`findLinkClusters` groups links by their nearest real-joint ancestor). At a branch point it's one physical piece that starts in the joint's distal connection and splits into rigidly connected tubes, which may overlap each other anywhere. A real joint isn't checked against the stretch of its incoming path within π·r of its proximal end (`incomingWindowCheckedLengths`), where the links before an empty or short incoming link meet it; the rest of that path, and its siblings, still are.

Two consequences:
- The penalty gives PSO no gradient to follow out of a collision. It only makes colliding placements look bad.
- With the default `power = 1`, the visits minimize plain link length, which is what the experiments report (total length). `power = 2` minimizes **squared** link lengths instead, which pushes waypoints toward the middle of their links to even out incoming and outgoing lengths. That helped before parent-locking; whether it still does is open. At power 2 a visit can lower its own objective and still raise total length.

## Options

### `optimizeTree` arguments

| Argument | Values | Notes |
|---|---|---|
| `traversal` | `"dfs"`, `"bfs"`, `"randomized"` | See [Traversals](#traversals) |
| `direction` | `"outward"`, `"inward"` | dfs/bfs only. Outward is root → leaves; inward is the reverse of the outward order |
| `orderBy` | `"default"`, `"longest"`, `"shortest"` | dfs/bfs only. The order in which to visit siblings at each branch. **Has no effect on chains** |
| `repeatTraversal` | int, or `"n"` | Number of passes. `"n"` means one pass per entry in `tree.Joints` (joints and waypoints together) |
| `configurations` | list of state lists | Each state list has an entry for every joint and waypoint (0 for waypoints). The result must be collision-free in all of them. `None` means a single all-zeros configuration |
| `childFraction` | float ≥ 0 | Weight on outgoing links in each visit's objective, and in the randomized traversal's weights (default 1) |
| `power` | float ≥ 0 | Power on link lengths in each visit's objective: 1 = plain length, 2 = squared lengths (default 1) |
| `traversalPower` | float ≥ 0 | Power on link lengths in the randomized traversal's weights only (default 2) |
| `waypointSampling` | `"box"`, `"multiscale"` | How waypoint visits sample random initial particles (default `"box"`); see [One waypoint visit](#one-waypoint-visit-optimizewaypointplacementonce) |
| `straightThroughGuess` | bool | Add two waypoint particles pointing along the line between the neighboring parts (default `True` since 2026-11-22; was `False`) |
| `alignmentWeight` | float ≥ 0 | Weight of the waypoint objective's alignment penalty, per radian (default 10; 0 turns it off) |
| `parentLocking` | bool | Let waypoints be parent-locked onto their parent's distal frame, as if deleted but kept available (default: on when `guarantee=True`, which it needs). See [Parent-locking](#parent-locking) |
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

The experiment scripts use 9 variants: {DFS, BFS} × {outward, inward} × {longest, shortest}, plus randomized with `traversalPower=3`.

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
python randomChain.py --note "why"   # jointCount, seed, cubeSize, chainCount; the traversals are the `versions` list in test()
python randomTree.py --note "why"    # jointCount, probabilityOfBranching, seed, treeCount; the `optimizations` list in testRandomTrees()
```

**Warning:** both scripts create their results folder as soon as they're imported, not just when run. So don't import from them; copy the function you need instead.

Each script:
1. Generates a random tree (seeded).
2. Picks collision-free configurations at random (`find_collision_free_configs`).
3. Runs each traversal variant on that tree, resetting the random seed before each one.
4. Plots time vs. length for every variant against `totalLengthLowerBound()`.

`randomChain.py` uses `repeatTraversal="n"`. `randomTree.py` uses the default of 1 pass.

**How long it takes:** an 8-joint gimbal chain has 23 entries in `tree.Joints` (joints plus waypoints). At `"n"` passes that's 529 visits, about 65–70 minutes per traversal. A visit that needs the penalized retry takes roughly 20 times longer than one that doesn't.

### Parameter sweeps

[parameterSweep.py](parameterSweep.py) runs every combination of `waypointSampling`, `straightThroughGuess`, `alignmentWeight`, visit `power` and `parentLocking` (defaults first) on small random gimbal chains (default: 4 joints, seeds 0–7, 5 configurations each), with DFS outward for one pass per entry in `tree.Joints`. Runs go one seed at a time, so stopping early still leaves complete comparisons for the finished seeds.

```
python parameterSweep.py --note "why" [--joints 4] [--seeds 8] [--first-seed 0] [--settings PATTERN,...] [--grid KEY=V1,V2 ...]
python parameterSweep.py --resume "Trials Before Experiments/<sweep folder>" --note "why" [--settings PATTERN,...]
python parameterSweep.py --summarize "Trials Before Experiments/<sweep folder>"
```

`--note` saves why the sweep (or resume) is being run in `run_info.json` (or `resume_<timestamp>.json`); every script warns if it's missing. Log the experiment in [EXPERIMENTS.md](EXPERIMENTS.md) too. `--settings` takes settings labels or wildcard patterns, such as `'*_power2'` or `'multiscale_straight_*'`. With `--resume`, it narrows which of the sweep's unfinished runs to do this time. Each resume writes `resume_<timestamp>.json` with the git commit and the settings it ran. A run that was interrupted partway is cleared and rerun from the start.

`--grid` replaces one setting's values for a new sweep (repeatable), for example `--grid power=1 --grid parentLocking=False,True`. The first combination is the baseline `--summarize` compares against. The grid is saved in `run_info.json`, and `--resume` and `--summarize` use it. Labels end in `_parentlock` for runs with `parentLocking`. `summary.csv` also counts each run's parent-locks and unlocks and its parent-locked waypoints at the end.

It writes `run_info.json`, `summary.csv` (one row per run) and `errors.log` to the sweep folder, and each run's `visits.csv`, final tree and checkpoints to `seed<k>/<settings label>/`. A 4-joint run takes about 4–6 minutes, so the full 16-setting, 8-seed grid takes about 11 hours.

**Overnight runs on a Windows laptop:**
- **Pause Windows Update** (Settings → Windows Update → Pause updates). An update restart at 06:29 killed the first sweep; sleep settings don't prevent update restarts.
- **Set the screen timeout to Never while plugged in.** On laptops with Modern Standby, the screen timing out puts the system in a low-power state even when sleep is set to Never. The first sweep kept running in it, but runs took about 15–25% longer, which affects the timing columns.
- **Keep the lid open.**

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
| `parent_lock_action` | With `parentLocking`: `parentLocked`, `unlocked` or `keptLocked` (a locked waypoint stayed locked); blank otherwise, or if the visit was reverted |
| `parent_locked_count` | The number of parent-locked waypoints after the visit |

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

This section is a snapshot of where things stand. The experiments behind it, with their questions, commands and results, are logged in [EXPERIMENTS.md](EXPERIMENTS.md); add an entry there for each new experiment.

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
  - **Joint sampling.** Particles were clipped to ±2× the incoming link length, which can be about 0 once optimization has collapsed a link. That collapsed the whole swarm onto one point.
  - **Waypoint moves.** These used to be in the world frame, so rotations swung waypoints around the world origin.
  - **Checkpoints.** Names now sort in visit order.
  - **Logging.** Per-visit `visits.csv` and per-run `run_info.json`.
  - **Visit objective.** `power` and `childFraction` now reach each visit's objective; the randomized traversal's weights have their own `traversalPower`.
  - **Bunching.** Investigated waypoints left at orientations where both adjacent links bend, when turning (and maybe slightly moving) the waypoint would clearly straighten them. Added three waypoint search options (`waypointSampling`, `straightThroughGuess`, `alignmentWeight`), five saved examples with a test harness, and a parameter sweep.

### Current status

**Defaults (decided 2026-11-24): power 1 with parent-locking**, together with box sampling, the straight-through guess and alignment weight 10. All of these will be revisited in the final tests.

**Results before 2026-11-23 are unreliable.** Until then, the collision checker never compared links entering real joints against each other, carried parent-locked waypoints weren't checked, and penalized joint retries scored their reversed-axis option on the wrong geometry (see [Known issues](#known-issues)). Their history is in [EXPERIMENTS.md](EXPERIMENTS.md); the waypoint-option results there (the straight-through guess, multiscale sampling, alignment weight) and the 8-joint seed-44 runs need rerunning.

**Parent-locking sweeps with the fixes (`c1b79125`; DFS outward, one pass per entry in `tree.Joints`, box sampling, straight-through guess, alignment weight 10).** No errors, no final tree collides, 2 reverts in 38 runs.

4-joint chains, 8 seeds (`2026.11.23_19.01.14_Sweep_Joints4_Seeds0-7`), final total length relative to power 2 without parent-locking on the same chain:

| Settings | Geometric mean ratio | Better on |
|---|---|---|
| Power 1 + parent-locking | **0.467** | 8 of 8 |
| Power 2 + parent-locking | 0.954 | 5 of 8 |
| Power 1 | 1.788 | 0 of 8 |

6-joint chains, 3 seeds, parent-locking only (`2026.11.23_23.04.06_Sweep_Joints6_Seeds0-2`):

| Chain | Power 1 + parent-locking | Power 2 + parent-locking | Lower bound |
|---|---|---|---|
| 0 | **53.1** | 90.5 | 15.4 |
| 1 | **44.8** | 108.1 | 22.3 |
| 2 | **59.1** | 148.5 | 17.8 |

- **Power 1 + parent-locking is shortest on all 11 chains:** 1.8–2.8× the lower bound on 4 joints, 2.0–3.4× on 6 joints. It ends with essentially every waypoint parent-locked (7 of 7; 10–11 of 11).
- **Power 1 without parent-locking gets stuck:** its waypoints collapse onto their parents too (5–7 of 7), but joints can't move with them.
- **At power 2, parent-locking waypoints lock and unlock repeatedly** (3–13 unlocks per run) and end with 1–4 locked.
- **Run time:** 6.7–9.2 minutes per 4-joint run; 27.5 (power 2) and 38.9 (power 1) minutes per 6-joint run with parent-locking, mostly penalized retries.
- **Not yet rerun with the fixes:** power 2 without parent-locking on 6 joints (earlier: 93.1 / 103.1 / 121.9, which pass the fixed checker but ran with the old checks).

### Known issues

- **Fixed 2026-11-23: links entering real joints were never checked against each other.** `findLinkClusters` left the link from the last waypoint into each real joint out of every cluster, and `buildCollisionMatrices` treated all such links as one cluster. Each cluster now includes the link into the next joint, and any link left out of every cluster gets its own. Results from before the fix are unreliable.
- **Fixed 2026-11-23: `linkLoss` with the collision penalty left the tree in the last configuration.** The joint visit's objective then scored the reversed-axis option on that posed tree, so penalized joint retries could pick a placement that collides once rebuilt in the neutral pose (the revert guard caught it). `linkLoss` now restores the tree's configuration.
- **Fixed 2026-11-23: with parent-locking, carried waypoints' links weren't collision-checked.** The penalty and the revert guard now check the visited node plus the parent-locked waypoints carried with it.
- **Fixed 2026-11-24: collision count and penalty disagreed about near-contacts.** Pairs closer than 0.001 counted as collisions, but the penalty was centered at distance 0, so a near-contact between 0 and 0.001 added almost nothing (1.4 × 10⁻⁶ in a traced case at 0.000976), and the post-visit check then reverted the visit. First fixed by centering the penalty at 0.001; then (same day) replaced by counting only real overlaps, with refined distances, and centering the penalty there. See [Objective](#objective).
- **Fixed 2026-11-24: the tubes of a branching link were checked against each other.** Every pair of links leaving the same branch point counted as colliding, since they start from the same frame in the same direction. Links are now grouped into rigid parts by their nearest real-joint ancestor (see [Objective](#objective)). The GUI branch's fix (`af6b0359`) merged only the waypoint siblings' clusters, which missed the direct link to a branch point's first child.
- **Added 2026-11-24: a window along each real joint's incoming path.** A real joint isn't checked against the links within π·r before its proximal end, since with an empty or short incoming link the link before it meets the joint there.
- **Open question: end-cap false positives with only one closest point at an end.** The collision check models tubes with rounded end caps and filters out overlaps that lie only between two caps (the disc test), but not an overlap between one tube's cap and the middle of another. Such a false positive set off the failure cascade in the first 6-joint parent-locking run (EXPERIMENTS.md, 2026-11-23), and parent-locking makes them likelier, since links pass closer to joints. Filtering them without much extra computation or false negatives in edge cases is hard (whiteboard discussions with Samhitha and others); left open for now.
- **Empty links.** See [Empty links](#empty-links).
- **Fixed 2026-11-24: the root joint had a broken initial guess.** The root has no parent (`Parents[0] = -1`), so the joint visit computed its line-up guess from `Joints[-1]`, the last joint. The root now gets no line-up guess: its particles are the current placement and random ones. The randomized traversal never visits the root.
- **Resuming can change the traversal order.** `optimizeTree` computes the order from the starting tree. `resumeFromCheckpoint` recomputes it from the checkpoint's tree, so `longest`/`shortest` orderings can differ after a resume.
- **The collision penalty is all-or-nothing** (see [Objective](#objective)), so PSO has no gradient to follow out of a collision.
- **Runs aren't exactly reproducible.** Earlier 3-joint runs with the same seed ended between 18.1 and 50.6, partly because the code changed between runs. Comparing variants needs several trials per setting.
- **Merged `chainsGUIdevelopment` (2026-11-25, branch `merge-gui-into-pso`).** The optimizer now runs on the GUI branch's version of the shared classes (end effectors, tip variants, link splitting and waypoint insertion, `PrintedTube` changes), with this branch's collision methods kept. Leftovers from the GUI branch:
  - `KinematicTree.export3DKinematicTree` still uses the OpenSCAD files in `scad/`, which the GUI branch deleted (the manifold-based export replaces it), so it no longer works.
  - `style.py` now makes surfaces opaque by default (`surfaceOpacityDefault` 1.0, was 0.5), so matplotlib plots of trees, including `viewSweepSetting.py` and `viewSweepChain.py`, show solid tubes.
- **Open question: the pair test looks only at the closest pair of points.** If that pair is dismissed by the end-cap filter (where two parts meet), a real contact elsewhere along the same two parts is hidden. The incoming-path window removes this for a joint and the links before it (where it hid a real contact in one saved tree); it remains for a joint and its outgoing links, and for neighboring links in different rigid parts. Related to the end-cap question above.
- **Hard-coded paths.** [manageCheckpoints.py](manageCheckpoints.py) and [view.py](view.py) contain `/home/samhitha/...` paths.

### Empty links

Links of length 0, where a waypoint sits exactly on a neighbor's frame, may fundamentally break both the sampling and the one-joint-at-a-time approach:
- **The optimization creates them.** The 8-joint seed-44 start tree has none. After 23 passes with the default settings, 2 of its 15 waypoints sit on a neighbor (W1 and W22); April's run had 5 by pass 19, plus several more within 0.001. With power 1, about every other link collapses.
- **Visits can't move a waypoint out of that position.** Every rotation in place is infeasible: in the bunching survey, 59% of W1's random placements were infeasible. Escaping needs a move away combined with a large turn, which none of the sampling options produce reliably. W1 isn't fixed by any setting.
- **Power 1 drives waypoints into it,** because sliding a waypoint onto its neighbor costs nothing in plain length. It's the main reason power 1 does worse.
- **The objective sees no penalty for it.** A zero-length link and its waypoint add nothing to `linkLoss`, but the waypoint then can't bend the path anymore.

Ideas to investigate:
- **Delete the waypoint** when it gets close enough to a neighbor, simplifying the tree. *Implemented as [parent-locking](#parent-locking), which treats the waypoint as deleted but keeps it available.*
- **Re-center the waypoint** to keep the tree's structure: move it to exactly the middle of the CSC path outward from it, so both of its links have length again. It's unclear what that should mean at a branch, where a waypoint has several outward paths. *Implemented as parent-locking's unlocking guess.* With it, W1 closes 100% of its gap in the bunching examples (linkLoss 22.61 → 11.31, against 11.32 for the heavy search).
- **Avoid creating them:** sample or guess positions away from the neighbor, or penalize very short links.

### Before systematic experiments

1. **Rerun the waypoint options** (`straightThroughGuess`, `waypointSampling`, `alignmentWeight`) with the collision fixes and the new defaults; the earlier results are unreliable.
2. ~~Test parent-locking~~ Done: power 1 + parent-locking is now the default (see [Current status](#current-status)). Still to do: parent-locking at branch points (the unlock guess assumes one outgoing path).
3. ~~Decide the near-contact question~~ Done: the penalty is centered at the collision threshold (see [Known issues](#known-issues)).
4. **Metric: decided, total length.** Sum of squares is only an option for the visits' objective (power 2).
5. **Rerun the 8-joint seed-44 chain** with the best settings from the sweep, on several traversals. So far only one traversal has been rerun, with default settings.
6. ~~Port the clustering fix from the GUI branch, then smoke-test `randomTree.py`.~~ Done 2026-11-24 (rigid parts; `randomTree.py` updated to this branch's joint classes). Still to do: the full merge of the GUI branch (see [Known issues](#known-issues)).
7. Profile:
   - the split between PSO and Nelder-Mead time;
   - joint visits vs. waypoint visits;
   - how collision checking scales with the number of configurations.
8. Fix the experiment list. The parameter space from the 5/7 meeting notes:
   - **Inputs:** number of joints, number of configurations.
   - **Our parameters:** traversal, number of passes (or visits), `childFraction`, `power` for visits (default 1), `parentLocking` (default on), `waypointSampling` (default box), `straightThroughGuess` (default on), `alignmentWeight` (default 10).
   - **PSO parameters:** global vs. local best; number of particles; sampling bounds; initial particles (uniform random vs. evenly spaced, from [TODO_tests.md](TODO_tests.md)); k; w.

   The planned plots compare time per joint against the objective divided by the best objective found by a reference method.
