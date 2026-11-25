# Experiment log

A running record of the experiments on the PSO tree optimizer: why each was run, what was run, what came out, and what we concluded. [PSO.md](PSO.md) describes the code and its *current* status; this file keeps the history behind it.

## How to use this file

- **Add an entry before or right after each experiment,** including small diagnostic ones. Add new entries at the bottom, and don't rewrite old ones. If a later experiment changes a conclusion, say so in the new entry, and add a one-line "Update:" pointer to the old one.
- **Record the question first.** An experiment without a stated question is hard to interpret later.
- **Pass `--note "<the question>"`** to `randomChain.py`, `randomTree.py` and `parameterSweep.py` (including resumes). It's saved in the results folder's `run_info.json`, so the folder records its own purpose even apart from this log.
- **Results folders** are in `Trials Before Experiments/` and mostly aren't in git, so they exist only on the machine that ran them. Each entry names its folder, the commit it ran on, and the command.
- **Planned experiments** can be logged ahead of time with `Status: planned`, then filled in.

Template:

```markdown
## YYYY-MM-DD: <short title>

**Question:** what we wanted to find out, and why.
**Run:** command, commit, results folder, settings.
**Result:** the numbers that matter.
**Learned / decided:** conclusions, and what changes because of them.
```

---

## Spring 2026: debugging the 8-joint crashes (reconstructed)

*Reconstructed in November 2026 from the meeting notes ("Meeting with Samhitha notes.docx"), the commit history and the saved results. It isn't a record kept at the time; Samhitha's own log of the May tests is in Google Drive.*

**Question:** why do longer chains crash with "Collision detected after retrying with penalty"?

**Run:**
- **January:** `randomChain.py` on 2- and 3-joint chains, seed 42, while the SDF collision checking was being integrated.
- **February 13:** a 7-joint chain, one traversal (`2026.02.13_18.40.31_Joints7_Chains1_Seed42`).
- **April 6:** the 8-joint seed-44 chain on several traversals, with the new per-visit checkpoints (`2026.04.06_21.38.31_Joints8_Chains1_Seed44`, committed as `b40e7303`).
- **May:** `LocalBestPSO` vs. `GlobalBestPSO`, and k and w values, on an 8-joint chain with 4 configurations. This code was never pushed.

**Result:**
- **January:** seven 3-joint runs on the same seed ended between 18.1 and 50.6, partly because the code changed between runs.
- **February:** the 7-joint run went from 540.6 to 212.4 in 2,639 s.
- **April:** DFS Outward Longest and Shortest both finished 23 passes at **243.9** (from 470.3, about 67 minutes each). They were identical, because on a chain the two orderings visit the same order. DFS Inward Longest crashed at joint 9 on pass 19, at length 287.2. A BFS Outward Shortest attempt stopped after one checkpoint, for unknown reasons. The crash tree was saved in `viewBrokenExample.py`.
- **May, from the notes:**
  - GlobalBest found a better result than LocalBest with reverting.
  - With LocalBest, lowering k from 4 to 2 helped a little. Raising w from 0.5 to 0.9 helped a lot but caused more reverts, all on one joint.
  - Even with n passes, the results weren't "smoothing out across waypoints."

**Learned / decided (at the time):**
- **3/23:** `collisionErrorWeight` was defaulting to 1; the default was removed (`2f2dc2cf`). A remaining puzzle: the retry seemed to run without its penalty.
- **4/15:** concluded that PSO wasn't keeping the best particle.
- **5/7:** added a revert when a retry still collides, and planned to revert whenever the objective rises.

Update: the 2026-11-20 entry below finds a different cause for the crash, and explains both the 3/23 and 4/15 puzzles.

---

## 2026-11-20: Why did the April 8-joint run crash at joint 9?

**Question:** is the joint-9 collision real, and why couldn't the penalized retry avoid it?

**Run:** replayed the crash tree (`viewBrokenExample.py`) and the April checkpoint before it (`iter_18_10.dill`), at commit `b40e7303` with diagnostic scripts:
- compared this branch's `collision()` with the GUI branch's version;
- checked the tree in each of the 5 configurations;
- checked whether pyswarms keeps each particle's best position.

**Result:**
- **The collision is real but only exists in the tree's current state.** Root joint 0 and link 9 come within 1.60 of each other at their centerlines, against the 2.0 they need. But the tree is collision-free in all 5 configurations.
- **The current state isn't one of the configurations.** The gimbal construction leaves prismatic joints fully extended (state 2.0), and no configuration has that.
- **The checks disagreed.** The penalized retry worked: its result is collision-free in every configuration, and its penalized loss of 146.0 equals the unpenalized one. The final checks then tested the current state and raised.
- **The GUI branch's `collision()` gives the same result,** so its end-cap fixes aren't involved.
- **pyswarms 1.3.0 does keep each particle's best,** and the current placement is always a starting particle, so a visit can't come back worse.
- **Explains the 3/23 puzzle:** the loss of 146.049 had no penalty because nothing collided in the optimization's configurations.
- **Explains the 4/15 puzzle:** resuming from `iter_18_22` replayed joints 21 → 10, because checkpoint names sorted by joint index and not by visit order. For an inward traversal, `iter_18_22` was the *first* visit of pass 18.

**Learned / decided:**
- Only the given configurations need to be collision-free (decided by Daniel). All collision checks now use the configuration list.
- Visits that still collide, or that raise their own objective, are reverted with a warning instead of raising.
- Checkpoint names now include a visit counter.
- Also fixed in the same commit:
  - joint sampling: particles were clipped to ±2× the incoming link length instead of the bounds;
  - waypoint moves: they were in the world frame, so rotations pivoted about the world origin.

All in `3570f996`.

---

## 2026-11-20: Small-chain regression baseline

**Question:** do the fixes run from start to finish, including resuming? This also gives a quick baseline to check later code changes against.

**Run:** a 4-joint random chain (`np.random.seed(42)`, same generator as `randomChain.py`, 4 configurations), 2 passes of DFS Inward Longest, then a resume from the middle checkpoint. Scratch script, not in the repo; its recipe is in [PSO.md](PSO.md) under Tests.

**Result:** length 77.26 → **66.53**, no collisions in any configuration, and the resume completed. Rerun after each later code change: with default settings it reproduced 66.53 exactly, with every visit's losses and lengths identical. With the new waypoint options (multiscale, straight-through guesses, alignment weight 0), it reached 56.60.

**Learned / decided:** use it as a regression check; defaults must reproduce 66.53.

---

## 2026-11-20: Rerun the 8-joint seed-44 chain with the fixes

**Question:** does the fixed code get through the chain that crashed in April, and how does it compare with April's best traversal?

**Run:** `python randomChain.py` with one traversal, BFS Outward Shortest, default settings, at commit `d3951070`. Folder `2026.11.20_16.15.57_Joints8_Chains1_Seed44`. It was stopped after visit 228 and resumed on 2026-11-21 with `resumeFromCheckpoint`. On a chain, BFS Outward visits the same order as April's DFS Outward runs.

**Result:**
- **Total length 470.3 → 242.8** after 23 passes, vs. April's 243.9. Sum of squared link lengths 3,329, vs. 4,497 for April.
- **Much faster early on,** for example 294.5 vs. 359.3 after pass 6. It then slowed to 254.4 by pass 14 and 242.8 by pass 23, while April dropped sharply after pass 15.
- **No penalized retries or reverts in 529 visits.** Time 4,050 s, vs. 4,042 s in April on a different machine.
- **The final tree is still visibly bunched,** for example in a loop near the base. Collisions don't explain it: with `viewConfigurations.py`, no two unattached parts come within 4.4 of each other in any configuration.

**Learned / decided:** the crash is fixed, and the result matches April's best. The bunching is an optimizer problem, not a collision constraint, so investigate it next.

---

## 2026-11-21: How much better could each waypoint be placed? (bunching survey)

**Question:** in the final tree above, how much could each waypoint's adjacent links shrink if it were placed better? And is the visits' failure to find that due to sampling, or to the alignment penalty?

**Run:** for each of the 15 waypoints, 3,000 random moves, then Nelder-Mead from the 4 best. A third of the moves rotate in place, a third make small moves, and a third use the optimizer's full range. Each is scored under the waypoint objective (alignment weight 10) and under plain `linkLoss`. Script: [debugExamples/bunching/surveyBunching.py](debugExamples/bunching/surveyBunching.py), at commit `d3951070` plus the then-uncommitted `power`/`childFraction` wiring.

**Result:**
- **Under the optimizer's own objective, 13 of the 14 non-leaf waypoints could be improved,** and every improvement is collision-free. For example, W4 157.7 → 114.2, W8 121.0 → 80.5, W16 721 → 545, W19 813 → 654.
- **The better placements are mostly 60–140° turns with moves of 1–9 units.**
- **For some waypoints, the alignment penalty blocks the straightening.** W1's `linkLoss` could halve (22.6 → 11.3) with a 119° turn, but the penalty makes that worse overall. W2, W7, W13 and W14 show smaller versions of the same thing.

**Learned / decided:**
- **Sampling hypothesis:** box sampling almost never pairs a large rotation with a small move, about 1 in 10,000 per particle.
- **Penalty hypothesis:** the alignment penalty can hold a waypoint in place.
- Saved five waypoints (W1, W4, W7, W8, W19) with their best-found placements as test cases in [debugExamples/bunching/](debugExamples/bunching/).
- Implemented three options: `waypointSampling="multiscale"`, `straightThroughGuess`, and `alignmentWeight`.

---

## 2026-11-21: Which waypoint options fix the bunching examples?

**Question:** on the five saved examples, how much does each option, and each combination, close the gap to the best-found placement?

**Run:** `python bunchingExamples.py [--sampling multiscale] [--straight] [--alignment 0]`, 3 seeds per example. On the code committed as `9282d653`. The measure is the share of the gap in local squared link length closed, toward the best placement found for plain `linkLoss`.

**Result:**

| Settings | W1 | W4 | W7 | W8 | W19 | Overall |
|---|---|---|---|---|---|---|
| Default | 0% | 14% | 0% | 9% | 0% | 5% |
| Alignment 0 | 0% | 17% | 0% | 9% | 0% | 5% |
| Multiscale | 0% | 50% | 0% | 31% | 49% | 26% |
| Straight-through | 0% | 82% | 81% | 66% | 79% | 62% |
| Multiscale + straight-through | 0% | 79% | 81% | 61% | 82% | 61% |
| All three | 0% | 83% | 89% | 70% | 82% | 65% |

**Learned / decided:**
- **The straight-through guess does most of the work,** and nearly identically across seeds, since the guesses are deterministic.
- **Multiscale helps alone but adds little on top of the straight-through guess.**
- **Alignment 0 matters only together with a good guess.**
- **W1 isn't fixed by any setting.** It sits on its parent's frame (incoming link length 0), and every turn in place is infeasible there.
- Next: test the options on whole-tree optimizations, in the sweep below.

---

## 2026-11-22: Overnight sweep over the waypoint options and visit power

**Question:** on whole optimizations of several small chains, which combination of `waypointSampling`, `straightThroughGuess`, `alignmentWeight` and visit `power` gives the shortest trees? In particular, does power 1, which targets total length directly, beat power 2?

**Run:** `python parameterSweep.py` at commit `667e7808`. Folder `2026.11.22_02.00.54_Sweep_Joints4_Seeds0-7`. The grid is all 16 combinations, on 4-joint random gimbal chains (seeds 0–7), 5 configurations each, DFS outward for one pass per entry in `tree.Joints` (121 visits per run).

**Result:**
- **Stopped after 52 of 128 runs** (seeds 0–2 complete): a Windows Update restart at 06:29 killed WSL. Earlier, at 03:47, the laptop had entered Modern Standby when its screen timed out. That didn't stop the sweep, but later runs were about 15–25% slower.
- **Final total length relative to the default settings on the same chain (3 seeds):**

| Settings (power 2) | Mean ratio | Range |
|---|---|---|
| multiscale + straight-through + alignment 0 | 0.747 | 0.64–0.89 |
| multiscale + straight-through + alignment 10 | 0.771 | 0.64–0.88 |
| box + straight-through + alignment 10 | 0.819 | 0.68–0.94 |
| multiscale + alignment 0 | 0.867 | 0.84–0.89 |
| box + straight-through + alignment 0 | 0.869 | 0.81–0.93 |
| multiscale alone | 0.943 | 0.88–0.98 |
| alignment 0 alone | 0.999 | 0.84–1.13 |

- **All 26 power-1 runs finished 12–76% longer than the default settings.** On seed 0, all 8 power-1 settings gave identical results, with link lengths `[0, 0, 12.2, 0, 0, 16.2, 0, 0, 24.2, 0]`.
- **No final tree collides.**
- **10 visits were reverted, each with a collision flagged and the loss unchanged.** Traced one: links 5 and 9 at distance 0.000976. The collision count uses a threshold of 0.001, but the smooth penalty is almost 0 at that distance (1.4 × 10⁻⁶).

**Learned / decided:**
- **Same ranking as on the examples:** the straight-through guess helps most, multiscale adds to it, and alignment 0 helps only in combination. The best combination gives trees about 25% shorter.
- **Power 1 fails because of empty links:** visits slide waypoints onto their neighbors, and the waypoints then can't move. A new [Empty links](PSO.md#empty-links) investigation was added to PSO.md, with Daniel's ideas: delete a waypoint that gets close enough to a neighbor, or re-center it on the CSC path outward from it.
- **Decision (Daniel):** finish the sweep with the 8 power-2 settings only, and come back to power 1 after the empty-links investigation. More power-1 runs wouldn't say much about the waypoint settings.
- **The near-contact mismatch is logged in PSO.md's known issues;** whether to count only true overlaps or make the penalty stricter is undecided.
- **For overnight runs:** pause Windows Update, set the screen timeout to Never while plugged in, and keep the lid open.

---

## 2026-11-22: Finish the sweep with power 2 only

**Status:** done (resume ran 17:38–21:59, 4.4 hours; commit `e202b4d6`).

**Question:** with 8 chains instead of 3, does the ranking of the power-2 settings hold, and do the top three separate?

**Run:** `python parameterSweep.py --resume "Trials Before Experiments/2026.11.22_02.00.54_Sweep_Joints4_Seeds0-7" --settings '*_power2' --note "Finish the sweep with power 2 only: does the ranking hold on 8 chains?"`: 38 runs on seeds 3–7.

**Result:** all 64 power-2 runs finished, none failed, and no final tree collides. Final total length relative to the default settings (box, no straight-through, alignment 10) on the same chain:

| Settings (power 2) | Mean ratio | Range | Beat default |
|---|---|---|---|
| multiscale + straight-through + alignment 0 | 0.732 | 0.59–1.05 | 7/8 |
| multiscale + straight-through + alignment 10 | 0.751 | 0.59–1.00 | 8/8 (seed 3 by under 0.5%) |
| box + straight-through + alignment 10 | 0.755 | 0.56–0.94 | 8/8 |
| box + straight-through + alignment 0 | 0.786 | 0.53–1.13 | 7/8 |
| multiscale + alignment 0 | 0.892 | 0.71–1.22 | 7/8 |
| multiscale + alignment 10 | 0.913 | 0.79–1.06 | 7/8 |
| alignment 0 alone | 0.979 | 0.80–1.13 | 5/8 |

Each setting's effect, comparing runs that differ only in that setting (32 pairs each):

| Setting | Geometric mean ratio | Better in |
|---|---|---|
| `straightThroughGuess=True` vs. False | 0.789 | 31/32 |
| `waypointSampling="multiscale"` vs. box | 0.938 | 20/32 |
| `alignmentWeight=0` vs. 10 | 0.987 | 18/32 |

- The best setting varies by seed, but it always uses the straight-through guess (4 multiscale, 4 box).
- Seed 3 is the exception: nothing beat the default by more than 6%.
- The best result on each chain is still 3.1–5.9 times its lower bound.
- 13 visits were reverted, all with the near-contact signature ("did not raise the objective").

**Learned / decided:**
- The ranking held for the straight-through guess: it's the one clear win.
- The top three don't separate. Their means are within 0.025 of each other, and their ranges overlap widely.
- Multiscale sampling helps slightly but not consistently.
- The alignment weight makes no measurable difference on these chains. On the bunching examples it mattered only together with multiscale sampling.
- Decided: `straightThroughGuess=True` is now the default. The sampling and alignment weight stay at box and 10 until they're tested on bigger chains, where there are more waypoints for them to act on.

## 2026-11-22: Pilot of parent-locking at power 1

**Status:** done (ran 23:15–23:39 on uncommitted changes on top of `70870712`; the code is otherwise as committed with parent-locking).

**Question:** does [parent-locking](PSO.md#parent-locking) stop power 1 from getting stuck on empty links, and does it run cleanly?

**Run:** `python parameterSweep.py --seeds 2 --grid power=1 --grid parentLocking=False,True --grid waypointSampling=box --grid alignmentWeight=10 --grid straightThroughGuess=True --note "..."`: folder `2026.11.22_23.15.45_Sweep_Joints4_Seeds0-1`, 4 runs.

**Result:**

| Chain | Power 1 | Power 1 + parent-locking | Best power-2 run in the earlier sweep | Lower bound |
|---|---|---|---|---|
| Seed 0 | 61.5 | **23.6** | 43.5 | 9.2 |
| Seed 1 | 88.6 | **19.6** | 45.2 | 7.7 |

- No errors, and no final tree collides. Seed 1's parent-locking run had 2 reverts, both with the near-contact signature.
- The sum of squared link lengths also fell far below the power-2 runs (80.7 and 47.7, against 178.8 and 167.7), even though those runs minimize it directly.
- In both runs, all 7 waypoints became parent-locked during the first pass, exactly as without parent-locking, and none were unlocked afterward. The first pass ends at the same length as without parent-locking. The gains all come from later passes, where joints can move with their locked waypoints carried along.
- The swarm locked the waypoints by itself, landing on particle 0 (the line-up guess, which is exactly the locked placement). Those visits weren't logged as `parentLocked`; that's fixed now.

**Learned / decided:**
- On these two chains, power 1 with parent-locking ends at about half the best power-2 result. In effect it removes the waypoints and lets the joints place themselves, with direct links between them.
- Whether waypoints help at all here is now an open question: on these chains, none were worth keeping.
- Next: run parent-locking on all 8 chains at both powers, to check the result holds and whether it helps power 2.

## 2026-11-23: Parent-locking at both powers, on 4-joint and 6-joint chains

**Status:** done (4-joint sweep 02:36–06:10, 6-joint sweep 06:10–13:30). Both ran on `5c5b8669`, before the collision-check fix below.

**Question:** does the pilot's power-1 gain hold on all 8 chains, does parent-locking help power 2, and does it carry over to bigger chains?

**Run:** `python parameterSweep.py [--joints 6 --seeds 3] --grid power=2,1 --grid parentLocking=False,True --grid waypointSampling=box --grid alignmentWeight=10 --grid straightThroughGuess=True --note "..."`: folders `2026.11.23_02.36.20_Sweep_Joints4_Seeds0-7` (32 runs) and `2026.11.23_06.10.21_Sweep_Joints6_Seeds0-2` (12 runs).

**Result, 4 joints** (final total length; no errors, and no final tree collides):

| Chain | Power 1 | Power 1 + parent-locking | Power 2 | Power 2 + parent-locking | Lower bound |
|---|---|---|---|---|---|
| 0 | 61.5 | **23.6** | 45.9 | 44.9 | 9.2 |
| 1 | 88.6 | **19.6** | 59.7 | 54.6 | 7.7 |
| 2 | 111.0 | **23.1** | 57.4 | 50.3 | 10.0 |
| 3 | 86.1 | **27.4** | 46.2 | 45.2 | 12.0 |
| 4 | 82.5 | **25.0** | 44.5 | 36.6 | 13.6 |
| 5 | 88.1 | **23.4** | 45.9 | 50.6 | 11.5 |
| 6 | 76.3 | **20.7** | 55.0 | 43.9 | 7.8 |
| 7 | 133.5 | **32.4** | 54.5 | 54.9 | 11.6 |

- Power 1 + parent-locking wins on every chain, at about half the best power-2 length, and ends with all 7 waypoints parent-locked every time. It had 13 reverts in total.
- At power 2, parent-locking is shorter on 6 of 8 chains, but keeps only 1–2 waypoints locked at the end, as the power-2 math predicts.
- Runs took 5.5–7.5 minutes each.

**Result, 6 joints** (final total length; ! marks a final tree the checker calls colliding, from the bug below; 11 waypoints per chain):

| Chain | Power 1 | Power 1 + parent-locking | Power 2 | Power 2 + parent-locking | Lower bound |
|---|---|---|---|---|---|
| 0 | 175.0 | 79.4 ! | 93.1 | 90.5 | 15.4 |
| 1 | 210.0 | 51.6 ! | 103.1 | 103.2 | 22.3 |
| 2 | 296.7 | **42.7** | 121.9 | 99.4 | 17.8 |

- Chain 2's power-1 + parent-locking run is collision-free (0 reverts), so its 42.7 stands: less than half the best power-2 result, as on 4 joints. Chains 0 and 1 ended colliding after 46 and 34 reverts.
- Power 1 without parent-locking gets stuck again, at 1.7–2.4× power 2, with 9–11 of 11 waypoints collapsed onto their parents but not carried.
- Power 2 + parent-locking: shorter on chain 2 (99.4 vs 121.9), about even on chains 0 and 1. It unlocked waypoints 8–11 times per run (5–9 per run on 4 joints) and ended with 1–2 locked, so at power 2 waypoints lock and unlock repeatedly.
- Mean run time: 19.5 / 45.0 / 26.1 / 29.5 minutes; about 7.3 hours for the 6-joint sweep.

**Bug found: the 6-joint parent-locking results are invalid.** Chain 0's power-1 + parent-locking run (79.4) ended with a tree the checker calls colliding, after 46 reverts in 69 minutes. Traced through its checkpoints:
- At visit 61, joint 9 moved with W10 and W11 carried along, which moved link 12 (leaving W11) into joint 0 in configuration 4. Both the collision penalty and the revert guard only checked pairs involving the visited node and its own links, so the pair joint 0 / link 12 was never checked. A deep collision of a carried link would have slipped through the same way.
- The collision itself is an end-cap false positive of the kind accepted so far: the closest point on joint 0 is the end of its centerline, and the overlap is −0.187 against its rounded cap but +0.335 against its flat end. It isn't filtered because only one of the two closest points is at an end.
- From then on, the revert guard rejected every change that still collided, even those that reduced the collisions, so the tree stayed in that state.

**Fixed in `6ed3f12c`:** collision checks and the collision penalty now cover the visited node plus its carried parent-locked waypoints. The revert guard accepts a change that still collides only in configurations the tree already collided in, when the objective didn't increase. Replaying visit 61 from the checkpoint before it with the fix gives a collision-free result.

**Learned / decided:**
- The 4-joint lengths stand, since none of those final trees collide, though a rerun with the fix may differ slightly. The 6-joint parent-locking runs, at both powers, need rerunning with the fix.
- Open: with parent-locking, links pass closer to joints, so end-cap false positives with only one closest point at an end may come up more often. Worth considering extending the end-cap filter to that case.
- Next: rerun the parent-locking settings on 6 joints (and possibly 4) with the fix.

## 2026-11-23: Collision checker never compared links entering real joints: all earlier results are unreliable

**Status:** fixed in the commit after `b142a961`.

**Question:** why does the final tree of chain 2's 6-joint power-1 + parent-locking run visibly collide, when the checker reported it collision-free?

**Found:** `KinematicTree.findLinkClusters` groups each run of waypoints with its links, and link/link pairs within a cluster aren't checked. It only followed waypoint children, so the link from the last waypoint into each real joint (links 3, 6, 9, 12, 15 on a 6-joint chain) was in no cluster. `buildCollisionMatrices` then saw those links as all in the same cluster (`None == None`), so **no two links entering real joints were ever checked against each other**, in any tree, in any run, by the penalty, the revert guard or the final check. On chain 2, links 6 and 15 overlap by 0.79 tube radii. The intended clusters (as in [sets.py](sets.py)'s docstring) include that link: `[1, 2, 3], [4, 5, 6], …`.

This is a different bug from the GUI branch's clustering fix `af6b0359` (sibling waypoints at a branch point), which still needs porting before running trees.

**Fix:** follow real-joint children too, so each cluster includes the link into the next joint; any link still left out of every cluster gets a cluster of its own.

**Recheck of all 138 final trees in the four sweeps with the fixed checker** (all configurations of each chain):

| Settings | Final trees that collide |
|---|---|
| Power 2, all settings, all sweeps | 1 of 99 |
| Power 1, no parent-locking | 18 of 37 |
| Power 1 + parent-locking, 4 joints | 4 of 10 |
| Power 1 + parent-locking, 6 joints | 3 of 3 |

**Learned / decided:**
- **All results recorded before this fix are unreliable**, every entry above included. The optimizer never saw these collisions, so every run, at either power, may have taken moves it shouldn't have, even where its final tree happens to pass the fixed check. Power 1 is the worst affected: once waypoints collapse, the links into joints are the only links with length, and none of them were compared. The power-1 + parent-locking gains in particular (about half the power-2 length) may come largely from trees that pass through themselves.
- The decisions based on those results (`straightThroughGuess=True` as the default; power 2 only) stand for now as working choices, but need confirming.
- Next: small runs to check the fixed code runs smoothly, then rerun the sweeps.

## 2026-11-23: Small runs after the collision fixes

**Status:** done.

**Question:** does the optimizer still run smoothly after the link-clustering and carried-waypoint fixes?

**Run:** `python parameterSweep.py --joints 3 --seeds 2 --grid power=2,1 --grid parentLocking=False,True --grid waypointSampling=box --grid alignmentWeight=10 --grid straightThroughGuess=True --note "..."`, twice: folder `2026.11.23_14.00.49_Sweep_Joints3_Seeds0-1` on `be0e9422`, then `2026.11.23_14.26.14_Sweep_Joints3_Seeds0-1` with the `linkLoss` fix below (uncommitted at the time, committed right after).

**Result:**

| Chain | Power 2 | Power 2 + parent-locking | Power 1 | Power 1 + parent-locking | Lower bound |
|---|---|---|---|---|---|
| 0 | 28.1 | 25.1 | 36.9 | 15.3 | 6.5 |
| 1 | 31.8 | 30.2 | 46.5 | 14.3 → 14.4 | 6.5 |

- Both times: no errors, no final tree collides under the fixed checker, and runs took 2.2–4.9 minutes.
- **The first run had 2 reverts** (chain 1, power 1 + parent-locking, joint 3), where the objective with the collision penalty went down although the result really collided (overlap 0.37). Cause, older than parent-locking: `linkLoss` posed the tree in each configuration for the penalty and left it in the last one. The joint visit's objective then reversed the joint's axis on that posed tree and scored it there, but the returned tree is rebuilt in the neutral pose, where it collided. Fixed by restoring the configuration in `linkLoss`. Replaying that visit with the fix: accepted, loss 10.12 → 9.64, collision-free.
- **The second run had no reverts**; chain 1's run kept all 5 penalized retries. The other 7 runs gave identical lengths (none had a penalized joint retry that chose the reversed axis).

**Learned / decided:**
- The fixed code runs cleanly. Every penalized joint retry before this fix could have chosen between its two axis options on the wrong geometry, one more reason the earlier results are unreliable.
- On 3-joint chains, power 1 + parent-locking stays at about half the power-2 length with the fixed checks. Too small to conclude anything; next is rerunning the 4-joint sweep and the 6-joint parent-locking settings.

## 2026-11-23: Rerun of the parent-locking sweeps with all three collision fixes

**Status:** done (4-joint sweep 19:01–23:04; 6-joint sweep started 23:04, cut off by the computer sleeping before its first run finished, resumed 23:17). Both on `c1b79125`.

**Question:** with the collision checks fixed, does power 1 + parent-locking still beat power 2, and does it hold on 6 joints?

**Run:** as in the earlier 2026-11-23 entry, except the 6-joint sweep ran only `parentLocking=True` (power 2 and power 1): folders `2026.11.23_19.01.14_Sweep_Joints4_Seeds0-7` (32 runs) and `2026.11.23_23.04.06_Sweep_Joints6_Seeds0-2` (6 runs).

**Result, 4 joints** (final total length; no errors, no final tree collides, 1 revert in all 32 runs):

| Chain | Power 1 | Power 1 + parent-locking | Power 2 | Power 2 + parent-locking | Lower bound |
|---|---|---|---|---|---|
| 0 | 61.5 | **22.2** | 45.9 | 44.9 | 9.2 |
| 1 | 93.9 | **19.0** | 59.7 | 55.2 | 7.7 |
| 2 | 119.9 | **23.6** | 57.2 | 60.7 | 10.0 |
| 3 | 86.1 | **28.0** | 46.2 | 45.2 | 12.0 |
| 4 | 82.5 | **25.1** | 44.5 | 36.6 | 13.6 |
| 5 | 92.1 | **21.9** | 45.9 | 50.6 | 11.5 |
| 6 | 76.3 | **20.3** | 55.0 | 43.9 | 7.8 |
| 7 | 133.5 | **32.4** | 54.5 | 54.9 | 11.6 |

Geometric mean ratio to power 2: power 1 + parent-locking 0.467 (better on 8 of 8), power 2 + parent-locking 0.954 (5 of 8), power 1 1.788 (0 of 8).

**Result, 6 joints** (no errors, no final tree collides, 1 revert):

| Chain | Power 1 + parent-locking | Power 2 + parent-locking | Lower bound |
|---|---|---|---|
| 0 | **53.1** | 90.5 | 15.4 |
| 1 | **44.8** | 108.1 | 22.3 |
| 2 | **59.1** | 148.5 | 17.8 |

- Before the fixes these were 79.4 / 51.6 / 42.7 (all colliding under the fixed checker) and 90.5 / 103.2 / 99.4.
- Power 1 + parent-locking ended with 10–11 of 11 waypoints locked; power 2 + parent-locking unlocked 8–13 times per run and ended with 2–3 locked.
- Mean run time 38.9 minutes for power 1 + parent-locking (149 penalized retries over 3 runs) and 27.5 for power 2 + parent-locking.

**Learned / decided:**
- **The power-1 + parent-locking result survives the collision fixes.** On 4 joints the lengths barely changed (within about 1.4 either way) and are now verified collision-free: about 47% of the power-2 length, better on every chain. On 6 joints it is 2.0–3.4× the lower bound, against 4.8–8.3× for power 2 + parent-locking.
- With parent-locking at power 1, essentially every waypoint ends locked: on these chains, the best trees found have no waypoints in use.
- Power 2 + parent-locking is unconvincing (5 of 8 on 4 joints, and chain 2 on 6 joints got much worse, 99.4 → 148.5), with waypoints locking and unlocking repeatedly.
- No power-2 baseline without parent-locking on 6 joints was rerun with the fixes; the earlier ones (93.1 / 103.1 / 121.9) passed the fixed checker but ran with the old checks.
