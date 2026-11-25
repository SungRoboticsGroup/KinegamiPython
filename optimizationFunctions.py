from KinematicTree import *
from treeTraversals import *
import pyswarms as ps
import os
import dill
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Any, Tuple
from datetime import datetime
import csv
import json

# Global file path for collision penalty logging (can be set externally)
collision_penalty_log_file = None

@dataclass
class CheckpointState:
    # Core tree state
    tree: Any  # KinematicTree
    
    # Random state for reproducibility
    random_state: Tuple
    
    # Traversal information
    traversal_info: Dict[str, Any] = field(default_factory=dict)
    # Contains: algorithm, direction, orderBy, current_index, repeat_index, 
    #           completed_indices (list of already-processed joint indices)
    
    # Collision detection state
    collision_matrices: Optional[Tuple] = None
    
    # Optimization history
    optimization_history: Dict[str, List] = field(default_factory=dict)
    # Contains: times, lengths, losses
    
    # Configuration
    configurations: Optional[List] = None
    failure_penalty: float = 0.0
    
    # Metadata
    metadata: Dict[str, Any] = field(default_factory=dict)
    # Contains: timestamp, iteration_count, total_joints, version

def saveCheckpoint(state: CheckpointState, filepath: str, compress: bool = False) -> str:
    """
    Save a checkpoint state to disk using dill serialization.
    
    Args:
        state: CheckpointState object containing all optimization state
        filepath: Path to save the checkpoint (without extension)
        compress: If True, use gzip compression (slower but smaller files)
    
    Returns:
        The full path of the saved checkpoint file
    """
    # Ensure directory exists
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    ext = '.dill.gz' if compress else '.dill'
    full_path = filepath + ext
    
    if compress:
        import gzip
        with gzip.open(full_path, 'wb') as f:
            dill.dump(state, f)
    else:
        with open(full_path, 'wb') as f:
            dill.dump(state, f)
    
    return full_path

def loadCheckpoint(filepath: str, restore_random_state: bool = True) -> CheckpointState:
    """
    Load a checkpoint state from disk.
    
    Args:
        filepath: Path to the checkpoint file
        restore_random_state: If True, restore numpy random state from checkpoint
    
    Returns:
        CheckpointState object with all saved optimization state
    """
    if filepath.endswith('.gz'):
        import gzip
        with gzip.open(filepath, 'rb') as f:
            state = dill.load(f)
    else:
        with open(filepath, 'rb') as f:
            state = dill.load(f)
    
    if restore_random_state and state.random_state is not None:
        np.random.set_state(state.random_state)
    
    return state

def inspectCheckpoint(filepath: str, load_tree: bool = False) -> Dict[str, Any]:
    """
    Inspect a checkpoint file and print summary without loading full tree.
    
    Args:
        filepath: Path to the checkpoint file
        load_tree: If True, also load and return the tree object
    
    Returns:
        Dict with checkpoint metadata and optionally the tree
    """
    state = loadCheckpoint(filepath, restore_random_state=False)
    
    info = {
        'traversal_info': state.traversal_info,
        'metadata': state.metadata,
        'failure_penalty': state.failure_penalty,
        'num_configurations': len(state.configurations) if state.configurations else 0,
        'has_collision_matrices': state.collision_matrices is not None,
        'optimization_history': {
            'num_iterations': len(state.optimization_history.get('times', [])),
            'latest_time': state.optimization_history.get('times', [None])[-1],
            'latest_length': state.optimization_history.get('lengths', [None])[-1],
        }
    }
    
    if load_tree:
        info['tree'] = state.tree
    
    # Print summary
    print(f"Checkpoint: {filepath}")
    print(f"  Timestamp: {state.metadata.get('timestamp', 'N/A')}")
    print(f"  Traversal: {state.traversal_info.get('algorithm', 'N/A')} "
          f"({state.traversal_info.get('direction', '')}, {state.traversal_info.get('orderBy', '')})")
    print(f"  Repeat: {state.traversal_info.get('repeat_index', 'N/A')} / "
          f"{state.traversal_info.get('total_repeats', 'N/A')}")
    print(f"  Current joint: {state.traversal_info.get('current_index', 'N/A')}")
    print(f"  Completed joints: {len(state.traversal_info.get('completed_indices', []))}")
    print(f"  Iterations: {info['optimization_history']['num_iterations']}")
    print(f"  Latest length: {info['optimization_history']['latest_length']:.4f}" 
          if info['optimization_history']['latest_length'] else "  Latest length: N/A")
    
    return info

def checkpointName(visit: int, repeat: int, index: int) -> str:
    # The visit counter comes first so names sort in the order the visits happened
    # (joint index order differs from visit order for most traversals, e.g. inward DFS)
    return f'visit{visit:05d}_iter_{repeat}_{index}'

def listCheckpoints(checkpoint_dir: str) -> list:
    """
    List all checkpoint files in a directory, sorted by iteration order.

    Args:
        checkpoint_dir: Path to directory containing checkpoint files

    Returns:
        List of dicts with 'path', 'visit', 'repeat', 'joint_idx' for each checkpoint
    """
    import re

    checkpoints = []

    if not os.path.exists(checkpoint_dir):
        return checkpoints

    for filename in os.listdir(checkpoint_dir):
        if filename.endswith('.dill'):
            match = re.match(r'(?:visit(\d+)_)?iter_(\d+)_(\d+)\.dill', filename)
            if match:
                checkpoints.append({
                    'path': os.path.join(checkpoint_dir, filename),
                    'visit': int(match.group(1)) if match.group(1) else None,
                    'repeat': int(match.group(2)),
                    'joint_idx': int(match.group(3)),
                    'filename': filename
                })

    if any(cp['visit'] is None for cp in checkpoints):
        # Older checkpoints have no visit counter, so within a repeat they can only be sorted
        # by joint index, which is NOT the order they were visited in (e.g. for inward traversals)
        print("Warning: some checkpoints predate visit numbering; their order within a repeat is by joint index, not visit order.")
        checkpoints.sort(key=lambda x: (x['repeat'], x['joint_idx']))
    else:
        checkpoints.sort(key=lambda x: x['visit'])

    print(f"Found {len(checkpoints)} checkpoints in {checkpoint_dir}")
    print(f"Checkpoints:")
    for cp in checkpoints:
        print(f"  {cp['filename']} (repeat {cp['repeat']}, joint {cp['joint_idx']})")
    
    return checkpoints

def resumeFromCheckpoint(checkpoint_path: str, 
                           showSteps: bool = True,
                           verbose: bool = True,
                           directory: Optional[str] = None) -> Any:
    """
    Resume optimization from a saved checkpoint.
    
    Args:
        checkpoint_path: Path to the checkpoint file to resume from
        showSteps: Whether to show visualization after each step
        verbose: Whether to print progress information
        directory: Output directory (if None, uses original from checkpoint)
    
    Returns:
        The optimized KinematicTree (or tuple with times/lengths if evaluate=True in original)
    """
    state = loadCheckpoint(checkpoint_path, restore_random_state=True)
    
    tree = state.tree
    traversal_info = state.traversal_info
    times = state.optimization_history.get('times', [])
    lengths = state.optimization_history.get('lengths', [])
    losses = state.optimization_history.get('losses', [])
    configurations = state.configurations
    failurePenalty = state.failure_penalty
    
    # Extract traversal parameters
    traversal_algo = traversal_info.get('algorithm', 'dfs')
    direction = traversal_info.get('direction', 'outward')
    orderBy = traversal_info.get('orderBy', 'longest')
    repeat_index = traversal_info.get('repeat_index', 0)
    total_repeats = traversal_info.get('total_repeats', 1)
    completed_indices = set(traversal_info.get('completed_indices', []))
    current_index = traversal_info.get('current_index')
    
    # Get other parameters from metadata
    metadata = state.metadata
    guarantee = metadata.get('guarantee', False)
    parallelize = metadata.get('parallelize', False)
    childFraction = metadata.get('childFraction', 1)
    if 'traversalPower' in metadata:
        power = metadata.get('power', 2)
        traversalPower = metadata['traversalPower']
    else:
        # Older checkpoints only stored the randomized traversal's weight power as 'power';
        # their visits always used power 2
        power = 2
        traversalPower = metadata.get('power', 2)
    evaluate = metadata.get('evaluate', False)
    waypointSampling = metadata.get('waypointSampling', "box")
    straightThroughGuess = metadata.get('straightThroughGuess', False)  # checkpoints from before the option existed ran without it
    alignmentWeight = metadata.get('alignmentWeight', 10)
    parentLocking = metadata.get('parentLocking', False)
    save_checkpoints = metadata.get('save_checkpoints', True)
    
    if directory is None:
        directory = metadata.get('directory')
    if directory is not None:
        os.makedirs(directory, exist_ok=True)

    if verbose:
        print(f"Resuming from checkpoint: {checkpoint_path}")
        print(f"  Repeat: {repeat_index + 1}/{total_repeats}")
        print(f"  Completed joints: {len(completed_indices)}")
        print(f"  Current tree length: {tree.totalLength():.4f}")
    
    # Reconstruct traversal
    treeTraversals = {
        "dfs": partial(dfs, direction=direction, orderBy=orderBy),
        "bfs": partial(bfs, direction=direction, orderBy=orderBy),
        "randomized": partial(randomized, 
                            power=traversalPower, 
                            count=len(tree.Joints), 
                            childFraction=childFraction,
                            isWeighted=True)
    }
    
    start = time.time() - (times[-1] if times else 0)  # Adjust start time
    
    def log(t, idx):
        diff = time.time() - start
        times.append(diff)
        lengths.append(t.totalLength())
        if directory is not None:
            save_path = os.path.join(directory, f"{diff}_{idx}")
            t.save(save_path, saveDir=False)
    
    checkpoint_dir = os.path.join(directory, 'checkpoints') if directory else None
    iteration_count = metadata.get('iteration_count', len(completed_indices))
    
    # Continue from where we left off
    for rep in range(repeat_index, total_repeats):
        full_traversal = list(treeTraversals[traversal_algo](tree))
        
        for index in full_traversal:
            # Skip already completed indices in current repeat
            if rep == repeat_index and index in completed_indices:
                continue
            
            if verbose:
                print(f"Optimizing joint index: {index} {tree.Joints[index]}")
            
            iters = 50
            tolerance = tree.r / 10
            
            visitRecord = {}
            visitStart = time.time()
            if isWaypoint(tree.Joints[index]):
                tree, loss = optimizeWaypointPlacement(
                    tree, index, maxiter=iters, tol=tolerance,
                    failurePenalty=failurePenalty, childFraction=childFraction, power=power,
                    ignoreLater=(not guarantee), parallelize=parallelize,
                    verbose=verbose, configurations=configurations, visitRecord=visitRecord,
                    waypointSampling=waypointSampling, straightThroughGuess=straightThroughGuess,
                    alignmentWeight=alignmentWeight, parentLocking=parentLocking)
            else:
                tree, loss = optimizeJointPlacement(
                    tree, index, maxiter=iters, tol=tolerance,
                    failurePenalty=failurePenalty, childFraction=childFraction, power=power,
                    ignoreLater=(not guarantee), parallelize=parallelize,
                    verbose=verbose, configurations=configurations, visitRecord=visitRecord,
                    parentLocking=parentLocking)
            visitTime = time.time() - visitStart

            log(tree, index)
            losses.append(loss)
            completed_indices.add(index)
            iteration_count += 1
            if directory is not None:
                visitRecord.update(visit=iteration_count, repeat=rep, index=index, visit_s=visitTime,
                                   elapsed_s=times[-1], total_length=lengths[-1],
                                   sum_squared_length=sumSquaredLinkLength(tree),
                                   parent_locked_count=countParentLocked(tree))
                appendVisitRecord(directory, visitRecord)
            
            # Save checkpoint
            if save_checkpoints and checkpoint_dir:
                checkpoint_state = CheckpointState(
                    tree=copy.deepcopy(tree),
                    random_state=np.random.get_state(),
                    traversal_info={
                        'algorithm': traversal_algo,
                        'direction': direction,
                        'orderBy': orderBy,
                        'current_index': index,
                        'repeat_index': rep,
                        'total_repeats': total_repeats,
                        'completed_indices': list(completed_indices),
                    },
                    collision_matrices=tree.buildCollisionMatrices(),
                    optimization_history={'times': times.copy(), 'lengths': lengths.copy(), 'losses': losses.copy()},
                    configurations=configurations,
                    failure_penalty=failurePenalty,
                    metadata={
                        'timestamp': datetime.now().isoformat(),
                        'iteration_count': iteration_count,
                        'total_joints': len(tree.Joints),
                        'directory': directory,
                        'guarantee': guarantee,
                        'parallelize': parallelize,
                        'childFraction': childFraction,
                        'power': power,
                        'traversalPower': traversalPower,
                        'waypointSampling': waypointSampling,
                        'straightThroughGuess': straightThroughGuess,
                        'alignmentWeight': alignmentWeight,
                        'parentLocking': parentLocking,
                        'evaluate': evaluate,
                        'save_checkpoints': save_checkpoints,
                        'version': '1.0',
                    }
                )
                saveCheckpoint(checkpoint_state, os.path.join(checkpoint_dir, checkpointName(iteration_count, rep, index)))
        
        # Clear completed indices for next repeat
        completed_indices.clear()
    
    if showSteps:
        tree.show()
    
    if verbose:
        print(f"TOTAL OPTIMIZATION TIME: {time.time() - start}")
    
    if directory is not None:
        tree.save(os.path.join(directory, "final"), saveDir=False)
        
        # Plot and save losses
        plt.figure(figsize=(8, 5))
        plt.plot(times, lengths, marker='o', linestyle='-', label='Length')
        plt.xlabel('Time')
        plt.ylabel('Length')
        plt.title('Optimization Progress (Resumed)')
        plt.legend()
        plt.grid(True)
        
        # Save to same directory as plot_0.png (results_dir, which is 2 levels up from directory)
        results_dir = os.path.dirname(os.path.dirname(directory))
        plot_save_path = os.path.join(results_dir, "lossPlotAfterResuming.png")
        try:
            plt.savefig(plot_save_path, dpi=300, bbox_inches='tight')
            if verbose:
                print(f"Saved plot to: {plot_save_path}")
        except Exception as e:
            print(f"Error saving plot: {e}")
        plt.close()
    
    if evaluate:
        return tree, times, lengths
    
    return tree

def setCollisionPenaltyLogFile(filepath):
    global collision_penalty_log_file
    collision_penalty_log_file = filepath

def logCollisionPenalty(message):
    print(message)
    if collision_penalty_log_file:
        with open(collision_penalty_log_file, 'a') as f:
            f.write(message + '\n')

def collidingConfigurations(tree, index, configurations):
    # Returns the configurations in which the moved joint/waypoint at index (or any of a list of indices)
    # collides with something.
    # Only the given configurations matter: the tree's current state is not checked.
    colliding = []
    for configuration in configurations:
        configured = copy.deepcopy(tree)
        for i in range(len(configured.Joints)):
            if not isWaypoint(configured.Joints[i]):
                configured.setJointState(i, configuration[i])
                configured.Joints[i].recomputeCollisionCapsules()
        if configured.detectCollisions(specificJointIndex=index, debug=False) > 0:
            colliding.append(configuration)
    return colliding

def warnIfInitialCollisions(subject, configurations, caller):
    for configuration in configurations:
        configured = copy.deepcopy(subject)
        for i in range(len(configured.Joints)):
            if not isWaypoint(configured.Joints[i]):
                configured.setJointState(i, configuration[i])
                configured.Joints[i].recomputeCollisionCapsules()
        if configured.detectCollisions(debug=True) > 0:
            print(f"Warning: Initial tree in {caller} contains collisions in configuration {configuration}.")

def localLinkScale(subject, index):
    # Length scale for sampling moves of the joint/waypoint at index: its longest adjacent link,
    # floored at the tube radius so near-zero-length links (which optimization can create by sliding
    # a waypoint onto its neighbor) don't collapse the whole swarm onto one point
    lengths = [subject.Links[index].path.length] + [subject.Links[c].path.length for c in subject.Children[index]]
    return max(max(lengths), subject.r)

def wrapToRange(values, low, high):
    return low + np.mod(values - low, high - low)

def rotationAngle(R):
    # Geodesic angle of a rotation matrix
    return np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))

# Parent-locking (option parentLocking): a waypoint is parent-locked when its proximal Dubins frame coincides
# with its parent's distal Dubins frame, so its incoming link is empty and the tree behaves as if the waypoint
# were deleted, while keeping it available to be moved out again later. Visits carry parent-locked descendants
# along with the joint or waypoint they're locked to, waypoint visits consider locking to the parent, and visits
# to a parent-locked waypoint consider unlocking it. Whether a waypoint is parent-locked is read from the
# geometry, so it needs no extra saved state.
PARENT_LOCK_TOLERANCE = 1e-6  # in units of the tube radius for positions, and directly for rotation matrix entries
PARENT_LOCK_MARGIN = 1e-3     # relative improvement in linkLoss needed to lock or unlock (keeps waypoints from flip-flopping)

def isParentLocked(tree, index):
    if index == 0 or not isWaypoint(tree.Joints[index]):
        return False
    parentFrame = tree.Joints[tree.Parents[index]].DistalDubinsFrame()
    frame = tree.Joints[index].ProximalDubinsFrame()
    return np.linalg.norm(frame.t - parentFrame.t) < PARENT_LOCK_TOLERANCE * tree.r and \
           np.max(np.abs(frame.R - parentFrame.R)) < PARENT_LOCK_TOLERANCE

def parentLockedDescendants(tree, index):
    # The waypoints locked to index, directly or through a chain of other locked waypoints, parents first
    locked = []
    frontier = [index]
    while frontier:
        node = frontier.pop(0)
        for c in tree.Children[node]:
            if isParentLocked(tree, c):
                locked.append(c)
                frontier.append(c)
    return locked

def countParentLocked(tree):
    return sum(isParentLocked(tree, i) for i in range(len(tree.Joints)))

def lockedOutgoingLinks(tree, index, locked):
    # The links leaving index and its parent-locked descendants (other than the empty links to those descendants)
    return [c for node in [index] + locked for c in tree.Children[node] if c not in locked]

def moveWithParentLocked(tree, index, transformation, locked, relative=True):
    # Like tree.transformJoint(index, transformation, propogate=False, safe=True, relative=relative,
    # recomputeBoundingBall=False), but the parent-locked descendants in locked move with the joint, snapped back
    # onto their parents' distal frames, so a move is never judged by the path to a waypoint left behind.
    # Returns False (leaving the tree unchanged) if a link can't be made.
    if not locked:
        return tree.transformJoint(index, transformation, propogate=False, safe=True, relative=relative,
                                   recomputeBoundingBall=False)
    backup = tree.dataDeepCopy()
    try:
        if relative:
            transformation = tree.Joints[index].Pose @ transformation @ tree.Joints[index].Pose.inv()
        tree.Joints[index].transformPoseBy(transformation)
        for d in locked:
            snap = tree.Joints[tree.Parents[d]].DistalDubinsFrame() * tree.Joints[d].ProximalDubinsFrame().inv()
            tree.Joints[d].transformPoseBy(snap)
        moved = [index] + locked
        if index > 0:
            tree.Links[index] = LinkCSC(tree.r, tree.Joints[tree.Parents[index]].DistalDubinsFrame(),
                                        tree.Joints[index].ProximalDubinsFrame(), tree.maxAnglePerElbow)
        else:
            tree.Links[0] = tree.Links[0].newLinkTransformedBy(transformation)
            tree.boundingBall = tree.boundingBall.newBallTransformedBy(transformation)
        for node in moved:
            for c in tree.Children[node]:
                tree.Links[c] = LinkCSC(tree.r, tree.Joints[node].DistalDubinsFrame(),
                                        tree.Joints[c].ProximalDubinsFrame(), tree.maxAnglePerElbow)
            tree.Joints[node].recomputeCollisionCapsules()
    except (ValueError, AssertionError):
        tree.setTo(backup)
        return False
    return True

def pathMidpointFrame(link):
    # Position and unit tangent at the middle (by arc length) of a link's CSC path
    path = link.path
    r = path.r
    lengthC1, lengthS, lengthC2 = r * path.theta1, path.tMag, r * path.theta2
    s = (lengthC1 + lengthS + lengthC2) / 2
    if s <= lengthC1:
        # along the first turn: the tangent turns from startDir toward w1, the direction to the circle's center
        a = s / r
        d, w = path.startDir, path.w1
        return path.circleCenter1 - r * np.cos(a) * w + r * np.sin(a) * d, np.cos(a) * d + np.sin(a) * w
    if s <= lengthC1 + lengthS:
        return path.turn1end + (s - lengthC1) * path.tUnit, path.tUnit
    # along the second turn, measured back from the end
    a = (lengthC1 + lengthS + lengthC2 - s) / r
    d, w = path.endDir, path.w2
    return path.circleCenter2 - r * np.cos(a) * w - r * np.sin(a) * d, np.cos(a) * d - np.sin(a) * w

def unlockGuess(subject, index):
    # For a parent-locked waypoint: the move (in its own frame, as [t, ZYZ angles]) to the middle of the path
    # leaving it, pointing along that path, so both of its links have length again. Splitting a path this way
    # leaves both halves feasible, so it's a safe way out of the locked position.
    if len(subject.Children[index]) == 0:
        return []
    waypoint = subject.Joints[index]
    position, tangent = pathMidpointFrame(subject.Links[subject.Children[index][0]])
    frame = waypoint.ProximalDubinsFrame()
    # Dubins frames point along x; turn the frame as little as possible to point along the tangent
    newFrame = SE3.Rt(rotationBetween(frame.R[:, 0], tangent) @ frame.R, position, check=False)
    move = waypoint.Pose.inv() * (newFrame * frame.inv()) * waypoint.Pose
    return [list(move.t) + list(move.eul())]

def optimizeWithCollisionRetry(subject, index, label, optimizeOnce, includeCollisionPenalty, configurations,
                               visitRecord=None, movedIndices=None):
    # optimizeOnce(includeCollisionPenalty) -> (tree, loss, initialLoss, timing), where initialLoss is the
    # same objective evaluated at the current placement and timing has 'pso_s', 'nm_s' (and 'notApplied'
    # if the best particle couldn't be applied).
    # If visitRecord is a dict, it's filled in with what happened on this visit (see VISIT_FIELDS).
    # movedIndices: the joints/waypoints the visit moves (default [index]; with parent-locking, also the waypoints
    # carried along), whose collisions are checked.
    if movedIndices is None:
        movedIndices = [index]
    if visitRecord is None:
        visitRecord = {}
    visitRecord.update(kind=label, retried=False, reverted=False, revert_reason="", collision_check_s=0.0)

    def timedCollidingConfigurations(tree):
        start = time.time()
        colliding = collidingConfigurations(tree, movedIndices, configurations)
        visitRecord['collision_check_s'] += time.time() - start
        return colliding

    # First try without the collision penalty (much cheaper), and only retry with it if that collides.
    tree, loss, initialLoss, timing = optimizeOnce(includeCollisionPenalty=False)
    visitRecord.update(nopenalty_pso_s=timing['pso_s'], nopenalty_nm_s=timing['nm_s'])
    if includeCollisionPenalty and len(timedCollidingConfigurations(tree)) > 0:
        logCollisionPenalty(f"Collision detected for {label} {index}, retrying with collision penalty")
        visitRecord['retried'] = True
        tree, loss, initialLoss, timing = optimizeOnce(includeCollisionPenalty=True)
        visitRecord.update(penalty_pso_s=timing['pso_s'], penalty_nm_s=timing['nm_s'])
    visitRecord['initial_loss'] = initialLoss
    visitRecord['parent_lock_action'] = timing.get('parentLockAction', '')

    # The current placement is always a candidate, so a visit should never make things worse.
    # If it does (or still collides), keep the current placement and move on in the traversal.
    colliding = timedCollidingConfigurations(tree) if includeCollisionPenalty else []
    increased = loss > initialLoss + 1e-9 * max(1, abs(initialLoss))
    if len(colliding) > 0 and not increased:
        # If the current placement already collided in each of these configurations, reverting can't help:
        # accept the change, which the objective (collision penalty included) says is no worse
        before = timedCollidingConfigurations(subject)
        if all(c in before for c in colliding):
            reason = (f"still colliding in configurations {colliding}, but the placement before this visit collided "
                      f"in {before}; accepted since the objective did not increase")
            logCollisionPenalty(f"WARNING: {label} {index}: {reason}. Loss {initialLoss} -> {loss}.")
            visitRecord.update(revert_reason=reason, final_loss=loss)
            return tree, loss
    if len(colliding) > 0 or increased:
        if len(colliding) > 0 and not increased:
            reason = f"collision in configurations {colliding} did not raise the objective (possible bug in the objective)"
        elif len(colliding) > 0:
            reason = f"still colliding in configurations {colliding} after retrying with penalty"
        else:
            reason = "objective increased"
        logCollisionPenalty(f"WARNING: reverting {label} {index} ({reason}). Loss {initialLoss} -> {loss}; keeping previous placement.")
        visitRecord.update(reverted=True, revert_reason=reason, final_loss=initialLoss, parent_lock_action='')
        return subject.copyAbbreviatedSelf(), initialLoss
    if timing.get('notApplied'):
        visitRecord.update(reverted=True, revert_reason="best particle could not be applied")
    visitRecord['final_loss'] = loss
    return tree, loss

# Columns of visits.csv, one row per visit
VISIT_FIELDS = ['visit', 'repeat', 'index', 'kind', 'initial_loss', 'final_loss',
                'retried', 'reverted', 'revert_reason',
                'nopenalty_pso_s', 'nopenalty_nm_s', 'penalty_pso_s', 'penalty_nm_s', 'collision_check_s',
                'visit_s', 'elapsed_s', 'total_length', 'sum_squared_length',
                'parent_lock_action', 'parent_locked_count']

def appendVisitRecord(directory, visitRecord):
    path = os.path.join(directory, 'visits.csv')
    isNew = not os.path.exists(path)
    with open(path, 'a', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=VISIT_FIELDS, restval='', extrasaction='ignore')
        if isNew:
            writer.writeheader()
        writer.writerow(visitRecord)

def sumSquaredLinkLength(tree):
    # Tree-wide version of the visits' objective at the default power 2 (no collision penalty), for comparing with total length
    return sum(link.path.length ** 2 for link in tree.Links)

def writeRunInfo(directory, settings, filename='run_info.json', note=None):
    # Records which code and settings produced a run, and why it was run (note), so results can be
    # traced back to them later
    if not note:
        print("Warning: no --note given, so this run's run info won't record why it was run.")
    import subprocess, platform
    repo = os.path.dirname(os.path.abspath(__file__))
    def git(*args):
        try:
            return subprocess.run(['git', *args], cwd=repo, capture_output=True, text=True, check=True).stdout.rstrip()
        except Exception as e:
            return f"unavailable ({e})"
    info = {
        'note': note,
        'start_time': datetime.now().isoformat(),
        'git_commit': git('rev-parse', 'HEAD'),
        'git_branch': git('rev-parse', '--abbrev-ref', 'HEAD'),
        'uncommitted_changes_to_tracked_files': git('status', '--porcelain', '--untracked-files=no').splitlines(),
        'python_version': platform.python_version(),
        'numpy_version': np.__version__,
        'pyswarms_version': ps.__version__,
        'settings': settings,
    }
    with open(os.path.join(directory, filename), 'w') as f:
        json.dump(info, f, indent=2, default=str)

def optimizeJointPlacement(subject, index, maxiter, tol, failurePenalty,
                           childFraction = 1, ignoreLater = False, parallelize = False,
                           verbose=True, power=1, includeCollisionPenalty=True,
                           configurations=None, visitRecord=None, parentLocking=None):
    if configurations is None:
        configurations = [[0] * len(subject.Joints)]
    if parentLocking is None:
        parentLocking = not ignoreLater  # default: on whenever the whole tree is considered
    if parentLocking and ignoreLater:
        raise ValueError("parentLocking needs the whole tree (ignoreLater=False, i.e. guarantee=True)")
    warnIfInitialCollisions(subject, configurations, "optimizeJointPlacement")

    optimizeOnce = partial(optimizeJointPlacementOnce, subject, index, maxiter, tol, failurePenalty,
                           childFraction, ignoreLater, parallelize, verbose, power, configurations, parentLocking)
    movedIndices = [index] + parentLockedDescendants(subject, index) if parentLocking else [index]
    return optimizeWithCollisionRetry(subject, index, "joint", optimizeOnce, includeCollisionPenalty, configurations,
                                      visitRecord, movedIndices)

def optimizeJointPlacementOnce(subject, index, maxiter, tol, failurePenalty, childFraction, ignoreLater,
                               parallelize, verbose, power, configurations, parentLocking, includeCollisionPenalty):
    collisionMatrices = subject.buildCollisionMatrices()
    # parent-locked waypoints move with the joint, and the links leaving them count as the joint's outgoing links
    locked = parentLockedDescendants(subject, index) if parentLocking else []
    outgoing = lockedOutgoingLinks(subject, index, locked) if locked else None
    # the collision penalty covers everything that moves, carried waypoints and their links included
    movedJointIndex = [index] + locked if locked else index
    def move(tree, transformation):
        return moveWithParentLocked(tree, index, transformation, locked)

    def objective(params, returnWhich : bool = False) -> float:
        tree = subject.copyAbbreviatedSelf(ignoreLater, index)

        translation = params[0]
        rotation = params[1]

        transform = SE3.Trans([0,0,translation]) @ SE3.Rz(rotation)

        pathNonExistancePenalty = failurePenalty * (len(subject.Joints) ** 4) * (len(subject.Children[index]) + 1)
        linkLossReversedZhat = pathNonExistancePenalty

        #try just moving it
        if move(tree, transform):
            linkLossSameZhat = linkLoss(tree,
                                        index,
                                        power=power, childFraction=childFraction,
                                        collisionMatrices=collisionMatrices,
                                        movedJointIndex=movedJointIndex,
                                        includeCollisionPenalty=includeCollisionPenalty,
                                        configurations=configurations,
                                        collisionErrorWeight=failurePenalty, outgoing=outgoing)
        else:
            linkLossSameZhat = pathNonExistancePenalty

        #try switching zhat
        tree.Joints[index].reverseZhat()
        if (linkLossSameZhat == pathNonExistancePenalty or
             not move(tree, SE3())):
            tree2 = subject.copyAbbreviatedSelf(ignoreLater, index)
            tree2.Joints[index].reverseZhat()
            if move(tree2, SE3.Trans([0,0,-translation]) @ SE3.Rz(-rotation)):
                linkLossReversedZhat = linkLoss(tree2,
                                                index,
                                                power=power, childFraction=childFraction,
                                                collisionMatrices=collisionMatrices,
                                                movedJointIndex=movedJointIndex,
                                                includeCollisionPenalty=includeCollisionPenalty,
                                                configurations=configurations,
                                                collisionErrorWeight=failurePenalty, outgoing=outgoing)
            else:
                linkLossReversedZhat = pathNonExistancePenalty
        else:
            linkLossReversedZhat = linkLoss(tree,
                                            index,
                                            power=power, childFraction=childFraction,
                                            collisionMatrices=collisionMatrices,
                                            movedJointIndex=movedJointIndex,
                                            includeCollisionPenalty=includeCollisionPenalty,
                                            configurations=configurations,
                                            collisionErrorWeight=failurePenalty, outgoing=outgoing)

        if returnWhich:
            if linkLossSameZhat <= linkLossReversedZhat:
                return 1
            else:
                return 2
        else:
            return min(linkLossSameZhat,linkLossReversedZhat)

    start = time.time()

    joint = subject.Joints[index]
    # The root has no parent (Parents[0] == -1), so it gets no line-up guess: only the current placement
    # and random particles
    hasParent = subject.Parents[index] >= 0
    if hasParent:
        parent = subject.Joints[subject.Parents[index]]

        frame1 = joint.Pose
        frame2 = parent.DistalDubinsFrame()

        transformation = frame1.inv() * frame2

        #try to make initial guess right next to each other
        initialPosition = transformation.t[2]
        initialRotation = np.arctan2(transformation.R[1, 0], transformation.R[0, 0])
        initialGuess = [initialPosition,initialRotation]
        initialTree = subject.copyAbbreviatedSelf(ignoreLater, index)

        if not move(initialTree, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation)):
            initialTree = subject.copyAbbreviatedSelf(ignoreLater, index)

        if not move(initialTree, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation)):
            initialGuess = [0,0]
    else:
        initialPosition = 0
        initialGuess = [0,0]

    initialLoss = objective([0,0])

    #calculate dist
    dist = localLinkScale(subject, index)
    dist2 = dist*2
    positionBound = [min(-dist2, initialPosition - dist2), max(dist2, initialPosition + dist2)]
    angleBound = [-np.pi*2, np.pi*2]
    bounds = [positionBound, angleBound]

    #initial swarm
    global joint_batch_objective_function
    def joint_batch_objective_function(X):
        return np.array([objective(x) for x in X])
    n_particles = 16

    # particle 0 is the initial guess and particle 1 is the current placement, followed by random particles;
    # the root has no initial guess, so its particle 0 is the current placement
    init_pos = np.tile(np.array(initialGuess, dtype='float64'), (n_particles,1))
    init_pos[1] = np.array([0,0])
    first = 2 if hasParent else 1
    #add random noise
    noise = np.zeros_like(init_pos[first:])
    noise[:, 0] = np.random.uniform(-dist2,dist2, n_particles - first)
    noise[:, 1] = np.random.uniform(-np.pi*2,np.pi*2, n_particles - first)
    init_pos[first:] += noise
    init_pos[first:, 0] = np.clip(init_pos[first:, 0], positionBound[0], positionBound[1])
    init_pos[first:, 1] = wrapToRange(init_pos[first:, 1], angleBound[0], angleBound[1])

    min_bound = np.array([b[0] for b in bounds])
    max_bound = np.array([b[1] for b in bounds])
    dimensions = 2
    if not min_bound.shape == (dimensions,) or not max_bound.shape == (dimensions,):
        raise ValueError(f"Bounds arrays must be of shape ({dimensions},)")

    optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,dimensions=2,options={'c1':0.6, 'c2':0.7, 'w':0.5},
                                        bounds=(min_bound, max_bound),init_pos=init_pos,ftol=tol)
    psoStart = time.time()
    minSwarmLoss, minSwarmResult = optimizer.optimize(joint_batch_objective_function, iters=int((maxiter + 1)/2),verbose=False,
                                                        n_processes=n_particles if parallelize else None)
    nmStart = time.time()
    nelderMead = minimize(objective, minSwarmResult, method="Nelder-Mead", bounds=bounds, tol=tol, options={
        'maxiter':int(maxiter/2),
        'fatol':tol,
    })
    result = nelderMead.x
    loss = nelderMead.fun
    timing = {'pso_s': nmStart - psoStart, 'nm_s': time.time() - nmStart}

    if verbose:
        msg = f"Optimized joint {index} in {time.time() - start}s -- Old loss: {initialLoss}, Improved Loss: {loss}"
        if includeCollisionPenalty:
            logCollisionPenalty(msg)
        else:
            print(msg)

    which = objective(result, returnWhich=True)

    tree = subject.copyAbbreviatedSelf()
    if which == 1:
        if not move(tree, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1])):
            raise Exception()
        return tree, loss, initialLoss, timing
    else:
        try:
            if not move(tree, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1])):
                raise Exception()

            tree.Joints[index].reverseZhat()
            if not move(tree, SE3.Trans([0,0,0])):
                raise Exception()

            return tree, loss, initialLoss, timing
        except:
            tree2 = subject.copyAbbreviatedSelf()
            tree2.Joints[index].reverseZhat()

            if not move(tree2, SE3.Trans([0,0,-result[0]]) @ SE3.Rz(-result[1])):
                raise Exception()
            return tree2, loss, initialLoss, timing

WAYPOINT_SAMPLINGS = ("box", "multiscale")

def randomRotation():
    # Uniformly distributed rotation matrix (from a uniformly distributed unit quaternion), using np.random
    w, x, y, z = (lambda q: q / np.linalg.norm(q))(np.random.normal(size=4))
    return np.array([[1 - 2*(y*y + z*z), 2*(x*y - z*w),     2*(x*z + y*w)],
                     [2*(x*y + z*w),     1 - 2*(x*x + z*z), 2*(y*z - x*w)],
                     [2*(x*z - y*w),     2*(y*z + x*w),     1 - 2*(x*x + y*y)]])

def rotationBetween(a, b):
    # Smallest rotation matrix taking unit vector a to unit vector b
    v = np.cross(a, b)
    c = np.dot(a, b)
    if np.linalg.norm(v) < 1e-12:
        if c > 0:
            return np.eye(3)
        # opposite directions: rotate pi about any axis perpendicular to a
        axis = np.cross(a, [1, 0, 0] if abs(a[0]) < 0.9 else [0, 1, 0])
        axis /= np.linalg.norm(axis)
        return 2 * np.outer(axis, axis) - np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx / (1 + c)

def straightThroughGuesses(subject, index):
    # Waypoint moves (in its own frame, as [t, ZYZ angles]) that point it along the line from the end of
    # the previous part to the start of the next one, the direction that would straighten both adjacent
    # links: one turning it in place, and one also moving it onto that line.
    waypoint = subject.Joints[index]
    start = subject.Joints[subject.Parents[index]].DistalDubinsFrame().t
    if len(subject.Children[index]) > 0:
        end = subject.Joints[subject.Children[index][0]].ProximalDubinsFrame().t
    else:
        end = waypoint.Pose.t
    direction = end - start
    if np.linalg.norm(direction) < 1e-9:
        return []
    direction = direction / np.linalg.norm(direction)
    orientation = rotationBetween(waypoint.Pose.R[:, 2], direction) @ waypoint.Pose.R
    onLine = start + np.dot(waypoint.Pose.t - start, direction) * direction
    guesses = []
    for position in (waypoint.Pose.t, onLine):
        move = waypoint.Pose.inv() * SE3.Rt(orientation, position, check=False)
        guesses.append(list(move.t) + list(move.eul()))
    return guesses

def optimizeWaypointPlacement(subject, index, maxiter, tol,
                              failurePenalty, childFraction = 1,
                              ignoreLater=False, parallelize=False, verbose = True, configurations=None,
                              includeCollisionPenalty=True, visitRecord=None, power=1,
                              waypointSampling="box", straightThroughGuess=True, alignmentWeight=10,
                              parentLocking=None):
    if configurations is None:
        configurations = [[0] * len(subject.Joints)]
    if waypointSampling not in WAYPOINT_SAMPLINGS:
        raise ValueError(f"waypointSampling must be one of {WAYPOINT_SAMPLINGS}")
    if parentLocking is None:
        parentLocking = not ignoreLater  # default: on whenever the whole tree is considered
    if parentLocking and ignoreLater:
        raise ValueError("parentLocking needs the whole tree (ignoreLater=False, i.e. guarantee=True)")
    warnIfInitialCollisions(subject, configurations, "optimizeWaypointPlacement")

    optimizeOnce = partial(optimizeWaypointPlacementOnce, subject, index, maxiter, tol, failurePenalty,
                           childFraction, ignoreLater, parallelize, verbose, power, configurations,
                           waypointSampling, straightThroughGuess, alignmentWeight, parentLocking)
    movedIndices = [index] + parentLockedDescendants(subject, index) if parentLocking else [index]
    return optimizeWithCollisionRetry(subject, index, "waypoint", optimizeOnce, includeCollisionPenalty, configurations,
                                      visitRecord, movedIndices)

def optimizeWaypointPlacementOnce(subject, index, maxiter, tol, failurePenalty, childFraction, ignoreLater,
                                  parallelize, verbose, power, configurations,
                                  waypointSampling, straightThroughGuess, alignmentWeight, parentLocking,
                                  includeCollisionPenalty):
    start = time.time()

    parent = subject.Joints[subject.Parents[index]]
    waypoint = subject.Joints[index]

    # Params are a move in the waypoint's own frame (relative=True): a translation along its own axes,
    # then a ZYZ rotation about its own origin. (In the world frame, rotations pivot about the world
    # origin, so even rotating in place would need a precisely matched counter-translation.)
    def paramsToTransform(params):
        return SE3.Trans(params[0:3]) @ SE3.Rz(params[3]) @ SE3.Ry(params[4]) @ SE3.Rz(params[5])

    # World-frame transform that would put the waypoint's proximal frame on the parent's distal frame,
    # expressed as a move in the waypoint's own frame for the initial guess
    alignTransform = parent.DistalDubinsFrame() * waypoint.ProximalDubinsFrame().inv()
    alignedOrientation = alignTransform.R @ waypoint.Pose.R
    localAlign = waypoint.Pose.inv() * alignTransform * waypoint.Pose
    initialGuess = list(localAlign.t) + list(localAlign.eul())

    collisionMatrices = subject.buildCollisionMatrices()
    # parent-locked waypoints move with this one, and the links leaving them count as its outgoing links
    locked = parentLockedDescendants(subject, index) if parentLocking else []
    outgoing = lockedOutgoingLinks(subject, index, locked) if locked else None
    # the collision penalty covers everything that moves, carried waypoints and their links included
    movedJointIndex = [index] + locked if locked else index
    selfLocked = parentLocking and isParentLocked(subject, index)

    def waypointLinkLoss(tree):
        # on a copy, since the collision penalty leaves the tree in the last configuration
        return linkLoss(tree.copyAbbreviatedSelf(), index, power=power, childFraction=childFraction,
                        includeCollisionPenalty=includeCollisionPenalty, configurations=configurations,
                        collisionMatrices=collisionMatrices, movedJointIndex=movedJointIndex,
                        collisionErrorWeight=failurePenalty, outgoing=outgoing)

    def objective(params):
        tree = subject.copyAbbreviatedSelf(ignoreLater, index)

        if not moveWithParentLocked(tree, index, paramsToTransform(params), locked):
            return failurePenalty * len(subject.Joints) * (len(subject.Children) + 1)

        # regularize toward aligning with the parent's distal frame
        misalignment = rotationAngle(tree.Joints[index].Pose.R.T @ alignedOrientation)
        return linkLoss(tree, index, power=power, childFraction=childFraction,
                        includeCollisionPenalty=includeCollisionPenalty, configurations=configurations,
                        collisionMatrices=collisionMatrices, movedJointIndex=movedJointIndex,
                        collisionErrorWeight=failurePenalty, outgoing=outgoing) + misalignment * alignmentWeight

    initialTree = subject.copyAbbreviatedSelf(ignoreLater, index)
    if not moveWithParentLocked(initialTree, index, paramsToTransform(initialGuess), locked):
        initialGuess = [0]*6

    initialLoss = objective([0]*6)

    dist = localLinkScale(subject, index)
    bounds = [(min(-dist*2, initialGuess[i] - dist*2), max(dist*2, initialGuess[i] + dist*2)) for i in range(0,3)] \
                + [(-np.pi*2, np.pi*2)] * 3

    global waypoint_batch_objective_function
    def waypoint_batch_objective_function(X):
        return np.array([objective(x) for x in X])

    n_particles = 24

    # particle 0 is the initial guess and particle 1 is the current placement,
    # followed by the straight-through guesses (if enabled), the unlocking guess (for a parent-locked
    # waypoint), then random particles
    init_pos = np.tile(np.array(initialGuess, dtype='float64'), (n_particles,1))
    init_pos[1] = np.array([0]*6, dtype='float64')
    guesses = straightThroughGuesses(subject, index) if straightThroughGuess else []
    if selfLocked:
        guesses = guesses + unlockGuess(subject, index)
    for k, guess in enumerate(guesses):
        init_pos[2 + k] = guess
    first = 2 + len(guesses)
    nRandom = n_particles - first
    if waypointSampling == "box":
        #add random noise
        noise = np.zeros_like(init_pos[first:])
        for i in range(0,3):
            noise[:, i] = np.random.uniform(-dist*2,dist*2, nRandom)
        for i in range(3, 6):
            noise[:, i] = np.random.uniform(-np.pi*2,np.pi*2, nRandom)
        init_pos[first:] += noise
    else: # "multiscale"
        # Thirds: rotate in place, small moves, and moves over the full range, all with rotations
        # sampled uniformly over orientations. Fixing a bunched waypoint usually needs a large
        # rotation with little or no translation, which box sampling almost never produces.
        scales = np.array_split(np.arange(nRandom), 3)
        translations = np.zeros((nRandom, 3))
        translations[scales[1]] = np.random.uniform(-dist/2, dist/2, (len(scales[1]), 3))
        translations[scales[2]] = np.random.uniform([b[0] for b in bounds[:3]], [b[1] for b in bounds[:3]],
                                                    (len(scales[2]), 3))
        init_pos[first:, 0:3] = translations
        init_pos[first:, 3:6] = [SE3.Rt(randomRotation(), np.zeros(3)).eul() for _ in range(nRandom)]
    for i in range(0,3):
        init_pos[2:, i] = np.clip(init_pos[2:, i], bounds[i][0], bounds[i][1])
    for i in range(3, 6):
        init_pos[2:, i] = wrapToRange(init_pos[2:, i], bounds[i][0], bounds[i][1])
    optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,dimensions=6,options={'c1':0.7, 'c2':0.5, 'w':0.5},bounds=(np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])),init_pos=init_pos,ftol=tol)
    psoStart = time.time()
    minSwarmLoss, minSwarmResult = optimizer.optimize(waypoint_batch_objective_function, iters=maxiter,verbose=False, n_processes=n_particles if parallelize else None)

    timing = {'pso_s': time.time() - psoStart, 'nm_s': 0.0}

    tree = subject.copyAbbreviatedSelf()
    if not moveWithParentLocked(tree, index, paramsToTransform(minSwarmResult), locked):
        logCollisionPenalty(f"WARNING: best particle for waypoint {index} could not be applied; keeping previous placement.")
        return subject.copyAbbreviatedSelf(), initialLoss, initialLoss, dict(timing, notApplied=True)

    if parentLocking:
        # Compare locked against unlocked on linkLoss alone (the alignment penalty is always 0 when locked, so
        # including it would favor locking for reasons unrelated to length), keeping the current state unless
        # the other is better by PARENT_LOCK_MARGIN.
        psoLinkLoss = waypointLinkLoss(tree)
        if selfLocked:
            currentLinkLoss = waypointLinkLoss(subject)
            if not isParentLocked(tree, index) and psoLinkLoss < currentLinkLoss * (1 - PARENT_LOCK_MARGIN):
                timing['parentLockAction'] = 'unlocked'
                return tree, minSwarmLoss, initialLoss, timing
            # The swarm stayed locked (or barely improved). Its objective includes the alignment penalty,
            # which is 0 when locked, so it rarely picks the unlocking guess even when that's much shorter;
            # judge that guess on linkLoss directly, as locking is judged.
            for guess in unlockGuess(subject, index):
                unlockedTree = subject.copyAbbreviatedSelf()
                if moveWithParentLocked(unlockedTree, index, paramsToTransform(guess), locked):
                    unlockedLinkLoss = waypointLinkLoss(unlockedTree)
                    if unlockedLinkLoss < currentLinkLoss * (1 - PARENT_LOCK_MARGIN):
                        timing['parentLockAction'] = 'unlocked'
                        if verbose:
                            print(f"Unlocked waypoint {index}: linkLoss {currentLinkLoss} locked vs. {unlockedLinkLoss} unlocked")
                        # Reported on linkLoss, the basis of the decision (the locked objective equals its linkLoss)
                        return unlockedTree, unlockedLinkLoss, initialLoss, timing
            timing['parentLockAction'] = 'keptLocked'
            return subject.copyAbbreviatedSelf(), initialLoss, initialLoss, timing
        elif isParentLocked(tree, index):
            # the swarm itself landed on the locked placement (particle 0, the line-up guess, is exactly it)
            timing['parentLockAction'] = 'parentLocked'
        else:
            lockedTree = subject.copyAbbreviatedSelf()
            lockMove = parent.DistalDubinsFrame() * waypoint.ProximalDubinsFrame().inv()
            if moveWithParentLocked(lockedTree, index, lockMove, locked, relative=False):
                lockedLinkLoss = waypointLinkLoss(lockedTree)
                if lockedLinkLoss < psoLinkLoss * (1 - PARENT_LOCK_MARGIN):
                    timing['parentLockAction'] = 'parentLocked'
                    misalignment = rotationAngle(lockedTree.Joints[index].Pose.R.T @ alignedOrientation)
                    lockedLoss = lockedLinkLoss + misalignment * alignmentWeight
                    if verbose:
                        print(f"Parent-locked waypoint {index}: linkLoss {psoLinkLoss} unlocked vs. {lockedLinkLoss} locked")
                    return lockedTree, lockedLoss, initialLoss, timing

    if verbose:
        msg = f"Optimized waypoint {index} in {time.time() - start}s -- Old Loss: {initialLoss}, Improved Loss: {minSwarmLoss}"
        if includeCollisionPenalty:
            logCollisionPenalty(msg)
        else:
            print(msg)
    return tree, minSwarmLoss, initialLoss, timing

def optimizeTree(subject, showSteps=False, childFraction=1, guarantee=False, parallelize=False, 
                 evaluate=False, verbose=True, directory=None, resetOnFail=False,
                 traversal="dfs", direction="outward", orderBy="longest", power=1, configurations=None, 
                 repeatTraversal=1, save_checkpoints=True, traversalPower=2,
                 waypointSampling="box", straightThroughGuess=True, alignmentWeight=10, parentLocking=None):
    """
    Optimize the placement of joints in a kinematic tree.
    
    Args:
        subject: KinematicTree to optimize
        showSteps: Show visualization after optimization
        childFraction: Weight on outgoing links in each visit's objective (and the randomized traversal's weights)
        guarantee: If True, don't ignore later joints
        parallelize: Enable parallel computation
        evaluate: If True, return (tree, times, lengths) tuple
        verbose: Print progress information
        directory: Output directory for saving results and checkpoints
        resetOnFail: Reset on optimization failure
        traversal: Traversal algorithm ('dfs', 'bfs', 'randomized')
        direction: Traversal direction ('outward', 'inward')
        orderBy: Traversal ordering ('default', 'longest', 'shortest')
        power: Power on link lengths in each visit's objective (1 = plain length, 2 = squared lengths)
        configurations: List of joint configurations to test
        repeatTraversal: Number of times to repeat the traversal
        save_checkpoints: If True, save dill checkpoints after each iteration
        traversalPower: Power on link lengths in the randomized traversal's weights
        waypointSampling: How waypoint visits sample their random initial particles: "box" (uniform over the
            bounds, the original behavior) or "multiscale" (thirds that rotate in place, move a little, or move
            over the full range, with uniformly distributed rotations)
        straightThroughGuess: If True, waypoint visits also start particles that point the waypoint along the
            line from the previous part to the next one, in place and moved onto that line (default True;
            the parameter sweep found it about 21% shorter on 4-joint chains)
        alignmentWeight: Weight (per radian) of the waypoint objective's penalty for misalignment with the
            parent's distal frame; 0 turns it off
        parentLocking: If True, waypoints can be parent-locked onto their parent's distal frame (as if deleted,
            but kept available): visits carry parent-locked descendants along, waypoint visits consider
            locking, and visits to parent-locked waypoints consider unlocking them. Needs guarantee=True.
            Default (None): on when guarantee=True, off otherwise.
    
    Returns:
        Optimized KinematicTree, or (tree, times, lengths) if evaluate=True
    """
    if parentLocking is None:
        parentLocking = guarantee
    if repeatTraversal == "n":
        repeatTraversal = len(subject.Joints)

    times = []
    lengths = []
    losses = []  # Track optimization losses
    completed_indices = set()  # Track completed joint indices
    iteration_count = 0
    
    # Setup checkpoint directory
    checkpoint_dir = None
    if save_checkpoints and directory is not None:
        checkpoint_dir = os.path.join(directory, 'checkpoints')
        os.makedirs(checkpoint_dir, exist_ok=True)
        if verbose:
            print(f"Checkpoints will be saved to: {checkpoint_dir}")

    if configurations == None:
        num_joints = len(subject.Joints)
        random_config = [0] * num_joints
        
        # Set non-waypoint joints to random states within their range limits
        # Add a tolerance margin to avoid getting too close to limits
        tolerance = 0.02  # 2% margin from limits
        for i in range(num_joints):
            if not isWaypoint(subject.Joints[i]):
                joint = subject.Joints[i]
                min_state, max_state = joint.stateRange()
                range_size = max_state - min_state
                margin = range_size * tolerance
                random_config[i] = np.random.uniform(min_state + margin, max_state - margin)
        
        neutral_config = [0] * num_joints  
        configurations = [neutral_config]

        print("Using random configuration for optimization:")
        print("Joint states and their limits:")
        for i in range(num_joints):
            if not isWaypoint(subject.Joints[i]):
                min_state, max_state = subject.Joints[i].stateRange()
                print(f"  Joint {i}: config = {random_config[i]:.3f} (limits: [{min_state:.3f}, {max_state:.3f}])")

    for i in range(0, len(subject.Joints)):
        subject.Joints[i].recomputeCollisionCapsules()

    if subject.detectCollisions(debug=True) > 0:
        print("Warning: initial tree in optimizeTree contains collisions.")
    if showSteps and isinstance(subject.Joints[0], OrigamiJoint):
        subject.show()

    failurePenalty = 0
    for i in range(0, len(subject.Joints)):
        if len(subject.Children[i]) > 0:
            continue
        j = i
        length = 0
        while j != 0:
            length += subject.Links[j].path.length ** 2
            j = subject.Parents[j]
        if length > failurePenalty:
            failurePenalty = length

    print(f"Collision penalty scale is {failurePenalty}")

    start = time.time()

    def log(t, idx):
        diff = time.time() - start
        times.append(diff)
        lengths.append(t.totalLength())
        if directory != None:
            save_path = os.path.join(directory, f"{diff}_{idx}")
            t.save(save_path, saveDir=False)

    tree = subject.copyAbbreviatedSelf()
    log(tree, -1)
    if directory != None:
        appendVisitRecord(directory, dict(visit=0, index=-1, kind='start', elapsed_s=times[-1],
                                          total_length=lengths[-1], sum_squared_length=sumSquaredLinkLength(tree),
                                          parent_locked_count=countParentLocked(tree)))

    # create the traversal calls
    treeTraversals = {
        "dfs": partial(dfs, direction=direction, orderBy=orderBy),
        "bfs": partial(bfs, direction=direction, orderBy=orderBy),
        "randomized": partial(randomized, 
                            power=traversalPower, 
                            count=len(subject.Joints), 
                            childFraction=childFraction,
                            isWeighted=True)
    }

    print("Doing the optimization:")
    for repeat_idx in range(repeatTraversal):
        for index in treeTraversals[traversal](subject):
            print("Optimizing joint index:", index, subject.Joints[index])
            iters = 50
            tolerance = subject.r/10

            visitRecord = {}
            visitStart = time.time()
            if isWaypoint(subject.Joints[index]):
                tree, loss = optimizeWaypointPlacement(tree,index, maxiter=iters, tol=tolerance, 
                                                    failurePenalty=failurePenalty, childFraction=childFraction, power=power, 
                                                    ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose, 
                                                    configurations=configurations, visitRecord=visitRecord,
                                                    waypointSampling=waypointSampling, straightThroughGuess=straightThroughGuess,
                                                    alignmentWeight=alignmentWeight, parentLocking=parentLocking)
            else:
                tree, loss = optimizeJointPlacement(tree,index, maxiter=iters, tol=tolerance, 
                                                    failurePenalty=failurePenalty, childFraction=childFraction, power=power, 
                                                    ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose, 
                                                    configurations=configurations, visitRecord=visitRecord,
                                                    parentLocking=parentLocking)
            visitTime = time.time() - visitStart

            log(tree, index)
            losses.append(loss)
            completed_indices.add(index)
            iteration_count += 1
            if directory != None:
                visitRecord.update(visit=iteration_count, repeat=repeat_idx, index=index, visit_s=visitTime,
                                   elapsed_s=times[-1], total_length=lengths[-1],
                                   sum_squared_length=sumSquaredLinkLength(tree),
                                   parent_locked_count=countParentLocked(tree))
                appendVisitRecord(directory, visitRecord)
            
            # Save checkpoint after each iteration
            if save_checkpoints and checkpoint_dir is not None:
                checkpoint_state = CheckpointState(
                    tree=copy.deepcopy(tree),
                    random_state=np.random.get_state(),
                    traversal_info={
                        'algorithm': traversal,
                        'direction': direction,
                        'orderBy': orderBy,
                        'current_index': index,
                        'repeat_index': repeat_idx,
                        'total_repeats': repeatTraversal,
                        'completed_indices': list(completed_indices),
                    },
                    collision_matrices=tree.buildCollisionMatrices(),
                    optimization_history={
                        'times': times.copy(), 
                        'lengths': lengths.copy(), 
                        'losses': losses.copy()
                    },
                    configurations=configurations,
                    failure_penalty=failurePenalty,
                    metadata={
                        'timestamp': datetime.now().isoformat(),
                        'iteration_count': iteration_count,
                        'total_joints': len(tree.Joints),
                        'directory': directory,
                        'guarantee': guarantee,
                        'parallelize': parallelize,
                        'childFraction': childFraction,
                        'power': power,
                        'traversalPower': traversalPower,
                        'waypointSampling': waypointSampling,
                        'straightThroughGuess': straightThroughGuess,
                        'alignmentWeight': alignmentWeight,
                        'parentLocking': parentLocking,
                        'evaluate': evaluate,
                        'save_checkpoints': save_checkpoints,
                        'version': '1.0',
                    }
                )
                checkpoint_path = saveCheckpoint(
                    checkpoint_state, 
                    os.path.join(checkpoint_dir, checkpointName(iteration_count, repeat_idx, index))
                )
                if verbose:
                    print(f"Checkpoint saved: {checkpoint_path}")
        
        # Clear completed indices for next repeat
        completed_indices.clear()

    if showSteps:
        tree.show()

    print(f"TOTAL OPTIMIZATION TIME: {time.time() - start}")

    if directory != None:
        tree.save(os.path.join(directory, "final"), saveDir=False)
    if (evaluate):
        return tree, times, lengths
    
    return tree