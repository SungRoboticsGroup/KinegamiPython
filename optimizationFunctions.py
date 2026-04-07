from KinematicTree import *
from treeTraversals import *
import pyswarms as ps
import os
import dill
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Any, Tuple
from datetime import datetime

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

def listCheckpoints(checkpoint_dir: str) -> list:
    """
    List all checkpoint files in a directory, sorted by iteration order.
    
    Args:
        checkpoint_dir: Path to directory containing checkpoint files
        
    Returns:
        List of dicts with 'path', 'repeat', 'joint_idx' for each checkpoint
    """
    import re
    
    checkpoints = []
    
    if not os.path.exists(checkpoint_dir):
        return checkpoints
    
    for filename in os.listdir(checkpoint_dir):
        if filename.endswith('.dill'):
            match = re.match(r'iter_(\d+)_(\d+)\.dill', filename)
            if match:
                checkpoints.append({
                    'path': os.path.join(checkpoint_dir, filename),
                    'repeat': int(match.group(1)),
                    'joint_idx': int(match.group(2)),
                    'filename': filename
                })
    
    # Sort by repeat, then by joint index
    checkpoints.sort(key=lambda x: (x['repeat'], x['joint_idx']))

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
    power = metadata.get('power', 2)
    evaluate = metadata.get('evaluate', False)
    save_checkpoints = metadata.get('save_checkpoints', True)
    
    if directory is None:
        directory = metadata.get('directory')
    
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
                            power=power, 
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
    iteration_count = len(completed_indices)
    
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
            
            if isWaypoint(tree.Joints[index]):
                tree, loss = optimizeWaypointPlacement(
                    tree, index, maxiter=iters, tol=tolerance,
                    failurePenalty=failurePenalty, childFraction=childFraction,
                    ignoreLater=(not guarantee), parallelize=parallelize,
                    verbose=verbose, configurations=configurations)
            else:
                tree, loss = optimizeJointPlacement(
                    tree, index, maxiter=iters, tol=tolerance,
                    failurePenalty=failurePenalty, childFraction=childFraction,
                    ignoreLater=(not guarantee), parallelize=parallelize,
                    verbose=verbose, configurations=configurations)
            
            if tree.detectCollisions(specificJointIndex=index, debug=True) > 0:
                print(f"Post-optimization collision detected at joint {index}.")
                raise Exception("Post-optimization collision detected.")
            
            log(tree, index)
            losses.append(loss)
            completed_indices.add(index)
            iteration_count += 1
            
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
                        'evaluate': evaluate,
                        'save_checkpoints': save_checkpoints,
                        'version': '1.0',
                    }
                )
                saveCheckpoint(checkpoint_state, os.path.join(checkpoint_dir, f'iter_{rep}_{index}'))
        
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

def optimizeJointPlacement(subject, index, maxiter, tol, failurePenalty, 
                           childFraction = 1, ignoreLater = False, parallelize = False, 
                           verbose=True, power=2, includeCollisionPenalty=True, 
                           configurations=None, retryingWithPenalty=False):

    # If includeCollisionPenalty is True, first try without it
    original_includeCollisionPenalty = includeCollisionPenalty
    verbose = True
    if original_includeCollisionPenalty and not retryingWithPenalty:
        # Temporarily disable collision penalty for first attempt
        includeCollisionPenalty = False

    if configurations is None:
        configurations = [[0] * len(subject.Joints)]

    for configuration in configurations:
        copied_subject = copy.deepcopy(subject)
        for i in range(0, len(copied_subject.Joints)):
            if (not isWaypoint(copied_subject.Joints[i])):
                copied_subject.setJointState(i,configuration[i])
                copied_subject.Joints[i].recomputeCollisionCapsules()
        if copied_subject.detectCollisions(debug=True) > 0:
            print(f"Warning: Initial tree in optimizeJointPlacement contains collisions in configuration {configuration}.")

    subjects = [copy.deepcopy(subject) for _ in configurations]
    for i in range(0,len(configurations)):
        subjects[i].setConfiguration(configurations[i])

    movedJointIndex = index
    collisionMatrices = subject.buildCollisionMatrices()
    
    def objective(params, returnWhich : bool = False) -> float:
        tree = subject.copyAbbreviatedSelf(ignoreLater, index)

        translation = params[0]
        rotation = params[1]

        transform = SE3.Trans([0,0,translation]) @ SE3.Rz(rotation)

        pathNonExistancePenalty = failurePenalty * (len(subject.Joints) ** 4) * (len(subject.Children[index]) + 1)
        linkLossReversedZhat = pathNonExistancePenalty

        #try just moving it
        if tree.transformJoint(index, transform, propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
            # print(f"Before linkLossSameZhat, includeCollisionPenalty = {includeCollisionPenalty} and collisionErrorWeight = {failurePenalty}")
            linkLossSameZhat = linkLoss(tree, 
                                        index, 
                                        power=power,
                                        collisionMatrices=collisionMatrices,
                                        movedJointIndex=movedJointIndex,
                                        includeCollisionPenalty=includeCollisionPenalty, 
                                        configurations=configurations,
                                        collisionErrorWeight=failurePenalty)
        else:
            linkLossSameZhat = pathNonExistancePenalty
        
        #try switching zhat
        tree.Joints[index].reverseZhat()
        if (linkLossSameZhat == pathNonExistancePenalty or
             not tree.transformJoint(index, SE3(), safe=True, relative=True, propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False)):
            tree2 = subject.copyAbbreviatedSelf(ignoreLater, index)
            tree2.Joints[index].reverseZhat()
            if tree2.transformJoint(index, SE3.Trans([0,0,-translation]) @ SE3.Rz(-rotation), safe=True, relative=True, 
                                    propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False):
                # print(f"Before linkLossReversedZhat (1), includeCollisionPenalty = {includeCollisionPenalty} and collisionErrorWeight = {failurePenalty}")
                linkLossReversedZhat = linkLoss(tree2, 
                                                index, 
                                                power=power,
                                                collisionMatrices=collisionMatrices,
                                                movedJointIndex=movedJointIndex,
                                                includeCollisionPenalty=includeCollisionPenalty,
                                                configurations=configurations,
                                                collisionErrorWeight=failurePenalty)
            else:
                linkLossReversedZhat = pathNonExistancePenalty
        else:
            # print(f"Before linkLossReversedZhat (2), includeCollisionPenalty = {includeCollisionPenalty} and collisionErrorWeight = {failurePenalty}")
            linkLossReversedZhat = linkLoss(tree, 
                                            index, 
                                            power=power,
                                            collisionMatrices=collisionMatrices,
                                            movedJointIndex=movedJointIndex,
                                            includeCollisionPenalty=includeCollisionPenalty, 
                                            configurations=configurations,
                                            collisionErrorWeight=failurePenalty)

        if returnWhich:
            if linkLossSameZhat <= linkLossReversedZhat:
                return 1
            else:
                return 2
        else:
            return min(linkLossSameZhat,linkLossReversedZhat)

    start = time.time()

    joint = subject.Joints[index]
    parent = subject.Joints[subject.Parents[index]]

    frame1 = joint.Pose
    frame2 = parent.DistalDubinsFrame()

    transformation = frame1.inv() * frame2

    #try to make initial guess right next to each other
    initialPosition = transformation.t[2]
    initialRotation = np.arctan2(transformation.R[1, 0], transformation.R[0, 0])
    initialGuess = [initialPosition,initialRotation]
    initialTree = subject.copyAbbreviatedSelf(ignoreLater, index)
    
    if not initialTree.transformJoint(index, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation), 
                                      propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
        initialTree = subject.copyAbbreviatedSelf(ignoreLater, index)
    
    if not initialTree.transformJoint(index, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation), 
                                      propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
        initialGuess = [0,0]

    initialLoss = objective([0,0])

    #calculate dist
    dist = subject.Links[index].path.length
    dist2 = dist*2
    #bounds = [(-dist*2, dist*2), (-np.pi*2, np.pi*2)]
    positionBound = [min(-dist2, initialPosition - dist2), max(dist2, initialPosition + dist2)]
    angleBound = [-np.pi*2, np.pi*2]
    bounds = [positionBound, angleBound]

    #initial swarm
    global joint_batch_objective_function
    def joint_batch_objective_function(X):
        return np.array([objective(x) for x in X])
    n_particles = 16

    init_pos = np.tile(np.array(initialGuess, dtype='float64'), (n_particles,1))
    init_pos[1] = np.array([0,0])
    #add random noise
    noise = np.zeros_like(init_pos[2:])
    noise[:, 0] = np.random.uniform(-dist2,dist2, n_particles - 2)
    noise[:, 1] = np.random.uniform(-np.pi*2,np.pi*2, n_particles - 2)
    init_pos[2:] += noise
    init_pos[2:, 0] = np.clip(init_pos[2:, 0], -dist2, dist2)
    init_pos[2:, 1] = np.clip(init_pos[2:, 1], -np.pi*2, np.pi*2)

    min_bound = np.array([b[0] for b in bounds])
    max_bound = np.array([b[1] for b in bounds])
    dimensions = 2
    if not min_bound.shape == (dimensions,) or not max_bound.shape == (dimensions,):
        raise ValueError(f"Bounds arrays must be of shape ({dimensions},)")

    optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,dimensions=2,options={'c1':0.6, 'c2':0.7, 'w':0.5},
                                        bounds=(min_bound, max_bound),init_pos=init_pos,ftol=tol)
    minSwarmLoss, minSwarmResult = optimizer.optimize(joint_batch_objective_function, iters=int((maxiter + 1)/2),verbose=False, 
                                                        n_processes=n_particles if parallelize else None)

    nelderMead = minimize(objective, minSwarmResult, method="Nelder-Mead", bounds=bounds, tol=tol, options={
        'maxiter':int(maxiter/2),
        'fatol':tol,
    })
    result = nelderMead.x
    loss = nelderMead.fun

    #print(minSwarmLoss, loss)

    if verbose:
        msg = f"Optimized joint {index} in {time.time() - start}s -- Old loss: {initialLoss}, Improved Loss: {loss}"
        if retryingWithPenalty:
            logCollisionPenalty(msg)
        else:
            print(msg)

    which = objective(result, returnWhich=True)

    tree = subject.copyAbbreviatedSelf()
    if which == 1:
        if not tree.transformJoint(index, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1]), 
                                   propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
            raise Exception()
        final_tree, final_loss = tree, loss
    else:
        try:
            if not tree.transformJoint(index, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1]), 
                                       propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
                raise Exception()
            
            tree.Joints[index].reverseZhat()
            if not tree.transformJoint(index, SE3.Trans([0,0,0]), safe=True, relative=True, 
                                       propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False):
                raise Exception()

            final_tree, final_loss = tree, loss
        except:
            tree2 = subject.copyAbbreviatedSelf()
            tree2.Joints[index].reverseZhat()

            if not tree2.transformJoint(index, SE3.Trans([0,0,-result[0]]) @ SE3.Rz(-result[1]), 
                                        safe=True, relative=True, propogate=False, 
                                        recomputeLinkPath=True, recomputeBoundingBall=False):
                raise Exception()
            final_tree, final_loss = tree2, loss

    # Check for collisions if we tried without penalty first
    if original_includeCollisionPenalty and not includeCollisionPenalty:
        # Check for collisions in all configurations
        has_collisions = False
        if configurations is not None:
            for configuration in configurations:
                test_tree = copy.deepcopy(final_tree)
                for i in range(len(test_tree.Joints)):
                    if not isWaypoint(test_tree.Joints[i]):
                        test_tree.setJointState(i, configuration[i])
                        test_tree.Joints[i].recomputeCollisionCapsules()  # Fixed: recompute for each joint
                if test_tree.detectCollisions(specificJointIndex=index, debug=False) > 0:
                    has_collisions = True
                    break
        
        if has_collisions:
            if verbose:
                logCollisionPenalty(f"Collision detected for joint {index}, retrying with collision penalty")
            # Retry with collision penalty enabled
            includeCollisionPenalty = True  # Ensure penalty is enabled for retry
            retryingWithPenalty = True
            final_tree, final_loss = optimizeJointPlacement(subject, index, maxiter, tol, failurePenalty,
                                        childFraction, ignoreLater,
                                        parallelize, verbose, power,
                                        includeCollisionPenalty=includeCollisionPenalty, 
                                        configurations=configurations,
                                        retryingWithPenalty=retryingWithPenalty)
        
            if final_tree.detectCollisions(specificJointIndex=index, debug=True) > 0:
                logCollisionPenalty(f"Collision still detected for joint {index} after retrying with penalty. Final tree: {final_tree}")
                logCollisionPenalty(f"Optimization loss for this joint was {final_loss}.")
                raise Exception("Collision detected after retrying with penalty.")
            
        elif verbose:
            print("No collisions detected, using result without penalty")
    
    return final_tree, final_loss

def optimizeWaypointPlacement(subject, index, maxiter, tol, 
                              failurePenalty, childFraction = 1, 
                              ignoreLater=False, parallelize=False, verbose = True, configurations=None, 
                              includeCollisionPenalty=True, retryingWithPenalty=False):

    # If includeCollisionPenalty is True, first try without it
    original_includeCollisionPenalty = includeCollisionPenalty
    if original_includeCollisionPenalty and not retryingWithPenalty:
        try:
            if verbose:
                print("Trying waypoint optimization without collision penalty...")
            
            # Temporarily disable collision penalty for first attempt
            includeCollisionPenalty = False
            
            # Continue with the optimization logic below (will use includeCollisionPenalty=False)
            # We'll check for collisions after optimization and potentially retry
            
        except Exception as e:
            if verbose:
                print(f"Setup for no-penalty waypoint optimization failed: {str(e)}")
            includeCollisionPenalty = True
    start = time.time()
    initialTree = subject.copyAbbreviatedSelf(ignoreLater, index)

    initialGuess = [0]*6

    parent = initialTree.Joints[initialTree.Parents[index]]
    waypoint = initialTree.Joints[index]

    transform = parent.DistalDubinsFrame() * waypoint.ProximalDubinsFrame().inv()
    initialGuess[0:3] = transform.t
    initialGuess[3:6] = SE3.Rt(transform.R, np.zeros(3)).eul()

    for configuration in configurations:
        copied_subject = copy.deepcopy(subject)
        for i in range(0, len(copied_subject.Joints)):
            if (not isWaypoint(copied_subject.Joints[i])):
                copied_subject.setJointState(i,configuration[i])
                copied_subject.Joints[i].recomputeCollisionCapsules()
        if copied_subject.detectCollisions(debug=True) > 0:
            print(f"Warning: Initial tree in optimizeWaypointPlacement contains collisions in configuration {configuration}.")

    subjects = [copy.deepcopy(subject) for _ in configurations]
    for i in range(0,len(configurations)):
        subjects[i].setConfiguration(configurations[i])

    movedJointIndex = index
    collisionMatrices = subject.buildCollisionMatrices()

    def objective(params):
        tree = subject.copyAbbreviatedSelf(ignoreLater, index)

        if not tree.transformJoint(index, SE3.Trans(params[0:3]) @ SE3.Rz(params[3]) @ SE3.Ry(params[4]) @ SE3.Rz(params[5]),  
                                   propogate=False, safe=True, relative=False, recomputeBoundingBall=False):
            return failurePenalty * len(subject.Joints) * (len(subject.Children) + 1)
        
        # print(f"Waypoint optimization, includeCollisionPenalty = {includeCollisionPenalty} and collisionErrorWeight = {failurePenalty}")
        return linkLoss(tree, index, includeCollisionPenalty=includeCollisionPenalty, configurations=configurations, 
                        collisionMatrices=collisionMatrices, movedJointIndex=movedJointIndex,
                        collisionErrorWeight=failurePenalty) + \
            np.linalg.norm(np.array(params[3:6]) - SE3.Rt(transform.R, np.zeros(3)).eul()) * 10

    if not initialTree.transformJoint(index, SE3.Trans(initialGuess[0:3]) @ SE3.Rz(initialGuess[3]) @ SE3.Ry(initialGuess[4]) @ SE3.Rz(initialGuess[5]),  
                                      propogate=False, safe=True, relative=False, recomputeBoundingBall=False):
        initialGuess = [0]*6

    #print(f"INITAL WAYPOINT GUESS LOSS: {objective(initialGuess)}")

    #initialTree.detectCollisions(debug=True)        
        
    initialLoss = objective([0]*6)

    dist = subject.Links[index].path.length + max(np.amax(np.abs(initialGuess)), np.amax(np.abs(subject.Joints[index].Pose.t)))
    bounds = [(-dist*2, dist*2)]*3 + [(-np.pi*2, np.pi*2)] * 3

    global waypoint_batch_objective_function
    def waypoint_batch_objective_function(X):
        return np.array([objective(x) for x in X])

    n_particles = 24

    init_pos = np.tile(np.array(initialGuess, dtype='float64'), (n_particles,1))
    init_pos[1] = np.array([0]*6, dtype='float64')
    #add random noise
    noise = np.zeros_like(init_pos[2:])
    for i in range(0,3):
        noise[:, i] = np.random.uniform(-dist*2,dist*2, n_particles - 2)
    for i in range(3, 6):
        noise[:, i] = np.random.uniform(-np.pi*2,np.pi*2, n_particles - 2)
    init_pos[2:] += noise
    for i in range(0,3):
        init_pos[2:, i] = np.clip(init_pos[2:, i], -dist*2, dist*2)
    for i in range(3, 6):
        init_pos[2:, i] = np.clip(init_pos[2:, i], -np.pi*2, np.pi*2)
    optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,dimensions=6,options={'c1':0.7, 'c2':0.5, 'w':0.5},bounds=(np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])),init_pos=init_pos,ftol=tol)
    minSwarmLoss, minSwarmResult = optimizer.optimize(waypoint_batch_objective_function, iters=maxiter,verbose=False, n_processes=n_particles if parallelize else None)
    
    tree = subject.copyAbbreviatedSelf()
    if tree.transformJoint(index, SE3.Trans(minSwarmResult[0:3]) @ SE3.Rz(minSwarmResult[3]) @ SE3.Ry(minSwarmResult[4]) @ SE3.Rz(minSwarmResult[5]),  
                           propogate=False, safe=True, relative=False, recomputeBoundingBall=False):
        if verbose:
            msg = f"Optimized waypoint {index} in {time.time() - start}s -- Old Loss: {initialLoss}, Improved Loss: {minSwarmLoss}"
            if retryingWithPenalty:
                logCollisionPenalty(msg)
            else:
                print(msg)
        
        # Check for collisions if we tried without penalty first
        if original_includeCollisionPenalty and not includeCollisionPenalty:
            # Check for collisions in all configurations
            has_collisions = False
            if configurations is not None:
                for configuration in configurations:
                    test_tree = copy.deepcopy(tree)
                    for i in range(len(test_tree.Joints)):
                        if not isWaypoint(test_tree.Joints[i]):
                            test_tree.setJointState(i, configuration[i])
                            test_tree.Joints[i].recomputeCollisionCapsules() 
                    if test_tree.detectCollisions(specificJointIndex=index, debug=False) > 0:
                        has_collisions = True
                        break
            
            if has_collisions:
                if verbose:
                    logCollisionPenalty(f"Collision detected for waypoint {index}, retrying with collision penalty")
                # Retry with collision penalty enabled
                tree, minSwarmResult = optimizeWaypointPlacement(subject, index, maxiter, tol, failurePenalty,
                                               childFraction, ignoreLater,
                                               parallelize, verbose, configurations,
                                               includeCollisionPenalty=True, retryingWithPenalty=True)
                
                if tree.detectCollisions(specificJointIndex=index, debug=True) > 0:
                    logCollisionPenalty(f"Collision still detected for waypoint {index} after retrying with penalty. Final tree: {repr(tree)}")
                    logCollisionPenalty(f"Optimization loss for this waypoint was {minSwarmResult}.")
                    raise Exception("Collision detected after retrying with penalty.")  
                
            elif verbose:
                print("No collisions detected, using waypoint result without penalty")
        
        return tree, minSwarmResult
    else:
        raise Exception("Optimization failed dramatically")

def optimizeTree(subject, showSteps=False, childFraction=1, guarantee=False, parallelize=False, 
                 evaluate=False, verbose=True, directory=None, resetOnFail=False,
                 traversal="dfs", direction="outward", orderBy="longest", power=2, configurations=None, 
                 repeatTraversal=1, save_checkpoints=True):
    """
    Optimize the placement of joints in a kinematic tree.
    
    Args:
        subject: KinematicTree to optimize
        showSteps: Show visualization after optimization
        childFraction: Weight for child link costs
        guarantee: If True, don't ignore later joints
        parallelize: Enable parallel computation
        evaluate: If True, return (tree, times, lengths) tuple
        verbose: Print progress information
        directory: Output directory for saving results and checkpoints
        resetOnFail: Reset on optimization failure
        traversal: Traversal algorithm ('dfs', 'bfs', 'randomized')
        direction: Traversal direction ('outward', 'inward')
        orderBy: Traversal ordering ('default', 'longest', 'shortest')
        power: Power for length cost calculation
        configurations: List of joint configurations to test
        repeatTraversal: Number of times to repeat the traversal
        save_checkpoints: If True, save dill checkpoints after each iteration
    
    Returns:
        Optimized KinematicTree, or (tree, times, lengths) if evaluate=True
    """
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

    # create the traversal calls
    treeTraversals = {
        "dfs": partial(dfs, direction=direction, orderBy=orderBy),
        "bfs": partial(bfs, direction=direction, orderBy=orderBy),
        "randomized": partial(randomized, 
                            power=power, 
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

            if isWaypoint(subject.Joints[index]):
                tree, loss = optimizeWaypointPlacement(tree,index, maxiter=iters, tol=tolerance, 
                                                    failurePenalty=failurePenalty, childFraction=childFraction, 
                                                    ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose, 
                                                    configurations=configurations)
            else:
                tree, loss = optimizeJointPlacement(tree,index, maxiter=iters, tol=tolerance, 
                                                    failurePenalty=failurePenalty, childFraction=childFraction, 
                                                    ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose, 
                                                    configurations=configurations)
                
            if tree.detectCollisions(specificJointIndex=index, debug=True) > 0:
                print(f"Post-optimization collision detected at joint {index}. Tree: {repr(tree)}")
                print(f"Optimization loss for this joint was {loss}.")
                raise Exception("Post-optimization collision detected.")
            log(tree, index)
            losses.append(loss)
            completed_indices.add(index)
            iteration_count += 1
            
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
                        'evaluate': evaluate,
                        'save_checkpoints': save_checkpoints,
                        'version': '1.0',
                    }
                )
                checkpoint_path = saveCheckpoint(
                    checkpoint_state, 
                    os.path.join(checkpoint_dir, f'iter_{repeat_idx}_{index}')
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