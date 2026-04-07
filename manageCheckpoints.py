from optimizationFunctions import (
    loadCheckpoint,
    inspectCheckpoint,
    resumeFromCheckpoint,
    listCheckpoints
)
import matplotlib.pyplot as plt

CHECKPOINT_DIR = '/home/samhitha/code/Trials Before Experiments/2026.04.06_21.12.11_Joints2_Chains1_Seed44/0/BFS Outward Shortest n repetitions/checkpoints'
CHECKPOINT_FILE = '/home/samhitha/code/Trials Before Experiments/2026.04.06_21.12.11_Joints2_Chains1_Seed44/0/BFS Outward Shortest n repetitions/checkpoints/iter_2_2.dill'

def exampleInspectAndVisualize(checkpoint_path=CHECKPOINT_FILE):
    info = inspectCheckpoint(checkpoint_path, load_tree=True)

    if 'tree' in info:
        info['tree'].show()
    
    return info

def exampleLoadCheckpoint(checkpoint_path=CHECKPOINT_FILE):

    state = loadCheckpoint(checkpoint_path)
    
    # Access the tree
    tree = state.tree
    print(f"Tree total length: {tree.totalLength():.4f}")
    print(f"Number of joints: {len(tree.Joints)}")
    
    # Access traversal info
    trav = state.traversal_info
    print(f"\nTraversal: {trav.get('algorithm')} ({trav.get('direction')}, {trav.get('orderBy')})")
    print(f"Current joint index: {trav.get('current_index')}")
    print(f"Repeat: {trav.get('repeat_index')} / {trav.get('total_repeats')}")
    print(f"Completed joints: {trav.get('completed_indices')}")
    
    # Access optimization history
    hist = state.optimization_history
    print(f"\nOptimization history:")
    print(f"  Iterations: {len(hist.get('times', []))}")
    print(f"  Latest time: {hist.get('times', [None])[-1]:.2f}s" if hist.get('times') else "  No times recorded")
    print(f"  Latest length: {hist.get('lengths', [None])[-1]:.4f}" if hist.get('lengths') else "  No lengths recorded")
    latest_loss = hist.get('losses', [None])[-1]
    latest_loss = latest_loss.item() if latest_loss.size == 1 else latest_loss.tolist()
    print(f"  Latest loss: {latest_loss}")
    
    # Access configurations
    print(f"\nConfigurations tested: {len(state.configurations) if state.configurations else 0}")
    print(f"Failure penalty: {state.failure_penalty}")
    
    # Access metadata
    meta = state.metadata
    print(f"\nMetadata:")
    print(f"  Timestamp: {meta.get('timestamp')}")
    print(f"  Iteration count: {meta.get('iteration_count')}")
    print(f"  Version: {meta.get('version')}")
    
    return state

def exampleCheckCollisions(checkpoint_path=CHECKPOINT_FILE):
    state = loadCheckpoint(checkpoint_path)
    
    # Use the saved collision matrices
    if state.collision_matrices is not None:
        count, error = state.tree.collisionsCountAndError(None, state.collision_matrices)
        print(f"Collisions at checkpoint: count={count}, error={error:.6f}")
    else:
        # Recompute if not saved
        matrices = state.tree.buildCollisionMatrices()
        count, error = state.tree.collisionsCountAndError(None, matrices)
        print(f"Collisions (recomputed): count={count}, error={error:.6f}")
    
    return count, error

def exampleLoadWithoutRestoringRandom(checkpoint_path=CHECKPOINT_FILE):
    state = loadCheckpoint(checkpoint_path, restore_random_state=False)
    print("Loaded checkpoint without restoring random state")
    return state

def exampleResumeFromCheckpoint(checkpoint_path=CHECKPOINT_FILE):
    tree = resumeFromCheckpoint(checkpoint_path, verbose=True)
    return tree

def findWhenCollisionAppeared():
    checkpoints = listCheckpoints(CHECKPOINT_DIR)
    
    for cp in checkpoints:
        state = loadCheckpoint(cp['path'], restore_random_state=False)
        matrices = state.collision_matrices or state.tree.buildCollisionMatrices()
        count, _ = state.tree.collisionsCountAndError(None, matrices)
        
        if count > 0:
            print(f"First collision at repeat {cp['repeat']}, joint {cp['joint_idx']}")
            print(f"Checkpoint: {cp['path']}")
            return cp
    
    print("No collisions found in any checkpoint")
    return None

def plotLengthProgression():
    
    checkpoints = listCheckpoints(CHECKPOINT_DIR)
    
    iterations = []
    lengths = []
    
    for i, cp in enumerate(checkpoints):
        state = loadCheckpoint(cp['path'], restore_random_state=False)
        iterations.append(i)
        lengths.append(state.tree.totalLength())
    
    plt.figure(figsize=(10, 6))
    plt.plot(iterations, lengths, 'b-o')
    plt.xlabel('Iteration')
    plt.ylabel('Total Tree Length')
    plt.title('Optimization Progress')
    plt.grid(True)
    plt.savefig('optimization_progress.png')
    plt.show()
    
    return iterations, lengths

if __name__ == "__main__":
    # Example usage:
    # list = listCheckpoints(CHECKPOINT_DIR)
    # info = exampleInspectAndVisualize()
    # state = exampleLoadCheckpoint()
    # count, error = exampleCheckCollisions()
    # state_no_random = exampleLoadWithoutRestoringRandom()
    # tree = exampleResumeFromCheckpoint()
    # collision_cp = findWhenCollisionAppeared()
    iterations, lengths = plotLengthProgression()