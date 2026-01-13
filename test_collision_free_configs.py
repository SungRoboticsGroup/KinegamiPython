#!/usr/bin/env python3
"""
Test script to generate configurations close to neutral and filter out
those with initial collisions.
"""

from KinematicTree import *
from makeKinematicTree import *
import random
import numpy as np
import copy

def generate_random_configs(tree, num_configs=100):
    """
    Generate random configurations using tree.randomConfiguration().
    Uses the built-in method that properly handles joint ranges.
    
    Args:
        tree: KinematicTree to generate configurations for
        num_configs: Number of configurations to generate
        
    Returns:
        List of configuration arrays
    """
    configs = []
    num_joints = len(tree.Joints)
    
    # Always include the neutral configuration
    neutral = [0] * num_joints
    configs.append(neutral)
    
    print(f"Generating {num_configs} random configurations using tree.randomConfiguration()...\n")
    
    for config_idx in range(num_configs - 1):
        # Use the tree's built-in randomConfiguration method
        config = tree.randomConfiguration(realJointsOnly=False).tolist()
        configs.append(config)
    
    return configs

def check_configuration_collisions(tree, config):
    """
    Check if a configuration has collisions.
    
    Args:
        tree: KinematicTree
        config: Configuration array
        
    Returns:
        Number of collisions detected
    """
    test_tree = copy.deepcopy(tree)
    
    # Set the configuration
    for i in range(len(test_tree.Joints)):
        if not isWaypoint(test_tree.Joints[i]):
            test_tree.setJointState(i, config[i])
        test_tree.Joints[i].recomputeCollisionCapsules()
    
    # Detect collisions
    collision_count = test_tree.detectCollisions(debug=True)
    return collision_count

def filter_collision_free_configs(tree, configs, verbose=True):
    """
    Filter out configurations with collisions.
    
    Args:
        tree: KinematicTree
        configs: List of configuration arrays
        verbose: Print detailed information
        
    Returns:
        List of collision-free configurations
    """
    collision_free = []
    
    print("=" * 60)
    print("TESTING CONFIGURATIONS FOR COLLISIONS")
    print("=" * 60)
    
    for idx, config in enumerate(configs):
        collision_count = check_configuration_collisions(tree, config)
        
        if collision_count == 0:
            collision_free.append(config)
            status = "✓ COLLISION-FREE"
        else:
            status = f"✗ HAS {collision_count} COLLISION(S)"
        
        if verbose:
            print(f"\nConfig {idx}: {status}")
            
            # Show non-zero joint values
            non_zero_joints = []
            for i in range(len(config)):
                if not isWaypoint(tree.Joints[i]) and abs(config[i]) > 1e-6:
                    min_state, max_state = tree.Joints[i].stateRange()
                    non_zero_joints.append(f"J{i}={config[i]:.3f} (range: [{min_state:.3f}, {max_state:.3f}])")
            
            if non_zero_joints:
                print(f"  Non-zero joints: {', '.join(non_zero_joints)}")
            else:
                print(f"  All joints at neutral (0)")
    
    print("\n" + "=" * 60)
    print(f"RESULTS: {len(collision_free)}/{len(configs)} configurations are collision-free")
    print("=" * 60)
    
    return collision_free

def print_collision_free_configs(configs):
    """
    Print collision-free configurations in a format that can be copied.
    """
    print("\n" + "=" * 60)
    print("COLLISION-FREE CONFIGURATIONS (Python list format)")
    print("=" * 60)
    print("\nconfigurations = [")
    for idx, config in enumerate(configs):
        config_str = "[" + ", ".join(f"{val:.6f}" for val in config) + "]"
        if idx < len(configs) - 1:
            print(f"    {config_str},")
        else:
            print(f"    {config_str}")
    print("]")
    print()

def generate_test_tree(nJoints=4, cubeSize=10, probabilityOfBranching=0.5, seed=43):
    """
    Generate a random kinematic tree for testing.
    Based on the generateTree function from randomTree.py
    
    Args:
        nJoints: Number of joints in the tree
        cubeSize: Size of the cube to place joints randomly
        probabilityOfBranching: Probability of creating a branch vs extending a chain
        seed: Random seed for reproducibility
    """
    np.random.seed(seed)
    
    bounds = (-cubeSize/2, cubeSize/2)
    poses = [ SE3.Rand(xrange=bounds, yrange=bounds, zrange=bounds)
                for _ in range(nJoints) ]

    r = 1
    numSides = 4
    neutralLength = 3

    root = RevoluteJoint(numSides, r, np.pi, poses[0]) if np.random.rand() > 0.5 \
            else PrismaticJoint(numSides, r, neutralLength, 3, np.pi/5, poses[0])

    specTree = JointSpecificationTree(root)

    for i in range(1, nJoints):
        branching = np.random.rand() < probabilityOfBranching and i > 1
        if branching:
            # randomly select a parent from among the non-leaves
            nonLeaves = specTree.nonLeaves()
            parent = nonLeaves[np.random.randint(0, len(nonLeaves))]
        else:
            # randomly select a parent from among the leaves
            leaves = specTree.leaves()
            parent = leaves[np.random.randint(0, len(leaves))]
        
        newJoint = RevoluteJoint(numSides, r, np.pi, poses[i]) if np.random.rand() > 0.5 \
            else PrismaticJoint(numSides, r, neutralLength, 3, np.pi/5, poses[i])
        specTree.addJoint(parent, newJoint)

    initialTree = makeTubularKinematicTree(specTree)
    return initialTree

def main():
    """
    Main test function - generates a tree and tests configurations.
    """
    print("=" * 60)
    print("COLLISION-FREE CONFIGURATION GENERATOR")
    print("=" * 60)
    
    # Generate a tree
    print("\nGenerating tree...")
    nJoints = 4  # Number of joints - adjust as needed
    probabilityOfBranching = 0.3  # Probability of branching
    seed = 43  # Random seed for reproducibility
    
    tree = generate_test_tree(nJoints=nJoints, 
                              cubeSize=10, 
                              probabilityOfBranching=probabilityOfBranching,
                              seed=seed)
    
    print(f"Generated tree with {len(tree.Joints)} joints")
    print(f"Random seed: {seed}")
    print()
    
    # Print joint information
    print("Joint types and ranges:")
    for i in range(len(tree.Joints)):
        if not isWaypoint(tree.Joints[i]):
            joint = tree.Joints[i]
            min_state, max_state = joint.stateRange()
            joint_type = type(joint).__name__
            print(f"  Joint {i} ({joint_type}): range [{min_state:.3f}, {max_state:.3f}]")
    print()
    
    # Generate random configurations using tree.randomConfiguration()
    num_configs = 100  # Number of configurations to test
    
    configs = generate_random_configs(tree, num_configs)
    
    # Filter out configurations with collisions
    collision_free_configs = filter_collision_free_configs(tree, configs, verbose=False)
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"Found {len(collision_free_configs)} collision-free configs out of {len(configs)} tested")
    print(f"Success rate: {len(collision_free_configs)/len(configs)*100:.1f}%")
    print(f"{'='*60}\n")
    
    # Print the collision-free configurations in a copyable format
    if collision_free_configs:
        print_collision_free_configs(collision_free_configs)
        
        # Save to file
        with open("collision_free_configs.txt", "w") as f:
            f.write(f"# Collision-free configurations\n")
            f.write(f"# Generated using tree.randomConfiguration()\n")
            f.write(f"# Success rate: {len(collision_free_configs)/len(configs)*100:.1f}%\n\n")
            f.write("configurations = [\n")
            for idx, config in enumerate(collision_free_configs):
                config_str = "    [" + ", ".join(f"{val:.6f}" for val in config) + "]"
                if idx < len(collision_free_configs) - 1:
                    f.write(config_str + ",\n")
                else:
                    f.write(config_str + "\n")
            f.write("]\n")
        print("✓ Saved to collision_free_configs.txt")
    else:
        print("\n⚠ WARNING: No collision-free configurations found!")
        print("Try testing more configurations or using a different tree.")

if __name__ == "__main__":
    main()
