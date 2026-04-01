"""
Debug script to display what collision pairs are checked for each joint.
This helps visualize why some collisions may be missed during optimization.

Usage:
    python debug_collision_pairs.py <path_to_tree_file>
    
Example:
    python debug_collision_pairs.py "save/my_tree.tree"
"""

from KinematicTree import *
from makeKinematicTree import *
import numpy as np
import sys


def display_collision_pairs(tree, joint_index=None):
    """
    Display collision pairs that would be checked for a specific joint index,
    or for all joints if joint_index is None.
    """
    collisionMatrices = tree.buildCollisionMatrices()
    
    print(f"Tree has {len(tree.Joints)} joints")
    print(f"Children: {tree.Children}")
    print()
    
    if joint_index is not None:
        joints_to_check = [joint_index]
    else:
        joints_to_check = range(len(tree.Joints))
    
    print("=" * 70)
    print("COLLISION PAIRS CHECKED PER JOINT (when movedJointIndex is set)")
    print("=" * 70)
    print()

    joints_to_check = [12]
    
    for idx in joints_to_check:
        jj_pairs, jl_pairs, ll_pairs = tree.collisionPairsFromMovingJoint(idx, collisionMatrices)
        
        print(f"Joint {idx} (children: {tree.Children[idx]}):")
        print(f"  Joint-Joint pairs ({len(jj_pairs.T) if jj_pairs.size > 0 else 0}):")
        if jj_pairs.size > 0:
            for pair in jj_pairs.T.tolist():
                print(f"    Joint {pair[0]} <-> Joint {pair[1]}")
        else:
            print("    (none)")
            
        print(f"  Joint-Link pairs ({len(jl_pairs.T) if jl_pairs.size > 0 else 0}):")
        if jl_pairs.size > 0:
            for pair in jl_pairs.T.tolist():
                print(f"    Joint {pair[0]} <-> Link {pair[1]}")
        else:
            print("    (none)")
            
        print(f"  Link-Link pairs ({len(ll_pairs.T) if ll_pairs.size > 0 else 0}):")
        if ll_pairs.size > 0:
            for pair in ll_pairs.T.tolist():
                print(f"    Link {pair[0]} <-> Link {pair[1]}")
        else:
            print("    (none)")
        print()
    
    # Show all pairs for comparison
    print("=" * 70)
    print("ALL COLLISION PAIRS (when movedJointIndex=None)")
    print("=" * 70)
    jj_all, jl_all, ll_all = tree.getAllCollisionPairs(collisionMatrices)
    
    print(f"\nAll Joint-Joint pairs ({len(jj_all.T) if jj_all.size > 0 else 0}):")
    if jj_all.size > 0:
        for pair in jj_all.T.tolist():
            print(f"  Joint {pair[0]} <-> Joint {pair[1]}")
    
    print(f"\nAll Joint-Link pairs ({len(jl_all.T) if jl_all.size > 0 else 0}):")
    if jl_all.size > 0:
        for pair in jl_all.T.tolist():
            print(f"  Joint {pair[0]} <-> Link {pair[1]}")
    
    print(f"\nAll Link-Link pairs ({len(ll_all.T) if ll_all.size > 0 else 0}):")
    if ll_all.size > 0:
        for pair in ll_all.T.tolist():
            print(f"  Link {pair[0]} <-> Link {pair[1]}")


def find_missing_pairs(tree, joint_index):
    """
    Find collision pairs that exist in ALL pairs but are NOT checked 
    when optimizing a specific joint.
    """
    collisionMatrices = tree.buildCollisionMatrices()
    
    jj_joint, jl_joint, ll_joint = tree.collisionPairsFromMovingJoint(joint_index, collisionMatrices)
    jj_all, jl_all, ll_all = tree.getAllCollisionPairs(collisionMatrices)
    
    def pairs_to_set(pairs):
        if pairs.size == 0:
            return set()
        return set(tuple(sorted(p)) for p in pairs.T.tolist())
    
    jj_joint_set = pairs_to_set(jj_joint)
    jl_joint_set = set(tuple(p) for p in jl_joint.T.tolist()) if jl_joint.size > 0 else set()
    ll_joint_set = pairs_to_set(ll_joint)
    
    jj_all_set = pairs_to_set(jj_all)
    jl_all_set = set(tuple(p) for p in jl_all.T.tolist()) if jl_all.size > 0 else set()
    ll_all_set = pairs_to_set(ll_all)
    
    missing_jj = jj_all_set - jj_joint_set
    missing_jl = jl_all_set - jl_joint_set
    missing_ll = ll_all_set - ll_joint_set
    
    print(f"\n{'=' * 70}")
    print(f"MISSING PAIRS when optimizing Joint {joint_index}")
    print(f"(These collisions will NOT be penalized!)")
    print(f"{'=' * 70}")
    
    print(f"\nMissing Joint-Joint pairs ({len(missing_jj)}):")
    for pair in sorted(missing_jj):
        print(f"  Joint {pair[0]} <-> Joint {pair[1]}")
    
    print(f"\nMissing Joint-Link pairs ({len(missing_jl)}):")
    for pair in sorted(missing_jl):
        print(f"  Joint {pair[0]} <-> Link {pair[1]}")
    
    print(f"\nMissing Link-Link pairs ({len(missing_ll)}):")
    for pair in sorted(missing_ll):
        print(f"  Link {pair[0]} <-> Link {pair[1]}")
    
    return missing_jj, missing_jl, missing_ll


def check_specific_collision(tree, obj1_type, obj1_idx, obj2_type, obj2_idx):
    """
    Check which joints would detect a specific collision pair.
    
    obj1_type/obj2_type: 'joint' or 'link'
    """
    collisionMatrices = tree.buildCollisionMatrices()
    
    print(f"\n{'=' * 70}")
    print(f"Which joints would detect: {obj1_type} {obj1_idx} <-> {obj2_type} {obj2_idx}?")
    print(f"{'=' * 70}")
    
    detecting_joints = []
    
    for joint_idx in range(len(tree.Joints)):
        jj_pairs, jl_pairs, ll_pairs = tree.collisionPairsFromMovingJoint(joint_idx, collisionMatrices)
        
        found = False
        
        if obj1_type == 'joint' and obj2_type == 'joint':
            if jj_pairs.size > 0:
                for pair in jj_pairs.T.tolist():
                    if set(pair) == {obj1_idx, obj2_idx}:
                        found = True
                        break
        elif obj1_type == 'joint' and obj2_type == 'link':
            if jl_pairs.size > 0:
                for pair in jl_pairs.T.tolist():
                    if pair[0] == obj1_idx and pair[1] == obj2_idx:
                        found = True
                        break
        elif obj1_type == 'link' and obj2_type == 'joint':
            if jl_pairs.size > 0:
                for pair in jl_pairs.T.tolist():
                    if pair[0] == obj2_idx and pair[1] == obj1_idx:
                        found = True
                        break
        elif obj1_type == 'link' and obj2_type == 'link':
            if ll_pairs.size > 0:
                for pair in ll_pairs.T.tolist():
                    if set(pair) == {obj1_idx, obj2_idx}:
                        found = True
                        break
        
        if found:
            detecting_joints.append(joint_idx)
    
    if detecting_joints:
        print(f"\nJoints that would detect this collision: {detecting_joints}")
    else:
        print(f"\nWARNING: NO joints would detect this collision!")
    
    print(f"\nJoints that would NOT detect this collision: {[j for j in range(len(tree.Joints)) if j not in detecting_joints]}")
    
    return detecting_joints


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        print("\nNo tree file specified. Using default example.")
        # Try to load a default tree
        try:
            tree = loadKinematicTree('/home/samhitha/code/Trials Before Experiments/2026.03.23_08.52.18_Joints8_Chains1_Seed44/0/DFS Inward Shortest n repetitions/10143.410218715668_10.tree')
        except:
            print("Could not load default tree. Please specify a tree file.")
            sys.exit(1)
    else:
        tree_path = sys.argv[1]
        print(f"Loading tree from: {tree_path}")
        tree = loadKinematicTree(tree_path)
    
    # Display all collision pairs for each joint
    display_collision_pairs(tree)
    
    # Example: Check which joints would detect a specific collision
    # Uncomment and modify as needed:
    # check_specific_collision(tree, 'joint', 0, 'link', 9)
    
    # Example: Find missing pairs for a specific joint
    # find_missing_pairs(tree, 3)
