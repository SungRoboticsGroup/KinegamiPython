from optimizationFunctions import *
from testqtgraph import *
from makeKinematicTree import *
from KinematicTree import *

import random
import copy
from collections import defaultdict, deque
import json
from functools import partial
import threading
import os 
import shutil
import re
from datetime import datetime

def findWaypointSets(tree : KinematicTree):
    """
    Find sets of waypoints that are adjacent to each other without any real joints in between,
    the corresponding sets of links adjacent to those waypoints (both incoming and outgoing), 
    and a dictionary mapping each link index to the set it belongs to.
    
    Example: a tree has a joint at index 3 which branches out into two chains:
        waypoint at index 4, waypoint at index 5, real joint at index 6
        waypoint at index 7, real joint at index 8
    This would yield waypoint set {4, 5, 7} and link set {4, 5, 6, 7, 8}
    
    :param tree: KinematicTree to analyze
    :return: (two lists, one dictionary)
             - List of waypoint sets (each set is a set of joint indices)
             - List of link sets (each set is a set of link indices)
             - Dictionary mapping each link index to the set it belongs to
    """
    # visited = [False] * len(tree.Joints)
    waypoint_sets = []
    link_sets = []
    link_dictionary = {}
    

    def addToWaypointCluster(joint_idx : int, set_idx : int):
        if set_idx == -1:
            set_idx = len(waypoint_sets)
            waypoint_sets.append(set())
            link_sets.append(set())
        waypoint_sets[set_idx].add(joint_idx)
        return set_idx
    
    def addToLinkCluster(link_idx : int, set_idx : int):
        if set_idx == -1:
            set_idx = len(link_sets)
            waypoint_sets.append(set())
            link_sets.append(set())
        link_sets[set_idx].add(link_idx)
        link_dictionary[link_idx] = set_idx
        return set_idx

    def findSetsFromJoint(joint_idx : int, set_idx : int = -1):
        if isWaypoint(tree.Joints[joint_idx]):
            set_idx = addToWaypointCluster(joint_idx, set_idx)

        for c in tree.Children[joint_idx]:
            set_idx = addToLinkCluster(c, set_idx)
            findSetsFromJoint(c, set_idx if isWaypoint(tree.Joints[c]) else -1)
            
    findSetsFromJoint(0)

    return waypoint_sets, link_sets, link_dictionary


    

def buildCollisionPairDictionary(tree):
    """
    Build a dictionary mapping each joint/waypoint index to the set of collision pairs
    that need to be checked when that joint/waypoint moves.
    
    Rules:
    - Waypoints don't collision check directly
    - Link/Link: check iff links are separated by a real joint (not in same waypoint cluster)
    - Real Joint / Real Joint: always check
    - Real Joint / Link: always check
    
    Returns:
        dict: Keys are joint indices, values are lists of tuples:
              [((idx1, type1), (idx2, type2)), ...]
              where type is 'joint' or 'link'
    """
    waypoint_sets, link_sets, whichLinkSet = findWaypointSets(tree)
    
    # Pre-classify all joints
    real_joints = [i for i in range(len(tree.Joints)) if not isWaypoint(tree.Joints[i])]
    
    collision_pairs = {}
    
    for node_idx in range(len(tree.Joints)):

        # Waypoint
        if isWaypoint(tree.Joints[node_idx]):

            # Doesn't need collision checking with anything
            collision_pairs[node_idx] = []
            continue
        
        pairs = []
        
        # Real joint
        for other_idx in real_joints:

            # Against other real joints
            if other_idx != node_idx:
                pairs.append(((node_idx, 'joint'), (other_idx, 'joint')))

            # Against links, including the ones it's adjacent to
            pairs.append(((node_idx, 'joint'), (other_idx, 'link')))
        
        # Link
        set_idx = whichLinkSet.get(node_idx)
        for other_idx in range(len(tree.Links)):

            # Skip itself
            if other_idx == node_idx:
                continue

            # If comparing with a different link, check which set that one belongs to
            other_set = whichLinkSet.get(other_idx)

            # Check if both links are in different sets or the one in question is not in a set at all
            if set_idx != other_set or set_idx is None:
                pairs.append(((node_idx, 'link'), (other_idx, 'link')))
        
        collision_pairs[node_idx] = pairs
    
    return collision_pairs


def generate_colors(x, cmap_name="rainbow"):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(i / (x - 1))[:3] for i in range(x)]

def plotColoredTrees(directory, collection = []):
    filenames = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".tree")]

    final_file = [f for f in filenames if 'final.tree' in f]
    other_files = [f for f in filenames if 'final.tree' not in f]
    other_files.sort(key=lambda x: float(re.search(r'/([\d.]+)_', x).group(1)))
    sorted_files = other_files + final_file
    
    trees =  [loadKinematicTree(f) for f in sorted_files]

    plottingColors = generate_colors(len(collection) if (len(collection) > 0) else len(trees))
    
    def plotCollection(col):
        labels = []
        ax = plt.figure().add_subplot(projection='3d')
        for ct, i in enumerate(col):
            color = plottingColors[ct]
            trees[i].addToPlot(ax, jointColor=color, jointEdgeColor=color, linkColor=color, surfaceOpacity=0.5,showLinkPath=False,showJointPoses=False)
            labels.append(plt.Line2D([0], [0], color=color, lw=4, label=str(i)))
    
        ax.set_aspect('equal')
        ax.legend(handles=labels, loc='upper right', fontsize='small')
        plt.axis('off')
        #plt.savefig(directory + "graph_output.png", dpi=800)
        plt.show(block=True)

    if (collection == None):
        for i in range(0,len(trees)):
            plotCollection([i, i+1])
        
        trees[-1].show()
    else:
        plotCollection(collection)

"""
# Setup for testing
jointCount = 12
probabilityOfBranching = 0.5
sparse = False
cubeSize = 100 if sparse else 10
title = str(jointCount)+" Joint " + ("Chains" if probabilityOfBranching == 0 else "Trees") + (" Sparse" if sparse else " Dense")
treeCount = 1
restartFrom = 0
multipleIterations=False

seed = 42
np.random.seed(seed)
saved_state = np.random.get_state()

timestamp = datetime.now().strftime("%Y.%m.%d_%H.%M.%S")
base_dir = "Trials Before Experiments"
experiment_name = f"{timestamp}_Joints{jointCount}_Trees{treeCount}_Seed{seed}"
results_dir = os.path.join(base_dir, experiment_name)
os.makedirs(results_dir, exist_ok=True)


def generateTree(nJoints):
    bounds = (-cubeSize/2, cubeSize/2)
    poses = [ SE3.Rand(xrange=bounds, yrange=bounds, zrange=bounds)
                for _ in range(nJoints) ]

    r = 1
    numSides = 4
    neutralLength = 3

    root = RevoluteJoint(numSides,r,np.pi,poses[0]) if np.random.rand() > 0.5 \
            else PrismaticJoint(numSides,r,neutralLength,3,np.pi/5,poses[0])

    specTree = JointSpecificationTree(root)

    for i in range(1,nJoints):
        #parent = np.random.randint(int((i - 1) * (1 - branchingRatio)), i)
        branching = np.random.rand() < probabilityOfBranching and i > 1
        if branching:
            # randomly select a parent from among the non-leaves
            nonLeaves = specTree.nonLeaves()
            parent = nonLeaves[np.random.randint(0,len(nonLeaves))]
        else:
            # randomly select a parent from among the leaves
            leaves = specTree.leaves()
            parent = leaves[np.random.randint(0,len(leaves))]
        
        newJoint = RevoluteJoint(numSides,r,np.pi,poses[i]) if np.random.rand() > 0.5 \
            else PrismaticJoint(numSides,r,neutralLength,3,np.pi/5,poses[i])
        specTree.addJoint(parent, newJoint)

    initialTree = makeTubularKinematicTree(specTree)
    assert(abs(specTree.totalLengthLowerBound() - initialTree.totalLengthLowerBound()) < 1e-5)
    return initialTree


# Main execution for testing
print(f"\n\nConstructing tree...")
np.random.set_state(saved_state)
construct = generateTree(jointCount)

# Print structure with link lengths
print("\n=== Tree Structure with Link Lengths ===\n")
# Function to recursively print the tree structure,
# broken down into sub-chains between wherever it branches
def printTreeStructure(tree, joint_idx, indent=""):
    joint = tree.Joints[joint_idx]
    children = tree.Children[joint_idx]
    
    print(f"{indent}- Joint {joint_idx} ({'Waypoint' if isWaypoint(joint) else 'Joint'})")
    
    for child_idx in children:
        link = tree.Links[child_idx]
        path = link.path
        print(f"{indent}   ---> Link to Joint {child_idx} | Length: {path.length:.2f}")
        printTreeStructure(tree, child_idx, indent + "\t")

printTreeStructure(construct, 0)

# Find sets of adjacent waypoints
waypoint_sets, link_sets, link_to_set = findWaypointSets(construct)

# Each link set should have a corresponding color for visualization
colors = generate_colors(len(link_sets))
linkColorList = [colors[link_to_set[i]] if i in link_to_set else (0.5,0.5,0.5) for i in range(len(construct.Links))]


print(f"Waypoint Sets: {waypoint_sets}")
print(f"Link Sets: {link_sets}")
print(f"Link to Set Mapping: {link_to_set}")


# Visualize the tree
print("\nVisualizing tree...")
construct.show(block=True, linkColor=linkColorList)
"""