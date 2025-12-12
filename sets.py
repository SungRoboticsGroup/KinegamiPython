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

jointCount = 20
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


def findWaypointSets(tree):

    visited = set()
    waypoint_sets = []
    
    def buildWaypointSet(joint_idx, current_set):

        if joint_idx in visited or joint_idx is None or joint_idx == -1:
            return
        
        # Check if this joint is a waypoint
        if not isWaypoint(tree.Joints[joint_idx]):
            return
        
        # Add to current set and mark as visited
        current_set.append(joint_idx)
        visited.add(joint_idx)
        
        # Check parent
        parent_idx = tree.Parents[joint_idx]
        if parent_idx is not None and parent_idx != -1 and parent_idx not in visited:
            if isWaypoint(tree.Joints[parent_idx]):
                buildWaypointSet(parent_idx, current_set)
        
        # Check children
        children_indices = tree.Children[joint_idx]
        for child_idx in children_indices:
            if child_idx not in visited and isWaypoint(tree.Joints[child_idx]):
                buildWaypointSet(child_idx, current_set)
    
    # Go through all joints and find waypoint sets
    for joint_idx in range(len(tree.Joints)):
        if joint_idx not in visited and isWaypoint(tree.Joints[joint_idx]):
            current_set = []
            buildWaypointSet(joint_idx, current_set)
            if current_set:
                waypoint_sets.append(current_set)
    
    return waypoint_sets


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


# Main execution
print(f"\n\nConstructing tree...")
np.random.set_state(saved_state)
construct = generateTree(jointCount)

# Print structure with link lengths
print("\n=== Tree Structure with Link Lengths ===\n")
for parent in range(len(construct.Children)):
    children = construct.Children[parent]
    if not children:
        continue

    print(f"- Joint {parent} -> {"Waypoint" if isWaypoint(construct.Joints[parent]) else "Joint"}")
    for child in children:
        link = construct.Links[child]
        path = link.path
        print(f"   ---> Joint {child} | Length: {path.length:.2f}")
    print("")

# Find sets of adjacent waypoints
waypoint_sets = findWaypointSets(construct)

print(f"\nFound {len(waypoint_sets)} waypoint sets:")
for set_idx, waypoint_set in enumerate(waypoint_sets):
    print(f"  Set {set_idx}: {waypoint_set}")

# Visualize the tree
print("\nVisualizing tree...")
construct.show(block=True)