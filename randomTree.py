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

jointCount = 4
probabilityOfBranching = 0.5
sparse = False
cubeSize = 100 if sparse else 10
title = str(jointCount)+" Joint " + ("Chains" if probabilityOfBranching == 0 else "Trees") + (" Sparse" if sparse else " Dense")
treeCount = 1
restartFrom = 0
multipleIterations=False

seed = 44
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

def find_collision_free_configs(tree, max_attempts=100, num_configs_needed=2):
    configs = [[0] * len(tree.Joints)] 
    
    for attempt in range(max_attempts):
        if len(configs) >= num_configs_needed+1:
            break
            
        # Generate random config
        config = tree.randomConfiguration(realJointsOnly=False).tolist()
        
        # Test for collisions
        test_tree = copy.deepcopy(tree)
        for i in range(len(test_tree.Joints)):
            if not isWaypoint(test_tree.Joints[i]):
                test_tree.setJointState(i, config[i])
            test_tree.Joints[i].recomputeCollisionCapsules()
        
        if test_tree.detectCollisions(debug=False) == 0:
            configs.append(config)
    
    print(f"Using {len(configs)} collision-free configs")
    return configs

def testRandomTrees(show: bool = False):
    optimizations = [
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="dfs", direction="outward", orderBy="longest"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="dfs", direction="outward", orderBy="shortest"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="dfs", direction="inward", orderBy="longest"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="dfs", direction="inward", orderBy="shortest"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="bfs", direction="outward", orderBy="longest"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="bfs", direction="outward", orderBy="shortest"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="bfs", direction="inward", orderBy="longest"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="bfs", direction="inward", orderBy="shortest"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="randomized", power=3)
    ]

    labels = [
            "DFS - Outward Longest",
            "DFS - Outward Shortest",
            "DFS - Inward Longest",
            "DFS - Inward Shortest",
            "BFS - Outward Longest",
            "BFS - Outward Shortest",
            "BFS - Inward Longest",
            "BFS - Inward Shortest",
            "Randomized"
    ]

    lowerBounds = []
    results = []
    outputs = []
    constructs = []

    restartDir = os.path.join(results_dir, f"trial_{restartFrom}")
    if os.path.exists(restartDir):
        shutil.rmtree(restartDir)


    for i in range(restartFrom,treeCount):
        print(f"\n\nConstructing tree {i}")
        # Reset random state before each trial to ensure deterministic behavior
        np.random.set_state(saved_state)
        construct = generateTree(jointCount)
        tree_save_path = os.path.join(results_dir, f"{i}.txt")
        
        # Save using repr() representation
        with open(tree_save_path, 'w') as f:
            f.write(repr(construct))
            
        construct.save(os.path.join(results_dir, str(i)), saveDir=False)
        lowerBounds.append(construct.totalLengthLowerBound())
        results.append([])
        outputs.append([])
        constructs.append(construct)
        
        # Generate collision-free configs for this specific tree
        collision_free_configs = find_collision_free_configs(construct, max_attempts=100, num_configs_needed=2)
        
        for no, f in enumerate(optimizations):
            print(f"\nTrying loss function {no}")
            # Create a subdirectory for this trial's results
            trial_dir = os.path.join(results_dir, f"trial_{i}", labels[no])
            os.makedirs(trial_dir, exist_ok=True)
            
            # Reset random state before each optimization to ensure deterministic behavior
            np.random.set_state(saved_state)
            optimized, times, losses = f(construct, 
                                        configurations=collision_free_configs,
                                        showSteps=False, 
                                        parallelize=True, 
                                        evaluate=True, 
                                        verbose=False, 
                                        directory=trial_dir)
            
            if multipleIterations:
                count = 2
                while (losses[0] - losses[-1] > 100):
                    print(f"Trying loss function {no} for the {count}th time")
                    # Reset random state before each retry
                    np.random.set_state(saved_state)
                    optimized, times, losses = f(construct, showSteps=False, parallelize=True, evaluate=True, verbose=False, directory=None)
                    count += 1
            #print(optimized.detectCollisions(plot=True, includeEnds=False, debug=True))
            results[i].append((times, losses))
            outputs[i].append(optimized)
        
        # Save checkpoint to results directory
        checkpoint_file = os.path.join(results_dir, f"random_results_chkpt{i}.json")
        with open(checkpoint_file, "w") as file:
            json.dump(results, file)

    # Save final results to results directory
    final_results_file = os.path.join(results_dir, "random_results.json")
    with open(final_results_file, "w") as file:
        json.dump(results, file)


    for index, result in enumerate(results):
        if show:
            constructs[index].show(block=False)
            for output in outputs[index]:
                output.show(block=False)


        plt.figure(figsize=(8, 5)) 
        
        idx = 0
        for x, y in result:
            plt.plot(x, y, marker='o', linestyle='-', label=labels[idx])
            idx += 1

        # plot the lower bound as a horizontal line
        plt.axhline(y=lowerBounds[index], color='r', linestyle='-', label='Lower Bound')

        plt.xlabel('Time')
        plt.ylabel('Loss')
        plt.title('Loss vs Time: ' + str(title))
        plt.legend()

        plt.grid(True)
        # save to an image file
        plt.savefig(os.path.join(results_dir, f"example{index}.png"))
        plt.show()

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


testRandomTrees(show=True)