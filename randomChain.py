from optimizationFunctions import *
from testqtgraph import *
from makeKinematicTree import *
from KinematicTree import *
from KinematicChain import *

import random
import copy
from collections import defaultdict, deque
import json
from functools import partial
import threading
import os 
import shutil
from datetime import datetime

jointCount = 3
sparse = False
cubeSize = 10
title = str(jointCount)+" Joint Generalized Gimbal Chains Cube Size " + str(cubeSize)
chainCount = 1  
restartFrom = 0
multipleIterations=False

seed = 42
np.random.seed(seed)
saved_state = np.random.get_state()

timestamp = datetime.now().strftime("%Y.%m.%d_%H.%M.%S")
base_dir = "Trials Before Experiments"
experiment_name = f"{timestamp}_Joints{jointCount}_Chains{chainCount}_Seed{seed}"
results_dir = os.path.join(os.getcwd(), base_dir, experiment_name)
os.makedirs(results_dir, exist_ok=True)

def generateRandomChain(nJoints):
    bounds = (-cubeSize/2, cubeSize/2)
    poses = [ SE3.Rand(xrange=bounds, yrange=bounds, zrange=bounds)
                for _ in range(nJoints) ]

    r = 1
    numSides = 4
    neutralLength = 3

    root = RevoluteJoint(numSides,r,np.pi,poses[0]) if np.random.rand() > 0.5 \
            else PrismaticJoint(numSides,r,neutralLength,3,np.pi/5,poses[0])

    chain = KinematicChain(root, gimbal=True)
    for i in range(1,nJoints):
        newJoint = RevoluteJoint(numSides,r,np.pi,poses[i]) if np.random.rand() > 0.5 \
            else PrismaticJoint(numSides,r,neutralLength,3,np.pi/5,poses[i])
        chain.appendGeneralizedGimbal(newJoint)

    return chain

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

def test():
    optimizations = [
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="dfs", direction="outward", orderBy="longest", repeatTraversal="n"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="dfs", direction="outward", orderBy="shortest"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="dfs", direction="inward", orderBy="longest", repeatTraversal="n"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="dfs", direction="inward", orderBy="shortest"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="bfs", direction="outward", orderBy="longest", repeatTraversal="n"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="bfs", direction="outward", orderBy="shortest"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="bfs", direction="inward", orderBy="longest", repeatTraversal="n"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="bfs", direction="inward", orderBy="shortest"),
                    partial(optimizeTree, childFraction=1, guarantee=True, traversal="randomized", power=3, repeatTraversal="n")
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
            "Randomized Weighted"
    ]

    lowerBounds = []
    results = []

    restartDir = "sim_results/" + title + "/" + str(restartFrom)
    if os.path.exists(restartDir):
        shutil.rmtree(restartDir)


    for i in range(restartFrom,chainCount):
        print(f"\n\nConstructing tree {i}")
        # Reset random state before each trial to ensure deterministic behavior
        np.random.set_state(saved_state)
        construct = generateRandomChain(jointCount)
        tree_save_path = os.path.join(results_dir, f"{i}.txt")
        
        # Save using repr() representation
        with open(tree_save_path, 'w') as f:
            f.write(repr(construct))

        lowerBounds.append(construct.totalLengthLowerBound())
        results.append([])
        
        # Generate collision-free configs for this specific chain
        # This shouldn't be necessary since we're using the gimbal construction...
        collision_free_configs = find_collision_free_configs(construct, max_attempts=100, num_configs_needed=2)
        
        for no, f in enumerate(optimizations):
            print(f"\nTrying loss function {no}")
            direc = os.path.join(results_dir, str(i), labels[no])
            os.makedirs(direc, exist_ok=True)
            # Reset random state before each optimization to ensure deterministic behavior
            np.random.set_state(saved_state)
            optimized, times, losses = f(construct, 
                                        configurations=collision_free_configs,
                                        showSteps=False, 
                                        parallelize=True, 
                                        evaluate=True, 
                                        verbose=False, 
                                        directory=direc)
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
        
        # Save checkpoint to results directory
        checkpoint_file = os.path.join(results_dir, f"random_results_chkpt{i}.json")
        with open(checkpoint_file, "w") as file:
            json.dump(results, file)

    # Save final results to results directory
    final_results_file = os.path.join(results_dir, "random_results.json")
    with open(final_results_file, "w") as file:
        json.dump(results, file)


    for index, result in enumerate(results):
        plt.figure(figsize=(8, 5)) 
        
        idx = 0
        for x, y in result:
            plt.plot(x, y, marker='o', linestyle='-', label=labels[idx])
            idx += 1

        # plot the lower bound as a horizontal line
        plt.axhline(y=lowerBounds[index], color='r', linestyle='-', label='Lower Bound')

        plt.xlabel('Time')
        plt.ylabel('Length')
        #plt.title('Optimization Traversals: ' + str(title))
        plt.legend()

        plt.grid(True)
        # save to an image file in the results directory
        plot_save_path = os.path.join(results_dir, f"plot_{index}.png")
        try:
            plt.savefig(plot_save_path, dpi=300, bbox_inches='tight')
            print(f"Saved plot to: {plot_save_path}")
        except Exception as e:
            print(f"Error saving plot: {e}")

        plt.show()
        plt.close()

def generate_colors(x, cmap_name="rainbow"):
    cmap = plt.get_cmap(cmap_name)
    return [cmap(i / (x - 1))[:3] for i in range(x)]

def plotColoredTrees(directory):
    files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".tree")]
    trees =  [loadKinematicTree(f) for f in files]
    plottingColors = generate_colors(len(trees))
    
    ax = plt.figure().add_subplot(projection='3d')

    labels = []
    for i in range(0,len(trees)):
        color = plottingColors[i]
        trees[i].addToPlot(ax, jointColor=color, jointEdgeColor=color, linkColor=color, surfaceOpacity=0.07,showLinkPath=False,showJointPoses=False)
        labels.append(plt.Line2D([0], [0], color=color, lw=4, label=str(i)))
    
    ax.set_aspect('equal')
    ax.legend(handles=labels, loc='upper right', fontsize='small')
    plt.axis('off')
    #plt.savefig(directory + "graph_output.png", dpi=800)
    plt.show(block=True)


test()