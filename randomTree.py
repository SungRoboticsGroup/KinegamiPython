from optimizationFunctions import *
from testqtgraph import *
from makeKinematicTree import *
from KinematicTree import *
from treeTraversals import dfs, bfs

import random
from collections import defaultdict, deque
import json
from functools import partial
import threading
import os 
import shutil
import re

jointCount = 12
probabilityOfBranching = 0.5
sparse = False
cubeSize = 100 if sparse else 10
title = str(jointCount)+" Joint " + ("Chains" if probabilityOfBranching == 0 else "Trees") + (" Sparse" if sparse else " Dense")
treeCount = 1
restartFrom = 0
multipleIterations=False

os.makedirs("sim_results/" + title, exist_ok=True)


def generateTree(nJoints):
    bounds = (-cubeSize/2, cubeSize/2)
    poses = [ SE3.Rand(xrange=bounds, yrange=bounds, zrange=bounds)
                for _ in range(nJoints) ]

    r = 1
    numSides = 4
    neutralLength = 3

    root = RevoluteJoint(numSides,
                         r,
                         np.pi,
                         poses[0]) if np.random.rand() > 0.5 else PrismaticJoint(numSides,
                                                                                 r,
                                                                                 neutralLength,
                                                                                 3,
                                                                                 np.pi/5,
                                                                                 poses[0])

    specTree = JointSpecificationTree(root)

    for i in range(1, nJoints):
        #parent = np.random.randint(int((i - 1) * (1 - branchingRatio)), i)
        branching = np.random.rand() < probabilityOfBranching and i > 1
        if branching:
            # randomly select a parent from among the non-leaves
            nonLeaves = specTree.nonLeaves()
            parent = nonLeaves[np.random.randint(0, len(nonLeaves))]
        else:
            # randomly select a parent from among the leaves
            leaves = specTree.leaves()
            parent = leaves[np.random.randint(0, len(leaves))]
        
        newJoint = RevoluteJoint(numSides,
                                 r,
                                 np.pi,
                                 poses[i]) if np.random.rand() > 0.5 else PrismaticJoint(numSides,
                                                                                         r,
                                                                                         neutralLength,
                                                                                         3,
                                                                                         np.pi/5,
                                                                                         poses[i])
        specTree.addJoint(parent, newJoint)

    initialTree = makeTubularKinematicTree(specTree)
    assert(abs(specTree.totalLengthLowerBound() - initialTree.totalLengthLowerBound()) < 1e-5)
    return initialTree

def testRandomTrees():
    # optimizations = [partial(squaredOptimize, childFraction=0,streamline=True,guarantee=True),
    #                 partial(squaredOptimize, childFraction=0,streamline=True,guarantee=False),
    #                 partial(squaredOptimize, childFraction=0,streamline=False,guarantee=False),
    #                 partial(squaredOptimize, childFraction=0,streamline=False,guarantee=False,resetOnFail=False),
    #                 partial(squaredOptimize, childFraction=0,streamline=False,guarantee=True),
    #                 partial(linearOptimize, childFraction=0, streamline=False,guarantee=False),
    #                 partial(squaredOptimize, childFraction=1,streamline=True,guarantee=False),
    #                 partial(squaredOptimize, childFraction=1,streamline=True,guarantee=True),
    #                 partial(perpetualOptimize, iterations=jointCount * 2, childFraction=1)]

    # optimizations = [partial(linearOptimize, childFraction=1, streamline=False, guarantee=False, traversal=dfs, direction="outward", orderBy="longest"),
    #                  partial(linearOptimize, childFraction=1, streamline=False, guarantee=False, traversal=dfs, direction="outward", orderBy="shortest"),
    #                  partial(linearOptimize, childFraction=1, streamline=False, guarantee=False, traversal=bfs, direction="outward", orderBy="longest"),
    #                  partial(linearOptimize, childFraction=1, streamline=False, guarantee=False, traversal=bfs, direction="outward", orderBy="shortest"),
    #                  partial(linearOptimize, childFraction=1, streamline=False, guarantee=False, traversal=dfs, direction="inward", orderBy="longest"),
    #                  partial(linearOptimize, childFraction=1, streamline=False, guarantee=False, traversal=dfs, direction="inward", orderBy="shortest"),
    #                  partial(linearOptimize, childFraction=1, streamline=False, guarantee=False, traversal=bfs, direction="inward", orderBy="longest"),
    #                  partial(linearOptimize, childFraction=1, streamline=False, guarantee=False, traversal=bfs, direction="inward", orderBy="shortest")]

    optimizations = [partial(squaredOptimize, childFraction=1, streamline=False, guarantee=False, traversal=dfs, direction="outward", orderBy="longest"),
                     partial(squaredOptimize, childFraction=1, streamline=False, guarantee=False, traversal=dfs, direction="outward", orderBy="shortest"),
                     partial(squaredOptimize, childFraction=1, streamline=False, guarantee=False, traversal=bfs, direction="outward", orderBy="longest"),
                     partial(squaredOptimize, childFraction=1, streamline=False, guarantee=False, traversal=bfs, direction="outward", orderBy="shortest"),
                     partial(squaredOptimize, childFraction=1, streamline=False, guarantee=False, traversal=dfs, direction="inward", orderBy="longest"),
                     partial(squaredOptimize, childFraction=1, streamline=False, guarantee=False, traversal=dfs, direction="inward", orderBy="shortest"),
                     partial(squaredOptimize, childFraction=1, streamline=False, guarantee=False, traversal=bfs, direction="inward", orderBy="longest"),
                     partial(squaredOptimize, childFraction=1, streamline=False, guarantee=False, traversal=bfs, direction="inward", orderBy="shortest")]

    # optimizations = [partial(perpetualOptimize, iterations=jointCount * 2, weighted=True, showSteps=True, childFraction=1),
    #                  partial(perpetualOptimize, iterations=jointCount * 2, weighted=False, showSteps=True, childFraction=1)]

    # labels = ["Streamline + Guarantee (SG)", 
    #         "Streamline No Guarantee (SNG)",
    #         "No Streamline No Guarantee (NSNG)",
    #         "NSNG, No Reset on Fail",
    #         "No Streamline Guarantee (NSG)",
    #         "Linear (L)",
    #         "Equal Child No Guarantee (ECNG)",
    #         "Equal Child Guarantee (ECG)",
    #         "Perpetual (P)"]

    labels = ["DFS Outward Longest",
              "DFS Outward Shortest",
              "BFS Outward Longest",
              "BFS Outward Shortest",
              "DFS Inward Longest",
              "DFS Inward Shortest",
              "BFS Inward Longest",
              "BFS Inward Shortest"]

    # labels = ["Perpetual Weighted", "Perpetual Unweighted"]

    lowerBounds = []
    results = []

    restartDir = "sim_results/" + title + "/" + str(restartFrom)
    if os.path.exists(restartDir):
        shutil.rmtree(restartDir)


    for i in range(restartFrom, treeCount):
        print(f"\n\nConstructing tree {i}")
        
        construct = generateTree(jointCount)
        construct.save("sim_results/" + title + "/" + str(i), saveDir=False)
        lowerBounds.append(construct.totalLengthLowerBound())
        results.append([])
        
        for no, f in enumerate(optimizations):
            print(f"\nTrying loss function {no}")
            direc = "sim_results/" + title + "/" + str(i) + "/" + labels[no]+ "/"
            os.makedirs(direc, exist_ok=True)
            
            optimized, times, losses = f(construct, showSteps=False, parallelize=True, evaluate=True, verbose=False, directory=direc)
            
            if multipleIterations:
                count = 2
                while (losses[0] - losses[-1] > 100):
                    print(f"Trying loss function {no} for the {count}th time")
                    optimized, times, losses = f(construct, showSteps=False, parallelize=True, evaluate=True, verbose=False, directory=None)
                    count += 1
            #print(optimized.detectCollisions(plot=True, includeEnds=False, debug=True))
            results[i].append((times, losses))
        
        with open("sim_results/" + title + "/random_results_chkpt" + str(i) + ".json", "w") as file:
            json.dump(results, file)

    with open("sim_results/" + title + "/random_results.json", "w") as file:
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
        plt.ylabel('Loss')
        plt.title('Loss vs Time: ' + str(title))
        plt.legend()

        plt.grid(True)
        # save to an image file
        plt.savefig("sim_results/" + str(title) + "/example" + str(index) + ".png")
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


testRandomTrees()