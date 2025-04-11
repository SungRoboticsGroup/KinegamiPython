from optimizationFunctions import *
from testqtgraph import *
from makeKinematicTree import *
from KinematicTree import *
from KinematicChain import *

import random
from collections import defaultdict, deque
import json
from functools import partial
import threading
import os 
import shutil

jointCount = 5
sparse = False
cubeSize = 20
title = str(jointCount)+" Joint Generalized Gimbal Chains Cube Size " + str(cubeSize)
chainCount = 1
restartFrom = 0
multipleIterations=False

os.makedirs("sim_results/" + title, exist_ok=True)

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

def test():
    optimizations = [partial(squaredOptimize, childParentRatio=1, streamline=True, guarantee=True),
                    partial(linearOptimize, childParentRatio=1, streamline=False, guarantee=True),
                    partial(linearOptimize, childParentRatio=1, streamline=True, guarantee=True),
                    partial(squaredOptimize, childParentRatio=1, streamline=False, guarantee=True)]
    labels = ["Quadratic+Streamline", 
            "Quadratic",
            "Linear+Streamline",
            "Linear"]

    lowerBounds = []
    results = []

    restartDir = "sim_results/" + title + "/" + str(restartFrom)
    if os.path.exists(restartDir):
        shutil.rmtree(restartDir)


    for i in range(restartFrom,chainCount):
        print(f"\n\nConstructing tree {i}")
        construct = generateRandomChain(jointCount)
        construct.save("sim_results/" + title + "/" + str(i), saveDir=False)
        lowerBounds.append(construct.totalLengthLowerBound())
        results.append([])
        for no, f in enumerate(optimizations):
            print(f"\nTrying loss function {no}")
            direc = "sim_results/" + title + "/" + str(i) + "/" + labels[no] + "/"
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
        plt.ylabel('Length')
        #plt.title('Optimization Traversals: ' + str(title))
        plt.legend()

        plt.grid(True)
        # save to an image file
        plt.savefig("sim_results/" + str(title) + "/example" + str(index) + ".png")
        plt.show()

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