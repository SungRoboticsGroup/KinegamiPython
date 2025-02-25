from optimizationFunctions import *
from testqtgraph import *
from makeKinematicTree import *

import random
from collections import defaultdict, deque
import json
from functools import partial
import threading
import os 
import shutil

title = "30 Joint Trees"
jointCount = 30
treeCount = 10
restartFrom = 0
multipleIterations=False

os.makedirs("sim_results/" + title, exist_ok=True)

branchingRatio = 0
def generateTree(nJoints):
    poses = [
        SE3.Rand(xrange=(-100,100),yrange=(-100,100),zrange=(-100,100))
        for _ in range(nJoints)
    ]

    r = 1
    numSides = 4
    neutralLength = 3
    tree = JointSpecificationTree(Waypoint(numSides,r,poses[0]))

    for i in range(1,nJoints):
        parent = np.random.randint(int((i - 1) * (1 - branchingRatio)), i)
        match np.random.randint(1,2):
            case 1:
                tree.addJoint(parent, RevoluteJoint(numSides,r,np.pi,poses[i]))
            case _:
                tree.addJoint(parent, PrismaticJoint(numSides,r,neutralLength,3,np.pi/5,poses[i]))

    try:
        return makeTubularKinematicTree(tree)
    except:
        return generateTree(nJoints)

optimizations = [partial(squaredOptimize, childParentRatio=0,streamline=True,guarantee=True),
                partial(squaredOptimize, childParentRatio=0,streamline=True,guarantee=False),
                partial(squaredOptimize, childParentRatio=0,streamline=False,guarantee=False),
                partial(squaredOptimize, childParentRatio=0,streamline=False,guarantee=False,resetOnFail=False),
                partial(squaredOptimize, childParentRatio=0,streamline=False,guarantee=True),
                partial(linearOptimize, childParentRatio=0, streamline=False,guarantee=False)]
labels = ["Streamline + Guarantee (SG)", 
          "Streamline No Guarantee (SNG)",
          "No Streamline No Guarantee (NSNG)",
          "NSNG, No Reset on Fail",
          "No Streamline Guarantee (NSG)",
          "Linear (L)"]

results = []

restartDir = "sim_results/" + title + "/" + str(restartFrom)
if os.path.exists(restartDir):
    shutil.rmtree(restartDir)

for i in range(restartFrom,treeCount):
    print(f"\n\nConstructing tree {i}")
    construct = generateTree(jointCount)
    construct.save("sim_results/" + title + "/" + str(i), saveDir=False)
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

for result in results:
    plt.figure(figsize=(8, 5)) 
    
    idx = 0
    for x, y in result:
        plt.plot(x, y, marker='o', linestyle='-', label=labels[idx])
        idx += 1

    plt.xlabel('Time')
    plt.ylabel('Loss')
    plt.title('Loss vs Time: ' + str(title))
    plt.legend()

    plt.grid(True)
    plt.show()
