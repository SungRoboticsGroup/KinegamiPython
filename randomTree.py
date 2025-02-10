from optimizationFunctions import *
from testqtgraph import *
from makeKinematicTree import *

import random
from collections import defaultdict, deque
import json
from functools import partial
import threading

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

    return makeTubularKinematicTree(tree)

optimizations = [partial(squaredOptimize, childParentRatio=0,streamline=True,guarantee=True),
                partial(squaredOptimize, childParentRatio=0,streamline=True,guarantee=False),
                partial(squaredOptimize, childParentRatio=0,streamline=False,guarantee=False),
                partial(squaredOptimize, childParentRatio=0,streamline=False,guarantee=True),
                partial(linearOptimize, childParentRatio=0, streamline=False,guarantee=False)]
labels = ["Streamline + Guarantee", 
          "Streamline No Guarantee",
          "No Streamline No Guarantee",
          "No Streamline Guarantee",
          "Linear"]

results = []

for i in range(0,2):
    construct = generateTree(6)
    results.append([])
    for f in optimizations:
        optimized, times, losses = f(construct, showSteps=False, parallelize=True, evaluate=True)
        #print(optimized.detectCollisions(plot=True, includeEnds=False, debug=True))
        results[i].append((times, losses))

with open("random_results.json", "w") as file:
    json.dump(results, file)

for result in results:
    plt.figure(figsize=(8, 5)) 
    
    idx = 0
    for x, y in result:
        plt.plot(x, y, marker='o', linestyle='-', label=labels[idx])
        idx += 1

    plt.xlabel('Time')
    plt.ylabel('Loss')
    plt.title('Loss vs Time')
    plt.legend()

    plt.grid(True)
    plt.show()
