from optimizationFunctions import *
from testqtgraph import *
from makeKinematicTree import *

import random
from collections import defaultdict, deque
import json
from functools import partial

def generateTree(nJoints):
    poses = [
        SE3(
            random.uniform(-100, 100),  # Random x-coordinate
            random.uniform(-100, 100),  # Random y-coordinate
            random.uniform(-100, 100),  # Random z-coordinate
        ) * SE3.RPY(
            random.uniform(0, 360),  # Random roll in degrees
            random.uniform(0, 360),  # Random pitch in degrees
            random.uniform(0, 360),  # Random yaw in degrees
            unit='deg'
        )
        for _ in range(nJoints)
    ]

    r = 1
    numSides = 4
    neutralLength = 3
    tree = JointSpecificationTree(Waypoint(numSides,r,poses[0]))

    for i in range(1,nJoints):
        match np.random.randint(1,3):
            case 1:
                tree.addJoint(np.random.randint(0, i), RevoluteJoint(numSides,r,np.pi,poses[i]))
            case 2:
                tree.addJoint(np.random.randint(0, i), PrismaticJoint(numSides,r,neutralLength,3,np.pi/5,poses[i]))
            case _:
                tree.addJoint(np.random.randint(0, i), Waypoint(numSides,r,poses[i]))

    return makeTubularKinematicTree(tree)

optimizations = [partial(squaredOptimize, childParentRatio=0,streamline=True,guarantee=True),
                partial(squaredOptimize, childParentRatio=0,streamline=True,guarantee=False),
                partial(squaredOptimize, childParentRatio=0,streamline=False,guarantee=False),
                partial(squaredOptimize, childParentRatio=0,streamline=False,guarantee=True)]

results = []
for i in range(0,3):
    construct = generateTree(12)
    results.append([])
    for f in optimizations:
        optimized, times, losses = f(construct, showSteps=False, parallelize=True, evaluate=True)
        #print(optimized.detectCollisions(plot=True, includeEnds=False, debug=True))
        results[i].append((times, losses))

with open("random_results.json", "w") as file:
    json.dump(results, file)

for result in results:
    plt.figure(figsize=(8, 5)) 

    for x, y in result:
        plt.plot(x, y, marker='o', linestyle='-', label='Loss over Time')

    plt.xlabel('Time')
    plt.ylabel('Loss')
    plt.title('Loss vs Time')
    plt.legend()

    plt.grid(True)
    plt.show()
