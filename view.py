from KinematicTree import *
import numpy as np
from numpy import array

# load and visualize results from randomChain.py

#init_path="/home/daniel/KinegamiPython/Trials Before Experiments/2025.12.05_12.03.23_Joints3_Chains1_Seed42/0.txt"
final_path="/home/daniel/KinegamiPython/Trials Before Experiments/2025.12.05_12.03.23_Joints3_Chains1_Seed42/0/DFS - Inward Longest184.94166469573975_7.tree"

"""
# read initial_path contents to string init_repr
with open(init_path, 'r') as f:
    init_repr = f.read()

init = eval(init_repr)
init.show(block=False)
"""

final = loadKinematicTree(final_path)
#final.show()

configs = np.array([[0, 0, 0, 0, 0, 0, 0, 0], 
           [-0.09206587897246177, 0.0, 0.0, 0.7158213272740908, 0.0, 0.0, 0.04471880714817855, 0.0], 
           [-0.32950833155822457, 0.0, 0.0, -0.2820390044713327, 0.0, 0.0, 0.9688588432292482, 0.0]])
for config in configs:
    final.setConfiguration(config, realJointsOnly=False)
    final.recursivelyRecomputeCollisionCapsules(0)
    print("collisions:", final.detectCollisions(debug=True))
    final.show(block=False)


final.show(block=True)