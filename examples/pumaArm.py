# -*- coding: utf-8 -*-
# Setting up path to be able to import from the parent directory
# https://tutorpython.com/tutorial/python-import-from-parent-directory
import os
import sys
child_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(child_dir, '..'))
sys.path.append(parent_dir)

from KinematicChain import *
from numpy import sin, cos
r = 100 #units: mm
numSides = 6

"""
An example PUMA arm, based on Denavit-Hartenberg parameters from:
C. S. G. Lee and M. Ziegler, 
"Geometric Approach in Solving Inverse Kinematics of PUMA Robots," 
in IEEE Transactions on Aerospace and Electronic Systems, 
vol. AES-20, no. 6, pp. 695-706, Nov. 1984, doi: 10.1109/TAES.1984.310452.
"""

# Function to convert Denavit-Hartenberg parameters to relative pose matrices
def relativePoseFromDHOriginal(theta, alpha, a, d):
    return SE3([[cos(theta), -cos(alpha)*sin(theta), sin(alpha)*sin(theta), a*cos(theta)],
                [sin(theta), cos(alpha)*cos(theta), -sin(alpha)*cos(theta), a*sin(theta)],
                [0,          sin(alpha),             cos(alpha),            d],
                [0,          0,                      0,                     1]])

# Construct pose matrices from the Denavit-Hartenberg parameters for a PUMA arm
# from Lee and Ziegler 1984
RelativePoses = [relativePoseFromDHOriginal(np.pi/2, -np.pi/2, 0, 0),
         relativePoseFromDHOriginal(0, 0, 431.8, 149.09),
         relativePoseFromDHOriginal(np.pi/2, np.pi/2, -20.32, 0),
         relativePoseFromDHOriginal(0, -np.pi/2, 0, 433.07),
         relativePoseFromDHOriginal(0, np.pi/2, 0, 0),
         relativePoseFromDHOriginal(0, 0, 0, 56.25)]

GlobalPoses = [RelativePoses[0]]
for Pose in RelativePoses[1:]:
    GlobalPoses.append(GlobalPoses[-1] @ Pose)



# Initialize the chain with a base waypoint
chain = KinematicChain(Waypoint(numSides, r, SE3.Trans(0,0,-660.4)@SE3.Rz(np.pi)))

# The first joint axis is horizontal, 660.4 mm above the start base
chain.append(RevoluteJoint(numSides, r, np.pi, SE3()), safe=False)
# Add the remaining revolute joints using the compact joint placement algorithm
for Pose in GlobalPoses[:-1]:
    chain.append(RevoluteJoint(numSides, r, np.pi, Pose), safe=False, relative=False)
# The last DH parameters (and thus pose matrix) is for the end effector
chain.append(EndTip(numSides, r, GlobalPoses[-1], 50), safe=False, relative=False)

# Plot the resulting chain
chain.show(block=False, showLinkPath=False, showJointPoses=False, showLinkPoses=False, showAxisGrids=False)


# Adjust the resulting chain to shorten links
chain.translateJointAlongAxisOfMotion(1, -200, propogate=False)
chain.translateJointAlongAxisOfMotion(2, -400, propogate=False)
chain.translateJointAlongAxisOfMotion(3, 1400, propogate=False)
chain.translateJointAlongAxisOfMotion(4, -900, propogate=False)
chain.translateJointAlongAxisOfMotion(5, -2100, propogate=False)
chain.translateJointAlongAxisOfMotion(6, -1500, propogate=True)

# Plot the chain structure
chain.show(block=False, showLinkPath=False, showJointPoses=False, showLinkPoses=False, showAxisGrids=False)
# Plot the crease pattern
chain.creasePattern().show()

