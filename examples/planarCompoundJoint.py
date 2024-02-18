# -*- coding: utf-8 -*-
# Setting up path to be able to import from the parent directory
# https://tutorpython.com/tutorial/python-import-from-parent-directory
import os
import sys
child_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(child_dir, '..'))
sys.path.append(parent_dir)

from KinematicChain import *
r = 1
numSides = 6
spacing = 4
numJoints = 3 
jointAxisScale = 3

# Initialize chain with a revolute joint
planar = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3()))

# Append more parallel revolute joints at distance spacing forward (in x axis)
for i in range(numJoints-1):
    planar.append(RevoluteJoint(numSides, r, np.pi, SE3.Trans(spacing,0,0)),
                  safe=False, fixedPosition=True, fixedOrientation=True)

# Plot the chain
planar.show(showLinkPath=False, showJointPoses=False, showLinkPoses=False, 
            jointAxisScale=jointAxisScale, showAxisGrids=False)




"""
# Alternatively, here is a convoluted way to construct the same example 
# demonstrating some additional methods:

# Initialize chain with a revolute joint
chain = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3()))

# Use safe joint placement algorithm to append numJoints-1 additional revolute
# joints parallel to the first with axes 
for i in range(numJoints-1):
    chain.append(RevoluteJoint(numSides, r, np.pi, SE3.Trans(spacing,0,0)))

chain.show(block=False, showLinkPath=False, showJointPoses=False, showLinkPoses=False, jointAxisScale=jointAxisScale)

# Delete waypoints
i = 0
while i < len(chain.Joints):
    if isinstance(chain.Joints[i], Waypoint):
        chain.delete(i)
    else:
        i += 1

chain.show(block=False)

# Move joints to x axis
for i in range(len(chain.Joints)):
    chain.translateJointAlongAxisOfMotion(i,-chain.Joints[i].Pose.t[2])
chain.show(showLinkPath=False, showJointPoses=False, showLinkPoses=False, 
           jointAxisScale=jointAxisScale, showAxisGrids=False)
"""