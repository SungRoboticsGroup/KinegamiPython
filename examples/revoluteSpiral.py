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
spacing = 10
numJoints = 6
angleDiff = np.pi/3

# Initialize chain with a revolute joint at the global origin
chain = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3()))

# Add other revolute joints on axes of motion in the XZ plane, 
# successively rotated by pi/3, at exact locations chosen by the safe
# joint placement algorithm
for i in range(numJoints-1):
    chain.append(RevoluteJoint(numSides, r, np.pi, SE3.Ry(angleDiff)))

chain.show(showLinkPath=False, showLinkPoses=False, showJointPoses=False)
