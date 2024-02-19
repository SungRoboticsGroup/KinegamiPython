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
numJoints = 5
angleDiff = np.pi/4

# Initialize chain with a revolute joint at the global origin
chain = KinematicChain(PrismaticJoint(numSides, r, 3, 5, np.pi/4, SE3()))

# Add other revolute joints on axes of motion in the XZ plane, 
# successively rotated by pi/4, at exact locations chosen by the compact
# joint placement algorithm
for i in range(numJoints-1):
    chain.append(PrismaticJoint(numSides, r, 3, 5, np.pi/4, SE3.Ry(angleDiff)), 
                 safe=False)

chain.show(showLinkPath=False, showLinkPoses=False, showJointPoses=False)
