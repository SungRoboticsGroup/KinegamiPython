# -*- coding: utf-8 -*-
# Setting up path to be able to import from the parent directory
# https://tutorpython.com/tutorial/python-import-from-parent-directory
import os
import sys
this_dir = os.path.dirname(__file__)
main_dir = os.path.abspath(os.path.join(this_dir, '../..'))
sys.path.append(main_dir)

from KinematicChain import *
r = 1
numSides = 6
numJoints = 3

# Construct a spherical compound joint
chain = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3()), gimbal=True)
for i in range(1, numJoints):
    chain.appendGeneralizedGimbal(RevoluteJoint(numSides, r, np.pi, SE3()), relative=False)

chain.show(showLinkPath=False, showLinkPoses=False, showJointPoses=False, 
           jointAxisScale=50, showSpheres=False, showScaleBar=False)
