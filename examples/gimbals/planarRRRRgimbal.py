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
numSides = 4
spacing = 10
numJoints = 4

gg = KinematicChain(RevoluteJoint(numSides, r, 3*np.pi/2, SE3()), gimbal=True)
for i in range(1, numJoints):
        gg.appendGeneralizedGimbal(RevoluteJoint(numSides, r, 3*np.pi/2, SE3(0, spacing*i, 0)), relative=False)

gg.show(showSpheres=False, showScaleBar=False, showJointPoses=False, 
        showLinkPoses=False, showLinkPath=False, showAxisGrids=False, jointAxisScale=30)
