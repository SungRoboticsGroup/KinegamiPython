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
spacing = 50
numSphericalJoints = 3

# Construct a spherical compound joint
chain = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3()), gimbal=True)
chain.appendGeneralizedGimbal(RevoluteJoint(numSides, r, np.pi, SE3.Rx(np.pi/2)), relative=False)
chain.appendGeneralizedGimbal(RevoluteJoint(numSides, r, np.pi, SE3.Ry(-np.pi/2)), relative=False)

for i in range(1, numSphericalJoints):
    T = SE3(spacing*i, 0, 0)
    chain.appendGeneralizedGimbal(RevoluteJoint(numSides, r, np.pi, T), relative=False)
    chain.appendGeneralizedGimbal(RevoluteJoint(numSides, r, np.pi, SE3.Rx(np.pi/2)@T), relative=False)
    chain.appendGeneralizedGimbal(RevoluteJoint(numSides, r, np.pi, SE3.Ry(-np.pi/2)@T), relative=False)

chain.show(showLinkPath=False, showLinkPoses=False, showJointPoses=False, 
           jointAxisScale=50, showSpheres=False, showScaleBar=False)
