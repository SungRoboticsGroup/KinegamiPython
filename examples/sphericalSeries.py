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
numSphericalJoints = 4

# Construct a spherical compound joint
spherical = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3()))
i = spherical.append(RevoluteJoint(numSides, r, np.pi, SE3.Rx(np.pi/2)), relative=False)
j = spherical.append(RevoluteJoint(numSides, r, np.pi, SE3.Ry(-np.pi/2)), relative=False)
spherical.translateJointAlongAxisOfMotion(j, 5, propogate=False)

# Initialize the chain with a waypoint
chain = KinematicChain(Waypoint(numSides, r, SE3.Trans(0,0,-5)))
# Repeatedly add compound spherical joints by copying the joints (and waypoints)
# from spherical above, translating them up to the appropriate level, and
# then adding them to the chain at that exact pose
for i in range(numSphericalJoints):
    for Joint in spherical.Joints:
        newJoint = copy.deepcopy(Joint)
        newJoint.transformPoseIntoFrame(SE3(0,0,i*spacing))
        chain.append(newJoint, relative=False, safe=False, fixedPosition=True, 
                     fixedOrientation=True)

chain.show(showLinkPath=False, showLinkPoses=False, showJointPoses=False)
