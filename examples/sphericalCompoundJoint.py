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

# A spherical compound joint is 3 revolute joints with linearly independent
# axes of motion intersecting at a single point (the rotational origin of
# the compound joint)

# Initialize chain with first joint at the origin, facing towards y
chain = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3.Rz(np.pi/2)))
# Append the other two joints with axes of motion along the other two global axes
# and exact poses found by the safe joint placement algorithm
i = chain.append(RevoluteJoint(numSides, r, np.pi, SE3.Rx(np.pi/2)), relative=False)
j = chain.append(RevoluteJoint(numSides, r, np.pi, SE3.Ry(-np.pi/2)), relative=False)
chain.show(block=False, showLinkPath=False, showJointPoses=False, 
           showLinkPoses=False, showGlobalFrame=True)

# Move the joints along their axes of motion to make the compound structure
# more compact
chain.translateJointAlongAxisOfMotion(i, 2, propogate=False)
chain.translateJointAlongAxisOfMotion(j, 5, propogate=False)
# Move a waypoint to make that link more compact
chain.transformJoint(len(chain.Joints)-2, SE3.Trans(-1,-2,0), propogate=False, relative=False)
chain.show(showLinkPath=False, showJointPoses=False, showLinkPoses=False)
