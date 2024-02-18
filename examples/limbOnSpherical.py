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
startRotation = SE3.Rx(np.pi/2)

# Construct spherical compound joint
chain = KinematicChain(RevoluteJoint(numSides, r, np.pi, startRotation))
i = chain.append(RevoluteJoint(numSides, r, np.pi, 
                               startRotation@SE3.Rx(np.pi/2)), relative=False)
j = chain.append(RevoluteJoint(numSides, r, np.pi, 
                               startRotation@SE3.Ry(-np.pi/2)), relative=False)
# Move last joint in spherical compound joint to make it more compact
chain.translateJointAlongAxisOfMotion(j, 5, propogate=False) 

# Attach two more joints on long links, then a tip at the end
k = chain.append(RevoluteJoint(numSides, r, np.pi, 
                            startRotation@SE3.Trans(0,0,20)@SE3.Ry(np.pi/2)), 
                            relative=False, safe=False, fixedPosition=True)
chain.append(RevoluteJoint(numSides, r, np.pi, SE3.Trans(20,0,0)), 
             safe=False, fixedPosition=True)
chain.append(EndTip(numSides, r, SE3.Trans(5,0,0)@SE3.Ry(np.pi/2), 5), 
             safe=False, fixedPosition=True)

# Adjust joint states to pull the limb slightly inwards
chain.setJointState(j, -np.pi/4)
chain.setJointState(k, -np.pi/3)

chain.show(showLinkPath=False, showLinkPoses=False, showJointPoses=False)
