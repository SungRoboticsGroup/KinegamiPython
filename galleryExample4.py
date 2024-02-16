# -*- coding: utf-8 -*-
from KinematicChain import *
r = 1
numSides = 6
spacing = 10
numJoints = 6
angleDiff = np.pi/3

chain = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3()))
for i in range(numJoints-1):
    chain.append(RevoluteJoint(numSides, r, np.pi, SE3.Ry(angleDiff)))

chain.show(showLinkPath=False, showLinkPoses=False, showJointPoses=False)
