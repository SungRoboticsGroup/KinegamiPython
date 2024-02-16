# -*- coding: utf-8 -*-
from KinematicChain import *
r = 1
numSides = 6
spacing = 10
numJoints = 5
angleDiff = np.pi/4

chain = KinematicChain(PrismaticJoint(numSides, r, 3, 5, np.pi/4, SE3()))
for i in range(numJoints-1):
    chain.append(PrismaticJoint(numSides, r, 3, 5, np.pi/4, SE3.Ry(angleDiff)), safe=False)

chain.show(showLinkPath=False, showLinkPoses=False, showJointPoses=False)
