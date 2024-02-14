# -*- coding: utf-8 -*-
from KinematicChain import *
import inspect
r = 1
numSides = 6

chain = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3.Rz(np.pi/2)))
i = chain.append(RevoluteJoint(numSides, r, np.pi, SE3.Rx(np.pi/2)), relative=False)
j = chain.append(RevoluteJoint(numSides, r, np.pi, SE3.Ry(-np.pi/2)), relative=False)
chain.show(block=False, showLinkPath=False, showJointPoses=False, showLinkPoses=False)
chain.translateJointAlongKinematicAxis(i, 2, propogate=False)
chain.translateJointAlongKinematicAxis(j, 5, propogate=False)
chain.transformJoint(len(chain.Joints)-2, SE3.Trans(-1,-2,0), propogate=False)
chain.show(showLinkPath=False, showJointPoses=False, showLinkPoses=False)
