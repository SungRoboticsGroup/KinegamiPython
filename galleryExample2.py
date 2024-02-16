# -*- coding: utf-8 -*-
from KinematicChain import *
r = 1
numSides = 6
startRotation = SE3.Rx(np.pi/2)

chain = KinematicChain(RevoluteJoint(numSides, r, np.pi, startRotation))
i = chain.append(RevoluteJoint(numSides, r, np.pi, startRotation@SE3.Rx(np.pi/2)), relative=False)
j = chain.append(RevoluteJoint(numSides, r, np.pi, startRotation@SE3.Ry(-np.pi/2)), relative=False)
chain.translateJointAlongKinematicAxis(j, 5, propogate=False)

k = chain.append(RevoluteJoint(numSides, r, np.pi, startRotation@SE3.Trans(0,0,20)@SE3.Ry(np.pi/2)), 
             relative=False, safe=False, fixedPosition=True)

chain.append(RevoluteJoint(numSides, r, np.pi, SE3.Trans(20,0,0)), 
             safe=False, fixedPosition=True)

chain.append(EndTip(numSides, r, SE3.Trans(5,0,0)@SE3.Ry(np.pi/2), 5), 
             safe=False, fixedPosition=True)

chain.setJointState(j, -np.pi/4)
chain.setJointState(k, -np.pi/3)

chain.show(showLinkPath=False, showLinkPoses=False, showJointPoses=False)
