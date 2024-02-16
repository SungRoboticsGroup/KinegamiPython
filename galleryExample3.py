# -*- coding: utf-8 -*-
from KinematicChain import *
r = 1
numSides = 6
spacing = 10
numSphericalJoints = 4

spherical = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3()))
i = spherical.append(RevoluteJoint(numSides, r, np.pi, SE3.Rx(np.pi/2)), relative=False)
j = spherical.append(RevoluteJoint(numSides, r, np.pi, SE3.Ry(-np.pi/2)), relative=False)
spherical.translateJointAlongKinematicAxis(j, 5, propogate=False)

chain = KinematicChain(Waypoint(numSides, r, SE3.Trans(0,0,-5)))
for i in range(numSphericalJoints):
    for Joint in spherical.Joints:
        newJoint = copy.deepcopy(Joint)
        newJoint.transformPoseBy(SE3(0,0,i*spacing))
        chain.append(newJoint, relative=False, safe=False, fixedPosition=True, fixedOrientation=True)

chain.show(showLinkPath=False, showLinkPoses=False, showJointPoses=False)
