# -*- coding: utf-8 -*-
from KinematicChain import *
import inspect
r = 1
numSides = 6
spacing = 4
numJoints = 3
jointAxisScale = 3

chain = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3()))
for i in range(numJoints-1):
    #chain.show(showSpheres=True, block=False)
    chain.append(RevoluteJoint(numSides, r, np.pi, SE3.Trans(spacing,0,0)))

chain.show(block=False, showLinkPath=False, showJointPoses=False, showLinkPoses=False, jointAxisScale=jointAxisScale)

# Delete waypoints
i = 0
while i < len(chain.Joints):
    if isinstance(chain.Joints[i], Waypoint):
        chain.delete(i)
    else:
        i += 1

chain.show(block=False)

# Move joints to x axis
for i in range(len(chain.Joints)):
    chain.translateJointAlongKinematicAxis(i,-chain.Joints[i].Pose.t[2])
chain.show(showLinkPath=False, showJointPoses=False, showLinkPoses=False, jointAxisScale=jointAxisScale, showAxisGrids=False)
