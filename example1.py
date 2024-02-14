# -*- coding: utf-8 -*-
from KinematicChain import *
r = 1
numSides = 4
chain = KinematicChain(StartTip(numSides, r, Pose=SE3.Trans(3,3,0), length=0.5)) 
# TODO: switch to 3 layers, replace all figures accordingly
prismatic = PrismaticJoint(numSides, r, neutralLength=3, numLayers=3, coneAngle=np.pi/4, Pose= SE3.Trans([5,5,0]))
prismaticIndex = chain.append(prismatic)
revolute = RevoluteJoint(numSides, r, np.pi, SE3.Ry(np.pi/4))
revoluteIndex = chain.append(revolute)
end = EndTip(numSides, r, Pose=SE3.Ry(np.pi/2), length=0.5)
endIndex = chain.append(end)
chain.show(showGlobalFrame=True, block=False)

chain.translateJointAlongKinematicAxis(revoluteIndex, -7)
chain.show(showGlobalFrame=True, block=False)
chain.translateJointAlongKinematicAxis(endIndex, -10, applyToPreviousWaypoint=True)
chain.show(showGlobalFrame=True, block=False)
chain.rotateJointAboutKinematicAxis(revoluteIndex, -np.pi/3)
chain.show(showGlobalFrame=True, block=False)
chain.translateJointAlongKinematicAxis(prismaticIndex, -2, propogate=False)
chain.show(showGlobalFrame=True, block=False)
chain.show(showGlobalFrame=False, showJointPoses=False, showLinkPath=False)
pattern = chain.creasePattern()
pattern.save(dxfName="examplePatterns/example1.dxf")

minPrismaticState, maxPrismaticState = chain.Joints[prismaticIndex].stateRange()
chain.setJointState(prismaticIndex, maxPrismaticState)
chain.setJointState(revoluteIndex, np.pi/2)
chain.show(showGlobalFrame=True)
