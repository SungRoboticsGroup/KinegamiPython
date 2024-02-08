# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:32:19 2023

@author: Daniel Feshbach
"""
from KinematicChain import *
r = 1
numSides = 6

# chain whose root is a waypoint at the global origin
chain = KinematicChain(StartFingertip(numSides, r, Pose=SE3(), length=0.5)) 
#chain.show(block=False, showJointSurface=False, showLinkSurface=False, showSpheres=True)

prismaticIndex = chain.append(PrismaticJoint(numSides, r, neutralLength=3, numLayers=6, 
                        coneAngle=np.pi/4, Pose= SE3() ) )
#chain.show(block=False, showJointSurface=False, showLinkSurface=False, showSpheres=True)

revoluteIndex = chain.append(RevoluteJoint(numSides, r, np.pi, 
                                           SE3.Ry(np.pi/4)))
#chain.show(block=False, showJointSurface=False, showLinkSurface=False, showSpheres=True)
end = chain.append(EndFingertip(numSides, r, Pose=SE3.Ry(np.pi/2), length=0.5))
#chain.show(block=False, showJointSurface=False, showLinkSurface=False, showSpheres=True)
chain.show(block=False)
#chain.show(block=True, showLinkPath=False, showJointPoses=False, showLinkPoses=False)

chain.translateJointAlongKinematicAxis(prismaticIndex, -5)
chain.show(block=False)
chain.translateJointAlongKinematicAxis(revoluteIndex, -7)
chain.show(block=False)
chain.translateJointAlongKinematicAxis(end, -10, applyToPreviousWaypoint=True)
chain.show(block=False)
chain.rotateJointAboutKinematicAxis(revoluteIndex, -np.pi/2)
chain.show(block=False)

pattern = chain.creasePattern()
pattern.show(show=True, block=False)

minPrismaticState, maxPrismaticState = chain.Joints[prismaticIndex].stateRange()
chain.setJointState(prismaticIndex, maxPrismaticState)

chain.setJointState(revoluteIndex, -np.pi/2)
chain.show(block=True)
