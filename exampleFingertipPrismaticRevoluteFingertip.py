# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:32:19 2023

@author: Daniel Feshbach
"""
from KinematicChain import *
r = 1
numSides = 6
chain = KinematicChain(StartFingertip(numSides, r, Pose=SE3(), length=0.5)) 
# TODO: switch to 3 layers, replace all figures accordingly
prismatic = PrismaticJoint(numSides, r, neutralLength=3, numLayers=6, coneAngle=np.pi/4, Pose= SE3.Trans([5,5,0]))
prismaticIndex = chain.append(prismatic)
revolute = RevoluteJoint(numSides, r, np.pi, SE3.Ry(np.pi/4))
revoluteIndex = chain.append(revolute)
end = EndFingertip(numSides, r, Pose=SE3.Ry(np.pi/2), length=0.5)
endIndex = chain.append(end)
chain.show()

chain.translateJointAlongKinematicAxis(revoluteIndex, -7)
chain.translateJointAlongKinematicAxis(endIndex, -10, applyToPreviousWaypoint=True)
chain.show()
chain.translateJointAlongKinematicAxis(endIndex, -5)
chain.rotateJointAboutKinematicAxis(revoluteIndex, -np.pi/3)
chain.show()
chain.translateJointAlongKinematicAxis(prismaticIndex, -5, propogate=False)
chain.show()

chain.creasePattern().show()

minPrismaticState, maxPrismaticState = chain.Joints[prismaticIndex].stateRange()
chain.setJointState(prismaticIndex, maxPrismaticState)
chain.setJointState(prismaticIndex, maxPrismaticState+1)
chain.setJointState(revoluteIndex, np.pi/2)
chain.show(block=True)
