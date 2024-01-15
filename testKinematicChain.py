# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:32:19 2023

@author: Daniel Feshbach
"""
from KinematicChain import *
r = 1
numSides = 6
guarantee = True

# chain whose root is a waypoint at the global origin
KC = KinematicChain(Fingertip(numSides, r, SE3.Ry(np.pi/4), 2, forward=False),
                    maxAnglePerElbow = np.pi/2) 

prismaticIndex = KC.addJointToEnd(PrismaticJoint(numSides, r, neutralLength=3, 
                    numLayers=6, coneAngle=np.pi/4,
                    Pose= SE3.Ry(np.pi/4) @ SE3([1,-3,0])),
                relative=True, fixedPosition=False, fixedOrientation=False, 
                guarantee=guarantee)

minPrismaticState, maxPrismaticState = KC.Joints[prismaticIndex].stateRange()

revoluteIndex = KC.addJointToEnd(RevoluteJoint(numSides, r, np.pi,
                                 SE3.Rx(np.pi/4) @ SE3([3,1,0])),
                relative=True, fixedPosition=False, fixedOrientation=False, 
                guarantee=guarantee)

KC.addJointToEnd(Fingertip(numSides, r, SE3(4,0,0)@SE3.Ry(np.pi/2), 2, forward=True),
                relative=True, fixedPosition=True, fixedOrientation=True, 
                guarantee=False)

KC.show(linkColor='orange', jointColor='blue', showJointPoses=False, 
        showJointAxes=False, jointAxisScale=5, showLinkPath=False,
        surfaceOpacity=1, showSpheres=True, block=False)

KC.translateJointAlongAxis(prismaticIndex, -5)
KC.translateJointAlongAxis(revoluteIndex, -5)
KC.rotateJointAboutAxis(revoluteIndex, -np.pi/4)

KC.show(block=False)


pattern = KC.tubularOrigamiPattern(numSides)
pattern.show(show=True)

KC.setJointState(prismaticIndex, maxPrismaticState)

KC.setJointState(revoluteIndex, -np.pi/2)
KC.show(block=True)