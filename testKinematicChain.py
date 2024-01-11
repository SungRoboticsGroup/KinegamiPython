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
KC = KinematicChain(Fingertip(numSides, r, SE3.Ry(np.pi/4), 2, forward=False)) 

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

KC.plot(showSpheres=False)
pattern = KC.tubularOrigamiPattern(numSides)
pattern.makeDXF(show=True)

KC.setJointState(prismaticIndex, maxPrismaticState)

KC.setJointState(revoluteIndex, -np.pi/2)
KC.plot(showSpheres=False)
