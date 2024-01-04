# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:32:19 2023

@author: Daniel Feshbach
"""
from KinematicChain import *

r = 1
numSides = 6

# chain whose root is a waypoint at the global origin
KC = KinematicChain(WayPoint(numSides, r, SE3())) 

KC.addJointToEnd(PrismaticJoint(numSides, r, neutralLength=3, 
                    numLayers=6, coneAngle=np.pi/4,
                    Pose= SE3.Ry(np.pi/4) @ SE3([1,-3,0])),
                relative=True, fixedPosition=False, fixedOrientation=False)

KC.addJointToEnd(RevoluteJoint(numSides, r, np.pi,
                                 SE3.Rx(np.pi/4) @ SE3([3,1,0])),
                relative=True, fixedPosition=False, fixedOrientation=False)

KC.addJointToEnd(WayPoint(numSides, r, SE3([1,1,1])),
                relative=True, fixedPosition=False, fixedOrientation=False) 

KC.addJointToEnd(WayPoint(numSides, r, SE3([3,1,1])),
                relative=True, fixedPosition=True, fixedOrientation=True)

KC.plot()
pattern = KC.tubularOrigamiPattern(numSides, splitLongElbowsInto=2)
pattern.makeDXF(show=True)
