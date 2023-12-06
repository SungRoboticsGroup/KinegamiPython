# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:32:19 2023

@author: Daniel Feshbach
"""
from KinematicChain import *

r = 1
numSides = 6

# tree whose root is a waypoint at the global frame
KC = KinematicChain(WayPoint(numSides, r, SE3())) 

KC.addJoint(PrismaticJoint(numSides, r, neutralLength=3, 
                    numLayers=6, coneAngle=np.pi/4,
                    Pose= SE3.Ry(np.pi/4) @ SE3([1,-3,0])),
                relative=True, fixedPosition=False, fixedOrientation=False)

KC.addJoint(RevoluteJoint(numSides, r, np.pi,
                                 SE3.Rx(np.pi/4) @ SE3([3,1,0])),
                relative=True, fixedPosition=False, fixedOrientation=False)

pattern = KC.tubularOrigamiPattern(numSides, splitLongElbowsInto=2)
pattern.makeDXF(show=True)
