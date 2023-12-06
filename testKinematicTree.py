# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 23:44:41 2023

@author: DanielFeshbach

This file exists to test KinematicTree and by extension Joint and dubinsPath
"""
import KinematicTree
from KinematicTree import *

r = 1
numSides = 6

# tree whose root is a waypoint at the global frame
KT = KinematicTree(WayPoint(numSides, r, SE3())) 

j = KT.addJoint(0, PrismaticJoint(numSides, r, neutralLength=3, 
                    numLayers=6, coneAngle=np.pi/4,
                    Pose= SE3.Ry(np.pi/4) @ SE3([1,-3,0])),
                relative=True, fixedPosition=False, fixedOrientation=False)

k = KT.addJoint(j, WayPoint(numSides, r, SE3.Rz(7*np.pi/6) @ SE3.Trans(4,4,4)), 
                relative=True, fixedPosition=False, fixedOrientation=False)

"""
m = KT.addJoint(j, WayPoint(numSides, r, SE3.Ry(np.pi/2) @ SE3.Trans(1,-2,-3)),
                relative=True, fixedPosition=False, fixedOrientation=False)
"""

i = KT.addJoint(0, RevoluteJoint(numSides, r, np.pi,
                                 SE3.Rx(np.pi/4) @ SE3([3,1,0])),
                relative=True, fixedPosition=False, fixedOrientation=False)

"""
n = KT.addJoint(i, WayPoint(numSides, r, SE3.Ry(np.pi/2) @ SE3.Trans(1,-2,-3)),
                relative=True, fixedPosition=False, fixedOrientation=False)
"""


c = KT.addJoint(j, RevoluteJoint(numSides, r, np.pi,
                             SE3.Rx(np.pi/2) @ SE3([0,-2,0])), guarantee=True)

"""


a = KT.addJoint(0, RevoluteJoint(numSides, r, np.pi, 0, 
                             SE3.Ry(np.pi/4) @ SE3([-2,-2,0])), guarantee=True)


# SE3.Rx(np.pi/4)@SE3([1,1,0]) needs turns >theta to reach from the root
b = KT.addJoint(a, RevoluteJoint(numSides, r, np.pi, 0, 
                             SE3.Rx(np.pi/4) @ SE3([1,1,0])), 
                relative=True, guarantee=True)

c = KT.addJoint(0, RevoluteJoint(numSides, r, np.pi, 0, 
                             SE3.Rx(np.pi/4) @ SE3([1,1,0])), 
                relative=True, guarantee=True)
"""
