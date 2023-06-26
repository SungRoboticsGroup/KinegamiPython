# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 23:44:41 2023

@author: dfesh

This file exists to test KinematicTree and by extension Joint and dubinsPath
"""
import KinematicTree
from KinematicTree import *

r = 1

# tree whose root is a waypoint at the global frame
KT = KinematicTree(WayPoint(r, SE3())) 


# SE3.Rx(np.pi/4)@SE3([1,1,0]) needs turns >theta to reach from the root
pose1 = SE3.Rx(np.pi/4) @ SE3([3,1,0])
i = KT.addJoint(0, RevoluteJoint(r, pose1, 1))

pose2 = SE3.Ry(np.pi/4) @ SE3([1,-3,0])
j = KT.addJoint(0, PrismaticJoint(r, pose2, 1))

pose3 = SE3.Rz(7*np.pi/6) @ SE3.Trans(4,4,4) @ pose2
k = KT.addJoint(j, WayPoint(r, pose3))

pose4 = SE3.Ry(np.pi/2) @ SE3.Trans(1,-2,-3) @ pose2
m = KT.addJoint(j, WayPoint(r, pose4))

KT.plot()

