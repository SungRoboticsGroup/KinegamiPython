# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 19:02:33 2023

@author: Daniel Feshbach
"""
from geometryHelpers import *

Elbow(1, SE3(), np.pi/2, np.pi/3).show(showFrames=True, block=False)
Elbow(2, SE3(10,3,0), np.pi/3, np.pi/2).show(showFrames=True, block=False)
Elbow(1, SE3(10,3,0), 3*np.pi/4, 3*np.pi/4).show(showFrames=True, block=False)
Elbow(1, SE3(2,3,0)@SE3.Rx(np.pi/3), np.pi/2, np.pi/4).show(showFrames=True, block=False)
Elbow(1, SE3.Ry(np.pi/3) * SE3(2,3,0), bendingAngle=3*np.pi/4, 
      rotationalAxisAngle=3*np.pi/4).show(showFrames=True, block=False)
CompoundElbow(1, SE3.Ry(np.pi/3) * SE3(2,3,0), bendingAngle=3*np.pi/4, 
              rotationalAxisAngle=3*np.pi/4, maxAnglePerElbow=np.pi/4).show(showFrames=True, block=blockDefault)
