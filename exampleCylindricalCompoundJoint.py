# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:32:19 2023

@author: Daniel Feshbach
"""
from KinematicChain import *
r = 1
numSides = 6

# chain whose root is a waypoint at the global origin
chain = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3()))
prismaticIndex = chain.append(PrismaticJoint(numSides, r, neutralLength=3, 
                                    numLayers=3, coneAngle=np.pi/4, Pose=SE3()))
chain.show()
chain.translateJointAlongKinematicAxis(prismaticIndex, -2)
chain.show()
