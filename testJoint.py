# -*- coding: utf-8 -*-
"""
Created on Wed Jan  3 21:34:38 2024

@author: Daniel Feshbach
"""
from Joint import *

numSides = 6
r = 1

testRevolute = RevoluteJoint(numSides, r, totalBendingAngle=np.pi, Pose=SE3.Rx(np.pi/3))
#testRevolute.plot(showSphere = True)
testRevolute.transformPoseBy(SE3(1,1,1) @ SE3.Ry(np.pi/4))
testRevolute.transformStateTo(-np.pi/3)
testRevolute.show(showSphere = True, block=False)

testPrismatic = PrismaticJoint(numSides, r, neutralLength=4, numLayers=4, 
                        coneAngle=np.pi/3, Pose=SE3.Ry(np.pi/4))
#testPrismatic.plot(showSphere=True)
testPrismatic.transformPoseBy(SE3(-2,3,5) @ SE3.Ry(np.pi/4))
testPrismatic.show(showSphere=True, block=False)
testPrismatic.transformStateTo(-1)
testPrismatic.show(showSphere=True, block=False)
testPrismatic.transformStateTo(0.6)
testPrismatic.show(showSphere=True, block=True)
