# -*- coding: utf-8 -*-
# Setting up path to be able to import from the parent directory
# https://tutorpython.com/tutorial/python-import-from-parent-directory
import os
import sys
child_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(child_dir, '..'))
sys.path.append(parent_dir)

from KinematicChain import *
r = 1
numSides = 4
chain = KinematicChain(StartTip(numSides, r, Pose=SE3(), length=0.5)) 
# TODO: switch to 3 layers, replace all figures accordingly
prismatic = PrismaticJoint(numSides, r, neutralLength=3, numLayers=3, 
                           coneAngle=np.pi/4, Pose= SE3.Rx(np.pi/2))
prismaticIndex = chain.append(prismatic)
revolute = RevoluteJoint(numSides, r, np.pi, SE3.Ry(np.pi/2))
revoluteIndex = chain.append(revolute)
end = EndTip(numSides, r, Pose=SE3.Trans(4,0,0)@SE3.Ry(np.pi/2), length=0.5)
endIndex = chain.append(end, fixedPosition=True, safe=False)
chain.show(showJointPoses=False, showLinkPath=False)