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
numSides = 6
jointAxisScale=3

# Initialize the chain with a revolute joint at the global origin
chain = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3()))

# Add a prismatic joint, with the same axis of motion as the revolute joint
# (the Z axis): The safe joint placement algorithm will find its position along
# and orientation about this axis
prismaticIndex = chain.append(PrismaticJoint(numSides, r, neutralLength=3, 
                                    numLayers=3, coneAngle=np.pi/4, Pose=SE3()))

chain.show(block=False, showLinkPath=False, showJointPoses=False, 
           showLinkPoses=False, jointAxisScale=jointAxisScale)

# Manual adjustments to make more compact (optional)
chain.translateJointAlongAxisOfMotion(prismaticIndex, -2)
chain.show(showLinkPath=False, showJointPoses=False, showLinkPoses=False, jointAxisScale=jointAxisScale, showAxisGrids=False)
