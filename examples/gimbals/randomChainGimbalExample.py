# -*- coding: utf-8 -*-
# Setting up path to be able to import from the parent directory
# https://tutorpython.com/tutorial/python-import-from-parent-directory
import os
import sys
this_dir = os.path.dirname(__file__)
main_dir = os.path.abspath(os.path.join(this_dir, '../..'))
sys.path.append(main_dir)

from optimizationFunctions import *
from testqtgraph import *
from makeKinematicTree import *
from KinematicTree import *
from KinematicChain import *

import random
from collections import defaultdict, deque
import json
from functools import partial
import threading
import os 
import shutil

def generateRandomChain(nJoints, cubeSize):
    bounds = (-cubeSize/2, cubeSize/2)
    poses = [ SE3.Rand(xrange=bounds, yrange=bounds, zrange=bounds)
                for _ in range(nJoints) ]

    r = 1
    numSides = 4
    neutralLength = 3

    root = RevoluteJoint(numSides,r,np.pi,poses[0]) if np.random.rand() > 0.5 \
            else PrismaticJoint(numSides,r,neutralLength,3,np.pi/5,poses[0])

    chain = KinematicChain(root, gimbal=True)
    for i in range(1,nJoints):
        newJoint = RevoluteJoint(numSides,r,np.pi,poses[i]) if np.random.rand() > 0.5 \
            else PrismaticJoint(numSides,r,neutralLength,3,np.pi/5,poses[i])
        chain.appendGeneralizedGimbal(newJoint)

    return chain

chain = generateRandomChain(nJoints=6, cubeSize=100)
chain.show(showSpheres=False, showScaleBar=False, showJointPoses=False, 
        showLinkPoses=False, showLinkPath=False, showAxisGrids=False, jointAxisScale=100)