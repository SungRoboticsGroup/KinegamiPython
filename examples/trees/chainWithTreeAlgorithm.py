# -*- coding: utf-8 -*-
# Setting up path to be able to import from the parent directory
# https://tutorpython.com/tutorial/python-import-from-parent-directory
import os
import sys
this_dir = os.path.dirname(__file__)
main_dir = os.path.abspath(os.path.join(this_dir, '../..'))
sys.path.append(main_dir)
from makeKinematicTree import *

r = 1
numSides = 4

spec = JointSpecificationTree(RevoluteJoint(numSides, r, np.pi, SE3()))
i = spec.addJoint(0, PrismaticJoint(numSides,r,3,3,np.pi/5,SE3(4,0,0)@SE3.Ry(np.pi/2)), relative=True)
i = spec.addJoint(i, RevoluteJoint(numSides,r,np.pi,SE3.Trans(0,0,4)@SE3.Ry(2*np.pi/3)), relative=True)
tree = makeTubularKinematicTree(spec, plotSteps=True)
tree.show()