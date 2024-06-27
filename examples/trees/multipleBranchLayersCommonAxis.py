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

realJoints = [0]
spec = JointSpecificationTree(RevoluteJoint(numSides, r, np.pi, SE3()))
realJoints.append(spec.addJoint(0, PrismaticJoint(numSides,r,3,3,np.pi/5,SE3()), relative=True))
realJoints.append(spec.addJoint(realJoints[1], RevoluteJoint(numSides,r,np.pi,SE3()), relative=True))
realJoints.append(spec.addJoint(realJoints[1], RevoluteJoint(numSides,r,np.pi,SE3()), relative=True))
realJoints.append(spec.addJoint(realJoints[2], RevoluteJoint(numSides,r,np.pi,SE3()), relative=True))
realJoints.append(spec.addJoint(realJoints[2], PrismaticJoint(numSides,r,3,3,np.pi/5,SE3()), relative=True))
tree = makeTubularKinematicTree(spec, plotSteps=False)
tree.show()