# -*- coding: utf-8 -*-
# Setting up path to be able to import from the parent directory
# https://tutorpython.com/tutorial/python-import-from-parent-directory
import os
import sys
this_dir = os.path.dirname(__file__)
main_dir = os.path.abspath(os.path.join(this_dir, '../..'))
sys.path.append(main_dir)
from makeKinematicTree import *
from testqtgraph import *

r = 1
numSides = 4

realJoints = [0]
spec = JointSpecificationTree(RevoluteJoint(numSides, r, np.pi, SE3()))
i = spec.addJoint(0, PrismaticJoint(numSides,r,3,3,np.pi/5,SE3.Ry(np.pi/3)), relative=False)
j = spec.addJoint(i, RevoluteJoint(numSides,r,np.pi,SE3()), relative=False)
spec.addJoint(i, RevoluteJoint(numSides,r,np.pi,SE3(10,0,0)@SE3.Ry(-np.pi/6)), relative=False)
spec.addJoint(j, RevoluteJoint(numSides,r,np.pi,SE3.Ry(-np.pi/4)), relative=False)
spec.addJoint(j, PrismaticJoint(numSides,r,3,3,np.pi/5,SE3()), relative=False)
tree = makeTubularKinematicTree(spec, plotSteps=False, orientUp=True)
tree.show()

printed = origamiToPrinted(tree, 0.05)
plotPrintedTree(printed, "multipleBranchLayersPrinted")