# -*- coding: utf-8 -*-
# Setting up path to be able to import from the parent directory
# https://tutorpython.com/tutorial/python-import-from-parent-directory
import os
import sys
this_dir = os.path.dirname(__file__)
main_dir = os.path.abspath(os.path.join(this_dir, '../..'))
sys.path.append(main_dir)

# Example 1A
from KinematicTree import *
from PrintedJoint import *
from testqtgraph import *

os.chdir(main_dir)

r = 1
multiplier = 1.5
screwRadius = 0.05*multiplier

tree = KinematicTree(PrintedWaypoint(r, SE3(), screwRadius))
palmJoint = PrintedOrthogonalRevoluteJoint(r, -np.pi/2, np.pi/2, SE3.Trans(2.5*multiplier,0,3*multiplier)@SE3.Ry(-np.pi/4)@SE3.Rx(np.pi/4), screwRadius)
palm = tree.addJoint(0, palmJoint,
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

thumb0Joint = PrintedOrthogonalRevoluteJoint(r, -np.pi/2, np.pi/2, SE3.Trans(palmJoint.neutralLength,0,0)@SE3.Rx(np.pi/2), screwRadius)
thumb0 = tree.addJoint(palm, thumb0Joint,
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

thumb1 = tree.addJoint(thumb0, PrintedOrthogonalRevoluteJoint(r, -np.pi/2, np.pi/2,
                                           SE3.Trans(4*multiplier,0,0), screwRadius),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

thumbEnd = tree.addJoint(thumb1, PrintedTip(r, 
        SE3.Trans(3*multiplier,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), screwRadius),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)


pointer1 = tree.addJoint(0, PrintedOrthogonalRevoluteJoint(r, -np.pi/2, np.pi/2,
            SE3.Trans(0,0,8*multiplier)@SE3.Ry(-np.pi/2)@SE3.Rx(np.pi), screwRadius),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

pointer2 = tree.addJoint(pointer1, PrintedOrthogonalRevoluteJoint(r, -np.pi/2, np.pi/2,
                        SE3.Trans(4*multiplier,0,0), screwRadius),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

pointerEnd = tree.addJoint(pointer2, PrintedTip(r, 
        SE3.Trans(4*multiplier,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), screwRadius),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

middle0 = tree.addJoint(0, PrintedOrthogonalRevoluteJoint(r, -np.pi/2, np.pi/2,
        SE3.Trans(-2.25*multiplier,0,4*multiplier)@SE3.Ry(-np.pi/2)@SE3.Rx(np.pi/2), screwRadius),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

middle1 = tree.addJoint(middle0, PrintedOrthogonalRevoluteJoint(r, -np.pi/2, np.pi/2,
                        SE3.Trans(4*multiplier,0,0)@SE3.Rx(np.pi/2), screwRadius),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

middle2 = tree.addJoint(middle1, PrintedOrthogonalRevoluteJoint(r, -np.pi/2, np.pi/2,
                        SE3.Trans(4*multiplier,0,0), screwRadius),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

middleEnd = tree.addJoint(middle2, PrintedTip(r, 
        SE3.Trans(4*multiplier,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), screwRadius),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

# Spherical grasp
tree.setJointState(middle0, np.pi/2)
tree.setJointState(palm, -np.pi/2)
for i in [middle1, middle2, pointer1, pointer2, thumb0, thumb1]:
    tree.setJointState(i, np.pi/3)

# Cylindrical grasp
# tree.setJointState(middle0, 0)
# tree.setJointState(palm, -np.pi/2)
# for i in [middle1, middle2, pointer1, pointer2, thumb0, thumb1]:
#     tree.setJointState(i, np.pi/3)

# # Pinch grasp
# tree.setJointState(middle0, 0)
# tree.setJointState(palm, -np.pi/2)
# for i in [middle2, pointer2, thumb1]:
#     tree.setJointState(i, 0)
# for i in [middle1, pointer1, thumb0]:
#     tree.setJointState(i, np.pi/2)

#tree.export3DKinematicTree("manualHand")
plotPrintedTree(tree, "manualHand")