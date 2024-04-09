# -*- coding: utf-8 -*-
# Setting up path to be able to import from the parent directory
# https://tutorpython.com/tutorial/python-import-from-parent-directory
import os
import sys
child_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(child_dir, '..'))
sys.path.append(parent_dir)

# Example 1A
from KinematicTree import *
r = 1
numSides = 4
tree = KinematicTree(Waypoint(numSides, r, Pose=SE3()))
palm = tree.addJoint(0, RevoluteJoint(numSides, r, np.pi, 
                        SE3.Trans(2.25,0,3)@SE3.Ry(-np.pi/4)@SE3.Rx(np.pi/4)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

thumb0 = tree.addJoint(palm, RevoluteJoint(numSides, r, np.pi,
                                           SE3.Trans(1.5,0,0)@SE3.Rx(np.pi/2)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

thumb1 = tree.addJoint(thumb0, RevoluteJoint(numSides, r, np.pi,
                                           SE3.Trans(4,0,0)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

thumbEnd = tree.addJoint(thumb1, EndTip(numSides, r, 
        SE3.Trans(3,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), 1),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)


pointer1 = tree.addJoint(0, RevoluteJoint(numSides, r, np.pi, 
                        SE3.Trans(0,0,8)@SE3.Ry(-np.pi/2)@SE3.Rx(np.pi)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

pointer2 = tree.addJoint(pointer1, RevoluteJoint(numSides, r, np.pi, 
                        SE3.Trans(4,0,0)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

pointerEnd = tree.addJoint(pointer2, EndTip(numSides, r, 
        SE3.Trans(4,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), 1),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

middle0 = tree.addJoint(0, RevoluteJoint(numSides, r, np.pi, 
                        SE3.Trans(-2.25,0,4)@SE3.Ry(-np.pi/2)@SE3.Rx(np.pi/2)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

middle1 = tree.addJoint(middle0, RevoluteJoint(numSides, r, np.pi, 
                        SE3.Trans(4,0,0)@SE3.Rx(np.pi/2)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

middle2 = tree.addJoint(middle1, RevoluteJoint(numSides, r, np.pi, 
                        SE3.Trans(4,0,0)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

middleEnd = tree.addJoint(middle2, EndTip(numSides, r, 
        SE3.Trans(4,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), 1),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

tree.show(block=False)

# Spherical grasp
tree.setJointState(middle0, np.pi/2)
tree.setJointState(palm, -np.pi/2)
for i in [middle1, middle2, pointer1, pointer2, thumb0, thumb1]:
    tree.setJointState(i, np.pi/3)
tree.show(block=False)

# Cylindrical grasp
tree.setJointState(middle0, 0)
tree.setJointState(palm, -np.pi/2)
for i in [middle1, middle2, pointer1, pointer2, thumb0, thumb1]:
    tree.setJointState(i, np.pi/3)
tree.show(block=False)

# Pinch grasp
tree.setJointState(middle0, 0)
tree.setJointState(palm, -np.pi/2)
for i in [middle2, pointer2, thumb1]:
    tree.setJointState(i, 0)
for i in [middle1, pointer1, thumb0]:
    tree.setJointState(i, np.pi/2)
tree.show()
