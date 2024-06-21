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
from testqtgraph import *

r = 0.3     #in inches - currently not working with r != 1, needs debugging
numSides = 4
screwHoleRadius = 0.05 #in inches for 2mm diameter screw with some tolerance
jointLength = 3#0.45/0.245
unextendedRevoluteJointLength = RevoluteJoint(numSides, r, np.pi, SE3()).neutralLength
extensionLength = (jointLength - unextendedRevoluteJointLength)/2
print(extensionLength)

tree = KinematicTree(Waypoint(numSides, r, Pose=SE3()))
palmJoint = ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, SE3.Trans(4,-1,1.5)@SE3.Ry(0)@SE3.Rx(np.pi/5.5) @ SE3.Rz(-np.pi/6))
palm = tree.addJoint(0, palmJoint,
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

thumb0Joint = ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, SE3.Trans(jointLength,0,0)@SE3.Rx(1.95*np.pi/3))
thumb0 = tree.addJoint(palm, thumb0Joint,
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

thumb1 = tree.addJoint(thumb0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
                                           SE3.Trans(5,0,0)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

thumbEnd = tree.addJoint(thumb1, EndTip(numSides, r, 
        SE3.Trans(3.75,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), 1),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)


pointer1 = tree.addJoint(0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
            SE3.Trans(0,0,8)@SE3.Ry(-np.pi/2)@SE3.Rx(np.pi)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

pointer2 = tree.addJoint(pointer1, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
                        SE3.Trans(4,0,0)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

pointerEnd = tree.addJoint(pointer2, EndTip(numSides, r, 
        SE3.Trans(4,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), 1),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

middle0 = tree.addJoint(0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
        SE3.Trans(-2.25,0,4)@SE3.Ry(-np.pi/2)@SE3.Rx(np.pi/2)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

middle1 = tree.addJoint(middle0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
                        SE3.Trans(4,0,0)@SE3.Rx(np.pi/2)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

middle2 = tree.addJoint(middle1, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
                        SE3.Trans(4,0,0)),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

middleEnd = tree.addJoint(middle2, EndTip(numSides, r, 
        SE3.Trans(4,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), 1),
              relative=True, safe=False, 
              fixedPosition=True, fixedOrientation=True)

# print("Branching parameters from root:")
# print(np.round(tree.branchingParametersFrom(0), decimals=3))
# tree.branchModuleFrom(0, "test", 0.045, 0.03, 0.3, 0.25)
# tree.show(block=False, showJointPoses=False)
# print("Branching parameters from root:")
# print(np.round(tree.branchingParametersFrom(0), decimals=3))
# tree.show(block=False, showJointPoses=False)

# print("Straight link lengths:")


# for i in [thumb1, thumbEnd, pointer2, pointerEnd, middle1, middle2, middleEnd]:
#     print(np.round(tree.Links[i].path.tMag, 3))

printedTree = origamiToPrinted(tree, screwHoleRadius)
# # Spherical grasp
printedTree.setJointState(middle0, np.pi/2)
printedTree.setJointState(palm, -np.pi/2.5)
for i in [middle1, middle2, pointer1, pointer2, thumb0, thumb1]:
    printedTree.setJointState(i, np.pi/4)
    printedTree.setJointState(pointer1, np.pi/2.5)
    printedTree.setJointState(thumb1, np.pi/3)
    printedTree.setJointState(middle1, np.pi/2.25)
    printedTree.setJointState(middle2, np.pi/4)
#tree.show(block=False, showJointPoses=False)
plotPrintedTree(printedTree, "origamiToPrintedTest")

# # # Cylindrical grasp
# printedTree.setJointState(middle0, 0)
# printedTree.setJointState(palm, -np.pi/2.25)
# for i in [middle1, middle2, pointer1, pointer2, thumb0, thumb1]:
#     printedTree.setJointState(i, np.pi/3)
#     printedTree.setJointState(thumb1, np.pi/3.75)
#     printedTree.setJointState(thumb0, np.pi/3)
# #tree.show(block=False, showJointPoses=False)
# #plotPrintedTree(printedTree, "origamiToPrintedTest")

# # Pinch grasp
# printedTree.setJointState(middle0, 0)
# printedTree.setJointState(palm, -np.pi/2)
# for i in [middle2, pointer2, thumb1]:
#     printedTree.setJointState(i, 0)
# for i in [middle1, pointer1, thumb0]:
#     printedTree.setJointState(i, np.pi/2)
# printedTree.setJointState(thumb0, np.pi/2.25)
# # # #tree.show(showJointPoses=False)
# plotPrintedTree(printedTree, "origamiToPrintedTest")

# printedTree.export3DKinematicTree("manualHand")
