import os
import sys
this_dir = os.path.dirname(__file__)
main_dir = os.path.abspath(os.path.join(this_dir, '../..'))
sys.path.append(main_dir)

# Example 1A
from KinematicTree import *
from testqtgraph import *

r = 1

numSides = 4

jointLength = 4
unextendedRevoluteJointLength = RevoluteJoint(numSides, r, np.pi, SE3()).neutralLength
extensionLength = (jointLength - unextendedRevoluteJointLength)/2

twistJointLength = 2.5

tree = KinematicTree[OrigamiJoint](Waypoint(numSides, r, SE3()))

spine1 = tree.addJoint(0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, SE3.Trans([0,jointLength,jointLength/2]) @ SE3.Rz(np.pi/2)), fixedOrientation=True, fixedPosition=True, safe=False)

topRightInAxis = tree.addJoint(spine1, PrismaticJoint(numSides, r, twistJointLength,3,np.pi/5, SE3.Trans([jointLength,jointLength,0]) @ SE3.Rx(-np.pi/2)), fixedOrientation=True, fixedPosition=True, safe=False)
topRightOrthogonal = tree.addJoint(topRightInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(np.pi/2) @ SE3.Trans([jointLength,0,jointLength])), fixedOrientation=True, fixedPosition=True, safe=False)
topRightFoot = tree.addJoint(topRightOrthogonal, EndTip(numSides, r, SE3.Trans([jointLength,0,0]), 0.5, pathIndex=0), fixedOrientation=True, fixedPosition=True, safe=False)

topLeftInAxis = tree.addJoint(spine1, PrismaticJoint(numSides, r, twistJointLength,3,np.pi/5, SE3.Trans([jointLength,-jointLength,0]) @ SE3.Rx(np.pi/2)), fixedOrientation=True, fixedPosition=True, safe=False)
topLeftOrthogonal = tree.addJoint(topLeftInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(-np.pi/2) @ SE3.Trans([jointLength,0,jointLength])), fixedOrientation=True, fixedPosition=True, safe=False)
topLeftFoot = tree.addJoint(topLeftOrthogonal, EndTip(numSides, r, SE3.Trans([jointLength,0,0]) @ SE3.Ry(np.pi/2), 0.5), fixedOrientation=True, fixedPosition=True, safe=False)


middleRightInAxis = tree.addJoint(0, PrismaticJoint(numSides, r, twistJointLength,3,np.pi/5, SE3.Trans([jointLength,0,jointLength/2]) @ SE3.Ry(np.pi/2)), fixedOrientation=True, fixedPosition=True, safe=False)
middleRightOrthogonal = tree.addJoint(middleRightInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Trans([jointLength,0,jointLength])), fixedOrientation=True, fixedPosition=True, safe=False)
middleRightFoot = tree.addJoint(middleRightOrthogonal, EndTip(numSides, r, SE3.Trans([jointLength,0,0]) @ SE3.Ry(np.pi/2), 0.5), fixedOrientation=True, fixedPosition=True, safe=False)

middleLeftInAxis = tree.addJoint(0, PrismaticJoint(numSides, r, twistJointLength,3,np.pi/5, SE3.Trans([-jointLength,0,jointLength/2]) @ SE3.Ry(-np.pi/2)), fixedOrientation=True, fixedPosition=True, safe=False)
middleLeftOrthogonal = tree.addJoint(middleLeftInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(np.pi) @ SE3.Trans([jointLength,0,jointLength])), fixedOrientation=True, fixedPosition=True, safe=False)
middleLeftFoot = tree.addJoint(middleLeftOrthogonal, EndTip(numSides, r, SE3.Trans([jointLength,0,0]) @ SE3.Ry(np.pi/2), 0.5), fixedOrientation=True, fixedPosition=True, safe=False)

spine2 = tree.addJoint(0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, SE3.Trans([0,-jointLength,jointLength/2]) @ SE3.Rz(-np.pi/2)), fixedOrientation=True, fixedPosition=True, safe=False)

bottomRightInAxis = tree.addJoint(spine2, PrismaticJoint(numSides, r, twistJointLength,3,np.pi/5, SE3.Trans([jointLength,-jointLength,0]) @ SE3.Rx(np.pi/2)), fixedOrientation=True, fixedPosition=True, safe=False)
bottomRightOrthogonal = tree.addJoint(bottomRightInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(-np.pi/2) @ SE3.Trans([jointLength,0,jointLength])), fixedOrientation=True, fixedPosition=True, safe=False)
bottomRightFoot = tree.addJoint(bottomRightOrthogonal, EndTip(numSides, r, SE3.Trans([jointLength,0,0]) @ SE3.Ry(np.pi/2), 0.5), fixedOrientation=True, fixedPosition=True, safe=False)

bottomLeftInAxis = tree.addJoint(spine2, PrismaticJoint(numSides, r, twistJointLength,3,np.pi/5, SE3.Trans([jointLength,jointLength,0]) @ SE3.Rx(-np.pi/2)), fixedOrientation=True, fixedPosition=True, safe=False)
bottomLeftOrthogonal = tree.addJoint(bottomLeftInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(np.pi/2) @ SE3.Trans([jointLength,0,jointLength])), fixedOrientation=True, fixedPosition=True, safe=False)
bottomLeftFoot = tree.addJoint(bottomLeftOrthogonal, EndTip(numSides, r, SE3.Trans([jointLength,0,0]) @ SE3.Ry(np.pi/2), 0.5), fixedOrientation=True, fixedPosition=True, safe=False)

tree.show()

# printedTree = origamiToPrinted(tree, 0.05)

# plotPrintedTree(printedTree, "manualHexapod")