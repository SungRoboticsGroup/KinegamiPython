import os
import sys
this_dir = os.path.dirname(__file__)
main_dir = os.path.abspath(os.path.join(this_dir, '../..'))
sys.path.append(main_dir)

# Example 1A
from KinematicTree import *
from testqtgraph import *
from makeKinematicTree import *

r = 1

numSides = 4

jointLength = 4
unextendedRevoluteJointLength = RevoluteJoint(numSides, r, np.pi, SE3()).neutralLength
extensionLength = (jointLength - unextendedRevoluteJointLength)/2

twistJointLength = 2.5

tree = JointSpecificationTree(Waypoint(numSides, r, SE3()))

spine1 = tree.addJoint(0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, SE3.Trans([0,jointLength,jointLength/2]) @ SE3.Rz(np.pi/2)))

topRightInAxis = tree.addJoint(spine1, PrismaticJoint(numSides, r, twistJointLength,3,np.pi/5, SE3.Trans([jointLength,jointLength,0]) @ SE3.Rx(-np.pi/2)))
topRightOrthogonal = tree.addJoint(topRightInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(np.pi/2) @ SE3.Trans([jointLength,0,jointLength])))
topRightFoot = tree.addJoint(topRightOrthogonal, EndTip(numSides, r, SE3.Trans([jointLength,0,0]), 0.5, pathIndex=0))

topLeftInAxis = tree.addJoint(spine1, PrismaticJoint(numSides, r, twistJointLength,3,np.pi/5, SE3.Trans([jointLength,-jointLength,0]) @ SE3.Rx(np.pi/2)))
topLeftOrthogonal = tree.addJoint(topLeftInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(-np.pi/2) @ SE3.Trans([jointLength,0,jointLength])))
topLeftFoot = tree.addJoint(topLeftOrthogonal, EndTip(numSides, r, SE3.Trans([jointLength,0,0]), 0.5, pathIndex=0))


middleRightInAxis = tree.addJoint(0, PrismaticJoint(numSides, r, twistJointLength,3,np.pi/5, SE3.Trans([jointLength,0,jointLength/2]) @ SE3.Ry(np.pi/2)))
middleRightOrthogonal = tree.addJoint(middleRightInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Trans([jointLength,0,jointLength])))
middleRightFoot = tree.addJoint(middleRightOrthogonal, EndTip(numSides, r, SE3.Trans([jointLength,0,0]), 0.5, pathIndex=0))

middleLeftInAxis = tree.addJoint(0, PrismaticJoint(numSides, r, twistJointLength,3,np.pi/5, SE3.Trans([-jointLength,0,jointLength/2]) @ SE3.Ry(-np.pi/2)))
middleLeftOrthogonal = tree.addJoint(middleLeftInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(np.pi) @ SE3.Trans([jointLength,0,jointLength])))
middleLeftFoot = tree.addJoint(middleLeftOrthogonal, EndTip(numSides, r, SE3.Trans([jointLength,0,0]), 0.5, pathIndex=0))

spine2 = tree.addJoint(0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, SE3.Trans([0,-jointLength,jointLength/2]) @ SE3.Rz(-np.pi/2)))

bottomRightInAxis = tree.addJoint(spine2, PrismaticJoint(numSides, r, twistJointLength,3,np.pi/5, SE3.Trans([jointLength,-jointLength,0]) @ SE3.Rx(np.pi/2)))
bottomRightOrthogonal = tree.addJoint(bottomRightInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(-np.pi/2) @ SE3.Trans([jointLength,0,jointLength])))
bottomRightFoot = tree.addJoint(bottomRightOrthogonal, EndTip(numSides, r, SE3.Trans([jointLength,0,0]), 0.5, pathIndex=0))

bottomLeftInAxis = tree.addJoint(spine2, PrismaticJoint(numSides, r, twistJointLength,3,np.pi/5, SE3.Trans([jointLength,jointLength,0]) @ SE3.Rx(-np.pi/2)))
bottomLeftOrthogonal = tree.addJoint(bottomLeftInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(np.pi/2) @ SE3.Trans([jointLength,0,jointLength])))
bottomLeftFoot = tree.addJoint(bottomLeftOrthogonal, EndTip(numSides, r, SE3.Trans([jointLength,0,0]), 0.5, pathIndex=0))

hexapod = makeTubularKinematicTree(tree, plotSteps=False, optimize=False)

printedTree = origamiToPrinted(hexapod, 0.01)

optimizedHexapod = printedTree.squaredOptimize(showSteps=False)

plotPrintedTree(optimizedHexapod, "test")

optimizedHexapod.save("optimizedHexapodPrinted")

print(optimizedHexapod.detectCollisions(plot=True, includeEnds=True, debug=True))