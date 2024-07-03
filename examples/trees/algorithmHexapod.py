import os
import sys
this_dir = os.path.dirname(__file__)
main_dir = os.path.abspath(os.path.join(this_dir, '../..'))
sys.path.append(main_dir)

# Example 1A
from KinematicTree import *
from testqtgraph import *
from makeKinematicTree import *

r = 0.5
numSides = 4

tree = JointSpecificationTree(Waypoint(numSides, r, SE3()))

spine1 = tree.addJoint(0, RevoluteJoint(numSides, r, np.pi, SE3.Trans([0,3,1]) @ SE3.Rz(np.pi/2)))

topRightInAxis = tree.addJoint(spine1, PrismaticJoint(numSides, r, 1,3,np.pi/5, SE3.Trans([3,3,0]) @ SE3.Rx(-np.pi/2)))
topRightOrthogonal = tree.addJoint(topRightInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(np.pi/2) @ SE3.Trans([2,0,2])))
topRightFoot = tree.addJoint(topRightOrthogonal, EndTip(numSides, r, SE3.Trans([3,0,0]), 0.5, pathIndex=0))

topLeftInAxis = tree.addJoint(spine1, PrismaticJoint(numSides, r, 1,3,np.pi/5, SE3.Trans([3,-3,0]) @ SE3.Rx(np.pi/2)))
topLeftOrthogonal = tree.addJoint(topLeftInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(-np.pi/2) @ SE3.Trans([2,0,2])))
topLeftFoot = tree.addJoint(topLeftOrthogonal, EndTip(numSides, r, SE3.Trans([3,0,0]), 0.5, pathIndex=0))


middleRightInAxis = tree.addJoint(0, PrismaticJoint(numSides, r, 1,3,np.pi/5, SE3.Trans([3,0,1]) @ SE3.Ry(np.pi/2)))
middleRightOrthogonal = tree.addJoint(middleRightInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Trans([2,0,2])))
middleRightFoot = tree.addJoint(middleRightOrthogonal, EndTip(numSides, r, SE3.Trans([3,0,0]), 0.5, pathIndex=0))

middleLeftInAxis = tree.addJoint(0, PrismaticJoint(numSides, r, 1,3,np.pi/5, SE3.Trans([-3,0,1]) @ SE3.Ry(-np.pi/2)))
middleLeftOrthogonal = tree.addJoint(middleLeftInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(np.pi) @ SE3.Trans([2,0,2])))
middleLeftFoot = tree.addJoint(middleLeftOrthogonal, EndTip(numSides, r, SE3.Trans([3,0,0]), 0.5, pathIndex=0))

spine2 = tree.addJoint(0, RevoluteJoint(numSides, r, np.pi, SE3.Trans([0,-3,1]) @ SE3.Rz(-np.pi/2)))

bottomRightInAxis = tree.addJoint(spine2, PrismaticJoint(numSides, r, 1,3,np.pi/5, SE3.Trans([3,-3,0]) @ SE3.Rx(np.pi/2)))
bottomRightOrthogonal = tree.addJoint(bottomRightInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(-np.pi/2) @ SE3.Trans([2,0,2])))
bottomRightFoot = tree.addJoint(bottomRightOrthogonal, EndTip(numSides, r, SE3.Trans([3,0,0]), 0.5, pathIndex=0))

bottomLeftInAxis = tree.addJoint(spine2, PrismaticJoint(numSides, r, 1,3,np.pi/5, SE3.Trans([3,3,0]) @ SE3.Rx(-np.pi/2)))
bottomLeftOrthogonal = tree.addJoint(bottomLeftInAxis, RevoluteJoint(numSides, r, np.pi, SE3.Rz(np.pi/2) @ SE3.Trans([2,0,2])))
bottomLeftFoot = tree.addJoint(bottomLeftOrthogonal, EndTip(numSides, r, SE3.Trans([3,0,0]), 0.5, pathIndex=0))

hexapod = makeTubularKinematicTree(tree, plotSteps=False, optimize=False)

hexapod.show()

optimizedHexapod = hexapod.postOptimize()

optimizedHexapod.show()

print(hexapod.detectCollisions(plot=True))