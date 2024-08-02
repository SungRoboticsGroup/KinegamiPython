from KinematicTree import *
from testqtgraph import *
import time
from makeKinematicTree import *

tree = loadKinematicTree("optimizedHexapod")
tree = origamiToPrinted(tree, 0.0625)
plotPrintedTree(tree, "optimizedHexapodPrinted")
tree.export3DKinematicTree("optimizedHexapodPrinted")
print(tree.detectCollisions(plot=False, debug=True))

# tree = loadKinematicTree("algorithmHand")

# r = 1
# numSides = 4
# jointLength = 0.45/0.245
# unextendedRevoluteJointLength = RevoluteJoint(numSides, r, np.pi, SE3()).neutralLength
# extensionLength = (jointLength - unextendedRevoluteJointLength)/2


# spec = JointSpecificationTree(Waypoint(numSides, r, Pose=SE3()))
# palmJoint = ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, SE3.Trans(4,0,3)@SE3.Ry(-np.pi/4)@SE3.Rx(np.pi/4))
# palm = spec.addJoint(0, palmJoint, relative=True)


# thumb0Joint = ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, SE3.Trans(jointLength,0,0)@SE3.Rx(np.pi/2))
# thumb0 = spec.addJoint(palm, thumb0Joint, relative=True)

# thumb1 = spec.addJoint(thumb0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
#                                            SE3.Trans(4,0,0)), relative=True)

# thumbEnd = spec.addJoint(thumb1, EndTip(numSides, r, 
#         SE3.Trans(3,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), 1), relative=True)

# tree = makeTubularKinematicTree(spec, plotSteps=False, optimize=False)

# tree2 = tree.optimizeJointsDifferentiable()

# tree2.show()
