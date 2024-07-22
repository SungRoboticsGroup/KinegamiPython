import os
import sys
this_dir = os.path.dirname(__file__)
main_dir = os.path.abspath(os.path.join(this_dir, '../..'))
sys.path.append(main_dir)

# Example 1A
from KinematicTree import *
from testqtgraph import *

scale = 30 #in mm
r = 1

angle = np.pi/3
screwR = 1.5/scale

jointLength = 75/scale
unextendedRevoluteJointLength = PrintedOrthogonalRevoluteJoint(r, -angle, angle, SE3(), screwR).neutralLength
extensionLength = jointLength - unextendedRevoluteJointLength

print(extensionLength)

twistJointLength = 60/scale
unextendedPrismaticJointLength = PrintedPrismaticJoint(r, 0, SE3(), screwR).neutralLength
prismaticExtensionLength = twistJointLength - unextendedPrismaticJointLength 

print(prismaticExtensionLength)

tree = KinematicTree[PrintedJoint](PrintedWaypoint(r, SE3(), screwR))

spine1Joint = PrintedOrthogonalRevoluteJoint(r, -angle, angle, SE3.Trans([0,jointLength, jointLength/2]) @ SE3.Rz(np.pi/2), screwR)
spine1Joint.extendSegment(extensionLength)
spine1 = tree.addJoint(0, spine1Joint, fixedOrientation=True, fixedPosition=True, safe=False)

topRightJoint1 = PrintedPrismaticJoint(r, 0, SE3.Trans([jointLength, jointLength, 0]) @ SE3.Rx(-np.pi/2), screwR)
topRightJoint1.extendSegment(prismaticExtensionLength)
topRightInAxis = tree.addJoint(spine1, topRightJoint1, fixedOrientation=True, fixedPosition=True, safe=False)

topRightJoint2 = PrintedOrthogonalRevoluteJoint(r, -angle, angle, SE3.Rz(np.pi/2) @ SE3.Trans([jointLength,0,jointLength]), screwR)
topRightJoint2.extendSegment(extensionLength)
topRightOrthogonal = tree.addJoint(topRightInAxis, topRightJoint2, fixedOrientation=True, fixedPosition=True, safe=False)

topRightFoot = tree.addJoint(topRightOrthogonal, PrintedTip(r, SE3.Trans([jointLength,0,0]), screwR, pathIndex=0), fixedOrientation=True, fixedPosition=True, safe=False)

topLeftJoint1 = PrintedPrismaticJoint(r, 0, SE3.Trans([jointLength, -jointLength, 0]) @ SE3.Rx(np.pi/2), screwR)
topLeftJoint1.extendSegment(prismaticExtensionLength)
topLeftInAxis = tree.addJoint(spine1, topLeftJoint1, fixedOrientation=True, fixedPosition=True, safe=False)

topLeftJoint2 = PrintedOrthogonalRevoluteJoint(r, -angle, angle, SE3.Rz(-np.pi/2) @ SE3.Trans([jointLength,0,jointLength]), screwR)
topLeftJoint2.extendSegment(extensionLength)
topLeftOrthogonal = tree.addJoint(topLeftInAxis, topLeftJoint2, fixedOrientation=True, fixedPosition=True, safe=False)

topLeftFoot = tree.addJoint(topLeftOrthogonal, PrintedTip(r, SE3.Trans([jointLength,0,0]), screwR, pathIndex=0), fixedOrientation=True, fixedPosition=True, safe=False)

middleRightJoint1 = PrintedPrismaticJoint(r, 0, SE3.Trans([jointLength, 0, jointLength/2]) @ SE3.Ry(np.pi/2), screwR)
middleRightJoint1.extendSegment(prismaticExtensionLength)
middleRightInAxis = tree.addJoint(0, middleRightJoint1, fixedOrientation=True, fixedPosition=True, safe=False)

middleRightJoint2 = PrintedOrthogonalRevoluteJoint(r, -angle, angle, SE3.Trans([jointLength,0,jointLength]), screwR)
middleRightJoint2.extendSegment(extensionLength)
middleRightOrthogonal = tree.addJoint(middleRightInAxis, middleRightJoint2, fixedOrientation=True, fixedPosition=True, safe=False)

middleRightFoot = tree.addJoint(middleRightOrthogonal, PrintedTip(r, SE3.Trans([jointLength,0,0]), screwR, pathIndex=0), fixedOrientation=True, fixedPosition=True, safe=False)


middleLeftJoint1 = PrintedPrismaticJoint(r, 0, SE3.Trans([-jointLength, 0, jointLength/2]) @ SE3.Ry(-np.pi/2), screwR)
middleLeftJoint1.extendSegment(prismaticExtensionLength)
middleLeftInAxis = tree.addJoint(0, middleLeftJoint1, fixedOrientation=True, fixedPosition=True, safe=False)

middleLeftJoint2 = PrintedOrthogonalRevoluteJoint(r, -angle, angle, SE3.Rz(np.pi) @ SE3.Trans([jointLength,0,jointLength]), screwR)
middleLeftJoint2.extendSegment(extensionLength)
middleLeftOrthogonal = tree.addJoint(middleLeftInAxis, middleLeftJoint2, fixedOrientation=True, fixedPosition=True, safe=False)

middleLeftFoot = tree.addJoint(middleLeftOrthogonal, PrintedTip(r, SE3.Trans([jointLength,0,0]), screwR, pathIndex=0), fixedOrientation=True, fixedPosition=True, safe=False)

spine2Joint = PrintedOrthogonalRevoluteJoint(r, -angle, angle, SE3.Trans([0,-jointLength, jointLength/2]) @ SE3.Rz(-np.pi/2), screwR)
spine2Joint.extendSegment(extensionLength)
spine2 = tree.addJoint(0, spine2Joint, fixedOrientation=True, fixedPosition=True, safe=False)

bottomRightJoint1 = PrintedPrismaticJoint(r, 0, SE3.Trans([jointLength, -jointLength, 0]) @ SE3.Rx(np.pi/2), screwR)
bottomRightJoint1.extendSegment(prismaticExtensionLength)
bottomRightInAxis = tree.addJoint(spine2, bottomRightJoint1, fixedOrientation=True, fixedPosition=True, safe=False)

bottomRightJoint2 = PrintedOrthogonalRevoluteJoint(r, -angle, angle, SE3.Rz(-np.pi/2) @ SE3.Trans([jointLength,0,jointLength]), screwR)
bottomRightJoint2.extendSegment(extensionLength)
bottomRightOrthogonal = tree.addJoint(bottomRightInAxis, bottomRightJoint2, fixedOrientation=True, fixedPosition=True, safe=False)

bottomRightFoot = tree.addJoint(bottomRightOrthogonal, PrintedTip(r, SE3.Trans([jointLength,0,0]), screwR, pathIndex=0), fixedOrientation=True, fixedPosition=True, safe=False)

bottomLeftJoint1 = PrintedPrismaticJoint(r, 0, SE3.Trans([jointLength, jointLength, 0]) @ SE3.Rx(-np.pi/2), screwR)
bottomLeftJoint1.extendSegment(prismaticExtensionLength)
bottomLeftInAxis = tree.addJoint(spine2, bottomLeftJoint1, fixedOrientation=True, fixedPosition=True, safe=False)

bottomLeftJoint2 = PrintedOrthogonalRevoluteJoint(r, -angle, angle, SE3.Rz(np.pi/2) @ SE3.Trans([jointLength,0,jointLength]), screwR)
bottomLeftJoint2.extendSegment(extensionLength)
bottomLeftOrthogonal = tree.addJoint(bottomLeftInAxis, bottomLeftJoint2, fixedOrientation=True, fixedPosition=True, safe=False)

bottomLeftFoot = tree.addJoint(bottomLeftOrthogonal, PrintedTip(r, SE3.Trans([jointLength,0,0]), screwR, pathIndex=0), fixedOrientation=True, fixedPosition=True, safe=False)

plotPrintedTree(tree, "test")
#tree.export3DKinematicTree("manualHexapodPrinted")