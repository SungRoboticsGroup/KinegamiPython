import os
import sys
this_dir = os.path.dirname(__file__)
main_dir = os.path.abspath(os.path.join(this_dir, '../..'))
sys.path.append(main_dir)

# Example 1A
from KinematicTree import *
from testqtgraph import *
from makeKinematicTree import *
from optimizationFunctions import *
from randomTree import *
from multiObjectiveOptimizationFunctions import *
from kinematicsOptimizationFunctions import *

from numpy import sin, cos
r = 100 #units: mm
numSides = 6

"""
An example PUMA arm, based on Denavit-Hartenberg parameters from:
C. S. G. Lee and M. Ziegler, 
"Geometric Approach in Solving Inverse Kinematics of PUMA Robots," 
in IEEE Transactions on Aerospace and Electronic Systems, 
vol. AES-20, no. 6, pp. 695-706, Nov. 1984, doi: 10.1109/TAES.1984.310452.
"""

# Function to convert Denavit-Hartenberg parameters to relative pose matrices
def relativePoseFromDHOriginal(theta, alpha, a, d):
    return SE3([[cos(theta), -cos(alpha)*sin(theta), sin(alpha)*sin(theta), a*cos(theta)],
                [sin(theta), cos(alpha)*cos(theta), -sin(alpha)*cos(theta), a*sin(theta)],
                [0,          sin(alpha),             cos(alpha),            d],
                [0,          0,                      0,                     1]])

# Construct pose matrices from the Denavit-Hartenberg parameters for a PUMA arm
# from Lee and Ziegler 1984
RelativePoses = [relativePoseFromDHOriginal(np.pi/2, -np.pi/2, 0, 0),
         relativePoseFromDHOriginal(0, 0, 431.8, 149.09),
         relativePoseFromDHOriginal(np.pi/2, np.pi/2, -20.32, 0),
         relativePoseFromDHOriginal(0, -np.pi/2, 0, 433.07),
         relativePoseFromDHOriginal(0, np.pi/2, 0, 0),
         relativePoseFromDHOriginal(0, 0, 0, 56.25)]

GlobalPoses = [RelativePoses[0]]
for Pose in RelativePoses[1:]:
    GlobalPoses.append(GlobalPoses[-1] @ Pose)



# Initialize the chain with a base waypoint
chain = KinematicTree(Waypoint(numSides, r, SE3.Trans(0,0,-660.4)@SE3.Rz(np.pi)))

# The first joint axis is horizontal, 660.4 mm above the start base
currentJoint = chain.addJoint(0, RevoluteJoint(numSides, r, np.pi, SE3()), safe=False)
# Add the remaining revolute joints using the compact joint placement algorithm
for Pose in GlobalPoses[:-1]:
    currentJoint = chain.addJoint(currentJoint, RevoluteJoint(numSides, r, np.pi, Pose), safe=False, relative=False)
# The last DH parameters (and thus pose matrix) is for the end effector
chain.addJoint(currentJoint, EndTip(numSides, r, GlobalPoses[-1], 50), safe=False, relative=False)

chain.transformAll(SE3.Trans(0,0,660.4))

# Plot the resulting chain
# chain.show()

chain.setJointState(3, np.pi)
pose1 = chain.Joints[-1].Pose

chain.setJointState(4, np.pi/3)
chain.setJointState(6, np.pi*2/3)
pose2 = chain.Joints[-1].Pose

newChain = squaredOptimizeForEndEffectorPoses(chain, [len(chain.Joints) - 1] * 2, [pose1, pose2], parallelize=True, guarantee=True)

# newChain.show()

newChain.save("testChain")
# newChain = loadKinematicTree("save/testChain")
robotChain = robotFromChain(newChain)
for pose in [pose1, pose2]:
    sol = robotChain.ikine_LM(pose, q0=[0]*len(robotChain.joints()))
    if not sol.success:
        print("Could not fit pose in the end")
    else:
        newChain.setConfiguration(sol.q, realJointsOnly=True)
        newChain.show()
# r = 1
# numSides = 4
# jointLength = 0.45/0.245
# unextendedRevoluteJointLength = RevoluteJoint(numSides, r, np.pi, SE3()).neutralLength
# extensionLength = (jointLength - unextendedRevoluteJointLength)/2

# # spec = JointSpecificationTree(RevoluteJoint(numSides, r, np.pi, SE3()))
# # i = spec.addJoint(0, PrismaticJoint(numSides,r,3,3,np.pi/5,SE3(4,0,0)@SE3.Ry(np.pi/2)), relative=True)
# # i = spec.addJoint(i, RevoluteJoint(numSides,r,np.pi,SE3.Trans(0,0,4)@SE3.Ry(2*np.pi/3)), relative=True)
# # i = spec.addJoint(i, RevoluteJoint(numSides,r,np.pi,SE3.Trans(1,0,0)@SE3.Rx(np.pi/4)), relative=True)
# # tree = makeTubularKinematicTree(spec, plotSteps=False, orientUp=True)

# spec = JointSpecificationTree(Waypoint(numSides, r, Pose=SE3()))
# palmJoint = ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, SE3.Trans(4,0,3)@SE3.Ry(-np.pi/4)@SE3.Rx(np.pi/4))
# palm = spec.addJoint(0, palmJoint, relative=True)


# thumb0Joint = ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, SE3.Trans(jointLength,0,0)@SE3.Rx(np.pi/2))
# thumb0 = spec.addJoint(palm, thumb0Joint, relative=True)

# thumb1 = spec.addJoint(thumb0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
#                                            SE3.Trans(4,0,0)), relative=True)

# thumbEnd = spec.addJoint(thumb1, EndTip(numSides, r, 
#         SE3.Trans(3,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), 1), relative=True)


# pointer1 = spec.addJoint(0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
#             SE3.Trans(0,0,8)@SE3.Ry(-np.pi/2)@SE3.Rx(np.pi)), relative=True)

# pointer2 = spec.addJoint(pointer1, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
#                         SE3.Trans(4,0,0)), relative=True)


# pointerEnd = spec.addJoint(pointer2, EndTip(numSides, r, 
#         SE3.Trans(4,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), 1), relative=True)

        

# middle0 = spec.addJoint(0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
#         SE3.Trans(-2.25,0,4)@SE3.Ry(-np.pi/2)@SE3.Rx(np.pi/2)), relative=True)

# middle1 = spec.addJoint(middle0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
#                         SE3.Trans(4,0,0)@SE3.Rx(np.pi/2)), relative=True)

# middle2 = spec.addJoint(middle1, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
#                         SE3.Trans(4,0,0)), relative=True)

# middleEnd = spec.addJoint(middle2, EndTip(numSides, r, 
#         SE3.Trans(4,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), 1), relative=True)

        
# tree = makeTubularKinematicTree(spec, plotSteps=False)

# print(len(tree.Joints))

# states = [[0]* len(tree.Joints),[np.pi/2,np.pi/4,np.pi/2,np.pi/3,np.pi/3,np.pi/3,np.pi/3,np.pi/3,np.pi/2,np.pi/4,np.pi/2,np.pi/3,np.pi/3,np.pi/3,np.pi/3,np.pi/3]]
# for t in method1(tree, states, parallelize=True, guarantee=True):
#     t.show()

