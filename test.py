from KinematicChain import *


r = 1
numSides = 4
gg = KinematicChain(RevoluteJoint(numSides, r, 3*np.pi/2, SE3()), gimbal=True)
gg.appendGeneralizedGimbal(RevoluteJoint(numSides, r, 3*np.pi/2, SE3.Rx(np.pi/2)@SE3.Trans(-10,0,0)))
#gg.show(showSpheres=True, block=False)
gg.appendGeneralizedGimbal(RevoluteJoint(numSides, r, 3*np.pi/2, SE3.Rx(-3*np.pi/2)@SE3.Ry(np.pi/3)@SE3.Trans(0,12,0)))
#gg.show(showSpheres=True, block=False)
gg.appendGeneralizedGimbal(RevoluteJoint(numSides, r, 3*np.pi/2, SE3.Ry(-np.pi/4)@SE3.Trans(0,-15,30)))
gg.show(showSpheres=False, showScaleBar=False, showJointPoses=False, 
        showLinkPoses=False, showLinkPath=False, showAxisGrids=False, jointAxisScale=30)


"""
from numpy import sin, cos
r = 100 #units: mm
numSides = 6
# An example PUMA arm, based on Denavit-Hartenberg parameters from:
# C. S. G. Lee and M. Ziegler, 
# "Geometric Approach in Solving Inverse Kinematics of PUMA Robots," 
# in IEEE Transactions on Aerospace and Electronic Systems, 
# vol. AES-20, no. 6, pp. 695-706, Nov. 1984, doi: 10.1109/TAES.1984.310452.


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


realJointIndices = []

# Initialize the chain with a base waypoint
chain = KinematicChain(Waypoint(numSides, r, SE3.Trans(0,0,-660.4)@SE3.Rz(np.pi)), gimbal=True)

# The first joint axis is horizontal, 660.4 mm above the start base
realJointIndices.append(chain.appendGeneralizedGimbal(RevoluteJoint(numSides, r, np.pi, SE3())))
# Add the remaining revolute joints using the compact joint placement algorithm
for Pose in GlobalPoses[:-1]:
    realJointIndices.append(chain.appendGeneralizedGimbal(RevoluteJoint(numSides, r, np.pi, Pose)))
# The last DH parameters (and thus pose matrix) is for the end effector
#chain.append(EndTip(numSides, r, GlobalPoses[-1], 50), safe=False, relative=False)

chain.transformAll(SE3.Trans(0,0,660.4), recomputeBoundingBall=False)
chain.checkBallsAreNested()

# Plot the resulting chain
chain.show(block=False, showGroundPlane=True, groundPlaneScale=800, showSpheres=False, showScaleBar=False, showJointPoses=False)

numRandomConfigs = 3
randomConfigs = [chain.randomConfiguration() for _ in range(numRandomConfigs)]
for i in range(numRandomConfigs):
    chain.setConfiguration(randomConfigs[i])
    chain.checkBallsAreNested()
    chain.show(block=(i==numRandomConfigs-1), showGroundPlane=True, groundPlaneScale=800, showSpheres=False, showScaleBar=False, showJointPoses=False)

"""