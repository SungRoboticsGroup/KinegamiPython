# -*- coding: utf-8 -*-
from KinematicChain import *
from numpy import sin, cos
r = 100 #units: mm
numSides = 6
guarantee = False


"""
DH parameters from:
C. s. g. Lee and M. Ziegler, 
"Geometric Approach in Solving Inverse Kinematics of PUMA Robots," 
in IEEE Transactions on Aerospace and Electronic Systems, 
vol. AES-20, no. 6, pp. 695-706, Nov. 1984, doi: 10.1109/TAES.1984.310452.
"""
def relativePoseFromDHOriginal(theta, alpha, a, d):
    return SE3([[cos(theta), -cos(alpha)*sin(theta), sin(alpha)*sin(theta), a*cos(theta)],
                [sin(theta), cos(alpha)*cos(theta), -sin(alpha)*cos(theta), a*sin(theta)],
                [0,          sin(alpha),             cos(alpha),            d],
                [0,          0,                      0,                     1]])


Poses = [relativePoseFromDHOriginal(np.pi/2, -np.pi/2, 0, 0),
         relativePoseFromDHOriginal(0, 0, 431.8, 149.09),
         relativePoseFromDHOriginal(np.pi/2, np.pi/2, -20.32, 0),
         relativePoseFromDHOriginal(0, -np.pi/2, 0, 433.07),
         relativePoseFromDHOriginal(0, np.pi/2, 0, 0),
         relativePoseFromDHOriginal(0, 0, 0, 56.25)]

chain = KinematicChain(Waypoint(numSides, r, SE3.Rz(np.pi)))
chain.append(RevoluteJoint(numSides, r, np.pi, SE3.Trans(0,0,660.4)), guarantee=guarantee)
for Pose in Poses[:-1]:
    chain.append(RevoluteJoint(numSides, r, np.pi, Pose), guarantee=guarantee)
chain.append(EndTip(numSides, r, Poses[-1], 50), guarantee=guarantee)
chain.show(block=False, showLinkPath=False, showJointPoses=False, showLinkPoses=False, showAxisGrids=False)

chain.translateJointAlongKinematicAxis(2, 400)
chain.translateJointAlongKinematicAxis(3, 600)
chain.translateJointAlongKinematicAxis(4, -550)
chain.rotateJointAboutKinematicAxis(4, np.pi/2, propogate=False)
chain.translateJointAlongKinematicAxis(5, -550)
chain.translateJointAlongKinematicAxis(6, -400)
chain.translateJointAlongKinematicAxis(7, -300)
chain.show(block=False, showLinkPath=False, showJointPoses=False, showLinkPoses=False, showAxisGrids=False)
chain.creasePattern().show()


