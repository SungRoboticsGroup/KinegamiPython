import os
import sys
this_dir = os.path.dirname(__file__)
main_dir = os.path.abspath(os.path.join(this_dir, '../..'))
sys.path.append(main_dir)

from KinematicTree import *
from KinematicChain import *
from testqtgraph import *
from makeKinematicTree import *

r = 0.02
numSides = 4
extensionLength = 0.035

# Lengths
l1 = 0.1
l2 = 0.1
hip_to_hip = 2.5*l1 # version PET
first_waypoint = 0.09 # version PET
second_waypoint = hip_to_hip + 0.014 +0.035 # version PET
indent = -0.08

# Create the chain version
chain = KinematicChain(Waypoint(numSides, r, SE3(), pathIndex=0))
knee1 = chain.append(ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, SE3.Tx(l2)), 
             safe=False, fixedOrientation=True)
hip1 = chain.append(ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, 
                           SE3.Rz(-np.pi/3) @ SE3.Tx(l1) @ SE3.Rz(-np.pi/3)), 
                           safe=False, fixedOrientation=True)
midWaypoint1 = chain.append(Waypoint(numSides, r, SE3.Trans(first_waypoint, 0, -indent), pathIndex=0), 
             safe=False, fixedOrientation=True)
midWaypoint2 = chain.append(Waypoint(numSides, r, SE3.Trans(second_waypoint-first_waypoint, 0, 0), 
                      pathIndex=0), 
             safe=False, fixedOrientation=True)
hip2 = chain.append(ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, 
                SE3.Trans(hip_to_hip-second_waypoint,0,indent) @ SE3.Rz(np.pi)), 
             safe=False, fixedOrientation=True)
knee2 = chain.append(ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, 
                SE3.Rz(np.pi/3) @ SE3.Tx(l1) @ SE3.Rz(np.pi/3)), 
             safe=False, fixedOrientation=True)
foot2 = chain.append(Waypoint(numSides, r, SE3.Tx(l2), pathIndex=0), 
             safe=False, fixedOrientation=True)
for j in range(midWaypoint1):
    z = chain.Joints[j].Pose.t[2]
    chain.transformJoint(j, SE3.Tz(-z), propogate=False)
mwz1 = chain.Joints[midWaypoint1].Pose.t[2]
chain.transformJoint(midWaypoint1, SE3.Tz(-indent-mwz1), propogate=False)
mwz2 = chain.Joints[midWaypoint2].Pose.t[2]
chain.transformJoint(midWaypoint2, SE3.Tz(-indent-mwz2), propogate=False)
for j in range(len(chain.Joints)-1, midWaypoint2, -1):
    z = chain.Joints[j].Pose.t[2]
    chain.transformJoint(j, SE3.Tz(-z), propogate=False)
chain.creasePattern().show(block=False)
chain.show(showScaleBar=False)

# Create the tree specification from the global poses from the chain
hip1leftJoint = copy.deepcopy(chain.Joints[hip1])
hip1rightJoint = copy.deepcopy(chain.Joints[hip1])
knee1leftJoint = copy.deepcopy(chain.Joints[knee1])
knee1rightJoint = copy.deepcopy(chain.Joints[knee1])
foot1leftJoint = copy.deepcopy(chain.Joints[0])
foot1rightJoint = copy.deepcopy(chain.Joints[0])
hip2leftJoint = copy.deepcopy(chain.Joints[hip2])
hip2rightJoint = copy.deepcopy(chain.Joints[hip2])
knee2leftJoint = copy.deepcopy(chain.Joints[knee2])
knee2rightJoint = copy.deepcopy(chain.Joints[knee2])
foot2leftJoint = copy.deepcopy(chain.Joints[foot2])
foot2rightJoint = copy.deepcopy(chain.Joints[foot2])

treeSpec = JointSpecificationTree(Waypoint(numSides, r, chain.Joints[midWaypoint1].Pose, pathIndex=0))
hip1leftIndex = treeSpec.addJoint(0, hip1leftJoint, relative=False)
knee1leftIndex = treeSpec.addJoint(hip1leftIndex, knee1leftJoint, relative=False)
foot1leftIndex = treeSpec.addJoint(knee1leftIndex, foot1leftJoint, relative=False)
hip2leftIndex = treeSpec.addJoint(0, hip2leftJoint, relative=False)
knee2leftIndex = treeSpec.addJoint(hip2leftIndex, knee2leftJoint, relative=False)
foot2leftIndex = treeSpec.addJoint(knee2leftIndex, foot2leftJoint, relative=False)
hip1rightIndex = treeSpec.addJoint(0, hip1rightJoint, relative=False)
knee1rightIndex = treeSpec.addJoint(hip1rightIndex, knee1rightJoint, relative=False)
foot1rightIndex = treeSpec.addJoint(knee1rightIndex, foot1rightJoint, relative=False)
hip2rightIndex = treeSpec.addJoint(0, hip2rightJoint, relative=False)
knee2rightIndex = treeSpec.addJoint(hip2rightIndex, knee2rightJoint, relative=False)
foot2rightIndex = treeSpec.addJoint(knee2rightIndex, foot2rightJoint, relative=False)

treeFromAlgorithm = makeTubularKinematicTree(treeSpec, plotSteps=False)
#treeFromAlgorithm.show(showScaleBar=False)

optimizedTree = treeFromAlgorithm.squaredOptimize(showSteps=False, childParentRatio = 0, streamline=True, guarantee=True, parallelize=True)
optimizedTree.show(showScaleBar=False)


