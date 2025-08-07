import numpy as np
import roboticstoolbox as rtb
from spatialmath import SE3
from makeKinematicTree import *
from Kinematics import *

r = 1
numSides = 4

spec = JointSpecificationTree(RevoluteJoint(numSides, r, np.pi, SE3()))
i = spec.addJoint(0, PrismaticJoint(numSides,r,3,3,np.pi/5,SE3(4,0,0)@SE3.Ry(np.pi/2)), relative=True)
i = spec.addJoint(i, RevoluteJoint(numSides,r,np.pi,SE3.Trans(0,0,4)@SE3.Ry(2*np.pi/3)), relative=True)
i = spec.addJoint(i, RevoluteJoint(numSides,r,np.pi,SE3.Trans(1,0,0)@SE3.Rx(np.pi/4)), relative=True)
t = makeTubularKinematicTree(spec, plotSteps=False, orientUp=True)
t.show(block=False)

robot = robotFromTree(t)

q = [0,0,0,0]  # Two joint angles (radians)

# Plot the robot
# robot.plot(q, block=False)
# plt.show()

# Desired end-effector position (x, y)

t2 = copy.deepcopy(t)
t2.setJointState(0,np.pi/2)
t2.setJointState(2,np.pi/3)

target_position = SE3(2, 3, 1)  # t2.Joints[-1].Pose

# Solve inverse kinematics
solution = robot.ikine_LM(target_position, q0=q)

# Display results
if solution.success:
    print(f"IK Solution Found: {solution.q}")
    # Verify by computing forward kinematics
    fk_solution = robot.fkine(solution.q)
    print(f"FK of Solution: {fk_solution}")
else:
    print("IK Solver Failed")
    tree2, loss = makeChainFitPoses(t, [target_position], maxiter=100)
    print(f"Loss : {loss}")

    robot2 = robotFromTree(tree2)
    solution = robot2.ikine_LM(target_position, q0=q)
    if solution.success:
        print("worked")
    else:
        print("no worked")


