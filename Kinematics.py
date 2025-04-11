import numpy as np
import roboticstoolbox as rtb
from spatialmath import SE3
from scipy.spatial.transform import Rotation as R
from KinematicTree import *

def robotFromTree(tree):
    arr = []
    qlim = [[],[]]
    angles = tree.Joints[0].Pose.eul()
    arr.append(rtb.ET.Rz(angles[0]))
    arr.append(rtb.ET.Ry(angles[1]))
    arr.append(rtb.ET.Rz(angles[2]))
    for i in range(0,len(tree.Joints)):
        pose = tree.Joints[i].Pose
        if (isinstance(tree.Joints[i], RevoluteJoint)):
            arr.append(rtb.ET.Rz())
            qlim[0].append(-np.pi)
            qlim[1].append(np.pi)
        elif (isinstance(tree.Joints[i], PrismaticJoint)):
            arr.append(rtb.ET.tz())
            qlim[0].append(0)
            qlim[1].append(2)
        
        if (i < len(tree.Joints) - 1):
            pose2 =  tree.Joints[i + 1].Pose
            transform = pose.inv() * pose2
            angles = transform.eul()
            arr.append(rtb.ET.tx(transform.t[0]))
            arr.append(rtb.ET.ty(transform.t[1]))
            arr.append(rtb.ET.tz(transform.t[2]))
            arr.append(rtb.ET.Rz(angles[0]))
            arr.append(rtb.ET.Ry(angles[1]))
            arr.append(rtb.ET.Rz(angles[2]))
    
    robot = rtb.ETS(arr)
    robot.qlim = np.array(qlim)

    return robot

def makeChainFitPoses(subject, poses, maxiter, parallelize=False):
    def poseLoss(t):
        loss = 0
        robot = robotFromTree(t)
        for pose in poses:
            sol = robot.ikine_LM(pose, q0=[0]*len(subject.Joints))
            if not sol.success:
                result = robot.fkine(sol.q)
                trans = np.linalg.norm(result.t - pose.t)
                # rot = 2 * np.arccos(np.clip(np.dot(R.from_matrix(pose.R).as_quat(), R.from_matrix(result.R).as_quat()), -1.0, 1.0))
                loss += trans# + t.r * rot
        return loss
    
    def objective(params):
        tree = subject.copyAbbreviatedSelf()

        for i in range(0, len(params),2):
            tree.Joints[int(i/2)].applyTransformationToPose(SE3(params[i],params[i+1],0))
        
        return poseLoss(tree)
    
    global ctr
    ctr = 0

    global batch_objective_function
    def batch_objective_function(X):
        losses = np.array([objective(x) for x in X])
        global ctr
        ctr += 1
        print(f"{ctr}: {np.min(losses)}")
        return losses

    start = time.time()


    initialGuess = [0]*len(subject.Joints) * 2

    n_particles = 10

    dist = np.max([np.linalg.norm(subject.Joints[0].Pose.t - pose.t) for pose in poses])
    bounds = [(-dist*2, dist*2)]*len(initialGuess)

    init_pos = np.tile(np.array(initialGuess, dtype='float64'), (n_particles,1))
    #add random noise
    noise = np.zeros_like(init_pos[1:])
    for i in range(0,3):
        noise[:, i] = np.random.uniform(-dist,dist, n_particles - 1)
    for i in range(3, 6):
        noise[:, i] = np.random.uniform(-dist,dist, n_particles - 1)
    init_pos[1:] += noise

#   bounds=(np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds]))
    optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,dimensions=len(initialGuess),options={'c1':0.7, 'c2':0.5, 'w':0.2},bounds=None,init_pos=None)
    loss, result = optimizer.optimize(batch_objective_function, iters=maxiter,verbose=False, n_processes=n_particles if parallelize else None)

    endTree = subject.copyAbbreviatedSelf()
    for i in range(0, len(result),2):
        endTree.Joints[int(i/2)].applyTransformationToPose(SE3(result[i],result[i+1],0))

    if (loss != 0):
        print("Unable to fit all poses")

    print(f"Time taken: {time.time() - start}")
    return endTree, loss

