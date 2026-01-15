from KinematicTree import *
from treeTraversals import *
import pyswarms as ps


def optimizeJointPlacement(subject, index, maxiter, tol, penaltyScale, 
                           childFraction = 1, ignoreLater = False, parallelize = False, 
                           verbose=True, power=2, includeCollisionPenalty=True, 
                           configurations=None, retryingWithPenalty=False):

    # If includeCollisionPenalty is True, first try without it
    original_includeCollisionPenalty = includeCollisionPenalty
    verbose = True
    if original_includeCollisionPenalty and not retryingWithPenalty:
        try:
            if verbose:
                print("Trying optimization without collision penalty...")
            
            # Temporarily disable collision penalty for first attempt
            includeCollisionPenalty = False
            
        except Exception as e:
            if verbose:
                print(f"Setup for no-penalty optimization failed: {str(e)}")
            includeCollisionPenalty = True

    if configurations is None:
        configurations = [[0] * len(subject.Joints)]

    for configuration in configurations:
        copied_subject = copy.deepcopy(subject)
        for i in range(0, len(copied_subject.Joints)):
            if (not isWaypoint(copied_subject.Joints[i])):
                copied_subject.setJointState(i,configuration[i])
                copied_subject.Joints[i].recomputeCollisionCapsules()
        if copied_subject.detectCollisions(debug=True) > 0:
            print(f"Warning: Initial tree in optimizeJointPlacement contains collisions in configuration {configuration}.")

    subjects = [copy.deepcopy(subject) for _ in configurations]
    for i in range(0,len(configurations)):
        subjects[i].setConfiguration(configurations[i])

    movedJointIndex = index
    collisionMatrices = subject.buildCollisionMatrices()
    
    def objective(params, returnWhich : bool = False) -> float:
        tree = subject.copyAbbreviatedSelf(ignoreLater, index)

        translation = params[0]
        rotation = params[1]

        transform = SE3.Trans([0,0,translation]) @ SE3.Rz(rotation)

        pathNonExistancePenalty = penaltyScale * (len(subject.Joints) ** 4) * (len(subject.Children[index]) + 1)
        linkLossReversedZhat = pathNonExistancePenalty

        #try just moving it
        if tree.transformJoint(index, transform, safe=True, relative=True, recomputeBoundingBall=False):
            linkLossSameZhat = linkLoss(tree, 
                                        index, 
                                        power=power,
                                        collisionMatrices=collisionMatrices,
                                        movedJointIndex=movedJointIndex,
                                        includeCollisionPenalty=includeCollisionPenalty, 
                                        configurations=configurations)
        else:
            linkLossSameZhat = pathNonExistancePenalty
        
        #try switching zhat
        tree.Joints[index].reverseZhat()
        if (linkLossSameZhat == pathNonExistancePenalty or
             not tree.transformJoint(index, SE3(), safe=True, relative=True, propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False)):
            tree2 = subject.copyAbbreviatedSelf(ignoreLater, index)
            tree2.Joints[index].reverseZhat()
            if tree2.transformJoint(index, SE3.Trans([0,0,-translation]) @ SE3.Rz(-rotation), safe=True, relative=True, 
                                    propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False):
                linkLossReversedZhat = linkLoss(tree2, 
                                                index, 
                                                power=power,
                                                collisionMatrices=collisionMatrices,
                                                movedJointIndex=movedJointIndex,
                                                includeCollisionPenalty=includeCollisionPenalty,
                                                configurations=configurations)
            else:
                linkLossReversedZhat = pathNonExistancePenalty
        else:
            linkLossReversedZhat = linkLoss(tree, 
                                            index, 
                                            power=power,
                                            collisionMatrices=collisionMatrices,
                                            movedJointIndex=movedJointIndex,
                                            includeCollisionPenalty=includeCollisionPenalty, 
                                            configurations=configurations)

        if returnWhich:
            if linkLossSameZhat <= linkLossReversedZhat:
                return 1
            else:
                return 2
        else:
            return min(linkLossSameZhat,linkLossReversedZhat)

    start = time.time()

    joint = subject.Joints[index]
    parent = subject.Joints[subject.Parents[index]]

    frame1 = joint.Pose
    frame2 = parent.DistalDubinsFrame()

    transformation = frame1.inv() * frame2

    #try to make initial guess right next to each other
    initialPosition = transformation.t[2]
    initialRotation = np.arctan2(transformation.R[1, 0], transformation.R[0, 0])
    initialGuess = [initialPosition,initialRotation]
    initialTree = subject.copyAbbreviatedSelf(ignoreLater, index)
    
    if not initialTree.transformJoint(index, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation), 
                                      propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
        initialTree = subject.copyAbbreviatedSelf(ignoreLater, index)
    
    if not initialTree.transformJoint(index, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation), 
                                      propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
        initialGuess = [0,0]

    initialLoss = objective([0,0])

    #calculate dist
    dist = subject.Links[index].path.length
    dist2 = dist*2
    #bounds = [(-dist*2, dist*2), (-np.pi*2, np.pi*2)]
    positionBound = [min(-dist2, initialPosition - dist2), max(dist2, initialPosition + dist2)]
    angleBound = [-np.pi*2, np.pi*2]
    bounds = [positionBound, angleBound]

    #initial swarm
    global joint_batch_objective_function
    def joint_batch_objective_function(X):
        return np.array([objective(x) for x in X])
    n_particles = 16

    init_pos = np.tile(np.array(initialGuess, dtype='float64'), (n_particles,1))
    init_pos[1] = np.array([0,0])
    #add random noise
    noise = np.zeros_like(init_pos[2:])
    noise[:, 0] = np.random.uniform(-dist2,dist2, n_particles - 2)
    noise[:, 1] = np.random.uniform(-np.pi*2,np.pi*2, n_particles - 2)
    init_pos[2:] += noise
    init_pos[2:, 0] = np.clip(init_pos[2:, 0], -dist2, dist2)
    init_pos[2:, 1] = np.clip(init_pos[2:, 1], -np.pi*2, np.pi*2)

    min_bound = np.array([b[0] for b in bounds])
    max_bound = np.array([b[1] for b in bounds])
    dimensions = 2
    if not min_bound.shape == (dimensions,) or not max_bound.shape == (dimensions,):
        raise ValueError(f"Bounds arrays must be of shape ({dimensions},)")

    optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,dimensions=2,options={'c1':0.6, 'c2':0.7, 'w':0.5},
                                        bounds=(min_bound, max_bound),init_pos=init_pos,ftol=tol)
    minSwarmLoss, minSwarmResult = optimizer.optimize(joint_batch_objective_function, iters=int((maxiter + 1)/2),verbose=False, 
                                                        n_processes=n_particles if parallelize else None)

    nelderMead = minimize(objective, minSwarmResult, method="Nelder-Mead", bounds=bounds, tol=tol, options={
        'maxiter':int(maxiter/2),
        'fatol':tol,
    })
    result = nelderMead.x
    loss = nelderMead.fun

    #print(minSwarmLoss, loss)

    if verbose:
        print(f"Optimized joint {index} in {time.time() - start}s -- Old loss: {initialLoss}, Improved Loss: {loss}")

    which = objective(result, returnWhich=True)

    tree = subject.copyAbbreviatedSelf()
    if which == 1:
        if not tree.transformJoint(index, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1]), 
                                   propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
            raise Exception()
        final_tree, final_loss = tree, loss
    else:
        try:
            if not tree.transformJoint(index, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1]), 
                                       propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
                raise Exception()
            
            tree.Joints[index].reverseZhat()
            if not tree.transformJoint(index, SE3.Trans([0,0,0]), safe=True, relative=True, 
                                       propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False):
                raise Exception()

            final_tree, final_loss = tree, loss
        except:
            tree2 = subject.copyAbbreviatedSelf()
            tree2.Joints[index].reverseZhat()

            if not tree2.transformJoint(index, SE3.Trans([0,0,-result[0]]) @ SE3.Rz(-result[1]), 
                                        safe=True, relative=True, propogate=False, 
                                        recomputeLinkPath=True, recomputeBoundingBall=False):
                raise Exception()
            final_tree, final_loss = tree2, loss

    # Check for collisions if we tried without penalty first
    if original_includeCollisionPenalty and not includeCollisionPenalty:
        # Check for collisions in all configurations
        has_collisions = False
        if configurations is not None:
            for configuration in configurations:
                test_tree = copy.deepcopy(final_tree)
                for i in range(len(test_tree.Joints)):
                    if not isWaypoint(test_tree.Joints[i]):
                        test_tree.setJointState(i, configuration[i])
                        test_tree.Joints[i].recomputeCollisionCapsules()  # Fixed: recompute for each joint
                if test_tree.detectCollisions(specificJointIndex=index, debug=False) > 0:
                    has_collisions = True
                    break
        
        if has_collisions:
            if verbose:
                print(f"Collision detected for joint {index}, retrying with collision penalty")
            # Retry with collision penalty enabled
            return optimizeJointPlacement(subject, index, maxiter, tol, penaltyScale,
                                        childFraction, ignoreLater,
                                        parallelize, verbose, power,
                                        includeCollisionPenalty=True, configurations=configurations,
                                        retryingWithPenalty=True)
        elif verbose:
            print("No collisions detected, using result without penalty")
    
    return final_tree, final_loss

def optimizeWaypointPlacement(subject, index, maxiter, tol, 
                              collisionPenaltyScale, childFraction = 1, 
                              ignoreLater=False, parallelize=False, verbose = True, configurations=None, 
                              includeCollisionPenalty=True, retryingWithPenalty=False):

    # If includeCollisionPenalty is True, first try without it
    original_includeCollisionPenalty = includeCollisionPenalty
    if original_includeCollisionPenalty and not retryingWithPenalty:
        try:
            if verbose:
                print("Trying waypoint optimization without collision penalty...")
            
            # Temporarily disable collision penalty for first attempt
            includeCollisionPenalty = False
            
            # Continue with the optimization logic below (will use includeCollisionPenalty=False)
            # We'll check for collisions after optimization and potentially retry
            
        except Exception as e:
            if verbose:
                print(f"Setup for no-penalty waypoint optimization failed: {str(e)}")
            includeCollisionPenalty = True
    start = time.time()
    initialTree = subject.copyAbbreviatedSelf(ignoreLater, index)

    initialGuess = [0]*6

    parent = initialTree.Joints[initialTree.Parents[index]]
    waypoint = initialTree.Joints[index]

    transform = parent.DistalDubinsFrame() * waypoint.ProximalDubinsFrame().inv()
    initialGuess[0:3] = transform.t
    initialGuess[3:6] = SE3.Rt(transform.R, np.zeros(3)).eul()

    for configuration in configurations:
        copied_subject = copy.deepcopy(subject)
        for i in range(0, len(copied_subject.Joints)):
            if (not isWaypoint(copied_subject.Joints[i])):
                copied_subject.setJointState(i,configuration[i])
                copied_subject.Joints[i].recomputeCollisionCapsules()
        if copied_subject.detectCollisions(debug=True) > 0:
            print(f"Warning: Initial tree in optimizeWaypointPlacement contains collisions in configuration {configuration}.")

    subjects = [copy.deepcopy(subject) for _ in configurations]
    for i in range(0,len(configurations)):
        subjects[i].setConfiguration(configurations[i])

    movedJointIndex = index
    collisionMatrices = subject.buildCollisionMatrices()

    def objective(params):
        tree = subject.copyAbbreviatedSelf(ignoreLater, index)

        if not tree.transformJoint(index, SE3.Trans(params[0:3]) @ SE3.Rz(params[3]) @ SE3.Ry(params[4]) @ SE3.Rz(params[5]),  
                                   propogate=False, safe=True, relative=False, recomputeBoundingBall=False):
            return collisionPenaltyScale * len(subject.Joints) * (len(subject.Children) + 1)
        
        return linkLoss(tree, index, includeCollisionPenalty=True, configurations=configurations, 
                        collisionMatrices=collisionMatrices, movedJointIndex=movedJointIndex) + \
            np.linalg.norm(np.array(params[3:6]) - SE3.Rt(transform.R, np.zeros(3)).eul()) * 10

    if not initialTree.transformJoint(index, SE3.Trans(initialGuess[0:3]) @ SE3.Rz(initialGuess[3]) @ SE3.Ry(initialGuess[4]) @ SE3.Rz(initialGuess[5]),  
                                      propogate=False, safe=True, relative=False, recomputeBoundingBall=False):
        initialGuess = [0]*6

    #print(f"INITAL WAYPOINT GUESS LOSS: {objective(initialGuess)}")

    #initialTree.detectCollisions(debug=True)        
        
    initialLoss = objective([0]*6)

    dist = subject.Links[index].path.length + max(np.amax(np.abs(initialGuess)), np.amax(np.abs(subject.Joints[index].Pose.t)))
    bounds = [(-dist*2, dist*2)]*3 + [(-np.pi*2, np.pi*2)] * 3

    global waypoint_batch_objective_function
    def waypoint_batch_objective_function(X):
        return np.array([objective(x) for x in X])

    n_particles = 24

    init_pos = np.tile(np.array(initialGuess, dtype='float64'), (n_particles,1))
    init_pos[1] = np.array([0]*6, dtype='float64')
    #add random noise
    noise = np.zeros_like(init_pos[2:])
    for i in range(0,3):
        noise[:, i] = np.random.uniform(-dist*2,dist*2, n_particles - 2)
    for i in range(3, 6):
        noise[:, i] = np.random.uniform(-np.pi*2,np.pi*2, n_particles - 2)
    init_pos[2:] += noise
    for i in range(0,3):
        init_pos[2:, i] = np.clip(init_pos[2:, i], -dist*2, dist*2)
    for i in range(3, 6):
        init_pos[2:, i] = np.clip(init_pos[2:, i], -np.pi*2, np.pi*2)
    optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,dimensions=6,options={'c1':0.7, 'c2':0.5, 'w':0.5},bounds=(np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])),init_pos=init_pos,ftol=tol)
    minSwarmLoss, minSwarmResult = optimizer.optimize(waypoint_batch_objective_function, iters=maxiter,verbose=False, n_processes=n_particles if parallelize else None)
    
    tree = subject.copyAbbreviatedSelf()
    if tree.transformJoint(index, SE3.Trans(minSwarmResult[0:3]) @ SE3.Rz(minSwarmResult[3]) @ SE3.Ry(minSwarmResult[4]) @ SE3.Rz(minSwarmResult[5]),  propogate=False, safe=True, relative=False, recomputeBoundingBall=False):
        if verbose:
            print(f"Optimized waypoint {index} in {time.time() - start}s -- Old Loss: {initialLoss}, Improved Loss: {minSwarmLoss}")
        
        # Check for collisions if we tried without penalty first
        if original_includeCollisionPenalty and not includeCollisionPenalty:
            # Check for collisions in all configurations
            has_collisions = False
            if configurations is not None:
                for configuration in configurations:
                    test_tree = copy.deepcopy(tree)
                    for i in range(len(test_tree.Joints)):
                        if not isWaypoint(test_tree.Joints[i]):
                            test_tree.setJointState(i, configuration[i])
                            test_tree.Joints[i].recomputeCollisionCapsules() 
                    if test_tree.detectCollisions(specificJointIndex=index, debug=False) > 0:
                        has_collisions = True
                        break
            
            if has_collisions:
                if verbose:
                    print(f"Collision detected for waypoint {index}, retrying with collision penalty")
                # Retry with collision penalty enabled
                return optimizeWaypointPlacement(subject, index, maxiter, tol, collisionPenaltyScale,
                                               childFraction, ignoreLater,
                                               parallelize, verbose, configurations,
                                               includeCollisionPenalty=True, retryingWithPenalty=True)
            elif verbose:
                print("No collisions detected, using waypoint result without penalty")
        
        return tree, minSwarmResult
    else:
        raise Exception("Optimization failed dramatically")


def optimizeTree(subject, showSteps=False, childFraction=1, guarantee=False, parallelize=False, 
                 evaluate=False, verbose=True, directory=None, resetOnFail=False,
                 traversal="dfs", direction="outward", orderBy="longest", power=2, configurations=None, 
                 repeatTraversal=1):
    if repeatTraversal == "n":
        repeatTraversal = len(subject.Joints)

    times = []
    lengths = []

    if configurations == None:
        num_joints = len(subject.Joints)
        random_config = [0] * num_joints
        
        # Set non-waypoint joints to random states within their range limits
        # Add a tolerance margin to avoid getting too close to limits
        tolerance = 0.02  # 2% margin from limits
        for i in range(num_joints):
            if not isWaypoint(subject.Joints[i]):
                joint = subject.Joints[i]
                min_state, max_state = joint.stateRange()
                range_size = max_state - min_state
                margin = range_size * tolerance
                random_config[i] = np.random.uniform(min_state + margin, max_state - margin)
        
        neutral_config = [0] * num_joints  
        configurations = [neutral_config]

        print("Using random configuration for optimization:")
        print("Joint states and their limits:")
        for i in range(num_joints):
            if not isWaypoint(subject.Joints[i]):
                min_state, max_state = subject.Joints[i].stateRange()
                print(f"  Joint {i}: config = {random_config[i]:.3f} (limits: [{min_state:.3f}, {max_state:.3f}])")

    for i in range(0, len(subject.Joints)):
        subject.Joints[i].recomputeCollisionCapsules()

    if subject.detectCollisions(debug=True) > 0:
        print("Warning: initial tree in optimizeTree contains collisions.")
    if showSteps and isinstance(subject.Joints[0], OrigamiJoint):
        subject.show()

    collisionPenaltyScale = 0
    for i in range(0, len(subject.Joints)):
        if len(subject.Children[i]) > 0:
            continue
        j = i
        length = 0
        while j != 0:
            length += subject.Links[j].path.length ** 2
            j = subject.Parents[j]
        if length > collisionPenaltyScale:
            collisionPenaltyScale = length

    print(f"Collision penalty scale is {collisionPenaltyScale}")

    start = time.time()

    def log(t, idx):
        diff = time.time() - start
        times.append(diff)
        lengths.append(t.totalLength())
        if directory != None:
            t.save(directory + str(diff) + "_" + str(idx), saveDir=False)

    tree = subject.copyAbbreviatedSelf()
    log(tree, -1)

    # create the traversal calls
    treeTraversals = {
        "dfs": partial(dfs, direction=direction, orderBy=orderBy),
        "bfs": partial(bfs, direction=direction, orderBy=orderBy),
        "randomized": partial(randomized, 
                            power=power, 
                            count=len(subject.Joints), 
                            childFraction=childFraction,
                            isWeighted=True)
    }

    print("Doing the optimization:")
    for _ in range(repeatTraversal):
        for index in treeTraversals[traversal](subject):
            print("Optimizing joint index:", index, subject.Joints[index])
            iters = 50
            tolerance = subject.r/10

            if isWaypoint(subject.Joints[index]):
                tree, loss = optimizeWaypointPlacement(tree,index, maxiter=iters, tol=tolerance, 
                                                    collisionPenaltyScale=collisionPenaltyScale, childFraction=childFraction, 
                                                    ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose, 
                                                    configurations=configurations)
            else:
                tree, loss = optimizeJointPlacement(tree,index, maxiter=iters, tol=tolerance, 
                                                    penaltyScale=collisionPenaltyScale, childFraction=childFraction, 
                                                    ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose, 
                                                    configurations=configurations)
            if tree.detectCollisions(specificJointIndex=index, debug=True) > 0:
                print(f"Post-optimization collision detected at joint {index}.")
                print(repr(tree))
                raise Exception("Post-optimization collision detected.")
            log(tree, index)

    if showSteps:
        tree.show()

    print(f"TOTAL OPTIMIZATION TIME: {time.time() - start}")

    if directory != None:
        tree.save(directory + "final", saveDir=False)
    if (evaluate):
        return tree, times, lengths
    
    return tree