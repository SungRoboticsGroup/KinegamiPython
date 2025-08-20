from KinematicTree import *
from numpy import random
import random 
from treeTraversals import *

def optimizeNodePlacement(subject, index, waypoint, power, maxiter, tol, collisionErrorWeight,
                           penaltyScale=1, childFraction=1, ignorePlacement=False, ignoreOtherBranchCollisions=False, 
                           includeCollisionPenalty=False, parallelize=False, verbose=True, states=None):

    def computeAndMinimizeLoss(params, returnWhich=False):
        subjects = [copy.deepcopy(subject) for _ in states]
        for i in range(0,len(states)):
            subjects[i].setConfiguration(states[i])

        selectedIndices = [index] if ignorePlacement else ([index] + subjects[0].Children[index])
        
        capsuleSelections = [subject.selectCollisionCapsules(specificJointIndices=selectedIndices, 
                                                            ignoreLater=ignoreOtherBranchCollisions) 
                                                            for subject in subjects]
        if waypoint:
            tree = subject.copyAbbreviatedSelf(ignoreOtherBranchCollisions, index)
            transform = SE3.Trans(params[0:3]) @ SE3.Rz(params[3]) @ SE3.Ry(params[4]) @ SE3.Rz(params[5])
        
            if not tree.transformJoint(index, 
                                       transform=transform,
                                       propogate=ignorePlacement, 
                                       safe=True, 
                                       relative=False, 
                                       recomputeBoundingBall=False):
                return collisionErrorWeight * len(subject.Joints) * (len(subject.Children) + 1)
            
            loss = linkLoss(tree, 
                            index, 
                            power,
                            childFraction=childFraction,
                            collisionPenaltyScale=1,
                            capsuleSelections=capsuleSelections,
                            includeCollisionPenalty=includeCollisionPenalty, 
                            ignorePlacement=ignorePlacement,
                            ignoreOtherBranchCollisions=ignoreOtherBranchCollisions,
                            states=states,
                            collisionErrorWeight=collisionErrorWeight) + \
                            np.linalg.norm(np.array(params[3:6]) - \
                            SE3.Rt(transform.R, np.zeros(3)).eul()) * 10
            
            return loss
        
        else:
            # lightweight copy of the tree 
            tree = subject.copyAbbreviatedSelf(ignoreOtherBranchCollisions, index)

            translation = params[0]
            rotation = params[1]

            transform = SE3.Trans([0,0,translation]) @ SE3.Rz(rotation)

            # if the path cannot be constructed, return a large penalty
            pathNonExistancePenalty = penaltyScale * (len(subject.Joints) ** 4) * (len(subject.Children[index]) + 1)
            linkLossReversedZhat = pathNonExistancePenalty

            # try just moving it
            # find the cost if valid, otherwise penalize
            if tree.transformJoint(index, 
                                   transform, 
                                   propogate=ignorePlacement, 
                                   safe=True, 
                                   relative=True, 
                                   recomputeBoundingBall=False):
                linkLossSameZhat = linkLoss(tree, 
                                            index, 
                                            power,
                                            childFraction=childFraction,
                                            collisionPenaltyScale=1,
                                            capsuleSelections=capsuleSelections,
                                            includeCollisionPenalty=includeCollisionPenalty, 
                                            ignorePlacement=ignorePlacement,
                                            ignoreOtherBranchCollisions=ignoreOtherBranchCollisions,
                                            states=states,
                                            collisionErrorWeight=collisionErrorWeight)

            else:
                linkLossSameZhat = pathNonExistancePenalty
            
            # try switching zhat
            if linkLossSameZhat < pathNonExistancePenalty:
                tree.Joints[index].reverseZhat()

            # if the original zhat didn't work or the new zhat resulted in an invalid transform
            if (linkLossSameZhat == pathNonExistancePenalty or
                not tree.transformJoint(index, SE3(), safe=True, relative=True, propogate=ignorePlacement, recomputeLinkPath=True, recomputeBoundingBall=False)):
                
                # make a new copy of the tree 
                tree2 = subject.copyAbbreviatedSelf(ignoreOtherBranchCollisions, index)

                # reverse the zhat again and try the inverse transform
                # if it's valid, find the cost, otherwise penalize 
                tree2.Joints[index].reverseZhat()
                if tree2.transformJoint(index, SE3.Trans([0,0,-translation]) @ SE3.Rz(-rotation), safe=True, relative=True, propogate=ignorePlacement, recomputeLinkPath=True, recomputeBoundingBall=False):
                    linkLossReversedZhat = linkLoss(tree, 
                                                    index, 
                                                    power,
                                                    childFraction=childFraction,
                                                    collisionPenaltyScale=1,
                                                    capsuleSelections=capsuleSelections,
                                                    includeCollisionPenalty=includeCollisionPenalty, 
                                                    ignorePlacement=ignorePlacement,
                                                    ignoreOtherBranchCollisions=ignoreOtherBranchCollisions,
                                                    states=states,
                                                    collisionErrorWeight=collisionErrorWeight)
                else:
                    linkLossReversedZhat = pathNonExistancePenalty
            
            # if the original zhat worked or the new one resulted in a valid transform,  find the cost
            else:
                linkLossReversedZhat = linkLoss(tree2, 
                                                index, 
                                                power,
                                                childFraction=childFraction,
                                                collisionPenaltyScale=1,
                                                capsuleSelections=capsuleSelections,
                                                includeCollisionPenalty=includeCollisionPenalty, 
                                                ignorePlacement=ignorePlacement,
                                                ignoreOtherBranchCollisions=ignoreOtherBranchCollisions,
                                                states=states,
                                                collisionErrorWeight=collisionErrorWeight)

            if (linkLossReversedZhat < pathNonExistancePenalty):
                linkLossReversedZhat = linkLoss(tree2, 
                                                index, 
                                                power,
                                                childFraction=childFraction,
                                                collisionPenaltyScale=1,
                                                capsuleSelections=capsuleSelections,
                                                includeCollisionPenalty=includeCollisionPenalty, 
                                                ignorePlacement=ignorePlacement,
                                                ignoreOtherBranchCollisions=ignoreOtherBranchCollisions,
                                                states=states,
                                                collisionErrorWeight=collisionErrorWeight)

            # if returnWhich is True, return 1 if the original zhat worked better, otherwise return 2
            # if returnWhich is False, return the minimum of the two costs
            if returnWhich:
                if linkLossSameZhat <= linkLossReversedZhat:
                    return 1
                else:
                    return 2
            else:
                return min(linkLossSameZhat,linkLossReversedZhat)

    def addNoiseAndClip(init_pos, dist, n_particles):
        if waypoint:
            noise = np.zeros_like(init_pos[2:])
            for i in range(0,3):
                noise[:, i] = np.random.uniform(-dist*2,dist*2, n_particles - 2)
            for i in range(3, 6):
                noise[:, i] = np.random.uniform(-np.pi*2,np.pi*2, n_particles - 2)
            init_pos[2:] += noise
            
            for i in range(0,3):
                init_pos[2:, i] = np.clip(init_pos[2:, 0], -dist*2, dist*2)
            for i in range(3, 6):
                init_pos[2:, i] = np.clip(init_pos[2:, 1], -np.pi*2, np.pi*2)

        else:
            noise = np.zeros_like(init_pos[2:])
            noise[:, 0] = np.random.uniform(-dist*2,dist*2, n_particles - 2)
            noise[:, 1] = np.random.uniform(-np.pi*2,np.pi*2, n_particles - 2)
            init_pos[2:] += noise

            init_pos[:, 0] = np.clip(init_pos[:, 0], -dist*2, dist*2)
            init_pos[:, 1] = np.clip(init_pos[:, 1], -np.pi*2, np.pi*2)   

    start = time.time()

    # preliminary steps before creating the optimizer
    if waypoint:
        initialTree = subject.copyAbbreviatedSelf(ignoreOtherBranchCollisions, index)

        initialGuess = [0]*6
        parent = initialTree.Joints[initialTree.Parents[index]]
        child = initialTree.Joints[index]

        # relative transformation from waypoint proximal frame to parent distal frame
        transform = parent.DistalDubinsFrame() * child.ProximalDubinsFrame().inv()
        initialGuess[0:3] = transform.t
        initialGuess[3:6] = SE3.Rt(transform.R, np.zeros(3)).eul()

        # try to transform the waypoint with the initial guess
        # if it fails, use the initial guess as [0]*6 and try to optimize from there
        if not initialTree.transformJoint(index, SE3.Trans(initialGuess[0:3]) @ SE3.Rz(initialGuess[3]) @ SE3.Ry(initialGuess[4]) @ SE3.Rz(initialGuess[5]),  propogate=ignorePlacement, safe=True, relative=False, recomputeBoundingBall=False):
            initialGuess = [0]*6

        # print(f"INITAL WAYPOINT GUESS LOSS: {objective(initialGuess)}")
        # initialTree.detectCollisions(debug=True)        
            
        # compute the loss at the initial guess (baseline)
        initialLoss = computeAndMinimizeLoss([0]*6)

        # calculate the distance and the bounds for the optimization
        dist = subject.Links[index].path.length + max(np.amax(np.abs(initialGuess)), np.amax(np.abs(subject.Joints[index].Pose.t)))
        bounds = [(-dist*2, dist*2)]*3 + [(-np.pi*2, np.pi*2)] * 3

    else:
        joint = subject.Joints[index]
        parent = subject.Joints[subject.Parents[index]]

        frame1 = joint.Pose
        frame2 = parent.DistalDubinsFrame()

        # relative transformation joint frame to parent distal frame
        transformation = frame1.inv() * frame2

        # try to make initial guess right next to each other
        # assume that adjusting the joint along its own z-axis is a good initial guess
        initialPosition = transformation.t[2]
        initialRotation = np.arctan2(transformation.R[1, 0], transformation.R[0, 0]) # rotation around z-axis
        initialGuess = [initialPosition,initialRotation] 
        initialTree = subject.copyAbbreviatedSelf(ignoreOtherBranchCollisions, index)

        # try the transform and rebuild the tree if it fails
        # if it fails again, use the initial guess as [0,0] and try to optimize from there
        if not initialTree.transformJoint(index, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation), propogate=ignorePlacement, safe=True, relative=True, recomputeBoundingBall=False):
            initialTree = subject.copyAbbreviatedSelf(ignoreOtherBranchCollisions, index)
        
        if not initialTree.transformJoint(index, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation), propogate=ignorePlacement, safe=True, relative=True, recomputeBoundingBall=False):
            initialGuess = [0,0]

        # compute the loss at the fallback pose
        initialLoss = computeAndMinimizeLoss([0,0])

        # calculate dist and the bounds for the optimization
        dist = subject.Links[index].path.length
        bounds = [(-dist*2, dist*2), (-np.pi*2, np.pi*2)]
    
    # initialize the swarm
    global joint_batch_objective_function
    def joint_batch_objective_function(X):
        return np.array([computeAndMinimizeLoss(x) for x in X])
    
    n_particles = 16 if waypoint else n_particles = 24

    init_pos = np.tile(np.array(initialGuess, dtype='float64'), (n_particles, 1))
    init_pos[1] = [0]*6 if waypoint else init_pos[1] = np.array([0,0]) 
    
    # add noise to the initial positions and clip them to the bounds
    addNoiseAndClip(init_pos=init_pos,
                    dist=dist,
                    n_particles=n_particles)

    # set up the PSO optimizer
    optimizer = ps.single.GlobalBestPSO(n_particles=n_particles, 
                                        dimensions=2, 
                                        options={'c1':0.6, 'c2':0.7, 'w':0.5},
                                        bounds=(np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])),
                                        init_pos=init_pos,
                                        ftol=tol)
    
    # run the optimizer
    minSwarmLoss, minSwarmResult = optimizer.optimize(joint_batch_objective_function, 
                                                      iters=int((maxiter + 1)/2),
                                                      verbose=False, 
                                                      n_processes=n_particles if parallelize else None)

    # final transformations and refinement
    if waypoint:
        # try to apply the optimized transform to the waypoint
        tree = subject.copyAbbreviatedSelf()
        if tree.transformJoint(index, SE3.Trans(minSwarmResult[0:3]) @ SE3.Rz(minSwarmResult[3]) @ SE3.Ry(minSwarmResult[4]) @ SE3.Rz(minSwarmResult[5]),  propogate=False, safe=True, relative=False, recomputeBoundingBall=False):
            if verbose:
                print(f"Optimized waypoint {index} in {time.time() - start}s -- Old Loss: {initialLoss}, Improved Loss: {minSwarmLoss}")
            return tree, minSwarmResult
        else:
            raise Exception("Optimization failed dramatically")

    else:
        # refine the search using Nelder-Mead
        nelderMead = minimize(computeAndMinimizeLoss, 
                            minSwarmResult, 
                            method="Nelder-Mead", 
                            bounds=bounds, 
                            tol=tol, 
                            options={'maxiter':int(maxiter/2),'fatol':tol,})
        
        # store the result and the loss
        result = nelderMead.x
        loss = nelderMead.fun

        if verbose:
            print(f"Optimized joint {index} in {time.time() - start}s -- Old loss: {initialLoss}, Improved Loss: {loss}")

        # decide which zhat to use based on the objective function
        which = computeAndMinimizeLoss(result, returnWhich=True)

        tree = subject.copyAbbreviatedSelf()

        # if it's the original zhat, try to transform the joint and return the tree and loss
        if which == 1:
            if not tree.transformJoint(index, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1]), propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
                raise Exception()
            return tree, loss
        
        # if it's the reversed zhat, revese zhat and try to transform the joint
        else:
            try:
                if not tree.transformJoint(index, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1]), propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
                    raise Exception()
                
                tree.Joints[index].reverseZhat()
                if not tree.transformJoint(index, SE3.Trans([0,0,0]), safe=True, relative=True, propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False):
                    raise Exception()

                return tree, loss
            
            # if the reversed zhat + original transform fails, reverse zhat and try the inverse transform
            except:
                tree2 = subject.copyAbbreviatedSelf()
                tree2.Joints[index].reverseZhat()
                if not tree2.transformJoint(index, SE3.Trans([0,0,-result[0]]) @ SE3.Rz(-result[1]), safe=True, relative=True, propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False):
                    raise Exception()
                return tree2, loss

def optimizeTree(subject, showSteps=False, childFraction=1, streamline = False, ignoreDescendantCollisions=False, 
                 parallelize=False, evaluate=False, verbose=True, directory=None, power = 2, 
                 resetOnFail=True, traversal=dfs, direction="outward", orderBy="longest"):
                 
    times = []
    lengths = []

    # if states is None, initialize it to a list of zeros
    if states == None:
        states = [[0] * len(subject.Joints)]

    for state in states:
        copied_subject = copy.deepcopy(subject)
        for i in range(0, len(copied_subject.Joints)):
            if (not isWaypoint(copied_subject.Joints[i])):
                copied_subject.setJointState(i,state[i])
            copied_subject.Joints[i].recomputeCollisionCapsules()
            if copied_subject.detectCollisions(debug=True) > 0:
                print("Warning: Initial tree contains collisions.")

    # recompute collision capsules for all joints
    for i in range(0, len(subject.Joints)):
        subject.Joints[i].recomputeCollisionCapsules()

    # find the collision error weight with a power
    collisionErrorWeight = np.sum(link.length ** power for link in subject.Links)
    print(f"Collision error is: {collisionErrorWeight}\n")

    # check for initial collisions
    if subject.detectCollisions(debug=True) > 0:
        print("Warning: Initial tree contains collisions.")

    # if showSteps is True, show the initial tree
    if showSteps and isinstance(subject.Joints[0], OrigamiJoint):
        subject.show()

    start = time.time()

    # function to log the time and length of the tree at a given index
    def log(t, idx):
        diff = time.time() - start
        times.append(diff)
        lengths.append(t.totalLength())
        if directory != None:
            t.save(directory + str(diff) + "_" + str(idx), saveDir=False)

    tree = subject.copyAbbreviatedSelf()
    log(tree, -1)

    # isOptimized[i] is True if joint i is optimized
    isOptimized = [True] + [False] * (len(subject.Joints) - 1)
    numOptimized = 1

    # function to traverse the tree based on the specified traversal method
    for index in traversal(subject, isOptimized, direction, orderBy):
        iters = 50
        tolerance = subject.r/10

        optimizedThisPass = []

        # optimize the node
        try:
            tree, loss = optimizeNodePlacement(tree, 
                                                index,
                                                waypoint = isWaypoint(subject.Joints[index]),
                                                power = 2,
                                                maxiter = iters, 
                                                tol = tolerance, 
                                                collisionErrorWeight = collisionErrorWeight, 
                                                childFraction = childFraction, 
                                                ignoreOtherBranchCollisions = (not ignoreDescendantCollisions), 
                                                includeCollisionPenalty = True,
                                                parallelize = parallelize, 
                                                verbose = verbose,
                                                states = states)
            
            # check for collisions after the optimization
            if tree.detectCollisions(specificJointIndices=[index], 
                                    ignoreOtherBranchCollisions=(not ignoreDescendantCollisions), 
                                    debug=True) > 0:
                raise Exception("Moving all children caused collision.")
            
            # if the optimization was successful, mark the joint as optimized
            if isOptimized[subject.Parents[index]]:
                isOptimized[index] = True
                numOptimized += 1
                optimizeStreak += 1
                optimizedThisPass.append(index)

        # handle the exception if the optimization fails        
        except Exception as e:
            if verbose:
                print(f"COULD NOT IGNORE CHILDREN PLACEMENT {index}: {e}")
            
            # if the joint is not the root, reset the parent to not optimized
            if index != 0:
                if isOptimized[tree.Parents[index]]:
                    isOptimized[tree.Parents[index]] = False
                    numOptimized -= 1
            
            # if resetOnFail is True, reset all optimized joints to not optimized
            if resetOnFail:
                print(f"RESET ON FAIL OCCURRED: JOINT {index} TIME: {time.time() - start}")
                for idx in optimizedThisPass:
                    if isOptimized[idx]:
                        isOptimized[idx] = False
                        numOptimized -= 1

            # try to optimize the joint placement again while ignoring the children
            tree, loss = optimizeNodePlacement(tree, 
                                                index, 
                                                maxiter = iters, 
                                                tol = tolerance, 
                                                penaltyScale = collisionErrorWeight, 
                                                childFraction = childFraction, 
                                                ignorePlacement = True, 
                                                ignoreOtherBranchCollisions = (not ignoreDescendantCollisions), 
                                                includeCollisionPenalty = True,
                                                parallelize = parallelize, 
                                                verbose = verbose,
                                                states = states)
            if verbose:
                print(tree.detectCollisions(specificJointIndices=[index], 
                                            ignoreOtherBranchCollisions=(not ignoreDescendantCollisions), 
                                            plot=False, 
                                            debug=True))
            break
        
        # log the time and length of the tree at the given index
        log(tree, index)

    if verbose:
        print("CURRENT COLLISIONS")
        if tree.detectCollisions(debug=True) == 0:
            print("NONE")

    if showSteps:
        tree.show()

    print(f"TOTAL OPTIMIZATION TIME: {time.time() - start}")

    # save the final tree if a directory is specified
    if directory != None:
        tree.save(directory + "final", saveDir=False)
    if (evaluate):
        return tree, times, lengths
    
    return tree