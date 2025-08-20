from KinematicTree import *
import random 

from collections import deque

from treeTraversals import dfs, bfs



# add a power argument and change the loss and collisionError to be dependent on it
# combine whatever is possible in the optimizeWaypoint/Joint functions and the final optimize function (optimizeTree)

def optimizeJointPlacement(subject, index, maxiter, tol, penaltyScale, childFraction = 1, 
                           ignorePlacement=False, ignoreLater = False, parallelize = False, 
                           verbose=True):
    
    parentIndex = subject.Parents[index]

    # if ignorePlacement is True, optimize only the joint placement and not its children
    # otherwise, optimize the joint placement and its children
    selectedIndices = [index] if ignorePlacement else ([index] + subject.Children[index])
    selectedCapsules = subject.selectCollisionCapsules(specificJointIndices=selectedIndices, 
                                                       ignoreLater=ignoreLater)
    
    def linkLoss(t, link, curveLossFactor = 2):
        #d = (np.abs(t.Links[index].path.theta1 * curveLossFactor) ** 3 + np.abs(t.Links[index].path.theta2 * curveLossFactor) ** 3) + (np.arccos(np.clip((np.trace(t.Joints[index].ProximalDubinsFrame().R.T @ t.Joints[t.Parents[index]].DistalDubinsFrame().R) - 1) / 2, -1.0, 1.0))) * (1/t.Links[index].path.length + 1)
        
        # penalize large curve angles in the Dubins path and favor straighter paths
        # d = t.Links[index].path.theta1 ** 2 + t.Links[index].path.theta2 ** 2

        # childrenLoss = 0 if len(t.Children[index]) == 0 else np.mean([t.Links[idx].path.length ** 2 for idx in t.Children[index]]) * childFraction
        childrenLoss = 0 if len(t.Children[index]) == 0 else np.sum([t.Links[idx].path.length ** 2 
                                                                     for idx in t.Children[index]]) * childFraction
        
        # final cost = own squared link length + collision error + children loss 
        # optional: curve angle penalty 
        return t.Links[index].path.length ** 2 + \
                t.getCollisionError(selectedIndices, selectedCapsules) * penaltyScale + \
                childrenLoss# + d * t.r
    
    def objective(params, returnWhich = False):
        # lightweight copy of the tree 
        tree = subject.copyAbbreviatedSelf(ignoreLater, index)

        translation = params[0]
        rotation = params[1]

        transform = SE3.Trans([0,0,translation]) @ SE3.Rz(rotation)

        # if the path cannot be constructed, return a large penalty
        pathNonExistancePenalty = penaltyScale * (len(subject.Joints) ** 4) * (len(subject.Children[index]) + 1)
        linkLossReversedZhat = pathNonExistancePenalty

        # try just moving it
        # find the cost if valid, otherwise penalize
        if tree.transformJoint(index, transform, propogate=ignorePlacement, safe=True, relative=True, recomputeBoundingBall=False):
            linkLossSameZhat = linkLoss(tree, index)
        else:
            linkLossSameZhat = pathNonExistancePenalty
        
        # try switching zhat
        tree.Joints[index].reverseZhat()

        # if the original zhat didn't work or the new zhat resulted in an invalid transform
        if (linkLossSameZhat == pathNonExistancePenalty or
             not tree.transformJoint(index, SE3(), safe=True, relative=True, propogate=ignorePlacement, recomputeLinkPath=True, recomputeBoundingBall=False)):
            
            # make a new copy of the tree 
            tree2 = subject.copyAbbreviatedSelf(ignoreLater, index)

            # reverse the zhat again and try the inverse transform
            # if it's valid, find the cost, otherwise penalize 
            tree2.Joints[index].reverseZhat()
            if tree2.transformJoint(index, SE3.Trans([0,0,-translation]) @ SE3.Rz(-rotation), safe=True, relative=True, propogate=ignorePlacement, recomputeLinkPath=True, recomputeBoundingBall=False):
                linkLossReversedZhat = linkLoss(tree2, index)
            else:
                linkLossReversedZhat = pathNonExistancePenalty
        
        # if the original zhat worked or the new one resulted in a valid transform,  find the cost
        else:
            linkLossReversedZhat = linkLoss(tree, index)

        # if returnWhich is True, return 1 if the original zhat worked better, otherwise return 2
        # if returnWhich is False, return the minimum of the two costs
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

    # relative transformation joint frame to parent distal frame
    transformation = frame1.inv() * frame2

    # try to make initial guess right next to each other
    # assume that adjusting the joint along its own z-axis is a good initial guess
    initialPosition = transformation.t[2]
    initialRotation = np.arctan2(transformation.R[1, 0], transformation.R[0, 0]) # rotation around z-axis
    initialGuess = [initialPosition,initialRotation] 
    initialTree = subject.copyAbbreviatedSelf(ignoreLater, index)
    
    # try the transform and rebuild the tree if it fails
    # if it fails again, use the initial guess as [0,0] and try to optimize from there
    if not initialTree.transformJoint(index, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation), propogate=ignorePlacement, safe=True, relative=True, recomputeBoundingBall=False):
        initialTree = subject.copyAbbreviatedSelf(ignoreLater, index)
    
    if not initialTree.transformJoint(index, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation), propogate=ignorePlacement, safe=True, relative=True, recomputeBoundingBall=False):
        initialGuess = [0,0]

    # compute the loss at the fallback pose
    initialLoss = objective([0,0])

    # calculate dist and the bounds for the optimization
    dist = subject.Links[index].path.length
    bounds = [(-dist*2, dist*2), (-np.pi*2, np.pi*2)]

    # initial swarm
    global joint_batch_objective_function
    def joint_batch_objective_function(X):
        return np.array([objective(x) for x in X])
    n_particles = 16

    # initialize the swarm with the same initial guess
    init_pos = np.tile(np.array(initialGuess, dtype='float64'), (n_particles,1))
    init_pos[1] = np.array([0,0])
    
    # add random noise
    noise = np.zeros_like(init_pos[2:])
    noise[:, 0] = np.random.uniform(-dist*2,dist*2, n_particles - 2)
    noise[:, 1] = np.random.uniform(-np.pi*2,np.pi*2, n_particles - 2)
    init_pos[2:] += noise

    # keep the values within the bounds
    init_pos[:, 0] = np.clip(init_pos[:, 0], -dist*2, dist*2)
    init_pos[:, 1] = np.clip(init_pos[:, 1], -np.pi*2, np.pi*2)

    # set up the PSO optimizer
    optimizer = ps.single.GlobalBestPSO(n_particles=n_particles, 
                                        dimensions=2, 
                                        options={'c1':0.6, 'c2':0.7, 'w':0.5},
                                        bounds=(np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])),
                                        init_pos=init_pos,
                                        ftol=tol)
    
    # run the optimizer
    minSwarmLoss, minSwarmResult = optimizer.optimize(joint_batch_objective_function, 
                                                      iters=int((maxiter + 1)/2),verbose=False, 
                                                      n_processes=n_particles if parallelize else None)

    # refine the search using Nelder-Mead
    nelderMead = minimize(objective, 
                          minSwarmResult, 
                          method="Nelder-Mead", 
                          bounds=bounds, 
                          tol=tol, 
                          options={'maxiter':int(maxiter/2),'fatol':tol,})
    
    # store the result and the loss
    result = nelderMead.x
    loss = nelderMead.fun

    #print(minSwarmLoss, loss)

    if verbose:
        print(f"Optimized joint {index} in {time.time() - start}s -- Old loss: {initialLoss}, Improved Loss: {loss}")

    # decide which zhat to use based on the objective function
    which = objective(result, returnWhich=True)

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

def optimizeWaypointPlacement(subject, index, maxiter, tol, collisionError, childFraction = 1, 
                              ignorePlacement=False, ignoreLater=False, parallelize=False, verbose = True):

    # make initial guess as close as possible to previous
    start = time.time()
    initialTree = subject.copyAbbreviatedSelf(ignoreLater, index)

    initialGuess = [0]*6

    parent = initialTree.Joints[initialTree.Parents[index]]
    waypoint = initialTree.Joints[index]

    # relative transformation from waypoint proximal frame to parent distal frame
    transform = parent.DistalDubinsFrame() * waypoint.ProximalDubinsFrame().inv()
    initialGuess[0:3] = transform.t
    initialGuess[3:6] = SE3.Rt(transform.R, np.zeros(3)).eul()

    # check only the joint if ignorePlacement is True, otherwise check the joint and its children
    selectedIndices = [index] if ignorePlacement else ([index] + subject.Children[index])
    selectedCapsules = subject.selectCollisionCapsules(specificJointIndices=selectedIndices, ignoreLater=ignoreLater)

    def linkLoss(t, link, curveLossFactor = np.pi):
        # penalize large curve angles in the Dubins path and favor straighter paths
        d = t.Links[index].path.theta1 ** 2 + t.Links[index].path.theta2 ** 2

        # compute the average length of the child links and multiply it by a tunable weight
        childrenLength = 0 if len(t.Children[index]) == 0 \
                        else np.sum([t.Links[idx].path.length ** 2 
                                     for idx in t.Children[index]]) * childFraction
        
        # final cost = own squared link length + collision error + children loss
        # optional: curve angle penalty
        return t.Links[index].path.length ** 2  + t.getCollisionError(selectedIndices, selectedCapsules) * \
                                                collisionError + childrenLength# + d * t.r

    def objective(params):
        tree = subject.copyAbbreviatedSelf(ignoreLater, index)

        # try to apply the full 6D transform to the waypoint
        # xyz translation + yaw, pitch, roll rotation
        # if the transform fails, return a large penalty
        if not tree.transformJoint(index, SE3.Trans(params[0:3]) @ SE3.Rz(params[3]) @ SE3.Ry(params[4]) @ SE3.Rz(params[5]),  propogate=ignorePlacement, safe=True, relative=False, recomputeBoundingBall=False):
            return collisionError * len(subject.Joints) * (len(subject.Children) + 1)
        
        # if the transform is valid, return the link loss and the euler deviation from the original transform
        return linkLoss(tree, index) + np.linalg.norm(np.array(params[3:6]) - SE3.Rt(transform.R, np.zeros(3)).eul()) * 10

    # try to transform the waypoint with the initial guess
    # if it fails, use the initial guess as [0]*6 and try to optimize from there
    if not initialTree.transformJoint(index, SE3.Trans(initialGuess[0:3]) @ SE3.Rz(initialGuess[3]) @ SE3.Ry(initialGuess[4]) @ SE3.Rz(initialGuess[5]),  propogate=ignorePlacement, safe=True, relative=False, recomputeBoundingBall=False):
        initialGuess = [0]*6

    # print(f"INITAL WAYPOINT GUESS LOSS: {objective(initialGuess)}")
    # initialTree.detectCollisions(debug=True)        
        
    # compute the loss at the initial guess (baseline)
    initialLoss = objective([0]*6)

    # calculate the distance and the bounds for the optimization
    dist = subject.Links[index].path.length + max(np.amax(np.abs(initialGuess)), np.amax(np.abs(subject.Joints[index].Pose.t)))
    bounds = [(-dist*2, dist*2)]*3 + [(-np.pi*2, np.pi*2)] * 3

    # initial swarm
    global waypoint_batch_objective_function
    def waypoint_batch_objective_function(X):
        return np.array([objective(x) for x in X])

    n_particles = 24

    # initialize the swarm with the same initial guess
    init_pos = np.tile(np.array(initialGuess, dtype='float64'), (n_particles,1))
    init_pos[1] = [0]*6

    # add random noise
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

    # set up the PSO optimizer
    optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,
                                        dimensions=6,options={'c1':0.7, 'c2':0.5, 'w':0.5},
                                        bounds=(np.array([b[0] for b in bounds]), 
                                        np.array([b[1] for b in bounds])),
                                        init_pos=init_pos,
                                        ftol=tol)

    # run the optimizer
    minSwarmLoss, minSwarmResult = optimizer.optimize(waypoint_batch_objective_function, 
                                                      iters=maxiter,
                                                      verbose=False, 
                                                      n_processes=n_particles if parallelize else None)
    
    # try to apply the optimized transform to the waypoint
    tree = subject.copyAbbreviatedSelf()
    if tree.transformJoint(index, SE3.Trans(minSwarmResult[0:3]) @ SE3.Rz(minSwarmResult[3]) @ SE3.Ry(minSwarmResult[4]) @ SE3.Rz(minSwarmResult[5]),  propogate=False, safe=True, relative=False, recomputeBoundingBall=False):
        if verbose:
            print(f"Optimized waypoint {index} in {time.time() - start}s -- Old Loss: {initialLoss}, Improved Loss: {minSwarmLoss}")
        return tree, minSwarmResult
    else:
        raise Exception("Optimization failed dramatically")

def squaredOptimize(subject, showSteps=False, childFraction=1, streamline = False, resetOnFail = False, guarantee=False, parallelize=False, evaluate=False, verbose = True, directory = None, traversal=dfs, direction="outward", orderBy="longest"):
    times = []
    lengths = []

    # recompute collision capsules for all joints
    for i in range(0, len(subject.Joints)):
        subject.Joints[i].recomputeCollisionCapsules()

    # check for initial collisions
    if subject.detectCollisions(debug=True) > 0:
        print("Warning: Initial tree contains collisions.")
    if showSteps and isinstance(subject.Joints[0], OrigamiJoint):
        subject.show()

    # CHANGE TO COLLISION ERROR WEIGHT
    # compute the collision error, which is the maximum squared length of the path from a joint to the root
    collisionError = 0
    for i in range(0, len(subject.Joints)):
        # check if the joint is a leaf node or not
        if len(subject.Children[i]) > 0:
            continue
        
        # if it is a leaf node, start there and traverse up to the root and find the squared length of the path
        j = i
        length = 0
        while j != 0:
            length += subject.Links[j].path.length ** 2
            j = subject.Parents[j]
        if length > collisionError:
            collisionError = length

    print(f"Collision error is {collisionError}")

    start = time.time()

    # function to log the time and length of the tree at a given index
    def log(t, idx):
        diff = time.time() - start
        times.append(diff)
        lengths.append(t.totalLength())

        if directory != None:
            t.save(directory + str(diff) + "_" + str(idx), saveDir=False)

    # make an abbreviated copy of the subject tree to work with
    tree = subject.copyAbbreviatedSelf()
    log(tree, -1)

    # isOptimized[i] is True if joint i is optimized
    isOptimized = [True] + [False] * (len(subject.Joints) - 1) 
    numOptimized = 1

    def optimizeFromIndex(index):       
        nonlocal tree
        nonlocal numOptimized

        # set the number of iterations and tolerance for the optimization
        iters = 20
        tolerance = subject.r/10

        # if the parent is optimized, use more iterations and a larger tolerance
        if isOptimized[subject.Parents[index]]:
            iters = 50
            tolerance = subject.r/10

        # if the joint is a waypoint, optimize its placement, otherwise optimize the joint placement
        if isWaypoint(subject.Joints[index]):
            tree, loss = optimizeWaypointPlacement(tree,
                                                   index, 
                                                   maxiter=iters, 
                                                   tol=tolerance, 
                                                   collisionError=collisionError, 
                                                   childFraction=childFraction, 
                                                   ignoreLater = (not guarantee), 
                                                   parallelize=parallelize, 
                                                   verbose=verbose)
        else:
            tree, loss = optimizeJointPlacement(tree,
                                                index, 
                                                maxiter=iters, 
                                                tol=tolerance, 
                                                penaltyScale=collisionError, 
                                                childFraction=childFraction, 
                                                ignoreLater = (not guarantee), 
                                                parallelize=parallelize, 
                                                verbose=verbose)

        # log the time and length of the tree at the given index
        log(tree, index)

        # make sure the parent is optimized
        if isOptimized[subject.Parents[index]]:
            isOptimized[index] = True
            numOptimized += 1
        else:
            optimizeFromIndex(subject.Parents[index])


    while numOptimized < len(subject.Joints):
        i = None

        # find the first unoptimized leaf joint
        for j in range(len(subject.Joints) - 1, 0, -1):
            if not isOptimized[j] and len(subject.Children[j]) == 0:
                i = j

        # start optimizing from that joint
        print(f"OPTIMIZING CHAIN ENDING AT {i}:")
        start2 = time.time()
        order = []

        parent = i
        
        # traverse up the tree to find the order of joints to optimize
        while not isOptimized[parent]:
            order.append(parent)
            parent = subject.Parents[parent]
        order.reverse()

        optimizeStreak = 0

        for index in order:
            optimizedThisPass = []

            # set the number of iterations and tolerance for the optimization
            iters = 50
            tolerance = subject.r/10
            
            # if it's a waypoint
            if isWaypoint(subject.Joints[index]):
                try:
                    # try to optimize the waypoint placement while ignoring the children
                    tree2, loss = optimizeWaypointPlacement(tree,
                                                            index, 
                                                            maxiter=iters, 
                                                            tol=tolerance, 
                                                            collisionError=collisionError, 
                                                            childFraction=childFraction, 
                                                            ignorePlacement=True, 
                                                            ignoreLater = (not guarantee), 
                                                            parallelize=parallelize, 
                                                            verbose=verbose)

                    # check for collisions after the optimization
                    if tree2.detectCollisions(specificJointIndices=[index], ignoreLater=(not guarantee), debug=True) > 0:
                        raise Exception("Moving all children caused collision.")

                    tree = tree2

                    # if the optimization was successful, mark the joint as optimized
                    if optimizeStreak == j:
                        isOptimized[index] = True
                        numOptimized += 1
                        optimizeStreak += 1
                        optimizedThisPass.append(index)
                
                # handle the exception if the optimization fails
                except Exception as e:
                    if verbose:
                        print(f"COULD NOT IGNORE CHILDREN PLACEMENT {index}: {e}")

                    # if the joint is not the root, reset the parent to not optimized
                    if j != 0:
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
                    
                    # try to optimize the waypoint placement again while ignoring the children
                    tree, loss = optimizeWaypointPlacement(tree,
                                                           index, 
                                                           maxiter=iters, 
                                                           tol=tolerance, 
                                                           collisionError=collisionError, 
                                                           childFraction=childFraction, 
                                                           ignorePlacement=False, 
                                                           ignoreLater = (not guarantee), 
                                                           parallelize=parallelize, 
                                                           verbose=verbose)
                    
                    # check for collisions after the optimization
                    if verbose:
                        print(tree.detectCollisions(specificJointIndices=[index], ignoreLater=(not guarantee), plot=False, debug=True))
                    #tree.show()
                    break
            
            # if it's a joint
            else:
                try:
                    # try to optimize the joint placement while ignoring the children
                    tree2, loss = optimizeJointPlacement(tree,
                                                         index, 
                                                         maxiter=iters, 
                                                         tol=tolerance, 
                                                         penaltyScale=collisionError, 
                                                         childFraction=childFraction, 
                                                         ignorePlacement=True, 
                                                         ignoreLater = (not guarantee), 
                                                         parallelize=parallelize, 
                                                         verbose=verbose)

                    # check for collisions after the optimization
                    if tree2.detectCollisions(specificJointIndices=[index], ignoreLater=(not guarantee), plot=False, debug=True) > 0:
                        raise Exception("Moving all children caused collision.")
                    
                    tree = tree2

                    # if the optimization was successful, mark the joint as optimized
                    if optimizeStreak == j:
                        isOptimized[index] = True
                        numOptimized += 1
                        optimizeStreak += 1
                        optimizedThisPass.append(index)

                # handle the exception if the optimization fails        
                except Exception as e:
                    if verbose:
                        print(f"COULD NOT IGNORE CHILDREN PLACEMENT {index}: {e}")
                    
                    # if the joint is not the root, reset the parent to not optimized
                    if j != 0:
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
                    tree, loss = optimizeJointPlacement(tree, index, maxiter=iters, tol=tolerance, penaltyScale=collisionError, childFraction=childFraction, ignorePlacement=True, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)
                    if verbose:
                        print(tree.detectCollisions(specificJointIndices=[index], ignoreLater=(not guarantee), plot=False, debug=True))
                    break
            
            # log the time and length of the tree at the given joint
            log(tree, index)
            
        # optimize all joints
        while not isOptimized[i]:
            optimizeFromIndex(i)

        # another pass, not totally necessary but helps streamline shape, commented out to improve runtime
        if streamline:
            for index in traversal(subject, direction=direction, orderBy=orderBy):
                iters = 50
                tolerance = subject.r/10
                
                if isWaypoint(subject.Joints[index]):
                    tree, loss = optimizeWaypointPlacement(tree, index, maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)
                else:
                    tree, loss = optimizeJointPlacement(tree, index, maxiter=iters, tol=tolerance, penaltyScale=collisionError, childFraction=childFraction, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)
                
                log(tree, index)
        
        if verbose:
            print("CURRENT COLLISIONS")
            if tree.detectCollisions(debug=True) == 0:
                print("NONE")

        if showSteps:
            tree.show()

        print(f"Optimized chain ending at {i} in {time.time() - start2}s \n")

    print(f"TOTAL OPTIMIZATION TIME: {time.time() - start}")
    
    # save the final tree if a directory is specified
    if directory != None:
        tree.save(directory + "final", saveDir=False)
    if (evaluate):
        return tree, times, lengths
    
    return tree

def linearOptimize(subject, showSteps=False, childFraction=1, streamline = False, guarantee=False, parallelize=False, evaluate=False, verbose=True, directory=None, traversal=dfs, direction="outward", orderBy="longest"):
    times = []
    lengths = []

    # recompute collision capsules for all joints
    for i in range(0, len(subject.Joints)):
        subject.Joints[i].recomputeCollisionCapsules()

    # check for initial collisions
    if subject.detectCollisions(debug=True) > 0:
        print("Warning: Initial tree contains collisions.")

    # if showSteps is True, show the initial tree
    if showSteps and isinstance(subject.Joints[0], OrigamiJoint):
        subject.show()

    # compute the collision error, which is the maximum squared length of the path from a joint to the root
    collisionError = 0
    for i in range(0, len(subject.Joints)):
        # check if the joint is a leaf node or not
        if len(subject.Children[i]) > 0:
            continue

        # if it is a leaf node, start there and traverse up to the root and find the squared length of the path
        j = i
        length = 0
        while j != 0:
            length += subject.Links[j].path.length ** 2
            j = subject.Parents[j]
        if length > collisionError:
            collisionError = length

    print(f"Collision error is {collisionError}")

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
    # using precomputed traversal orders
    for index in traversal(subject, direction, orderBy):
        # set the number of iterations and tolerance for the optimization
        iters = 50
        tolerance = subject.r/10

        # if it's a waypoint, optimize its placement ignoring the children
        if isWaypoint(subject.Joints[index]):
            tree, loss = optimizeWaypointPlacement(tree, index, maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)
        
        # if it's a joint, optimize the joint placement ignoring the children
        else:
            tree, loss = optimizeJointPlacement(tree, index, maxiter=iters, tol=tolerance, penaltyScale=collisionError, childFraction=childFraction, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)

        # log the time and length of the tree at the given index
        log(tree, index)


    # def optimizeFromIndex(index):       
    #     nonlocal tree
    #     nonlocal numOptimized

    #     # set the number of iterations and tolerance for the optimization
    #     iters = 50
    #     tolerance = subject.r/10

    #     # if it's a waypoint, optimize its placement ignoring the children
    #     if isWaypoint(subject.Joints[index]):
    #         tree, loss = optimizeWaypointPlacement(tree,index, maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)
        
    #     # if it's a joint, optimize the joint placement ignoring the children
    #     else:
    #         tree, loss = optimizeJointPlacement(tree,index, maxiter=iters, tol=tolerance, penaltyScale=collisionError, childFraction=childFraction, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)

    #     # log the time and length of the tree at the given index
    #     log(tree, index)
        
    #     # make sure every child is optimized
    #     for child in subject.Children[index]:
    #         optimizeFromIndex(child)
    
    # # run the optimization for all joints
    # # go down the tree and optimize each joint
    # for child in subject.Children[0]:
    #     optimizeFromIndex(child)

    print(f"TOTAL OPTIMIZATION TIME: {time.time() - start}")

    # save the final tree if a directory is specified
    if directory != None:
        tree.save(directory + "final", saveDir=False)
    if (evaluate):
        return tree, times, lengths
    
    return tree

def perpetualOptimize(subject, iterations, weighted = True, showSteps=False, childFraction=1, parallelize=False, evaluate=False, verbose = True, directory = None):
    # log the time and length of the tree at a given index
    def log(t, idx):
        diff = time.time() - start
        times.append(diff)
        lengths.append(t.totalLength())
        if directory != None:
            t.save(directory + str(diff) + "_" + str(idx), saveDir=False)
            
    times = []
    lengths = []

    # recompute collision capsules for all joints
    for i in range(0, len(subject.Joints)):
        subject.Joints[i].recomputeCollisionCapsules()

    # check for initial collisions
    if subject.detectCollisions(debug=True) > 0:
        print("Warning: Initial tree contains collisions.")

    # if showSteps is True, show the initial tree
    if showSteps and isinstance(subject.Joints[0], OrigamiJoint):
        subject.show()

    # compute the collision error, which is the maximum squared length of the path from a joint to the root
    collisionError = 0
    for i in range(0, len(subject.Joints)):
        if len(subject.Children[i]) > 0:
            continue
        j = i
        length = 0
        while j != 0:
            length += subject.Links[j].path.length ** 2
            j = subject.Parents[j]
        if length > collisionError:
            collisionError = length

    print(f"Collision error is {collisionError}")
 
    start = time.time()

    tree = subject.copyAbbreviatedSelf()
    log(tree, -1)

    subIters = 50
    tolerance = subject.r/10

    for _ in range(0, iterations):
        # randomly select a joint to optimize
        index = np.random.randint(1, len(tree.Joints))

        # if weighted is True, select a joint based on the length of its path
        if weighted:
            weights = weights = [tree.Links[i].path.length for i in range(1,len(tree.Joints))]
            index = random.choices(range(1,len(tree.Joints)), weights=weights, k=1)[0]

        # if it's a waypoint, optimize its placement, otherwise optimize the joint placement
        if isWaypoint(subject.Joints[index]):
            tree, loss = optimizeWaypointPlacement(tree,index, maxiter=subIters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignoreLater = False, parallelize=parallelize, verbose=verbose)
        else:
            tree, loss = optimizeJointPlacement(tree,index, maxiter=subIters, tol=tolerance, penaltyScale=collisionError, childFraction=childFraction, ignoreLater = False, parallelize=parallelize, verbose=verbose)
        
        log(tree, index)
        

    print(f"TOTAL OPTIMIZATION TIME: {time.time() - start}")

    if directory != None:
        tree.save(directory + "final", saveDir=False)
    if (evaluate):
        return tree, times, lengths

    return tree

    