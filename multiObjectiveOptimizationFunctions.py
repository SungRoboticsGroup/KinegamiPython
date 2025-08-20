from KinematicTree import *

def optimizeJointPlacementMulti(subject, states, index, maxiter, tol, collisionError, childFraction = 1, ignorePlacement=False, ignoreLater = False, parallelize = False, verbose=True):
    parentIndex = subject.Parents[index]

    
    subjects = [copy.deepcopy(subject) for _ in states]
    for i in range(0,len(states)):
        subjects[i].setConfiguration(states[i])

    selectedIndices = [index] if ignorePlacement else ([index] + subjects[0].Children[index])
    
    capsuleSelections = [subject.selectCollisionCapsules(specificJointIndices=selectedIndices, 
                                                         ignoreLater=ignoreLater) 
                                                         for subject in subjects]
    
    def linkLoss(t, link, curveLossFactor = 2):
        #d = (np.abs(t.Links[index].path.theta1 * curveLossFactor) ** 3 + np.abs(t.Links[index].path.theta2 * curveLossFactor) ** 3) + (np.arccos(np.clip((np.trace(t.Joints[index].ProximalDubinsFrame().R.T @ t.Joints[t.Parents[index]].DistalDubinsFrame().R) - 1) / 2, -1.0, 1.0))) * (1/t.Links[index].path.length + 1)
        d = t.Links[index].path.theta1 ** 2 + t.Links[index].path.theta2 ** 2
        childrenLength = 0 if len(t.Children[index]) == 0 else np.mean([t.Links[idx].path.length ** 2 for idx in t.Children[index]]) * childFraction
        return t.Links[index].path.length ** 2 + childrenLength# + d * t.r

    defaultError = collisionError * len(subjects[0].Joints) * (len(subjects[0].Children) + 1) * len(states)    
    
    stateSubsets = [([states[i][j] for j in range(0,len(states[i])) 
                      if ((j <= index) or (j in subject.Children[index]))] 
                      if ignoreLater else states[i]) for i in range(0,len(states))]

    def objective(params, returnWhich = False):
        tree = subject.copyAbbreviatedSelf(ignoreLater, index)

        linkLossUnchanged = 0
        linkLossReversedZhat = 0

        translation = params[0]
        rotation = params[1]

        transform = SE3.Trans([0,0,translation]) @ SE3.Rz(rotation)

        #try just moving it
        if tree.transformJoint(index, transform, propogate=ignorePlacement, 
                               safe=True, relative=True, recomputeBoundingBall=False):
            
            linkLossUnchanged = linkLoss(tree, index)

            for i in range(0,len(states)):
                selectedCapsules = capsuleSelections[i]
                try:
                    tree.setConfiguration(stateSubsets[i])
                except:
                    print("unwanted")
                    linkLossUnchanged = defaultError
                    break
                linkLossUnchanged += tree.getCollisionError(selectedIndices, selectedCapsules) * \
                    collisionError
        else:
            linkLossUnchanged = defaultError
        
        tree.setConfiguration([0]*len(tree.Joints))
        #try switching zhat

        if linkLossUnchanged < defaultError:
            tree.Joints[index].reverseZhat()

        if (not tree.transformJoint(index, SE3.Trans([0,0,0]), 
                                    safe=True, relative=True, propogate=ignorePlacement, 
                                    recomputeLinkPath=True, recomputeBoundingBall=False)):
            
            tree2 = tree.copyAbbreviatedSelf(ignoreLater, index)
            tree2.Joints[index].reverseZhat()
            
            if tree2.transformJoint(index, SE3.Trans([0,0,-translation]) @ SE3.Rz(-rotation), safe=True, relative=True, propogate=ignorePlacement, recomputeLinkPath=True, recomputeBoundingBall=False):
                linkLossReversedZhat = linkLoss(tree2, index)
            else:
                linkLossReversedZhat = defaultError
        else:
            linkLossReversedZhat = linkLoss(tree, index)
            tree2 = tree

        if (linkLossReversedZhat < defaultError):
            for i in range(0,len(states)):
                selectedCapsules = capsuleSelections[i]
                try:
                    tree2.setConfiguration(stateSubsets[i])
                except:
                    print("unwanted")
                    linkLossReversedZhat = defaultError
                    break
                linkLossReversedZhat += tree2.getCollisionError(selectedIndices, selectedCapsules) * collisionError

        if returnWhich:
            if linkLossUnchanged <= linkLossReversedZhat:
                return 1
            else:
                return 2
        else:
            return min(linkLossUnchanged,linkLossReversedZhat)
            

    start = time.time()

    joint = subjects[0].Joints[index]
    parent = subjects[0].Joints[subjects[0].Parents[index]]

    frame1 = joint.Pose
    frame2 = parent.DistalDubinsFrame()

    transformation = frame1.inv() * frame2

    #try to make initial guess right next to each other
    initialPosition = transformation.t[2]
    initialRotation = np.arctan2(transformation.R[1, 0], transformation.R[0, 0])
    initialGuess = [initialPosition,initialRotation]
    initialTree = subjects[0].copyAbbreviatedSelf(ignoreLater, index)
    
    if not initialTree.transformJoint(index, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation), propogate=ignorePlacement, safe=True, relative=True, recomputeBoundingBall=False):
        initialGuess = [0,0]

    initialLoss = objective([0,0])

    #calculate dist
    dist = subjects[0].Links[index].path.length
    bounds = [(-dist*2, dist*2), (-np.pi*2, np.pi*2)]

    #initial swarm
    global joint_batch_objective_function
    def joint_batch_objective_function(X):
        return np.array([objective(x) for x in X])
    n_particles = 16

    init_pos = np.tile(np.array(initialGuess, dtype='float64'), (n_particles,1))
    init_pos[1] = np.array([0,0])
    #add random noise
    noise = np.zeros_like(init_pos[2:])
    noise[:, 0] = np.random.uniform(-dist*2,dist*2, n_particles - 2)
    noise[:, 1] = np.random.uniform(-np.pi*2,np.pi*2, n_particles - 2)
    init_pos[2:] += noise
    init_pos[2:, 0] = np.clip(init_pos[2:, 0], -dist*2, dist*2)
    init_pos[2:, 1] = np.clip(init_pos[2:, 1], -np.pi*2, np.pi*2)

    optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,dimensions=2,options={'c1':0.6, 'c2':0.7, 'w':0.5},bounds=(np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])),init_pos=init_pos,ftol=tol)
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
        if not tree.transformJoint(index, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1]), propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
            raise Exception()
    else:
        try:
            if not tree.transformJoint(index, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1]), propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
                raise Exception()
            
            tree.Joints[index].reverseZhat()
            if not tree.transformJoint(index, SE3.Trans([0,0,0]), safe=True, relative=True, propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False):
                raise Exception()
        except:
            tree2 = subject.copyAbbreviatedSelf()
            tree2.Joints[index].reverseZhat()
            if not tree2.transformJoint(index, SE3.Trans([0,0,-result[0]]) @ SE3.Rz(-result[1]), safe=True, relative=True, propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False):
                raise Exception()
            tree = tree2
    
    return tree, loss

def optimizeWaypointPlacementMulti(subject, states, index, maxiter, tol, collisionError, childFraction = 1, ignorePlacement=False, ignoreLater=False, parallelize=False, verbose = True):

    subjects = [copy.deepcopy(subject) for _ in states]
    for i in range(0,len(states)):
        subjects[i].setConfiguration(states[i])

    #make initial guess as close as possible to previous
    start = time.time()
    initialTree = subjects[0].copyAbbreviatedSelf(ignoreLater, index)

    initialGuess = [0]*6

    parent = initialTree.Joints[initialTree.Parents[index]]
    waypoint = initialTree.Joints[index]

    transform = parent.DistalDubinsFrame() * waypoint.ProximalDubinsFrame().inv()
    initialGuess[0:3] = transform.t
    initialGuess[3:6] = SE3.Rt(transform.R, np.zeros(3)).eul()

    selectedIndices = [index] if ignorePlacement else ([index] + subjects[0].Children[index])
    capsuleSelections = [subject.selectCollisionCapsules(specificJointIndices=selectedIndices, 
                                                         ignoreLater=ignoreLater) 
                                                         for subject in subjects]

    stateSubsets = [([states[i][j] for j in range(0,len(states[i])) 
                      if ((j <= index) or (j in subject.Children[index]))] 
                      if ignoreLater else states[i]) for i in range(0,len(states))]

    def linkLoss(t, link, curveLossFactor = np.pi):
        d = t.Links[index].path.theta1 ** 2 + t.Links[index].path.theta2 ** 2
        childrenLength = 0 if len(t.Children[index]) == 0 else np.mean([t.Links[idx].path.length ** 2 for idx in t.Children[index]]) * childFraction
        return t.Links[index].path.length ** 2 + childrenLength# + d * t.r

    defaultError = collisionError * len(subject.Joints) * (len(subject.Children) + 1) * len(states)

    def objective(params):
        tree = subject.copyAbbreviatedSelf(ignoreLater, index)
        transform = SE3.Trans(params[0:3]) @ SE3.Rz(params[3]) @ SE3.Ry(params[4]) @ \
            SE3.Rz(params[5])

        if not tree.transformJoint(index, transform,  propogate=ignorePlacement, safe=True, relative=False, recomputeBoundingBall=False):
            return defaultError
            
        loss = linkLoss(tree, index) + np.linalg.norm(np.array(params[3:6]) - SE3.Rt(transform.R, np.zeros(3)).eul()) * 10
        
        for i in range(0, len(states)):
            selectedCapsules = capsuleSelections[i]
            #TODO: Why are there exceptions here, should just transform link
            try:
                tree.setConfiguration(stateSubsets[i])
            except:
                return defaultError
            loss += tree.getCollisionError(selectedIndices, selectedCapsules) * collisionError

        return loss

    if not initialTree.transformJoint(index, SE3.Trans(initialGuess[0:3]) @ SE3.Rz(initialGuess[3]) @ SE3.Ry(initialGuess[4]) @ SE3.Rz(initialGuess[5]),  propogate=ignorePlacement, safe=True, relative=False, recomputeBoundingBall=False):
        initialGuess = [0]*6

    #print(f"INITAL WAYPOINT GUESS LOSS: {objective(initialGuess)}")

    #initialTree.detectCollisions(debug=True)
        
        
    initialLoss = objective([0]*6)

    dist = subjects[0].Links[index].path.length + max(np.amax(np.abs(initialGuess)), np.amax(np.abs(subjects[0].Joints[index].Pose.t)))
    bounds = [(-dist*2, dist*2)]*3 + [(-np.pi*2, np.pi*2)] * 3

    global waypoint_batch_objective_function
    def waypoint_batch_objective_function(X):
        return np.array([objective(x) for x in X])

    n_particles = 24

    init_pos = np.tile(np.array(initialGuess, dtype='float64'), (n_particles,1))
    init_pos[1] = [0]*6
    #add random noise
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
    optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,dimensions=6,options={'c1':0.7, 'c2':0.5, 'w':0.5},bounds=(np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])),init_pos=init_pos,ftol=tol)
    minSwarmLoss, minSwarmResult = optimizer.optimize(waypoint_batch_objective_function, iters=maxiter,verbose=False, n_processes=n_particles if parallelize else None)
    
    tree = subject.copyAbbreviatedSelf()
    if tree.transformJoint(index, SE3.Trans(minSwarmResult[0:3]) @ SE3.Rz(minSwarmResult[3]) @ SE3.Ry(minSwarmResult[4]) @ SE3.Rz(minSwarmResult[5]),  propogate=False, safe=True, relative=False, recomputeBoundingBall=False):
        if verbose:
            print(f"Optimized waypoint {index} in {time.time() - start}s -- Old Loss: {initialLoss}, Improved Loss: {minSwarmLoss}")
    else:
        raise Exception("Optimization failed dramatically")

    return tree, minSwarmResult

def method1(subject, states = None, showSteps=False, childFraction=1, streamline = False, resetOnFail = True, guarantee=True, parallelize=False, evaulate=False, verbose=True, directory=None):
    if states == None:
        states = [[0] * len(subject.Joints)]
    
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

    for state in states:
        tree = copy.deepcopy(subject)
        for i in range(0, len(tree.Joints)):
            if (not isWaypoint(tree.Joints[i])):
                tree.setJointState(i,state[i])
            tree.Joints[i].recomputeCollisionCapsules()
            if tree.detectCollisions(debug=True) > 0:
                print("Warning: Initial tree contains collisions.")
    
    tree = copy.deepcopy(subject)

    times = []
    lengths = []

    start = time.time()

    def log(t, idx):
        diff = time.time() - start
        times.append(diff)
        lengths.append(t.totalLength())

        if directory != None:
            t.save(directory + str(diff) + "_" + str(idx), saveDir=False)

    # log(tree, -1)

    isOptimized = [True] + [False] * (len(subject.Joints) - 1) #isOptimized[i] is True if joint i is optimized
    numOptimized = 1

    def optimizeFromIndex(index):       
        nonlocal tree
        nonlocal numOptimized

        iters = 20
        tolerance = subject.r/10
        if isOptimized[subject.Parents[index]]:
            iters = 50
            tolerance = subject.r/10

        if isWaypoint(subject.Joints[index]):
            tree, loss = optimizeWaypointPlacementMulti(tree, states, index, maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)
        else:
            tree, loss = optimizeJointPlacementMulti(tree, states, index, maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)

        # log(tree, index)

        if isOptimized[subject.Parents[index]]:
            isOptimized[index] = True
            numOptimized += 1
        else:
            optimizeFromIndex(subject.Parents[index])


    while numOptimized < len(subject.Joints):
        i = None
        for j in range(len(subject.Joints) - 1, 0, -1):
            if not isOptimized[j] and len(subject.Children[j]) == 0:
                i = j

        print(f"OPTIMIZING CHAIN ENDING AT {i}:")
        start2 = time.time()
        order = []

        parent = i
        while not isOptimized[parent]:
            order.append(parent)
            parent = subject.Parents[parent]
        order.reverse()

        optimizeStreak = 0

        for j in range(0, len(order)):
            iters = 50
            tolerance = subject.r/10
            
            if isWaypoint(subject.Joints[order[j]]):
                try:
                    tree2, loss = optimizeWaypointPlacementMulti(tree, states, order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignorePlacement=True, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)

                    if tree2.detectCollisions(specificJointIndices=[order[j]], ignoreLater=(not guarantee), debug=True) > 0:
                        raise Exception("Moving all children caused collision.")

                    tree = tree2

                    if optimizeStreak == j:
                        isOptimized[order[j]] = True
                        numOptimized += 1
                        optimizeStreak += 1
                except Exception as e:
                    if verbose:
                        print(f"COULD NOT IGNORE CHILDREN PLACEMENT {order[j]}: {e}")
                    if j != 0:
                        if isOptimized[tree.Parents[order[j]]]:
                            isOptimized[tree.Parents[order[j]]] = False
                            numOptimized -= 1
                    
                    if resetOnFail:
                        print(f"RESET ON FAIL OCCURRED: JOINT {order[j]} TIME: {time.time() - start}")
                        for idx in order:
                            if isOptimized[idx]:
                                isOptimized[idx] = False
                                numOptimized -= 1

                    tree, loss = optimizeWaypointPlacementMulti(tree, states, order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignorePlacement=False, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)
                    break
            else:
                try:
                    tree2, loss = optimizeJointPlacementMulti(tree, states, order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignorePlacement=True, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)

                    if tree2.detectCollisions(specificJointIndices=[order[j]], ignoreLater=(not guarantee), plot=False, debug=True) > 0:
                        raise Exception("Moving all children caused collision.")
                
                    tree = tree2

                    if optimizeStreak == j:
                        isOptimized[order[j]] = True
                        numOptimized += 1
                        optimizeStreak += 1
                except Exception as e:
                    if verbose:
                        print(f"COULD NOT IGNORE CHILDREN PLACEMENT {order[j]}: {e}")
                    if j != 0:
                        if isOptimized[tree.Parents[order[j]]]:
                            isOptimized[tree.Parents[order[j]]] = False
                            numOptimized -= 1
                    if resetOnFail:
                        print(f"RESET ON FAIL OCCURRED: JOINT {order[j]} TIME: {time.time() - start}")
                        for idx in order:
                            if isOptimized[idx]:
                                isOptimized[idx] = False
                                numOptimized -= 1

                    tree, loss = optimizeJointPlacementMulti(tree, states,order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignorePlacement=True, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)
                    break

            # log(tree, order[j])

        while not isOptimized[i]:
            optimizeFromIndex(i)

        # # Another pass, not totally necessary but helps streamline shape, commented out to improve runtime
        # if streamline:
        #     for j in range(0, len(order)):
        #         iters = 50
        #         tolerance = subject.r/10
                
        #         if isWaypoint(subject.Joints[order[j]]):
        #             tree, loss = optimizeWaypointPlacement(tree,order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)
        #         else:
        #             tree, loss = optimizeJointPlacement(tree,order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)
                
        #         log(tree, order[j])
        if verbose:
            print(f"CURRENT COLLISIONS: {i}")
            if tree.detectCollisions(debug=True) == 0:
                print("NONE")

        print(f"Optimized chain ending at {i} in {time.time() - start2}s \n")

    print(f"TOTAL OPTIMIZATION TIME: {time.time() - start}")
    
    #TODO:
    # if directory != None:
    #     tree.save(directory + "final", saveDir=False)
    # if (evaluate):
    #     return tree, times, lengths
    trees = []
    for state in states:
        t = copy.deepcopy(tree)
        for i in range(0, len(t.Joints)):
            if (not isWaypoint(t.Joints[i])):
                t.setJointState(i,state[i])
        trees.append(t)
    return trees

