from KinematicTree import *

def optimizeJointPlacementMulti(subjects, index, maxiter, tol, collisionError, childFraction = 1, ignorePlacement=False, ignoreLater = False, parallelize = False, verbose=True):
    parentIndex = subjects[0].Parents[index]

    selectedIndices = [index] if ignorePlacement else ([index] + subjects[0].Children[index])
    
    capsuleSelections = [subject.selectCollisionCapsules(specificJointIndices=selectedIndices, ignoreLater=ignoreLater) for subject in subjects]
    def linkLoss(t, link, treeIdx, curveLossFactor = 2):
        selectedCapsules = capsuleSelections[treeIdx]
        #d = (np.abs(t.Links[index].path.theta1 * curveLossFactor) ** 3 + np.abs(t.Links[index].path.theta2 * curveLossFactor) ** 3) + (np.arccos(np.clip((np.trace(t.Joints[index].ProximalDubinsFrame().R.T @ t.Joints[t.Parents[index]].DistalDubinsFrame().R) - 1) / 2, -1.0, 1.0))) * (1/t.Links[index].path.length + 1)
        d = t.Links[index].path.theta1 ** 2 + t.Links[index].path.theta2 ** 2
        childrenLength = 0 if len(t.Children[index]) == 0 else np.mean([t.Links[idx].path.length ** 2 for idx in t.Children[index]]) * childFraction
        return t.Links[index].path.length ** 2 + t.getCollisionError(selectedIndices, selectedCapsules) * collisionError + childrenLength# + d * t.r

    defaultError = collisionError**2 * len(subjects[0].Joints) * (len(subjects[0].Children) + 1)
    def single_tree_objective(tree, treeIdx, params, which):
        translation = params[0]
        rotation = params[1]

        transform = SE3.Trans([0,0,translation]) @ SE3.Rz(rotation)


        if (which == 1):
            #try just moving it
            if tree.transformJoint(index, transform, propogate=ignorePlacement, safe=True, relative=True, recomputeBoundingBall=False):
                return linkLoss(tree, index, treeIdx)
            else:
                return defaultError
        
        #try switching zhat
        tree2 = tree.copyAbbreviatedSelf(ignoreLater, index)

        if tree.transformJoint(index, transform, propogate=ignorePlacement, safe=True, relative=True, recomputeBoundingBall=False):
            tree.Joints[index].reverseZhat()

        if (not tree.transformJoint(index, SE3.Trans([0,0,0]), safe=True, relative=True, propogate=ignorePlacement, recomputeLinkPath=True, recomputeBoundingBall=False)):
            tree2.Joints[index].reverseZhat()
            if tree2.transformJoint(index, SE3.Trans([0,0,-translation]) @ SE3.Rz(-rotation), safe=True, relative=True, propogate=ignorePlacement, recomputeLinkPath=True, recomputeBoundingBall=False):
                return linkLoss(tree2, index, treeIdx)
            else:
                return defaultError
        else:
            return linkLoss(tree, index, treeIdx)
    
    def objective(params, returnWhich = False):
        defaultError = (collisionError**2 * len(subjects[0].Joints) * (len(subjects[0].Children) + 1)) * len(subjects)
        linkLossUnchanged = 0
        linkLossReversedZhat = 0

        for i, subject in enumerate(subjects):
            tree = subject.copyAbbreviatedSelf(ignoreLater, index)
            linkLossUnchanged += single_tree_objective(tree, i, params, 1)
        
        for i, subject in enumerate(subjects):
            tree = subject.copyAbbreviatedSelf(ignoreLater, index)
            linkLossReversedZhat += single_tree_objective(tree, i, params, 2)

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

    returnTrees = []
    for subject in subjects:
        tree = subject.copyAbbreviatedSelf()
        if which == 1:
            if not tree.transformJoint(index, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1]), propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
                raise Exception()
            returnTrees.append(tree)
        else:
            try:
                if not tree.transformJoint(index, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1]), propogate=False, safe=True, relative=True, recomputeBoundingBall=False):
                    raise Exception()
                
                tree.Joints[index].reverseZhat()
                if not tree.transformJoint(index, SE3.Trans([0,0,0]), safe=True, relative=True, propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False):
                    raise Exception()

                returnTrees.append(tree)
            except:
                tree2 = subject.copyAbbreviatedSelf()
                tree2.Joints[index].reverseZhat()
                if not tree2.transformJoint(index, SE3.Trans([0,0,-result[0]]) @ SE3.Rz(-result[1]), safe=True, relative=True, propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False):
                    raise Exception()
                returnTrees.append(tree2)
    
    return returnTrees, loss

def optimizeWaypointPlacementMulti(subjects, index, maxiter, tol, collisionError, childFraction = 1, ignorePlacement=False, ignoreLater=False, parallelize=False, verbose = True):

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
    capsuleSelections = [subject.selectCollisionCapsules(specificJointIndices=selectedIndices, ignoreLater=ignoreLater) for subject in subjects]

    def linkLoss(t, treeIdx, curveLossFactor = np.pi):
        selectedCapsules = capsuleSelections[treeIdx]
        d = t.Links[index].path.theta1 ** 2 + t.Links[index].path.theta2 ** 2
        childrenLength = 0 if len(t.Children[index]) == 0 else np.mean([t.Links[idx].path.length ** 2 for idx in t.Children[index]]) * childFraction
        return t.Links[index].path.length ** 2  + t.getCollisionError(selectedIndices, selectedCapsules) * collisionError + childrenLength# + d * t.r

    defaultError = collisionError**2 * len(subjects[0].Joints) * (len(subjects[0].Children[index]) + 1)
    def single_tree_objective(subject, params, treeIdx):
        tree = subject.copyAbbreviatedSelf(ignoreLater, index)

        if not tree.transformJoint(index, SE3.Trans(params[0:3]) @ SE3.Rz(params[3]) @ SE3.Ry(params[4]) @ SE3.Rz(params[5]),  propogate=ignorePlacement, safe=True, relative=False, recomputeBoundingBall=False):
            return defaultError
        
        return linkLoss(tree, treeIdx) + np.linalg.norm(np.array(params[3:6]) - SE3.Rt(transform.R, np.zeros(3)).eul()) * 10

    def objective(params):
        loss = 0
        for i, subject in enumerate(subjects):
            tree = subject.copyAbbreviatedSelf(ignoreLater, index)
            loss += single_tree_objective(tree, params, i)
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
    
    returnTrees = []
    for subject in subjects:
        tree = subject.copyAbbreviatedSelf()
        if tree.transformJoint(index, SE3.Trans(minSwarmResult[0:3]) @ SE3.Rz(minSwarmResult[3]) @ SE3.Ry(minSwarmResult[4]) @ SE3.Rz(minSwarmResult[5]),  propogate=False, safe=True, relative=False, recomputeBoundingBall=False):
            if verbose:
                print(f"Optimized waypoint {index} in {time.time() - start}s -- Old Loss: {initialLoss}, Improved Loss: {minSwarmLoss}")
            returnTrees.append(tree)
        else:
            raise Exception("Optimization failed dramatically")

    return returnTrees, minSwarmResult

def method1(subject, states = None, showSteps=False, childFraction=1, streamline = False, resetOnFail = True, guarantee=False, parallelize=False, evaluate=False, verbose=True, directory=None):
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

    trees = []
    for state in states:
        tree = copy.deepcopy(subject)
        for i in range(0, len(tree.Joints)):
            if (not isWaypoint(tree.Joints[i])):
                tree.setJointState(i,state[i])
            tree.Joints[i].recomputeCollisionCapsules()
            if tree.detectCollisions(debug=True) > 0:
                print("Warning: Initial tree contains collisions.")
        trees.append(tree)
    
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
        nonlocal trees
        nonlocal numOptimized

        iters = 20
        tolerance = subject.r/10
        if isOptimized[subject.Parents[index]]:
            iters = 50
            tolerance = subject.r/10

        if isWaypoint(subject.Joints[index]):
            trees, loss = optimizeWaypointPlacementMulti(trees,index, maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)
        else:
            trees, loss = optimizeJointPlacementMulti(trees,index, maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)

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
                    trees2, loss = optimizeWaypointPlacementMulti(trees,order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignorePlacement=True, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)

                    for tree in trees2:
                        if tree.detectCollisions(specificJointIndices=[order[j]], ignoreLater=(not guarantee), debug=True) > 0:
                            raise Exception("Moving all children caused collision.")

                    trees = trees2

                    if optimizeStreak == j:
                        isOptimized[order[j]] = True
                        numOptimized += 1
                        optimizeStreak += 1
                except Exception as e:
                    if verbose:
                        print(f"COULD NOT IGNORE CHILDREN PLACEMENT {order[j]}: {e}")
                    if j != 0:
                        if isOptimized[trees[0].Parents[order[j]]]:
                            isOptimized[trees[0].Parents[order[j]]] = False
                            numOptimized -= 1
                    
                    if resetOnFail:
                        print(f"RESET ON FAIL OCCURRED: JOINT {order[j]} TIME: {time.time() - start}")
                        for idx in order:
                            if isOptimized[idx]:
                                isOptimized[idx] = False
                                numOptimized -= 1

                    trees, loss = optimizeWaypointPlacementMulti(trees,order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignorePlacement=False, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)
                    break
            else:
                try:
                    trees2, loss = optimizeJointPlacementMulti(trees,order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignorePlacement=True, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)

                    for tree in trees2:
                        if tree.detectCollisions(specificJointIndices=[order[j]], ignoreLater=(not guarantee), plot=False, debug=True) > 0:
                            raise Exception("Moving all children caused collision.")
                    
                    trees = trees2

                    if optimizeStreak == j:
                        isOptimized[order[j]] = True
                        numOptimized += 1
                        optimizeStreak += 1
                except Exception as e:
                    if verbose:
                        print(f"COULD NOT IGNORE CHILDREN PLACEMENT {order[j]}: {e}")
                    if j != 0:
                        if isOptimized[trees[0].Parents[order[j]]]:
                            isOptimized[trees[0].Parents[order[j]]] = False
                            numOptimized -= 1
                    if resetOnFail:
                        print(f"RESET ON FAIL OCCURRED: JOINT {order[j]} TIME: {time.time() - start}")
                        for idx in order:
                            if isOptimized[idx]:
                                isOptimized[idx] = False
                                numOptimized -= 1

                    trees, loss = optimizeJointPlacementMulti(trees,order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childFraction=childFraction, ignorePlacement=True, ignoreLater = (not guarantee), parallelize=parallelize, verbose=verbose)
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
            for i, tree in enumerate(trees):
                print(f"CURRENT COLLISIONS: {i}")
                if tree.detectCollisions(debug=True) == 0:
                    print("NONE")

        print(f"Optimized chain ending at {i} in {time.time() - start2}s \n")

    print(f"TOTAL OPTIMIZATION TIME: {time.time() - start}")
    
    if directory != None:
        neutralTree = trees[0]
        neutralTree.save(directory + "final", saveDir=False)
    if evaluate:
        return trees, times, lengths
    
    return trees

