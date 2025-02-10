from KinematicTree import *

def optimizeJointPlacement(subject, index, maxiter, tol, collisionError, childParentRatio = 1, ignorePlacement=False, ignoreLater = False, parallelize = False):
    parentIndex = subject.Parents[index]

    selectedIndices = [index] if ignorePlacement else ([index] + subject.Children[index])
    selectedCapsules = subject.selectCollisionCapsules(specificJointIndices=selectedIndices, ignoreLater=ignoreLater)
    def linkLoss(t, link, curveLossFactor = 2):
        #d = (np.abs(t.Links[index].path.theta1 * curveLossFactor) ** 3 + np.abs(t.Links[index].path.theta2 * curveLossFactor) ** 3) + (np.arccos(np.clip((np.trace(t.Joints[index].ProximalDubinsFrame().R.T @ t.Joints[t.Parents[index]].DistalDubinsFrame().R) - 1) / 2, -1.0, 1.0))) * (1/t.Links[index].path.length + 1)
        d = t.Links[index].path.theta1 ** 2 + t.Links[index].path.theta2 ** 2
        childrenLength = 0 if len(t.Children[index]) == 0 else np.mean([t.Links[idx].path.length ** 2 for idx in t.Children[index]]) * childParentRatio
        return t.Links[index].path.length ** 2 + t.detectCollisionsWithCapsules(selectedIndices, selectedCapsules) * collisionError + childrenLength# + d * t.r
    def objective(params, returnWhich = False):
        tree = subject.copyAbbreviatedSelf(ignoreLater, index)

        translation = params[0]
        rotation = params[1]

        transform = SE3.Trans([0,0,translation]) @ SE3.Rz(rotation)

        linkLossUnchanged = collisionError * len(subject.Joints) * (len(subject.Children) + 1)
        linkLossReversedZhat = collisionError * len(subject.Joints) * (len(subject.Children) + 1)
        linkLossReversedParent = collisionError * len(subject.Joints) * (len(subject.Children) + 1)

        try:
            if not tree.transformJoint(index, transform, propogate=ignorePlacement, safe=False, relative=True, recomputeBoundingBall=False):
                raise Exception()
            
            linkLossUnchanged = linkLoss(tree, index)
        except:
            pass
        
        try:
            if linkLossUnchanged == linkLossReversedZhat:
                #only do this if already transformed
                raise Exception()

            tree.Joints[index].reverseZhat()
            if not tree.transformJoint(index, SE3.Trans([0,0,0]), safe=False, relative=True, propogate=ignorePlacement, recomputeLinkPath=True, recomputeBoundingBall=False):
                raise Exception()

            linkLossReversedZhat = linkLoss(tree, index)
        except:
            try:
                tree2 = subject.copyAbbreviatedSelf(ignoreLater, index)
                tree2.Joints[index].reverseZhat()
                if not tree2.transformJoint(index, SE3.Trans([0,0,-translation]) @ SE3.Rz(-rotation), safe=False, relative=True, propogate=ignorePlacement, recomputeLinkPath=True, recomputeBoundingBall=False):
                    raise Exception()

                linkLossReversedZhat = linkLoss(tree2, index)
            except:
                pass

        if returnWhich:
            if linkLossUnchanged <= linkLossReversedParent and linkLossUnchanged <= linkLossReversedZhat:
                return 1
            elif linkLossReversedZhat <= linkLossUnchanged and linkLossReversedZhat <= linkLossReversedParent:
                return 2
            else:
                return 3
        else:
            return min(min(linkLossUnchanged,linkLossReversedZhat), linkLossReversedParent)

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
    initialTree = tree = subject.copyAbbreviatedSelf(ignoreLater, index)
    try:
        if not initialTree.transformJoint(index, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation), propogate=ignorePlacement, safe=False, relative=True, recomputeBoundingBall=False):
            raise Exception()
    except:
        initialGuess = [0,0]

    initialLoss = objective([0,0])

    #calculate dist
    dist = subject.Links[index].path.length
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

    print(f"Optimized joint {index} in {time.time() - start}s -- Old loss: {initialLoss}, Improved Loss: {loss}")

    which = objective(result, returnWhich=True)

    tree = subject.copyAbbreviatedSelf()
    if which == 1:
        if not tree.transformJoint(index, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1]), propogate=False, safe=False, relative=True, recomputeBoundingBall=False):
            raise Exception()
        return tree, loss
    elif which == 2:
        try:
            if not tree.transformJoint(index, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1]), propogate=False, safe=False, relative=True, recomputeBoundingBall=False):
                raise Exception()
            
            tree.Joints[index].reverseZhat()
            if not tree.transformJoint(index, SE3.Trans([0,0,0]), safe=False, relative=True, propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False):
                raise Exception()

            return tree, loss
        except:
            tree2 = subject.copyAbbreviatedSelf()
            tree2.Joints[index].reverseZhat()
            if not tree2.transformJoint(index, SE3.Trans([0,0,-result[0]]) @ SE3.Rz(-result[1]), safe=False, relative=True, propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False):
                raise Exception()
            return tree2, loss
    else:
        tree.Joints[parentIndex].reverseZhat()
        if not tree.transformJoint(parentIndex, SE3.Trans([0,0,0]), safe=False, relative=True, propogate=False, recomputeLinkPath=True, recomputeBoundingBall=False):
            raise Exception()
        if not tree.transformJoint(index, SE3.Trans([0,0,result[0]]) @ SE3.Rz(result[1]), propogate=False, safe=False, relative=True, recomputeBoundingBall=False):
            raise Exception()
        return tree, loss

def optimizeWaypointPlacement(subject, index, maxiter, tol, collisionError, childParentRatio = 1, ignorePlacement=False, ignoreLater=False, parallelize=False):

    #make initial guess as close as possible to previous
    start = time.time()
    initialTree = subject.copyAbbreviatedSelf(ignoreLater, index)

    initialGuess = [0]*6

    parent = initialTree.Joints[initialTree.Parents[index]]
    waypoint = initialTree.Joints[index]

    transform = parent.DistalDubinsFrame() * waypoint.ProximalDubinsFrame().inv()
    initialGuess[0:3] = transform.t
    initialGuess[3:6] = SE3.Rt(transform.R, np.zeros(3)).eul()

    selectedIndices = [index] if ignorePlacement else ([index] + subject.Children[index])
    selectedCapsules = subject.selectCollisionCapsules(specificJointIndices=selectedIndices, ignoreLater=ignoreLater)

    def linkLoss(t, link, curveLossFactor = np.pi):
        d = t.Links[index].path.theta1 ** 2 + t.Links[index].path.theta2 ** 2
        childrenLength = 0 if len(t.Children[index]) == 0 else np.mean([t.Links[idx].path.length ** 2 for idx in t.Children[index]]) * childParentRatio
        return t.Links[index].path.length ** 2  + t.detectCollisionsWithCapsules(selectedIndices, selectedCapsules) * collisionError + childrenLength# + d * t.r

    def objective(params):
        tree = subject.copyAbbreviatedSelf(ignoreLater, index)

        try:
            if not tree.transformJoint(index, SE3.Trans(params[0:3]) @ SE3.Rz(params[3]) @ SE3.Ry(params[4]) @ SE3.Rz(params[5]),  propogate=ignorePlacement, safe=False, relative=False, recomputeBoundingBall=False):
                raise Exception()
        except:
            # tb = traceback.format_exc()
            # print(f"Traceback details:\n{tb}")
            return collisionError * len(subject.Joints) * (len(subject.Children) + 1)
        
        return linkLoss(tree, index) + np.linalg.norm(np.array(params[3:6]) - SE3.Rt(transform.R, np.zeros(3)).eul()) * 10

    try:
        if not initialTree.transformJoint(index, SE3.Trans(initialGuess[0:3]) @ SE3.Rz(initialGuess[3]) @ SE3.Ry(initialGuess[4]) @ SE3.Rz(initialGuess[5]),  propogate=ignorePlacement, safe=False, relative=False, recomputeBoundingBall=False):
            raise Exception("failed")

        #print(f"INITAL WAYPOINT GUESS LOSS: {objective(initialGuess)}")

        #initialTree.detectCollisions(debug=True)
    except Exception as e:
        initialGuess = [0]*6
        
    initialLoss = objective([0]*6)

    dist = subject.Links[index].path.length + max(np.amax(np.abs(initialGuess)), np.amax(np.abs(subject.Joints[index].Pose.t)))
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
    if tree.transformJoint(index, SE3.Trans(minSwarmResult[0:3]) @ SE3.Rz(minSwarmResult[3]) @ SE3.Ry(minSwarmResult[4]) @ SE3.Rz(minSwarmResult[5]),  propogate=False, safe=False, relative=False, recomputeBoundingBall=False):
        print(f"Optimized waypoint {index} in {time.time() - start}s -- Old Loss: {initialLoss}, Improved Loss: {minSwarmLoss}")
        return tree, minSwarmResult
    else:
        raise Exception("Optimization failed dramatically")

def squaredOptimize(subject, showSteps=False, childParentRatio=1, streamline = False, guarantee=False, parallelize=False, evaluate=False):
    times = []
    losses = []

    for i in range(0, len(subject.Joints)):
        subject.Joints[i].recomputeCollisionCapsules()

    if subject.detectCollisions(debug=True) > 0:
        print("Warning: Initial tree contains collisions.")
    if showSteps and isinstance(subject.Joints[0], OrigamiJoint):
        subject.show()

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
            tree, loss = optimizeWaypointPlacement(tree,index, maxiter=iters, tol=tolerance, collisionError=collisionError, childParentRatio=childParentRatio, ignoreLater = (not guarantee), parallelize=parallelize)
        else:
            tree, loss = optimizeJointPlacement(tree,index, maxiter=iters, tol=tolerance, collisionError=collisionError, childParentRatio=childParentRatio, ignoreLater = (not guarantee), parallelize=parallelize)

        if isOptimized[subject.Parents[index]]:
            isOptimized[index] = True
            numOptimized += 1
        else:
            optimizeFromIndex(subject.Parents[index])
        
        times.append(time.time() - start)
        losses.append(optimizationLoss(tree))


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
                    tree2, loss = optimizeWaypointPlacement(tree,order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childParentRatio=childParentRatio, ignorePlacement=True, ignoreLater = (not guarantee), parallelize=parallelize)

                    if tree2.detectCollisions(specificJointIndices=[order[j]], ignoreLater=(not guarantee), debug=True) > 0:
                        raise Exception("Moving all children caused collision.")

                    tree = tree2

                    if optimizeStreak == j:
                        isOptimized[order[j]] = True
                        numOptimized += 1
                        optimizeStreak += 1
                except Exception as e:
                    print(f"COULD NOT IGNORE CHILDREN PLACEMENT {order[j]}: {e}")
                    if j != 0:
                        if isOptimized[tree.Parents[order[j]]]:
                            isOptimized[tree.Parents[order[j]]] = False
                            numOptimized -= 1
                    for idx in order:
                        if isOptimized[idx]:
                            isOptimized[idx] = False
                            numOptimized -= 1

                    tree, loss = optimizeWaypointPlacement(tree,order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childParentRatio=childParentRatio, ignorePlacement=False, ignoreLater = (not guarantee), parallelize=parallelize)
                    print(tree.detectCollisions(specificJointIndices=[order[j]], ignoreLater=(not guarantee), plot=False, debug=True))
                    #tree.show()
                    break
            else:
                try:
                    tree2, loss = optimizeJointPlacement(tree,order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childParentRatio=childParentRatio, ignorePlacement=True, ignoreLater = (not guarantee), parallelize=parallelize)

                    if tree2.detectCollisions(specificJointIndices=[order[j]], ignoreLater=(not guarantee), plot=False, debug=True) > 0:
                        raise Exception("Moving all children caused collision.")
                    
                    tree = tree2

                    if optimizeStreak == j:
                        isOptimized[order[j]] = True
                        numOptimized += 1
                        optimizeStreak += 1
                except Exception as e:
                    print(f"COULD NOT IGNORE CHILDREN PLACEMENT {order[j]}: {e}")
                    if j != 0:
                        if isOptimized[tree.Parents[order[j]]]:
                            isOptimized[tree.Parents[order[j]]] = False
                            numOptimized -= 1
                    for idx in order:
                        if isOptimized[idx]:
                            isOptimized[idx] = False
                            numOptimized -= 1

                    tree, loss = optimizeJointPlacement(tree,order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childParentRatio=childParentRatio, ignorePlacement=True, ignoreLater = (not guarantee), parallelize=parallelize)
                    print(tree.detectCollisions(specificJointIndices=[order[j]], ignoreLater=(not guarantee), plot=False, debug=True))
                    break

            times.append(time.time() - start)
            losses.append(optimizationLoss(tree))
            

        while not isOptimized[i]:
            optimizeFromIndex(i)

        # Another pass, not totally necessary but helps streamline shape, commented out to improve runtime
        if streamline:
            for j in range(0, len(order)):
                iters = 50
                tolerance = subject.r/10
                
                if isWaypoint(subject.Joints[order[j]]):
                    tree, loss = optimizeWaypointPlacement(tree,order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childParentRatio=childParentRatio, ignoreLater = (not guarantee), parallelize=parallelize)
                else:
                    tree, loss = optimizeJointPlacement(tree,order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, childParentRatio=childParentRatio, ignoreLater = (not guarantee), parallelize=parallelize)
                
                times.append(time.time() - start)
                losses.append(optimizationLoss(tree))

        print("CURRENT COLLISIONS")
        if tree.detectCollisions(debug=True) == 0:
            print("NONE")

        if showSteps:
            tree.show()

        print(f"Optimized chain ending at {i} in {time.time() - start2}s \n")

    print(f"TOTAL OPTIMIZATION TIME: {time.time() - start}")
    
    if (evaluate):
        return tree, times, losses
        
    return tree

def linearOptimize(subject, showSteps=False, childParentRatio=1, streamline = False, guarantee=False, parallelize=False, evaluate=False):
    times = []
    losses = []
    start = time.time()

    for i in range(0, len(subject.Joints)):
        subject.Joints[i].recomputeCollisionCapsules()

    if subject.detectCollisions(debug=True) > 0:
        print("Warning: Initial tree contains collisions.")
    if showSteps and isinstance(subject.Joints[0], OrigamiJoint):
        subject.show()

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
    isOptimized = [True] + [False] * (len(subject.Joints) - 1) #isOptimized[i] is True if joint i is optimized
    numOptimized = 1

    def optimizeFromIndex(index):       
        nonlocal tree
        nonlocal numOptimized

        iters = 50
        tolerance = subject.r/10

        if isWaypoint(subject.Joints[index]):
            tree, loss = optimizeWaypointPlacement(tree,index, maxiter=iters, tol=tolerance, collisionError=collisionError, childParentRatio=childParentRatio, ignoreLater = (not guarantee), parallelize=parallelize)
        else:
            tree, loss = optimizeJointPlacement(tree,index, maxiter=iters, tol=tolerance, collisionError=collisionError, childParentRatio=childParentRatio, ignoreLater = (not guarantee), parallelize=parallelize)

        for child in subject.Children[index]:
            optimizeFromIndex(child)
        
        times.append(time.time() - start)
        losses.append(optimizationLoss(tree))
    
    for child in subject.Children[0]:
        optimizeFromIndex(child)

    print(f"TOTAL OPTIMIZATION TIME: {time.time() - start}")

    if (evaluate):
        return tree, times, losses
    
    return tree