from testqtgraph import *
from makeKinematicTree import *
from KinematicTree import *
from OrigamiJoint import *

import random
from collections import defaultdict, deque
import json
from functools import partial
import threading
import os 
import shutil
import re
import numpy as np
from numpy import array
from spatialmath import SE3

jointCount = 12
probabilityOfBranching = 0.5
sparse = False
cubeSize = 100 if sparse else 10

def linkLoss(tree : KinematicTree, index : int, power : float, childFraction : float, 
             collisionPenaltyScale : float, capsuleSelections=None, 
             includeCollisionPenalty=False, ignorePlacement=False, 
             ignoreOtherBranchCollisions=False, states=None, collisionErrorWeight=1):
    
    assert(childFraction>=0 and power>=0 and index>=0 and collisionPenaltyScale>=0)

    # incoming link
    loss = tree.Links[index].path.length**power 

    # outgoing links
    if len(tree.Children[index]) > 0 and childFraction > 0: 
        loss += np.sum([tree.Links[idx].path.length ** 2 for idx in tree.Children[index]]) * childFraction

    # collision penalty (optional)
    if includeCollisionPenalty:
        for i in range(0, len(states)):
            selectedCapsules = capsuleSelections[i]
            tree.setConfiguration(states[i])
            loss += tree.collisionError(selectedCapsules,
                                        ignoreLater=ignoreOtherBranchCollisions, 
                                        ignorePlacement=ignorePlacement) * collisionErrorWeight
    
    return loss 

def dfs(subject, isOptimized, direction="outward", orderBy="default"):
    stack = [0]
    visited = set()
    order = []
    reverse = False
    branchLengths = {}

    def subtreeLengthToLeaf(node):
        if node in branchLengths:
            return branchLengths[node]
        
        if not subject.Children[node]:
            branchLengths[node] = 0
            return 0
        
        if orderBy == "shortest":
            minDepth = float('inf')
            for child in subject.Children[node]:
                depth = subject.Links[child].path.length + subtreeLengthToLeaf(child)
                minDepth = min(minDepth, depth)
            branchLengths[node] = minDepth

            return minDepth
        
        else:
            maxDepth = 0
            for child in subject.Children[node]:
                depth = subject.Links[child].path.length + subtreeLengthToLeaf(child)
                maxDepth = max(maxDepth, depth)
            branchLengths[node] = maxDepth

            return maxDepth
    
    while stack:
        node = stack.pop()
        order.append(node)
        visited.add(node)

        children = subject.Children[node]

        children = reversed(children)

        if (direction == "outward" and orderBy == "shortest") or (direction == "inward" and orderBy == "longest"):
            reverse = True

        if orderBy in ["longest", "shortest"]:
            children = sorted(children, key=lambda c: subject.Links[c].path.length + subtreeLengthToLeaf(c), reverse=reverse)

        for child in children:
            if child not in visited:
                stack.append(child)

    yield (order if direction == "outward" else list(reversed(order)))

def bfs(subject, isOptimized, direction="outward", orderBy="default"):
    queue = deque([0])
    visited = set()
    order = []
    reverse = False
    level_no = 0

    while queue:
        level_size = len(queue)
        current_level = []

        for _ in range(level_size):
            node = queue.popleft()
            order.append(node)

            for child in subject.Children[node]:
                if child not in visited:
                    current_level.append(child)
                    visited.add(child)

        if (direction == "outward" and orderBy == "longest") or (direction == "inward" and orderBy == "shortest"):
            reverse = True

        if orderBy in ["longest", "shortest"]:
            current_level.sort(key=lambda c: subject.Links[c].path.length, reverse=reverse)

        queue.extend(current_level)

    yield (order if direction == "outward" else list(reversed(order)))

def randomized(subject, isOptimized, power, count, childFraction, isWeighted=True):
    for _ in range(count):
        if isWeighted:
            weights = [linkLoss(subject, i, power, childFraction, 0) for i in range(1, len(subject.Joints))]
            print(weights)
            yield random.choices(range(1, len(subject.Joints)), weights=weights, k=1)[0]
        else:
            yield np.random.randint(1, len(subject.Joints))

def squared(subject, isOptimized, direction="outward", orderBy="default"):
    visited = set()

    while not all(isOptimized):
        dfs_output = next(dfs(subject, direction=direction, orderBy=orderBy))

        leaf = None
        # find an unoptimized leaf node
        for node in reversed(dfs_output):
            if not subject.Children[node] and not isOptimized[node]:
                leaf = node
                break

        if leaf is None:
            return
        
        # get the unoptimized chain
        unoptimized_chain = []
        j = leaf

        while not isOptimized[j]:
            unoptimized_chain.append(j)
            j = subject.Parents[j]
        unoptimized_chain.reverse()

        print(f"\nNext chain to optimize (leaf {leaf}): {unoptimized_chain}")

        # yield one joint at a time
        for node in unoptimized_chain:
            if not isOptimized[node]:
                visited.add(node)
                yield node

def testMultipleArguments(construct):
    directions = ["outward", "inward"]
    orders = ["default", "longest", "shortest"]

    for direction in directions:
        for order in orders:
            print(f"\nDFS (direction={direction} and orderBy={order}):")
            for node in dfs(construct, direction=direction, orderBy=order):
                print(node, end=' ')
            print("\n")

            print(f"BFS (direction={direction} and orderBy={order}):")
            for node in bfs(construct, direction=direction, orderBy=order):
                print(node, end=' ')
            print("\n")

def generateTree(nJoints):
    bounds = (-cubeSize/2, cubeSize/2)
    poses = [ SE3.Rand(xrange=bounds, yrange=bounds, zrange=bounds)
                for _ in range(nJoints) ]

    r = 1
    numSides = 4
    neutralLength = 3

    root = RevoluteJoint(numSides,
                         r,
                         np.pi,
                         poses[0]) if np.random.rand() > 0.5 else PrismaticJoint(numSides,
                                                                                 r,
                                                                                 neutralLength,
                                                                                 3,
                                                                                 np.pi/5,
                                                                                 poses[0])

    specTree = JointSpecificationTree(root)

    for i in range(1, nJoints):
        #parent = np.random.randint(int((i - 1) * (1 - branchingRatio)), i)
        branching = np.random.rand() < probabilityOfBranching and i > 1
        if branching:
            # randomly select a parent from among the non-leaves
            nonLeaves = specTree.nonLeaves()
            parent = nonLeaves[np.random.randint(0, len(nonLeaves))]
        else:
            # randomly select a parent from among the leaves
            leaves = specTree.leaves()
            parent = leaves[np.random.randint(0, len(leaves))]
        
        newJoint = RevoluteJoint(numSides,
                                 r,
                                 np.pi,
                                 poses[i]) if np.random.rand() > 0.5 else PrismaticJoint(numSides,
                                                                                         r,
                                                                                         neutralLength,
                                                                                         3,
                                                                                         np.pi/5,
                                                                                         poses[i])
        
        specTree.addJoint(parent, newJoint)

    initialTree = makeTubularKinematicTree(specTree)
    assert(abs(specTree.totalLengthLowerBound() - initialTree.totalLengthLowerBound()) < 1e-5)
    return initialTree

def testRandomTrees():
    construct = generateTree(jointCount)

    # Open file for writing outputs
    with open("exampleTreeForTraversal.txt", "w") as f:

        # Save the tree structure as repr
        f.write("=== Serialized Tree ===\n")
        f.write(repr(construct) + "\n\n")

        # Print structure with link lengths
        print("\n=== Tree Structure with Link Lengths ===\n", file=f)
        for parent in range(len(construct.Children)):
            children = construct.Children[parent]
            if not children:
                continue

            print(f"- Joint {parent} ->", file=f)
            for child in children:
                link = construct.Links[child]
                path = link.path
                print(f"   ---> Joint {child} | Length: {path.length:.2f}", file=f)
            print("", file=f)

        # DFS and BFS traversal orders
        directions = ["outward", "inward"]
        orders = ["default", "longest", "shortest"]

        print("\n=== Traversal Orders ===", file=f)
        for direction in directions:
            for order in orders:
                print(f"\n\nDFS (direction='{direction}', orderBy='{order}'):", file=f)
                for node in dfs(construct, direction=direction, orderBy=order):
                    print(" ", node, file=f)

                print(f"\nBFS (direction='{direction}', orderBy='{order}'):", file=f)
                for node in bfs(construct, direction=direction, orderBy=order):
                    print(" ", node, file=f)

    print("Tree and traversal outputs saved to exampleTreeForTraversal.txt")



# testTree1 = generateTree(jointCount)
# print("Original Tree:" + repr(testTree1)+ "\n")

# print("\n=======================================\n")

# print(testTree1.__repr__)

# print("\n=======================================\n")

# testTree2 = eval(repr(testTree1))
# print("Reconstructed Tree:" + repr(testTree2) + "\n")

# testRandomTrees()

tree = generateTree(jointCount)

# Initially only root is optimized
isOptimized = [True] + [False] * (len(tree.Joints) - 1)

print("\n=== Tree Structure with Link Lengths ===\n")
for parent in range(len(tree.Children)):
    children = tree.Children[parent]
    if not children:
        continue

    print(f"- Joint {parent} ->")
    for child in children:
        link = tree.Links[child]
        path = link.path
        print(f"   --> Joint {child} | Length: {path.length:.2f}")
    print("")

print("\nOptimization order:")
for node in squared(tree, isOptimized, direction="outward", orderBy="longest"):
    print(f"Optimizing joint {node}...")

    # Simulate success or failure of optimization
    success = random.random() > 0.2
    if success:
        print(f"Optimized joint {node}")
        isOptimized[node] = True
    else:
        print(f"Optimization failed for joint {node}")

print("\nFinal optimization status:")
for i, status in enumerate(isOptimized):
    print(f"  Joint {i}: {'YES' if status else 'NO'}")

