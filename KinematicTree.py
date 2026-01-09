# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 23:13:27 2023

@author: Daniel Feshbach
"""
from ast import Raise
import Joint
from Joint import *
from OrigamiJoint import *
import PathCSC
from PathCSC import *
import scipy
from scipy.optimize import NonlinearConstraint, minimize
import queue
import TubularPattern
from TubularPattern import *
from LinkCSC import LinkCSC
from PrintedJoint import *
import os
import time
from typing import Generic, TypeVar
from functools import partial
from geometryHelpers import *
import pyswarms as ps
import logging
import collections
import traceback
import style

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('pyswarms')
logger.setLevel(logging.DEBUG)

J = TypeVar("J", bound=Joint)

class KinematicTree(Generic[J]):
    """
    Nodes are Joint objects
    Edges are Dubins linkages from parent distal frame to child proximal frame    
    Attributes (GLOBAL COORDINATES):
        r               tubular radius
        Joints          array of Joint objects (nodes)
        Parents         array of parent indices in self.Joints
        Paths           array of CSC Dubins paths to each joint from its parent
        boundingBall    ball bounding all proximal, central, and distal origins
        Children        array of arrays of child indices of each joint
    """
    def __init__(self, root : J, maxAnglePerElbow : float = np.pi/2, joints : list[Joint] = None,
                    links : list[LinkCSC] = None, parents : list[int] = None, children : list[list[int]] = None, boundingBall : Ball = None):
            self.r = root.r

            saveParamsAreNone = np.array([joints is None, links is None, parents is None, children is None, boundingBall is None])
            assert(np.all(saveParamsAreNone) or np.all(np.logical_not(saveParamsAreNone)))

            try:
                self.numSides = root.numSides
            except:
                self.numSides = 4
            if joints:
                self.Joints = joints
            else:
                self.Joints = [root]
            if parents:
                self.Parents = parents
            else:
                self.Parents = [-1]     # root has no parent
            if links:
                self.Links = links
            else:
                self.Links = [LinkCSC(self.r, root.ProximalDubinsFrame(),
                                        root.ProximalDubinsFrame(),
                                        maxAnglePerElbow)]
            assert(maxAnglePerElbow >= 0 and maxAnglePerElbow <= np.pi)
            self.maxAnglePerElbow = maxAnglePerElbow

            if boundingBall:
                self.boundingBall = boundingBall
            else:
                self.boundingBall = root.boundingBall()
                if self.boundingBall.r < self.r:
                    self.boundingBall = Ball(root.Pose.t, self.r)
                
            if children:
                self.Children = children
            else:
                self.Children = [[]]

    def __repr__(self):
        numpy_precision = np.get_printoptions()['precision']
        if numpy_precision < 16:
            np.set_printoptions(precision=16)
        output = (
            f"KinematicTree(root={repr(self.Joints[0])}, "
            f"maxAnglePerElbow={repr(self.maxAnglePerElbow)}, "
            f"joints={repr(self.Joints)}, "
            f"links={repr(self.Links)}, "
            f"parents={repr(self.Parents)}, "
            f"children={repr(self.Children)}, "
            f"boundingBall={repr(self.boundingBall)})"
        )
        np.set_printoptions(precision=numpy_precision)
        return output
    
    def dataDeepCopy(self):
        return copy.deepcopy([self.r, self.Joints, self.Parents, 
                              self.Links, self.maxAnglePerElbow, 
                              self.boundingBall, self.Children])
    
    def setTo(self, data : list):
        self.r, self.Joints, self.Parents, \
            self.Links, self.maxAnglePerElbow, \
            self.boundingBall, self.Children = data

    """
    Returns the new Joint's index. 
    relative - boolean: is newJoint input in parent-relative coordinates (True)
                        or global coordiantes (False)?
    fixedPosition - boolean: should the joint be located
                    (True) exactly at its given position, or
                    (False) somewhere kinematically equivalent (i.e., on the
                            same z axis) chosen by the placement algorithm
    fixedOrientation - boolean: should the joint be oriented
                    (True) exactly at its given orientation, or
                    (False) something kinematically equivalent (i.e., with the
                            same z axis) with x axis constructed as the common 
                            normal from the parent
    safe - boolean: if this is True, allow the algorithm 
            to insert intermediate waypoints to route the path from the 
            parent to guarantee it avoids local self-intersection (i.e., run 
            Algorithm 9 from the Kinegami paper instead of Algorithm 8).
            Not compatible with fixedPosition or fixedOrientation.
    endPlane - Plane: defaults to None, but if this is specified and 
            fixedPosition is False, the algorithm will place the new joint such
            that its whole bounding sphere is >= 4r from this plane.
    """
    def addJoint(self, parentIndex : int, newJoint : J, 
                 relative : bool = True, fixedPosition : bool = False, 
                 fixedOrientation : bool = False, 
                 safe : bool = True, endPlane : Plane = None,
                 chooseXhatToMinPath : bool = False) -> int:
        
        if isinstance(newJoint, OrigamiJoint):
            if newJoint.r != self.r:
                raise ValueError("ERROR: newJoint.r != self.r")
            if newJoint.numSides != self.Joints[parentIndex].numSides:
                raise ValueError("ERROR: newJoint.numSides != self.Joints[parentIndex].numSides")
        if safe and fixedPosition:
            raise ValueError("ERROR: trying to call addJoint with \
                safe and fixedPosition both True")
        if safe and fixedOrientation:
            raise ValueError("ERROR: trying to call addJoint with \
                safe and fixedOrientation both True")
        
        newJoint = copy.deepcopy(newJoint)
        parent = self.Joints[parentIndex]
        if relative:
            newJoint.transformPoseIntoFrame(parent.Pose)

        if safe: # Algorithm 9 from [Chen et al. 2023]
            jointsToAdd = placeJointAndWayPoints(newJoint, parent,
                                                          self.boundingBall)
            i = parentIndex
            for joint in jointsToAdd:
                i = self.addJoint(parentIndex=i, newJoint=joint, 
                                  safe=False, fixedPosition=True, 
                                  fixedOrientation=True, relative=False)
            return i
        
        if not fixedPosition: 
            if endPlane is None: #Algorithm 8 from [Chen et al. 2023]
                newJoint = moveJointNearNeighborBut4rFromBall(newJoint, parent,
                                                          self.boundingBall)
            else: # Tree algorithm for WAFR
                newJoint = moveJointNearNeighborBut4rPastPlane(newJoint, parent,
                                                            endPlane)

        if not fixedOrientation:
            if chooseXhatToMinPath:
                if not newJoint.pathIndex() == 2:
                    def objective(angleToRotateAboutZ):
                        newJointCopy = copy.deepcopy(newJoint)
                        newJointCopy.applyTransformationToPose(SE3.Rz(angleToRotateAboutZ))
                        endDubinsFrame = newJointCopy.ProximalDubinsFrame()
                        startDubinsFrame = parent.DistalDubinsFrame()
                        path = shortestCSC(newJoint.r, startDubinsFrame.t, startDubinsFrame.R[:,0], 
                                        endDubinsFrame.t, endDubinsFrame.R[:,0])
                        if path is None or norm(path.error) > 0.001 * path.r:
                            return np.inf
                        else:
                            return path.length
                        
                    result = minimize(objective, 0)
                    newJoint.applyTransformationToPose(SE3.Rz(result.x[0]))
            elif endPlane is None:
                xhat = commonNormal(parent.Pose.t, parent.Pose.R[:,2],
                                    newJoint.Pose.t, newJoint.Pose.R[:,2],
                                    undefined=newJoint.Pose.R[:,0])
                newJoint.setXhatAboutZhat(xhat)
                outwardDirection = newJoint.Pose.t - parent.Pose.t
                if np.dot(newJoint.pathDirection(), outwardDirection) < 0:
                    newJoint.reversePathDirection()
            else: # make xhat point as forward as possible
                def newXhat(angleToRotateAboutZ):
                    return (SE3.Rz(angleToRotateAboutZ) * newJoint.Pose.R[:,0]).flatten()
                def objective(angleToRotateAboutZ):
                    return -np.dot(newXhat(angleToRotateAboutZ), endPlane.nhat)
                result = minimize(objective, 0)
                # xhat = newXhat(result.x[0])
                newJoint.applyTransformationToPose(SE3.Rz(result.x[0]))



        newLink = LinkCSC(self.r, parent.DistalDubinsFrame(), 
                                newJoint.ProximalDubinsFrame(),
                                self.maxAnglePerElbow)
        if newLink is None:
            print("WARNING: no valid path found to newJoint, chain not changed.")
            return None
        
        self.boundingBall = minBoundingBall(self.boundingBall, 
                                            newLink.elbow2BoundingBall)
        self.boundingBall = minBoundingBall(self.boundingBall, 
                                            newLink.elbow1BoundingBall)
        self.boundingBall = minBoundingBall(self.boundingBall, 
                                            newJoint.boundingBall())
        
        newIndex = len(self.Joints)
        self.Joints.append(newJoint)
        self.Children[parentIndex].append(newIndex)
        self.Children.append([])
        self.Parents.append(parentIndex)
        self.Links.append(newLink)
        
        #set twist angle for newly added joint
        if (isinstance(newJoint, PrintedJoint)):
            jointProximalFrame = newJoint.ProximalDubinsFrame()
            prevDistalFrame = self.Joints[parentIndex].DistalDubinsFrame()
            twistAngle = signedAngle(jointProximalFrame.R[:,1], prevDistalFrame.R[:,1], jointProximalFrame.R[:,0])
            newJoint.setTwistAngle(twistAngle)

        return newIndex
    
    
    def recomputeBoundingBall(self):
        self.boundingBall = self.Joints[0].boundingBall()
        for joint in self.Joints[1:]:
            if joint is None:
                continue
            self.boundingBall = minBoundingBall(self.boundingBall,
                                                joint.boundingBall())
        for link in self.Links:
            if link is None:
                continue
            self.boundingBall = minBoundingBall(self.boundingBall,
                                                link.elbow1BoundingBall)
            self.boundingBall = minBoundingBall(self.boundingBall,
                                                link.elbow2BoundingBall)
            
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                  proximalColor='c', centerColor='m', distalColor='y',
                  showJointSurface=True, jointColor=jointColorDefault, jointEdgeColor=jointEdgeColorDefault,
                  jointAxisScale=jointAxisScaleDefault, showJointPoses=True,
                  linkColor=linkColorDefault, surfaceOpacity=surfaceOpacityDefault, showLinkSurface=True, 
                  showLinkPoses=False, showLinkPath=True, pathColor=pathColorDefault,
                  showPathCircles=False, sphereColor=sphereColorDefault,
                  showSpheres=False, showGlobalFrame=False, globalAxisScale=globalAxisScaleDefault, 
                  showCollisionBoxes=False, showSpecificCapsules = ([],[]), plotPoint=None, addCapsules = []):
        xyzHandles = []
        abcHandles = []
        
        if showGlobalFrame:
            handles = addPosesToPlot(SE3(), ax, globalAxisScale, xColor, yColor, zColor)
            if not handles is None:
                xyzHandles.append(handles)
        
        for joint in self.Joints:
            if joint is None:
                continue
            handles = joint.addToPlot(ax, xColor, yColor, zColor, 
                                    proximalColor, centerColor, distalColor, 
                                    sphereColor=sphereColor, showSphere=showSpheres, 
                                    surfaceColor=jointColor, edgeColor=jointEdgeColor, surfaceOpacity=surfaceOpacity,
                                    showSurface=showJointSurface, axisScale=jointAxisScale,
                                    showPoses=showJointPoses)
            if not handles is None:
                xyzHandles.append(handles)


        for i, link in enumerate(self.Links):
            if link is None:
                continue

            color = linkColor[i] if isinstance(linkColor, list) or isinstance(linkColor, np.ndarray) else linkColor
            handles = link.addToPlot(ax, color=color, 
                                   alpha=surfaceOpacity, 
                                   showPath=showLinkPath, 
                                   pathColor=pathColor,
                                   showPathCircles=showPathCircles, 
                                   showFrames=showLinkPoses,
                                   showBoundary=showLinkSurface)
            if showLinkPoses:
                for elbowHandles in handles:
                    abcHandles.append(elbowHandles)

        
        for capsule in addCapsules:
            capsule.addToPlot(ax)

        for jointIndex, capsuleIndex in showSpecificCapsules[0]:
            self.Joints[jointIndex].collisionCapsules[capsuleIndex].addToPlot(ax)

        for linkIndex, capsuleIndex in showSpecificCapsules[1]:
            self.Links[linkIndex].collisionCapsules[capsuleIndex].addToPlot(ax)
        
        # for jointIndex, capsuleIndex in showCollisionBoxes[0]:
        #     self.Joints[jointIndex].collisionCapsules[capsuleIndex].box.addToPlot(ax)

        # for linkIndex, capsuleIndex in showCollisionBoxes[1]:
        #     self.Links[linkIndex].collisionCapsules[capsuleIndex].box.addToPlot(ax)
        
        if showSpheres:
            if self.nestedBallsRelativeToJoints is None:
                self.boundingBall.addToPlot(ax, color=sphereColor, 
                                        alpha=0.05, frame=True)
            else:
                for i, ball in enumerate(self.nestedBallsGlobal()):
                    # ball is expressed relative to joint i, but we need it in global coordinates
                    ball.addToPlot(ax, color=sphereColor, alpha=0.05, frame=True)
        
        if not plotPoint is None:
            ax.scatter(plotPoint[0], plotPoint[1], plotPoint[2], color='black', s=50)
        return np.array(xyzHandles), np.array(abcHandles)
    
    def copyAbbreviatedSelf(self, isolate=False, isolateJoint = 0):
        # tree = KinematicTree(Waypoint(self.numSides, self.r, SE3()), self.maxAnglePerElbow)
        # # for i in range(1, len(self.Joints)):
        # #     tree.addJoint(self.Parents[i], self.Joints[i].copy(), relative=False, fixedPosition=True, fixedOrientation=True, safe = False)
        # tree.Joints = [copy.deepcopy(joint) for joint in self.Joints]
        # tree.Links = [copy.deepcopy(link) for link in self.Links]
        # tree.Parents = self.Parents
        # tree.Children = self.Children
        # return tree
        try:
            newTree = KinematicTree(copy.deepcopy(self.Joints[self.Parents[self.Parents[isolateJoint]]]), self.maxAnglePerElbow)
        except:
            newTree = KinematicTree(copy.deepcopy(self.Joints[0]), self.maxAnglePerElbow)

        if isolate:

            numChildren = len(self.Children[isolateJoint])
            childIdx = 0
            
            i = 1
            while (childIdx < numChildren or i <= isolateJoint):
                if i == self.Parents[isolateJoint]:
                    newTree.Joints.append(copy.deepcopy(self.Joints[i]))
                    newTree.Children.append([isolateJoint])
                    newTree.Parents.append(None)
                    newTree.Links.append(copy.deepcopy(self.Links[i]))
                elif isolateJoint == i:
                    newTree.Joints.append(copy.deepcopy(self.Joints[i]))
                    newTree.Children.append(self.Children[isolateJoint].copy())
                    newTree.Parents.append(self.Parents[i])
                    newTree.Links.append(copy.deepcopy(self.Links[i]))
                elif numChildren > 0 and self.Children[isolateJoint][childIdx] == i:
                    newTree.Joints.append(copy.deepcopy(self.Joints[i]))
                    newTree.Links.append(copy.deepcopy(self.Links[i]))
                    newTree.Parents.append(self.Parents[i])
                    newTree.Children.append([])
                    childIdx += 1
                else:
                    newTree.Links.append(None)
                    newTree.Joints.append(None)
                    newTree.Parents.append(None)
                    newTree.Children.append(None)
                i += 1
                                
        else:
            newTree.Joints = [copy.deepcopy(joint) for joint in self.Joints]
            newTree.Links = [copy.deepcopy(link) for link in self.Links]
            newTree.Parents = self.Parents.copy()
            newTree.Children = self.Children.copy()

        return newTree

    def detectCollisions(self, specificJointIndices = None, plot=False, includeEnds=False, debug=False, ignoreLater=False, ignoreWaypoints=True):
        toCheck = list(range(len(self.Joints))) if specificJointIndices is None else specificJointIndices
        numCollisions = 0
        # for index in toCheck:
        #     capsules = self.selectCollisionCapsules(specificJointIndices=[index], ignoreLater=ignoreLater)
        #     numCollisions += self.detectCollisionsWithCapsules([index], capsules, debug=debug)

        collisionPairDict = self.buildCollisionPairDictionary()
        numCollisions = self.detectCollisionsWithPairs(toCheck, collisionPairDict, debug=debug)
        
        return numCollisions

    def selectCollisionCapsules(self, specificJointIndices = None, ignoreLater = False, ignoreWaypoints=True):
        allCapsules = [[],[]]
        EPSILON = 0.001

        others = [x for x in list(range(0,len(self.Joints))) if not x in specificJointIndices]
        if ignoreLater:
            latest = max(specificJointIndices)
            others = [x for x in others if x <= latest]

        def posesAreSame(pose1, pose2):
            return np.allclose(pose1.t, pose2.t, rtol=1e-05, atol=1e-08)

        def caseWaypointTooClose(t1, j1, t2, j2):
            joint1 = t1.Joints[j1]
            joint2 = t2.Joints[j2]

            if (not isWaypoint(joint1)) or (not isWaypoint(joint2)):
                return False

            return posesAreSame(joint1.DistalDubinsFrame(), joint2.ProximalDubinsFrame()) or posesAreSame(joint1.ProximalDubinsFrame(), joint2.DistalDubinsFrame())# or (t1.Parents[j1] != None and posesAreSame(joint1.ProximalDubinsFrame(), t1.Joints[t1.Parents[j1]].DistalDubinsFrame()) and caseWaypointTooClose(t1, t1.Parents[j1], t2, j2)) or (t2.Parents[j2] != None and posesAreSame(joint2.ProximalDubinsFrame(), t2.Joints[t2.Parents[j2]].DistalDubinsFrame()) and caseWaypointTooClose(t1, j1, t2, t2.Parents[j2]))
        
        def caseWaypointOnTopOfBranch(t1, l1, t2, l2):
            while t1.Parents[t1.Parents[l1]] != None and t1.Parents[t1.Parents[l1]] != -1 and isWaypoint(t1.Joints[t1.Parents[l1]]) and posesAreSame(t1.Joints[t1.Parents[l1]].ProximalDubinsFrame(), t1.Joints[t1.Parents[t1.Parents[l1]]].DistalDubinsFrame()):
                l1 = t1.Parents[l1]
            while t2.Parents[t2.Parents[l2]] != None and t2.Parents[t2.Parents[l2]] != -1 and isWaypoint(t2.Joints[t2.Parents[l2]]) and posesAreSame(t2.Joints[t2.Parents[l2]].ProximalDubinsFrame(), t2.Joints[t2.Parents[t2.Parents[l2]]].DistalDubinsFrame()):
                l2 = t2.Parents[l2]

            if isWaypoint(t1.Joints[l1]):
                if linksInSameBranch(t1, l1, t2, l2):
                    return True
            if isWaypoint(t2.Joints[l2]):
                if linksInSameBranch(t1, l1, t2, l2):
                    return True
            
            return False

        def waypointOnTopOfJoint(t1, j1, t2, j2):
            joint1 = t1.Joints[j1]
            joint2 = t2.Joints[j2]
            
            return t1.Parents[j1] == j2 or t2.Parents[j2] == j1 or posesAreSame(t1.Links[j1].StartDubinsPose, joint2.DistalDubinsFrame()) or posesAreSame(joint1.DistalDubinsFrame(), t2.Links[j2].StartDubinsPose)# or (isWaypoint(t1.Parents[j1]) and waypointOnTopOfJoint(t1, t1.Parents[j1], t2, j2)) or (isWaypoint(t2.Parents[j2]) and waypointOnTopOfJoint(t1, j1, t2, t2.Parents[j2]))
        
        def linksInSameBranch(t1, l1, t2, l2):
            while t1.Parents[t1.Parents[l1]] != None and t1.Parents[t1.Parents[l1]] != -1 and isWaypoint(t1.Joints[t1.Parents[l1]]) and posesAreSame(t1.Joints[t1.Parents[l1]].ProximalDubinsFrame(), t1.Joints[t1.Parents[t1.Parents[l1]]].DistalDubinsFrame()):
                l1 = t1.Parents[l1]
            while t2.Parents[t2.Parents[l2]] != None and t2.Parents[t2.Parents[l2]] != -1 and isWaypoint(t2.Joints[t2.Parents[l2]]) and posesAreSame(t2.Joints[t2.Parents[l2]].ProximalDubinsFrame(), t2.Joints[t2.Parents[t2.Parents[l2]]].DistalDubinsFrame()):
                l2 = t2.Parents[l2]

            return posesAreSame(t1.Links[l1].StartDubinsPose, t2.Links[l2].StartDubinsPose)

        for i in others:
            #NEED TO SET i=i IN LAMBDA BECAUSE CAPTURE BY REFERENCE
            joint = self.Joints[i]

            allCapsules[0].append((i, joint.collisionCapsules, lambda tree, jointIdx: True))
            
            #link-joint

            allCapsules[1].append((i, joint.collisionCapsules, lambda tree, linkIdx, i=i: not (tree.Parents[linkIdx] == i or caseWaypointTooClose(tree, tree.Parents[linkIdx], self, i) or waypointOnTopOfJoint(self, i, tree, linkIdx) or linksInSameBranch(self, i, tree, linkIdx) or caseWaypointOnTopOfBranch(self, i, tree, linkIdx))))

            link = self.Links[i]
            
            #joint-link
            allCapsules[0].append((i, link.collisionCapsules, lambda tree, jointIdx, i=i: not (self.Parents[i] == jointIdx or caseWaypointTooClose(self, self.Parents[i], tree, jointIdx) or waypointOnTopOfJoint(self, i, tree, jointIdx))))
            
            #link-link

            #base of branch
            allCapsules[1].append((i, link.collisionCapsules, lambda tree, linkIdx, i=i: not (self.Parents[i] == tree.Parents[linkIdx] or caseWaypointTooClose(self, i, tree, linkIdx) or waypointOnTopOfJoint(self, i, tree, linkIdx) or linksInSameBranch(self, i, tree, linkIdx) or caseWaypointOnTopOfBranch(self, i, tree, linkIdx))))

            #end cap (taken care of by joints)
            # allCapsules[1].append((i, [link.collisionCapsules[-1]], lambda tree, linkIdx, i=i: not caseWaypointTooClose(self, i, tree, tree.Parents[linkIdx]) and not waypointOnTopOfJoint(self, i, tree, linkIdx)))

        return allCapsules

    def detectCollisionsWithCapsules(self, indices, capsulesToCheck, show=False,debug=False):
        # vectorized doesn't improve speed
        # start = time.time()
        # def process_capsules(capsule_source, check_list):
        #     all_capsules = np.concatenate([getattr(self, capsule_source)[idx].collisionCapsules for idx in indices])
        #     all_check_capsules = np.concatenate([capsuleList for capsuleList, _ in check_list])
        #     all_funcs = np.concatenate([np.concatenate([[func(self, idx)] * len(getattr(self, capsule_source)[idx].collisionCapsules) for idx in indices]).repeat(len(capsuleList)) for capsuleList, func in check_list])
            
        #     # matrix of all possible capsule pairs
        #     capsule_pairs = np.array(np.meshgrid(all_capsules, all_check_capsules)).T.reshape(-1, 2)
            
        #     # collision check
        #     collisions = np.array([all_funcs[i] and capsule_pairs[i][1].collidesWith(capsule_pairs[i][0])[0] 
        #                         for i in range(0,len(capsule_pairs))])

        #     return np.sum(collisions)

        # joint_collisions = process_capsules('Joints', capsulesToCheck[0])
        # link_collisions = process_capsules('Links', capsulesToCheck[1])

        # numCollisions = joint_collisions + link_collisions + self.detectCollisions(specificJointIndices = indices)
        # #print(time.time() - start)
        # return numCollisions

        numCollisions = 0
        angles = list(range(0,90,15))
        angles1 = angles * len(angles)
        angles2 = np.concatenate([[a] * len(angles) for a in angles])

        for idx in indices:
            capsules = self.Joints[idx].collisionCapsules
            for (i, capsuleList, func) in capsulesToCheck[0]:
                if func(self, idx):
                    for capsule2 in capsuleList:
                        for capsule1 in capsules:
                            if separatingAxisTheorem(capsule1.box, capsule2.box):
                                #if np.all(vectorized_box_collision(np.array([capsule1.box.rotate(angle) for angle in angles1]), np.array([capsule2.box.rotate(angle) for angle in angles2]))):
                                #if capsule_box_collision(capsule1.start, capsule1.end, capsule1.radius, capsule2.box.points) and capsule_box_collision(capsule2.start, capsule2.end, capsule2.radius, capsule1.box.points):
                                didCollide1, pt = capsule2.collidesWith(capsule1)
                                
                                if didCollide1:
                                    #if capsule1.collidesWith(capsule2)[0]:
                                    numCollisions += 1
                                    if debug:
                                        print("joint", idx, i)
                                    if show:
                                        self.show(addCapsules=[capsule1, capsule2], plotPoint=pt)

            
            capsules = self.Links[idx].collisionCapsules
            for (i, capsuleList, func) in capsulesToCheck[1]:
                if func(self, idx):
                    for capsule2 in capsuleList:
                        for capsule1 in capsules:
                            if separatingAxisTheorem(capsule1.box, capsule2.box):
                                #if np.all(vectorized_box_collision(np.array([capsule1.box.rotate(angle) for angle in angles1]), np.array([capsule2.box.rotate(angle) for angle in angles2]))):
                                #if capsule_box_collision(capsule1.start, capsule1.end, capsule1.radius, capsule2.box.points) and capsule_box_collision(capsule2.start, capsule2.end, capsule2.radius, capsule1.box.points):
                                didCollide1, pt = capsule2.collidesWith(capsule1)

                                if didCollide1:
                                    #if capsule1.collidesWith(capsule2)[0]:
                                    numCollisions += 1
                                    if debug:
                                        print("link", idx, i)
                                    if show:
                                        self.show(addCapsules=[capsule1, capsule2], plotPoint=pt)

        #collision between indices
        #numCollisions += self.detectCollisions(specificJointIndices = indices, debug=True, plot=True)

        #print(time.time() - start)
        return numCollisions

    def getCollisionError(self, indices, capsulesToCheck):
        totalError = 0

        for idx in indices:
            capsules = self.Joints[idx].collisionCapsules
            for (i, capsuleList, func) in capsulesToCheck[0]:
                if func(self, idx):
                    for capsule2 in capsuleList:
                        for capsule1 in capsules:
                            if separatingAxisTheorem(capsule1.box, capsule2.box):
                                totalError += capsule2.collisionErrorWith(capsule1)
                                
            
            capsules = self.Links[idx].collisionCapsules
            for (i, capsuleList, func) in capsulesToCheck[1]:
                if func(self, idx):
                    for capsule2 in capsuleList:
                        for capsule1 in capsules:
                            if separatingAxisTheorem(capsule1.box, capsule2.box):
                                totalError += capsule2.collisionErrorWith(capsule1)

        return totalError

    def detectCollisionsWithPairs(self, indices, collisionPairDict, show=True, debug=True,
                                  coarseDistanceThreshold=0.5, fineDistanceThreshold=0.001):
        numCollisions = 0
        coarseDensity = 1 / coarseDistanceThreshold
        fineDensity = 1 / fineDistanceThreshold
        
        for idx in indices:
            pairs = collisionPairDict.get(idx, [])
            
            for (obj1, obj2) in pairs:
                idx1, type1 = obj1
                idx2, type2 = obj2

                tube1 = self.Joints[idx1] if type1 == 'joint' else self.Links[idx1]
                tube2 = self.Joints[idx2] if type2 == 'joint' else self.Links[idx2]
                
                # Filter out if either tube is empty (no length)
                epsilon = 1e-2 * fineDistanceThreshold
                if tube1.length() > epsilon and tube2.length() > epsilon:
                    # Filter based on distance between bounding boxes
                    min1, max1 = tube1.boundingBox()
                    min2, max2 = tube2.boundingBox()
                    if sdf_aabb(min1, max1, min2, max2) < 0: # bounding boxes overlap
                        # Check SDF values at coarse density to see if close enough to consider collision
                        if np.min(tube2.sdf(tube1.interpolate(density = coarseDensity))) < coarseDistanceThreshold or \
                        np.min(tube1.sdf(tube2.interpolate(density = coarseDensity))) < coarseDistanceThreshold:
                            # Finer check
                            points1 = tube1.interpolate(density = fineDensity)
                            dists1to2 = tube2.sdf(points1)
                            minIdx1to2 = np.argmin(dists1to2)
                            minDist1to2 = dists1to2[minIdx1to2]
                            points2 = tube2.interpolate(density = fineDensity)
                            dists2to1 = tube1.sdf(points2)
                            minIdx2to1 = np.argmin(dists2to1)
                            minDist2to1 = dists2to1[minIdx2to1]
                            if minDist1to2 < tube1.r and minDist2to1 < tube2.r:
                                # Filter out collisions that are solely in the hemispherical end caps
                                # Check if both closest points are at endpoints
                                isEndpoint1 = (minIdx1to2 == 0 or minIdx1to2 == len(points1) - 1)
                                isEndpoint2 = (minIdx2to1 == 0 or minIdx2to1 == len(points2) - 1)
                                
                                isTrueCollision = True
                                if isEndpoint1 and isEndpoint2:
                                    # Both closest points are endpoints - need to check disc intersection
                                    # Get the circles at the relevant endpoints
                                    circle1 = tube1.startCircle(forward=False) if (minIdx1to2 == 0) else tube1.endCircle(forward=True)
                                    circle2 = tube2.startCircle(forward=False) if (minIdx2to1 == 0) else tube2.endCircle(forward=True)
                                    
                                    # It's definitely a true collision if it's on the inside side of either disc's plane
                                    # Otherwise, check if the discs cross
                                    isTrueCollision = Plane(circle1.c, circle1.n).signedDistanceToPoint(points2[minIdx2to1]) < 0\
                                                or Plane(circle2.c, circle2.n).signedDistanceToPoint(points1[minIdx1to2]) < 0\
                                                or discs_cross(circle1, circle2)
                                
                                if isTrueCollision:
                                    numCollisions += 1
                                    if show:
                                        if minDist1to2 < minDist2to1:
                                            pt = points1[minIdx1to2]
                                        else:
                                            pt = points2[minIdx2to1]
                                        self.show(plotPoint=pt, block=False)
                                        placeholder = 42
                                    if debug:
                                        print(f"{type1} {idx1} vs {type2} {idx2}")



                


                """
                # Get capsules for each object
                capsules1 = self.Joints[idx1].collisionCapsules if type1 == 'joint' else self.Links[idx1].collisionCapsules
                capsules2 = self.Joints[idx2].collisionCapsules if type2 == 'joint' else self.Links[idx2].collisionCapsules
                
                # Check all capsule pairs
                for capsule1 in capsules1:
                    for capsule2 in capsules2:
                        if separatingAxisTheorem(capsule1.box, capsule2.box):
                            didCollide, pt = capsule1.collidesWith(capsule2)
                            
                            if didCollide:
                                numCollisions += 1
                                if show:
                                    self.show(addCapsules=[capsule1, capsule2], plotPoint=pt, block=False)
                                    placeholder = 42
                                if debug:
                                    print(f"{type1} {idx1} vs {type2} {idx2}")
                """
        
        return numCollisions

    def getCollisionErrorFromDict(self, indices, collisionPairDict):
        totalError = 0
        
        for idx in indices:
            pairs = collisionPairDict.get(idx, [])
            
            for (obj1, obj2) in pairs:
                idx1, type1 = obj1
                idx2, type2 = obj2
                
                # Get capsules for each object
                capsules1 = self.Joints[idx1].collisionCapsules if type1 == 'joint' else self.Links[idx1].collisionCapsules
                capsules2 = self.Joints[idx2].collisionCapsules if type2 == 'joint' else self.Links[idx2].collisionCapsules
                
                # Sum collision errors for all capsule pairs
                for capsule1 in capsules1:
                    for capsule2 in capsules2:
                        if separatingAxisTheorem(capsule1.box, capsule2.box):
                            totalError += capsule1.collisionErrorWith(capsule2)
        
        return totalError

    def findWaypointSets(self):
        # Find clusters of waypoints and their incoming + outgoing links, that are connected without real joints between them

        visited = set()
        waypoint_sets = []
        link_sets = []
        
        def buildWaypointSet(joint_idx, current_wp_set, current_link_set):
            if joint_idx in visited or joint_idx is None or joint_idx == -1:
                return
            
            # If not a waypoint, include only the link and return
            if not isWaypoint(self.Joints[joint_idx]):
                current_link_set.append(joint_idx)
                return
            
            # Add waypoint to current set and mark as visited
            current_wp_set.append(joint_idx)
            visited.add(joint_idx)

            # Add the incoming link to the link set
            incoming_link_idx = joint_idx
            if incoming_link_idx is not None and incoming_link_idx != -1:
                current_link_set.append(incoming_link_idx)
            
            # Check parent
            parent_idx = self.Parents[joint_idx]
            if parent_idx is not None and parent_idx != -1 and parent_idx not in visited:
                if isWaypoint(self.Joints[parent_idx]):
                    buildWaypointSet(parent_idx, current_wp_set, current_link_set)
            
            # Check children
            children_indices = self.Children[joint_idx]
            for child_idx in children_indices:
                if child_idx not in visited and isWaypoint(self.Joints[child_idx]):
                    buildWaypointSet(child_idx, current_wp_set, current_link_set)
        
        # Go through all joints and find waypoint sets
        for joint_idx in range(len(self.Joints)):
            if joint_idx not in visited and isWaypoint(self.Joints[joint_idx]):
                current_wp_set = []
                current_link_set = []
                buildWaypointSet(joint_idx, current_wp_set, current_link_set)
                if current_wp_set:
                    waypoint_sets.append(current_wp_set)
                    link_sets.append(current_link_set)
        
        return link_sets
    
    def buildCollisionPairDictionary(self):
        """
        Build a dictionary mapping each joint index to collision pairs that need checking.
        
        Rules:
        - Waypoints don't collision check directly
        - Link/Link: check iff separated by real joint (different waypoint cluster)
        - Real Joint/Real Joint: always check
        - Real Joint/Link: always check
        
        Returns:
            dict: Keys are joint indices, values are lists of tuples:
                  [((idx1, type1), (idx2, type2)), ...]
                  where type is 'joint' or 'link'
        """
        link_sets = self.findWaypointSets()
        
        # Pre-classify all joints
        real_joints = [i for i in range(len(self.Joints)) if not isWaypoint(self.Joints[i])]
        
        # Create link-to-set mapping 
        whichLinkSet = {}
        for set_idx, link_set in enumerate(link_sets):
            for link_idx in link_set:
                whichLinkSet[link_idx] = set_idx
        
        collision_pairs = {}
        
        for node_idx in range(len(self.Joints)):
            # Waypoint 
            if isWaypoint(self.Joints[node_idx]):

                # Don't check for collisions
                collision_pairs[node_idx] = []
                continue
            
            pairs = []
            
            # Real joint
            for other_idx in real_joints:

                # Against other real joints
                if other_idx != node_idx:
                    pairs.append(((node_idx, 'joint'), (other_idx, 'joint')))

                # Against links
                pairs.append(((node_idx, 'joint'), (other_idx, 'link')))
            
            # Link 
            set_idx = whichLinkSet.get(node_idx)
            for other_idx in range(len(self.Links)):

                # Skips itself
                if other_idx == node_idx:
                    continue

                # If comparing with a different link, check which set that one belongs to
                other_set = whichLinkSet.get(other_idx)

                # Check if both links are in different sets or the one in question is not in a set at all
                if set_idx != other_set or set_idx is None:
                    pairs.append(((node_idx, 'link'), (other_idx, 'link')))
            
            collision_pairs[node_idx] = pairs
        
        return collision_pairs

    
    def branchingParametersFrom(self, parentIndex : int):
        linksToChildren = [self.Links[childIndex] for childIndex in self.Children[parentIndex]]
        return [link.branchingParameters() for link in linksToChildren]
    
    def exportLink3DFile(self, parentIndex : int, folder : str, fileFormat = "stl", pose=False, manifold=False):

        name = f"linkfrom_{parentIndex}_to_"
        for endpointIndex in self.Children[parentIndex]:
            name += f"{endpointIndex}_"
        name += "." + fileFormat
        
        source = self.Joints[parentIndex]

        children = self.Children[parentIndex]
        params = np.round(self.branchingParametersFrom(parentIndex), 4)

        #check link min length        
        for i in self.Children[parentIndex]:
            if self.Links[i].path.length < (self.Joints[i].screwRadius*2 + self.Joints[i].printParameters.holeMargin*2 + self.Joints[i].printParameters.gridHoleRadius):
                if (np.abs(self.Links[i].path.theta1) > 0.0001 or np.abs(self.Links[i].path.theta2) > 0.0001): #replace with some epsilon
                    raise InvalidDimensionException(f"The link between joint {parentIndex} and joint {i} are too close together for a link to be created between them.")
                else:
                    print(f"Skipping link between joint {parentIndex} and joint {i}. Extending joint {parentIndex} by {self.Links[i].path.tMag} instead")
                    source.extendSegment(self.Links[i].path.tMag)
                    #return None

        sourceParameters = source.printParameters

        defs = [f"tolerance={sourceParameters.tolerance};\n",f"hole_radius={source.screwRadius};\n",
                f"grid_hole_radius={sourceParameters.gridHoleRadius};\n",f"outer_radius={source.r};\n",
                f"thickness={sourceParameters.thickness};\n",f"hole_attach_height={sourceParameters.holeMargin};\n",
                f"attach_thickness={sourceParameters.attachThickness};\n"]
        
        if pose:
            name = "link_"
            for parameter in defs:
                name += parameter[parameter.index("=") + 1:-2] + "_"

        branch_scad = "scad/branch.scad"
        if pose:
            branch_scad = "scad/poses/branch_pose.scad"
        with open(branch_scad, "r") as file:
            lines = file.readlines()
        truncated = lines[7:239] #first 7 are parameter definitions, first 239 lines are function definitions

        linkEndpoints = [self.Joints[childIndex] for childIndex in self.Children[parentIndex]]
        new_lines = ["branch([\n"]
        for i in range(0,len(params)):
            path = params[i]
            nextInnerRadius = linkEndpoints[i].r - linkEndpoints[i].printParameters.thickness
            nextHoleMargin = linkEndpoints[i].printParameters.holeMargin
            nextScrewRadius = linkEndpoints[i].screwRadius

            if pose:
                name += f"{nextInnerRadius}_{nextHoleMargin}_{nextScrewRadius}_"
                for p in path:
                    name += f"{p}_"
            
            new_lines.append(f"[ {path[0]}, {path[1]}, {path[2]}*outer_radius, {path[3]}*outer_radius, {path[4]}, {path[5]}, {path[6]}*outer_radius, {nextScrewRadius}, {nextHoleMargin}, {nextInnerRadius}],\n")
        new_lines.append("],outer_radius,inner_radius);")

        if pose:
            if os.path.isfile(f"3d_output/{folder}/{name}.stl"):
                return f"3d_output/{folder}/{name}.stl"
            else:
                name += ".stl"
        
        try:
            with open(f"scad_output/{folder}/{name}.scad", "w+") as file:
                truncated.extend(new_lines)
                defs.extend(truncated)
                file.writelines(defs)
        except:
            #file name is too long
            name = f"linkfrom_{parentIndex}_to_"
            for endpointIndex in self.Children[parentIndex]:
                name += f"{endpointIndex}_"
            name += ".stl"

            with open(f"scad_output/{folder}/{name}.scad", "w+") as file:
                truncated.extend(new_lines)
                defs.extend(truncated)
                file.writelines(defs)

        if manifold:
            os.system(f"openscad --backend Manifold -q -o 3d_output/{folder}/{name} -m scad_output/{folder}/{name}.scad")
        else:
            os.system(f"openscad -q -o 3d_output/{folder}/{name} scad_output/{folder}/{name}.scad")

        return f"3d_output/{folder}/{name}"

    def export3DKinematicTree(self, folder = "", fileFormat = "stl", manifold=False):
        if (folder != ""):
            os.makedirs(f"scad_output/{folder}", exist_ok=True)
            os.makedirs(f"3d_output/{folder}", exist_ok=True)

        print(f"Printing modules for {folder[:-1]}...")

        tree = self.copyAbbreviatedSelf()
        #TODO: fix issue with overlapping waypoints
        # for i in range(0, len(self.Joints)):
        #     currentFrame = tree.Joints[i].DistalDubinsFrame().t

        #     while True:
        #         children = tree.Children[i].copy()
        #         for j in range(0, len(children)):
        #             if isWaypoint(tree.Joints[children[j]]) and np.allclose(tree.Joints[children[j]].ProximalDubinsFrame().t, currentFrame, rtol=1e-05, atol=1e-08):
        #                 tree.Children[i].remove(children[j])
        #                 for child in tree.Children[children[j]]:
        #                     tree.Parents[child] = i
        #                 tree.Children[i] += tree.Children[children[j]]
        #                 tree.Children[children[j]] = []
                
        #         if collections.Counter(children) == collections.Counter(tree.Children[i]):
        #             break

        # for i in range(0, len(tree.Joints)):
        #     tree.transformJoint(i, SE3(), recomputeLinkPath=True, safe=False)

        # print("Done shrinking tree.")

        #export all the links
        for i in range(0,len(tree.Children)):
            start = time.time()
            if len(tree.Children[i]) > 0:
                tree.exportLink3DFile(i,folder,fileFormat,manifold=manifold)
            print(f"Finished link {i}/{len(tree.Children) - 1}, Time: {time.time() - start} \r")
        
        #export all the joints
        for i in range(0,len(tree.Joints)):
            start = time.time()
            if not isinstance(tree.Joints[i],PrintedWaypoint):
                tree.Joints[i].export3DFile(i,folder,fileFormat,manifold=manifold)
            print(f"Finished joint {i}/{len(tree.Joints) - 1}, Time: {time.time() - start} \r")
        

    def show(self, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             showJointSurface=True, jointColor=jointColorDefault, jointEdgeColor = jointEdgeColorDefault, 
             jointAxisScale=jointAxisScaleDefault, showJointPoses=True,
             linkColor=linkColorDefault, surfaceOpacity=surfaceOpacityDefault, showLinkSurface=True, 
             showLinkPoses=False, showLinkPath=True, pathColor=pathColorDefault,
             showPathCircles=False, sphereColor=sphereColorDefault,
             showSpheres=False, block=blockDefault, showAxisGrids=False, 
             showGlobalFrame=False, globalAxisScale=globalAxisScaleDefault,
             showGroundPlane=False, groundPlaneScale=groundPlaneScaleDefault,
             groundPlaneColor=groundPlaneColorDefault, showCollisionBoxes=([],[]), 
             showSpecificCapsules=([],[]), plotPoint = None, addCapsules=[],
             showScaleBar=True):
        ax = plt.figure().add_subplot(projection='3d')
        if showGroundPlane: #https://stackoverflow.com/questions/36060933/plot-a-plane-and-points-in-3d-simultaneously
            xx, yy = np.meshgrid(range(groundPlaneScale), range(groundPlaneScale))
            xx = xx - groundPlaneScale/2
            yy = yy - groundPlaneScale/2
            z = 0*xx #(9 - xx - yy) / 2 

            # plot the plane
            ax.plot_surface(xx, yy, z, alpha=surfaceOpacity/4, color=groundPlaneColor)
                
        xyzHandles, abcHandles = self.addToPlot(ax, xColor=xColor, yColor=yColor, zColor=zColor,
                                                proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor,
                                                showJointSurface=showJointSurface, jointColor=jointColor,
                                                jointAxisScale=jointAxisScale, showJointPoses=showJointPoses,
                                                linkColor=linkColor, surfaceOpacity=surfaceOpacity, showLinkSurface=showLinkSurface,
                                                showLinkPoses=showLinkPoses, showLinkPath=showLinkPath, pathColor=pathColor,
                                                showPathCircles=showPathCircles, sphereColor=sphereColor,
                                                showSpheres=showSpheres, showGlobalFrame=showGlobalFrame, globalAxisScale=globalAxisScale,
                                                showCollisionBoxes=showCollisionBoxes, showSpecificCapsules=showSpecificCapsules, 
                                                plotPoint=plotPoint, addCapsules=addCapsules)


        # Get the current limits of the axes
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        zlim = ax.get_zlim()

        # Define the scale bar length
        scale_bar_length = 10

        scale_bar_x = [xlim[1] - scale_bar_length - 0.05 * (xlim[1] - xlim[0]), xlim[1] - 0.05 * (xlim[1] - xlim[0])]
        scale_bar_y = [ylim[0] + 0.05 * (ylim[1] - ylim[0]), ylim[0] + 0.05 * (ylim[1] - ylim[0])]
        scale_bar_z = [zlim[0] + 0.05 * (zlim[1] - zlim[0]), zlim[0] + 0.05 * (zlim[1] - zlim[0])]

        if showScaleBar:
            # Plot the scale bar
            ax.plot(scale_bar_x, scale_bar_y, scale_bar_z, color='black', linewidth=2)

            # Label the scale bar
            ax.text(scale_bar_x[0] + scale_bar_length / 2, scale_bar_y[0], scale_bar_z[0], 
                    f'{scale_bar_length} units', color='black', fontsize=10)

        handleGroups = []
        labels = []
        if showJointPoses or showGlobalFrame:
            xHats = xyzHandles[:,0]
            yHats = xyzHandles[:,1]
            zHats = xyzHandles[:,2]
            origins = xyzHandles[:,3]
            handleGroups += [tuple(xHats), tuple(yHats), tuple(zHats)]
            labels += [r'$\^x$', r'$\^y$', r'$\^z$']
        if showLinkPoses:
            aHats = abcHandles[:,0]
            bHats = abcHandles[:,1]
            cHats = abcHandles[:,2]
            handleGroups += [tuple(aHats), tuple(bHats), tuple(cHats)]
            labels += [r'$\^a$', r'$\^b$', r'$\^c$']
        if not handleGroups==[]:
            ax.legend(handleGroups, labels)
        
        ax.set_aspect('equal')
        if not showAxisGrids:
            plt.axis('off')

        plt.savefig('img_output.png', dpi=800)
        plt.show(block=block)


    def transformAll(self, Transformation : SE3, recomputeBoundingBall : bool = True):
        self.transformJoint(0, Transformation, safe=False, recomputeBoundingBall=recomputeBoundingBall)

    """ 
    Apply given transformation (SE3() object) to given joint (index), 
    and to its descendants if propogate (defaults to True).

    Returns True if it succeeds (the transformation gives valid links).
    In safe=True mode (default), if it fails it will leave the chain unchanged,
    print a warning, and return False rather than throwing an error.   
    """
    def transformJoint(self, jointIndex : int, Transformation : SE3, 
                       propogate : bool = True, recomputeBoundingBall : bool = True,
                       recomputeLinkPath : bool = True, 
                       safe : bool = True, relative : bool = False, printErrors=False) -> bool:
        if relative:
            Transformation = self.Joints[jointIndex].Pose @ Transformation @ self.Joints[jointIndex].Pose.inv()
        
        if safe:
            backup = self.dataDeepCopy()
            try:
                self.transformJoint(jointIndex, Transformation, 
                       propogate, recomputeBoundingBall,
                       recomputeLinkPath, safe=False, relative=False)
            except ValueError as err:
                if printErrors:
                    print("WARNING: Value Error in transformJoint:")
                    print(err)
                    print("Reverting chain to before outer call.")
                self.setTo(backup)
                return False
            except AssertionError as err:
                if printErrors:
                    print("WARNING: Assertion Error in transformJoint:")
                    print(err)
                    print("Reverting chain to before outer call.")
                self.setTo(backup)
                return False

        else:
            self.Joints[jointIndex].transformPoseBy(Transformation)
            joint = self.Joints[jointIndex]
            if recomputeLinkPath and jointIndex > 0:
                parent = self.Joints[self.Parents[jointIndex]]

                self.Links[jointIndex] = LinkCSC(self.r, parent.DistalDubinsFrame(), 
                                        joint.ProximalDubinsFrame(),
                                        self.maxAnglePerElbow)
            else:
                self.Links[jointIndex] = self.Links[jointIndex].newLinkTransformedBy(Transformation)
            if propogate:
                for c in self.Children[jointIndex]:
                    self.transformJoint(c, Transformation, propogate=True, 
                                        recomputeBoundingBall=False,
                                        recomputeLinkPath=False,
                                        safe=False, relative=False)
                
                self.recursivelyRecomputeCollisionCapsules(jointIndex)
            else:
                for c in self.Children[jointIndex]:
                    child = self.Joints[c]
                    self.Links[c] = LinkCSC(self.r, joint.DistalDubinsFrame(), 
                                            child.ProximalDubinsFrame(),
                                            self.maxAnglePerElbow)
                self.Joints[jointIndex].recomputeCollisionCapsules()
            
            if recomputeBoundingBall:
                self.recomputeBoundingBall()
            elif jointIndex == 0:
                self.boundingBall = self.boundingBall.newBallTransformedBy(Transformation)

        return True
                
    def setJointState(self, jointIndex : int, newState : float) -> bool:
        joint = self.Joints[jointIndex]

        if joint is None:
            return False
            
        minState, maxState = joint.stateRange()

        # no change required, return true
        if newState == joint.state or np.abs(newState - joint.state) < 1e-8 * (maxState - minState):
            return True

        if newState < minState or newState > maxState:
            # print("WARNING: state out of range in setJointRange, "+
            #         "state unchanged.")
            return False
        
        Transformation = joint.TransformStateTo(newState)
        for c in self.Children[jointIndex]:
            # TODO: it is possible, and would be more efficient, to make this 
            # transform the existing links rather than recompute them
            self.transformJoint(c, Transformation, propogate=True, recomputeLinkPath=False,
                                recomputeBoundingBall=False, safe=False)
        self.recomputeBoundingBall()

        self.recursivelyRecomputeCollisionCapsules(jointIndex)
        return True
    
    def recursivelyRecomputeCollisionCapsules(self, index):
        self.Joints[index].recomputeCollisionCapsules()
        for child in self.Children[index]:
            self.recursivelyRecomputeCollisionCapsules(child)
    
    # Returns True if it succeeds (the transformation gives valid links).
    # In safe=True mode (default), if it fails it will leave the chain unchanged,
    # print a warning, and return False rather than throwing an error.    
    def translateJointAlongAxisOfMotion(self, jointIndex : int, 
                                         distance : float, 
                                         propogate : bool = True, 
                                         applyToPreviousWaypoint : bool = False, 
                                         safe : bool = True) -> bool:
        if safe:
            backup = self.dataDeepCopy()
            try:
                self.translateJointAlongAxisOfMotion(jointIndex, distance, 
                            propogate, applyToPreviousWaypoint, safe = False)
            except ValueError as err:
                print("WARNING: something went wrong in translateJointAlongAxisOfMotion:")
                print(err)
                print("Reverting chain to before outer call.")
                self.setTo(backup)
                return False
        else:
            Translation = SE3(distance * self.Joints[jointIndex].Pose.R[:,2])
            if applyToPreviousWaypoint and type(self.Joints[jointIndex-1])==Waypoint:
                if propogate:
                    self.transformJoint(jointIndex-1, Translation, True, safe=False)
                else:
                    self.transformJoint(jointIndex-1, Translation, False, safe=False)
                    self.transformJoint(jointIndex, Translation, False, safe=False)
            else:
                self.transformJoint(jointIndex, Translation, propogate, safe=False)
        return True
    
    # Returns True if it succeeds (the transformation gives valid links).
    # In safe=True mode (default), if it fails it will leave the chain unchanged,
    # print a warning, and return False rather than throwing an error.   
    def rotateJointAboutAxisOfMotion(self, jointIndex : int, angle : float,
                             propogate : bool = True, 
                             applyToPreviousWaypoint : bool = False,
                             safe : bool = True) -> bool:
        if safe:
            backup = self.dataDeepCopy()
            try:
                self.rotateJointAboutAxisOfMotion(jointIndex, angle, propogate, 
                             applyToPreviousWaypoint, safe = False)
            except ValueError as err:
                print("WARNING: something went wrong in rotateJointAboutAxisOfMotion:")
                print(err)
                print("Reverting chain to before outer call.")
                self.setTo(backup)
                return False
        else:
            Pose = self.Joints[jointIndex].Pose
            Rotation = RotationAboutLine(Pose.R[:,2], Pose.t, angle)
            if applyToPreviousWaypoint and type(self.Joints[jointIndex-1])==Waypoint:
                if propogate:
                    self.transformJoint(jointIndex-1, Rotation, True, safe=False)
                else:
                    self.transformJoint(jointIndex-1, Rotation, False, safe=False)
                    self.transformJoint(jointIndex, Rotation, False, safe=False)
            else:
                self.transformJoint(jointIndex, Rotation, propogate, safe=False)
        return True
    
    def save(self, filename: str, saveDir = True):
        #TODO: ADD EXTENSIONS FOR PRINTED JOINTS
        name = filename + ".tree"
        if saveDir:
            name = "save/" + filename
        
        with open(name, "w") as f:
            save = str(self.maxAnglePerElbow) + "\n"
            for i in range(0, len(self.Joints)):
                joint = self.Joints[i]
                
                save += str(self.Parents[i]) + " "
                if isinstance(joint, Waypoint):
                    save += "Waypoint " + str(joint.numSides) + " " + str(joint.r) + " " + str(joint.pidx) + " "
                elif isinstance(joint, RevoluteJoint):
                    save += "RevoluteJoint " + str(joint.numSides) + " " + str(joint.r) + " " + str(joint.totalBendingAngle) + " " + str(joint.numSinkLayers) + " " + str(joint.initialState) + " "
                elif isinstance(joint, ExtendedRevoluteJoint):
                    save += "ExtendedRevoluteJoint " + str(joint.numSides) + " " + str(joint.r) + " " + str(joint.totalBendingAngle) + " " + str(joint.tubeLength) + " " + str(joint.numSinkLayers) + " " + str(joint.initialState) + " "
                elif isinstance(joint, PrismaticJoint):
                    save += "PrismaticJoint " + str(joint.numSides) + " " + str(joint.r) + " " + str(joint.neutralLength) + " " + str(joint.numLayers) + " " + str(joint.coneAngle) + " " + str(joint.initialState) + " "
                elif isinstance(joint, Tip):
                    save += "Tip " + str(joint.numSides) + " " + str(joint.r) + " " + str(joint.neutralLength) + " " + str(joint.forward) + " " + str(joint.pidx) + " "
                elif isinstance(joint, PrintedWaypoint):
                    save += "PrintedWaypoint " + str(joint.r) + " " + str(joint.screwRadius) + " " + str(joint.pidx) + " " + joint.printParameters.toString() + " "
                elif isinstance(joint, PrintedPrismaticJoint):
                    save += "PrintedPrismaticJoint " + str(joint.r) + " " + str(joint.extensionLength) + " " + str(joint.screwRadius) + " " + str(joint.initialState) + " " + str(joint.minLength) + " " + joint.printParameters.toString() + " "
                elif isinstance(joint, PrintedOrthogonalRevoluteJoint):
                    save += "PrintedOrthogonalRevoluteJoint " + str(joint.r) + " " + str(joint.startBendingAngle) + " " + str(joint.endBendingAngle) + " " + str(joint.screwRadius) + " " + str(joint.initialState) + " " + str(joint.bottomLength) + " " + str(joint.topLength) + " " + joint.printParameters.toString() + " "
                elif isinstance(joint, PrintedInAxisRevoluteJoint):
                    save += "PrintedInAxisRevoluteJoint " + str(joint.r) + " " + str(joint.neutralLength) + " " + str(joint.screwRadius) + " " + str(joint.initialState) + " " + joint.printParameters.toString() + " "
                elif isinstance(joint, PrintedTip):
                    save += "PrintedTip " + str(joint.r) + " " + str(joint.screwRadius) + " " + str(joint.pidx) + " " + joint.printParameters.toString() + " "
                else:
                    raise Exception("Not Implemented")
                save += "[" + ''.join([str(x) + "," for x in joint.Pose.A.reshape((16,)).tolist()])
                save += "\n"
            
            f.write(save)
            f.close()

    def totalLengthLowerBound(self):
        sum = self.Joints[0].neutralLength/2
        for i in range(1, len(self.Joints)):
            # if this isn't a waypoint
            if not isinstance(self.Joints[i], Waypoint):
                # backtrack until finding the parent that's not a waypoint, or reach the root
                p = self.Parents[i]
                while isinstance(self.Joints[p], Waypoint) and p > 0:
                    p = self.Parents[p]
                
                parent = self.Joints[p]
                child = self.Joints[i]
                
                jointLength = parent.neutralLength/2 + child.neutralLength/2
                distanceBetweenZaxes = shortestDistanceBetweenLines(parent.Pose.t,
                                                                parent.Pose.R[:,2],
                                                                child.Pose.t,
                                                                child.Pose.R[:,2])
                sum += max(jointLength, distanceBetweenZaxes)
        return sum

    def totalLength(self):
        sum = 0
        for i in range(len(self.Joints)):
            sum += self.Joints[i].neutralLength
        for i in range(1, len(self.Links)):
            sum += self.Links[i].path.length
        return sum

    def leaves(self):
        return [i for i in range(len(self.Joints)) if len(self.Children[i]) == 0]
    
    def nonLeaves(self):
        return [i for i in range(len(self.Joints)) if len(self.Children[i]) > 0]

    def realJointIndices(self) -> np.ndarray:
        """
        Returns a numpy array of the real (non-waypoint) joints in the chain.
        """
        return np.array([i for i in range(len(self.Joints)) if not isinstance(self.Joints[i], Waypoint) and not isinstance(self.Joints[i], PrintedWaypoint) and not isinstance(self.Joints[i], EndTip)])

    # Returns the current configuration of the chain as a numpy array of joint states.
    def configuration(self, realJointsOnly : bool = False) -> np.ndarray:
        if realJointsOnly:
            return np.array([self.Joints[i].state for i in self.realJointIndices()])
        else:
            return np.array([joint.state for joint in self.Joints])
    
    def setConfiguration(self, newConfig, realJointsOnly : bool = False):
        if realJointsOnly:
            if len(newConfig) != len(self.realJointIndices()):
                print(len(newConfig), len(self.realJointIndices()))
                raise ValueError("Length mismatch between new configuration and real joints")
            for i, jointIndex in enumerate(self.realJointIndices()):
                self.setJointState(i, newConfig[i])
        else:
            if len(newConfig) != len(self.Joints):
                print(len(newConfig), len(self.Joints))
                raise ValueError("Length mismatch between new configuration and all joints")
            for i, joint in enumerate(self.Joints):
                self.setJointState(i, newConfig[i])
    
    # Returns the state ranges of the joints in the chain as a 2D numpy array.
    # Each row represents a joint, and each column represents the min and max state of that joint.
    def stateRanges(self, realJointsOnly : bool = False) -> np.ndarray:
        if realJointsOnly:
            return np.array([self.Joints[i].stateRange() for i in self.realJointIndices()])
        else:
            return np.array([joint.stateRange() for joint in self.Joints])
    
    def randomConfiguration(self, realJointsOnly : bool = False) -> np.ndarray:
        if realJointsOnly:
            return np.array([np.random.uniform(*self.Joints[i].stateRange()) for i in self.realJointIndices()])
        else:
            return np.array([np.random.uniform(*joint.stateRange()) for joint in self.Joints])

def optimizationLoss(tree):
    loss = 0
    for i in range(1, len(tree.Joints)):
        loss += tree.Links[i].path.length
    return loss

def loadKinematicTree(filename : str):
    def getJoint(line):
        first = line.split(' ')
        pose = SE3(np.array([float(x) for x in line.split('[')[1].split(",")[:-1]]).reshape(4,4))
        match first[1]:
            case "Waypoint":
                numSides = int(first[2])
                r = float(first[3])
                pathIndex = int(first[4])
                return Waypoint(numSides, r, pose, pathIndex)
            case "RevoluteJoint":
                numSides = int(first[2])
                r = float(first[3])
                totalBendingAngle = float(first[4])
                numSinkLayers = int(first[5])
                initialState = float(first[6])
                return RevoluteJoint(numSides, r, totalBendingAngle, pose, numSinkLayers, initialState)
            case "ExtendedRevoluteJoint":
                numSides = int(first[2])
                r = float(first[3])
                totalBendingAngle = float(first[4])
                tubeLength = float(first[5])
                numSinkLayers = int(first[6])
                initialState = float(first[7])
                return ExtendedRevoluteJoint(numSides, r, totalBendingAngle, tubeLength, pose, numSinkLayers, initialState)
            case "PrismaticJoint":
                numSides = int(first[2])
                r = float(first[3])
                neutralLength = float(first[4])
                numLayers = int(first[5])
                coneAngle = float(first[6])
                initialState = float(first[7])
                return PrismaticJoint(numSides, r, neutralLength, numLayers, coneAngle, pose, initialState)
            case "Tip":
                numSides = int(first[2])
                r = float(first[3])
                length = float(first[4])
                closesForward = bool(first[5])
                pidx = int(first[6])
                return Tip(numSides, r, pose, length, closesForward, pathIndex=pidx)
            case "PrintedWaypoint":
                r = float(first[2])
                screwRadius = float(first[3])
                pathIndex = int(first[4])
                printParameters = PrintParameters.fromString(first[5])
                return PrintedWaypoint(r, pose, screwRadius, pathIndex, printParameters)
            case "PrintedPrismaticJoint":
                r = float(first[2])
                extensionLength = float(first[3])
                screwRadius = float(first[4])
                initialState = float(first[5])
                minLength = float(first[6])
                printParameters = PrintParameters.fromString(first[7])
                newJoint = PrintedPrismaticJoint(r, extensionLength, pose, screwRadius, printParameters, initialState)
                newJoint.extendSegment(minLength - newJoint.minLength)
                return newJoint
            case "PrintedOrthogonalRevoluteJoint":
                r = float(first[2])
                startAngle = float(first[3])
                endAngle = float(first[4])
                screwRadius = float(first[5])
                initialState = float(first[6])
                bottomLength = float(first[7])
                topLength = float(first[8])
                printParameters = PrintParameters.fromString(first[9])
                newJoint = PrintedOrthogonalRevoluteJoint(r, startAngle, endAngle, pose, screwRadius, printParameters, initialState)
                newJoint.extendSegment(topLength - newJoint.topLength)
                newJoint.extendBottomSegment(bottomLength - newJoint.bottomLength)
                return newJoint
            case "PrintedInAxisRevoluteJoint":
                r = float(first[2])
                neutralLength = float(first[3])
                screwRadius = float(first[4])
                initialState = float(first[5])
                printParameters = PrintParameters.fromString(first[6])
                return PrintedInAxisRevoluteJoint(r, neutralLength, pose, screwRadius, printParameters, initialState)
            case "PrintedTip":
                r = float(first[2])
                screwRadius = float(first[3])
                pathIndex = int(first[4])
                printParameters = PrintParameters.fromString(first[5])
                return PrintedTip(r, pose, screwRadius, pathIndex, printParameters)


        raise Exception(f"{first[1]} not implemented in save")
            
    try:
        with open(filename) as f:
            lines = f.readlines()
            rootJoint = getJoint(lines[1])
            if isinstance(rootJoint, OrigamiJoint):
                tree = KinematicTree[OrigamiJoint](rootJoint, float(lines[0]))
            else:
                tree = KinematicTree[PrintedJoint](rootJoint, float(lines[0]))
            for i in range(2, len(lines)):
                parent = int(lines[i].split(" ")[0])
                tree.addJoint(parent, getJoint(lines[i]), relative=False, fixedPosition=True, fixedOrientation=True, safe=False)
            
            return tree
    except Exception as e:
        print(e)
        raise Exception(f"file {filename} doesnt exist")

def isWaypoint(joint):
    if joint is None:
        return False
    return isinstance(joint, Waypoint) or isinstance(joint, PrintedWaypoint)

def curvinessOfLink(link : LinkCSC):
    return link.path.theta1 ** 1.5 * link.path.r + link.path.theta2 ** 1.5 * link.path.r

def origamiToPrinted(t : KinematicTree[OrigamiJoint], screwRadius: float):
    tree = copy.deepcopy(t)

    newTree = KinematicTree[PrintedJoint](tree.Joints[0].toPrinted(screwRadius), tree.maxAnglePerElbow)
    for i in range(1, len(tree.Joints)):
        try:
            tree.setJointState(i, tree.Joints[i].initialState)
            newJoint = tree.Joints[i].toPrinted(screwRadius)
            parent = newTree.Joints[tree.Parents[i]]
            
            dist_to_prox = tree.Joints[tree.Parents[i]].DistalDubinsFrame().inv() * tree.Joints[i].ProximalDubinsFrame()
            prox_to_pose = newJoint.ProximalDubinsFrame().inv() * newJoint.Pose
            newJoint.Pose = parent.DistalDubinsFrame() @ dist_to_prox @ prox_to_pose
            newTree.addJoint(tree.Parents[i], newJoint, relative=False, safe=False, 
            fixedPosition=True, fixedOrientation=True)
        except Exception as e:
            #EDGE CASE WHERE WAYPOINT POSE IS SAME AS PARENT, has some rounding error
            if (tree.Joints[tree.Parents[i]].Pose == tree.Joints[i].Pose):
                newTree.addJoint(tree.Parents[i], copy.deepcopy(newTree.Joints[tree.Parents[i]]), relative=False, safe = False, fixedPosition=True, fixedOrientation=True)
            else:
                print(f"Unable to convert tree to 3D print because of joint {i} (parent is joint {tree.Parents[i]}): {e}\n(Try increasing placing joints further apart)")
                return None
    return newTree

def printedToOrigami(tree : KinematicTree[PrintedJoint], numSides: int, numLayers : int = 1):
    newTree = KinematicTree[OrigamiJoint](tree.Joints[0].toOrigami(numSides, numLayers), tree.maxAnglePerElbow)
    for i in range(1, len(tree.Joints)):
        try:
            tree.setJointState(i, tree.Joints[i].initialState)
            newTree.addJoint(tree.Parents[i], tree.Joints[i].toOrigami(numSides, numLayers), relative=False, safe=False, 
            fixedPosition=True, fixedOrientation=True)
        except Exception as e:
            print(f"Unable to convert tree to origami because of joint {i} (parent is joint {tree.Parents[i]}): {e}\n(Try adjusting parameters)")
            return None
    return newTree

""" 
Places joint along its joint axis, as close as possible to the given neighbor 
joint while >= 4r from the given ball. Modifies and returns joint.
 - joint is interpreted in global coordinates.
 - neighbor is usually for a parent, but could be for a child if 
   building a chain backwards (as in the 2022 paper).
 - ball is intended to enclose other parts of the tree whose location is 
   already fixed. If building a chain backwards, for example, this would be 
   all the descendant joints.

Placing the joint 4r outside of the bounding sphere is the approach from 
Algorithm 8 from the 2022 Kinegami paper: it does NOT guarantee local 
collision-avoidance, and does not insert waypoints.
"""
def moveJointNearNeighborBut4rFromBall(jointToPlace, neighbor, ball):
    def newPosition(zChange):
        return jointToPlace.Pose.t + zChange*jointToPlace.Pose.R[:,2]

    # optimization objective
    def distanceFromNeighbor(zChange): 
        return norm(neighbor.Pose.t - newPosition(zChange))
    
    # for constraint 
    def distanceBetweenBallCenters(zChange):
        return norm(ball.c - newPosition(zChange))
        
    farEnough = NonlinearConstraint(distanceBetweenBallCenters, 
        lb= 4*jointToPlace.r + ball.r + jointToPlace.boundingRadius(), 
        ub= np.inf)
        
    result = minimize(distanceFromNeighbor, 0, constraints=(farEnough))
    zChange = result.x[0]
    jointToPlace.translateAlongZ(zChange)
    return jointToPlace

"""
Places joint along its joint axis, as close as possible to the given neighbor
(parent) while ensuring its whole bounding sphere is >= 4r from the given plane.
Modifies and returns joint.
 - joint is interpreted in global coordinates.
 - plane is intended to represent the end plane of the bounding cylinder
    in the tree construction algorithm.
"""
def moveJointNearNeighborBut4rPastPlane(jointToPlace, neighbor, plane : Plane):
    def newPosition(zChange):
        return jointToPlace.Pose.t + zChange*jointToPlace.Pose.R[:,2]

    # optimization objective
    def distanceFromNeighbor(zChange): 
        return norm(neighbor.Pose.t - newPosition(zChange))
    
    # for constraint 
    def distanceInFront(zChange):
        return plane.signedDistanceToPoint(newPosition(zChange))
        
    farEnough = NonlinearConstraint(distanceInFront, 
        lb= 4*jointToPlace.r + jointToPlace.boundingRadius(), 
        ub= np.inf)
        
    result = minimize(distanceFromNeighbor, 0, constraints=(farEnough))
    zChange = result.x[0]
    newPos = newPosition(zChange)
    jointToPlace.translateAlongZ(zChange)
    return jointToPlace



"""
Places jointToPlace along its joint axis, as close as possible to the given 
neighbor while inserting intermediate waypoints as appropriate to guarantee
no local collision avoidance, based on the parallel-plane strategy from
Conjecture 7 and Algorithm 9 of the 2022 Kinegami paper.
- jointToPlace should be given in global coordinates.
- neighbor is usually a parent, but could be for a child if 
  building a chain backwards (as in the 2022 paper).
- ball is intended to enclose other parts of the tree whose location is 
  already fixed. If building a chain backwards, for example, this would be 
  all the descendant joints.
- backwards (optional, defaults to False) specifies whether we're inserting
  joints backwards (i.e., neighbor is jointToPlace's child) or forwards 
  (neighbor is jointToPlace's parent). If backwards, we need to take the 
  neighbor path direction and jointToPlace z direction in reverse.
Returns the list of joints to insert, including any intermediate waypoints.
"""
def placeJointAndWayPoints(jointToPlace, neighbor, ball, backwards=False):
    """
    This algorithm (Algorithm 9 from the 2022 Kinegami paper) 
    avoids local self-intersection (Dubins paths with turning angles >pi) 
    using Conjecture 7 from that paper.               
    The central idea is to insert 1 or 2 intermediate waypoints 
    that route the path from parallel planes >=4r apart, then place 
    the new joint along its axis based on proximity to the last 
    waypoint.
    """    
    assert(neighbor.r == jointToPlace.r)
    r = jointToPlace.r
    toReturn = []
    """ First, construct the plane tangent to the bounding sphere
    and normal to the neighbors's path direction. """
    if backwards:
        nhat1 = -neighbor.pathDirection()
    else:
        nhat1 = neighbor.pathDirection()
    # point on bounding sphere in direction nhat1
    s1 = ball.c + ball.r * nhat1
    tangentPlane1 = Plane(s1, nhat1)
    
    """ Construct a waypoint where this plane intersects the 
    neighbor's path axis, with orientation matching neighbor. """
    neighborPathAxis = Line(neighbor.Pose.t, nhat1)
    originW1 = tangentPlane1.intersectionWithLine(neighborPathAxis)
    # guaranteed to be a point because line is normal to plane
    PoseW1 = SE3.Rt(neighbor.Pose.R, originW1)
    W1 = Waypoint(jointToPlace.numSides, jointToPlace.r, PoseW1, 
                  neighbor.pathIndex())
    toReturn.append(W1)
    
    """ Translate the tangent plane forward by 4r + the new joint's 
    bounding radius, to check if we need a second waypoint. """
    farPoint1 = s1 + nhat1 * (4*jointToPlace.r + jointToPlace.boundingRadius())
    farPlane1 = Plane(farPoint1, nhat1)
    zhatNew = jointToPlace.Pose.R[:,2]
    if backwards:
        nhat2 = -zhatNew
    else:
        nhat2 = zhatNew
    zAxisNew = Line(jointToPlace.Pose.t, zhatNew)
    if farPlane1.sidesOfLine(zAxisNew) == [-1]:
        """ The new joint's Z axis is entirely on the near side of 
        farPlane, so we need to route through a second waypoint. 
        We construct this waypoint similarly to the above, but
        constructing points and planes in the direction of the
        new joint's Z rather than the neighbor's path direction,
        and the waypoint path direction is also along zhatNew. """
        
                
        s2 = ball.c + ball.r*nhat2
        tangentPlane2 = Plane(s2, nhat2)
        # TODO: think through what this is doing and if it's correct given 
        # that I'm not using the separate [a b c] frames
        R_NeighborToW2 = SO3.AngleAxis(np.pi/2, cross(nhat1, nhat2))
        RotationW2 = R_NeighborToW2.R @ neighbor.Pose.R
        originW2 = tangentPlane2.intersectionWithLine(
                                        Line(originW1 + r*nhat1, nhat2))
        PoseW2 = SE3.Rt(RotationW2, originW2)
        W2 = Waypoint(jointToPlace.numSides, r, PoseW2, neighbor.pathIndex())
        toReturn.append(W2)
        
        farPoint2 = s2 + nhat2 * (4*r + jointToPlace.boundingRadius())
        farPlane2 = Plane(farPoint2, nhat2)
        
    else:
        """
        We only need 1 waypoint, but some values of waypoint 2 are used 
        to place+orient jointToPlace, so we need to define those as the
        corresponding values of waypoint 1
        """
        farPlane2 = farPlane1 
        originW2 = originW1
    
    def newPosition(zChange):
        return jointToPlace.Pose.t + zChange*jointToPlace.Pose.R[:,2]
    
    def distanceFromW2(zChange):
        return norm(originW2 - newPosition(zChange))
    
    def signedDistanceToFarPlane2(zChange):
        return farPlane2.signedDistanceToPoint(newPosition(zChange))
    
    beyondFarPlane2 = NonlinearConstraint(signedDistanceToFarPlane2, 
        lb= 0, ub= np.inf)
    
    result = minimize(distanceFromW2, 0, constraints=(beyondFarPlane2))
    zChange = result.x[0]
    jointToPlace.translateAlongZ(zChange)

    if np.dot(jointToPlace.pathDirection(), nhat2) < 0:
        jointToPlace.reversePathDirection()

    toReturn.append(jointToPlace)
    return toReturn
    
