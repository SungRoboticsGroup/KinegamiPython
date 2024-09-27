# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 23:13:27 2023

@author: Daniel Feshbach
"""
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
from pyswarm import pso
import pyswarms as ps
import logging
import collections
import traceback

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
    def __init__(self, root : J, maxAnglePerElbow : float = np.pi/2):
        self.r = root.r
        try:
            self.numSides = root.numSides
        except:
            self.numSides = 4

        self.Joints = [root]
        self.Parents = [-1]     # root has no parent
        self.Links = [LinkCSC(self.r, root.ProximalDubinsFrame(),
                                      root.ProximalDubinsFrame(),
                                      maxAnglePerElbow)]
        assert(maxAnglePerElbow >= 0 and maxAnglePerElbow <= np.pi)
        self.maxAnglePerElbow = maxAnglePerElbow 
        self.boundingBall = root.boundingBall()
        if self.boundingBall.r < self.r:
            self.boundingBall = Ball(root.Pose.t, self.r)
        self.Children = [[]]
    
    def dataDeepCopy(self):
        return copy.deepcopy([self.r, self.numSides, self.Joints, self.Parents, 
                              self.Links, self.maxAnglePerElbow, 
                              self.boundingBall, self.Children])
    
    def setTo(self, data : list):
        self.r, self.numSides, self.Joints, self.Parents, \
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
                 safe : bool = True, endPlane : Plane = None) -> int:
        
        if isinstance(newJoint, OrigamiJoint):
            if newJoint.r != self.r:
                raise ValueError("ERROR: newJoint.r != self.r")
            if newJoint.numSides != self.numSides:
                raise ValueError("ERROR: newJoint.numSides != self.numSides")
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
            if endPlane is None:
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
            self.boundingBall = minBoundingBall(self.boundingBall,
                                                joint.boundingBall())
        for link in self.Links:
            self.boundingBall = minBoundingBall(self.boundingBall,
                                                link.elbow1BoundingBall)
            self.boundingBall = minBoundingBall(self.boundingBall,
                                                link.elbow2BoundingBall)
            
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                  proximalColor='c', centerColor='m', distalColor='y',
                  showJointSurface=True, jointColor=jointColorDefault,
                  jointAxisScale=jointAxisScaleDefault, showJointPoses=True,
                  linkColor=linkColorDefault, surfaceOpacity=surfaceOpacityDefault, showLinkSurface=True, 
                  showLinkPoses=False, showLinkPath=True, pathColor=pathColorDefault,
                  showPathCircles=False, sphereColor=sphereColorDefault,
                  showSpheres=False, showGlobalFrame=False, globalAxisScale=globalAxisScaleDefault, showCollisionBoxes=False, showSpecificCapsules = ([],[]), plotPoint=None, addCapsules = []):
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
                                    surfaceColor=jointColor, surfaceOpacity=surfaceOpacity,
                                    showSurface=showJointSurface, axisScale=jointAxisScale,
                                    showPoses=showJointPoses)
            if not handles is None:
                xyzHandles.append(handles)


        for link in self.Links:
            if link is None:
                continue
            handles = link.addToPlot(ax, color=linkColor, 
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
        
        for jointIndex, capsuleIndex in showCollisionBoxes[0]:
            self.Joints[jointIndex].collisionCapsules[capsuleIndex].box.addToPlot(ax)

        for linkIndex, capsuleIndex in showCollisionBoxes[1]:
            self.Links[linkIndex].collisionCapsules[capsuleIndex].box.addToPlot(ax)
        
        if showSpheres:
            self.boundingBall.addToPlot(ax, color=sphereColor, 
                                        alpha=0.05, frame=True)
        
        if not plotPoint is None:
            ax.scatter(plotPoint[0], plotPoint[1], plotPoint[2], color='black', s=50)
        return np.array(xyzHandles), np.array(abcHandles)
    
    def copyAbbreviatedSelf(self, isolate=False, isolateJoint = None):
        # tree = KinematicTree(Waypoint(self.numSides, self.r, SE3()), self.maxAnglePerElbow)
        # # for i in range(1, len(self.Joints)):
        # #     tree.addJoint(self.Parents[i], self.Joints[i].copy(), relative=False, fixedPosition=True, fixedOrientation=True, safe = False)
        # tree.Joints = [copy.deepcopy(joint) for joint in self.Joints]
        # tree.Links = [copy.deepcopy(link) for link in self.Links]
        # tree.Parents = self.Parents
        # tree.Children = self.Children
        # return tree
        newTree = KinematicTree(copy.deepcopy(self.Joints[0]), self.maxAnglePerElbow)
        if isolate:
            assert(isolateJoint is not None)

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
        for index in toCheck:
            capsules = self.selectCollisionCapsules(specificJointIndices=[index], ignoreLater=ignoreLater)
            numCollisions += self.detectCollisionsWithCapsules([index], capsules, debug=debug)
        
        return numCollisions

        def posesAreSame(pose1, pose2):
            return np.allclose(pose1.t, pose2.t, rtol=1e-05, atol=1e-08)# and np.allclose(pose1.n, pose2.n, rtol=1e-05, atol=1e-08)

        def waypointEdgeCase(l, j):
            joint = self.Joints[j]
            link = self.Links[l]

            closeWaypoints = False

            if isWaypoint(joint) or isWaypoint(self.Joints[l]):
                if self.Parents[l] is not None:
                    parent = self.Joints[self.Parents[l]]
                    if isWaypoint(parent) and posesAreSame(parent.DistalDubinsFrame(), self.Joints[l].ProximalDubinsFrame()) and (waypointEdgeCase(self.Parents[l], j) or (self.Parents[self.Parents[l]] == self.Parents[j])):
                        closeWaypoints = True
                        #print("CLOSE WAYPOINTS link")
                
                if self.Parents[j] is not None:
                    parent = self.Joints[self.Parents[j]]
                    if isWaypoint(parent) and posesAreSame(parent.DistalDubinsFrame(), self.Joints[j].ProximalDubinsFrame()) and (waypointEdgeCase(l, self.Parents[j]) or (self.Parents[self.Parents[j]] == self.Parents[l])):
                        closeWaypoints = True
                        #print("CLOSE WAYPOINTS joint")

                return np.allclose(joint.DistalDubinsFrame().t, link.StartDubinsPose.t, rtol=1e-03, atol=1e-05) or np.allclose(joint.ProximalDubinsFrame().t, link.EndDubinsPose.t, rtol=1e-03, atol=1e-05) or closeWaypoints
            
            return False
            
        def waypointEdgeCase2(wp1, wp2):
            return np.allclose(wp1.DistalDubinsFrame().t, wp2.ProximalDubinsFrame().t, rtol=1e-05, atol=1e-08) or np.allclose(wp1.ProximalDubinsFrame().t, wp2.DistalDubinsFrame().t, rtol=1e-05, atol=1e-08)

        numCollisions = 0
        EPSILON = 0.001

        jointsToCheck = list(range(len(self.Joints))) if specificJointIndices is None else specificJointIndices
        linksToCheck = list(range(len(self.Links))) if specificJointIndices is None else specificJointIndices
        
        others = list(range(0,len(self.Joints)))
        if ignoreLater:
            latest = max(specificJointIndices)
            others = [x for x in others if (x <= latest and self.Joints[x] != None)]
        else:
            others = [x for x in others if self.Joints[x] != None]
        
        start1 = time.time()
        checked = []
        #joint to joint collision:
        for i in jointsToCheck:
            checked.append(i)
            joint = self.Joints[i]
            for j in others:#range(0,len(self.Joints)):
                if j in checked:
                    continue
                joint2 = self.Joints[j]
                #check not joint-link-joint
                if joint != joint2 and not posesAreSame(self.Links[i].StartDubinsPose, joint2.DistalDubinsFrame()) and not posesAreSame(self.Links[j].StartDubinsPose, joint.DistalDubinsFrame()) and not posesAreSame(joint.DistalDubinsFrame(), joint2.ProximalDubinsFrame()) and not posesAreSame(joint2.DistalDubinsFrame(), joint.ProximalDubinsFrame()):
                    jointsCollided = False
                    idx1 = 0
                    for capsule1 in joint.collisionCapsules:
                        if not jointsCollided:
                            idx2 = 0
                            for capsule2 in joint2.collisionCapsules:
                                didCollide, collisionPoint = capsule2.collidesWith(capsule1, includeEnds=includeEnds)
                                if didCollide:
                                    if plot:
                                        self.show(showSpecificCapsules=([(i, idx1),(j, idx2)],[]), plotPoint=collisionPoint)

                                    if debug:
                                        print(f"JOINT JOINT: {i} {j}")
                                    numCollisions += 1
                                    jointsCollided = True
                                    break
                                idx2 += 1
                        idx1 += 1

        #print(f"joint joint time: {time.time() - start1}")
        
        start2 = time.time()
        #joint to link collision
        for i in jointsToCheck:
            joint = self.Joints[i]
            for j in others:#range(0,len(self.Links)):
                link = self.Links[j]
                #check link and joint not connected
                if not posesAreSame(link.StartDubinsPose, joint.DistalDubinsFrame()) and not posesAreSame(link.EndDubinsPose, joint.ProximalDubinsFrame()) and not waypointEdgeCase(j, i):
                    collided = False
                    idx1 = 0
                    for capsule1 in joint.collisionCapsules:
                        if not collided:
                            idx2 = 0
                            for capsule2 in link.collisionCapsules:
                                didCollide, collisionPoint = capsule2.collidesWith(capsule1, includeEnds=includeEnds)
                                if didCollide:
                                    if plot:
                                        self.show(showSpecificCapsules=([(i, idx1)], [(j, idx2)]), showCollisionBoxes = False, plotPoint = collisionPoint)
                                    numCollisions += 1
                                    if debug:
                                        print(f"JOINT LINK {i} {j}")
                                    collided = True
                                    break
                                idx2 += 1
                        idx1 += 1
        
        #print(f"joint link time: {time.time() - start2}")

        start3 = time.time()
        #link to joint collision
        for j in linksToCheck:
            link = self.Links[j]
            for i in others:#range(0,len(self.Joints)):
                if i in jointsToCheck:
                    continue
                joint = self.Joints[i]
                #check link and joint not connected
                if not posesAreSame(link.StartDubinsPose, joint.DistalDubinsFrame()) and not posesAreSame(link.EndDubinsPose, joint.ProximalDubinsFrame()) and not waypointEdgeCase(j, i):
                    collided = False
                    idx1 = 0
                    for capsule1 in joint.collisionCapsules:
                        if not collided:
                            idx2 = 0
                            for capsule2 in link.collisionCapsules:
                                didCollide, collisionPoint = capsule2.collidesWith(capsule1, includeEnds=includeEnds)
                                if didCollide:
                                    if plot:
                                        self.show(showSpecificCapsules=([(i, idx1)], [(j, idx2)]), showCollisionBoxes = False, plotPoint = collisionPoint)
                                    numCollisions += 1
                                    if debug:
                                        print(f"LINK JOINT {i} {j}")
                                    collided = True
                                    break
                                idx2 += 1
                        idx1 += 1
        #print(f"link joint time: {time.time() - start3}")

        start4 = time.time()
        checked = []
        #link to link collision
        for i in linksToCheck:
            checked.append(i)
            link = self.Links[i]
            for j in others:#range(0,len(self.Links)):
                if j in checked:
                    continue
                link2 = self.Links[j]
                #check that not part of the same branch and not link-joint-link
                if not posesAreSame(link.StartDubinsPose, link2.StartDubinsPose) and not posesAreSame(link.EndDubinsPose, link2.StartDubinsPose) and not posesAreSame(link2.EndDubinsPose, link.StartDubinsPose) and not posesAreSame(self.Joints[i].DistalDubinsFrame(), link2.StartDubinsPose) and not posesAreSame(self.Joints[j].DistalDubinsFrame(), link.StartDubinsPose):
                    linksCollided = False
                    idx1 = 0
                    for capsule1 in link.collisionCapsules:
                        if not linksCollided:
                            idx2 = 0
                            for capsule2 in link2.collisionCapsules:
                                didCollide, collisionPoint = capsule2.collidesWith(capsule1, includeEnds=includeEnds)
                                if didCollide and ((not isWaypoint(self.Joints[i]) and not isWaypoint(self.Joints[j]) and not waypointEdgeCase2(self.Joints[i], self.Joints[j])) or not ignoreWaypoints):
                                    if plot:
                                        self.show(showSpecificCapsules=([], [(i, idx1), (j, idx2)]), showCollisionBoxes = False, plotPoint = collisionPoint)
                                    numCollisions += 1

                                    if debug:
                                        print(f"LINK LINK {i} {j}")
                                    linksCollided = True
                                    break
                                idx2 += 1
                        idx1 += 1
                
                elif link2 != link and not isWaypoint(self.Joints[i]) and not isWaypoint(self.Joints[j]) and posesAreSame(link.StartDubinsPose, link2.StartDubinsPose):
                    idx1 = 0
                    for capsule in link.collisionCapsules:
                        didOverlap, collisionPoint = capsule.frameOverlap(link2.EndDubinsPose, self.r)
                        if didOverlap:
                            numCollisions += 1
                            if debug:
                                print(f"LINK FRAME OVERLAP {i} {j}")
                            if plot:
                                self.show(showSpecificCapsules=([], [(i, idx1), (j, len(link2.collisionCapsules) - 1)]), plotPoint = collisionPoint)
                            break
                        idx1 += 1
        #print(f"link link time: {time.time() - start4}")
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

    
    def branchingParametersFrom(self, parentIndex : int):
        linksToChildren = [self.Links[childIndex] for childIndex in self.Children[parentIndex]]
        return [link.branchingParameters() for link in linksToChildren]
    
    def exportLink3DFile(self, parentIndex : int, folder : str, fileFormat = "stl", pose=False):

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

        os.system(f"openscad -q -o 3d_output/{folder}/{name} scad_output/{folder}/{name}.scad")

        return f"3d_output/{folder}/{name}"

    def export3DKinematicTree(self, folder = "", fileFormat = "stl"):
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
                tree.exportLink3DFile(i,folder,fileFormat)
            print(f"Finished link {i}/{len(tree.Children) - 1}, Time: {time.time() - start} \r")
        
        #export all the joints
        for i in range(0,len(tree.Joints)):
            start = time.time()
            if not isinstance(tree.Joints[i],PrintedWaypoint):
                tree.Joints[i].export3DFile(i,folder,fileFormat)
            print(f"Finished joint {i}/{len(tree.Joints) - 1}, Time: {time.time() - start} \r")
        

    def show(self, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             showJointSurface=True, jointColor=jointColorDefault, 
             jointAxisScale=jointAxisScaleDefault, showJointPoses=True,
             linkColor=linkColorDefault, surfaceOpacity=surfaceOpacityDefault, showLinkSurface=True, 
             showLinkPoses=False, showLinkPath=True, pathColor=pathColorDefault,
             showPathCircles=False, sphereColor=sphereColorDefault,
             showSpheres=False, block=blockDefault, showAxisGrids=False, 
             showGlobalFrame=False, globalAxisScale=globalAxisScaleDefault,
             showGroundPlane=False, groundPlaneScale=groundPlaneScaleDefault,
             groundPlaneColor=groundPlaneColorDefault, showCollisionBoxes=([],[]), showSpecificCapsules=([],[]), plotPoint = None, addCapsules=[]):
        ax = plt.figure().add_subplot(projection='3d')
        if showGroundPlane: #https://stackoverflow.com/questions/36060933/plot-a-plane-and-points-in-3d-simultaneously
            xx, yy = np.meshgrid(range(groundPlaneScale), range(groundPlaneScale))
            xx = xx - groundPlaneScale/2
            yy = yy - groundPlaneScale/2
            z = 0*xx #(9 - xx - yy) / 2 

            # plot the plane
            ax.plot_surface(xx, yy, z, alpha=surfaceOpacity/4, color=groundPlaneColor)
        xyzHandles, abcHandles = self.addToPlot(ax, xColor, yColor, zColor, 
                                    proximalColor, centerColor, distalColor,
                                    showJointSurface, jointColor, 
                                    jointAxisScale, showJointPoses,
                                    linkColor, surfaceOpacity, showLinkSurface, 
                                    showLinkPoses, showLinkPath, pathColor,
                                    showPathCircles, sphereColor,
                                    showSpheres, showGlobalFrame, globalAxisScale, showCollisionBoxes=showCollisionBoxes, showSpecificCapsules=showSpecificCapsules, plotPoint=plotPoint, addCapsules=addCapsules)
        
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
        plt.show(block=block)


    def transformAll(self, Transformation : SE3):
        self.transformJoint(0, Transformation, safe=False)

    """ 
    Apply given transformation (SE3() object) to given joint (index), 
    and to its descendants if propogate (defaults to True).

    Returns True if it succeeds (the transformation gives valid links).
    In safe=True mode (default), if it fails it will leave the chain unchanged,
    print a warning, and return False rather than throwing an error.   
    """
    def transformJoint(self, jointIndex : int, Transformation : SE3, 
                       propogate : bool = True, recomputeBoundingBall=True,
                       recomputeLinkPath : bool = True, 
                       safe : bool = True, relative : bool = False) -> bool:
        if relative:
            Transformation = self.Joints[jointIndex].Pose @ Transformation @ self.Joints[jointIndex].Pose.inv()
        
        if safe:
            backup = self.dataDeepCopy()
            try:
                self.transformJoint(jointIndex, Transformation, 
                       propogate, recomputeBoundingBall,
                       recomputeLinkPath, safe=False, relative=False)
            except ValueError as err:
                print("WARNING: something went wrong in transformJoint:")
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

        return True
                
    def setJointState(self, jointIndex : int, newState : float) -> bool:
        joint = self.Joints[jointIndex]
        minState, maxState = joint.stateRange()
        if newState < minState or newState > maxState:
            print("WARNING: state out of range in setJointRange, "+
                    "state unchanged.")
            return False
        Transformation = joint.TransformStateTo(newState)
        for c in self.Children[jointIndex]:
            # TODO: it is possible, and would be more efficient, to make this 
            # transform the existing links rather than recompute them
            self.transformJoint(c, Transformation, propogate=True, recomputeLinkPath=True,
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

    def optimizeJointPlacement(self, index, maxiter, tol, collisionError, showGuess = False, ignorePlacement=False, ignoreLater = False):
        parentIndex = self.Parents[index]

        selectedIndices = [index] if ignorePlacement else ([index] + self.Children[index])
        selectedCapsules = self.selectCollisionCapsules(specificJointIndices=selectedIndices, ignoreLater=ignoreLater)
        def linkLoss(t, link, curveLossFactor = 2):
            #d = (np.abs(t.Links[index].path.theta1 * curveLossFactor) ** 3 + np.abs(t.Links[index].path.theta2 * curveLossFactor) ** 3) + (np.arccos(np.clip((np.trace(t.Joints[index].ProximalDubinsFrame().R.T @ t.Joints[t.Parents[index]].DistalDubinsFrame().R) - 1) / 2, -1.0, 1.0))) * (1/t.Links[index].path.length + 1)
            #return (t.Links[index].path.length - t.r*2) ** 2  + d + t.detectCollisionsWithCapsules(selectedIndices, selectedCapsules) * collisionError #t.detectCollisions(specificJointIndices=[index] if ignorePlacement else ([index] + self.Children[index]), includeEnds=ignorePlacement, ignoreLater=ignoreLater) * collisionError
            d = t.Links[index].path.theta1 ** 2 + t.Links[index].path.theta2 ** 2
            return t.Links[index].path.length ** 2 + t.detectCollisionsWithCapsules(selectedIndices, selectedCapsules) * collisionError + d * t.r
        def objective(params, returnWhich = False):
            tree = self.copyAbbreviatedSelf(ignoreLater, index)

            translation = params[0]
            rotation = params[1]

            transform = SE3.Trans([0,0,translation]) @ SE3.Rz(rotation)

            linkLossUnchanged = collisionError * len(self.Joints) * (len(self.Children) + 1)
            linkLossReversedZhat = collisionError * len(self.Joints) * (len(self.Children) + 1)
            linkLossReversedParent = collisionError * len(self.Joints) * (len(self.Children) + 1)

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
                    tree2 = self.copyAbbreviatedSelf(ignoreLater, index)
                    tree2.Joints[index].reverseZhat()
                    if not tree2.transformJoint(index, SE3.Trans([0,0,-translation]) @ SE3.Rz(-rotation), safe=False, relative=True, propogate=ignorePlacement, recomputeLinkPath=True, recomputeBoundingBall=False):
                        raise Exception()

                    linkLossReversedZhat = linkLoss(tree2, index)
                except:
                    pass
            # TODO: does this need to happen? affects runtime significantly
            # try:
            #     tree3 = self.copyAbbreviatedSelf()
            #     tree3.Joints[parentIndex].reverseZhat()
            #     if not tree3.transformJoint(parentIndex, SE3.Trans([0,0,0]), safe=False, relative=True, propogate=False, recomputeLinkPath=True):
            #         raise Exception()
            #     if not tree3.transformJoint(index, SE3.Trans([0,0,translation]) @ SE3.Rz(rotation), propogate=ignorePlacement, safe=False, relative=True) or tree3.detectCollisions(specificJointIndices=[parentIndex, index], ignoreLater=ignoreLater):
            #         raise Exception()
                
            #     linkLossReversedParent = linkLoss(tree3, index)
            # except:
            #     pass

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

        joint = self.Joints[index]
        parent = self.Joints[self.Parents[index]]

        frame1 = joint.Pose
        frame2 = parent.DistalDubinsFrame()

        transformation = frame1.inv() * frame2

        #try to make initial guess right next to each other
        initialPosition = transformation.t[2]
        initialRotation = np.arctan2(transformation.R[1, 0], transformation.R[0, 0])
        initialGuess = [initialPosition,initialRotation]
        initialTree = tree = self.copyAbbreviatedSelf(ignoreLater, index)
        try:
            if not initialTree.transformJoint(index, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation), propogate=ignorePlacement, safe=False, relative=True, recomputeBoundingBall=False):
                raise Exception()
        except:
            initialGuess = [0,0]

        initialLoss = objective([0,0])

        #calculate dist
        current = index
        dist = 0
        while current != 0:
            dist += self.Links[current].path.length
            current = self.Parents[current]
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
        noise[:, 0] = np.linspace(-dist*2,dist*2, num=n_particles - 2)
        noise[:, 1] = np.random.uniform(-np.pi*2,np.pi*2, n_particles - 2)
        init_pos[2:] += noise
        init_pos[2:, 0] = np.clip(init_pos[2:, 0], -dist*2, dist*2)
        init_pos[2:, 1] = np.clip(init_pos[2:, 1], -np.pi*2, np.pi*2)

        optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,dimensions=2,options={'c1':0.6, 'c2':0.7, 'w':0.5},bounds=(np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])),init_pos=init_pos,ftol=tol)
        minSwarmLoss, minSwarmResult = optimizer.optimize(joint_batch_objective_function, iters=int((maxiter + 1)/2),verbose=False, n_processes=n_particles)

        nelderMead = minimize(objective, minSwarmResult, method="Nelder-Mead", bounds=bounds, tol=tol, options={
            'maxiter':int(maxiter/2),
            'fatol':tol,
        })
        result = nelderMead.x
        loss = nelderMead.fun

        #print(minSwarmLoss, loss)

        print(f"Optimized joint {index} in {time.time() - start}s -- Old loss: {initialLoss}, Improved Loss: {loss}")

        which = objective(result, returnWhich=True)

        tree = self.copyAbbreviatedSelf()
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
                tree2 = self.copyAbbreviatedSelf()
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

    def optimizeWaypointPlacement(self, index, maxiter, tol, collisionError, ignorePlacement=False, ignoreLater=False):
        current = index
        dist = 0
        while current != 0:
            dist += self.Links[current].path.length
            current = self.Parents[current]
        bounds = [(-dist*2, dist*2)]*3 + [(-np.pi*2, np.pi*2)] * 3

        #make initial guess as close as possible to previous
        start = time.time()
        initialTree = self.copyAbbreviatedSelf(ignoreLater, index)

        initialGuess = [0]*6

        parent = initialTree.Joints[initialTree.Parents[index]]
        waypoint = initialTree.Joints[index]

        transform = parent.DistalDubinsFrame() * waypoint.ProximalDubinsFrame().inv()
        initialGuess[0:3] = transform.t
        initialGuess[3:6] = SE3.Rt(transform.R, np.zeros(3)).eul()

        selectedIndices = [index] if ignorePlacement else ([index] + self.Children[index])
        selectedCapsules = self.selectCollisionCapsules(specificJointIndices=selectedIndices, ignoreLater=ignoreLater)

        def linkLoss(t, link, curveLossFactor = np.pi):
            d = t.Links[index].path.theta1 ** 2 + t.Links[index].path.theta2 ** 2
            return t.Links[index].path.length ** 2  + t.detectCollisionsWithCapsules(selectedIndices, selectedCapsules) * collisionError + d

        def objective(params):
            tree = self.copyAbbreviatedSelf(ignoreLater, index)

            try:
                if not tree.transformJoint(index, SE3.Trans(params[0:3]) @ SE3.Rz(params[3]) @ SE3.Ry(params[4]) @ SE3.Rz(params[5]),  propogate=ignorePlacement, safe=False, relative=False, recomputeBoundingBall=False):
                    raise Exception()
            except:
                # tb = traceback.format_exc()
                # print(f"Traceback details:\n{tb}")
                return collisionError * len(self.Joints) * (len(self.Children) + 1)
            
            return linkLoss(tree, index) + np.linalg.norm(np.array(params[3:6]) - SE3.Rt(transform.R, np.zeros(3)).eul()) * 10

        try:
            if not initialTree.transformJoint(index, SE3.Trans(initialGuess[0:3]) @ SE3.Rz(initialGuess[3]) @ SE3.Ry(initialGuess[4]) @ SE3.Rz(initialGuess[5]),  propogate=ignorePlacement, safe=False, relative=False, recomputeBoundingBall=False):
                raise Exception("failed")

            print(f"INITAL WAYPOINT GUESS LOSS: {objective(initialGuess)}")

            #initialTree.detectCollisions(debug=True)
        except Exception as e:
            try:
                success = False
                for i in range(0, len(self.Children[index])):
                    if success:
                        break
                    initialTree = self.copyAbbreviatedSelf(ignoreLater, index)
                    rots = SE3.Rt(initialTree.Joints[self.Children[index][i]].ProximalDubinsFrame().R, np.zeros(3)).eul()
                    frame1 = SE3.Trans(initialTree.Joints[self.Parents[index]].DistalDubinsFrame().t) @ SE3.Rz(rots[0]) @ SE3.Ry(rots[1]) @ SE3.Rz(rots[2])
                    transform2 = frame1 * initialTree.Joints[index].ProximalDubinsFrame().inv()
                    initialGuess[0:3] = transform2.t
                    initialGuess[3:6] = SE3.Rt(transform2.R, np.zeros(3)).eul()
                    try:
                        if initialTree.transformJoint(index, SE3.Trans(initialGuess[0:3]) @ SE3.Rz(initialGuess[3]) @ SE3.Ry(initialGuess[4]) @ SE3.Rz(initialGuess[5]),  propogate=ignorePlacement, safe=False, relative=False, recomputeBoundingBall=False):
                            #print(initialTree.Joints[index].ProximalDubinsFrame().t, initialTree.Joints[initialTree.Parents[index]].DistalDubinsFrame().t)
                            success = True
                    except:
                        pass

                if not success:
                    raise Exception("failed")

                print(f"NVM: {e}, INITAL WAYPOINT GUESS LOSS: {objective(initialGuess)}")
            except Exception as e:
                print(f"INITIAL WAYPOINT GUESS NOT POSSIBLE : {e}")
                initialGuess = [0]*6
            
        initialLoss = objective([0]*6)

        global waypoint_batch_objective_function
        def waypoint_batch_objective_function(X):
            return np.array([objective(x) for x in X])

        n_particles = 24

        init_pos = np.tile(np.array(initialGuess, dtype='float64'), (n_particles,1))
        init_pos[1] = np.array([0]*6)
        #add random noise
        noise = np.zeros_like(init_pos[2:])
        for i in range(0,3):
            noise[:, i] = np.linspace(-dist*2,dist*2, num=n_particles - 2)
        for i in range(3, 6):
            noise[:, i] = np.random.uniform(-np.pi*2,np.pi*2, n_particles - 2)
        init_pos[2:] += noise
        for i in range(0,3):
            init_pos[2:, i] = np.clip(init_pos[2:, 0], -dist*2, dist*2)
        for i in range(3, 6):
            init_pos[2:, i] = np.clip(init_pos[2:, 1], -np.pi*2, np.pi*2)
        optimizer = ps.single.GlobalBestPSO(n_particles=n_particles,dimensions=6,options={'c1':0.7, 'c2':0.5, 'w':0.5},bounds=(np.array([b[0] for b in bounds]), np.array([b[1] for b in bounds])),init_pos=init_pos,ftol=tol)
        minSwarmLoss, minSwarmResult = optimizer.optimize(waypoint_batch_objective_function, iters=maxiter,verbose=False, n_processes=n_particles)
        
        tree = self.copyAbbreviatedSelf()
        if tree.transformJoint(index, SE3.Trans(minSwarmResult[0:3]) @ SE3.Rz(minSwarmResult[3]) @ SE3.Ry(minSwarmResult[4]) @ SE3.Rz(minSwarmResult[5]),  propogate=False, safe=False, relative=False, recomputeBoundingBall=False):
            print(f"Optimized waypoint {index} in {time.time() - start}s -- Old Loss: {initialLoss}, Improved Loss: {minSwarmLoss}")
            return tree, minSwarmResult
        else:
            raise Exception("Optimization failed dramatically")

    def squaredOptimize(self, showSteps=False, guarantee=False):
        for i in range(0, len(self.Joints)):
            self.Joints[i].recomputeCollisionCapsules()

        if self.detectCollisions(debug=True) > 0:
            print("Warning: Initial tree contains collisions.")
        if showSteps and isinstance(self.Joints[0], OrigamiJoint):
            self.show()

        collisionError = 0
        for i in range(0, len(self.Joints)):
            if len(self.Children[i]) > 0:
                continue
            j = i
            length = 0
            while j != 0:
                length += self.Links[j].path.length ** 2
                j = self.Parents[j]
            if length > collisionError:
                collisionError = length

        print(f"Collision error is {collisionError}")

        start = time.time()

        tree = self.copyAbbreviatedSelf()
        isOptimized = [True] + [False] * (len(self.Joints) - 1) #isOptimized[i] is True if joint i is optimized
        numOptimized = 1

        def optimizeFromIndex(index):       
            nonlocal tree
            nonlocal numOptimized

            iters = 20
            tolerance = self.r/10
            if isOptimized[self.Parents[index]]:
                iters = 50
                tolerance = self.r/10

            if isWaypoint(self.Joints[index]):
                tree, loss = tree.optimizeWaypointPlacement(index, maxiter=iters, tol=tolerance, collisionError=collisionError, ignoreLater = (not guarantee))
            else:
                tree, loss = tree.optimizeJointPlacement(index, maxiter=iters, tol=tolerance, collisionError=collisionError, ignoreLater = (not guarantee))

            if isOptimized[self.Parents[index]]:
                isOptimized[index] = True
                numOptimized += 1
            else:
                optimizeFromIndex(self.Parents[index])


        while numOptimized < len(self.Joints):
            i = None
            for j in range(len(self.Joints) - 1, 0, -1):
                if not isOptimized[j] and len(self.Children[j]) == 0:
                    i = j

            if i == None:
                print("SOMETHING TERRIBLY WRONG")
                break

            print(f"OPTIMIZING CHAIN ENDING AT {i}:")
            start2 = time.time()
            order = []

            parent = i
            while not isOptimized[parent]:
                order.append(parent)
                parent = self.Parents[parent]
            order.reverse()

            optimizeStreak = 0

            for j in range(0, len(order)):
                iters = 50
                tolerance = self.r/10
                
                if isWaypoint(self.Joints[order[j]]):
                    try:
                        tree2, loss = tree.optimizeWaypointPlacement(order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, ignorePlacement=True, ignoreLater = (not guarantee))

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
                        # for idx in order:
                        #     if isOptimized[idx]:
                        #         isOptimized[idx] = False
                        #         numOptimized -= 1

                        tree, loss = tree.optimizeWaypointPlacement(order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, ignorePlacement=False, ignoreLater = (not guarantee))
                        print(tree.detectCollisions(specificJointIndices=[order[j]], ignoreLater=(not guarantee), plot=False, debug=True))
                        #tree.show()
                        break
                else:
                    try:
                        tree2, loss = tree.optimizeJointPlacement(order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, showGuess=False, ignorePlacement=True, ignoreLater = (not guarantee))

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
                        # for idx in order:
                        #     if isOptimized[idx]:
                        #         isOptimized[idx] = False
                        #         numOptimized -= 1

                        tree, loss = tree.optimizeJointPlacement(order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, showGuess=False, ignorePlacement=True, ignoreLater = (not guarantee))
                        print(tree.detectCollisions(specificJointIndices=[order[j]], ignoreLater=(not guarantee), plot=False, debug=True))
                        break
                
            #tree.show()

            while not isOptimized[i]:
                optimizeFromIndex(i)
            
            print("CURRENT COLLISIONS")
            if tree.detectCollisions(debug=True) == 0:
                print("NONE")
            else:
                pass
                #tree.show()

            # print("")
            # Another pass, not totally necessary but helps streamline shape, commented out to improve runtime
            # for j in range(0, len(order)):
            #     iters = 50
            #     tolerance = self.r/10
                
            #     if isWaypoint(self.Joints[order[j]]):
            #         tree, loss = tree.optimizeWaypointPlacement(order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, ignoreLater = (not guarantee))
            #     else:
            #         tree, loss = tree.optimizeJointPlacement(order[j], maxiter=iters, tol=tolerance, collisionError=collisionError, showGuess=False, ignoreLater = (not guarantee))

            # print("CURRENT COLLISIONS (p2):")
            # if tree.detectCollisions(debug=True) == 0:
            #     print("NONE")

            print(f"Optimized chain ending at {i} in {time.time() - start2}s \n")

        print(f"TOTAL OPTIMIZATION TIME: {time.time() - start}")

        return tree
    
    def save(self, filename: str):
        #TODO: ADD EXTENSIONS FOR PRINTED JOINTS
        with open(f"save/{filename}.tree", "w") as f:
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
        with open(f"save/{filename}.tree") as f:
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
        raise Exception(f"file save/{filename}.tree doesnt exist")

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
    
