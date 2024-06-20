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
from scipy.optimize import dual_annealing

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
                  showSpheres=False, showGlobalFrame=False, globalAxisScale=globalAxisScaleDefault, showCollisionBoxes=False, showSpecificCapsules = ([],[]), plotPoint=None):
        xyzHandles = []
        abcHandles = []
        
        if showGlobalFrame:
            handles = addPosesToPlot(SE3(), ax, globalAxisScale, xColor, yColor, zColor)
            if not handles is None:
                xyzHandles.append(handles)
        
        for joint in self.Joints:
            handles = joint.addToPlot(ax, xColor, yColor, zColor, 
                                    proximalColor, centerColor, distalColor, 
                                    sphereColor=sphereColor, showSphere=showSpheres, 
                                    surfaceColor=jointColor, surfaceOpacity=surfaceOpacity,
                                    showSurface=showJointSurface, axisScale=jointAxisScale,
                                    showPoses=showJointPoses)
            if not handles is None:
                xyzHandles.append(handles)

            if showCollisionBoxes:
                for capsule in joint.collisionCapsules():
                    capsule.addToPlot(ax)

        for link in self.Links:
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
            
            if showCollisionBoxes:
                for capsule in link.collisionCapsules():
                    capsule.addToPlot(ax)
        
        for jointIndex, capsuleIndex in showSpecificCapsules[0]:
            self.Joints[jointIndex].collisionCapsules()[capsuleIndex].addToPlot(ax)

        for linkIndex, capsuleIndex in showSpecificCapsules[1]:
            self.Links[linkIndex].collisionCapsules()[capsuleIndex].addToPlot(ax)
        
        if showSpheres:
            self.boundingBall.addToPlot(ax, color=sphereColor, 
                                        alpha=0.05, frame=True)
        
        if not plotPoint is None:
            ax.scatter(plotPoint[0], plotPoint[1], plotPoint[2], color='black', s=50)
        return np.array(xyzHandles), np.array(abcHandles)
    
    def detectCollisions(self, specificJointIndex = None, plot=False):
        #TODO: check if base circles intersect within the same branch
        def posesAreSame(pose1, pose2):
            return np.allclose(pose1.t, pose2.t, rtol=1e-05, atol=1e-08) and np.allclose(pose1.n, pose2.n, rtol=1e-05, atol=1e-08)

        numCollisions = 0
        EPSILON = 0.001

        jointsToCheck = list(range(len(self.Joints))) if specificJointIndex is None else [specificJointIndex]

        linksToCheck = list(range(len(self.Links))) if specificJointIndex is None else [specificJointIndex] + [child for child in self.Children[specificJointIndex]]
        
        #joint to joint collision:
        for i in jointsToCheck:
            joint = self.Joints[i]
            for j in range(0, len(self.Joints)):
                joint2 = self.Joints[j]
                #check not joint-link-joint
                if joint != joint2 and not posesAreSame(self.Links[i].StartDubinsPose, joint2.DistalDubinsFrame()) and not posesAreSame(self.Links[j].StartDubinsPose, joint.DistalDubinsFrame()) and not posesAreSame(joint.DistalDubinsFrame(), joint2.ProximalDubinsFrame()) and not posesAreSame(joint2.DistalDubinsFrame(), joint.ProximalDubinsFrame()):
                    jointsCollided = False
                    idx1 = 0
                    for capsule1 in joint.collisionCapsules():
                        if not jointsCollided:
                            idx2 = 0
                            for capsule2 in joint2.collisionCapsules():
                                didCollide, collisionPoint = capsule2.collidesWith(capsule1)
                                if didCollide:
                                    if plot:
                                        self.show(showSpecificCapsules=([(i, idx1),(j, idx2)],[]), showCollisionBoxes=False, plotPoint=collisionPoint)
                                    numCollisions += 1
                                    jointsCollided = True
                                    break
                                idx2 += 1
                        idx1 += 1

        #joint to link collision
        for i in jointsToCheck:
            joint = self.Joints[i]
            for j in range(0,len(self.Links)):
                link = self.Links[j]
                #check link and joint not connected
                if not posesAreSame(link.StartDubinsPose, joint.DistalDubinsFrame()) and not posesAreSame(link.EndDubinsPose, joint.ProximalDubinsFrame()):
                    collided = False
                    idx1 = 0
                    for capsule1 in joint.collisionCapsules():
                        if not collided:
                            idx2 = 0
                            for capsule2 in link.collisionCapsules():
                                didCollide, collisionPoint = capsule2.collidesWith(capsule1)
                                if didCollide:
                                    if plot:
                                        self.show(showSpecificCapsules=([(i, idx1)], [(j, idx2)]), showCollisionBoxes = False, plotPoint = collisionPoint)
                                    numCollisions += 1
                                    collided = True
                                    break
                                idx2 += 1
                        idx1 += 1

        #link to link collision
        for i in linksToCheck:
            link = self.Links[i]
            for j in range(0,len(self.Links)):
                link2 = self.Links[j]
                #check that not part of the same branch and not link-joint-link
                if not posesAreSame(link.StartDubinsPose, link2.StartDubinsPose) and not posesAreSame(link.EndDubinsPose, link2.StartDubinsPose) and not posesAreSame(link2.EndDubinsPose, link.StartDubinsPose) and not posesAreSame(self.Joints[i].DistalDubinsFrame(), link2.StartDubinsPose) and not posesAreSame(self.Joints[j].DistalDubinsFrame(), link.StartDubinsPose):
                    linksCollided = False
                    idx1 = 0
                    for capsule1 in link.collisionCapsules():
                        if not linksCollided:
                            idx2 = 0
                            for capsule2 in link2.collisionCapsules():
                                didCollide, collisionPoint = capsule2.collidesWith(capsule1)
                                if didCollide:
                                    if plot:
                                        self.show(showSpecificCapsules=([], [(i, idx1), (j, idx2)]), showCollisionBoxes = False, plotPoint = collisionPoint)
                                    numCollisions += 1
                                    linksCollided = True
                                    break
                                idx2 += 1
                        idx1 += 1
                
                elif link2 != link and not isinstance(self.Joints[i], Waypoint) and not isinstance(self.Joints[j], Waypoint) and posesAreSame(link.StartDubinsPose, link2.StartDubinsPose):
                    idx1 = 0
                    for capsule in link.collisionCapsules():
                        didOverlap, collisionPoint = capsule.frameOverlap(link2.EndDubinsPose, self.r)
                        if didOverlap:
                            numCollisions += 1
                            if plot:
                                self.show(showSpecificCapsules=([], [(i, idx1), (j, len(link2.collisionCapsules()) - 1)]), showCollisionBoxes=False, plotPoint = collisionPoint)
                            break
                        idx1 += 1
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
        params = np.round(self.branchingParametersFrom(parentIndex), 4)

        #check link min length        
        for i in self.Children[parentIndex]:
            if self.Links[i].path.length < (self.Joints[i].screwRadius*2 + self.Joints[i].printParameters.holeMargin*2 + self.Joints[i].printParameters.gridHoleRadius):
                if (np.abs(self.Links[i].path.theta1) > 0.0001 or np.abs(self.Links[i].path.theta2) > 0.0001): #replace with some epsilon
                    raise InvalidDimensionException(f"The link between joint {parentIndex} and joint {i} are too close together for a link to be created between them.")
                else:
                    print(f"Skipping link between joint {parentIndex} and joint {i}. Extending joint {parentIndex} by {self.Links[i].path.tMag} instead")
                    source.extendSegment(self.Links[i].path.tMag)

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

        #export all the links
        for i in range(0,len(self.Children)):
            start = time.time()
            if len(self.Children[i]) > 0:
                self.exportLink3DFile(i,folder,fileFormat)
            print(f"Finished link {i}/{len(self.Children) - 1}, Time: {time.time() - start} \r")
        
        #export all the joints
        for i in range(0,len(self.Joints)):
            start = time.time()
            if not isinstance(self.Joints[i],PrintedWaypoint):
                self.Joints[i].export3DFile(i,folder,fileFormat)
            print(f"Finished joint {i}/{len(self.Joints) - 1}, Time: {time.time() - start} \r")
        

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
             groundPlaneColor=groundPlaneColorDefault, showCollisionBoxes=False, showSpecificCapsules=([],[]), plotPoint = None):
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
                                    showSpheres, showGlobalFrame, globalAxisScale, showCollisionBoxes=showCollisionBoxes, showSpecificCapsules=showSpecificCapsules, plotPoint=plotPoint)
        
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
            else:
                for c in self.Children[jointIndex]:
                    child = self.Joints[c]
                    self.Links[c] = LinkCSC(self.r, joint.DistalDubinsFrame(), 
                                            child.ProximalDubinsFrame(),
                                            self.maxAnglePerElbow)
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
        return True
    
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

    def optimizeWaypointsAndJointPlacement(self, waypoint1Index, waypoint2Index, jointIndex, maxiter=3):
        start = time.time()

        dist = self.Links[waypoint1Index].path.length + self.Links[waypoint2Index].path.length + self.Links[jointIndex].path.length

        #waypoint 1 t, waypoint 2 t, joint translation, waypoint 1 R, waypoint 2 R, joint rotation
        bounds = [(-dist*2, dist*2)]*7 + [(-np.pi*2, np.pi*2)] * 7

        def transformJointToPose(tree, index, frame2):
            joint = tree.Joints[index]
            frame1 = joint.ProximalDubinsFrame()
            
            transformation = frame2 * frame1.inv()

            try:
                if not tree.transformJoint(index, transformation, propogate=False, safe=False, relative=False):
                    transformJointToPose(index, frame1)
                    return False

                if tree.detectCollisions() > 0:
                    transformJointToPose(index, frame1)
                    return False
            except:
                transformJointToPose(index, frame1)
                return False
            
            return True

        def objective(params):            
            tree = copy.deepcopy(self)
            
            t1 = params[0:3]
            t2 = params[3:6]
            t3 = params[6]
            rx1 = params[7]
            ry1 = params[8]
            rz1 = params[9]
            rx2 = params[10]
            ry2 = params[11]
            rz2 = params[12]
            r3 = params[13]

            start2 = time.time()

            try:
                if tree.transformJoint(waypoint1Index, SE3.Trans(t1) @ SE3.Rz(rx1) @ SE3.Ry(ry1) @ SE3.Rz(rz1),  propogate=False, safe=False, relative=False) and tree.transformJoint(waypoint2Index, SE3.Trans(t2) @ SE3.Rz(rx2) @ SE3.Ry(ry2) @ SE3.Rz(rz2), propogate=False, relative=False) and tree.transformJoint(jointIndex, SE3.Trans([0,0,t3]) @ SE3.Rz(r3), propogate=False, relative=True):
                    distance = tree.Links[waypoint1Index].path.length + tree.Links[waypoint2Index].path.length + tree.Links[jointIndex].path.length
                    #np.linalg.norm(tree.Joints[jointIndex].ProximalDubinsFrame().t - tree.Joints[tree.Parents[waypoint1Index]].DistalDubinsFrame().t)
                    ans = distance + 1000*(tree.detectCollisions(specificJointIndex=waypoint1Index) + tree.detectCollisions(specificJointIndex=waypoint2Index) + tree.detectCollisions(specificJointIndex=jointIndex))
                else:
                    ans = 100000
            except:
                ans = 100000
            
            #print(time.time() - start2)
            return ans

        #try to make as close as possible
        initialTree = copy.deepcopy(self)

        initialGuess = [0]*14

        parent = initialTree.Joints[initialTree.Parents[waypoint1Index]]
        joint = initialTree.Joints[jointIndex]
        waypoint1 = initialTree.Joints[waypoint1Index]
        waypoint2 = initialTree.Joints[waypoint2Index]
        transform1 = parent.DistalDubinsFrame() * waypoint1.ProximalDubinsFrame().inv()
        initialGuess[0:3] = transform1.t
        initialGuess[7:10] = SE3.Rt(transform1.R, np.zeros(3)).eul()
            
        try:
            initialTree.transformJoint(waypoint1Index, transform1, propogate=False, safe=False, relative=False)
            if (initialTree.detectCollisions() > 0):
                initialTree.transformJoint(waypoint1Index, transform1.inv(), propogate=False, safe=False, relative=False)
                raise Exception()
        except:
            initialGuess[0:3] = [0]*3
            initialGuess[7:10] = [0]*3

        transform2 = waypoint1.DistalDubinsFrame() * waypoint2.ProximalDubinsFrame().inv()
        initialGuess[3:6] = transform2.t
        initialGuess[10:13] = SE3.Rt(transform2.R, np.zeros(3)).eul()

        try:
            initialTree.transformJoint(waypoint2Index, SE3.Trans(initialGuess[3:6]) @ SE3.Rz(initialGuess[10]) @ SE3.Ry(initialGuess[11]) @ SE3.Rz(initialGuess[12]), propogate=False, safe=False, relative=False)
            if (initialTree.detectCollisions() > 0):
                initialTree.transformJoint(waypoint2Index, transform2.inv(), propogate=False, safe=False, relative=False)
                raise Exception()
        except:
            initialGuess[3:6] = [0]*3
            initialGuess[10:13] = [0]*3

        transform3 = joint.ProximalDubinsFrame().inv() * waypoint2.DistalDubinsFrame()

        initialGuess[6] = transform3.t[2]
        initialGuess[13] = np.arctan2(transform3.R[1, 0], transform3.R[0, 0])
        try:
            jointTransform = SE3.Trans([0,0,initialGuess[6]]) @ SE3.Rz(initialGuess[13])
            initialTree.transformJoint(jointIndex, jointTransform, propogate=False, safe=False, relative=True)
            if (initialTree.detectCollisions() > 0):
                initialTree.transformJoint(jointIndex, jointTransform.inv(), propogate=False, safe=False, relative=False)
                raise Exception()
        except:
            initialGuess[6] = 0
            initialGuess[13] = 0
        
        initialLoss = objective(initialGuess)

        result = minimize(objective, initialGuess, method="L-BFGS-B", bounds=None,
            options={
            'maxiter': maxiter, 
            'ftol': 1e-2,
            #'disp': True,
        })
        
        #for some reason minimize sometimes returns value greater than initial loss
        if (result.fun > initialLoss):
            print(f"Optimized new branch {jointIndex} with waypoints:\nInitial Loss: {initialLoss}, Improved Loss: {initialLoss}, TIME: {time.time() - start}")
            return initialTree, initialLoss

        tree = copy.deepcopy(self)

        tree.transformJoint(waypoint1Index, SE3.Trans(result.x[0:3]) @ SE3.Rx(result.x[7]) @ SE3.Ry(result.x[8]) @ SE3.Rz(result.x[9]), propogate=False, safe=False, relative=True)
        tree.transformJoint(waypoint2Index, SE3.Trans(result.x[3:6]) @ SE3.Rx(result.x[10]) @ SE3.Ry(result.x[11]) @ SE3.Rz(result.x[12]), propogate=False, relative=True)
        tree.transformJoint(jointIndex, SE3.Trans([0,0,result.x[6]]) @ SE3.Rz(result.x[13]), propogate=False, relative=True)

        print(f"Optimized new branch {jointIndex} with waypoints:\nInitial Loss: {initialLoss}, Improved Loss: {result.fun}, TIME: {time.time() - start}")

        return tree, result.fun

    def optimizeJointPlacement(self, index, maxiter=1000):
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
        try:
            tree = copy.deepcopy(self)
            tree.transformJoint(index, SE3.Trans([0,0,initialPosition]) @ SE3.Rz(initialRotation), propogate=False, safe=False, relative=True)
        except:
            initialGuess = [0,0]

        initialLoss = self.calculateJointFitness(initialGuess, index)

        # dist = np.linalg.norm(self.Joints[index].ProximalDubinsFrame().t - self.Joints[self.Parents[index]].DistalDubinsFrame().t) * 2
        # bounds = [(-dist*2, dist*2), (-np.pi*2, np.pi*2)]

        result = minimize(partial(self.calculateJointFitness, index=index), initialGuess, method="Nelder-Mead",
            options={
            'maxiter': maxiter, 
            #'disp': True,
        })

        print(f"Optimized joint {index} in {time.time() - start}s")
        print(f"Old loss: {initialLoss}, Improved Loss: {result.fun}")

        tree = copy.deepcopy(self)
        tree.transformJoint(index, SE3.Trans([0,0,result.x[0]]) @ SE3.Rz(result.x[1]), propogate=False, safe=False, relative=True)

        return tree, result.fun

    def calculateJointFitness(self, params, index):
        #TODO: exponential growth (in detectCollisions?)
        translation = params[0]
        rotation = params[1]
        tree = copy.deepcopy(self)
        try:
            if not tree.transformJoint(index, SE3.Trans([0,0,translation]) @ SE3.Rz(rotation), propogate=False, safe=False, relative=True):
                return 100000
        except Exception as e:
            return 100000

        path = tree.Links[index].path

        pathCurviness = path.theta1 ** 2 * path.r + path.theta2 ** 2 * path.r

        jointDistance = np.linalg.norm(tree.Joints[index].ProximalDubinsFrame().t - tree.Joints[tree.Parents[index]].DistalDubinsFrame().t)

        #start = time.time()
        collisions = tree.detectCollisions(specificJointIndex=index) * 1000
        #print(time.time() - start)

        loss = path.length + pathCurviness * 4 + jointDistance + collisions
        #print(f"Evaluating fitness_function at x={params}, : {loss}")
        return loss

    def postOptimize(self):
        #TODO
        pass

    def save(self, filename):
        with open(f"save/{filename}.tree", "w") as f:
            save = str(self.maxAnglePerElbow) + "\n"
            for i in range(0, len(self.Joints)):
                joint = self.Joints[i]
                
                save += str(self.Parents[i]) + " "
                if isinstance(joint, Waypoint):
                    save += "Waypoint " + str(joint.numSides) + " " + str(joint.r) + " " + str(joint.pidx) + " "
                elif isinstance(joint, ExtendedRevoluteJoint):
                    save += "ExtendedRevoluteJoint " + str(joint.numSides) + " " + str(joint.r) + " " + str(joint.totalBendingAngle) + " " + str(joint.tubeLength) + " " + str(joint.numSinkLayers) + " " + str(joint.initialState) + " "
                else:
                    raise Exception("Not Implemented")
                save += "[" + ''.join([str(x) + "," for x in joint.Pose.A.reshape((16,)).tolist()])
                save += "\n"
            
            f.write(save)
            f.close()
                

    # global optimization algorithm
    def globalOptimize(self, iterations=10):

        start = time.time()

        bounds = [(-self.r*100, self.r*100)] * (len(self.Joints) - 1) + [(-np.pi*2, np.pi*2)] * (len(self.Joints) - 1)
        x0 = np.zeros(len(self.Joints)*2 - 2)
        print(f"INITIAL LOSS: {self.calculateTreeLoss(x0)}")
        result = dual_annealing(self.calculateTreeLoss, bounds=bounds, maxiter=iterations, maxfun=100, x0=x0)

        print(result.x)
        print(result.fun)
        print(f"TOOK TIME: {time.time() - start}, LOSS: {result.fun}")

        return self.transformedTree(result.x)

    def transformedTree(self, params):
        tree = copy.deepcopy(self)
        
        jointsToTransform = list(range(1,len(tree.Joints)))
        
        ctr = 0
        while len(jointsToTransform) > 0:
            idx = jointsToTransform[0]
            try:
                if not tree.transformJoint(idx, SE3.Trans([0,0,params[idx]]) @ SE3.Rz(params[idx + len(tree.Joints) - 2]), propogate=False, relative=True, safe=False):
                    raise Exception("transform didnt work")

                jointsToTransform.remove(idx)
                ctr = 0
            except:
                return None
                # jointsToTransform.remove(idx)
                # jointsToTransform.append(idx)
                # ctr += 1
                # if ctr > len(jointsToTransform):
                #     return None
        
        return tree

    def calculateTreeLoss(self, params):
        start = time.time()

        tree = self.transformedTree(params)

        if tree == None:
            print(time.time() - start)
            return 100000

        #calculate total path length
        totalLength = 0
        for link in tree.Links:
            totalLength += link.path.length

        jointDistance = 0
        for i in range(1,len(tree.Joints)):
            jointDistance += np.linalg.norm(tree.Joints[i].ProximalDubinsFrame().t - tree.Joints[tree.Parents[i]].DistalDubinsFrame().t)
        
        #ideas:
            # add incentive for joints being super close together
            # add disincentive for individual super long paths

        loss = totalLength + tree.detectCollisions()*1000

        print(time.time() - start)
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
            case "ExtendedRevoluteJoint":
                numSides = int(first[2])
                r = float(first[3])
                totalBendingAngle = float(first[4])
                tubeLength = float(first[5])
                numSinkLayers = int(first[6])
                initialState = float(first[7])
                return ExtendedRevoluteJoint(numSides, r, totalBendingAngle, tubeLength, pose, numSinkLayers, initialState)
    try:
        with open(f"save/{filename}.tree") as f:
            lines = f.readlines()
            tree = KinematicTree[OrigamiJoint](getJoint(lines[1]), float(lines[0]))
            for i in range(2, len(lines)):
                parent = int(lines[i].split(" ")[0])
                tree.addJoint(parent, getJoint(lines[i]), relative=False, fixedPosition=True, fixedOrientation=True, safe=False)
            
            return tree
    except Exception as e:
        print(e)
        raise Exception(f"file save/{filename}.tree doesnt exist")

def origamiToPrinted(tree : KinematicTree[OrigamiJoint], screwRadius: float):
    newTree = KinematicTree[PrintedJoint](tree.Joints[0].toPrinted(screwRadius), tree.maxAnglePerElbow)
    for i in range(1, len(tree.Joints)):
        try:
            tree.setJointState(i, tree.Joints[i].initialState)
            newJoint = tree.Joints[i].toPrinted(screwRadius)
            newTree.addJoint(tree.Parents[i], newJoint, relative=False, safe=False, 
            fixedPosition=True, fixedOrientation=True)
        except Exception as e:
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
    
