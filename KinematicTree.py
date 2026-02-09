# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 23:13:27 2023

@author: Daniel Feshbach
"""
from ast import Raise
import Joint
from Joint import *
import PathCSC
from PathCSC import *
import scipy
from scipy.optimize import NonlinearConstraint, minimize
import queue
from LinkCSC import LinkCSC
import os
import time
from typing import Generic, TypeVar, get_args, Union, Optional, Tuple
from functools import partial
from geometryHelpers import *
import collections
import traceback
import style
from Tube import Tube

# import logging
# logging.basicConfig(level=logging.INFO)
# logger = logging.getLogger('pyswarms')
# logger.setLevel(logging.DEBUG)

F = TypeVar("F", bound="Tube")
class KinematicTree(Generic[F]):
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
    def __init__(self, root : Joint, maxAnglePerElbow : float = np.pi/12, 
                 joints : Optional[list[Joint]] = None, links : Optional[list[LinkCSC]] = None, 
                 parents : Optional[list[int]] = None, children : Optional[list[list[int]]] = None, 
                 boundingBall : Optional[Ball] = None, units : str = "Centimeter (cm)"):
            self.r = root.r
            self.units = units

            #saveParamsAreNone = np.array([joints is None, links is None, parents is None, children is None, boundingBall is None])
            if not joints is None and not links is None and not parents is None and not children is None and not boundingBall is None:
                self.Joints = joints
                self.Parents = parents
                self.Links = links
            elif joints is None and links is None and parents is None and children is None and boundingBall is None:
                self.Joints = [root]
                self.Parents = [-1]     # root has no parent
                # Use fabrication-specific link constructor if available
                link_constructor = self._get_link_constructor()
                self.Links = [link_constructor(self.r, root.ProximalDubinsFrame(),
                                        root.ProximalDubinsFrame(),
                                        maxAnglePerElbow)]
            if maxAnglePerElbow > 0 and maxAnglePerElbow < np.pi:
                self.maxAnglePerElbow = maxAnglePerElbow
            else:
                raise ValueError("ERROR: maxAnglePerElbow must be in (0, pi)")

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

    def _get_link_constructor(self):
        """Get the link constructor for this tree, defaulting to LinkCSC"""
        constructor = getattr(type(self), '_link_constructor', None)
        if constructor is None:
            return LinkCSC
        # Check if it's a method (for parametric constructors)
        if callable(constructor) and hasattr(constructor, '__self__'):
            return constructor
        return constructor
    
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
            f"boundingBall={repr(self.boundingBall)}, "
            f"units={repr(self.units)})"
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
    endPlane - Optional[Plane]: defaults to None, but if this is specified and 
            fixedPosition is False, the algorithm will place the new joint such
            that its whole bounding sphere is >= 4r from this plane.
    """
    def addJoint(self, parentIndex : int, newJoint : Joint, 
                 relative : bool = True, fixedPosition : bool = False, 
                 fixedOrientation : bool = False, 
                 safe : bool = True, endPlane : Optional[Plane] = None,
                 chooseXhatToMinPath : bool = False, 
                 relativeToDistalDubins : bool = False,
                 cachedLink : Optional[LinkCSC] = None) -> int:
        # Validate fabrication type if this tree has a type constraint
        # Skip validation for waypoints as they are fabrication-agnostic
        from Joint import Waypoint
        if not isinstance(newJoint, Waypoint):
            # Check for class-level _fabrication_type attribute first
            fabrication_type = getattr(type(self), '_fabrication_type', None)
            if fabrication_type is None:
                # Fall back to checking __orig_class__ for generic instantiation
                orig_class = getattr(self, '__orig_class__', None)
                if orig_class is not None:
                    type_args = get_args(orig_class)
                    if type_args and type_args[0] is not type(None):
                        fabrication_type = type_args[0]
            
            if fabrication_type is not None:
                if not isinstance(newJoint, fabrication_type):
                    raise TypeError(f"Joint must inherit from {fabrication_type.__name__}, got {type(newJoint).__name__}")
        
        # A forward-closing Tip (end tip) must remain a leaf — no children allowed
        parent = self.Joints[parentIndex]
        if isinstance(parent, Tip) and parent.forward:
            raise ValueError("Cannot add a child joint to an end tip (forward-closing Tip). "
                             "End tips must remain leaf nodes.")
        
        if safe and fixedPosition:
            raise ValueError("ERROR: trying to call addJoint with \
                safe and fixedPosition both True")
        if safe and fixedOrientation:
            raise ValueError("ERROR: trying to call addJoint with \
                safe and fixedOrientation both True")
        
        newJoint = copy.deepcopy(newJoint)
        parent = self.Joints[parentIndex]
        if relative and not relativeToDistalDubins:
            newJoint.transformPoseIntoFrame(parent.Pose)
        elif relativeToDistalDubins:
            newJoint.transformPoseIntoFrame(parent.DistalDubinsFrame())

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



        # Use fabrication-specific link constructor if available
        link_constructor = self._get_link_constructor()
        if not safe and not cachedLink is None:
            newLink = cachedLink
        else:
            newLink = link_constructor(self.r, parent.DistalDubinsFrame(), 
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

        return newIndex
    
    
    def deleteJoint(self, jointIndex : int, recursive : bool = False) -> None:
        """
        Delete a joint from the tree.
        
        Parameters:
            jointIndex : int - Index of the joint to delete
            recursive : bool - If True, delete all descendants. If False, re-parent 
                              children to the deleted joint's parent (default: False)
        
        When recursive=True:
        1. Finds all descendants of the joint (recursively)
        2. Removes them from Joints, Links, Parents, and Children
        3. Adjusts all remaining joint indices in Parents and Children
        
        When recursive=False:
        1. Re-parents all children to the deleted joint's parent
        2. Creates new links from parent to each child
        3. Removes only the specified joint
        4. Adjusts remaining joint indices
        """
        if jointIndex < 0 or jointIndex >= len(self.Joints):
            raise ValueError(f"Joint index {jointIndex} out of range [0, {len(self.Joints)-1}]")
        
        if jointIndex == 0:
            raise ValueError("Cannot delete the root joint (index 0)")
        
        parent_index = self.Parents[jointIndex]
        
        if recursive:
            # Original recursive deletion behavior
            # Step 1: Find all descendants (including the joint itself)
            descendants = set()
            to_process = [jointIndex]
            
            while to_process:
                current = to_process.pop()
                descendants.add(current)
                # Add all children of current joint to the processing queue
                to_process.extend(self.Children[current])
            
            # Sort descendants in descending order so we can delete from the end
            # This prevents index shifting issues during deletion
            descendants_sorted = sorted(descendants, reverse=True)
            
            # Step 2: Remove the joint from its parent's children list
            if parent_index >= 0:
                self.Children[parent_index] = [c for c in self.Children[parent_index] if c not in descendants]
            
            # Step 3: Delete all descendants from the lists (from end to start)
            for idx in descendants_sorted:
                del self.Joints[idx]
                del self.Links[idx]
                del self.Parents[idx]
                del self.Children[idx]
            
            # Step 4: Adjust all indices in Parents and Children
            # Create a mapping from old indices to new indices
            index_mapping = {}
            offset = 0
            for old_idx in range(len(self.Joints) + len(descendants_sorted)):
                if old_idx in descendants:
                    offset += 1
                else:
                    index_mapping[old_idx] = old_idx - offset
            
            # Update Parents array
            for i in range(len(self.Parents)):
                if self.Parents[i] >= 0:
                    self.Parents[i] = index_mapping[self.Parents[i]]
            
            # Update Children array
            for i in range(len(self.Children)):
                self.Children[i] = [index_mapping[c] for c in self.Children[i]]
        else:
            # Non-recursive: re-parent children to the deleted joint's parent
            children_to_reparent = self.Children[jointIndex].copy()
            
            # Step 1: Update parent's children list - remove deleted joint, add its children
            if parent_index >= 0:
                self.Children[parent_index] = [c for c in self.Children[parent_index] if c != jointIndex]
                self.Children[parent_index].extend(children_to_reparent)
            
            # Step 2: Create new links from parent to each child and update their parents
            link_constructor = self._get_link_constructor()
            parent_joint = self.Joints[parent_index]
            
            for child_idx in children_to_reparent:
                child_joint = self.Joints[child_idx]
                # Create new link from parent's distal frame to child's proximal frame
                new_link = link_constructor(self.r, parent_joint.DistalDubinsFrame(),
                                           child_joint.ProximalDubinsFrame(),
                                           self.maxAnglePerElbow)
                if new_link is None:
                    print(f"WARNING: Could not create valid link from joint {parent_index} to joint {child_idx}")
                    # Keep the old link as a fallback
                else:
                    self.Links[child_idx] = new_link
                
                # Update parent reference
                self.Parents[child_idx] = parent_index
            
            # Step 3: Delete the joint
            del self.Joints[jointIndex]
            del self.Links[jointIndex]
            del self.Parents[jointIndex]
            del self.Children[jointIndex]
            
            # Step 4: Adjust all indices in Parents and Children
            # All indices > jointIndex need to be decremented by 1
            for i in range(len(self.Parents)):
                if self.Parents[i] > jointIndex:
                    self.Parents[i] -= 1
                elif self.Parents[i] == jointIndex:
                    # This shouldn't happen since we already re-parented children
                    raise RuntimeError(f"Found orphaned joint at index {i}")
            
            for i in range(len(self.Children)):
                self.Children[i] = [c - 1 if c > jointIndex else c for c in self.Children[i]]
        
        # Recompute bounding ball since we removed joints
        self.recomputeBoundingBall()
    
    
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
            
    
    def addSubtree(self, parentIndex : int, subtree : 'KinematicTree') -> list[int]:
        """
        Copy all joints and links from another tree into this tree, 
        attaching the subtree's root as a child of parentIndex.
        
        Uses the subtree's existing links (cachedLink) to avoid recomputing paths.
        
        Parameters:
            parentIndex : int - Index of the joint in this tree to attach the subtree's root to
            subtree : KinematicTree - The tree whose joints/links to copy in
            
        Returns:
            list[int] - The new indices of all added joints (in the order they appeared in the subtree)
        """
        if parentIndex < 0 or parentIndex >= len(self.Joints):
            raise ValueError(f"Parent index {parentIndex} out of range [0, {len(self.Joints)-1}]")
        
        # A backward-closing Tip (start tip) must remain a root — it cannot become a non-root
        # node by being grafted as a subtree under another joint
        subtree_root = subtree.Joints[0]
        if isinstance(subtree_root, Tip) and not subtree_root.forward:
            raise ValueError("Cannot attach a subtree whose root is a start tip (backward-closing Tip). "
                             "Start tips must remain root nodes.")
        
        # Map from old subtree indices to new indices in this tree
        index_mapping = {}
        
        # BFS to add joints in parent-first order so parents exist before children
        queue_items = collections.deque()
        # Start with root of subtree (index 0), attaching to parentIndex in this tree
        queue_items.append((0, parentIndex))
        
        new_indices = []
        
        while queue_items:
            subtree_idx, new_parent_idx = queue_items.popleft()
            subtree_joint = subtree.Joints[subtree_idx]
            
            # Use cached link for non-root subtree joints (root's link is a 
            # degenerate self-link, so we need a fresh link from the new parent)
            cachedLink = subtree.Links[subtree_idx] if subtree_idx != 0 else None
            
            new_idx = self.addJoint(
                parentIndex=new_parent_idx,
                newJoint=subtree_joint,
                relative=False,
                fixedPosition=True,
                fixedOrientation=True,
                safe=False,
                cachedLink=cachedLink
            )
            
            if new_idx is None:
                print(f"WARNING: Failed to add subtree joint {subtree_idx}")
                continue
            
            index_mapping[subtree_idx] = new_idx
            new_indices.append(new_idx)
            
            # Enqueue all children of this subtree joint
            for child_idx in subtree.Children[subtree_idx]:
                queue_items.append((child_idx, new_idx))
        
        return new_indices


    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                  proximalColor='c', centerColor='m', distalColor='y',
                  showJointSurface=True, jointAxisScale=jointAxisScaleDefault, showJointPoses=True,
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
            handles = joint.addToPlot(ax=ax, xColor=xColor, yColor=yColor, zColor=zColor, 
                                    proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor, 
                                    sphereColor=sphereColor, showSphere=showSpheres, 
                                    surfaceOpacity=surfaceOpacity,
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
            plotPoint = plotPoint.reshape((-1,3))
            ax.scatter(plotPoint[0,:], plotPoint[1,:], plotPoint[2,:], color='red', s=50)
        return np.array(xyzHandles), np.array(abcHandles)
    
    def addToWidget(self, widget, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                  proximalColor=proximalColorDefault, centerColor=centerColorDefault, distalColor=distalColorDefault,
                  showJointSurface=True, jointColor=jointColorDefault,
                  jointAxisScale=jointAxisScaleDefault, showJointPoses=False,
                  linkColor=linkColorDefault, surfaceOpacity=surfaceOpacityDefault, showLinkSurface=True, 
                  showLinkPoses=False, showLinkPath=True, pathColor=pathColorDefault,
                  showPathCircles=False, sphereColor=sphereColorDefault,
                  showSpheres=False, showGlobalFrame=False, globalAxisScale=globalAxisScaleDefault, lastJoint=None, 
                  selectedJoint=None, selectedLink=None):
        # TODO: IMPLEMENT showGlobalFrame
        if showSpheres:
            self.boundingBall.addToWidget(widget, color=sphereColor)
            
        for index, joint in enumerate(self.Joints):
            if index == selectedJoint:
                joint.addToWidget(widget, xColor, yColor, zColor, 
                        proximalColor, centerColor, distalColor, 
                        sphereColor=selectedJointColor, showSphere=True,
                        surfaceColor=jointColor, showSurface=showJointSurface, 
                        axisScale=jointAxisScale, showPoses=showJointPoses, poseAxisScaleMultipler=2)
            else:
                joint.addToWidget(widget, xColor, yColor, zColor, 
                        proximalColor, centerColor, distalColor, 
                        sphereColor, showSphere=showSpheres, 
                        surfaceColor=jointColor, showSurface=showJointSurface, 
                        axisScale=jointAxisScale, showPoses=showJointPoses)
                
        for index, link in enumerate(self.Links):
            if index == selectedLink:
                link.addToWidget(widget, color=selectedLinkColor, 
                                alpha=linkOpacityDefault,
                                showPath=showLinkPath, 
                                pathColor=pathColor,
                                showPathCircles=showPathCircles, 
                                showFrames=showLinkPoses,
                                showBoundary=showLinkSurface,
                                linkID=index)
            else:
                link.addToWidget(widget, color=linkColorDefault, 
                                alpha=linkOpacityDefault,
                                showPath=showLinkPath, 
                                pathColor=pathColor,
                                showPathCircles=showPathCircles, 
                                showFrames=showLinkPoses,
                                showBoundary=showLinkSurface,
                                linkID=index)

        widget.add_chain(self)

    def copyAbbreviatedSelf(self, isolate=False, isolateJoint = 0):
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

    def detectCollisions(self, specificJointIndex : Optional[int] = None, plot: bool = False, debug: bool = False) -> int:
        collisionMatrices = self.buildCollisionMatrices()
        numCollisions, _ = self.collisionsCountAndError(specificJointIndex, collisionMatrices, debug=debug, show=plot)
        return numCollisions

    
    def collisionPairsFromMovingJoint(self, movingJointIndex : int, 
            collisionMatrices : Tuple[np.ndarray, np.ndarray, np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        jointJointCollisionMatrix, linkLinkCollisionMatrix, jointLinkCollisionMatrix = collisionMatrices
        i = movingJointIndex
        
        # Vectorized: Find all collidable objects for this joint
        # Joint-Joint collisions where row i is True
        joints_paired_with_joint_i = np.where(jointJointCollisionMatrix[i, :])[0]
        jj_pairs = np.vstack((np.full(len(joints_paired_with_joint_i), i), joints_paired_with_joint_i)) if len(joints_paired_with_joint_i) > 0 else np.empty((2, 0), dtype=int)
        
        # Joint-Link collisions where row i is True
        links_paired_with_joint_i = np.where(jointLinkCollisionMatrix[i, :])[0]
        jl_pairs_list = [np.vstack((np.full(len(links_paired_with_joint_i), i), links_paired_with_joint_i))] if len(links_paired_with_joint_i) > 0 else []
        
        # Find all collidable objects for this joint's incoming and outgoing links
        adjacentLinks = np.array([i] + self.Children[i])
        ll_pairs_list = []
        for linkIdx in adjacentLinks:
            links_paired_with_link_idx = np.where(linkLinkCollisionMatrix[linkIdx, :])[0]
            if len(links_paired_with_link_idx) > 0:
                ll_pairs_list.append(np.vstack((np.full(len(links_paired_with_link_idx), linkIdx), links_paired_with_link_idx)))
            joints_paired_with_link_idx = np.where(jointLinkCollisionMatrix[:, linkIdx])[0]
            if len(joints_paired_with_link_idx) > 0:
                jl_pairs_list.append(np.vstack((joints_paired_with_link_idx, np.full(len(joints_paired_with_link_idx), linkIdx))))
        ll_pairs = np.hstack(ll_pairs_list) if ll_pairs_list else np.empty((2,0), dtype=int)
        jl_pairs = np.hstack(jl_pairs_list) if jl_pairs_list else np.empty((2,0), dtype=int)
        
        # Remove duplicates (e.g., (i,j) and (j,i))
        jj_pairs = np.unique(np.sort(jj_pairs, axis=0), axis=1) if jj_pairs.size > 0 else jj_pairs
        ll_pairs = np.unique(np.sort(ll_pairs, axis=0), axis=1) if ll_pairs.size > 0 else ll_pairs

        return jj_pairs, jl_pairs, ll_pairs

    
    def getAllCollisionPairs(self, collisionMatrices : Tuple[np.ndarray, np.ndarray, np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        jointJointCollisionMatrix, linkLinkCollisionMatrix, jointLinkCollisionMatrix = collisionMatrices
        jj_pairs = np.array(np.where(np.triu(jointJointCollisionMatrix, k=1)))
        jl_pairs = np.array(np.where(jointLinkCollisionMatrix))
        ll_pairs = np.array(np.where(np.triu(linkLinkCollisionMatrix, k=1)))
        return jj_pairs, jl_pairs, ll_pairs
    
    
    def collisionsCountAndError(self, movingJointIndex: Optional[int], 
                                  collisionMatrices : Tuple[np.ndarray, np.ndarray, np.ndarray], 
                                  show : bool = False, debug: bool = False,
                                  coarseDistanceThreshold : float = 0.5, fineDistanceThreshold: float = 0.001) -> Tuple[int, float]:
        jointJointCollisionMatrix, linkLinkCollisionMatrix, jointLinkCollisionMatrix = collisionMatrices
        numCollisions = 0
        pairs = []

        if movingJointIndex is None:
            jj_pairs, jl_pairs, ll_pairs = self.getAllCollisionPairs(collisionMatrices)
        else:
            jj_pairs, jl_pairs, ll_pairs = self.collisionPairsFromMovingJoint(movingJointIndex, collisionMatrices)          
        
        collisions = []
        totalError = 0.0
        for j1, j2 in jj_pairs.T:
            collisionResult = self.collision(self.Joints[j1], self.Joints[j2], 
                                                coarseDistanceThreshold, fineDistanceThreshold)
            if collisionResult is not None:
                collisions.append(( (j1, 'Joint'), (j2, 'Joint'), collisionResult))
        for j, l in jl_pairs.T:
            collisionResult = self.collision(self.Joints[j], self.Links[l], 
                                                coarseDistanceThreshold, fineDistanceThreshold)
            if collisionResult is not None:
                collisions.append(( (j, 'Joint'), (l, 'Link'), collisionResult))
        for l1, l2 in ll_pairs.T:
            collisionResult = self.collision(self.Links[l1], self.Links[l2], 
                                                coarseDistanceThreshold, fineDistanceThreshold)
            if collisionResult is not None:
                collisions.append(( (l1, 'Link'), (l2, 'Link'), collisionResult))
        
        for collision in collisions:
            (idx1, type1), (idx2, type2), collisionResult = collision
            minPoint1to2, minDist1to2, minPoint2to1, minDist2to1 = collisionResult
            distance = min(minDist1to2, minDist2to1)
            # Smooth error function based on logistic function
            k = 50 #steepness of the transition
            totalError += 1 / (1 + np.exp(-k * distance))
            if show:
                self.show(block=False)
            if debug:
                print(f"{type1} {idx1} vs {type2} {idx2}")
            
        return len(collisions), totalError

    # def getCollisionErrorFromDict(self, indices, collisionPairDict):
    #     totalError = 0
        
    #     for idx in indices:
    #         pairs = collisionPairDict.get(idx, [])
            
    #         for (obj1, obj2) in pairs:
    #             idx1, type1 = obj1
    #             idx2, type2 = obj2
                
    #             collisionResult = self.collision(idx1, type1, idx2, type2)
    #             if collisionResult is not None:
    #                 minPoint1to2, minDist1to2, minPoint2to1, minDist2to1 = collisionResult
    #                 distance = min(minDist1to2, minDist2to1)
    #                 # Smooth error function based on logistic function
    #                 k = 50 #steepness of the transition
    #                 error = 1 / (1 + np.exp(-k * distance))
    #                 totalError += error
        
    #     return totalError
    
    
    # Returns None for no collision or the pair of closest points for a collision
    def collision(self, tube1, tube2, coarseDistanceThreshold=0.5, 
                  fineDistanceThreshold=0.001) -> Optional[Tuple[np.ndarray, float, np.ndarray, float]]:
        coarseDensity = 1 / coarseDistanceThreshold
        fineDensity = 1 / fineDistanceThreshold        
        # Filter out if either tube is empty (no length)
        epsilon = 1e-2 * fineDistanceThreshold
        if tube1.length() > epsilon and tube2.length() > epsilon:
            # Filter based on distance between bounding boxes
            min1, max1 = tube1.boundingBox()
            min2, max2 = tube2.boundingBox()
            if sdf_aabb(min1, max1, min2, max2) < 0: # bounding boxes overlap
                # Check SDF values at coarse density to see if close enough to consider collision
                if np.min(tube2.sdf(tube1.interpolate(density = coarseDensity))) < tube1.r + coarseDistanceThreshold or \
                np.min(tube1.sdf(tube2.interpolate(density = coarseDensity))) < tube2.r + coarseDistanceThreshold:
                    # Finer check
                    points1 = tube1.interpolate(density = fineDensity)
                    dists1to2 = tube2.sdf(points1)
                    minIdx1to2 = np.argmin(dists1to2)
                    minDist1to2 = dists1to2[minIdx1to2]
                    points2 = tube2.interpolate(density = fineDensity)
                    dists2to1 = tube1.sdf(points2)
                    minIdx2to1 = np.argmin(dists2to1)
                    minDist2to1 = dists2to1[minIdx2to1]
                    if minDist1to2 < tube1.r + fineDistanceThreshold and minDist2to1 < tube2.r + fineDistanceThreshold:
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
                            return points1[minIdx1to2], minDist1to2, points2[minIdx2to1], minDist2to1                
        return None

    def findLinkClusters(self):
        # Find clusters of waypoints and their incoming + outgoing links, that are connected without real joints between them
        # Every link except link 0 should be in a cluster

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
    
    
    def buildCollisionMatrices(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Build matrices indicating which joint/joint, link/link, and joint/link pairs need collision checking.
        
        Returns:
            jointJointCollisionMatrix (np.ndarray): Boolean matrix for joint/joint collisions.
            linkLinkCollisionMatrix (np.ndarray): Boolean matrix for link/link collisions.
            jointLinkCollisionMatrix (np.ndarray): Boolean matrix for joint/link collisions.
        
        Rules:
        - Waypoints don't collision check directly (all False)
        - Link 0 is empty, not checked against anything (all False)
        - Link/Link: check iff separated by real joint (different waypoint cluster)
        - Real Joint/Real Joint: always check
        - Real Joint/Link: always check
        """
        link_sets = self.findLinkClusters()
        # Create link-to-set mapping 
        whichLinkSet = {}
        for set_idx, link_set in enumerate(link_sets):
            for link_idx in link_set:
                whichLinkSet[link_idx] = set_idx
        
        # Pre-classify all joints
        real_joints = [i for i in range(len(self.Joints)) if not isWaypoint(self.Joints[i])]

        # Initialize matrices
        jointJointCollisionMatrix = np.zeros((len(self.Joints), len(self.Joints)), dtype=bool)
        linkLinkCollisionMatrix = np.zeros((len(self.Links), len(self.Links)), dtype=bool)
        jointLinkCollisionMatrix = np.zeros((len(self.Joints), len(self.Links)), dtype=bool)

        for i in range(len(self.Joints)):
            # Joint i
            if not isWaypoint(self.Joints[i]):
                # Check real joint against all other real joints and against all links
                for j in range(len(self.Joints)):
                    jointJointCollisionMatrix[i, j] = not isWaypoint(self.Joints[j]) and i != j  
                    jointLinkCollisionMatrix[i, j] = True
            # Link i
            if i > 0:
                set_idx = whichLinkSet.get(i)
                for j in range(len(self.Links)):
                    linkLinkCollisionMatrix[i, j] = whichLinkSet.get(j) != set_idx
        
        return jointJointCollisionMatrix, linkLinkCollisionMatrix, jointLinkCollisionMatrix
    
    
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
        link_sets = self.findLinkClusters()
        
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
             showJointSurface=True, jointAxisScale=jointAxisScaleDefault, showJointPoses=True,
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
                                                showJointSurface=showJointSurface, jointAxisScale=jointAxisScale, showJointPoses=showJointPoses,
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
                       safe : bool = True, relative : bool = False, localOrient : bool = True, printErrors=False) -> bool:
        if relative:
            if localOrient:
                # Use the joint's full orientation for the transformation
                Transformation = self.Joints[jointIndex].Pose @ Transformation @ self.Joints[jointIndex].Pose.inv()
            else:
                # Only use the joint's position, not orientation
                translate_transform = SE3(self.Joints[jointIndex].Pose.t)
                Transformation = translate_transform @ Transformation @ translate_transform.inv()
        
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
                link_constructor = self._get_link_constructor()
                self.Links[jointIndex] = link_constructor(self.r, parent.DistalDubinsFrame(), 
                                        joint.ProximalDubinsFrame(),
                                        self.maxAnglePerElbow)
            else:
                self.Links[jointIndex] = self.Links[jointIndex].newLinkTransformedBy(Transformation)
            if jointIndex == 0 and not self.Links[0].length() == 0:
                raise ValueError("Error in transformJoint: Link 0 is supposed to stay empty (length 0).")
            
            if propogate:
                for c in self.Children[jointIndex]:
                    self.transformJoint(c, Transformation, propogate=True, 
                                        recomputeBoundingBall=False,
                                        recomputeLinkPath=False,
                                        safe=False, relative=False)
                
                self.recursivelyRecomputeCollisionCapsules(jointIndex)
            else:
                link_constructor = self._get_link_constructor()
                for c in self.Children[jointIndex]:
                    child = self.Joints[c]
                    self.Links[c] = link_constructor(self.r, joint.DistalDubinsFrame(), 
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
        # Don't add extension if filename already has one
        if not (filename.endswith('.tree') or filename.endswith('.chain')):
            name = filename + ".tree"
        else:
            name = filename
            
        if saveDir and not os.path.isabs(filename):
            name = os.path.join("save", name)
        
        with open(name, "w") as f:
            save = self.__repr__()
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
    from numpy import array
    from geometryHelpers import Ball, Cylinder, Plane, Circle3D, Arc3D
    try:
        with open(filename) as f:
            data = f.read()
            # Provide necessary imports for eval
            eval_namespace = {
                'KinematicTree': KinematicTree,
                'LinkCSC': LinkCSC,
                'PathCSC': PathCSC,
                'SE3': SE3,
                'array': array,
                'Ball': Ball,
                'Cylinder': Cylinder,
                'Plane': Plane,
                'Circle3D': Circle3D,
                'Arc3D': Arc3D,
                'Prismatic': Prismatic,
                'TransverseRevolute': TransverseRevolute,
                'CoaxialRevolute': CoaxialRevolute,
                'Waypoint': Waypoint,
                'Tip': Tip,
            }
            tree = eval(data, eval_namespace)
            f.close()
            return tree
    except Exception as e:
        print(e)
        raise Exception(f"Could not load file {filename}: {e}")

def isWaypoint(joint):
    if joint is None:
        return False
    return isinstance(joint, Waypoint)

def curvinessOfLink(link : LinkCSC):
    return link.path.theta1 ** 1.5 * link.path.r + link.path.theta2 ** 1.5 * link.path.r


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
    W1 = Waypoint(jointToPlace.r, PoseW1, neighbor.pathIndex())
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
        W2 = Waypoint(r, PoseW2, neighbor.pathIndex())
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
    
