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
import dill
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

class JointHalf:
    """
    One half of a real joint, as its own collision object. The proximal half runs from the joint's proximal end
    to its center, and the distal half from its center to its distal end (for a prismatic joint, which doesn't
    bend, the split is at the midpoint of its current ends). The two halves belong to different rigid bodies:
    the proximal half moves with the links coming into the joint, the distal half with the links leaving it.
    Each half is a flat-ended cylinder: flat at the split and at its outer end (for prismatic joints too,
    whose own sdf still has round ends: see the TODO in Prismatic.sdf). Provides what KinematicTree.collision
    needs from a tube.
    """
    def __init__(self, joint, which : str):
        assert which in ('proximal', 'distal')
        self.joint = joint
        self.which = which
        self.r = joint.r

    def __repr__(self):
        return f"JointHalf({self.which} half of {type(self.joint).__name__})"

    def splitPoint(self) -> np.ndarray:
        if isinstance(self.joint, Revolute):
            return np.asarray(self.joint.Pose.t, dtype=float)
        return (np.asarray(self.joint.ProximalDubinsFrame().t, dtype=float) +
                np.asarray(self.joint.DistalDubinsFrame().t, dtype=float)) / 2

    def endpoints(self) -> Tuple[np.ndarray, np.ndarray]:
        # start and end, in path order (root to leaves)
        if self.which == 'proximal':
            return np.asarray(self.joint.ProximalDubinsFrame().t, dtype=float), self.splitPoint()
        return self.splitPoint(), np.asarray(self.joint.DistalDubinsFrame().t, dtype=float)

    def direction(self) -> np.ndarray:
        return np.asarray(self.joint.pathDirection() if self.which == 'proximal'
                          else self.joint.distalPathDirection(), dtype=float)

    def length(self) -> float:
        a, b = self.endpoints()
        return float(np.linalg.norm(b - a))

    def interpolate(self, density : float, xp : ModuleType = np) -> ArrayLike:
        a, b = self.endpoints()
        n = max(2, int(np.ceil(density * np.linalg.norm(b - a))))
        return xp.asarray(a.reshape(1, 3) + np.linspace(0, 1, n).reshape(-1, 1) * (b - a).reshape(1, 3))

    def sdf(self, point : ArrayLike, xp : ModuleType = np) -> ArrayLike:
        point = xp.asarray(point).reshape(-1, 3)
        a, b = self.endpoints()
        n = self.direction()
        d = sdf_capsule(xp, point, a=xp.asarray(a).reshape(1, 3), b=xp.asarray(b).reshape(1, 3), r=self.r).flatten()
        # flat at both ends: the start (normal pointing back along the half) and the end
        d = xp.maximum(d, sdf_plane(xp, point, a, -n))
        return xp.maximum(d, sdf_plane(xp, point, b, n))

    def boundingBox(self, xp : ModuleType = np) -> ArrayLike:
        a, b = self.endpoints()
        return xp.vstack((xp.minimum(a, b) - self.r, xp.maximum(a, b) + self.r))

    def _radialVector(self) -> np.ndarray:
        frame = self.joint.ProximalFrame() if self.which == 'proximal' else self.joint.DistalFrame()
        return frame.R[:, 1]

    def startCircle(self, forward : bool = True) -> Circle3D:
        n = self.direction()
        return Circle3D(radius=self.r, center=self.endpoints()[0], normal=n if forward else -n,
                        radialVector=self._radialVector())

    def endCircle(self, forward : bool = True) -> Circle3D:
        n = self.direction()
        return Circle3D(radius=self.r, center=self.endpoints()[1], normal=n if forward else -n,
                        radialVector=self._radialVector())


class CollisionStructure:
    """
    The collision objects of a tree and the pairs of them that are checked (KinematicTree.buildCollisionMatrices).
    It depends only on the tree's structure, not its geometry, so it can be reused while joints move;
    KinematicTree.collisionCandidates rebuilds it if the structure has changed (signature).
    objects: list of (kind, node, which): kind 'link' (the link into node), 'half' (which is 'proximal' or
        'distal'), or 'joint' (a real joint checked whole: one that doesn't split into halves, e.g. a tip or an
        end effector)
    pairs: (2, P) array of indices into objects
    ancestors: for each node, the set of its proper ancestors
    """
    # position of each kind of object within its node, along the path from the root
    ORDER = {('link', None): 0, ('half', 'proximal'): 1, ('joint', None): 1, ('half', 'distal'): 2}

    def __init__(self, objects, pairs, ancestors, signature):
        self.objects = objects
        self.pairs = pairs
        self.ancestors = ancestors
        self.signature = signature

    def before(self, x : int, y : int) -> bool:
        # whether object x lies on the path from the root to object y, before it
        kx, nx, wx = self.objects[x]
        ky, ny, wy = self.objects[y]
        if nx == ny:
            return self.ORDER[(kx, wx)] < self.ORDER[(ky, wy)]
        return nx in self.ancestors[ny]

    def label(self, k : int) -> str:
        kind, node, which = self.objects[k]
        return f"Link {node}" if kind == 'link' else f"Joint {node}" + (f" ({which} half)" if which else "")

    def key(self, k : int) -> Tuple[int, str]:
        # (index, 'Joint' or 'Link'), as getCollidingPairs reports pairs
        kind, node, _ = self.objects[k]
        return (node, 'Link' if kind == 'link' else 'Joint')


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

    def isConsistent(self, tolerance : float = 0.1) -> bool:
        """Fast check that all links connect to their parent's distal
        and their own joint's proximal Dubins frames within tolerance."""
        for i in range(1, len(self.Joints)):
            parent_distal = self.Joints[self.Parents[i]].DistalDubinsFrame().t
            link_start = self.Links[i].StartDubinsPose.t
            if np.linalg.norm(parent_distal - link_start) > tolerance:
                return False
            joint_proximal = self.Joints[i].ProximalDubinsFrame().t
            link_end = self.Links[i].EndDubinsPose.t
            if np.linalg.norm(joint_proximal - link_end) > tolerance:
                return False
        return True

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
        # Likewise an end effector (e.g. a gripper) must remain a leaf
        if isinstance(parent, EndEffector):
            raise ValueError("Cannot add a child joint to an end effector. "
                             "End effectors must remain leaf nodes.")
        
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


    def insertWaypointsIntoLink(self, linkIndex: int, numWaypoints: int = 1) -> None:
        """
        Insert numWaypoints evenly-spaced Waypoint joints along self.Links[linkIndex].
        Mutates the tree in place. linkIndex must be > 0.
        """
        assert linkIndex > 0, "Use the GUI root-case handler for the root joint"
        assert numWaypoints >= 1

        fractions = [i / (numWaypoints + 1) for i in range(1, numWaypoints + 1)]
        sub_links = self.Links[linkIndex].splitAtFractions(fractions)

        parentIndex = self.Parents[linkIndex]
        first_new_idx = len(self.Joints)

        new_waypoints = []
        for i in range(numWaypoints):
            pose = _waypoint_pose_at_split(sub_links[i].path, sub_links[i + 1].path)
            new_waypoints.append(Waypoint(self.r, pose))

        # Rewire existing child: new parent is the last inserted waypoint
        self.Parents[linkIndex] = first_new_idx + numWaypoints - 1
        self.Links[linkIndex] = sub_links[-1]

        # Replace linkIndex with first new waypoint in parent's children list
        parent_children = self.Children[parentIndex]
        idx_in_parent = parent_children.index(linkIndex)
        parent_children[idx_in_parent] = first_new_idx

        # Append new waypoints and their incoming sub-links
        for i, (wp, sub_link) in enumerate(zip(new_waypoints, sub_links)):
            wp_parent = parentIndex if i == 0 else first_new_idx + i - 1
            wp_children = ([first_new_idx + i + 1] if i < numWaypoints - 1
                           else [linkIndex])
            self.Joints.append(wp)
            self.Links.append(sub_link)
            self.Parents.append(wp_parent)
            self.Children.append(wp_children)
            self.boundingBall = minBoundingBall(self.boundingBall, wp.boundingBall())

        # Recompute collision capsules for all affected links
        for i in range(numWaypoints):
            self.Links[first_new_idx + i].recomputeCollisionCapsules()
        self.Links[linkIndex].recomputeCollisionCapsules()


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
    
    def _addOneJoint(self, widget, index, selectedJoint, collidingJoints,
                     xColor, yColor, zColor, proximalColor, centerColor, distalColor,
                     sphereColor, showSpheres, jointColor, showJointSurface,
                     showJointAxis, jointAxisScale, showJointPoses):
        """Add GL items for a single joint (used by addToWidget and incremental updates)."""
        joint = self.Joints[index]
        if joint is None:
            return
        isColliding = index in collidingJoints
        if index == selectedJoint:
            joint.addToWidget(widget, xColor, yColor, zColor,
                    proximalColor, centerColor, distalColor,
                    sphereColor=selectedJointColor, showSphere=True,
                    surfaceColor=jointColor,
                    showSurface=showJointSurface, showAxis=showJointAxis,
                    axisScale=jointAxisScale, showPoses=showJointPoses, poseAxisScaleMultipler=2)
        else:
            joint.addToWidget(widget, xColor, yColor, zColor,
                    proximalColor, centerColor, distalColor,
                    sphereColor, showSphere=showSpheres,
                    surfaceColor=jointColor, showSurface=showJointSurface, showAxis=showJointAxis,
                    axisScale=jointAxisScale, showPoses=showJointPoses)
        if isColliding and (isinstance(joint, Revolute) or isinstance(joint, Prismatic)):
            joint.proximalCylinder().addToWidget(widget, color_list=collisionJointColor, is_joint=True)
            joint.distalCylinder().addToWidget(widget, color_list=collisionJointColor, is_joint=True)
            joint.centerSphere().addToWidget(widget, color=collisionJointColor)
        elif isColliding and isinstance(joint, EndEffector):
            # Highlight the geometry actually used for collision detection
            # (an under-approximating tube of radius r; see EndEffector.sdf)
            joint.collisionCylinder().addToWidget(widget, color_list=collisionJointColor, is_joint=True)

    def _addOneLink(self, widget, index, selectedLink, collidingLinks,
                    linkColor, linkOpacity, showLinkPath, pathColor,
                    showPathCircles, showLinkPoses, showLinkSurface):
        """Add GL items for a single link (used by addToWidget and incremental updates)."""
        link = self.Links[index]
        if link is None:
            return
        isColliding = index in collidingLinks
        if index == selectedLink:
            color = collisionLinkColor if isColliding else selectedLinkColor
        elif isColliding:
            color = collisionLinkColor
        else:
            color = linkColor
        link.addToWidget(widget, color=color,
                        alpha=linkOpacityDefault,
                        showPath=showLinkPath,
                        pathColor=pathColor,
                        showPathCircles=showPathCircles,
                        showFrames=showLinkPoses,
                        showBoundary=showLinkSurface,
                        linkID=index)

    def addToWidget(self, widget, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault,
                  proximalColor=proximalColorDefault, centerColor=centerColorDefault, distalColor=distalColorDefault,
                  showJointSurface=True, jointColor=jointColorDefault,
                  jointAxisScale=jointAxisScaleDefault, showJointAxis=True, showJointPoses=False,
                  linkColor=linkColorDefault, linkOpacity=linkOpacityDefault, surfaceOpacity=surfaceOpacityDefault, showLinkSurface=True,
                  showLinkPoses=False, showLinkPath=True, pathColor=pathColorDefault,
                  showPathCircles=False, sphereColor=sphereColorDefault,
                  showSpheres=False, showGlobalFrame=False, globalAxisScale=globalAxisScaleDefault, lastJoint=None,
                  selectedJoint=None, selectedLink=None,
                  collidingJoints=None, collidingLinks=None,
                  only_rebuild_index=None):
        # TODO: IMPLEMENT showGlobalFrame
        if collidingJoints is None:
            collidingJoints = set()
        if collidingLinks is None:
            collidingLinks = set()

        if only_rebuild_index is not None:
            # Incremental mode: remove old GL items for only_rebuild_index, rebuild those,
            # and apply model-matrix updates to everything else.
            joint = self.Joints[only_rebuild_index]
            link = self.Links[only_rebuild_index]
            if joint is not None:
                for item in (getattr(joint, '_gl_items_proximal', []) +
                             getattr(joint, '_gl_items_center', []) +
                             getattr(joint, '_gl_items_distal', [])):
                    try:
                        widget.plot_widget.removeItem(item)
                    except Exception:
                        pass
                joint.clearGLCache()
                self._addOneJoint(widget, only_rebuild_index, selectedJoint, collidingJoints,
                                   xColor, yColor, zColor, proximalColor, centerColor, distalColor,
                                   sphereColor, showSpheres, jointColor, showJointSurface,
                                   showJointAxis, jointAxisScale, showJointPoses)
            if link is not None:
                for item in getattr(link, '_gl_items', []):
                    try:
                        widget.plot_widget.removeItem(item)
                    except Exception:
                        pass
                link.clearGLCache()
                self._addOneLink(widget, only_rebuild_index, selectedLink, collidingLinks,
                                  linkColor, linkOpacity, showLinkPath, pathColor,
                                  showPathCircles, showLinkPoses, showLinkSurface)
            # Model-matrix update for all other joints and links
            for i, j in enumerate(self.Joints):
                if j is not None and i != only_rebuild_index:
                    j.updateCachedGLTransforms()
            for i, lnk in enumerate(self.Links):
                if lnk is not None and i != only_rebuild_index:
                    lnk.updateCachedGLTransforms()
            return

        if showSpheres:
            self.boundingBall.addToWidget(widget, color=sphereColor)

        for index, joint in enumerate(self.Joints):
            self._addOneJoint(widget, index, selectedJoint, collidingJoints,
                               xColor, yColor, zColor, proximalColor, centerColor, distalColor,
                               sphereColor, showSpheres, jointColor, showJointSurface,
                               showJointAxis, jointAxisScale, showJointPoses)

        for index, link in enumerate(self.Links):
            self._addOneLink(widget, index, selectedLink, collidingLinks,
                             linkColor, linkOpacity, showLinkPath, pathColor,
                             showPathCircles, showLinkPoses, showLinkSurface)

    def updateGLTransforms(self):
        """Fast-update all cached GL items by applying model-matrix transforms.
        Only works when all joints and links have a valid GL cache
        (i.e. addToWidget was called earlier and shapes haven't changed).
        Returns True if the fast path succeeded, False if a full redraw is needed."""
        for joint in self.Joints:
            if joint is None:
                continue
            if not joint.hasGLCache():
                return False
            joint.updateCachedGLTransforms()
        for link in self.Links:
            if link is None:
                continue
            if not link.hasGLCache():
                return False
            link.updateCachedGLTransforms()
        return True

    def clearAllGLCaches(self):
        """Clear GL cache on all joints and links (before a full rebuild)."""
        for joint in self.Joints:
            if joint is not None:
                joint.clearGLCache()
        for link in self.Links:
            if link is not None:
                link.clearGLCache()

    def resyncFromLightweight(self):
        """Rebuild link geometry from current joint poses after lightweight
        animation.  Call this when switching from animation back to editing
        or export so that elbows, cylinders, path objects, bounding balls,
        and collision capsules are all consistent again."""
        link_constructor = self._get_link_constructor()
        for i, link in enumerate(self.Links):
            if link is None or i == 0:
                continue
            parent_idx = self.Parents[i]
            parent = self.Joints[parent_idx]
            joint = self.Joints[i]
            self.Links[i] = link_constructor(
                self.r, parent.DistalDubinsFrame(),
                joint.ProximalDubinsFrame(),
                self.maxAnglePerElbow)
            # Transfer GL cache from old link to new one
            if not getattr(link, '_gl_shape_dirty', True):
                self.Links[i]._gl_items = link._gl_items
                self.Links[i]._gl_ref_pose = link._gl_ref_pose
                self.Links[i]._gl_shape_dirty = False
        self.recomputeBoundingBall()

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

    # ---------------------------------------------------------------------------------------------------
    # Collision checking
    #
    # Objects: every link (except the empty link 0), and every real joint, split into a proximal and a distal
    # half (JointHalf) if it's a revolute or prismatic joint, else checked whole (tips, end effectors).
    # Waypoints aren't objects. A rigid body is then one real joint's distal half, the links of its rigid part
    # (see findLinkClusters), and the proximal halves of the real joints where that part ends.
    #
    # Pairs (buildCollisionMatrices):
    # - two links: checked iff they're in different rigid parts. A part's tubes may overlap each other
    #   anywhere: at a branch point it's one piece that splits into several tubes.
    # - two joint objects: checked unless they're the two halves of the same joint (each joint's range is
    #   designed so its own halves don't collide).
    # - a joint object and a link: always checked, including the links of the joint's own bodies.
    #
    # Window (collisionCandidates): where one object of a pair lies on the path from the root to the other, the
    # stretch of one of them within PATH_WINDOW_RADII * r of the other, measured along that path, isn't checked.
    # Tubes have bend radius at least r, so a path can't come back within 2r of where it was in less than pi*r
    # (it only touches after exactly a U-turn), while parts next to each other along the path always meet
    # where they connect, which the closest-pair test would otherwise always pick (and dismiss as an end-cap
    # contact, hiding anything else). This one rule covers a joint against its incoming and outgoing links,
    # the links on either side of a short joint, and adjacent joints touching end to end (e.g. with only
    # parent-locked waypoints between them, where a proximal half and the next joint's distal half can
    # overlap if both bend far toward the same side), however many joints are in between. The window doesn't
    # pass between branches: a joint and the tubes of a sibling branch are checked everywhere.
    # Which side is cut: a link rather than a joint object, else the object farther from the root, at its end
    # nearer the other one.
    # ---------------------------------------------------------------------------------------------------

    # Length of the window along the path, in units of the tube radius
    PATH_WINDOW_RADII = np.pi

    @staticmethod
    def splitsIntoHalves(joint) -> bool:
        return isinstance(joint, (Revolute, Prismatic))

    def collisionStructureSignature(self):
        return (tuple(self.Parents), tuple(isWaypoint(j) for j in self.Joints),
                tuple(self.splitsIntoHalves(j) for j in self.Joints))

    def detectCollisions(self, specificJointIndex = None, plot: bool = False, debug: bool = False) -> int:
        # specificJointIndex: None for the whole tree, or an index or list of indices (see collisionCandidates)
        collisionMatrices = self.buildCollisionMatrices()
        numCollisions, _ = self.collisionsCountAndError(specificJointIndex, collisionMatrices, debug=debug, show=plot)
        return numCollisions

    def buildCollisionMatrices(self) -> CollisionStructure:
        """
        Build the collision objects and the pairs of them to check (see the comment above PATH_WINDOW_RADII).
        Kept under its old name: it used to return joint/joint, link/link and joint/link matrices.
        """
        partOf = {link: p for p, links in enumerate(self.findLinkClusters()) for link in links}
        objects = []
        for i in range(len(self.Joints)):
            if i > 0:
                objects.append(('link', i, None))
            if not isWaypoint(self.Joints[i]):
                if self.splitsIntoHalves(self.Joints[i]):
                    objects += [('half', i, 'proximal'), ('half', i, 'distal')]
                else:
                    objects.append(('joint', i, None))
        pairs = []
        for a in range(len(objects)):
            kindA, nodeA, _ = objects[a]
            for b in range(a + 1, len(objects)):
                kindB, nodeB, _ = objects[b]
                if kindA == 'link' and kindB == 'link':
                    if partOf.get(nodeA) == partOf.get(nodeB):
                        continue
                elif kindA != 'link' and kindB != 'link':
                    if nodeA == nodeB:
                        continue
                pairs.append((a, b))
        ancestors = []
        for i in range(len(self.Joints)):
            found, node = set(), self.Parents[i]
            while node != -1:
                found.add(node)
                node = self.Parents[node]
            ancestors.append(found)
        pairs = np.array(pairs, dtype=int).T if pairs else np.empty((2, 0), dtype=int)
        return CollisionStructure(objects, pairs, ancestors, self.collisionStructureSignature())

    def collisionObjectTube(self, structure : CollisionStructure, k : int):
        kind, node, which = structure.objects[k]
        if kind == 'link':
            return self.Links[node]
        if kind == 'half':
            return JointHalf(self.Joints[node], which)
        return self.Joints[node]

    def collisionObjectIntervals(self, structure : CollisionStructure, tubes) -> list:
        # [(start, end)] of each object along the path from the root (the root joint's proximal end is at 0)
        n = len(self.Joints)
        sProx, sSplit, sDist = np.zeros(n), np.zeros(n), np.zeros(n)
        stack = [0]
        while stack:
            i = stack.pop()
            parent = self.Parents[i]
            sProx[i] = 0.0 if parent == -1 else sDist[parent] + self.Links[i].length()
            joint = self.Joints[i]
            if not isWaypoint(joint) and self.splitsIntoHalves(joint):
                sSplit[i] = sProx[i] + JointHalf(joint, 'proximal').length()
                sDist[i] = sSplit[i] + JointHalf(joint, 'distal').length()
            else:
                prox, dist = np.asarray(joint.ProximalDubinsFrame().t), np.asarray(joint.DistalDubinsFrame().t)
                center = np.asarray(joint.Pose.t)
                sSplit[i] = sProx[i]
                sDist[i] = sProx[i] + np.linalg.norm(center - prox) + np.linalg.norm(dist - center)
            stack.extend(self.Children[i])
        intervals = []
        for kind, node, which in structure.objects:
            if kind == 'link':
                intervals.append((sDist[self.Parents[node]], sProx[node]))
            elif which == 'proximal':
                intervals.append((sProx[node], sSplit[node]))
            elif which == 'distal':
                intervals.append((sSplit[node], sDist[node]))
            else:
                intervals.append((sProx[node], sDist[node]))
        return intervals

    def collisionCandidates(self, movingJointIndex=None, collisionMatrices=None):
        """
        The pairs to check, with the window applied. Returns (candidates, structure): candidates is
        [(x, y, tubeX, tubeY, keepArcY)], where x and y index structure.objects, and keepArcY is None (check all of tubeY) or (lo, hi), the stretch of tubeY
        by arc length from its start that is checked (the cut object is always second). Pairs that lie
        entirely in the window are left out.
        movingJointIndex: None for every pair, or an index (or a list of indices, for joints that move
        together, e.g. a joint and the parent-locked waypoints carried with it): the pairs involving any of
        their objects or the links into them and out of them.
        collisionMatrices: from buildCollisionMatrices; rebuilt if missing, of another kind (e.g. from an old
        checkpoint), or from a tree of a different structure.
        """
        s = collisionMatrices
        if not isinstance(s, CollisionStructure) or s.signature != self.collisionStructureSignature():
            s = self.buildCollisionMatrices()
        pairs = s.pairs
        if movingJointIndex is not None:
            nodes = [int(movingJointIndex)] if isinstance(movingJointIndex, (int, np.integer)) else [int(i) for i in movingJointIndex]
            links = set(nodes) | {c for i in nodes for c in self.Children[i]}
            involved = np.array([(kind == 'link' and node in links) or (kind != 'link' and node in nodes)
                                 for kind, node, _ in s.objects], dtype=bool)
            pairs = pairs[:, involved[pairs[0]] | involved[pairs[1]]] if pairs.size > 0 else pairs
        tubes = [self.collisionObjectTube(s, k) for k in range(len(s.objects))]
        intervals = self.collisionObjectIntervals(s, tubes)
        window = self.PATH_WINDOW_RADII * self.r
        candidates = []
        for x, y in pairs.T:
            x, y = int(x), int(y)
            keep = None
            if s.before(x, y):
                up, down = x, y
            elif s.before(y, x):
                up, down = y, x
            else:
                up = None
            if up is not None:
                cutLength = window - (intervals[down][0] - intervals[up][1])
                if cutLength > 0:
                    if s.objects[up][0] == 'link' and s.objects[down][0] != 'link':
                        cut, other = up, down
                        keep = (0.0, tubes[up].length() - cutLength)  # the link's end, nearer the joint, is cut
                    else:
                        cut, other = down, up
                        keep = (cutLength, tubes[down].length())  # the start, nearer the other one, is cut
                    if keep[1] <= keep[0]:
                        continue
                    x, y = other, cut
            candidates.append((x, y, tubes[x], tubes[y], keep))
        return candidates, s

    def getCollidingPairs(self, coarseDistanceThresholdRatio : float = 0.5,
                          fineDistanceThresholdRatio : float = 0.001) -> list[Tuple[Tuple[int, str], Tuple[int, str]]]:
        """
        Return a list of all currently colliding pairs of objects.

        Each element is ((idx1, type1), (idx2, type2)) where type is 'Joint' or 'Link' (a joint half is
        reported as its joint), without repeats.
        """
        coarseDistanceThreshold = self.r * coarseDistanceThresholdRatio
        fineDistanceThreshold = self.r * fineDistanceThresholdRatio
        candidates, s = self.collisionCandidates()
        collidingPairs = []
        for x, y, tubeX, tubeY, keep in candidates:
            if self.collision(tubeX, tubeY, coarseDistanceThreshold, fineDistanceThreshold, keepArc2=keep) is not None:
                pair = (s.key(x), s.key(y))
                if pair not in collidingPairs:
                    collidingPairs.append(pair)
        return collidingPairs

    def collisionsCountAndError(self, movingJointIndex,
                                  collisionMatrices = None,
                                  show : bool = False, debug: bool = False,
                                  coarseDistanceThreshold : float = 0.5, fineDistanceThreshold: float = 0.001) -> Tuple[int, float]:
        # Returns (the number of colliding pairs of collision objects, the collision penalty). A pair of joints
        # can count more than once, through their halves.
        candidates, s = self.collisionCandidates(movingJointIndex, collisionMatrices)
        collisions = []
        for x, y, tubeX, tubeY, keep in candidates:
            collisionResult = self.collision(tubeX, tubeY, coarseDistanceThreshold, fineDistanceThreshold, keepArc2=keep)
            if collisionResult is not None:
                collisions.append((x, y, collisionResult))

        totalError = 0.0
        for x, y, collisionResult in collisions:
            minPoint1to2, minDist1to2, minPoint2to1, minDist2to1 = collisionResult
            distance = min(minDist1to2, minDist2to1)
            # Smooth error function based on logistic function, centered on the collision threshold (an overlap
            # of CONTACT_TOLERANCE, about 0): 0.5 there, rising to 1 - epsilon at an overlap of
            # fineDistanceThreshold (k = steepness of the transition). Every pair counted as a collision above
            # overlaps by more than the threshold, so each one adds at least 0.5, consistent with the count.
            epsilon = 1e-6
            k = np.log(1/epsilon - 1) / fineDistanceThreshold # transition over [threshold - fineDistanceThreshold, threshold]
            contact = -self.CONTACT_TOLERANCE * self.r
            totalError += 1 / (1 + np.exp(k * (distance - contact)))
            if show:
                self.show(block=False)
            if debug:
                print(f"{s.label(x)} vs {s.label(y)}")

        return len(collisions), totalError

    # Returns None for no collision or the pair of closest points for a collision
    # Every object is treated as a centerline (interpolate()) with radius .r.
    # TODO(collision): this means an EndEffector (e.g. a gripper) is
    # UNDER-approximated as a tube of radius r (the incoming tube radius): its
    # wider parts (carriage, jaws) are not collision-checked. See
    # EndEffector.sdf and KNOWN_ISSUES.md.

    # A pair collides when the tubes overlap by more than this, in units of the tube radius: only real overlaps
    # count, not tubes that just touch (e.g. after a U-turn of radius r, the tube before the turn and the part
    # after it are exactly 2r apart and touch along a line). It's slightly below 0 to absorb the error of the
    # geometry in an exact touch: the CSC path solver places links about 1e-5 r off (7.3e-6 r in a constructed
    # U-turn alongside a joint), so it was raised from 1e-9 on 2026-11-26. Sampling error is removed by
    # refining each sampled minimum (refinedMinimum).
    CONTACT_TOLERANCE = 1e-4

    @staticmethod
    def refinedMinimum(tube, points, dists, other, spacing, fullCount=None, refineFrom=-np.inf, maxCandidates=5, start=0):
        # The minimum over tube's whole center line of other.sdf, from its values dists at the sample points
        # (spaced at most `spacing` apart). Samples can miss the true minimum by up to spacing/2, so each sampled
        # local minimum within that of the best is refined by a bounded search between its neighboring samples,
        # along the true center line: the CSC path for a link (its samples are evenly spaced in its path
        # parameter), and straight segments between samples otherwise (a joint's center line is straight
        # between its proximal, center and distal points, all of which are samples, and a joint half's is
        # straight).
        # fullCount: the number of samples over the whole tube, if points are only some of them; start: the
        # index of points[0] among them.
        # refineFrom: refine only if the sampled minimum is in [refineFrom, refineFrom + spacing): below that
        # the true minimum is certainly below refineFrom as well, and from refineFrom + spacing/2 up it certainly
        # isn't, so refining can't change a decision made at refineFrom (it's only needed near the threshold,
        # which keeps it cheap).
        # Returns (point, value).
        n = len(points)
        best = int(np.argmin(dists))
        if n < 2 or not (refineFrom <= dists[best] < refineFrom + spacing):
            return points[best], dists[best]
        isLocalMin = np.ones(n, dtype=bool)
        isLocalMin[1:] &= dists[1:] <= dists[:-1]
        isLocalMin[:-1] &= dists[:-1] <= dists[1:]
        candidates = np.where(isLocalMin & (dists <= dists[best] + spacing / 2))[0]
        candidates = candidates[np.argsort(dists[candidates])][:maxCandidates]
        path = getattr(tube, 'path', None)
        if fullCount is None:
            fullCount = start + n
        def pointAt(u):  # u: fractional sample index within points
            if path is not None:
                return path.interpolate_vectorized(np.array([(u + start) / (fullCount - 1)]))[0]
            k = min(int(np.floor(u)), n - 2)
            return points[k] + (u - k) * (points[k + 1] - points[k])
        bestPoint, bestValue = points[best], dists[best]
        for c in candidates:
            lo, hi = max(0, c - 1), min(n - 1, c + 1)
            result = scipy.optimize.minimize_scalar(lambda u: float(other.sdf(pointAt(u).reshape(1, 3))[0]),
                                                    bounds=(lo, hi), method='bounded', options={'xatol': 1e-6})
            if result.fun < bestValue:
                bestPoint, bestValue = pointAt(result.x), result.fun
        return bestPoint, bestValue

    def collision(self, tube1, tube2, coarseDistanceThreshold=0.5,
                  fineDistanceThreshold=0.001, keepArc2=None) -> Optional[Tuple[np.ndarray, float, np.ndarray, float]]:
        # fineDistanceThreshold: the spacing of the center-line samples. Overlaps are judged against
        # CONTACT_TOLERANCE, from refined minimum distances (see refinedMinimum).
        # keepArc2: if given, (lo, hi): only the stretch of tube2 from arc length lo to hi (from its start)
        # counts (the rest is in the window, see collisionCandidates). The overlap is then measured from that
        # stretch's center line to tube1's surface, and its closest point is paired with the nearest point of
        # tube1's center line. (Measuring from tube1's center line to all of tube2 would find the place where
        # the two meet, inside the window, every time.) Otherwise both directions are measured and must agree.
        coarseDensity = 1 / coarseDistanceThreshold
        fineDensity = 1 / fineDistanceThreshold
        # Filter out if either tube is empty (no length)
        epsilon = 1e-2 * fineDistanceThreshold
        if not (tube1.length() > epsilon and tube2.length() > epsilon):
            return None
        # Filter based on distance between bounding boxes
        min1, max1 = tube1.boundingBox()
        min2, max2 = tube2.boundingBox()
        if sdf_aabb(min1, max1, min2, max2) >= 0:
            return None
        # Check SDF values at coarse density to see if close enough to consider collision
        if not (np.min(tube2.sdf(tube1.interpolate(density = coarseDensity))) < tube1.r + coarseDistanceThreshold or
                np.min(tube1.sdf(tube2.interpolate(density = coarseDensity))) < tube2.r + coarseDistanceThreshold):
            return None
        # Finer check
        points1 = tube1.interpolate(density = fineDensity)
        allPoints2 = tube2.interpolate(density = fineDensity)
        contact = -self.CONTACT_TOLERANCE * min(tube1.r, tube2.r)
        if keepArc2 is None:
            closest1, minDist1to2 = self.refinedMinimum(tube1, points1, tube2.sdf(points1), tube2,
                                                        fineDistanceThreshold, refineFrom=tube1.r + contact)
            minDist1to2 = minDist1to2 - tube1.r
            closest2, minDist2to1 = self.refinedMinimum(tube2, allPoints2, tube1.sdf(allPoints2), tube1,
                                                        fineDistanceThreshold, refineFrom=tube2.r + contact)
            minDist2to1 = minDist2to1 - tube2.r
            if not (minDist1to2 < contact and minDist2to1 < contact):
                return None
        else:
            # samples are evenly spaced by arc length
            n2 = len(allPoints2)
            scale = (n2 - 1) / tube2.length()
            i0 = max(0, int(np.ceil(keepArc2[0] * scale - 1e-9)))
            i1 = min(n2, int(np.floor(keepArc2[1] * scale + 1e-9)) + 1)
            if i1 <= i0:
                return None
            points2 = allPoints2[i0:i1]
            closest2, minDist2to1 = self.refinedMinimum(tube2, points2, tube1.sdf(points2), tube1,
                                                        fineDistanceThreshold, fullCount=n2,
                                                        refineFrom=tube2.r + contact, start=i0)
            minDist2to1 = minDist2to1 - tube2.r
            if not minDist2to1 < contact:
                return None
            closest1 = points1[int(np.argmin(norm(points1 - closest2, axis=1)))]
            minDist1to2 = minDist2to1
        # Filter out collisions that are solely in the hemispherical end caps.
        # The SDF is point-to-surface, so a query point inside a tube's
        # hemispherical cap appears as within r of the tube surface.
        # This happens where adjacent tubes connect. We filter these by
        # checking end-circle disc planes.
        # Use distance-based endpoint check (within r of start/end)
        # rather than exact index, since the closest point may be
        # slightly inward from the exact endpoint.
        # (Ported from chainsGUIdevelopment b38c016a and 392f2061.)
        startDist1 = norm(closest1 - points1[0])
        endDist1 = norm(closest1 - points1[-1])
        nearStart1 = startDist1 < tube1.r
        nearEnd1 = endDist1 < tube1.r
        isEndpoint1 = nearStart1 or nearEnd1

        startDist2 = norm(closest2 - allPoints2[0])  # tube2's real ends, not a cut
        endDist2 = norm(closest2 - allPoints2[-1])
        nearStart2 = startDist2 < tube2.r
        nearEnd2 = endDist2 < tube2.r
        isEndpoint2 = nearStart2 or nearEnd2

        if isEndpoint1 and isEndpoint2:
            # Both closest points are at/near endpoints - need to check disc intersection
            # Get the circles at the relevant endpoints.
            # When both nearStart and nearEnd are true (tube shorter than r),
            # pick the closer endpoint to avoid choosing the wrong disc plane.
            if nearStart1 and nearEnd1:
                circle1 = tube1.startCircle(forward=False) if startDist1 <= endDist1 else tube1.endCircle(forward=True)
            else:
                circle1 = tube1.startCircle(forward=False) if nearStart1 else tube1.endCircle(forward=True)
            if nearStart2 and nearEnd2:
                circle2 = tube2.startCircle(forward=False) if startDist2 <= endDist2 else tube2.endCircle(forward=True)
            else:
                circle2 = tube2.startCircle(forward=False) if nearStart2 else tube2.endCircle(forward=True)

            # The overlap above is between balls around the center lines, which bulge past flat ends. It's real
            # only if the solids themselves overlap: if the end discs cross, or a cross-section of either tube
            # near its closest point reaches inside the other. (This replaced a test of whether either closest
            # point is behind the other's end plane, which is true anywhere in a half-space, e.g. beside a
            # joint folded far enough that its end comes back near the part before it.)
            isTrueCollision = discs_cross(circle1, circle2) \
                        or self.crossSectionsReach(tube1, points1, closest1, tube2, contact) \
                        or self.crossSectionsReach(tube2, allPoints2, closest2, tube1, contact)
            if not isTrueCollision:
                return None
        return closest1, minDist1to2, closest2, minDist2to1

    @staticmethod
    def crossSectionsReach(tube, points, closestPoint, other, contact, spacingRadii=0.125, rings=4, spokes=24) -> bool:
        # Whether any cross-section disc of tube (radius tube.r, normal to its center line), at the center-line
        # samples points within tube.r of closestPoint, reaches inside other (other.sdf below contact). Every such
        # disc is part of the tube's solid, so this tests the solids near an end directly, flat ends included.
        # The discs are spaced about spacingRadii * r apart, each sampled on `rings` circles of `spokes` points.
        near = np.where(np.linalg.norm(points - closestPoint, axis=1) <= tube.r)[0]
        if len(near) == 0:
            return False
        spacing = np.linalg.norm(points[-1] - points[0]) / max(1, len(points) - 1) if len(points) > 1 else tube.r
        step = max(1, int(round(spacingRadii * tube.r / max(spacing, 1e-12))))
        indices = np.unique(np.concatenate([near[::step], near[-1:]]))
        angles = np.linspace(0, 2 * np.pi, spokes, endpoint=False)
        radii = tube.r * (1 - 1e-9) * np.arange(1, rings + 1) / rings
        samples = []
        for i in indices:
            tangent = points[min(i + 1, len(points) - 1)] - points[max(i - 1, 0)]
            if np.linalg.norm(tangent) < 1e-12:
                continue
            tangent = tangent / np.linalg.norm(tangent)
            u = np.cross(tangent, [1.0, 0, 0] if abs(tangent[0]) < 0.9 else [0, 1.0, 0])
            u /= np.linalg.norm(u)
            v = np.cross(tangent, u)
            ring = np.cos(angles)[:, None] * u + np.sin(angles)[:, None] * v
            samples.append(points[i] + (radii[:, None, None] * ring[None, :, :]).reshape(-1, 3))
            samples.append(points[i].reshape(1, 3))
        if not samples:
            return False
        return bool(np.min(other.sdf(np.vstack(samples))) < contact)

    def findLinkClusters(self):
        # Groups the links into rigid parts, each list being one part. Link i (the link into node i) belongs to
        # the part that starts at its nearest real-joint ancestor: following Parents from node i's parent through
        # any waypoints. So a part is everything from a real joint's distal end up to the next real joints. At a
        # branch point it's one physical piece that starts in the joint's distal connection and splits into
        # several rigidly connected tubes (direct ones and ones routed through waypoints), which may overlap
        # each other anywhere. Links in the same part aren't collision-checked against each other; everything
        # else still is, including the next real joints and their descendants against all of the part's tubes.
        # Parent-locked waypoints are still waypoints, so locking doesn't change the parts.
        parts = {}
        for i in range(1, len(self.Links)):
            source = self.Parents[i]
            while source != -1 and isWaypoint(self.Joints[source]):
                source = self.Parents[source]
            parts.setdefault(source, []).append(i)
        return list(parts.values())


    def buildCollisionPairDictionary(self):
        """
        Build a dictionary mapping each joint index to collision pairs that need checking.
        
        Rules:
        - Waypoints don't collision check directly
        - Link/Link: check iff in different rigid parts (see findLinkClusters), i.e. separated by a real joint
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
                       safe : bool = True, relative : bool = False, localOrient : bool = True,
                       printErrors=False, lightweight : bool = False) -> bool:
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
            except Exception as err:
                if printErrors:
                    print("WARNING: Unexpected error in transformJoint:")
                    print(err)
                    print("Reverting chain to before outer call.")
                self.setTo(backup)
                return False

        else:
            joint = self.Joints[jointIndex]
            if lightweight:
                # Fast path: only update link poses for GL transform, skip
                # full link reconstruction (elbows, cylinders, path objects)
                self.Joints[jointIndex].transformPoseBy(Transformation)
                self.Links[jointIndex].transformPosesInPlace(Transformation)
            elif recomputeLinkPath and jointIndex > 0:
                parent = self.Joints[self.Parents[jointIndex]]
                link_constructor = self._get_link_constructor()
                # Build the new link BEFORE moving the joint so that if link_constructor
                # raises (e.g. CSC solver fails), the joint hasn't moved yet and the tree
                # remains consistent. After transformPoseBy(T), ProximalDubinsFrame() ==
                # T @ current_ProximalDubinsFrame(), so we pass that directly.
                new_link = link_constructor(self.r, parent.DistalDubinsFrame(),
                                            Transformation @ joint.ProximalDubinsFrame(),
                                            self.maxAnglePerElbow)
                self.Joints[jointIndex].transformPoseBy(Transformation)
                self.Links[jointIndex] = new_link
            else:
                self.Joints[jointIndex].transformPoseBy(Transformation)
                self.Links[jointIndex] = self.Links[jointIndex].newLinkTransformedBy(Transformation)
            joint = self.Joints[jointIndex]  # reference is valid; transformPoseBy modifies in-place
            if not lightweight and jointIndex == 0 and not self.Links[0].length() == 0:
                raise ValueError("Error in transformJoint: Link 0 is supposed to stay empty (length 0).")
            
            if propogate:
                for c in self.Children[jointIndex]:
                    self.transformJoint(c, Transformation, propogate=True, 
                                        recomputeBoundingBall=False,
                                        recomputeLinkPath=False,
                                        safe=False, relative=False,
                                        lightweight=lightweight)
                
                if not lightweight:
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
                
    def setJointState(self, jointIndex : int, newState : float, lightweight : bool = False) -> bool:
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
                                recomputeBoundingBall=False, safe=False,
                                lightweight=lightweight)
        if not lightweight:
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
            except Exception as err:
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
            except Exception as err:
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
        # Don't add extension if filename already has one
        if not (filename.endswith('.tree') or filename.endswith('.chain')):
            name = filename + ".tree"
        else:
            name = filename
            
        if saveDir and not os.path.isabs(filename):
            name = os.path.join("save", name)
        
        os.makedirs(os.path.dirname(os.path.abspath(name)), exist_ok=True)
        with open(name, "wb") as f:
            dill.dump(self, f)

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

def isPickledTreeFile(filename : str) -> bool:
    # KinematicTree.save writes a dill pickle (since the GUI-branch merge, 2026-11-25), whose first byte is the
    # pickle protocol marker 0x80; older saves (and parameterSweep's start.tree) are text, a tree's repr
    with open(filename, "rb") as f:
        return f.read(1) == b"\x80"

def loadKinematicTree(filename : str): 
    # Loads either format: a dill pickle (KinematicTree.save) or a text repr
    if isPickledTreeFile(filename):
        return loadTree(filename)
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
                'ParallelJawGripper': ParallelJawGripper,
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


def _waypoint_pose_at_split(sub_path_before, sub_path_after) -> SE3:
    """Build an SE3 pose for a Waypoint(pathIndex=2) at the junction of two sub-paths."""
    pos  = sub_path_before.endPosition
    zhat = sub_path_before.endDir / np.linalg.norm(sub_path_before.endDir)
    xhat = sub_path_after.circleNormal1   # guaranteed ⊥ zhat
    yhat = np.cross(zhat, xhat)
    return SE3.Rt(SO3(np.column_stack([xhat, yhat, zhat])), pos)

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
    

def loadTree(filename: str):
    """Load any KinematicTree (or subclass) from a dill-serialized file.
    Also supports legacy text-based files saved with __repr__/eval."""
    if not isPickledTreeFile(filename):
        return loadKinematicTree(filename)
    try:
        with open(filename, "rb") as f:
            return dill.load(f)
    except Exception as e:
        raise Exception(f"Could not load file {filename}: {e}")
