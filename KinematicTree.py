# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 23:13:27 2023

@author: dfesh
"""
import Joint
from Joint import *
import dubinsPath
from dubinsPath import *
import scipy
from scipy.optimize import NonlinearConstraint, minimize
import queue

class KinematicTree:
    """
    Nodes are Joint objects
    Edges are Dubins linkages from parent distal frame to child proximal frame    
    Attributes (GLOBAL COORDINATES):
        r               tubular radius
        root            root joint
        Joints          array of Joint objects (nodes)
        Parents         array of parent indices in self.Joints
        Paths           array of CSC Dubins paths to each joint from its parent
        boundingBall    ball bounding all proximal, central, and distal origins
        Children        array of arrays of child indices of each joint
    """
    def __init__(self, root):
        self.root = root
        self.r = root.r
        self.numSides = root.numSides
        self.Joints = [root]
        self.Parents = [-1]     # root has no parent
        self.Paths = [emptyCSC(self.r, root.proximalPosition(), \
                                       root.pathDirection())]
        self.boundingBall = root.boundingBall()
        if self.boundingBall.r < self.r:
            self.boundingBall = Ball(self.root.Pose.t, self.r)
        self.Children = [[]]

    """
    Returns the new Joint's index. 
    relative - boolean: is newJoint input in parent-relative coordinates (True)
                        or global coordiantes (False)?
    fixedPosition - boolean: should the joint be located
                    (True) exactly at its given position, or
                    (False) somewhere kinematically equivalent (i.e., on the
                            same z axis) chosen by the placement algorithm?
    fixedOrientation - boolean: should the joint be oriented
                    (True) exactly at its given orientation, or
                    (False) something kinematically equivalent (i.e., with the
                            same z axis) with x axis constructed as the common 
                            normal from the parent?
    guarantee - boolean: if this is True, allow the algorithm to insert 
                intermediate waypoints to route the path from the parent to 
                guarantee it avoids local self-intersection (i.e., run 
                Algorithm 9 from the Kinegami paper instead of Algorithm 8).
                This setting takes priority over fixedPosition AND 
                fixedOrientation, i.e., if guarantee is True then it doesn't 
                matter what you set those to, the algorithm will adjust the
                position and orientation of newJoint.
    """
    def addJoint(self, parentIndex, newJoint, relative=False, 
                 fixedPosition=True, fixedOrientation=True, guarantee=False):
        assert(newJoint.r == self.r)
        assert(newJoint.numSides == self.numSides)
        parent = self.Joints[parentIndex]
        if relative:
            newJoint.transformFromFrame(parent.Pose)
        
        if guarantee: # Algorithm 9
            #return self.addJointWithWayPoints(parentIndex, newJoint)
            jointsToAdd = placeJointAndWayPoints(newJoint, parent,
                                                          self.boundingBall)
            i = parentIndex
            for joint in jointsToAdd:
                i = self.addJoint(parentIndex=i, newJoint=joint)
            return i
        
        if not fixedPosition: #Algorithm 8
            newJoint = moveJointNearNeighborBut4rFromBall(newJoint, parent,
                                                          self.boundingBall)
        
        if not fixedOrientation:
            xhat = commonNormal(parent.Pose.t, parent.Pose.R[:,2],
                                newJoint.Pose.t, newJoint.Pose.R[:,2],
                                undefined=parent.Pose.R[:,0])
            newJoint.setXhatAboutZhat(xhat)
        self.boundingBall = minBoundingBall(self.boundingBall, 
                                            newJoint.boundingBall())
        newIndex = len(self.Joints)
        self.Joints.append(newJoint)
        self.Children[parentIndex].append(newIndex)
        self.Children.append([])
        self.Parents.append(parentIndex)        
        self.Paths.append(shortestCSC(self.r, 
                    parent.distalPosition(), parent.pathDirection(),
                    newJoint.proximalPosition(), newJoint.pathDirection()))
        self.plot()
        
        return newIndex
    
    """
    IN PROCESS OF IMPLEMENTING, COMMENTING OUT FOR THIS COMMIT
    
    Input:
    childrenToPlace is an array of joints
    parentIndex is the index of the parent joint
    """
    """
    def addChildrenWithWayPoints(childrenToPlace, parentIndex):
        parent = self.Joints[parentIndex]
        
        # WAYPOINT 1 (common to all children)
        nhat1 = parent.pathDirection()
        # point on bounding sphere in direction nhat1
        s1 = ball.c + ball.r * nhat1
        tangentPlane1 = Plane(s1, nhat1)
        # Construct a waypoint where this plane intersects the 
        # parent's path axis, with orientation matching parent.
        parentPathAxis = Line(parent.Pose.t, nhat1)
        originW1 = tangentPlane1.intersectionWithLine(parentPathAxis)
        # guaranteed to be a point because line is normal to plane
        PoseW1 = SE3.Rt(parent.Pose.R, originW1)
        W1 = WayPoint(jointToPlace.numSides, jointToPlace.r, 
                                   PoseW1, parent.pathIndex())
        # don't add in W1 yet because don't want to update bounding sphere?
    """            
        
    
    
    def addToPlot(self, ax, xColor='r', yColor='b', zColor='g', 
                  proximalColor='c', centerColor='m', distalColor='y',
                  pathColor='black', showCircles=False, sphereColor='black',
                  showSpheres=True):
        jointPlotHandles = []
        for joint in self.Joints:
            jointPlotHandles.append(joint.addToPlot(ax, xColor, yColor, zColor, 
                                    proximalColor, centerColor, distalColor, 
                                    sphereColor, showSpheres))
        
        for path in self.Paths:
            path.addToPlot(ax, showCircles, False, pathColor=pathColor, 
                           cscBoundaryMarker=None)
        
        if showSpheres:
            self.boundingBall.addToPlot(ax, color=sphereColor, 
                                        alpha=0.05, frame=True)
        
        return np.array(jointPlotHandles)
        
    
    def plot(self, xColor='r', yColor='b', zColor='g', 
             proximalColor='c', centerColor='m', distalColor='y',
             pathColor='black', showCircles=False):
        ax = plt.figure().add_subplot(projection='3d')
        jointPlotHandles = self.addToPlot(ax, xColor, yColor, zColor, 
                                    proximalColor, centerColor, distalColor,
                                    pathColor, showCircles)
        xHats = jointPlotHandles[:,0]
        yHats = jointPlotHandles[:,1]
        zHats = jointPlotHandles[:,2]
        origins = jointPlotHandles[:,3]
        ax.legend([tuple(xHats), tuple(yHats), tuple(zHats)], 
                  [r'$\^x$', r'$\^y$', r'$\^z$'])
        ax.set_aspect('equal')
        plt.show()
    
    
    def isLeaf(self, jointIndex):
        return len(Children[jointIndex]) == 0
    
    """ Returns list of indices of descendants of a given joint index, 
        in breadth-first order"""
    def DescendantsBreadthFirst(self, jointIndex):
        descendants = Children[jointIndex].copy()
        i = 0
        while i < len(descendants):
            d = descendants[i]
            descendants += self.DescendantsBreadthFirst(d)
            i += 0
        return descendants
    
    """ Apply given transformation (SE3() object) to given joint (index), 
        and to its descendants if recursive (defaults to False) """
    def transformJoint(self, jointIndex, Transformation, recursive=False):
        self.Joints[jointIndex].transformBy(Transformation)
        if recursive:
            for c in self.Children[jointIndex]:
                self.transformJoint(c, Transformation, recursive=True)


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
    W1 = WayPoint(jointToPlace.numSides, jointToPlace.r, PoseW1, 
                  neighbor.pathIndex())
    toReturn.append(W1)
    
    """ Translate the tangent plane forward by 4r + the new joint's 
    bounding radius, to check if we need a second waypoint. """
    farPoint1 = s1 + nhat1 * (4*jointToPlace.r + jointToPlace.boundingRadius())
    farPlane1 = Plane(farPoint1, nhat1)
    zhatNew = jointToPlace.Pose.R[:,2]
    zAxisNew = Line(jointToPlace.Pose.t, zhatNew)
    if farPlane1.sidesOfLine(zAxisNew) == [-1]:
        """ The new joint's Z axis is entirely on the near side of 
        farPlane, so we need to route through a second waypoint. 
        We construct this waypoint similarly to the above, but
        constructing points and planes in the direction of the
        new joint's Z rather than the neighbor's path direction,
        and the waypoint path direction is also along zhatNew. """
        if backwards:
            nhat2 = -zhatNew
        else:
            nhat2 = zhatNew
                
        s2 = ball.c + ball.r*nhat2
        tangentPlane2 = Plane(s2, nhat2)
        
        # TODO: think through what this is doing and if it's correct given 
        # that I'm not using the separate [a b c] frames
        R_NeighborToW2 = SO3.AngleAxis(np.pi/2, cross(nhat1, nhat2))
        RotationW2 = R_NeighborToW2.R @ neighbor.Pose.R
        originW2 = tangentPlane2.intersectionWithLine(
                                        Line(originW1 + r*nhat1, nhat2))
        PoseW2 = SE3.Rt(RotationW2, originW2)
        W2 = WayPoint(jointToPlace.numSides, r, PoseW2, neighbor.pathIndex())
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
    toReturn.append(jointToPlace)
    return toReturn
    
