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

class KinematicTree:
    """
    Nodes are Joint objects
    Edges are Dubins linkages from parent distal frame to child proximal frame    
    Attributes (GLOBAL COORDINATES):
        r           tubular radius
        root        root joint
        Joints      array of Joint objects (nodes)
        Parents     array of parent indices in self.Joints
        Paths       array of CSC Dubins paths to each joint from its parent
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
    
    
    # Algorithm 9
    def addJointWithWayPoints(self, parentIndex, newJoint, relative=False):
        assert(newJoint.r == self.r)
        assert(newJoint.numSides == self.numSides)
        parent = self.Joints[parentIndex]
        if relative:
            newJoint.transformIntoFrame(parent.Pose)
        
        """
        This algorithm (Algorithm 9 from the 2022 Kinegami paper) 
        avoids local self-intersection (Dubins paths with turning 
        angles >pi) using Conjecture 7 from that paper.               
        The central idea is to insert 1 or 2 intermediate waypoints 
        that route the path from parallel planes >=4r apart, then place 
        the new joint along its axis based on proximity to the last 
        waypoint.
        """
        
        """ First, construct the plane tangent to the bounding sphere
        and normal to the parent's path direction. """
        ppd = parent.pathDirection()
        # point on bounding sphere in direction ppd
        s = self.boundingBall.c + self.boundingBall.r * ppd 
        tangentPlane = Plane(s, ppd)
        
        """ Construct a waypoint where this plane intersects the 
        parent's path axis, with orientation matching parent. """
        parentPathAxis = Line(parent.Pose.t, ppd)
        originW1 = tangentPlane.intersectionWithLine(parentPathAxis)
        # guaranteed to be a point because line is normal to plane
        PoseW1 = SE3.Rt(parent.Pose.R, originW1)
        ppidx = parent.pathIndex() #path direction matches parent's
        W1 = WayPoint(self.numSides, self.r, PoseW1, ppidx)
        w1idx = self.addJoint(parentIndex, W1)
                                     
        """ Translate the tangent plane forward by 4r + the new joint's 
        bounding radius, to check if we need a second waypoint. """
        njbr = newJoint.boundingRadius()
        farPoint = s + ppd * (4*self.r + njbr)
        farPlane = Plane(farPoint, ppd)
        zhatNew = newJoint.Pose.R[:,2]
        zAxisNew = Line(newJoint.Pose.t, zhatNew)
        if farPlane.sidesOfLine(zAxisNew) == [-1]:
            """ The new joint's Z axis is entirely on the near side of 
            farPlane, so we need to route through a second waypoint. 
            We construct this waypoint similarly to the above, but
            constructing points and planes in the direction of the
            new joint's Z rather than the parent's path direction,
            and the waypoint path direction is also along zhatNew. """
                                
            s2 = self.boundingBall.c + self.boundingBall.r*zhatNew
            tangentPlane2 = Plane(s2, zhatNew)
            
            # TODO: think through what this is doing and if it's
            # correct given that I'm not using the separate [a b c] 
            # frames
            R_ParentToW2 = SO3.AngleAxis(np.pi/2, 
                                cross(parent.pathDirection(), zhatNew))
            RotationW2 = R_ParentToW2 @ parent.Pose.R
            originW2 = tangentPlane2.intersectionWithLine(
                                        Line(parent.Pose.t, zhatNew))
            PoseW2 = SE3.Rt(RotationW2, originW2)
            waypoint2 = WayPoint(self.numSides, self.r, PoseW2, ppidx)
            w2idx = self.addJoint(w1idx, waypoint2)
            
            farPoint2 = s2 + (4*self.r + njbr)*zhatNew
            farPlane2 = Plane(farPoint, zhatNew)
            
        else:
            """
            We only need 1 waypoint, but some values of waypoint 2 are used 
            to place+orient newJoint, so we need to define those as the
            corresponding values of waypoint 1
            """
            farPlane2 = farPlane 
            originW2 = originW1
            w2idx = w1idx
        
        def newPosition(zChange):
            return newJoint.Pose.t + zChange*newJoint.Pose.R[:,2]
        
        def distanceFromW2(zChange):
            return norm(originW2 - newPosition(zChange))
        
        def signedDistanceToFarPlane2(zChange):
            return farPlane2.signedDistanceToPoint(newPosition(zChange))
        
        farEnough = NonlinearConstraint(signedDistanceToFarPlane2, 
            lb= 0, ub= np.inf)
        
        result = minimize(distanceFromW2, 0, constraints=(farEnough))
        zChange = result.x[0]
        newJoint.translateAlongZ(zChange)
        return self.addJoint(w2idx, newJoint)
    


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
            newJoint.transformIntoFrame(parent.Pose)
        
        if guarantee: # Algorithm 9
            return self.addJointWithWayPoints(parentIndex, newJoint)
        
        
        ######################################################################
        if not fixedPosition: #Algorithm 8
            def newPosition(zChange):
                return newJoint.Pose.t + zChange*newJoint.Pose.R[:,2]
        
            # optimization objective
            def distanceFromParent(zChange): 
                return norm(parent.Pose.t - newPosition(zChange))
            
            # for constraint 
            def distanceBetweenBallCenters(zChange):
                return norm(self.boundingBall.c - newPosition(zChange))    
                
            farEnough = NonlinearConstraint(distanceBetweenBallCenters, 
                lb= 4*self.r + self.boundingBall.r + newJoint.boundingRadius(), 
                ub= np.inf)
                
            result = minimize(distanceFromParent, 0, constraints=(farEnough))
            zChange = result.x[0]
            newJoint.translateAlongZ(zChange)
        ######################################################################
        if not fixedOrientation:
            xhat = commonNormal(parent.Pose.t, parent.Pose.R[:,2],
                                newJoint.Pose.t, newJoint.Pose.R[:,2],
                                undefined=parent.Pose.R[:,0])
            newJoint.setXhatAboutZhat(xhat)
        ######################################################################
        self.boundingBall = minBoundingBall(self.boundingBall, 
                                            newJoint.boundingBall())
        self.Joints.append(newJoint)
        self.Parents.append(parentIndex)        
        self.Paths.append(shortestCSC(self.r, 
                    parent.distalPosition(), parent.pathDirection(),
                    newJoint.proximalPosition(), newJoint.pathDirection()))
        
        self.plot()
        
        return len(self.Joints)-1
    
    
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