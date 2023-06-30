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
        self.Joints = [root]
        self.Parents = [-1]     # root has no parent
        self.Paths = [emptyCSC(self.r, root.proximalPosition(), \
                                       root.pathDirection())]
        self.boundingBall = root.boundingBall()
    
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
    """
    def addJoint(self, parentIndex, newJoint, relative=False, 
                 fixedPosition=True, fixedOrientation=True):
        assert(newJoint.r == self.r)
        parent = self.Joints[parentIndex]
        if relative:
            newJoint.transformIntoFrame(parent.Pose)
            
        ######################################################################
        if not fixedPosition:
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
            
            result = minimize(distanceFromParent, 0, 
                               constraints=(farEnough))
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