# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 21:54:54 2023

@author: dfesh
"""
from spatialmath import SE3
from abc import ABC, abstractmethod
from geometryHelpers import *
import matplotlib.pyplot as plt
from tubularOrigami import *

class Joint(ABC):
    """
    r is the tubular radius
    Pose is an SE3 object where the Z axis is the joint axis
    length is from proximal to distal position (in relaxed configuration)
    """
    def __init__(self, r : float, length : float, Pose : SE3):
        self.r = r
        self.Pose = Pose
        self.length = length
    
    @abstractmethod #0 for xhat, 2 for zhat
    def pathIndex(self):
        pass
    
    def boundingRadius(self):
        return max(self.r, self.length / 2)
    
    def boundingBall(self):
        return Ball(self.Pose.t, self.boundingRadius())
    
    def pathDirection(self):
        return self.Pose.R[:,self.pathIndex()]
    
    def distalPose(self):
        return SE3.Trans((self.length/2) * self.pathDirection()) @ self.Pose
    
    def proximalPose(self):
        return SE3.Trans(-(self.length/2) * self.pathDirection()) @ self.Pose
    
    # Indices 0,1,2,3 with 0,1,2 cycled to begin with pathDirection
    def dubinsColumnOrder(self):
        return np.hstack((np.roll(np.arange(3), -self.pathIndex()),[3]))
    
    # Pose with axes cycled so that the first axis direction is pathDirection
    def dubinsFrame(self):
        return SE3(self.Pose.A[:,self.dubinsColumnOrder()])
    
    def proximalDubinsFrame(self):
        return SE3(self.proximalPose().A[:,self.dubinsColumnOrder()])
    
    def distalDubinsFrame(self):
        return SE3(self.distalPose().A[:,self.dubinsColumnOrder()])
    
    def proximalPosition(self):
        return self.proximalPose().t
    
    def distalPosition(self):
        return self.distalPose().t    
    
    # Frame should be an SE3 object
    def transformFromFrame(self, Frame):
        self.Pose = self.Pose @ Frame
    
    def translateAlongZ(self, zChange):
        self.Pose = self.Pose @ SE3.Trans([0,0,zChange])
    
    def setXhatAboutZhat(self, xhatNew):
        xhatNew = xhatNew / norm(xhatNew)
        zhat = self.Pose.R[:,2]
        assert(dot(zhat, xhatNew)==0) #input must be orthogonal to Z axis
        yhatNew = cross(zhat, xhatNew)
        Transform = np.eye(4)
        Transform[0:3,0] = xhatNew
        Transform[0:3,1] = yhatNew
        Transform[0:3,2] = zhat
        Transform[0:3,3] = self.Pose.t
        self.Pose = SE3(Transform)
    
    def addToPlot(self, ax, xColor='r', yColor='b', zColor='g', 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor='black', showSphere=False):
        Poses = np.array([self.proximalPose(), self.distalPose(), self.Pose])
        oColors = np.array([proximalColor, distalColor, centerColor])
        plotHandles = addPosesToPlot(Poses, ax, self.r, 
                                     xColor, yColor, zColor, oColors)
        scale = 5
        zhat = self.Pose.R[:,2]
        JointAxis = np.array([self.Pose.t - scale*self.r*zhat,
                              self.Pose.t + scale*self.r*zhat])
        ax.plot(JointAxis[:,0], JointAxis[:,1], JointAxis[:,2], 
                linestyle='--', color='silver')
        if showSphere:
            self.boundingBall().addToPlot(ax, color=sphereColor, alpha=0.05)
        return plotHandles
        
    
    def plot(self, xColor='r', yColor='b', zColor='g', 
             proximalColor='c', centerColor='m', distalColor='y'):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandles = self.addToPlot(ax, xColor, yColor, zColor)
        xHats, yHats, zHats, origins = plotHandles
        ax.set_aspect('equal')
        ax.legend([xHats, yHats, zHats], [r'$\^x$', r'$\^y$', r'$\^z$'])
    
class OrigamiJoint(Joint):
    def __init__(self, numSides, r, length, Pose):
        self.numSides = numSides
        self.polygonInnerAngle = np.pi * (numSides-2)/(2*numSides)
        super().__init__(r, length, Pose)
    
class RevoluteJoint(OrigamiJoint):
    """
    Origami revolute joint with rotation range [-totalBendingAngle/2, totalBendingAngle/2]
    and numSinkLayers recursive sink gadget layers 
    """
    def __init__(self, numSides : int, r : float, totalBendingAngle : float, 
                 Pose : SE3, numSinkLayers : int = 1):
        polygonInnerAngle = np.pi * (numSides-2)/(2*numSides)
        length = 2*r*np.sin(polygonInnerAngle)*np.tan(totalBendingAngle/4) #2*delta from paper
        super().__init__(numSides, r, length, Pose)
        self.pattern = RevoluteJointPattern(self.numSides, self.r, 
                                            totalBendingAngle, numSinkLayers)
    
    def pathIndex(self):
        return 0 # xhat
    
    def boundingRadius(self):
        return self.length / 2
    
class PrismaticJoint(OrigamiJoint):
    def __init__(self, numSides : int, r : float, neutralLength : float, 
                 numLayers : int, coneAngle : float, Pose : SE3):
        super().__init__(numSides, r, neutralLength, Pose)
        self.pattern = PrismaticJointPattern(numSides, r, neutralLength, 
                                             numLayers, coneAngle)
        
    def pathIndex(self):
        return 2 # zhat
    
    
class WayPoint(OrigamiJoint):
    # path direction through a waypoint defaults to zhat
    def __init__(self, numSides, r, Pose, pathIndex=2):
        super().__init__(numSides, r, 0, Pose)
        self.pidx = pathIndex
        self.pattern = TubularPattern(numSides, r)
    
    def pathIndex(self):
        return self.pidx
