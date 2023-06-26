# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 21:54:54 2023

@author: dfesh
"""
from spatialmath import SE3
from abc import ABC, abstractmethod
from geometryHelpers import *
import matplotlib.pyplot as plt

class Joint(ABC):
    """
    r is the tubular radius
    Pose is an SE3 object where the Z axis is the joint axis
    length is from proximal to distal position (in relaxed configuration)
    """
    def __init__(self, r, Pose, length):
        self.r = r
        self.Pose = Pose
        self.length = length
    
    @abstractmethod #0 for xhat, 2 for zhat
    def pathIndex(self):
        pass
    
    def pathDirection(self):
        return self.Pose.R[:,self.pathIndex()]
    
    def distalPose(self):
        return SE3.Trans((self.length/2) * self.pathDirection()) @ self.Pose
    
    def proximalPose(self):
        return SE3.Trans(-(self.length/2) * self.pathDirection()) @ self.Pose
    
    def proximalPosition(self):
        return self.proximalPose().t
    
    def distalPosition(self):
        return self.distalPose().t    
    
    def addToPlot(self, ax, xColor='r', yColor='b', zColor='g', 
             proximalColor='c', centerColor='m', distalColor='y'):
        Poses = np.array([self.proximalPose(), self.distalPose(), self.Pose])
        oColors = np.array([proximalColor, distalColor, centerColor])
        plotHandles = addPosesToPlot(Poses, ax, self.r, 
                                     xColor, yColor, zColor, oColors)
        scale = 4
        JointAxis = np.array([self.Pose.t - scale*self.r*self.pathDirection(),
                              self.Pose.t + scale*self.r*self.pathDirection()])
        ax.plot(JointAxis[:,0], JointAxis[:,1], JointAxis[:,2], 
                linestyle='--', color='silver')
        return plotHandles
        
    
    def plot(self, xColor='r', yColor='b', zColor='g', 
             proximalColor='c', centerColor='m', distalColor='y'):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandles = self.addToPlot(ax, xColor, yColor, zColor)
        xHats, yHats, zHats, origins = plotHandles
        ax.set_aspect('equal')
        ax.legend([xHats, yHats, zHats], [r'$\^x$', r'$\^y$', r'$\^z$'])
    
    
class RevoluteJoint(Joint):
    def __init__(self, r, Pose, length):
        super().__init__(r, Pose, length)
    
    def pathIndex(self):
        return 0 # xhat
    
    
class PrismaticJoint(Joint):
    def __init__(self, r, Pose, length):
        super().__init__(r, Pose, length)
    
    def pathIndex(self):
        return 2 # zhat
    
    
class WayPoint(Joint):
    def __init__(self, r, Pose):
        super().__init__(r, Pose, 0)
    
    def pathIndex(self):
        return 2 # zhat