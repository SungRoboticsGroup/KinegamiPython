# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 21:54:54 2023

@author: dfesh
"""
from spatialmath import SE3
from abc import ABC, abstractmethod

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
    
    @abstractmethod #unit vector
    def pathDirection(self):
        pass
    
    def proximalPosition(self):
        return self.Pose.t - (self.length/2) * self.pathDirection()
    
    def distalPosition(self):
        return self.Pose.t + (self.length/2) * self.pathDirection()
    
    
class RevoluteJoint(Joint):
    def __init__(self, r, Pose, length):
        super().__init__(r, Pose, length)
    
    def pathDirection(self):
        return self.Pose.R[:,0] # xhat
    
    
class PrismaticJoint(Joint):
    def __init__(self, r, Pose, length):
        super().__init__(r, Pose, length)
    
    def pathDirection(self):
        return self.Pose.R[:,2] # xhat
    
    
class WayPoint(Joint):
    def __init__(self, r, Pose):
        super().__init__(r, Pose, 0)
    
    def pathDirection(self):
        return self.Pose.R[:,2] # xhat