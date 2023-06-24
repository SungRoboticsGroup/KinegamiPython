# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 23:13:27 2023

@author: dfesh
"""
import Joint
from Joint import *
import dubinsPath
from dubinsPath import *

class KinematicTree:
    """
    Nodes are Joint objects in parent-relative coordinates
    Edges are Dubins linkages from parent distal frame to child proximal frame
    
    Attributes:
        r           tubular radius
        root        root joint, in global coordinates
        Joints      array of Joint objects (nodes) in parent coordinates
        Parents     array of parent indices in self.Joints
        Paths       array of CSC Dubins paths (in parent coordiantes)
                    to each joint from its parent
    """
    def __init__(self, root):
        self.root = root
        self.r = root.r
        self.Joints = [root]
        self.Parents = [-1]     # root has no parent
        self.Paths = [emptyCSC(self.r, root.proximalPosition(), \
                                       root.pathDirection())]
    
    # return its index
    def addJoint(self, parentIndex, NewJoint):
        assert(NewJoint.r == self.r)
        self.Joints.append(NewJoint)
        self.Parents.append(parentIndex)
        parent = self.Joints[parentIndex]
        self.Paths.append(shortestCSC(self.r, 
                    parent.distalPosition(), parent.pathDirection(), 
                    NewJoint.proximalPosition(), NewJoint.pathDirection()))
        return len(self.Joints)-1