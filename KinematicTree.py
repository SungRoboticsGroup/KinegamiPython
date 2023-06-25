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
    
    
    def addToPlot(self, ax, xColor='r', yColor='b', zColor='g', 
                  proximalColor='c', centerColor='m', distalColor='y',
                  pathColor='black', showCircles=False):
        for joint in self.Joints:
            joint.addToPlot(ax, xColor, yColor, zColor, 
                            proximalColor, centerColor, distalColor)
        for path in self.Paths:
            path.addToPlot(ax, showCircles, False, pathColor=pathColor)
        
    
    def plot(self, xColor='r', yColor='b', zColor='g', 
             proximalColor='c', centerColor='m', distalColor='y',
             pathColor='black', showCircles=False):
        ax = plt.figure().add_subplot(projection='3d')
        self.addToPlot(ax, xColor, yColor, zColor, 
                       proximalColor, centerColor, distalColor,
                       pathColor, showCircles)
        ax.set_aspect('equal')
        ax.legend()
        plt.show()