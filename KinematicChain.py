# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 14:07:27 2023

@author: Daniel Feshbach
"""
from Joint import *
from KinematicTree import KinematicTree
from TubularPattern import *
from LinkCSC import LinkCSC


"""
A KinematicTree with no branching.
"""
class KinematicChain(KinematicTree):
    def __init__(self, root : Joint, maxAnglePerElbow : float = np.pi/2):
        super().__init__(root, maxAnglePerElbow)
    
    """ Add the given joint to the end of the chain, return its index """
    def addJointToEnd(self, newJoint : Joint, relative : bool = False, 
                 fixedPosition : bool = False, fixedOrientation : bool = False, 
                 guarantee : bool = False) -> int:
        parentIndex = len(self.Joints) - 1
        return super().addJoint(parentIndex, newJoint, relative, fixedPosition,
                                fixedOrientation, guarantee)
    
    def TubularPatternPattern(self, numSides : int,
                              twistPortion : float = 0.2) -> TubularPattern:
        pattern = self.Joints[0].pattern
        for j in range(1, len(self.Joints)):
            prevJoint = self.Joints[j-1]
            thisJoint = self.Joints[j]
            link = LinkCSC(self.r, prevJoint.distalDubinsFrame(), 
                                   thisJoint.proximalDubinsFrame(),
                                   self.maxAnglePerElbow)
            linkPattern = link.creasePattern(numSides, twistPortion)
            pattern.append(linkPattern)
            pattern.append(thisJoint.pattern)
            
        return pattern
