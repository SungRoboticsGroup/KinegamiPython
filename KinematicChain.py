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
    def __init__(self, startJoint : Joint, maxAnglePerElbow : float = np.pi/2):
        super().__init__(startJoint, maxAnglePerElbow)
    
    """ Add the given joint to the end of the chain, return its index """
    def append(self, newJoint : Joint, relative : bool = True, 
                 fixedPosition : bool = False, fixedOrientation : bool = False, 
                 safe : bool = True) -> int:
        parentIndex = len(self.Joints) - 1
        return super().addJoint(parentIndex, newJoint, relative, fixedPosition,
                                fixedOrientation, safe)
    
    def creasePattern(self, twistPortion : float = 0.2) -> TubularPattern:
        pattern = self.Joints[0].pattern
        for j in range(1, len(self.Joints)):
            pattern.append(self.Links[j].creasePattern(self.numSides, twistPortion))
            pattern.append(self.Joints[j].pattern)
        return pattern
    
    def delete(self, jointIndex : int, safe : bool = True) -> bool:
        assert(jointIndex>=0)
        if safe:
            backup = self.dataDeepCopy()
            try:
                self.delete(jointIndex, safe=False)
            except ValueError as err:
                print("WARNING: something went wrong in delete:")
                print(err)
                print("Deletion canceled.")
                self.setTo(backup)
                return False
        else:
            if jointIndex > 0 and jointIndex < len(self.Joints) - 1:
                prevJoint = self.Joints[jointIndex - 1]
                nextJoint = self.Joints[jointIndex + 1]
                newLink = LinkCSC(self.r, prevJoint.DistalDubinsFrame(),
                                nextJoint.ProximalDubinsFrame(), self.maxAnglePerElbow)
                self.Links[jointIndex - 1] = newLink
                self.Links.pop(jointIndex)
            elif jointIndex == 0 and len(self.Joints) > 1:
                self.Links.pop(0)
            elif jointIndex == len(self.Joints) - 1:
                self.Links.pop(-1)

            self.Joints.pop(jointIndex)
            if hasattr(self, 'Children'):
                for i in range(len(self.Children)):
                    self.Children[i] = [child - 1 if child > jointIndex else child for child in self.Children[i]]
                self.Children.pop(jointIndex)

            if hasattr(self, 'Parents'):
                self.Parents = [parent - 1 if parent > jointIndex else parent for parent in self.Parents]
                self.Parents.pop(jointIndex)
            
            if len(self.Joints) > 0:
                self.recomputeBoundingBall()
            """"
            nextJoint = self.Joints[jointIndex+1]
            prevJoint = self.Joints[jointIndex-1] if jointIndex>0 else nextJoint
            newLink = LinkCSC(self.r, prevJoint.DistalDubinsFrame(), 
                                    nextJoint.ProximalDubinsFrame(),
                                    self.maxAnglePerElbow)
            linksBefore = self.Links[:jointIndex]
            linksAfter = self.Links[jointIndex+2:]
            self.Links = linksBefore + [newLink] + linksAfter
            self.Joints = self.Joints[:jointIndex] + self.Joints[jointIndex+1:]
            self.recomputeBoundingBall()
            
            self.Children = []
            for i in range(len(self.Joints)-1):
                self.Children.append([i+1])
            self.Children.append([])
            
            self.Parents = []
            for i in range(len(self.Joints)):
                self.Parents.append(i-1)
            """
        return True