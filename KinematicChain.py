# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 14:07:27 2023

@author: Daniel Feshbach
"""
from Joint import *
from OrigamiJoint import *
from KinematicTree import KinematicTree
from TubularPattern import *
from LinkCSC import LinkCSC

def remove_duplicates(arr):
    seen = []
    unique_arr = []
    for item in arr:
        if item not in seen:
            unique_arr.append(item)
            seen.append(item)
    return unique_arr

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
        chainPattern = copy.deepcopy(self.Joints[0].pattern)
        for j in range(1, len(self.Joints)):
            chainPattern.append(self.Links[j].creasePattern(self.numSides, twistPortion))
            chainPattern.append(self.Joints[j].pattern)
        return chainPattern
    
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
            if not all(len(children) <= 1 for children in self.Children):
                # parent of the selected joint
                parentIndex = self.Parents[jointIndex]
                if parentIndex >= 0:
                    self.Children[parentIndex].remove(jointIndex)
                        

                # children of the selected joint (represented as indices)
                children_to_reassign = self.Children[jointIndex]

                # assign every child of the joint to the new parent
                # and assign the new parent to every child
                for childIndex in children_to_reassign:
                    self.Parents[childIndex] = parentIndex
                    if parentIndex != -1:
                        self.Children[parentIndex].append(childIndex)

                # creates new links between children and parent
                if parentIndex != -1:
                     parentJoint = self.Joints[parentIndex]
                     for childIndex in children_to_reassign:
                         childJoint = self.Joints[childIndex]
                         newLink = LinkCSC(self.r, parentJoint.DistalDubinsFrame(), 
                                         childJoint.ProximalDubinsFrame(), self.maxAnglePerElbow)
                         self.Links[childIndex] = newLink

                # delete the selected joint
                self.Joints.pop(jointIndex)
                self.Children.pop(jointIndex)
                self.Parents.pop(jointIndex)
                self.Links.pop(jointIndex)

                self.Parents = [p - 1 if p >= jointIndex else p for p in self.Parents]

                self.Children = [[child - 1 if child >= jointIndex else child for child in children] 
                                    for children in self.Children]
                
                # creates new links between children and parent
                if parentIndex != -1:
                    parentJoint = self.Joints[parentIndex]
                    for childIndex in children_to_reassign:
                        childIndex -= 1
                        childJoint = self.Joints[childIndex]
                        newLink = LinkCSC(self.r, parentJoint.DistalDubinsFrame(), 
                                        childJoint.ProximalDubinsFrame(), self.maxAnglePerElbow)
                        self.Links[childIndex] = newLink

                if len(self.Joints) > 0:
                    self.recomputeBoundingBall()
            else:
                parentIndex = self.Parents[jointIndex]
                prevJoint = self.Joints[parentIndex] if jointIndex>0 else nextJoint

                children = self.Children[jointIndex]
                for child in children:
                    nextJoint = self.Joints[child]

                    newLink = LinkCSC(self.r, prevJoint.DistalDubinsFrame(), 
                                            nextJoint.ProximalDubinsFrame(),
                                            self.maxAnglePerElbow) 

                    linksBefore = self.Links[:jointIndex]
                    linksAfter = self.Links[jointIndex+2:]
                    self.Links = linksBefore + [newLink] + linksAfter
                    self.Joints = self.Joints[:jointIndex] + self.Joints[jointIndex+1:]

                self.recomputeBoundingBall()
                
                # recompute children and parents

                self.Children = []
                for i in range(len(self.Joints)-1):
                    self.Children.append([i+1])
                self.Children.append([])
                
                self.Parents = []
                for i in range(len(self.Joints)):
                    self.Parents.append(i-1)
            

        return True