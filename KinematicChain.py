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
    
    """ WARNING: this might have bugs. Or the bugs might be in the GUI.
    I don't have time to figure it out right now, so I'm implementing the
    (less efficient) deleteJoint function instead of using this method. """
    def delete(self, jointIndex : int, safe : bool = True) -> bool:
        assert(jointIndex>=0 and jointIndex<len(self.Joints) and len(self.Joints)>1)
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
        elif jointIndex == len(self.Joints)-1:
            self.Links = self.Links[:-1]
            self.Joints = self.Joints[:-1]
            self.Children = self.Children[:-1]
            self.Children[-1] = []
            self.Parents = self.Parents[:-1]
            self.recomputeBoundingBall()
        else:
            nextJoint = self.Joints[jointIndex+1]
            if jointIndex>0:
                prevJoint = self.Joints[jointIndex-1]
                newLink = LinkCSC(self.r, prevJoint.DistalDubinsFrame(), 
                                        nextJoint.ProximalDubinsFrame(),
                                        self.maxAnglePerElbow)
            else:
                newLink = LinkCSC(self.r, nextJoint.ProximalDubinsFrame(), 
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
            return True
    
    def save(self, filename: str):
        with open(f"save/{filename}.chain", "w") as f:
            save = str(self.maxAnglePerElbow) + "\n"
            for i in range(0, len(self.Joints)):
                joint = self.Joints[i]
                
                save += str(self.Parents[i]) + " "
                if isinstance(joint, Waypoint):
                    save += "Waypoint " + str(joint.numSides) + " " + str(joint.r) + " " + str(joint.pidx) + " "
                elif isinstance(joint, RevoluteJoint):
                    save += "RevoluteJoint " + str(joint.numSides) + " " + str(joint.r) + " " + str(joint.totalBendingAngle) + " " + str(joint.numSinkLayers) + " " + str(joint.initialState) + " "
                elif isinstance(joint, ExtendedRevoluteJoint):
                    save += "ExtendedRevoluteJoint " + str(joint.numSides) + " " + str(joint.r) + " " + str(joint.totalBendingAngle) + " " + str(joint.tubeLength) + " " + str(joint.numSinkLayers) + " " + str(joint.initialState) + " "
                elif isinstance(joint, PrismaticJoint):
                    save += "PrismaticJoint " + str(joint.numSides) + " " + str(joint.r) + " " + str(joint.neutralLength) + " " + str(joint.numLayers) + " " + str(joint.coneAngle) + " " + str(joint.initialState) + " "
                elif isinstance(joint, Tip):
                    save += "Tip " + str(joint.numSides) + " " + str(joint.r) + " " + str(joint.neutralLength) + " " + str(joint.forward) + " "
                else:
                    raise Exception("Not Implemented")
                save += "[" + ''.join([str(x) + "," for x in joint.Pose.A.reshape((16,)).tolist()])
                save += "\n"
            
            f.write(save)
            f.close()
    
    
def loadKinematicChain(filename : str):
    def getJoint(line):
        first = line.split(' ')
        pose = SE3(np.array([float(x) for x in line.split('[')[1].split(",")[:-1]]).reshape(4,4))
        match first[1]:
            case "Waypoint":
                numSides = int(first[2])
                r = float(first[3])
                pathIndex = int(first[4])
                return Waypoint(numSides, r, pose, pathIndex)
            case "RevoluteJoint":
                numSides = int(first[2])
                r = float(first[3])
                totalBendingAngle = float(first[4])
                numSinkLayers = int(first[5])
                initialState = float(first[6])
                return RevoluteJoint(numSides, r, totalBendingAngle, pose, numSinkLayers, initialState)
            case "ExtendedRevoluteJoint":
                numSides = int(first[2])
                r = float(first[3])
                totalBendingAngle = float(first[4])
                tubeLength = float(first[5])
                numSinkLayers = int(first[6])
                initialState = float(first[7])
                return ExtendedRevoluteJoint(numSides, r, totalBendingAngle, tubeLength, pose, numSinkLayers, initialState)
            case "PrismaticJoint":
                numSides = int(first[2])
                r = float(first[3])
                neutralLength = float(first[4])
                numLayers = int(first[5])
                coneAngle = float(first[6])
                initialState = float(first[7])
                return PrismaticJoint(numSides, r, neutralLength, numLayers, coneAngle, pose, initialState)
            case "Tip":
                numSides = int(first[2])
                r = float(first[3])
                length = float(first[4])
                closesForward = bool(first[5])
                return Tip(numSides, r, pose, length, closesForward)

        raise Exception(f"{first[1]} not implemented in save")
        
    try:
        with open(f"save/{filename}.chain") as f:
            lines = f.readlines()
            chain = KinematicChain(getJoint(lines[1]), float(lines[0]))
            for i in range(2, len(lines)):
                parent = int(lines[i].split(" ")[0])
                chain.addJoint(parent, getJoint(lines[i]), relative=False, fixedPosition=True, fixedOrientation=True, safe=False)
            
            return chain
    except Exception as e:
        print(e)
        raise Exception(f"file save/{filename}.chain does not exist")


def chainWithJointDeleted(chain : KinematicChain, jointIndex : int) -> KinematicChain:
    assert(len(chain.Joints) > 1)
    if jointIndex == 0:
        newChain = KinematicChain(chain.Joints[1])
        for i in range(2, len(chain.Joints)):
            newChain.append(chain.Joints[i], relative=False, fixedPosition=True, 
                            fixedOrientation=True, safe=False)
        return newChain
    else:
        newChain = KinematicChain(chain.Joints[0])
        for i in range(1, jointIndex):
            newChain.append(chain.Joints[i], relative=False, fixedPosition=True, 
                            fixedOrientation=True, safe=False)
        for i in range(jointIndex+1, len(chain.Joints)):
            newChain.append(chain.Joints[i], relative=False, fixedPosition=True, 
                            fixedOrientation=True, safe=False)
        return newChain