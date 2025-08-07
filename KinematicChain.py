# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 14:07:27 2023

@author: Daniel Feshbach
"""
from Joint import *
from OrigamiJoint import *
from KinematicTree import *
from TubularPattern import *
from LinkCSC import LinkCSC
from PrintedJoint import *

"""
A KinematicTree with no branching.
"""
class KinematicChain(KinematicTree):
    def __init__(self, startJoint : Joint, maxAnglePerElbow : float = np.pi/2, gimbal : bool = False):
        super().__init__(startJoint, maxAnglePerElbow)
        if gimbal:
            self.boundingBall.expandToCenterOnLine(Line(startJoint.Pose.t, startJoint.Pose.R[:,2]))
            self.nestedBallsRelativeToJoints = []
            self.addNestedBall()
            self.appendOuterWaypoint()
            self.addNestedBall()
    
    """ Add the given joint to the end of the chain, return its index """
    def append(self, newJoint : Joint, relative : bool = True, 
                 fixedPosition : bool = False, fixedOrientation : bool = False, 
                 safe : bool = True, chooseXhatToMinPath : bool = False) -> int:
        parentIndex = len(self.Joints) - 1
        return super().addJoint(parentIndex, newJoint, relative, fixedPosition,
                                fixedOrientation, safe, None, chooseXhatToMinPath)
    
    def appendGlobalFixed(self, newJoint : Joint) -> int:
        parentIndex = len(self.Joints) - 1
        return super().addJoint(parentIndex, newJoint, relative=False, 
                                fixedPosition=True, fixedOrientation=True, safe=False)
    
    def creasePattern(self, twistPortion : float = 0.2) -> TubularPattern:
        chainPattern = copy.deepcopy(self.Joints[0].pattern)
        numSides = self.Joints[0].numSides
        for j in range(1, len(self.Joints)):
            chainPattern.append(self.Links[j].creasePattern(numSides, twistPortion))
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
    

    # returns the index of the last joint that is not a waypoint, 
    # or 0 (root) if there are no non-waypoint joints
    def lastRealJointIndex(self) -> Joint:
        # loop over the joint indices in reverse order
        for i in range(len(self.Joints)-1, -1, -1):
            if not isinstance(self.Joints[i], Waypoint) and not isinstance(self.Joints[i], PrintedWaypoint):
                return i
        return 0


    # add a waypoint on the intersection of the bounding ball and 
    # the parent's z axis, facing outwards
    def appendOuterWaypoint(self, parent = None):
        if parent is None:
            parent = self.Joints[self.lastRealJointIndex()]
        
        parentZhat = parent.Pose.R[:,2]

        if not Line(parent.Pose.t, parentZhat).contains(self.boundingBall.c):
            raise ValueError("Bounding ball is not centered on the parent's Z axis")
        
        if not Ray(parent.Pose.t, -parentZhat).contains(self.boundingBall.c):
            raise ValueError("Parent zhat is facing inwards, not outwards")

        outerPoint = self.boundingBall.c + parentZhat * self.boundingBall.r
        # make sure it's at least 4r in the parentZhat direction from the bounding sphere of the parent
        parentBB = parent.boundingBall()
        if norm(outerPoint - parentBB.c) < parentBB.r + 4*parent.r:
            outerPoint = parentBB.c + (parentBB.r + 4*parent.r) * parentZhat
        
        outwardPoseOnBoundary = SE3.Rt(parent.Pose.R, outerPoint)

        # create a new waypoint at the intersection
        if isinstance(parent, PrintedJoint):
            newWaypoint = PrintedWaypoint(parent.r, outwardPoseOnBoundary, screwRadius=parent.screwRadius, pathIndex=2)
        elif isinstance(parent, OrigamiJoint):
            newWaypoint = Waypoint(parent.numSides, parent.r, outwardPoseOnBoundary, pathIndex=2)

        self.appendGlobalFixed(newWaypoint)
    
    def addNestedBall(self, verify : bool = True):
        self.nestedBallsRelativeToJoints.append(self.boundingBall.newBallTransformedBy(self.Joints[-1].Pose.inv()))
        if verify:
            self.checkBallsAreNested()
    
    
    def nestedBallsGlobal(self):
        return [self.nestedBallsRelativeToJoints[i].newBallTransformedBy(self.Joints[i].Pose) for i in range(len(self.Joints))]
    
    def checkBallsAreNested(self):
        if len(self.nestedBallsRelativeToJoints) != len(self.Joints):
            raise ValueError("Length mismatch between joints and nested balls")
        globalBalls = self.nestedBallsGlobal()
        for i in range(len(globalBalls)-1):
            if not globalBalls[i+1].containsBall(globalBalls[i]):
                raise ValueError("Nested balls are not nested")
            

    def appendGeneralizedGimbal(self, newJoint : Joint, relative : bool = False, addOutwardWaypoint : bool = True) -> int:
        lrji = self.lastRealJointIndex()
        lastRealJoint = self.Joints[lrji]
        
        """
        ballsGlobal = self.nestedBallsGlobal()
        if not Line(lastRealJoint.Pose.t, lastRealJoint.Pose.R[:,2]).contains(ballsGlobal[lrji].c):
            raise ValueError("Bounding ball is not centered on the parent's Z axis")
        """
                
        if relative:
            newJoint.transformPoseIntoFrame(lastRealJoint.Pose)
        
        # make sure the chain goes right up to the bounding sphere and faces outwards
        endJoint = self.Joints[-1]
        if not endJoint.boundingBall().isTangentToBall(self.boundingBall) and \
                Ray(endJoint.Pose.t, -endJoint.pathDirection()).contains(self.boundingBall.c):
            raise ValueError("Bounding ball is not tangent to the end joint's bounding sphere")

        # reverse the new joint's zhat if that makes it more aligned with the previous joint's path direction
        endDir = endJoint.pathDirection()
        if np.dot(endDir, newJoint.Pose.R[:,2]) < np.dot(endDir, -newJoint.Pose.R[:,2]):
            newJoint.reverseZhat()
        newZhat = newJoint.Pose.R[:,2]

        tangentPoint = self.boundingBall.c + newZhat * self.boundingBall.r
        tangentPlane = Plane(tangentPoint, newZhat)
        arcToZhat = arcToDirection(startPoint=endJoint.DistalFrame().t, startDir=endDir, endDir=newZhat, r=endJoint.r)
        assert(norm(arcToZhat.endTangent - newZhat) < 1e-8)

        intersect = tangentPlane.intersectionWithLine(Line(arcToZhat.endPoint, newZhat))
        waypointPose = SE3.Rt(newJoint.Pose.R, intersect)

        # create a new waypoint at the intersection
        if isinstance(endJoint, PrintedJoint):
            newWaypoint = PrintedWaypoint(endJoint.r, waypointPose, screwRadius=endJoint.screwRadius, pathIndex=2)
        elif isinstance(endJoint, OrigamiJoint):
            newWaypoint = Waypoint(endJoint.numSides, endJoint.r, waypointPose, pathIndex=2)

        self.appendGlobalFixed(newWaypoint)
        self.addNestedBall()

        # if the new joint is prismatic, make sure it's fully expanded
        if isinstance(newJoint, PrismaticJoint) or isinstance(newJoint, PrintedPrismaticJoint):
            newJoint.state = newJoint.stateRange()[1]

        newJoint = moveJointNearNeighborBut4rPastPlane(newJoint, newWaypoint, tangentPlane)
        self.append(newJoint, relative=False, fixedPosition=True, fixedOrientation=False, safe=False,
                    chooseXhatToMinPath=True)
                    
        self.boundingBall.expandToCenterOnLine(Line(newJoint.Pose.t, newZhat))
        newJointIndex = len(self.Joints)-1
        self.addNestedBall()

        if addOutwardWaypoint:
            # add a waypoint on the intersection of the bounding ball and 
            # the new joint's z axis, facing outwards
            self.appendOuterWaypoint(parent=newJoint)
            self.addNestedBall()
        
        return newJointIndex
        





            
