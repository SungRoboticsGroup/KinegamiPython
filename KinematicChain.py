# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 14:07:27 2023

@author: Daniel Feshbach
"""
from KinematicTree import *

"""
A KinematicTree with no branching.
"""
class KinematicChain(KinematicTree[F]):
    def __init__(self, startJoint : Joint, maxAnglePerElbow : float = np.pi/2, gimbal : bool = False,
                 joints : Optional[list[Joint]] = None, links : Optional[list[LinkCSC]] = None, 
                 parents : Optional[list[int]] = None, children : Optional[list[list[int]]] = None, 
                 boundingBall : Optional[Ball] = None, units : str = "Centimeter (cm)"):
        super().__init__(startJoint, maxAnglePerElbow, joints, links, parents, children, boundingBall, units)
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
    
    def delete(self, jointIndex : int, safe : bool = True) -> bool:
        assert(jointIndex>=0)
        if safe:
            backup = self.dataDeepCopy()
            try:
                self.delete(jointIndex, safe=False)
                return True
            except ValueError as err:
                print("WARNING: something went wrong in delete:")
                print(err)
                print("Deletion canceled.")
                self.setTo(backup)
                return False
        else:
            nextJoint = self.Joints[jointIndex+1]
            prevJoint = self.Joints[jointIndex-1] if jointIndex>0 else nextJoint
            link_constructor = self._get_link_constructor()
            newLink = link_constructor(self.r, prevJoint.DistalDubinsFrame(), 
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
    

    # returns the index of the last joint that is not a waypoint, 
    # or 0 (root) if there are no non-waypoint joints
    def lastRealJointIndex(self) -> int:
        # loop over the joint indices in reverse order
        for i in range(len(self.Joints)-1, -1, -1):
            if not isinstance(self.Joints[i], Waypoint):
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
        newWaypoint = Waypoint(parent.r, outwardPoseOnBoundary, pathIndex=2)
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
        newWaypoint = Waypoint(endJoint.r, waypointPose, pathIndex=2)

        self.appendGlobalFixed(newWaypoint)
        self.addNestedBall()

        # if the new joint is prismatic, make sure it's fully expanded
        if isinstance(newJoint, Prismatic):
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
        





            
