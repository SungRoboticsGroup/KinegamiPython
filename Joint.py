# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 21:54:54 2023

@author: dfesh
"""
from spatialmath import SE3
from abc import ABC, abstractmethod
from geometryHelpers import *
import matplotlib.pyplot as plt
from tubularOrigami import *

class Joint(ABC):
    """
    r is the tubular radius
    Pose is an SE3 object where the Z axis is the joint axis, representing
        the center in the relaxed configuration
    neutralLength is from proximal to distal position (in state=0 configuration)
    initialState (defaults to 0) is the initial value for the joint state
    """
    def __init__(self, r : float, neutralLength : float, Pose : SE3, 
                 initialState : float = 0):
        self.r = r
        self.Pose = Pose
        self.neutralLength = neutralLength
        self.state = 0
        self.transformStateTo(initialState)
    
    @abstractmethod #0 for xhat, 2 for zhat
    def pathIndex(self) -> int:
        pass
    
    @abstractmethod
    def stateChangeTransformation(self, stateChange : float) -> SE3:
        pass
    
    @abstractmethod
    def stateRange(self) -> list:
        pass
    
    @abstractmethod
    def boundingRadius(self) -> float:
        pass
    
    @abstractmethod
    def boundingBall(self) -> Ball:
        pass
    
    def proximalPose(self) -> SE3:
        return SE3.Trans(-(self.neutralLength/2) * self.pathDirection()) @ self.Pose
    
    def neutralDistalPose(self) -> SE3:
        return SE3.Trans((self.neutralLength/2) * self.pathDirection()) @ self.Pose
    
    def stateTransformationFromNeutral(self) -> SE3:
        return self.stateChangeTransformation(self.state)
    
    def distalPose(self) -> SE3:
        return self.stateTransformationFromNeutral() @ self.neutralDistalPose()        
    
    def transformStateTo(self, state : float) -> SE3:
        minState, maxState = self.stateRange()
        assert(minState <= state and state <= maxState)
        stateChange = state - self.state
        self.state = state
        return self.stateChangeTransformation(stateChange)
    
    def pathDirection(self) -> np.ndarray:
        return self.Pose.R[:,self.pathIndex()]
    
    # Indices 0,1,2,3 with 0,1,2 cycled to begin with pathDirection
    def dubinsColumnOrder(self) -> np.ndarray:
        return np.hstack((np.roll(np.arange(3), -self.pathIndex()),[3]))
    
    # Pose with axes cycled so that the first axis direction is pathDirection
    def dubinsFrame(self) -> SE3:
        return SE3(self.Pose.A[:,self.dubinsColumnOrder()])
    
    def proximalDubinsFrame(self) -> SE3():
        return SE3(self.proximalPose().A[:,self.dubinsColumnOrder()])
    
    def distalDubinsFrame(self) -> SE3():
        return SE3(self.distalPose().A[:,self.dubinsColumnOrder()])
    
    def proximalPosition(self) -> np.ndarray:
        return self.proximalPose().t
    
    def distalPosition(self) -> np.ndarray:
        return self.distalPose().t    
    
    # Frame should be an SE3 object
    def transformPoseFromFrameToGlobal(self, Frame : SE3):
        self.Pose = self.Pose @ Frame
    
    def transformPoseBy(self, Transformation : SE3):
        self.Pose = Transformation @ self.Pose
    
    def translateAlongZ(self, zChange : float):
        self.Pose = self.Pose @ SE3.Trans([0,0,zChange])
    
    def setXhatAboutZhat(self, xhatNew):
        xhatNew = xhatNew / norm(xhatNew)
        zhat = self.Pose.R[:,2]
        assert(dot(zhat, xhatNew) < 0.0001) #input must be orthogonal to Z axis
        yhatNew = cross(zhat, xhatNew)
        Transform = np.eye(4)
        Transform[0:3,0] = xhatNew
        Transform[0:3,1] = yhatNew
        Transform[0:3,2] = zhat
        Transform[0:3,3] = self.Pose.t
        self.Pose = SE3(Transform)
    
    def addToPlot(self, ax, xColor='r', yColor='b', zColor='g', 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor='black', showSphere=False):
        Poses = np.array([self.proximalPose(), self.distalPose(), self.Pose])
        oColors = np.array([proximalColor, distalColor, centerColor])
        plotHandles = addPosesToPlot(Poses, ax, self.r, 
                                     xColor, yColor, zColor, oColors)
        scale = 5
        zhat = self.Pose.R[:,2]
        JointAxis = np.array([self.Pose.t - scale*self.r*zhat,
                              self.Pose.t + scale*self.r*zhat])
        ax.plot(JointAxis[:,0], JointAxis[:,1], JointAxis[:,2], 
                linestyle='--', color='silver')
        if showSphere:
            self.boundingBall().addToPlot(ax, color=sphereColor, alpha=0.05)
        return plotHandles
        
    
    def plot(self, xColor='r', yColor='b', zColor='g', 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor='black', showSphere=False):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandles = self.addToPlot(ax, xColor, yColor, zColor,
                                     proximalColor, centerColor, distalColor,
                                     sphereColor, showSphere)
        xHats, yHats, zHats, origins = plotHandles
        ax.set_aspect('equal')
        ax.legend([xHats, yHats, zHats], [r'$\^x$', r'$\^y$', r'$\^z$'])
    
class OrigamiJoint(Joint):
    def __init__(self, numSides : int, r : float, neutralLength : float, Pose : SE3(), 
                 initialState : float = 0):
        self.numSides = numSides
        self.polygonInnerAngle = np.pi * (numSides-2)/(2*numSides)
        super().__init__(r, neutralLength, Pose)
    
class RevoluteJoint(OrigamiJoint):
    """
    Origami revolute joint with rotation range [-totalBendingAngle/2, totalBendingAngle/2]
    and numSinkLayers recursive sink gadget layers.
    """
    def __init__(self, numSides : int, r : float, totalBendingAngle : float, 
                 Pose : SE3, numSinkLayers : int = 1,
                 initialState : float = 0):
        polygonInnerAngle = np.pi * (numSides-2)/(2*numSides)
        neutralLength = 2*r*np.sin(polygonInnerAngle)*np.tan(totalBendingAngle/4) #2*delta from paper
        self.totalBendingAngle = totalBendingAngle
        super().__init__(numSides, r, neutralLength, Pose, initialState)
        self.pattern = RevoluteJointPattern(self.numSides, self.r, 
                                            totalBendingAngle, numSinkLayers)
    
    def pathIndex(self) -> int:
        return 0 # xhat
    
    def stateRange(self) -> list:
        return [-self.totalBendingAngle/2, self.totalBendingAngle/2]
    
    def stateChangeTransformation(self, stateChange : float) -> SE3:
        return RotationAboutLine(rotAxisDir=self.Pose.R[:,2],
                              rotAxisPoint=self.Pose.t,
                              angle=stateChange)
    
    def boundingRadius(self) -> float:
        return norm([self.r, self.neutralLength / 2])
    
    def boundingBall(self) -> Ball:
        return Ball(self.Pose.t, self.boundingRadius())
        
class PrismaticJoint(OrigamiJoint):
    def __init__(self, numSides : int, r : float, neutralLength : float, 
                 numLayers : int, coneAngle : float, Pose : SE3, 
                 initialState : float = 0):
        neutralLayerHeight = neutralLength / numLayers
        flatLayerHalfHeight = neutralLayerHeight / (2*np.sin(coneAngle))
        self.minLength = numLayers*flatLayerHalfHeight
        self.maxLength = 2*self.minLength
        super().__init__(numSides, r, neutralLength, Pose, initialState)
        self.pattern = PrismaticJointPattern(numSides, r, neutralLength, 
                                             numLayers, coneAngle)
        
    def pathIndex(self) -> int:
        return 2 # zhat
    
    def stateChangeTransformation(self, stateChange : float) -> SE3:
        return SE3.Trans(stateChange * self.pathDirection())
    
    def stateRange(self) -> list:
        return [self.minLength - self.neutralLength, self.maxLength - self.neutralLength]
    
    def length(self) -> float:
        return self.neutralLength + self.state
    
    def boundingRadius(self) -> float:
        return norm([self.r, self.length() / 2])
    
    def center(self) -> np.ndarray:
        return self.Pose.t + (self.state/2) * self.pathDirection()
    
    def boundingBall(self) -> Ball:
        return Ball(self.center(), self.boundingRadius())
    
class WayPoint(OrigamiJoint):
    # path direction through a waypoint defaults to zhat
    def __init__(self, numSides : int, r : float, Pose : SE3, pathIndex : int = 2):
        super().__init__(numSides, r, 0, Pose)
        self.pidx = pathIndex
        self.pattern = TubularPattern(numSides, r)
    
    def pathIndex(self) -> int:
        return self.pidx
    
    def stateChangeTransformation(self, stateChange : float) -> SE3:
        return SE3()
    
    def stateRange(self) -> list:
        return [0,0]
    
    def boundingRadius(self) -> float:
        return self.r
    
    def boundingBall(self) -> Ball:
        return Ball(self.Pose.t, self.boundingRadius())
    
