# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 21:54:54 2023

@author: dfesh
"""
from spatialmath import SE3
from abc import ABC, abstractmethod
from geometryHelpers import *
import matplotlib.pyplot as plt
from TubularPattern import *
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
from typing import Union, Any
from numpy.typing import ArrayLike, NDArray
from types import ModuleType

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
        self.initialState = initialState
        self.TransformStateTo(initialState)

        self.collisionCapsules = self.getCapsules()
        self.proximalDubins = self.ProximalDubinsFrame()
        self.distalDubins = self.DistalDubinsFrame()
    
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
    
    @abstractmethod
    def sdf(self, point: ArrayLike, xp: ModuleType = np) -> Union[float, ArrayLike]:
        pass
    
    def boundingBox(self, xp: ModuleType = np) -> ArrayLike:
        """
        Potentially overestimated axis-aligned bounding box of the joint,
        computed based on the spheres at the proximal, center, and distal positions.
        
        Parameters:
        -----------
        xp : ModuleType
            Numerical module (e.g., numpy or cupy)
            
        Returns:
        --------
        np.ndarray
            Bounding box as [(min_x, min_y, min_z), 
                             (max_x, max_y, max_z)]
        """
        # The number of points sampled along each arc is chosen to ensure they are spaced by at most tolerance*r
        points = xp.vstack((xp.asarray(self.proximalPosition()).reshape(1,3),
                            xp.asarray(self.Pose.t).reshape(1,3),
                            xp.asarray(self.distalPosition()).reshape(1,3)))
        min_corner = xp.min(points, axis=0) - self.r
        max_corner = xp.max(points, axis=0) + self.r
        return xp.vstack((min_corner, max_corner))

    
    def copy(self):
        return Joint(self.r, self.neutralLength, self.Pose, self.initialState)
    
    def ProximalFrame(self) -> SE3:
        return SE3.Trans(-(self.neutralLength/2) * self.pathDirection()) @ self.Pose
    
    def NeutralDistalFrame(self) -> SE3:
        return SE3.Trans((self.neutralLength/2) * self.pathDirection()) @ self.Pose
    
    def StateTransformationFromNeutral(self) -> SE3:
        return self.stateChangeTransformation(self.state)
    
    def DistalFrame(self) -> SE3:
        return self.StateTransformationFromNeutral() @ self.NeutralDistalFrame()        
    
    def TransformStateTo(self, state : float) -> SE3:
        minState, maxState = self.stateRange()
        assert(minState <= state and state <= maxState)
        stateChange = state - self.state
        transformation = self.stateChangeTransformation(stateChange)
        self.state = state
        return transformation
    
    def pathDirection(self) -> np.ndarray:
        return self.Pose.R[:,self.pathIndex()]
    
    def reverseZhat(self):
        self.Pose = self.Pose @ SE3.Rx(np.pi)
        self.proximalDubins = self.ProximalDubinsFrame()
        self.distalDubins = self.DistalDubinsFrame()
    
    def reversePathDirection(self):
        if self.pathIndex() == 2:
            self.Pose = self.Pose @ SE3.Rx(np.pi)
        else:
            self.Pose = self.Pose @ SE3.Rz(np.pi)
        
        self.proximalDubins = self.ProximalDubinsFrame()
        self.distalDubins = self.DistalDubinsFrame()

    # Indices 0,1,2,3 with 0,1,2 cycled to begin with pathDirection
    def dubinsColumnOrder(self) -> np.ndarray:
        return np.hstack((np.roll(np.arange(3), -self.pathIndex()),[3]))
    
    # Pose with axes cycled so that the first axis direction is pathDirection
    def DubinsFrame(self) -> SE3:
        return SE3(self.Pose.A[:,self.dubinsColumnOrder()])
    
    def ProximalDubinsFrame(self) -> SE3:
        return SE3(self.ProximalFrame().A[:,self.dubinsColumnOrder()])
    
    def DistalDubinsFrame(self) -> SE3:
        return SE3(self.DistalFrame().A[:,self.dubinsColumnOrder()])
    
    def proximalPosition(self) -> np.ndarray:
        return self.ProximalFrame().t
    
    def distalPosition(self) -> np.ndarray:
        return self.DistalFrame().t
    
    def transformPoseIntoFrame(self, Frame : SE3):
        self.Pose = Frame @ self.Pose
        self.proximalDubins = self.ProximalDubinsFrame()
        self.distalDubins = self.DistalDubinsFrame()

    def transformPoseBy(self, Transformation: SE3):
        self.Pose = Transformation @ self.Pose
        self.proximalDubins = self.ProximalDubinsFrame()
        self.distalDubins = self.DistalDubinsFrame()
    
    def applyTransformationToPose(self, Transformation : SE3):
        self.Pose = self.Pose @ Transformation
        self.proximalDubins = self.ProximalDubinsFrame()
        self.distalDubins = self.DistalDubinsFrame()

    def translateAlongZ(self, zChange : float):
        self.Pose = self.Pose @ SE3.Trans([0,0,zChange])
        self.proximalDubins = self.ProximalDubinsFrame()
        self.distalDubins = self.DistalDubinsFrame()

    def rotateAboutZ(self, angleToRotateAboutZ):
        self.applyTransformationToPose(SE3.Rz(angleToRotateAboutZ))

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
        self.proximalDubins = self.ProximalDubinsFrame()
        self.distalDubins = self.DistalDubinsFrame()

    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=jointColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=False, 
             axisScale=10, showPoses=True):
        if showAxis:
            zhat = self.Pose.R[:,2]
            JointAxis = np.array([self.Pose.t - axisScale*self.r*zhat,
                                  self.Pose.t + axisScale*self.r*zhat])
            ax.plot(JointAxis[:,0], JointAxis[:,1], JointAxis[:,2], 
                    linestyle='--', color='silver')
        if showSphere:
            self.boundingBall().addToPlot(ax, color=sphereColor, alpha=0.05)
        if showPoses:
            Poses = np.array([self.ProximalFrame(), self.DistalFrame(), self.Pose])
            oColors = np.array([proximalColor, distalColor, centerColor])
            plotHandles = addPosesToPlot(Poses, ax, self.r, 
                                         xColor, yColor, zColor, oColors)
        else:
            plotHandles = None
        return plotHandles
        
    
    def show(self, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, surfaceColor='m',
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=False,
             axisScale=jointAxisScaleDefault, showPoses=True, block=blockDefault):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandles = self.addToPlot(ax, xColor, yColor, zColor,
                                     proximalColor, centerColor, distalColor,
                                     sphereColor, showSphere, surfaceColor,
                                     surfaceOpacity, showSurface, showAxis,
                                     axisScale, showPoses)
        if showPoses:
            xHats, yHats, zHats, origins = plotHandles
            ax.legend([xHats, yHats, zHats], [r'$\^x$', r'$\^y$', r'$\^z$'])
        ax.set_aspect('equal')
        plt.show(block=block)

    def getCapsules(self):
        return []

    def recomputeCollisionCapsules(self):
        self.collisionCapsules = self.getCapsules()
    
    def interpolate(self, density: float, xp : ModuleType = np) -> ArrayLike:
        """
        Sample points along the joint centerline.
        
        Parameters:
        -----------
        density : float
            Number of points per unit length
            
        Returns:
        --------
        np.ndarray
            Sampled points, shape (N, 3)
        """
        v1 = xp.asarray(self.Pose.t - self.proximalPosition())
        v2 = xp.asarray(self.distalPosition() - self.Pose.t)
        n1 = max(2, int(xp.ceil(density * xp.linalg.norm(v1))))
        n2 = max(2, int(xp.ceil(density * xp.linalg.norm(v2))))
        result = xp.empty((n1 + n2 - 1, 3))
        result[:n1] = xp.linspace(0, 1, n1).reshape(-1, 1) * v1.reshape(1, 3) + xp.asarray(self.proximalPosition()).reshape(1, 3)
        result[n1:] = (xp.linspace(0, 1, n2).reshape(-1, 1) * v2.reshape(1, 3) + xp.asarray(self.Pose.t).reshape(1, 3))[1:]
        return result
    
    def startCircle(self, forward : bool = True) -> Circle3D:
        """
        Return the circular face at the proximal end of the joint.
        
        Returns:
        --------
        Circle3D
            The circular disc at the proximal position, oriented with normal
            pointing along the path direction
        """
        return Circle3D(
            radius=self.r,
            center=self.proximalPosition(),
            normal=self.pathDirection() if forward else -self.pathDirection(),
            radialVector=self.ProximalFrame().R[:,1]
        )
    
    def endCircle(self, forward : bool = True) -> Circle3D:
        """
        Return the circular face at the distal end of the joint.
        
        Returns:
        --------
        Circle3D
            The circular disc at the distal position, oriented with normal
            pointing along the path direction
        """
        return Circle3D(
            radius=self.r,
            center=self.distalPosition(),
            normal=self.pathDirection() if forward else -self.pathDirection(),
            radialVector=self.DistalFrame().R[:,1]
        )

    def length(self) -> np.floating[Any]:
        """Return the current length of the joint from proximal to center to distal."""
        return np.linalg.norm(self.Pose.t - self.proximalPosition()) +\
               np.linalg.norm(self.distalPosition() - self.Pose.t)