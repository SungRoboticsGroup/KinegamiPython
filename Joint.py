# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 21:54:54 2023

@author: dfesh
"""
from spatialmath import SE3
from abc import ABC, abstractmethod
from TubularPattern import jointAxisScaleDefault, jointColorDefault, jointEdgeColorDefault, sphereColorDefault, surfaceOpacityDefault, xColorDefault, yColorDefault, zColorDefault
from geometryHelpers import *
import matplotlib.pyplot as plt
from TubularPattern import *
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
from typing import Union, Any
from numpy.typing import ArrayLike, NDArray
from types import ModuleType
from geometryHelpers import jointAxisScaleDefault, jointColorDefault, jointEdgeColorDefault, sphereColorDefault, surfaceOpacityDefault, xColorDefault, yColorDefault, zColorDefault
import pyqtgraph.opengl as gl

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
    def boundingRadius(self) -> float | np.floating[Any]:
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
        return Joint(self.r, self.neutralLength, self.Pose, self.state)
    
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

    def distalPathDirection(self) -> np.ndarray:
        return self.DistalFrame().R[:,self.pathIndex()]
    
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
             surfaceColor=jointColorDefault, edgeColor=jointEdgeColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=False, 
             axisScale=jointAxisScaleDefault, showPoses=True):
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
    
    def addToWidget(self, widget, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                    proximalColor=proximalColorDefault, centerColor=centerColorDefault, distalColor=distalColorDefault,
                    sphereColor=sphereColorDefault, showSphere=False, surfaceColor=jointColorDefault, 
                    showSurface=True, showAxis=False, axisScale=jointAxisScaleDefault, showPoses=True, poseAxisScaleMultipler=None):
        if showAxis:
            zhat = self.Pose.R[:, 2] 
            jointAxis = np.array([self.Pose.t - 10 * self.r * zhat,
                                self.Pose.t + 10 * self.r * zhat])
            line_item = gl.GLLinePlotItem(pos=jointAxis, color=showAxisColor, width=2, antialias=True)  # Using a silver color
            widget.plot_widget.addItem(line_item)

        if showPoses:
            for pose, color in zip([self.ProximalFrame(), self.DistalFrame(), self.Pose], [proximalColor, distalColor, centerColor]):
                for i, axis_color in enumerate([xColor, yColor, zColor]):
                    poseAxisScale = self.r
                    if poseAxisScaleMultipler:
                        poseAxisScale *= poseAxisScaleMultipler
                    start_point = pose.t
                    end_point = start_point + poseAxisScale * pose.R[:, i]
                    points = np.array([start_point, end_point])
                    line = gl.GLLinePlotItem(pos=points, color=axis_color, width=2, antialias=True)
                    widget.plot_widget.addItem(line)

        if showSphere:
            self.boundingBall().addToWidget(widget, sphereColor)
    
    def addArrows(self, widget, selectedArrow=-1, local=True, frame: SE3 = None, mode=""):
        import math
        desired_px = 80
        thickness_px = 10
        rad = widget.plot_widget.world_length_for_pixel_length(desired_px)
        arrow_thick = widget.plot_widget.world_length_for_pixel_length(thickness_px)
        colors = rotateArrowColors
        center = self.Pose.t

        # choose axes
        if local:
            axes = [self.Pose.R[:,i] for i in range(3)]
        else:
            axes = np.eye(3)
        if frame:
            axes = [frame.R[:,i] for i in range(3)]

        if mode == "Translate":
            for i, a in enumerate(axes):
                col = selectedArrowColor if i==selectedArrow else colors[i]
                start = center
                end = center + rad * a
                from printedGUI import OverlayLine
                widget.plot_widget.addItem(
                    OverlayLine(pos=np.array([start,end]), color=col, width=8, antialias=True)
                )

        elif mode == "Rotate":
            for i, axis in enumerate(axes):
                helper = np.array([1,0,0])
                if abs(np.dot(axis,helper))>0.9:
                    helper = np.array([0,1,0])
                u = np.cross(axis, helper); u/=np.linalg.norm(u)
                v = np.cross(axis, u)
                pts = np.array([
                    center + (rad+arrow_thick)*(u*np.cos(t)+v*np.sin(t))
                    for t in np.linspace(0,2*math.pi,64)
                ])
                col = selectedArrowColor if i==selectedArrow else colors[i]
                from printedGUI import OverlayLine
                widget.plot_widget.addItem(
                    OverlayLine(pos=pts, color=col, width=8, antialias=True)
                )

    def addTranslateArrows(self, widget, selectedArrow=-1, local=True, frame : SE3 = None):
        self.addArrows(widget, selectedArrow, local, frame, mode="Translate")
    
    def addRotateArrows(self, widget, selectedArrow=-1, local=True, frame : SE3 = None):
        self.addArrows(widget, selectedArrow, local, frame, mode="Rotate")
    
    def show(self, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=jointColorDefault, edgeColor=jointEdgeColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=False,
             axisScale=jointAxisScaleDefault, showPoses=True, block=blockDefault):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandles = self.addToPlot(ax=ax, xColor=xColor, yColor=yColor, zColor=zColor,
                                     proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor,
                                     sphereColor=sphereColor, showSphere=showSphere, 
                                     surfaceColor=surfaceColor, edgeColor=edgeColor,
                                     surfaceOpacity=surfaceOpacity, showSurface=showSurface, showAxis=showAxis,
                                     axisScale=axisScale, showPoses=showPoses)
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
            normal=self.distalPathDirection() if forward else -self.distalPathDirection(),
            radialVector=self.DistalFrame().R[:,1]
        )

    def length(self) -> float | np.floating[Any]:
        """Return the current length of the joint from proximal to center to distal."""
        return np.linalg.norm(self.Pose.t - self.proximalPosition()) +\
               np.linalg.norm(self.distalPosition() - self.Pose.t)


class Prismatic(Joint):
    def __init__(self, r : float, neutralLength : float, minLength : float, maxLength : float, Pose : SE3, initialState : float = 0):        
        # Validity checks
        if minLength < 0 or maxLength < 0:
            raise ValueError("Prismatic joint lengths must be non-negative")
        if minLength >= maxLength:
            raise ValueError("Prismatic joint minLength must be less than maxLength")
        if neutralLength + initialState < minLength or neutralLength + initialState > maxLength:
            raise ValueError("Prismatic joint neutralLength + initialState must be within [minLength, maxLength]")
        if neutralLength < minLength or neutralLength > maxLength:
            raise ValueError("Prismatic joint neutralLength must be within [minLength, maxLength]")

        self.minLength = minLength
        self.maxLength = maxLength
        super().__init__(r, neutralLength, Pose, initialState)
    
    def __repr__(self):
        return f"Prismatic(r={repr(self.r)}, neutralLength={repr(self.neutralLength)}, minLength={repr(self.minLength)}, maxLength={repr(self.maxLength)}, Pose={repr(self.Pose)}, initialState={repr(self.state)})"

    def pathIndex(self) -> int:
        return 2 # zhat
    
    def stateChangeTransformation(self, stateChange : float) -> SE3:
        return SE3.Trans(stateChange * self.pathDirection())
    
    def stateRange(self) -> list:
        return [self.minLength - self.neutralLength, self.maxLength - self.neutralLength]
    
    def length(self) -> float | np.floating[Any]:
        return self.neutralLength + self.state
    
    def boundingRadius(self) -> float | np.floating[Any]:
        return np.linalg.norm([self.r, self.length() / 2])
    
    def center(self) -> np.ndarray:
        return self.Pose.t + (self.state/2) * self.pathDirection()
    
    def boundingBall(self) -> Ball:
        return Ball(self.center(), self.boundingRadius())
    
    def boundingCylinder(self) -> Cylinder:
        uhat = (self.Pose @ SE3.Rz(np.pi/2)).R[:,1]
        return Cylinder(self.r, self.ProximalFrame().t, self.pathDirection(), 
                        self.length(), uhat)
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=prismaticColorDefault, edgeColor=prismaticEdgeColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=True, 
             axisScale=10, showPoses=True):
        plotHandles = super().addToPlot(ax=ax, xColor=xColor, yColor=yColor, zColor=zColor, 
                          proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor, 
                          sphereColor=sphereColor, showSphere=showSphere,
                          surfaceColor=surfaceColor, edgeColor=edgeColor,
                          surfaceOpacity=surfaceOpacity, showSurface=False, showAxis=showAxis,
                          axisScale=axisScale, showPoses=showPoses)
        if showSurface:
            self.boundingCylinder().addToPlot(ax, color=surfaceColor, 
                                              alpha=surfaceOpacity, 
                                              edgeColor=edgeColor)            
        return plotHandles
    
    def addToWidget(self, widget, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                    proximalColor=proximalColorDefault, centerColor=centerColorDefault, distalColor=distalColorDefault,
                    sphereColor=sphereColorDefault, showSphere=False, 
                    surfaceColor=prismaticColorDefault, 
                    showSurface=True, showAxis=True, axisScale=10, showPoses=True, poseAxisScaleMultipler=None):
        import pyqtgraph.opengl as gl
        from style import prismaticColorList
        
        # Call parent's addToWidget for poses and axis
        super().addToWidget(widget=widget, xColor=xColor, yColor=yColor, zColor=zColor, 
                          proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor, 
                          sphereColor=sphereColor, showSphere=showSphere,
                          surfaceColor=surfaceColor, showSurface=False, showAxis=showAxis,
                          axisScale=axisScale, showPoses=showPoses, poseAxisScaleMultipler=poseAxisScaleMultipler)
        if showSurface:
            # Use surfaceColor override (e.g. collision highlighting) if provided, otherwise default color list
            color = surfaceColor if surfaceColor != prismaticColorDefault else prismaticColorList
            self.boundingCylinder().addToWidget(widget, color_list=color, is_joint=True)
    
    def sdf(self, point: ArrayLike, xp: ModuleType = np) -> Union[float, ArrayLike]:
        # TODO: ADD PLANAR CUTOFFS
        point = xp.asarray(point).reshape(-1,3)
        return sdf_capsule(xp, point, 
                            a = xp.asarray([self.ProximalDubinsFrame().t]).reshape(1,3),
                            b = xp.asarray([self.DistalDubinsFrame().t]).reshape(1,3),
                            r = self.r).flatten()

class Revolute(Joint):
    """
    Docstring for Revolute
    A revolute joint that rotates about a specified axis.
    Attributes:
        r (float): Radius of the joint.
        Pose (SE3): The pose of the joint.
        pathIndex (int): The index of the axis of rotation (0 for x, 1 for y, 2 for z).
        neutralLength (float): The length of the joint.
        minAngle (float): Minimum angle of rotation in radians. Can be None or -np.inf for no limit.
        maxAngle (float): Maximum angle of rotation in radians. Can be None or np.inf for no limit.
        initialState (float): Initial angle of rotation in radians.
    """
    def __init__(self, r : float, Pose : SE3, pathIndex : int, neutralLength : float,
                 minAngle : Optional[float], maxAngle : Optional[float], initialState : float = 0.0):
        if minAngle is None:
            minAngle = -np.inf
        if maxAngle is None:
            maxAngle = np.inf
        if initialState < -np.pi or initialState > np.pi:
            raise ValueError("Revolute joint initialState must be within [-pi, pi]")
        if not pathIndex in [0, 1, 2]:
            raise ValueError("pathIndex must be 0 (x), 1 (y), or 2 (z)")
        if initialState < minAngle or initialState > maxAngle:
            raise ValueError("Revolute joint initialState must be within [minAngle, maxAngle]")
        
        self.minAngle = minAngle
        self.maxAngle = maxAngle
        self.pidx = pathIndex
        super().__init__(r, neutralLength, Pose, initialState)
        

    def pathIndex(self) -> int:
        return self.pidx
    
    def stateRange(self) -> list:
        return [self.minAngle, self.maxAngle]
    
    def stateChangeTransformation(self, stateChange : float) -> SE3:
        return RotationAboutLine(rotAxisDir=self.Pose.R[:,2],
                              rotAxisPoint=self.Pose.t,
                              angle=stateChange)
    
    def boundingRadius(self) -> float | np.floating[Any]:
        return np.linalg.norm([self.r, self.neutralLength / 2])
    
    def boundingBall(self) -> Ball:
        return Ball(self.Pose.t, self.boundingRadius())
    
    def proximalCylinder(self) -> Cylinder:
        return Cylinder(self.r, self.proximalPosition(), self.pathDirection(), 
                            self.neutralLength/2, self.ProximalDubinsFrame().R[:,1])
    
    def distalCylinder(self) -> Cylinder:
        return Cylinder(self.r, self.distalPosition(), -self.distalPathDirection(), 
                            self.neutralLength/2, self.DistalDubinsFrame().R[:,1])
    
    def centerSphere(self) -> Ball:
        return Ball(self.Pose.t, self.r)
    
    def sdf(self, point: ArrayLike, xp: ModuleType = np) -> Union[float, ArrayLike]:
        point = xp.asarray(point).reshape(-1,3)
        # Compute the min of the SDFs of the proximal and distal capsules
        # Cut off the outer hemisphere caps by maxing with signed distances to planes
        # a and b are the capsule starts and ends respectively, so each are (2,3) arrays
        capsule_sdfs = sdf_capsule(xp, point, 
                a = xp.vstack((self.proximalPosition(), self.distalPosition())), 
                b = xp.vstack((self.Pose.t, self.Pose.t)), 
                r = self.r)  # (N, 2)
        plane_sdfs = xp.stack([
            sdf_plane(xp, point, self.proximalPosition(), -self.pathDirection()),
            sdf_plane(xp, point, self.distalPosition(), self.distalPathDirection())
        ], axis=1)  # (N, 2)
        return xp.min(xp.maximum(capsule_sdfs, plane_sdfs), axis=1)  # (N,)

    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=revoluteColorDefault, edgeColor=revoluteEdgeColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=True,
             axisScale=10, showPoses=True):
        plotHandles = super().addToPlot(ax=ax, xColor=xColor, yColor=yColor, zColor=zColor, 
                          proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor, 
                          sphereColor=sphereColor, showSphere=showSphere,
                          surfaceColor=surfaceColor, edgeColor=edgeColor,
                          surfaceOpacity=surfaceOpacity, showSurface=False, showAxis=showAxis,
                          axisScale=axisScale, showPoses=showPoses)
        if showSurface:
            self.proximalCylinder().addToPlot(ax, color=surfaceColor, edgeColor=edgeColor, alpha=surfaceOpacity)
            self.distalCylinder().addToPlot(ax, color=surfaceColor, edgeColor=edgeColor, alpha=surfaceOpacity)
            self.centerSphere().addToPlot(ax, color=surfaceColor, alpha=surfaceOpacity)
        return plotHandles
    
    def addToWidget(self, widget, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                    proximalColor=proximalColorDefault, centerColor=centerColorDefault, distalColor=distalColorDefault,
                    sphereColor=sphereColorDefault, showSphere=False, 
                    surfaceColor=revoluteColorDefault, 
                    showSurface=True, showAxis=True, axisScale=10, showPoses=True, poseAxisScaleMultipler=None):
        import pyqtgraph.opengl as gl
        from style import revoluteColorList
        
        # Call parent's addToWidget for poses and axis
        super().addToWidget(widget=widget, xColor=xColor, yColor=yColor, zColor=zColor, 
                          proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor, 
                          sphereColor=sphereColor, showSphere=showSphere,
                          surfaceColor=surfaceColor, showSurface=False, showAxis=showAxis,
                          axisScale=axisScale, showPoses=showPoses, poseAxisScaleMultipler=poseAxisScaleMultipler)
        if showSurface:
            # Use surfaceColor override (e.g. collision highlighting) if provided, otherwise default color list
            color = surfaceColor if surfaceColor != revoluteColorDefault else revoluteColorList
            self.proximalCylinder().addToWidget(widget, color_list=color, is_joint=True)
            self.distalCylinder().addToWidget(widget, color_list=color, is_joint=True)
            self.centerSphere().addToWidget(widget, color=color)

class TransverseRevolute(Revolute):
    def __init__(self, r : float, Pose : SE3, minAngle : Optional[float] = None, maxAngle : Optional[float] = None, 
                 neutralLength : Optional[float] = None, initialState : float = 0.0,
                 checkCircleOverlap : bool = True):
        if neutralLength is None: 
            if minAngle is None or maxAngle is None:
                raise ValueError("TransverseRevolute joints must have either neutralLength or both minAngle and maxAngle specified")
            # compute neutralLength to achieve min and max angles
            # without the end circles overlapping
            largerAngle = max(abs(minAngle), abs(maxAngle))
            neutralLength = 2*r*np.tan(largerAngle/2)
        elif checkCircleOverlap:
            # check that the provided neutralLength is sufficient
            angleAtWhichCirclesTouch = 2 * np.arctan(neutralLength / (2*r))
            if minAngle is None:
                minAngle = -angleAtWhichCirclesTouch
            if maxAngle is None:
                maxAngle = angleAtWhichCirclesTouch
            if minAngle < -angleAtWhichCirclesTouch-1e-6 or maxAngle > angleAtWhichCirclesTouch+1e-6:
                raise Warning("Provided neutralLength is too small to prevent end circle overlap for the given minAngle and maxAngle")
        Revolute.__init__(self, r, Pose, pathIndex=0, neutralLength=neutralLength, 
                         minAngle=minAngle, maxAngle=maxAngle, initialState=initialState)
    
    def __repr__(self):
        return f"TransverseRevolute(r={repr(self.r)}, Pose={repr(self.Pose)}, minAngle={repr(self.minAngle)}, "+\
                f"maxAngle={repr(self.maxAngle)}, neutralLength={repr(self.neutralLength)}, initialState={repr(self.state)})"

class CoaxialRevolute(Revolute):
    def __init__(self, r : float, Pose : SE3, neutralLength : float, minAngle : Optional[float], 
                 maxAngle : Optional[float], initialState : float = 0.0):
        super().__init__(r, Pose, pathIndex=2, neutralLength=neutralLength, 
                         minAngle=minAngle, maxAngle=maxAngle, initialState=initialState)
    
    def __repr__(self):
        return f"CoaxialRevolute(r={repr(self.r)}, Pose={repr(self.Pose)}, neutralLength={repr(self.neutralLength)}, "+\
                f"minAngle={repr(self.minAngle)}, maxAngle={repr(self.maxAngle)}, initialState={repr(self.state)})"

class Waypoint(Joint):
    # path direction through a waypoint defaults to zhat
    def __init__(self, r : float, Pose : SE3, pathIndex : int = 2):
        assert(pathIndex in [0,1,2])
        self.pidx = pathIndex
        super().__init__(r, 0, Pose, 0)
    
    def __repr__(self):
        return f"Waypoint(r={repr(self.r)}, Pose={repr(self.Pose)}, pathIndex={repr(self.pidx)})"
    
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

    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=linkColorDefault, edgeColor=linkColorDefault,
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
            Poses = np.array([self.Pose])
            oColors = np.array([centerColor])
            plotHandles = addPosesToPlot(Poses, ax, self.r, 
                                         xColor, yColor, zColor, oColors)
        else:
            plotHandles = None
        return plotHandles

    def sdf(self, point: ArrayLike, xp: ModuleType = np) -> Union[float, ArrayLike]:
        point = xp.asarray(point).reshape(-1,3)
        center = xp.asarray(self.Pose.t).reshape(3)
        # SDF of a sphere centered at Pose.t with radius r
        return xp.linalg.norm(point - center, axis=1) - self.r

    def addToWidget(self, widget, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                    proximalColor=proximalColorDefault, centerColor=centerColorDefault, distalColor=distalColorDefault,
                    sphereColor=sphereColorDefault, showSphere=False, surfaceColor=jointColorDefault, 
                    showSurface=True, showAxis=False, axisScale=jointAxisScaleDefault, showPoses=True, poseAxisScaleMultipler=None):
        """Draw the waypoint as a circle with a small dot at the center."""
        circle = Circle3D(self.r, self.Pose.t, self.Pose.R[:, self.pidx])
        circle.addToWidget(widget, color=surfaceColor, width=0.05*self.r)
        # Small dot at the waypoint center
        dot = Ball(self.Pose.t, 0.05 * self.r)
        dot.addToWidget(widget, color=surfaceColor)
        # Call parent's addToWidget for poses and axis
        super().addToWidget(widget=widget, xColor=xColor, yColor=yColor, zColor=zColor, 
                          proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor, 
                          sphereColor=sphereColor, showSphere=showSphere,
                          surfaceColor=surfaceColor, showSurface=False, showAxis=showAxis,
                          axisScale=axisScale, showPoses=showPoses, poseAxisScaleMultipler=poseAxisScaleMultipler)
        

class Tip(Joint):
    def __init__(self, r : float, Pose : SE3, length : float, 
                 closesForward : bool = True, pathIndex : int = 2):
        if not pathIndex in [0,1,2]:
            raise ValueError("pathIndex must be 0 (x), 1 (y), or 2 (z)")
        self.pidx = pathIndex
        super().__init__(r, length, Pose, 0)
        self.forward = closesForward
    
    def __repr__(self):
        return f"Tip(r={repr(self.r)}, Pose={repr(self.Pose)}, length={repr(self.neutralLength)}, closesForward={repr(self.forward)}, pathIndex={repr(self.pidx)})"
    
    def pathIndex(self) -> int:
        return self.pidx
    
    def stateChangeTransformation(self, stateChange : float) -> SE3:
        return SE3()
    
    def stateRange(self) -> list:
        return [0,0]
    
    def boundingRadius(self) -> float | np.floating[Any]:
        return np.linalg.norm([self.r, self.neutralLength/2])
    
    def boundingBall(self) -> Ball:
        return Ball(self.Pose.t, self.boundingRadius())
    
    def sdf(self, point: ArrayLike, xp: ModuleType = np) -> Union[float, ArrayLike]:
        point = xp.asarray(point).reshape(-1,3)
        return xp.maximum(
            sdf_capsule(xp, point, 
                            a = xp.asarray([self.ProximalDubinsFrame().t]).reshape(1,3),
                            b = xp.asarray([self.DistalDubinsFrame().t]).reshape(1,3),
                            r = self.r).flatten(),
            sdf_plane(xp, point, 
                            self.DistalDubinsFrame().t if self.forward else self.ProximalDubinsFrame().t,
                            self.pathDirection() if self.forward else -self.distalPathDirection())
        )

    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                  proximalColor='c', centerColor='m', distalColor='y', 
                  sphereColor=sphereColorDefault, showSphere=False, 
                  surfaceColor=jointColorDefault, edgeColor=jointEdgeColorDefault, 
                  surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=False, 
                  axisScale=jointAxisScaleDefault, showPoses=True):
        plotHandles = super().addToPlot(ax=ax, xColor=xColor, yColor=yColor, zColor=zColor, 
                                 proximalColor=proximalColor, centerColor=centerColor, 
                                 distalColor=distalColor, sphereColor=sphereColor, showSphere=showSphere, 
                                 surfaceColor=surfaceColor, edgeColor=edgeColor, 
                                 surfaceOpacity=surfaceOpacity, showSurface=False, 
                                 showAxis=showAxis, axisScale=axisScale, showPoses=showPoses)
        
        if showSurface:
            # Plot the tip as a stretched hemisphere
            # Generate hemisphere on unit sphere
            n_lat = 20  # latitude divisions
            n_lon = 32  # longitude divisions
            
            # Hemisphere goes from 0 to pi/2 in latitude (or -pi/2 to 0 if closing backward)
            if self.forward:
                theta = np.linspace(0, np.pi/2, n_lat)  # 0 at pole (distal end)
            else:
                theta = np.linspace(np.pi/2, np.pi, n_lat)  # pi at pole (proximal end)
            
            phi = np.linspace(0, 2*np.pi, n_lon)
            theta_grid, phi_grid = np.meshgrid(theta, phi)
            
            # Unit hemisphere points
            x_sphere = np.sin(theta_grid) * np.cos(phi_grid)
            y_sphere = np.sin(theta_grid) * np.sin(phi_grid)
            z_sphere = np.cos(theta_grid)
            
            # Stretch hemisphere: radial (x,y) by self.r, axial (z) by self.neutralLength
            # The hemisphere naturally goes from z=0 to z=1 (or z=0 to z=-1)
            # We need it to span from proximal to distal
            if self.forward:
                # z goes from 1 (at theta=0) to 0 (at theta=pi/2)
                # Map to go from proximalPosition (z=1) to distalPosition (z=0)
                scale_z = self.neutralLength
            else:
                # z goes from 0 (at theta=pi/2) to -1 (at theta=pi)
                # Map to go from distalPosition (z=0) to proximalPosition (z=-1)
                scale_z = self.neutralLength
            
            # Scale the hemisphere
            x_scaled = self.r * x_sphere
            y_scaled = self.r * y_sphere
            z_scaled = scale_z * z_sphere
            
            # Transform to world coordinates based on path direction
            if self.forward:
                base_frame = self.ProximalFrame()
            else:
                base_frame = self.DistalFrame()
            
            # Get basis vectors - need to align z_scaled with path direction
            path_dir = self.pathDirection()
            # Create rotation matrix with path_dir as z-axis
            R = base_frame.R
            
            # Reorder to match pathIndex
            if self.pidx == 0:  # xhat is path direction
                # x_scaled -> path, y_scaled -> y, z_scaled -> z
                x_world = base_frame.t[0] + z_scaled * R[0, 0] + x_scaled * R[0, 1] + y_scaled * R[0, 2]
                y_world = base_frame.t[1] + z_scaled * R[1, 0] + x_scaled * R[1, 1] + y_scaled * R[1, 2]
                z_world = base_frame.t[2] + z_scaled * R[2, 0] + x_scaled * R[2, 1] + y_scaled * R[2, 2]
            elif self.pidx == 1:  # yhat is path direction
                # x_scaled -> x, y_scaled -> path, z_scaled -> z
                x_world = base_frame.t[0] + x_scaled * R[0, 0] + z_scaled * R[0, 1] + y_scaled * R[0, 2]
                y_world = base_frame.t[1] + x_scaled * R[1, 0] + z_scaled * R[1, 1] + y_scaled * R[1, 2]
                z_world = base_frame.t[2] + x_scaled * R[2, 0] + z_scaled * R[2, 1] + y_scaled * R[2, 2]
            else:  # pidx == 2, zhat is path direction
                # x_scaled -> x, y_scaled -> y, z_scaled -> path
                x_world = base_frame.t[0] + x_scaled * R[0, 0] + y_scaled * R[0, 1] + z_scaled * R[0, 2]
                y_world = base_frame.t[1] + x_scaled * R[1, 0] + y_scaled * R[1, 1] + z_scaled * R[1, 2]
                z_world = base_frame.t[2] + x_scaled * R[2, 0] + y_scaled * R[2, 1] + z_scaled * R[2, 2]
            
            ax.plot_surface(x_world, y_world, z_world, color=surfaceColor, 
                          alpha=surfaceOpacity, edgecolor=edgeColor if edgeColor else None,
                          linewidth=0.5 if edgeColor else 0)
        
        return plotHandles
    
    def addToWidget(self, widget, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                    proximalColor=proximalColorDefault, centerColor=centerColorDefault, distalColor=distalColorDefault,
                    sphereColor=sphereColorDefault, showSphere=False, 
                    surfaceColor=linkColorDefault, 
                    showSurface=True, showAxis=False, axisScale=jointAxisScaleDefault, showPoses=True, poseAxisScaleMultipler=None):
        import pyqtgraph.opengl as gl
        from style import linkColorList
        
        # Call parent's addToWidget for poses and axis  
        super().addToWidget(widget=widget, xColor=xColor, yColor=yColor, zColor=zColor, 
                          proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor, 
                          sphereColor=sphereColor, showSphere=showSphere,
                          surfaceColor=surfaceColor, showSurface=False, showAxis=showAxis,
                          axisScale=axisScale, showPoses=showPoses, poseAxisScaleMultipler=poseAxisScaleMultipler)
        
        if showSurface:
            # Render the tip as a stretched hemisphere, matching addToPlot logic
            n_lat = 6    # latitude divisions
            n_lon = 8    # longitude divisions

            # Hemisphere theta range
            if self.forward:
                theta = np.linspace(0, np.pi/2, n_lat)
            else:
                theta = np.linspace(np.pi/2, np.pi, n_lat)
            phi = np.linspace(0, 2*np.pi, n_lon, endpoint=False)

            # Build vertex grid (n_lat x n_lon) + 1 pole vertex
            vertices = []
            for i in range(n_lat):
                for j in range(n_lon):
                    x_s = np.sin(theta[i]) * np.cos(phi[j])
                    y_s = np.sin(theta[i]) * np.sin(phi[j])
                    z_s = np.cos(theta[i])
                    # Scale
                    x_sc = self.r * x_s
                    y_sc = self.r * y_s
                    z_sc = self.neutralLength * z_s
                    vertices.append([x_sc, y_sc, z_sc])

            # Pole vertex (the closed tip of the hemisphere)
            if self.forward:
                # theta=0 pole: z_sphere=cos(0)=1
                vertices.append([0, 0, self.neutralLength])
            else:
                # theta=pi pole: z_sphere=cos(pi)=-1
                vertices.append([0, 0, -self.neutralLength])

            pole_idx = n_lat * n_lon
            vertices = np.array(vertices, dtype=np.float32)

            # Build triangle faces
            faces = []
            for i in range(n_lat - 1):
                for j in range(n_lon):
                    j_next = (j + 1) % n_lon
                    v00 = i * n_lon + j
                    v01 = i * n_lon + j_next
                    v10 = (i + 1) * n_lon + j
                    v11 = (i + 1) * n_lon + j_next
                    faces.append([v00, v01, v11])
                    faces.append([v00, v11, v10])

            # Fan triangles connecting pole to the first ring (theta[0] row)
            if self.forward:
                for j in range(n_lon):
                    j_next = (j + 1) % n_lon
                    faces.append([pole_idx, j_next, j])
            else:
                last_row = (n_lat - 1) * n_lon
                for j in range(n_lon):
                    j_next = (j + 1) % n_lon
                    faces.append([pole_idx, last_row + j, last_row + j_next])

            faces = np.array(faces, dtype=np.int32)

            # Transform vertices to world coordinates (same logic as addToPlot)
            if self.forward:
                base_frame = self.ProximalFrame()
            else:
                base_frame = self.DistalFrame()
            R = base_frame.R

            x_sc = vertices[:, 0]
            y_sc = vertices[:, 1]
            z_sc = vertices[:, 2]

            if self.pidx == 0:
                xw = base_frame.t[0] + z_sc * R[0, 0] + x_sc * R[0, 1] + y_sc * R[0, 2]
                yw = base_frame.t[1] + z_sc * R[1, 0] + x_sc * R[1, 1] + y_sc * R[1, 2]
                zw = base_frame.t[2] + z_sc * R[2, 0] + x_sc * R[2, 1] + y_sc * R[2, 2]
            elif self.pidx == 1:
                xw = base_frame.t[0] + x_sc * R[0, 0] + z_sc * R[0, 1] + y_sc * R[0, 2]
                yw = base_frame.t[1] + x_sc * R[1, 0] + z_sc * R[1, 1] + y_sc * R[1, 2]
                zw = base_frame.t[2] + x_sc * R[2, 0] + z_sc * R[2, 1] + y_sc * R[2, 2]
            else:  # pidx == 2
                xw = base_frame.t[0] + x_sc * R[0, 0] + y_sc * R[0, 1] + z_sc * R[0, 2]
                yw = base_frame.t[1] + x_sc * R[1, 0] + y_sc * R[1, 1] + z_sc * R[1, 2]
                zw = base_frame.t[2] + x_sc * R[2, 0] + y_sc * R[2, 1] + z_sc * R[2, 2]

            world_verts = np.column_stack([xw, yw, zw]).astype(np.float32)

            meshdata = gl.MeshData(vertexes=world_verts, faces=faces)
            # Use surfaceColor override (e.g. collision highlighting) if provided, otherwise default color list
            tip_color = surfaceColor if surfaceColor != linkColorDefault else linkColorList
            meshitem = gl.GLMeshItem(meshdata=meshdata, color=tuple(tip_color),
                                     shader='shaded', smooth=True)
            meshitem.setGLOptions('translucent')
            meshitem.setObjectName("Joint")
            widget.plot_widget.addItem(meshitem)
        

class StartTip(Tip):
    def __init__(self, r : float, Pose : SE3, length : float, pathIndex : int = 2):
        super().__init__(r, Pose, length, closesForward=False, pathIndex=pathIndex)

class EndTip(Tip):
    def __init__(self, r : float, Pose : SE3, length : float, pathIndex : int = 2):
        super().__init__(r, Pose, length, closesForward=True, pathIndex=pathIndex)
