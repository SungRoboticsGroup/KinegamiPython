# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 21:54:54 2023

@author: dfesh
"""
from spatialmath import SE3, SO3
from abc import ABC, abstractmethod
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
import pyqtgraph.opengl as gl
from pyqtgraph import Transform3D
from PyQt5.QtGui import QMatrix4x4, QVector4D, QVector3D
from OpenGL.GL import glDisable, glEnable, GL_DEPTH_TEST
import math

from TubularPattern import *
from geometryHelpers import *
from style import *

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
        # self.id = 0
        self.collisionCapsules = self.getCapsules()
    
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
    def cloneWithNewRadius(self, new_r: float):
        pass
    
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
    
    def pathDirectionLocal(self) -> np.ndarray:
        direction = np.array([0,0,0])
        direction[self.pathIndex()] = 1
        return direction

    def reverseZhat(self):
        self.Pose = self.Pose @ SE3.Rx(np.pi)
    
    def reversePathDirection(self):
        if self.pathIndex() == 2:
            self.Pose = self.Pose @ SE3.Rx(np.pi)
        else:
            self.Pose = self.Pose @ SE3.Rz(np.pi)

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

    def transformPoseBy(self, Transformation: SE3):
        self.Pose = Transformation @ self.Pose
    
    def applyTransformationToPose(self, Transformation : SE3):
        self.Pose = self.Pose @ Transformation

    def translateAlongZ(self, zChange : float):
        self.Pose = self.Pose @ SE3.Trans([0,0,zChange])

    def placeInFrontOf(self, otherJoint : 'Joint', distance : float):
        translation = SE3.Trans((distance + self.neutralLength/2)*otherJoint.pathDirectionLocal())
        position = otherJoint.distalPosition() + (distance + self.neutralLength/2)*otherJoint.pathDirection()
        rotation = SO3()
        if self.pathIndex() == 0 and otherJoint.pathIndex() == 2:
            rotation = SO3.Ry(np.pi/2)
        elif self.pathIndex() == 2 and otherJoint.pathIndex() == 0:
            rotation = SO3.Ry(-np.pi/2)
        self.Pose = SE3.Rt(rotation @ SO3(otherJoint.Pose.R), position)
    
    def translateAlongX(self, xChange : float):
        self.Pose = self.Pose @ SE3.Trans([xChange,0,0])

    def translateAlongY(self, yChange : float):
        self.Pose = self.Pose @ SE3.Trans([0,yChange,0])

    def rotateAboutZ(self, angleToRotateAboutZ):
        self.applyTransformationToPose(SE3.Rz(angleToRotateAboutZ))
    
    def getCapsules(self):
        return []

    def recomputeCollisionCapsules(self):
        self.collisionCapsules = self.getCapsules()

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

    def generate_extended_axis(self, point1, point2, length):
        direction_vector = point2 - point1
    
        midpoint = (point1 + point2) / 2

        p1 = midpoint - direction_vector * length
        p2 = midpoint + direction_vector * length

        return np.array([p1, p2])
    
    def get_transform3D(self, translate=False):
        m = self.Pose.R
        t = self.Pose.t

        trans = [0, 0, 0]
        if (translate):
            trans = t

        transform = Transform3D()
        transform.setRow(0, QVector4D(m[0][0], m[0][1], m[0][2], trans[0]))
        transform.setRow(1, QVector4D(m[1][0], m[1][1], m[1][2], trans[1]))
        transform.setRow(2, QVector4D(m[2][0], m[2][1], m[2][2], trans[2]))
        transform.setRow(3, QVector4D(0, 0, 0, 1))

        return transform
    
    def create_torus_mesh(self, radius, tube_radius, radial_segments, tubular_segments):
        theta = np.linspace(0, 2 * np.pi, radial_segments)
        phi = np.linspace(0, 2 * np.pi, tubular_segments)
        theta, phi = np.meshgrid(theta, phi)
        theta, phi = theta.flatten(), phi.flatten()

        x = (radius + tube_radius * np.cos(phi)) * np.cos(theta)
        y = (radius + tube_radius * np.cos(phi)) * np.sin(theta)
        z = tube_radius * np.sin(phi)

        vertices = np.vstack([x, y, z]).T
        faces = []

        for i in range(radial_segments):
            for j in range(tubular_segments):
                next_i = (i + 1) % radial_segments
                next_j = (j + 1) % tubular_segments

                faces.append([i * tubular_segments + j,
                            next_i * tubular_segments + j,
                            i * tubular_segments + next_j])
                faces.append([next_i * tubular_segments + j,
                            next_i * tubular_segments + next_j,
                            i * tubular_segments + next_j])

        faces = np.array(faces)
        return gl.MeshData(vertexes=vertices, faces=faces)
    
    def addArrows(self, widget, selectedArrow=-1, local=True, frame: SE3 = None, mode=""):
            desired_axis_pixels = 80
            rad = widget.plot_widget.world_length_for_pixel_length(desired_axis_pixels)
            colors = rotateArrowColors
            center = self.Pose.t
            extended_axis_color = [(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)]

            if local:
                axes = [self.Pose.R[:, i] for i in range(3)]
            else:
                axes = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]
            if frame:
                axes = [frame.R[:, i] for i in range(3)]

            if mode == "Translate":
                for i, axis_vec in enumerate(axes):
                    col = selectedArrowColor if i == selectedArrow else colors[i]
                    start = center
                    end = center + (rad + 1.0) * axis_vec
                    line_item = OverlayLine(pos=np.array([start, end]), color=col, width=8, antialias=True)
                    widget.plot_widget.addItem(line_item)

                if selectedArrow != -1:
                    dir_pt = center + rad * axes[selectedArrow]
                    ext = self.generate_extended_axis(center, dir_pt, 1000)
                    ext_line = gl.GLLinePlotItem(pos=ext,
                                                color=extended_axis_color[selectedArrow],
                                                width=8,
                                                antialias=True)
                    widget.plot_widget.addItem(ext_line)

            elif mode == "Rotate":
                for i, axis in enumerate(axes):
                    helper = np.array([1,0,0])
                    if abs(np.dot(axis, helper)) > 0.9:
                        helper = np.array([0,1,0])
                    u = np.cross(axis, helper); u /= np.linalg.norm(u)
                    v = np.cross(axis, u)

                    pts = np.array([
                        center + (rad + 0.2)*(u*np.cos(t) + v*np.sin(t))
                        for t in np.linspace(0, 2*math.pi, 64)
                    ])

                    col = selectedArrowColor if i == selectedArrow else colors[i]
                    width = 8
                    circle = OverlayLine(pos=pts, color=col, width=width, antialias=True)
                    widget.plot_widget.addItem(circle)

    def addTranslateArrows(self, widget, selectedArrow=-1, local=True, frame : SE3 = None):
        self.addArrows(widget, selectedArrow, local, frame, mode="Translate")
    
    def addRotateArrows(self, widget, selectedArrow=-1, local=True, frame : SE3 = None):
        self.addArrows(widget, selectedArrow, local, frame, mode="Rotate")

class OverlayLine(gl.GLLinePlotItem):
    def paint(self):
        glDisable(GL_DEPTH_TEST)
        super().paint()
        glEnable(GL_DEPTH_TEST)