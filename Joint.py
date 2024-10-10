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

from TubularPattern import *
from geometryHelpers import *
from style import *

class LineItemWithID(gl.GLLinePlotItem):
    def __init__(self, id : int = -1, **kwds):
        """All keyword arguments are passed to setData()"""
        super().__init__(**kwds)
        self.id = id

class LineSphere(gl.GLMeshItem):
    def __init__(self, position = [], rotation=0.0, **kwds):
        super().__init__(**kwds)
        self.position = position
        self.rotation = rotation

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
        self.id = 0
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
    
    def changeRadius(self, r : float):
        self.r = r
    
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

    def generate_extended_line_points(self, point1, point2, gap):
        point1 = np.array(point1)
        point2 = np.array(point2)
        
        direction = point2 - point1
        distance = np.linalg.norm(direction)
        direction = direction / distance

        point1 = point2 - direction * distance * 4
        
        extended_length = 8 * distance
        num_points = int(extended_length / gap) + 1
        start_point = point1 - direction * distance

        points = [start_point + i * gap * direction for i in range(num_points)]
    
        return np.array(points)

    def generate_extended_axis(self, point1, point2, length):
        direction_vector = point2 - point1
    
        midpoint = (point1 + point2) / 2

        p1 = midpoint - direction_vector * length
        p2 = midpoint + direction_vector * length

        return np.array([p1, p2])
    
    def rotation_matrix(self, axis, theta):
        # rodrigues rotation formula
        axis = np.asarray(axis)
        axis = axis / np.linalg.norm(axis)
        a = np.cos(theta / 2.0)
        b, c, d = -axis * np.sin(theta / 2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                        [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                        [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
    
    def generate_circle_points(self, axis, center, rad=1.0, num_points=10, rotation=0.0):
        #angles where the points are placed
        angles = self.generate_angles(num_points)

        axis = np.array(axis)
        axis = axis / np.linalg.norm(axis)
        
        #vector not parallel to the axis
        if (axis == [1, 0, 0]).all() or (axis == [-1, 0, 0]).all():
            not_parallel = np.array([0, 1, 0])
        else:
            not_parallel = np.array([1, 0, 0])

        v1 = np.cross(axis, not_parallel)
        v1 = v1 / np.linalg.norm(v1)
        v1 = v1 * rad

        points = []

        for angle in angles:
            R = self.rotation_matrix(axis, angle + rotation)

            point = center + v1

            line_point = np.array(point.tolist()) - center

            rotated_point = np.dot(R, line_point) + center
            points.append(rotated_point)
        
        return points

    def generate_angles(self, num_points=10):
        angles = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
        return angles
    
    def rotate_vector(self, vector, axis='x'):
        if axis == 'x':
            rotation_matrix = np.array([
                [1, 0, 0],
                [0, 0, -1],
                [0, 1, 0]
            ])
        elif axis == 'y':
            rotation_matrix = np.array([
                [0, 0, 1],
                [0, 1, 0],
                [-1, 0, 0]
            ])

        rotated_vector = np.dot(rotation_matrix, vector)
        return rotated_vector

    def addTranslateArrows(self, widget, selectedArrow=-1, local=True, frame : SE3 = None):
        rad = self.boundingBall().r
        colors = [(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)]

        if selectedArrow != -1:
            colors[selectedArrow] = (1, 1, 1, 1)

        extended_axis_color = [(1, 0, 0, 0.2), (0, 1, 0, 0.2), (0, 0, 1, 0.2)]
        point1 = self.Pose.t

        if local:
            axes = [self.Pose.R[:, i] for i in range(3)]
        else:
            axes = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]

        if frame:
            axes = [frame.R[:, i] for i in range(3)]

        if selectedArrow != -1:
            point2 = point1 + rad * self.r * axes[selectedArrow]
            extended_line_points = self.generate_extended_line_points(point1, point2, 0.1)
            extended_axis = self.generate_extended_axis(point1, point2, 20)
            extended_axis_line = gl.GLLinePlotItem(pos=extended_axis, color=extended_axis_color[selectedArrow], width=5, antialias=True)
            widget.plot_widget.addItem(extended_axis_line)

            for line in extended_line_points:
                md = gl.MeshData.sphere(rows=3, cols=2)
                sphere = LineSphere(meshdata=md, color=lineSphereColor, shader='shaded', smooth=True, position=line)
                sphere.setObjectName("line_sphere")
                sphere.setGLOptions('translucent')
                sphere.scale(2.0, 2.0, 0.1)
                sphere.translate(*line)
                widget.plot_widget.addItem(sphere)

        for i, axis in enumerate(axes):
            pos = np.array([point1, point1 + rad * self.r * axis])
            line = LineItemWithID(pos=pos, color=colors[i], width=10, antialias=True, id=i)
            line.setObjectName("Arrow")
            widget.plot_widget.addItem(line)
    
    def addRotateArrows(self, widget, selectedArrow=-1, local=True, frame : SE3 = None):
        if local:
            axes = [self.Pose.R[:, i] for i in range(3)]
        else: 
            axes = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]

        if frame:
            axes = [frame.R[:, i] for i in range(3)]

        rad = self.boundingBall().r
        colors = rotateArrowColors

        if selectedArrow != -1:
            colors[selectedArrow] = selectedArrowColor

        center = self.Pose.t
        extended_circle_color = extendedCircleColors
        num_points = 20

        R = self.Pose.R
        theta = [np.arctan2(R[2, 1], R[2, 2]), np.arctan2(R[0, 2], R[0, 0]), np.arctan2(R[1, 0], R[1, 1])]

        # generate the marker points around the arrows
        if selectedArrow != -1:
            points = self.generate_circle_points(axes[selectedArrow], center, rad, num_points, theta[selectedArrow])
            angles = self.generate_angles(num_points)

            for i, point in enumerate(points[:num_points]):
                md = gl.MeshData.sphere(rows=3, cols=3)
                sphere = LineSphere(meshdata=md, color=lineSphereColor, shader='shaded', smooth=True, position=point, rotation=angles[i])
                sphere.setObjectName("rotate_sphere")
                sphere.setGLOptions('translucent')
                sphere.scale(0.1, 0.1, 0.1)
                sphere.translate(*point)
                widget.plot_widget.addItem(sphere)

        #swap the axes so that the selected axis gets rendered last -> appears above the other axes
        numbered_axes = [[0, axes[0]], [1, axes[1]], [2, axes[2]]]

        if selectedArrow != -1:
            numbered_axes[selectedArrow], numbered_axes[2] = numbered_axes[2], numbered_axes[selectedArrow]

        for tul in numbered_axes:    
            i = tul[0]
            axis = tul[1]

            points = self.generate_circle_points(axis, center, rad, num_points, theta[selectedArrow])
            points.append(points[0])

            #generate the lighter circles
            circle = LineItemWithID(pos=points, color=extended_circle_color[i], width=5, antialias=True, id=i)
            circle.setObjectName("circle")
            widget.plot_widget.addItem(circle)

            arrow = np.array(points[-num_points//4:])
            arrow2 = np.array(points[:num_points//4])
            arrows = np.append(arrow, arrow2, axis=0)

            # generate the arrows
            line = LineItemWithID(pos=arrows, color=colors[i], width=15, antialias=True, id=i)
            line.setObjectName("Arrow")
            widget.plot_widget.addItem(line)
