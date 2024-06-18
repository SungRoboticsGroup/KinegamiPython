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
import pyqtgraph.opengl as gl
from meshHelpers import *

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
        self.id = 0
    
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
    
    def ProximalDubinsFrame(self) -> SE3():
        return SE3(self.ProximalFrame().A[:,self.dubinsColumnOrder()])
    
    def DistalDubinsFrame(self) -> SE3():
        return SE3(self.DistalFrame().A[:,self.dubinsColumnOrder()])
    
    def proximalPosition(self) -> np.ndarray:
        return self.ProximalFrame().t
    
    def distalPosition(self) -> np.ndarray:
        return self.DistalFrame().t
    
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

    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor=proximalColorDefault, centerColor=centerColorDefault, distalColor=distalColorDefault,
             sphereColor=sphereColorDefault, showSphere=False, surfaceColor=jointColorDefault, 
             showSurface=True, showAxis=False, axisScale=10, showPoses=True):
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
            line_item = gl.GLLinePlotItem(pos=jointAxis, color=(0.75, 0.75, 0.75, 1), width=2, antialias=True)  # Using a silver color
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
    
    def generate_circle_points(self, axis, center, r=1.0, num_points=10, rotation=0.0):
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
        v2 = np.cross(axis, v1)
        v2 = v2 / np.linalg.norm(v2)

        line1_point1 = center + v1
        line1_point2 = center - v1
        line2_point1 = center + v2
        line2_point2 = center - v2
        
        line1 = (line1_point1.tolist(), line1_point2.tolist())
        line2 = (line2_point1.tolist(), line2_point2.tolist())

        points = []

        for angle in angles:
            R = self.rotation_matrix(axis, angle + rotation)

            point = center + v1

            line_point = np.array(point.tolist()) - center

            rotated_point = np.dot(R, line_point) + center
            points.append(rotated_point)

        line1 = (line1_point1.tolist(), line1_point2.tolist())
        line2 = (line2_point1.tolist(), line2_point2.tolist())

        lines = [line1, line2]

        line_results = []

        for line in lines:
            # Center the line points at the origin
            line_point1 = np.array(line[0]) - center
            line_point2 = np.array(line[1]) - center
            
            # Rotate the points
            rotated_point1 = np.dot(R, line_point1) + center
            rotated_point2 = np.dot(R, line_point2) + center
            line_results.append([rotated_point1, rotated_point2])
        
        return points, line_results

    def generate_angles(self, num_points=10):
        angles = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
        return angles

    def addTranslateArrows(self, widget, selectedArrow=-1):
        rad = self.boundingBall().r
        colors = [(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)]

        if selectedArrow != -1:
            colors[selectedArrow] = (1, 1, 1, 1)

        axes = [self.Pose.R[:, i] for i in range(3)]
        extended_axis_color = [(1, 0, 0, 0.2), (0, 1, 0, 0.2), (0, 0, 1, 0.2)]
        point1 = self.Pose.t

        if selectedArrow != -1:
            point2 = point1 + rad * self.r * axes[selectedArrow]
            extended_line_points = self.generate_extended_line_points(point1, point2, 0.1)
            extended_axis = self.generate_extended_axis(point1, point2, 20)
            extended_axis_line = gl.GLLinePlotItem(pos=extended_axis, color=extended_axis_color[selectedArrow], width=5, antialias=True)
            widget.plot_widget.addItem(extended_axis_line)

            for line in extended_line_points:
                md = gl.MeshData.sphere(rows=2, cols=2)
                sphere = LineSphere(meshdata=md, color=[0, 0, 0, 0], shader='shaded', smooth=True, position=line)
                sphere.setObjectName("line_sphere")
                sphere.setGLOptions('translucent')
                sphere.scale(0.02, 0.02, 0.02)
                sphere.translate(*line)
                widget.plot_widget.addItem(sphere)

        for i, axis in enumerate(axes):
            pos = np.array([point1, point1 + rad * self.r * axis])
            line = LineItemWithID(pos=pos, color=colors[i], width=10, antialias=True, id=i)
            line.setObjectName("Arrow")
            widget.plot_widget.addItem(line)
    
    def addRotateArrows(self, widget, selectedArrow=-1):
        rad = self.boundingBall().r
        colors = [(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)]

        if selectedArrow != -1:
            colors[selectedArrow] = (1, 1, 1, 1)

        axes = [self.Pose.R[:, i] for i in range(3)]
        center = self.Pose.t
        extended_circle_color = [(1, 0, 0, 0.5), (0, 1, 0, 0.5), (0, 0, 1, 0.5)]
        num_points = 20

        R = self.Pose.R
        theta = [np.arctan2(R[2, 1], R[2, 2]), np.arctan2(R[0, 2], R[0, 0]), np.arctan2(R[1, 0], R[1, 1])]

        # generate the marker points around the arrows
        if selectedArrow != -1:
            points, lines = self.generate_circle_points(axes[selectedArrow], center, rad, num_points, theta[selectedArrow])
            angles = self.generate_angles(num_points)

            for i, point in enumerate(points[:num_points]):
                md = gl.MeshData.sphere(rows=3, cols=3)
                sphere = LineSphere(meshdata=md, color=[0, 0, 0, 0], shader='shaded', smooth=True, position=point, rotation=angles[i])
                sphere.setObjectName("rotate_sphere")
                sphere.setGLOptions('translucent')
                sphere.scale(0.1, 0.1, 0.1)
                sphere.translate(*point)
                widget.plot_widget.addItem(sphere)
            
            #marker lines for debugging
            for line in lines:
                l = LineItemWithID(pos=line, color=(0,0,0,1), width=10, antialias=True)
                widget.plot_widget.addItem(l)

        for i, axis in enumerate(axes):    
            points, lines = self.generate_circle_points(axis, center, rad, num_points, theta[selectedArrow])
            points.append(points[0])

            #generate the lighter circles
            circle = LineItemWithID(pos=points, color=extended_circle_color[i], width=5, antialias=True, id=i)
            circle.setObjectName("circle")
            widget.plot_widget.addItem(circle)

            # generate the arrows
            line = LineItemWithID(pos=points[:len(points)//2 + 1], color=colors[i], width=10, antialias=True, id=i)
            line.setObjectName("Arrow")
            widget.plot_widget.addItem(line)

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
    
class OrigamiJoint(Joint):
    def __init__(self, numSides : int, r : float, neutralLength : float, Pose : SE3(), 
                 initialState : float = 0, id : int = 0):
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
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, surfaceColor='m',
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=True,
             axisScale=jointAxisScaleDefault, showPoses=True):
        plotHandles = super().addToPlot(ax, xColor, yColor, zColor, proximalColor,
                          centerColor, distalColor, sphereColor, showSphere,
                          surfaceColor, surfaceOpacity, showSurface, showAxis,
                          axisScale, showPoses)
        if showSurface:
            scale = self.pattern.baseSideLength / 2
            CenterSegment = np.array([self.Pose.t - scale * self.Pose.R[:,2],
                                      self.Pose.t + scale * self.Pose.R[:,2]])
            #https://stackoverflow.com/questions/63207496/how-to-visualize-polyhedrons-defined-by-their-vertices-in-3d-with-matplotlib-or
            
            radialCount = self.numSides + 1
            angle = np.linspace(0, 2*np.pi, radialCount) + np.pi/self.numSides
            u = self.r * np.cos(angle)
            v = self.r * np.sin(angle)
            ProximalPose = self.ProximalFrame()
            uhatProximal = ProximalPose.R[:,1]
            vhatProximal = ProximalPose.R[:,2]
            ProximalBase = ProximalPose.t + u.reshape(-1,1) @ uhatProximal.reshape(1,3) + v.reshape(-1,1) @ vhatProximal.reshape(1,3)
            
            ProximalPoints = np.vstack((ProximalBase, CenterSegment))
            ProximalHull = ConvexHull(ProximalPoints)
            for s in ProximalHull.simplices:
                tri = Poly3DCollection([ProximalPoints[s]])
                tri.set_color(surfaceColor)
                tri.set_alpha(surfaceOpacity)
                ax.add_collection3d(tri)
            
            DistalPose = self.DistalFrame()
            uhatDistal = DistalPose.R[:,1]
            vhatDistal = DistalPose.R[:,2]
            DistalBase = DistalPose.t + u.reshape(-1,1) @ uhatDistal.reshape(1,3) + v.reshape(-1,1) @ vhatDistal.reshape(1,3)
            DistalPoints = np.vstack((DistalBase, CenterSegment))
            DistalHull = ConvexHull(DistalPoints)
            for s in DistalHull.simplices:
                tri = Poly3DCollection([DistalPoints[s]])
                tri.set_color(surfaceColor)
                tri.set_alpha(surfaceOpacity)
                ax.add_collection3d(tri)
            
        return plotHandles
        
    def addToWidget(self, widget, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                    proximalColor=proximalColorDefault, centerColor=centerColorDefault, distalColor=distalColorDefault,
                    sphereColor=sphereColorDefault, showSphere=False, surfaceColor=jointColorDefault, 
                    showSurface=True, showAxis=True, axisScale=jointAxisScaleDefault, showPoses=True, poseAxisScaleMultipler=None):
        if showSurface:
            scale = self.pattern.baseSideLength / 2
            centerSegment = np.array([self.Pose.t - scale * self.Pose.R[:, 2],
                                    self.Pose.t + scale * self.Pose.R[:, 2]])

            radialCount = self.numSides + 1
            angle = np.linspace(0, 2 * np.pi, radialCount) + np.pi / self.numSides
            u = self.r * np.cos(angle)
            v = self.r * np.sin(angle)

            pose = self.ProximalFrame()
            uhat = pose.R[:, 1]
            vhat = pose.R[:, 2]
            basePoints = np.array([pose.t + u[i] * uhat + v[i] * vhat for i in range(radialCount - 1)])
            allPoints = np.vstack([basePoints, centerSegment])

            hull = ConvexHull(allPoints)
            vertices = allPoints[hull.vertices]
            faces = hull.simplices

            pose2 = self.DistalFrame()
            uhat2 = pose2.R[:, 1]
            vhat2 = pose2.R[:, 2]
            basePoints2 = np.array([pose2.t + u[i] * uhat2 + v[i] * vhat2 for i in range(radialCount - 1)])
            allPoints2 = np.vstack([basePoints2, centerSegment])

            hull2 = ConvexHull(allPoints2)
            vertices = np.append(vertices, allPoints2[hull2.vertices], axis=0)
            faces = np.append(faces, hull2.simplices + 6, axis=0)

            meshdata = gl.MeshData(vertexes=vertices, faces=faces)
            item = MeshItemWithID(meshdata=meshdata, color=surfaceColor, shader='shaded', smooth=False, drawEdges=True, id=self.id)
            item.setGLOptions('translucent')
            item.setObjectName("Joint")
            widget.plot_widget.addItem(item)

        super().addToWidget(widget, xColor, yColor, zColor, proximalColor,
                            centerColor, distalColor, sphereColor, showSphere, 
                            surfaceColor, showSurface, showAxis,
                            axisScale, showPoses, poseAxisScaleMultipler)

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
    
    def boundingCylinder(self) -> Cylinder:
        uhat = (self.Pose @ SE3.Rz(np.pi/self.numSides)).R[:,1]
        return Cylinder(self.r, self.ProximalFrame().t, self.pathDirection(), 
                        self.length(), uhat)
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=jointColorDefault, edgeColor=jointColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=True, 
             axisScale=10, showPoses=True):
        plotHandles = super().addToPlot(ax, xColor, yColor, zColor, proximalColor,
                          centerColor, distalColor, sphereColor, showSphere,
                          surfaceColor, surfaceOpacity, showSurface, showAxis,
                          axisScale, showPoses)
        if showSurface:
            self.boundingCylinder().addToPlot(ax, color=surfaceColor, 
                                              alpha=surfaceOpacity, 
                                              edgeColor=edgeColor,
                                              numPointsPerCircle=self.numSides)
            
        return plotHandles
    
    def addToWidget(self, widget, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                    proximalColor=proximalColorDefault, centerColor=centerColorDefault, distalColor=distalColorDefault,
                    sphereColor=sphereColorDefault, showSphere=False, surfaceColor=jointColorDefault,
                    showSurface=True, showAxis=True, axisScale=jointAxisScaleDefault, showPoses=True, poseAxisScaleMultipler=None):
        
        if showSurface:
            self.boundingCylinder().addToWidget(widget, numPointsPerCircle=self.numSides, numCircles=2, color_list=surfaceColor, 
                                                is_joint=True, id=self.id)

        super().addToWidget(widget, xColor, yColor, zColor, proximalColor,
                            centerColor, distalColor, sphereColor, showSphere,
                            surfaceColor, showSurface, showAxis,
                            axisScale, showPoses, poseAxisScaleMultipler)
                        
class Waypoint(OrigamiJoint):
    # path direction through a waypoint defaults to zhat
    def __init__(self, numSides : int, r : float, Pose : SE3, pathIndex : int = 2):
        assert(pathIndex in [0,1,2])
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
    
    def _setView(self, v):
        self.__view = v
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, surfaceColor='m',
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
            Poses = np.array([self.Pose])
            oColors = np.array([centerColor])
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
            jointAxis = np.array([self.Pose.t - axisScale * self.r * zhat,
                                self.Pose.t + axisScale * self.r * zhat])
            line_item = gl.GLLinePlotItem(pos=jointAxis, color=(0.75, 0.75, 0.75, 1), width=2, antialias=True)  # Using a silver color
            widget.plot_widget.addItem(line_item)

        if showPoses:
            for i, axis_color in enumerate([xColor, yColor, zColor]):
                poseAxisScale = self.r
                if poseAxisScaleMultipler:
                    poseAxisScale *= poseAxisScaleMultipler
                start_point = self.Pose.t
                end_point = start_point + poseAxisScale * self.Pose.R[:, i]
                points = np.array([start_point, end_point])
                line = gl.GLLinePlotItem(pos=points, color=axis_color, width=2, antialias=True)
                widget.plot_widget.addItem(line)
        
        if showSphere:
            self.boundingBall().addToWidget(widget, sphereColor, is_waypoint=True, id=self.id)
        else:
            self.boundingBall().addToWidget(widget, (0,0,0,0.5), is_waypoint=True, id=self.id)

class Tip(OrigamiJoint):
    def __init__(self, numSides : int, r : float, Pose : SE3, length : float, 
                 closesForward : bool = True):
        super().__init__(numSides, r, length, Pose)
        self.pattern = TipPattern(numSides, r, length, closesForward)
        self.forward = closesForward
    
    def pathIndex(self) -> int:
        return 2 # zhat
    
    def stateChangeTransformation(self, stateChange : float) -> SE3:
        return SE3()
    
    def stateRange(self) -> list:
        return [0,0]
    
    def boundingRadius(self) -> float:
        return norm([self.r, self.neutralLength/2])
    
    def boundingBall(self) -> Ball:
        return Ball(self.Pose.t, self.boundingRadius())
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, surfaceColor='m',
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=False,
             axisScale=jointAxisScaleDefault, showPoses=True):
        plotHandles = super().addToPlot(ax, xColor, yColor, zColor, proximalColor,
                          centerColor, distalColor, sphereColor, showSphere, 
                          surfaceColor, surfaceOpacity, showSurface, showAxis,
                          axisScale, showPoses)
        if showSurface:
            #https://stackoverflow.com/questions/63207496/how-to-visualize-polyhedrons-defined-by-their-vertices-in-3d-with-matplotlib-or
            radialCount = self.numSides + 1
            angle = np.linspace(0, 2*np.pi, radialCount) + np.pi/self.numSides
            u = self.r * np.cos(angle)
            v = self.r * np.sin(angle)
            scale = self.pattern.baseSideLength / 2
            
            tipSegmentIndex, uhatIndex, vhatIndex = 0,0,1
            
            if self.forward:
                DistalPose = self.DistalFrame()
                TipSegment = np.array([DistalPose.t - scale * DistalPose.R[:,tipSegmentIndex],
                                          DistalPose.t + scale * DistalPose.R[:,tipSegmentIndex]])
                
                ProximalPose = self.ProximalFrame()
                uhatProximal = ProximalPose.R[:,uhatIndex]
                vhatProximal = ProximalPose.R[:,vhatIndex]
                ProximalBase = ProximalPose.t + u.reshape(-1,1) @ uhatProximal.reshape(1,3) + v.reshape(-1,1) @ vhatProximal.reshape(1,3)
                ProximalPoints = np.vstack((ProximalBase, TipSegment))
                ProximalHull = ConvexHull(ProximalPoints)
                for s in ProximalHull.simplices:
                    tri = Poly3DCollection([ProximalPoints[s]])
                    tri.set_color(surfaceColor)
                    tri.set_alpha(surfaceOpacity)
                    ax.add_collection3d(tri)
            else:
                ProximalPose = self.ProximalFrame()
                TipSegment = np.array([ProximalPose.t - scale * ProximalPose.R[:,tipSegmentIndex],
                                          ProximalPose.t + scale * ProximalPose.R[:,tipSegmentIndex]])
                
                DistalPose = self.DistalFrame()
                uhatDistal = DistalPose.R[:,uhatIndex]
                vhatDistal = DistalPose.R[:,vhatIndex]
                DistalBase = DistalPose.t + u.reshape(-1,1) @ uhatDistal.reshape(1,3) + v.reshape(-1,1) @ vhatDistal.reshape(1,3)
                DistalPoints = np.vstack((DistalBase, TipSegment))
                DistalHull = ConvexHull(DistalPoints)
                for s in DistalHull.simplices:
                    tri = Poly3DCollection([DistalPoints[s]])
                    tri.set_color(surfaceColor)
                    tri.set_alpha(surfaceOpacity)
                    ax.add_collection3d(tri)
            
        return plotHandles
    
    def addToWidget(self, widget, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                    proximalColor=proximalColorDefault, centerColor=centerColorDefault, distalColor=distalColorDefault,
                    sphereColor=sphereColorDefault, showSphere=False, surfaceColor=jointColorDefault,
                    showSurface=True, showAxis=False, axisScale=jointAxisScaleDefault, showPoses=True, poseAxisScaleMultipler=None):
        
        if showSurface:
            radialCount = self.numSides + 1
            angle = np.linspace(0, 2 * np.pi, radialCount) + np.pi / self.numSides
            u = self.r * np.cos(angle)
            v = self.r * np.sin(angle)
            scale = self.pattern.baseSideLength / 2

            proximalPose = self.ProximalFrame()
            distalPose = self.DistalFrame()
            proximalBase = proximalPose.t + np.outer(u, proximalPose.R[:, 0]) + np.outer(v, proximalPose.R[:, 1])
            distalBase = distalPose.t + np.outer(u, distalPose.R[:, 0]) + np.outer(v, distalPose.R[:, 1])

            if self.forward:
                tipSegment = np.array([distalPose.t - scale * distalPose.R[:, 0],
                                    distalPose.t + scale * distalPose.R[:, 0]])
                tipPoints = np.vstack((proximalBase, tipSegment))
            else:
                tipSegment = np.array([proximalPose.t - scale * proximalPose.R[:, 0],
                                    proximalPose.t + scale * proximalPose.R[:, 0]])
                tipPoints = np.vstack((distalBase, tipSegment))

            hull = ConvexHull(tipPoints)
            for s in hull.simplices:
                vertices = tipPoints[s]
                meshdata = gl.MeshData(vertexes=vertices, faces=[np.arange(len(vertices))])
                item = MeshItemWithID(meshdata=meshdata, color=surfaceColor, smooth=False, drawEdges=True, shader='shaded', glOptions='translucent', id=self.id)
                widget.plot_widget.addItem(item)
        
        super().addToWidget(widget, xColor, yColor, zColor, proximalColor,
                            centerColor, distalColor, sphereColor, showSphere,
                            surfaceColor, showSurface, showAxis,
                            axisScale, showPoses, poseAxisScaleMultipler)

class StartTip(Tip):
    def __init__(self, numSides : int, r : float, Pose : SE3, length : float):
        super().__init__(numSides, r, Pose, length, closesForward=False)

class EndTip(Tip):
    def __init__(self, numSides : int, r : float, Pose : SE3, length : float):
        super().__init__(numSides, r, Pose, length, closesForward=True)