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
    
    def ProximalPose(self) -> SE3:
        return SE3.Trans(-(self.neutralLength/2) * self.pathDirection()) @ self.Pose
    
    def NeutralDistalPose(self) -> SE3:
        return SE3.Trans((self.neutralLength/2) * self.pathDirection()) @ self.Pose
    
    def StateTransformationFromNeutral(self) -> SE3:
        return self.stateChangeTransformation(self.state)
    
    def DistalPose(self) -> SE3:
        return self.StateTransformationFromNeutral() @ self.NeutralDistalPose()        
    
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
        return SE3(self.ProximalPose().A[:,self.dubinsColumnOrder()])
    
    def DistalDubinsFrame(self) -> SE3():
        return SE3(self.DistalPose().A[:,self.dubinsColumnOrder()])
    
    def proximalPosition(self) -> np.ndarray:
        return self.ProximalPose().t
    
    def distalPosition(self) -> np.ndarray:
        return self.DistalPose().t
    
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
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, surfaceColor='m',
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
            Poses = np.array([self.ProximalPose(), self.DistalPose(), self.Pose])
            oColors = np.array([proximalColor, distalColor, centerColor])
            plotHandles = addPosesToPlot(Poses, ax, self.r, 
                                         xColor, yColor, zColor, oColors)
        else:
            plotHandles = None
        return plotHandles

    def addToWidget(self, widget, xColor=(1, 0, 0, 1), yColor=(0, 1, 0, 1), zColor=(0, 0, 1, 1), 
                    proximalColor=(0, 1, 1, 1), centerColor=(1, 0, 1, 1), distalColor=(1, 1, 0, 1),
                    sphereColor=(1, 0, 0, 0.3), showSphere=False, surfaceColor=(1, 0, 1, 0.5),
                    surfaceOpacity=0.5, showSurface=True, showAxis=False, 
                    axisScale=10, showPoses=True):
        if showAxis:
            for i, color in enumerate([xColor, yColor, zColor]):
                start_point = self.Pose.t
                end_point = start_point + axisScale * self.Pose.R[:, i]
                points = np.array([start_point, end_point])
                line = gl.GLLinePlotItem(pos=points, color=color, width=2, antialias=True)
                widget.plot_widget.addItem(line)

        if showPoses:
            for pose, color in zip([self.ProximalPose(), self.DistalPose(), self.Pose], [proximalColor, distalColor, centerColor]):
                for i, axis_color in enumerate([xColor, yColor, zColor]):
                    start_point = pose.t
                    end_point = start_point + axisScale * pose.R[:, i]
                    points = np.array([start_point, end_point])
                    line = gl.GLLinePlotItem(pos=points, color=axis_color, width=2, antialias=True)
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
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, surfaceColor='m',
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=True,
             axisScale=10, showPoses=True):
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
            ProximalPose = self.ProximalPose()
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
            
            DistalPose = self.DistalPose()
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
        
    def addToWidget(self, widget, xColor=(1, 0, 0, 1), yColor=(0, 1, 0, 1), zColor=(0, 0, 1, 1),
                    proximalColor=(0, 1, 1, 1), centerColor=(1, 0, 1, 1), distalColor=(1, 1, 0, 1),
                    sphereColor=(1, 0, 0, 0.3), showSphere=False, surfaceColor=(1, 0, 1, 0.5),
                    surfaceOpacity=0.5, showSurface=True, showAxis=True,
                    axisScale=10, showPoses=True):

        super().addToWidget(widget, xColor, yColor, zColor, proximalColor,
                            centerColor, distalColor, sphereColor, showSphere,
                            surfaceColor, surfaceOpacity, showSurface, showAxis,
                            axisScale, showPoses)

        if showSurface:
            scale = self.pattern.baseSideLength / 2
            centerSegment = np.array([self.Pose.t - scale * self.Pose.R[:, 2],
                                    self.Pose.t + scale * self.Pose.R[:, 2]])

            radialCount = self.numSides + 1
            angle = np.linspace(0, 2 * np.pi, radialCount) + np.pi / self.numSides
            u = self.r * np.cos(angle)
            v = self.r * np.sin(angle)

            for pose in [self.ProximalPose(), self.DistalPose()]:
                uhat = pose.R[:, 1]
                vhat = pose.R[:, 2]
                basePoints = np.array([pose.t + u[i] * uhat + v[i] * vhat for i in range(radialCount - 1)])
                allPoints = np.vstack([basePoints, centerSegment])

                hull = ConvexHull(allPoints)
                vertices = allPoints[hull.vertices]
                faces = hull.simplices

                color_list = (surfaceColor[0], surfaceColor[1], surfaceColor[2], surfaceOpacity)

                meshdata = gl.MeshData(vertexes=vertices, faces=faces)
                item = gl.GLMeshItem(meshdata=meshdata, color=tuple(color_list), shader='shaded', smooth=False, drawEdges=True)
                item.setGLOptions('translucent')
                widget.plot_widget.addItem(item)

        if showSphere:
            color_list = (sphereColor[0], sphereColor[1], sphereColor[2], sphereColor[3])
            self.boundingBall().addToWidget(widget, color_list)

        
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
        return Cylinder(self.r, self.ProximalPose().t, self.pathDirection(), self.length())
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, surfaceColor='m',
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=True, 
             axisScale=10, showPoses=True):
        plotHandles = super().addToPlot(ax, xColor, yColor, zColor, proximalColor,
                          centerColor, distalColor, sphereColor, showSphere,
                          surfaceColor, surfaceOpacity, showSurface, showAxis,
                          axisScale, showPoses)
        if showSurface:
            self.boundingCylinder().addToPlot(ax, color=surfaceColor, alpha=surfaceOpacity)
            
        return plotHandles
    
    def addToWidget(self, widget, xColor=(1, 0, 0, 1), yColor=(0, 1, 0, 1), zColor=(0, 0, 1, 1),
                    proximalColor=(0, 1, 1, 1), centerColor=(1, 0, 1, 1), distalColor=(1, 1, 0, 1),
                    sphereColor=(1, 0, 0, 0.3), showSphere=False, surfaceColor=(1, 0, 1, 0.5),
                    surfaceOpacity=0.5, showSurface=True, showAxis=True,
                    axisScale=10, showPoses=True):
        
        super().addToWidget(widget, xColor, yColor, zColor, proximalColor,
                            centerColor, distalColor, sphereColor, showSphere,
                            surfaceColor, surfaceOpacity, showSurface, showAxis,
                            axisScale, showPoses)
        if showSurface:
            color_list = (surfaceColor[0], surfaceColor[1], surfaceColor[2], surfaceOpacity)
            self.boundingCylinder().addToWidget(widget, numPointsPerCircle=32, numCircles=10, color_list=tuple(color_list))

        if showSphere:
            color_list = (sphereColor[0], sphereColor[1], sphereColor[2], sphereColor[3])
            self.boundingBall().addToWidget(widget, color_list)
                        
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
    
    def addToWidget(self, widget, xColor=(1, 0, 0, 1), yColor=(0, 1, 0, 1), zColor=(0, 0, 1, 1),
                    proximalColor=(0, 1, 1, 1), centerColor=(1, 0, 1, 1), distalColor=(1, 1, 0, 1),
                    sphereColor=(1, 0, 0, 0.3), showSphere=False, surfaceColor=(1, 0, 1, 0.5),
                    surfaceOpacity=0.5, showSurface=True, showAxis=True,
                    axisScale=10, showPoses=True):
        
        if showAxis:
            for i, color in enumerate([xColor, yColor, zColor]):
                start_point = self.Pose.t
                end_point = start_point + axisScale * self.Pose.R[:, i]
                points = np.array([start_point, end_point])
                line = gl.GLLinePlotItem(pos=points, color=color, width=2, antialias=True)
                widget.plot_widget.addItem(line)

        if showPoses:
            for i, axis_color in enumerate([xColor, yColor, zColor]):
                start_point = self.Pose.t
                end_point = start_point + axisScale * self.Pose.R[:, i]
                points = np.array([start_point, end_point])
                line = gl.GLLinePlotItem(pos=points, color=axis_color, width=2, antialias=True)
                widget.plot_widget.addItem(line)

        if showSphere:
            color_list = (sphereColor[0], sphereColor[1], sphereColor[2], sphereColor[3])
            self.boundingBall().addToWidget(widget, color_list)

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
             axisScale=10, showPoses=True):
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
                DistalPose = self.DistalPose()
                TipSegment = np.array([DistalPose.t - scale * DistalPose.R[:,tipSegmentIndex],
                                          DistalPose.t + scale * DistalPose.R[:,tipSegmentIndex]])
                
                ProximalPose = self.ProximalPose()
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
                ProximalPose = self.ProximalPose()
                TipSegment = np.array([ProximalPose.t - scale * ProximalPose.R[:,tipSegmentIndex],
                                          ProximalPose.t + scale * ProximalPose.R[:,tipSegmentIndex]])
                
                DistalPose = self.DistalPose()
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
    
    def addToWidget(self, widget, xColor=(1, 0, 0, 1), yColor=(0, 1, 0, 1), zColor=(0, 0, 1, 1),
                    proximalColor=(0, 1, 1, 1), centerColor=(1, 0, 1, 1), distalColor=(1, 1, 0, 1),
                    sphereColor=(1, 0, 0, 0.3), showSphere=False, surfaceColor=(1, 0, 1, 0.5),
                    surfaceOpacity=0.5, showSurface=True, showAxis=True,
                    axisScale=10, showPoses=True):
        
        super().addToWidget(widget, xColor, yColor, zColor, proximalColor,
                            centerColor, distalColor, sphereColor, showSphere,
                            surfaceColor, surfaceOpacity, showSurface, showAxis,
                            axisScale, showPoses)

        if showSurface:
            color_list = (surfaceColor[0], surfaceColor[1], surfaceColor[2], surfaceOpacity)
            radialCount = self.numSides + 1
            angle = np.linspace(0, 2 * np.pi, radialCount) + np.pi / self.numSides
            u = self.r * np.cos(angle)
            v = self.r * np.sin(angle)
            scale = self.pattern.baseSideLength / 2

            proximalPose = self.ProximalPose()
            distalPose = self.DistalPose()
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
                item = gl.GLMeshItem(meshdata=meshdata, color=color_list, smooth=False, drawEdges=True, shader='shaded', glOptions='translucent')
                widget.plot_widget.addItem(item)
        
        if showSphere:
            color_list = (sphereColor[0], sphereColor[1], sphereColor[2], sphereColor[3])
            self.boundingBall().addToWidget(widget, color_list)

class StartTip(Tip):
    def __init__(self, numSides : int, r : float, Pose : SE3, length : float):
        super().__init__(numSides, r, Pose, length, closesForward=False)

class EndTip(Tip):
    def __init__(self, numSides : int, r : float, Pose : SE3, length : float):
        super().__init__(numSides, r, Pose, length, closesForward=True)