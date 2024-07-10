from spatialmath import SE3
from abc import ABC, abstractmethod
from geometryHelpers import *
import matplotlib.pyplot as plt
from TubularPattern import *
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
import Joint
from Joint import *
from CollisionDetection import *

class OrigamiJoint(Joint):
    def __init__(self, numSides : int, r : float, neutralLength : float, Pose : SE3(), 
                 initialState : float = 0):
        self.numSides = numSides
        self.polygonInnerAngle = np.pi * (numSides-2)/(2*numSides)
        super().__init__(r, neutralLength, Pose, initialState)
    
    def toPrinted(self, screwRadius):
        pass

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
    
    def toPrinted(self, screwRadius):
        from PrintedJoint import PrintedOrthogonalRevoluteJoint
        return PrintedOrthogonalRevoluteJoint(self.r, -self.totalBendingAngle/2, self.totalBendingAngle/2, self.Pose, screwRadius, initialState=self.initialState)

    def getCapsules(self):
        return [CollisionCapsule(base=self.ProximalDubinsFrame(), radius=self.r, height=self.neutralLength/2),
                CollisionCapsule(base=self.DistalDubinsFrame(), radius=self.r, height = -self.neutralLength/2)]
    
class ExtendedRevoluteJoint(OrigamiJoint):
    def __init__(self, numSides : int, r : float, totalBendingAngle : float, 
                 tubeLength: float, Pose : SE3, numSinkLayers : int = 1,
                 initialState : float = 0):
        polygonInnerAngle = np.pi * (numSides-2)/(2*numSides)
        self.revoluteLength = 2*r*np.sin(polygonInnerAngle)*np.tan(totalBendingAngle/4) #2*delta from paper
        self.tubeLength = tubeLength
        neutralLength = self.revoluteLength + 2*tubeLength
        self.totalBendingAngle = totalBendingAngle
        self.numSinkLayers = numSinkLayers
        super().__init__(numSides, r, neutralLength, Pose, initialState)
        revolutePattern = RevoluteJointPattern(self.numSides, self.r, 
                                            totalBendingAngle, numSinkLayers)
        self.pattern = TubeFittingPattern(numSides, r, tubeLength).append(revolutePattern).append(TubeFittingPattern(numSides, r, tubeLength))
    
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
    
    def RevoluteProximalFrame(self) -> SE3:
        PF = self.ProximalFrame()
        return SE3.Trans(self.tubeLength*PF.R[:,0]) @ PF
    
    def RevoluteDistalFrame(self) -> SE3:
        DF = self.DistalFrame()
        return SE3.Trans(-self.tubeLength*DF.R[:,0]) @ DF

    def proximalCylinder(self) -> Cylinder:
        PF = self.DistalFrame()
        uhat = (PF @ SE3.Rx(np.pi/self.numSides)).R[:,1]
        return Cylinder(self.r, self.ProximalFrame().t, PF.R[:,0], 
                        self.tubeLength, uhat)
    
    def distalCylinder(self) -> Cylinder:
        DF = self.DistalFrame()
        uhat = (DF @ SE3.Rx(np.pi/self.numSides)).R[:,1]
        return Cylinder(self.r, self.DistalFrame().t, -DF.R[:,0], 
                        self.tubeLength, uhat)

    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=jointColorDefault, edgeColor=jointEdgeColorDefault,
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
            RevoluteProximalPose = self.RevoluteProximalFrame()
            uhatProximal = RevoluteProximalPose.R[:,1]
            vhatProximal = RevoluteProximalPose.R[:,2]
            RevoluteProximalBase = RevoluteProximalPose.t + u.reshape(-1,1) @ uhatProximal.reshape(1,3) + v.reshape(-1,1) @ vhatProximal.reshape(1,3)
            
            ProximalPoints = np.vstack((RevoluteProximalBase, CenterSegment))
            ProximalHull = ConvexHull(ProximalPoints)
            for s in ProximalHull.simplices:
                tri = Poly3DCollection([ProximalPoints[s]])
                tri.set_color(surfaceColor)
                tri.set_alpha(surfaceOpacity)
                ax.add_collection3d(tri)
            
            RevoluteDistalPose = self.RevoluteDistalFrame()
            uhatDistal = RevoluteDistalPose.R[:,1]
            vhatDistal = RevoluteDistalPose.R[:,2]
            DistalBase = RevoluteDistalPose.t + u.reshape(-1,1) @ uhatDistal.reshape(1,3) + v.reshape(-1,1) @ vhatDistal.reshape(1,3)
            DistalPoints = np.vstack((DistalBase, CenterSegment))
            DistalHull = ConvexHull(DistalPoints)
            for s in DistalHull.simplices:
                tri = Poly3DCollection([DistalPoints[s]])
                tri.set_facecolor(surfaceColor)
                tri.set_edgecolor(edgeColor)
                tri.set_alpha(surfaceOpacity)
                ax.add_collection3d(tri)
            
            self.proximalCylinder().addToPlot(ax, color=surfaceColor, 
                                              alpha=surfaceOpacity, 
                                              edgeColor=edgeColor,
                                              numPointsPerCircle=self.numSides)
            self.distalCylinder().addToPlot(ax, color=surfaceColor, 
                                              alpha=surfaceOpacity, 
                                              edgeColor=edgeColor,
                                              numPointsPerCircle=self.numSides)
        return plotHandles
    
    def show(self, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=jointColorDefault, edgeColor=jointEdgeColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=True,
             axisScale=10, showPoses=True, block=blockDefault):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandles = self.addToPlot(ax, xColor, yColor, zColor, 
                                        proximalColor, centerColor, distalColor,
                                        sphereColor, showSphere, 
                                        surfaceColor, edgeColor,
                                        surfaceOpacity, showSurface, showAxis,
                                        axisScale, showPoses)
        if showPoses:
            xHats, yHats, zHats, origins = plotHandles
            ax.legend([xHats, yHats, zHats], [r'$\^x$', r'$\^y$', r'$\^z$'])
        ax.set_aspect('equal')
        plt.show(block=block)

    def toPrinted(self, screwRadius):
        from PrintedJoint import PrintedOrthogonalRevoluteJoint

        joint = PrintedOrthogonalRevoluteJoint(self.r, -self.totalBendingAngle/2, self.totalBendingAngle/2, self.Pose, screwRadius, initialState=self.initialState)
        if joint.neutralLength < self.neutralLength:
            amount = (self.neutralLength - joint.neutralLength)/2
            joint.extendBottomSegment(amount)
            joint.extendSegment(amount)
        return joint

    def collisionBox(self):
        prox = self.ProximalDubinsFrame()
        points = np.array([
            prox.t + self.r*prox.a + self.revoluteLength*prox.o,
            prox.t - self.r*prox.a + self.revoluteLength*prox.o,
            prox.t - self.r*prox.a - self.revoluteLength*prox.o,
            prox.t + self.r*prox.a - self.revoluteLength*prox.o,
            prox.t + self.r*prox.a + self.revoluteLength*prox.o + self.revoluteLength*prox.n,
            prox.t - self.r*prox.a + self.revoluteLength*prox.o + self.revoluteLength*prox.n,
            prox.t - self.r*prox.a - self.revoluteLength*prox.o + self.revoluteLength*prox.n,
            prox.t + self.r*prox.a - self.revoluteLength*prox.o + self.revoluteLength*prox.n,
        ])

        return CollisionBox(points)
    
    def getCapsules(self):
        return [CollisionCapsule(base=self.ProximalDubinsFrame(), radius=self.r, height=self.neutralLength/2),
                CollisionCapsule(base=self.DistalDubinsFrame(), radius=self.r, height = -self.neutralLength/2)]


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
             surfaceColor=jointColorDefault, edgeColor=jointEdgeColorDefault,
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
    
    def toPrinted(self, screwRadius):
        from PrintedJoint import PrintedPrismaticJoint
        return PrintedPrismaticJoint(self.r, self.maxLength - self.minLength, self.Pose, screwRadius, initialState=self.initialState)
    
    def getCapsules(self):
        return [CollisionCapsule(base=self.ProximalDubinsFrame(), radius=self.r, height=self.minLength)]
        
    
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
        return math.sqrt(self.r)
    
    def boundingBall(self) -> Ball:
        return Ball(self.Pose.t, self.boundingRadius())
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=jointColorDefault, edgeColor=jointEdgeColorDefault,
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
    
    def toPrinted(self, screwRadius):
        from PrintedJoint import PrintedWaypoint
        return PrintedWaypoint(self.r, self.Pose, screwRadius, pathIndex = self.pidx)
    
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
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=linkColorDefault, edgeColor=jointEdgeColorDefault,
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
            
            tipSegmentIndex, uhatIndex, vhatIndex = 1,0,1
            
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
                    tri.set_facecolor(linkColorDefault)
                    tri.set_edgecolor(None)
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
                    tri.set_facecolor(surfaceColor)
                    tri.set_edgecolor(edgeColor)
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
                # facing up
                tipSegment = np.array([distalPose.t - scale * distalPose.R[:, 0],
                                    distalPose.t + scale * distalPose.R[:, 0]])
                tipPoints = np.vstack((proximalBase, tipSegment))
            else:
                # facing down
                tipSegment = np.array([proximalPose.t - scale * proximalPose.R[:, 0],
                                    proximalPose.t + scale * proximalPose.R[:, 0]])
                tipPoints = np.vstack((distalBase, tipSegment))

            hull = ConvexHull(tipPoints)
            for s in hull.simplices:
                vertices = tipPoints[s]
                meshdata = gl.MeshData(vertexes=vertices, faces=[np.arange(len(vertices))])
                item = MeshItemWithID(meshdata=meshdata, color=surfaceColor, smooth=False, drawEdges=True, shader='shaded', glOptions='translucent', id=self.id)
                item.setObjectName("Joint")
                widget.plot_widget.addItem(item)
        
        super().addToWidget(widget, xColor, yColor, zColor, proximalColor,
                            centerColor, distalColor, sphereColor, showSphere,
                            surfaceColor, showSurface, showAxis,
                            axisScale, showPoses, poseAxisScaleMultipler)

    def toPrinted(self, screwRadius):
        from PrintedJoint import PrintedTip
        return PrintedTip(self.r, self.Pose, screwRadius)
    
class StartTip(Tip):
    def __init__(self, numSides : int, r : float, Pose : SE3, length : float):
        super().__init__(numSides, r, Pose, length, closesForward=False)
        print("start")

class EndTip(Tip):
    def __init__(self, numSides : int, r : float, Pose : SE3, length : float):
        super().__init__(numSides, r, Pose, length, closesForward=True)
        print("end")