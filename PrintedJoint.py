from spatialmath import SE3
from abc import ABC, abstractmethod
from geometryHelpers import *
import matplotlib.pyplot as plt
from TubularPattern import *
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
import Joint
from Joint import *
import os

class InvalidDimensionException(Exception):
    def __init__(self, message="Invalid Dimensions for 3D printed joint."):
        self.message = message
        super().__init__(self.message)

class PrintParameters:
    def __init__(self, thickness : float, gridHoleRadius : float, holeMargin : float, holeGridMargin: float, attachThickness : float, tolerance: float,
                nextR, nextScrewRadius, nextThickness, nextHoleMargin):
        self.thickness = thickness
        self.gridHoleRadius = gridHoleRadius
        self.holeMargin = holeMargin
        self.attachThickness = attachThickness
        self.tolerance = tolerance
        self.nextR = nextR
        self.nextScrewRadius = nextScrewRadius
        self.nextThickness = nextThickness
        self.nextHoleMargin = nextHoleMargin
        self.holeGridMargin = holeGridMargin
    
    def calculateNumHoles(self, outerRadius : float):
        #TODO
        return 0
    
    @staticmethod
    def default(r: float, screwRadius : float):
        thickness = r * 0.2
        gridHoleMargin = screwRadius*2.5#printParameters.gridHoleRadius*2 + printParameters.holeMargin
        return PrintParameters(thickness, screwRadius*2, screwRadius, gridHoleMargin, thickness, 0.01, r, screwRadius, thickness, screwRadius)

class PrintedJoint(Joint):
    def __init__(self, r : float, neutralLength : float, Pose : SE3(), screwRadius : float, printParameters: PrintParameters, initialState : float = 0):
        self.screwRadius = screwRadius
        self.printParameters = printParameters
        super().__init__(r, neutralLength, Pose, initialState)

class PrintedOrthogonalRevoluteJoint(PrintedJoint):
    def __init__(self, r : float, startBendingAngle : float, endBendingAngle : float, Pose : SE3, screwRadius : float, printParameters: PrintParameters = None, initialState : float = 0):
        if printParameters == None:
            printParameters = PrintParameters.default(r, screwRadius)

        self.startBendingAngle = startBendingAngle
        self.endBendingAngle = endBendingAngle
        self.maxBendingAngle = max(-startBendingAngle, endBendingAngle)
        self.minBendingAngle = min(-startBendingAngle, endBendingAngle)
        
        if self.maxBendingAngle >= np.pi or self.minBendingAngle <= -np.pi:
            raise InvalidDimensionException(f"Tried to construct an orthogonal revolute joint which can turn {self.maxBendingAngle} degrees in a singular direction. However, an orthogonal revolute joint can turn less than 180 degrees in each direction")

        maxTurnDist = r / np.tan((np.pi - self.maxBendingAngle) / 2)
        
        #subject to change, check scad
        self.bottomLength = printParameters.holeMargin*4 + screwRadius*4 + maxTurnDist + printParameters.attachThickness + printParameters.holeGridMargin
        self.topLength = printParameters.holeMargin*3 + printParameters.gridHoleRadius*4 + maxTurnDist
        neutralLength = self.bottomLength + self.topLength

        super().__init__(r, neutralLength, Pose, screwRadius, printParameters, initialState)

    def extendSegment(self, amount):
        self.topLength += amount
    
    def extendBottomSegment(self, amount):
        self.bottomLength += amount

    def recomputeDimensions(self):
        self.bottomLength = self.printParameters.holeMargin*4 + self.screwRadius*4 + maxTurnDist + self.printParameters.attachThickness + self.printParameters.holeGridMargin
        self.topLength = self.printParameters.holeMargin*3 + self.printParameters.gridHoleRadius*4 + maxTurnDist

    def export3DFile(self, index: int, folder = "", fileFormat = "stl"):
        name = f"joint_{index}." + fileFormat

        defs = [f"tolerance={self.printParameters.tolerance};\n",f"hole_radius={self.screwRadius};\n",
                f"grid_hole_radius={self.printParameters.gridHoleRadius};\n",f"outer_radius={self.r};\n",
                f"thickness={self.printParameters.thickness};\n",f"hole_attach_height={self.printParameters.holeMargin};\n",
                f"attach_thickness={self.printParameters.attachThickness};\n", f"bottom_rotation_height={self.bottomLength};\n",
                f"top_rotation_height={self.topLength};\n", f"next_hole_radius={self.printParameters.nextScrewRadius};\n",
                f"next_hole_attach_height={self.printParameters.nextHoleMargin};\n", f"next_inner={self.printParameters.nextR - self.printParameters.nextThickness};\n",
                f"max_turn={np.rad2deg(self.maxBendingAngle)};\n", f"min_turn={np.rad2deg(self.minBendingAngle*180/np.pi)};\n"]

        with open("scad/orthogonal_revolute.scad", "r") as file:
            lines = file.readlines()
        truncated = lines[14:] #first 14 are parameter definitions

        with open(f"scad_output/{folder}{name}.scad", "w+") as file:
            defs.extend(truncated)
            file.writelines(defs)
        
        os.system(f"openscad -q -o 3d_output/{folder}{name} scad_output/{folder}/{name}.scad > /dev/null")
    
    def pathIndex(self) -> int:
        return 0 # xhat
    
    def stateRange(self) -> list:
        return [self.startBendingAngle, self.endBendingAngle]
    
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
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=jointColorDefault, edgeColor=jointEdgeColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=True,
             axisScale=10, showPoses=True):
        plotHandles = super().addToPlot(ax, xColor, yColor, zColor, proximalColor,
                          centerColor, distalColor, sphereColor, showSphere,
                          surfaceColor, surfaceOpacity, showSurface, showAxis,
                          axisScale, showPoses)
        if showSurface:
            scale = self.r/2
            CenterSegment = np.array([self.Pose.t - scale * self.Pose.R[:,2],
                                      self.Pose.t + scale * self.Pose.R[:,2]])
            #https://stackoverflow.com/questions/63207496/how-to-visualize-polyhedrons-defined-by-their-vertices-in-3d-with-matplotlib-or
            
            radialCount = 5
            angle = np.linspace(0, 2*np.pi, radialCount) + np.pi/4
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
                tri.set_facecolor(surfaceColor)
                tri.set_edgecolor(edgeColor)
                tri.set_alpha(surfaceOpacity)
                ax.add_collection3d(tri)
            
        return plotHandles

        

class PrintedPrismaticJoint(PrintedJoint):
    def __init__(self, r : float, minLength : float, maxLength: float, Pose : SE3, screwRadius : float, printParameters: PrintParameters = None, initialState : float = 0):
        if printParameters == None:
            printParameters = PrintParameters.default(r, screwRadius)
        
        #subject to change, check scad
        outHoleAngle = 60
        guideHeight = printParameters.holeMargin*2 + screwRadius*2 + printParameters.gridHoleRadius*2 + printParameters.grid_hole_radius*2/sin(outHoleAngle) + printParameters.thickness/tan(outHoleAngle)
        maxPossibleLength = minLength*2 - guideHeight - printParameters.holeMargin*2 - printParameters.gridHoleRadius
        if (maxLength > maxPossibleLength):
            raise InvalidDimensionException(f"Tried to construct a prismatic joint with minimum length of {minLength} and a maximum length of {maxLength}, yet with current parameters a maximum length of at most {maxPossibleLength} is possible")
        
        self.minLength = minLength
        self.maxLength = maxLength

        super().__init__(r, minLength, Pose, screwRadius, printParameters, initialState)

    def extendSegment(self, amount):
        self.minLength += amount

    def export3DFile(self, index: int, folder = "", fileFormat = "stl"):
        name = f"joint_{index}." + fileFormat

        defs = [f"tolerance={self.printParameters.tolerance};\n",f"hole_radius={self.screwRadius};\n",
                f"grid_hole_radius={self.printParameters.gridHoleRadius};\n",f"outer_radius={self.r};\n",
                f"thickness={self.printParameters.thickness};\n",f"hole_attach_height={self.printParameters.holeMargin};\n",
                f"compressed_height={self.minLength};\n", f"extended_height={self.maxLength};\n",
                f"next_hole_radius={self.printParameters.nextScrewRadius};\n", f"next_hole_attach_height={self.printParameters.nextHoleMargin};\n", 
                f"next_inner={self.printParameters.nextR - self.printParameters.nextThickness};\n",]

        with open("scad/prismatic.scad", "r") as file:
            lines = file.readlines()
        truncated = lines[11:] #first 11 are parameter definitions

        with open(f"scad_output/{folder}{name}.scad", "w+") as file:
            defs.extend(truncated)
            file.writelines(defs)
        
        os.system(f"openscad -q -o 3d_output/{folder}{name} scad_output/{folder}{name}.scad > /dev/null")

class PrintedTip(PrintedJoint):
    def __init__(self, r : float, Pose : SE3(), screwRadius : float, printParameters: PrintParameters = None):
        if printParameters == None:
            printParameters = PrintParameters.default(r, screwRadius)

        neutralLength = printParameters.holeMargin + screwRadius*2 + r
        super().__init__(r, neutralLength, Pose, screwRadius, printParameters)

    def export3DFile(self, index: float, folder = "", fileFormat = "stl"):
        name = f"tip_{index}." + fileFormat

        defs = [f"eps={self.printParameters.tolerance/100};\n",f"hole_attach_height={self.printParameters.holeMargin};\n",
                f"hole_radius={self.screwRadius};\n",f"outer_radius={self.r};\n",
                f"thickness={self.printParameters.thickness};\n"]

        with open("scad/tip.scad", "r") as file:
            lines = file.readlines()
        truncated = lines[5:] #first 5 are parameter definitions

        with open(f"scad_output/{folder}{name}.scad", "w+") as file:
            defs.extend(truncated)
            file.writelines(defs)
        
        os.system(f"openscad -q -o 3d_output/{folder}{name} scad_output/{folder}{name}.scad > /dev/null")
    
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
             surfaceColor=jointColorDefault, edgeColor=jointEdgeColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=False,
             axisScale=10, showPoses=True):
        plotHandles = super().addToPlot(ax, xColor, yColor, zColor, proximalColor,
                          centerColor, distalColor, sphereColor, showSphere,
                          surfaceColor, surfaceOpacity, showSurface, showAxis,
                          axisScale, showPoses)
        if showSurface:
            #https://stackoverflow.com/questions/63207496/how-to-visualize-polyhedrons-defined-by-their-vertices-in-3d-with-matplotlib-or
            radialCount = 5
            angle = np.linspace(0, 2*np.pi, radialCount) + np.pi/4
            u = self.r * np.cos(angle)
            v = self.r * np.sin(angle)
            scale = self.r / 2
            
            tipSegmentIndex, uhatIndex, vhatIndex = 1,0,1

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
                tri.set_facecolor(surfaceColor)
                tri.set_edgecolor(edgeColor)
                tri.set_alpha(surfaceOpacity)
                ax.add_collection3d(tri)
            
        return plotHandles


class PrintedWaypoint(PrintedJoint):
    # path direction through a waypoint defaults to zhat
    def __init__(self, r : float, Pose : SE3, screwRadius : float, pathIndex : int = 2, printParameters: PrintParameters = None):
        if printParameters == None:
            printParameters = PrintParameters.default(r, screwRadius)

        assert(pathIndex in [0,1,2])
        self.pidx = pathIndex
        super().__init__(r, 0, Pose, screwRadius, printParameters, 0)

        #temporary
        self.numSides = 4

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