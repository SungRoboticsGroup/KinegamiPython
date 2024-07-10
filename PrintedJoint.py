from spatialmath import SE3
from abc import ABC, abstractmethod
from geometryHelpers import *
import matplotlib.pyplot as plt
from TubularPattern import *
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
from scipy.spatial.transform import Rotation
import Joint
from Joint import *
from CollisionDetection import *
import os

class InvalidDimensionException(Exception):
    def __init__(self, message="Invalid Dimensions for 3D printed joint."):
        self.message = message
        super().__init__(self.message)

class PrintParameters:
    def __init__(self, thickness : float, gridHoleRadius : float, holeMargin : float, holeGridMargin: float, attachThickness : float, tolerance: float,
                nextR, nextScrewRadius, nextThickness, nextHoleMargin):
        self.thickness = np.round(thickness, 4)
        self.gridHoleRadius = np.round(gridHoleRadius, 4)
        self.holeMargin = np.round(holeMargin, 4)
        self.attachThickness = np.round(attachThickness, 4)
        self.tolerance = np.round(tolerance, 4)
        self.nextR = np.round(nextR, 4)
        self.nextScrewRadius = np.round(nextScrewRadius, 4)
        self.nextThickness = np.round(nextThickness, 4)
        self.nextHoleMargin = np.round(nextHoleMargin, 4)
        self.holeGridMargin = np.round(holeGridMargin, 4)
    
    def toString(self):
        return str(self.thickness) + ";" + str(self.gridHoleRadius) + ";" + str(self.holeMargin) + ";" + str(self.attachThickness) + ";" + str(self.tolerance) + ";" + str(self.nextR) + ";" + str(self.nextScrewRadius) + ";" + str(self.nextThickness) + ";" + str(self.nextHoleMargin)
    
    @staticmethod
    def default(r: float, screwRadius : float):
        thickness = r * 0.2
        gridHoleMargin = screwRadius*2.5
        return PrintParameters(thickness, screwRadius, screwRadius, gridHoleMargin, thickness/2, 0.075, r, screwRadius, thickness, screwRadius)

    @staticmethod
    def fromString(string):
        x = [float(param) for param in string.split(';')]
        return PrintParameters(x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9])

class PrintedJoint(Joint):
    def __init__(self, r : float, neutralLength : float, Pose : SE3(), screwRadius : float, printParameters: PrintParameters, initialState : float = 0):
        self.screwRadius = np.round(screwRadius, 4)
        self.printParameters = printParameters
        self.twistAngle = 0
        super().__init__(r, neutralLength, Pose, initialState)
    
    def setTwistAngle(self, angle):
        self.twistAngle = angle

    def toOrigami(self, numSides : int, numLayers : int = 1):
        pass

    def extendSegment(self, amount):
        pass

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
        name = f"orthogonal_revolute_{index}." + fileFormat

        defs = [f"tolerance={self.printParameters.tolerance};\n",f"hole_radius={self.screwRadius};\n",
                f"grid_hole_radius={self.printParameters.gridHoleRadius};\n",f"outer_radius={self.r};\n",
                f"thickness={self.printParameters.thickness};\n",f"hole_attach_height={self.printParameters.holeMargin};\n",
                f"attach_thickness={self.printParameters.attachThickness};\n", f"bottom_rotation_height={np.round(self.bottomLength, 4)};\n",
                f"top_rotation_height={np.round(self.topLength, 4)};\n", f"next_hole_radius={self.printParameters.nextScrewRadius};\n",
                f"next_hole_attach_height={self.printParameters.nextHoleMargin};\n", f"next_inner={self.printParameters.nextR - self.printParameters.nextThickness};\n",
                f"max_turn={np.round(np.rad2deg(self.maxBendingAngle), 4)};\n", f"min_turn={np.round(np.rad2deg(self.minBendingAngle), 4)};\n",
                f"hole_twist={np.round(np.rad2deg(self.twistAngle) % 90, 4)};\n"]

        with open("scad/orthogonal_revolute.scad", "r") as file:
            lines = file.readlines()
        truncated = lines[15:] #first 15 are parameter definitions

        with open(f"scad_output/{folder}/{name}.scad", "w+") as file:
            defs.extend(truncated)
            file.writelines(defs)
        
        os.system(f"openscad -q -o 3d_output/{folder}/{name} scad_output/{folder}/{name}.scad")

    def renderPose(self, folder):
        rot = SE3.Ry(np.pi/2)
        name = f"orthogonal_revolute_"

        defs = [f"tolerance={self.printParameters.tolerance};\n",f"hole_radius={self.screwRadius};\n",
                f"grid_hole_radius={self.printParameters.gridHoleRadius};\n",f"outer_radius={self.r};\n",
                f"thickness={self.printParameters.thickness};\n",f"hole_attach_height={self.printParameters.holeMargin};\n",
                f"attach_thickness={self.printParameters.attachThickness};\n", f"bottom_rotation_height={np.round(self.bottomLength, 4)};\n",
                f"top_rotation_height={np.round(self.topLength, 4)};\n", f"next_hole_radius={self.printParameters.nextScrewRadius};\n",
                f"next_hole_attach_height={self.printParameters.nextHoleMargin};\n", f"next_inner={self.printParameters.nextR - self.printParameters.nextThickness};\n",
                f"max_turn={np.round(np.rad2deg(self.maxBendingAngle), 4)};\n", f"min_turn={np.round(np.rad2deg(self.minBendingAngle), 4)};\n",
                f"hole_twist={np.round(np.rad2deg(self.twistAngle) % 90, 4)};\n"]
        defs2 = defs.copy()

        for parameter in defs:
            name += parameter[parameter.index("=") + 1:-2] + "_"

        #orient bend the right way
        if self.maxBendingAngle == self.endBendingAngle:
            rot = rot @ SE3.Rz(np.pi)
        #check if already has rendered somewhere
        if os.path.isfile(f"3d_output/{folder}/poses/{name}1.stl") and os.path.isfile(f"3d_output/{folder}/poses/{name}2.stl"):
            return f"3d_output/{folder}/poses/{name}1.stl", rot, f"3d_output/{folder}/poses/{name}2.stl", rot      
        
        with open("scad/poses/orthogonal_revolute_pose1.scad", "r") as file:
            lines = file.readlines()
            truncated = lines[15:] #first 15 are parameter definitions

        with open(f"scad_output/{folder}/poses/{name}1.scad", "w+") as file:
            defs.extend(truncated)
            file.writelines(defs)

        os.system(f"openscad -q -o 3d_output/{folder}/poses/{name}1.stl scad_output/{folder}/poses/{name}1.scad")

        with open("scad/poses/orthogonal_revolute_pose2.scad", "r") as file:
            lines2 = file.readlines()
            truncated2 = lines2[15:] #first 15 are parameter definitions

        with open(f"scad_output/{folder}/poses/{name}2.scad", "w+") as file:
            defs2.extend(truncated2)
            file.writelines(defs2)
        
        os.system(f"openscad -q -o 3d_output/{folder}/poses/{name}2.stl scad_output/{folder}/poses/{name}2.scad")

        return f"3d_output/{folder}/poses/{name}1.stl", rot, f"3d_output/{folder}/poses/{name}2.stl", rot
    
    def pathIndex(self) -> int:
        return 0 # xhat
    
    def stateRange(self) -> list:
        return [self.startBendingAngle, self.endBendingAngle]
    
    def stateChangeTransformation(self, stateChange : float) -> SE3:
        return RotationAboutLine(rotAxisDir=self.Pose.R[:,2],
                              rotAxisPoint=(self.ProximalDubinsFrame() @ SE3.Trans(self.bottomLength, 0, 0)).t,#rotate about axis, not necessarily center
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

    def toOrigami(self, numSides : int, numLayers : int):
        from OrigamiJoint import ExtendedRevoluteJoint
        return ExtendedRevoluteJoint(numSides, self.r, max(-self.startBendingAngle, self.endBendingAngle)*2, self.bottomLength, self.Pose, numSinkLayers=numLayers, initialState=self.initialState)

    def getCapsules(self):
        return [CollisionCapsule(base=self.ProximalDubinsFrame(), radius=self.r, height=self.bottomLength),
                CollisionCapsule(base=self.DistalDubinsFrame(), radius=self.r, height = -self.topLength)]

class PrintedInAxisRevoluteJoint(PrintedJoint):
    def __init__(self, r: float, neutralLength : float, Pose : SE3, screwRadius : float, printParameters : PrintParameters = None, initialState : float = 0):
        if printParameters == None:
            printParameters = PrintParameters.default(r, screwRadius)
        
        super().__init__(r, neutralLength, Pose, screwRadius, printParameters, initialState)
    
    def pathIndex(self) -> int:
        return 2 # zhat
    
    def stateChangeTransformation(self, stateChange : float) -> SE3:
        return RotationAboutLine(rotAxisDir=self.Pose.R[:,0],
                              rotAxisPoint=self.ProximalDubinsFrame().t,
                              angle=stateChange)
    
    def stateRange(self) -> list:
        return [0, np.pi*2]
    
    def boundingRadius(self) -> float:
        return norm([self.r, self.neutralLength / 2])
    
    def center(self) -> np.ndarray:
        return self.Pose.t + (self.state/2) * self.pathDirection()
    
    def boundingBall(self) -> Ball:
        return Ball(self.center(), self.boundingRadius())
    
    def boundingCylinder(self) -> Cylinder:
        uhat = (self.Pose @ SE3.Rz(np.pi/4)).R[:,1]
        return Cylinder(self.r, self.ProximalFrame().t, self.pathDirection(), 
                        self.neutralLength, uhat)
    
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
                                              numPointsPerCircle=4)
            
        return plotHandles
    
    def getCapsules(self):
        return [CollisionCapsule(base=self.ProximalDubinsFrame(), radius=self.r, height=self.neutralLength)]

class PrintedPrismaticJoint(PrintedJoint):
    def __init__(self, r : float, extensionLength : float, Pose : SE3, screwRadius : float,  printParameters: PrintParameters = None, initialState : float = 0):
        if printParameters == None:
            printParameters = PrintParameters.default(r, screwRadius)
        
        #subject to change, check scad
        outHoleAngle = np.rad2deg(60)
        guideHeight = printParameters.holeMargin*2 + screwRadius*2 + printParameters.gridHoleRadius*2 + printParameters.gridHoleRadius*2/np.sin(outHoleAngle) + printParameters.thickness/np.tan(outHoleAngle)

        self.minLength = extensionLength + guideHeight + printParameters.holeMargin*2 + printParameters.gridHoleRadius
        self.maxLength = self.minLength*2 - guideHeight - printParameters.holeMargin*2 - printParameters.gridHoleRadius

        # maxPossibleLength = minLength*2 - guideHeight - printParameters.holeMargin*2 - printParameters.gridHoleRadius
        # if (maxLength > maxPossibleLength):
        #     raise InvalidDimensionException(f"Tried to construct a prismatic joint with minimum length of {minLength} and a maximum length of {maxLength}, yet with current parameters a maximum length of at most {maxPossibleLength} is possible")

        # self.minLength = minLength
        # self.maxLength = maxLength

        super().__init__(r, self.minLength, Pose, screwRadius, printParameters, initialState)

    def extendSegment(self, amount):
        self.minLength += amount

    def export3DFile(self, index: int, folder : str, fileFormat = "stl"):
        name = f"prismatic_{index}." + fileFormat

        if os.path.isfile(f"3d_output/{folder}/{name}"):
            return

        defs = [f"tolerance={self.printParameters.tolerance};\n",f"hole_radius={self.screwRadius};\n",
                f"grid_hole_radius={self.printParameters.gridHoleRadius};\n",f"outer_radius={self.r};\n",
                f"thickness={self.printParameters.thickness};\n",f"hole_attach_height={self.printParameters.holeMargin};\n",
                f"compressed_height={np.round(self.minLength, 4)};\n", f"extended_height={np.round(self.maxLength, 4)};\n",
                f"next_hole_radius={self.printParameters.nextScrewRadius};\n", f"next_hole_attach_height={self.printParameters.nextHoleMargin};\n", 
                f"next_inner={self.printParameters.nextR - self.printParameters.nextThickness};\n",f"hole_twist={np.round(np.rad2deg(self.twistAngle) % 90, 4)};"]

        with open("scad/prismatic.scad", "r") as file:
            lines = file.readlines()
        truncated = lines[12:] #first 12 are parameter definitions

        with open(f"scad_output/{folder}{name}.scad", "w+") as file:
            defs.extend(truncated)
            file.writelines(defs)
        
        os.system(f"openscad -q -o 3d_output/{folder}/{name} scad_output/{folder}/{name}.scad")
    
    def renderPose(self, folder):
        rot = SE3.Ry(np.pi/2)
        name = f"prismatic_"

        defs = [f"tolerance={self.printParameters.tolerance};\n",f"hole_radius={self.screwRadius};\n",
                f"grid_hole_radius={self.printParameters.gridHoleRadius};\n",f"outer_radius={self.r};\n",
                f"thickness={self.printParameters.thickness};\n",f"hole_attach_height={self.printParameters.holeMargin};\n",
                f"compressed_height={np.round(self.minLength, 4)};\n", f"extended_height={np.round(self.maxLength, 4)};\n",
                f"next_hole_radius={self.printParameters.nextScrewRadius};\n", f"next_hole_attach_height={self.printParameters.nextHoleMargin};\n", 
                f"next_inner={self.printParameters.nextR - self.printParameters.nextThickness};\n",f"hole_twist={np.round(np.rad2deg(self.twistAngle) % 90, 4)};"]
        defs2 = defs.copy()

        for parameter in defs:
            name += parameter[parameter.index("=") + 1:-2] + "_"

        #check if already has rendered somewhere
        if os.path.isfile(f"3d_output/{folder}/poses/{name}1.stl") and os.path.isfile(f"3d_output/{folder}/poses/{name}2.stl"):
            return f"3d_output/{folder}/poses/{name}1.stl", rot, f"3d_output/{folder}/poses/{name}2.stl", rot      
        
        with open("scad/poses/prismatic_pose1.scad", "r") as file:
            lines = file.readlines()
            truncated = lines[12:] #first 12 are parameter definitions

        with open(f"scad_output/{folder}/poses/{name}1.scad", "w+") as file:
            defs.extend(truncated)
            file.writelines(defs)

        os.system(f"openscad -q -o 3d_output/{folder}/poses/{name}1.stl scad_output/{folder}/poses/{name}1.scad")

        with open("scad/poses/prismatic_pose2.scad", "r") as file:
            lines2 = file.readlines()
            truncated2 = lines2[12:] #first 12 are parameter definitions

        with open(f"scad_output/{folder}/poses/{name}2.scad", "w+") as file:
            defs2.extend(truncated2)
            file.writelines(defs2)
        
        os.system(f"openscad -q -o 3d_output/{folder}/poses/{name}2.stl scad_output/{folder}/poses/{name}2.scad")

        return f"3d_output/{folder}/poses/{name}1.stl", rot, f"3d_output/{folder}/poses/{name}2.stl", rot
    
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
        uhat = (self.Pose @ SE3.Rz(np.pi/4)).R[:,1]
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
                                              numPointsPerCircle=4)
            
        return plotHandles

    def toOrigami(self, numSides : float, numLayers : int):
        from OrigamiJoint import PrismaticJoint
        coneAngle = np.arcsin(self.neutralLength/self.maxLength)
        return PrismaticJoint(numSides, self.r, self.neutralLength, numLayers, coneAngle, self.Pose, initialState=self.initialState)

    def getCapsules(self):
        return [CollisionCapsule(base=self.ProximalDubinsFrame(), radius=self.r, height=self.minLength)]

class PrintedTip(PrintedJoint):
    def __init__(self, r : float, Pose : SE3(), screwRadius : float, pathIndex = 2, printParameters: PrintParameters = None):
        if printParameters == None:
            printParameters = PrintParameters.default(r, screwRadius)

        neutralLength = printParameters.holeMargin + screwRadius*2 + r

        self.pidx = pathIndex

        super().__init__(r, neutralLength, Pose, screwRadius, printParameters)

    def export3DFile(self, index: float, folder = "", fileFormat = "stl"):
        name = f"tip_{index}." + fileFormat

        defs = [f"eps={self.printParameters.tolerance/100};\n",f"hole_attach_height={self.printParameters.holeMargin};\n",
                f"hole_radius={self.screwRadius};\n",f"outer_radius={self.r};\n",
                f"thickness={self.printParameters.thickness};\n", f"hole_twist={np.round(np.rad2deg(self.twistAngle) % 90, 4)};\n"]

        with open("scad/tip.scad", "r") as file:
            lines = file.readlines()
        truncated = lines[6:] #first 6 are parameter definitions

        with open(f"scad_output/{folder}/{name}.scad", "w+") as file:
            defs.extend(truncated)
            file.writelines(defs)
        
        os.system(f"openscad -q -o 3d_output/{folder}/{name} scad_output/{folder}/{name}.scad")
    
    def renderPose(self, folder):
        rot1 = SE3.Ry(np.pi/2)

        name = f"tip_"

        defs = [f"eps={self.printParameters.tolerance/100};\n",f"hole_attach_height={self.printParameters.holeMargin};\n",
                f"hole_radius={self.screwRadius};\n",f"outer_radius={self.r};\n",
                f"thickness={self.printParameters.thickness};\n", f"hole_twist={np.round(np.rad2deg(self.twistAngle) % 90, 4)};\n"]
        
        for parameter in defs:
            name += parameter[parameter.index("=") + 1:-2] + "_"
        
        filepath = f"3d_output/{folder}/poses/{name}.stl"
        #check if already has rendered somewhere
        if os.path.isfile(filepath):
            return filepath, rot1, None, None

        

        with open("scad/poses/tip_pose.scad", "r") as file:
            lines = file.readlines()
        truncated = lines[6:] #first 6 are parameter definitions

        with open(f"scad_output/{folder}/poses/{name}.scad", "w+") as file:
            defs.extend(truncated)
            file.writelines(defs)
        
        os.system(f"openscad -q -o {filepath} scad_output/{folder}/poses/{name}.scad")

        return filepath, rot1, None, None 

    def pathIndex(self) -> int:
        return self.pidx
    
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
            
            match self.pidx:
                case 0: tipSegmentIndex, uhatIndex, vhatIndex = 2,1,2
                case 1: tipSegmentIndex, uhatIndex, vhatIndex = 0,2,0
                case 2: tipSegmentIndex, uhatIndex, vhatIndex = 1,0,1

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

    def toOrigami(self, numSides : int, numLayers: int = 1):
        from OrigamiJoint import Tip
        return Tip(numSides, self.r, self.Pose, self.neutralLength)

    def getCapsules(self):
        return [CollisionCapsule(base=self.ProximalDubinsFrame(), radius=self.r, height=self.neutralLength)]

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

    def toOrigami(self, numSides : float, numLayers: int = 1):
        from OrigamiJoint import Waypoint
        return Waypoint(numSides, self.r, self.Pose, self.pidx)