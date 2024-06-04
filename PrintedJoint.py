from spatialmath import SE3
from abc import ABC, abstractmethod
from geometryHelpers import *
import matplotlib.pyplot as plt
from TubularPattern import *
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
import Joint
from Joint import *

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
        thickness = screwRadius
        gridHoleMargin = screwRadius*2.5#printParameters.gridHoleRadius*2 + printParameters.holeMargin
        return PrintParameters(thickness, screwRadius, screwRadius/2, gridHoleMargin, r/2, 0.01, r, screwRadius, thickness, screwRadius/2)

class PrintedJoint(Joint):
    def __init__(self, r : float, neutralLength : float, Pose : SE3(), screwRadius : float, printParameters: PrintParameters, nextLink : PrintedJoint, initialState : float = 0):
        self.screwRadius = screwRadius
        self.printParameters = printParameters
        self.nextLink = nextLink
        super().__init__(r, neutralLength, Pose, initialState)

class PrintedOrthogonalRevoluteJoint(PrintedJoint):
    def __init__(self, r : float, totalBendingAngle : float, Pose : SE3, screwRadius : float, printParameters: PrintParameters = None, nextLink = None, initialState : float = 0):
        if printParameters == None:
            printParameters = PrintParameters.default(r, screwRadius)

        leftBendingAngle = -1*(initialState - totalBendingAngle/2)
        rightBendingAngle = initialState + totalBendingAngle/2
        self.maxBendingAngle = max(leftBendingAngle, rightBendingAngle)
        self.minBendingAngle = min(leftBendingAngle, rightBendingAngle)

        if (leftBendingAngle >= 180 | rightBendingAngle >= 180):
            raise InvalidDimensionException(f"Tried to construct an orthogonal revolute joint which can turn {self.maxBendingAngle} degrees in a singular direction. However, an orthogonal revolute joint can turn less than 180 degrees in each direction").

        maxTurnDist = r / tan((180 - self.maxBendingAngle) / 2)
        
        #subject to change, check scad
        self.bottomLength = printParameters.holeMargin*4 + screwRadius*4 + maxTurnDist + printParameters.attachThickness + printParameters.holeGridMargin
        self.topLength = printParameters.holeMargin*3 + printParameters.gridHoleRadius*4 + maxTurnDist
        neutralLength = bottomLength + topLength

        super().__init__(r, neutralLength, Pose, screwRadius, printParameters, nextLink, initialState)

    def extendTopSegment(self, amount):
        self.topLength += amount
    
    def extendSegment(self, amount):
        self.bottomLength += amount

    def recomputeDimensions(self):
        self.bottomLength = self.printParameters.holeMargin*4 + self.screwRadius*4 + maxTurnDist + self.printParameters.attachThickness + self.printParameters.holeGridMargin
        self.topLength = self.printParameters.holeMargin*3 + self.printParameters.gridHoleRadius*4 + maxTurnDist

    def export3DFile(self, index: int, fileFormat = "stl"):
        name = f"joint_{index}." + fileFormat

        defs = [f"tolerance={self.printParameters.tolerance};\n",f"hole_radius={self.screwRadius};\n",
                f"grid_hole_radius={self.printParameters.gridHoleRadius};\n",f"outer_radius={self.r};\n",
                f"thickness={self.printParameters.thickness};\n",f"hole_attach_height={self.printParameters.holeMargin};\n",
                f"attach_thickness={self.printParameters.attachThickness};\n", f"bottom_rotation_height={self.bottomLength};\n",
                f"top_rotation_height={self.topLength};\n", f"next_hole_radius={self.printParameters.nextScrewRadius};\n",
                f"next_hole_attach_height={self.printParameters.nextHoleMargin};\n", f"next_inner={self.printParameters.nextInnerRadius};\n",
                f"max_turn={self.maxBendingAngle};\n", f"min_turn={self.minBendingAngle};\n"]

        with open("scad/orthogonal_revolute.scad", "r") as file:
            lines = file.readlines()
        truncated = lines[14:] #first 14 are parameter definitions

        with open(f"scad_output/{name}", "w") as file:
            defs.extend(truncated)
            file.writelines(defs)
        
        os.system(f"openscad -o scad_output/{name} scad_output/{name}.scad")

        

class PrintedPrismaticJoint(PrintedJoint):
    def __init__(self, r : float, minLength : float, maxLength: float, Pose : SE3, screwRadius : float, printParameters: PrintParameters = None, nextLink = None, initialState : float = 0):
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

        super().__init__(r, minLength, Pose, screwRadius, printParameters, nextLink, initialState)

    def extendSegment(self, amount):
        self.minLength += amount

    def export3DFile(self, index: int, fileFormat = "stl"):
        name = f"joint_{index}." + fileFormat

        defs = [f"tolerance={self.printParameters.tolerance};\n",f"hole_radius={self.screwRadius};\n",
                f"grid_hole_radius={self.printParameters.gridHoleRadius};\n",f"outer_radius={self.r};\n",
                f"thickness={self.printParameters.thickness};\n",f"hole_attach_height={self.printParameters.holeMargin};\n",
                f"compressed_height={self.minLength};\n", f"extended_height={self.maxLength};\n",
                f"next_hole_radius={self.printParameters.nextScrewRadius};\n", f"next_hole_attach_height={self.printParameters.nextHoleMargin};\n", 
                f"next_inner={self.printParameters.nextInnerRadius};\n",]

        with open("scad/prismatic.scad", "r") as file:
            lines = file.readlines()
        truncated = lines[11:] #first 11 are parameter definitions

        with open(f"scad_output/{name}", "w") as file:
            defs.extend(truncated)
            file.writelines(defs)
        
        os.system(f"openscad -o scad_output/{name} scad_output/{name}.scad")

class PrintedTip(PrintedJoint):
    def __init__(self, r : float, Pose : SE3(), screwRadius : float, printParameters: PrintParameters = None):
        if printParameters == None:
            printParameters = PrintParameters.default(r, screwRadius)

        neutralLength = printParameters.holeMargin + screwRadius*2 + r
        super().__init__(r, neutralLength, Pose, screwRadius, printParameters, nextLink=None)

    def export3DFile(self, index: float, fileFormat = "stl"):
        name = f"tip_{index}." + fileFormat

        defs = [f"eps={self.printParameters.tolerance/100};\n",f"hole_attach_height={self.printParameters.holeMargin};\n",
                f"hole_radius={self.screwRadius};\n",f"outer_radius={self.r};\n",
                f"thickness={self.printParameters.thickness};\n"]

        with open("scad/tip.scad", "r") as file:
            lines = file.readlines()
        truncated = lines[5:] #first 5 are parameter definitions

        with open(f"scad_output/{name}", "w") as file:
            defs.extend(truncated)
            file.writelines(defs)
        
        os.system(f"openscad -o scad_output/{name} scad_output/{name}.scad")