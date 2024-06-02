from spatialmath import SE3
from abc import ABC, abstractmethod
from geometryHelpers import *
import matplotlib.pyplot as plt
from TubularPattern import *
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
import Joint
from Joint import *

class PrintParameters:
    def __init__(self, thickness : float, gridHoleRadius : float):
        self.thickness = thickness
        self.gridHoleRadius = gridHoleRadius
    
    def calculateNumHoles(self, outerRadius : float):
        #TODO
        return 0

class PrintedJoint(Joint):
    def __init__(self, r : float, neutralLength : float, Pose : SE3(), screwRadius : float, initialState : float = 0):
        super().__init__(r, neutralLength, Pose, initialState)


