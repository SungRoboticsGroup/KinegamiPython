# -*- coding: utf-8 -*-
import numpy as np
from numpy import cross, dot, arctan2
import scipy
from scipy.linalg import null_space
from numpy.linalg import norm
import matplotlib.pyplot as plt
import geometryHelpers
from geometryHelpers import *
from scipy.optimize import minimize

class PathCS:
    def __init__(self, startPose : SE3, endPosition : np.ndarray, r : float):
        assert(r>0)
        self.r = r
        self.startPose = startPose
        self.startPosition = startPose.t
        self.startDir = startPose.R[:,0]
        self.endPosition = endPosition

        StartNormalPlane = Plane(self.startPosition, self.startDir)
        circleDirection = StartNormalPlane.projectionOfPoint(endPosition) - self.startPosition

        if norm(circleDirection) < 1e-6:
            self.arc = None
            self.segment = np.array([self.startPosition, endPosition])
        else:
            circleDirection = circleDirection / norm(circleDirection)
            circleCenter = self.startPosition + r*circleDirection
            
            def distanceFromTangent(theta):
                theta = math.remainder(theta[0], 2*np.pi-1e-6)
                arc = Arc3D(circleCenter, self.startPosition, self.startDir, theta)
                tangentRay = Ray(arc.endPoint, arc.endTangent)
                return tangentRay.distanceToPoint(endPosition)
            
            solution = minimize(distanceFromTangent, np.pi/4, method='Nelder-Mead')
            assert(solution.success)
            theta = math.remainder(solution.x[0], 2*np.pi-1e-6)

            arc = Arc3D(circleCenter, self.startPosition, self.startDir, theta)
            segment = np.array([arc.endPoint, endPosition])
            self.arc = arc
            self.segment = segment
        
        segmentDir = self.segment[1] - self.segment[0]
        if norm(segmentDir) < 1e-6:
            assert(not self.arc is None)
            segmentDir = self.arc.endTangent
        self.segmentDir = segmentDir / norm(segmentDir)
