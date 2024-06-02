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


    def addToPlot(self, ax, showWholeCircle=True, showStart=True, 
                    startColor='r', pathColor=pathColorDefault,
                    cscBoundaryMarker='*', alpha=1):
        if showStart:
            x1,y1,z1 = self.startPosition
            u1,v1,w1 = self.startDir
            ax.quiver(x1,y1,z1,u1,v1,w1, length=self.r, color=startColor, label='start')
        
        if not self.arc is None:
            if showWholeCircle:
                cx, cy, cz = Circle3D(self.r, self.arc.circleCenter, 
                                    self.arc.binormal).interpolate().T
                ax.plot(cx, cy, cz, color = startColor)
                
            # path arc (C component)
            self.arc.addToPlot(ax, pathColor, alpha=alpha)

        # path segment (S component)
        sx,sy,sz = self.segment.T
        ax.plot(sx, sy, sz, color = pathColor, 
                marker=cscBoundaryMarker, alpha=alpha)

    def show(self, showWholeCircle=True, showStart=True, 
                    startColor='r', pathColor=pathColorDefault,
                    cscBoundaryMarker='*', alpha=1, block=blockDefault):
        ax = plt.figure().add_subplot(projection='3d')
        self.addToPlot(ax, showWholeCircle, showPoses, startColor, 
                       pathColor, cscBoundaryMarker, alpha)
        ax.set_aspect('equal')
        ax.legend()
        plt.show(block=block)


"""
# Tests
PathCS(SE3.Ry(np.pi/3)@SE3.Trans([1,1,-2]), np.array([3,0,1]), 1).show()
PathCS(SE3.Ry(np.pi/3)@SE3.Trans([1,1,-2]), np.array([-3,-1,1]), 1).show()
PathCS(SE3.Rz(np.pi/4), np.array([-3,4,2]), 1).show()
PathCS(SE3.Ry(np.pi/3)@SE3.Trans([1,1,-2]), np.array([3,0,4]), 1).show()
PathCS(SE3.Ry(np.pi/3)@SE3.Trans([1,1,-2]), np.array([-3,-1,4]), 1).show()
PathCS(SE3.Rz(np.pi/4), np.array([-6,-3,-7]), 1).show()
"""