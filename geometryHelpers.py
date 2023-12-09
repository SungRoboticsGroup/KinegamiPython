# -*- coding: utf-8 -*-
"""
Assorted geometry-related helper functions and classes
"""

import numpy as np
from numpy import cross, dot, arctan2
import scipy
from scipy.spatial.transform import Rotation
from scipy.linalg import null_space
from numpy.linalg import norm
from spatialmath import SO3, SE3
import matplotlib.pyplot as plt


def unit(v):
    return v / np.linalg.norm(v)


"""
Output n points spaced evenly from start to end.
"""
def lerpN(start, end, n):
    assert(start.shape == end.shape)
    d = start.shape[0]
    if n==0:
        return np.empty((0,d))
    elif n==1:
        return np.array([start + (end-start)/2])
    else:
        ts = (np.arange(n)/(n-1)).reshape((-1,1))
        return start + ts*(end-start)


"""
Given a numpy array whose rows are points, and an array of the same points
cycled by a roll of k, find k. If this fails (i.e. if the points don't
actually match via a roll) return None.
"""
def rollToMatch(points, pointsCycled, EPSILON=0.00001):
    n = points.shape[0]
    assert(points.shape == pointsCycled.shape)
    for k in range(n):
        if np.all(np.abs(np.roll(points, k, axis=0) - pointsCycled) < EPSILON):
            return k
    return None


""" 
Given 3D np vectors u and v, return a unit vector orthogonal to both, 
obeying the right-hand-rule if applicable 
(i.e., if neither is 0 and they aren't colinear).
"""
def unitNormalToBoth(u, v):
    # ensure correct dimensions
    u = u.flatten()
    v = v.flatten()
    assert(u.shape==(3,))
    assert(v.shape==(3,))
    
    cp = cross(u,v)
    # if u and v are not colinear and both nonzero
    if norm(cp) > 0:
        return cp / norm(cp)
    else: 
        # return any vector orthogonal to both
        nullSpaceBasis = null_space(np.array([u, v]))
        # any nullspace basis vector will work, 
        # we'll use the first one since it may be the only one
        return nullSpaceBasis[:,0].flatten()


"""
For A and B matrices (of the same shape) storing vectors as rows,
return the unsigned angle between corresponding rows
"""
def unsignedAngles(A,B):
    assert(A.shape == B.shape)
    AnormalizedRows = (A.T * (1/norm(A, axis=1))).T
    BnormalizedRows = (B.T * (1/norm(B, axis=1))).T
    dotProductsOfRows = np.sum(AnormalizedRows*BnormalizedRows, axis=1)
    return np.arccos(np.clip(dotProductsOfRows, -1.0, 1.0))

"""
Unsigned angle between two vectors, based on
https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
"""
def unsignedAngle(a,b):
    return np.arccos(np.clip(np.dot(a/norm(a), b/norm(b)), -1.0, 1.0))
    
    
"""
Signed angle (radians) from vector a to vector b, around the normal vector n.
All inputs should be numpy arrays of shape (3,)
"""
def signedAngle(a, b, n):
    a = a / norm(a)
    b = b / norm(b)
    if norm(n)>0:
        n = n / norm(n)
        
    return arctan2(dot(cross(a,b),n), dot(a,b))


"""
For A and B matrices (of the same shape) storing 2D vectors as rows,
return the signed angle between corresponding rows
"""
def signedAngles2D(A, B):
    assert(A.shape == B.shape)
    k = A.shape[0]
    A3 = np.hstack((A, np.zeros((k,1))))
    B3 = np.hstack((B, np.zeros((k,1))))
    return np.array([signedAngle(A3[i,:], B3[i,:], [0,0,1])
                     for i in range(k)])

def signedAngle2D(a, b):
    return signedAngle(np.hstack((a,[0])), np.hstack((b,[0])), [0,0,1])

"""
For A and B matrices (of the same shape) storing 3D vectors as rows,
return the signed angle about n between corresponding rows
"""
def signedAngles3D(A, B, n):
    assert(A.shape == B.shape)
    k = A.shape[0]
    return np.array([signedAngle(A[i], B[i], n) for i in range(k)]).flatten()

# wrap angles to [0,2pi)
def wrapAngle(angle):
    return angle % (2*np.pi)

class Circle3D:
    def __init__(self, radius, center, normal):
        assert(norm(normal)>0)
        self.r = radius
        self.c = center
        self.n = normal / norm(normal)
    
    def interpolate(self, count=50):
        angle = np.linspace(0, 2*np.pi, count).reshape(-1,1)
        u = self.r * np.cos(angle)
        v = self.r * np.sin(angle)
        
        # construct basis for circle plane
        uhat = unitNormalToBoth(self.n, self.n).reshape(1,3)
        vhat = cross(self.n, uhat).reshape(1,3)
        
        # 3d circle points
        return self.c + u @ uhat + v @ vhat

class Cylinder:
    def __init__(self, radius : float, start : np.ndarray, 
                 direction : np.ndarray, length : float):
        assert(length > 0)
        assert(norm(direction)>0)
        self.start = start
        self.direction = direction / norm(direction)
        self.length = length
        self.end = start + self.length * self.direction
        self.r = radius
    
    def interpolateCircles(self, numPointsPerCircle=12, numCircles=2):
        radialCount = numPointsPerCircle+1 #the +1 is because the first equals the last
        angle = np.linspace(0, 2*np.pi, radialCount) 
        u = self.r * np.cos(angle)
        v = self.r * np.sin(angle)
        circlePlaneBasis = null_space([self.direction])
        uhat = circlePlaneBasis[:,0]
        vhat = circlePlaneBasis[:,1]
        circle = u.reshape(-1,1) @ uhat.reshape(1,3) + v.reshape(-1,1) @ vhat.reshape(1,3)
        
        segment = np.linspace(self.start, self.end, numCircles)
        circlePoints = np.tile(circle, (numCircles,1)) + np.repeat(segment, radialCount, axis=0)
        return circlePoints.reshape((numCircles, radialCount, 3))
    
    def addToPlot(self, ax, numPointsPerCircle=12, numCircles=2, color='black', alpha=0.5, frame=False):
        circles = self.interpolateCircles(numPointsPerCircle, numCircles)
        X = circles[:,:,0]
        Y = circles[:,:,1]
        Z = circles[:,:,2]
        if frame:
            return ax.plot_wireframe(X, Y, Z, color=color, alpha=alpha)
        else:
            return ax.plot_surface(X, Y, Z, color=color, alpha=alpha)
    
    def plot(self, numPointsPerCircle=12, numCircles=2, color='black', alpha=0.5, frame=False):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandles = self.addToPlot(ax, numPointsPerCircle, numCircles, color, alpha, frame)
        ax.set_aspect('equal')

    
class Arc3D:
    def __init__(self, circleCenter, startPoint, startDir, theta):
        assert(norm(startPoint-circleCenter)>0)
        # verify orthogonality up to numerical stability
        assert(abs(np.dot(startDir, startPoint-circleCenter)) < 0.000001)
        
        self.circleCenter = circleCenter
        self.startPoint = startPoint
        self.startTangent = startDir / norm(startDir)
        self.theta = theta
        
        self.centerToStart = self.startPoint - self.circleCenter
        self.r = norm(self.centerToStart)
        self.startNormal = - self.centerToStart / self.r
                
        self.binormal = cross(self.startTangent, self.startNormal)
        self.rot = Rotation.from_rotvec(self.theta * self.binormal)
        self.centerToEnd = self.rot.apply(self.centerToStart)
        self.endPoint = self.circleCenter + self.centerToEnd
        self.endNormal = - self.centerToEnd / self.r
        self.endTangent = cross(self.endNormal, self.binormal)
    
    def interpolate(self, count=50):
        angle = np.linspace(0, self.theta, count).reshape(-1,1)
        u = self.r * np.cos(angle)
        v = self.r * np.sin(angle)
        
        # construct basis for circle plane
        uhat = -self.startNormal.reshape(1,3)
        vhat = cross(self.binormal, uhat).reshape(1,3)
        
        # 3d circle points
        return self.circleCenter + u @ uhat + v @ vhat

# add given reference frames to matplotlib figure ax with a 3d subplot
# pose is a matrix of SE3() objects
# returns the plot handles for the xHats, yHats, zHats, origins
def addPosesToPlot(Poses, ax, axisLength, xColor='r', yColor='b', zColor='g', oColors='black'):
    if Poses.shape == (4,4): # so it can plot a single frame
        Poses = np.array([Poses])
    
    ux, vx, wx = Poses[:,0:3,0].T # frame xhat coordinates
    uy, vy, wy = Poses[:,0:3,1].T # frame yhat coordinates
    uz, vz, wz = Poses[:,0:3,2].T # frame zhat coordinates
    ox, oy, oz = Poses[:,0:3,3].T # frame origin coordinates
    
    # https://matplotlib.org/stable/gallery/mplot3d/quiver3d.html
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.quiver.html#matplotlib.axes.Axes.quiver
    xHats = ax.quiver(ox, oy, oz, ux, vx, wx, length=axisLength, color=xColor) #plot xhat vectors
    yHats = ax.quiver(ox, oy, oz, uy, vy, wy, length=axisLength, color=yColor) #plot yhat vectors
    zHats = ax.quiver(ox, oy, oz, uz, vz, wz, length=axisLength, color=zColor) #plot zhat vectors
    origins = ax.scatter(ox, oy, oz, c=oColors)
    
    return (xHats, yHats, zHats, origins)


class Ball:
    # closed ball centered at self.c of radius self.r
    def __init__(self, center, radius):
        self.c = center
        self.r = radius
    
    def containsPoint(self, point):
        return norm(self.c - point) <= self.r
    
    def addToPlot(self, ax, color='black', alpha=0.1, frame=False):
        #https://www.tutorialspoint.com/plotting-a-3d-cube-a-sphere-and-a-vector-in-matplotlib
        u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
        x = self.c[0] + self.r*np.cos(u)*np.sin(v)
        y = self.c[1] + self.r*np.sin(u)*np.sin(v)
        z = self.c[2] + self.r*np.cos(v)
        if frame:
            return ax.plot_wireframe(x, y, z, color=color, alpha=alpha)
        else:
            return ax.plot_surface(x, y, z, color=color, alpha=alpha)
    
    def plot(self, color='black', alpha=1, frame=False):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandles = self.addToPlot(ax, color, alpha, frame)
        ax.set_aspect('equal')


"""
I'm pretty sure the below algorithm works but I haven't proven correctness yet.
The minimum bounding balls of balls is much more complicated for n balls, but
fortunately we only need to deal with the n=2 case.
For the n case, see:
Fischer, Kaspar, and Bernd Gartner. "The smallest enclosing ball of balls: 
combinatorial structure and algorithms." Proceedings of the nineteenth annual 
symposium on Computational geometry. 2003.
"""
# minimum bounding ball of 2 balls
def minBoundingBall(ball1, ball2):
    if ball1.r <= ball2.r:
        smaller, larger = ball1, ball2
    else:
        smaller, larger = ball2, ball1
    
    v = larger.c - smaller.c
    if norm(v) <= larger.r: # if larger contains smaller
        return larger
    else:
        vhat = v / norm(v)
        p = smaller.c - smaller.r*vhat
        q = larger.c + larger.r*vhat
        pq = q-p
        return Ball(p + pq/2, norm(pq)/2)

def distanceBetweenBalls(ball1, ball2):
    centerDistance = norm(ball2.c - ball1.c)
    return max(0, centerDistance - ball1.r - ball2.r)


"""
Returns the common normal from line 1 to line 2, input in point-direction form.
If the lines intersect, return the given value representing undefined.
"""
def commonNormal(point1, direction1, point2, direction2, undefined="undefined"):
    direction1 = direction1 / norm(direction1)
    direction2 = direction2 / norm(direction2)
    
    v12 = point2 - point1
    
    cp = cross(direction1, direction2)
    if norm(cp) > 0:
        nhat = cp / norm(cp)
    else: # z axes are parallel
        normalInPlane = cross(direction1, cross(v12, direction1))
        if norm(normalInPlane) == 0: #axes are coincident
            return undefined
        else:
            nhat = normalInPlane / norm(normalInPlane)
            
    projectionAlongNormal = dot(v12, nhat)
    if projectionAlongNormal > 0: # nhat is oriented correctly
        return nhat
    elif projectionAlongNormal < 0: # nhat is opposite orientation
        return -nhat
    else: # axes instersect
        return undefined


class Line:
    def __init__(self, point, direction):
        self.p = point
        self.dhat = direction / norm(direction)
    
class Plane:
    def __init__(self, point, normal):
        self.p = point
        self.nhat = normal / norm(normal)
        
    # sign is + if the point is on the +nhat side, 0 if on plane, - otherwise
    def signedDistanceToPoint(self, point, epsilon=0.00000001):
        projectionLength = dot(point - self.p, self.nhat)
        if abs(projectionLength) < epsilon: #for numerical stability
            return 0
        else:
            return projectionLength
    
    # return 0 if point is on the plane, +1 if on the +nhat side, -1 otherwise
    def sideOfPoint(self, point, epsilon=0.00000001):
        return np.sign(self.signedDistanceToPoint(point, epsilon))
    
    def containsPoint(self, point, epsilon=0.00000001):
        # plane contains point iff self.p - point is orthogonal to self.nhat
        return self.sideOfPoint(point, epsilon) == 0
    
    def intersectionWithLine(self, line, epsilon=0.00000001):
        if abs(dot(self.nhat, line.dhat)) < epsilon: 
            # line is either on the plane or parallel to it
            if self.containsPoint(line.p): 
                # line is on the plane
                return line
            else: 
                # line is parallel to the plane
                return "empty"
        else:
            # line intersects the plane at a single point
            # so the linear system constructed per wikipedia has point solution
            # https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
            # what they call cross(p01, p02) is just our nhat
            t = dot(self.nhat, line.p - self.p) / dot(-line.dhat, self.nhat)
            intersect = line.p + t * line.dhat
            assert(self.containsPoint(intersect))
            return intersect
        
    
    """ 
    Return an array containing the sides of the plane that the given line
    intersects with: 0 for on the plane, 1 for +nhat side, -1 for -nhat side
    """
    def sidesOfLine(self, line):
        intersection = self.intersectionWithLine(line)
        if intersection == 'empty':
            # line is parallel to plane, need to check on which side
            return [self.sideOfPoint(line.p)]
        elif type(intersection) == Line:
            # whole line is on plane
            return [0]
        else:
            # line crosses the plane at a point
            return [-1,0,1]
    
