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
from spatialmath import SE3
import matplotlib.pyplot as plt

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
Signed angle (radians) from vector a to vector b, around the normal vector n.
All inputs should be numpy arrays of shape (3,)
"""
def signedAngle(a, b, n):
    a = a / norm(a)
    b = b / norm(b)
    if norm(n)>0:
        n = n / norm(n)
        
    return arctan2(dot(cross(a,b),n), dot(a,b))

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
    
    
    
    