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