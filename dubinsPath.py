# -*- cendPosing: utf-8 -*-
"""
Created on Thu Jun 15 23:21:57 2023

@author: dfesh
"""

import numpy as np
from numpy import cross, dot, arctan2
import scipy
from scipy.linalg import null_space
from numpy.linalg import norm
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

"""
Given tDirMag as a 4D vector containing [tDir, tMag] representing t,
and the other Dubins path parameters, construct the PathCSC and return
its error, a 4D vector composed of:
    
    First 3 indices: 
        t - (turn2start - turn1end),
            correctness of T based on the Hota and Ghose 2010 construction
    
    Last index:
        norm(tDir) - 1
            enforcing that tDir should become a unit vector 

This function exists to be input in scipy.optimize.fsolve
https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.fsolve.html
so it should be of the form func(x, *args) where we solve func(x)=0
for x given other parameters in args. 

In this case, x is tDirMag (we're solving for t and representing it as 
direction-magnitude to ensure tDir becomes the common tangent between circles 
even if t=[0,0,0]) and args are the Dubins start/end frames and the signs 
controlling the circle directions.
"""
def pathErrorCSC(tDirMag, r, startPos, startDir, endPos, endDir, 
                 circle1sign, circle2sign):
    return PathCSC(tDirMag, r, startPos, startDir, endPos, endDir, 
                   circle1sign, circle2sign).error

def shortestCSC(r, startPos, startDir, endPos, endDir):
    t0 = endPos - startPos
    if norm(t0)==0:
        t0 = r * startDir
    t0mag = norm(t0)
    t0unit = t0 / t0mag
    t0DirMag = np.append(t0unit, t0mag)
    
    
    # solve for the path under all 4 sign combinations
    # solutions are values for tDirMag
    solPP = scipy.optimize.fsolve(pathErrorCSC, x0=t0DirMag, args=
                            (r, startPos, startDir, endPos, endDir, 1, 1))
    solPM = scipy.optimize.fsolve(pathErrorCSC, x0=t0DirMag, args=
                            (r, startPos, startDir, endPos, endDir, 1, -1))
    solMP = scipy.optimize.fsolve(pathErrorCSC, x0=t0DirMag, args=
                            (r, startPos, startDir, endPos, endDir, -1, 1))
    solMM = scipy.optimize.fsolve(pathErrorCSC, x0=t0DirMag, args=
                            (r, startPos, startDir, endPos, endDir, -1, -1))
    
    # Construct the paths found based on the solutions found for tDirMag
    pathPP = PathCSC(solPP, r, startPos, startDir, endPos, endDir, 1, 1)
    pathPM = PathCSC(solPM, r, startPos, startDir, endPos, endDir, 1, -1)
    pathMP = PathCSC(solMP, r, startPos, startDir, endPos, endDir, -1, 1)
    pathMM = PathCSC(solMM, r, startPos, startDir, endPos, endDir, -1, -1)
    
    paths = [pathPP, pathPM, pathMP, pathMM]
    errorNorms = np.array([norm(path.error) for path in paths])
    lengths = np.array([path.length for path in paths])
    print(errorNorms)
    print(lengths)
    # exclude invalid paths from length-min selection
    lengths[errorNorms > 0.001*r] = np.inf
    
    return paths[np.argmin(lengths)]
    
    
    
# wrap angles to [0,2pi)
def wrapAngle(angle):
    return angle % (2*np.pi)

class Circle3D:
    def __init__(self, radius, center, normal):
        assert(norm(normal)>0)
        self.r = radius
        self.c = center
        self.n = normal / norm(normal)
    
    def points(self, count=50):
        angle = np.linspace(0, 2*np.pi, count).reshape(-1,1)
        u = self.r * np.cos(angle)
        v = self.r * np.sin(angle)
        
        # construct basis for circle plane
        uhat = unitNormalToBoth(self.n, self.n).reshape(1,3)
        vhat = cross(self.n, uhat).reshape(1,3)
        
        # 3d circle points
        return self.c + u @ uhat + v @ vhat
        # TODO: test this (especially vectorization dimensions)
        
        
        


class PathCSC:
    """
    Representation for CSC Dubins paths based on the S section vector.
    Can't handle when one of the turn angles is exactly pi.
    Based on the construction in S. Hota and D. Ghose, 
    "Optimal Geometrical Path in 3D with Curvature Constraint", IROS 2010.
    """
    
    """
    r is the turn radius (a scalar)
    startPos, startDir, endPos, endDir define the Dubins frames
                                       and are np arrays of shape (3,)
    tDirMag defines the t vector in direction-magnitude form, an np array of 
        shape (4,) with direction first 3 indices and magnitdue last
    circle1sign and circle2sign are each 1 or -1, and control the direction of 
        the circles (relative to startPos, endPos) in the planes defined by
        startDir/endDir and tDir.
    """
    def __init__(self, tDirMag, r, startPos, startDir, endPos, endDir, 
                 circle1sign, circle2sign):
        assert(r>=0 and abs(circle1sign)==1 and abs(circle2sign)==1)
        self.r = r
        self.startPos = startPos
        self.startDir = startDir
        self.endPos = endPos
        self.endDir = endDir
        self.circle1sign = circle1sign  # direction of circle 1 in plane of startDir, tUnit
        self.circle2sign = circle2sign  # direction of circle 2 in plane of endDir, tUnit
        
        """
        The straight component of the path, T, is tangent to both circles
        in a direction agreeing with startDir and endDir.
        The optimization finds the path by solving for this component
        in direction-magnitude representation.
        """
        self.tUnit = tDirMag[0:3] / norm(tDirMag[0:3]) # common unit tangent between circles
        self.tMag = abs(tDirMag[3])   # S section length (along tUnit)
        self.t = self.tMag * self.tUnit
        
        # Measure turning angles along each circle
        self.theta1 = wrapAngle(signedAngle(self.startDir, self.tUnit, 
                                unitNormalToBoth(self.startDir, self.tUnit)))
        self.theta2 = wrapAngle(signedAngle(self.tUnit, self.endDir, 
                                unitNormalToBoth(self.tUnit, self.endDir)))
        
        """ 
        The Hota and Ghose vector construction solves for where the path 
        departs circle1 and enters circle2 given the common tangent direction
        tUnit. 
        """
        self.circleNormal1 = unitNormalToBoth(self.tUnit, self.startDir) # normal to circle 1
        self.w1 = self.circle1sign * cross(self.startDir, self.circleNormal1) # unit vector from startPos to circle center
        self.y1 = self.circle1sign * cross(self.tUnit, self.circleNormal1) # unit vector from circle 1 departure point to circle center
        self.circleCenter1 = self.startPos + self.r * self.w1
        self.turn1end = self.circleCenter1 - self.r * self.y1
        
        self.circleNormal2 = unitNormalToBoth(self.tUnit, self.endDir) # normal to circle 2        
        self.w2 = circle2sign * cross(self.endDir, self.circleNormal2) # unit vector from startPos to circle center
        self.y2 = circle2sign * cross(self.tUnit, self.circleNormal2) # unit vector from circle 2 entering point to circle center
        self.circleCenter2 = self.endPos + self.r * self.w2
        self.turn2start = self.circleCenter2 - self.r * self.y2
        
        """
        If T is a valid solution, it should equal the vector from where 
        the path departs circle 1 to where it enters circle 2
        (according to the Hota and Ghose construction).
        """ 
        tError = self.t - (self.turn2start - self.turn1end)
        # We also need a part of the error to enforce that tDir becomes unit.
        self.error = np.append(tError, norm(tDirMag[0:3]) - 1)
        
        self.length = self.r*self.theta1 + self.tMag + self.r*self.theta2
    
    def __str__(self):
        return str(np.append(self.tUnit, self.tMag))
    
    def __repr__(self):
        tDirMag = np.append(self.tUnit, self.tMag)
        return "PathCSC(tDirMag="+str(tDirMag)+\
                        ", r="+str(self.r)+\
                        ", startPos="+str(self.startPos)+\
                        ", startDir="+str(self.startDir)+\
                        ", endPos="+str(self.endPos)+\
                        ", endDir="+str(self.endDir)+\
                        ", circle1sign="+str(self.circle1sign)+\
                        ", circle2sign="+str(self.circle2sign)+")"
    
    # add to existing matplotlib axis ax
    def addToPlot(self, ax, startColor='red', endColor='blue', pathColor='green'):
        # start pose
        x1,y1,z1 = self.startPos
        u1,v1,w1 = self.startDir
        ax.quiver(x1,y1,z1,u1,v1,w1,
                  length=self.r, color=startColor, label='start')
        
        # end pose
        x2,y2,z2 = self.endPos
        u2,v2,w2 = self.endDir
        ax.quiver(x2,y2,z2,u2,v2,w2,
                  length=self.r, color=endColor, label='end')
        
        # circles
        c1x, c1y, c1z = Circle3D(self.r, 
                            self.circleCenter1, self.circleNormal1).points().T
        ax.plot(c1x, c1y, c1z, color = startColor)
        c2x, c2y, c2z = Circle3D(self.r, 
                            self.circleCenter2, self.circleNormal2).points().T
        ax.plot(c2x, c2y, c2z, color = endColor)
        
        # path S component
        sx,sy,sz = np.array([self.turn1end, self.turn2start]).T
        ax.plot(sx, sy, sz, color = pathColor, marker='*')
        
    
    def plot(self, startColor='red', endColor='blue', pathColor='green'):
        ax = plt.figure().add_subplot(projection='3d')
        self.addToPlot(ax, startColor, endColor, pathColor)
        ax.set_aspect('equal')
        ax.legend()
        plt.show()
        
        
        
        
        
        
        
        