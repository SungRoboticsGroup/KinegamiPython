# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 23:21:57 2023

@author: Daniel Feshbach

This is more or less a python version of a MATLAB function from:
Wei-Hsi Chen, Woohyeok Yang, Lucien Peach, Daniel E. Koditschek, and Cynthia R. Sung,
"Kinegami: Algorithmic Design of Compliant Kinematic Chains From Tubular Origami",
IEEE Transaction on Robotics 2022.

This particular part of the Kinegami code is based on a construction from:
Sikha Hota and Debasish Ghose,
"Optimal Geometrical Path in 3D with Curvature Constraint", 
IROS 2010.
"""

import numpy as np
from numpy import cross, dot, arctan2
import scipy
from scipy.linalg import null_space
from numpy.linalg import norm
import matplotlib.pyplot as plt
from spatialmath import SE3
import pyqtgraph.opengl as gl

import warnings

from geometryHelpers import *

from style import *

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
def pathErrorCSC(tDirMag, r, startPosition, startDir, endPosition, endDir, 
                 circle1sign, circle2sign):
    return PathCSC(tDirMag, r, startPosition, startDir, endPosition, endDir, 
                   circle1sign, circle2sign).error

def shortestCSC(r, startPosition, startDir, endPosition, endDir):
    startDir = startDir / norm(startDir)
    endDir = endDir / norm(endDir)
    
    t0 = endPosition - startPosition
    if norm(t0)==0:
        t0 = r * startDir
    t0mag = norm(t0)
    t0unit = t0 / t0mag
    t0DirMag = np.append(t0unit, t0mag)
    
    
    # solve for the path under all 4 sign combinations
    # solutions are values for tDirMag
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        solPP = scipy.optimize.fsolve(pathErrorCSC, x0=t0DirMag, args=
                                (r, startPosition, startDir, endPosition, endDir, 1, 1))
        solPM = scipy.optimize.fsolve(pathErrorCSC, x0=t0DirMag, args=
                                (r, startPosition, startDir, endPosition, endDir, 1, -1))
        solMP = scipy.optimize.fsolve(pathErrorCSC, x0=t0DirMag, args=
                                (r, startPosition, startDir, endPosition, endDir, -1, 1))
        solMM = scipy.optimize.fsolve(pathErrorCSC, x0=t0DirMag, args=
                                (r, startPosition, startDir, endPosition, endDir, -1, -1))
    
    # Construct the paths found based on the solutions found for tDirMag
    pathPP = PathCSC(solPP, r, startPosition, startDir, endPosition, endDir, 1, 1)
    pathPM = PathCSC(solPM, r, startPosition, startDir, endPosition, endDir, 1, -1)
    pathMP = PathCSC(solMP, r, startPosition, startDir, endPosition, endDir, -1, 1)
    pathMM = PathCSC(solMM, r, startPosition, startDir, endPosition, endDir, -1, -1)
    
    paths = [pathPP, pathPM, pathMP, pathMM]
    errorNorms = np.array([norm(path.error) for path in paths])
    lengths = np.array([path.length for path in paths])
    # exclude invalid paths from length-min selection
    lengths[errorNorms > 0.001*r] = np.inf
    
    return paths[np.argmin(lengths)]
        
# empty path of radius r at position p in direction d
def emptyCSC(r, p, d):
    d = d / norm(d)
    return PathCSC(np.append(d,0), r, p, d, p, d, 1, 1)

class PathCSC:
    """
    Representation for CSC Dubins paths based on the S section vector.
    Can't handle when one of the turn angles is exactly pi.
    Based on the construction in S. Hota and D. Ghose, 
    "Optimal Geometrical Path in 3D with Curvature Constraint", IROS 2010.
    """
                       
    """
    r is the turn radius (a scalar)
    startPosition, startDir, endPosition, endDir define the Dubins frames
                                       and are np arrays of shape (3,)
    tDirMag defines the t vector in direction-magnitude form, an np array of 
        shape (4,) with direction first 3 indices and magnitdue last
    circle1sign and circle2sign are each 1 or -1, and control the direction of 
        the circles (relative to startPosition, endPosition) in the planes defined by
        startDir/endDir and tDir.
    """
    def __init__(self, tDirMag, r, startPosition, startDir, endPosition, endDir, 
                 circle1sign, circle2sign):
        assert(r>=0 and abs(circle1sign)==1 and abs(circle2sign)==1)
        self.r = r
        self.startPosition = startPosition
        self.startDir = startDir
        self.endPosition = endPosition
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
        
        """ 
        The Hota and Ghose vector construction solves for where the path 
        departs circle1 and enters circle2 given the common tangent direction
        tUnit. 
        """
        self.circleNormal1 = unitNormalToBoth(self.tUnit, self.startDir) # normal to circle 1
        self.w1 = self.circle1sign * cross(self.startDir, self.circleNormal1) # unit vector from startPosition to circle center
        self.y1 = self.circle1sign * cross(self.tUnit, self.circleNormal1) # unit vector from circle 1 departure point to circle center
        self.circleCenter1 = self.startPosition + self.r * self.w1
        self.turn1end = self.circleCenter1 - self.r * self.y1
        
        self.circleNormal2 = unitNormalToBoth(self.tUnit, self.endDir) # normal to circle 2        
        self.w2 = circle2sign * cross(self.endDir, self.circleNormal2) # unit vector from endPosition to circle center
        self.y2 = circle2sign * cross(self.tUnit, self.circleNormal2) # unit vector from circle 2 entering point to circle center
        self.circleCenter2 = self.endPosition + self.r * self.w2
        self.turn2start = self.circleCenter2 - self.r * self.y2
        
        """
        If T is a valid solution, it should equal the vector from where 
        the path departs circle 1 to where it enters circle 2
        (according to the Hota and Ghose construction).
        """ 
        tError = self.t - (self.turn2start - self.turn1end)
        # We also need a part of the error to enforce that tDir becomes unit.
        self.error = np.append(tError, norm(tDirMag[0:3]) - 1)
        
        
        # Measure turning angles along each circle
        """
        self.theta1 = wrapAngle(signedAngle(self.startDir, self.tUnit, 
                                unitNormalToBoth(self.startDir, self.tUnit)))
        self.theta2 = wrapAngle(signedAngle(self.tUnit, self.endDir, 
                                unitNormalToBoth(self.tUnit, self.endDir)))
        """
        self.theta1 = wrapAngle(signedAngle(self.startDir, self.tUnit, 
                                cross(self.startDir, self.w1)))
        self.theta2 = wrapAngle(signedAngle(self.tUnit, self.endDir, 
                                cross(self.endDir, self.w2)))
        self.length = self.r*self.theta1 + self.tMag + self.r*self.theta2

    def newPathTransformedBy(self, Transformation : SE3):
        new_turn1end = Transformation * self.turn1end
        new_turn1endPlusTunit = Transformation * (self.turn1end + self.tUnit)
        new_tUnit = new_turn1endPlusTunit - new_turn1end
        #new_tUnit = Transformation * self.tUnit
        new_tDirMag = np.vstack((new_tUnit, [self.tMag])).flatten()
        
        new_startPosition = (Transformation * self.startPosition).flatten()
        new_startPosPlusStartDir = (Transformation*(self.startPosition + self.startDir)).flatten()
        new_startDir = new_startPosPlusStartDir - new_startPosition
        
        new_endPosition = (Transformation * self.endPosition).flatten()
        new_endPosPlusEndDir = (Transformation*(self.endPosition + self.endDir)).flatten()
        new_endDir = new_endPosPlusEndDir - new_endPosition
        
        return PathCSC(new_tDirMag, self.r, new_startPosition, new_startDir,
                       new_endPosition, new_endDir, 
                       self.circle1sign, self.circle2sign)
    
    def __str__(self):
        return str(np.append(self.tUnit, self.tMag))
    
    def __repr__(self):
        tDirMag = np.append(self.tUnit, self.tMag)
        return "PathCSC(tDirMag="+repr(tDirMag)+\
                        ", r="+repr(self.r)+\
                        ", startPosition="+repr(self.startPosition)+\
                        ", startDir="+repr(self.startDir)+\
                        ", endPosition="+repr(self.endPosition)+\
                        ", endDir="+repr(self.endDir)+\
                        ", circle1sign="+repr(self.circle1sign)+\
                        ", circle2sign="+repr(self.circle2sign)+")"
    
    # add to existing GLViewWidget ax
    def add(self, widget, showCircles=True, showPoses=True, 
                  startColor=pathStartColor, endColor=pathEndColor, pathColor=pathColorDefault,
                  cscBoundaryMarker='*', showTunit=False):
        if showPoses:
            # start pose
            lineStart = gl.GLLinePlotItem(pos=(self.startPosition, self.startPosition + self.startDir), color=startColor, width=2) 
            widget.addItem(lineStart)

            # end pose
            lineEnd = gl.GLLinePlotItem(pos=(self.endPosition, self.endPosition + self.endDir), color=endColor, width=2) 
            widget.addItem(lineEnd)

        if showCircles:
            # start circle
            c1x, c1y, c1z = Circle3D(self.r, 
                            self.circleCenter1, self.circleNormal1).interpolate().T

            pointsStart = np.vstack([c1x, c1y, c1z]).transpose()
            circleStart = gl.GLLinePlotItem(pos=pointsStart, color=startColor, width=1)
            widget.addItem(circleStart)

            # end circle
            c2x, c2y, c2z = Circle3D(self.r, 
                            self.circleCenter2, self.circleNormal2).interpolate().T
            pointsEnd = np.vstack([c2x, c2y, c2z]).transpose()
            circleEnd = gl.GLLinePlotItem(pos=pointsEnd, color=endColor, width=1)
            widget.addItem(circleEnd)
        
        # path arcs (C components)
        a1x, a1y, a1z = Arc3D(self.circleCenter1, 
                self.startPosition, self.startDir, self.theta1).interpolate().T
        arcPoints1 = np.vstack([a1x, a1y, a1z]).transpose()
        arc1 = gl.GLLinePlotItem(pos=arcPoints1, color=pathColor, width=2)
        widget.addItem(arc1)

        a2x, a2y, a2z = Arc3D(self.circleCenter2, 
                self.turn2start, self.tUnit, self.theta2).interpolate().T
        arcPoints2 = np.vstack([a2x, a2y, a2z]).transpose()
        arc2 = gl.GLLinePlotItem(pos=arcPoints2, color=pathColor, width=2)
        widget.addItem(arc2)
        
        # path S component
        sx,sy,sz = np.array([self.turn1end, self.turn2start]).T
        pathSPoints = np.vstack([sx, sy, sz]).transpose()
        pathS = gl.GLLinePlotItem(pos=pathSPoints, color=pathColor, width=2)
        widget.addItem(pathS)

        point1 = gl.GLScatterPlotItem(pos=pathSPoints[0], color=pathPointColor, size=10)
        widget.addItem(point1)
        point2 = gl.GLScatterPlotItem(pos=pathSPoints[1], color=pathPointColor, size=10)
        widget.addItem(point2)

        if showTunit:
            tUnit = gl.GLLinePlotItem(pos=(self.turn1end, self.turn1end + self.tUnit), color=tUnitColor, width=2) 
            widget.addItem(tUnit)
    
    def addToWidget(self, widget, showCircles=True, showPoses=True, 
                  startColor=pathStartColor, endColor=pathEndColor, pathColor=pathWidgetColor,
                  cscBoundaryMarker='*', showTunit=False, block=blockDefault):

        #make a new plot and then call add to plot, then display that plot
        self.add(widget, showCircles, showPoses, startColor, endColor, 
                       pathColor, cscBoundaryMarker, showTunit)
        ax.set_aspect('equal')
        ax.legend()
        plt.show(block=block)