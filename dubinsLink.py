# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 15:42:09 2023

@author: Daniel Feshbach
"""
import numpy as np
from numpy import cross, dot, arctan2
from numpy.linalg import norm
from spatialmath import SO3, SE3
from dubinsPath import *
from geometryHelpers import *
import matplotlib.pyplot as plt
import tubularOrigami
from tubularOrigami import TubularPattern, TubeFittingPattern, \
                            ElbowFittingPattern, TwistFittingPattern


class LinkCSC:
    def __init__(self, r : float, StartDubinsFrame : SE3, EndDubinsFrame : SE3,
                 splitLongElbowsInto : int = 2, EPSILON : float = 0.0001):
        assert(r>0)
        assert(splitLongElbowsInto >= 1)
        self.r = r
        self.splitLongElbowsInto = splitLongElbowsInto
        self.EPSILON = EPSILON
        self.StartDubinsFrame = StartDubinsFrame
        self.EndDubinsFrame = EndDubinsFrame
        self.path = shortestCSC(r, self.StartDubinsFrame.t, self.StartDubinsFrame.R[:,0],
                                self.EndDubinsFrame.t, self.EndDubinsFrame.R[:,0])
        assert(norm(self.path.error) < self.EPSILON)
        assert(self.path.theta1 >= 0)
        assert(self.path.theta1 <= np.pi)
        assert(self.path.theta2 >= 0)
        assert(self.path.theta2 <= np.pi)
        self.rot1AxisDir = np.cross(self.StartDubinsFrame.R[:,0], self.path.w1)
        self.rot1AxisAngle = signedAngle(self.StartDubinsFrame.R[:,1],
                                          self.rot1AxisDir,
                                          self.StartDubinsFrame.R[:,0])
        self.rot2AxisDir = np.cross(self.EndDubinsFrame.R[:,0], self.path.w2)
        self.rot2AxisAngle = signedAngle(self.EndDubinsFrame.R[:,1],
                                          self.rot2AxisDir,
                                          self.EndDubinsFrame.R[:,0])
        
        
        if self.path.theta1 > self.EPSILON:
            self.elbow1 = CompoundElbow(self.r, self.StartDubinsFrame, 
                                       self.path.theta1, self.rot1AxisAngle, 
                                       self.splitLongElbowsInto, self.EPSILON)
            self.Elbow1EndFrame = self.elbow1.EndFrame
            self.elbow1BoundingBall = self.elbow1.boundingBall()
        else:
            self.Elbow1EndFrame = self.StartDubinsFrame
            self.elbow1BoundingBall = Ball(self.StartDubinsFrame.t, self.r)
        assert(norm(self.Elbow1EndFrame.R[:,0] - self.path.tUnit) < self.EPSILON) # verify alignment
        assert(norm(self.Elbow1EndFrame.t - self.path.turn1end) < self.EPSILON)
        
        if self.path.tMag > self.EPSILON:
            self.cylinder = Cylinder(self.r, self.Elbow1EndFrame.t, 
                            self.Elbow1EndFrame.R[:,0], self.path.tMag)
        
        if self.path.theta2 > self.EPSILON:
            Elbow2StartOrientation = SO3.AngleAxis(-self.path.theta2, self.rot2AxisDir) * self.EndDubinsFrame.R 
            self.Elbow2StartFrame = SE3.Rt(Elbow2StartOrientation, self.path.turn2start)
            assert(norm(self.Elbow2StartFrame.t - (self.Elbow1EndFrame * SE3.Tx(self.path.tMag)).t) < self.EPSILON)
            assert(norm(self.Elbow2StartFrame.R[:,0] - self.path.tUnit) < self.EPSILON) # verify alignment
            
            self.elbow2 = CompoundElbow(self.r, self.Elbow2StartFrame, 
                                       self.path.theta2, self.rot2AxisAngle, 
                                       self.splitLongElbowsInto, self.EPSILON)
            assert(norm(self.elbow2.EndFrame - self.EndDubinsFrame) < self.EPSILON)
            self.elbow2BoundingBall = self.elbow2.boundingBall()
        else:
            self.Elbow2StartFrame = self.EndDubinsFrame
            self.elbow2BoundingBall = Ball(self.EndDubinsFrame.t, self.r)
        
        
    def addToPlot(self, ax, numSides : int = 32, color : str = 'black', 
                  alpha : float = 0.5, wireFrame : bool = False, 
                  showFrames : bool = False, showPath : bool = True, 
                  showPathCircles : bool = False, showBoundary : bool = True,
                  showElbowBoundingBalls : bool = False):
        allElbowHandleSets = []
        if showBoundary:
            if self.path.theta1 > self.EPSILON:
                elbow1HandleSets = self.elbow1.addToPlot(ax, numSides, color, 
                                                alpha, wireFrame, showFrames,
                                                showElbowBoundingBalls)
                if showFrames:
                    allElbowHandleSets += elbow1HandleSets
            
            if self.path.tMag > self.EPSILON:
                self.cylinder.addToPlot(ax, numSides, color, alpha, 
                                                      wireFrame)
            
            if self.path.theta2 > self.EPSILON:
                elbow2HandleSets = self.elbow2.addToPlot(ax, numSides, color, 
                                                    alpha, wireFrame, showFrames,
                                                    showElbowBoundingBalls)
                if showFrames:
                    allElbowHandleSets += elbow2HandleSets
        
        if showPath:
            self.path.addToPlot(ax, showCircles=showPathCircles, showPoses=showFrames)
        
        return allElbowHandleSets
    
    
    def plot(self, numSides : int = 32, color : str = 'black', 
                  alpha : float = 0.5, wireFrame : bool = False, 
                  showFrames : bool = False, showPath : bool = True, 
                  showPathCircles : bool = False, showBoundary : bool = True,
                  showElbowBoundingBalls : bool = False):
        ax = plt.figure().add_subplot(projection='3d')
        allElbowHandleSets = self.addToPlot(ax, numSides, color, alpha, wireFrame, 
                                     showFrames, showPath, showPathCircles,
                                     showBoundary, showElbowBoundingBalls)
        ax.set_aspect('equal')
    
    def creasePattern(self, numSides : int, twistPortion : float = 0.2) -> TubularPattern:
        assert(numSides >= 4 and numSides%2==0)
        assert(twistPortion > 0)
        
        composed = TubularPattern(numSides, self.r)
        
        if self.path.theta1 > self.EPSILON:
            #w1 is the unit vector from startPosition to circle1center
            if self.path.theta1 > np.pi/self.splitLongElbowsInto:
                elbow1PartPattern = ElbowFittingPattern(numSides, self.r, 
                                    bendingAngle=self.path.theta1/self.splitLongElbowsInto, 
                                    rotationalAxisAngle=self.rot1AxisAngle)
                for i in range(self.splitLongElbowsInto):
                    composed.append(elbow1PartPattern)
            else:
                elbow1Pattern = ElbowFittingPattern(numSides, self.r, 
                                    bendingAngle=self.path.theta1, 
                                    rotationalAxisAngle=self.rot1AxisAngle)
                composed.append(elbow1Pattern)
        
        twistAngle = signedAngle(self.Elbow1EndFrame.R[:,1], 
                                 self.Elbow2StartFrame.R[:,1], 
                                 self.path.tUnit)
        twistLen = 0
        if abs(twistAngle) > self.EPSILON:
            twistLen = twistPortion * self.path.tMag
            twistPattern = TwistFittingPattern(numSides, self.r, twistAngle, twistLen)
            composed.append(twistPattern)
        
        tubeLen = self.path.tMag - twistLen
        tubePattern = TubeFittingPattern(numSides, self.r, tubeLen)
        composed.append(tubePattern)
        
        if self.path.theta2 > self.EPSILON:
            #w1 is the unit vector from startPosition to circle1center
            if self.path.theta2 > np.pi/self.splitLongElbowsInto:
                elbow2PartPattern = ElbowFittingPattern(numSides, self.r, 
                                    bendingAngle=self.path.theta2/self.splitLongElbowsInto, 
                                    rotationalAxisAngle=self.rot2AxisAngle)
                for i in range(self.splitLongElbowsInto):
                    composed.append(elbow2PartPattern)
            else:
                elbow2Pattern = ElbowFittingPattern(numSides, self.r, 
                                    bendingAngle=self.path.theta2, 
                                    rotationalAxisAngle=self.rot2AxisAngle)
                composed.append(elbow2Pattern)
        
        return composed
