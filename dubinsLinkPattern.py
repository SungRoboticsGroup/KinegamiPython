# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 14:39:26 2023

@author: Daniel Feshbach
"""
import tubularOrigami
from tubularOrigami import TubularPattern, TubeFittingPattern, \
                            ElbowFittingPattern, TwistFittingPattern
import dubinsPath
from dubinsPath import *

# applies to DUBINS FRAMES, where the dubins direction is column 0 (ahat)
def elbowTranformation(rotAxis : np.ndarray,
                       bendingAngle : float,
                       r : float,
                       EPSILON : float = 0.0001) -> SE3:
    if bendingAngle < EPSILON:
        return SE3()
    else:
        Forward = SE3.Tx(r*np.tan(bendingAngle/2))
        Rotate = SE3.AngleAxis(bendingAngle, rotAxis)
        return Forward @ Rotate @ Forward
    
def dubinsLinkPattern(numSides : int,
                      r : float,
                      ProximalDubinsFrame : SE3,
                      DistalDubinsFrame : SE3,
                      path : PathCSC = None, 
                      splitLongElbowsInto : int = 2,
                      twistPortion : float = 0.2,
                      EPSILON : float = 0.001) -> TubularPattern:
    
    assert(splitLongElbowsInto >= 1)
    proximalPathDir = ProximalDubinsFrame.R[:,0]
    bProximal = ProximalDubinsFrame.R[:,1]
    distalPathDir = DistalDubinsFrame.R[:,0]
    bDistal = DistalDubinsFrame.R[:,1]
    
    if path is None:
        path = shortestCSC(r, ProximalDubinsFrame.t, proximalPathDir,
                              DistalDubinsFrame.t, distalPathDir)
    else:
        assert(path.r == r)
        assert(norm(path.startPosition - ProximalDubinsFrame.t) < EPSILON)
        assert(norm(path.endPosition - DistalDubinsFrame.t) < EPSILON)
        assert(norm(path.startDir - ProximalDubinsFrame.R[:,0]) < EPSILON)
        assert(norm(path.endDir - DistalDubinsFrame.R[:,0]) < EPSILON)
    
    assert(path.theta1 >= 0)
    assert(path.theta1 <= np.pi)
    assert(path.theta2 >= 0)
    assert(path.theta2 <= np.pi)
        
    composed = TubularPattern(numSides, r)
    
    rot1AxisDir = np.cross(proximalPathDir, path.w1)
    Elbow1Transform = elbowTranformation(rot1AxisDir, path.theta1, r, EPSILON)
    Elbow1EndFrame = Elbow1Transform @ ProximalDubinsFrame
    assert(np.linalg.norm(Elbow1EndFrame.R[:,0] - path.tUnit) < EPSILON) # verify alignment
    
    rot2AxisDir = np.cross(distalPathDir, path.w2)
    Elbow2Transform = elbowTranformation(rot2AxisDir, path.theta2, r, EPSILON)
    Elbow2StartFrame = Elbow2Transform.inv() @ DistalDubinsFrame
    assert(np.linalg.norm(Elbow2StartFrame.R[:,0] - path.tUnit) < EPSILON) # verify alignment
    
    if path.theta1 > EPSILON:
        #w1 is the unit vector from startPosition to circle1center
        elbow1AxisAngle = signedAngle(bProximal, rot1AxisDir, proximalPathDir)
        if path.theta1 > np.pi/splitLongElbowsInto:
            elbow1PartPattern = ElbowFittingPattern(numSides, r, 
                                bendingAngle=path.theta1/splitLongElbowsInto, 
                                rotationalAxisAngle=elbow1AxisAngle)
            for i in range(splitLongElbowsInto):
                composed.append(elbow1PartPattern)
        else:
            elbow1Pattern = ElbowFittingPattern(numSides, r, 
                                bendingAngle=path.theta1, 
                                rotationalAxisAngle=elbow1AxisAngle)
            composed.append(elbow1Pattern)
    
    twistAngle = signedAngle(Elbow1EndFrame.R[:,1], Elbow2StartFrame.R[:,1], 
                             path.tUnit)
    twistLen = 0
    if abs(twistAngle) > EPSILON:
        twistLen = twistPortion * path.tMag
        twistPattern = TwistFittingPattern(numSides, r, twistAngle, twistLen)
        composed.append(twistPattern)
    
    tubeLen = path.tMag - twistLen
    tubePattern = TubeFittingPattern(numSides, r, tubeLen)
    composed.append(tubePattern)
    
    if path.theta2 > EPSILON:
        #w1 is the unit vector from startPosition to circle1center
        elbow2AxisAngle = signedAngle(bDistal, rot2AxisDir, distalPathDir)
        if path.theta2 > np.pi/splitLongElbowsInto:
            elbow2PartPattern = ElbowFittingPattern(numSides, r, 
                                bendingAngle=path.theta2/splitLongElbowsInto, 
                                rotationalAxisAngle=elbow2AxisAngle)
            for i in range(splitLongElbowsInto):
                composed.append(elbow2PartPattern)
        else:
            elbow2Pattern = ElbowFittingPattern(numSides, r, 
                                bendingAngle=path.theta2, 
                                rotationalAxisAngle=elbow2AxisAngle)
            composed.append(elbow2Pattern)
    
    return composed
