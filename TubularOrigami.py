# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 17:11:58 2023

@author: dfesh

Implementing the catalog of tubular origami crease patterns from the original 
Kinegami paper (Chen et al. 2022).

These functions construct 2D crease patterns which fold into tubular modules
of radius r whose ends are identical regular polygons where the number of sides
(called numSides here and n_s in the original paper) is even.
Let x coordinates be around the tube and y coordinates be along it, so the
tube wraparound is represented by modular wraparound in x. 
The base polygon's side length can be computed from r and numSides.
The proximal base polygon's vertices have y=0 in the crease pattern, and the
distal base has constant y>0. 
The proximal base vertices have x coordinates
    0, baseSideLength, 2*baseSideLength, ..., (numSides-1)*baseSideLength
and x coordinates wrap modulo numSides*baseSideLength.
We store the pattern as a graph where vertices are points and edges are
partitioned into mountain and valley lists (each encoded by vertex pairs).
Since edges are within their panel of width baseSideLength, we can detect
edges with ostensible x component > baseSideLength as those that wrap around.
?????????????????????????
When outputting to a dxf pattern, we can find all edges in x range
[0, baseSideLength] with ???????????????????????


Can I use integer coordinates for the proximal base, and then just scale
based on r at the end????
Should we just manually maintain a marking of which edges are wraparound?
Or manually construct the duplicate at the end?
"""
import numpy as np
from abc import ABC, abstractmethod
import ezdxf

class TubularPattern(ABC):
    def __init__(self, numSides, r, proximalMarker=[0,0]):
        assert(r>0)
        assert(numSides%2 == 0)
        self.numSides = numSides                                   #n_s
        self.polygonInnerAngle = np.pi * (numSides-2)/(2*numSides) #\delta
        self.r = r
        self.baseSideLength = 2*r*np.cos(self.polygonInnerAngle)
        self.proximalMarker = np.array(proximalMarker).reshape((1,2))
        proximalBaseX = np.arange(self.numSides) + 0.5
        nzeros = np.zeros(self.numSides)
        relative = np.vstack((proximalBaseX, nzeros)).T
        self.proximalBase = self.proximalMarker + relative
        self.width = self.width = self.baseSideLength * self.numSides
    
    def wrapVertices(self):
        self.Vertices[:,0] %= self.width
    
    def MountainFolds(self):
        return self.Vertices[self.MountainEdges]
    
    def ValleyFolds(self):
        return self.Vertices[self.ValleyEdges]
    
    def saveDXF(self, name):
        doc = ezdxf.new()

        # add new entities to the modelspace
        msp = doc.modelspace()
        doc.layers.add(name="Mountain", color=5)
        doc.layers.add(name="Valley", color=1)
        doc.layers.add(name="Boundary", color=2)
        
        ymin = np.min(self.Vertices[:,1])
        ymax = np.max(self.Vertices[:,1])
        xmin = 0
        xmax = self.width
        
        msp.add_line((xmin,ymin),(xmin,ymax), dxfattribs={"layer": "Boundary"})
        msp.add_line((xmin,ymin),(xmax,ymin), dxfattribs={"layer": "Boundary"})
        msp.add_line((xmax,ymax),(xmin,ymax), dxfattribs={"layer": "Boundary"})
        msp.add_line((xmax,ymax),(xmax,ymin), dxfattribs={"layer": "Boundary"})
        
        if self.MountainEdges.shape[1] > 0:  
            MF = self.MountainFolds()
            WrapAround = abs(MF[:,1,0]-MF[:,0,0]) > self.baseSideLength
            for seg in MF[np.logical_not(WrapAround)]:
                msp.add_line(seg[0,:],seg[1,:],dxfattribs={"layer":"Mountain"})
            for seg in MF[WrapAround]:
                seg[seg[:,0] <= self.baseSideLength, 0] += self.width
                msp.add_line(seg[0,:],seg[1,:],dxfattribs={"layer":"Mountain"})
        
        if self.ValleyEdges.shape[1] > 0:
            VF = self.ValleyFolds()
            WrapAround = abs(VF[:,1,0]-VF[:,0,0]) > self.baseSideLength
            for seg in VF[np.logical_not(WrapAround)]:
                msp.add_line(seg[0,:],seg[1,:],dxfattribs={"layer":"Valley"})
            for seg in VF[WrapAround]:
                seg[seg[:,0] <= self.baseSideLength, 0] += self.width
                msp.add_line(seg[0,:],seg[1,:],dxfattribs={"layer":"Valley"})
        
        doc.saveas(name+".dxf")

class Tube(TubularPattern):
    def __init__(self, numSides, r, height, proximalMarker=[0,0]):
        super().__init__(numSides, r, proximalMarker)
        self.height = height
        self.distalBase = self.proximalBase + [0, self.height]
        self.Vertices = np.vstack((self.proximalBase, self.distalBase))
        proxInds = np.arange(self.numSides)
        distInds = self.numSides + proxInds
        self.MountainEdges = np.vstack((proxInds, distInds)).T
        self.ValleyEdges = np.array([[]]) # What shape to give this???
        self.distalMarker = self.proximalMarker + [0, self.height]
        self.distalMarker[0] %= self.width
        self.wrapVertices()
        
        
class Twist(TubularPattern):
    def __init__(self, numSides, r, twistAngle, moduleHeight, 
                 proximalMarker=[0,0]):
        super().__init__(numSides, r, proximalMarker)
        twistAngle = twistAngle % (2*np.pi)
        """ 
        Twist in increments of 2pi/numSides is be accomplished by rotating the
        next modules itself (shift crease pattern columns modularly).
        So the twist module only needs to instantiate the remaining portion:
        """
        baseTwist = twistAngle % (2*np.pi / self.numSides) #\overline{\alpha}
        baseXshift = (self.baseSideLength/2) * ( 1 - np.cos(baseTwist) + 
              (1/np.tan(np.pi/numSides))*np.sin(baseTwist) )
        patternHeight = np.sqrt(moduleHeight**2 + 
                        (self.baseSideLength * 
                         (1/np.sin(np.pi/numSides)) *
                         np.sin((np.pi/numSides) - (baseTwist/2)) *
                         np.sin(baseTwist/2) )**2 )
        self.distalBase = self.proximalBase + [baseXshift, patternHeight]
        self.distalBase[:,0] %= numSides * self.baseSideLength
        self.Vertices = np.vstack((self.proximalBase, self.distalBase))
        
        proxInds = np.arange(self.numSides)
        distInds = self.numSides + proxInds
        rightDiagonals = np.vstack((proxInds, distInds)).T
        leftDiagonals = np.vstack((proxInds, np.roll(distInds,1))).T
        self.MountainEdges = np.vstack((rightDiagonals, leftDiagonals))
        
        proximalValleys = np.vstack((proxInds, np.roll(proxInds,1))).T
        distalValleys = np.vstack((distInds, np.roll(distInds,1))).T
        self.ValleyEdges = np.vstack((proximalValleys,distalValleys))
        
        referenceXshift = (np.floor(numSides * twistAngle / (2*np.pi)) * \
                            self.baseSideLength) + baseXshift
        self.distalMarker = self.proximalMarker + \
                            np.array([referenceXshift, patternHeight])
        self.distalMarker[0] %= self.width
        self.wrapVertices()
        
    
