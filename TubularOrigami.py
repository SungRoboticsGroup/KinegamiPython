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

TODO: Cut off the segments after they leave the cut boundary
"""
import numpy as np
import matplotlib.pyplot as plt
import ezdxf
from ezdxf.addons.drawing import RenderContext, Frontend
from ezdxf.addons.drawing.matplotlib import MatplotlibBackend
from ezdxf.addons.drawing.config import Configuration
from ezdxf.addons.drawing.properties import LayoutProperties


class TubularPattern():
    def __init__(self, numSides, r, proximalMarker=[0,0]):
        assert(r>0)
        assert(numSides%2 == 0)
        self.numSides = numSides                                   #n_s
        self.r = r
        self.polygonInnerAngle = np.pi * (numSides-2)/(2*numSides) #\delta
        self.baseSideLength = 2*r*np.cos(self.polygonInnerAngle)
        self.width = self.baseSideLength * self.numSides
        
        self.proximalMarker = np.array(proximalMarker).reshape((1,2))
        proximalBaseX = np.arange(self.numSides) + 0.5
        nzeros = np.zeros(self.numSides)
        relative = np.vstack((proximalBaseX, nzeros)).T
        proximalBase = self.proximalMarker + relative
        self.proximalBaseIndices = np.arange(self.numSides)
        
              
        """
        Instance variables initialized here for an empty pattern:
        SHOULD ALL BE RESET OR ADJUSTED IN SUBCLASS CONSTRUCTORS
        """
        self.Vertices = proximalBase
        self.distalBaseIndices = self.proximalBaseIndices
        self.MountainEdges = np.empty(shape=(0,2), dtype=np.int32)
        self.ValleyEdges = np.empty(shape=(0,2), dtype=np.int32)
        self.patternHeight = 0
        self.distalMarker = self.proximalMarker
    
    def hasMountainFolds(self):
        return self.MountainEdges.shape[0] > 0
    
    def hasValleyFolds(self):
        return self.ValleyEdges.shape[0] > 0
    
    def isEmpty(self):
        return self.hasMountainFolds() or self.hasValleyFolds()
    
    def wrapToWidth(self):
        self.Vertices[:,0] %= self.width
        self.proximalMarker[0,0] %= self.width
        self.distalMarker[0,0] %= self.width
            
    def shift(self, v):
        self.Vertices += v
        self.proximalMarker = self.proximalMarker + v
        self.distalMarker = self.distalMarker + v
        self.wrapToWidth()
    
    def proximalBase(self):
        return self.Vertices[self.proximalBaseIndices,:]
    
    def distalBase(self):
        return self.Vertices[self.distalBaseIndices,:]
    
    """ index to roll other.proximalBase() by to align with self.distalBase(),
        or None if this is impossible """
    def rollToAlign(self, other : 'TubularPattern'):
        for i in range(self.numSides):
            if np.array_equal(self.distalBase(), 
                        np.roll(other.proximalBase(), i, axis=0)):
                return i
        return None
    
    def append(self, other : 'TubularPattern'):
        assert(other.numSides == self.numSides)
        assert(other.r == self.r)
        assert(self.MountainEdges.dtype.kind == 'i')
        assert(self.ValleyEdges.dtype.kind == 'i')
        other.shift(self.distalMarker - other.proximalMarker)
        currentNumVertices = self.Vertices.shape[0]
        mountainEdgesToAdd = other.MountainEdges + currentNumVertices
        self.MountainEdges = np.vstack((self.MountainEdges,mountainEdgesToAdd))
        valleyEdgesToAdd = other.ValleyEdges + currentNumVertices
        self.ValleyEdges = np.vstack((self.ValleyEdges, valleyEdgesToAdd))
        self.distalBaseIndices = other.distalBaseIndices + currentNumVertices
        self.Vertices = np.vstack((self.Vertices, other.Vertices))
        self.patternHeight += other.patternHeight
        self.distalMarker = other.distalMarker
    
    def MountainFolds(self):
        if (self.hasMountainFolds()):
            return self.Vertices[self.MountainEdges]
        else:
            return np.array([]).reshape((0,2,2))
    
    def ValleyFolds(self):
        if (self.hasValleyFolds()):
            return self.Vertices[self.ValleyEdges]
        else:
            return np.array([]).reshape((0,2,2))
    
    def makeDXF(self, saveas="", show=False):
        doc = ezdxf.new()

        # add new entities to the modelspace
        msp = doc.modelspace()
        doc.layers.add(name="Mountain", color=5)
        doc.layers.add(name="Valley", color=1)
        doc.layers.add(name="Boundary", color=2)
        doc.layers.add(name="Cut", color=3)
        
        ymin = self.proximalMarker[0,1]
        ymax = ymin + self.patternHeight
        xmin = 0
        xmax = self.width
        msp.add_line((xmin,ymin),(xmin,ymax), dxfattribs={"layer": "Boundary"})
        msp.add_line((xmin,ymin),(xmax,ymin), dxfattribs={"layer": "Boundary"})
        msp.add_line((xmax,ymax),(xmin,ymax), dxfattribs={"layer": "Boundary"})
        msp.add_line((xmax,ymax),(xmax,ymin), dxfattribs={"layer": "Boundary"})
        
        xmin = -0.5*self.baseSideLength
        xmax = self.width + 0.5*self.baseSideLength
        msp.add_line((xmin,ymin),(xmin,ymax), dxfattribs={"layer": "Cut"})
        msp.add_line((xmin,ymin),(xmax,ymin), dxfattribs={"layer": "Cut"})
        msp.add_line((xmax,ymax),(xmin,ymax), dxfattribs={"layer": "Cut"})
        msp.add_line((xmax,ymax),(xmax,ymin), dxfattribs={"layer": "Cut"})
        
        if self.hasMountainFolds():
            MF = self.MountainFolds()
            WrappingEdges = abs(MF[:,1,0]-MF[:,0,0]) > self.baseSideLength
            NonWrappingFolds = MF[np.logical_not(WrappingEdges)]
            WrappingFolds = MF[WrappingEdges]
            lastColStart = self.width - self.baseSideLength
            
            for seg in NonWrappingFolds:
                msp.add_line(seg[0,:],seg[1,:],dxfattribs={"layer":"Mountain"})
            for seg in WrappingFolds:
                rightcopy = np.copy(seg)
                leftcopy = np.copy(seg)
                rightcopy[seg[:,0] <= self.baseSideLength, 0] += self.width
                msp.add_line(rightcopy[0,:],rightcopy[1,:],
                             dxfattribs={"layer":"Mountain"})
                leftcopy[seg[:,0] >= lastColStart, 0] -= self.width
                msp.add_line(leftcopy[0,:],leftcopy[1,:],
                             dxfattribs={"layer":"Mountain"})
            FirstColFolds = NonWrappingFolds[np.logical_or(
                NonWrappingFolds[:,0,0] <= composed.baseSideLength, 
                NonWrappingFolds[:,1,0] <= composed.baseSideLength)]
            for seg in FirstColFolds:
                seg[:,0] += self.width
                msp.add_line(seg[0,:],seg[1,:],
                             dxfattribs={"layer":"Mountain"})
            
            LastColFolds = NonWrappingFolds[np.logical_or(
                NonWrappingFolds[:,0,0] >= lastColStart, 
                NonWrappingFolds[:,1,0] >= lastColStart)]
            for seg in LastColFolds:
                seg[:,0] -= self.width
                msp.add_line(seg[0,:],seg[1,:],
                             dxfattribs={"layer":"Mountain"})
            
        
        if self.hasValleyFolds():
            VF = self.ValleyFolds()
            WrappingEdges = abs(VF[:,1,0]-VF[:,0,0]) > self.baseSideLength
            NonWrappingFolds = VF[np.logical_not(WrappingEdges)]
            WrappingFolds = VF[WrappingEdges]
            lastColStart = self.width - self.baseSideLength
            
            for seg in NonWrappingFolds:
                msp.add_line(seg[0,:],seg[1,:],dxfattribs={"layer":"Valley"})
            for seg in WrappingFolds:
                rightcopy = np.copy(seg)
                leftcopy = np.copy(seg)
                rightcopy[seg[:,0] <= self.baseSideLength, 0] += self.width
                msp.add_line(rightcopy[0,:],rightcopy[1,:],
                             dxfattribs={"layer":"Valley"})
                leftcopy[seg[:,0] >= lastColStart, 0] -= self.width
                msp.add_line(leftcopy[0,:],leftcopy[1,:],
                             dxfattribs={"layer":"Valley"})
            FirstColFolds = NonWrappingFolds[np.logical_or(
                NonWrappingFolds[:,0,0] <= composed.baseSideLength, 
                NonWrappingFolds[:,1,0] <= composed.baseSideLength)]
            for seg in FirstColFolds:
                seg[:,0] += self.width
                msp.add_line(seg[0,:],seg[1,:],
                             dxfattribs={"layer":"Valley"})
            
            LastColFolds = NonWrappingFolds[np.logical_or(
                NonWrappingFolds[:,0,0] >= lastColStart, 
                NonWrappingFolds[:,1,0] >= lastColStart)]
            for seg in LastColFolds:
                seg[:,0] -= self.width
                msp.add_line(seg[0,:],seg[1,:],
                             dxfattribs={"layer":"Valley"})
        
        if saveas:
            doc.saveas(saveas+".dxf")
        
        if show:
            config = Configuration()
            fig: plt.Figure = plt.figure(figsize=(xmax-xmin, ymax-ymin))
            ax: plt.Axes = fig.add_axes([0, 0, 1, 1])
            ctx = RenderContext(doc)
            out = MatplotlibBackend(ax)
            # get the modelspace properties
            msp_properties = LayoutProperties.from_layout(msp)
            # set light gray background color and black foreground color
            msp_properties.set_colors("#eaeaea")
            Frontend(ctx, out, config=config).draw_layout(msp, finalize=False,
                                            layout_properties=msp_properties)
            ax.set_ylim(ymin,ymax)
            ax.set_xlim(xmin,xmax)
            ax.show()
        
        return doc
        

class Tube(TubularPattern):
    def __init__(self, numSides, r, height, proximalMarker=[0,0]):
        super().__init__(numSides, r, proximalMarker)
        self.patternHeight = height
        proximalBase = self.Vertices
        distalBase = proximalBase + [0, self.patternHeight]
        self.Vertices = np.vstack((proximalBase, distalBase))
        self.distalBaseIndices = self.numSides + self.proximalBaseIndices
        self.MountainEdges = np.vstack((self.proximalBaseIndices, 
                                        self.distalBaseIndices)).T
        self.distalMarker = self.proximalMarker + [0, self.patternHeight]
        self.wrapToWidth()
        
        
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
        self.patternHeight = np.sqrt(moduleHeight**2 + 
                        (self.baseSideLength * 
                         (1/np.sin(np.pi/numSides)) *
                         np.sin((np.pi/numSides) - (baseTwist/2)) *
                         np.sin(baseTwist/2) )**2 )
        proximalBase = self.Vertices
        distalBase = proximalBase + [baseXshift, self.patternHeight]
        distalBase[:,0] %= self.width
        self.Vertices = np.vstack((proximalBase, distalBase))
        
        self.distalBaseIndices = self.numSides + self.proximalBaseIndices
        rightDiagonals = np.vstack((self.proximalBaseIndices, 
                                    self.distalBaseIndices)).T
        leftDiagonals = np.vstack((self.proximalBaseIndices, 
                                   np.roll(self.distalBaseIndices,1))).T
        self.MountainEdges = np.vstack((rightDiagonals, leftDiagonals))
        
        proximalValleys = np.vstack((self.proximalBaseIndices, 
                                     np.roll(self.proximalBaseIndices,1))).T
        distalValleys = np.vstack((self.distalBaseIndices, 
                                   np.roll(self.distalBaseIndices,1))).T
        self.ValleyEdges = np.vstack((proximalValleys,distalValleys))
        
        referenceXshift = (np.floor(numSides * twistAngle / (2*np.pi)) * \
                            self.baseSideLength) + baseXshift
        self.distalMarker = self.proximalMarker + \
                            np.array([referenceXshift, self.patternHeight])
        self.wrapToWidth()
        

""" Testing """
r = 1
numSides = 6
composed = TubularPattern(numSides, r)
tube = Tube(numSides, r, 3)
tube.makeDXF(saveas="tube")
composed.append(Tube(numSides, r, 3))
twist = Twist(numSides, r, 0.45*np.pi, 1)
twist.makeDXF(saveas="twist")
composed.append(Twist(numSides, r, 0.45*np.pi, 1))
composed.append(Tube(numSides, r, 3))
composed.makeDXF(saveas="tube_twist_tube", show=True)
