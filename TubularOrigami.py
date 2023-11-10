# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 17:11:58 2023

@author: dfesh

Implementing the catalog of tubular origami crease patterns from the original 
Kinegami paper (Chen et al. 2022).

These functions construct 2D crease patterns which fold into tubular modules
of radius r whose ends are identical regular polygons where the number of sides
is even.

Each module stores a proximal reference point along the proximal base,
and a distal reference point along the distal base, to encode how they 
align together when composed.

We represent them with a graph structure that encodes the geometry and topology
in a graph structure where vertices (points on the surface) are parameterized
in 2D by:
    x - distance about the tube modulo numSides*baseSideLength, i.e., the
        horizontal coordinate when cut and unwrapped to flat
    y - distance about the tube, i.e., the vertical coordiante when flat

In tubular topology, a non-vertical edge can either proceed clockwise or 
counterclockwise from one vertex to the other. To disambiguate and enable us
to detect which edges wrap around in the x coordinate, we encode all edges in
the +x direction. This means an edge (i,j) proceeds from vertex i in the 
+x direction (possibly with wraparound) to vertex j. Therefore, if 
vertex j is to the left of vertex i, we know the fold has wraparound. 
(When detecting vertical edges, an epsilon term provides numerical stability.)
Obeying this orientation is important when implementing new fold patterns. 

To make this manufacturable in a flat topology (e.g. by a laser etching),
we consider it to be cut along the x=0 axis, with the geometry in 
x=[0,baseSideLength] duplicated at the end so that the wraparound can adhere
together along a 2D surface. This duplication, along with other difficulties of
converting wraparound edges to a flat pattern, is dealt with in the method
that constructs a DXF file. (We organize in this way to separate fabrication 
details from the underlying pattern: the construction and representation of 
TubularPattern objects is in the tubular topology).


TODO:
    Make append copy, not modify, the other pattern
    
    Make variable naming consistent:
        Matrices/arrays of >1 axis: Capitalize
        Vectors/1-axis arrays: plural, don't capitalize
        Vcalars: singular, don't capitalize
    
    Review all uses of %, np.mod, and math.remainder to make sure the correct
    one is being used for the situation
    
    Figure out how to append revolute joints (mountain folds at the bases)
    and twist fittings (valley folds at the bases)
"""
import numpy as np
import matplotlib.pyplot as plt
import ezdxf
from ezdxf.addons.drawing import RenderContext, Frontend
from ezdxf.addons.drawing.matplotlib import MatplotlibBackend
from ezdxf.addons.drawing.config import Configuration
from ezdxf.addons.drawing.properties import LayoutProperties
import math
from math import remainder
import geometryHelpers
from geometryHelpers import rollToMatch, unsignedAngles, signedAngle, \
    signedAngles2D, signedAngles3D
from spatialmath import SO3, SE3

class TubularPattern():
    def __init__(self, numSides, r, proximalMarker=[0, 0]):
        assert (r > 0)
        assert (numSides % 2 == 0)
        self.numSides = numSides  # n_s
        self.r = r
        self.polygonInnerAngle = np.pi * (numSides-2)/(2*numSides)  # \delta
        self.baseSideLength = 2*r*np.cos(self.polygonInnerAngle)
        self.width = self.baseSideLength * self.numSides

        self.proximalMarker = np.array(proximalMarker).reshape((1, 2))
        ProximalBaseX = np.arange(self.numSides)*self.baseSideLength\
                                            +0.5*self.baseSideLength
        nzeros = np.zeros(self.numSides)
        relative = np.vstack((ProximalBaseX, nzeros)).T
        ProximalBase = self.proximalMarker + relative
        self.proximalBaseIndices = np.arange(self.numSides)
        self.EPSILON = self.baseSideLength * 0.00001

        """
        Instance variables initialized here for an empty pattern:
        SHOULD ALL BE RESET OR ADJUSTED IN SUBCLASS CONSTRUCTORS
        """
        self.Vertices = ProximalBase
        self.distalBaseIndices = self.proximalBaseIndices
        self.Edges = np.empty(shape=(0, 2), dtype=np.int32)
        self.EdgeLabelsMV = np.empty(0, dtype=bool) #True for M, False for V
        self.patternHeight = 0
        self.distalMarker = self.proximalMarker

    """ Add the given edges, with the given mountain/valley labels, to the 
        crease pattern. Filters our self-edges (i,i). If deDuplicate is set
        to True, also filters out duplicate edges. """
    def addEdges(self, EdgesToAdd, LabelsToAdd, deDuplicate=False):
        assert(EdgesToAdd.shape[0] == LabelsToAdd.shape[0])
        assert(EdgesToAdd.shape[1] == 2)
        nonEmpty = EdgesToAdd[:,0] != EdgesToAdd[:,1]
        self.Edges = np.vstack((self.Edges, EdgesToAdd[nonEmpty]))
        self.EdgeLabelsMV = np.hstack((self.EdgeLabelsMV, 
                                       LabelsToAdd[nonEmpty]))
        
        if deDuplicate:
            self.Edges, indices = np.unique(self.Edges, return_index=True, 
                                            axis=0)
            self.EdgeLabelsMV = self.EdgeLabelsMV[indices]
            # TODO: is there a nice way to verify that duplicated edges have 
            # the same MV label?

    def addMountainEdges(self, EdgesToAdd, deDuplicate=False):
        EdgesToAdd = np.array(EdgesToAdd)
        newTrues = np.ones(EdgesToAdd.shape[0], dtype=bool)
        self.addEdges(EdgesToAdd, newTrues, deDuplicate)
    
    def addValleyEdges(self, EdgesToAdd, deDuplicate=False):
        EdgesToAdd = np.array(EdgesToAdd)
        newFalses = np.zeros(EdgesToAdd.shape[0], dtype=bool)
        self.addEdges(EdgesToAdd, newFalses, deDuplicate)

    def wrapToWidth(self):
        self.Vertices[:, 0] %= self.width
        self.proximalMarker[0, 0] %= self.width
        self.distalMarker[0, 0] %= self.width

    def shift(self, v):
        self.Vertices += v
        self.proximalMarker = self.proximalMarker + v
        self.distalMarker = self.distalMarker + v
        self.wrapToWidth()

    def ProximalBase(self):
        return self.Vertices[self.proximalBaseIndices, :]

    def DistalBase(self):
        return self.Vertices[self.distalBaseIndices, :]

    def rollProximalIndices(self, k):
        PI = self.proximalBaseIndices
        self.Vertices[PI, :] = np.roll(self.Vertices[PI, :], k, axis=0)
        Mask = self.Edges < self.numSides
        self.Edges[Mask] = (self.Edges[Mask] + k) % self.numSides
        self.proximalBaseIndices = np.roll(self.proximalBaseIndices, k)
    
    def append(self, other: 'TubularPattern'):
        assert (other.numSides == self.numSides)
        assert (other.r == self.r)
        assert (self.Edges.dtype.kind == 'i')
        other.shift(self.distalMarker - other.proximalMarker)

        """ Merge other's proximal base vertices with current distal base """
        indexShift = self.Vertices.shape[0] - other.numSides
        otherNonProximalVertices = np.delete(other.Vertices,
                                             other.proximalBaseIndices, axis=0)
        self.Vertices = np.vstack((self.Vertices, otherNonProximalVertices))

        """ Reindex other's proximal base vertices such that they match
            the current distal base vertices they correspond to """
        k = rollToMatch(other.ProximalBase(), self.DistalBase(),
                        EPSILON=self.EPSILON)
        assert (not k is None)
        other.rollProximalIndices(k)
        EdgesToAdd = other.Edges
        EdgesToAdd[EdgesToAdd>=self.numSides] += indexShift #non-proximal
        EdgesToAdd[EdgesToAdd<self.numSides] += np.min(self.distalBaseIndices)
        self.addEdges(EdgesToAdd, other.EdgeLabelsMV, deDuplicate=True)
        
        self.distalBaseIndices = other.distalBaseIndices + indexShift
        self.patternHeight += other.patternHeight
        self.distalMarker = other.distalMarker
    
    """ Plot the graph directly, without accounting for wraparound, boundary,
        or other considerations for actual manufacturing.
        Mainly useful for debugging - comparing it to the result of makeDXF
        helps separate mistakes in the actual graph from bugs in makeDXF. """
    def plotRawGraph(self, saveas="", show=True, directed=False):
        doc = ezdxf.new()

        # add new entities to the modelspace
        msp = doc.modelspace()
        if directed:
            doc.layers.add(name="MountainCW", color=5)
            doc.layers.add(name="MountainCCW", color=4)
            doc.layers.add(name="ValleyCW", color=1)
            doc.layers.add(name="ValleyCCW", color=6)
        else:
            doc.layers.add(name="Mountain", color=5)
            doc.layers.add(name="Valley", color=1)

        ymin = self.proximalMarker[0, 1]
        ymax = ymin + self.patternHeight
        xmin = -0.5*self.baseSideLength
        xmax = self.width + 0.5*self.baseSideLength
        
        Folds = self.Vertices[self.Edges]
        for i in range(self.Edges.shape[0]):
            seg = Folds[i,:,:]
            layer = "Mountain" if self.EdgeLabelsMV[i] else "Valley"
            if directed:
                layer += "CCW" if seg[0,0] > seg[1,0]+self.EPSILON else "CW"
            msp.add_line(seg[0, :], seg[1, :], dxfattribs={"layer": layer})
            
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
            ax.set_ylim(ymin, ymax)
            ax.set_xlim(xmin, xmax)
            plt.show()

        return doc
    
    def makeDXF(self, saveas="", show=False, debug=False):
        doc = ezdxf.new()

        # add new entities to the modelspace
        msp = doc.modelspace()
        doc.layers.add(name="Mountain", color=5)
        doc.layers.add(name="Valley", color=1)
        doc.layers.add(name="Boundary", color=2)
        doc.layers.add(name="Cut", color=3)

        ymin = self.proximalMarker[0, 1]
        ymax = ymin + self.patternHeight
        msp.add_line((0, ymin), (0, ymax),
                     dxfattribs={"layer": "Boundary"})
        msp.add_line((0, ymin), (self.width, ymin),
                     dxfattribs={"layer": "Boundary"})
        msp.add_line((self.width, ymax), (0, ymax),
                     dxfattribs={"layer": "Boundary"})
        msp.add_line((self.width, ymax), (self.width, ymin),
                     dxfattribs={"layer": "Boundary"})

        xmin = -0.5*self.baseSideLength
        xmax = self.width + 0.5*self.baseSideLength
        msp.add_line((xmin, ymin), (xmin, ymax), dxfattribs={"layer": "Cut"})
        msp.add_line((xmin, ymin), (xmax, ymin), dxfattribs={"layer": "Cut"})
        msp.add_line((xmax, ymax), (xmin, ymax), dxfattribs={"layer": "Cut"})
        msp.add_line((xmax, ymax), (xmax, ymin), dxfattribs={"layer": "Cut"})
        
        Folds = self.Vertices[self.Edges]
        lastColStart = self.width - self.baseSideLength
        if debug:
            print("Vertices:")
            print(self.Vertices)
            print("Edges:")
            print(self.Edges)
            print("Folds:")
            print(Folds)
        for i in range(self.Edges.shape[0]):
            seg = Folds[i,:,:]
            isWrapping = seg[0,0] > seg[1,0]+self.EPSILON
            layer = "Mountain" if self.EdgeLabelsMV[i] else "Valley"
            
            if isWrapping:
                rightcopy = np.copy(seg)
                leftcopy = np.copy(seg)
                rightcopy[seg[:, 0] <= self.baseSideLength, 0] += self.width
                msp.add_line(rightcopy[0, :], rightcopy[1, :],
                             dxfattribs={"layer": layer})
                leftcopy[seg[:, 0] >= lastColStart, 0] -= self.width
                msp.add_line(leftcopy[0, :], leftcopy[1, :],
                             dxfattribs={"layer": layer})
            else:
                msp.add_line(seg[0, :], seg[1, :], dxfattribs={"layer": layer})
            
                # If it intersects with the first column, duplicate on right
                firstColEnd = self.baseSideLength + self.EPSILON
                if (seg[0,0] <= firstColEnd) or (seg[1,0] <= firstColEnd):
                    rightcopy = np.copy(seg)    
                    rightcopy[:,0] += self.width
                    msp.add_line(rightcopy[0, :], rightcopy[1, :],
                                 dxfattribs={"layer": layer})
                
                # If it intersects with the last column, duplicate on left
                lastColStart = self.width - self.baseSideLength - self.EPSILON
                if (seg[0,0]>=lastColStart) or (seg[1,0]>=lastColStart):
                    leftcopy = np.copy(seg)
                    leftcopy[:,0] -= self.width
                    msp.add_line(leftcopy[0, :], leftcopy[1, :],
                                 dxfattribs={"layer": layer})
            
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
            ax.set_ylim(ymin, ymax)
            ax.set_xlim(xmin, xmax)
            plt.show()

        return doc


class Tube(TubularPattern):
    def __init__(self, numSides, r, height, proximalMarker=[0, 0]):
        super().__init__(numSides, r, proximalMarker)
        self.patternHeight = height
        ProximalBase = self.Vertices
        DistalBase = ProximalBase + [0, self.patternHeight]
        self.Vertices = np.vstack((ProximalBase, DistalBase))
        self.distalBaseIndices = self.numSides + self.proximalBaseIndices
        self.addMountainEdges(np.vstack((self.proximalBaseIndices,
                                        self.distalBaseIndices)).T)
        self.distalMarker = self.proximalMarker + [0, self.patternHeight]
        self.wrapToWidth()

class RevoluteJointPattern(TubularPattern):
    def __init__(self, numSides, r, totalBendingAngle, numSinkLayers, 
                 proximalMarker=[0, 0]):
        super().__init__(numSides, r, proximalMarker)
        rotAxisHeight = r*np.sin(self.polygonInnerAngle)*\
                            np.tan(totalBendingAngle/4)
        flatRotAxisHeight = r*np.sin(self.polygonInnerAngle)/\
                                np.cos(totalBendingAngle/4)
        self.patternHeight = 2*flatRotAxisHeight
        ProximalBase = self.Vertices
        DistalBase = ProximalBase + [0, self.patternHeight]
        self.Vertices = np.vstack((ProximalBase, DistalBase))
        self.distalBaseIndices = self.numSides + self.proximalBaseIndices
        MidAlignedVertices = ProximalBase + [0, flatRotAxisHeight]
        midAlignedIndices = self.Vertices.shape[0] + np.arange(numSides)
        self.Vertices = np.vstack((self.Vertices, MidAlignedVertices))
        MidHalfVertices = MidAlignedVertices + [self.baseSideLength/2, 0]
        midHalfIndices = self.Vertices.shape[0] + np.arange(numSides)
        self.Vertices = np.vstack((self.Vertices, MidHalfVertices))
        
        # Base mountain folds
        self.addMountainEdges(np.vstack((self.proximalBaseIndices, 
                                np.roll(self.proximalBaseIndices, -1))).T)
        self.addMountainEdges(np.vstack((self.distalBaseIndices,
                                np.roll(self.distalBaseIndices, -1))).T)
        
        # Folds of the flapping rectangles
        j = numSides-1
        k = (numSides-1)//2
        self.addMountainEdges([
                    [self.proximalBaseIndices[0], midAlignedIndices[0]],
                    [midAlignedIndices[0], midHalfIndices[0]],
                    [midAlignedIndices[0], self.distalBaseIndices[0]], 
                    #######
                    [self.proximalBaseIndices[k], midAlignedIndices[k]],
                    [midHalfIndices[k-1], midAlignedIndices[k]],
                    [midAlignedIndices[k], self.distalBaseIndices[k]], 
                    #######
                    [self.proximalBaseIndices[k+1], midAlignedIndices[k+1]],
                    [midAlignedIndices[k+1], midHalfIndices[k+1]],
                    [midAlignedIndices[k+1], self.distalBaseIndices[k+1]],
                    #######
                    [self.proximalBaseIndices[j], midAlignedIndices[j]],
                    [midHalfIndices[j-1], midAlignedIndices[j]],
                    [midAlignedIndices[j], self.distalBaseIndices[j]]
                ])
        self.addValleyEdges([    [midAlignedIndices[k], midHalfIndices[k]],
                                 [midHalfIndices[k], midAlignedIndices[k+1]],
                                 [midAlignedIndices[j], midHalfIndices[j]],
                                 [midHalfIndices[j], midAlignedIndices[0]]  ])
        
        """
        Why are these things experessed in Chen et al. 2022 as dependent on i?
        It seems like they should be the same for all the indices on which
        they're actually defined. But for now I'll write it this way to match.
        """
        psi = np.arctan2(self.baseSideLength, self.patternHeight)
        deltas = np.pi * (numSides + 2 - 4*(np.arange(numSides)+1)) / (2*numSides)
        frontLengths = np.sqrt(rotAxisHeight**2 + (r*np.sin(deltas))**2)
        gammas = np.arctan(frontLengths / (r*np.abs(np.cos(deltas))))
        tuckYs = frontLengths / np.cos((np.pi/2) - gammas - psi)
        
        for i in range(numSides):
            # Crimping vertices and folds
            if not i in [0, k, k+1, j]:
                lowerCrimpVertex = ProximalBase[i] + [0,tuckYs[i]]
                lowerCrimpIndex = self.Vertices.shape[0]
                upperCrimpVertex = DistalBase[i] - [0,tuckYs[i]]
                upperCrimpIndex = self.Vertices.shape[0]+1
                self.Vertices = np.vstack((self.Vertices, 
                                          [lowerCrimpVertex,upperCrimpVertex]))
                self.addValleyEdges([[midHalfIndices[i-1], upperCrimpIndex],
                                [midHalfIndices[i-1], midAlignedIndices[i]],
                                [midAlignedIndices[i], midHalfIndices[i]],
                                [lowerCrimpIndex, midAlignedIndices[i]]])
                self.addMountainEdges([[midHalfIndices[i-1], lowerCrimpIndex],
                                       [lowerCrimpIndex, midHalfIndices[i]],
                                       [upperCrimpIndex, midHalfIndices[i]],
                                       [midAlignedIndices[i],upperCrimpIndex]])
        
        for i in range(numSides):    
            # Tucking X folds    
            if not i in [(numSides-1)//2, numSides-1]:
                mhi = midHalfIndices[i]
                self.addValleyEdges([[self.proximalBaseIndices[i], mhi],
                                     [self.distalBaseIndices[i], mhi],
                                     [mhi, self.proximalBaseIndices[i+1]],
                                     [mhi, self.distalBaseIndices[i+1]]])    
        
        #TODO: RECURSIVE SINK GADGET
        
        self.distalMarker = self.proximalMarker + [0, self.patternHeight]
        self.wrapToWidth()

class PrismaticJointPattern(TubularPattern):
    def __init__(self, numSides, r, neutralLength, numLayers, coneAngle,
                 proximalMarker=[0, 0]):
        super().__init__(numSides, r, proximalMarker)
        assert(coneAngle > 0 and coneAngle < np.pi/2)
        assert(numLayers > 0)
        neutralLayerHeight = neutralLength / numLayers
        flatLayerHalfHeight = neutralLayerHeight / (2*np.sin(coneAngle))
        maxLength = 2*numLayers*flatLayerHalfHeight
        self.patternHeight = 4*maxLength
        flatConeAngle = (np.pi/2)*(1 - 2*np.cos(coneAngle)/numSides)
        
        # Standard distal base setup
        ProximalBase = self.Vertices
        DistalBase = ProximalBase + [0, self.patternHeight]
        self.Vertices = np.vstack((ProximalBase, DistalBase))
        self.distalBaseIndices = numSides + self.proximalBaseIndices
        
        Mid1 = ProximalBase + [0, maxLength]
        mid1Indices = numSides + self.distalBaseIndices
        self.Vertices = np.vstack((self.Vertices,Mid1))
        
        xsTiled = np.tile(ProximalBase[:,0], numLayers+1)
        reboY0 = proximalMarker[1] + 2*maxLength
        layerYs = reboY0 + np.arange(numLayers+1)*2*flatLayerHalfHeight
        ysRepeated = np.repeat(layerYs, numSides)
        ReboLayerBoundariesVertices = np.vstack((xsTiled, ysRepeated)).T
        ReboLayerBoundariesIndices = self.Vertices.shape[0] + \
            np.arange((numLayers+1)*numSides).reshape((numLayers+1, numSides))
        self.Vertices = np.vstack((self.Vertices, ReboLayerBoundariesVertices))
        # TODO: can this be written more elegantly with np.meshgrid?
        
        ReboLayerMidAlignedVertices = [0, flatLayerHalfHeight] +\
                                    ReboLayerBoundariesVertices[:-numSides,:]
        ReboLayerMidAlignedIndices = self.Vertices.shape[0] + \
            np.arange(numLayers*numSides).reshape((numLayers, numSides))
        self.Vertices = np.vstack((self.Vertices, ReboLayerMidAlignedVertices))
        
        offset = [-flatLayerHalfHeight/np.tan(flatConeAngle), 0]
        ReboLayerMidOffsetVertices = ReboLayerMidAlignedVertices + offset
        ReboLayerMidOffsetIndices = self.Vertices.shape[0] + \
            np.arange(numLayers*numSides).reshape((numLayers, numSides))
        self.Vertices = np.vstack((self.Vertices, ReboLayerMidOffsetVertices))
        
        self.addMountainEdges(np.vstack((self.proximalBaseIndices,
                                         mid1Indices)).T)
        self.addMountainEdges(np.vstack((mid1Indices, 
                                         np.roll(mid1Indices, -1, axis=0))).T)
        
        reboBottom = ReboLayerBoundariesIndices[0,:]
        reboTop = ReboLayerBoundariesIndices[numLayers,:]
        self.addMountainEdges(np.vstack((mid1Indices, reboBottom)).T)
        self.addMountainEdges(np.vstack((reboTop, self.distalBaseIndices)).T)
        
        # Construct REBO folds
        self.addValleyEdges(np.vstack((reboTop, np.roll(reboTop,-1))).T)
        # TODO: vectortize?
        for i in range(numLayers):
            bottom = ReboLayerBoundariesIndices[i,:]
            nextBottom = np.roll(bottom, -1)
            midAligned = ReboLayerMidAlignedIndices[i,:]
            midOffset = ReboLayerMidOffsetIndices[i,:]
            nextMidOffset = np.roll(midOffset, -1)
            top = ReboLayerBoundariesIndices[i+1,:]
            
            self.addValleyEdges(np.vstack((bottom, nextBottom)).T)
            self.addValleyEdges(np.vstack((bottom, midAligned)).T)
            self.addValleyEdges(np.vstack((midAligned, top)).T)
            self.addValleyEdges(np.vstack((midOffset, midAligned)).T)
            self.addMountainEdges(np.vstack((midAligned, nextMidOffset)).T)
            self.addMountainEdges(np.vstack((midOffset, bottom)).T)
            self.addMountainEdges(np.vstack((midOffset, top)).T)
            
        self.distalMarker = self.proximalMarker + [0, self.patternHeight]
        self.wrapToWidth()  
        
        
        
    
class ElbowFitting(TubularPattern):
    def __init__(self, numSides, r, bendingAngle, rotationalAxisAngle,
                 proximalMarker=[0,0]):
        super().__init__(numSides, r, proximalMarker)
        rotationalAxisAngle = np.mod(rotationalAxisAngle, 2*np.pi)
        bendingAngle = math.remainder(bendingAngle, 2*np.pi) #wrap to [-pi,pi]
        assert(abs(bendingAngle) < np.pi)
        assert(abs(bendingAngle) > self.EPSILON)
        if bendingAngle < 0:
            bendingAngle = abs(bendingAngle)
            rotationalAxisAngle = np.mod(rotationalAxisAngle+np.pi, 2*np.pi)
        
        dw = self.r * np.tan(bendingAngle / 2)
        baseAngles = (2*np.pi/self.numSides)*(np.arange(self.numSides)+0.5)
        midPolygonHeights = dw * (np.sin(baseAngles - rotationalAxisAngle) + 1)
        maxMidPolygonHeight = np.max(midPolygonHeights)
        self.patternHeight = 2*maxMidPolygonHeight
        
        ProximalBase = self.Vertices
        DistalBase = ProximalBase + [0, self.patternHeight]
        self.Vertices = np.vstack((ProximalBase, DistalBase))
        self.distalBaseIndices = self.numSides + self.proximalBaseIndices
        self.distalMarker = self.proximalMarker + [0, self.patternHeight]
        
        Heights0Y = np.vstack((np.zeros(self.numSides), midPolygonHeights)).T
        LowerTuckBoundary = ProximalBase + Heights0Y
        UpperTuckBoundary = DistalBase - Heights0Y
        
        Prev2D = np.roll(LowerTuckBoundary, 1, axis=0)
        Prev2D[0,0] -= self.width
        Next2D = np.roll(LowerTuckBoundary, -1, axis=0)
        Next2D[-1,0] += self.width
        Back2D = Prev2D - LowerTuckBoundary
        Forward2D = Next2D - LowerTuckBoundary
        Angles2D = np.mod(signedAngles2D(Forward2D, Back2D), 2*np.pi)
        
        MidPolygon3D = np.vstack((self.r * np.cos(baseAngles),
                                  self.r * np.sin(baseAngles),
                                  midPolygonHeights)).T
        rotAxis = np.array([np.cos(rotationalAxisAngle), 
                            np.sin(rotationalAxisAngle), 
                            0])
        midNormal3D = SO3.AngVec(bendingAngle/2, rotAxis) * np.array([0,0,1])
        Prev3D = np.roll(MidPolygon3D, 1, axis=0)
        Next3D = np.roll(MidPolygon3D, -1, axis=0)
        Back3D = Prev3D-MidPolygon3D
        Forward3D = Next3D-MidPolygon3D
        Angles3D = signedAngles3D(Forward3D, Back3D, midNormal3D)
        
        TuckAngles = (Angles2D - Angles3D)/2 #\varepsilon_i for each i
        mMPHcopies = maxMidPolygonHeight * np.ones(self.numSides) 
        MidAligned = np.vstack((ProximalBase[:,0], mMPHcopies)).T
        MidOffset = np.vstack((ProximalBase[:,0] +\
                        (mMPHcopies-midPolygonHeights)*np.tan(TuckAngles), 
                        mMPHcopies)).T
            
        lowerTuckBoundaryIndices = []
        upperTuckBoundaryIndices = []
        midAlignedIndices = []
        midOffsetIndices = []
        VerticesToAdd = []
        nv = self.Vertices.shape[0] #number of vertices / next vertex index
        for i in range(self.numSides):
            height = midPolygonHeights[i]
            if abs(height-maxMidPolygonHeight) < self.EPSILON:
                # all four are the same vertex, don't want to duplicate
                VerticesToAdd.append(MidAligned[i])
                lowerTuckBoundaryIndices.append(nv)
                upperTuckBoundaryIndices.append(nv)
                midAlignedIndices.append(nv)
                midOffsetIndices.append(nv)
                nv += 1
            else:
                # add the four vertices separately
                lowerTuckBoundaryIndices.append(nv)
                VerticesToAdd.append(LowerTuckBoundary[i]) #lower
                midAlignedIndices.append(nv+1)
                VerticesToAdd.append(MidAligned[i])
                midOffsetIndices.append(nv+2)
                VerticesToAdd.append(MidOffset[i])
                upperTuckBoundaryIndices.append(nv+3)
                VerticesToAdd.append(UpperTuckBoundary[i]) #upper
                nv += 4
        
        self.Vertices = np.vstack((self.Vertices, VerticesToAdd))
        
        lowerBoundaryEdges = np.vstack((np.roll(lowerTuckBoundaryIndices, 1), 
                                       lowerTuckBoundaryIndices)).T
        self.addMountainEdges(lowerBoundaryEdges)
        lowerVerticalEdges = np.vstack((lowerTuckBoundaryIndices,
                                        self.proximalBaseIndices)).T
        self.addMountainEdges(lowerVerticalEdges)
        
        upperBoundaryEdges = np.vstack((np.roll(upperTuckBoundaryIndices, 1), 
                                               upperTuckBoundaryIndices)).T
        self.addMountainEdges(upperBoundaryEdges, deDuplicate=True)
        upperVerticalEdges = np.vstack((upperTuckBoundaryIndices,
                                        self.distalBaseIndices)).T
        self.addMountainEdges(upperVerticalEdges, deDuplicate=True)
        
        # Horizontal tucking edges (across the middle)
        self.addValleyEdges(np.vstack((midAlignedIndices, midOffsetIndices)).T)
        nextMidAligned = np.roll(midAlignedIndices, -1)
        self.addValleyEdges(np.vstack((midOffsetIndices, nextMidAligned)).T)
        
        # Vertical tucking edges (aligned with the base vertices)
        self.addValleyEdges(np.vstack((lowerTuckBoundaryIndices, 
                                       midAlignedIndices)).T)
        self.addMountainEdges(np.vstack((midAlignedIndices, 
                                         upperTuckBoundaryIndices)).T)
        
        # Diagonal tucking edges (boundary vertices to mid offset vertices)
        self.addMountainEdges(np.vstack((lowerTuckBoundaryIndices,
                                        midOffsetIndices)).T, deDuplicate=True)
        self.addValleyEdges(np.vstack((upperTuckBoundaryIndices,
                                       midOffsetIndices)).T, deDuplicate=True)
        self.distalMarker = self.proximalMarker + [0, self.patternHeight]
        self.wrapToWidth()


class Twist(TubularPattern):
    def __init__(self, numSides, r, twistAngle, moduleHeight,
                 proximalMarker=[0, 0]):
        super().__init__(numSides, r, proximalMarker)
        twistAngle = twistAngle % (2*np.pi)
        """ 
        Twist in increments of 2pi/numSides is be accomplished by rotating the
        next modules itself (shift crease pattern columns modularly).
        So the twist module only needs to instantiate the remaining portion:
        """
        baseTwist = twistAngle % (2*np.pi / self.numSides)  # \overline{\alpha}
        baseXshift = (self.baseSideLength/2) * (1 - np.cos(baseTwist) +
                                (1/np.tan(np.pi/numSides))*np.sin(baseTwist))
        self.patternHeight = np.sqrt(moduleHeight**2 +
                                     (self.baseSideLength *
                                     (1/np.sin(np.pi/numSides)) *
                                     np.sin((np.pi/numSides) - (baseTwist/2)) *
                                     np.sin(baseTwist/2))**2)
        ProximalBase = self.Vertices
        DistalBase = ProximalBase + [baseXshift, self.patternHeight]
        DistalBase[:, 0] %= self.width
        self.Vertices = np.vstack((ProximalBase, DistalBase))

        self.distalBaseIndices = self.numSides + self.proximalBaseIndices
        rightDiagonals = np.vstack((self.proximalBaseIndices,
                                    self.distalBaseIndices)).T
        leftDiagonals = np.vstack((np.roll(self.distalBaseIndices, 1),
                                   self.proximalBaseIndices)).T
        self.addMountainEdges(np.vstack((rightDiagonals, leftDiagonals)))

        
        self.addValleyEdges(np.vstack((np.roll(self.proximalBaseIndices, 1),
                                       self.proximalBaseIndices)).T)
        self.addValleyEdges(np.vstack((np.roll(self.distalBaseIndices, 1),
                                       self.distalBaseIndices)).T)
        
        referenceXshift = (np.floor(numSides * twistAngle / (2*np.pi)) *
                           self.baseSideLength) + baseXshift
        self.distalMarker = self.proximalMarker + \
            np.array([referenceXshift, self.patternHeight])
        self.wrapToWidth()


""" Testing """
r = 1
numSides = 8

composed = TubularPattern(numSides, r)
tube = Tube(numSides, r, 2)
#tube.plotRawGraph(saveas="tube_raw_graph", directed=True)
tube.makeDXF(saveas="tube", show=False)


#tube.plotRawGraph(directed=True)
composed.append(Tube(numSides, r, 3))
twist = Twist(numSides, r, 0.45*np.pi, 1)
#twist.plotRawGraph(saveas="twist_raw_graph", directed=True)
twist.makeDXF(saveas="twist", show=False)
composed.append(Twist(numSides, r, 0.9*np.pi, 1))
composed.append(Tube(numSides, r, 3))
composed.makeDXF(saveas="tube_twist_tube", show=False)
#composed.plotRawGraph(saveas="tube_twist_tube_raw_graph", directed=True)


doubleTwist = TubularPattern(numSides, r)
doubleTwist.append(Twist(numSides, r, 0.9*np.pi, 1))
doubleTwist.append(Twist(numSides, r, 0.3*np.pi, 1))
doubleTwist.makeDXF(saveas="twist_twist", show=False)
#doubleTwist.plotRawGraph(saveas="twist_twist_raw_graph", directed=True)

elbow = ElbowFitting(numSides, r, -np.pi/4, np.pi/3)
#elbow.plotRawGraph(saveas="elbow_raw_graph", directed=True)
elbow.makeDXF(saveas="elbow", show=False)
composed.append(elbow)
composed.makeDXF(saveas="tube_twist_tube_elbow", show=False)


prismatic = PrismaticJointPattern(numSides, r, 1, 3, np.pi/3)
#prismatic.plotRawGraph(saveas="prismatic_raw_graph", directed=True)
prismatic.makeDXF(saveas="prismatic", show=False)
composed.append(prismatic)
composed.makeDXF(saveas="tube_twist_tube_elbow_prismatic", show=False)
composed.append(Tube(numSides, r, 1))
composed.makeDXF(saveas="tube_twist_tube_elbow_prismatic_tube", show=True)

revolute = RevoluteJointPattern(numSides, r, np.pi, 0)
revolute.plotRawGraph(saveas="revolute_raw_graph", directed=True)
revolute.makeDXF(saveas="revolute", show=True)
