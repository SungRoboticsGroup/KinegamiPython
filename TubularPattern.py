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
    x - distance about the tube modulo width=numSides*baseSideLength, i.e., the
        horizontal coordinate when cut and unwrapped to flat
    y - distance along the tube, i.e., the vertical coordinate when flat

The fold segment connecting vertices (x0,y0),(x1,y1) is a shortest path 
connecting them in the tubular topology. This is either the direct
line segment within x range [0,width], or the linear path wrapping around
the vertical edge. We can check which is shorter by comparing abs(x0-x1) to
width/2. 

If abs(x0-x1) is exactly width/2, then both paths must be creases. However, 
they may not have the same mountain-valley assignment, so they have to be 
encoded separately, distingushed by their ordering of the vertices.
In this case, we say the edge where x0<x1 is the direct path and the one where
x0>x1 is the wraparound path.

(To see why both paths must be creases if abs(x0-x1) is exactly width/2, note 
that in this case the paths have the same length d. If one is a crease, it maps 
to a straight segment in the folded 3D structure, so the vertices map to 
exactly d apart, since folding is an isometry. Since the other path must also
map to a path of length d connecting these vertices in 3D which are d apart,
it must also map to the segment connecting them, i.e., it is also a crease.)

To make this manufacturable from a sheet (e.g. by a laser cutter),
we consider it to be cut along the x=0 axis, with the geometry in 
x=[0,baseSideLength] duplicated at the end so that the wraparound can adhere
together along a 2D surface. This duplication, along with other difficulties of
converting wraparound edges to a flat pattern, is dealt with in the method
that constructs a DXF file. (We organize in this way to separate fabrication 
details from the underlying pattern: the construction and representation of 
TubularPattern objects is in the tubular topology).

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
    signedAngles2D, signedAngles3D, lerpN, unit, signedAngle2D, unsignedAngle
from spatialmath import SO3, SE3
import copy
from style import *

class TubularPattern():
    def __init__(self, numSides, r, proximalMarker=[0, 0]):
        assert (r > 0)
        assert (numSides % 2 == 0)
        assert (numSides >= 4)
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
        self.EPSILON = self.baseSideLength * 0.0001

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
        self.assertValidationChecks()

    def reverse(self):
        self.Vertices[:,1] = self.patternHeight - self.Vertices[:,1] + self.proximalMarker[:,1]
        self.proximalBaseIndices, self.distalBaseIndices = self.distalBaseIndices, self.proximalBaseIndices
        self.proximalMarker[0,0], self.distalMarker[0,0] = self.distalMarker[0,0], self.proximalMarker[0,0]
        return self

    """ Add the given vertices (a numpy nx2), return their indices """    
    def addVertices(self, VerticesToAdd):
        newIndices = self.Vertices.shape[0] + np.arange(VerticesToAdd.shape[0])
        self.Vertices = np.vstack((self.Vertices, VerticesToAdd))
        return newIndices

    """ Add the given edges, with the given mountain/valley labels, to the 
        crease pattern. Filters our self-edges (i,i). If deDuplicate is set
        to True, also filters out duplicate edges. """
    def addEdges(self, EdgesToAdd, LabelsToAdd, deDuplicate=False):
        assert(EdgesToAdd.shape[0] == LabelsToAdd.shape[0])
        assert(EdgesToAdd.shape[1] == 2)
        assert(EdgesToAdd.dtype == int)
        nonEmpty = EdgesToAdd[:,0] != EdgesToAdd[:,1]
        self.Edges = np.vstack((self.Edges, EdgesToAdd[nonEmpty]))
        self.EdgeLabelsMV = np.hstack((self.EdgeLabelsMV, 
                                       LabelsToAdd[nonEmpty]))
        
        if deDuplicate:
            self.Edges, indices = np.unique(self.Edges, return_index=True, 
                                            axis=0)
            self.EdgeLabelsMV = self.EdgeLabelsMV[indices]

    def addMountainEdges(self, EdgesToAdd, deDuplicate=False):
        EdgesToAdd = np.array(EdgesToAdd)
        assert(EdgesToAdd.shape[0]>0)
        EdgesToAdd = np.array(EdgesToAdd)
        newTrues = np.ones(EdgesToAdd.shape[0], dtype=bool)
        self.addEdges(EdgesToAdd, newTrues, deDuplicate)
    
    def addValleyEdges(self, EdgesToAdd, deDuplicate=False):
        EdgesToAdd = np.array(EdgesToAdd)
        assert(EdgesToAdd.shape[0]>0)
        EdgesToAdd = np.array(EdgesToAdd)
        newFalses = np.zeros(EdgesToAdd.shape[0], dtype=bool)
        self.addEdges(EdgesToAdd, newFalses, deDuplicate)

    def wrapToWidth(self):
        # For numerical stability, need to ensure that x coordinates 
        # within EPSILON of 0 or width are treated as 0 or width
        basicallyZero = abs(self.Vertices[:, 0]) < self.EPSILON
        self.Vertices[basicallyZero, 0] = 0
        basicallyWidth = abs(self.Vertices[:, 0] - self.width) < self.EPSILON
        self.Vertices[basicallyWidth, 0] = self.width
        
        self.Vertices[:, 0] %= self.width
        self.proximalMarker[0, 0] %= self.width
        self.distalMarker[0, 0] %= self.width

    def translate(self, vector):
        self.Vertices += vector
        self.proximalMarker = self.proximalMarker + vector
        self.distalMarker = self.distalMarker + vector
        self.wrapToWidth()

    def ProximalBase(self):
        return self.Vertices[self.proximalBaseIndices, :]

    def DistalBase(self):
        return self.Vertices[self.distalBaseIndices, :]
    
    def isEmpty(self):
        return self.Vertices.shape[0] == self.numSides
    
    def hasProximalBase(self):
        return not self.proximalBaseIndices is None
    
    def hasDistalBase(self):
        return not self.distalBaseIndices is None

    def rollProximalIndices(self, k):
        PI = self.proximalBaseIndices
        self.Vertices[PI, :] = np.roll(self.Vertices[PI, :], k, axis=0)
        Mask = self.Edges < self.numSides
        self.Edges[Mask] = (self.Edges[Mask] + k) % self.numSides
        self.proximalBaseIndices = np.roll(self.proximalBaseIndices, k)
    
    def assertValidationChecks(self):
        if self.hasDistalBase():
            DistalBase = self.DistalBase()
            assert(np.all(DistalBase[:,1]==DistalBase[0,1]))
        if self.hasProximalBase():
            ProximalBase = self.ProximalBase()
            assert(np.all(ProximalBase[:,1]==ProximalBase[0,1]))
        assert(np.all(self.Vertices[:, 0] >= 0))
        assert(np.all(self.Vertices[:, 0] <= self.width))
    
    def append(self, other: 'TubularPattern'):
        assert(self.hasDistalBase() and other.hasProximalBase())
        self.assertValidationChecks()
        other.assertValidationChecks()
        other = copy.deepcopy(other)
        assert (other.numSides == self.numSides)
        assert (other.r == self.r)
        assert (self.Edges.dtype.kind == 'i')
        other.translate(self.distalMarker - other.proximalMarker)
        
        if not other.isEmpty():
            """ Merge other's proximal base vertices with current distal base """
            indexShift = self.Vertices.shape[0] - other.numSides
            otherNonProximalVertices = np.delete(other.Vertices,
                                                 other.proximalBaseIndices, axis=0)
            self.Vertices = np.vstack((self.Vertices, otherNonProximalVertices))
    
            """ Reindex other's proximal base vertices such that they match
                the current distal base vertices they correspond to """
            k = rollToMatch(other.ProximalBase(), self.DistalBase(),
                            EPSILON=self.EPSILON)
            if k is None:
                print(self.ProximalBase())
                print(self.DistalBase())
                print(other.ProximalBase())
                print("ERROR: can't align self.DistalBase() with other.ProximalBase()")
            
            assert (not k is None)
            other.rollProximalIndices(k)
            EdgesToAdd = other.Edges
            EdgesToAdd[EdgesToAdd>=self.numSides] += indexShift #non-proximal
            EdgesToAdd[EdgesToAdd<self.numSides] += np.min(self.distalBaseIndices)
            self.addEdges(EdgesToAdd, other.EdgeLabelsMV, deDuplicate=True)
            
            if other.hasDistalBase():
                self.distalBaseIndices = other.distalBaseIndices + indexShift
            else:
                self.distalBaseIndices = None
            self.patternHeight += other.patternHeight
            self.distalMarker = other.distalMarker
            self.assertValidationChecks()
        
        return self
    
    """ Plot the graph directly, without accounting for wraparound, boundary,
        or other considerations for actual manufacturing.
        Mainly useful for debugging - comparing it to the result of makeDXF
        helps separate mistakes in the actual graph from bugs in makeDXF. """
    def showRawGraph(self, saveas="", show=True, directed=False):
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
            plt.show(block=False)

        return doc

    def save(self, dxfName, mountainColor=mountainColorDefault, 
             valleyColor=valleyColorDefault, boundaryColor=boundaryColorDefault, 
             cutColor=cutColorDefault):
        return self.show(dxfName=dxfName, show=False, printGraph=False,
                         block=False, mountainColor=mountainColor,
                         valleyColor=valleyColor, boundaryColor=boundaryColor,
                         cutColor=cutColor) 

    def show(self, dxfName="", show=True, printGraph=False, block=blockDefault,
             mountainColor=mountainColorDefault, valleyColor=valleyColorDefault,
             boundaryColor=boundaryColorDefault, cutColor=cutColorDefault):
        if len(dxfName) < 4 or dxfName[-4:]!=".dxf":
            dxfName += ".dxf"
        
        doc = ezdxf.new()

        # add new entities to the modelspace
        msp = doc.modelspace()
        doc.layers.add(name="Mountain", color=mountainColor)
        doc.layers.add(name="Valley", color=valleyColor)
        doc.layers.add(name="Boundary", color=boundaryColor)
        doc.layers.add(name="Cut", color=cutColor)

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
        if printGraph:
            print("Vertices:")
            print(self.Vertices)
            print("Edges:")
            print(self.Edges)
            print("Folds:")
            print(Folds)
        for i in range(self.Edges.shape[0]):
            seg = Folds[i,:,:]
            # The normal wraparound case is abs(x0-x1) > width/2.
            # If abs(x0-x1)==width/2, both directions should exist as separate
            # edges, and the one with x0>x1 is the wraparound version.
            # Note that (abs(x0-x1)==width/2 and x0>x1) iff x0-x1==width/2.
            # We implement both of these checks within epsilon bounds for 
            # numerical stability.
            dx = seg[0,0] - seg[1,0] # x0 - x1
            isWrapping = (abs(dx) > (self.width/2 - self.EPSILON)) or \
                            (abs(dx - self.width/2) < self.EPSILON)
            layer = "Mountain" if self.EdgeLabelsMV[i] else "Valley"
            
            if isWrapping:
                rightcopy = np.copy(seg)
                leftcopy = np.copy(seg)
                #rightcopy[seg[:, 0] <= self.baseSideLength, 0] += self.width
                rightcopy[np.argmin(seg[:,0]), 0] += self.width
                msp.add_line(rightcopy[0, :], rightcopy[1, :],
                             dxfattribs={"layer": layer})
                #leftcopy[seg[:, 0] >= lastColStart, 0] -= self.width
                leftcopy[np.argmax(seg[:,0]), 0] -= self.width
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
            
        if dxfName:
            doc.saveas(dxfName)

        if show:
            config = Configuration()
            fig: plt.Figure = plt.figure(figsize=(xmax-xmin, ymax-ymin))
            ax: plt.Axes = fig.add_axes([0, 0, 1, 1])
            ctx = RenderContext(doc)
            out = MatplotlibBackend(ax)
            # get the modelspace properties
            msp_properties = LayoutProperties.from_layout(msp)
            # set light gray background color and black foreground color
            msp_properties.set_colors("#ffffff") # light gray #eaeaea
            Frontend(ctx, out, config=config).draw_layout(msp, finalize=False,
                                             layout_properties=msp_properties)
            ax.set_ylim(ymin, ymax)
            ax.set_xlim(xmin, xmax)
            plt.show(block=block)

        return doc

class TubeFittingPattern(TubularPattern):
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
        self.assertValidationChecks()

class TipPattern(TubularPattern):
    def __init__(self, numSides, r, length, forward=True, proximalMarker=[0, 0]):
        super().__init__(numSides, r, proximalMarker)
        assert(length > 0)
        totalBendingAngle = 4 * np.arctan(length / r*np.sin(self.polygonInnerAngle))
        self.patternHeight = length / np.sin(totalBendingAngle/4)
        
        self.patternHeight = self.patternHeight
        ProximalBase = self.Vertices
        
        self.distalBaseIndices = None #Can't append anything after a Tip
        
        MidAlignedVertices = ProximalBase + [0, self.patternHeight]
        midAlignedIndices = self.Vertices.shape[0] + np.arange(numSides)
        self.Vertices = np.vstack((self.Vertices, MidAlignedVertices))
        MidHalfVertices = MidAlignedVertices + [self.baseSideLength/2, 0]
        midHalfIndices = self.Vertices.shape[0] + np.arange(numSides)
        self.Vertices = np.vstack((self.Vertices, MidHalfVertices))
        
        # Base mountain folds
        self.addMountainEdges(np.vstack((self.proximalBaseIndices, 
                                np.roll(self.proximalBaseIndices, -1))).T)
        
        # Folds of the flapping rectangles
        
        ep = numSides-1 # end panel index
        mp = (numSides-1)//2 # mid panel index
        nonSinkPanelIndices = [mp, ep]
        
        self.addMountainEdges([
                    [self.proximalBaseIndices[0], midAlignedIndices[0]],
                    [midAlignedIndices[0], midHalfIndices[0]],
                    #######
                    [self.proximalBaseIndices[mp], midAlignedIndices[mp]],
                    [midHalfIndices[mp-1], midAlignedIndices[mp]],
                    #######
                    [self.proximalBaseIndices[mp+1], midAlignedIndices[mp+1]],
                    [midAlignedIndices[mp+1], midHalfIndices[mp+1]],
                    #######
                    [self.proximalBaseIndices[ep], midAlignedIndices[ep]],
                    [midHalfIndices[ep-1], midAlignedIndices[ep]]
                ])
        self.addValleyEdges([    [midAlignedIndices[mp], midHalfIndices[mp]],
                                 [midHalfIndices[mp], midAlignedIndices[mp+1]],
                                 [midAlignedIndices[ep], midHalfIndices[ep]],
                                 [midHalfIndices[ep], midAlignedIndices[0]]  ])
        
        psi = np.arctan2(self.baseSideLength, self.patternHeight)
        deltas = np.pi * (numSides + 2 - 4*(np.arange(numSides)+1)) / (2*numSides)
        frontLengths = np.sqrt(length**2 + (r*np.sin(deltas))**2)
        gammas = np.arctan(frontLengths / (r*np.abs(np.cos(deltas))))
        tuckYs = frontLengths / np.cos((np.pi/2) - gammas - psi)
        
        lowerCrimpIndices = {}
        
        for i in range(numSides):
            # Crimping vertices and folds
            if not i in [0, mp, mp+1, ep]:
                lowerCrimpVertex = ProximalBase[i] + [0,tuckYs[i]]
                lowerCrimpIndices[i] = self.Vertices.shape[0]
                self.Vertices = np.vstack((self.Vertices, 
                                          [lowerCrimpVertex]))
                self.addValleyEdges([[lowerCrimpIndices[i], midAlignedIndices[i]]])
                self.addMountainEdges([[self.proximalBaseIndices[i], lowerCrimpIndices[i]]])

        for i in range(numSides):    
            # Check if this is a sink panel
            # (rather than one of the 2 side flaps that stay rectangular)
            if not i in nonSinkPanelIndices: # Tucking X folds
                mhi = midHalfIndices[i]
                self.addValleyEdges([[self.proximalBaseIndices[i], mhi],
                                     [mhi, self.proximalBaseIndices[i+1]]])    
            if not i in [0, mp, mp+1, ep]: 
                self.addValleyEdges([[midHalfIndices[i-1], midAlignedIndices[i]],
                                     [midAlignedIndices[i], midHalfIndices[i]]])
                self.addMountainEdges([[midHalfIndices[i-1], lowerCrimpIndices[i]],
                                       [lowerCrimpIndices[i], midHalfIndices[i]]])
        self.distalMarker = self.proximalMarker + [0, self.patternHeight]
        self.wrapToWidth()
        self.assertValidationChecks()

        if not forward:
            self.reverse()
    
class RevoluteJointPattern(TubularPattern):
    def __init__(self, numSides, r, totalBendingAngle, numSinkLayers=1, 
                 proximalMarker=[0, 0]):
        assert(numSinkLayers >= 1)
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
        
        # store for Tip construction
        self.midAlignedIndices = midAlignedIndices 
        self.flatRotAxisHeight = flatRotAxisHeight 
        
        # Base mountain folds
        self.addMountainEdges(np.vstack((self.proximalBaseIndices, 
                                np.roll(self.proximalBaseIndices, -1))).T)
        self.addMountainEdges(np.vstack((self.distalBaseIndices,
                                np.roll(self.distalBaseIndices, -1))).T)
        
        # Folds of the flapping rectangles
        
        ep = numSides-1 # end panel index
        mp = (numSides-1)//2 # mid panel index
        nonSinkPanelIndices = [mp, ep]
        
        self.addMountainEdges([
                    [self.proximalBaseIndices[0], midAlignedIndices[0]],
                    [midAlignedIndices[0], midHalfIndices[0]],
                    [midAlignedIndices[0], self.distalBaseIndices[0]], 
                    #######
                    [self.proximalBaseIndices[mp], midAlignedIndices[mp]],
                    [midHalfIndices[mp-1], midAlignedIndices[mp]],
                    [midAlignedIndices[mp], self.distalBaseIndices[mp]], 
                    #######
                    [self.proximalBaseIndices[mp+1], midAlignedIndices[mp+1]],
                    [midAlignedIndices[mp+1], midHalfIndices[mp+1]],
                    [midAlignedIndices[mp+1], self.distalBaseIndices[mp+1]],
                    #######
                    [self.proximalBaseIndices[ep], midAlignedIndices[ep]],
                    [midHalfIndices[ep-1], midAlignedIndices[ep]],
                    [midAlignedIndices[ep], self.distalBaseIndices[ep]]
                ])
        self.addValleyEdges([    [midAlignedIndices[mp], midHalfIndices[mp]],
                                 [midHalfIndices[mp], midAlignedIndices[mp+1]],
                                 [midAlignedIndices[ep], midHalfIndices[ep]],
                                 [midHalfIndices[ep], midAlignedIndices[0]]  ])
        
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
        
        lowerCrimpIndices = {}
        upperCrimpIndices = {}
        
        for i in range(numSides):
            # Crimping vertices and folds
            if not i in [0, mp, mp+1, ep]:
                lowerCrimpVertex = ProximalBase[i] + [0,tuckYs[i]]
                lowerCrimpIndices[i] = self.Vertices.shape[0]
                upperCrimpVertex = DistalBase[i] - [0,tuckYs[i]]
                upperCrimpIndices[i] = self.Vertices.shape[0]+1
                self.Vertices = np.vstack((self.Vertices, 
                                          [lowerCrimpVertex,upperCrimpVertex]))
                self.addValleyEdges([[lowerCrimpIndices[i], midAlignedIndices[i]]])
                self.addMountainEdges([
                    [midAlignedIndices[i], upperCrimpIndices[i]],
                    [self.proximalBaseIndices[i], lowerCrimpIndices[i]],
                    [upperCrimpIndices[i], self.distalBaseIndices[i]]
                    ])

        if numSinkLayers==1:                
            for i in range(numSides):    
                # Check if this is a sink panel
                # (rather than one of the 2 side flaps that stay rectangular)
                if not i in nonSinkPanelIndices: # Tucking X folds
                    mhi = midHalfIndices[i]
                    self.addValleyEdges([[self.proximalBaseIndices[i], mhi],
                                         [self.distalBaseIndices[i], mhi],
                                         [mhi, self.proximalBaseIndices[i+1]],
                                         [mhi, self.distalBaseIndices[i+1]]])    
                if not i in [0, mp, mp+1, ep]:
                    self.addValleyEdges([
                            [midHalfIndices[i-1], upperCrimpIndices[i]],
                            [midHalfIndices[i-1], midAlignedIndices[i]],
                            [midAlignedIndices[i], midHalfIndices[i]] ])
                    self.addMountainEdges([
                            [midHalfIndices[i-1], lowerCrimpIndices[i]],
                            [lowerCrimpIndices[i], midHalfIndices[i]],
                            [upperCrimpIndices[i], midHalfIndices[i]]   ])
        else: #RECURSIVE SINK GADGET
            for p in range(numSides): # panel index
                if not p in nonSinkPanelIndices:        
                    bottomLeftIndex = self.proximalBaseIndices[p]
                    bottomRightIndex = self.proximalBaseIndices[p+1]
                    midLeftIndex = midAlignedIndices[p]
                    midRightIndex = midAlignedIndices[p+1]
                    topLeftIndex = self.distalBaseIndices[p]
                    topRightIndex = self.distalBaseIndices[p+1]
                    centerIndex = midHalfIndices[p]
                    center = self.Vertices[centerIndex]
                    
                    left, bottom = self.Vertices[bottomLeftIndex]
                    right, top = self.Vertices[topRightIndex]
                    
                    pLeft = (p-1)%numSides
                    pRight = (p+1)%numSides
                    crimpLeft = not pLeft in nonSinkPanelIndices
                    crimpRight = not pRight in nonSinkPanelIndices
                    
                    
                    
                    # Helper function, returns indices used
                    def addInteriorVertices(startIndex):
                        start = self.Vertices[startIndex]
                        end = self.Vertices[centerIndex]
                        newIndices = self.addVertices(
                            lerpN(start, end, numSinkLayers+1)[1:-1,:])
                        return np.hstack(([startIndex],newIndices,[centerIndex]))
                   
                    
                    def addInteriorCrimpVertices(crimpIndex):
                        crimp = self.Vertices[crimpIndex]
                        d = self.baseSideLength / (2*numSinkLayers)
                        interiorStart = center + (numSinkLayers-1)*d*unit(crimp-center)
                        newIndices = self.addVertices(lerpN(interiorStart, center, 
                                                      numSinkLayers)[:-1,:])
                        return interiorStart, np.hstack(([crimpIndex],newIndices,[centerIndex]))
                    
                    def addInteriorDiagVertices(crimpInteriorStart, outerIndex):
                        outer = self.Vertices[outerIndex]
                        dirFromCenter = unit(outer - center)
                        crimpInteriorFromCenter = crimpInteriorStart - center 
                        mag = np.linalg.norm(crimpInteriorFromCenter) /\
                                    np.cos(unsignedAngle(dirFromCenter, 
                                                    crimpInteriorFromCenter))
                        interiorStart = center + mag*dirFromCenter
                        newIndices = self.addVertices(lerpN(interiorStart, 
                                                center, numSinkLayers)[:-1,:])
                        return np.hstack(([outerIndex],newIndices,[centerIndex]))
                    
                    def EdgeSequence(vertexIndices, wrap):
                        Edges = np.vstack((vertexIndices, 
                                           np.roll(vertexIndices,-1))).T
                        return Edges if wrap else Edges[:-1]
                    
                    def EdgesFromVertexIndexSequences(VertexIndexSequenceRows):
                        Next = np.roll(VertexIndexSequenceRows, -1, axis=1)
                        Stacked = np.stack((VertexIndexSequenceRows, Next))
                        return Stacked[:,:,:-1].reshape((2,-1)).T
                    
                    def addFoldSequence(vertexIndices, startWithValley=True,
                                            alternating=True, flipEdges=False):
                        if flipEdges:
                            Edges = np.vstack((np.roll(vertexIndices,-1),
                                               vertexIndices)).T[:-1]
                        else:
                            Edges = np.vstack((vertexIndices, 
                                           np.roll(vertexIndices,-1))).T[:-1]
                        
                        if alternating:
                            Labels = np.arange(vertexIndices.shape[0]-1)%2 == startWithValley
                        else:
                            Labels = np.zeros(vertexIndices.shape[0]-1) == startWithValley
                        
                        self.addEdges(Edges, Labels)
                        
                    
                    def addInnerLayersEdges(TopIndices, BottomIndices, 
                                        diamondLeft=False, diamondRight=False):
                        # RING FOLDS
                        ringLayers = np.arange(numSinkLayers-1)
                        layerIsMountain = ringLayers % 2 == 1 # V, M, V, M, ...
                        layerIsValley = np.logical_not(layerIsMountain)
                        
                        # When we have diamonds on a side, the bottom edge 
                        # within the diamond has M/V assignment opposite that
                        # of its layer.
                        # This edge is the first (for a left diamond) or last
                        # (for a right diamond) in the bottom part of the layer
                        
                        if np.any(layerIsValley):
                            self.addValleyEdges(EdgesFromVertexIndexSequences(
                                TopIndices[layerIsValley, :]))
                            if diamondLeft and diamondRight:
                                self.addValleyEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsValley, 1:-1]))
                                self.addMountainEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsValley, :2]))
                                self.addMountainEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsValley, -2:]))
                            elif diamondRight:
                                self.addValleyEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsValley, :-1]))
                                self.addMountainEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsValley, -2:]))
                            elif diamondLeft:
                                self.addValleyEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsValley, 1:]))
                                self.addMountainEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsValley, :2]))
                            else:
                                self.addValleyEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsValley, :]))
                                
                        if np.any(layerIsMountain):
                            self.addMountainEdges(EdgesFromVertexIndexSequences(
                                TopIndices[layerIsMountain, :]))
                            if diamondLeft and diamondRight:
                                self.addMountainEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsMountain, 1:-1]))
                                self.addValleyEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsMountain, :2]))
                                self.addValleyEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsMountain, -2:]))
                            elif diamondRight:
                                self.addMountainEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsMountain, :-1]))
                                self.addValleyEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsMountain, -2:]))
                            elif diamondLeft:
                                self.addMountainEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsMountain, 1:]))
                                self.addValleyEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsMountain, :2]))
                            else:
                                self.addMountainEdges(EdgesFromVertexIndexSequences(
                                    BottomIndices[layerIsMountain, :]))
                    
                    if not crimpLeft and not crimpRight: # crimp neither side
                        leftIndices = addInteriorVertices(midLeftIndex)
                        rightIndices = addInteriorVertices(midRightIndex)
                        topLeftDiagIndices = addInteriorVertices(topLeftIndex)
                        topRightDiagIndices = addInteriorVertices(topRightIndex)
                        bottomLeftDiagIndices = addInteriorVertices(bottomLeftIndex)
                        bottomRightDiagIndices = addInteriorVertices(bottomRightIndex)
                        
                        # INWARDS-OUTWARDS FOLDS
                        # the index lists are oriented inwards,
                        # so the ones on the rigth need to flip to have
                        # folds oriented in increasing x
                        addFoldSequence(leftIndices, startWithValley=False)
                        addFoldSequence(rightIndices, flipEdges=True, alternating=False)
                        addFoldSequence(topLeftDiagIndices)
                        addFoldSequence(topRightDiagIndices, flipEdges=True)
                        addFoldSequence(bottomLeftDiagIndices)
                        addFoldSequence(bottomRightDiagIndices, flipEdges=True)
                        
                        # Each row is the indices of that recursive sink layer,
                        # counting from outwards in.
                        # It's split into top and bottom components to ensure
                        # they're oriented left-to-right.
                        InnerLayersTopIndices = np.vstack((leftIndices,
                                                    topLeftDiagIndices,
                                                    topRightDiagIndices,
                                                    rightIndices)).T[1:-1]
                        InnerLayersBottomIndices = np.vstack((leftIndices,
                                                    bottomLeftDiagIndices,
                                                    bottomRightDiagIndices,
                                                    rightIndices)).T[1:-1]
                        
                        addInnerLayersEdges(InnerLayersTopIndices, 
                                            InnerLayersBottomIndices)

                    elif not crimpLeft:
                        # crimp right only
                        leftIndices = addInteriorVertices(midLeftIndex)
                        topLeftDiagIndices = addInteriorVertices(topLeftIndex)
                        bottomLeftDiagIndices = addInteriorVertices(bottomLeftIndex)
                        
                        upperRightCrimpInteriorStart, upperRightCrimpDiagIndices = \
                            addInteriorCrimpVertices(upperCrimpIndices[pRight])
                        topRightDiagIndices = addInteriorDiagVertices(
                            upperRightCrimpInteriorStart, topRightIndex)
                        
                        lowerRightCrimpInteriorStart, lowerRightCrimpDiagIndices = \
                            addInteriorCrimpVertices(lowerCrimpIndices[pRight])
                        bottomRightDiagIndices = addInteriorDiagVertices(
                            lowerRightCrimpInteriorStart, bottomRightIndex)
                        
                        rightIndices = addInteriorDiagVertices(
                            lowerRightCrimpInteriorStart, midRightIndex)
                        
                        # INWARDS-OUTWARDS FOLDS
                        # the index lists are oriented inwards,
                        # so the ones on the rigth need to flip to have
                        # folds oriented in increasing x
                        addFoldSequence(leftIndices, startWithValley=False)
                        addFoldSequence(topLeftDiagIndices)
                        addFoldSequence(topRightDiagIndices, flipEdges=True)
                        addFoldSequence(upperRightCrimpDiagIndices, flipEdges=True)
                        addFoldSequence(rightIndices, flipEdges=True, alternating=False)
                        addFoldSequence(lowerRightCrimpDiagIndices, 
                                        startWithValley=False, flipEdges=True)
                        addFoldSequence(bottomRightDiagIndices, flipEdges=True)
                        addFoldSequence(bottomLeftDiagIndices)
                        
                        # Each row is the indices of that recursive sink layer,
                        # counting from outwards in.
                        # It's split into top and bottom components to ensure
                        # they're oriented left-to-right.
                        InnerLayersTopIndices = np.vstack((leftIndices,
                                                    topLeftDiagIndices,
                                                    topRightDiagIndices,
                                                    upperRightCrimpDiagIndices,
                                                    rightIndices)).T[1:-1]
                        InnerLayersBottomIndices = np.vstack((leftIndices,
                                                    bottomLeftDiagIndices,
                                                    bottomRightDiagIndices,
                                                    lowerRightCrimpDiagIndices,
                                                    rightIndices)).T[1:-1]
                        
                        addInnerLayersEdges(InnerLayersTopIndices, 
                                            InnerLayersBottomIndices, 
                                            diamondRight=True)
                        
                    
                
                    elif not crimpRight:
                        # crimp left only
                        rightIndices = addInteriorVertices(midRightIndex)
                        topRightDiagIndices = addInteriorVertices(topRightIndex)
                        bottomRightDiagIndices = addInteriorVertices(bottomRightIndex)
                        
                        upperLeftCrimpInteriorStart, upperLeftCrimpDiagIndices = \
                            addInteriorCrimpVertices(upperCrimpIndices[p])
                        topLeftDiagIndices = addInteriorDiagVertices(
                            upperLeftCrimpInteriorStart, topLeftIndex)
                        
                        lowerLeftCrimpInteriorStart, lowerLeftCrimpDiagIndices = \
                            addInteriorCrimpVertices(lowerCrimpIndices[p])
                        bottomLeftDiagIndices = addInteriorDiagVertices(
                            lowerLeftCrimpInteriorStart, bottomLeftIndex)
                        
                        leftIndices = addInteriorDiagVertices(
                            lowerLeftCrimpInteriorStart, midLeftIndex)
                        
                        # INWARDS-OUTWARDS FOLDS
                        # the index lists are oriented inwards,
                        # so the ones on the rigth need to flip to have
                        # folds oriented in increasing x
                        addFoldSequence(leftIndices, alternating=False)
                        addFoldSequence(upperLeftCrimpDiagIndices, startWithValley=False)
                        addFoldSequence(topLeftDiagIndices)
                        addFoldSequence(topRightDiagIndices, flipEdges=True)
                        addFoldSequence(rightIndices, flipEdges=True, startWithValley=False)
                        addFoldSequence(bottomRightDiagIndices, flipEdges=True)
                        addFoldSequence(bottomLeftDiagIndices)
                        addFoldSequence(lowerLeftCrimpDiagIndices, startWithValley=False)
                        
                        # Each row is the indices of that recursive sink layer,
                        # counting from outwards in.
                        # It's split into top and bottom components to ensure
                        # they're oriented left-to-right.
                        InnerLayersTopIndices = np.vstack((leftIndices,
                                                    upperLeftCrimpDiagIndices,
                                                    topLeftDiagIndices,
                                                    topRightDiagIndices,
                                                    rightIndices)).T[1:-1]
                        InnerLayersBottomIndices = np.vstack((leftIndices,
                                                    lowerLeftCrimpDiagIndices,
                                                    bottomLeftDiagIndices,
                                                    bottomRightDiagIndices,
                                                    rightIndices)).T[1:-1]
                        
                        addInnerLayersEdges(InnerLayersTopIndices, 
                                            InnerLayersBottomIndices,
                                            diamondLeft=True)
                        
                    
                    else:
                        upperLeftCrimpInteriorStart, upperLeftCrimpDiagIndices = \
                            addInteriorCrimpVertices(upperCrimpIndices[p])
                        topLeftDiagIndices = addInteriorDiagVertices(
                            upperLeftCrimpInteriorStart, topLeftIndex)
                        
                        lowerLeftCrimpInteriorStart, lowerLeftCrimpDiagIndices = \
                            addInteriorCrimpVertices(lowerCrimpIndices[p])
                        bottomLeftDiagIndices = addInteriorDiagVertices(
                            lowerLeftCrimpInteriorStart, bottomLeftIndex)
                        
                        leftIndices = addInteriorDiagVertices(
                            lowerLeftCrimpInteriorStart, midLeftIndex)
                        
                        # crimp on both sides
                        upperRightCrimpInteriorStart, upperRightCrimpDiagIndices = \
                            addInteriorCrimpVertices(upperCrimpIndices[pRight])
                        topRightDiagIndices = addInteriorDiagVertices(
                            upperRightCrimpInteriorStart, topRightIndex)
                        
                        lowerRightCrimpInteriorStart, lowerRightCrimpDiagIndices = \
                            addInteriorCrimpVertices(lowerCrimpIndices[pRight])
                        bottomRightDiagIndices = addInteriorDiagVertices(
                            lowerRightCrimpInteriorStart, bottomRightIndex)
                        
                        rightIndices = addInteriorDiagVertices(
                            lowerRightCrimpInteriorStart, midRightIndex)
                        
                        # INWARDS-OUTWARDS FOLDS
                        addFoldSequence(bottomLeftDiagIndices)
                        addFoldSequence(lowerLeftCrimpDiagIndices, startWithValley=False)
                        addFoldSequence(leftIndices, alternating=False)
                        addFoldSequence(upperLeftCrimpDiagIndices, startWithValley=False)
                        addFoldSequence(topLeftDiagIndices)                        
                        # the index lists are oriented inwards,
                        # so the ones on the right need to flip to have
                        # folds oriented in increasing x
                        addFoldSequence(topRightDiagIndices, flipEdges=True)
                        addFoldSequence(upperRightCrimpDiagIndices, flipEdges=True)
                        addFoldSequence(rightIndices, flipEdges=True, alternating=False)
                        addFoldSequence(lowerRightCrimpDiagIndices, 
                                        startWithValley=False, flipEdges=True)
                        addFoldSequence(bottomRightDiagIndices, flipEdges=True)
                        
                        # Each row is the indices of that recursive sink layer,
                        # counting from outwards in.
                        # It's split into top and bottom components to ensure
                        # they're oriented left-to-right.
                        InnerLayersTopIndices = np.vstack((leftIndices,
                                                    upperLeftCrimpDiagIndices,
                                                    topLeftDiagIndices,
                                                    topRightDiagIndices,
                                                    upperRightCrimpDiagIndices,
                                                    rightIndices)).T[1:-1]
                        InnerLayersBottomIndices = np.vstack((leftIndices,
                                                    lowerLeftCrimpDiagIndices,
                                                    bottomLeftDiagIndices,
                                                    bottomRightDiagIndices,
                                                    lowerRightCrimpDiagIndices,
                                                    rightIndices)).T[1:-1]
                        
                        addInnerLayersEdges(InnerLayersTopIndices, 
                                            InnerLayersBottomIndices,
                                            diamondLeft=True, 
                                            diamondRight=True)
        
        self.distalMarker = self.proximalMarker + [0, self.patternHeight]
        self.wrapToWidth()
        self.assertValidationChecks()

class ReboPattern(TubularPattern):
    def __init__(self, numSides, r, neutralLength, numLayers, coneAngle,
                 proximalValleys = False, proximalMarker=[0, 0]):
        super().__init__(numSides, r, proximalMarker)
        assert(coneAngle > 0 and coneAngle < np.pi/2)
        assert(numLayers > 0)
        neutralLayerHeight = neutralLength / numLayers
        flatLayerHalfHeight = neutralLayerHeight / (2*np.sin(coneAngle))
        maxLength = 2*numLayers*flatLayerHalfHeight
        self.patternHeight = maxLength
        flatConeAngle = (np.pi/2)*(1 - 2*np.cos(coneAngle)/numSides)
        
        ProximalBase = self.Vertices
        xsTiled = np.tile(ProximalBase[:,0], numLayers+1)
        layerYs = np.arange(numLayers+1)*2*flatLayerHalfHeight
        ysRepeated = np.repeat(layerYs, numSides)
        ReboLayerBoundariesVertices = np.vstack((xsTiled, ysRepeated)).T
        ReboLayerBoundariesIndices = np.arange((numLayers+1)*numSides).reshape((numLayers+1, numSides))
        self.Vertices = ReboLayerBoundariesVertices
        assert(np.linalg.norm(self.Vertices[:numSides,:] - ProximalBase) < self.EPSILON)
        self.distalBaseIndices = ReboLayerBoundariesIndices[-1,:]
        
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
        
        # Construct folds
        for i in range(numLayers):
            bottom = ReboLayerBoundariesIndices[i,:]
            nextBottom = np.roll(bottom, -1)
            midAligned = ReboLayerMidAlignedIndices[i,:]
            midOffset = ReboLayerMidOffsetIndices[i,:]
            nextMidOffset = np.roll(midOffset, -1)
            top = ReboLayerBoundariesIndices[i+1,:]
            
            if proximalValleys and i==0:
                self.addValleyEdges(np.vstack((bottom, nextBottom)).T)
            else:    
                self.addMountainEdges(np.vstack((bottom, nextBottom)).T)
                
            self.addMountainEdges(np.vstack((bottom, midAligned)).T)
            self.addMountainEdges(np.vstack((midAligned, top)).T)
            self.addMountainEdges(np.vstack((midOffset, midAligned)).T)

            self.addValleyEdges(np.vstack((midAligned, nextMidOffset)).T)
            self.addValleyEdges(np.vstack((midOffset, bottom)).T)
            self.addValleyEdges(np.vstack((midOffset, top)).T)
        reboTop = ReboLayerBoundariesIndices[numLayers,:]
        self.addMountainEdges(np.vstack((reboTop, np.roll(reboTop,-1))).T)
            
        self.distalMarker = self.proximalMarker + [0, self.patternHeight]
        self.wrapToWidth()  
        self.assertValidationChecks()

class PrismaticJointPattern(TubularPattern):
    def __init__(self, numSides, r, neutralLength, numLayers, coneAngle,
                 proximalMarker=[0, 0]):
        super().__init__(numSides, r, proximalMarker)
        assert(coneAngle > 0 and coneAngle < np.pi/2)
        assert(numLayers > 0)
        neutralLayerHeight = neutralLength / numLayers
        flatLayerHalfHeight = neutralLayerHeight / (2*np.sin(coneAngle))
        maxLength = 2*numLayers*flatLayerHalfHeight
        
        tube = TubeFittingPattern(numSides, r, height=maxLength)
        self.append(tube)
        mid1Indices = self.distalBaseIndices
        self.addMountainEdges(np.vstack((mid1Indices, 
                                         np.roll(mid1Indices, -1, axis=0))).T)
        valleyTube = copy.deepcopy(tube)
        valleyTube.EdgeLabelsMV = np.logical_not(valleyTube.EdgeLabelsMV)
        self.append(valleyTube)
        self.append(ReboPattern(numSides, r, neutralLength, numLayers, coneAngle, True))
        self.append(tube)
           
class ElbowFittingPattern(TubularPattern):
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
        self.assertValidationChecks()

class TwistFittingPattern(TubularPattern):
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
        self.assertValidationChecks()


class KreslingJointPattern(TubularPattern):
    """ Construction reference:
        Priyanka Bhovad, Joshua Kaufmann, Suyi Li. 
        "Peristaltic locomotion without digital controllers: 
        Exploiting multi-stability in origami to coordinate robotic motion".
        Extreme Mechanics Letters, 2019.
        
        Stiffness scales with angleRatio, which must be in (0,1).
        Pattern is bistable if angleRatio > 0.5, and larger angleRatio
        gives longer extended stable length.
    """
    def __init__(self, numSides, r, compressedStableLength, angleRatio,
                 proximalMarker=[0, 0]):
        assert(angleRatio>0 and angleRatio<1) 
        assert(compressedStableLength > 0)
        super().__init__(numSides, r, proximalMarker)
        gamma = (np.pi/2) - (np.pi/numSides)
        D = 2*r*np.cos(gamma*(1-angleRatio))
        P = self.baseSideLength
        assert(abs(r - 0.5*P/np.sin(np.pi/numSides)) < self.EPSILON) #sanity check that their R matches ours
        L = compressedStableLength
        valleyLength = np.sqrt(D**2 + L**2)
        mountainLength = np.sqrt(P**2 + D**2 - 2*P*D*np.cos(gamma*angleRatio) + L**2)
        valleyAngle = np.arccos((P**2 + valleyLength**2 - mountainLength**2) / (2*P*valleyLength))
        
        valleyDirection = np.array([np.cos(valleyAngle), np.sin(valleyAngle)])
        valleyVector = valleyLength * valleyDirection
        self.patternHeight = valleyVector[1]
        DistalBase = self.ProximalBase() + valleyVector
        self.distalBaseIndices = self.addVertices(DistalBase)
        
        self.addValleyEdges(np.vstack((self.proximalBaseIndices,
                                       self.distalBaseIndices)).T)
        self.addMountainEdges(np.vstack((self.proximalBaseIndices,
                                    np.roll(self.distalBaseIndices, 1))).T)
        
        self.distalMarker = self.proximalMarker + valleyVector
        self.wrapToWidth()
        self.assertValidationChecks()
        
