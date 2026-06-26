# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 15:42:09 2023

@author: Daniel Feshbach
"""
import numpy as np
from numpy import cross, dot, arctan2
from numpy.linalg import norm
from spatialmath import SO3, SE3
from PathCSC import *
from geometryHelpers import *
import matplotlib.pyplot as plt
import TubularPattern
from TubularPattern import TubularPattern, TubeFittingPattern, \
                            ElbowFittingPattern, TwistFittingPattern
from CollisionDetection import *


def _dubins_frame_from_pos_and_dir(pos: np.ndarray, path_dir: np.ndarray) -> SE3:
    """Build a minimal valid SE3 Dubins frame with path_dir as column 0 (x-axis)."""
    xhat = path_dir / norm(path_dir)
    yhat = unitNormalToBoth(xhat, np.array([0., 0., 1.]))
    zhat = np.cross(xhat, yhat)
    return SE3.Rt(SO3(np.column_stack([xhat, yhat, zhat])), pos)


class LinkCSC:
    def __init__(self, r : float, StartDubinsPose : SE3, EndDubinsPose : SE3,
                 maxAnglePerElbow : float = np.pi/2, 
                 path : Optional[PathCSC] = None, EPSILON : float = 0.01,
                 startRadius : Optional[float] = None, endRadius : Optional[float] = None):
        if r <= 0:
            raise ValueError("ERROR: Tried to create a LinkCSC with non-positive radius")
        if not (maxAnglePerElbow >= 0 and maxAnglePerElbow <= np.pi):
            raise ValueError("ERROR: maxAnglePerElbow must be in [0, pi]")
        if not startRadius is None and startRadius <= 0:
            raise ValueError("ERROR: Tried to create a LinkCSC with non-positive startRadius")
        if not endRadius is None and endRadius <= 0:
            raise ValueError("ERROR: Tried to create a LinkCSC with non-positive endRadius")

        self.r = r
        self.startRadius = startRadius if startRadius is not None else r
        self.endRadius = endRadius if endRadius is not None else r
        self.maxAnglePerElbow = maxAnglePerElbow
        self.EPSILON = EPSILON
        self.DISTANCE_EPSILON = self.r * self.EPSILON
        self.StartDubinsPose = StartDubinsPose
        self.EndDubinsPose = EndDubinsPose
        if path is None:
            self.path = shortestCSC(r, self.StartDubinsPose.t, self.StartDubinsPose.R[:,0],
                                self.EndDubinsPose.t, self.EndDubinsPose.R[:,0])
        else:
            self.path = path

        if norm(self.path.error) > self.DISTANCE_EPSILON:
            raise ValueError(f"ERROR: Tried to generate a link for an invalid path\n---\nPath details: {repr(self.path)}")
        if self.path.theta1 < -self.EPSILON:
            raise ValueError(f"ERROR: Tried to generate a link for a path with theta1 < 0\n---\nPath details: {repr(self.path)}")
        if self.path.theta1 >= np.pi:
            raise ValueError(f"ERROR: Tried to generate a link for a path with theta1 >= pi\n---\nPath details: {repr(self.path)}")
        if self.path.theta2 < -self.EPSILON:
            raise ValueError(f"ERROR: Tried to generate a link for a path with theta2 < 0\n---\nPath details: {repr(self.path)}")
        if self.path.theta2 >= np.pi:
            raise ValueError(f"ERROR: Tried to generate a link for a path with theta2 >= pi\n---\nPath details: {repr(self.path)}")

        self.rot1AxisDir = np.cross(self.StartDubinsPose.R[:,0], self.path.w1)
        self.rot1AxisAngle = signedAngle(self.StartDubinsPose.R[:,1],
                                          self.rot1AxisDir,
                                          self.StartDubinsPose.R[:,0])
        self.rot2AxisDir = np.cross(self.EndDubinsPose.R[:,0], self.path.w2)
        self.rot2AxisAngle = signedAngle(self.EndDubinsPose.R[:,1],
                                          self.rot2AxisDir,
                                          self.EndDubinsPose.R[:,0])
        
        
        if self.path.theta1 > self.EPSILON:
            self.elbow1 = CompoundElbow(self.r, self.StartDubinsPose, 
                                       self.path.theta1, self.rot1AxisAngle, 
                                       self.maxAnglePerElbow, self.EPSILON)
            self.Elbow1EndFrame = self.elbow1.EndFrame
            self.elbow1BoundingBall = self.elbow1.boundingBall()
        else:
            self.elbow1 = None
            self.Elbow1EndFrame = self.StartDubinsPose
            self.elbow1BoundingBall = Ball(self.StartDubinsPose.t, self.r)
        assert(norm(self.Elbow1EndFrame.R[:,0] - self.path.tUnit) < self.DISTANCE_EPSILON) # verify alignment
        assert(norm(self.Elbow1EndFrame.t - self.path.turn1end) < self.DISTANCE_EPSILON)
        
        if self.path.tMag > self.DISTANCE_EPSILON:
            self.cylinder = Cylinder(self.r, self.Elbow1EndFrame.t, 
                            self.Elbow1EndFrame.R[:,0], self.path.tMag)
        
        if self.path.theta2 > self.EPSILON:
            Elbow2StartOrientation = SO3.AngleAxis(-self.path.theta2, self.rot2AxisDir) * self.EndDubinsPose.R 
            self.Elbow2StartFrame = SE3.Rt(Elbow2StartOrientation, self.path.turn2start)
            assert(norm(self.Elbow2StartFrame.t - (self.Elbow1EndFrame * SE3.Tx(self.path.tMag)).t) < self.EPSILON)
            assert(norm(self.Elbow2StartFrame.R[:,0] - self.path.tUnit) < self.DISTANCE_EPSILON) # verify alignment
            
            self.elbow2 = CompoundElbow(self.r, self.Elbow2StartFrame, 
                                       self.path.theta2, self.rot2AxisAngle, 
                                       self.maxAnglePerElbow, self.EPSILON)
            assert(norm(self.elbow2.EndFrame - self.EndDubinsPose) < self.DISTANCE_EPSILON)
            self.elbow2BoundingBall = self.elbow2.boundingBall()
        else:
            self.elbow2 = None
            self.Elbow2StartFrame = self.EndDubinsPose
            self.elbow2BoundingBall = Ball(self.EndDubinsPose.t, self.r)

        # Pre-compute lengths and arcs for interpolation
        self.lengthC1 = self.r * self.path.theta1 if self.elbow1 else 0
        self.lengthS = self.path.tMag
        self.lengthC2 = self.r * self.path.theta2 if self.elbow2 else 0
        
        if self.lengthC1 > 0:
            self.arc1 = Arc3D(self.path.circleCenter1, self.StartDubinsPose.t,
                             self.StartDubinsPose.R[:,0], self.path.theta1)
        else:
            self.arc1 = None
            
        if self.lengthC2 > 0:
            self.arc2 = Arc3D(self.path.circleCenter2, self.path.turn2start,
                             self.path.tUnit, self.path.theta2)
        else:
            self.arc2 = None

        self.collisionCapsules = self.getCapsules()

        # GL mesh cache for fast animation (populated by addToWidget)
        self._gl_items: list = []          # cached GLMeshItem / GLLinePlotItem refs
        self._gl_ref_pose: SE3 | None = None  # StartDubinsPose when meshes were built
        self._gl_shape_dirty: bool = True     # True => must rebuild meshes

    # GL cache attributes that hold pyqtgraph GL objects (cannot be pickled)
    _GL_CACHE_KEYS = frozenset({'_gl_items', '_gl_ref_pose', '_gl_shape_dirty'})

    def __getstate__(self):
        """Exclude GL cache from pickling / deepcopy."""
        return {k: v for k, v in self.__dict__.items()
                if k not in LinkCSC._GL_CACHE_KEYS}

    def __setstate__(self, state):
        self.__dict__.update(state)
        self._gl_items = []
        self._gl_ref_pose = None
        self._gl_shape_dirty = True

    def __repr__(self):
        return (
            "LinkCSC("
            f"r={repr(self.r)},"
            f"StartDubinsPose={repr(self.StartDubinsPose)},"
            f"EndDubinsPose={repr(self.EndDubinsPose)},"
            f"maxAnglePerElbow={repr(self.maxAnglePerElbow)},"
            f"path={repr(self.path)},"
            f"EPSILON={repr(self.EPSILON)},"
            ")"
        )
    
    def branchingParameters(self):
        dir1 = self.StartDubinsPose.R[:,0]
        b1 = self.StartDubinsPose.R[:,1]

        twist1 = signedAngle(b1, self.path.w1, dir1)
        if self.path.theta1 < self.EPSILON:
            twist1 = 0

        twist1Deg = np.rad2deg(twist1)
        theta1Deg = np.rad2deg(self.path.theta1)
        y1 = self.path.y1 if self.path.theta1 > self.EPSILON else b1

        twist2 = signedAngle(y1, self.path.y2, self.path.tUnit)
        if self.path.theta2 < self.EPSILON:
            twist2 = 0

        twist2Deg = np.rad2deg(twist2)
        theta2Deg = np.rad2deg(self.path.theta2)
        return [twist1Deg, theta1Deg, self.r, self.path.tMag, twist2Deg, theta2Deg, self.r]

    def newLinkTransformedBy(self, Transformation : SE3):
        newLink = LinkCSC(self.r, Transformation @ self.StartDubinsPose, 
                       Transformation @ self.EndDubinsPose, 
                       maxAnglePerElbow = self.maxAnglePerElbow, 
                       path = self.path.newPathTransformedBy(Transformation),
                       EPSILON = self.EPSILON)
        # Transfer GL cache: the shape hasn't changed, only the pose
        if not getattr(self, '_gl_shape_dirty', True):
            newLink._gl_items = self._gl_items
            newLink._gl_ref_pose = self._gl_ref_pose
            newLink._gl_shape_dirty = False
        return newLink
        
    def addToPlot(self, ax, numSides : int = 32, color : str = linkColorDefault, 
                  alpha : float = 0.5, wireFrame : bool = False, 
                  showFrames : bool = False, showPath : bool = True, 
                  pathColor : str = pathColorDefault,
                  showPathCircles : bool = False, showBoundary : bool = True,
                  showElbowBoundingBalls : bool = False):
        allElbowHandleSets = []
        if showBoundary:
            if not self.elbow1 is None:
                elbow1HandleSets = self.elbow1.addToPlot(ax, numSides, color, 
                                                alpha, wireFrame, showFrames,
                                                showElbowBoundingBalls)
                if showFrames:
                    allElbowHandleSets += elbow1HandleSets
            
            if self.path.tMag > self.DISTANCE_EPSILON:
                self.cylinder.addToPlot(ax, numSides, color, alpha, wireFrame)
            
            if not self.elbow2 is None:
                elbow2HandleSets = self.elbow2.addToPlot(ax, numSides, color, 
                                                    alpha, wireFrame, showFrames,
                                                    showElbowBoundingBalls)
                if showFrames:
                    allElbowHandleSets += elbow2HandleSets
        elif showFrames:
            # show the start and end frames
            startFrameHandles = addPosesToPlot(np.array([self.StartDubinsPose]), ax, 
                                              axisLength=self.r, xColor='darkred', 
                                              yColor='darkblue', zColor='darkgreen')
            endFrameHandles = addPosesToPlot(np.array([self.EndDubinsPose]), ax, 
                                            axisLength=self.r, xColor='darkred', 
                                            yColor='darkblue', zColor='darkgreen')
            allElbowHandleSets.append(startFrameHandles)
            allElbowHandleSets.append(endFrameHandles)

        if showPath:
            self.path.addToPlot(ax, showCircles=showPathCircles, 
                                showPoses=False, pathColor=pathColor)
        
        return allElbowHandleSets
    
    def addToWidget(self, widget, numSides : int = 8, color : str = linkColorDefault, 
                  alpha : float = 0.5, wireFrame : bool = False, 
                  showFrames : bool = False, showPath : bool = True, 
                  pathColor : str = pathColorDefault,
                  showPathCircles : bool = False, showBoundary : bool = True,
                  showElbowBoundingBalls : bool = False, linkID : int = None):
        import pyqtgraph.opengl as gl
        from meshHelpers import LinkMesh
        # Reset GL cache for this rebuild
        self._gl_items = []
        self._gl_ref_pose = SE3(self.StartDubinsPose.A.copy())
        self._gl_shape_dirty = False
        allElbowHandleSets = []
        if showBoundary:
            vertices = []
            faces = []
            vertex_offset = 0

            if not self.elbow1 is None:
                v, f = self.elbow1.circleEllipseCircleQT(numSides)

                vertices.extend(v)
                faces.extend(f)
                vertex_offset += len(v)
            
            if self.path.tMag > self.DISTANCE_EPSILON:
                v, f = self.cylinder.interpolateQtCircles(numSides, 2)

                vertices.extend(v)
                f2 = [[item + vertex_offset for item in sublist] for sublist in f]
                faces.extend(f2)
                vertex_offset += len(v)

            if not self.elbow2 is None:
                v, f = self.elbow2.circleEllipseCircleQT(numSides)

                vertices.extend(v)
                f2 = [[item + vertex_offset for item in sublist] for sublist in f]
                faces.extend(f2)
                vertex_offset += len(v)

            if (len(vertices) > 0 and len(faces) > 0):
                meshdata = gl.MeshData(vertexes=np.array(vertices), faces=np.array(faces))
                # Apply alpha to color if it's an RGBA tuple
                draw_color = color
                if isinstance(color, (tuple, list)) and len(color) >= 4:
                    draw_color = (color[0], color[1], color[2], alpha)
                link_shader = None if alpha >= 1.0 else 'shaded'
                draw_edges = True if alpha >= 1.0 else wireFrame
                meshitem = LinkMesh(id=linkID, meshdata=meshdata, color=draw_color, drawEdges=draw_edges, shader=link_shader, smooth=True)
                meshitem.setObjectName("Link")
                if alpha >= 1.0:
                    meshitem.setGLOptions('opaque')
                else:
                    meshitem.setGLOptions('translucent')
                widget.plot_widget.addItem(meshitem)
                self._gl_items.append(meshitem)
        elif showFrames:
            # show the start and end frames
            startFrameHandles = addPosesToWidget(np.array([self.StartDubinsPose]), widget, 
                                              axisLength=self.r, xColor='darkred', 
                                              yColor='darkblue', zColor='darkgreen')
            endFrameHandles = addPosesToWidget(np.array([self.EndDubinsPose]), widget, 
                                            axisLength=self.r, xColor='darkred', 
                                            yColor='darkblue', zColor='darkgreen')
            allElbowHandleSets.append(startFrameHandles)
            allElbowHandleSets.append(endFrameHandles)

        if showPath:
            self.path.addToWidget(widget, showCircles=showPathCircles, 
                                showPoses=False, pathColor=pathColor)
        
        return allElbowHandleSets
    
    def updateCachedGLTransforms(self):
        """Apply the rigid-body delta between cached and current pose
        to all cached GL items via setTransform (GPU model matrix, zero-cost).
        Call this instead of removing + re-adding items when only the pose changed."""
        if not self._gl_items or self._gl_ref_pose is None:
            return
        from pyqtgraph import Transform3D
        delta = self.StartDubinsPose @ self._gl_ref_pose.inv()
        mat = Transform3D(delta.A)
        for item in self._gl_items:
            item.setTransform(mat)

    def clearGLCache(self):
        """Discard cached GL items (e.g. before a full rebuild)."""
        self._gl_items = []
        self._gl_ref_pose = None
        self._gl_shape_dirty = True

    def hasGLCache(self) -> bool:
        return not self._gl_shape_dirty

    def transformPosesInPlace(self, Transformation : SE3):
        """Lightweight pose update: only move start/end poses without
        rebuilding path, elbows, or cylinder geometry.  Sufficient for
        updateCachedGLTransforms() which only reads StartDubinsPose."""
        self.StartDubinsPose = Transformation @ self.StartDubinsPose
        self.EndDubinsPose = Transformation @ self.EndDubinsPose

    def endBoundingBalls2r(self):
        startBall = Ball(self.path.circleCenter1, 2*self.r) if self.elbow1 else Ball(self.StartDubinsPose.t, self.r)
        endBall = Ball(self.path.circleCenter2, 2*self.r) if self.elbow2 else Ball(self.EndDubinsPose.t, self.r)
        return startBall, endBall

    def show(self, numSides : int = 32, color : str = linkColorDefault, 
                  alpha : float = 0.5, wireFrame : bool = False, 
                  showFrames : bool = False, showPath : bool = True, 
                  pathColor : str = pathColorDefault,
                  showPathCircles : bool = False, showBoundary : bool = True,
                  showElbowBoundingBalls : bool = False, block : bool = False):
        ax = plt.figure().add_subplot(projection='3d')
        allElbowHandleSets = self.addToPlot(ax, numSides, color, alpha, wireFrame, 
                                     showFrames, showPath, pathColor, showPathCircles,
                                     showBoundary, showElbowBoundingBalls)
        ax.set_aspect('equal')
        plt.show(block=block)    

    def collisionBoxes(self):
        return [CollisionBox(startDubinsFrame=self.StartDubinsPose, endDubinsFrame=self.Elbow1EndFrame, r=self.r),
                CollisionBox(startDubinsFrame=self.Elbow1EndFrame, endDubinsFrame=self.Elbow2StartFrame, r=self.r),
                CollisionBox(startDubinsFrame=self.Elbow2StartFrame, endDubinsFrame=self.EndDubinsPose, r=self.r)]
    
    def getCapsules(self):
        capsules = [CollisionCapsule(base=self.Elbow1EndFrame, radius=self.r, height=np.linalg.norm(self.Elbow2StartFrame.t - self.Elbow1EndFrame.t))]

        lastFrame = [float('inf')] * 3
        if self.elbow1:
            for elbow in self.elbow1.elbows:
                if np.linalg.norm(elbow.StartFrame.t - lastFrame) > self.r/2:
                    capsules.append(CollisionCapsule(base=elbow.StartFrame, radius = self.r, height=np.linalg.norm(elbow.StartFrame.t - elbow.midPoint)))
                    capsules.append(CollisionCapsule(base=elbow.EndFrame, radius = self.r, height=-np.linalg.norm(elbow.EndFrame.t - elbow.midPoint)))
                    lastFrame = elbow.EndFrame.t
        
        lastFrame = [float('inf')] * 3
        if self.elbow2:
            for elbow in self.elbow2.elbows:
                if np.linalg.norm(elbow.StartFrame.t - lastFrame) > self.r/2:
                    capsules.append(CollisionCapsule(base=elbow.StartFrame, radius = self.r, height=np.linalg.norm(elbow.StartFrame.t - elbow.midPoint)))
                    capsules.append(CollisionCapsule(base=elbow.EndFrame, radius = self.r, height=-np.linalg.norm(elbow.EndFrame.t - elbow.midPoint)))
                    lastFrame = elbow.EndFrame.t
        return capsules

    def recomputeCollisionCapsules(self):
        self.collisionCapsules = self.getCapsules()

    def length(self) -> float:
        """Return the total length of the CSC link."""
        return self.lengthC1 + self.lengthS + self.lengthC2
    
    
    def interpolateAt(self, t : float) -> np.ndarray:
        """
        Return the 3D position at parameter t in [0, 1] along the CSC path.
        """
        return self.path.interpolateAt(t)
            
    def interpolate(self, count : Optional[int] = None, density: Optional[float] = None) -> np.ndarray:
        return self.path.interpolate(count=count, density=density)
    
    def interpolate_vectorized(self, t_array: np.ndarray) -> np.ndarray:
        """
        Vectorized interpolation: compute 3D positions for multiple t values at once.
        
        Parameters:
        -----------
        t_array : np.ndarray
            Array of parameter values in [0, 1] (shape (n,))
            
        Returns:
        --------
        np.ndarray
            Array of 3D positions (shape (n, 3))
        """
        return self.path.interpolate_vectorized(t_array)
    
    def sdf(self, point: ArrayLike, radius: Optional[float] = None, xp: ModuleType = np) -> Union[float, ArrayLike]:
        """
        Compute the signed distance from a 3D point to this link.
        
        The SDF is constructed piecewise from the arc and straight segments
        of the CSC path, returning the minimum distance to any segment.
        
        Parameters:
        -----------
        point : np.ndarray
            3D points (shape (3,) or (N, 3)) at which to evaluate the SDF
        radius : float, optional
            Tube radius. If None, uses self.r
            
        Returns:
        --------
        float
            Signed distance (negative inside, positive outside)
        """
        if radius is None:
            radius = self.r
        
        point = xp.asarray(point).reshape(-1, 3)
        N = xp.shape(point)[0]
        distances = xp.full((3,N), xp.inf, dtype=xp.float64)
        
        # 1. First arc (elbow1)
        if self.arc1 is not None and self.path.theta1 > self.EPSILON:
            distances[0,:] = sdf_torus_section_flat_ended(xp, point, 
                                    center=xp.asarray(self.arc1.circleCenter).reshape(1, 3),
                                    R_w2l=xp.asarray(self.arc1.localOrientation()).reshape(1, 3, 3),
                                    sc=xp.asarray(self.arc1.sdfSinCos()).reshape(1, 2),
                                    ra=xp.asarray([self.arc1.r]), # arc radius
                                    rb=xp.asarray([radius])).flatten() # tube radius
            
        # 2. Straight section
        if self.path.tMag > self.DISTANCE_EPSILON:
            # Capsule from turn1end to turn2start
            distances[1,:] = sdf_capsule(xp, point,
                                xp.asarray(self.path.turn1end).reshape(1,3), # (L,3)
                                xp.asarray(self.path.turn2start).reshape(1,3), # (L,3)
                                radius).flatten() # r: float or (L,)
        
        # 3. Second arc (elbow2)
        if self.arc2 is not None and self.path.theta2 > self.EPSILON:
            distances[2,:] = sdf_torus_section_flat_ended(xp, point, 
                                    center=xp.asarray(self.arc2.circleCenter).reshape(1, 3),
                                    R_w2l=xp.asarray(self.arc2.localOrientation()).reshape(1, 3, 3),
                                    sc=xp.asarray(self.arc2.sdfSinCos()).reshape(1, 2),
                                    ra=xp.asarray([self.arc2.r]), # arc radius
                                    rb=xp.asarray([radius])).flatten() # tube radius
        
        # Return minimum distance to any section
        return xp.min(distances, axis=0)

    def boundingBox(self, xp: ModuleType = np, tolerance: float = 0.1) -> ArrayLike:
        """
        Compute axis-aligned bounding box of the link, possibly overestimated tolerance*r.
        
        Parameters:
        -----------
        xp : ModuleType
            Numerical module (e.g., numpy or cupy)
            
        Returns:
        --------
        np.ndarray
            Bounding box as [(min_x, min_y, min_z), 
                             (max_x, max_y, max_z)]
        """
        #points = self.interpolate(density=20.0)

        # The number of points sampled along each arc is chosen to ensure they are spaced by at most tolerance*r
        points = xp.vstack((
            self.arc1.interpolate(int(xp.ceil(self.arc1.theta / tolerance)) + 1, xp=xp) if self.arc1 else xp.asarray(self.StartDubinsPose.t).reshape(1,3),
            self.arc2.interpolate(int(xp.ceil(self.arc2.theta / tolerance)) + 1, xp=xp) if self.arc2 else xp.asarray(self.EndDubinsPose.t).reshape(1,3)
            ))
        
        min_corner = xp.min(points, axis=0) - (1 + tolerance)*self.r
        max_corner = xp.max(points, axis=0) + (1 + tolerance)*self.r
        return xp.vstack((min_corner, max_corner))
    
    def startCircle(self, forward : bool = True) -> Circle3D:
        """
        Return the circular face at the start of the link.
        
        Returns:
        --------
        Circle3D
            The circular disc at the starting position (proximal end of the link),
            oriented with normal pointing along the link direction
        """
        return Circle3D(
            radius=self.r,
            center=self.StartDubinsPose.t,
            normal=self.StartDubinsPose.R[:,0] if forward else -self.StartDubinsPose.R[:,0],
            radialVector=self.StartDubinsPose.R[:,1]
        )
    
    def endCircle(self, forward : bool = True) -> Circle3D:
        """
        Return the circular face at the end of the link.
        
        Returns:
        --------
        Circle3D
            The circular disc at the ending position (distal end of the link),
            oriented with normal pointing along the link direction
        """
        return Circle3D(
            radius=self.r,
            center=self.EndDubinsPose.t,
            normal=self.EndDubinsPose.R[:,0] if forward else -self.EndDubinsPose.R[:,0],
            radialVector=self.EndDubinsPose.R[:,1]
        )

    def splitAtFractions(self, fractions: list) -> list:
        """
        Split this link at arc-length fractions (each in (0,1)).
        Returns len(fractions)+1 sub-links. Returns [self] for empty fractions.
        """
        sub_paths = self.path.splitAtFractions(fractions)
        if len(sub_paths) == 1:
            return [self]
        n = len(sub_paths)
        result = []
        for i, sub_path in enumerate(sub_paths):
            start_frame = (self.StartDubinsPose if i == 0
                           else _dubins_frame_from_pos_and_dir(
                               sub_path.startPosition, sub_path.startDir))
            end_frame   = (self.EndDubinsPose if i == n - 1
                           else _dubins_frame_from_pos_and_dir(
                               sub_path.endPosition, sub_path.endDir))
            result.append(LinkCSC(self.r, start_frame, end_frame,
                                  self.maxAnglePerElbow, path=sub_path,
                                  EPSILON=self.EPSILON))
        return result
