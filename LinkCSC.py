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



class LinkCSC:
    def __init__(self, r : float, StartDubinsPose : SE3, EndDubinsPose : SE3,
                 maxAnglePerElbow : float = np.pi/2, 
                 path : PathCSC = None, EPSILON : float = 0.01):
        assert(r>0)
        assert(maxAnglePerElbow >= 0 and maxAnglePerElbow <= np.pi)
        self.r = r
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
            raise ValueError("ERROR: Tried to generate a link for an invalid path")
        if self.path.theta1 < -EPSILON:
            raise ValueError("ERROR: Tried to generate a link for a path with theta1 < 0")
        if self.path.theta1 >= np.pi:
            raise ValueError("ERROR: Tried to generate a link for a path with theta1 >= pi")
        if self.path.theta2 < -EPSILON:
            raise ValueError("ERROR: Tried to generate a link for a path with theta2 < 0")
        if self.path.theta2 >= np.pi:
            raise ValueError("ERROR: Tried to generate a link for a path with theta2 >= pi")

        
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
        return LinkCSC(self.r, Transformation @ self.StartDubinsPose, 
                       Transformation @ self.EndDubinsPose, 
                       maxAnglePerElbow = self.maxAnglePerElbow, 
                       path = self.path.newPathTransformedBy(Transformation),
                       EPSILON = self.EPSILON)
        
    def addToPlot(self, ax, numSides : int = 32, color : str = linkColorDefault, 
                  alpha : float = 0.5, wireFrame : bool = False, 
                  showFrames : bool = False, showPath : bool = True, 
                  pathColor : str = pathColorDefault,
                  showPathCircles : bool = False, showBoundary : bool = True,
                  showElbowBoundingBalls : bool = False,
                  showManifold : bool = False, startRadius : float = None, 
                  endRadius : float = None, hullBends : bool = False,
                  wallThickness : float = None, 
                  extendBackward : float = 0, extendForward : float = 0):
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
        
        if showManifold:
            mesh = self.manifold(startRadius, endRadius, numSides, hullBends, True, 
                                 wallThickness, extendBackward, extendForward).to_mesh()
            vertices = mesh.vert_properties[:, :3]
            triangles = mesh.tri_verts

            # Create a list of triangle vertex coordinates
            faces = [vertices[tri] for tri in triangles]
            mesh_collection = Poly3DCollection(faces, alpha=0.7, edgecolor='k')
            ax.add_collection3d(mesh_collection)
        
        
        if showPath:
            self.path.addToPlot(ax, showCircles=showPathCircles, 
                                showPoses=showFrames, pathColor=pathColor)
        
        return allElbowHandleSets
    
    def endBoundingBalls2r(self):
        startBall = Ball(self.path.circleCenter1, 2*self.r) if self.elbow1 else Ball(self.StartDubinsPose.t, self.r)
        endBall = Ball(self.path.circleCenter2, 2*self.r) if self.elbow2 else Ball(self.EndDubinsPose.t, self.r)
        return startBall, endBall

    def show(self, numSides : int = 32, color : str = linkColorDefault, 
                  alpha : float = 0.5, wireFrame : bool = False, 
                  showFrames : bool = False, showPath : bool = True, 
                  pathColor : str = pathColorDefault,
                  showPathCircles : bool = False, showBoundary : bool = True,
                  showElbowBoundingBalls : bool = False, block : bool = True,
                  showManifold : bool = False, startRadius : float = None, 
                  endRadius : float = None, hullBends : bool = False, 
                  wallThickness : float = None, 
                  extendBackward : float = 0, extendForward : float = 0):
        ax = plt.figure().add_subplot(projection='3d')
        allElbowHandleSets = self.addToPlot(ax, numSides, color, alpha, wireFrame, 
                                     showFrames, showPath, pathColor, showPathCircles,
                                     showBoundary, showElbowBoundingBalls,
                                     showManifold, startRadius, endRadius, hullBends,
                                     wallThickness, extendBackward, extendForward)
        ax.set_aspect('equal')
        ax.axis('off')
        plt.show(block=block)
    
    def creasePattern(self, numSides : int, twistPortion : float = 0.2) -> TubularPattern:
        assert(numSides >= 4 and numSides%2==0)
        assert(twistPortion > 0)
        
        composed = TubularPattern(numSides, self.r)
        
        if self.path.theta1 > self.EPSILON:
            numElbows1 = (int)(np.ceil(self.path.theta1 / self.maxAnglePerElbow)) 
            elbow1PartPattern = ElbowFittingPattern(numSides, self.r, 
                                    bendingAngle=self.path.theta1 / numElbows1, 
                                    rotationalAxisAngle=self.rot1AxisAngle)
            for i in range(numElbows1):
                composed.append(elbow1PartPattern)
        
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
            numElbows2 = (int)(np.ceil(self.path.theta2 / self.maxAnglePerElbow)) 
            elbow2PartPattern = ElbowFittingPattern(numSides, self.r, 
                                    bendingAngle=self.path.theta2 / numElbows2, 
                                    rotationalAxisAngle=self.rot2AxisAngle)
            for i in range(numElbows2):
                composed.append(elbow2PartPattern)
        
        return composed
    

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
    

    def manifold(self, startRadius : float = None, endRadius : float = None, numSides : int = 20, 
                 hullBends : bool = False, stabilize : bool = True, wallThickness : float = None,
                 extendBackward : float = 0, extendForward : float = 0) -> m3d.Manifold:
        if startRadius is None:
            startRadius = self.r
        if endRadius is None:
            endRadius = self.r
        assert(startRadius > 0 and startRadius <= self.r and endRadius > 0 and endRadius <= self.r)
        assert(numSides >= 3)

        cumulativeLengths = np.cumsum([0, self.r * self.path.theta1 if self.elbow1 else 0,
                            self.path.tMag, self.r * self.path.theta2 if self.elbow2 else 0])
        totalLength = cumulativeLengths[-1]
        ts = cumulativeLengths / totalLength if totalLength > 0 else np.zeros(cumulativeLengths.shape)
        rs = (1-ts) * startRadius + ts * endRadius

        output = m3d.Manifold() # empty manifold
        innerExtendLength = self.DISTANCE_EPSILON if stabilize else 0

        if self.elbow1:
            bend1 = Bend(arcRadius=self.r, StartFrame=self.StartDubinsPose,
                         bendingAngle=self.path.theta1, rotationalAxisAngle=self.rot1AxisAngle,
                         startRadius=rs[0], endRadius=rs[1], 
                         numSides=numSides, maxSectionAngle=self.maxAnglePerElbow)
            C1 = bend1.manifold(hullBends, extendBackward=extendBackward, extendForward=innerExtendLength)
            output = C1
        
        if self.path.tMag > self.DISTANCE_EPSILON:
            distanceBackward = innerExtendLength if self.elbow1 else extendBackward
            distanceForward = innerExtendLength if self.elbow2 else extendForward
            length = self.path.tMag + distanceBackward + distanceForward
            S = m3d.Manifold.cylinder(height=length, radius_low=rs[1], radius_high=rs[2], 
                                      circular_segments=numSides)
            S = S.translate((0,0,-distanceBackward)).rotate((0,90,0)).transform(self.Elbow1EndFrame.A[:3,:])
            output += S
        
        
        if self.elbow2:
            bend2 = Bend(arcRadius=self.r, StartFrame=self.Elbow2StartFrame,
                         bendingAngle=self.path.theta2, rotationalAxisAngle=self.rot2AxisAngle,
                         startRadius=rs[2], endRadius=rs[3], 
                         numSides=numSides, maxSectionAngle=self.maxAnglePerElbow)
            C2 = bend2.manifold(hullBends, extendBackward=innerExtendLength, extendForward=extendForward)
            output += C2

        
        if not wallThickness is None and wallThickness > 0:
            inner = self.manifold(startRadius = startRadius - wallThickness, 
                                  endRadius = endRadius - wallThickness, 
                                  numSides = numSides, hullBends = hullBends, 
                                  stabilize = stabilize, wallThickness = None,
                                  extendBackward = extendBackward + self.DISTANCE_EPSILON,
                                  extendForward = extendForward + self.DISTANCE_EPSILON)
            output -= inner

        return output

    def connectableModule(self, wallThickness : float, holeDiameter : float, numHoles : int = 4,
                           startRadius : float = None, endRadius : float = None) -> m3d.Manifold:
        numSides = 100
        connectionLength = 2 * holeDiameter
        tube = self.manifold(startRadius, endRadius, numSides, stabilize=True,
                             wallThickness=wallThickness)
        
        holeSlicer = m3d.Manifold()
        holeAnglesDegrees = np.linspace(0, 360, numHoles, endpoint=False)
        for angle in holeAnglesDegrees:
            hole = m3d.Manifold.cylinder(height=self.r+self.DISTANCE_EPSILON, 
                                         radius_low=holeDiameter/2, 
                                         radius_high=holeDiameter/2, 
                                         circular_segments=numSides)
            hole = hole.rotate((0,90,0)).rotate((0,0,angle))
            holeSlicer += hole
    
        inset = m3d.Manifold.cylinder(height=2*connectionLength, 
                                       radius_low=startRadius-wallThickness+self.DISTANCE_EPSILON, 
                                       radius_high=startRadius-wallThickness+self.DISTANCE_EPSILON,
                                       circular_segments=numSides)
        tube -= inset.translate((0,0,-connectionLength)).rotate((0,90,0)).transform(self.StartDubinsPose.A[:3,:])
        tube -= holeSlicer.translate((0,0,holeDiameter)).rotate((0,90,0)).transform(self.StartDubinsPose.A[:3,:])

        outset = m3d.Manifold.cylinder(height=2*connectionLength, 
                                       radius_low=endRadius-wallThickness+self.DISTANCE_EPSILON, 
                                       radius_high=endRadius-wallThickness+self.DISTANCE_EPSILON,
                                       circular_segments=numSides)
        outset -= m3d.Manifold.cylinder(height=2*connectionLength, 
                                       radius_low=endRadius-2*wallThickness, 
                                       radius_high=endRadius-2*wallThickness,
                                       circular_segments=numSides)
        outset -= holeSlicer.translate((0,0,3*holeDiameter))
        outset = outset.translate((0,0,-connectionLength)).rotate((0,90,0)).transform(self.EndDubinsPose.A[:3,:])

        tube += outset
        return tube
    

    def saveModule(self, filename : str, wallThickness : float, holeDiameter : float, numHoles : int = 4,
                           startRadius : float = None, endRadius : float = None) -> None:
        module = self.connectableModule(wallThickness, holeDiameter, numHoles, startRadius, endRadius)
        mesh_data = module.to_mesh()
        vertices = mesh_data.vert_properties[:, :3]  # Get XYZ coordinates
        faces = mesh_data.tri_verts
        tri_mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
        tri_mesh.export(filename)


    def interpolateAt(self, t : float) -> np.ndarray:
        """
        Return the 3D position at parameter t in [0, 1] along the CSC path.
        """
        t = max(0.0, min(t, 1.0))

        totalLength = self.lengthC1 + self.lengthS + self.lengthC2

        if totalLength < self.DISTANCE_EPSILON:
            return self.StartDubinsPose.t
        
        s = t * totalLength

        if s <= self.lengthC1:
            if self.arc1:
                t1 = s / self.lengthC1
                return self.arc1.interpolateAt(t1)
            else:
                return self.StartDubinsPose.t
        elif s <= self.lengthC1 + self.lengthS:
            t2 = s - self.lengthC1
            return self.path.turn1end + t2 * self.path.tUnit
        else:
            if self.arc2:
                localS = s - self.lengthC1 - self.lengthS
                t3 = localS / self.lengthC2
                return self.arc2.interpolateAt(t3)
            else:
                return self.EndDubinsPose.t
            
    def interpolate(self, count : int = 10, density: float = None) -> np.ndarray:
        # density is points per unit length
        if density:
            assert density > 0, "Density must be positive"
            totalLength = self.lengthC1 + self.lengthS + self.lengthC2
            count = max(2, int(np.ceil(totalLength * density)) + 1)
        else:
            assert count >= 2, "Count must be at least 2"
        return self.interpolate_vectorized(np.linspace(0, 1, count))
    
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
        # Clamp t values
        t_array = np.clip(t_array, 0.0, 1.0)
        
        totalLength = self.lengthC1 + self.lengthS + self.lengthC2
        if totalLength < self.DISTANCE_EPSILON:
            return np.tile(self.StartDubinsPose.t, (len(t_array), 1))
        
        # Convert t to arc length s
        s_array = t_array * totalLength
        
        # Partition t_array indices by segment
        arc1_mask = s_array <= self.lengthC1
        straight_mask = (s_array > self.lengthC1) & (s_array <= self.lengthC1 + self.lengthS)
        arc2_mask = s_array > self.lengthC1 + self.lengthS
        
        # Initialize output array
        points = np.zeros((len(t_array), 3), dtype=np.float64)
        
        # Arc 1 segment
        if np.any(arc1_mask):
            if self.arc1:
                arc1_indices = np.where(arc1_mask)[0]
                arc1_t = s_array[arc1_indices] / self.lengthC1
                points[arc1_indices] = self.arc1.interpolate_vectorized(arc1_t)
            else:
                arc1_indices = np.where(arc1_mask)[0]
                points[arc1_indices] = np.tile(self.StartDubinsPose.t, (len(arc1_indices), 1))
        
        # Straight segment
        if np.any(straight_mask):
            straight_indices = np.where(straight_mask)[0]
            t2_array = s_array[straight_indices] - self.lengthC1
            # Straight line: start + t2 * direction
            points[straight_indices] = (self.path.turn1end[np.newaxis, :] + 
                                       t2_array[:, np.newaxis] * self.path.tUnit[np.newaxis, :])
        
        # Arc 2 segment
        if np.any(arc2_mask):
            if self.arc2:
                arc2_indices = np.where(arc2_mask)[0]
                localS = s_array[arc2_indices] - self.lengthC1 - self.lengthS
                arc2_t = localS / self.lengthC2
                points[arc2_indices] = self.arc2.interpolate_vectorized(arc2_t)
            else:
                arc2_indices = np.where(arc2_mask)[0]
                points[arc2_indices] = np.tile(self.EndDubinsPose.t, (len(arc2_indices), 1))
        
        return points
    
    def _sdCapsule(self, p: np.ndarray, a: np.ndarray, b: np.ndarray, r: float) -> float:
        """
        Signed Distance Function for a capsule between points a and b.
        
        Parameters:
        -----------
        p : np.ndarray
            3D query point
        a : np.ndarray
            Start point of capsule
        b : np.ndarray
            End point of capsule
        r : float
            Radius of capsule
        
        Returns:
        --------
        float
            Signed distance (negative inside, positive outside)
        """
        pa = p - a
        ba = b - a
        h = np.clip(np.dot(pa, ba) / np.dot(ba, ba), 0.0, 1.0)
        return np.linalg.norm(pa - ba * h) - r
    
    def sdf(self, point: np.ndarray, radius: float = None) -> float:
        """
        Compute the signed distance from a 3D point to this link.
        
        The SDF is constructed piecewise from the arc and straight segments
        of the CSC path, returning the minimum distance to any segment.
        
        Parameters:
        -----------
        point : np.ndarray
            3D point to query (shape (3,))
        radius : float, optional
            Tube radius. If None, uses self.r
            
        Returns:
        --------
        float
            Signed distance (negative inside, positive outside)
        """
        if radius is None:
            radius = self.r
        
        distances = []
        
        # 1. First arc (elbow1)
        if self.arc1 is not None and self.path.theta1 > self.EPSILON:
            dist1 = self.arc1.sdf(point, radius)
            distances.append(dist1)
        else:
            # Just a sphere at the start if no arc
            dist1 = np.linalg.norm(point - self.StartDubinsPose.t) - radius
            distances.append(dist1)
        
        # 2. Straight section
        if self.path.tMag > self.DISTANCE_EPSILON:
            # Capsule from turn1end to turn2start
            dist2 = self._sdCapsule(point, self.path.turn1end, 
                                   self.path.turn1end + self.path.tMag * self.path.tUnit, 
                                   radius)
            distances.append(dist2)
        
        # 3. Second arc (elbow2)
        if self.arc2 is not None and self.path.theta2 > self.EPSILON:
            dist3 = self.arc2.sdf(point, radius)
            distances.append(dist3)
        else:
            # Just a sphere at the end if no arc
            dist3 = np.linalg.norm(point - self.EndDubinsPose.t) - radius
            distances.append(dist3)
        
        # Return minimum distance (union of all segments)
        return min(distances)
