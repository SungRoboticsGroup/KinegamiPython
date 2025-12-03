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
                 extendBackward : float = 0, extendForward : float = 0, maxSectionAngle : float = None) -> m3d.Manifold:
        if maxSectionAngle is None:
            maxSectionAngle = self.maxAnglePerElbow
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

        output = m3d.Manifold() # empty manifold

        print("computing manifold")
        if self.elbow1:
            bend1 = Bend(arcRadius=self.r, StartFrame=self.StartDubinsPose,
                         bendingAngle=self.path.theta1, rotationalAxisAngle=self.rot1AxisAngle,
                         startRadius=rs[0], endRadius=rs[1], 
                         numSides=numSides, maxSectionAngle=maxSectionAngle)
            C1 = bend1.manifold(hullBends, extendBackward=extendBackward, extendForward=innerExtendLength)
            output += C1
            startS = bend1.circles[-1]
        else:
            startS = Circle3D(radius=rs[1], center=self.StartDubinsPose.t, normal=self.StartDubinsPose.R[:,0], 
                              radialVector=self.StartDubinsPose.R[:,1]).interpolate(numSides+1)[:-1]
        if self.elbow2:
            bend2 = Bend(arcRadius=self.r, StartFrame=self.Elbow2StartFrame,
                         bendingAngle=self.path.theta2, rotationalAxisAngle=self.rot2AxisAngle,
                         startRadius=rs[2], endRadius=rs[3], 
                         numSides=numSides, maxSectionAngle=maxSectionAngle)
            C2 = bend2.manifold(hullBends, extendBackward=innerExtendLength, extendForward=extendForward)
            output += C2
            endS = bend2.circles[0]
        else:
            endS = Circle3D(radius=rs[2], center=self.EndDubinsPose.t, normal=self.EndDubinsPose.R[:,0], 
                           radialVector=self.EndDubinsPose.R[:,1]).interpolate(numSides+1)[:-1]
            
        if self.path.tMag > self.DISTANCE_EPSILON:
            S = m3d.Manifold.hull_points(np.vstack((startS, endS)))
            output += S

            """
            distanceBackward = innerExtendLength if self.elbow1 else extendBackward
            distanceForward = innerExtendLength if self.elbow2 else extendForward
            length = self.path.tMag + distanceBackward + distanceForward
            S = m3d.Manifold.cylinder(height=length, radius_low=rs[1], radius_high=rs[2], 
                                      circular_segments=numSides)
            S = S.translate((0,0,-distanceBackward)).rotate((0,90,0)).transform(self.Elbow1EndFrame.A[:3,:])
            output += S
            """

        
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
                           startRadius : float = None, endRadius : float = None, numSides : int = 50,
                           hullBends : bool = False, trussify : bool = False, trussCenterline : bool = False, 
                           trussNumSides : int = 6, trussMaxSectionAngle : float = np.pi/2) -> m3d.Manifold:
        connectionLength = 2 * holeDiameter
        if trussify:
            insetStartRadius = startRadius-wallThickness/2
            insetEndRadius = endRadius-wallThickness/2
            insetSolid = self.manifold(insetStartRadius, insetEndRadius,
                                  trussNumSides, stabilize=True, wallThickness=None, 
                                  hullBends=hullBends, maxSectionAngle=trussMaxSectionAngle)
            insetSolid = insetSolid.refine_to_length(self.r)
            if trussCenterline:
                """insetSolid -= self.manifold(0.4*insetStartRadius, 0.4*insetEndRadius,
                                  trussNumSides, stabilize=True, wallThickness=None, 
                                  hullBends=hullBends, maxSectionAngle=trussMaxSectionAngle)"""
                        
                numCenterlinePoints = max(3, int(self.path.length/self.r))
                centerlineVertices = self.path.interpolate(numPoints=numCenterlinePoints)
                centerlineEdges = np.hstack((np.arange(numCenterlinePoints-1).reshape(-1,1), 
                                             np.arange(1, numCenterlinePoints).reshape(-1,1)))
                insetSolidVertices, insetSolidEdges = manifoldToGraph(insetSolid)
                cV, cE = connectOuterToInner(insetSolidVertices, insetSolidEdges, 
                                          centerlineVertices, centerlineEdges,
                                          nearestCount=2)
                tube = trussManifold(cV, cE, diameter=wallThickness)
            else:
                tube = manifoldToTruss(insetSolid, wallThickness)
            baseHeight = 1.5*connectionLength
            baseTopRadius = (startRadius**2 - baseHeight**2)**0.5
            base = m3d.Manifold.cylinder(height=baseHeight, radius_low=startRadius, 
                                         radius_high=baseTopRadius, 
                                         circular_segments=numSides)
            if not trussCenterline:
                base -= m3d.Manifold.cylinder(height=baseHeight+self.DISTANCE_EPSILON, 
                                          radius_low=startRadius-wallThickness,
                                          radius_high=baseTopRadius-wallThickness, 
                                          circular_segments=numSides)
            tube += base.rotate((0,90,0)).transform(self.StartDubinsPose.A[:3,:])
            trimStart = m3d.Manifold.cylinder(height=connectionLength, radius_low=startRadius+wallThickness/2, 
                                         radius_high=startRadius+wallThickness/2, circular_segments=numSides)
            tube -= trimStart.translate((0,0,-connectionLength)).rotate((0,90,0)).transform(self.StartDubinsPose.A[:3,:])
            trimEnd = m3d.Manifold.cylinder(height=connectionLength, radius_low=endRadius+wallThickness/2, 
                                         radius_high=endRadius+wallThickness/2, circular_segments=numSides)
            tube -= trimEnd.translate((0,0,0)).rotate((0,90,0)).transform(self.EndDubinsPose.A[:3,:])

        else:
            if trussCenterline:
                raise ValueError("trussInfill==True only works with truss==True")
            tube = self.manifold(startRadius, endRadius, numSides, stabilize=True,
                                wallThickness=wallThickness, hullBends=hullBends)
        
        holeSlicer = m3d.Manifold()
        holeAnglesDegrees = np.linspace(0, 360, numHoles, endpoint=False)
        for angle in holeAnglesDegrees:
            hole = m3d.Manifold.cylinder(height=self.r+self.DISTANCE_EPSILON, 
                                         radius_low=holeDiameter/2, 
                                         radius_high=holeDiameter/2, 
                                         circular_segments=numSides)
            hole = hole.rotate((0,90,0)).rotate((0,0,angle))
            holeSlicer += hole
    
        inset = m3d.Manifold.cylinder(height=2*connectionLength+2*self.DISTANCE_EPSILON, 
                                       radius_low=startRadius-wallThickness+self.DISTANCE_EPSILON, 
                                       radius_high=startRadius-wallThickness+self.DISTANCE_EPSILON,
                                       circular_segments=numSides)
        tube -= inset.translate((0,0,-connectionLength-self.DISTANCE_EPSILON)).rotate((0,90,0)).transform(self.StartDubinsPose.A[:3,:])
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