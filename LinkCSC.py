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
                  showFrames : bool = False, showPath : bool = False, 
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
        
        if showPath:
            self.path.addToPlot(ax, showCircles=showPathCircles, 
                                showPoses=showFrames, pathColor=pathColor)
        
        return allElbowHandleSets
    
    def addToWidget(self, widget, numSides : int = 32, color = linkColorDefault, 
                  alpha : float = 0.5, wireFrame : bool = False, 
                  showFrames : bool = False, showPath : bool = False, 
                  pathColor : str = pathColorDefault,
                  showPathCircles : bool = False, showBoundary : bool = True,
                  showElbowBoundingBalls : bool = False, 
                  linkID : int = -1):
        allElbowHandleSets = []
        if showBoundary:

            vertices = []
            faces = []

            if not self.elbow1 is None:
                self.elbow1.addToWidget(widget, numSides, color, 
                  alpha, wireFrame, 
                  showFrames)
            
            if self.path.tMag > self.EPSILON:
                v, f = self.cylinder.addToWidget(widget, numSides, 2, color)
                #print(v)
                #print(f)

            if not self.elbow2 is None:
                self.elbow2.addToWidget(widget, numSides, color, 
                  alpha, wireFrame, showFrames)

                # v, f = self.elbow2.circleEllipseCircleQT(numSides, debug=True)
                
                # v = np.array(v)
                # f = np.array(f)

                # meshdata = gl.MeshData(vertexes=v, faces=f)
                # meshitem = gl.GLMeshItem(meshdata=meshdata, color=color, drawEdges=wireFrame, shader='shaded', smooth=True)
                # meshitem.setObjectName("Link")
                # widget.plot_widget.addItem(meshitem)

            """
            if not self.elbow1 is None:
                v, f = self.elbow1.circleEllipseCircleQT(numSides)

                v = np.array(v)[0]
                f = np.array(f)[0]

                meshdata = gl.MeshData(vertexes=v, faces=f)
                meshitem = gl.GLMeshItem(meshdata=meshdata, color=(1,0,0,1), drawEdges=wireFrame, shader='shaded', smooth=True)
                meshitem.setObjectName("Link")
                widget.plot_widget.addItem(meshitem)
            
            if self.path.tMag > self.EPSILON:
                v, f = self.cylinder.interpolateQtCircles(numSides, 2)

                meshdata = gl.MeshData(vertexes=v, faces=f)
                meshitem = gl.GLMeshItem(meshdata=meshdata, color=(0,0,0,1), drawEdges=wireFrame, shader='shaded', smooth=True)
                meshitem.setObjectName("Link")
                widget.plot_widget.addItem(meshitem)
                
            if not self.elbow2 is None:
                v, f = self.elbow2.circleEllipseCircleQT(numSides)
                
                v = np.array(v)[0]
                f = np.array(f)[0]

                meshdata = gl.MeshData(vertexes=v, faces=f)
                meshitem = gl.GLMeshItem(meshdata=meshdata, color=color, drawEdges=wireFrame, shader='shaded', smooth=True)
                meshitem.setObjectName("Link")
                widget.plot_widget.addItem(meshitem)
            """
                
            # if (len(vertices) > 0 and len(faces) > 0):
            # vertices = np.array(vertices)[0]
            # faces = np.array(faces)[0]
                
            # meshdata = gl.MeshData(vertexes=vertices, faces=faces)
            # meshitem = gl.GLMeshItem(meshdata=meshdata, color=(1,0,0,1), drawEdges=wireFrame, shader='shaded', smooth=True)
            # meshitem.setObjectName("Link")

            # print(vertices)
            # print(faces)
        
        # if showPath:
        #    self.path.addToPlot(ax, showCircles=showPathCircles, 
        #                        showPoses=showFrames, pathColor=pathColor)
        
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
                  showElbowBoundingBalls : bool = False, block : bool = False):
        ax = plt.figure().add_subplot(projection='3d')
        allElbowHandleSets = self.addToPlot(ax, numSides, color, alpha, wireFrame, 
                                     showFrames, showPath, pathColor, showPathCircles,
                                     showBoundary, showElbowBoundingBalls)
        ax.set_aspect('equal')
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

        if self.elbow1:
            for elbow in self.elbow1.elbows:
                capsules.append(CollisionCapsule(base=elbow.StartFrame, radius = self.r, height=np.linalg.norm(elbow.StartFrame.t - elbow.midPoint)))
                capsules.append(CollisionCapsule(base=elbow.EndFrame, radius = self.r, height=-np.linalg.norm(elbow.EndFrame.t - elbow.midPoint)))
        
        if self.elbow2:
            for elbow in self.elbow2.elbows:
                capsules.append(CollisionCapsule(base=elbow.StartFrame, radius = self.r, height=np.linalg.norm(elbow.StartFrame.t - elbow.midPoint)))
                capsules.append(CollisionCapsule(base=elbow.EndFrame, radius = self.r, height=-np.linalg.norm(elbow.EndFrame.t - elbow.midPoint)))
        return capsules

    def recomputeCollisionCapsules(self):
        self.collisionCapsules = self.getCapsules()
