from Joint import *
from TubularPattern import *
from LinkCSC import *
from KinematicChain import *
from Tube import Tube

class OrigamiTube(Tube):
    def __init__(self, numSides : int):
        if not (numSides >= 4 and numSides % 2 == 0):
            raise ValueError("OrigamiTube requires numSides to be an even integer >= 4")
        self.numSides = numSides
        self.polygonInnerAngle = np.pi * (numSides-2)/(2*numSides)
        self._pattern = None
    
    @abstractmethod
    def creasePattern(self) -> TubularPattern:
        pass

class OrigamiExtendedRevolute(OrigamiTube, TransverseRevolute):
    """
    An origami revolute joint with tubular extensions on both sides.
    The revolute joint is defined by a total bending angle and the tube length on each side.
    If outerLength is set to None, it is computed to ensure that the end circles do not intersect
    at the specified bending angle range.

    Attributes:
        numSides (int): Number of sides of the origami tube (must be even and >= 4).
        r (float): Radius of the tube.
        totalBendingAngle (float): Total bending angle of the revolute joint in radians.
        outerLength (Optional[float]): Length of the tubular extension on each side of the revolute joint. 
                            If None, it is computed automatically to avoid intersection of end circles.
        numSinkLayers (int): Number of sink layers in the revolute joint pattern.
        initialState (float): Initial angle of the revolute joint in radians. Defaults to 0. 
    """
    def __init__(self, numSides : int, r : float, totalBendingAngle : float, 
                 Pose : SE3, outerLength: Optional[float] = None, numSinkLayers : int = 1,
                 initialState : float = 0, checkCircleOverlap : bool = True):
        if not (totalBendingAngle > 0 and totalBendingAngle < 2*np.pi):
            raise ValueError("OrigamiExtendedRevolute requires totalBendingAngle to be in (0, 2*pi) radians")
        self.totalBendingAngle = totalBendingAngle
        self.numSinkLayers = numSinkLayers
        OrigamiTube.__init__(self, numSides)
        self.revoluteLength = 2*r*np.sin(self.polygonInnerAngle)*np.tan(totalBendingAngle/4) #2*delta from paper
        TransverseRevolute.__init__(self, r, Pose, minAngle=-totalBendingAngle/2, 
                                    maxAngle=totalBendingAngle/2, 
                                    neutralLength=None if outerLength is None else self.revoluteLength + 2*outerLength, 
                                    initialState=initialState, checkCircleOverlap=checkCircleOverlap)
        self.outerLength = (self.neutralLength - self.revoluteLength) / 2
    
    def creasePattern(self) -> TubularPattern:
        if self._pattern is None:
            revolutePattern = RevoluteJointPattern(self.numSides, self.r, 
                                                self.totalBendingAngle, self.numSinkLayers)
            if self.outerLength > 0:
                self._pattern = TubeFittingPattern(self.numSides, self.r, self.outerLength).append(revolutePattern).append(TubeFittingPattern(self.numSides, self.r, self.outerLength))
            else:
                self._pattern = revolutePattern
        return self._pattern
    
    def RevoluteProximalFrame(self) -> SE3:
        PF = self.ProximalFrame()
        return SE3.Trans(self.outerLength*PF.R[:,0]) @ PF
    
    def RevoluteDistalFrame(self) -> SE3:
        DF = self.DistalFrame()
        return SE3.Trans(-self.outerLength*DF.R[:,0]) @ DF
    
    def proximalExtension(self) -> Cylinder:
        PF = self.ProximalFrame()
        uhat = (PF @ SE3.Rx(np.pi/self.numSides)).R[:,1]
        return Cylinder(self.r, PF.t, PF.R[:,0], 
                        self.outerLength, uhat)
    
    def distalExtension(self) -> Cylinder:
        DF = self.DistalFrame()
        uhat = (DF @ SE3.Rx(np.pi/self.numSides)).R[:,1]
        return Cylinder(self.r, DF.t, -DF.R[:,0], 
                        self.outerLength, uhat)
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=revoluteColorDefault, edgeColor=revoluteEdgeColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=True,
             axisScale=10, showPoses=True):
        plotHandles = super().addToPlot(ax, xColor, yColor, zColor, proximalColor,
                          centerColor, distalColor, sphereColor, showSphere,
                          showAxis, axisScale, showPoses)
        if showSurface:
            scale = self.creasePattern().baseSideLength / 2
            CenterSegment = np.array([self.Pose.t - scale * self.Pose.R[:,2],
                                      self.Pose.t + scale * self.Pose.R[:,2]])
            #https://stackoverflow.com/questions/63207496/how-to-visualize-polyhedrons-defined-by-their-vertices-in-3d-with-matplotlib-or
            
            radialCount = self.numSides + 1
            angle = np.linspace(0, 2*np.pi, radialCount) + np.pi/self.numSides
            u = self.r * np.cos(angle)
            v = self.r * np.sin(angle)
            RevoluteProximalPose = self.RevoluteProximalFrame()
            uhatProximal = RevoluteProximalPose.R[:,1]
            vhatProximal = RevoluteProximalPose.R[:,2]
            RevoluteProximalBase = RevoluteProximalPose.t + u.reshape(-1,1) @ uhatProximal.reshape(1,3) + v.reshape(-1,1) @ vhatProximal.reshape(1,3)
            
            ProximalPoints = np.vstack((RevoluteProximalBase, CenterSegment))
            ProximalHull = ConvexHull(ProximalPoints)
            for s in ProximalHull.simplices:
                tri = Poly3DCollection([ProximalPoints[s]])
                tri.set_color(surfaceColor)
                tri.set_alpha(surfaceOpacity)
                ax.add_collection3d(tri)
            
            RevoluteDistalPose = self.RevoluteDistalFrame()
            uhatDistal = RevoluteDistalPose.R[:,1]
            vhatDistal = RevoluteDistalPose.R[:,2]
            DistalBase = RevoluteDistalPose.t + u.reshape(-1,1) @ uhatDistal.reshape(1,3) + v.reshape(-1,1) @ vhatDistal.reshape(1,3)
            DistalPoints = np.vstack((DistalBase, CenterSegment))
            DistalHull = ConvexHull(DistalPoints)
            for s in DistalHull.simplices:
                tri = Poly3DCollection([DistalPoints[s]])
                tri.set_facecolor(surfaceColor)
                tri.set_edgecolor(edgeColor)
                tri.set_alpha(surfaceOpacity)
                ax.add_collection3d(tri)
            
            self.proximalExtension().addToPlot(ax, color=surfaceColor, 
                                              alpha=surfaceOpacity, 
                                              edgeColor=edgeColor,
                                              numPointsPerCircle=self.numSides)
            self.distalExtension().addToPlot(ax, color=surfaceColor, 
                                              alpha=surfaceOpacity, 
                                              edgeColor=edgeColor,
                                              numPointsPerCircle=self.numSides)
        return plotHandles


class OrigamiRevolute(OrigamiExtendedRevolute):
    def __init__(self, numSides : int, r : float, totalBendingAngle : float, 
                 Pose : SE3, numSinkLayers : int = 1, initialState : float = 0):
        super().__init__(numSides, r, totalBendingAngle, Pose, 0, numSinkLayers, initialState, 
                         checkCircleOverlap=False)


class OrigamiPrismatic(OrigamiTube, Prismatic):
    def __init__(self, numSides : int, r : float, neutralLength : float, 
                 numLayers : int, coneAngle : float, Pose : SE3, 
                 initialState : float = 0):
        OrigamiTube.__init__(self, numSides)
        neutralLayerHeight = neutralLength / numLayers
        flatLayerHalfHeight = neutralLayerHeight / (2*np.sin(coneAngle))
        minLength = numLayers*flatLayerHalfHeight
        maxLength = 2*minLength
        self.numLayers = numLayers
        self.coneAngle = coneAngle
        Prismatic.__init__(self, r, neutralLength, minLength, maxLength, Pose, initialState)
    
    def creasePattern(self) -> TubularPattern:
        if self._pattern is None:
            self._pattern = PrismaticJointPattern(self.numSides, self.r, self.neutralLength, 
                                                 self.numLayers, self.coneAngle)
        return self._pattern
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=prismaticColorDefault, edgeColor=prismaticEdgeColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=True, 
             axisScale=10, showPoses=True):
        plotHandles = super().addToPlot(ax=ax, xColor=xColor, yColor=yColor, zColor=zColor, 
                          proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor, 
                          sphereColor=sphereColor, showSphere=showSphere,
                          surfaceColor=surfaceColor, edgeColor=edgeColor,
                          surfaceOpacity=surfaceOpacity, showSurface=False, showAxis=showAxis,
                          axisScale=axisScale, showPoses=showPoses)
        if showSurface:
            self.boundingCylinder().addToPlot(ax, color=surfaceColor, 
                                              alpha=surfaceOpacity, 
                                              edgeColor=edgeColor,
                                              numPointsPerCircle=self.numSides)            
        return plotHandles
    
class OrigamiTip(OrigamiTube, Tip):
    def __init__(self, numSides : int, r : float, Pose : SE3, length : float, 
                 closesForward : bool = True, pathIndex : int = 2):
        OrigamiTube.__init__(self, numSides)
        Tip.__init__(self, r, Pose, length, closesForward, pathIndex)
    
    def creasePattern(self) -> TubularPattern:
        if self._pattern is None:
            self._pattern = TipPattern(self.numSides, self.r, self.length, self.forward)
        return self._pattern
    
    def pathIndex(self) -> int:
        return self.pidx
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=linkColorDefault, edgeColor=linkColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=True,
             axisScale=10, showPoses=True):
        plotHandles = super().addToPlot(ax, xColor, yColor, zColor, proximalColor,
                          centerColor, distalColor, sphereColor, showSphere,
                          surfaceColor, surfaceOpacity, showSurface, showAxis,
                          axisScale, showPoses)

        if showSurface:
            #https://stackoverflow.com/questions/63207496/how-to-visualize-polyhedrons-defined-by-their-vertices-in-3d-with-matplotlib-or
            radialCount = self.numSides + 1
            angle = np.linspace(0, 2*np.pi, radialCount) + np.pi/self.numSides
            u = self.r * np.cos(angle)
            v = self.r * np.sin(angle)
            scale = self.creasePattern().baseSideLength / 2
            
            match self.pidx:
                case 0: tipSegmentIndex, uhatIndex, vhatIndex = 2,1,2
                case 1: tipSegmentIndex, uhatIndex, vhatIndex = 0,2,0
                case 2: tipSegmentIndex, uhatIndex, vhatIndex = 1,0,1

            
            if self.forward:
                DistalPose = self.DistalFrame()
                TipSegment = np.array([DistalPose.t - scale * DistalPose.R[:,tipSegmentIndex],
                                          DistalPose.t + scale * DistalPose.R[:,tipSegmentIndex]])
                
                ProximalPose = self.ProximalFrame()
                uhatProximal = ProximalPose.R[:,uhatIndex]
                vhatProximal = ProximalPose.R[:,vhatIndex]
                ProximalBase = ProximalPose.t + u.reshape(-1,1) @ uhatProximal.reshape(1,3) + v.reshape(-1,1) @ vhatProximal.reshape(1,3)
                ProximalPoints = np.vstack((ProximalBase, TipSegment))
                ProximalHull = ConvexHull(ProximalPoints)
                for s in ProximalHull.simplices:
                    tri = Poly3DCollection([ProximalPoints[s]])
                    tri.set_facecolor(linkColorDefault)
                    tri.set_edgecolor(None)
                    tri.set_alpha(surfaceOpacity)
                    ax.add_collection3d(tri)
            else:
                ProximalPose = self.ProximalFrame()
                TipSegment = np.array([ProximalPose.t - scale * ProximalPose.R[:,tipSegmentIndex],
                                          ProximalPose.t + scale * ProximalPose.R[:,tipSegmentIndex]])
                
                DistalPose = self.DistalFrame()
                uhatDistal = DistalPose.R[:,uhatIndex]
                vhatDistal = DistalPose.R[:,vhatIndex]
                DistalBase = DistalPose.t + u.reshape(-1,1) @ uhatDistal.reshape(1,3) + v.reshape(-1,1) @ vhatDistal.reshape(1,3)
                DistalPoints = np.vstack((DistalBase, TipSegment))
                DistalHull = ConvexHull(DistalPoints)
                for s in DistalHull.simplices:
                    tri = Poly3DCollection([DistalPoints[s]])
                    tri.set_facecolor(surfaceColor)
                    tri.set_edgecolor(edgeColor)
                    tri.set_alpha(surfaceOpacity)
                    ax.add_collection3d(tri)
            
        return plotHandles

class OrigamiStartTip(OrigamiTip):
    def __init__(self, numSides : int, r : float, Pose : SE3, length : float):
        super().__init__(numSides, r, Pose, length, closesForward=False, pathIndex=2)

class OrigamiEndTip(OrigamiTip):
    def __init__(self, numSides : int, r : float, Pose : SE3, length : float):
        super().__init__(numSides, r, Pose, length, closesForward=True, pathIndex=2)


class OrigamiLinkCSC(OrigamiTube, LinkCSC):
    def __init__(self, numSides : int, r : float, StartDubinsPose : SE3, EndDubinsPose : SE3,
                 maxAnglePerElbow : float = np.pi/2, path : Optional[PathCSC] = None, 
                 EPSILON : float = 0.01):
        OrigamiTube.__init__(self, numSides)
        LinkCSC.__init__(self, r, StartDubinsPose, EndDubinsPose, maxAnglePerElbow, path, EPSILON)
    
    @classmethod
    def from_link_params(cls, r: float, StartDubinsPose: SE3, EndDubinsPose: SE3, 
                         maxAnglePerElbow: float = np.pi/2, path: Optional[PathCSC] = None,
                         EPSILON: float = 0.01, **kwargs):
        """Factory method that infers numSides from context or uses default"""
        # For now, use a default numSides. In practice, this could be inferred from the tree's joints
        numSides = kwargs.get('numSides', 6)
        return cls(numSides, r, StartDubinsPose, EndDubinsPose, maxAnglePerElbow, path, EPSILON)
    
    def creasePattern(self, twistPortion : float = 0.2) -> TubularPattern:
        assert(self.numSides >= 4 and self.numSides%2==0)
        assert(twistPortion > 0)
        
        composed = TubularPattern(self.numSides, self.r)
        
        if self.path.theta1 > self.EPSILON:
            numElbows1 = (int)(np.ceil(self.path.theta1 / self.maxAnglePerElbow)) 
            elbow1PartPattern = ElbowFittingPattern(self.numSides, self.r, 
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
            twistPattern = TwistFittingPattern(self.numSides, self.r, twistAngle, twistLen)
            composed.append(twistPattern)
        
        tubeLen = self.path.tMag - twistLen
        tubePattern = TubeFittingPattern(self.numSides, self.r, tubeLen)
        composed.append(tubePattern)
        
        if self.path.theta2 > self.EPSILON:
            numElbows2 = (int)(np.ceil(self.path.theta2 / self.maxAnglePerElbow)) 
            elbow2PartPattern = ElbowFittingPattern(self.numSides, self.r, 
                                    bendingAngle=self.path.theta2 / numElbows2, 
                                    rotationalAxisAngle=self.rot2AxisAngle)
            for i in range(numElbows2):
                composed.append(elbow2PartPattern)
        
        return composed
    
    def newLinkTransformedBy(self, Transformation : SE3):
        """Override to preserve OrigamiLinkCSC type when transforming"""
        return OrigamiLinkCSC(self.numSides, self.r, 
                             Transformation @ self.StartDubinsPose, 
                             Transformation @ self.EndDubinsPose, 
                             maxAnglePerElbow = self.maxAnglePerElbow, 
                             path = self.path.newPathTransformedBy(Transformation),
                             EPSILON = self.EPSILON)


class OrigamiKinematicChain(KinematicChain[OrigamiTube]):
    """KinematicChain constrained to OrigamiTube fabrication with a fixed numSides"""
    _fabrication_type = OrigamiTube  # Class-level fabrication type constraint
    
    def __init__(self, startJoint : Joint, numSides : int, maxAnglePerElbow : float = np.pi/2, gimbal : bool = False,
                 joints : Optional[list[Joint]] = None, links : Optional[list[LinkCSC]] = None, 
                 parents : Optional[list[int]] = None, children : Optional[list[list[int]]] = None, 
                 boundingBall : Optional[Ball] = None, units : str = "Centimeter (cm)"):
        # Validate numSides
        if not (numSides >= 4 and numSides % 2 == 0):
            raise ValueError("OrigamiKinematicChain requires numSides to be an even integer >= 4")
        
        # Check that startJoint has correct numSides if it's an OrigamiTube
        if isinstance(startJoint, OrigamiTube):
            if startJoint.numSides != numSides:
                raise ValueError(f"startJoint has numSides={startJoint.numSides} but chain requires numSides={numSides}")
        elif not isinstance(startJoint, Waypoint):
            raise ValueError("startJoint must be an OrigamiTube or Waypoint")
        
        self.numSides = numSides
        super().__init__(startJoint=startJoint, maxAnglePerElbow=maxAnglePerElbow, gimbal=gimbal,
                         joints=joints, links=links, parents=parents, children=children, 
                         boundingBall=boundingBall, units=units)
    
    def append(self, newJoint : Joint, relative : bool = True, 
                 fixedPosition : bool = False, fixedOrientation : bool = False, 
                 safe : bool = True, chooseXhatToMinPath : bool = False, 
                 relativeToDistalDubins : bool = False) -> int:
        """Override append to enforce numSides constraint on OrigamiTube joints"""
        # Check that newJoint has correct numSides if it's an OrigamiTube
        if isinstance(newJoint, OrigamiTube):
            if newJoint.numSides != self.numSides:
                raise ValueError(f"Cannot append joint with numSides={newJoint.numSides} to chain with numSides={self.numSides}")
        elif not isinstance(newJoint, Waypoint):
            raise ValueError("newJoint must be an OrigamiTube or Waypoint")
        
        return super().append(newJoint, relative, fixedPosition, fixedOrientation, safe, 
                              chooseXhatToMinPath, relativeToDistalDubins=relativeToDistalDubins)
    
    def _get_link_constructor(self):
        """Return a callable that creates OrigamiLinkCSC with the chain's numSides"""
        def make_origami_link(r, start_pose, end_pose, max_angle_per_elbow, path=None, epsilon=0.01):
            return OrigamiLinkCSC(self.numSides, r, start_pose, end_pose, max_angle_per_elbow, path, epsilon)
        return make_origami_link
    
    
    def creasePattern(self, twistPortion : float = 0.2) -> TubularPattern:
        chainPattern = copy.deepcopy(self.Joints[0].creasePattern()) if not isinstance(self.Joints[0], Waypoint) else TubularPattern(self.numSides, self.Joints[0].r)
        for j in range(1, len(self.Joints)):
            chainPattern.append(self.Links[j].creasePattern(twistPortion))
            if not isinstance(self.Joints[j], Waypoint):
                chainPattern.append(self.Joints[j].creasePattern())
        return chainPattern

def loadOrigamiChain(filename : str):
    """Load an OrigamiKinematicChain from a file with origami-specific classes"""
    from numpy import array
    from geometryHelpers import Ball, Cylinder, Plane, Circle3D, Arc3D
    from KinematicTree import KinematicTree
    from LinkCSC import LinkCSC
    from PathCSC import PathCSC
    from spatialmath import SE3
    from Joint import Prismatic, TransverseRevolute, CoaxialRevolute, Waypoint, Tip
    
    try:
        with open(filename) as f:
            data = f.read()
            # Provide necessary imports for eval, including origami-specific classes
            eval_namespace = {
                'KinematicTree': KinematicTree,
                'LinkCSC': LinkCSC,
                'OrigamiLinkCSC': OrigamiLinkCSC,
                'PathCSC': PathCSC,
                'SE3': SE3,
                'array': array,
                'Ball': Ball,
                'Cylinder': Cylinder,
                'Plane': Plane,
                'Circle3D': Circle3D,
                'Arc3D': Arc3D,
                'Prismatic': Prismatic,
                'TransverseRevolute': TransverseRevolute,
                'CoaxialRevolute': CoaxialRevolute,
                'Waypoint': Waypoint,
                'Tip': Tip,
                'OrigamiRevolute': OrigamiRevolute,
                'OrigamiPrismatic': OrigamiPrismatic,
                'OrigamiTip': OrigamiTip,
                'StartTip': StartTip,
                'EndTip': EndTip,
                'OrigamiKinematicChain': OrigamiKinematicChain,
            }
            chain = eval(data, eval_namespace)
            f.close()
            return chain
    except Exception as e:
        print(e)
        raise Exception(f"Could not load file {filename}: {e}")
