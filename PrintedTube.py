from Joint import *
from TubularPattern import *
from TubularPattern import revoluteColorDefault, revoluteEdgeColorDefault, sphereColorDefault, surfaceOpacityDefault, xColorDefault, yColorDefault, zColorDefault
from geometryHelpers import revoluteColorDefault, revoluteEdgeColorDefault, sphereColorDefault, surfaceOpacityDefault, xColorDefault, yColorDefault, zColorDefault
from LinkCSC import *
from Tube import Tube
from KinematicTree import KinematicTree
from mpl_toolkits.mplot3d import Axes3D

class PrintedTube(Tube):
    def __init__(self, wallThickness : float, holeDiameter : float, numHoles : int):
        self.wallThickness = wallThickness
        self.holeDiameter = holeDiameter
        self.numHoles = numHoles
    
class TransverseRDS3225(PrintedTube, TransverseRevolute):
    """
    RDS3225 Servo Motor with brackets and 3D-printed parts attached to make it attach transversely to tubes.
    Works for 180 degree or 270 degree versions. Dimensions are from our CAD model of the 3D printed joint, in mm.
    In the future this parameterized model could be written in the manifold library so wallThickness, holeDiameter, 
    and numHoles could be adjusted as needed. Radius r is as small as possible to bound the servo.

    Angles are in radians where zero is facing forward (x axis of the Pose, i.e., aligning the incoming and outgoing tubes).
    The states are [-np.pi/2, np.pi/2] for 180 degree version, or [-3*np.pi/4, 3*np.pi/4] for 270 degree version.
    When using this in a robot, attach the servo horn and brackets to face forward at 90 or 135 degrees respectively, 
    and then subtract 90 or 135 degrees from the joint state when commanding the servo.

    # TODO: edit hole placement in CAD model, update dimensions here
    """
    # Dimensions from CAD model (mm) — update these when CAD changes
    R = 33.0              # outer radius of tube in mm
    WALL_THICKNESS = 3.0  # wall thickness in mm
    HOLE_DIAMETER = 3.0   # We use M3 bolts
    NUM_HOLES = 4
    NEUTRAL_LENGTH = 76.791  # length in mm from opposite ends of CAD model, excluding the protrusion to inset into the next tube

    def __init__(self, Pose : SE3, version : int | float | str, initialState : float = 0.0):
        if version == 180 or version == "180" or version == 180.0 or version == np.pi:
            maxBendingAngle = np.pi
        elif version == 270 or version == "270" or version == 270.0 or version == 3*np.pi/2:
            maxBendingAngle = 3*np.pi/2
        else:
            raise ValueError("version must be 180 or 270")
        
        PrintedTube.__init__(self, wallThickness=self.WALL_THICKNESS, holeDiameter=self.HOLE_DIAMETER, numHoles=self.NUM_HOLES)
        TransverseRevolute.__init__(self, r=self.R, Pose=Pose, 
                                    minAngle=-maxBendingAngle/2, maxAngle=maxBendingAngle/2, 
                                    neutralLength=self.NEUTRAL_LENGTH, initialState=initialState,
                                    checkCircleOverlap=False)
    
class CoaxialRDS3225(PrintedTube, CoaxialRevolute):
    """
    RDS3225 Servo Motor with brackets and 3D-printed parts attached to make it attach coaxially to tubes.
    Works for 180 degree or 270 degree versions. Dimensions are from our CAD model of the 3D printed joint, in mm.
    In the future this parameterized model could be written in the manifold library so wallThickness, holeDiameter, 
    and numHoles could be adjusted as needed. Radius r is as small as possible to bound the servo.

    Angles are in radians where zero is facing forward (x axis of the Pose, i.e., aligning the incoming and outgoing tubes).
    The states are [-np.pi/2, np.pi/2] for 180 degree version, or [-3*np.pi/4, 3*np.pi/4] for 270 degree version.
    When using this in a robot, attach the servo horn and brackets to face forward at 90 or 135 degrees respectively, 
    and then subtract 90 or 135 degrees from the joint state when commanding the servo.

    # TODO: edit hole placement in CAD model, update dimensions here
    """
    # Dimensions from CAD model (mm) — update these when CAD changes
    R = 33.0              # outer radius of tube in mm
    WALL_THICKNESS = 3.0  # wall thickness in mm
    HOLE_DIAMETER = 3.0   # We use M3 bolts
    NUM_HOLES = 4
    NEUTRAL_LENGTH = 62.4  # length in mm

    def __init__(self, Pose : SE3, version : [int, float, str], initialState : float = 0.0):
        if version == 180 or version == "180" or version == 180.0 or version == np.pi:
            maxBendingAngle = np.pi
        elif version == 270 or version == "270" or version == 270.0 or version == 3*np.pi/2:
            maxBendingAngle = 3*np.pi/2
        else:
            raise ValueError("version must be 180 or 270")
        
        PrintedTube.__init__(self, wallThickness=self.WALL_THICKNESS, holeDiameter=self.HOLE_DIAMETER, numHoles=self.NUM_HOLES)
        CoaxialRevolute.__init__(self, r=self.R, Pose=Pose, neutralLength=self.NEUTRAL_LENGTH,
                                 minAngle=-maxBendingAngle/2, maxAngle=maxBendingAngle/2, 
                                 initialState=initialState)
    
class PrintedHemisphere(PrintedTube, Tip):
    def __init__(self, r : float, Pose : SE3, closesForward : bool, pathIndex : int = 2):
        PrintedTube.__init__(self, wallThickness=3.0, holeDiameter=3.0, numHoles=4)
        Tip.__init__(self, r, Pose, length=r, closesForward=closesForward, pathIndex=pathIndex)

class PrintedStartHemisphere(PrintedHemisphere):
    def __init__(self, r : float, Pose : SE3, pathIndex : int = 2):
        super().__init__(r, Pose, closesForward=False, pathIndex=pathIndex)

class PrintedEndHemisphere(PrintedHemisphere):
    def __init__(self, r : float, Pose : SE3, pathIndex : int = 2):
        super().__init__(r, Pose, closesForward=True, pathIndex=pathIndex)


class PrintedLinkCSC(PrintedTube, LinkCSC):
    def __init__(self, r : float, StartDubinsPose : SE3, EndDubinsPose : SE3,
                 wallThickness : float, holeDiameter : float, numHoles : int,
                 maxAnglePerElbow : float = np.pi/10, path : Optional[PathCSC] = None, 
                 EPSILON : float = 0.01, startRadius : Optional[float] = None, 
                 endRadius : Optional[float] = None):
        PrintedTube.__init__(self, wallThickness, holeDiameter, numHoles)
        LinkCSC.__init__(self, r, StartDubinsPose, EndDubinsPose,
                            maxAnglePerElbow=maxAnglePerElbow, path=path, EPSILON=EPSILON)
        self.startRadius = startRadius if startRadius is not None else r
        self.endRadius = endRadius if endRadius is not None else r
        
    def newLinkTransformedBy(self, Transformation : SE3):
        """Override to preserve PrintedLinkCSC type when transforming"""
        return PrintedLinkCSC(self.r, Transformation @ self.StartDubinsPose, 
                              Transformation @ self.EndDubinsPose,
                              wallThickness = self.wallThickness,
                              holeDiameter = self.holeDiameter,
                              numHoles = self.numHoles,
                              maxAnglePerElbow = self.maxAnglePerElbow, 
                              path = self.path.newPathTransformedBy(Transformation),
                              EPSILON = self.EPSILON,
                              startRadius = self.startRadius,
                              endRadius = self.endRadius)

    def manifold(self, startRadius : Optional[float] = None, endRadius : Optional[float] = None, 
                 numSides : int = 20,  hullBends : bool = False, stabilize : bool = True, 
                 wallThickness : Optional[float] = None, extendBackward : float = 0, 
                 extendForward : float = 0, maxSectionAngle : Optional[float] = None) -> m3d.Manifold:
        if maxSectionAngle is None:
            maxSectionAngle = self.maxAnglePerElbow
        if startRadius is None:
            startRadius = self.startRadius
        if endRadius is None:
            endRadius = self.endRadius
        assert(startRadius > 0 and startRadius <= self.r and endRadius > 0 and endRadius <= self.r)
        assert(numSides >= 3)

        cumulativeLengths = np.cumsum([0, self.r * self.path.theta1 if self.elbow1 else 0,
                            self.path.tMag, self.r * self.path.theta2 if self.elbow2 else 0])
        totalLength = cumulativeLengths[-1]
        ts = cumulativeLengths / totalLength if totalLength > 0 else np.zeros(cumulativeLengths.shape)
        rs = (1-ts) * startRadius + ts * endRadius
        innerExtendLength = self.DISTANCE_EPSILON if stabilize else 0
        output = m3d.Manifold() # empty manifold

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
        
        if not wallThickness is None and wallThickness > 0:
            inner = self.manifold(startRadius = startRadius - wallThickness, 
                                  endRadius = endRadius - wallThickness, 
                                  numSides = numSides, hullBends = hullBends, 
                                  stabilize = stabilize, wallThickness = None,
                                  extendBackward = extendBackward + self.DISTANCE_EPSILON,
                                  extendForward = extendForward + self.DISTANCE_EPSILON)
            output -= inner

        return output

    def connectableModule(self, numSides : int = 50,
                           hullBends : bool = False, trussify : bool = False, trussCenterline : bool = False, 
                           trussNumSides : int = 6, trussMaxSectionAngle : float = np.pi/2) -> m3d.Manifold:
        connectionLength = 2 * self.holeDiameter
        if trussify:
            insetStartRadius = self.startRadius-self.wallThickness/2
            insetEndRadius = self.endRadius-self.wallThickness/2
            insetSolid = self.manifold(insetStartRadius, insetEndRadius,
                                  trussNumSides, stabilize=True, wallThickness=None, 
                                  hullBends=hullBends, maxSectionAngle=trussMaxSectionAngle)
            insetSolid = insetSolid.refine_to_length(self.r)
            if trussCenterline:
                numCenterlinePoints = max(3, int(self.path.length/self.r))
                centerlineVertices = self.path.interpolate(count=numCenterlinePoints)
                centerlineEdges = np.hstack((np.arange(numCenterlinePoints-1).reshape(-1,1), 
                                             np.arange(1, numCenterlinePoints).reshape(-1,1)))
                insetSolidVertices, insetSolidEdges = manifoldToGraph(insetSolid)
                cV, cE = connectOuterToInner(insetSolidVertices, insetSolidEdges, 
                                          centerlineVertices, centerlineEdges,
                                          nearestCount=2)
                tube = trussManifold(cV, cE, diameter=self.wallThickness)
            else:
                tube = manifoldToTruss(insetSolid, self.wallThickness)
            baseHeight = 1.5*connectionLength
            baseTopRadius = (self.startRadius**2 - baseHeight**2)**0.5
            base = m3d.Manifold.cylinder(height=baseHeight, radius_low=self.startRadius, 
                                         radius_high=baseTopRadius, 
                                         circular_segments=numSides)
            if not trussCenterline:
                base -= m3d.Manifold.cylinder(height=baseHeight+self.DISTANCE_EPSILON, 
                                          radius_low=self.startRadius-self.wallThickness,
                                          radius_high=baseTopRadius-self.wallThickness, 
                                          circular_segments=numSides)
            tube += base.rotate((0,90,0)).transform(self.StartDubinsPose.A[:3,:])
            trimStart = m3d.Manifold.cylinder(height=connectionLength, 
                                              radius_low=self.startRadius+self.wallThickness/2, 
                                              radius_high=self.startRadius+self.wallThickness/2, 
                                              circular_segments=numSides)
            tube -= trimStart.translate((0,0,-connectionLength)).rotate((0,90,0)).transform(self.StartDubinsPose.A[:3,:])
            trimEnd = m3d.Manifold.cylinder(height=connectionLength, radius_low=self.endRadius+self.wallThickness/2, 
                                         radius_high=self.endRadius+self.wallThickness/2, circular_segments=numSides)
            tube -= trimEnd.translate((0,0,0)).rotate((0,90,0)).transform(self.EndDubinsPose.A[:3,:])

        else:
            if trussCenterline:
                raise ValueError("trussInfill==True only works with truss==True")
            tube = self.manifold(self.startRadius, self.endRadius, numSides, stabilize=True,
                                wallThickness=self.wallThickness, hullBends=hullBends)
        
        holeSlicer = m3d.Manifold()
        holeAnglesDegrees = np.linspace(0, 360, self.numHoles, endpoint=False)
        for angle in holeAnglesDegrees:
            hole = m3d.Manifold.cylinder(height=self.r+self.DISTANCE_EPSILON, 
                                         radius_low=self.holeDiameter/2, 
                                         radius_high=self.holeDiameter/2, 
                                         circular_segments=numSides)
            hole = hole.rotate((0,90,0)).rotate((0,0,angle))
            holeSlicer += hole
    
        inset = m3d.Manifold.cylinder(height=2*connectionLength+2*self.DISTANCE_EPSILON, 
                                       radius_low=self.startRadius-self.wallThickness+self.DISTANCE_EPSILON, 
                                       radius_high=self.startRadius-self.wallThickness+self.DISTANCE_EPSILON,
                                       circular_segments=numSides)
        tube -= inset.translate((0,0,-connectionLength-self.DISTANCE_EPSILON)).rotate((0,90,0)).transform(self.StartDubinsPose.A[:3,:])
        tube -= holeSlicer.translate((0,0,self.holeDiameter)).rotate((0,90,0)).transform(self.StartDubinsPose.A[:3,:])
        outset = m3d.Manifold.cylinder(height=2*connectionLength, 
                                       radius_low=self.endRadius-self.wallThickness+self.DISTANCE_EPSILON, 
                                       radius_high=self.endRadius-self.wallThickness+self.DISTANCE_EPSILON,
                                       circular_segments=numSides)
        outset -= m3d.Manifold.cylinder(height=2*connectionLength, 
                                       radius_low=self.endRadius-2*self.wallThickness, 
                                       radius_high=self.endRadius-2*self.wallThickness,
                                       circular_segments=numSides)
        outset -= holeSlicer.translate((0,0,3*self.holeDiameter))
        outset = outset.translate((0,0,-connectionLength)).rotate((0,90,0)).transform(self.EndDubinsPose.A[:3,:])

        tube += outset
        return tube
    

    def saveModule(self, filename : str, startRadius : float = None, endRadius : float = None) -> None:
        module = self.connectableModule()
        mesh_data = module.to_mesh()
        vertices = mesh_data.vert_properties[:, :3]  # Get XYZ coordinates
        faces = mesh_data.tri_verts
        tri_mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
        tri_mesh.export(filename)

    def addToPlot(self, ax: Axes3D, numSides : int = 32, color : str = linkColorDefault, 
                  alpha : float = 0.5, wireFrame : bool = False, 
                  showFrames : bool = False, showPath : bool = True, 
                  pathColor : str = pathColorDefault,
                  showPathCircles : bool = False, showBoundary : bool = True,
                  showElbowBoundingBalls : bool = False, showModule : bool = False):
        """
        Add this link to a 3D plot, with optional module visualization.
        
        Parameters
        ----------
        showModule : bool, default=False
            If True, displays the connectable module instead of the usual surface.
            When True, showBoundary is automatically set to False.
        
        Other parameters are inherited from LinkCSC.addToPlot()
        """
        if showModule:
            # Generate and display the connectable module
            module = self.connectableModule()
            mesh = module.to_mesh()
            vertices = mesh.vert_properties[:, :3]
            triangles = mesh.tri_verts
            
            # Create a list of triangle vertex coordinates
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
            faces = [vertices[tri] for tri in triangles]
            mesh_collection = Poly3DCollection(faces, alpha=alpha, 
                                              edgecolor='k' if wireFrame else None,
                                              facecolors=color)
            ax.add_collection3d(mesh_collection)
            
        return super().addToPlot(ax, numSides=numSides, color=color, alpha=alpha,
                                   wireFrame=wireFrame, showFrames=showFrames,
                                   showPath=showPath, pathColor=pathColor,
                                   showPathCircles=showPathCircles, 
                                   showBoundary=showBoundary and not showModule,
                                   showElbowBoundingBalls=showElbowBoundingBalls)

    def show(self, numSides : int = 32, color : str = linkColorDefault, 
             alpha : float = 0.5, wireFrame : bool = False, 
             showFrames : bool = False, showPath : bool = True, 
             pathColor : str = pathColorDefault,
             showPathCircles : bool = False, showBoundary : bool = True,
             showElbowBoundingBalls : bool = False, block : bool = False,
             showModule : bool = False):
        """
        Display this link in a 3D plot window.
        
        Parameters
        ----------
        showModule : bool, default=False
            If True, displays the connectable module instead of the usual surface.
            When True, showBoundary is automatically set to False.
        block : bool, default=False
            If True, blocks execution until the plot window is closed.
        
        Other parameters are inherited from LinkCSC.show()
        """
        ax: Axes3D = plt.figure().add_subplot(projection='3d')
        allElbowHandleSets = self.addToPlot(ax, numSides, color, alpha, wireFrame, 
                                     showFrames, showPath, pathColor, showPathCircles,
                                     showBoundary, showElbowBoundingBalls, showModule)
        ax.set_aspect('equal')
        plt.show(block=block)

def branchingModule(links : list[PrintedLinkCSC], numSides : int = 50, hullBends : bool = False,
                     maxSectionAngle : float = np.pi/10) -> m3d.Manifold:
    """Generate a connectable branching module from multiple PrintedLinkCSC links sharing the same start pose"""
    if len(links) < 1:
        raise ValueError("At least one link is required to create a branching module.")
    # Verify that all links share the same start pose
    startPose = links[0].StartDubinsPose
    DISTANCE_EPSILON = links[0].DISTANCE_EPSILON
    wallThickness = links[0].wallThickness
    holeDiameter = links[0].holeDiameter
    numHoles = links[0].numHoles
    r = links[0].r
    startRadius = links[0].startRadius
    
    outer = m3d.Manifold()
    inner = m3d.Manifold()
    for link in links:
        if not np.all(np.isclose(link.StartDubinsPose.A, startPose.A)):
            raise ValueError("All links must share the same start pose to create a branching module.")
        if not abs(r - link.r) < DISTANCE_EPSILON:
            raise ValueError("All links must have the same radius to create a branching module.")
        if not abs(wallThickness - link.wallThickness) < DISTANCE_EPSILON:
            raise ValueError("All links must have the same wall thickness to create a branching module.")
        if not abs(holeDiameter - link.holeDiameter) < DISTANCE_EPSILON:
            raise ValueError("All links must have the same hole diameter to create a branching module.")
        if not numHoles == link.numHoles:
            raise ValueError("All links must have the same number of holes to create a branching module.")
        if not abs(startRadius - link.startRadius) < DISTANCE_EPSILON:
            raise ValueError("All links must have the same start radius to create a branching module.")
        
        outer += link.manifold(numSides=numSides, hullBends=hullBends, stabilize=True, wallThickness=None,
                               maxSectionAngle=maxSectionAngle)
        inner += link.manifold(startRadius=startRadius - wallThickness, 
                               endRadius=link.endRadius - link.wallThickness,
                               numSides=numSides, hullBends=hullBends, stabilize=True, wallThickness=None,
                               maxSectionAngle=maxSectionAngle, extendBackward=link.DISTANCE_EPSILON,
                               extendForward=link.DISTANCE_EPSILON)
    tube = outer - inner

    # Create the inset at the base
    holeSlicer = m3d.Manifold()
    holeAnglesDegrees = np.linspace(0, 360, numHoles, endpoint=False)
    connectionLength = 2 * holeDiameter
    for angle in holeAnglesDegrees:
        hole = m3d.Manifold.cylinder(height=r+DISTANCE_EPSILON, 
                                         radius_low=holeDiameter/2, 
                                         radius_high=holeDiameter/2, 
                                         circular_segments=numSides)
        hole = hole.rotate((0,90,0)).rotate((0,0,angle))
        holeSlicer += hole
    
    inset = m3d.Manifold.cylinder(height=2*connectionLength+2*DISTANCE_EPSILON, 
                                       radius_low=startRadius-wallThickness+DISTANCE_EPSILON, 
                                       radius_high=startRadius-wallThickness+DISTANCE_EPSILON,
                                       circular_segments=numSides)
    tube -= inset.translate((0,0,-connectionLength-DISTANCE_EPSILON)).rotate((0,90,0)).transform(startPose.A[:3,:])
    tube -= holeSlicer.translate((0,0,holeDiameter)).rotate((0,90,0)).transform(startPose.A[:3,:])

    # Create the outsets at each link end
    for link in links:
        outset = m3d.Manifold.cylinder(height=2*connectionLength, 
                                       radius_low=link.endRadius-wallThickness+DISTANCE_EPSILON, 
                                       radius_high=link.endRadius-wallThickness+DISTANCE_EPSILON,
                                       circular_segments=numSides)
        outset -= m3d.Manifold.cylinder(height=2*connectionLength, 
                                       radius_low=link.endRadius-2*wallThickness, 
                                       radius_high=link.endRadius-2*wallThickness,
                                       circular_segments=numSides)
        outset -= holeSlicer.translate((0,0,3*holeDiameter))
        outset = outset.translate((0,0,-connectionLength)).rotate((0,90,0)).transform(link.EndDubinsPose.A[:3,:])
        tube += outset
    
    return tube
        



class PrintedKinematicTree(KinematicTree):
    """KinematicTree constrained to PrintedTube fabrication"""
    _fabrication_type = PrintedTube  # Class-level fabrication type constraint
    Links: list[PrintedLinkCSC]  # Type annotation override for proper type checking
    
    def __init__(self, root : Joint, maxAnglePerElbow : float = np.pi/10,
                 joints : Optional[list[Joint]] = None, links : Optional[list[LinkCSC]] = None, 
                 parents : Optional[list[int]] = None, children : Optional[list[list[int]]] = None, 
                 boundingBall : Optional[Ball] = None):
        super().__init__(root=root, maxAnglePerElbow=maxAnglePerElbow,
                         joints=joints, links=links, parents=parents, children=children, boundingBall=boundingBall)
    
    def _get_link_constructor(self):
        """Return callable that creates PrintedLinkCSC with proper fabrication parameters"""
        def make_printed_link(r, start_pose, end_pose, max_angle_per_elbow, path=None, epsilon=0.01):
            # Infer fabrication parameters from existing joints
            wallThickness = 5
            holeDiameter = 2
            numHoles = 4
            for joint in self.Joints:
                if isinstance(joint, PrintedTube):
                    wallThickness = joint.wallThickness
                    holeDiameter = joint.holeDiameter
                    numHoles = joint.numHoles
                    break
            return PrintedLinkCSC(r, start_pose, end_pose, wallThickness, 
                                  holeDiameter, numHoles, max_angle_per_elbow, path, epsilon)
        return make_printed_link
    
    def getLinkModules(self, numSides : int = 50, hullBends : bool = False,
                          maxSectionAngle : float = np.pi/10) -> dict[int, m3d.Manifold]:
        """Generate connectable branching modules for all joints with children"""
        branchingModules = {}
        for jointIndex, childIndices in enumerate(self.Children):
            if len(childIndices) >= 1:
                links = [link for childIndex in childIndices if (link := self.Links[childIndex]).path.length > link.DISTANCE_EPSILON]
                branchingModuleManifold = branchingModule(links, numSides, hullBends, maxSectionAngle)
                branchingModules[jointIndex] = branchingModuleManifold
        return branchingModules
    
    def saveLinkModules(self, baseFilename : str, numSides : int = 50, hullBends : bool = False,
                          maxSectionAngle : float = np.pi/10) -> None:
        """Save connectable branching modules for all joints with children to files"""
        branchingModules = self.getLinkModules(numSides, hullBends, maxSectionAngle)
        for linkIndex, module in branchingModules.items():
            filename = f"{baseFilename}_link{linkIndex}.stl"
            mesh_data = module.to_mesh()
            vertices = mesh_data.vert_properties[:, :3]  # Get XYZ coordinates
            faces = mesh_data.tri_verts
            tri_mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
            tri_mesh.export(filename)
    
    def showLinkModules(self, numSides : int = 50, hullBends : bool = False,
                          maxSectionAngle : float = np.pi/10, block : bool = False) -> None:
        """Display connectable branching modules for all joints with children in a 3D plot window"""
        branchingModules = self.getLinkModules(numSides, hullBends, maxSectionAngle)
        fig = plt.figure()
        ax: Axes3D = fig.add_subplot(projection='3d')
        for linkIndex, module in branchingModules.items():
            mesh_data = module.to_mesh()
            vertices = mesh_data.vert_properties[:, :3]
            triangles = mesh_data.tri_verts
            
            # Create a list of triangle vertex coordinates
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
            faces = [vertices[tri] for tri in triangles]
            mesh_collection = Poly3DCollection(faces, alpha=0.5, 
                                              edgecolor='k',
                                              facecolors='cyan')
            ax.add_collection3d(mesh_collection)
        ax.set_aspect('equal')
        plt.show(block=block)
