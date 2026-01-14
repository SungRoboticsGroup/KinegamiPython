from Joint import *
from TubularPattern import *
from TubularPattern import revoluteColorDefault, revoluteEdgeColorDefault, sphereColorDefault, surfaceOpacityDefault, xColorDefault, yColorDefault, zColorDefault
from geometryHelpers import revoluteColorDefault, revoluteEdgeColorDefault, sphereColorDefault, surfaceOpacityDefault, xColorDefault, yColorDefault, zColorDefault
from LinkCSC import *
from Tube import Tube
from KinematicTree import KinematicTree

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
    Dimensions: 
    r = 30.0  # outer radius of tube in mm
    wallThickness = 5.0  # wall thickness in mm
    holeDiameter = 2.0 # We use M2 bolts
    numHoles = 4
    neutralLength = 76.791 # length in mm from opposite ends of CAD model, excluding the protrusion to inset into the next tube
    """
    def __init__(self, Pose : SE3, version : int | float | str, initialState : float = 0.0):
        if version == 180 or version == "180" or version == 180.0 or version == np.pi:
            maxBendingAngle = np.pi
        elif version == 270 or version == "270" or version == 270.0 or version == 3*np.pi/2:
            maxBendingAngle = 3*np.pi/2
        else:
            raise ValueError("version must be 180 or 270")
        
        PrintedTube.__init__(self, wallThickness=5, holeDiameter=2, numHoles=4)
        TransverseRevolute.__init__(self, r=30, Pose=Pose, 
                                    minAngle=-maxBendingAngle/2, maxAngle=maxBendingAngle/2, 
                                    neutralLength=76.791, initialState=initialState)
    
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
    Dimensions:
    r = 30.0  # outer radius of tube in mm
    wallThickness = 5.0  # wall thickness in mm
    holeDiameter = 2.0 # We use M2 bolts
    numHoles = 4
    length = 64.5 
    """
    def __init__(self, Pose : SE3, version : [int, float, str], initialState : float = 0.0):
        if version == 180 or version == "180" or version == 180.0 or version == np.pi:
            maxBendingAngle = np.pi
        elif version == 270 or version == "270" or version == 270.0 or version == 3*np.pi/2:
            maxBendingAngle = 3*np.pi/2
        else:
            raise ValueError("version must be 180 or 270")
        
        
        PrintedTube.__init__(self, wallThickness=5, holeDiameter=2, numHoles=4)
        CoaxialRevolute.__init__(self, r=30, Pose=Pose, neutralLength=64.5,
                                 minAngle=-maxBendingAngle/2, maxAngle=maxBendingAngle/2, 
                                 initialState=initialState)
    

class PrintedLinkCSC(PrintedTube, LinkCSC):
    def __init__(self, r : float, StartDubinsPose : SE3, EndDubinsPose : SE3,
                 wallThickness : float, holeDiameter : float, numHoles : int,
                 maxAnglePerElbow : float = np.pi/2, path : Optional[PathCSC] = None, 
                 EPSILON : float = 0.01, startRadius : Optional[float] = None, 
                 endRadius : Optional[float] = None):
        PrintedTube.__init__(self, wallThickness, holeDiameter, numHoles)
        LinkCSC.__init__(self, r, StartDubinsPose, EndDubinsPose,
                            maxAnglePerElbow=maxAnglePerElbow, path=path, EPSILON=EPSILON)
        self.startRadius = startRadius if startRadius is not None else r
        self.endRadius = endRadius if endRadius is not None else r
        

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
                """insetSolid -= self.manifold(0.4*insetStartRadius, 0.4*insetEndRadius,
                                  trussNumSides, stabilize=True, wallThickness=None, 
                                  hullBends=hullBends, maxSectionAngle=trussMaxSectionAngle)"""
                        
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


class PrintedKinematicTree(KinematicTree):
    """KinematicTree constrained to PrintedTube fabrication"""
    _fabrication_type = PrintedTube  # Class-level fabrication type constraint
    
    def __init__(self, root : Joint, maxAnglePerElbow : float = np.pi/2,
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
    
    # TODO: make method to export meshes for all links in the tree as connectable modules