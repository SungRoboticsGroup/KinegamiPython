import os
from Joint import *
from TubularPattern import *
from TubularPattern import revoluteColorDefault, revoluteEdgeColorDefault, sphereColorDefault, surfaceOpacityDefault, xColorDefault, yColorDefault, zColorDefault
from geometryHelpers import revoluteColorDefault, revoluteEdgeColorDefault, sphereColorDefault, surfaceOpacityDefault, xColorDefault, yColorDefault, zColorDefault
from LinkCSC import *
from Tube import Tube
from KinematicTree import KinematicTree
from mpl_toolkits.mplot3d import Axes3D

def _circle_points(center, normal, radius, n=32):
    """Generate n points equally spaced around a circle in 3D."""
    if abs(normal[0]) < 0.9:
        u = np.cross(normal, [1, 0, 0])
    else:
        u = np.cross(normal, [0, 1, 0])
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return center + radius * (np.cos(angles)[:, None] * u + np.sin(angles)[:, None] * v)

def _circle_verts_and_tris(center, normal, radius, n=32):
    """Generate vertices and triangle indices for a filled circle disc in 3D (for pyqtgraph)."""
    if abs(normal[0]) < 0.9:
        u = np.cross(normal, [1, 0, 0])
    else:
        u = np.cross(normal, [0, 1, 0])
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)
    angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
    ring = center + radius * (np.cos(angles)[:, None] * u + np.sin(angles)[:, None] * v)
    verts = np.vstack([center.reshape(1, 3), ring])  # index 0 = center
    tris = np.array([[0, i + 1, (i % n) + 2] for i in range(n - 1)] + [[0, n, 1]])
    return verts, tris

class PrintedTube(Tube):
    # Fabrication parameters — edit these in one place for all PrintedTube subclasses
    WALL_THICKNESS = 3.0   # wall thickness in mm
    HOLE_DIAMETER = 3.0    # We use M3 bolts
    NUM_HOLES = 4
    LOOSE_FIT_TOLERANCE = 0.1   # amount in mm to increase radius for loose fit (e.g. clearance holes)
    TIGHT_FIT_TOLERANCE = 0.05  # amount in mm to reduce radius for tight fit (e.g. friction-fit inserts)

    def __init__(self):
        pass

    @property
    def wallThickness(self):
        return self.WALL_THICKNESS

    @property
    def holeDiameter(self):
        return self.HOLE_DIAMETER

    @property
    def numHoles(self):
        return self.NUM_HOLES

    @property
    def looseFitTolerance(self):
        return self.LOOSE_FIT_TOLERANCE

    @property
    def tightFitTolerance(self):
        return self.TIGHT_FIT_TOLERANCE
    
class TransverseRDS3225(PrintedTube, TransverseRevolute):
    """
    RDS3225 Servo Motor with brackets and 3D-printed parts attached to make it attach transversely to tubes.
    Can be set up for the 180 or 270 degree servo, but the 270 version actually has smaller 
    range of motion (214.5 degrees) due to collisions of the proximal and distal attachments. 
    Dimensions are from our CAD model of the 3D printed joint, in mm.
    In the future this parameterized model could be written in the manifold library so wallThickness, holeDiameter, 
    and numHoles could be adjusted as needed. Radius r is as small as possible to bound the servo.

    Angles are in radians where zero is facing forward (x axis of the Pose, i.e., aligning the incoming and outgoing tubes).
    The states are [-np.pi/2, np.pi/2] for 180 degree version, or [-3*np.pi/4, 3*np.pi/4] for 270 degree version.
    When using this in a robot, attach the servo horn and brackets to face forward at 90 or 135 degrees respectively, 
    and then subtract 90 or 135 degrees from the joint state when commanding the servo.

    
    """
    # Dimensions from CAD model (mm) — update these when CAD changes
    R = 33.0              # outer radius of tube in mm
    NEUTRAL_LENGTH = 93.783  # length in mm from opposite ends of CAD model, excluding the protrusion to inset into the next tube

    def __init__(self, Pose : SE3, version : int | float | str, initialState : float = 0.0):
        if version == 180 or version == "180" or version == 180.0 or version == np.pi:
            maxBendingAngle = np.pi
        elif version == 270 or version == "270" or version == 270.0 or version == 3*np.pi/2:
            # 2*170.25 degrees is the largest range we can achieve without the proximal and distal attachments 
            # colliding with each other (based on CAD model measurements) — this is less than 270 degrees but still more than 180 degrees
            maxBendingAngle = 2 * np.deg2rad(107.25)
        else:
            raise ValueError("version must be 180 or 270")
        
        PrintedTube.__init__(self)
        TransverseRevolute.__init__(self, r=self.R, Pose=Pose, neutralLength=self.NEUTRAL_LENGTH, 
                                    minAngle=-maxBendingAngle/2, maxAngle=maxBendingAngle/2, checkCircleOverlap=True,
                                    initialState=initialState)
    
    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=revoluteColorDefault, edgeColor=revoluteEdgeColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=True,
             axisScale=jointAxisScaleDefault, showPoses=True):
        """Override to display a simplified servo icon instead of the default revolute visualization.
        
        Draws a 40x40x18 box (narrow along z-axis of Pose) offset -11 in x from Pose center,
        plus 9mm cylinders extending from the proximal and distal frames toward the Pose center.
        Calls Joint.addToPlot (skipping TransverseRevolute/Revolute surface drawing).
        """
        # Call Joint.addToPlot directly to get frames/axis without revolute surface
        plotHandles = Joint.addToPlot(self, ax=ax, xColor=xColor, yColor=yColor, zColor=zColor, 
                          proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor, 
                          sphereColor=sphereColor, showSphere=showSphere,
                          surfaceColor=surfaceColor, edgeColor=edgeColor,
                          surfaceOpacity=surfaceOpacity, showSurface=False, showAxis=showAxis,
                          axisScale=axisScale, showPoses=showPoses)
        
        if showSurface:
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
            
            xhat = self.Pose.R[:, 0]  # path direction
            yhat = self.Pose.R[:, 1]
            zhat = self.Pose.R[:, 2]  # rotation axis
            
            # Box dimensions: 40 along x, 20 along y, 40 along z
            bx, by, bz = 40.0, 20.0, 40.0
            box_center = self.Pose.t + (-11.0) * xhat
            
            # 8 corners of the box in local frame, then transform
            hx, hy, hz = bx / 2, by / 2, bz / 2
            local_corners = np.array([
                [-hx, -hy, -hz], [+hx, -hy, -hz], [+hx, +hy, -hz], [-hx, +hy, -hz],
                [-hx, -hy, +hz], [+hx, -hy, +hz], [+hx, +hy, +hz], [-hx, +hy, +hz],
            ])
            # Transform to world: columns of R are [xhat, yhat, zhat]
            R = np.column_stack([xhat, yhat, zhat])
            corners = (R @ local_corners.T).T + box_center
            
            # 6 faces of the box (indices into corners)
            face_indices = [
                [0, 1, 2, 3],  # -z face
                [4, 5, 6, 7],  # +z face
                [0, 1, 5, 4],  # -y face
                [2, 3, 7, 6],  # +y face
                [0, 3, 7, 4],  # -x face
                [1, 2, 6, 5],  # +x face
            ]
            faces = [[corners[i] for i in face] for face in face_indices]
            box_collection = Poly3DCollection(faces, alpha=surfaceOpacity,
                                              facecolors='black',
                                              edgecolors=edgeColor)
            ax.add_collection3d(box_collection)
            
            # Cylinders: 9mm long, radius self.r, with solid end caps
            cyl_length = 9.0
            numCapPoints = 32
            proximal_pos = self.ProximalFrame().t
            proximal_cyl = Cylinder(self.r, proximal_pos, xhat, cyl_length)
            proximal_cyl.addToPlot(ax, color=surfaceColor, alpha=surfaceOpacity, edgeColor=edgeColor)
            
            distalFrame = self.DistalFrame()
            distal_pos = distalFrame.t
            distal_xhat = distalFrame.R[:, 0]
            distal_cyl = Cylinder(self.r, distal_pos, -distal_xhat, cyl_length)
            distal_cyl.addToPlot(ax, color=surfaceColor, alpha=surfaceOpacity, edgeColor=edgeColor)
            
            # 4 vertical lines along each cylinder at ±y/±z in the respective frames
            prox_inner = proximal_pos + cyl_length * xhat
            for direction in [yhat, -yhat, zhat, -zhat]:
                line_pts = np.array([proximal_pos + self.r * direction, prox_inner + self.r * direction])
                ax.plot(line_pts[:, 0], line_pts[:, 1], line_pts[:, 2], color=edgeColor, linewidth=1)
            distal_inner = distal_pos - cyl_length * distal_xhat
            distal_yhat = distalFrame.R[:, 1]
            distal_zhat = distalFrame.R[:, 2]
            for direction in [distal_yhat, -distal_yhat, distal_zhat, -distal_zhat]:
                line_pts = np.array([distal_pos + self.r * direction, distal_inner + self.r * direction])
                ax.plot(line_pts[:, 0], line_pts[:, 1], line_pts[:, 2], color=edgeColor, linewidth=1)
            
            # End caps for solid cylinders
            for cap_center, cap_normal in [
                (proximal_pos + cyl_length * xhat, xhat),       # inner cap of proximal cyl
                (distal_pos - cyl_length * distal_xhat, -distal_xhat),  # inner cap of distal cyl
            ]:
                cap_pts = _circle_points(cap_center, cap_normal, self.r, numCapPoints)
                cap_face = Poly3DCollection([cap_pts], alpha=surfaceOpacity,
                                            facecolors=surfaceColor, edgecolors=edgeColor)
                ax.add_collection3d(cap_face)
            
            # Bracket dimensions
            bracket_color = 'gray'
            bracket_lw = 2
            bracket_hy = 10.0  # half of 20 units wide in y
            
            # Proximal [-bracket: connects ±z edges of box's -x face to inner end of proximal cylinder
            box_prox_face = box_center - hx * xhat  # center of box's -x face
            prox_inner = proximal_pos + cyl_length * xhat  # inner end of proximal cylinder
            # Draw two [ outlines at y = ±bracket_hy, plus connecting edges
            for sign_y in [+1, -1]:
                y_offset = sign_y * bracket_hy * yhat
                prox_bracket = np.array([
                    box_prox_face + hz * zhat + y_offset,
                    prox_inner + hz * zhat + y_offset,
                    prox_inner - hz * zhat + y_offset,
                    box_prox_face - hz * zhat + y_offset,
                ])
                ax.plot(prox_bracket[:, 0], prox_bracket[:, 1], prox_bracket[:, 2],
                        color=bracket_color, linewidth=bracket_lw)
            # Connecting edges between the two y sides at the 4 bracket corners
            for z_sign in [+1, -1]:
                for x_pos in [box_prox_face, prox_inner]:
                    corner = x_pos + z_sign * hz * zhat
                    edge = np.array([corner + bracket_hy * yhat, corner - bracket_hy * yhat])
                    ax.plot(edge[:, 0], edge[:, 1], edge[:, 2],
                            color=bracket_color, linewidth=bracket_lw)
            # Rectangular faces connecting +y and -y bracket outlines (3 segments of the [)
            prox_bracket_faces = [
                # Top horizontal: box_prox+z to prox_inner+z
                [box_prox_face + hz*zhat + bracket_hy*yhat, prox_inner + hz*zhat + bracket_hy*yhat,
                 prox_inner + hz*zhat - bracket_hy*yhat, box_prox_face + hz*zhat - bracket_hy*yhat],
                # Vertical: prox_inner+z to prox_inner-z
                [prox_inner + hz*zhat + bracket_hy*yhat, prox_inner - hz*zhat + bracket_hy*yhat,
                 prox_inner - hz*zhat - bracket_hy*yhat, prox_inner + hz*zhat - bracket_hy*yhat],
                # Bottom horizontal: prox_inner-z to box_prox-z
                [prox_inner - hz*zhat + bracket_hy*yhat, box_prox_face - hz*zhat + bracket_hy*yhat,
                 box_prox_face - hz*zhat - bracket_hy*yhat, prox_inner - hz*zhat - bracket_hy*yhat],
            ]
            prox_face_collection = Poly3DCollection(prox_bracket_faces, alpha=surfaceOpacity*0.5,
                                                     facecolors=bracket_color, edgecolors=bracket_color)
            ax.add_collection3d(prox_face_collection)
            
            # Distal [-bracket: connects from axis of motion (Pose.t) to inner end of distal cylinder, in distal frame
            distal_inner = distal_pos - cyl_length * distal_xhat  # inner end of distal cylinder
            distal_yhat = distalFrame.R[:, 1]
            for sign_y in [+1, -1]:
                y_offset = sign_y * bracket_hy * distal_yhat
                distal_bracket = np.array([
                    self.Pose.t + hz * zhat + y_offset,
                    distal_inner + hz * zhat + y_offset,
                    distal_inner - hz * zhat + y_offset,
                    self.Pose.t - hz * zhat + y_offset,
                ])
                ax.plot(distal_bracket[:, 0], distal_bracket[:, 1], distal_bracket[:, 2],
                        color=bracket_color, linewidth=bracket_lw)
            # Connecting edges between the two y sides at the 4 bracket corners
            for z_sign in [+1, -1]:
                for x_pos in [self.Pose.t, distal_inner]:
                    corner = x_pos + z_sign * hz * zhat
                    edge = np.array([corner + bracket_hy * distal_yhat, corner - bracket_hy * distal_yhat])
                    ax.plot(edge[:, 0], edge[:, 1], edge[:, 2],
                            color=bracket_color, linewidth=bracket_lw)
            # Rectangular faces connecting +y and -y bracket outlines (3 segments of the [)
            distal_bracket_faces = [
                # Top horizontal: Pose.t+z to distal_inner+z
                [self.Pose.t + hz*zhat + bracket_hy*distal_yhat, distal_inner + hz*zhat + bracket_hy*distal_yhat,
                 distal_inner + hz*zhat - bracket_hy*distal_yhat, self.Pose.t + hz*zhat - bracket_hy*distal_yhat],
                # Vertical: distal_inner+z to distal_inner-z
                [distal_inner + hz*zhat + bracket_hy*distal_yhat, distal_inner - hz*zhat + bracket_hy*distal_yhat,
                 distal_inner - hz*zhat - bracket_hy*distal_yhat, distal_inner + hz*zhat - bracket_hy*distal_yhat],
                # Bottom horizontal: distal_inner-z to Pose.t-z
                [distal_inner - hz*zhat + bracket_hy*distal_yhat, self.Pose.t - hz*zhat + bracket_hy*distal_yhat,
                 self.Pose.t - hz*zhat - bracket_hy*distal_yhat, distal_inner - hz*zhat - bracket_hy*distal_yhat],
            ]
            distal_face_collection = Poly3DCollection(distal_bracket_faces, alpha=surfaceOpacity*0.5,
                                                       facecolors=bracket_color, edgecolors=bracket_color)
            ax.add_collection3d(distal_face_collection)
        
        return plotHandles

    
    def addToWidget(self, widget, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                    proximalColor=proximalColorDefault, centerColor=centerColorDefault, distalColor=distalColorDefault,
                    sphereColor=sphereColorDefault, showSphere=False, 
                    surfaceColor=revoluteColorDefault, 
                    showSurface=True, showAxis=True, axisScale=jointAxisScaleDefault, showPoses=True, poseAxisScaleMultipler=None):
        #Override to display a simplified servo icon in the pyqtgraph widget.
        
        #Mirrors the addToPlot method: box, cylinders with inner caps, and [-brackets.
        #Calls Joint.addToWidget directly (skipping Revolute/TransverseRevolute surface).
        
        import pyqtgraph.opengl as gl
        from style import revoluteColorList
        
        # Call Joint.addToWidget directly for poses/axis without revolute surface
        Joint.addToWidget(self, widget=widget, xColor=xColor, yColor=yColor, zColor=zColor, 
                          proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor, 
                          sphereColor=sphereColor, showSphere=showSphere,
                          surfaceColor=surfaceColor, showSurface=False, showAxis=showAxis,
                          axisScale=axisScale, showPoses=showPoses, poseAxisScaleMultipler=poseAxisScaleMultipler)
        
        if showSurface:
            xhat = self.Pose.R[:, 0]
            yhat = self.Pose.R[:, 1]
            zhat = self.Pose.R[:, 2]
            
            # Box dimensions: 40 along x, 20 along y, 40 along z
            bx, by, bz = 40.0, 20.0, 40.0
            box_center = self.Pose.t + (-11.0) * xhat
            hx, hy, hz = bx / 2, by / 2, bz / 2
            
            # 8 corners of the box
            local_corners = np.array([
                [-hx, -hy, -hz], [+hx, -hy, -hz], [+hx, +hy, -hz], [-hx, +hy, -hz],
                [-hx, -hy, +hz], [+hx, -hy, +hz], [+hx, +hy, +hz], [-hx, +hy, +hz],
            ])
            R = np.column_stack([xhat, yhat, zhat])
            corners = (R @ local_corners.T).T + box_center
            
            # 12 triangles for 6 box faces
            box_tris = np.array([
                [0,1,2], [0,2,3],  # -z
                [4,5,6], [4,6,7],  # +z
                [0,1,5], [0,5,4],  # -y
                [2,3,7], [2,7,6],  # +y
                [0,3,7], [0,7,4],  # -x
                [1,2,6], [1,6,5],  # +x
            ])
            box_mesh = gl.GLMeshItem(vertexes=corners, faces=box_tris,
                                      color=(0, 0, 0, 0.6), smooth=False, 
                                      drawEdges=True, edgeColor=(0.3, 0.3, 0.3, 1))
            box_mesh.setGLOptions('opaque')
            box_mesh.setObjectName("Joint")
            widget.plot_widget.addItem(box_mesh)
            
            # Cylinders: 9mm long, radius self.r
            cyl_length = 9.0
            numCylPoints = 32
            proximal_pos = self.ProximalFrame().t
            proximal_cyl = Cylinder(self.r, proximal_pos, xhat, cyl_length)
            proximal_cyl.addToWidget(widget, color_list=revoluteColorList, is_joint=True, opaque=True)
            
            distalFrame = self.DistalFrame()
            distal_pos = distalFrame.t
            distal_xhat = distalFrame.R[:, 0]
            distal_cyl = Cylinder(self.r, distal_pos, -distal_xhat, cyl_length)
            distal_cyl.addToWidget(widget, color_list=revoluteColorList, is_joint=True, opaque=True)
            
            # 4 vertical lines along each cylinder at ±y/±z in the respective frames
            line_color = (0.3, 0.3, 0.3, 1.0)
            prox_inner = proximal_pos + cyl_length * xhat
            for direction in [yhat, -yhat, zhat, -zhat]:
                line_pts = np.array([proximal_pos + self.r * direction, prox_inner + self.r * direction])
                line = gl.GLLinePlotItem(pos=line_pts, color=line_color, width=2, antialias=True)
                widget.plot_widget.addItem(line)
            distal_inner_pos = distal_pos - cyl_length * distal_xhat
            distal_yhat = distalFrame.R[:, 1]
            distal_zhat = distalFrame.R[:, 2]
            for direction in [distal_yhat, -distal_yhat, distal_zhat, -distal_zhat]:
                line_pts = np.array([distal_pos + self.r * direction, distal_inner_pos + self.r * direction])
                line = gl.GLLinePlotItem(pos=line_pts, color=line_color, width=2, antialias=True)
                widget.plot_widget.addItem(line)
            
            # Inner end caps
            for cap_center, cap_normal in [
                (proximal_pos + cyl_length * xhat, xhat),
                (distal_pos - cyl_length * distal_xhat, -distal_xhat),
            ]:
                cap_verts, cap_tris = _circle_verts_and_tris(cap_center, cap_normal, self.r, numCylPoints)
                cap_mesh = gl.GLMeshItem(vertexes=cap_verts, faces=cap_tris,
                                          color=tuple(revoluteColorList), smooth=True)
                cap_mesh.setGLOptions('opaque')
                cap_mesh.setObjectName("Joint")
                widget.plot_widget.addItem(cap_mesh)
            
            # Bracket dimensions
            bracket_color = (0.5, 0.5, 0.5, 0.8)
            bracket_line_color = (0.5, 0.5, 0.5, 1.0)
            bracket_hy = 10.0
            
            # Proximal [-bracket
            box_prox_face = box_center - hx * xhat
            prox_inner = proximal_pos + cyl_length * xhat
            
            # Bracket outlines at y = ±bracket_hy
            for sign_y in [+1, -1]:
                y_offset = sign_y * bracket_hy * yhat
                pts = np.array([
                    box_prox_face + hz * zhat + y_offset,
                    prox_inner + hz * zhat + y_offset,
                    prox_inner - hz * zhat + y_offset,
                    box_prox_face - hz * zhat + y_offset,
                ])
                line = gl.GLLinePlotItem(pos=pts, color=bracket_line_color, width=2, antialias=True)
                widget.plot_widget.addItem(line)
            # Connecting edges
            for z_sign in [+1, -1]:
                for x_pos in [box_prox_face, prox_inner]:
                    corner = x_pos + z_sign * hz * zhat
                    edge = np.array([corner + bracket_hy * yhat, corner - bracket_hy * yhat])
                    line = gl.GLLinePlotItem(pos=edge, color=bracket_line_color, width=2, antialias=True)
                    widget.plot_widget.addItem(line)
            # Bracket face quads (as triangulated meshes)
            prox_bracket_quads = [
                [box_prox_face + hz*zhat + bracket_hy*yhat, prox_inner + hz*zhat + bracket_hy*yhat,
                 prox_inner + hz*zhat - bracket_hy*yhat, box_prox_face + hz*zhat - bracket_hy*yhat],
                [prox_inner + hz*zhat + bracket_hy*yhat, prox_inner - hz*zhat + bracket_hy*yhat,
                 prox_inner - hz*zhat - bracket_hy*yhat, prox_inner + hz*zhat - bracket_hy*yhat],
                [prox_inner - hz*zhat + bracket_hy*yhat, box_prox_face - hz*zhat + bracket_hy*yhat,
                 box_prox_face - hz*zhat - bracket_hy*yhat, prox_inner - hz*zhat - bracket_hy*yhat],
            ]
            for quad in prox_bracket_quads:
                verts = np.array(quad)
                tris = np.array([[0, 1, 2], [0, 2, 3]])
                mesh = gl.GLMeshItem(vertexes=verts, faces=tris, color=bracket_color, smooth=False)
                mesh.setGLOptions('translucent')
                mesh.setObjectName("Joint")
                widget.plot_widget.addItem(mesh)
            
            # Distal [-bracket
            distal_inner = distal_pos - cyl_length * distal_xhat
            distal_yhat = distalFrame.R[:, 1]
            
            for sign_y in [+1, -1]:
                y_offset = sign_y * bracket_hy * distal_yhat
                pts = np.array([
                    self.Pose.t + hz * zhat + y_offset,
                    distal_inner + hz * zhat + y_offset,
                    distal_inner - hz * zhat + y_offset,
                    self.Pose.t - hz * zhat + y_offset,
                ])
                line = gl.GLLinePlotItem(pos=pts, color=bracket_line_color, width=2, antialias=True)
                widget.plot_widget.addItem(line)
            for z_sign in [+1, -1]:
                for x_pos in [self.Pose.t, distal_inner]:
                    corner = x_pos + z_sign * hz * zhat
                    edge = np.array([corner + bracket_hy * distal_yhat, corner - bracket_hy * distal_yhat])
                    line = gl.GLLinePlotItem(pos=edge, color=bracket_line_color, width=2, antialias=True)
                    widget.plot_widget.addItem(line)
            distal_bracket_quads = [
                [self.Pose.t + hz*zhat + bracket_hy*distal_yhat, distal_inner + hz*zhat + bracket_hy*distal_yhat,
                 distal_inner + hz*zhat - bracket_hy*distal_yhat, self.Pose.t + hz*zhat - bracket_hy*distal_yhat],
                [distal_inner + hz*zhat + bracket_hy*distal_yhat, distal_inner - hz*zhat + bracket_hy*distal_yhat,
                 distal_inner - hz*zhat - bracket_hy*distal_yhat, distal_inner + hz*zhat - bracket_hy*distal_yhat],
                [distal_inner - hz*zhat + bracket_hy*distal_yhat, self.Pose.t - hz*zhat + bracket_hy*distal_yhat,
                 self.Pose.t - hz*zhat - bracket_hy*distal_yhat, distal_inner - hz*zhat - bracket_hy*distal_yhat],
            ]
            for quad in distal_bracket_quads:
                verts = np.array(quad)
                tris = np.array([[0, 1, 2], [0, 2, 3]])
                mesh = gl.GLMeshItem(vertexes=verts, faces=tris, color=bracket_color, smooth=False)
                mesh.setGLOptions('translucent')
                mesh.setObjectName("Joint")
                widget.plot_widget.addItem(mesh)
        
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
    NEUTRAL_LENGTH = 62.4  # length in mm

    def __init__(self, Pose : SE3, version : [int, float, str], initialState : float = 0.0):
        if version == 180 or version == "180" or version == 180.0 or version == np.pi:
            maxBendingAngle = np.pi
        elif version == 270 or version == "270" or version == 270.0 or version == 3*np.pi/2:
            maxBendingAngle = 3*np.pi/2
        else:
            raise ValueError("version must be 180 or 270")
        
        PrintedTube.__init__(self)
        CoaxialRevolute.__init__(self, r=self.R, Pose=Pose, neutralLength=self.NEUTRAL_LENGTH,
                                 minAngle=-maxBendingAngle/2, maxAngle=maxBendingAngle/2, 
                                 initialState=initialState)

    def addToPlot(self, ax, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
             proximalColor='c', centerColor='m', distalColor='y',
             sphereColor=sphereColorDefault, showSphere=False, 
             surfaceColor=revoluteColorDefault, edgeColor=revoluteEdgeColorDefault,
             surfaceOpacity=surfaceOpacityDefault, showSurface=True, showAxis=True,
             axisScale=jointAxisScaleDefault, showPoses=True):
        """Override to display a simplified coaxial servo icon instead of the default revolute visualization.
        
        Draws a 40x20x40 black box (the servo) centered at (11, 0, 1.95) in the Pose frame,
        a proximal cylinder of radius self.r and length 9mm,
        and a distal cylinder of radius (self.r - self.wallThickness) and length 12.9mm 
        oriented in the distal frame so it rotates with the joint state.
        Calls Joint.addToPlot directly (skipping Revolute/CoaxialRevolute surface drawing).
        """
        from mpl_toolkits.mplot3d.art3d import Poly3DCollection

        # Call Joint.addToPlot directly to get frames/axis without revolute surface
        plotHandles = Joint.addToPlot(self, ax=ax, xColor=xColor, yColor=yColor, zColor=zColor, 
                          proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor, 
                          sphereColor=sphereColor, showSphere=showSphere,
                          surfaceColor=surfaceColor, edgeColor=edgeColor,
                          surfaceOpacity=surfaceOpacity, showSurface=False, showAxis=showAxis,
                          axisScale=axisScale, showPoses=showPoses)
        
        if showSurface:
            xhat = self.Pose.R[:, 0]
            yhat = self.Pose.R[:, 1]
            zhat = self.Pose.R[:, 2]  # path direction for coaxial
            
            # --- Servo box: 40x20x40 centered at (11, 0, -1.95) in Pose frame ---
            bx, by, bz = 40.0, 20.0, 40.0
            box_center = self.Pose.t + 11.0 * xhat + (-1.95) * zhat
            hx, hy, hz = bx / 2, by / 2, bz / 2
            
            local_corners = np.array([
                [-hx, -hy, -hz], [+hx, -hy, -hz], [+hx, +hy, -hz], [-hx, +hy, -hz],
                [-hx, -hy, +hz], [+hx, -hy, +hz], [+hx, +hy, +hz], [-hx, +hy, +hz],
            ])
            R = np.column_stack([xhat, yhat, zhat])
            corners = (R @ local_corners.T).T + box_center
            
            face_indices = [
                [0, 1, 2, 3],  # -z face
                [4, 5, 6, 7],  # +z face
                [0, 1, 5, 4],  # -y face
                [2, 3, 7, 6],  # +y face
                [0, 3, 7, 4],  # -x face
                [1, 2, 6, 5],  # +x face
            ]
            faces = [[corners[i] for i in face] for face in face_indices]
            box_collection = Poly3DCollection(faces, alpha=surfaceOpacity,
                                              facecolors='black',
                                              edgecolors=edgeColor)
            ax.add_collection3d(box_collection)
            
            # --- Proximal cylinder: radius self.r, 9mm long, from ProximalFrame along zhat ---
            prox_cyl_length = 9.0
            numCapPoints = 32
            proximal_pos = self.ProximalFrame().t
            proximal_cyl = Cylinder(self.r, proximal_pos, zhat, prox_cyl_length)
            proximal_cyl.addToPlot(ax, color=surfaceColor, alpha=surfaceOpacity, edgeColor=edgeColor)
            
            # Inner cap for proximal cylinder
            prox_cap_center = proximal_pos + prox_cyl_length * zhat
            prox_cap_pts = _circle_points(prox_cap_center, zhat, self.r, numCapPoints)
            prox_cap_face = Poly3DCollection([prox_cap_pts], alpha=surfaceOpacity,
                                              facecolors=surfaceColor, edgecolors=edgeColor)
            ax.add_collection3d(prox_cap_face)
            
            # 4 vertical lines along the proximal cylinder at ±x/±y in the Pose frame
            for direction in [xhat, -xhat, yhat, -yhat]:
                line_start = proximal_pos + self.r * direction
                line_end = prox_cap_center + self.r * direction
                line_pts = np.array([line_start, line_end])
                ax.plot(line_pts[:, 0], line_pts[:, 1], line_pts[:, 2],
                        color=edgeColor, linewidth=1)
            
            # --- Distal cylinder: radius (r - wallThickness), 12.9mm long, in distal frame ---
            dist_cyl_length = 12.9
            distalFrame = self.DistalFrame()
            distal_pos = distalFrame.t
            distal_zhat = distalFrame.R[:, 2]  # path direction in distal frame
            distal_r = self.r - self.wallThickness
            distal_cyl = Cylinder(distal_r, distal_pos, -distal_zhat, dist_cyl_length)
            distal_cyl.addToPlot(ax, color=surfaceColor, alpha=surfaceOpacity, edgeColor=edgeColor)
            
            # Inner cap for distal cylinder
            dist_cap_center = distal_pos - dist_cyl_length * distal_zhat
            dist_cap_pts = _circle_points(dist_cap_center, -distal_zhat, distal_r, numCapPoints)
            dist_cap_face = Poly3DCollection([dist_cap_pts], alpha=surfaceOpacity,
                                              facecolors=surfaceColor, edgecolors=edgeColor)
            ax.add_collection3d(dist_cap_face)
            
            # 4 vertical lines along the distal cylinder at ±x/±y in the distal frame
            distal_xhat = distalFrame.R[:, 0]
            distal_yhat = distalFrame.R[:, 1]
            for direction in [distal_xhat, -distal_xhat, distal_yhat, -distal_yhat]:
                line_start = distal_pos + distal_r * direction
                line_end = dist_cap_center + distal_r * direction
                line_pts = np.array([line_start, line_end])
                ax.plot(line_pts[:, 0], line_pts[:, 1], line_pts[:, 2],
                        color=edgeColor, linewidth=1)
            
            # Mark proximal and distal frame origins with labeled scatter points
            ax.scatter(*proximal_pos, color='cyan', s=60, zorder=5)
            ax.text(*proximal_pos, '  Prox', color='cyan', fontsize=8)
            ax.scatter(*distal_pos, color='yellow', s=60, zorder=5)
            ax.text(*distal_pos, '  Dist', color='yellow', fontsize=8)
        
        return plotHandles

    
    def addToWidget(self, widget, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, 
                    proximalColor=proximalColorDefault, centerColor=centerColorDefault, distalColor=distalColorDefault,
                    sphereColor=sphereColorDefault, showSphere=False, 
                    surfaceColor=revoluteColorDefault, 
                    showSurface=True, showAxis=True, axisScale=jointAxisScaleDefault, showPoses=True, poseAxisScaleMultipler=None):
        #Override to display a simplified coaxial servo icon in the pyqtgraph widget.
        
        #Mirrors the addToPlot method: black box, proximal cylinder (r), distal cylinder (r-wallThickness),
        #each with inner end caps.
        #Calls Joint.addToWidget directly (skipping Revolute/CoaxialRevolute surface).
        import pyqtgraph.opengl as gl
        from style import revoluteColorList
        
        # Call Joint.addToWidget directly for poses/axis without revolute surface
        Joint.addToWidget(self, widget=widget, xColor=xColor, yColor=yColor, zColor=zColor, 
                          proximalColor=proximalColor, centerColor=centerColor, distalColor=distalColor, 
                          sphereColor=sphereColor, showSphere=showSphere,
                          surfaceColor=surfaceColor, showSurface=False, showAxis=showAxis,
                          axisScale=axisScale, showPoses=showPoses, poseAxisScaleMultipler=poseAxisScaleMultipler)
        
        if showSurface:
            xhat = self.Pose.R[:, 0]
            yhat = self.Pose.R[:, 1]
            zhat = self.Pose.R[:, 2]  # path direction for coaxial
            
            # --- Servo box: 40x20x40 centered at (11, 0, -1.95) in Pose frame ---
            bx, by, bz = 40.0, 20.0, 40.0
            box_center = self.Pose.t + 11.0 * xhat + (-1.95) * zhat
            hx, hy, hz = bx / 2, by / 2, bz / 2
            
            local_corners = np.array([
                [-hx, -hy, -hz], [+hx, -hy, -hz], [+hx, +hy, -hz], [-hx, +hy, -hz],
                [-hx, -hy, +hz], [+hx, -hy, +hz], [+hx, +hy, +hz], [-hx, +hy, +hz],
            ])
            R = np.column_stack([xhat, yhat, zhat])
            corners = (R @ local_corners.T).T + box_center
            
            box_tris = np.array([
                [0,1,2], [0,2,3],  # -z
                [4,5,6], [4,6,7],  # +z
                [0,1,5], [0,5,4],  # -y
                [2,3,7], [2,7,6],  # +y
                [0,3,7], [0,7,4],  # -x
                [1,2,6], [1,6,5],  # +x
            ])
            box_mesh = gl.GLMeshItem(vertexes=corners, faces=box_tris,
                                      color=(0, 0, 0, 0.6), smooth=False, 
                                      drawEdges=True, edgeColor=(0.3, 0.3, 0.3, 1))
            box_mesh.setGLOptions('opaque')
            box_mesh.setObjectName("Joint")
            widget.plot_widget.addItem(box_mesh)
            
            # --- Proximal cylinder: radius self.r, 9mm long, from ProximalFrame along zhat ---
            prox_cyl_length = 9.0
            numCylPoints = 32
            proximal_pos = self.ProximalFrame().t
            proximal_cyl = Cylinder(self.r, proximal_pos, zhat, prox_cyl_length)
            proximal_cyl.addToWidget(widget, color_list=revoluteColorList, is_joint=True, opaque=True)
            
            # Inner cap for proximal cylinder
            prox_cap_center = proximal_pos + prox_cyl_length * zhat
            prox_cap_verts, prox_cap_tris = _circle_verts_and_tris(prox_cap_center, zhat, self.r, numCylPoints)
            prox_cap_mesh = gl.GLMeshItem(vertexes=prox_cap_verts, faces=prox_cap_tris,
                                           color=tuple(revoluteColorList), smooth=True)
            prox_cap_mesh.setGLOptions('opaque')
            prox_cap_mesh.setObjectName("Joint")
            widget.plot_widget.addItem(prox_cap_mesh)
            
            # 4 vertical lines along the proximal cylinder at ±x/±y in the Pose frame
            line_color = (0.3, 0.3, 0.3, 1.0)
            for direction in [xhat, -xhat, yhat, -yhat]:
                line_start = proximal_pos + self.r * direction
                line_end = prox_cap_center + self.r * direction
                line_pts = np.array([line_start, line_end])
                line = gl.GLLinePlotItem(pos=line_pts, color=line_color, width=2, antialias=True)
                widget.plot_widget.addItem(line)
            
            # --- Distal cylinder: radius (r - wallThickness), 12.9mm long, in distal frame ---
            dist_cyl_length = 12.9
            distalFrame = self.DistalFrame()
            distal_pos = distalFrame.t
            distal_zhat = distalFrame.R[:, 2]  # path direction in distal frame
            distal_r = self.r - self.wallThickness
            distal_cyl = Cylinder(distal_r, distal_pos, -distal_zhat, dist_cyl_length)
            distal_cyl.addToWidget(widget, color_list=revoluteColorList, is_joint=True, opaque=True)
            
            # Inner cap for distal cylinder
            dist_cap_center = distal_pos - dist_cyl_length * distal_zhat
            dist_cap_verts, dist_cap_tris = _circle_verts_and_tris(dist_cap_center, -distal_zhat, distal_r, numCylPoints)
            dist_cap_mesh = gl.GLMeshItem(vertexes=dist_cap_verts, faces=dist_cap_tris,
                                           color=tuple(revoluteColorList), smooth=True)
            dist_cap_mesh.setGLOptions('opaque')
            dist_cap_mesh.setObjectName("Joint")
            widget.plot_widget.addItem(dist_cap_mesh)
            
            # 4 vertical lines along the distal cylinder at ±x/±y in the distal frame
            distal_xhat = distalFrame.R[:, 0]
            distal_yhat = distalFrame.R[:, 1]
            line_color = (0.3, 0.3, 0.3, 1.0)
            for direction in [distal_xhat, -distal_xhat, distal_yhat, -distal_yhat]:
                line_start = distal_pos + distal_r * direction
                line_end = dist_cap_center + distal_r * direction
                line_pts = np.array([line_start, line_end])
                line = gl.GLLinePlotItem(pos=line_pts, color=line_color, width=2, antialias=True)
                widget.plot_widget.addItem(line)
    
class PrintedHemisphere(PrintedTube, Tip):
    def __init__(self, r : float, Pose : SE3, closesForward : bool, pathIndex : int = 2):
        PrintedTube.__init__(self)
        Tip.__init__(self, r, Pose, length=r, closesForward=closesForward, pathIndex=pathIndex)

    def connectableModule(self, numSides : int = 50) -> m3d.Manifold:
        """Generate a connectable hemisphere module for 3D printing.
        
        Creates a hollow hemisphere with an inset or outset connection at
        the open face, matching the link connection pattern.
        
        Forward tips (end caps) get an inset at the open (proximal) face.
        Backward tips (start caps) get an outset at the open (distal) face.
        """
        connectionLength = 2 * self.holeDiameter
        eps = 0.01 * self.r  # small epsilon for boolean clearance

        # Build a hollow sphere at the origin, then trim and transform
        outer = m3d.Manifold.sphere(self.r, circular_segments=numSides)
        inner = m3d.Manifold.sphere(self.r - self.wallThickness, circular_segments=numSides)
        shell = outer - inner

        # The tree/link connection uses proximal frame for end caps and distal
        # frame for start caps. Use that same open face frame for cap geometry.
        openFrame = self.ProximalDubinsFrame() if self.forward else self.DistalDubinsFrame()
        openCenter = openFrame.t

        # Trim: keep the closed half relative to the open face center.
        # Forward tip closes in +openFrame xhat; backward tip closes in -xhat.
        openAxis = openFrame.R[:,0]
        trim_normal = openAxis if self.forward else -openAxis
        # trim_by_plane keeps the side in the direction of the normal.
        # origin_offset = dot(normal, point_on_plane) positions the plane through openCenter.
        shell = shell.translate(tuple(openCenter))
        origin_offset = float(np.dot(trim_normal, openCenter))
        shell = shell.trim_by_plane(tuple(trim_normal), origin_offset)

        if self.forward:
            # Forward tip: open face is proximal → receives an outset from the
            # preceding link → needs an INSET at the open face.
            # Inset holes use loose fit (bolts clear through).
            holeSlicer = m3d.Manifold()
            holeAnglesDegrees = np.linspace(0, 360, self.numHoles, endpoint=False)
            for angle in holeAnglesDegrees:
                hole = m3d.Manifold.cylinder(height=self.r + eps,
                                             radius_low=self.holeDiameter/2 + self.looseFitTolerance,
                                             radius_high=self.holeDiameter/2 + self.looseFitTolerance,
                                             circular_segments=numSides)
                hole = hole.rotate((0,90,0)).rotate((0,0,angle))
                holeSlicer += hole

            inset = m3d.Manifold.cylinder(
                height=2*connectionLength + 2*eps,
                radius_low=self.r - self.wallThickness + eps,
                radius_high=self.r - self.wallThickness + eps,
                circular_segments=numSides)
            shell -= inset.translate((0,0,-connectionLength-eps)) \
                         .rotate((0,90,0)).transform(openFrame.A[:3,:])
            shell -= holeSlicer.translate((0,0,self.holeDiameter)) \
                               .rotate((0,90,0)).transform(openFrame.A[:3,:])
        else:
            # Backward tip: open face is distal → inserts into the following
            # link's inset → needs an OUTSET at the open face.
            # Outset holes use tight fit (bolts friction-fit in).
            holeSlicer = m3d.Manifold()
            holeAnglesDegrees = np.linspace(0, 360, self.numHoles, endpoint=False)
            for angle in holeAnglesDegrees:
                hole = m3d.Manifold.cylinder(height=self.r + eps,
                                             radius_low=self.holeDiameter/2 + self.tightFitTolerance,
                                             radius_high=self.holeDiameter/2 + self.tightFitTolerance,
                                             circular_segments=numSides)
                hole = hole.rotate((0,90,0)).rotate((0,0,angle))
                holeSlicer += hole

            outsetRadius = self.r - self.wallThickness
            outset = m3d.Manifold.cylinder(
                height=2*connectionLength,
                radius_low=outsetRadius - self.tightFitTolerance,
                radius_high=outsetRadius - self.tightFitTolerance,
                circular_segments=numSides)
            outset += m3d.Manifold.cylinder(
                height=connectionLength,
                radius_low=outsetRadius + eps,
                radius_high=outsetRadius + eps,
                circular_segments=numSides)
            outset -= m3d.Manifold.cylinder(
                height=2*connectionLength,
                radius_low=self.r - 2*self.wallThickness,
                radius_high=self.r - 2*self.wallThickness,
                circular_segments=numSides)
            outset -= holeSlicer.translate((0,0,3*self.holeDiameter))
            outset = outset.translate((0,0,-connectionLength)) \
                           .rotate((0,90,0)).transform(openFrame.A[:3,:])
            shell += outset

        return shell

class PrintedStartHemisphere(PrintedHemisphere):
    def __init__(self, r : float, Pose : SE3, pathIndex : int = 2):
        super().__init__(r, Pose, closesForward=False, pathIndex=pathIndex)

class PrintedEndHemisphere(PrintedHemisphere):
    def __init__(self, r : float, Pose : SE3, pathIndex : int = 2):
        super().__init__(r, Pose, closesForward=True, pathIndex=pathIndex)


class PrintedLinkCSC(PrintedTube, LinkCSC):
    def __init__(self, r : float, StartDubinsPose : SE3, EndDubinsPose : SE3,
                 maxAnglePerElbow : float = np.pi/10, path : Optional[PathCSC] = None, 
                 EPSILON : float = 0.01, startRadius : Optional[float] = None, 
                 endRadius : Optional[float] = None):
        PrintedTube.__init__(self)
        LinkCSC.__init__(self, r, StartDubinsPose, EndDubinsPose,
                            maxAnglePerElbow=maxAnglePerElbow, path=path, EPSILON=EPSILON)
        self.startRadius = startRadius if startRadius is not None else r
        self.endRadius = endRadius if endRadius is not None else r
        
    def newLinkTransformedBy(self, Transformation : SE3):
        """Override to preserve PrintedLinkCSC type when transforming"""
        return PrintedLinkCSC(self.r, Transformation @ self.StartDubinsPose, 
                              Transformation @ self.EndDubinsPose,
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
                                         radius_low=self.holeDiameter/2 + self.looseFitTolerance, 
                                         radius_high=self.holeDiameter/2 + self.looseFitTolerance, 
                                         circular_segments=numSides)
            hole = hole.rotate((0,90,0)).rotate((0,0,angle))
            holeSlicer += hole
    
        inset = m3d.Manifold.cylinder(height=2*connectionLength+2*self.DISTANCE_EPSILON, 
                                       radius_low=self.startRadius-self.wallThickness+self.DISTANCE_EPSILON, 
                                       radius_high=self.startRadius-self.wallThickness+self.DISTANCE_EPSILON,
                                       circular_segments=numSides)
        tube -= inset.translate((0,0,-connectionLength-self.DISTANCE_EPSILON)).rotate((0,90,0)).transform(self.StartDubinsPose.A[:3,:])
        tube -= holeSlicer.translate((0,0,self.holeDiameter)).rotate((0,90,0)).transform(self.StartDubinsPose.A[:3,:])
        outsetRadius = self.endRadius-self.wallThickness
        outset = m3d.Manifold.cylinder(height=2*connectionLength, 
                                       radius_low=outsetRadius-self.tightFitTolerance, 
                                       radius_high=outsetRadius-self.tightFitTolerance,
                                       circular_segments=numSides)
        outset += m3d.Manifold.cylinder(height=connectionLength, 
                                       radius_low=outsetRadius+self.DISTANCE_EPSILON, 
                                       radius_high=outsetRadius+self.DISTANCE_EPSILON,
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
    looseFitTolerance = links[0].looseFitTolerance
    tightFitTolerance = links[0].tightFitTolerance
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
    holeSlicerLoose = m3d.Manifold()
    holeAnglesDegrees = np.linspace(0, 360, numHoles, endpoint=False)
    connectionLength = 2 * holeDiameter
    for angle in holeAnglesDegrees:
        hole = m3d.Manifold.cylinder(height=1.5*r, 
                                         radius_low=holeDiameter/2 + looseFitTolerance, 
                                         radius_high=holeDiameter/2 + looseFitTolerance, 
                                         circular_segments=numSides)
        hole = hole.rotate((0,90,0)).rotate((0,0,angle))
        holeSlicerLoose += hole
    
    holeSlicerTight = m3d.Manifold()
    for angle in holeAnglesDegrees:
        hole = m3d.Manifold.cylinder(height=1.5*r, 
                                         radius_low=holeDiameter/2 + tightFitTolerance, 
                                         radius_high=holeDiameter/2 + tightFitTolerance, 
                                         circular_segments=numSides)
        hole = hole.rotate((0,90,0)).rotate((0,0,angle))
        holeSlicerTight += hole
    
    inset = m3d.Manifold.cylinder(height=2*connectionLength+2*DISTANCE_EPSILON, 
                                       radius_low=startRadius-wallThickness+DISTANCE_EPSILON, 
                                       radius_high=startRadius-wallThickness+DISTANCE_EPSILON,
                                       circular_segments=numSides)
    tube -= inset.translate((0,0,-connectionLength-DISTANCE_EPSILON)).rotate((0,90,0)).transform(startPose.A[:3,:])
    tube -= holeSlicerLoose.translate((0,0,holeDiameter)).rotate((0,90,0)).transform(startPose.A[:3,:])

    # Create the outsets at each link end
    for link in links:
        outsetRadius = link.endRadius-wallThickness
        outset = m3d.Manifold.cylinder(height=2*connectionLength, 
                                       radius_low=outsetRadius-tightFitTolerance, 
                                       radius_high=outsetRadius-tightFitTolerance,
                                       circular_segments=numSides)
        outset += m3d.Manifold.cylinder(height=connectionLength, 
                                       radius_low=outsetRadius+DISTANCE_EPSILON, 
                                       radius_high=outsetRadius+DISTANCE_EPSILON,
                                       circular_segments=numSides)
        outset -= m3d.Manifold.cylinder(height=2*connectionLength, 
                                       radius_low=link.endRadius-2*wallThickness, 
                                       radius_high=link.endRadius-2*wallThickness,
                                       circular_segments=numSides)
        outset -= holeSlicerTight.translate((0,0,3*holeDiameter))
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
            return PrintedLinkCSC(r, start_pose, end_pose, max_angle_per_elbow, path, epsilon)
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

    def getCapModules(self, numSides : int = 50) -> dict[int, m3d.Manifold]:
        """Generate connectable cap modules for all PrintedHemisphere joints"""
        capModules = {}
        for jointIndex, joint in enumerate(self.Joints):
            if isinstance(joint, PrintedHemisphere):
                capModules[jointIndex] = joint.connectableModule(numSides=numSides)
        return capModules
    
    def saveLinkModules(self, baseFilename : str, numSides : int = 50, hullBends : bool = False,
                          maxSectionAngle : float = np.pi/10) -> None:
        """Save connectable branching modules and cap modules to files"""
        branchingModules = self.getLinkModules(numSides, hullBends, maxSectionAngle)
        for linkIndex, module in branchingModules.items():
            filename = f"{baseFilename}_link{linkIndex}.stl"
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            mesh_data = module.to_mesh()
            vertices = mesh_data.vert_properties[:, :3]  # Get XYZ coordinates
            faces = mesh_data.tri_verts
            tri_mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
            tri_mesh.export(filename)
        capModules = self.getCapModules(numSides)
        for jointIndex, module in capModules.items():
            filename = f"{baseFilename}_cap{jointIndex}.stl"
            os.makedirs(os.path.dirname(filename), exist_ok=True)
            mesh_data = module.to_mesh()
            vertices = mesh_data.vert_properties[:, :3]  # Get XYZ coordinates
            faces = mesh_data.tri_verts
            tri_mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
            tri_mesh.export(filename)
    
    def showLinkModules(self, numSides : int = 50, hullBends : bool = False,
                          maxSectionAngle : float = np.pi/10, block : bool = False) -> None:
        """Display connectable branching and cap modules in a 3D plot window"""
        branchingModules = self.getLinkModules(numSides, hullBends, maxSectionAngle)
        capModules = self.getCapModules(numSides)
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
        for jointIndex, module in capModules.items():
            mesh_data = module.to_mesh()
            vertices = mesh_data.vert_properties[:, :3]
            triangles = mesh_data.tri_verts

            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
            faces = [vertices[tri] for tri in triangles]
            mesh_collection = Poly3DCollection(faces, alpha=0.5,
                                              edgecolor='k',
                                              facecolors='orange')
            ax.add_collection3d(mesh_collection)
        ax.set_aspect('equal')
        plt.show(block=block)
