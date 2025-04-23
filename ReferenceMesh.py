import pyqtgraph.opengl as gl
from spatialmath import SE3
import numpy as np
from style import *
from pyqtgraph import Transform3D
from PyQt5.QtGui import QVector4D 

class ReferenceMesh():

    def __init__(self, mesh : gl.GLMeshItem):
        self.mesh = mesh
        self.Pose = SE3()
        self.r = 10.0

    def updateScale(self, scale):
        prev = self.mesh.transform().matrix()
        prev[0][0] = scale
        prev[1][1] = scale
        prev[2][2] = scale
        new_mat = Transform3D(prev)
        self.mesh.setTransform(new_mat)

    def get_transform3D(self, translate=False):
        m = self.Pose.R
        t = self.Pose.t

        trans = [0, 0, 0]
        if (translate):
            trans = t

        transform = Transform3D()
        transform.setRow(0, QVector4D(m[0][0], m[0][1], m[0][2], trans[0]))
        transform.setRow(1, QVector4D(m[1][0], m[1][1], m[1][2], trans[1]))
        transform.setRow(2, QVector4D(m[2][0], m[2][1], m[2][2], trans[2]))
        transform.setRow(3, QVector4D(0, 0, 0, 1))

        return transform
    
    def create_torus_mesh(self, radius, tube_radius, radial_segments, tubular_segments):
        theta = np.linspace(0, 2 * np.pi, radial_segments)
        phi = np.linspace(0, 2 * np.pi, tubular_segments)
        theta, phi = np.meshgrid(theta, phi)
        theta, phi = theta.flatten(), phi.flatten()

        x = (radius + tube_radius * np.cos(phi)) * np.cos(theta)
        y = (radius + tube_radius * np.cos(phi)) * np.sin(theta)
        z = tube_radius * np.sin(phi)

        vertices = np.vstack([x, y, z]).T
        faces = []

        for i in range(radial_segments):
            for j in range(tubular_segments):
                next_i = (i + 1) % radial_segments
                next_j = (j + 1) % tubular_segments

                faces.append([i * tubular_segments + j,
                            next_i * tubular_segments + j,
                            i * tubular_segments + next_j])
                faces.append([next_i * tubular_segments + j,
                            next_i * tubular_segments + next_j,
                            i * tubular_segments + next_j])

        faces = np.array(faces)
        return gl.MeshData(vertexes=vertices, faces=faces)
    
    def generate_extended_axis(self, point1, point2, length):
        direction_vector = point2 - point1
    
        midpoint = (point1 + point2) / 2

        p1 = midpoint - direction_vector * length
        p2 = midpoint + direction_vector * length

        return np.array([p1, p2])
    
    def transform(self, transform: SE3):
        self.Pose = transform
        self.mesh.setTransform(self.get_transform3D())
    
    def addArrows(self, widget, selectedArrow=-1, local=True, frame : SE3 = None, mode=""):
        rad = self.r
        colors = rotateArrowColors
        opacity = [0.8, 0.8, 0.8]

        # if selectedArrow != -1:
        #     colors[selectedArrow] = selectedArrowColor
        #     opacity = [0.1, 0.1, 0.1]
        #     opacity[selectedArrow] = 0.8

        extended_axis_color = [(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)]
        center = self.Pose.t

        if local:
            axes = [self.Pose.R[:, i] for i in range(3)]
        else:
            axes = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]

        if frame:
            axes = [frame.R[:, i] for i in range(3)]

        transform = self.get_transform3D()

        if (mode == "Translate"):
            meshdata = gl.MeshData.cylinder(rows=2, cols=20, radius=[0.1,0.1], length=rad+1)
        elif (mode == "Rotate"):
            self.tube_rad=0.1
            self.tor_rad=rad + 0.2
            meshdata = self.create_torus_mesh(radius=self.tor_rad, tube_radius=self.tube_rad, radial_segments=20, tubular_segments=20)

        axis_x = gl.GLMeshItem(meshdata=meshdata, color=colors[0], shader='shaded', smooth=True)
        axis_x.rotate(90, 0, 1, 0, True)
        #widget.plot_widget.addItem(axis_x)

        axis_y = gl.GLMeshItem(meshdata=meshdata, color=colors[1], shader='shaded', smooth=True)
        axis_y.rotate(-90, 1, 0, 0, True)
        #widget.plot_widget.addItem(axis_y)

        axis_z = gl.GLMeshItem(meshdata=meshdata, color=colors[2], shader='shaded', smooth=True)
        #widget.plot_widget.addItem(axis_z)

        transform = self.get_transform3D(False)

        if (local):
            axis_x.applyTransform(transform, False)
            axis_y.applyTransform(transform, False)
            axis_z.applyTransform(transform, False)

        axis_x.translate(center[0], center[1], center[2])
        axis_y.translate(center[0], center[1], center[2])
        axis_z.translate(center[0], center[1], center[2])

        if selectedArrow != -1:
            # generate the line here
            dir = center + rad * self.r * axes[selectedArrow]
            extended_axis = self.generate_extended_axis(center, dir, 1)
            extended_axis_line = gl.GLLinePlotItem(pos=extended_axis, color=extended_axis_color[selectedArrow], width=3, antialias=True)
            widget.plot_widget.addItem(extended_axis_line)

        widget.plot_widget.addItem(axis_x)
        widget.plot_widget.addItem(axis_y)
        widget.plot_widget.addItem(axis_z)
        
    def addTranslateArrows(self, widget, selectedArrow=-1, local=True, frame : SE3 = None):
         self.addArrows(widget, selectedArrow, local, frame, mode="Translate")

    def addRotateArrows(self, widget, selectedArrow=-1, local=True, frame : SE3 = None):
        self.addArrows(widget, selectedArrow, local, frame, mode="Rotate")