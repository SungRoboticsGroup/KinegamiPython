import pyqtgraph.opengl as gl
from spatialmath import SE3
import numpy as np
from style import *
from pyqtgraph import Transform3D
from OpenGL.GL import glDisable, glEnable, GL_DEPTH_TEST
import math
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

    def addArrows(self, widget, selectedArrow=-1, local=True, frame: SE3 = None, mode=""):
        # match GUI sizing
        desired_px = 80
        thickness_px = 10
        rad = widget.plot_widget.world_length_for_pixel_length(desired_px)
        tube_rad = widget.plot_widget.world_length_for_pixel_length(thickness_px)

        center = self.Pose.t
        colors = rotateArrowColors

        # build a cylinder or torus mesh with those dims
        if mode=="Translate":
            meshdata = gl.MeshData.cylinder(
                rows=2, cols=20,
                radius=[tube_rad,tube_rad], length=rad
            )
        else:  # Rotate
            meshdata = self.create_torus_mesh(
                radius=rad,
                tube_radius=tube_rad,
                radial_segments=20,
                tubular_segments=20
            )

        # add all three axes
        for i in range(3):
            item = gl.GLMeshItem(
                meshdata=meshdata,
                color=colors[i],
                shader='shaded',
                smooth=True
            )
            if mode=="Translate":
                if i==0: item.rotate(90, 0,1,0, True)
                if i==1: item.rotate(-90,1,0,0, True)
            item.translate(*center)
            widget.plot_widget.addItem(item)
        
    def addTranslateArrows(self, widget, selectedArrow=-1, local=True, frame: SE3=None):
        axes = (self.Pose.R if local else np.eye(3))
        if frame is not None:
            axes = frame.R @ axes

        center = self.Pose.t
        colors = rotateArrowColors

        for i in range(3):
            vec = axes[:, i]
            start = center
            end = center + self.mesh.scale() * vec
            line = OverlayLine(pos=np.array([start, end]),
                               color=colors[i],
                               width=8,
                               antialias=True)
            widget.plot_widget.addItem(line)

    def addRotateArrows(self, widget, selectedArrow=-1, local=True, frame: SE3=None):
        axes = (self.Pose.R if local else np.eye(3))
        if frame is not None:
            axes = frame.R @ axes

        center = self.Pose.t
        vis_rad = self.mesh.scale()
        colors = rotateArrowColors

        for i in range(3):
            pts = self.generate_circle_points(
                axis=axes[:, i],
                center=center,
                rad=vis_rad + 0.2,
                num_points=64
            )
            circle = OverlayLine(pos=np.array(pts),
                                 color=colors[i],
                                 width=8,
                                 antialias=True)
            widget.plot_widget.addItem(circle)

class OverlayLine(gl.GLLinePlotItem):
    def paint(self):
        glDisable(GL_DEPTH_TEST)
        super().paint()
        glEnable(GL_DEPTH_TEST)
