import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QPixmap, QSurfaceFormat, QKeyEvent, QVector3D, QMatrix4x4, QVector4D
import sys
from style import *
from PyQt5 import QtCore as qc
from pyqtgraph.Qt import QtCore
import math
import numpy as np

import warnings
warnings.filterwarnings("ignore")

def create_torus_mesh(radius, tube_radius, radial_segments, tubular_segments):
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

def create_cube_mesh(size):
    d = size / 2.0
    vertices = np.array([
        [-d, -d, -d],
        [ d, -d, -d],
        [ d,  d, -d],
        [-d,  d, -d],
        [-d, -d,  d],
        [ d, -d,  d],
        [ d,  d,  d],
        [-d,  d,  d]
    ])

    faces = np.array([
        [0, 1, 2], [0, 2, 3],
        [4, 5, 6], [4, 6, 7],
        [0, 1, 5], [0, 5, 4],
        [2, 3, 7], [2, 7, 6],
        [0, 3, 7], [0, 7, 4],
        [1, 2, 6], [1, 6, 5]
    ])

    return gl.MeshData(vertexes=vertices, faces=faces)

class ClickableGLViewWidget(gl.GLViewWidget):
    def __init__(self, parent=None):
        super(ClickableGLViewWidget, self).__init__(parent)
        fmt = QSurfaceFormat()
        fmt.setDepthBufferSize(24)
        # fmt.setVersion(3, 3)
        fmt.setProfile(QSurfaceFormat.CompatibilityProfile)
        QSurfaceFormat.setDefaultFormat(fmt)
        self.setFormat(fmt)
        self.locked = False
        self.axis = None 
        self.is_dragging = False
        self.cube_start_pos = None
        self.drag_prev_vector = 0

        self.bounding_balls = []
        self.radius = 1

        self.axes = {
            'x': QVector3D(1, 0, 0),
            'y': QVector3D(0, 1, 0),
            'z': QVector3D(0, 0, 1)
        }

        dist = self.opts['distance']
        self.near_clip = dist * 0.001
        self.far_clip = dist * 1000.
        self.parent_window = parent

    key_pressed = qc.pyqtSignal(str)
    
    def get_world_coordinates(self, event):
        pos = event.localPos()
        ndc_x = (2.0 * pos.x()) / self.width() - 1.0
        ndc_y = 1.0 - (2.0 * pos.y()) / self.height()

        ndc = QVector3D(ndc_x, ndc_y, -1.0)
        view = self.viewMatrix()
        proj = self.projectionMatrix()

        inverted_matrix = (proj * view).inverted()[0]

        near_point = inverted_matrix.map(ndc)
        ndc.setZ(1.0)
        far_point = inverted_matrix.map(ndc)

        direction = far_point - near_point
        direction.normalize()

        return near_point, direction

    def compute_plane_intersection(self, org, dir, normal, point):
        threshold = 0.001
        direction = dir.normalized()
        n = normal.normalized()
        a = QVector3D.dotProduct(n, direction)

        if (abs(a) < threshold):
            return None
        
        b = point - org
        t = QVector3D.dotProduct(n, b) / a

        return org + t * direction
    
    def get_normalized_plane_vectors(self, event):
        if self.is_dragging and self.axis:
            origin, dir = self.get_world_coordinates(event)
            center = self.cube_start_pos
            qcenter = QVector3D(center[0], center[1], center[2])
            qaxis = self.axes[self.axis]
            npos = self.compute_plane_intersection(origin, dir, qaxis, qcenter)

            plane_vector = npos - qcenter
            plane_vector.normalize()

            normal = QVector3D.dotProduct(dir, qaxis) * qaxis
            normal.normalize()

            return plane_vector, normal
        

    def get_axis_angle_delta(self, event):
        if self.is_dragging and self.axis:
            origin, dir = self.get_world_coordinates(event)
            center = self.cube_start_pos
            qcenter = QVector3D(center[0], center[1], center[2])
            qaxis = self.axes[self.axis]
            
            plane_vector, normal = self.get_normalized_plane_vectors(event)

            prev_vector = self.drag_prev_vector

            d = QVector3D.dotProduct(plane_vector, prev_vector)
            angle = math.acos(d / (plane_vector.length() * prev_vector.length()))

            cross_product = QVector3D.crossProduct(prev_vector, plane_vector)
            
            if QVector3D.dotProduct(cross_product, normal) < 0:
                angle *= -1

            self.drag_prev_vector = plane_vector

            return math.degrees(angle), normal
    
    def get_zero_axis(self, axis: str):
        match axis:
            case 'x':
                return self.axes['x']
            case 'y':
                return self.axes['y']
            case 'z':
                return self.axes['z']
            case _:
                raise ValueError("Invalid axis")
    
    def compute_axis_intersection(self, org, dir, center, axis):
        threshold = 0.001
        # Normalize direction vectors
        direction1 = dir.normalized()
        direction2 = axis.normalized()
        
        # Compute the cross product of the direction vectors
        cross_directions = QVector3D.crossProduct(direction1, direction2)
        cross_directions_norm = cross_directions.length()
        
        # If the cross product is zero, the lines are parallel
        if cross_directions_norm < threshold:
            raise ValueError("The lines are parallel and do not intersect.")
        
        # Compute the vector between the origins
        origin_diff = center - org
        
        # Compute the determinants
        det1 = QVector3D.dotProduct(origin_diff, QVector3D.crossProduct(direction2, cross_directions))
        det2 = QVector3D.dotProduct(origin_diff, QVector3D.crossProduct(direction1, cross_directions))
        
        # Compute the parameters for the points of closest approach
        t1 = det1 / cross_directions_norm**2
        t2 = det2 / cross_directions_norm**2
        
        # Compute the points of closest approach
        point1 = org + t1 * direction1
        point2 = center + t2 * direction2
        
        return point2
    
    def get_axis_intersection(self, event):
        if self.is_dragging and self.axis:
            origin, dir = self.get_world_coordinates(event)
            center = self.cube_start_pos
            qcenter = QVector3D(center[0], center[1], center[2])
            axis = self.axis
            new_pos_3D = self.compute_axis_intersection(origin, dir, qcenter, axis)

            return new_pos_3D
    
    def compute_sphere_intersection(self, org, dir, cen, rad):
        a = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]
        b = 2 * (dir[0] * (org[0] - cen[0]) + dir[1] * (org[1] - cen[1]) + dir[2] * (org[2] - cen[2]))
        c = (org[0] - cen[0]) ** 2 +  (org[1] - cen[1]) ** 2 + (org[2] - cen[2]) ** 2 - rad * rad

        discrim = b*b - 4*a*c

        if (discrim < -0.00001) :
            return 2000
        else: 
            t0 = (-b - math.sqrt(discrim)) / (2*a)

            if (t0 > 0):
                return t0
            else:
                return (-b + math.sqrt(discrim)) / (2*a)
            
    def compute_cylinder_intersection(self, org: QVector3D, dir: QVector3D, start: QVector3D, axis: QVector3D, rad, len): 
        n = dir.normalized()
        a = axis.normalized()
        b = start - org

        n_cross_a = QVector3D.crossProduct(n, a)

        discrim = QVector3D.dotProduct(n_cross_a, n_cross_a) * rad * rad - QVector3D.dotProduct(a, a) * (QVector3D.dotProduct(b, n_cross_a) ** 2)

        if (discrim < -0.00001):
            return 2000
        else:
            d = (QVector3D.dotProduct(n_cross_a, QVector3D.crossProduct(b, a)) - math.sqrt(discrim)) / QVector3D.dotProduct(n_cross_a, n_cross_a)
            
            if (d < 0):
                d = (QVector3D.dotProduct(n_cross_a, QVector3D.crossProduct(b, a)) + math.sqrt(discrim)) / QVector3D.dotProduct(n_cross_a, n_cross_a)

            t = QVector3D.dotProduct(a, (n * d - b))
            if (t > 0 and t < len):
                return d
            else:
                return 2000

    def compute_torus_intersection(self, org: QVector3D, dir: QVector3D, center: QVector3D, normal: QVector3D, major_radius: float, minor_radius: float):
        n = dir.normalized()
        o = org - center
        n_dot_n = QVector3D.dotProduct(n, n)
        o_dot_n = QVector3D.dotProduct(o, n)
        o_dot_o = QVector3D.dotProduct(o, o)
        n_dot_o = QVector3D.dotProduct(n, o)
        n_dot_normal = QVector3D.dotProduct(n, normal)
        o_dot_normal = QVector3D.dotProduct(o, normal)
        
        R = major_radius
        r = minor_radius
        
        A = n_dot_n * n_dot_n
        B = 4 * n_dot_n * n_dot_o
        C = 2 * n_dot_n * (o_dot_o - R * R - r * r) + 4 * n_dot_o * n_dot_o + 4 * R * R * n_dot_normal * n_dot_normal
        D = 4 * (o_dot_o - R * R - r * r) * n_dot_o + 8 * R * R * o_dot_normal * n_dot_normal
        E = (o_dot_o - R * R - r * r) * (o_dot_o - R * R - r * r) - 4 * R * R * (r * r - o_dot_normal * o_dot_normal)
        
        coeffs = [A, B, C, D, E]
        roots = np.roots(coeffs)
        
        real_roots = [root.real for root in roots if np.isreal(root) and root.real > 0]
        
        if not real_roots:
            return 2000
        
        return min(real_roots)
    
    def mousePressEvent(self, event):
        origin, direction = self.get_world_coordinates(event)

        self.cube_start_pos = self.parent_window.cube_center.copy()

        closest = 1000
        closest_cylinder = None
        closest_axis = None

        for i in range(3):
            cylinder = self.parent_window.axes[i]
            trans = self.cube_start_pos.copy()

            if i == 0:
                axis = 'x'
            elif i == 1:
                axis = 'y'
            else:
                axis = 'z'

            center = QVector3D(trans[0], trans[1], trans[2])

            hit_location = self.compute_torus_intersection(origin, direction, center, self.axes[axis], major_radius=self.parent_window.tor_rad, minor_radius=self.parent_window.tube_rad+0.1)
            if (hit_location < closest):
                closest = hit_location
                closest_cylinder = cylinder
                closest_axis = axis

        if (hit_location < 1000):
            print("hit cylinder")

        if closest_cylinder:
            print(closest_cylinder.objectName())
            self.axis = closest_axis
            self.is_dragging = True

            self.drag_prev_vector, _ = self.get_normalized_plane_vectors(event)

            self.parent_window.update_visibility(self.axis)
            # self.parent_window.draw_axis_line(self.cube_start_pos, self.axis)

        event.setAccepted(True)

    def mouseMoveEvent(self, event):
        if self.is_dragging and self.axis:
            da, normal = self.get_axis_angle_delta(event)

            self.parent_window.update_mesh(da, normal)
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        self.axis = None
        self.is_dragging = False
        self.cube_start_pos = None
        self.drag_prev_vector = None

        # set axes
        obj = self.parent_window.cube
        cnt = self.parent_window.cube_center
        obj.translate(-cnt[0], -cnt[1], -cnt[2])
        rot = obj.viewTransform()
        obj.translate(cnt[0], cnt[1], cnt[2])

        self.axes['x'] = rot * QVector3D(1, 0, 0)
        self.axes['y'] = rot * QVector3D(0, 1, 0)
        self.axes['z'] = rot * QVector3D(0, 0, 1)

        print(self.axes)

        self.parent_window.update_visibility(None)
        
 
class PointEditorWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Point Editor")
        self.setGeometry(100, 100, 800, 600)

        self.plot_widget = ClickableGLViewWidget(self)
        self.setCentralWidget(self.plot_widget)
        self.plot_widget.setBackgroundColor(backgroundColorDefault)

        self.grid = gl.GLGridItem()
        self.plot_widget.addItem(self.grid)
        self.grid.setColor(gridColorDefault)

        sphere_md = gl.MeshData.sphere(rows=20, cols=20)
        self.c_rad = 0.5

        self.tube_rad=0.1
        self.tor_rad=1
        torus_md = create_torus_mesh(radius=self.tor_rad, tube_radius=self.tube_rad, radial_segments=20, tubular_segments=20)
        cylinder_md = gl.MeshData.cylinder(rows=2, cols=20, radius=[self.c_rad, self.c_rad], length=1)
        cube_md = create_cube_mesh(1)

        sphere_red = gl.GLMeshItem(meshdata=sphere_md, color=tuple((1, 0, 0, 0.5)), shader='shaded', glOptions='translucent', smooth=True)
        sphere_red.translate(3, 0, 0)
        sphere_red.setObjectName("Red Sphere")
        self.plot_widget.addItem(sphere_red)

        cube_grey = gl.GLMeshItem(meshdata=cube_md, color=tuple((1, 1, 1, 1)), shader='shaded', glOptions='translucent', smooth=False)
        cube_grey.translate(0, 3, 0)
        cube_grey.setObjectName("Grey Cube")
        self.plot_widget.addItem(cube_grey)

        self.spheres = [sphere_red]

        self.cube = cube_grey

        self.cube_center = [0,3,0]

        axis_x = gl.GLMeshItem(meshdata=torus_md, color=tuple((1, 0, 0, .5)), shader='shaded', glOptions='translucent', smooth=True)
        x_trans = self.cube_center.copy()
        axis_x.translate(x_trans[0], x_trans[1], x_trans[2])
        axis_x.rotate(90, 0, 1, 0, True)
        axis_x.setObjectName("X axis")
        self.plot_widget.addItem(axis_x)

        axis_y = gl.GLMeshItem(meshdata=torus_md, color=tuple((0, 1, 0, .5)), shader='shaded', glOptions='translucent', smooth=True)
        y_trans = self.cube_center.copy()
        axis_y.translate(y_trans[0], y_trans[1], y_trans[2])
        axis_y.rotate(90, 1, 0, 0, True)
        axis_y.setObjectName("Y axis")
        self.plot_widget.addItem(axis_y)

        axis_z = gl.GLMeshItem(meshdata=torus_md, color=tuple((0, 0, 1, .5)), shader='shaded', glOptions='translucent', smooth=True)
        z_trans = self.cube_center.copy()
        axis_z.translate(z_trans[0], z_trans[1], z_trans[2])
        axis_z.setObjectName("Z axis")
        self.plot_widget.addItem(axis_z)

        self.axes = [axis_x, axis_y, axis_z]

        self.axis_line = None

        # cylinder_test = gl.MeshData.cylinder(rows=2, cols=20, radius=[0.1, 0.1], length=4)
        # cylinder = gl.GLMeshItem(meshdata=cylinder_test, color=tuple((0, 0, 1, 1)), shader='shaded', smooth=True)
        # cylinder.rotate(90, 1, 0, 0)
        # cylinder.setObjectName("Cylinder")
        # self.plot_widget.addItem(cylinder)

        self.radius = 1

    def update_mesh(self, angle, axis):
        # Update the position of the spheres
        obj = self.cube
        cnt = self.cube_center


        obj.translate(-cnt[0], -cnt[1], -cnt[2])
        obj.rotate(angle, axis[0], axis[1], axis[2], local=False)
        obj.translate(cnt[0], cnt[1], cnt[2])

        for i, a in enumerate(self.axes):
            a.translate(-cnt[0], -cnt[1], -cnt[2])
            a.rotate(angle, axis[0], axis[1], axis[2], local=False)
            a.translate(cnt[0], cnt[1], cnt[2])

    def update_visibility(self, axis):
        if axis == 'x':
            c = self.axes[0].opts['color']
            new_color = (c[0], c[1], c[2], 1)
            self.axes[0].setColor(new_color)

            self.axes[1].setVisible(False)
            self.axes[2].setVisible(False)
        elif axis == 'y':
            c = self.axes[1].opts['color']
            new_color = (c[0], c[1], c[2], 1)
            self.axes[1].setColor(new_color)

            self.axes[0].setVisible(False)
            self.axes[2].setVisible(False)
        elif axis == 'z':
            c = self.axes[2].opts['color']
            new_color = (c[0], c[1], c[2], 1)
            self.axes[2].setColor(new_color)

            self.axes[0].setVisible(False)
            self.axes[1].setVisible(False)
        else:
            for a in self.axes:
                c = a.opts['color']
                new_color = (c[0], c[1], c[2], 0.5)
                a.setColor(new_color)
                a.setVisible(True)

    def draw_point(self, point):
        if not point:
            return
        sphere_md = gl.MeshData.sphere(rows=5, cols=5, radius=.1)
        sphere = gl.GLMeshItem(meshdata=sphere_md, color=tuple((0.1, 0.1, 0.1, 0.5)), shader='shaded', glOptions='translucent', smooth=True)
        sphere.translate(point[0], point[1], point[2])
        self.plot_widget.addItem(sphere)

    def draw_axis_line(self, origin, axis):
        if self.axis_line:
            self.plot_widget.removeItem(self.axis_line)
            self.axis_line = None
        if axis:
            length = 1000
            qstart_pt = origin + axis * length
            qend_pt = origin - axis * length
            start_pt = [qstart_pt[0], qstart_pt[1], qstart_pt[2]]
            end_pt = [qend_pt[0], qend_pt[1], qend_pt[2]]

            axis_line = np.array([start_pt, end_pt])
            line = gl.GLLinePlotItem(pos=axis_line, color=(axis[0], axis[1], axis[2], 1), width=3, antialias=True)
            self.axis_line = line
            self.plot_widget.addItem(line)
        

        
        
        


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PointEditorWindow()
    window.show()
    sys.exit(app.exec_())