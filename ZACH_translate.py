import pyqtgraph as pg
import pyqtgraph.opengl as gl
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QPixmap, QSurfaceFormat, QKeyEvent, QVector3D, QMatrix4x4
import sys
from style import *
from PyQt5 import QtCore as qc
from pyqtgraph.Qt import QtCore
import math
import numpy as np

import warnings
warnings.filterwarnings("ignore")

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
        self.sphere_start_pos = None
        self.drag_start_pos = None

        self.bounding_balls = []
        self.radius = 1

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
            center = self.sphere_start_pos
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
    
    def mousePressEvent(self, event):
        origin, direction = self.get_world_coordinates(event)

        self.sphere_start_pos = self.parent_window.sphere_center.copy()

        closest = 1000
        closest_cylinder = None
        closest_axis = None

        # for sphere in self.parent_window.spheres: 
        #     trans = sphere.transform().column(3)
        #     center = QVector3D(trans[0], trans[1], trans[2])
        #     hit_location = self.compute_sphere_intersection(origin, direction, center, self.parent_window.radius)
        #     if (hit_location < closest):
        #         closest = hit_location
        #         closest_sphere = sphere

        for i in range(3):
            cylinder = self.parent_window.axes[i]
            trans = self.sphere_start_pos.copy()
            rad = self.parent_window.radius

            axis = [0,0,0]
            axis[i] = 1
            qaxis = QVector3D(axis[0], axis[1], axis[2])
            start = QVector3D(trans[0] - qaxis[0]*rad, trans[1] - qaxis[1]*rad, trans[2] - qaxis[2]*rad)
            hit_location = self.compute_cylinder_intersection(origin, direction, start, qaxis, self.parent_window.c_rad, 4)
            if (hit_location < closest):
                closest = hit_location
                closest_cylinder = cylinder
                closest_axis = qaxis

        if (hit_location < 1000):
            print("hit cylinder")

        if closest_cylinder:
            print(closest_cylinder.objectName())
            self.axis = closest_axis
            self.is_dragging = True

            self.drag_start_pos = self.get_axis_intersection(event)
            self.parent_window.draw_axis_line(self.drag_start_pos, self.axis)

        event.setAccepted(True)

    def mouseMoveEvent(self, event):
        if self.is_dragging and self.axis:
            new_pos_3D = self.get_axis_intersection(event)
            qsphere_start = QVector3D(self.sphere_start_pos[0], self.sphere_start_pos[1], self.sphere_start_pos[2])
            new_pos_3D = qsphere_start + new_pos_3D - self.drag_start_pos
            # print(new_pos_3D)

            # draw point debug
            # self.parent_window.draw_point(new_pos_3D)

            # move the sphere
            self.parent_window.sphere_center = [new_pos_3D.x(), new_pos_3D.y(), new_pos_3D.z()]
            self.parent_window.update_mesh()
        else:
            super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        self.axis = None
        self.is_dragging = False
        self.sphere_start_pos = None
        self.parent_window.draw_axis_line(None, None)
        
 
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
        self.c_rad = 0.1
        cylinder_md = gl.MeshData.cylinder(rows=2, cols=20, radius=[self.c_rad, self.c_rad], length=1)

        sphere_red = gl.GLMeshItem(meshdata=sphere_md, color=tuple((1, 0, 0, 0.5)), shader='shaded', glOptions='translucent', smooth=True)
        sphere_red.translate(3, 0, 0)
        sphere_red.setObjectName("Red Sphere")
        self.plot_widget.addItem(sphere_red)

        sphere_green = gl.GLMeshItem(meshdata=sphere_md, color=tuple((0, 1, 0, 0.5)), shader='shaded', glOptions='translucent', smooth=True)
        sphere_green.translate(0, 3, 0)
        sphere_green.setObjectName("Green Sphere")
        self.plot_widget.addItem(sphere_green)

        self.spheres = [sphere_red, sphere_green]

        self.sphere_center = [0,3,0]

        axis_x = gl.GLMeshItem(meshdata=cylinder_md, color=tuple((1, 0, 0, 1)), shader='shaded', smooth=True)
        x_trans = self.sphere_center.copy()
        x_trans[0] += 1
        axis_x.translate(x_trans[0], x_trans[1], x_trans[2])
        axis_x.rotate(90, 0, 1, 0, True)
        axis_x.setObjectName("X axis")
        self.plot_widget.addItem(axis_x)

        axis_y = gl.GLMeshItem(meshdata=cylinder_md, color=tuple((0, 1, 0, 1)), shader='shaded', smooth=True)
        y_trans = self.sphere_center.copy()
        y_trans[1] += 2
        axis_y.translate(y_trans[0], y_trans[1], y_trans[2])
        axis_y.rotate(90, 1, 0, 0, True)
        axis_y.setObjectName("Y axis")
        self.plot_widget.addItem(axis_y)

        axis_z = gl.GLMeshItem(meshdata=cylinder_md, color=tuple((0, 0, 1, 1)), shader='shaded', smooth=True)
        z_trans = self.sphere_center.copy()
        z_trans[2] += 1
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

    def update_mesh(self):
        # Update the position of the spheres
        sphere = self.spheres[1]
        sphere.resetTransform()
        sphere.translate(self.sphere_center[0], self.sphere_center[1], self.sphere_center[2])

        # Update the position of the axes
        for i, axis in enumerate(self.axes):
            axis.resetTransform()
            trans = self.sphere_center.copy()
            if i == 0:
                trans[0] += 1
                axis.translate(trans[0], trans[1], trans[2])
                axis.rotate(90, 0, 1, 0, True)
            elif i == 1:
                trans[1] += 2
                axis.translate(trans[0], trans[1], trans[2])
                axis.rotate(90, 1, 0, 0, True)
            elif i == 2:
                trans[2] += 1
                axis.translate(trans[0], trans[1], trans[2])

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