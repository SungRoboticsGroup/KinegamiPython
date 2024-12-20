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
        self.mesh = None 
        self.is_dragging = False
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
        # not bothering with computing endcaps because they're too small to matter anyway

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

        closest = 1000
        closest_sphere = None

        for sphere in self.parent_window.spheres: 
            trans = sphere.transform().column(3)
            center = QVector3D(trans[0], trans[1], trans[2])
            hit_location = self.compute_sphere_intersection(origin, direction, center, self.parent_window.radius)
            if (hit_location < closest):
                closest = hit_location
                closest_sphere = sphere

        hit_location = self.compute_cylinder_intersection(origin, direction, QVector3D(0,0,0), QVector3D(0,-1,0), 0.1, 4)
        if (hit_location < 1000):
            print("hit cylinder")

        if closest_sphere:
            print(closest_sphere.objectName())

        event.setAccepted(True)

    def mouseReleaseEvent(self, event):
        pass
 
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
        cylinder_md = gl.MeshData.cylinder(rows=2, cols=20, radius=[0.1, 0.1], length=1)

        sphere_red = gl.GLMeshItem(meshdata=sphere_md, color=tuple((1, 0, 0, 0.5)), shader='shaded', glOptions='translucent', smooth=True)
        sphere_red.translate(3, 0, 0)
        sphere_red.setObjectName("Red Sphere")
        self.plot_widget.addItem(sphere_red)

        sphere_green = gl.GLMeshItem(meshdata=sphere_md, color=tuple((0, 1, 0, 0.5)), shader='shaded', glOptions='translucent', smooth=True)
        sphere_green.translate(0, 3, 0)
        sphere_green.setObjectName("Green Sphere")
        self.plot_widget.addItem(sphere_green)

        self.spheres = [sphere_red, sphere_green]

        axis_x = gl.GLMeshItem(meshdata=cylinder_md, color=tuple((1, 0, 0, 1)), shader='shaded', smooth=True)
        axis_x.rotate(90, 0, 1, 0)
        axis_x.translate(1, 3, 0)
        axis_x.setObjectName("X axis")
        # self.plot_widget.addItem(axis_x)

        axis_y = gl.GLMeshItem(meshdata=cylinder_md, color=tuple((0, 1, 0, 1)), shader='shaded', smooth=True)
        axis_y.rotate(90, 1, 0, 0)
        axis_y.translate(0, 5, 0)
        axis_y.setObjectName("Y axis")
        # self.plot_widget.addItem(axis_y)

        axis_z = gl.GLMeshItem(meshdata=cylinder_md, color=tuple((0, 0, 1, 1)), shader='shaded', smooth=True)
        axis_z.translate(0, 3, 1)
        axis_z.setObjectName("Z axis")
        # self.plot_widget.addItem(axis_z)

        cylinder_test = gl.MeshData.cylinder(rows=2, cols=20, radius=[0.1, 0.1], length=4)
        cylinder = gl.GLMeshItem(meshdata=cylinder_test, color=tuple((0, 0, 1, 1)), shader='shaded', smooth=True)
        cylinder.rotate(90, 1, 0, 0)
        cylinder.setObjectName("Cylinder")
        self.plot_widget.addItem(cylinder)

        self.radius = 1

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PointEditorWindow()
    window.show()
    sys.exit(app.exec_())