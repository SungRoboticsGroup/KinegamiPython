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
    
    def compute_intersection(self, org, dir, cen, rad):
        # does a ray-sphere intersection

        a = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]
        b = 2 * (dir[0] * (org[0] - cen[0]) + dir[1] * (org[1] - cen[1]) + dir[2] * (org[2] - cen[2]))
        c = (org[0] - cen[0]) ** 2 +  (org[1] - cen[1]) ** 2 + (org[2] - cen[2]) ** 2 - rad * rad

        discrim = b*b - 4*a*c

        if (discrim < 0) :
            print("miss")
        else: 
            print("hit")

        # t0 = (-b - math.sqrt(b*b - 4*a*c)) / (2*a)
        # t1 = (-b + math.sqrt(b*b - 4*a*c)) / (2*a)
    
    def mousePressEvent(self, event):
        self.selected_axis = (1, 0, 0)

        origin, direction = self.get_world_coordinates(event)
        self.start_pos_3D = self.compute_intersection(origin, direction, QVector3D(3,0,0), 1)
        event.setAccepted(True)

    def mouseMoveEvent(self, event):
        pass

    def mouseReleaseEvent(self, event):
        pass

    def keyPressEvent(self, event: QKeyEvent):
        if event.key() == Qt.Key_T:
            print("test")
        elif event.key() == Qt.Key_R:
            self.key_pressed.emit("Rotate")
        elif event.key() == Qt.Key_Return or event.key() == Qt.Key_Enter:
            self.key_pressed.emit("Enter")
        elif event.key() == Qt.Key_Escape:
            self.key_pressed.emit("Escape")
        elif event.key() == Qt.Key_Delete or event.key() == Qt.Key_Backspace:
            self.key_pressed.emit("Delete")
        elif event.key() == Qt.Key_X:
            self.key_pressed.emit("X")
        elif event.key() == Qt.Key_Y:
            self.key_pressed.emit("Y")
        elif event.key() == Qt.Key_Z:
            self.key_pressed.emit("Z")
        elif event.key() == Qt.Key_G:
            self.key_pressed.emit("G")
 
class PointEditorWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Point Editor")
        self.setGeometry(100, 100, 800, 600)

        self.plot_widget = ClickableGLViewWidget()
        self.setCentralWidget(self.plot_widget)
        self.plot_widget.setBackgroundColor(backgroundColorDefault)

        self.grid = gl.GLGridItem()

        self.plot_widget.addItem(self.grid)
        self.grid.setColor(gridColorDefault)

        md = gl.MeshData.sphere(rows=20, cols=20)
        sphere = gl.GLMeshItem(meshdata=md, color=tuple((0, 0, 0, 0.5)), shader='shaded', smooth=True, id=id)

        transform = QMatrix4x4()
        transform.translate(3, 0, 0)

        # Apply the transformation to the sphere
        sphere.setTransform(transform)

        self.plot_widget.addItem(sphere)
    
    @QtCore.pyqtSlot(str)
    def key_pressed(self, key):
        if key == "Translate":
            self.control_type = key
            print("test")
        elif key == "Rotate":
            self.control_type = key
            self.rotateJointRadioButton.setChecked(True)
        elif key == "Delete":
            if self.chain and self.selected_joint != -1:
                self.delete_selected_joint()
        elif key == "X":
            self.arrow_selection_changed(0)
        elif key == "Y":
            self.arrow_selection_changed(1)
        elif key == "Z":
            self.arrow_selection_changed(2)
        elif key == "G":
            self.toggle_grid_func()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PointEditorWindow()
    window.show()
    sys.exit(app.exec_())