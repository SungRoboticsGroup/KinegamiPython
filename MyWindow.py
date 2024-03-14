"""
@author: Raymond Feng
"""

import sys
import numpy as np
import pyqtgraph.opengl as gl
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton, QDockWidget, QComboBox, QHBoxLayout, QLabel
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap, QPainter
from spatialmath import SE3
import math
from PathCSC import *

class PointEditorWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Point Editor")
        self.setGeometry(100, 100, 800, 600)

        self.plot_widget = gl.GLViewWidget()
        self.setCentralWidget(self.plot_widget)
        self.plot_widget.setBackgroundColor(255,255,255, 255)

        grid = gl.GLGridItem()
        self.plot_widget.addItem(grid)
        grid.setColor((0,0,0,255))

        self.current_point = 0

        self.points = []
        self.poses = []

        self.add_point_button = QPushButton("Add Point")
        self.add_point_button.clicked.connect(self.add_point)
        self.select_point_options = QComboBox()
        self.select_point_options.activated.connect(self.index_changed)
        self.select_point_options.currentIndexChanged.connect(self.index_changed)

        self.move_x_pos = QPushButton("+X")
        self.move_x_pos.clicked.connect(self.move_point_x_pos)
        self.move_x_neg = QPushButton("-X")
        self.move_x_neg.clicked.connect(self.move_point_x_neg)

        self.move_y_pos = QPushButton("+Y")
        self.move_y_pos.clicked.connect(self.move_point_y_pos)
        self.move_y_neg = QPushButton("-Y")
        self.move_y_neg.clicked.connect(self.move_point_y_neg)

        self.move_z_pos = QPushButton("+Z")
        self.move_z_pos.clicked.connect(self.move_point_z_pos)
        self.move_z_neg = QPushButton("-Z")
        self.move_z_neg.clicked.connect(self.move_point_z_neg)

        self.rotate_button_x = QPushButton("Rotate X")
        self.rotate_button_x.clicked.connect(self.rotate_point_x)
        self.rotate_button_y = QPushButton("Rotate Y")
        self.rotate_button_y.clicked.connect(self.rotate_point_y)
        self.rotate_button_z = QPushButton("Rotate Z")
        self.rotate_button_z.clicked.connect(self.rotate_point_z)

        button_layout = QVBoxLayout()
        button_layout.addWidget(self.select_point_options)
        button_layout.addWidget(self.add_point_button)

        x_layout = QHBoxLayout()
        x_layout.addWidget(self.move_x_neg)
        x_layout.addWidget(self.move_x_pos) 

        y_layout = QHBoxLayout()
        y_layout.addWidget(self.move_y_neg)
        y_layout.addWidget(self.move_y_pos) 

        z_layout = QHBoxLayout()
        z_layout.addWidget(self.move_z_neg)
        z_layout.addWidget(self.move_z_pos) 

        rotate_layout = QHBoxLayout()
        rotate_layout.addWidget(self.rotate_button_x)
        rotate_layout.addWidget(self.rotate_button_y)
        rotate_layout.addWidget(self.rotate_button_z)

        button_layout.addLayout(x_layout)
        button_layout.addLayout(y_layout)
        button_layout.addLayout(z_layout)
        button_layout.addLayout(rotate_layout)

        button_widget = QWidget()
        button_widget.setLayout(button_layout)
        dock = QDockWidget("title")
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        dock.setWidget(button_widget)

        axis_key_layout = QVBoxLayout()
        self.axis_key_widget = QWidget()
        self.axis_key_widget.setLayout(axis_key_layout)

        self.x_axis_widget = self.create_axis_label('x̂', Qt.red)
        self.y_axis_widget = self.create_axis_label('ŷ', Qt.green)
        self.z_axis_widget = self.create_axis_label('ẑ', Qt.blue)

        axis_key_layout.addWidget(self.x_axis_widget)
        axis_key_layout.addWidget(self.y_axis_widget)
        axis_key_layout.addWidget(self.z_axis_widget)

        axis_key_dock = QDockWidget("Axis Key", self)
        axis_key_dock.setWidget(self.axis_key_widget)
        self.addDockWidget(Qt.RightDockWidgetArea, axis_key_dock)

    def add_point(self):

        #initialize pose
        pose = SE3()
        point = gl.GLScatterPlotItem(pos=pose.t, color=(1, 1, 1, 1), size=10)

        self.points.append(point)
        self.poses.append(pose)

        self.select_point_options.addItem(str(len(self.points)))

        self.update_plot()

    def move_point_x_pos(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Tx(.2)
        self.points[self.current_point].pos = self.poses[self.current_point].t 
        self.update_plot()

    def move_point_x_neg(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Tx(-.2)
        self.points[self.current_point].pos = self.poses[self.current_point].t 
        self.update_plot()

    def move_point_y_pos(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Ty(.2)
        self.points[self.current_point].pos = self.poses[self.current_point].t 
        self.update_plot()

    def move_point_y_neg(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Ty(-.2)
        self.points[self.current_point].pos = self.poses[self.current_point].t 
        self.update_plot()

    def move_point_z_pos(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Tz(.2)
        self.points[self.current_point].pos = self.poses[self.current_point].t 
        self.update_plot()

    def move_point_z_neg(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Tz(-.2)
        self.points[self.current_point].pos = self.poses[self.current_point].t 
        self.update_plot()

    def rotate_point_x(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Rx(math.pi/12)
        self.update_plot()

    def rotate_point_y(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Ry(math.pi/12)
        self.update_plot()

    def rotate_point_z(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Rz(math.pi/12)
        self.update_plot()

    def index_changed(self, index):
        self.points[self.current_point].color = (1,1,1,1)
        self.current_point = index
        self.points[self.current_point].color = (1,0,0,1)
        self.update_plot()

    def update_plot(self):
        self.plot_widget.clear()
        self.setCentralWidget(self.plot_widget)

        grid = gl.GLGridItem()
        self.plot_widget.addItem(grid)

        index = 0
        for point in self.points:
            md = gl.MeshData.sphere(rows=10, cols=20, radius=0.9)
            center = point.pos
            
            #plot sphere
            m1 = gl.GLMeshItem(
                meshdata=md,
                smooth=True,
                color=(0.4, 0.4, 0.4, 0.1),
                shader="balloon",
                glOptions="additive",
            )
            m1.translate(*center)
            self.plot_widget.addItem(m1)

            #plot point
            self.plot_widget.addItem(point)

            # -------- plot all the poses --------
            pose = self.poses[index]
            center = np.append(point.pos, 1)

            #endpoints
            p_x = np.append(point.pos + [1,0,0], 1)
            p_y = np.append(point.pos + [0,1,0], 1)
            p_z = np.append(point.pos + [0,0,1], 1)

            #apply transformation to a centered point
            rotated_point_x = np.dot(pose, p_x - center)
            rotated_point_y = np.dot(pose, p_y - center)
            rotated_point_z = np.dot(pose, p_z - center)
    
            #translate back from center
            rotated_point_x += center
            rotated_point_y += center
            rotated_point_z += center

            xLine = gl.GLLinePlotItem(pos=np.array([point.pos, rotated_point_x[:3]]), color=(1, 0, 0, 1)) 
            yLine = gl.GLLinePlotItem(pos=np.array([point.pos, rotated_point_y[:3]]), color=(0, 1, 0, 1)) 
            zLine = gl.GLLinePlotItem(pos=np.array([point.pos, rotated_point_z[:3]]), color=(0, 0, 1, 1)) 

            self.plot_widget.addItem(xLine)
            self.plot_widget.addItem(yLine)
            self.plot_widget.addItem(zLine)

            index = index + 1

    def create_axis_label(self, text, color):
        line_pixmap = QPixmap(20, 2)
        line_pixmap.fill(color)
        line_label = QLabel()
        line_label.setPixmap(line_pixmap)

        text_label = QLabel(text)

        layout = QHBoxLayout()
        layout.addWidget(line_label)
        layout.addWidget(text_label)
        layout.setContentsMargins(0, 0, 0, 0) 

        widget = QWidget()
        widget.setLayout(layout)
        return widget

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PointEditorWindow()
    window.show()
    sys.exit(app.exec_())