import sys
import numpy as np
import pyqtgraph.opengl as gl
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton, QDockWidget, QComboBox, QDoubleSpinBox, QHBoxLayout, QLabel
from PyQt5.QtGui import QColor
from PyQt5.QtCore import Qt

'''
- import SE3
- import numpy as np
- SE3.Rx() to rotate, @ to apply, * to apply to numpy vectors
'''

class PointEditorWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Point Editor")
        self.setGeometry(100, 100, 800, 600)

        self.plot_widget = gl.GLViewWidget()
        self.setCentralWidget(self.plot_widget)

        grid = gl.GLGridItem()
        self.plot_widget.addItem(grid)

        self.current_point = 0

        self.points = []

        self.add_point_button = QPushButton("Add Point")
        self.add_point_button.clicked.connect(self.add_point)
        self.select_point_options = QComboBox()
        self.select_point_options.activated.connect(self.index_changed)
        self.select_point_options.currentIndexChanged.connect(self.index_changed)

        self.move_x_spinbox = QDoubleSpinBox()
        self.move_x_spinbox.valueChanged.connect(self.move_point_x)
        self.move_x_spinbox.setRange(-20, 20)
        self.move_x_spinbox.setSingleStep(0.1)
        self.move_x_label = QLabel("X: ")

        self.move_y_spinbox = QDoubleSpinBox()
        self.move_y_spinbox.valueChanged.connect(self.move_point_y)
        self.move_y_spinbox.setRange(-20, 20)
        self.move_y_spinbox.setSingleStep(0.1)
        self.move_y_label = QLabel("Y: ")

        self.move_z_spinbox = QDoubleSpinBox()
        self.move_z_spinbox.valueChanged.connect(self.move_point_z)
        self.move_z_spinbox.setRange(-20, 20)
        self.move_z_spinbox.setSingleStep(0.1)
        self.move_z_label = QLabel("Z: ")

        self.rotate_button = QPushButton("Rotate")
        self.rotate_button.clicked.connect(self.rotate_points)

        button_layout = QVBoxLayout()
        button_layout.addWidget(self.select_point_options)
        button_layout.addWidget(self.add_point_button)

        x_layout = QHBoxLayout()
        x_layout.addWidget(self.move_x_label)
        x_layout.addWidget(self.move_x_spinbox) 

        y_layout = QHBoxLayout()
        y_layout.addWidget(self.move_y_label)
        y_layout.addWidget(self.move_y_spinbox) 

        z_layout = QHBoxLayout()
        z_layout.addWidget(self.move_z_label)
        z_layout.addWidget(self.move_z_spinbox) 

        button_layout.addLayout(x_layout)
        button_layout.addLayout(y_layout)
        button_layout.addLayout(z_layout)
        #button_layout.addWidget(self.rotate_button)

        button_widget = QWidget()
        button_widget.setLayout(button_layout)
        dock = QDockWidget("title")
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        dock.setWidget(button_widget)

    def add_point(self):
        point = gl.GLScatterPlotItem(pos=np.array([0,0,0]), color=(1, 1, 1, 1), size=10)
        self.points.append(point)
        self.select_point_options.addItem(str(len(self.points)))

        self.update_plot()

    def move_point_x(self, value):
        self.points[self.current_point].pos[0] = value
        self.update_plot()

    def move_point_y(self, value):
        self.points[self.current_point].pos[1] = value
        self.update_plot()

    def move_point_z(self, value):
        self.points[self.current_point].pos[2] = value
        self.update_plot()

    def rotate_points(self):
        angle = 10 if QApplication.keyboardModifiers() == Qt.ShiftModifier else 1
        point = self.points[self.current_point]
        rot = point.opts['rotate']
        point.opts['rotate'] = (rot[0], rot[1] + angle, rot[2])

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

        for point in self.points:
            md = gl.MeshData.sphere(rows=10, cols=20, radius=0.9)
            center = point.pos

            m1 = gl.GLMeshItem(
                meshdata=md,
                smooth=True,
                color=(0.4, 0.4, 0.4, 0.1),
                shader="balloon",
                glOptions="additive",
            )
            m1.translate(*center)
            self.plot_widget.addItem(m1)

            self.plot_widget.addItem(point)
            self.plot_widget.addItem(gl.GLLinePlotItem(pos=np.array([point.pos, [point.pos[0]+1, point.pos[1], point.pos[2]]]), color=(1, 0, 0, 1)))
            self.plot_widget.addItem(gl.GLLinePlotItem(pos=np.array([point.pos, [point.pos[0], point.pos[1]+1, point.pos[2]]]), color=(0, 1, 0, 1)))
            self.plot_widget.addItem(gl.GLLinePlotItem(pos=np.array([point.pos, [point.pos[0], point.pos[1], point.pos[2]+1]]), color=(0, 0, 1, 1)))

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PointEditorWindow()
    window.show()
    sys.exit(app.exec_())