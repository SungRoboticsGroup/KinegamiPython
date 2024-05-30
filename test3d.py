import sys
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton
from PyQt5.QtGui import QSurfaceFormat
import pyqtgraph.opengl as gl
from pyqtgraph.Qt import QtCore
from OpenGL.GL import *
from OpenGL.GLU import *

class ClickableGLViewWidget(gl.GLViewWidget):
    def __init__(self, parent=None):
        super(ClickableGLViewWidget, self).__init__(parent)
        fmt = QSurfaceFormat()
        fmt.setDepthBufferSize(24)
        fmt.setVersion(3, 3)
        fmt.setProfile(QSurfaceFormat.CompatibilityProfile)
        QSurfaceFormat.setDefaultFormat(fmt)
        self.setFormat(fmt)
        self.locked = False
        self.mesh = None 

    def toggle_lock(self):
        self.locked = not self.locked
        print("Screen lock toggled:", "Locked" if self.locked else "Unlocked")

    def mousePressEvent(self, event):
        if self.locked and event.button() == QtCore.Qt.LeftButton:
            self.project_click(event.pos())
        else:
            super(ClickableGLViewWidget, self).mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if not self.locked:
            super(ClickableGLViewWidget, self).mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if not self.locked:
            super(ClickableGLViewWidget, self).mouseReleaseEvent(event)

    def project_click(self, position):
        viewport = glGetIntegerv(GL_VIEWPORT)
        modelview = glGetDoublev(GL_MODELVIEW_MATRIX)
        projection = glGetDoublev(GL_PROJECTION_MATRIX)
        winX = position.x()
        winY = viewport[3] - position.y()
        winZ = glReadPixels(winX, winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT)
        winZ = winZ[0][0]

        if winZ == 1.0:
            print("Clicked on empty space")
        else:
            x, y, z = gluUnProject(winX, winY, winZ, modelview, projection, viewport)
            print(f"Clicked 3D coordinates: {x}, {y}, {z}")

            sphere_center = (-1.5, 0.5)  
            sphere_radius = 0.1  
            distance = ((x - sphere_center[0]) ** 2 + (y - sphere_center[1]) ** 2) ** 0.5

            if distance <= sphere_radius:
                print("Clicked inside the sphere")
                self.mesh.setColor((1, 1, 0, 1)) 
            else:
                print("Clicked outside the sphere")
                self.mesh.setColor((1, 0, 0, 1)) 

class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setWindowTitle('3D Click Position')
        self.setGeometry(100, 100, 800, 600)
        self.glWidget = ClickableGLViewWidget()
        self.setCentralWidget(self.glWidget)

        self.toggleButton = QPushButton('Toggle Lock', self)
        self.toggleButton.clicked.connect(self.glWidget.toggle_lock)
        self.toggleButton.setGeometry(650, 550, 140, 40)

        grid = gl.GLGridItem()
        self.glWidget.addItem(grid)

        object = gl.MeshData.sphere(rows=10, cols=10)
        mesh = gl.GLMeshItem(meshdata=object, color=(1, 0, 0, 1), shader='shaded')
        self.glWidget.addItem(mesh)
        self.glWidget.mesh = mesh 

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
