import sys
import numpy as np
import pyqtgraph.opengl as gl
from PyQt5 import QtCore as qc
from PyQt5.QtWidgets import QApplication, QMainWindow
from MyWindowOrigami import ClickableGLViewWidget

class PointEditorWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Point Editor")
        self.setGeometry(100, 100, 800, 600)

        self.plot_widget = ClickableGLViewWidget()
        self.setCentralWidget(self.plot_widget)
        self.plot_widget.setBackgroundColor(255,255,255, 255)

        num_points = 1000
        pos = np.random.random((num_points, 3)) * 50  # 3D positions in a range of 0 to 50

        # Create a GLScatterPlotItem
        scatter = gl.GLScatterPlotItem(pos=pos, size=2, color=(1,0,0,1), pxMode=False)

        # Add the scatter plot item to the GLViewWidget
        self.plot_widget.addItem(scatter)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PointEditorWindow()
    window.show()
    sys.exit(app.exec_())