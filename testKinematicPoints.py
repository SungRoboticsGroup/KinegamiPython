# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:32:19 2023

@author: Daniel Feshbach
"""
import sys
import matplotlib
matplotlib.use('Qt5Agg')

from KinematicChain import *
from PyQt5 import QtWidgets, QtGui, QtCore

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

r = 1
numSides = 6

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

class MainWindow(QtWidgets.QMainWindow):

    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        #creates the canvas
        sc = MplCanvas(self, width=5, height=4, dpi=100)

        #creates the toolbar
        toolbar = NavigationToolbar(sc, self)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(toolbar)
        layout.addWidget(sc)

        # Create a placeholder widget to hold our toolbar and canvas.
        widget = QtWidgets.QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)

        #self.show()

        # chain whose root is a waypoint at the global origin
        self = KinematicChain(StartFingertip(numSides, r, Pose=SE3(), length=0.5)) 
        prismaticIndex = self.append(PrismaticJoint(numSides, r, neutralLength=3, numLayers=6, 
                                coneAngle=np.pi/4, Pose= SE3() ) )
        revoluteIndex = self.append(RevoluteJoint(numSides, r, np.pi, 
                                                  SE3.Ry(np.pi/4)))
        self.show(block=True)

app = QtWidgets.QApplication(sys.argv)
w = MainWindow()
app.exec_()
