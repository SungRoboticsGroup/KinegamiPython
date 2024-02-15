# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 17:32:19 2023

@author: Daniel Feshbach
"""

import sys
import matplotlib
import pyqtgraph as pg
import pyqtgraph.opengl as gl

from PyQt5 import QtWidgets, QtGui, QtCore
from KinematicChain import *
from PyQt5.QtWidgets import QRadioButton, QVBoxLayout, QApplication, QWidget, QPushButton
from scipy.spatial.transform import Rotation as R

matplotlib.use('Qt5Agg')

numSides = 6
r = 1

joints = []

class Options(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        # Button with a parent widget

        self.joint1Option = QRadioButton('Joint 1')
        self.joint2Option = QRadioButton('Joint 2')

        self.joint1Option.setChecked(True)

        self.joint1Option.toggled.connect(self.radio_button_toggled)
        self.joint2Option.toggled.connect(self.radio_button_toggled)

        self.moveButton = QPushButton("transform")
        self.moveButton.clicked.connect(self.buttonClicked)
 
        layout = QVBoxLayout()
        layout.addWidget(self.joint1Option)
        layout.addWidget(self.joint2Option)
        layout.addWidget(self.moveButton)

        self.setLayout(layout)
    
    def radio_button_toggled(self):
        # Check which radio button is toggled and print its text
        if self.joint1Option.isChecked():
            super.selectedJointIndex = 1
        elif self.joint2Option.isChecked():
            super.selectedJointIndex = 2

    def buttonClicked(self): {
        print(super.selectedJointIndex)
    }

class JointPlot(QWidget):
    def __init__(self):
        self.app = None
        self.selected_joint_index = None

    def setup(self):
        self.app = gl.GLViewWidget()
        self.app.setWindowTitle('Joint Plot')
        self.app.setGeometry(100, 100, 800, 600)

        #3d grid
        grid = gl.GLGridItem()
        self.app.addItem(grid)

    def show(self):
        self.app.show()

    def add_joint(self, joint):
        joints.append(joint)

    def plot_joints(self):
        for joint in joints :
            pos = joint.distalPosition()
            self.app.addItem(gl.GLScatterPlotItem(pos=pos, color=(1.0, 0.0, 0.0, 1.0), size=10))
        
        #plot point
        

        #plot arrows
        #axes = joint.DubinsFrame().A

        #finish code here

if __name__ == '__main__':
    app = QApplication(sys.argv)

    selectedJointIndex = 0

    joint_plot = JointPlot()
    joint_plot.setup()
    joint_plot.show()

    buttons = Options()
    buttons.show()

    #plot joints
    joint1 = Waypoint(numSides, r, Pose=SE3(1,1,1), pathIndex=2)
    joint2 = Waypoint(numSides, r, Pose=SE3(3,3,3), pathIndex=2)
    #joint2.transformPoseBy(SE3(-2,3,5) @ SE3.Ry(np.pi/4))
    
    joint_plot.add_joint(joint1)
    joint_plot.add_joint(joint2)

    joint_plot.plot_joints()

    #button_list_widget = joint_plot.create_button_list()
    #button_list_widget.show()

    sys.exit(app.exec_())