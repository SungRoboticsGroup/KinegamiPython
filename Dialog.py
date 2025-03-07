"""
@author: Raymond Feng, Andy Wang, Daniel Feshbach
"""

import numpy as np
from PyQt5.QtWidgets import QVBoxLayout, QPushButton, QHBoxLayout, QLabel, QDialog, QLineEdit
from OpenGL.GL import *
from OpenGL.GLU import *
from spatialmath import SE3
import math
from PathCSC import *
from KinematicChain import *
from scipy.spatial.transform import Rotation as R
from testqtgraph import *
from style import *

class EditJointStateDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Edit Joint State')
        self.setGeometry(100, 100, 200, 100)

        layout = QVBoxLayout()

        self.state_input = QLineEdit(self)
        self.state_input.setPlaceholderText('Enter new joint state')
        layout.addWidget(QLabel('State:'))
        layout.addWidget(self.state_input)

        self.apply_button = QPushButton('Apply', self)
        self.apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(self.apply_button)

        self.cancel_button = QPushButton('Cancel', self)
        self.cancel_button.clicked.connect(self.reject)
        layout.addWidget(self.cancel_button)

        self.setLayout(layout)

    def onApplyClicked(self):
        self.accept()

    def get_state(self):
        try:
            state = float(self.state_input.text())
            return state
        except ValueError:
            self.show_error("Please enter a valid join state.")
            # QMessageBox.warning(self, "Invalid Input", "Please enter a valid join state.")
            self.exec_() 
            return None
        

class AddJointDialog(QDialog):
    jointToAdd = None

    def __init__(self, parent=None):
        super().__init__(parent)

    def getJoint(self):
        return self.jointToAdd
        
    def parse_angle(self, exp):
        try:
            result = eval(exp, {'np': np})
            return result
        except Exception as e:
            print("Error:", e)
            return None
    
    def parse_pose(self, exp):
        try:
            return eval(exp)
        except Exception as e:
            print("Error:", e)
            return None
        
class AddPrismaticDialog(AddJointDialog):
    def __init__(self, numSides, r, prevJoint : Joint = None, add_to_root = False):
        super().__init__()
        self.setWindowTitle('Add new prismatic joint')
        self.setGeometry(100, 100, 300, 100)

        layout = QVBoxLayout()
        
        length_layout = QHBoxLayout()
        length_label = QLabel("Neutral Length (default: 3r):")
        self.length_input = QLineEdit()
        length_layout.addWidget(length_label)
        length_layout.addWidget(self.length_input)
        layout.addLayout(length_layout)

        numLayers_layout = QHBoxLayout()
        numLayers_label = QLabel("Number of Layers (default: 3):")
        self.numLayers_input = QLineEdit()
        numLayers_layout.addWidget(numLayers_label)
        numLayers_layout.addWidget(self.numLayers_input)
        layout.addLayout(numLayers_layout)

        angle_layout = QHBoxLayout()
        angle_label = QLabel("Cone Angle (degrees, default: 60):")
        self.angle_input = QLineEdit()
        angle_layout.addWidget(angle_label)
        angle_layout.addWidget(self.angle_input)
        layout.addLayout(angle_layout)

        apply_button = QPushButton('Add Prismatic Joint')
        apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(apply_button)

        self.setLayout(layout)

        self.numSides = numSides
        self.r = r
        self.prevJoint = prevJoint
        self.add_to_root = add_to_root
                
    def onApplyClicked(self):
        try:            
            neutralLength = 3*self.r if self.length_input.text()=="" else float(self.length_input.text())
            numLayers = 3 if self.numLayers_input.text()=="" else int(self.numLayers_input.text())
            coneAngleText = 60 if self.angle_input.text()=="" else float(self.angle_input.text())

            if (self.prevJoint is None):
                pose = SE3()
            else:
                distance = 4 * self.r + norm(self.prevJoint.distalPosition()-self.prevJoint.Pose.t) + neutralLength/2
                if self.add_to_root: distance *= -1
                pose = SE3(0,0,distance)
                if self.prevJoint.pathIndex() == 0:
                    pose = SE3.Ry(np.pi/2) @ pose

            self.jointToAdd = PrismaticJoint(self.numSides, self.r, neutralLength, numLayers, math.radians(coneAngleText), pose)
            self.accept()
        except ValueError:
            self.show_error('Please enter valid numbers.')
            # error_dialog = ErrorDialog('Please enter valid integers.')
            # error_dialog.exec_()

class AddRevoluteDialog(AddJointDialog):
    def __init__(self, numSides, r, prevJoint : Joint = None, add_to_root = False):
        super().__init__()
        self.setWindowTitle('Add new revolute joint')
        self.setGeometry(100, 100, 300, 100)

        layout = QVBoxLayout()
        
        angle_layout = QHBoxLayout()
        angle_label = QLabel("Total Bending Angle (degrees, default: 180):")
        self.angle_input = QLineEdit()
        angle_layout.addWidget(angle_label)
        angle_layout.addWidget(self.angle_input)
        layout.addLayout(angle_layout)

        apply_button = QPushButton('Add Revolute Joint')
        apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(apply_button)

        self.setLayout(layout)

        self.numSides = numSides
        self.r = r
        self.prevJoint = prevJoint
        self.add_to_root = add_to_root
        #self.prevClass = prevClass
        
    def onApplyClicked(self):
        bendingAngleText = 180 if self.angle_input.text()=="" else float(self.angle_input.text())
        
        
        self.jointToAdd = RevoluteJoint(self.numSides, self.r, math.radians(bendingAngleText), SE3())

        if not self.prevJoint is None:
            distance = 4 * self.r + norm(self.prevJoint.distalPosition()-self.prevJoint.Pose.t) + self.jointToAdd.neutralLength/2
            if self.add_to_root: distance *= -1
            pose = SE3(distance,0,0)
            if self.prevJoint.pathIndex() == 2:
                pose = SE3.Ry(-np.pi/2) @ pose
            self.jointToAdd.Pose = pose
    
        self.accept()

class AddTipDialog(AddJointDialog):
    isStart = True

    def __init__(self, numSides, r, prevJoint : Joint = None, add_to_root = False):
        super().__init__()
        self.setWindowTitle('Add new joint')
        self.setGeometry(100, 100, 300, 100)

        self.numSides = numSides
        self.r = r
        self.isStart = True
        #self.prevClass = prevClass
        self.prevJoint = prevJoint

        layout = QVBoxLayout()

        length_layout = QHBoxLayout()
        length_label = QLabel("Length:")
        self.length_input = QLineEdit()
        length_layout.addWidget(length_label)
        length_layout.addWidget(self.length_input)
        layout.addLayout(length_layout)
            
        apply_button = QPushButton('Add Tip')
        apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(apply_button)

        self.setLayout(layout)
        self.add_to_root = add_to_root

    def onApplyClicked(self):
        try:
            length = float(self.length_input.text())
            if self.prevJoint is None:
                self.jointToAdd = StartTip(self.numSides, self.r, SE3(), length=length)
            else:
                distance = 4*self.r + norm(self.prevJoint.distalPosition()-self.prevJoint.Pose.t) + length/2
                if self.add_to_root: distance *= -1
                pose = SE3(0,0,distance)
                if self.prevJoint.pathIndex() == 0:
                    pose = SE3.Ry(np.pi/2) @ pose
                self.jointToAdd = EndTip(self.numSides, self.r, pose, length=length)

            self.accept()
        except ValueError:
            self.show_error('Please enter valid integers.')
            # error_dialog = ErrorDialog('Please enter valid integers.')
            # error_dialog.exec_()
    
    def updateVariable(self):
        if self.radio_start.isChecked():
            self.isStart = True
        elif self.radio_end.isChecked():
            self.isStart = False