import numpy as np
from PyQt5.QtWidgets import QVBoxLayout, QPushButton, QHBoxLayout, QLabel, QDockWidget, QWidget, QLineEdit
from OpenGL.GL import *
from OpenGL.GLU import *
from spatialmath import SE3
import math
from PathCSC import *
from KinematicChain import *
from scipy.spatial.transform import Rotation as R
from testqtgraph import *
from style import *

class AddJointMenu(QWidget):
    jointToAdd = None

    def __init__(self, parent=None):
        super().__init__(parent)
        self.numSides = 0
        self.r = 0
        self.prevJoint = None
        self.add_to_root = False

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
        
    def update(self):
        self.numSides = self.window().num_sides
        self.r = self.window().radius

        if (self.window().chain and len(self.window().chain.Joints) > 0):
            if not self.add_to_root:
                self.prevJoint = self.window().chain.Joints[-1]
                self.add_to_root = False
            else:
                self.prevJoint = self.window().chain.Joints[0]
                self.add_to_root = True
        else: 
            self.prevJoint = None
            self.add_to_root = False
    
class AddPrismaticMenu(AddJointMenu):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.add_joint = self.window().add_joint

        self.initUI()        

    def initUI(self):
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

        cancel_button = QPushButton('Cancel')
        cancel_button.clicked.connect(self.window().add_prismatic_toggle)
        layout.addWidget(cancel_button)

        self.setLayout(layout)

        self.update()

    def onApplyClicked(self):
        try:
            self.update()

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
            self.add_joint(self.jointToAdd)
            self.window().add_prismatic_toggle()
        except ValueError:
            self.show_error('Please enter valid numbers.')
            # error_dialog = ErrorDialog('Please enter valid integers.')
            # error_dialog.exec_()


class AddRevoluteMenu(AddJointMenu):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.add_joint = self.window().add_joint

        self.initUI()   

    def initUI(self):
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

        cancel_button = QPushButton('Cancel')
        cancel_button.clicked.connect(self.window().add_revolute_toggle)
        layout.addWidget(cancel_button)

        self.setLayout(layout)

        self.update()
        
    def onApplyClicked(self):
        self.update()

        bendingAngleText = 180 if self.angle_input.text()=="" else float(self.angle_input.text())
        
        
        self.jointToAdd = RevoluteJoint(self.numSides, self.r, math.radians(bendingAngleText), SE3())

        if not self.prevJoint is None:
            distance = 4 * self.r + norm(self.prevJoint.distalPosition()-self.prevJoint.Pose.t) + self.jointToAdd.neutralLength/2
            if self.add_to_root: distance *= -1
            pose = SE3(distance,0,0)
            if self.prevJoint.pathIndex() == 2:
                pose = SE3.Ry(-np.pi/2) @ pose
            self.jointToAdd.Pose = pose
    
        self.add_joint(self.jointToAdd)
        self.window().add_revolute_toggle()

class AddTipMenu(AddJointMenu):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.add_joint = self.window().add_joint

        self.isStart = True

        self.initUI()

    def initUI(self):
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

        cancel_button = QPushButton('Cancel')
        cancel_button.clicked.connect(self.window().add_tip_toggle)
        layout.addWidget(cancel_button)

        self.setLayout(layout)

    def onApplyClicked(self):
        try:
            self.update()

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

            self.add_joint(self.jointToAdd)
            self.window().add_tip_toggle()

        except ValueError:
            self.show_error('Please enter valid integers.')
            # error_dialog = ErrorDialog('Please enter valid integers.')
            # error_dialog.exec_()