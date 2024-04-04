"""
@author: Raymond Feng
"""

import sys
import numpy as np
import pyqtgraph.opengl as gl
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton, QDockWidget, QComboBox, QHBoxLayout, QLabel, QDialog, QLineEdit, QCheckBox, QMessageBox
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap, QPainter
from spatialmath import SE3
import math
from Joint import *
from PathCSC import *
from KinematicChain import *
import re

class TransformDialog(QDialog):
    propogateTransform = True

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Transform Joint')
        self.setGeometry(100, 100, 300, 250)  

        layout = QVBoxLayout()

        self.operationSelection = QComboBox(self)
        self.operationSelection.addItems(["Translate", "Rotate"])
        layout.addWidget(QLabel('Operation:'))
        layout.addWidget(self.operationSelection)

        self.axisSelection = QComboBox(self)
        self.axisSelection.addItems(["X", "Y", "Z"])
        self.axisSelectionLabel = QLabel('Axis:')
        layout.addWidget(self.axisSelectionLabel)
        layout.addWidget(self.axisSelection)
        self.axisSelection.setVisible(False) 
        self.axisSelectionLabel.setVisible(False)  

        self.operationSelection.currentTextChanged.connect(self.toggleAxisSelection)

        self.propogateTransformCheckbox = QCheckBox("Propagate")

        self.propogateTransformCheckbox.setChecked(TransformDialog.propogateTransform)

        layout.addWidget(self.propogateTransformCheckbox)

        self.setupTransformInputAndButtons(layout)
        self.setLayout(layout)

    def toggleAxisSelection(self):
        isRotateSelected = self.operationSelection.currentText() == "Rotate"
        self.axisSelection.setVisible(isRotateSelected)
        self.axisSelectionLabel.setVisible(isRotateSelected)

    def setupTransformInputAndButtons(self, layout):
        self.transform_input = QLineEdit(self)
        self.transform_input.setPlaceholderText('x, y, z for Translate or degree for Rotate:')
        layout.addWidget(QLabel('Values:'))
        layout.addWidget(self.transform_input)

        self.apply_button = QPushButton('Apply', self)
        self.apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(self.apply_button)

        self.cancel_button = QPushButton('Cancel', self)
        self.cancel_button.clicked.connect(self.reject)
        layout.addWidget(self.cancel_button)

    def onApplyClicked(self):
        TransformDialog.propogateTransform = self.propogateTransformCheckbox.isChecked()
        self.accept()


    def get_transform(self):
        while True:
            try:
                operation = self.operationSelection.currentText()
                input_text = self.transform_input.text()
                values = [float(val.strip()) for val in input_text.split(',')]
                if operation == "Translate":
                    if len(values) != 3:
                        raise ValueError("Invalid Input", "Please enter exactly 3 values for translation.")
                    transformation = SE3.Trans(values[0], values[1], values[2])
                else:
                    if len(values) != 1:
                        raise ValueError("Invalid Input", "Please enter exactly 1 value for rotation angle in degrees.")
                    angle_degrees = values[0]
                    angle_radians = math.radians(angle_degrees) 
                    axis = self.axisSelection.currentText()
                    if axis == "X":
                        transformation = SE3.Rx(angle_radians)
                    elif axis == "Y":
                        transformation = SE3.Ry(angle_radians)
                    else:
                        transformation = SE3.Rz(angle_radians)
                return transformation
            except ValueError as e:
                QMessageBox.warning(self, "Invalid Input", str(e))
                self.exec_()
                return None

class RotationDialog(QDialog):
    propogateRotation = True
    applyToPreviousWaypointRotation = False

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Rotate Joint Along Axis Of Motion')
        self.setGeometry(100, 100, 200, 100)

        layout = QVBoxLayout()

        self.propogateRotationCheckbox = QCheckBox("Propogate")
        self.applyToPreviousWaypointRotationCheckbox = QCheckBox("Apply to Previous Waypoint")

        self.propogateRotationCheckbox.setChecked(RotationDialog.propogateRotation)
        self.applyToPreviousWaypointRotationCheckbox.setChecked(RotationDialog.applyToPreviousWaypointRotation)

        layout.addWidget(self.propogateRotationCheckbox)
        layout.addWidget(self.applyToPreviousWaypointRotationCheckbox)

        self.setupAngleInputAndButtons(layout)
        self.setLayout(layout)

    def setupAngleInputAndButtons(self, layout):
        self.angle_input = QLineEdit(self)
        self.angle_input.setPlaceholderText('Enter angle in degrees')
        layout.addWidget(QLabel('Angle:'))
        layout.addWidget(self.angle_input)

        self.apply_button = QPushButton('Apply', self)
        self.apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(self.apply_button)

        self.cancel_button = QPushButton('Cancel', self)
        self.cancel_button.clicked.connect(self.reject)
        layout.addWidget(self.cancel_button)

    def onApplyClicked(self):
        RotationDialog.propogateRotation = self.propogateRotationCheckbox.isChecked()
        RotationDialog.applyToPreviousWaypointRotation = self.applyToPreviousWaypointRotationCheckbox.isChecked()
        self.accept()

    def get_angle(self):
        while True:
            try:
                angle = float(self.angle_input.text())
                return angle
            except ValueError:
                QMessageBox.warning(self, "Invalid Input", "Please enter a valid angle in degrees.")
                self.exec_() 
                return None
                
class TranslationDialog(QDialog):
    propogateTranslation = True
    applyToPreviousWaypointTranslation = False

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Translate Joint Along Axis of Motion')
        self.setGeometry(100, 100, 200, 100)

        layout = QVBoxLayout()

        self.propogateTranslationCheckbox = QCheckBox("Propogate")
        self.applyToPreviousWaypointTranslationCheckbox = QCheckBox("Apply to Previous Waypoint")

        self.propogateTranslationCheckbox.setChecked(TranslationDialog.propogateTranslation)
        self.applyToPreviousWaypointTranslationCheckbox.setChecked(TranslationDialog.applyToPreviousWaypointTranslation)

        layout.addWidget(self.propogateTranslationCheckbox)
        layout.addWidget(self.applyToPreviousWaypointTranslationCheckbox)

        self.setupDistanceInputAndButtons(layout)
        self.setLayout(layout)

    def setupDistanceInputAndButtons(self, layout):
        self.distance_input = QLineEdit()
        self.distance_input.setPlaceholderText('Enter translation distance')
        layout.addWidget(QLabel('Distance:'))
        layout.addWidget(self.distance_input)

        apply_button = QPushButton('Apply')
        apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(apply_button)

        cancel_button = QPushButton('Cancel')
        cancel_button.clicked.connect(self.reject)
        layout.addWidget(cancel_button)

    def onApplyClicked(self):
        TranslationDialog.propogateTranslation = self.propogateTranslationCheckbox.isChecked()
        TranslationDialog.applyToPreviousWaypointTranslation = self.applyToPreviousWaypointTranslationCheckbox.isChecked()
        self.accept()

    def get_distance(self):
        while True:
            try:
                distance = float(self.distance_input.text())
                return distance
            except ValueError:
                QMessageBox.warning(self, "Invalid Input", "Please enter a valid translation distance.")
                self.exec_() 
                return None

class DeleteDialog(QDialog):

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Confirm Delete')
        self.setGeometry(100, 100, 200, 100)

        layout = QVBoxLayout()
        layout.addWidget(QLabel('Are you sure you want to delete the joint?'))

        apply_button = QPushButton('Delete')
        apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(apply_button)

        cancel_button = QPushButton('Cancel')
        cancel_button.clicked.connect(self.reject)
        layout.addWidget(cancel_button)

        self.setLayout(layout)

    def onApplyClicked(self):
        self.accept()

class SuccessDialog(QDialog):
    def __init__(self, message, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Success')
        self.setGeometry(100, 100, 200, 100)

        layout = QVBoxLayout()
        
        self.message_label = QLabel(message)
        self.message_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.message_label)
        
        close_button = QPushButton('Close')
        close_button.clicked.connect(self.close)
        layout.addWidget(close_button)

        self.setLayout(layout)

class ErrorDialog(QDialog):
    def __init__(self, message, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Error')
        self.setGeometry(100, 100, 200, 100)

        layout = QVBoxLayout()
        
        self.message_label = QLabel(message)
        self.message_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(self.message_label)
        
        close_button = QPushButton('Close')
        close_button.clicked.connect(self.close)
        layout.addWidget(close_button)

        self.setLayout(layout)

class AddJointDialog(QDialog):
    jointToAdd = None

    def __init__(self, parent=None):
        super().__init__(parent)
    
    def getJoint(self):
        try:
            return self.jointToAdd
        except ValueError:
            QMessageBox.warning(self, "ERROR", "Invalid input")
            self.exec_() 
            return None
        
    def parse_angle(self, exp):
        try:
            result = eval(exp, {'np': np})
            return result
        except Exception as e:
            print("Error:", e)
            return None
    
    def parse_individual_pose(self, exp):
        try:
            if exp.startswith("SE3.Trans"):
                translation_str = exp.split('([')[1].split('])')[0]
                translation = [int(val) for val in translation_str.split(',')]
                return SE3.Trans(translation)
            
            if exp.startswith("SE3.Rx"):
                rot_str = exp.split('(')[1].split(')')[0]
                rotation = self.parse_angle(rot_str)
                return SE3.Rx(rotation)
            
            if exp.startswith("SE3.Ry"):
                rot_str = exp.split('(')[1].split(')')[0]
                rotation = self.parse_angle(rot_str)
                return SE3.Ry(rotation)
            
            if exp.startswith("SE3.Rz"):
                rot_str = exp.split('(')[1].split(')')[0]
                rotation = self.parse_angle(rot_str)
                return SE3.Rz(rotation)
            
        except Exception as e:
            print("Error:", e)
            return None
    
    def parse_pose(self, exp):
        pattern = r'\s*\*\s*'
        poses = re.split(pattern, exp)
        
        result = SE3()
        for pose in poses:
            parsedPose = self.parse_individual_pose(pose)
            result = result * parsedPose

        return result

class AddPrismaticDialog(AddJointDialog):
    def __init__(self, numSides, r, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Add new prismatic joint')
        self.setGeometry(100, 100, 300, 100)

        layout = QVBoxLayout()
        
        length_layout = QHBoxLayout()
        length_label = QLabel("Neutral Length:")
        self.length_input = QLineEdit()
        length_layout.addWidget(length_label)
        length_layout.addWidget(self.length_input)
        layout.addLayout(length_layout)

        numLayers_layout = QHBoxLayout()
        numLayers_label = QLabel("Number of Layers:")
        self.numLayers_input = QLineEdit()
        numLayers_layout.addWidget(numLayers_label)
        numLayers_layout.addWidget(self.numLayers_input)
        layout.addLayout(numLayers_layout)

        angle_layout = QHBoxLayout()
        angle_label = QLabel("Cone Angle (radians):")
        self.angle_input = QLineEdit()
        angle_layout.addWidget(angle_label)
        angle_layout.addWidget(self.angle_input)
        layout.addLayout(angle_layout)

        pose_layout = QHBoxLayout()
        pose_label = QLabel("Pose (SE3):")
        self.pose_input = QLineEdit()
        pose_layout.addWidget(pose_label)
        pose_layout.addWidget(self.pose_input)
        layout.addLayout(pose_layout)

        apply_button = QPushButton('Add')
        apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(apply_button)

        self.setLayout(layout)

        self.numSides = numSides
        self.r = r
    
    def onApplyClicked(self):
        try:
            neutralLength = int(self.length_input.text())
            numLayers = int(self.numLayers_input.text())
            coneAngleText = self.angle_input.text()
            poseText = self.pose_input.text()

            self.jointToAdd = PrismaticJoint(self.numSides, self.r, neutralLength, numLayers, self.parse_angle(coneAngleText), self.parse_pose(poseText))
            self.accept()
        except ValueError:
            error_dialog = ErrorDialog('Please enter valid integers.')
            error_dialog.exec_()
        
class AddRevoluteDialog(AddJointDialog):
    def __init__(self, numSides, r, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Add new joint')
        self.setGeometry(100, 100, 300, 100)

        layout = QVBoxLayout()
        
        angle_layout = QHBoxLayout()
        angle_label = QLabel("Total Bending Angle (radians):")
        self.angle_input = QLineEdit()
        angle_layout.addWidget(angle_label)
        angle_layout.addWidget(self.angle_input)
        layout.addLayout(angle_layout)

        pose_layout = QHBoxLayout()
        pose_label = QLabel("Pose (SE3):")
        self.pose_input = QLineEdit()
        pose_layout.addWidget(pose_label)
        pose_layout.addWidget(self.pose_input)
        layout.addLayout(pose_layout)
        
        apply_button = QPushButton('Add')
        apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(apply_button)

        self.setLayout(layout)

        self.numSides = numSides
        self.r = r

    def onApplyClicked(self):
        bendingAngleText = self.angle_input.text()
        poseText = self.pose_input.text()

        self.jointToAdd = RevoluteJoint(self.numSides, self.r, self.parse_angle(bendingAngleText), self.parse_pose(poseText))
        self.accept()
        
class AddWaypointDialog(AddJointDialog):
    def __init__(self, numSides, r, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Add new joint')
        self.setGeometry(100, 100, 300, 100)

        layout = QVBoxLayout()
        
        pose_layout = QHBoxLayout()
        pose_label = QLabel("Pose (SE3):")
        self.pose_input = QLineEdit()
        pose_layout.addWidget(pose_label)
        pose_layout.addWidget(self.pose_input)
        layout.addLayout(pose_layout)
        
        apply_button = QPushButton('Add')
        apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(apply_button)

        self.setLayout(layout)

        self.numSides = numSides
        self.r = r

    def onApplyClicked(self):
        poseText = self.pose_input.text()

        self.jointToAdd = Waypoint(self.numSides, self.r, Pose=self.parse_pose(poseText))
        self.accept()

class AddTipDialog(AddJointDialog):
    isStart = True

    def __init__(self, numSides, r, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Add new joint')
        self.setGeometry(100, 100, 300, 100)

        layout = QVBoxLayout()

        length_layout = QHBoxLayout()
        length_label = QLabel("Length:")
        self.length_input = QLineEdit()
        length_layout.addWidget(length_label)
        length_layout.addWidget(self.length_input)
        layout.addLayout(length_layout)
        
        pose_layout = QHBoxLayout()
        pose_label = QLabel("Pose (SE3):")
        self.pose_input = QLineEdit()
        pose_layout.addWidget(pose_label)
        pose_layout.addWidget(self.pose_input)
        layout.addLayout(pose_layout)
        
        apply_button = QPushButton('Add')
        apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(apply_button)

        self.setLayout(layout)

        self.numSides = numSides
        self.r = r

    def onApplyClicked(self):
        try:
            neutralLength = float(self.length_input.text())
            poseText = self.pose_input.text()

            if (self.isStart):
                self.jointToAdd = StartTip(self.numSides, self.r, self.parse_pose(poseText), length=neutralLength)
            else:
                self.jointToAdd = EndTip(self.numSides, self.r, self.parse_pose(poseText), length=neutralLength)

            self.accept()
        except ValueError:
            error_dialog = ErrorDialog('Please enter valid integers.')
            error_dialog.exec_()

class CreateNewChainDialog(QDialog):
    numSides = 4
    r = 1

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Create new chain')
        self.setGeometry(100, 100, 300, 100)

        layout = QVBoxLayout()

        numSides_layout = QHBoxLayout()
        numSides_label = QLabel("Number of Sides:")
        self.numSides_input = QLineEdit()
        numSides_layout.addWidget(numSides_label)
        numSides_layout.addWidget(self.numSides_input)
        layout.addLayout(numSides_layout)

        r_layout = QHBoxLayout()
        r_label = QLabel("Radius:")
        self.r_input = QLineEdit()
        r_layout.addWidget(r_label)
        r_layout.addWidget(self.r_input)
        layout.addLayout(r_layout)
        
        create_button = QPushButton('Create')
        create_button.clicked.connect(self.onCreateClicked)
        layout.addWidget(create_button)

        self.setLayout(layout)

    def onCreateClicked(self):
        try:
            self.numSides = int(self.numSides_input.text())
            self.r = int(self.r_input.text())
            self.accept()
        except ValueError:
            error_dialog = ErrorDialog('Please enter valid integers.')
            error_dialog.exec_()

    def getNumSides(self):
        return self.numSides
    
    def getR(self):
        return self.r

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
        self.chain = None

        self.points = []
        self.poses = []

        self.numSides = 4
        self.r = 1

        # //////////////////////////////////    ADD JOINTS    ///////////////////////////////////
        self.add_prismatic = QPushButton("Add Prismatic Joint")
        self.add_revolute = QPushButton("Add Revolute Joint")
        self.add_waypoint = QPushButton("Add Waypoint")
        self.add_tip = QPushButton("Add Tip")
        self.create_new_chain = QPushButton("Create New Chain")

        add_joints_layout = QVBoxLayout()
        add_joints_layout.addWidget(self.add_prismatic)
        add_joints_layout.addWidget(self.add_revolute)
        add_joints_layout.addWidget(self.add_waypoint)
        add_joints_layout.addWidget(self.add_tip)
        add_joints_layout.addWidget(self.create_new_chain)

        self.add_prismatic.clicked.connect(self.add_prismatic_func)
        self.add_revolute.clicked.connect(self.add_revolute_func)
        self.add_waypoint.clicked.connect(self.add_waypoint_func)
        self.add_tip.clicked.connect(self.add_tip_func)
        self.create_new_chain.clicked.connect(self.create_new_chain_func)

        add_joints_dock = QDockWidget("Add joints", self)
        self.add_joints_widget = QWidget()
        self.add_joints_widget.setLayout(add_joints_layout)
        add_joints_dock.setWidget(self.add_joints_widget)

        self.addDockWidget(Qt.RightDockWidgetArea, add_joints_dock)

        # //////////////////////////////////    AXIS KEY    ////////////////////////////////////
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

        # ////////////////////////////////    EDIT JOINTS    ///////////////////////////////////
        self.select_joint_options = QComboBox()
        self.transform_joint_button = QPushButton("Transform Joint")
        self.rotate_joint_button = QPushButton("Rotate Joint Along Axis of Motion")
        self.translate_joint_button = QPushButton("Translate Joint Along Axis of Motion")
        self.delete_joint_button = QPushButton("Delete Joint")

        joint_layout = QVBoxLayout()
        joint_layout.addWidget(self.select_joint_options)
        joint_layout.addWidget(self.transform_joint_button)
        joint_layout.addWidget(self.rotate_joint_button)
        joint_layout.addWidget(self.translate_joint_button)
        joint_layout.addWidget(self.delete_joint_button)
        main_layout = QVBoxLayout()
        main_layout.addLayout(joint_layout)

        button_widget = QWidget()
        button_widget.setLayout(main_layout)
        dock = QDockWidget("Edit Joints", self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        dock.setWidget(button_widget)

        self.transform_joint_button.clicked.connect(self.transform_joint)
        self.rotate_joint_button.clicked.connect(self.rotate_joint)
        self.translate_joint_button.clicked.connect(self.translate_joint)
        self.delete_joint_button.clicked.connect(self.delete_joint)

        self.selected_joint = -1

        self.select_joint_options.currentIndexChanged.connect(self.joint_selection_changed)

    def joint_selection_changed(self, index):
        if index != self.selected_joint:
            self.selected_joint = index
            self.update_joint()

    def add_chain(self, chain):
        self.chain = chain
        self.select_joint_options.blockSignals(True)
        self.select_joint_options.clear() 
    
        for index, joint in enumerate(self.chain.Joints):
            self.select_joint_options.addItem("Joint " + str(index) + " - " + joint.__class__.__name__)
    
        self.select_joint_options.blockSignals(False) 
        self.select_joint_options.setCurrentIndex(self.selected_joint)

    def transform_joint(self):
        dialog = TransformDialog(self) 
        if not self.chain:
            error_dialog = ErrorDialog('Please initialize a chain.')
            error_dialog.exec_()
        if self.selected_joint == -1:
            error_dialog = ErrorDialog('Please select a joint.')
            error_dialog.exec_()
        elif dialog.exec_() == QDialog.Accepted:
            transformation = dialog.get_transform()
            if transformation is not None:
                propagate = dialog.propogateTransformCheckbox.isChecked()

                if self.chain.transformJoint(self.selected_joint, transformation, propagate):
                    self.update_joint()
                    success_dialog = SuccessDialog('Joint successfully transformed!')
                    success_dialog.exec_()
                else:
                    error_dialog = ErrorDialog('Error transforming joint.')
                    error_dialog.exec_()

    def rotate_joint(self):
        dialog = RotationDialog(self)
        if not self.chain:
            error_dialog = ErrorDialog('Please initialize a chain.')
            error_dialog.exec_()
        if self.selected_joint == -1:
            error_dialog = ErrorDialog('Please select a joint.')
            error_dialog.exec_()
        elif dialog.exec_() == QDialog.Accepted:
            angle = dialog.get_angle()
            if angle is not None:
                propogate = dialog.propogateRotationCheckbox.isChecked()
                apply_to_previous_waypoint = dialog.applyToPreviousWaypointRotationCheckbox.isChecked()
                self.selected_joint = self.select_joint_options.currentIndex()

                if self.chain.rotateJointAboutAxisOfMotion(self.selected_joint, math.radians(angle), propogate, apply_to_previous_waypoint):
                    self.update_joint()
                    success_dialog = SuccessDialog('Joint successfully rotated!')
                    success_dialog.exec_()
                else:
                    error_dialog = ErrorDialog('Error rotating joint.')
                    error_dialog.exec_()

    def translate_joint(self):
        dialog = TranslationDialog(self)
        if not self.chain:
            error_dialog = ErrorDialog('Please initialize a chain.')
            error_dialog.exec_()
        if self.selected_joint == -1:
            error_dialog = ErrorDialog('Please select a joint.')
            error_dialog.exec_()
        elif dialog.exec_() == QDialog.Accepted:
            distance = dialog.get_distance()
            if distance is not None:
                propogate = dialog.propogateTranslationCheckbox.isChecked()
                apply_to_previous_waypoint = dialog.applyToPreviousWaypointTranslationCheckbox.isChecked()
                self.selected_joint = self.select_joint_options.currentIndex()

                if self.chain.translateJointAlongAxisOfMotion(self.selected_joint, distance, propogate, apply_to_previous_waypoint):
                    self.update_joint()
                    success_dialog = SuccessDialog('Joint successfully translated!')
                    success_dialog.exec_()
                else:
                    error_dialog = ErrorDialog('Error translating joint.')
                    error_dialog.exec_()

    def delete_joint(self):
        dialog = DeleteDialog(self)
        if not self.chain:
            error_dialog = ErrorDialog('Please initialize a chain.')
            error_dialog.exec_()
        if self.selected_joint == -1:
            error_dialog = ErrorDialog('Please select a joint.')
            error_dialog.exec_()
        elif dialog.exec_() == QDialog.Accepted:
            self.selected_joint = self.select_joint_options.currentIndex() 
            if self.chain.delete(self.selected_joint):
                self.update_joint()
                success_dialog = SuccessDialog('Joint successfully deleted!')
                success_dialog.exec_()
            else:
                error_dialog = ErrorDialog('Error deleting joint.')
                error_dialog.exec_()

    def update_joint(self):
        self.select_joint_options.blockSignals(True)

        self.plot_widget.clear()
        self.select_joint_options.clear()

        grid = gl.GLGridItem()
        self.plot_widget.addItem(grid)
        grid.setColor((0,0,0,255))
        
        if self.chain is not None:
            self.chain.addToWidget(self, selectedJoint=self.selected_joint)

        self.select_joint_options.blockSignals(False)

    def update_plot(self):
        self.plot_widget.clear()
        self.setCentralWidget(self.plot_widget)

        grid = gl.GLGridItem()
        self.plot_widget.addItem(grid)
        grid.setColor((0,0,0,255))

        self.chain.addToWidget(self)

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
    
    def add_joint(self, dialog):
        if dialog.exec_() == QDialog.Accepted:
            joint = dialog.getJoint()
            if (self.chain == None) :
                self.chain = KinematicChain(joint)
            else :
                self.chain.append(joint)
            
            self.update_plot()

    def add_prismatic_func(self):
        dialog = AddPrismaticDialog(self.numSides, self.r)
        self.add_joint(dialog)

    def add_revolute_func(self):
        dialog = AddRevoluteDialog(self.numSides, self.r)
        self.add_joint(dialog)

    def add_waypoint_func(self):
        dialog = AddWaypointDialog(self.numSides, self.r)
        self.add_joint(dialog)

    def add_tip_func(self):
        dialog = AddTipDialog(self.numSides, self.r)
        self.add_joint(dialog)

    def create_new_chain_func(self):
        dialog = CreateNewChainDialog()
        if dialog.exec_() == QDialog.Accepted:
            self.chain = None
            self.numSides = dialog.getNumSides()
            self.r = dialog.getR()
            self.update_plot()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PointEditorWindow()
    window.show()
    sys.exit(app.exec_())
