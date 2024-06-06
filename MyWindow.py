"""
@author: Raymond Feng
"""

import sys
import numpy as np
import pyqtgraph.opengl as gl
from PyQt5 import QtCore as qc
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton, QDockWidget, QComboBox, QHBoxLayout, QLabel, QDialog, QLineEdit, QCheckBox, QMessageBox, QButtonGroup, QRadioButton, QSlider
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap, QSurfaceFormat
from pyqtgraph.Qt import QtCore
import pyqtgraph as pg
from OpenGL.GL import *
from OpenGL.GLU import *
from spatialmath import SE3
import math
from Joint import *
from PathCSC import *
from KinematicChain import *
import re

class EditJointStateDialog(QDialog):

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Edit Joint State')
        self.setGeometry(100, 100, 200, 100)

        layout = QVBoxLayout()

        self.setupStateInput(layout)
        self.setLayout(layout)

    def setupStateInput(self, layout):
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

    def onApplyClicked(self):
        self.accept()

    def get_state(self):
        while True:
            try:
                state = float(self.state_input.text())
                return state
            except ValueError:
                QMessageBox.warning(self, "Invalid Input", "Please enter a valid join state.")
                self.exec_() 
                return None

class TransformDialog(QDialog):
    propogateTransform = True
    relativeTransform = True

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
        self.relativeTransformCheckbox = QCheckBox("Relative")

        self.propogateTransformCheckbox.setChecked(TransformDialog.propogateTransform)
        self.relativeTransformCheckbox.setChecked(TransformDialog.relativeTransform)

        layout.addWidget(self.propogateTransformCheckbox)
        layout.addWidget(self.relativeTransformCheckbox)

        self.setupTransformInputAndButtons(layout)
        self.setLayout(layout)

    def toggleAxisSelection(self):
        isRotateSelected = self.operationSelection.currentText() == "Rotate"
        self.axisSelection.setVisible(isRotateSelected)
        self.axisSelectionLabel.setVisible(isRotateSelected)

    def setupTransformInputAndButtons(self, layout):
        self.transform_input = QLineEdit(self)
        self.transform_input.setPlaceholderText('x, y, z for Translate or degree for Rotate')
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
        TransformDialog.relativeTransform = self.relativeTransformCheckbox.isChecked()
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
    relative = True
    fixedPosition = False
    fixedOrientation = False

    def __init__(self, parent=None):
        super().__init__(parent)
    
    def getJoint(self):
        return self.jointToAdd
        
    def getIsRelative(self):
        return self.relative
    
    def getFixedPosition(self):
        return self.fixedPosition
    
    def getFixedOrientation(self):
        return self.fixedOrientation
        
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
        
    def on_relative_clicked(self):
        self.relative = True

    def on_absolute_clicked(self):
        self.relative = False

    def on_fixed_position_toggled(self):
        self.fixedPosition = not self.fixedPosition

    def on_fixed_orientation_toggled(self):
        self.fixedOrientation = not self.fixedOrientation
        
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

        relative_layout = QHBoxLayout()
        self.radio_group = QButtonGroup()
        self.relative = QRadioButton('Relative')
        self.absolute = QRadioButton('Absolute')
        self.relative.setChecked(True)
        self.relative.clicked.connect(self.on_relative_clicked)
        self.absolute.clicked.connect(self.on_absolute_clicked)
        relative_layout.addWidget(self.relative)
        relative_layout.addWidget(self.absolute)
        layout.addLayout(relative_layout)

        self.position_checkbox = QCheckBox('Fixed Position')
        self.position_checkbox.stateChanged.connect(self.on_fixed_position_toggled)
        layout.addWidget(self.position_checkbox)

        self.orientation_checkbox = QCheckBox('Fixed Orientation')
        self.orientation_checkbox.stateChanged.connect(self.on_fixed_orientation_toggled)
        layout.addWidget(self.orientation_checkbox)

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

        relative_layout = QHBoxLayout()
        self.radio_group = QButtonGroup()
        self.relative = QRadioButton('Relative')
        self.absolute = QRadioButton('Absolute')
        self.relative.setChecked(True)
        self.relative.clicked.connect(self.on_relative_clicked)
        self.absolute.clicked.connect(self.on_absolute_clicked)
        relative_layout.addWidget(self.relative)
        relative_layout.addWidget(self.absolute)
        layout.addLayout(relative_layout)

        self.position_checkbox = QCheckBox('Fixed Position')
        self.position_checkbox.stateChanged.connect(self.on_fixed_position_toggled)
        layout.addWidget(self.position_checkbox)

        self.orientation_checkbox = QCheckBox('Fixed Orientation')
        self.orientation_checkbox.stateChanged.connect(self.on_fixed_orientation_toggled)
        layout.addWidget(self.orientation_checkbox)

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

        relative_layout = QHBoxLayout()
        self.radio_group = QButtonGroup()
        self.relative = QRadioButton('Relative')
        self.absolute = QRadioButton('Absolute')
        self.relative.setChecked(True)
        self.relative.clicked.connect(self.on_relative_clicked)
        self.absolute.clicked.connect(self.on_absolute_clicked)
        relative_layout.addWidget(self.relative)
        relative_layout.addWidget(self.absolute)
        layout.addLayout(relative_layout)

        self.position_checkbox = QCheckBox('Fixed Position')
        self.position_checkbox.stateChanged.connect(self.on_fixed_position_toggled)
        layout.addWidget(self.position_checkbox)

        self.orientation_checkbox = QCheckBox('Fixed Orientation')
        self.orientation_checkbox.stateChanged.connect(self.on_fixed_orientation_toggled)
        layout.addWidget(self.orientation_checkbox)
        
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

        relative_layout = QHBoxLayout()
        self.radio_group = QButtonGroup()
        self.relative = QRadioButton('Relative')
        self.absolute = QRadioButton('Absolute')
        self.relative.setChecked(True)
        self.relative.clicked.connect(self.on_relative_clicked)
        self.absolute.clicked.connect(self.on_absolute_clicked)
        relative_layout.addWidget(self.relative)
        relative_layout.addWidget(self.absolute)
        layout.addLayout(relative_layout)

        self.position_checkbox = QCheckBox('Fixed Position')
        self.position_checkbox.stateChanged.connect(self.on_fixed_position_toggled)
        layout.addWidget(self.position_checkbox)

        self.orientation_checkbox = QCheckBox('Fixed Orientation')
        self.orientation_checkbox.stateChanged.connect(self.on_fixed_orientation_toggled)
        layout.addWidget(self.orientation_checkbox)
        
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

        self.bounding_balls = []
        self.radius = 1

        dist = self.opts['distance']
        self.near_clip = dist * 0.001
        self.far_clip = dist * 1000.

    click_signal = qc.pyqtSignal(int)
    click_signal_arrow = qc.pyqtSignal(int)

    selected_index = -1

    def toggle_lock(self):
        self.locked = not self.locked
        print("Screen lock toggled:", "Locked" if self.locked else "Unlocked")

    def mousePressEvent(self, event):
        if not self.locked:
            super(ClickableGLViewWidget, self).mousePressEvent(event)

    def mouseMoveEvent(self, event):
        if not self.locked:
            super(ClickableGLViewWidget, self).mouseMoveEvent(event)

    def mouseReleaseEvent(self, event):
        if self.locked:
            lpos = event.position() if hasattr(event, 'position') else event.localPos()
            region = [lpos.x()-5, lpos.y()-5, 10, 10]
            # itemsAt seems to take in device pixels
            dpr = self.devicePixelRatioF()
            region = tuple([x * dpr for x in region])

            arrowIndex = -1
            jointSeen = False   # Basically just selects the first joint in distance order.

            joints = []
            arrows = []

            for item in self.itemsAt(region):
                if (item.objectName() == "Arrow"):
                    arrows.append(item)

                if (item.objectName() == "Joint" or item.objectName() == "Waypoint"):
                    joints.append(item)

            if (len(arrows) > 0 and self.selectedIndex != -1):
                arrowIndex = arrows[0].id
            else:
                if(len(joints) > 0):
                    self.selectedIndex = joints[0].id
                else:
                    self.selectedIndex = -1
            
            self.click_signal.emit(self.selectedIndex)
            self.click_signal_arrow.emit(arrowIndex)
        if not self.locked:
            super(ClickableGLViewWidget, self).mouseReleaseEvent(event)
    
    def get_ray(self, x_coord: int, y_coord: int) -> tuple[np.ndarray, np.ndarray]:
        """
        Method returns the ray origin (current camera position) and ray unit vector for
        selection of triangulated meshes in the GLViewWidget.

        :param x_coord: Mouse click local x coordinate within the widget.
        :param y_coord: Mouse click local y coordinate within the widget.

        Note: Mouse click coordinate system origin is top left of the GLViewWidget.

        from @gordankos on github
        https://github.com/pyqtgraph/pyqtgraph/issues/2647
        """
        x0, y0, width, height = self.getViewport()
        ray_origin = np.array(self.cameraPosition())

        projection_matrix = np.array(self.projectionMatrix().data()).reshape(4, 4)
        view_matrix = np.array(self.viewMatrix().data()).reshape(4, 4)
        view_matrix = np.transpose(view_matrix)

        ndc_x = (4.0 * x_coord / width) - 1.0                    
        ndc_y = (4.0 * y_coord) / height - 1.0                

        clip_coords = np.array([ndc_x, ndc_y, -1.0, 1.0])

        p = np.linalg.inv(view_matrix) @ np.linalg.inv(projection_matrix) @ clip_coords

        eye_coords = np.linalg.inv(projection_matrix) @ clip_coords
        eye_coords /= eye_coords[3]
        eye_coords = np.array([eye_coords[0], eye_coords[1], -1.0, 0.0])

        ray_direction = np.linalg.inv(view_matrix) @ eye_coords
        ray_direction = ray_direction[:3] / np.linalg.norm(ray_direction[:3])

        return ray_origin, ray_direction

    def project_click(self, pos, s, r):
        """
        pos: (local X coord within widget, local Y coord within widget)
        - The top left corner is the origin point

        s: center of the sphere that is being tested
        r: radius of the sphere that is being tested
        """
        o, d = self.get_ray(pos.x(), pos.y())

        print("origin:")
        print(o)
        print("direction:")
        print(d)

        a = d[0] ** 2 + d[1] ** 2 + d[2] ** 2
        b = 2 * (d[0] * (o[0] - s[0]) + d[1] * (o[1] - s[1]) + d[2] * (o[2] - s[2]))
        c = (o[0] - s[0]) ** 2 + (o[1] - s[1]) ** 2 + (o[2] - s[2]) ** 2 - r * r

        discrim = b * b - 4 * a * c

        if (discrim < 0): 
            root_discrim = math.sqrt(discrim)

        t = 0

        t0 = (-b - root_discrim) / (2 * a)
        if (t0 < 0) :
            t1 = (-b + root_discrim) / (2 * a)
            t = t1
        else:
            t = t0

        p = o + t * d
        return t
 
class PointEditorWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Point Editor")
        self.setGeometry(100, 100, 800, 600)

        self.plot_widget = ClickableGLViewWidget()
        self.setCentralWidget(self.plot_widget)
        self.plot_widget.setBackgroundColor(255,255,255, 255)

        self.toggleButton = QPushButton('Toggle Lock', self)
        self.toggleButton.clicked.connect(self.plot_widget.toggle_lock)
        self.toggleButton.setGeometry(20, 20, 140, 40)
        self.toggleButton.setStyleSheet('background-color: black; color: white;')

        grid = gl.GLGridItem()

        self.plot_widget.addItem(grid)
        grid.setColor((0,0,0,255))

        self.current_point = 0
        self.chain = None
        self.chain_created = False

        self.numSides = 4
        self.r = 1
        self.plot_widget.radius = 1
        self.crease_pattern = None
        self.selected_joint = -1
        self.selected_arrow = -1
        self.selected_axis_name = 'N/A'

        self.plot_widget.click_signal.connect(self.joint_selection_changed)
        self.plot_widget.click_signal_arrow.connect(self.arrow_selection_changed)

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

        add_joints_dock = QDockWidget("Add Joints", self)
        self.add_joints_widget = QWidget()
        self.add_joints_widget.setLayout(add_joints_layout)
        add_joints_dock.setWidget(self.add_joints_widget)

        self.addDockWidget(Qt.RightDockWidgetArea, add_joints_dock)

        # //////////////////////////////////    AXIS KEY    ////////////////////////////////////
        axis_key_layout = QVBoxLayout()
        self.axis_key_widget = QWidget()
        self.axis_key_widget.setLayout(axis_key_layout)

        self.x_axis_widget = self.create_axis_label('x', Qt.red)
        self.y_axis_widget = self.create_axis_label('y', Qt.green)
        self.z_axis_widget = self.create_axis_label('z', Qt.blue)

        axis_key_layout.addWidget(self.x_axis_widget)
        axis_key_layout.addWidget(self.y_axis_widget)
        axis_key_layout.addWidget(self.z_axis_widget)

        axis_key_dock = QDockWidget("Axis Key", self)
        axis_key_dock.setWidget(self.axis_key_widget)
        self.addDockWidget(Qt.RightDockWidgetArea, axis_key_dock)

        # ////////////////////////////////    EDIT JOINTS    ///////////////////////////////////
        self.select_joint_options = QComboBox()
        self.delete_joint_button = QPushButton("Delete Joint")
        self.edit_joint_state_button = QPushButton("Edit Joint State")
        self.current_state_label = QLabel('Min State ≤ Current State ≤ Max State')

        joint_layout = QVBoxLayout()
        joint_layout.addWidget(self.select_joint_options)
        joint_layout.addWidget(self.delete_joint_button)
        joint_layout.addWidget(self.edit_joint_state_button)
        joint_layout.addWidget(self.current_state_label)

        main_layout = QVBoxLayout()
        main_layout.addLayout(joint_layout)

        button_widget = QWidget()
        button_widget.setLayout(main_layout)
        dock = QDockWidget("Edit Joints", self)
        self.addDockWidget(Qt.RightDockWidgetArea, dock)
        dock.setWidget(button_widget)

        self.delete_joint_button.clicked.connect(self.delete_joint)
        self.edit_joint_state_button.clicked.connect(self.edit_joint_state)

        self.select_joint_options.currentIndexChanged.connect(self.joint_selection_changed)

        self.rotationLabel = QLabel("Rotation Angle: 0°", self)
        self.rotationSlider = QSlider(Qt.Horizontal, self)
        self.rotationSlider.setMinimum(-360)
        self.rotationSlider.setMaximum(360)
        self.rotationSlider.setValue(0)
        self.rotationSlider.valueChanged.connect(self.adjust_rotation)

        self.translate_label = QLabel('Transform: 0', self)
        self.translate_slider = QSlider(Qt.Horizontal)
        self.translate_slider.setMinimum(-100)
        self.translate_slider.setMaximum(100)
        self.translate_slider.setValue(0)
        self.translate_slider.valueChanged.connect(self.adjust_translation)

        sliderLayout = QVBoxLayout()
        sliderLayout.addWidget(self.rotationLabel) 
        sliderLayout.addWidget(self.rotationSlider) 
        sliderLayout.addWidget(self.translate_label)
        sliderLayout.addWidget(self.translate_slider)

        sliderWidget = QWidget(self)
        sliderWidget.setLayout(sliderLayout)
        sliderWidget.setGeometry(20, 70, 200, 50)  

        slider_joints_dock = QDockWidget("Edit Joints (Sliders)", self)
        self.slider_joints_widget = QWidget()
        self.slider_joints_widget.setLayout(sliderLayout)
        slider_joints_dock.setWidget(self.slider_joints_widget)

        self.addDockWidget(Qt.RightDockWidgetArea, slider_joints_dock)
        self.oldRotVal = 0
        self.oldTransVal = 0
        self.rotationSlider.setDisabled(True)

        # ////////////////////////////////    CREASE PATTERN   ///////////////////////////////////
        crease_dock_layout = QVBoxLayout()
        self.crease_pattern_name_input = QLineEdit()
        self.crease_pattern_name_input.setPlaceholderText('Name')
        crease_dock_layout.addWidget(self.crease_pattern_name_input)

        self.save_crease_pattern_button = QPushButton('Save Crease Pattern')
        self.save_crease_pattern_button.clicked.connect(self.save_crease_pattern)  
        crease_dock_layout.addWidget(self.save_crease_pattern_button)  

        crease_button_widget = QWidget()
        crease_button_widget.setLayout(crease_dock_layout) 

        crease_dock = QDockWidget("Crease Pattern", self)
        crease_dock.setWidget(crease_button_widget)
        self.addDockWidget(Qt.RightDockWidgetArea, crease_dock)

    def update_transform_label(self):


    def save_crease_pattern(self):
        crease_pattern_name = self.crease_pattern_name_input.text()
        if self.crease_pattern == None:
            self.crease_pattern = self.chain.creasePattern()
        if len(crease_pattern_name) > 0:
            self.crease_pattern.show(dxfName=crease_pattern_name)
        else:
            self.crease_pattern.show()

    @QtCore.pyqtSlot(int)
    def joint_selection_changed(self, index):
        if index != self.selected_joint:
            self.selected_joint = index
            self.selected_arrow = -1
            self.update_joint()

    @QtCore.pyqtSlot(int)
    def arrow_selection_changed(self, index):
        if index != self.selected_arrow:
            self.selected_arrow = index
            self.update_joint()
            if self.selected_arrow != 1:
                self.rotationSlider.setDisabled(False)
                self.transform_slider
            else:
                self.rotationSlider.setDisabled(True)

    def add_chain(self, chain):
        self.chain = chain
        self.select_joint_options.blockSignals(True)
        self.select_joint_options.clear() 
    
        for joint in self.chain.Joints:
            self.select_joint_options.addItem("Joint " + str(joint.id) + " - " + joint.__class__.__name__)
    
        self.select_joint_options.blockSignals(False) 
        self.select_joint_options.setCurrentIndex(self.selected_joint)

    def edit_joint_state(self):
        dialog = EditJointStateDialog(self) 
        if not self.chain:
            error_dialog = ErrorDialog('Please initialize a chain.')
            error_dialog.exec_()
        if self.selected_joint == -1:
            error_dialog = ErrorDialog('Please select a joint.')
            error_dialog.exec_()
        elif dialog.exec_() == QDialog.Accepted:
            edit = dialog.get_state()
            if edit is not None:

                if self.chain.setJointState(self.selected_joint, edit):
                    self.update_joint()
                    success_dialog = SuccessDialog('Joint state successfully edited!')
                    success_dialog.exec_()
                    min = self.chain.Joints[self.selected_joint].stateRange()[0]
                    max = self.chain.Joints[self.selected_joint].stateRange()[1]
                    current = self.chain.Joints[self.selected_joint].state
                    self.current_state_label.setText(f"{min} ≤ {current} ≤ {max}")
                else:
                    error_dialog = ErrorDialog('Error editing joint state.')
                    error_dialog.exec_()

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
                relative = dialog.relativeTransformCheckbox.isChecked()

                if self.chain.transformJoint(self.selected_joint, transformation, propagate, relative=relative):
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

    def adjust_rotation(self, value):
        self.rotationLabel.setText(f"Rotation Angle: {value}°")
        angle_radians = math.radians(value - self.oldRotVal)
        self.oldRotVal = value
        if self.chain and self.selected_joint != -1:
            transformation = SE3.Rz(angle_radians)
            if self.chain.transformJoint(self.selected_joint, transformation):
                self.update_joint()

    def adjust_translation(self, value):
        self.translate_label.setText(f'Transform: {float(value) / 10}')
        actualVal = float(value) / 10
        amount = actualVal - self.oldTransVal
        self.oldTransVal = actualVal
        if self.chain and self.selected_joint != -1:
            transformation = SE3()
            if (self.selected_arrow == 0):
                transformation = SE3.Tx(amount)
            if (self.selected_arrow == 1):
                transformation = SE3.Ty(amount)
            if (self.selected_arrow == 2):
                transformation = SE3.Tz(amount)
    
            if self.chain.transformJoint(self.selected_joint, transformation, relative=True):
                self.update_joint()

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

        for joint in self.chain.Joints:
            if (joint.id == self.selected_joint):
                joint.addArrows(self, selectedArrow=self.selected_arrow)

    def update_plot(self):
        self.plot_widget.clear()
        self.setCentralWidget(self.plot_widget)

        grid = gl.GLGridItem()
        self.plot_widget.addItem(grid)
        grid.setColor((0,0,0,255))

        for index, joint in enumerate(self.chain.Joints):
            joint.id = index

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
            joint : Joint = dialog.getJoint()

            if (self.chain == None):
                joint.id = 0
            else:
                joint.id = len(self.chain.Joints)
                
            isRelative = dialog.getIsRelative()
            isFixedPosition = dialog.getFixedPosition()
            isFixedOrientation = dialog.getFixedOrientation()
            isSafe = True

            if (isFixedPosition or isFixedOrientation):
                isSafe = False

            if (self.chain == None) :
                self.chain = KinematicChain(joint)
            else :
                self.chain.append(joint, relative=isRelative, fixedPosition=isFixedPosition, fixedOrientation=isFixedOrientation, safe=isSafe)
            
            self.update_plot()

    def chain_not_created(self):
        error_dialog = ErrorDialog('Please create a chain first.')
        error_dialog.exec_()
        self.create_new_chain_func()

    def add_prismatic_func(self):
        if (not self.chain_created) :
            self.chain_not_created()
        else:
            dialog = AddPrismaticDialog(self.numSides, self.r)
            self.add_joint(dialog)

    def add_revolute_func(self):
        if (not self.chain_created) :
            self.chain_not_created()
        else:
            dialog = AddRevoluteDialog(self.numSides, self.r)
            self.add_joint(dialog)

    def add_waypoint_func(self):
        if (not self.chain_created) :
            self.chain_not_created()
        else:
            dialog = AddWaypointDialog(self.numSides, self.r)
            self.add_joint(dialog)

    def add_tip_func(self):
        if (not self.chain_created) :
            self.chain_not_created()
        else:
            dialog = AddTipDialog(self.numSides, self.r)
            self.add_joint(dialog)

    def create_new_chain_func(self):
        dialog = CreateNewChainDialog()
        if dialog.exec_() == QDialog.Accepted:

            # clears everything
            self.chain = None
            self.chain_created = True
            self.numSides = dialog.getNumSides()
            self.selected_joint = -1
            self.r = dialog.getR()
            self.plot_widget.clear()
            self.plot_widget.radius = self.r
            self.setCentralWidget(self.plot_widget)

            grid = gl.GLGridItem()
            self.plot_widget.addItem(grid)
            grid.setColor((0,0,0,255))

            success_dialog = SuccessDialog('Chain created!')
            success_dialog.exec_()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PointEditorWindow()
    window.show()
    sys.exit(app.exec_())