"""
@author: Raymond Feng and Andy Wang
"""

import sys
import numpy as np
import pyqtgraph.opengl as gl
from PyQt5 import QtCore as qc
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton, QDockWidget, QComboBox, QHBoxLayout, QLabel, QDialog, QLineEdit, QCheckBox, QMessageBox, QButtonGroup, QRadioButton, QSlider, QSizePolicy
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QPixmap, QSurfaceFormat, QKeyEvent, QPixmap, QIcon
from pyqtgraph.Qt import QtCore
import pyqtgraph as pg
from OpenGL.GL import *
from OpenGL.GLU import *
from spatialmath import SE3
import math
from PathCSC import *
from KinematicChain import *
import re
from scipy.spatial.transform import Rotation as R
from testqtgraph import *

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
            QMessageBox.warning(self, "Invalid Input", "Please enter a valid join state.")
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
        angle_label = QLabel("Cone Angle (degrees):")
        self.angle_input = QLineEdit()
        angle_layout.addWidget(angle_label)
        angle_layout.addWidget(self.angle_input)
        layout.addLayout(angle_layout)

        apply_button = QPushButton('Add')
        apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(apply_button)

        self.setLayout(layout)

        self.numSides = numSides
        self.r = r

    def onApplyClicked(self):
        try:
            neutralLength = float(self.length_input.text())
            numLayers = int(self.numLayers_input.text())
            coneAngleText = float(self.angle_input.text())

            self.jointToAdd = PrismaticJoint(self.numSides, self.r, neutralLength, numLayers, math.radians(coneAngleText), SE3())
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
        angle_label = QLabel("Total Bending Angle (degrees):")
        self.angle_input = QLineEdit()
        angle_layout.addWidget(angle_label)
        angle_layout.addWidget(self.angle_input)
        layout.addLayout(angle_layout)

        apply_button = QPushButton('Add')
        apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(apply_button)

        self.setLayout(layout)

        self.numSides = numSides
        self.r = r

    def onApplyClicked(self):
        bendingAngleText = float(self.angle_input.text())

        self.jointToAdd = RevoluteJoint(self.numSides, self.r, math.radians(bendingAngleText), SE3())
        self.accept()

class AddTipDialog(AddJointDialog):
    isStart = True

    def __init__(self, numSides, r, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Add new joint')
        self.setGeometry(100, 100, 300, 100)

        self.numSides = numSides
        self.r = r
        self.isStart = True

        layout = QVBoxLayout()

        length_layout = QHBoxLayout()
        length_label = QLabel("Length:")
        self.length_input = QLineEdit()
        length_layout.addWidget(length_label)
        length_layout.addWidget(self.length_input)
        layout.addLayout(length_layout)

        radio_layout = QHBoxLayout()
        self.radio_start = QRadioButton('Start Tip')
        self.radio_end = QRadioButton('End Tip')
        self.radio_start.setChecked(True)
        self.radio_start.toggled.connect(self.updateVariable)
        self.radio_end.toggled.connect(self.updateVariable)
        radio_layout.addWidget(self.radio_start)
        radio_layout.addWidget(self.radio_end)
        layout.addLayout(radio_layout)
        
        apply_button = QPushButton('Add')
        apply_button.clicked.connect(self.onApplyClicked)
        layout.addWidget(apply_button)

        self.setLayout(layout)

    def onApplyClicked(self):
        try:
            neutralLength = float(self.length_input.text())

            if (self.isStart):
                self.jointToAdd = StartTip(self.numSides, self.r, SE3(), length=neutralLength)
            else:
                self.jointToAdd = EndTip(self.numSides, self.r, SE3(), length=neutralLength)

            self.accept()
        except ValueError:
            error_dialog = ErrorDialog('Please enter valid integers.')
            error_dialog.exec_()

    def updateVariable(self):
        if self.radio_start.isChecked():
            self.isStart = True
        elif self.radio_end.isChecked():
            self.isStart = False

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
    
class ImageRadioButton(QRadioButton):
    def __init__(self, unchecked_img, checked_img, tooltip_text, parent=None):
        super().__init__(parent)
        self.unchecked_img = QPixmap(unchecked_img)
        self.checked_img = QPixmap(checked_img)
        self.setIconSize(self.unchecked_img.size())
        self.update_icon()

        # Hide the default radio button indicator
        self.setStyleSheet("QRadioButton::indicator { width: 0px; height: 0px; }")

        # Connect the toggled signal to update the icon when the state changes
        self.toggled.connect(self.update_icon)

        self.setToolTip(tooltip_text)

    def update_icon(self):
        if self.isChecked():
            self.setIcon(QIcon(self.checked_img))
        else:
            self.setIcon(QIcon(self.unchecked_img))   

class ClickableGLViewWidget(gl.GLViewWidget):
    def __init__(self, parent=None):
        super(ClickableGLViewWidget, self).__init__(parent)
        fmt = QSurfaceFormat()
        fmt.setDepthBufferSize(24)
        # fmt.setVersion(3, 3)
        fmt.setProfile(QSurfaceFormat.CompatibilityProfile)
        QSurfaceFormat.setDefaultFormat(fmt)
        self.setFormat(fmt)
        self.locked = False
        self.mesh = None 
        self.is_dragging = False
        self.drag_start_pos = None

        self.bounding_balls = []
        self.radius = 1

        dist = self.opts['distance']
        self.near_clip = dist * 0.001
        self.far_clip = dist * 1000.
    
    lock_status_changed = pyqtSignal(bool)
    click_signal = qc.pyqtSignal(int)
    click_signal_arrow = qc.pyqtSignal(int)
    drag_change_position = qc.pyqtSignal(np.ndarray)
    drag_change_rotation = qc.pyqtSignal(float)
    key_pressed = qc.pyqtSignal(str)

    selected_index = -1
    selected_arrow = None

    camera_type = "Rotate"

    def toggle_lock(self):
        self.locked = not self.locked
        self.lock_status_changed.emit(self.locked)
        print("Screen lock toggled:", "Locked" if self.locked else "Unlocked")

    def mousePressEvent(self, event):
        if (event.buttons() and Qt.LeftButton and event.buttons() != QtCore.Qt.MouseButton.MiddleButton):
            self.is_dragging = False
            self.drag_start_pos = event.pos()

            lpos = event.position() if hasattr(event, 'position') else event.localPos()
            region = [lpos.x()-5, lpos.y()-5, 10, 10]
            # itemsAt seems to take in device pixels
            dpr = self.devicePixelRatioF()
            region = tuple([x * dpr for x in region])

            joints = []
            arrows = []

            for item in self.itemsAt(region):
                if (item.objectName() == "Arrow"):
                    arrows.append(item)

                if (item.objectName() == "Joint" or item.objectName() == "Waypoint"):
                    joints.append(item)

            if (len(arrows) > 0):
                self.selected_arrow = arrows[0]
                self.click_signal_arrow.emit(self.selected_arrow.id)
            else: 
                self.selected_arrow = None

    def mouseMoveEvent(self, event):
        if (event.buttons() and (Qt.LeftButton or Qt.MiddleButton)) and (event.pos() - self.drag_start_pos).manhattanLength() >= QApplication.startDragDistance():
            self.is_dragging = True

        if (self.selected_arrow):
            lpos = event.position() if hasattr(event, 'position') else event.localPos()
            region = [lpos.x()-10, lpos.y()-10, 20, 20]
            dpr = self.devicePixelRatioF()
            region = tuple([x * dpr for x in region])

            selected_translate_point = None
            selected_rotate_point = None

            for item in self.itemsAt(region):
                if (item.objectName() == "line_sphere"):
                    selected_translate_point = item
                    break
                if (item.objectName() == "rotate_sphere"):
                    selected_rotate_point = item
                    break
            
            if (selected_translate_point):
                self.drag_change_position.emit(selected_translate_point.position)
            if (selected_rotate_point):
                self.drag_change_rotation.emit(selected_rotate_point.rotation)

        else:
            if (self.is_dragging):
                lpos = event.position() if hasattr(event, 'position') else event.localPos()
                if not hasattr(self, 'mousePos'):
                    self.mousePos = lpos
                diff = lpos - self.mousePos
                self.mousePos = lpos

                if event.buttons() == QtCore.Qt.MouseButton.MiddleButton:
                    self.pan(diff.x(), diff.y(), 0, relative='view')
                elif event.buttons() == QtCore.Qt.MouseButton.LeftButton:
                    if (self.camera_type == "Rotate"):
                        self.orbit(-diff.x(), diff.y())
                    elif (self.camera_type == "Pan"):
                        self.pan(diff.x(), diff.y(), 0, relative='view')

    def mouseReleaseEvent(self, event):
        if not self.is_dragging:
            lpos = event.position() if hasattr(event, 'position') else event.localPos()
            region = [lpos.x()-5, lpos.y()-5, 10, 10]
            dpr = self.devicePixelRatioF()
            region = tuple([x * dpr for x in region])

            arrow_index = -1

            joints = []
            arrows = []

            for item in self.itemsAt(region):
                if (item.objectName() == "Arrow"):
                    arrows.append(item)

                if (item.objectName() == "Joint" or item.objectName() == "Waypoint"):
                    joints.append(item)

                if (item.objectName() == "Link"):
                    print("link clicked")

            if (len(arrows) > 0 and self.selected_index != -1):
                arrow_index = arrows[0].id
            else:
                if(len(joints) > 0):
                    self.selected_index = joints[0].id
                else:
                    self.selected_index = -1
            
            self.click_signal.emit(self.selected_index)
            self.click_signal_arrow.emit(arrow_index)

    def keyPressEvent(self, event: QKeyEvent):
        if event.key() == Qt.Key_W:
            self.key_pressed.emit("Translate")
        elif event.key() == Qt.Key_E:
            self.key_pressed.emit("Rotate")
        elif event.key() == Qt.Key_Return or event.key() == Qt.Key_Enter:
            self.key_pressed.emit("Enter")
        elif event.key() == Qt.Key_Escape:
            self.key_pressed.emit("Escape")
        elif event.key() == Qt.Key_Delete or event.key() == Qt.Key_Backspace:
            self.key_pressed.emit("Delete")
        elif event.key() == Qt.Key_X:
            self.key_pressed.emit("X")
        elif event.key() == Qt.Key_Y:
            self.key_pressed.emit("Y")
        elif event.key() == Qt.Key_Z:
            self.key_pressed.emit("Z")
 
class PointEditorWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Point Editor")
        self.setGeometry(100, 100, 800, 600)

        self.plot_widget = ClickableGLViewWidget()
        self.setCentralWidget(self.plot_widget)
        self.plot_widget.setBackgroundColor(255,255,255,255)

        grid = gl.GLGridItem()

        self.plot_widget.addItem(grid)
        grid.setColor((0,0,0,255))

        self.current_point = 0
        self.chain = None
        self.chain_created = False
        self.stl_generated = False

        self.numSides = 4
        self.r = 1
        self.plot_widget.radius = 1

        self.crease_pattern = None
        self.control_type = "Translate"
        self.is_local = True

        self.selected_joint = -1
        self.selected_arrow = -1
        self.selected_axis_name = 'N/A'
        self.last_joint = -1
        self.selected_frame = -1

        self.plot_widget.click_signal.connect(self.joint_selection_changed)
        self.plot_widget.click_signal_arrow.connect(self.arrow_selection_changed)
        self.plot_widget.drag_change_position.connect(self.drag_translate)
        self.plot_widget.drag_change_rotation.connect(self.drag_rotate)
        self.plot_widget.key_pressed.connect(self.key_pressed)

        # //////////////////////////////////    ADD TOP DOCK    ///////////////////////////////////
        top_dock_widget = QDockWidget("Top Dock", self)
        top_dock_widget.setAllowedAreas(Qt.TopDockWidgetArea)

        self.key_bar = QWidget()
        self.key_bar_layout = QHBoxLayout(self.key_bar)
        self.key_bar.setLayout(self.key_bar_layout)
        self.key_bar.setFixedHeight(40)
        self.init_key_bar()

        top_dock_widget.setWidget(self.key_bar)

        self.addDockWidget(Qt.TopDockWidgetArea, top_dock_widget)

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

        # ////////////////////////////////    SLIDERS    ///////////////////////////////////
        slider_joints_dock = QDockWidget("Edit Joints (Sliders)", self)
        slider_joints_widget = QWidget() 
        mainLayout = QVBoxLayout(slider_joints_widget)  

        self.rotationLabel = QLabel("Rotate N/A Axis: 0°")
        self.rotationSlider = QSlider(Qt.Horizontal)
        self.rotationSlider.setMinimum(-360)
        self.rotationSlider.setMaximum(360)
        self.rotationSlider.setValue(0)
        self.rotationSlider.setDisabled(True) 
        self.rotationSlider.valueChanged.connect(self.adjust_rotation)

        self.translate_label = QLabel('Translate N/A Axis: 0', self)
        self.translate_slider = QSlider(Qt.Horizontal, self)
        self.translate_slider.setMinimum(-100)
        self.translate_slider.setMaximum(100)
        self.translate_slider.setValue(0)
        self.translate_slider.valueChanged.connect(self.adjust_translation)

        self.state_label = QLabel('Edit Joint N/A State: 0', self)
        self.state_slider = QSlider(Qt.Horizontal, self)
        self.state_slider.setMinimum(-100)
        self.state_slider.setMaximum(100)
        self.state_slider.setValue(0)
        self.state_slider.valueChanged.connect(self.adjust_state)

        self.radius_label = QLabel('Edit Joint Radius: 0', self)
        self.radius_slider = QSlider(Qt.Horizontal, self)
        self.radius_slider.setMinimum(1)
        self.radius_slider.setMaximum(10)
        self.radius_slider.setValue(1)
        self.radius_slider.valueChanged.connect(self.adjust_radius)

        rotationLayout = QHBoxLayout()
        self.rotationInput = QLineEdit(self)
        self.rotationInput.setPlaceholderText("Enter angle in degrees")
        mainLayout.addWidget(self.rotationLabel)
        rotationLayout.addWidget(self.rotationSlider)
        rotationLayout.addWidget(self.rotationInput)  

        translationLayout = QHBoxLayout()
        self.translationInput = QLineEdit(self)
        self.translationInput.setPlaceholderText("Enter distance")
        mainLayout.addWidget(self.translate_label)
        translationLayout.addWidget(self.translate_slider)
        translationLayout.addWidget(self.translationInput) 

        stateLayout = QHBoxLayout()
        self.stateInput = QLineEdit(self)
        self.stateInput.setPlaceholderText("Enter state")
        mainLayout.addWidget(self.state_label)
        stateLayout.addWidget(self.state_slider)
        stateLayout.addWidget(self.stateInput) 

        radiusLayout = QHBoxLayout()
        self.radiusInput = QLineEdit(self)
        self.radiusInput.setPlaceholderText("Enter radius")
        mainLayout.addWidget(self.radius_label)
        radiusLayout.addWidget(self.radius_slider)
        radiusLayout.addWidget(self.radiusInput) 

        self.rotationInput.textChanged.connect(self.adjust_rotation)
        self.translationInput.textChanged.connect(self.adjust_translation)
        self.stateInput.textChanged.connect(self.adjust_state)
        self.radiusInput.textChanged.connect(self.adjust_radius)

        mainLayout.addLayout(rotationLayout)
        mainLayout.addLayout(translationLayout)
        mainLayout.addLayout(stateLayout)
        mainLayout.addLayout(radiusLayout)

        checkboxLayout = QHBoxLayout() 
        self.propogateSliderCheckbox = QCheckBox("Propagate")
        self.relativeSliderCheckbox = QCheckBox("Relative")
        self.propogateSliderCheckbox.setChecked(True)
        self.relativeSliderCheckbox.setChecked(True)
        checkboxLayout.addWidget(self.propogateSliderCheckbox)
        checkboxLayout.addWidget(self.relativeSliderCheckbox)

        self.relativeSliderCheckbox.stateChanged.connect(self.relative_clicked)

        mainLayout.addLayout(checkboxLayout)

        slider_joints_dock.setWidget(slider_joints_widget)
        self.addDockWidget(Qt.RightDockWidgetArea, slider_joints_dock)

        self.oldRotVal = 0
        self.oldTransVal = 0
        self.oldStateVal = 0
        self.oldRadiusVal = 1
        self.rotationSlider.setDisabled(True)
        self.translate_slider.setDisabled(True)
        self.state_slider.setDisabled(True)
        self.radius_slider.setDisabled(True)
        self.rotationInput.setDisabled(True)
        self.translationInput.setDisabled(True)
        self.stateInput.setDisabled(True)
        self.radiusInput.setDisabled(True)

        # ////////////////////////////////    CREASE PATTERN   ///////////////////////////////////
        crease_dock = QDockWidget("Save Crease Pattern", self)
        crease_dock.setAllowedAreas(Qt.RightDockWidgetArea)

        crease_dock_widget = QWidget()
        crease_dock_layout = QHBoxLayout(crease_dock_widget)

        self.crease_pattern_name_input = QLineEdit()
        self.crease_pattern_name_input.setPlaceholderText('Name')
        crease_dock_layout.addWidget(self.crease_pattern_name_input)

        self.save_crease_pattern_button = QPushButton('Save Crease Pattern')
        self.save_crease_pattern_button.clicked.connect(self.save_crease_pattern)  
        crease_dock_layout.addWidget(self.save_crease_pattern_button)  

        crease_dock_widget.setLayout(crease_dock_layout)
        crease_dock.setWidget(crease_dock_widget)
        self.addDockWidget(Qt.RightDockWidgetArea, crease_dock)

        # ////////////////////////////////    WIDGETS DOCK    ///////////////////////////////////
        self.controls_dock = QDockWidget("Control options", self)
        self.controls_options_widget = QWidget()
        self.controls_layout = QVBoxLayout()

        self.control1 = ImageRadioButton("ui/move_unchecked.png", "ui/move_checked.png", "Translate Joint")
        self.control2 = ImageRadioButton("ui/rotate_unchecked.png", "ui/rotate_checked.png", "Rotate Joint")
        self.control1.setChecked(True)

        self.control1.toggled.connect(self.change_control_type)
        self.control2.toggled.connect(self.change_control_type)
    
        self.controls_layout.addWidget(self.control1)
        self.controls_layout.addWidget(self.control2)
        
        self.controls_options_widget.setLayout(self.controls_layout)
        self.controls_dock.setWidget(self.controls_options_widget)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.controls_dock)

        # ////////////////////////////////    WIDGETS DOCK    ///////////////////////////////////
        self.camera_dock = QDockWidget("Camera controls", self)
        self.camera_options_widget = QWidget()
        self.camera_layout = QVBoxLayout()

        self.camera1 = QRadioButton("Rotate Camera")
        self.camera2 = QRadioButton("Pan Camera")
        self.camera1.setChecked(True)

        self.camera1.toggled.connect(self.change_camera_type)
        self.camera2.toggled.connect(self.change_camera_type)
    
        self.camera_layout.addWidget(self.camera1)
        self.camera_layout.addWidget(self.camera2)
        
        self.camera_options_widget.setLayout(self.camera_layout)
        self.camera_dock.setWidget(self.camera_options_widget)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.camera_dock)

        # ////////////////////////////////    SELECT JOINT AS FRAME    ///////////////////////////////////
        self.frame_dock = QDockWidget("Control options", self)
        self.frame_widget = QWidget()
        self.frame_layout = QVBoxLayout()

        self.set_frame_joint = QPushButton("Set Frame")
        self.set_frame_joint.clicked.connect(self.set_joint_as_frame)
        #self.frame_layout.addWidget(self.set_frame_joint)

        frame_string = None
        if (self.selected_frame == -1):
            frame_string = "N/A"
        else: 
            frame_string = str(frame_string)

        self.frame_label = QLabel("Frame: " + frame_string)
        #self.frame_layout.addWidget(self.frame_label)

        self.remove_frame = QPushButton("Cancel")
        self.remove_frame.clicked.connect(self.remove_frame_button)
        #self.frame_layout.addWidget(self.remove_frame)

        self.debug_btn = QPushButton("Debug")
        self.debug_btn.clicked.connect(self.debug)
        self.frame_layout.addWidget(self.debug_btn)

        self.insert_btn = QPushButton("Insert Waypoint")
        self.insert_btn.clicked.connect(self.insert_waypoint)
        self.frame_layout.addWidget(self.insert_btn)
        
        self.frame_widget.setLayout(self.frame_layout)
        self.frame_dock.setWidget(self.frame_widget)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.frame_dock)

        self.frame_dock.setMaximumSize(300, 200)

        # ////////////////////////////////    STL CONVERSION    ///////////////////////////////////
        self.random_btn_dock = QDockWidget("Control options", self)
        self.random_btn_widget = QWidget()
        self.random_btn_layout = QVBoxLayout()

        self.random_btn = QPushButton("Generate STL")
        self.random_btn.clicked.connect(self.generate_stl)
        #self.random_btn_layout.addWidget(self.random_btn)

        self.frame_layout.addWidget(self.random_btn)
        
        # self.random_btn_widget.setLayout(self.random_btn_layout)
        # self.random_btn_dock.setWidget(self.random_btn_widget)
        # self.addDockWidget(Qt.LeftDockWidgetArea, self.random_btn_dock)

        # self.random_btn_dock.setMaximumSize(300, 100)

    def insert_waypoint(self):
        print("test")

    def debug(self):
        for i in range (0, len(self.chain.Joints)):
            print("Joint:", i, "Children:", self.chain.Children[i], "Parent:", self.chain.Parents[i])

    def generate_stl(self):
        newTree = origamiToPrinted(self.chain, 0.05)

        self.plot_widget.clear()

        grid = gl.GLGridItem()
        self.plot_widget.addItem(grid)
        grid.setColor((0,0,0,255))

        for i in range(0,len(newTree.Children)):
            start = time.time()
            if len(newTree.Children[i]) > 0:
                filepath = newTree.exportLink3DFile(i, "test" + "/poses", pose=True)
                if filepath:
                    plotSTL(self.plot_widget, filepath, newTree.Joints[i].DistalDubinsFrame() @ SE3.Ry(np.pi/2) @SE3.Rz(-np.pi/2), color=(1,1,1,1))
                print(f"plotted links from {i}, Time: {time.time() - start}s")

        #export and plot all the joints
        for i in range(0,len(newTree.Joints)):
            start = time.time()
            if not isinstance(newTree.Joints[i],PrintedWaypoint):
                file1, rot1, file2, rot2 = newTree.Joints[i].renderPose("test")
                plotSTL(self.plot_widget, file1, newTree.Joints[i].ProximalDubinsFrame() @ rot1, color=(0,0,1,0))
                if (file2 != None):
                    plotSTL(self.plot_widget, file2, newTree.Joints[i].DistalDubinsFrame() @ rot2, color=(0,0,1,0))
            print(f"plotted joint {i}, Time: {time.time() - start}s")

        self.stl_generated = True
        #plotPrintedTree(newTree, "manualHandPrinted")

    def set_joint_as_frame(self):
        self.selected_frame = self.selected_joint
        self.frame_label.setText("Joint selected as frame: " + str(self.selected_frame))
        self.relativeSliderCheckbox.setDisabled(True)
        self.update_joint()

    def remove_frame_button(self):
        self.selected_frame = -1
        self.frame_label.setText("Joint selected as frame: N/A")
        self.relativeSliderCheckbox.setDisabled(False)
        self.update_joint()

    def init_key_bar(self):
        instructions = [
            "Middle Mouse Button: Pan Around",
            "W: Translate",
            "E: Rotate",
            "Delete: Delete Joint",
            "X: Select X Axis",
            "Y: Select Y Axis",
            "Z: Select Z Axis"
        ]
        for instruction in instructions:
            label = QLabel(instruction)
            label.setAlignment(Qt.AlignCenter) 
            self.key_bar_layout.addWidget(label)

        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.key_bar_layout.addWidget(spacer)

    def relative_clicked(self):
        self.is_local = self.relativeSliderCheckbox.isChecked()
        self.update_joint()

    def change_control_type(self):
        if self.control1.isChecked():
            self.control_type = "Translate"
        elif self.control2.isChecked():
            self.control_type = "Rotate"

        self.update_joint()

    def change_camera_type(self):
        if self.camera1.isChecked():
            self.plot_widget.camera_type = "Rotate"
        elif self.camera2.isChecked():
            self.plot_widget.camera_type = "Pan"

        self.update_joint()
    
    @QtCore.pyqtSlot(str)
    def key_pressed(self, key):
        if key == "Translate":
            self.control_type = key
            self.control1.setChecked(True)
        elif key == "Rotate":
            self.control_type = key
            self.control2.setChecked(True)
        elif key == "Delete":
            if (self.chain and self.selected_joint != -1):
                self.selected_joint = self.select_joint_options.currentIndex() 
                temp_last = len(self.chain.Joints) - 1
                if self.chain.delete(self.selected_joint):
                    self.reload_IDs()
                    self.update_joint()
                    success_dialog = SuccessDialog('Joint successfully deleted!')
                    success_dialog.exec_()
                    self.last_joint = temp_last
                else:
                    error_dialog = ErrorDialog('Error deleting joint.')
                    error_dialog.exec_()
        elif key == "X":
            self.arrow_selection_changed(0)
        elif key == "Y":
            self.arrow_selection_changed(1)
        elif key == "Z":
            self.arrow_selection_changed(2)

        self.update_joint()

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
            self.selected_axis_name = 'N/A'
            self.update_joint()
            self.update_rotation_slider()
            self.update_translate_slider()
            min = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[0])
            max = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[1])
            current = math.degrees(self.chain.Joints[self.selected_joint].state)
            self.current_state_label.setText(f"Min State: {int(min)} ≤ Current State: {int(current)} ≤ Max State: {int(max)}")
            self.update_state_slider()
            self.update_radius_slider()

    @QtCore.pyqtSlot(int)
    def arrow_selection_changed(self, index):
        if index != self.selected_arrow and self.selected_joint != -1:
            self.selected_arrow = index
                        
            if (index == 0): self.selected_axis_name = 'X'
            elif (index == 1): self.selected_axis_name = 'Y'
            elif (index == 2): self.selected_axis_name = 'Z'
            else: self.selected_axis_name = 'N/A'

            self.update_joint()
            self.update_rotation_slider()
            self.update_translate_slider()

    @QtCore.pyqtSlot(np.ndarray)
    def drag_translate(self, new_position):
        propogate = self.propogateSliderCheckbox.isChecked()
        relative = self.relativeSliderCheckbox.isChecked()
        old_position = self.chain.Joints[self.selected_joint].Pose.t
        trans = new_position - old_position
        transformation = SE3.Trans(trans[0], trans[1], trans[2])

        if self.chain.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=relative):
            self.update_joint()

        self.update_rotation_slider()
        self.update_translate_slider()

    @QtCore.pyqtSlot(float)
    def drag_rotate(self, new_rotation):
        transformation = SE3()
        propogate = self.propogateSliderCheckbox.isChecked()
        relative = self.relativeSliderCheckbox.isChecked()
        if (self.selected_axis_name == 'X'):
            transformation = SE3.Rx(new_rotation)
        elif (self.selected_axis_name == 'Y'):
            transformation = SE3.Ry(new_rotation)
        elif (self.selected_axis_name == 'Z'):
            transformation = SE3.Rz(new_rotation)

        if self.chain.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=relative):
            self.update_joint()
            self.update_rotation_slider()
            self.update_translate_slider()

    def update_rotation_slider(self):
        self.rotationSlider.setMinimum(-360)
        self.rotationSlider.setMaximum(360)
        if self.selected_arrow != -1:
            rotation_matrix = self.chain.Joints[self.selected_joint].Pose.R
            angle_degrees = self.rotation_angle_from_matrix(rotation_matrix, self.selected_arrow)
            self.rotationSlider.blockSignals(True)
            self.rotationLabel.setText(f"Rotate {self.selected_axis_name} Axis: {int(angle_degrees)}°")
            self.rotationSlider.setValue(int(angle_degrees))
            self.oldRotVal = angle_degrees
            self.rotationSlider.setDisabled(False)
            self.rotationSlider.blockSignals(False)
            self.rotationInput.blockSignals(True)
            self.rotationInput.setText(str(int(angle_degrees)))
            self.rotationInput.blockSignals(False)
            self.rotationInput.setDisabled(False)
        else:
            self.rotationSlider.blockSignals(True)
            self.rotationLabel.setText(f"Rotation {self.selected_axis_name} Axis: 0°")
            self.rotationSlider.setValue(0)
            self.rotationSlider.setDisabled(True)
            self.rotationSlider.blockSignals(False)
            self.rotationInput.blockSignals(True)
            self.rotationInput.setText("")
            self.rotationInput.blockSignals(False)
            self.rotationInput.setDisabled(True)

    def update_translate_slider(self):
        if self.selected_arrow != -1:
            amount = self.chain.Joints[self.selected_joint].Pose.t[self.selected_arrow]
            self.translate_slider.blockSignals(True)
            self.translate_label.setText(f"Translate {self.selected_axis_name} Axis: {int(amount * 10)}")
            self.translate_slider.setValue(int(amount * 10))
            self.oldTransVal = amount
            self.translate_slider.setDisabled(False)
            self.translate_slider.blockSignals(False)
            self.translationInput.blockSignals(True)
            self.translationInput.setText(str(int(amount * 10)))
            self.translationInput.blockSignals(False)
            self.translationInput.setDisabled(False)
        else:
            self.translate_slider.blockSignals(True)
            self.translate_label.setText(f"Translate {self.selected_axis_name} Axis: 0")
            self.translate_slider.setValue(0)
            self.translate_slider.setDisabled(True)
            self.translate_slider.blockSignals(False)
            self.translationInput.blockSignals(True)
            self.translationInput.setText("")
            self.translationInput.blockSignals(False)
            self.translationInput.setDisabled(True)

    def update_state_slider(self):
        min = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[0])
        max = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[1])
        if self.selected_joint != -1 and min != 0 and max != 0:
            self.state_slider.setMinimum(int(min))
            self.state_slider.setMaximum(int(max))
            current = int(math.degrees(self.chain.Joints[self.selected_joint].state))
            self.oldStateVal = current
            self.state_slider.blockSignals(True)
            self.state_label.setText(f"Edit Joint {self.selected_joint} State: {current}")
            self.state_slider.setValue(current)
            self.state_slider.setDisabled(False)
            self.state_slider.blockSignals(False)
            self.stateInput.blockSignals(True)
            self.stateInput.setText(str(current))
            self.stateInput.blockSignals(False)
            self.stateInput.setDisabled(False)
        else:
            self.state_slider.blockSignals(True)
            self.state_label.setText(f"Edit Joint N/A State: 0")
            self.state_slider.setValue(0)
            self.state_slider.setDisabled(True)
            self.state_slider.blockSignals(False)
            self.stateInput.blockSignals(True)
            self.stateInput.setText("")
            self.stateInput.blockSignals(False)
            self.stateInput.setDisabled(True)

    def update_radius_slider(self):
        if self.selected_joint != -1:
            radius = int(self.chain.Joints[self.selected_joint].r)
            self.oldRadiusVal = radius
            self.radius_slider.blockSignals(True)
            self.radius_label.setText(f"Edit Joint Radius: {radius}")
            self.radius_slider.setValue(radius)
            self.radius_slider.setDisabled(False)
            self.radius_slider.blockSignals(False)
            self.radiusInput.blockSignals(True)
            self.radiusInput.setText(str(radius))
            self.radiusInput.blockSignals(False)
            self.radiusInput.setDisabled(False)
        else:
            self.radius_slider.blockSignals(True)
            self.radius_label.setText(f"Edit Joint Radius: 0")
            self.radius_slider.setValue(1)
            self.radius_slider.setDisabled(True)
            self.radius_slider.blockSignals(False)
            self.radiusInput.blockSignals(True)
            self.radiusInput.setText("")
            self.radiusInput.blockSignals(False)
            self.radiusInput.setDisabled(True)

    def rotation_angle_from_matrix(self, rotation_matrix, axis):
        rot = R.from_matrix(rotation_matrix)
        euler_angles = rot.as_euler('xyz', degrees=True)
        return euler_angles[axis]
    
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

                if self.chain.setJointState(self.selected_joint, math.radians(edit)):
                    self.update_joint()
                    success_dialog = SuccessDialog('Joint state successfully edited!')
                    success_dialog.exec_()
                    min = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[0])
                    max = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[1])
                    current = math.degrees(self.chain.Joints[self.selected_joint].state)
                    self.current_state_label.setText(f"Min State: {int(min)} ≤ Current State: {int(current)} ≤ Max State: {int(max)}")
                else:
                    error_dialog = ErrorDialog('Error editing joint state.')
                    error_dialog.exec_()

    def adjust_rotation(self, value):
        # if not isinstance(value, float):
        #     value = value.strip()
        value = float(value) if len(str(value)) > 0 else 0
        angle_radians = math.radians(value - self.oldRotVal)
        self.rotationLabel.setText(f"Rotate {self.selected_axis_name} Axis: {int(value)}°")
        if self.chain and self.selected_joint != -1:
            if self.selected_arrow == 0:
                transformation = SE3.Rx(angle_radians)
            elif self.selected_arrow == 1:
                transformation = SE3.Ry(angle_radians)
            else:
                transformation = SE3.Rz(angle_radians)
            propogate = self.propogateSliderCheckbox.isChecked()
            relative = self.relativeSliderCheckbox.isChecked()
            if self.chain.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=relative):
                self.update_joint()
                self.oldRotVal = int(value)
                self.update_rotation_slider()
            else:
                self.rotationSlider.blockSignals(True)
                self.rotationSlider.setValue(int(self.oldRotVal))
                if value > self.oldRotVal:
                    self.rotationSlider.setMaximum(int(self.oldRotVal))
                else:
                    self.rotationSlider.setMinimum(int(self.oldRotVal))
                self.rotationLabel.setText(f"Rotate {self.selected_axis_name} Axis: {int(self.oldRotVal)}°")
                self.rotationSlider.blockSignals(False)

    def adjust_translation(self, value):
        if not isinstance(value, float) and not isinstance(value, int):
            value = value.strip()
        value = float(value) if value else 0
        actualVal = value / 10
        self.translate_label.setText(f'Translate {self.selected_axis_name} Axis: {int(value)}')
        amount = actualVal - self.oldTransVal
        if self.chain and self.selected_joint != -1:
            propogate = self.propogateSliderCheckbox.isChecked()
            relative = self.relativeSliderCheckbox.isChecked()
            transformation = SE3()
            if (self.selected_arrow == 0):
                transformation = SE3.Tx(amount)
            if (self.selected_arrow == 1):
                transformation = SE3.Ty(amount)
            if (self.selected_arrow == 2):
                transformation = SE3.Tz(amount)
            if self.chain.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=relative):
                self.update_joint()
                self.oldTransVal = actualVal
                self.update_translate_slider()
            else:
                self.translate_slider.blockSignals(True)
                self.translate_slider.setValue(int(self.oldTransVal * 10))
                self.translate_label.setText(f"Translate {self.selected_axis_name} Axis: {int(self.oldTransVal * 10)}")
                self.translate_label.blockSignals(False)

    def adjust_state(self, value):
        if not isinstance(value, float) and not isinstance(value, int):
            value = value.strip()
        value = float(value) if value else 0
        actualVal = math.radians(value)
        self.state_label.setText(f'Edit Joint {self.selected_joint} State: {int(value)}')
        if self.chain and self.selected_joint != -1:
            if self.chain.setJointState(self.selected_joint, actualVal):
                self.update_joint()
                self.OldStateVal = value
                min = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[0])
                max = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[1])
                current = math.degrees(self.chain.Joints[self.selected_joint].state)
                self.current_state_label.setText(f"Min State: {int(min)} ≤ Current State: {int(current)} ≤ Max State: {int(max)}")
                self.update_state_slider()
            else:
                self.state_slider.blockSignals(True)
                self.state_slider.setValue(int(self.oldStateVal))
                self.state_label.setText(f"Edit Joint {self.selected_joint} State: {int(self.oldStateVal)}")
                self.state_label.blockSignals(False)

    def adjust_radius(self, value):
        if not isinstance(value, float) and not isinstance(value, int):
            value = value.strip()
        value = float(value) if value else 0
        self.radius_label.setText(f'Edit Joint Radius: {int(value)}')
        if self.chain and self.selected_joint != -1:
            if self.chain.changeRadius(value):
                self.r = value
                self.update_joint()
                self.OldRadiusVal = value
                self.update_radius_slider()
            else:
                self.radius_slider.blockSignals(True)
                self.radius_slider.setValue(int(self.oldStateVal))
                self.radius_label.setText(f"Edit Joint Radius: {int(self.oldRadiusVal)}")
                self.radius_label.blockSignals(False)

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
            temp_last = len(self.chain.Joints) - 1
            if self.chain.delete(self.selected_joint):
                self.reload_IDs()
                self.update_joint()
                success_dialog = SuccessDialog('Joint successfully deleted!')
                success_dialog.exec_()
                self.last_joint = temp_last
            else:
                error_dialog = ErrorDialog('Error deleting joint.')
                error_dialog.exec_()

    def reload_IDs(self):
        if self.chain is not None:
            for index, joint in enumerate(self.chain.Joints):
                joint.id = index

    def update_joint(self):
        self.select_joint_options.blockSignals(True)

        if (not self.stl_generated):
            self.plot_widget.clear()
            self.select_joint_options.clear()

            grid = gl.GLGridItem()
            self.plot_widget.addItem(grid)
            grid.setColor((0,0,0,255))
            
            if self.chain is not None:
                self.chain.addToWidget(self, selectedJoint=self.selected_joint, lastJoint = self.last_joint)

                self.select_joint_options.blockSignals(False)

                if (self.control_type == "Translate"):
                    for joint in self.chain.Joints:
                        if (joint.id == self.selected_joint):
                            if (self.selected_frame == -1):
                                joint.addTranslateArrows(self, selectedArrow=self.selected_arrow, local=self.is_local)
                            else:
                                frame_joint = self.chain.Joints[self.selected_frame]
                                joint.addTranslateArrows(self, selectedArrow=self.selected_arrow, local=self.is_local, frame=frame_joint.Pose)
                elif (self.control_type == "Rotate"):
                    for joint in self.chain.Joints:
                        if (joint.id == self.selected_joint):
                            if (self.selected_frame == -1):
                                joint.addRotateArrows(self, selectedArrow=self.selected_arrow, local=self.is_local)
                            else: 
                                frame_joint = self.chain.Joints[self.selected_frame]
                                joint.addRotateArrows(self, selectedArrow=self.selected_arrow, local=self.is_local, frame=frame_joint.Pose)

            if self.selected_arrow != -1:
                self.rotationSlider.setDisabled(False)
                self.translate_slider.setDisabled(False)
            else:
                self.rotationSlider.setDisabled(True)
                self.translate_slider.setDisabled(True)

    def update_plot(self):
        self.plot_widget.clear()
        self.setCentralWidget(self.plot_widget)

        grid = gl.GLGridItem()
        self.plot_widget.addItem(grid)
        grid.setColor((0,0,0,255))

        for index, joint in enumerate(self.chain.Joints):
            joint.id = index

        self.chain.addToWidget(self, selectedJoint=self.selected_joint, lastJoint = self.last_joint)

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

            if (self.chain == None) or len(self.chain.Joints) == 0:
                self.chain = KinematicChain(joint)
            elif joint.id != 0:
                self.chain.addJoint(parentIndex = self.selected_joint, newJoint = joint)
            else:
                self.chain.append(joint)

            self.update_plot()
            self.select_joint_options.setCurrentIndex(len(self.chain.Joints) - 1)

    def chain_not_created(self):
        error_dialog = ErrorDialog('Please create a chain first.')
        error_dialog.exec_()
        self.create_new_chain_func()

    def is_parent_joint_selected(self):
        if self.selected_joint == -1 and self.chain is not None:
            QMessageBox.warning(self, "Selection Required", "Please select a parent joint before adding a new one.")
            return False
        return True

    def add_prismatic_func(self):
        if (not self.chain_created):
            self.chain_not_created()
        elif self.is_parent_joint_selected():
            dialog = AddPrismaticDialog(self.numSides, self.r)
            self.add_joint(dialog)

    def add_revolute_func(self):
        if (not self.chain_created):
            self.chain_not_created()
        elif self.is_parent_joint_selected():
            dialog = AddRevoluteDialog(self.numSides, self.r)
            self.add_joint(dialog)

    def add_waypoint_func(self):
        if (not self.chain_created):
            self.chain_not_created()
        elif self.is_parent_joint_selected():
            waypoint = Waypoint(self.numSides, self.r, SE3())
            
            if (self.chain == None):
                waypoint.id = 0
            else:
                waypoint.id = len(self.chain.Joints)

            if (self.chain == None) or len(self.chain.Joints) == 0:
                self.chain = KinematicChain(waypoint)
            elif waypoint.id != 0:
                self.chain.addJoint(parentIndex = self.selected_joint, newJoint = waypoint, safe=False)
            else:
                self.chain.append(waypoint, safe=False)

            self.update_plot()
            self.select_joint_options.setCurrentIndex(len(self.chain.Joints) - 1)

    def add_tip_func(self):
        if (not self.chain_created):
            self.chain_not_created()
        elif self.is_parent_joint_selected():
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

    def keyPressEvent(self, event: QKeyEvent):
        if event.key() == Qt.Key_W:
            self.control_type = "Translate"
            self.control1.setChecked(True)
            self.update_joint()
        elif event.key() == Qt.Key_E:
            self.control_type = "Rotate"
            self.control2.setChecked(True)
            self.update_joint()
        elif event.key() == Qt.Key_Delete or event.key() == Qt.Key_Backspace:
            if (self.chain and self.selected_joint != -1):
                self.selected_joint = self.select_joint_options.currentIndex() 
                temp_last = len(self.chain.Joints) - 1
                if self.chain.delete(self.selected_joint):
                    self.reload_IDs()
                    self.update_joint()
                    success_dialog = SuccessDialog('Joint successfully deleted!')
                    success_dialog.exec_()
                    self.last_joint = temp_last
                else:
                    error_dialog = ErrorDialog('Error deleting joint.')
                    error_dialog.exec_()
        elif event.key() == Qt.Key_X:
            self.arrow_selection_changed(0)
        elif event.key() == Qt.Key_Y:
            self.arrow_selection_changed(1)
        elif event.key() == Qt.Key_Z:
            self.arrow_selection_changed(2)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PointEditorWindow()
    window.show()
    sys.exit(app.exec_())