"""
@author: Raymond Feng, Andy Wang, Daniel Feshbach
"""

import sys, os
import dill
import time as _time_module
import threading as _threading

if sys.stdout is None:
    class DummyStream:
        def write(self, data):
            pass
        def flush(self):
            pass
        def isatty(self):
            return False

    sys.stdout = DummyStream()
if sys.stderr is None:
    sys.stderr = DummyStream()
    
import numpy as np
import pyqtgraph.opengl as gl
import PyQt5
from PyQt5 import QtWidgets
from PyQt5 import QtCore as qc
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton, QDockWidget, QComboBox, QHBoxLayout, QLabel, QDialog, QLineEdit, QCheckBox, QMessageBox, QButtonGroup, QRadioButton, QSlider, QSizePolicy, QFileDialog, QShortcut, QGridLayout, QSpinBox
from PyQt5.QtCore import Qt, pyqtSignal, QTimer, QTime, QEvent
from PyQt5.QtGui import QPixmap, QSurfaceFormat, QKeyEvent, QPixmap, QIcon, QMatrix4x4, QVector3D, QMatrix3x3, QKeySequence
from PyQt5 import sip
from pyqtgraph.Qt import QtCore
import pyqtgraph as pg
from OpenGL.GL import *
from OpenGL.GLU import *
from spatialmath import SE3
import math
from PathCSC import *
from KinematicChain import *
from OrigamiTube import *
from Joint import Prismatic, Revolute
import re
from scipy.spatial.transform import Rotation as R
from style import *
from ReferenceMesh import *
from Dialog import *
from origamiJointWidget import *
from IntersectionHelper import *
from measurementWidget import MeasurementWidget

import warnings
warnings.filterwarnings("ignore")

if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    PyQt5.QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    PyQt5.QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

from guiBase import DeleteWidget, AddMeshWidget, EditGridWidget, CollapsibleDockWidget, ImageRadioButton, OverlayLine, BaseClickableGLViewWidget, BaseWindowKinegamiGUI, EnvironmentWidget

class AddChainWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Create New Chain')

        layout = QVBoxLayout()

        # Input for the number of sides
        num_sides_label = QLabel("Number of Sides:")
        layout.addWidget(num_sides_label)
        self.num_sides_input = QLineEdit()
        self.num_sides_input.setPlaceholderText("Enter number of sides")
        layout.addWidget(self.num_sides_input)

        # Input for the radius
        radius_label = QLabel("Radius:")
        layout.addWidget(radius_label)
        self.radius_input = QLineEdit()
        self.radius_input.setPlaceholderText("Enter radius")
        layout.addWidget(self.radius_input)

        # Apply button to create the chain
        self.create_button = QPushButton('Create Chain', self)
        self.create_button.clicked.connect(self.on_create_clicked)
        layout.addWidget(self.create_button)

        # Cancel button to close the widget
        self.cancel_button = QPushButton('Cancel', self)
        self.cancel_button.clicked.connect(self.on_cancel_clicked)
        layout.addWidget(self.cancel_button)

        self.setLayout(layout)

    def on_create_clicked(self):
        try:
            # Get input values
            numSides = int(self.num_sides_input.text())
            radius = float(self.radius_input.text())
            
            # Call a function in the main window to create the new chain
            self.window().create_new_chains(numSides, radius)            
            # Optionally hide the widget after successful creation
            self.window().add_chain_popup_dock.setVisible(False)
            self.window().add_chain_dock.setVisible(True)
            self.window().radius_slider.setEnabled(True)
            #self.window().joint_range_slider.setEnabled(True)
        except ValueError:
            self.show_error("Please enter valid integers.")

    def on_cancel_clicked(self):
        # Hide the widget if the user cancels
        self.window().add_chain_popup_dock.setVisible(False)
        self.window().add_chain_dock.setVisible(True)

    def show_error(self, message):
        QMessageBox.warning(self, "Invalid Input", message)

class EditDimensionsWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Edit Dimensions')

        self.layout = QVBoxLayout()

        self.radio_layout = QVBoxLayout()

        self.radio1 = QRadioButton("Rescale Dimensions")
        self.radio1.setChecked(True)
        self.radio2 = QRadioButton("Preserve Dimensions")

        self.button_group = QButtonGroup(self)
        self.button_group.addButton(self.radio1)
        self.button_group.addButton(self.radio2)

        # Add widgets to layout
        self.radio_layout.addWidget(self.radio1)
        self.radio_layout.addWidget(self.radio2)

        self.units_layout = QHBoxLayout()
        self.label1 = QLabel("Target units:")
        self.units_layout.addWidget(self.label1)

        self.target_unit = QComboBox()
        self.target_unit.addItems(["Centimeter (cm)", "Inch (in)"])
        self.target_unit.currentIndexChanged.connect(self.crease_pattern_unit_changed)
        self.target_unit.setCurrentIndex(0)
        self.units_layout.addWidget(self.target_unit)

        self.button_layout = QHBoxLayout()

        # Apply button to create the chain
        self.apply_button = QPushButton('Apply Changes', self)
        self.apply_button.clicked.connect(self.on_apply_clicked)
        self.button_layout.addWidget(self.apply_button)

        # Cancel button to close the widget
        self.cancel_button = QPushButton('Cancel', self)
        self.cancel_button.clicked.connect(self.on_cancel_clicked)
        self.button_layout.addWidget(self.cancel_button)

        self.layout.addLayout(self.radio_layout)
        self.layout.addLayout(self.units_layout)
        self.layout.addLayout(self.button_layout)

        self.setLayout(self.layout)

    target_units = qc.pyqtSignal(str)

    def on_apply_clicked(self):
        units = self.target_unit.currentText()
        if units == "Centimeter (cm)":
            units_short = "cm"
        elif units == "Inch (in)":
            units_short = "in"
        
        self.window().edit_grid_widget.unit_label.setText(f"({units_short})")

        prev_units = self.window().units
        self.window().units = units
        self.target_units.emit(units)
        self.window().edit_dims_dock.setVisible(False)

        if self.radio1.isChecked():
            self.window().rescale_dimensions(prev_units, units)
        elif self.radio2.isChecked():
            self.window().preserve_dimensions(prev_units, units)

    def on_cancel_clicked(self):
        # Hide the widget if the user cancels
        self.window().edit_dims_dock.setVisible(False)

    def crease_pattern_unit_changed(self):
        print("Selected unit:", self.target_unit.currentText())
    
class ClickableGLViewWidget(BaseClickableGLViewWidget):
    def get_world_coordinates(self, event):
        pos = event.localPos()
        ndc_x = (2.0 * pos.x()) / self.width() - 1.0
        ndc_y = 1.0 - (2.0 * pos.y()) / self.height()

        ndc = QVector3D(ndc_x, ndc_y, -1.0)
        view = self.viewMatrix()
        proj = self.projectionMatrix()

        inverted_matrix = (proj * view).inverted()[0]

        near_point = inverted_matrix.map(ndc)
        ndc.setZ(1.0)
        far_point = inverted_matrix.map(ndc)

        direction = far_point - near_point
        direction.normalize()

        return near_point, direction

    def mouseMoveEvent(self, event):
        if self.parent_window.measure_select_mode:
            if (event.buttons() and (Qt.LeftButton or Qt.MiddleButton)) and (event.pos() - self.drag_start_pos).manhattanLength() >= QApplication.startDragDistance():
                self.is_dragging = True
            if not self.is_dragging:
                return
            self._handle_camera_drag(event)
            return

        if (event.buttons() and (Qt.LeftButton or Qt.MiddleButton)) and (event.pos() - self.drag_start_pos).manhattanLength() >= QApplication.startDragDistance():
            self.is_dragging = True

        if (self.is_dragging):
            if (self.selected_axis):
                new_pos_3D = self.get_closest_point(event)

                if (self.mesh_selected):
                    selected_joint = self.parent_window.referenceMesh
                else:
                    selected_joint = self.parent_window.chain.Joints[self.parent_window.selected_joint]
                joint_center = selected_joint.Pose.t
                qsphere_start = QVector3D(joint_center[0], joint_center[1], joint_center[2])

                trans = new_pos_3D - qsphere_start
                transformation = SE3.Trans(trans[0], trans[1], trans[2])

                propogate = self.parent_window.propogate_slider_checkbox.isChecked()

                if (self.mesh_selected):
                    self.parent_window.referenceMesh.transform(transformation)
                else:
                    # Throttle: skip stale queued events to prevent Qt event-queue backlog.
                    _now = _time_module.perf_counter()
                    if _now - getattr(self.parent_window, '_last_drag_frame_t', 0.0) < 0.200:
                        return
                    self.parent_window._last_drag_frame_t = _now

                    # Save backup before first drag transform
                    if not hasattr(self, '_drag_backup') or self._drag_backup is None:
                        self._drag_backup = self.parent_window.chain.dataDeepCopy()

                    self.parent_window.chain.transformJoint(self.parent_window.selected_joint, transformation, propogate=propogate, relative=False)
                    
                    # Check consistency — revert and stop drag if broken
                    if not self.parent_window.chain.isConsistent():
                        print("WARNING: Tree became inconsistent during drag translation. Reverting.")
                        self.parent_window.chain.setTo(self._drag_backup)
                        self._drag_backup = None
                        self.is_dragging = False
                        self.selected_axis = None
                        self.parent_window.update_joint()
                        return
                    # Keep a rolling backup of the last good state
                    self._drag_backup = self.parent_window.chain.dataDeepCopy()
                
                self.parent_window.update_joint()
            elif (self.selected_torus):
                # Throttle: return BEFORE get_axis_angle_delta so last_drag_pos is NOT
                # updated — the next processed frame then accumulates the full mouse motion.
                _now = _time_module.perf_counter()
                if _now - getattr(self.parent_window, '_last_drag_frame_t', 0.0) < 0.200:
                    return
                self.parent_window._last_drag_frame_t = _now

                da, normal = self.get_axis_angle_delta(event)
                if (not self.facing_same_dir):
                    da = -da

                if (self.mesh_selected):
                    pass
                else:
                    # Save backup before first drag transform
                    if not hasattr(self, '_drag_backup') or self._drag_backup is None:
                        self._drag_backup = self.parent_window.chain.dataDeepCopy()

                    self.parent_window.rotate_joint(da, self.selected_axis_orig)

                    # Check consistency — revert and stop drag if broken
                    if not self.parent_window.chain.isConsistent():
                        print("WARNING: Tree became inconsistent during drag rotation. Reverting.")
                        self.parent_window.chain.setTo(self._drag_backup)
                        self._drag_backup = None
                        self.is_dragging = False
                        self.selected_torus = None
                        self.parent_window.update_joint()
                        return
                    # Keep a rolling backup of the last good state
                    self._drag_backup = self.parent_window.chain.dataDeepCopy()

                self.parent_window.update_joint()
            else:
                self._handle_camera_drag(event)


class WindowKinegamiGUI(BaseWindowKinegamiGUI):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Kinegami Interactive Editor")
        # Default window to fill available screen space
        self.setWindowState(Qt.WindowMaximized)

        self.plot_widget = ClickableGLViewWidget(parent_window=self)
        self.setCentralWidget(self.plot_widget)

        self.plot_widget.setBackgroundColor(backgroundColorDefault)

        self.grid = gl.GLGridItem()
        self.plot_widget.addItem(self.grid)
        self.grid.setColor(gridColorDefault)
        self.grid_on = True

        self.units = "Centimeter (cm)"

        self.current_point = 0
        self.chain = None
        self.versions = []
        self.version_index = -1
        self.total_version_counter = 0
        self.chain_created = False
        self.stl_generated = False
        self.referenceMesh = None
        self._lightweight_dirty = False  # True when link geometry is stale from lightweight updates

        self.control_type = "Translate"
        self.is_local = True

        self.selected_joint = -1
        self.selected_arrow = -1
        self.selected_link = -1
        self.selected_axis_name = 'N/A'
        self.last_joint = -1
        self.selected_frame = -1
        self.mesh_selected = False
        self.mesh_visible = True
        self.mesh_scale = 1.0

        self.grid_color = gridColorDefault
        self.grid_spacing = 1.0
        self.grid_size = 10

        self.num_sides = 4
        self.radius = 1.0

        self.environment_items = []
        self.environment_visible = True
        self._last_drag_frame_t = 0.0

        self.measure_select_mode = False
        self.measure_active_endpoint = 0
        self.measure_points = [None, None]
        self._measurement_line_item = None

        self.plot_widget.click_signal.connect(self.joint_selection_changed)
        self.plot_widget.click_signal_arrow.connect(self.arrow_selection_changed)
        self.plot_widget.click_signal_link.connect(self.link_selection_changed)
        self.plot_widget.click_signal_mesh.connect(self.mesh_selected_slot)
        self.plot_widget.drag_change_position.connect(self.drag_translate)
        self.plot_widget.drag_change_rotation.connect(self.drag_rotate)
        self.plot_widget.done_transforming.connect(self.done_transforming)
        self.plot_widget.key_pressed.connect(self.key_pressed)
        
        # //////////////////////////////////    Keyboard Options    ///////////////////////////////////
        top_dock_widget = CollapsibleDockWidget("Keyboard Controls", self)
        top_dock_widget.setAllowedAreas(Qt.TopDockWidgetArea)

        self.key_bar = QWidget()
        self.key_bar_layout = QHBoxLayout(self.key_bar)  # Layout is initialized and set to the widget here
        self.key_bar_layout.setContentsMargins(0, 0, 0, 0)
        self.key_bar.setFixedHeight(20)
        self.init_key_bar()

        top_dock_widget.setWidget(self.key_bar)
        
        self.environment_widget = EnvironmentWidget(self)
        self.add_mesh_widget = AddMeshWidget(self)
        self.mesh_scale = 1.0
        self.add_mesh_widget.change_scale.connect(self.change_mesh_scale)

        # //////////////////////////////////    MESSAGE DISPLAY    ///////////////////////////////////
        self.message_display_dock = CollapsibleDockWidget("Messages", self)
        self.message_display_dock.setAllowedAreas(Qt.BottomDockWidgetArea)

        self.status_label = QLabel('')
        self.status_label.setWordWrap(True)
        self.message_display_dock.setWidget(self.status_label)
        self.message_display_dock.setVisible(False)

        # //////////////////////////////////    CONFIGURATIONS    ///////////////////////////////////
        self.configurations_widget = QWidget()
        self.configurations_widget.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        self.configurations_layout = QVBoxLayout(self.configurations_widget)
        self.configurations_layout.setSpacing(2)
        self.configurations_layout.setContentsMargins(0, 0, 0, 0)
        self.config_text_boxes = []  # List to store text boxes for joint states
        self.config_sliders = []  # List to store sliders for joint states
        self.config_joint_indices = []  # List to map text box index to actual joint index
        self.saved_configurations = []  # List of saved configurations (each is a list of joint states)
        self.config_durations = []  # List of durations (in seconds) between each pair of configs
        self.animation_timer = None  # QTimer for animation
        self.is_animating = False  # Animation state
        self.animation_start_time = 0  # Animation start timestamp
        self.animation_start_value = 0  # Slider value when animation started
        self.animation_loop = False  # Whether animation should loop
        
        self.configurations_dock = CollapsibleDockWidget("Configurations and Motion", self)
        self.configurations_dock.setWidget(self.configurations_widget)
        self.configurations_dock.setVisible(True)

        # //////////////////////////////////    ADD JOINTS    ///////////////////////////////////
        self.add_prismatic = QPushButton("Add Prismatic Joint")
        self.add_prismatic_menu = AddPrismaticMenu(self)
        self.add_prismatic_menu.setVisible(False)
        self.add_revolute = QPushButton("Add Revolute Joint")
        self.add_revolute_menu = AddRevoluteMenu(self)
        self.add_revolute_menu.setVisible(False)
        self.add_tip = QPushButton("Add Tip")
        self.add_tip_menu = AddTipMenu(self)
        self.add_tip_menu.setVisible(False)
        self.create_new_chain = QPushButton("Create New Chain")
        self.edit_dimension_menu = EditDimensionMenu(self)
        self.edit_dimension_menu.setVisible(False)
        self.edit_dimension_button = QPushButton("Edit Dimension")

        add_waypoints_layout = QVBoxLayout()
        self.add_waypoint = QPushButton("Add Waypoint")
        add_waypoints_layout.addWidget(self.add_waypoint)
        insert_wp_row = QHBoxLayout()
        self.insert_waypoint_button = QPushButton("Insert Intermediate Waypoint(s)")
        self.insert_waypoint_count = QSpinBox()
        self.insert_waypoint_count.setMinimum(1)
        self.insert_waypoint_count.setMaximum(99)
        self.insert_waypoint_count.setValue(1)
        self.insert_waypoint_count.setFixedWidth(50)
        insert_wp_row.addWidget(self.insert_waypoint_button)
        insert_wp_row.addWidget(self.insert_waypoint_count)
        add_waypoints_layout.addLayout(insert_wp_row)

        radius_layout = QVBoxLayout()
        radius_label = QLabel("Chain Radius:")
        self.radius_slider = QSlider(Qt.Horizontal, self)
        self.radius_slider.setMinimum(0)  # Minimum value
        self.radius_slider.setMaximum(100)  # Maximum value
        self.radius_slider.setValue(10)  # Initial value
        self.radius_slider.setEnabled(False)
        self.radius_slider.valueChanged.connect(self.onUpdateRadius)
        radius_layout.addWidget(radius_label)
        radius_layout.addWidget(self.radius_slider)

        add_joints_layout = QVBoxLayout()
        add_chain_layout = QVBoxLayout()
        # Add-to-root checkbox moved here from Edit Joints
        self.joint_add_root_button = QCheckBox("Add to Root")
        self.joint_add_root_button.setChecked(False)
        self.joint_add_root_button.stateChanged.connect(self.add_to_root_func)
        add_joints_layout.addWidget(self.add_prismatic)
        add_joints_layout.addWidget(self.add_prismatic_menu)
        add_joints_layout.addWidget(self.add_revolute)
        add_joints_layout.addWidget(self.add_revolute_menu)
        add_joints_layout.addLayout(add_waypoints_layout)
        add_joints_layout.addWidget(self.add_tip)
        add_joints_layout.addWidget(self.add_tip_menu)
        add_joints_layout.addWidget(self.edit_dimension_button)
        add_joints_layout.addWidget(self.edit_dimension_menu)
        # place Add-to-root checkbox at bottom
        add_joints_layout.addWidget(self.joint_add_root_button)

        add_chain_layout.addWidget(self.create_new_chain)

        self.add_prismatic.clicked.connect(self.add_prismatic_toggle)
        self.add_revolute.clicked.connect(self.add_revolute_toggle)
        self.add_waypoint.clicked.connect(self.add_waypoint_func)
        self.insert_waypoint_button.clicked.connect(self.insert_waypoint_func)
        self.add_tip.clicked.connect(self.add_tip_toggle)
        self.create_new_chain.clicked.connect(self.create_new_chain_func)
        self.edit_dimension_button.clicked.connect(self.edit_joint_dimension)

        self.add_chain_dock = CollapsibleDockWidget("New Chain", self)
        self.add_chain_button_widget = QWidget()
        self.add_chain_button_widget.setLayout(add_chain_layout)
        self.add_chain_dock.setWidget(self.add_chain_button_widget)

        add_joints_dock = CollapsibleDockWidget("Add Joints", self)
        self.add_joints_widget = QWidget()
        self.add_joints_widget.setLayout(add_joints_layout)
        add_joints_dock.setWidget(self.add_joints_widget)

        self.add_to_root = False

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
        self.editing_widget = QWidget()
        self.joint_editing_layout = QVBoxLayout()

        self.translate_joint_label = QLabel("Translate")
        self.rotate_joint_label = QLabel("Rotate")
        try:
            base_path = sys._MEIPASS
        except AttributeError:
            base_path = os.path.abspath(".")
        move_unchecked = os.path.join(base_path, "ui/move_unchecked.png")
        move_checked = os.path.join(base_path, "ui/move_checked.png")
        rotate_unchecked = os.path.join(base_path, "ui/rotate_unchecked.png")
        rotate_checked = os.path.join(base_path, "ui/rotate_checked.png")
        self.translate_joint_radio_button = ImageRadioButton(move_unchecked, move_checked, "Translate")
        self.rotate_joint_radio_button = ImageRadioButton(rotate_unchecked, rotate_checked, "Rotate")
        self.translate_joint_radio_button.setChecked(True)
        self.translate_joint_radio_button.toggled.connect(self.change_control_type)
        self.rotate_joint_radio_button.toggled.connect(self.change_control_type)
        
        self.editing_widget.setLayout(self.joint_editing_layout)

        # "Add to Root" checkbox (moved to Add Joints dock)
        self.select_joint_options = QComboBox()
        self.select_link_options = QComboBox()
        self.delete_joint_button = QPushButton("Delete Joint")
        self.current_state_label = QLabel('Min State ≤ Current State ≤ Max State')
        self.current_state_label.setWordWrap(True)

        #joint_layout = QVBoxLayout()
        self.joint_editing_layout.addWidget(self.select_joint_options)
        self.joint_editing_layout.addWidget(self.delete_joint_button)
        self.joint_editing_layout.addWidget(self.select_link_options)
        edit_joints_dock = CollapsibleDockWidget("Edit Joints", self)
        #self.joint_editing_layout.addWidget(button_widget)
        edit_joints_dock.setWidget(self.editing_widget)

        self.delete_joint_button.clicked.connect(self.delete_joint)
        self.select_joint_options.currentIndexChanged.connect(self.joint_selection_changed)
        self.select_link_options.currentIndexChanged.connect(self.link_selection_changed)

        # self.rotationLabel = QLabel("Rotate N/A Axis: 0°")
        self.rotation_slider = QSlider(Qt.Horizontal)
        self.rotation_slider.setMinimum(-360)
        self.rotation_slider.setMaximum(360)
        self.rotation_slider.setValue(0)
        self.rotation_slider.setDisabled(True) 
        self.rotation_slider.sliderMoved.connect(self.adjust_rotation)
        self.rotation_slider.sliderReleased.connect(self.rotation_slider_released)

        # self. = QLabel('Translate N/A Axis: 0', self)
        self.translation_slider = QSlider(Qt.Horizontal, self)
        self.translation_slider.sliderMoved.connect(self.adjust_translation)
        self.translation_slider.sliderReleased.connect(self.translation_slider_released)

        # self.state_label = QLabel('Edit Joint N/A State: 0', self)
        self.state_slider = QSlider(Qt.Horizontal, self)
        self.state_slider.setMinimum(-100)
        self.state_slider.setMaximum(100)
        self.state_slider.setValue(0)
        self.state_slider.sliderMoved.connect(self.state_slider_moved)
        self.state_slider.sliderReleased.connect(self.state_slider_released)

        translation_layout = QVBoxLayout()
        translation_header_layout = QHBoxLayout()
        translation_header_layout.addWidget(self.translate_joint_radio_button)
        translation_header_layout.addWidget(self.translate_joint_label)
        translation_layout.addLayout(translation_header_layout)
        translation_slider_layout = QHBoxLayout()
        self.translation_textbox = QLineEdit(self)
        self.translation_textbox.setPlaceholderText("Enter distance")
        self.translation_textbox.returnPressed.connect(self.translation_textbox_return)
        translation_slider_layout.addWidget(self.translation_slider)
        translation_slider_layout.addWidget(self.translation_textbox)
        translation_layout.addLayout(translation_slider_layout)

        rotation_layout = QVBoxLayout()
        rotation_header_layout = QHBoxLayout()
        rotation_header_layout.addWidget(self.rotate_joint_radio_button)
        rotation_header_layout.addWidget(self.rotate_joint_label)
        rotation_layout.addLayout(rotation_header_layout)
        rotation_slider_layout = QHBoxLayout()
        self.rotation_textbox = QLineEdit(self)
        self.rotation_textbox.setPlaceholderText("Enter angle in degrees")
        self.rotation_textbox.returnPressed.connect(self.rotation_textbox_return)
        rotation_slider_layout.addWidget(self.rotation_slider)
        rotation_slider_layout.addWidget(self.rotation_textbox)
        rotation_layout.addLayout(rotation_slider_layout)

        state_layout = QVBoxLayout()
        state_layout.addWidget(self.current_state_label)
        state_slider_layout = QHBoxLayout()
        self.state_textbox = QLineEdit(self)
        self.state_textbox.setPlaceholderText("Enter state")
        self.state_textbox.returnPressed.connect(self.state_textbox_return)
        state_slider_layout.addWidget(self.state_slider)
        state_slider_layout.addWidget(self.state_textbox)
        state_layout.addLayout(state_slider_layout)

        checkbox_layout = QVBoxLayout() 
        self.propogate_slider_checkbox = QCheckBox("Propagate")
        self.local_orient_slider_checkbox = QCheckBox("Local Orientation")
        self.propogate_slider_checkbox.setChecked(True)
        self.local_orient_slider_checkbox.setChecked(True)
        self.local_orient_slider_checkbox.stateChanged.connect(self.local_orient_clicked)
        checkbox_layout.addWidget(self.propogate_slider_checkbox)
        checkbox_layout.addWidget(self.local_orient_slider_checkbox)

        """
        joint_range_layout = QHBoxLayout()
        joint_range_label = QLabel("Joint Range of Motion:")
        self.joint_range_slider = QSlider(Qt.Horizontal, self)
        self.joint_range_slider.setMinimum(0)  # Minimum value
        self.joint_range_slider.setMaximum(180)  # Maximum value
        self.joint_range_slider.setValue(90)  # Initial value
        self.joint_range_slider.setEnabled(False)
        self.joint_range_slider.valueChanged.connect(self.onUpdateJointState)
        joint_range_layout.addWidget(joint_range_label)
        joint_range_layout.addWidget(self.joint_range_slider)
        """

        self.joint_editing_layout.addLayout(state_layout)
        self.joint_editing_layout.addLayout(checkbox_layout)
        self.joint_editing_layout.addLayout(rotation_layout)
        self.joint_editing_layout.addLayout(translation_layout)
        self.joint_editing_layout.addLayout(radius_layout)
        #self.joint_editing_layout.addLayout(joint_range_layout)

        self.old_rot_val = 0
        self.old_trans_val = 0
        self.old_state_slider_val = 0
        self.old_radius_val = 1
        self.rotation_slider.setDisabled(True)
        self.translation_slider.setDisabled(True)
        self.state_slider.setDisabled(True)
        self.rotation_textbox.setDisabled(True)
        self.translation_textbox.setDisabled(True)
        self.state_textbox.setDisabled(True)

        self.reset_translation_tools()
        self.reset_rotation_tools()
        self.set_state_tools()

        # ////////////////////////////////    OPTIONS    ///////////////////////////////////
        self.options_dock = CollapsibleDockWidget("Options", self)
        self.options_widget = QWidget()
        self.options_layout = QVBoxLayout()

        #self.debug_btn = QPushButton("Debug")
        #self.debug_btn.clicked.connect(self.debug)
        #self.options_layout.addWidget(self.debug_btn)

        # Keyboard shortcuts for undo/redo
        # Ctrl+Z / Ctrl+Y on Windows/Linux, Cmd+Z / Cmd+Shift+Z on Mac
        self.undo_shortcut = QShortcut(QKeySequence.Undo, self)
        self.undo_shortcut.activated.connect(self.undo)
        self.redo_shortcut = QShortcut(QKeySequence.Redo, self)
        self.redo_shortcut.activated.connect(self.redo)
        self.redo_shortcut2 = QShortcut(QKeySequence("Ctrl+Shift+Z"), self)
        self.redo_shortcut2.activated.connect(self.redo)

        # Ctrl+S: Save Chain, Ctrl+E: Export Crease Pattern
        self.save_shortcut = QShortcut(QKeySequence.Save, self)
        self.save_shortcut.activated.connect(self.save_chain)
        self.export_shortcut = QShortcut(QKeySequence("Ctrl+E"), self)
        self.export_shortcut.activated.connect(self.save_crease_pattern)
        
        self.edit_grid_button = QPushButton("Edit Grid")
        self.edit_grid_button.clicked.connect(self.edit_grid_func)
        self.options_layout.addWidget(self.edit_grid_button)

        self.toggle_grid = QCheckBox("Grid Visibility")
        self.toggle_grid.setChecked(True)
        self.toggle_grid.toggled.connect(self.toggle_grid_func)
        self.options_layout.addWidget(self.toggle_grid)

        self.units_label = QLabel(f"Current units: {self.units}")
        self.units_label.setWordWrap(True)
        self.options_layout.addWidget(self.units_label)

        self.edit_dims_button = QPushButton("Edit Dimensions")
        self.edit_dims_button.clicked.connect(self.edit_dims_func)
        self.options_layout.addWidget(self.edit_dims_button)

        self.options_layout.addWidget(self.environment_widget)
        self.options_layout.addWidget(self.add_mesh_widget)

        self.options_widget.setLayout(self.options_layout)
        self.options_dock.setWidget(self.options_widget)
        #self.options_dock.setMaximumSize(300, 150)

        # ////////////////////////////////    CAMERA CONTROLS DOCK    ///////////////////////////////////
        self.camera_controls_dock = CollapsibleDockWidget("Camera Controls", self)
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
        self.camera_controls_dock.setWidget(self.camera_options_widget)
        
        #self.camera_controls_dock.setMaximumSize(300, 150)

        # ////////////////////////////////    FILE   ///////////////////////////////////
        file_dock = CollapsibleDockWidget("File", self)
        #file_dock.setAllowedAreas(Qt.RightDockWidgetArea)

        file_dock_widget = QWidget()
        file_dock_layout = QVBoxLayout(file_dock_widget)

        self.save_chain_button = QPushButton('Save')
        self.save_chain_button.clicked.connect(self.save_chain)
        file_dock_layout.addWidget(self.save_chain_button)

        self.load_chain_button = QPushButton('Load')
        self.load_chain_button.clicked.connect(self.load_chain)
        file_dock_layout.addWidget(self.load_chain_button)

        self.save_crease_pattern_button = QPushButton('Export Crease Pattern')
        self.save_crease_pattern_button.clicked.connect(self.save_crease_pattern)  
        file_dock_layout.addWidget(self.save_crease_pattern_button) 

        self.undo_button = QPushButton("Undo")
        self.undo_button.clicked.connect(self.undo)
        file_dock_layout.addWidget(self.undo_button)

        self.redo_button = QPushButton("Redo")
        self.redo_button.clicked.connect(self.redo)
        file_dock_layout.addWidget(self.redo_button)

        file_dock_widget.setLayout(file_dock_layout)
        file_dock.setWidget(file_dock_widget)

        # ////////////////////////////////    STL CONVERSION    ///////////////////////////////////
        # self.random_btn_dock = CollapsibleDockWidget("Export Options", self)
        # self.random_btn_widget = QWidget()
        # self.random_btn_layout = QVBoxLayout()
   
        # self.random_btn = QPushButton("Generate STL")
        # self.random_btn.clicked.connect(self.generate_stl)
        # file_dock_layout.addWidget(self.random_btn)

        # file_dock_layout.addWidget(self.random_btn)
        
        # self.random_btn_widget.setLayout(self.random_btn_layout)
        # self.random_btn_dock.setWidget(self.random_btn_widget)
        # self.addDockWidget(Qt.LeftDockWidgetArea, self.random_btn_dock)

        # self.random_btn_dock.setMaximumSize(300, 100)

        self.delete_joint_widget = DeleteWidget(self)
        self.delete_joint_dock = CollapsibleDockWidget("Confirm Delete", self)
        self.delete_joint_dock.setWidget(self.delete_joint_widget)
        self.delete_joint_dock.setVisible(False)
        
        self.add_chain_button_widget = AddChainWidget(self)
        self.add_chain_popup_dock = CollapsibleDockWidget("Create New Chain", self)
        self.add_chain_popup_dock.setWidget(self.add_chain_button_widget)
        self.add_chain_popup_dock.setVisible(False)  # Initially hidden

        self.edit_grid_widget = EditGridWidget(self)
        self.edit_grid_dock = CollapsibleDockWidget("Edit Grid", self)
        self.edit_grid_dock.setWidget(self.edit_grid_widget)
        self.edit_grid_dock.setVisible(False)

        self.edit_dims_widget = EditDimensionsWidget(self)
        self.edit_dims_dock = CollapsibleDockWidget("Edit Dimensions", self)
        self.edit_dims_dock.setWidget(self.edit_dims_widget)
        self.edit_dims_dock.setVisible(False)
        self.edit_dims_widget.target_units.connect(self.change_units)
        
        # ////////////////////////////////    MEASUREMENT    ///////////////////////////////////
        self.measurement_widget = MeasurementWidget(self)
        self.measurement_dock = CollapsibleDockWidget("Measure", self)
        self.measurement_dock.setWidget(self.measurement_widget)

        self.measurement_widget.measure_select_toggled.connect(self._toggle_measure_select)
        self.measurement_widget.frame_type_changed.connect(self._measure_frame_type_changed)
        self.measurement_widget.set_endpoint_active.connect(self._set_measure_active_endpoint)
        self.measurement_widget.clear_requested.connect(self._clear_measurement)
        self.measurement_widget.axis_frame_changed.connect(lambda: self.update_joint())
        self.plot_widget.measure_joint_selected.connect(self._measure_joint_selected_slot)
        self.measurement_dock.visibilityChanged.connect(self._on_measurement_dock_visibility)

        self.addDockWidget(Qt.TopDockWidgetArea, top_dock_widget)

        self.addDockWidget(Qt.LeftDockWidgetArea, file_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.edit_dims_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.edit_grid_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.options_dock)

        self.addDockWidget(Qt.LeftDockWidgetArea, self.add_chain_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.add_chain_popup_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.camera_controls_dock)
        
        self.addDockWidget(Qt.RightDockWidgetArea, self.add_chain_dock)
        self.addDockWidget(Qt.RightDockWidgetArea, self.add_chain_popup_dock)
        self.addDockWidget(Qt.RightDockWidgetArea, add_joints_dock)
        self.addDockWidget(Qt.RightDockWidgetArea, edit_joints_dock)
        self.addDockWidget(Qt.RightDockWidgetArea, self.delete_joint_dock)
        
        self.addDockWidget(Qt.BottomDockWidgetArea, self.configurations_dock)
        self.tabifyDockWidget(self.configurations_dock, self.measurement_dock)
        self.tabifyDockWidget(self.configurations_dock, self.message_display_dock)
        self.configurations_dock.raise_()

        # Make side docks take precedence at corners (extend full height)
        self.setCorner(Qt.TopLeftCorner, Qt.LeftDockWidgetArea)
        self.setCorner(Qt.TopRightCorner, Qt.RightDockWidgetArea)
        self.setCorner(Qt.BottomLeftCorner, Qt.LeftDockWidgetArea)
        self.setCorner(Qt.BottomRightCorner, Qt.RightDockWidgetArea)

        self.factor = {
            'Centimeter (cm)': 1.0,
            'Meter (m)': 100.0,
            'Milimeter (mm)': 10.0,
            'Inch (in)': 2.54,
            'Feet (ft)': 30.48
        }

    @property
    def model(self):
        return self.chain

    @model.setter
    def model(self, v):
        self.chain = v

    @property
    def model_created(self) -> bool:
        return self.chain_created

    def _post_undo_hook(self, prev_model) -> None:
        prev_units = prev_model.units if prev_model else None
        if self.model is not None:
            new_units = self.model.units
            self.units = new_units
        else:
            new_units = prev_units
        if prev_units != new_units:
            self.rescale_dimensions(prev_units, new_units)
            self.change_units(new_units)

    def _save_model(self, autosave_id=None) -> None:
        self.save_chain(autosave_id=autosave_id)

    @QtCore.pyqtSlot(str)
    def change_units(self, key):
        self.units = key
        self.units_label.setText(f"Current units: {self.units}")
        self.chain.units = key
        # self.log_version()
        self.update_joint()

    def rescale_dimensions(self, prev, new):
        if (prev != new):
            conversion_factor = self.factor[prev] / self.factor[new]
            self.plot_widget.opts['distance'] *= conversion_factor
            self.update_joint()

    def preserve_dimensions(self, prev, new):
        if (prev != new):
            conversion_factor = self.factor[prev] / self.factor[new]
            
            self.change_chain_size(conversion_factor)
    
    def change_chain_size(self, factor):
        new_radius = self.radius * factor

        if (self.chain_created):
            self.chain.changeRadius(new_radius)

            for joint in self.chain.Joints:
                joint.Pose.t *= factor

            new_chain = None

            for joint in self.chain.Joints:
                if (new_chain == None):
                    new_chain = KinematicChain(self.chain.Joints[0], units=self.units)
                else:
                    new_chain.append(joint, relative=False, fixedPosition=True,
                                            fixedOrientation=True, safe=False)
                
            self.chain = new_chain

            self.plot_widget.opts['distance'] *= factor

            self.log_version()
            self.update_joint()

    def edit_dims_func(self):
        visibility = self.edit_dims_dock.isVisible()
        self.edit_dims_dock.setVisible(not visibility)

    def onUpdateRadius(self, value):
        value = value / 10.0
        self.radius = value
        self.chain.changeRadius(value)
        # self.update_joint() # now called in joint_selection_changed
        self.joint_selection_changed(self.selected_joint, force=True)

    def delete_selected_joint(self):
        self.selected_joint = self.select_joint_options.currentIndex()
        deleted_index = self.selected_joint
        temp_last = len(self.chain.Joints) - 1
        if len(self.chain.Joints) == 1:
            self.chain = None
            self.chain_created = False
            self.plot_widget.clear()
            self.setCentralWidget(self.plot_widget)
            self.show_success('Joint successfully deleted!')
            self.last_joint = -1
            # Clear saved configurations when chain is deleted
            self.saved_configurations.clear()
        else:
            backup = copy.deepcopy(self.chain)
            try:
                # Update saved configs before deleting
                self.update_saved_configs_for_joint_deleted(deleted_index)
                self.chain = chainWithJointDeleted(self.chain, self.selected_joint)
                # self.reload_IDs()
                self.update_joint()
                self.show_success('Joint successfully deleted!')
            except Exception as e:
                print(e)
                self.chain = backup
                self.show_error('Error deleting joint.')
        self.window().delete_joint_dock.setVisible(False)

    def create_new_chain_func(self):
        # Show the AddChainWidget dock when this function is called
        self.add_chain_popup_dock.setVisible(True)
        self.add_chain_dock.setVisible(False)

    def create_new_chains(self, num_sides, radius):
        # Implement chain creation logic here
        self.chain = None
        self.radius = radius
        self.num_sides = num_sides
        self.chain_created = True
        self.selected_joint = -1
        self.plot_widget.clear()
        self.setCentralWidget(self.plot_widget)

        self.grid = gl.GLGridItem()
        self.grid.setColor(gridColorDefault)

        if self.grid_on:
            self.plot_widget.addItem(self.grid)

        if self.referenceMesh is not None:
            self.plot_widget.addItem(self.referenceMesh.mesh)

        self.update_joint()
        self.log_version()
        self.show_success('Chain created!')
        self.edit_dimension_menu.setVisible(False)
        self.edit_dimension_button.setVisible(True)

    def generate_stl(self):
        newTree = origamiToPrinted(self.chain, 0.05)

        self.plot_widget.clear()

        self.grid = gl.GLGridItem()
        self.grid.setColor((0,0,0,255))

        if self.grid_on:
            self.plot_widget.addItem(self.grid)

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

    def save_crease_pattern(self):
        if self.chain:
            options = QFileDialog.Options()
            try:
                base_path = sys._MEIPASS 
            except AttributeError:
                base_path = os.path.abspath(".")
            file_path, _ = QFileDialog.getSaveFileName(
                self, "Save File", os.path.join(base_path, "save"), "DXF Files (*.dxf)", options=options
            )
            crease_pattern = self.chain.creasePattern()
            # Render crease pattern into a standalone Qt window using matplotlib
            if file_path:
                crease_pattern.show(dxfName=file_path, show=False)
            self._show_crease_pattern_window(crease_pattern)

    def _show_crease_pattern_window(self, crease_pattern):
        """Display the crease pattern in a dedicated Qt dialog with an embedded matplotlib canvas."""
        import matplotlib
        matplotlib.use('Agg')  # Use non-interactive backend for rendering
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
        from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
        import ezdxf
        from ezdxf.addons.drawing import RenderContext, Frontend
        from ezdxf.addons.drawing.matplotlib import MatplotlibBackend
        from ezdxf.addons.drawing.config import Configuration
        from ezdxf.addons.drawing.properties import LayoutProperties
        from style import mapPropertiesColor

        # Build the DXF doc from the crease pattern (reuse save() logic without writing file)
        doc = crease_pattern.show(show=False)

        xmin = -0.5 * crease_pattern.baseSideLength
        xmax = crease_pattern.width + 0.5 * crease_pattern.baseSideLength
        ymin = crease_pattern.proximalMarker[0, 1]
        ymax = ymin + crease_pattern.patternHeight
        msp = doc.modelspace()

        # Render into a matplotlib figure
        config = Configuration()
        fig = plt.Figure(figsize=(max(xmax - xmin, 4), max(ymax - ymin, 4)))
        ax = fig.add_axes([0, 0, 1, 1])
        ctx = RenderContext(doc)
        out = MatplotlibBackend(ax, adjust_figure=False)
        msp_properties = LayoutProperties.from_layout(msp)
        msp_properties.set_colors(mapPropertiesColor)
        Frontend(ctx, out, config=config).draw_layout(msp, finalize=False,
                                                       layout_properties=msp_properties)
        ax.set_aspect("equal", adjustable="box")
        ax.set_ylim(ymin, ymax)
        ax.set_xlim(xmin, xmax)

        # Create a standalone Qt dialog to host the figure
        dialog = QDialog(self)
        dialog.setWindowTitle("Crease Pattern")
        dialog.setAttribute(Qt.WA_DeleteOnClose)
        dialog.resize(900, 700)
        layout = QVBoxLayout(dialog)
        canvas = FigureCanvas(fig)
        toolbar = NavigationToolbar(canvas, dialog)
        layout.addWidget(toolbar)
        layout.addWidget(canvas)
        canvas.draw()
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

        # Reset backend so future matplotlib usage isn't affected
        matplotlib.use('QtAgg')

    def _gather_session_state(self, autosave=False):
        """Collect all serializable editor state into a dict for .session files.
        When autosave=True, undo/redo history is excluded to keep the file small."""
        state = {
            'chain': self.chain,
            'saved_configurations': self.saved_configurations,
            'config_durations': self.config_durations,
            'selected_joint': self.selected_joint,
            'selected_frame': self.selected_frame,
            'grid_color': self.grid_color,
            'grid_spacing': self.grid_spacing,
            'grid_size': self.grid_size,
            'grid_on': self.grid_on,
            'units': self.units,
            'mesh_scale': self.mesh_scale,
            'control_type': self.control_type,
            'is_local': self.is_local,
            'animation_loop': self.animation_loop,
            'versions': [] if autosave else self.versions,
            'version_index': self.version_index,
            'total_version_counter': self.total_version_counter,
            'num_sides': self.num_sides,
            'radius': self.radius,
            'environment_visible': self.environment_visible,
            # Camera state
            'camera_distance': self.plot_widget.opts.get('distance', 15),
            'camera_elevation': self.plot_widget.opts.get('elevation', 30),
            'camera_azimuth': self.plot_widget.opts.get('azimuth', 45),
        }
        # Convert camera center from pyqtgraph.Vector to plain list for serialization
        center = self.plot_widget.opts.get('center', None)
        if center is not None:
            state['camera_center'] = [float(center.x()), float(center.y()), float(center.z())]
        # Reference mesh: save vertex/face data and pose if present
        if self.referenceMesh is not None:
            try:
                md = self.referenceMesh.mesh.opts.get('meshdata', None)
                if md is not None:
                    state['reference_mesh'] = {
                        'vertexes': md.vertexes(),
                        'faces': md.faces(),
                        'pose': self.referenceMesh.Pose,
                        'r': self.referenceMesh.r,
                    }
            except Exception:
                pass
        return state

    def _restore_session_state(self, state):
        """Restore editor state from a session dict."""
        self.chain = state.get('chain', None)
        self.saved_configurations = state.get('saved_configurations', [])
        self.config_durations = state.get('config_durations', [])
        self.selected_joint = state.get('selected_joint', -1)
        self.selected_frame = state.get('selected_frame', -1)
        self.grid_color = state.get('grid_color', gridColorDefault)
        self.grid_spacing = state.get('grid_spacing', 1.0)
        self.grid_size = state.get('grid_size', 10)
        self.grid_on = state.get('grid_on', True)
        self.units = state.get('units', 'Centimeter (cm)')
        self.mesh_scale = state.get('mesh_scale', 1.0)
        self.control_type = state.get('control_type', 'Translate')
        self.is_local = state.get('is_local', True)
        self.animation_loop = state.get('animation_loop', False)
        self.versions = state.get('versions', [])
        self.version_index = state.get('version_index', -1)
        self.total_version_counter = state.get('total_version_counter', 0)
        self.num_sides = state.get('num_sides', 4)
        self.radius = state.get('radius', 1.0)
        self.environment_visible = state.get('environment_visible', True)
        self.environment_items = []  # environment items are not serialized; re-import after load

        if self.chain is not None:
            self.chain_created = True

        # Restore camera
        self.plot_widget.opts['distance'] = state.get('camera_distance', self.grid_size * 1.5)
        self.plot_widget.opts['elevation'] = state.get('camera_elevation', 30)
        self.plot_widget.opts['azimuth'] = state.get('camera_azimuth', 45)
        center = state.get('camera_center', None)
        if center is not None:
            self.plot_widget.opts['center'] = pg.Vector(*center)

        # Restore grid
        self.grid.setColor(self.grid_color)
        self.grid.setSize(self.grid_size, self.grid_size, self.grid_size)
        self.grid.setSpacing(self.grid_spacing, self.grid_spacing, self.grid_spacing)

        # Restore units label
        self.units_label.setText(f"Current units: {self.units}")

        # Restore control type radio buttons
        if self.control_type == 'Translate':
            self.translate_joint_radio_button.setChecked(True)
        else:
            self.rotate_joint_radio_button.setChecked(True)
        self.local_orient_slider_checkbox.setChecked(self.is_local)

        # Restore reference mesh if present
        if 'reference_mesh' in state:
            try:
                rm = state['reference_mesh']
                meshdata = gl.MeshData(vertexes=rm['vertexes'], faces=rm['faces'])
                mesh_item = gl.GLMeshItem(meshdata=meshdata, smooth=True, shader='shaded')
                mesh_item.setObjectName("Mesh")
                self.referenceMesh = ReferenceMesh(mesh=mesh_item)
                self.referenceMesh.Pose = rm['pose']
                self.referenceMesh.r = rm['r']
            except Exception:
                self.referenceMesh = None
        else:
            self.referenceMesh = None

    def save_chain(self, autosave_id=None):
        # confusing why autosave_id is sometimes False
        if autosave_id is False:
            autosave_id = None
        if self.chain:
            if autosave_id is None:
                print("save dialog")
                options = QFileDialog.Options()
                try:
                    base_path = sys._MEIPASS 
                except AttributeError:
                    base_path = os.path.abspath(".")  
                file_path, _ = QFileDialog.getSaveFileName(
                    self, "Save File", os.path.join(base_path, "save"),
                    "Session Files (*.session);;Tree Files (*.tree);;Chain Files (*.chain)",
                    options=options
                )
            else:
                try:
                    base_path = sys._MEIPASS
                except AttributeError:
                    base_path = os.path.abspath(".")
                file_path = os.path.join(base_path, "save", "autosave", f"autosave_{autosave_id}.session")

            if file_path:
                if file_path.endswith('.session'):
                    is_autosave = autosave_id is not None
                    session_state = self._gather_session_state(autosave=is_autosave)
                    if is_autosave:
                        def _do_autosave(state, path):
                            try:
                                with open(path, 'wb') as f:
                                    dill.dump(state, f)
                            except Exception as e:
                                print(f"[autosave] error: {e}")
                        _threading.Thread(
                            target=_do_autosave,
                            args=(session_state, file_path),
                            daemon=True,
                        ).start()
                    else:
                        try:
                            with open(file_path, 'wb') as f:
                                dill.dump(session_state, f)
                            self.show_success(f'Session saved to {os.path.basename(file_path)}')
                        except Exception as e:
                            self.show_error(f'Error saving session: {e}')
                else:
                    self.chain.save(file_path)
        
    def load_chain(self):
        options = QFileDialog.Options()
        try:
            base_path = sys._MEIPASS
        except AttributeError:
            base_path = os.path.abspath(".")
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open File", os.path.join(base_path, "save"),
            "Session Files (*.session);;Tree Files (*.tree);;Chain Files (*.chain);;All Files (*.*)",
            options=options
        )
        if file_path:
            if file_path.endswith('.session'):
                try:
                    with open(file_path, 'rb') as f:
                        state = dill.load(f)
                    self._restore_session_state(state)
                    self.show_success(f'Session loaded from {os.path.basename(file_path)}')
                except Exception as e:
                    self.show_error(f'Error loading session: {e}')
                    return
            else:
                self.chain = loadTree(file_path)
                self.radius = self.chain.r
                self.num_sides = self.chain.numSides
                self.chain_created = True
                
                self.units = self.chain.units
                self.units_label.setText(f"Current units: {self.units}")

            self.update_joint()
            if not file_path.endswith('.session'):
                self.log_version()

    def done_transforming(self, done):
        if done:
            self.log_version()

    """
    def update_slider(self, slider_type):
        if (slider_type == "rotation"):
            slider = self.rotation_slider
            textbox = self.rotation_textbox

            slider.setMinimum(-360)
            slider.setMaximum(360)

            rotation_matrix = self.chain.Joints[self.selected_joint].Pose.R
            angle_degrees = self.rotation_angle_from_matrix(rotation_matrix, self.selected_arrow)
            slider_value = int(angle_degrees)
            slider_text = str(int(angle_degrees))

            set_slider = self.selected_arrow != -1

            if (set_slider):
                self.old_rot_val = angle_degrees

        elif (slider_type == "translation"):
            slider = self.translation_slider
            textbox = self.translation_textbox

            amount = self.chain.Joints[self.selected_joint].Pose.t[self.selected_arrow]
            slider_value = int(amount * 10) #TODO: should this adjust for radius?
            slider_text = str(slider_value)

            set_slider = self.selected_arrow != -1

            if (set_slider):
                self.old_trans_val = amount

        elif (slider_type == "state"):
            slider = self.state_slider
            textbox = self.state_input

            min = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[0])
            max = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[1])

            current = int(math.degrees(self.chain.Joints[self.selected_joint].state))
            slider_value = current
            slider_text = str(current)

            set_slider = self.selected_joint != -1 and min != 0 and max != 0

            if (set_slider):
                self.state_slider.setMinimum(int(min))
                self.state_slider.setMaximum(int(max))
                self.old_state_slider_val = current

        if not set_slider:
            slider_value = 0
            slider_text = ""

        slider.blockSignals(True)
        slider.setValue(slider_value)
        slider.setDisabled(False)
        slider.blockSignals(False)
        textbox.blockSignals(True)
        textbox.setText(slider_text)
        textbox.blockSignals(False)
        textbox.setDisabled(False)    

    def update_state_slider(self):
        self.update_slider("state")
    """

    def add_chain(self, chain):
        self.chain = chain
        self.select_joint_options.blockSignals(True)
        self.select_joint_options.clear()
    
        for i, joint in enumerate(self.chain.Joints):
            self.select_joint_options.addItem("Joint " + str(i) + " - " + joint.__class__.__name__)
    
        self.select_joint_options.blockSignals(False) 
        self.select_joint_options.setCurrentIndex(self.selected_joint)

        self.select_link_options.blockSignals(True)
        self.select_link_options.clear()

        for link in enumerate(self.chain.Links): 
            self.select_link_options.addItem("Link " + str(i) + " - " + link.__class__.__name__)
        
        self.select_link_options.blockSignals(False)
        self.select_link_options.setCurrentIndex(self.selected_link)

    import copy

    def edit_joint_dimension(self):
        self._backup_chain = copy.deepcopy(self.chain)
        self._saved_states = []

        for idx, joint in enumerate(self._backup_chain.Joints):
            original_state = joint.state
            self._saved_states.append(original_state)

            self._backup_chain.setJointState(idx, 0)

        self._selected_joint = self.selected_joint

        target_joint = self._backup_chain.Joints[self._selected_joint]
        self.prev_joint = self._backup_chain.Joints[self._selected_joint - 1] if self._selected_joint > 0 else None

        if target_joint.__class__.__name__ in ["Prismatic", "OrigamiPrismatic"]:
            self.edit_dimension_menu.updatePrismatic()
        elif target_joint.__class__.__name__ in ["Revolute", "TransverseRevolute", "CoaxialRevolute", "OrigamiRevolute", "OrigamiExtendedRevolute"]:
            self.edit_dimension_menu.updateRevolute()
        elif target_joint.__class__.__name__ in ["StartTip", "EndTip", "Tip"]:
            self.edit_dimension_menu.updateTip()
        else:
            self.show_error("Uneditable joint type.")
            return

        self.edit_dimension_toggle()

    def finish_joint_edit(self, new_joint):
        try:
            new_chain = None
            for idx, joint in enumerate(self._backup_chain.Joints):
                if idx == self._selected_joint:
                    min_val, max_val = new_joint.stateRange()
                    original_state = self._saved_states[idx]
                    if original_state < min_val:
                        original_state = min_val
                    elif original_state > max_val:
                        original_state = max_val
                    self._saved_states[idx] = original_state

                    new_joint.Pose = joint.Pose
                    if new_chain is None:
                        new_chain = OrigamiKinematicChain(new_joint, numSides=self.num_sides, units=self.units)
                    else:
                        new_chain.append(new_joint, relative=False,
                                        fixedPosition=True, fixedOrientation=True, safe=False)
                else:
                    if new_chain is None:
                        new_chain = OrigamiKinematicChain(joint, numSides=self.num_sides, units=self.units)
                    else:
                        new_chain.append(joint, relative=False,
                                        fixedPosition=True, fixedOrientation=True, safe=False)
            
            for idx, joint in enumerate(new_chain.Joints):
                new_chain.setJointState(idx, self._saved_states[idx])
            
            self.chain = new_chain
            self.selected_joint = self._selected_joint
            self.window().log_version()
            self.show_success("Chain updated successfully!")
        except Exception as e:
            self.chain = self._backup_chain
            self.show_error("Error rebuilding chain: " + str(e))
        
        self.update_joint(force_recreate_config_widget=True)

    def edit_joint_state(self):
        dialog = EditJointStateDialog(self) 
        if not self.chain:
            self.show_error('Please initialize a chain.')
            # error_dialog = ErrorDialog('Please initialize a chain.')
            # error_dialog.exec_()
        if self.selected_joint == -1:
            self.show_error('Please select a joint.')
            # error_dialog = ErrorDialog('Please select a joint.')
            # error_dialog.exec_()
        elif dialog.exec_() == QDialog.Accepted:
            edit = dialog.get_state()
            if edit is not None:

                if self.chain.setJointState(self.selected_joint, math.radians(edit)):
                    self.update_joint()
                    self.show_success('Joint state successfully edited!')
                    # success_dialog = SuccessDialog('Joint state successfully edited!')
                    # success_dialog.exec_()
                    min = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[0])
                    max = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[1])
                    current = math.degrees(self.chain.Joints[self.selected_joint].state)
                    self.current_state_label.setText(f"Min State: {int(min)} ≤ Current State: {int(current)} ≤ Max State: {int(max)}")
                else:
                    self.show_error('Error editing joint state.')
                    # error_dialog = ErrorDialog('Error editing joint state.')
                    # error_dialog.exec_()

    def adjust_rotation(self, value):
        # if not isinstance(value, float):
        #     value = value.strip()
        value = float(value) if len(str(value)) > 0 else 0
        angle_radians = math.radians(value - self.old_rot_val)
        if self.chain and self.selected_joint != -1:
            if self.selected_arrow == 0:
                transformation = SE3.Rx(angle_radians)
            elif self.selected_arrow == 1:
                transformation = SE3.Ry(angle_radians)
            else:
                transformation = SE3.Rz(angle_radians)
            propogate = self.propogate_slider_checkbox.isChecked()
            localOrient = self.local_orient_slider_checkbox.isChecked()
            if self.chain.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=True, localOrient=localOrient):
                self.old_rot_val = int(value)
                self.update_joint()             
                self.rotation_slider.blockSignals(True)
                self.rotation_slider.setDisabled(False)
                self.rotation_slider.blockSignals(False)
                self.rotation_textbox.blockSignals(True)
                self.rotation_textbox.setText(str(value))
                self.rotation_textbox.blockSignals(False)
                self.rotation_textbox.setDisabled(False)
            else:
                self.rotation_slider.blockSignals(True)
                self.rotation_slider.setValue(int(self.old_rot_val))
                self.rotation_slider.blockSignals(False)

    def adjust_translation(self, value):
        if not isinstance(value, float) and not isinstance(value, int):
            value = value.strip()
        value = float(value) if value else 0
        actualVal = value / 10
        amount = actualVal - self.old_trans_val
        if self.chain and self.selected_joint != -1:
            propogate = self.propogate_slider_checkbox.isChecked()
            localOrient = self.local_orient_slider_checkbox.isChecked()
            transformation = SE3()
            if (self.selected_arrow == 0):
                transformation = SE3.Tx(amount)
            if (self.selected_arrow == 1):
                transformation = SE3.Ty(amount)
            if (self.selected_arrow == 2):
                transformation = SE3.Tz(amount)
            if self.chain.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=True, localOrient=localOrient):
                self.old_trans_val = actualVal
                self.update_joint()                
                self.translation_slider.blockSignals(True)
                self.translation_slider.setDisabled(False)
                self.translation_slider.blockSignals(False)
                self.translation_textbox.blockSignals(True)
                self.translation_textbox.setText(str(actualVal))
                self.translation_textbox.blockSignals(False)
                self.translation_textbox.setDisabled(False)
            else:
                self.translation_slider.blockSignals(True)
                self.translation_slider.setValue(int(self.old_trans_val * 10))
                self.translation_slider.blockSignals(False)
    
    """
    def adjust_state(self, stateTextboxValue):
        if not isinstance(stateTextboxValue, float) and not isinstance(stateTextboxValue, int):
            stateTextboxValue = stateTextboxValue.strip()
        stateFromTextbox = float(stateTextboxValue) if stateTextboxValue else 0
        scaledInfo = self.scaled_state_info(stateFromTextbox)
        if not scaledInfo is None:
            actual, slider, textbox = scaledInfo
            actualState = actual[2]
            if self.chain.setJointState(self.selected_joint, actualState):
                self.update_joint()
                self.current_state_label.setText(f"Min State: {textbox[0]} ≤ Current State: {textbox[2]} ≤ Max State: {textbox[1]}")
                self.old_state_slider_val = slider[2]
                self.state_slider.blockSignals(True)
                self.state_slider.setValue(slider[2])
                self.state_slider.setDisabled(False)
                self.state_slider.blockSignals(False)
                self.state_textbox.blockSignals(True)
                self.state_textbox.setText(textbox[2])
                self.state_textbox.blockSignals(False)
                self.state_textbox.setDisabled(False)
            else:
                self.state_slider.blockSignals(True)
                self.state_slider.setValue(int(self.old_state_slider_val))
    """

    def delete_joint(self):
        # dialog = DeleteDialog(self)
        if not self.chain:
            self.show_error('Please initialize a chain.')
            # error_dialog = ErrorDialog('Please initialize a chain.')
            # error_dialog.exec_()
        if self.selected_joint == -1:
            self.show_error('Please select a joint.')
            # error_dialog = ErrorDialog('Please select a joint.')
            # error_dialog.exec_()
        else:
            self.show_delete_widget()
    """
    def reload_IDs(self): 
        if self.chain is not None:
            for index, joint in enumerate(self.chain.Joints):
                joint.id = index
    """

    def update_joint(self, force_recreate_config_widget : bool = False):
        # If link geometry is stale from lightweight updates, rebuild first
        if self._lightweight_dirty and self.chain is not None:
            self.chain.resyncFromLightweight()
            self._lightweight_dirty = False

        self.edit_dimension_menu.setVisible(False)
        self.edit_dimension_button.setVisible(True)
        self.select_joint_options.blockSignals(True)
        self.select_link_options.blockSignals(True)

        self.units_label.setText(f"Current units: {self.units}")
        
        # Check if we need to recreate config widgets or just update values.
        # Compare actual real-joint indices, not just count, so that operations
        # that shift indices (add-to-root, waypoint insert, undo/redo, load)
        # are detected even when the number of real joints stays the same.
        need_recreate_config_widget = force_recreate_config_widget
        if not need_recreate_config_widget:
            if self.chain is None or not self.chain_created:
                need_recreate_config_widget = len(getattr(self, 'config_joint_indices', [])) > 0
            elif not hasattr(self, 'config_joint_indices'):
                need_recreate_config_widget = True
            else:
                current_real_indices = [i for i, j in enumerate(self.chain.Joints)
                                        if type(j).__name__ not in ['Waypoint', 'PrintedWaypoint']]
                if current_real_indices != self.config_joint_indices:
                    need_recreate_config_widget = True
        
        if need_recreate_config_widget:
            self.update_configurations()
        else:
            self.update_config_values()

        if (not self.stl_generated):
            self.plot_widget.clear()
            self._measurement_line_item = None
            self.select_joint_options.clear()
            self.select_link_options.clear()
            self.setCentralWidget(self.plot_widget)

            if self.grid_on:
                self.plot_widget.addItem(self.grid)

            if self.mesh_visible and self.referenceMesh is not None:
                self.plot_widget.addItem(self.referenceMesh.mesh)

            if self.environment_visible:
                for env_item in self.environment_items:
                    self.plot_widget.addItem(env_item)

            if self.chain is not None:
                self.chain.addToWidget(
                    self,
                    selectedJoint=self.selected_joint,
                    selectedLink=self.selected_link,
                    lastJoint=self.last_joint
                )
                self.add_chain(self.chain)

        self._refresh_measurement_overlay()
        # ─────────────────────────────────────────────────────────────────────

        if self.mesh_selected and self.referenceMesh is not None:
            if self.control_type == "Translate":
                self.referenceMesh.addTranslateArrows(
                    self,
                    selectedArrow=self.selected_arrow,
                    local=self.is_local
                )
            else:
                self.referenceMesh.addRotateArrows(
                    self,
                    selectedArrow=self.selected_arrow,
                    local=self.is_local
                )

        if self.chain is not None and self.selected_joint != -1:
            joint = self.chain.Joints[self.selected_joint]
            frame_pose = None
            if self.selected_frame >= 0:
                frame_pose = self.chain.Joints[self.selected_frame].Pose

            if self.control_type == "Translate":
                joint.addTranslateArrows(
                    self,
                    selectedArrow=self.selected_arrow,
                    local=self.is_local,
                    frame=frame_pose
                )
            else:
                joint.addRotateArrows(
                    self,
                    selectedArrow=self.selected_arrow,
                    local=self.is_local,
                    frame=frame_pose
                )

        if self.selected_arrow != -1:
            self.rotation_slider.setDisabled(False)
            self.translation_slider.setDisabled(False)
        else:
            self.rotation_slider.setDisabled(True)
            self.translation_slider.setDisabled(True)

        #print("current radius: " + str(self.chain.r))
                
    # def add_joint(self, dialog):
    #     if dialog.exec_() == QDialog.Accepted:
    #         joint : Joint = dialog.getJoint()
            # if not self.add_to_root:
            #     if (self.chain == None or len(self.chain.Joints) == 0) :
            #         self.chain = KinematicChain(joint)
            #     else :
            #         self.chain.append(joint, relative=True, fixedPosition=True, fixedOrientation=False, safe=False)
            #         #self.chain.addJoint(self.selected_joint, joint, relative=True, fixedPosition=True, fixedOrientation=False, safe=False)
            #     self.selected_joint = len(self.chain.Joints) - 1
            # else:

            #     if (self.chain == None or len(self.chain.Joints) == 0) :
            #         self.chain = KinematicChain(joint)
            #     else:
            #         old_root = self.chain.Joints[0]
            #         joint.Pose = old_root.Pose @ joint.Pose

            #         new_chain = KinematicChain(joint)

            #         for i, jt in enumerate(self.chain.Joints):
            #             if i == 0:
            #                 cachedLink = None
            #             else:
            #                 cachedLink = self.chain.Links[i]
            #             new_chain.append(jt, relative=False, fixedPosition=True, fixedOrientation=True, safe=False, cachedLink=cachedLink)
            #         self.chain = new_chain

            #     self.selected_joint = 0

            # self.update_joint()
            # self.log_version()
    def add_joint(self, joint : Joint):
        # Check if this is a real joint (not a waypoint)
        is_real_joint = type(joint).__name__ not in ['Waypoint', 'PrintedWaypoint']
        
        if (self.chain == None or len(self.chain.Joints) == 0):
            # No chain exists, create one with this joint as root
            self.chain = OrigamiKinematicChain(joint, numSides=self.num_sides, units=self.units)
            self.selected_joint = 0 if self.add_to_root else len(self.chain.Joints) - 1
        elif self.add_to_root:
            # Add as new root: create a new chain with this joint, then add old chain as subtree
            self.add_joint_as_new_root(joint)
        else:
            self.chain.append(joint, relativeToDistalDubins=True, fixedPosition=True, fixedOrientation=True, safe=False)
            self.selected_joint = len(self.chain.Joints) - 1

        # Update saved configurations if a real joint was added
        if is_real_joint:
            new_joint_index = self.selected_joint
            real_indices = [i for i, j in enumerate(self.chain.Joints)
                           if type(j).__name__ not in ['Waypoint', 'PrintedWaypoint']]
            config_pos = real_indices.index(new_joint_index) if new_joint_index in real_indices else len(real_indices) - 1
            self.update_saved_configs_for_joint_added(config_pos)

        self.update_joint()
        self.log_version()
        self.joint_selection_changed(self.selected_joint, force=True)
    
    def add_joint_as_new_root(self, joint : Joint):
        """Add a joint as the new root, making the old chain a subtree of the new root.
        
        The new joint's pose should already be in global coordinates.
        """
        old_chain = self.chain
        new_chain = OrigamiKinematicChain(joint, numSides=self.num_sides, units=self.units)
        new_chain.addSubtree(0, old_chain)
        self.chain = new_chain
        self.selected_joint = 0
    
    # def add_joint_func(self, joint_type):
    #     numSides = self.num_sides

    #     if (not self.chain_created):
    #         self.chain_not_created()
    #         return

    #     if (self.chain and len(self.chain.Joints) > 0):
    #         if not self.add_to_root:
    #             prevJoint = self.chain.Joints[-1]
    #             add_to_root = False
    #         else:
    #             prevJoint = self.chain.Joints[0]
    #             add_to_root = True
    #     else: 
    #         prevJoint = None
    #         add_to_root = False

    #     if joint_type == "prismatic":
    #         dialog = AddPrismaticDialog(numSides, self.radius, prevJoint=prevJoint, add_to_root=add_to_root)
    #     elif joint_type == "revolute":
    #         dialog = AddRevoluteDialog(numSides, self.radius, prevJoint=prevJoint, add_to_root=add_to_root)
    #     elif joint_type == "tip":
    #         dialog = AddTipDialog(numSides, self.radius, prevJoint=prevJoint, add_to_root=add_to_root)

    #     self.add_joint(dialog)

    # def add_prismatic_func(self):
    #     self.add_prismatic_menu.setEnabled(True)
    #     self.add_prismatic.setEnabled(False)

    # def add_revolute_func(self):
    #     self.add_joint_func("revolute")

    # def add_tip_func(self):
    #     self.add_joint_func("tip")

    def add_prismatic_toggle(self):
        if (not self.chain_created):
            self.chain_not_created()
            return
        self.add_prismatic_menu.setVisible(not self.add_prismatic_menu.isVisible())
        self.add_prismatic.setVisible(not self.add_prismatic.isVisible())

    def add_revolute_toggle(self):
        if (not self.chain_created):
            self.chain_not_created()
            return
        self.add_revolute_menu.setVisible(not self.add_revolute_menu.isVisible())
        self.add_revolute.setVisible(not self.add_revolute.isVisible())

    def add_tip_toggle(self):
        if (not self.chain_created):
            self.chain_not_created()
            return
        self.add_tip_menu.setVisible(not self.add_tip_menu.isVisible())
        self.add_tip.setVisible(not self.add_tip.isVisible())

    def edit_dimension_toggle(self):
        if (not self.chain_created):
            self.chain_not_created()
            return
        self.edit_dimension_menu.setVisible(not self.edit_dimension_menu.isVisible())
        self.edit_dimension_button.setVisible(not self.edit_dimension_button.isVisible())

    def add_waypoint_func(self):
        numSides = self.num_sides
        if (not self.chain_created):
            self.chain_not_created()
        elif (self.selected_link != -1):
            link = self.chain.Links[self.selected_link]
            nextJoint = self.chain.Links[self.selected_link]
            
            pos = link.cylinder.start + 0.5 * link.cylinder.length * link.cylinder.direction
            newPos = SE3().Trans(x=pos[0], y=pos[1], z=pos[2])

            rot = link.cylinder.orientation()
            newRot = SE3(rot)

            waypoint = Waypoint(self.radius, newPos * newRot)
            waypoint_index = len(self.chain.Joints)
            
            self.chain.append(newJoint = waypoint, 
                                relative=False, fixedPosition=True, fixedOrientation=False, safe=False)
            
            nextJoint_index = self.chain.Joints.index(nextJoint)

            self.chain.Links[nextJoint_index] = LinkCSC(self.chain.r, waypoint.DistalDubinsFrame(), 
                                            nextJoint.ProximalDubinsFrame(),
                                            self.chain.maxAnglePerElbow)
            self.chain.Parents[nextJoint_index] = waypoint_index
            self.chain.Children[waypoint_index].append(nextJoint_index)
            self.update_joint()
            self.log_version()

        else: 
            if (self.chain == None) or len(self.chain.Joints) == 0:
                waypoint = Waypoint(self.radius, SE3())
                self.chain = OrigamiKinematicChain(waypoint, numSides=self.num_sides, units=self.units)
            elif self.add_to_root:
                # Compute pose in global coordinates behind the old root
                root = self.chain.Joints[0]
                root_proximal_dubins = root.ProximalDubinsFrame()
                r = root.r
                distance = 4 * r  # Waypoint neutralLength is 0
                # Waypoint has pathIndex=2, so rotate by Ry(pi/2) so z-hat aligns with dubins x-hat
                pose = root_proximal_dubins @ SE3.Rt(SE3.Ry(np.pi/2).R, np.array([-distance, 0, 0]))
                waypoint = Waypoint(self.radius, pose)
                self.add_joint_as_new_root(waypoint)
            else:
                prevJoint = self.chain.Joints[-1]
                # Calculate pose relative to distal Dubins frame of previous joint
                distance = 4 * self.radius
                pose = SE3.Rt(SE3.Ry(np.pi/2).R, np.array([distance, 0, 0]))
                waypoint = Waypoint(self.radius, pose)
                self.chain.append(newJoint = waypoint, relativeToDistalDubins=True, 
                                  fixedPosition=True, fixedOrientation=True, safe=False)

            self.update_joint()
            self.log_version()
            if self.add_to_root:
                self.select_joint_options.setCurrentIndex(0)
            else:
                self.select_joint_options.setCurrentIndex(len(self.chain.Joints) - 1)

    def insert_waypoint_func(self):
        if self.chain is None or len(self.chain.Joints) == 0:
            self.show_error("No chain to insert into.")
            return
        if self.selected_joint == -1:
            self.show_error("Please select a joint first.")
            return

        num = self.insert_waypoint_count.value()

        try:
            if self.selected_joint == 0:
                self._insert_waypoints_before_root(num)
            else:
                self.chain.insertWaypointsIntoLink(self.selected_joint, num)
                self.update_joint()
                self.log_version()
                self.joint_selection_changed(self.selected_joint, force=True)
        except Exception as e:
            self.show_error(str(e))

    def _insert_waypoints_before_root(self, numWaypoints: int):
        root = self.chain.Joints[0]
        r = root.r
        original_root_frame = root.ProximalDubinsFrame()
        eps = r * 0.01
        for i in range(numWaypoints, 0, -1):
            d = i * eps
            pose = original_root_frame @ SE3.Rt(SE3.Ry(np.pi/2).R, np.array([-d, 0, 0]))
            waypoint = Waypoint(r, pose)
            self.add_joint_as_new_root(waypoint)
        self.update_joint()
        self.log_version()
        self.select_joint_options.setCurrentIndex(0)

    def chain_not_created(self):
        self.show_error('Please create a chain first.')
        # error_dialog = ErrorDialog('Please create a chain first.')
        # error_dialog.exec_()
        self.create_new_chain_func()

    def is_parent_joint_selected(self):
        if self.selected_joint == -1 and self.chain is not None:
            self.show_error("Please select a parent joint or link first.")
            # QMessageBox.warning(self, "Selection Required", "Please select a parent joint or link before adding a new one.")
            return False
        return True

    def keyPressEvent(self, event: QKeyEvent):
        if event.key() == Qt.Key_T:
            self.control_type = "Translate"
            self.translate_joint_radio_button.setChecked(True)
            self.update_joint()
        elif event.key() == Qt.Key_R:
            self.control_type = "Rotate"
            self.rotate_joint_radio_button.setChecked(True)
            self.update_joint()
        elif event.key() == Qt.Key_Delete or event.key() == Qt.Key_Backspace:
            if (self.chain and self.selected_joint != -1):
                self.delete_selected_joint()
        elif event.key() == Qt.Key_X:
            self.arrow_selection_changed(0)
        elif event.key() == Qt.Key_Y:
            self.arrow_selection_changed(1)
        elif event.key() == Qt.Key_Z:
            self.arrow_selection_changed(2)
        elif event.key() == Qt.Key_G:
            self.toggle_grid_func()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = WindowKinegamiGUI()
    window.show()
    sys.exit(app.exec_())