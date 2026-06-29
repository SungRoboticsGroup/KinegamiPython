import sys, os
import dill
import time as _time_module
import cProfile as _cProfile
import pstats as _pstats
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
from PrintedTube import *
from KinematicTree import loadTree
import re
from scipy.spatial.transform import Rotation as R
from style import *
from ReferenceMesh import *
from Dialog import *
from printedJointWidget import *
from IntersectionHelper import *
from measurementWidget import MeasurementWidget

from typing import Optional

import warnings
warnings.filterwarnings("ignore")

if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    PyQt5.QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    PyQt5.QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

from guiBase import DeleteWidget, AddMeshWidget, EditGridWidget, CollapsibleDockWidget, ImageRadioButton, OverlayLine, BaseClickableGLViewWidget, BaseWindowKinegamiGUI, EnvironmentWidget


class ClearTreeWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Clear Tree')

        layout = QVBoxLayout()

        # Confirmation message
        confirm_label = QLabel("Are you sure you want to clear the current tree?")
        confirm_label.setWordWrap(True)
        layout.addWidget(confirm_label)

        # Apply button to clear the tree
        self.clear_button = QPushButton('Clear Tree', self)
        self.clear_button.clicked.connect(self.on_clear_clicked)
        layout.addWidget(self.clear_button)

        # Cancel button to close the widget
        self.cancel_button = QPushButton('Cancel', self)
        self.cancel_button.clicked.connect(self.on_cancel_clicked)
        layout.addWidget(self.cancel_button)

        self.setLayout(layout)

    def on_clear_clicked(self):
        self.window().clear_tree()
        self.window().clear_tree_popup_dock.setVisible(False)
        self.window().clear_tree_dock.setVisible(True)

    def on_cancel_clicked(self):
        # Hide the widget if the user cancels
        self.window().clear_tree_popup_dock.setVisible(False)
        self.window().clear_tree_dock.setVisible(True)

    def show_error(self, message):
        QMessageBox.warning(self, "Invalid Input", message)

class ClickableGLViewWidget(BaseClickableGLViewWidget):
    def get_world_coordinates(self, event):
        pos = event.localPos()
        ndc_x = (2.0 * pos.x()) / self.width() - 1.0
        ndc_y = 1.0 - (2.0 * pos.y()) / self.height()

        view = self.viewMatrix()
        proj = self.projectionMatrix()
        inverted_matrix = (proj * view).inverted()[0]

        # Use camera position as ray origin for better precision when far away
        camera_pos = self.cameraPosition()

        ndc = QVector3D(ndc_x, ndc_y, 0.0)  # Use z=0 (middle of NDC depth range)
        world_point = inverted_matrix.map(ndc)

        direction = world_point - camera_pos
        direction.normalize()

        return camera_pos, direction

    def _should_pick_items(self) -> bool:
        # itemsAt does a full GL_SELECT re-render (~150ms) just to detect mesh clicks.
        # Skip it entirely when no reference mesh is loaded — there's nothing to detect.
        return self.parent_window.referenceMesh is not None

    def _mpe_profiling_context(self):
        _prof = getattr(self.parent_window, '_profiling_enabled', False)
        return _time_module.perf_counter() if _prof else None

    def _mpe_profiling_log(self, ctx, picked_items, skipped: bool) -> None:
        _prof = getattr(self.parent_window, '_profiling_enabled', False)
        if not _prof:
            return
        if skipped:
            print(f"[PROF] mousePressEvent itemsAt=SKIPPED (no mesh loaded)")
        else:
            _items_at_ms = (_time_module.perf_counter() - ctx) * 1000
            print(f"[PROF] mousePressEvent itemsAt={_items_at_ms:.1f}ms  items={len(picked_items)}")

    def _on_drag_release(self) -> None:
        self.parent_window._last_drag_frame_t = 0.0

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
                    selected_joint = self.parent_window.tree.Joints[self.parent_window.selected_joint]
                joint_center = selected_joint.Pose.t
                qsphere_start = QVector3D(joint_center[0], joint_center[1], joint_center[2])

                trans = new_pos_3D - qsphere_start
                transformation = SE3.Trans(trans[0], trans[1], trans[2])

                propogate = self.parent_window.propogate_slider_checkbox.isChecked()

                if (self.mesh_selected):
                    self.parent_window.referenceMesh.transform(transformation)
                else:
                    # Throttle: skip stale queued events to prevent Qt event-queue
                    # backlog.  If the user drags for 2 s at 100 events/s and each
                    # frame takes 25 ms, 200 events queue up → 5 s freeze on release.
                    # Skipped events cost ~0.001 ms (just a timestamp check), so the
                    # backlog drains in <1 ms instead of seconds.
                    _now = _time_module.perf_counter()
                    if _now - getattr(self.parent_window, '_last_drag_frame_t', 0.0) < 0.200:
                        return
                    self.parent_window._last_drag_frame_t = _now

                    # Start cProfile on first processed drag frame if armed.
                    if getattr(self.parent_window, '_drag_profile_armed', False):
                        if getattr(self.parent_window, '_drag_profiler', None) is None:
                            self.parent_window._drag_profiler = _cProfile.Profile()
                        self.parent_window._drag_profiler.enable()

                    _prof = getattr(self.parent_window, '_profiling_enabled', False)
                    _drag_t0 = _time_module.perf_counter() if _prof else None

                    # Save backup before first drag transform (for reversion on error)
                    if not hasattr(self, '_drag_backup') or self._drag_backup is None:
                        self._drag_backup = self.parent_window.tree.dataDeepCopy()
                    if _prof: print(f"[PROF drag-trans] backup1: {(_time_module.perf_counter()-_drag_t0)*1000:.1f}ms")

                    # Capture old link GL items NOW — transformJoint replaces the
                    # Link object entirely, so after the call _gl_items is gone.
                    _sidx = self.parent_window.selected_joint
                    _old_link_items = list(getattr(
                        self.parent_window.tree.Links[_sidx], '_gl_items', []))

                    # Real transform, no bounding-ball recompute per frame.
                    # Incremental GL update rebuilds only the moved joint+link.
                    try:
                        self.parent_window.tree.transformJoint(
                            _sidx, transformation,
                            propogate=propogate, relative=False, safe=False,
                            lightweight=False, recomputeBoundingBall=False)
                    except Exception as _drag_err:
                        if _prof: print(f"[PROF drag-trans] transformJoint failed: {_drag_err}")
                        if self._drag_backup is not None:
                            self.parent_window.tree.setTo(self._drag_backup)
                        self.parent_window.update_joint()
                        return
                    if _prof: print(f"[PROF drag-trans] transformJoint: {(_time_module.perf_counter()-_drag_t0)*1000:.1f}ms")

                if not self.parent_window._incremental_drag_update_gl(_sidx, _old_link_items):
                    self.parent_window.update_joint()
            elif (self.selected_torus):
                # Throttle: same backlog-prevention as translation drag.
                # Crucially, we return BEFORE calling get_axis_angle_delta so
                # that last_drag_pos is NOT updated on skipped frames — the next
                # processed frame then computes da = full accumulated mouse motion.
                _now = _time_module.perf_counter()
                if _now - getattr(self.parent_window, '_last_drag_frame_t', 0.0) < 0.200:
                    return
                self.parent_window._last_drag_frame_t = _now

                # Start cProfile on first processed drag frame if armed.
                if getattr(self.parent_window, '_drag_profile_armed', False):
                    if getattr(self.parent_window, '_drag_profiler', None) is None:
                        self.parent_window._drag_profiler = _cProfile.Profile()
                    self.parent_window._drag_profiler.enable()

                da, normal = self.get_axis_angle_delta(event)
                if (not self.facing_same_dir):
                    da = -da

                if (self.mesh_selected):
                    pass
                else:
                    _prof = getattr(self.parent_window, '_profiling_enabled', False)
                    _drag_t0 = _time_module.perf_counter() if _prof else None

                    # Save backup before first drag transform (for reversion on error)
                    if not hasattr(self, '_drag_backup') or self._drag_backup is None:
                        self._drag_backup = self.parent_window.tree.dataDeepCopy()
                    if _prof: print(f"[PROF drag-rot] backup1: {(_time_module.perf_counter()-_drag_t0)*1000:.1f}ms")

                    # Capture old link GL items before transformJoint replaces the object.
                    _sidx = self.parent_window.selected_joint
                    _old_link_items = list(getattr(
                        self.parent_window.tree.Links[_sidx], '_gl_items', []))

                    # Real transform, no bounding-ball recompute per frame.
                    # Incremental GL update rebuilds only the moved joint+link.
                    _axis = self.selected_axis_orig
                    _transformation_rot = SE3.AngleAxis(da, [_axis[0], _axis[1], _axis[2]], unit='deg')
                    _propogate = self.parent_window.propogate_slider_checkbox.isChecked()
                    _localOrient = self.parent_window.local_orient_slider_checkbox.isChecked()
                    try:
                        self.parent_window.tree.transformJoint(
                            _sidx, _transformation_rot,
                            propogate=_propogate, relative=True, localOrient=_localOrient,
                            safe=False, lightweight=False, recomputeBoundingBall=False)
                    except Exception as _drag_err:
                        if _prof: print(f"[PROF drag-rot] transformJoint failed: {_drag_err}")
                        if self._drag_backup is not None:
                            self.parent_window.tree.setTo(self._drag_backup)
                        self.parent_window.update_joint()
                        return
                    if _prof: print(f"[PROF drag-rot] transformJoint: {(_time_module.perf_counter()-_drag_t0)*1000:.1f}ms")

                if not self.parent_window._incremental_drag_update_gl(_sidx, _old_link_items):
                    self.parent_window.update_joint()
            else:
                self._handle_camera_drag(event)


class WindowKinegamiGUI(BaseWindowKinegamiGUI):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Kinematic Tree Interactive Editor")
        # Default window to fill available screen space
        self.setWindowState(Qt.WindowMaximized)

        self.plot_widget = ClickableGLViewWidget(parent_window=self)
        self.setCentralWidget(self.plot_widget)

        self.plot_widget.setBackgroundColor(backgroundColorDefault)

        # Initialize grid attributes first
        self.grid_color = gridColorDefault
        self.grid_spacing = 10.0
        self.grid_size = 300 # total size of grid
        
        # Set camera distance appropriate for the grid size
        # Set camera distance to view the whole grid comfortably
        self.plot_widget.opts['distance'] = self.grid_size * 1.5
        
        # Create and configure grid
        self.grid = gl.GLGridItem()
        self.grid.setColor(self.grid_color)
        self.grid.setSize(self.grid_size, self.grid_size, self.grid_size)
        self.grid.setSpacing(self.grid_spacing, self.grid_spacing, self.grid_spacing)
        self.plot_widget.addItem(self.grid)
        self.grid_on = True

        self.units = "Millimeter (mm)"
        self.default_radius = TransverseRDS3225.R

        self.current_point = 0
        self.tree = None
        self.versions = []
        self.version_index = -1
        self.total_version_counter = 0
        self.tree_created = False
        self.stl_generated = False
        self.referenceMesh = None

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
        self.environment_items = []
        self.environment_visible = True

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
        self.add_mesh_widget.change_scale.connect(self.change_mesh_scale)

        # //////////////////////////////////    MESSAGE DISPLAY    ///////////////////////////////////
        self.message_display_dock = CollapsibleDockWidget("Messages", self)
        self.message_display_dock.setAllowedAreas(Qt.BottomDockWidgetArea)

        # Create a label to display success/error messages
        self.status_label = QLabel('')
        self.status_label.setWordWrap(True)
        self.message_display_dock.setWidget(self.status_label)
        
        # Initially hidden - will show when messages appear
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
        self._lightweight_dirty = False  # True when link geometry is stale from lightweight updates
        self._collision_colors_cleared = False  # True when collision highlights have been removed for motion
        self._profiling_enabled = False  # Toggle with Ctrl+Shift+P to print update_joint timings
        self._prof_call_count = 0
        self._collision_dirty = True   # Recompute collisions on next update_joint
        self._cached_colliding_joints = set()
        self._cached_colliding_links = set()
        self._gizmo_gl_items = []      # GL items for the translate/rotate gizmo arrows
        
        self.configurations_dock = CollapsibleDockWidget("Configurations and Motion", self)
        self.configurations_dock.setWidget(self.configurations_widget)
        self.configurations_dock.setVisible(True)

        # //////////////////////////////////    ADD JOINTS    ///////////////////////////////////
        self.add_transverse_revolute = QPushButton("Add Transverse Revolute Joint")
        self.add_coaxial_revolute = QPushButton("Add Coaxial Revolute Joint")
        self.add_tip = QPushButton("Add Tip")
        self.clear_tree_button = QPushButton("Clear Tree")

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

        # Select parent prompt (shown when no parent is selected and user tries to add a joint)
        self.pending_add_joint_func = None
        self.select_parent_prompt = QWidget()
        select_parent_layout = QVBoxLayout()
        self.select_parent_label = QLabel("Select a parent joint:")
        self.select_parent_label.setWordWrap(True)
        self.select_parent_combo = QComboBox()
        self.select_parent_add_btn = QPushButton("Add")
        self.select_parent_cancel_btn = QPushButton("Cancel")
        select_parent_layout.addWidget(self.select_parent_label)
        select_parent_layout.addWidget(self.select_parent_combo)
        select_parent_layout.addWidget(self.select_parent_add_btn)
        select_parent_layout.addWidget(self.select_parent_cancel_btn)
        self.select_parent_prompt.setLayout(select_parent_layout)
        self.select_parent_prompt.setVisible(False)
        self.select_parent_add_btn.clicked.connect(self._on_select_parent_add)
        self.select_parent_cancel_btn.clicked.connect(self._on_select_parent_cancel)
        self.select_parent_combo.currentIndexChanged.connect(self._on_select_parent_combo_changed)

        add_joints_layout = QVBoxLayout()
        clear_tree_layout = QVBoxLayout()
        # Add-to-root checkbox moved here from Edit Joints
        self.joint_add_root_button = QCheckBox("Add to Root")
        self.joint_add_root_button.setChecked(False)
        self.joint_add_root_button.stateChanged.connect(self.add_to_root_func)
        add_joints_layout.addWidget(self.add_transverse_revolute)
        add_joints_layout.addWidget(self.add_coaxial_revolute)
        add_joints_layout.addLayout(add_waypoints_layout)
        add_joints_layout.addWidget(self.add_tip)
        add_joints_layout.addWidget(self.select_parent_prompt)
        # place Add-to-root checkbox at bottom
        add_joints_layout.addWidget(self.joint_add_root_button)

        clear_tree_layout.addWidget(self.clear_tree_button)

        self.add_transverse_revolute.clicked.connect(self.add_transverse_revolute_toggle)
        self.add_coaxial_revolute.clicked.connect(self.add_coaxial_revolute_toggle)
        self.add_waypoint.clicked.connect(self.add_waypoint_func)
        self.insert_waypoint_button.clicked.connect(self.insert_waypoint_func)
        self.add_tip.clicked.connect(self.add_tip_toggle)
        self.clear_tree_button.clicked.connect(self.clear_tree_func)

        self.clear_tree_dock = CollapsibleDockWidget("Clear Tree", self)
        self.clear_tree_button_widget = QWidget()
        self.clear_tree_button_widget.setLayout(clear_tree_layout)
        self.clear_tree_dock.setWidget(self.clear_tree_button_widget)

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
        #self.joint_editing_layout.addLayout(radius_layout)
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

        # ////////////////////////////////    VISUAL OPTIONS    ///////////////////////////////////
        self.options_dock = CollapsibleDockWidget("Visual Options", self)
        self.options_widget = QWidget()
        self.options_layout = QVBoxLayout()

        # Keyboard shortcuts for undo/redo
        # Ctrl+Z / Ctrl+Y on Windows/Linux, Cmd+Z / Cmd+Shift+Z on Mac
        self.undo_shortcut = QShortcut(QKeySequence.Undo, self)
        self.undo_shortcut.activated.connect(self.undo)
        self.redo_shortcut = QShortcut(QKeySequence.Redo, self)
        self.redo_shortcut.activated.connect(self.redo)
        self.redo_shortcut2 = QShortcut(QKeySequence("Ctrl+Shift+Z"), self)
        self.redo_shortcut2.activated.connect(self.redo)

        # Ctrl+S: Save Tree, Ctrl+E: Export Link Modules
        self.save_shortcut = QShortcut(QKeySequence.Save, self)
        self.save_shortcut.activated.connect(self.save_tree)
        self.export_shortcut = QShortcut(QKeySequence("Ctrl+E"), self)
        self.export_shortcut.activated.connect(self.export_link_modules)
        self.profile_shortcut = QShortcut(QKeySequence("Ctrl+Shift+P"), self)
        self.profile_shortcut.activated.connect(self.toggle_profiling)
        self.drag_profile_shortcut = QShortcut(QKeySequence("Ctrl+Shift+D"), self)
        self.drag_profile_shortcut.activated.connect(self._start_drag_profile)
        self._drag_profile_armed = False
        self._drag_profiler = None

        self.edit_grid_button = QPushButton("Edit Grid")
        self.edit_grid_button.clicked.connect(self.edit_grid_func)
        self.options_layout.addWidget(self.edit_grid_button)

        self.toggle_grid = QCheckBox("Grid Visibility")
        self.toggle_grid.setChecked(True)
        self.toggle_grid.toggled.connect(self.toggle_grid_func)
        self.options_layout.addWidget(self.toggle_grid)

        # Collision highlighting toggle
        self.show_collisions = True
        self.collision_toggle = QCheckBox("Collision Highlighting")
        self.collision_toggle.setChecked(True)
        self.collision_toggle.toggled.connect(self._toggle_collision_highlighting)
        self.options_layout.addWidget(self.collision_toggle)

        self.units_label = QLabel(f"Units: {self.units}")
        self.units_label.setWordWrap(True)
        self.options_layout.addWidget(self.units_label)

        self.options_layout.addWidget(self.environment_widget)
        self.options_layout.addWidget(self.add_mesh_widget)

        self.options_widget.setLayout(self.options_layout)
        self.options_dock.setWidget(self.options_widget)

        # ////////////////////////////////    CAMERA CONTROLS DOCK    ///////////////////////////////////
        self.camera_controls_dock = CollapsibleDockWidget("Camera Controls", self)
        self.camera_options_widget = QWidget()
        self.camera_layout = QVBoxLayout()

        self.camera1 = QRadioButton("Rotate Camera")
        self.camera2 = QRadioButton("Pan Camera\n(Shift+Drag or\nMiddle Mouse)")
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

        self.save_tree_button = QPushButton('Save')
        self.save_tree_button.clicked.connect(self.save_tree)
        file_dock_layout.addWidget(self.save_tree_button)

        self.load_tree_button = QPushButton('Load')
        self.load_tree_button.clicked.connect(self.load_tree)
        file_dock_layout.addWidget(self.load_tree_button)

        self.export_link_modules_button = QPushButton('Export Link Modules')
        self.export_link_modules_button.clicked.connect(self.export_link_modules)  
        file_dock_layout.addWidget(self.export_link_modules_button) 

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
        
        self.clear_tree_widget = ClearTreeWidget(self)
        self.clear_tree_popup_dock = CollapsibleDockWidget("Clear Tree", self)
        self.clear_tree_popup_dock.setWidget(self.clear_tree_widget)
        self.clear_tree_popup_dock.setVisible(False)  # Initially hidden

        self.edit_grid_widget = EditGridWidget(self)
        self.edit_grid_dock = CollapsibleDockWidget("Edit Grid", self)
        self.edit_grid_dock.setWidget(self.edit_grid_widget)
        self.edit_grid_dock.setVisible(False)
        
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
        # self.addDockWidget(Qt.LeftDockWidgetArea, self.edit_dims_dock)  # Removed - units locked to mm
        self.addDockWidget(Qt.LeftDockWidgetArea, self.edit_grid_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.options_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.clear_tree_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.clear_tree_popup_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.camera_controls_dock)
        
        self.addDockWidget(Qt.RightDockWidgetArea, self.clear_tree_dock)
        self.addDockWidget(Qt.RightDockWidgetArea, self.clear_tree_popup_dock)
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

    LOG_CAPACITY = 100
    AUTOSAVE_FREQUENCY = 1

    @property
    def model(self):
        return self.tree

    @model.setter
    def model(self, v):
        self.tree = v

    @property
    def model_created(self) -> bool:
        return self.tree is not None

    def _mark_collision_dirty(self) -> None:
        self._collision_dirty = True
        self._collision_colors_cleared = False

    def _on_joint_selection_changed_extra(self, index: int) -> None:
        if self.select_parent_prompt.isVisible() and index >= 0:
            self.select_parent_combo.blockSignals(True)
            self.select_parent_combo.setCurrentIndex(index)
            self.select_parent_combo.blockSignals(False)

    def _on_fast_update_gl_extra(self) -> None:
        if not self._collision_colors_cleared:
            if getattr(self, 'show_collisions', True):
                from style import linkColorDefault
                for link in self.tree.Links:
                    for item in getattr(link, '_gl_items', []):
                        if hasattr(item, 'setColor'):
                            item.setColor(linkColorDefault)
            self._collision_colors_cleared = True

    def _save_model(self, autosave_id=None) -> None:
        self.save_tree(autosave_id=autosave_id)

    def _should_skip_joint_in_config(self, joint) -> bool:
        return isinstance(joint, (Waypoint, Tip))

    @QtCore.pyqtSlot(str)
    def change_units(self, key):
        self.units = key
        self.units_label.setText(f"Units: {self.units}")
        self.tree.units = key
        # self.log_version()
        self.update_joint()

    # def onUpdateRadius(self, value):
    #     # Not used in printed trees
    #     pass
        value = value / 10.0
        self.radius = value
        self.tree.changeRadius(value)
        # self.update_joint() # now called in joint_selection_changed
        self.joint_selection_changed(self.selected_joint, force=True)

    def _toggle_collision_highlighting(self, checked):
        self.show_collisions = checked
        self.update_joint()

    def delete_selected_joint(self):
        self.selected_joint = self.select_joint_options.currentIndex()
        deleted_index = self.selected_joint
        
        if self.selected_joint == -1:
            self.show_error('Please select a joint to delete.')
            return
        
        if len(self.tree.Joints) == 1:
            self.tree = None
            self.plot_widget.clear()
            self.setCentralWidget(self.plot_widget)
            self.show_success('Joint successfully deleted!')
            self.last_joint = -1
            # Clear saved configurations when tree is deleted
            self.saved_configurations.clear()
        else:
            backup_tree = copy.deepcopy(self.tree)
            try:
                # Update saved configs before deleting
                self.update_saved_configs_for_joint_deleted(deleted_index)

                # Delete joint with recursive=False (re-parents children)
                self.tree.deleteJoint(self.selected_joint, recursive=False)

                # Update the selected joint
                if self.selected_joint >= len(self.tree.Joints):
                    self.selected_joint = len(self.tree.Joints) - 1

                self._mark_collision_dirty()
                self.update_joint()
                self.show_success('Joint successfully deleted!')
            except Exception as e:
                print(f"Error deleting joint: {e}")
                import traceback
                traceback.print_exc()
                self.tree = backup_tree
                self.show_error('Error deleting joint.')
        self.window().delete_joint_dock.setVisible(False)

    def clear_tree_func(self):
        # Show the ClearTreeWidget dock when this function is called
        self.clear_tree_popup_dock.setVisible(True)
        self.clear_tree_dock.setVisible(False)

    def clear_tree(self):
        # Clear the tree
        self.tree = None
        self.tree_created = True
        self.selected_joint = -1
        self.plot_widget.clear()
        self.setCentralWidget(self.plot_widget)

        self.grid = gl.GLGridItem()
        self.grid.setColor(self.grid_color)
        self.grid.setSize(self.grid_size, self.grid_size, self.grid_size)
        self.grid.setSpacing(self.grid_spacing, self.grid_spacing, self.grid_spacing)

        if self.grid_on:
            self.plot_widget.addItem(self.grid)

        if self.referenceMesh is not None:
            self.plot_widget.addItem(self.referenceMesh.mesh)

        self.update_joint()
        self.log_version()
        self.show_success('Tree cleared!')

    def export_link_modules(self):
        if self.tree:
            options = QFileDialog.Options()
            try:
                base_path = sys._MEIPASS 
            except AttributeError:
                base_path = os.path.abspath(".")
            folder_path = QFileDialog.getExistingDirectory(
                self, "Select Folder for Link Modules", os.path.join(base_path, "save"), options=options
            )
            if folder_path:
                base_filename = os.path.join(folder_path, "module")
                self.tree.saveLinkModules(base_filename)
            self.tree.showLinkModules()

    def _gather_session_state(self, autosave=False):
        """Collect all serializable editor state into a dict for .session files.
        When autosave=True, undo/redo history is excluded to keep the file small
        and avoid serializing up to 100 deep-copies of the tree on every action."""
        state = {
            'tree': self.tree,
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
            # Camera state
            'camera_distance': self.plot_widget.opts.get('distance', 450),
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
                if md is None:
                    md = self.referenceMesh.mesh.meshDataChanged()
                    md = self.referenceMesh.mesh.opts.get('meshdata', None)
                if md is not None:
                    state['reference_mesh'] = {
                        'vertexes': md.vertexes(),
                        'faces': md.faces(),
                        'pose': self.referenceMesh.Pose,
                        'r': self.referenceMesh.r,
                    }
            except Exception:
                pass  # Skip mesh if we can't extract data
        return state

    def _restore_session_state(self, state):
        """Restore editor state from a session dict."""
        self.tree = state.get('tree', None)
        self.saved_configurations = state.get('saved_configurations', [])
        self.config_durations = state.get('config_durations', [])
        self.selected_joint = state.get('selected_joint', -1)
        self.selected_frame = state.get('selected_frame', -1)
        self.grid_color = state.get('grid_color', gridColorDefault)
        self.grid_spacing = state.get('grid_spacing', 10.0)
        self.grid_size = state.get('grid_size', 300)
        self.grid_on = state.get('grid_on', True)
        self.units = state.get('units', 'Millimeter (mm)')
        self.mesh_scale = state.get('mesh_scale', 1.0)
        self.control_type = state.get('control_type', 'Translate')
        self.is_local = state.get('is_local', True)
        self.animation_loop = state.get('animation_loop', False)
        self.versions = state.get('versions', [])
        self.version_index = state.get('version_index', -1)
        self.total_version_counter = state.get('total_version_counter', 0)

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
        self.units_label.setText(f"Units: {self.units}")

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

    def save_tree(self, autosave_id=None):
        # confusing why autosave_id is sometimes False
        if autosave_id is False:
            autosave_id = None
        if self.tree:
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
                        # Run dill serialize+write on a background thread so the
                        # UI is never blocked.  Autosaves don't need a success dialog.
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
                    self.tree.save(file_path)
        
    def load_tree(self):
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
                self.tree = loadTree(file_path)
            self._mark_collision_dirty()
            self.update_joint()
            if not file_path.endswith('.session'):
                self.log_version()

    def done_transforming(self, done):
        if done:
            self.log_version()
            # Phase 1: fast visual update using cached collision result so the
            # joint snaps to its final position immediately (~50 ms).
            self.update_joint()
            # Phase 2: recompute collision detection, then redraw with highlights.
            # Deferred via QTimer so Phase 1 can render before the ~500 ms check.
            self._mark_collision_dirty()
            def _phase2():
                self.update_joint()
                # If cProfile drag capture was active, stop it now (after full redraw).
                if getattr(self, '_drag_profiler', None) is not None:
                    self._stop_drag_profile()
            QTimer.singleShot(0, _phase2)

    """
    def update_slider(self, slider_type):
        if (slider_type == "rotation"):
            slider = self.rotation_slider
            textbox = self.rotation_textbox

            slider.setMinimum(-360)
            slider.setMaximum(360)

            rotation_matrix = self.tree.Joints[self.selected_joint].Pose.R
            angle_degrees = self.rotation_angle_from_matrix(rotation_matrix, self.selected_arrow)
            slider_value = int(angle_degrees)
            slider_text = str(int(angle_degrees))

            set_slider = self.selected_arrow != -1

            if (set_slider):
                self.old_rot_val = angle_degrees

        elif (slider_type == "translation"):
            slider = self.translation_slider
            textbox = self.translation_textbox

            amount = self.tree.Joints[self.selected_joint].Pose.t[self.selected_arrow]
            slider_value = int(amount * 10) #TODO: should this adjust for radius?
            slider_text = str(slider_value)

            set_slider = self.selected_arrow != -1

            if (set_slider):
                self.old_trans_val = amount

        elif (slider_type == "state"):
            slider = self.state_slider
            textbox = self.state_input

            min = math.degrees(self.tree.Joints[self.selected_joint].stateRange()[0])
            max = math.degrees(self.tree.Joints[self.selected_joint].stateRange()[1])

            current = int(math.degrees(self.tree.Joints[self.selected_joint].state))
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

    def add_tree(self, tree):
        self.tree = tree
        self.select_joint_options.blockSignals(True)
        self.select_joint_options.clear()
    
        for i, joint in enumerate(self.tree.Joints):
            self.select_joint_options.addItem("Joint " + str(i) + " - " + joint.__class__.__name__)
    
        self.select_joint_options.blockSignals(False) 
        self.select_joint_options.setCurrentIndex(self.selected_joint)

        self.select_link_options.blockSignals(True)
        self.select_link_options.clear()

        for link in enumerate(self.tree.Links): 
            self.select_link_options.addItem("Link " + str(i) + " - " + link.__class__.__name__)
        
        self.select_link_options.blockSignals(False)
        self.select_link_options.setCurrentIndex(self.selected_link)

    import copy

    def edit_joint_dimension(self):
        self._backup_tree = copy.deepcopy(self.tree)
        self._saved_states = []

        for idx, joint in enumerate(self._backup_tree.Joints):
            original_state = joint.state
            self._saved_states.append(original_state)

            self._backup_tree.setJointState(idx, 0)

        self._selected_joint = self.selected_joint

        target_joint = self._backup_tree.Joints[self._selected_joint]
        self.prev_joint = self._backup_tree.Joints[self._selected_joint - 1] if self._selected_joint > 0 else None

        if target_joint.__class__.__name__ == "TransverseRevolute":
            self.edit_dimension_menu.updateTransverseRevolute()
        elif target_joint.__class__.__name__ == "CoaxialRevolute":
            self.edit_dimension_menu.updateCoaxialRevolute()
        elif target_joint.__class__.__name__ in ["StartTip", "EndTip", "Tip", "PrintedHemisphere"]:
            self.edit_dimension_menu.updateTip()
        else:
            self.show_error("Uneditable joint type.")
            return

        self.edit_dimension_toggle()

    def finish_joint_edit(self, new_joint):
        try:
            new_tree = None
            for idx, joint in enumerate(self._backup_tree.Joints):
                if idx == self._selected_joint:
                    min_val, max_val = new_joint.stateRange()
                    original_state = self._saved_states[idx]
                    if original_state < min_val:
                        original_state = min_val
                    elif original_state > max_val:
                        original_state = max_val
                    self._saved_states[idx] = original_state

                    new_joint.Pose = joint.Pose
                    if new_tree is None:
                        new_tree = PrintedKinematicTree(new_joint)
                    else:
                        new_tree.append(new_joint, relative=False,
                                        fixedPosition=True, fixedOrientation=True, safe=False)
                else:
                    if new_tree is None:
                        new_tree = PrintedKinematicTree(joint)
                    else:
                        new_tree.append(joint, relative=False,
                                        fixedPosition=True, fixedOrientation=True, safe=False)
            
            for idx, joint in enumerate(new_tree.Joints):
                new_tree.setJointState(idx, self._saved_states[idx])
            
            self.tree = new_tree
            self.selected_joint = self._selected_joint
            self.window().log_version()
            self.show_success("Tree updated successfully!")
        except Exception as e:
            self.tree = self._backup_tree
            self.show_error("Error rebuilding tree: " + str(e))

        self._mark_collision_dirty()
        self.update_joint(force_recreate_config_widget=True)

    def edit_joint_state(self):
        dialog = EditJointStateDialog(self) 
        if not self.tree:
            self.show_error('Please initialize a tree.')
            # error_dialog = ErrorDialog('Please initialize a tree.')
            # error_dialog.exec_()
        if self.selected_joint == -1:
            self.show_error('Please select a joint.')
            # error_dialog = ErrorDialog('Please select a joint.')
            # error_dialog.exec_()
        elif dialog.exec_() == QDialog.Accepted:
            edit = dialog.get_state()
            if edit is not None:

                if self.tree.setJointState(self.selected_joint, math.radians(edit)):
                    self._mark_collision_dirty()
                    self.update_joint()
                    self.show_success('Joint state successfully edited!')
                    # success_dialog = SuccessDialog('Joint state successfully edited!')
                    # success_dialog.exec_()
                    min = math.degrees(self.tree.Joints[self.selected_joint].stateRange()[0])
                    max = math.degrees(self.tree.Joints[self.selected_joint].stateRange()[1])
                    current = math.degrees(self.tree.Joints[self.selected_joint].state)
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
        if self.tree and self.selected_joint != -1:
            if self.selected_arrow == 0:
                transformation = SE3.Rx(angle_radians)
            elif self.selected_arrow == 1:
                transformation = SE3.Ry(angle_radians)
            else:
                transformation = SE3.Rz(angle_radians)
            propogate = self.propogate_slider_checkbox.isChecked()
            localOrient = self.local_orient_slider_checkbox.isChecked()
            _prof = getattr(self, '_profiling_enabled', False)
            if _prof:
                _t0 = _time_module.perf_counter()
            result = self.tree.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=True, localOrient=localOrient)
            if _prof:
                print(f"[PROF adjust_rotation] transformJoint: {(_time_module.perf_counter()-_t0)*1000:.1f}ms  result={result}")
            if result:
                self.old_rot_val = int(value)
                self._mark_collision_dirty()
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
        if self.tree and self.selected_joint != -1:
            propogate = self.propogate_slider_checkbox.isChecked()
            localOrient = self.local_orient_slider_checkbox.isChecked()
            transformation = SE3()
            if (self.selected_arrow == 0):
                transformation = SE3.Tx(amount)
            if (self.selected_arrow == 1):
                transformation = SE3.Ty(amount)
            if (self.selected_arrow == 2):
                transformation = SE3.Tz(amount)
            _prof = getattr(self, '_profiling_enabled', False)
            if _prof:
                _t0 = _time_module.perf_counter()
            result = self.tree.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=True, localOrient=localOrient)
            if _prof:
                print(f"[PROF adjust_translation] transformJoint: {(_time_module.perf_counter()-_t0)*1000:.1f}ms  result={result}")
            if result:
                self.old_trans_val = actualVal
                self._mark_collision_dirty()
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
            if self.tree.setJointState(self.selected_joint, actualState):
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
        if not self.tree:
            self.show_error('Please initialize a tree.')
            # error_dialog = ErrorDialog('Please initialize a tree.')
            # error_dialog.exec_()
        if self.selected_joint == -1:
            self.show_error('Please select a joint.')
            # error_dialog = ErrorDialog('Please select a joint.')
            # error_dialog.exec_()
        else:
            self.show_delete_widget()
    """
    def reload_IDs(self): 
        if self.tree is not None:
            for index, joint in enumerate(self.tree.Joints):
                joint.id = index
    """

    def _incremental_drag_update_gl(self, selected_joint: int, old_link_items=None) -> bool:
        """Incremental GL update during drag: remove/rebuild only the selected
        joint and its incoming link, then model-matrix update everything else.
        No scene clear, no collision detection, no sidebar update.
        Returns True on success, False if a full redraw is needed instead.

        old_link_items must be the GL items from the link object *before*
        transformJoint was called — transformJoint replaces the Link object
        entirely, so by call time self.tree.Links[selected_joint] is already
        a new object with an empty _gl_items list."""
        if self.tree is None:
            return False
        # All unaffected joints/links must have a valid GL cache for model-matrix update
        for i, j in enumerate(self.tree.Joints):
            if j is not None and i != selected_joint and not j.hasGLCache():
                return False
        for i, lnk in enumerate(self.tree.Links):
            if lnk is not None and i != selected_joint and not lnk.hasGLCache():
                return False
        # Remove the old link's GL items before addToWidget rebuilds them.
        # (The new link object has an empty _gl_items, so addToWidget won't find them.)
        if old_link_items:
            for item in old_link_items:
                try:
                    self.plot_widget.removeItem(item)
                except Exception:
                    pass
        self.tree.addToWidget(
            self,
            selectedJoint=self.selected_joint,
            selectedLink=self.selected_link,
            lastJoint=self.last_joint,
            showSpheres=False,
            collidingJoints=set(),
            collidingLinks=set(),
            only_rebuild_index=selected_joint,
        )

        # Remove old gizmo items and re-add at the joint's new position.
        for item in getattr(self, '_gizmo_gl_items', []):
            try:
                self.plot_widget.removeItem(item)
            except Exception:
                pass
        self._gizmo_gl_items = []
        if self.selected_joint != -1:
            joint = self.tree.Joints[self.selected_joint]
            frame_pose = self.tree.Joints[self.selected_frame].Pose if self.selected_frame >= 0 else None
            _items_before = list(self.plot_widget.items)
            if self.control_type == "Translate":
                joint.addTranslateArrows(self, selectedArrow=self.selected_arrow,
                                         local=self.is_local, frame=frame_pose)
            else:
                joint.addRotateArrows(self, selectedArrow=self.selected_arrow,
                                      local=self.is_local, frame=frame_pose)
            self._gizmo_gl_items = [it for it in self.plot_widget.items if it not in _items_before]

        self._refresh_measurement_overlay()
        self.plot_widget.update()
        return True

    def update_joint(self, force_recreate_config_widget : bool = False):
        # ── PROFILING ────────────────────────────────────────────────────────
        _prof = getattr(self, '_profiling_enabled', False)
        _t0 = _time_module.perf_counter() if _prof else None
        _timings = {} if _prof else None
        if _prof:
            import traceback as _tb
            _stack = _tb.extract_stack()
            # Show the 2 frames above update_joint (its direct caller and grandcaller)
            _caller_frames = [f"{f.name}:{f.lineno}" for f in _stack[-4:-1]]
            _caller_str = " <- ".join(reversed(_caller_frames))
        def _mark(label):
            if _prof:
                _timings[label] = _time_module.perf_counter() - _t0
        # ─────────────────────────────────────────────────────────────────────

        # If link geometry is stale from lightweight updates, rebuild first
        if self._lightweight_dirty and self.tree is not None:
            try:
                self.tree.resyncFromLightweight()
            except Exception as _resync_err:
                print(f"[WARNING] Drag resync failed ({_resync_err}), reverting to pre-drag state.")
                _drag_backup = getattr(self.plot_widget, '_drag_backup', None)
                if _drag_backup is not None:
                    self.tree.setTo(_drag_backup)
                    self.plot_widget._drag_backup = None
            self._lightweight_dirty = False
        _mark('resync')

        self.select_joint_options.blockSignals(True)
        self.select_link_options.blockSignals(True)

        self.units_label.setText(f"Units: {self.units}")

        # Check if we need to recreate config widgets or just update values.
        # Compare actual real-joint indices, not just count, so that operations
        # that shift indices (add-to-root, waypoint insert, undo/redo, load)
        # are detected even when the number of real joints stays the same.
        need_recreate_config_widget = force_recreate_config_widget
        if not need_recreate_config_widget:
            if self.tree is None:
                need_recreate_config_widget = len(getattr(self, 'config_joint_indices', [])) > 0
            elif not hasattr(self, 'config_joint_indices'):
                need_recreate_config_widget = True
            else:
                current_real_indices = [i for i, j in enumerate(self.tree.Joints)
                                        if not isinstance(j, (Waypoint, Tip))]
                if current_real_indices != self.config_joint_indices:
                    need_recreate_config_widget = True

        if need_recreate_config_widget:
            self.update_configurations()
        else:
            self.update_config_values()
        _mark('config_widgets')

        if (not self.stl_generated):
            self.plot_widget.clear()
            self._measurement_line_item = None
            _mark('gl_clear')
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

            if self.tree is not None:
                # Compute colliding pairs for highlighting — use cache when geometry hasn't
                # changed, and skip recompute during drag (fast visual feedback during drag,
                # fresh check once after drag ends via done_transforming).
                    collidingJoints = set()
                    collidingLinks = set()
                    if getattr(self, 'show_collisions', True):
                        try:
                            if self._collision_dirty:
                                collidingPairs = self.tree.getCollidingPairs()
                                for (idx1, type1), (idx2, type2) in collidingPairs:
                                    if type1 == 'Joint':
                                        collidingJoints.add(idx1)
                                    else:
                                        collidingLinks.add(idx1)
                                    if type2 == 'Joint':
                                        collidingJoints.add(idx2)
                                    else:
                                        collidingLinks.add(idx2)
                                self._cached_colliding_joints = collidingJoints
                                self._cached_colliding_links = collidingLinks
                                self._collision_dirty = False
                            else:
                                collidingJoints = self._cached_colliding_joints
                                collidingLinks = self._cached_colliding_links
                        except Exception as e:
                            print(f"Collision detection error: {e}")
                    _mark('collision_detection')
                    self.tree.addToWidget(
                        self,
                        selectedJoint=self.selected_joint,
                        selectedLink=self.selected_link,
                        lastJoint=self.last_joint,
                        showSpheres=False,  # Don't show bounding spheres for printed joints
                        collidingJoints=collidingJoints,
                        collidingLinks=collidingLinks
                    )
                    _mark('addToWidget')
                    self.add_tree(self.tree)
                    _mark('add_tree')

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

        if self.tree is not None and self.selected_joint != -1:
            joint = self.tree.Joints[self.selected_joint]
            frame_pose = None
            if self.selected_frame >= 0:
                frame_pose = self.tree.Joints[self.selected_frame].Pose

            _items_before_gizmo = list(self.plot_widget.items)
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
            self._gizmo_gl_items = [it for it in self.plot_widget.items
                                     if it not in _items_before_gizmo]
        else:
            self._gizmo_gl_items = []
        _mark('arrows')

        if self.selected_arrow != -1:
            self.rotation_slider.setDisabled(False)
            self.translation_slider.setDisabled(False)
        else:
            self.rotation_slider.setDisabled(True)
            self.translation_slider.setDisabled(True)

        # ── PROFILING: print timings ─────────────────────────────────────────
        if _prof:
            total = _time_module.perf_counter() - _t0
            prev = 0.0
            segments = []
            for label, t in _timings.items():
                segments.append(f"  {label}: {(t-prev)*1000:.1f}ms")
                prev = t
            self._prof_call_count += 1
            print(f"[PROF #{self._prof_call_count}] update_joint total={total*1000:.1f}ms  caller: {_caller_str}")
            print("\n".join(segments))
        # ─────────────────────────────────────────────────────────────────────

        #print("current radius: " + str(self.tree.r))
                
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
    

    def toggle_profiling(self):
        self._profiling_enabled = not self._profiling_enabled
        self._prof_call_count = 0
        state = "ON" if self._profiling_enabled else "OFF"
        print(f"[PROF] update_joint profiling {state}")

    # ── cProfile drag capture ─────────────────────────────────────────────────
    # Press Ctrl+Shift+D to arm.  Then do exactly ONE drag (press, move, release).
    # A snakeviz window opens automatically with the full call-tree for that drag.

    def _start_drag_profile(self):
        """Arm cProfile for the next drag cycle."""
        self._drag_profiler = _cProfile.Profile()
        self._drag_profiler.enable()
        self._drag_profile_armed = True
        print("[CPROF] Drag profiling armed — do one drag now.")

    def _stop_drag_profile(self):
        """Stop cProfile, save stats, and launch snakeviz."""
        if not getattr(self, '_drag_profiler', None):
            return
        self._drag_profiler.disable()
        self._drag_profile_armed = False
        _stats_path = os.path.join(os.path.dirname(__file__), 'drag_profile.prof')
        self._drag_profiler.dump_stats(_stats_path)
        self._drag_profiler = None
        print(f"[CPROF] Stats saved to {_stats_path}")
        print(f"[CPROF] Opening snakeviz...")
        import subprocess
        subprocess.Popen([sys.executable, '-m', 'snakeviz', _stats_path])
    # ─────────────────────────────────────────────────────────────────────────

    def _mark_collision_dirty(self):
        """Call whenever joint/link geometry changes so the next update_joint recomputes collisions."""
        self._collision_dirty = True

    def add_joint(self, joint : Joint):
        # Check if this is a real joint (not a waypoint)
        is_real_joint = type(joint).__name__ not in ['Waypoint', 'PrintedWaypoint']
        
        if (self.tree == None or len(self.tree.Joints) == 0):
            # No tree exists, create one with this joint as root
            self.tree = PrintedKinematicTree(joint)
            self.selected_joint = 0
        elif self.add_to_root:
            # Add as new root: create a new tree with this joint, then add old tree as subtree
            self.add_joint_as_new_root(joint)
        else:
            # Tree exists, need a selected joint to add to
            if self.selected_joint == -1:
                self.show_error('Please select a parent joint first.')
                return
            
            # Add joint as child of selected joint
            new_joint_index = self.tree.addJoint(
                parentIndex=self.selected_joint,
                newJoint=joint,
                relativeToDistalDubins=True,
                fixedPosition=True,
                fixedOrientation=True,
                safe=False
            )
            self.selected_joint = new_joint_index

        # Update saved configurations if a real joint was added
        if is_real_joint:
            # Determine the config-list position of the newly added joint.
            # config_joint_indices hasn't been refreshed yet, but we can
            # compute where the new joint will land among the real joints.
            new_joint_index = self.selected_joint  # set above for every path
            real_indices = [i for i, j in enumerate(self.tree.Joints)
                           if not isinstance(j, (Waypoint, Tip))]
            config_pos = real_indices.index(new_joint_index) if new_joint_index in real_indices else len(real_indices) - 1
            self.update_saved_configs_for_joint_added(config_pos)

        self._mark_collision_dirty()
        self.update_joint()
        self.log_version()
        self.joint_selection_changed(self.selected_joint, force=True)
    
    def add_joint_as_new_root(self, joint : Joint):
        """Add a joint as the new root, making the old tree a subtree of the new root.
        
        The new joint's pose should already be in global coordinates.
        """
        old_tree = self.tree
        new_tree = PrintedKinematicTree(joint)
        new_tree.addSubtree(0, old_tree)
        self.tree = new_tree
        self.selected_joint = 0

    # def add_joint_func(self, joint_type):
    #     numSides = self.num_sides

    #     if (not self.tree_created):
    #         self.tree_not_created()
    #         return

    #     if (self.tree and len(self.tree.Joints) > 0):
    #         if not self.add_to_root:
    #             prevJoint = self.tree.Joints[-1]
    #             add_to_root = False
    #         else:
    #             prevJoint = self.tree.Joints[0]
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

    def _can_add_joint_directly(self):
        """Check if conditions are met to add a joint without prompting for parent selection."""
        # No tree or empty tree -> will create new tree
        if self.tree is None or len(self.tree.Joints) == 0:
            return True
        # Add to root mode
        if self.add_to_root:
            return True
        # Check dropdown first
        dropdown_idx = self.select_joint_options.currentIndex()
        if dropdown_idx >= 0 and self.tree and dropdown_idx < len(self.tree.Joints):
            self.selected_joint = dropdown_idx
            return True
        # Check selected_joint attribute
        if self.selected_joint >= 0 and self.selected_joint < len(self.tree.Joints):
            return True
        return False

    def _show_select_parent_prompt(self, add_func, joint_type_name):
        """Show the select parent prompt when no parent is selected."""
        self.pending_add_joint_func = add_func
        self.select_parent_combo.clear()
        if self.tree and len(self.tree.Joints) > 0:
            for i, joint in enumerate(self.tree.Joints):
                self.select_parent_combo.addItem(f"Joint {i} - {joint.__class__.__name__}")
        self.select_parent_label.setText(f"Select a parent joint for {joint_type_name}:")
        self.select_parent_add_btn.setText(f"Add {joint_type_name}")
        self.select_parent_prompt.setVisible(True)

    def _on_select_parent_add(self):
        """Called when user clicks Add in the select parent prompt."""
        idx = self.select_parent_combo.currentIndex()
        if idx >= 0 and self.tree and idx < len(self.tree.Joints):
            self.selected_joint = idx
            self.select_joint_options.blockSignals(True)
            self.select_joint_options.setCurrentIndex(idx)
            self.select_joint_options.blockSignals(False)
            self.select_parent_prompt.setVisible(False)
            if self.pending_add_joint_func:
                func = self.pending_add_joint_func
                self.pending_add_joint_func = None
                func()
        else:
            self.show_error("Please select a valid parent joint.")

    def _on_select_parent_cancel(self):
        """Cancel the pending add joint operation."""
        self.select_parent_prompt.setVisible(False)
        self.pending_add_joint_func = None

    def _on_select_parent_combo_changed(self, index):
        """When the user picks a joint from the select-parent dropdown, visually highlight it."""
        if index >= 0 and self.tree and index < len(self.tree.Joints):
            self.joint_selection_changed(index)

    def _create_and_add_transverse_revolute(self):
        """Create and add a transverse revolute joint."""
        try:
            if self.tree is None or len(self.tree.Joints) == 0:
                pose = SE3.Ry(-math.pi/2)
            elif self.add_to_root:
                root = self.tree.Joints[0]
                root_proximal_dubins = root.ProximalDubinsFrame()
                r = root.r
                distance = 4 * r + TransverseRDS3225.NEUTRAL_LENGTH / 2
                pose = root_proximal_dubins @ SE3.Trans(-distance, 0, 0)
            else:
                distance = 4 * TransverseRDS3225.R + TransverseRDS3225.NEUTRAL_LENGTH / 2
                pose = SE3.Rt(SE3().R, np.array([distance, 0, 0]))
            joint = TransverseRDS3225(pose, version=270)
            self.add_joint(joint)
        except Exception as e:
            self.show_error(str(e))

    def _create_and_add_coaxial_revolute(self):
        """Create and add a coaxial revolute joint."""
        try:
            if self.tree is None or len(self.tree.Joints) == 0:
                pose = SE3()
            elif self.add_to_root:
                root = self.tree.Joints[0]
                root_proximal_dubins = root.ProximalDubinsFrame()
                r = root.r
                distance = 4 * r + CoaxialRDS3225.NEUTRAL_LENGTH / 2
                pose = root_proximal_dubins @ SE3.Rt(SE3.Ry(math.pi/2).R, np.array([-distance, 0, 0]))
            else:
                distance = 4 * CoaxialRDS3225.R + CoaxialRDS3225.NEUTRAL_LENGTH / 2
                pose = SE3.Rt(SE3.Ry(np.pi/2).R, np.array([distance, 0, 0]))
            joint = CoaxialRDS3225(pose, version=270)
            self.add_joint(joint)
        except Exception as e:
            self.show_error(str(e))

    def _create_and_add_tip(self):
        """Create and add a tip (PrintedHemisphere)."""
        try:
            if self.tree is None or len(self.tree.Joints) == 0:
                pose = SE3()
            elif self.add_to_root:
                root = self.tree.Joints[0]
                root_proximal_dubins = root.ProximalDubinsFrame()
                r = root.r
                distance = 4 * r + r / 2  # Tip length is r
                pose = root_proximal_dubins @ SE3.Rt(SE3.Ry(math.pi/2).R, np.array([-distance, 0, 0]))
            else:
                prevJoint = self.tree.Joints[self.selected_joint]
                distance = prevJoint.r * 4 + self.default_radius / 2
                pose = SE3.Rt(SE3.Ry(np.pi/2).R, np.array([distance, 0, 0]))
            if self.tree is None or len(self.tree.Joints) == 0 or self.add_to_root:
                joint = PrintedStartHemisphere(r=TransverseRDS3225.R, Pose=pose)
            else:
                joint = PrintedEndHemisphere(r=TransverseRDS3225.R, Pose=pose)
            self.add_joint(joint)
        except Exception as e:
            self.show_error(str(e))

    def add_transverse_revolute_toggle(self):
        if self._can_add_joint_directly():
            self._create_and_add_transverse_revolute()
        else:
            self._show_select_parent_prompt(
                self._create_and_add_transverse_revolute,
                "Transverse Revolute"
            )

    def add_coaxial_revolute_toggle(self):
        if self._can_add_joint_directly():
            self._create_and_add_coaxial_revolute()
        else:
            self._show_select_parent_prompt(
                self._create_and_add_coaxial_revolute,
                "Coaxial Revolute"
            )

    def add_tip_toggle(self):
        if self._can_add_joint_directly():
            self._create_and_add_tip()
        else:
            self._show_select_parent_prompt(
                self._create_and_add_tip,
                "Tip"
            )

    def edit_dimension_toggle(self):
        self.edit_dimension_menu.setVisible(not self.edit_dimension_menu.isVisible())
        self.edit_dimension_button.setVisible(not self.edit_dimension_button.isVisible())

    def _create_and_add_waypoint(self, distance: Optional[float] = None):
        """Create and add a waypoint."""
        try:
            if (self.tree is None) or len(self.tree.Joints) == 0:
                waypoint = Waypoint(self.default_radius, SE3())
                self.tree = PrintedKinematicTree(waypoint)
                self.selected_joint = 0
                self.update_joint()
                self.log_version()
                self.joint_selection_changed(self.selected_joint, force=True)
                return
            elif self.add_to_root:
                # Compute pose in global coordinates behind the old root
                root = self.tree.Joints[0]
                root_proximal_dubins = root.ProximalDubinsFrame()
                r = root.r
                d = distance if distance is not None else 4 * r
                # Waypoint has pathIndex=2, so rotate by Ry(pi/2) so z-hat aligns with dubins x-hat
                pose = root_proximal_dubins @ SE3.Rt(SE3.Ry(np.pi/2).R, np.array([-d, 0, 0]))
                waypoint = Waypoint(r, pose)
                self.add_joint_as_new_root(waypoint)
            else:
                prevJoint = self.tree.Joints[self.selected_joint]
                # Calculate pose relative to distal Dubins frame of previous joint
                d = distance if distance is not None else 4 * prevJoint.r
                pose = SE3.Rt(SE3.Ry(np.pi/2).R, np.array([d, 0, 0]))
                waypoint = Waypoint(prevJoint.r, pose)
                self.tree.addJoint(parentIndex=self.selected_joint,
                    newJoint=waypoint, relativeToDistalDubins=True,
                    fixedPosition=True, fixedOrientation=True, safe=False)
                self.selected_joint = len(self.tree.Joints) - 1

            self.update_joint()
            self.log_version()
            if self.add_to_root:
                self.select_joint_options.setCurrentIndex(0)
            else:
                self.select_joint_options.setCurrentIndex(len(self.tree.Joints) - 1)
        except Exception as e:
            self.show_error(str(e))

    def add_waypoint_func(self):
        if (self.selected_link != -1):
            link = self.tree.Links[self.selected_link]
            nextJoint = self.tree.Links[self.selected_link]
            
            pos = link.cylinder.start + 0.5 * link.cylinder.length * link.cylinder.direction
            newPos = SE3().Trans(x=pos[0], y=pos[1], z=pos[2])

            rot = link.cylinder.orientation()
            newRot = SE3(rot)

            waypoint = Waypoint(link.r, newPos * newRot)
            waypoint_index = len(self.tree.Joints)
            parentIndex = self.tree.Parents[self.selected_link]
            
            self.tree.addJoint(parentIndex=parentIndex, newJoint=waypoint, relative=False, 
                               fixedPosition=True, fixedOrientation=False, safe=False)
            
            for child_index in self.tree.Children[self.selected_link]:
                self.tree.Parents[child_index] = waypoint_index
                self.tree.Children[waypoint_index].append(child_index)
                self.tree.Links[child_index] = PrintedLinkCSC(self.tree.r, waypoint.DistalDubinsFrame(), 
                                                self.tree.Joints[child_index].ProximalDubinsFrame())
            self.update_joint()
            self.log_version()
        elif self._can_add_joint_directly():
            self._create_and_add_waypoint()
        else:
            self._show_select_parent_prompt(
                self._create_and_add_waypoint,
                "Waypoint"
            )

    def insert_waypoint_func(self):
        if self.tree is None or len(self.tree.Joints) == 0:
            self.show_error("No tree to insert into.")
            return
        if self.selected_joint == -1:
            self.show_error("Please select a joint first.")
            return

        num = self.insert_waypoint_count.value()

        try:
            if self.selected_joint == 0:
                self._insert_waypoints_before_root(num)
            else:
                self.tree.insertWaypointsIntoLink(self.selected_joint, num)
                self._mark_collision_dirty()
                self.update_joint()
                self.log_version()
                self.joint_selection_changed(self.selected_joint, force=True)
        except Exception as e:
            self.show_error(str(e))

    def _insert_waypoints_before_root(self, numWaypoints: int):
        root = self.tree.Joints[0]
        r = root.r
        original_root_frame = root.ProximalDubinsFrame()
        eps = r * 0.01
        for i in range(numWaypoints, 0, -1):
            d = i * eps
            pose = original_root_frame @ SE3.Rt(SE3.Ry(np.pi/2).R, np.array([-d, 0, 0]))
            waypoint = Waypoint(r, pose)
            self.add_joint_as_new_root(waypoint)
        self._mark_collision_dirty()
        self.update_joint()
        self.log_version()
        self.select_joint_options.setCurrentIndex(0)

    def is_parent_joint_selected(self):
        if self.selected_joint == -1 and self.tree is not None:
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
            if (self.tree and self.selected_joint != -1):
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