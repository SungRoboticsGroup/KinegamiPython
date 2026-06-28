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

import importlib.util
from typing import Optional

import warnings
warnings.filterwarnings("ignore")

if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
    PyQt5.QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    PyQt5.QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

class DeleteWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Confirm Delete')

        layout = QVBoxLayout()
        delete_label = QLabel('Are you sure you want to delete the joint?')
        delete_label.setWordWrap(True)
        layout.addWidget(delete_label)

        self.apply_button = QPushButton('Confirm')
        self.apply_button.clicked.connect(self.on_apply_clicked)
        layout.addWidget(self.apply_button)

        self.cancel_button = QPushButton('Cancel')
        self.cancel_button.clicked.connect(self.on_cancel_clicked)
        layout.addWidget(self.cancel_button)

        self.setLayout(layout)

    def on_apply_clicked(self):
        self.window().delete_selected_joint()

    def on_cancel_clicked(self):
        self.window().delete_joint_dock.setVisible(False)

# Widget to add a new mesh to the scene, imported from a file
class AddMeshWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Reference Mesh')
        layout = QVBoxLayout()
        self.scale = 1.0

        # Input for the file path
        self.file_input_button = QPushButton('Import Mesh', self)
        self.file_input_button.clicked.connect(self.on_import_stl)
        layout.addWidget(self.file_input_button)

        scale_layout = QHBoxLayout()
        scale_label = QLabel("Scale:")
        self.scale_slider = QSlider(Qt.Horizontal, self)
        self.scale_slider.setMinimum(0)  # Minimum value
        self.scale_slider.setMaximum(100)  # Maximum value
        self.scale_slider.setValue(50)  # Initial value
        self.scale_slider.setEnabled(False)
        self.scale_slider.valueChanged.connect(self.on_update_scale)
        scale_layout.addWidget(scale_label)
        scale_layout.addWidget(self.scale_slider)
        layout.addLayout(scale_layout)

        # Clear button to remove the mesh
        self.clear_button = QPushButton('Clear Mesh', self)
        self.clear_button.clicked.connect(self.on_clear_clicked)
        layout.addWidget(self.clear_button)

        self.visible_toggle = QCheckBox('Mesh Visibility')
        self.visible_toggle.setChecked(True)
        self.visible_toggle.toggled.connect(self.toggle_visibility)  
        layout.addWidget(self.visible_toggle)

        
        self.setLayout(layout)

    change_scale = qc.pyqtSignal(float)

    def on_import_stl(self):
        options = QFileDialog.Options()
        try:
            base_path = sys._MEIPASS 
        except AttributeError:
            base_path = os.path.abspath(".")
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Import STL", os.path.join(base_path, "referenceMeshes"), "STL Files (*.stl)", options=options
        )
        if file_path:
            mesh = stlToMeshItem(file_path, scale=1)
            mesh.setObjectName("Mesh")
            self.window().referenceMesh = ReferenceMesh(mesh=mesh)
            self.scale_slider.setEnabled(True)
            self.window().update_joint()

    def toggle_visibility(self):
        self.window().mesh_visible = self.visible_toggle.isChecked()
        print("Mesh visibility toggled:", "Visible" if self.window().mesh_visible else "Hidden")
        self.window().update_joint()

    def on_update_scale(self, value):
        if (self.window().referenceMesh is not None):
            self.scale = value / 50.0
            self.window().referenceMesh.updateScale(self.scale)
            self.change_scale.emit(self.scale)

    def on_add_clicked(self):
        try:
            # Get the file path and call a function in the main window to add the mesh
            file_path = self.file_input.text()
            scale_factor_string = None
            if scale_factor_string is None or scale_factor_string == "":
                scale_factor = 1
            else:
                scale_factor = float(scale_factor_string)
            if (scale_factor <= 0):
                raise ValueError
            mesh = stlToMeshItem(file_path, scale=scale_factor)
            mesh.setObjectName("Mesh")
            self.window().referenceMesh = ReferenceMesh(mesh=mesh)
            self.scale_slider.setEnabled(True)
            self.window().update_joint()
            #plotSTL(self.window().plot_widget, file_path, SE3(), scale=scale_factor)
            #self.window().add_mesh_dock.setVisible(False)
        except ValueError:
            self.show_error("Please enter a valid file path and scale factor.")

    def on_clear_clicked(self):
        self.window().referenceMesh = None
        #self.scale_slider.setEnabled(False)
        self.window().update_joint()

    def show_error(self, message):
        QMessageBox.warning(self, "Invalid Input", message)

class EnvironmentWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout()

        self.load_button = QPushButton('Import Environment', self)
        self.load_button.clicked.connect(self.on_load)
        layout.addWidget(self.load_button)

        self.clear_button = QPushButton('Clear Environment', self)
        self.clear_button.clicked.connect(self.on_clear)
        layout.addWidget(self.clear_button)

        self.visible_toggle = QCheckBox('Environment Visibility')
        self.visible_toggle.setChecked(True)
        self.visible_toggle.toggled.connect(self.on_toggle_visibility)
        layout.addWidget(self.visible_toggle)

        self.setLayout(layout)

    def on_load(self):
        try:
            base_path = sys._MEIPASS
        except AttributeError:
            base_path = os.path.abspath(".")
        environments_dir = os.path.join(base_path, "environments")
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Import Environment", environments_dir, "Python Files (*.py)"
        )
        if not file_path:
            return
        env_dir = os.path.dirname(os.path.abspath(file_path))
        repo_root = os.path.dirname(env_dir)
        try:
            if env_dir not in sys.path:
                sys.path.insert(0, env_dir)
            if repo_root not in sys.path:
                sys.path.insert(0, repo_root)
            spec = importlib.util.spec_from_file_location("_environment", file_path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            self.window().environment_items = module.get_items()
        except Exception as e:
            QMessageBox.warning(self, "Load Error", f"Failed to load environment:\n{e}")
            return
        self.window().update_joint()

    def on_toggle_visibility(self, checked):
        self.window().environment_visible = checked
        self.window().update_joint()

    def on_clear(self):
        self.window().environment_items = []
        self.window().update_joint()


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

class EditGridWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Edit Grid')

        layout = QVBoxLayout()

        # Input for line spacing
        spacing_header = QHBoxLayout()
        spacing_label = QLabel("Grid line spacing")
        self.unit_label = QLabel("(cm)")
        spacing_header.addWidget(spacing_label)
        spacing_header.addWidget(self.unit_label)
        spacing_header.addStretch()
        layout.addLayout(spacing_header)
        self.spacing_input = QLineEdit()
        self.spacing_input.setPlaceholderText("Enter spacing")
        layout.addWidget(self.spacing_input)

        # Input for the amount of lines
        line_amt_label = QLabel("Size (grid lines):")
        layout.addWidget(line_amt_label)
        self.line_amt_input = QLineEdit()
        self.line_amt_input.setPlaceholderText("Enter line amount")
        layout.addWidget(self.line_amt_input)

        button_layout = QHBoxLayout()

        # Apply button to create the tree
        apply_button = QPushButton('Apply Changes', self)
        apply_button.clicked.connect(self.on_apply_clicked)
        button_layout.addWidget(apply_button)

        # Cancel button to close the widget
        self.cancel_button = QPushButton('Cancel', self)
        self.cancel_button.clicked.connect(self.on_cancel_clicked)
        button_layout.addWidget(self.cancel_button)

        layout.addLayout(button_layout)

        self.setLayout(layout)

    def on_apply_clicked(self):
        try:
            # Get input values
            spacing = int(self.spacing_input.text())
            line_amt = int(self.line_amt_input.text())

            # Send signals to parent
            self.window().grid_size = line_amt * spacing
            self.window().grid_spacing = spacing
            self.window().initialize_grid()



        except ValueError:
            self.show_error("Please enter valid integers.")

    def on_cancel_clicked(self):
        # Hide the widget if the user cancels
        self.window().edit_grid_dock.setVisible(False)

    def show_error(self, message):
        QMessageBox.warning(self, "Invalid Input", message)
    
class CollapsibleDockWidget(QDockWidget):
    _TB_HEIGHT = 24

    def __init__(self, title, parent=None):
        super().__init__(title, parent)
        self._is_collapsed = False
        self.setFeatures(QDockWidget.DockWidgetMovable | QDockWidget.DockWidgetFloatable)
        self._setup_title_bar(title)

    def _setup_title_bar(self, title):
        bar = QWidget()
        bar.setFixedHeight(self._TB_HEIGHT)
        bar.setObjectName("dockTitleBar")
        bar.setStyleSheet("""
            #dockTitleBar { background-color: #4a4a4a; }
            #dockTitleBar QLabel { color: white; background: transparent; }
            #dockTitleBar QPushButton { color: white; background: transparent; border: none; }
            #dockTitleBar QPushButton:hover { background-color: #666666; border-radius: 3px; }
        """)
        layout = QHBoxLayout(bar)
        layout.setContentsMargins(4, 1, 4, 1)
        layout.setSpacing(4)

        self._collapse_btn = QPushButton("▼")
        self._collapse_btn.setFixedSize(self._TB_HEIGHT - 2, self._TB_HEIGHT - 2)
        self._collapse_btn.setFlat(True)
        self._collapse_btn.clicked.connect(self._toggle_collapse)
        layout.addWidget(self._collapse_btn)

        layout.addWidget(QLabel(title))
        layout.addStretch()

        float_btn = QPushButton("⧉")
        float_btn.setFixedSize(self._TB_HEIGHT - 2, self._TB_HEIGHT - 2)
        float_btn.setFlat(True)
        float_btn.setToolTip("Float / Dock")
        float_btn.clicked.connect(lambda: self.setFloating(not self.isFloating()))
        layout.addWidget(float_btn)

        self.setTitleBarWidget(bar)

    def _toggle_collapse(self):
        self._is_collapsed = not self._is_collapsed
        w = self.widget()
        if w is not None:
            w.setVisible(not self._is_collapsed)
        self._collapse_btn.setText("▶" if self._is_collapsed else "▼")
        if not self.isFloating():
            self.setMaximumHeight(self._TB_HEIGHT + 2 if self._is_collapsed else 16777215)


class ImageRadioButton(QRadioButton):
    def __init__(self, unchecked_img, checked_img, tooltip_text, parent=None):
        super().__init__(parent)
        self.unchecked_img = QPixmap(unchecked_img)
        self.checked_img = QPixmap(checked_img)
        
        # Only set icon size if pixmaps loaded successfully
        if not self.unchecked_img.isNull():
            self.setIconSize(self.unchecked_img.size())
        
        self.update_icon()

        # Hide the default radio button indicator
        self.setStyleSheet("QRadioButton::indicator { width: 0px; height: 0px; }")

        # Connect the toggled signal to update the icon when the state changes
        self.toggled.connect(self.update_icon)

        self.setToolTip(tooltip_text)

    def update_icon(self):
        # Only set icons if pixmaps are valid
        if self.isChecked() and not self.checked_img.isNull():
            self.setIcon(QIcon(self.checked_img))
        elif not self.isChecked() and not self.unchecked_img.isNull():
            self.setIcon(QIcon(self.unchecked_img))   

class OverlayLine(gl.GLLinePlotItem):
    def paint(self):
        glDisable(GL_DEPTH_TEST)
        super().paint()
        glEnable(GL_DEPTH_TEST)

class ClickableGLViewWidget(gl.GLViewWidget):
    def __init__(self, parent_window, parent=None):
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
        self.last_drag_pos = None
        self.parent_window = parent_window
        self.orbit_speed = 0.3

        dist = self.opts['distance']
        self.near_clip = dist * 0.001
        self.far_clip = dist * 1000.
        
        self.hit_cylinder = False
        self.selected_axis = None
        self.selected_axis_name = None
        self.selected_torus = None
        self.drag_prev_vector = 0
        self.facing_same_dir = False

        self.selected_joint_axes = [QVector3D(1,0,0), QVector3D(0,1,0), QVector3D(0,0,1)]
        self.original_axes = [QVector3D(1,0,0), QVector3D(0,1,0), QVector3D(0,0,1)]
    
    lock_status_changed = pyqtSignal(bool)
    click_signal = qc.pyqtSignal(int)
    click_signal_arrow = qc.pyqtSignal(int)
    click_signal_link = qc.pyqtSignal(int)
    click_signal_mesh = qc.pyqtSignal(bool)
    drag_change_position = qc.pyqtSignal(np.ndarray)
    drag_change_rotation = qc.pyqtSignal(float)
    key_pressed = qc.pyqtSignal(str)
    done_transforming = qc.pyqtSignal(bool)
    measure_joint_selected = qc.pyqtSignal(int)

    selected_index = -1
    selected_link_index = -1
    mesh_selected = False

    selected_joint_temp = None
    camera_type = "Rotate"

    def toggle_lock(self):
        self.locked = not self.locked
        self.lock_status_changed.emit(self.locked)
        print("Screen lock toggled:", "Locked" if self.locked else "Unlocked")
    
    def get_world_coordinates(self, event):
        pos = event.localPos()
        ndc_x = (2.0 * pos.x()) / self.width() - 1.0
        ndc_y = 1.0 - (2.0 * pos.y()) / self.height()

        view = self.viewMatrix()
        proj = self.projectionMatrix()
        inverted_matrix = (proj * view).inverted()[0]

        # Use camera position as ray origin for better precision when far away
        camera_pos = self.cameraPosition()
        
        # Get a point in world space from the click position
        ndc = QVector3D(ndc_x, ndc_y, 0.0)  # Use z=0 (middle of NDC depth range)
        world_point = inverted_matrix.map(ndc)

        # Direction from camera to the world point
        direction = world_point - camera_pos
        direction.normalize()

        return camera_pos, direction
    
    def get_normalized_plane_vectors(self, event):
        if self.is_dragging and self.selected_torus:
            origin, dir = self.get_world_coordinates(event)
            selected_joint = self.parent_window.tree.Joints[self.parent_window.selected_joint]
            center = selected_joint.Pose.t
            qcenter = QVector3D(center[0], center[1], center[2])
            qaxis = self.selected_torus
            npos = compute_plane_intersection(origin, dir, qaxis, qcenter)

            plane_vector = npos - qcenter
            plane_vector.normalize()

            normal = QVector3D.dotProduct(dir, qaxis) * qaxis
            normal.normalize()

            return plane_vector, normal
    
    def get_axis_angle_delta(self, event):
        if self.is_dragging and self.selected_torus:
            origin, dir = self.get_world_coordinates(event)
            selected_joint = self.parent_window.tree.Joints[self.parent_window.selected_joint]
            center = selected_joint.Pose.t
            qcenter = QVector3D(center[0], center[1], center[2])
            qaxis = self.selected_torus
            
            plane_vector, normal = self.get_normalized_plane_vectors(event)

            prev_vector = self.drag_prev_vector

            d = QVector3D.dotProduct(plane_vector, prev_vector)
            prod = plane_vector.length() * prev_vector.length()
            if prod == 0:
                return
            angle = math.acos(d / (plane_vector.length() * prev_vector.length()))

            cross_product = QVector3D.crossProduct(prev_vector, plane_vector)
            
            if QVector3D.dotProduct(cross_product, normal) < 0:
                angle *= -1

            self.drag_prev_vector = plane_vector

            return math.degrees(angle), normal
    
    def get_closest_point(self, event):
        if self.is_dragging and self.selected_axis:
            origin, dir = self.get_world_coordinates(event)

            if (self.mesh_selected):
                selected_joint = self.parent_window.referenceMesh
            else:
                selected_joint = self.parent_window.tree.Joints[self.parent_window.selected_joint]

            joint_center = selected_joint.Pose.t
            qcenter = QVector3D(joint_center[0], joint_center[1], joint_center[2])
            axis = self.selected_axis
            new_pos_3D = compute_closest_point_on_axis(origin, dir, qcenter, axis)

            return new_pos_3D
        
    def compute_joint_axes(self, selected_joint):
        axes = [selected_joint.Pose.R[:, i] for i in range(3)]
        axis_x = QVector3D(axes[0][0], axes[0][1], axes[0][2])
        axis_y = QVector3D(axes[1][0], axes[1][1], axes[1][2])
        axis_z = QVector3D(axes[2][0], axes[2][1], axes[2][2])

        self.selected_joint_axes[0] = axis_x
        self.selected_joint_axes[1] = axis_y
        self.selected_joint_axes[2] = axis_z
        
    def calculate_translation_isect(self, origin, direction, joint_center, selected_joint, event):
        desired_px = 80
        thickness_px = 10
        arrow_len = self.world_length_for_pixel_length(desired_px)
        arrow_thick = self.world_length_for_pixel_length(thickness_px)

        shortest_loc = float('inf')
        idx = -1
        hit = False

        for i in range(3):
            axis = (self.selected_joint_axes[i]
                    if self.is_local else self.original_axes[i])
            loc = compute_cylinder_intersection(
                origin, direction,
                joint_center, axis,
                arrow_thick,  
                arrow_len 
            )
            if loc < shortest_loc:
                shortest_loc = loc
                idx = i
                hit = True

        if hit:
            self.click_signal_arrow.emit(idx)
            self.selected_axis = (self.selected_joint_axes[idx]
                                  if self.is_local else self.original_axes[idx])
            self.is_dragging = True
        else:
            self.is_dragging = False
            self.selected_axis = None

    def calculate_rotation_isect(self, origin, direction, joint_center, selected_joint):
        desired_px = 80
        thickness_px = 10
        arrow_len = self.world_length_for_pixel_length(desired_px)
        arrow_thick = self.world_length_for_pixel_length(thickness_px)

        shortest_loc = float('inf')
        shortest_idx = -1
        hit_torus = False

        for i in range(3):
            axis = (self.selected_joint_axes[i]
                    if self.is_local else self.original_axes[i])
            loc = compute_torus_intersection(
                origin, direction,
                joint_center, axis,
                major_radius=arrow_len,   
                minor_radius=arrow_thick   
            )
            if loc < shortest_loc:
                shortest_loc = loc
                shortest_idx = i
                hit_torus = True

        if hit_torus:
            self.selected_torus = (self.selected_joint_axes[shortest_idx]
                                      if self.is_local else self.original_axes[shortest_idx])
            self.selected_axis_orig = self.original_axes[shortest_idx]

            dot = QVector3D.dotProduct(direction, self.selected_torus)
            self.facing_same_dir = (dot > 0)
            self.click_signal_arrow.emit(shortest_idx)
        else:
            self.is_dragging = False
            self.selected_torus = None

    def mousePressEvent(self, event):
        _prof = getattr(self.parent_window, '_profiling_enabled', False)
        _mpe_t0 = _time_module.perf_counter() if _prof else None

        # check to see if link or mesh is selected
        lpos = event.position() if hasattr(event, 'position') else event.localPos()
        region = [lpos.x()-5, lpos.y()-5, 10, 10]
        dpr = self.devicePixelRatioF()
        region = tuple([x * dpr for x in region])

        mesh = []

        # itemsAt does a full GL_SELECT re-render (~150ms) just to detect mesh clicks.
        # Skip it entirely when no reference mesh is loaded — there's nothing to detect.
        if not self.is_dragging and self.parent_window.referenceMesh is not None:
            import io
            _old_stdout = sys.stdout
            sys.stdout = io.StringIO()
            try:
                picked_items = list(self.itemsAt(region))
            finally:
                sys.stdout = _old_stdout

            if _prof:
                _items_at_ms = (_time_module.perf_counter() - _mpe_t0) * 1000
                print(f"[PROF] mousePressEvent itemsAt={_items_at_ms:.1f}ms  items={len(picked_items)}")

            for item in picked_items:
                if (item.objectName() == "Mesh"):
                    mesh.append(item)
        elif _prof:
            print(f"[PROF] mousePressEvent itemsAt=SKIPPED (no mesh loaded)")

        if (not self.is_dragging):
            if (len(mesh) == 0):
                self.mesh_selected = False
            else:
                self.mesh_selected = True
        
        #print("Selected:", self.mesh_selected)

        self.click_signal_link.emit(self.selected_link_index)
        self.click_signal_mesh.emit(self.mesh_selected)
        
        self.last_drag_pos = event.pos()

        if (event.buttons() and event.button() == Qt.MouseButton.MiddleButton):
            self.drag_start_pos = event.pos()

        if (event.buttons() and Qt.LeftButton and event.buttons() != QtCore.Qt.MouseButton.MiddleButton):
            self.is_dragging = False
            self.drag_start_pos = event.pos()

            if self.parent_window.measure_select_mode:
                origin, direction = self.get_world_coordinates(event)
                closest_joint_dist = float('inf')
                closest_joint = None
                if self.parent_window.tree:
                    for joint in self.parent_window.tree.Joints:
                        center = joint.Pose.t
                        radius = joint.r
                        hit_location = compute_sphere_intersection(origin, direction, center, 1.1 * radius)
                        if hit_location < closest_joint_dist:
                            closest_joint_dist = hit_location
                            closest_joint = joint
                self.selected_joint_temp = closest_joint
                self.last_drag_pos = event.pos()
                return

            origin, direction = self.get_world_coordinates(event)

            self.is_local = self.parent_window.is_local

            # raycasting for mesh's widgets
            if (self.mesh_selected and self.parent_window.selected_joint == -1):
                self.selected_axis = None
                self.selected_torus = None
                mesh = self.window().referenceMesh

                mesh_center = mesh.Pose.t
                mesh_center = QVector3D(mesh_center[0], mesh_center[1], mesh_center[2])
                self.compute_joint_axes(mesh)

                # translation widget - cylinder intersection
                if (self.parent_window.control_type == "Translate"):
                    self.calculate_translation_isect(origin, direction, mesh_center, mesh, event)

                # rotation widget - torus intersection
                elif (self.parent_window.control_type == "Rotate"):
                    self.calculate_rotation_isect(origin, direction, mesh_center, mesh)

                # print(self.selected_axis)

            # raycasting for widgets
            elif (self.parent_window.selected_joint != -1):
                self.selected_axis = None
                self.selected_torus = None
                selected_joint = self.parent_window.tree.Joints[self.parent_window.selected_joint]
                joint_center = selected_joint.Pose.t
                joint_center = QVector3D(joint_center[0], joint_center[1], joint_center[2])

                self.compute_joint_axes(selected_joint)

                # translation widget - cylinder intersection
                if (self.parent_window.control_type == "Translate"):
                    self.calculate_translation_isect(origin, direction, joint_center, selected_joint, event)

                # rotation widget - torus intersection
                elif (self.parent_window.control_type == "Rotate"):
                    self.calculate_rotation_isect(origin, direction, joint_center, selected_joint)

            if self.selected_torus:
                self.is_dragging = True
                self.drag_prev_vector, _ = self.get_normalized_plane_vectors(event)

            # check to see if a joint is clicked
            closest_joint_dist = float('inf')
            closest_joint = None

            if (self.parent_window.tree):
                for joint in self.parent_window.tree.Joints:
                    center = joint.Pose.t
                    radius = joint.r
                    hit_location = compute_sphere_intersection(origin, direction, center, 1.1*radius)
                    if (hit_location < closest_joint_dist):  # Will automatically exclude float('inf') sentinel
                        closest_joint_dist = hit_location
                        closest_joint = joint
            
            self.selected_joint_temp = closest_joint

    def mouseMoveEvent(self, event):
        if self.parent_window.measure_select_mode:
            if (event.buttons() and (Qt.LeftButton or Qt.MiddleButton)) and (event.pos() - self.drag_start_pos).manhattanLength() >= QApplication.startDragDistance():
                self.is_dragging = True
            if not self.is_dragging:
                return
            curr_pos = event.position() if hasattr(event, 'position') else event.localPos()
            diff = curr_pos - self.last_drag_pos
            self.last_drag_pos = curr_pos
            if event.buttons() == QtCore.Qt.MouseButton.MiddleButton:
                self.pan(diff.x(), diff.y(), 0, relative='view')
            elif event.buttons() == QtCore.Qt.MouseButton.LeftButton:
                if event.modifiers() & Qt.ShiftModifier:
                    self.pan(diff.x(), diff.y(), 0, relative='view')
                elif self.camera_type == "Rotate":
                    self.orbit(-diff.x() * self.orbit_speed, diff.y() * self.orbit_speed)
                elif self.camera_type == "Pan":
                    self.pan(diff.x(), diff.y(), 0, relative='view')
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
                curr_pos = event.position() if hasattr(event, 'position') else event.localPos()

                diff = curr_pos - self.last_drag_pos
                self.last_drag_pos = curr_pos

                # Check for Shift modifier (for trackpad users)
                is_shift_pressed = event.modifiers() & Qt.ShiftModifier

                if event.buttons() == QtCore.Qt.MouseButton.MiddleButton:
                    self.pan(diff.x(), diff.y(), 0, relative='view')
                elif event.buttons() == QtCore.Qt.MouseButton.LeftButton:
                    # Shift + Left Click/Drag = Pan (trackpad-friendly)
                    if is_shift_pressed:
                        self.pan(diff.x(), diff.y(), 0, relative='view')
                    elif (self.camera_type == "Rotate"):
                        self.orbit(-diff.x()*self.orbit_speed, diff.y()*self.orbit_speed)
                    elif (self.camera_type == "Pan"):
                        self.pan(diff.x(), diff.y(), 0, relative='view')

    def mouseReleaseEvent(self, event):
        if self.parent_window.measure_select_mode:
            was_dragging = self.is_dragging
            self.is_dragging = False
            if not was_dragging and self.selected_joint_temp is not None:
                joint_index = -1
                for i, j in enumerate(self.parent_window.tree.Joints):
                    if j is self.selected_joint_temp:
                        joint_index = i
                        break
                self.measure_joint_selected.emit(joint_index)
            self.selected_joint_temp = None
            return

        if self.is_dragging and (self.selected_axis or self.selected_torus):
            self._drag_backup = None  # Clear drag backup on release
            self.parent_window._last_drag_frame_t = 0.0  # Reset throttle for next drag
            self.done_transforming.emit(True)

        if (self.selected_joint_temp != None):
            # Find joint index by object identity (is) instead of equality
            joint_index = -1
            for i, joint in enumerate(self.parent_window.tree.Joints):
                if joint is self.selected_joint_temp:
                    joint_index = i
                    break
            self.click_signal.emit(joint_index)
        elif (not self.is_dragging):
            self.selected_axis = None
            self.selected_torus = None

            self.click_signal.emit(-1)
        self.is_dragging = False
        self.drag_prev_vector = None


    def keyPressEvent(self, event: QKeyEvent):
        if event.key() == Qt.Key_T:
            self.key_pressed.emit("Translate")
        elif event.key() == Qt.Key_R:
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
        elif event.key() == Qt.Key_G:
            self.key_pressed.emit("G")

    def world_length_for_pixel_length(self, pixel_length: float) -> float:
        dist = self.opts['distance']
        fov_deg = self.opts.get('fov', 60)
        fov = math.radians(fov_deg)

        h_pixels = self.height()
        if h_pixels == 0:
            return 1.0 
        angle_span = fov * (pixel_length / h_pixels)
        return dist * math.tan(angle_span)
    
    def wheelEvent(self, ev):
        super().wheelEvent(ev)
        self.parent_window.update_joint()
 
class WindowKinegamiGUI(QMainWindow):

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

    def edit_grid_func(self):
        visibility = self.edit_grid_dock.isVisible()
        self.edit_grid_dock.setVisible(not visibility)

    @QtCore.pyqtSlot(str)
    def change_units(self, key):
        self.units = key
        self.units_label.setText(f"Units: {self.units}")
        self.tree.units = key
        # self.log_version()
        self.update_joint()

    def update_configurations(self):
        """Update the configurations widget with text boxes for each real (non-Waypoint) joint"""
        # Clear existing text boxes and layouts completely
        while self.configurations_layout.count():
            item = self.configurations_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
            elif item.layout():
                # Clear the layout
                while item.layout().count():
                    child = item.layout().takeAt(0)
                    if child.widget():
                        child.widget().deleteLater()
                    elif child.layout():
                        while child.layout().count():
                            subchild = child.layout().takeAt(0)
                            if subchild.widget():
                                subchild.widget().deleteLater()
                        child.layout().deleteLater()
                item.layout().deleteLater()
        self.config_text_boxes.clear()
        self.config_sliders.clear()
        self.config_joint_indices.clear()
        
        if self.tree is None:
            return
        
        # Build individual joint column widgets (stored for grid-based layout)
        self.config_joint_widgets = []
        for actual_joint_index, joint in enumerate(self.tree.Joints):
            # Skip Waypoint and Tip joints (they have no movable state)
            if isinstance(joint, (Waypoint, Tip)):
                continue

            # Create a container widget for each joint
            joint_widget = QWidget()
            joint_widget.setFixedWidth(55)
            joint_layout = QVBoxLayout(joint_widget)
            joint_layout.setContentsMargins(0, 0, 0, 0)
            joint_layout.setSpacing(2)
            
            joint_label = QLabel(f"J{actual_joint_index}")
            joint_label.setAlignment(Qt.AlignCenter)
            
            text_box = QLineEdit()
            text_box.setPlaceholderText("0.0")
            text_box.setFixedWidth(50)
            text_box.setAlignment(Qt.AlignCenter)
            # Connect to handler with lambda to capture the actual joint index
            text_box.returnPressed.connect(lambda idx=actual_joint_index: self.config_textbox_return(idx))
            text_box.installEventFilter(self)
            
            # Create slider
            slider = QSlider(Qt.Horizontal)
            slider.setFixedWidth(50)
            
            # Set slider range and value based on joint type
            if isinstance(joint, Prismatic):
                display_state = joint.state
                state_range = joint.stateRange()
                # Scale by 100*r for prismatic joints
                scale = 100 * self.tree.r
                slider.setMinimum(int(state_range[0] * scale))
                slider.setMaximum(int(state_range[1] * scale))
                slider.setValue(int(display_state * scale))
            elif isinstance(joint, Revolute):
                display_state = math.degrees(joint.state)
                state_range = joint.stateRange()
                # Use degrees for revolute joints
                slider.setMinimum(int(math.degrees(state_range[0])))
                slider.setMaximum(int(math.degrees(state_range[1])))
                slider.setValue(int(display_state))
            else:
                display_state = 0.0
            
            # Connect slider to handler
            slider.sliderMoved.connect(lambda value, idx=actual_joint_index: self.config_slider_moved(idx, value))
            slider.sliderReleased.connect(lambda idx=actual_joint_index: self.config_slider_released(idx))
            
            text_box.setText(str(round(display_state, 2)))
            
            joint_layout.addWidget(joint_label)
            joint_layout.addWidget(text_box, 0, Qt.AlignHCenter)
            joint_layout.addWidget(slider, 0, Qt.AlignHCenter)
            joint_layout.addStretch()
            
            self.config_joint_widgets.append(joint_widget)
            self.config_text_boxes.append(text_box)
            self.config_sliders.append(slider)
            self.config_joint_indices.append(actual_joint_index)
        
        # Create "Define Configuration" button (stored for grid-based layout)
        self.config_save_button = QPushButton("Define Configuration")
        self.config_save_button.clicked.connect(self.save_current_configuration)
        self.config_save_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        
        # Create horizontal layout for saved configs and interpolation slider
        self.saved_configs_container = QHBoxLayout()
        self.configurations_layout.addLayout(self.saved_configs_container)
        
        # Display saved configurations
        self.display_saved_configurations()
        
        # Highlight the selected joint if any
        self.highlight_selected_config_box()
    
    def update_config_values(self):
        """Update the values in configuration text boxes and sliders without recreating widgets"""
        if not hasattr(self, 'config_text_boxes') or not self.config_text_boxes:
            return
        
        if self.tree is None:
            return
        
        # Update each text box and slider with current joint state
        for i, joint_index in enumerate(self.config_joint_indices):
            if joint_index >= len(self.tree.Joints):
                continue
            
            joint = self.tree.Joints[joint_index]
            text_box = self.config_text_boxes[i]
            slider = self.config_sliders[i]
            
            # Block signals to prevent triggering updates during value setting
            text_box.blockSignals(True)
            slider.blockSignals(True)
            
            if isinstance(joint, Prismatic):
                display_state = joint.state
                text_box.setText(str(round(display_state, 2)))
                slider.setValue(int(display_state * 100 * self.tree.r))
            elif isinstance(joint, Revolute):
                display_state = math.degrees(joint.state)
                text_box.setText(str(round(display_state, 2)))
                slider.setValue(int(display_state))
            
            text_box.blockSignals(False)
            slider.blockSignals(False)

    def highlight_selected_config_box(self):
        """Highlight the configuration text box for the currently selected joint"""
        if not hasattr(self, 'config_text_boxes') or not self.config_text_boxes:
            return
        
        # Reset all text boxes to default style
        for text_box in self.config_text_boxes:
            text_box.setStyleSheet("")
        
        # Highlight the selected joint's text box
        if self.selected_joint != -1 and self.selected_joint in self.config_joint_indices:
            try:
                text_box_index = self.config_joint_indices.index(self.selected_joint)
                self.config_text_boxes[text_box_index].setStyleSheet(
                    "background-color: #FFD700; border: 2px solid #FFA500;"
                )
            except (ValueError, IndexError):
                pass

    def save_current_configuration(self):
        """Save the current configuration from the text boxes"""
        if not self.config_text_boxes:
            return
        
        current_config = []
        for text_box in self.config_text_boxes:
            try:
                value = float(text_box.text())
                current_config.append(value)
            except ValueError:
                current_config.append(0.0)
        
        # Determine where to insert the configuration
        insert_index = len(self.saved_configurations)  # Default: append at end


        # if not hasattr(self, 'config_interp_slider'):
        #     print("Slider has never been created")
        # elif sip.isdeleted(self.config_interp_slider):
        #     print("Slider has been deleted")
        # elif self.config_interp_slider is None:
        #     print("Slider is None")
        # else:
        #     print(self.config_interp_slider.value())

        # Check if interpolation slider exists and is between configurations
        if hasattr(self, 'config_interp_slider') and \
            len(self.saved_configurations) > 1 and \
        not sip.isdeleted(self.config_interp_slider):
            slider_value = self.config_interp_slider.value()
            max_value = self.config_interp_slider.maximum()
            # Invert because vertical sliders have max at top
            inverted_value = max_value - slider_value
            config_value = inverted_value / 100.0
            
            # Round to nearest integer to find closest config position
            nearest_config = round(config_value)
            
            # Check if we're between two configurations (not exactly at one)
            # Use a small tolerance for floating point comparison
            if abs(config_value - nearest_config) > 0.01:
                # We're between configs - insert after the lower one
                insert_index = int(config_value) + 1
                # Make sure we don't exceed bounds
                insert_index = min(insert_index, len(self.saved_configurations))
        
        self.saved_configurations.insert(insert_index, current_config)
        self.display_saved_configurations()
        
        # Move the interpolation slider to the newly added configuration
        if hasattr(self, 'config_interp_slider') and \
        not sip.isdeleted(self.config_interp_slider):
            max_value = self.config_interp_slider.maximum()
            slider_value = max_value - (insert_index * 100)
            self.config_interp_slider.blockSignals(True)
            self.config_interp_slider.setValue(slider_value)
            self.config_interp_slider.blockSignals(False)
    
    def display_saved_configurations(self):
        """Display all saved configurations as rows below the current configuration"""
        # Save the current interpolation slider position before recreating
        saved_slider_value = None
        if hasattr(self, 'config_interp_slider') and \
            len(self.saved_configurations) > 1 and \
            not sip.isdeleted(self.config_interp_slider):
            saved_slider_value = self.config_interp_slider.value()
        
        # Clear the saved configs container
        # First, reparent joint widgets and save button so they're not destroyed
        if hasattr(self, 'config_joint_widgets'):
            for w in self.config_joint_widgets:
                w.setParent(None)
        if hasattr(self, 'config_save_button') and self.config_save_button is not None:
            self.config_save_button.setParent(None)
        while self.saved_configs_container.count():
            item = self.saved_configs_container.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
            elif item.layout():
                while item.layout().count():
                    child = item.layout().takeAt(0)
                    if child.widget():
                        child.widget().deleteLater()
                item.layout().deleteLater()
        
        num_joints = len(self.config_joint_widgets)
        
        if not self.saved_configurations:
            # Even with no configs, show the textbox row centered using a grid
            self.saved_configs_container.addStretch()
            combined_widget = QWidget()
            grid = QGridLayout(combined_widget)
            grid.setContentsMargins(0, 0, 0, 0)
            grid.setSpacing(0)
            # Row 0: empty cols 0-1, joint widgets, save button
            for i, jw in enumerate(self.config_joint_widgets):
                grid.addWidget(jw, 0, i + 2)
            grid.addWidget(self.config_save_button, 0, num_joints + 2, 1, 2)
            combined_widget.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
            self.saved_configs_container.addWidget(combined_widget)
            self.saved_configs_container.addStretch()
            return
        
        # Create a combined widget using QGridLayout for perfect column alignment
        combined_widget = QWidget()
        grid = QGridLayout(combined_widget)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setSpacing(0)
        
        # Row 0: header row with joint widgets
        # Cols 0-1 are empty (reserved for ▲▼ buttons in config rows)
        for i, jw in enumerate(self.config_joint_widgets):
            grid.addWidget(jw, 0, i + 2)
        grid.addWidget(self.config_save_button, 0, num_joints + 2, 1, 2)
        
        # Rows 1+: saved configuration rows
        for config_index, config in enumerate(self.saved_configurations):
            row = config_index + 1
            
            # ▲ button
            up_button = QPushButton("▲")
            up_button.setFixedSize(30, 22)
            up_button.setEnabled(config_index > 0)
            up_button.clicked.connect(lambda checked, idx=config_index: self.move_configuration_up(idx))
            grid.addWidget(up_button, row, 0)
            
            # ▼ button
            down_button = QPushButton("▼")
            down_button.setFixedSize(30, 22)
            down_button.setEnabled(config_index < len(self.saved_configurations) - 1)
            down_button.clicked.connect(lambda checked, idx=config_index: self.move_configuration_down(idx))
            grid.addWidget(down_button, row, 1)
            
            # Value labels aligned with joint widgets above
            for j, value in enumerate(config):
                value_label = QLabel(str(round(value, 2)))
                value_label.setAlignment(Qt.AlignCenter)
                value_label.setFixedWidth(55)
                grid.addWidget(value_label, row, j + 2)
            
            # "Set" button
            set_button = QPushButton("Set")
            set_button.setFixedSize(60, 22)
            set_button.clicked.connect(lambda checked, idx=config_index: self.set_configuration(idx))
            grid.addWidget(set_button, row, num_joints + 2)
            
            # Delete button (red X)
            delete_button = QPushButton("✗")
            delete_button.setFixedSize(30, 22)
            delete_button.setStyleSheet("background-color: #FF4444; color: white; font-weight: bold;")
            delete_button.clicked.connect(lambda checked, idx=config_index: self.delete_configuration(idx))
            grid.addWidget(delete_button, row, num_joints + 3)
        
        combined_widget.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
        
        # Center the combined widget, with interpolation/animation to the right
        self.saved_configs_container.addStretch()
        self.saved_configs_container.addWidget(combined_widget)
        
        # Add interpolation slider on the right
        if len(self.saved_configurations) > 1:
            slider_widget = QWidget()
            slider_layout = QVBoxLayout(slider_widget)
            slider_layout.setContentsMargins(10, 0, 10, 0)
            slider_layout.setAlignment(Qt.AlignHCenter)
            
            slider_label = QLabel("Interpolate")
            slider_label.setAlignment(Qt.AlignCenter)
            slider_label.setMinimumWidth(80)
            slider_layout.addWidget(slider_label, 0, Qt.AlignHCenter)
            
            self.config_interp_slider = QSlider(Qt.Vertical)
            self.config_interp_slider.setMinimum(0)
            self.config_interp_slider.setMaximum((len(self.saved_configurations) - 1) * 100)  # 100 steps per config
            
            # Restore previous slider position if it was saved, otherwise default to 0
            if saved_slider_value is not None:
                # Make sure the saved value is within the new range
                saved_slider_value = min(saved_slider_value, self.config_interp_slider.maximum())
                self.config_interp_slider.setValue(saved_slider_value)
            else:
                self.config_interp_slider.setValue(0)
            
            self.config_interp_slider.setTickPosition(QSlider.TicksBothSides)
            self.config_interp_slider.setTickInterval(100)  # Tick at each configuration
            self.config_interp_slider.sliderMoved.connect(self.interpolate_configurations)
            self.config_interp_slider.sliderReleased.connect(self.interpolation_slider_released)
            slider_layout.addWidget(self.config_interp_slider, 0, Qt.AlignHCenter)
            
            self.saved_configs_container.addWidget(slider_widget)
            
            # Add animation controls column on the right
            animation_widget = QWidget()
            animation_layout = QVBoxLayout(animation_widget)
            animation_layout.setContentsMargins(10, 0, 10, 0)
            animation_layout.setAlignment(Qt.AlignTop | Qt.AlignHCenter)
            
            # Add label above play/pause button
            animation_label = QLabel("Animate (s)")
            animation_label.setAlignment(Qt.AlignCenter)
            animation_label.setMinimumWidth(90)
            animation_layout.addWidget(animation_label, 0, Qt.AlignHCenter)
            
            # Create horizontal layout for play button and loop checkbox
            play_loop_layout = QHBoxLayout()
            play_loop_layout.setSpacing(5)
            
            # Play/Pause button
            if self.is_animating:
                self.play_pause_button = QPushButton("⏸")
                self.play_pause_button.setToolTip("Pause")
            else:
                self.play_pause_button = QPushButton("▶")
                self.play_pause_button.setToolTip("Play")
            self.play_pause_button.setFixedSize(30, 24)
            self.play_pause_button.clicked.connect(self.toggle_animation)
            play_loop_layout.addWidget(self.play_pause_button)
            
            # Loop checkbox next to play button
            self.animation_loop_checkbox = QCheckBox("Loop")
            self.animation_loop_checkbox.setMinimumWidth(70)
            self.animation_loop_checkbox.setChecked(self.animation_loop)
            self.animation_loop_checkbox.stateChanged.connect(self.toggle_animation_loop)
            play_loop_layout.addWidget(self.animation_loop_checkbox)
            
            animation_layout.addLayout(play_loop_layout)
            
            # Ensure we have durations list matching the number of segments
            num_segments = len(self.saved_configurations) - 1
            while len(self.config_durations) < num_segments:
                self.config_durations.append(1.0)
            while len(self.config_durations) > num_segments:
                self.config_durations.pop()
            
            # Add duration text boxes aligned with bottom of each interval (configuration row)
            for i in range(num_segments):
                # Add spacing to align with the configuration row at the END of this interval
                if i == 0:
                    animation_layout.addSpacing(0)
                else:
                    animation_layout.addSpacing(2)  # Match configs_column spacing
                
                duration_box = QLineEdit(str(self.config_durations[i]))
                duration_box.setFixedWidth(50)
                duration_box.setAlignment(Qt.AlignCenter)
                duration_box.setToolTip(f"Duration (seconds) from config {i} to {i+1}")
                duration_box.editingFinished.connect(lambda idx=i: self.update_duration(idx))
                animation_layout.addWidget(duration_box, 0, Qt.AlignHCenter)
            
            animation_layout.addStretch()
            self.saved_configs_container.addWidget(animation_widget)
        
        self.saved_configs_container.addStretch()
    
    def interpolate_configurations(self, slider_value):
        """Interpolate between saved configurations based on slider value"""
        if len(self.saved_configurations) < 2:
            return
        
        # Invert because vertical sliders have max at top
        max_value = self.config_interp_slider.maximum()
        inverted_value = max_value - slider_value
        
        # Convert slider value to config space (0.0 to len-1)
        config_value = inverted_value / 100.0
        
        # Determine which two configs to interpolate between
        config_index = int(config_value)
        if config_index >= len(self.saved_configurations) - 1:
            config_index = len(self.saved_configurations) - 2
        
        # Interpolation factor (0.0 to 1.0)
        t = config_value - config_index
        
        config_a = self.saved_configurations[config_index]
        config_b = self.saved_configurations[config_index + 1]
        
        # Interpolate each joint state
        for i in range(min(len(config_a), len(config_b))):
            if i >= len(self.config_joint_indices):
                break
            
            joint_index = self.config_joint_indices[i]
            if joint_index >= len(self.tree.Joints):
                continue
            
            # Linear interpolation
            interpolated_value = config_a[i] * (1 - t) + config_b[i] * t
            
            joint = self.tree.Joints[joint_index]
            if isinstance(joint, Prismatic):
                actualState = interpolated_value
            elif isinstance(joint, Revolute):
                actualState = math.radians(interpolated_value)
            else:
                continue
            
            self.tree.setJointState(joint_index, actualState, lightweight=True)
        
        # Fast path: just move existing GL items instead of full rebuild
        if not self._fast_update_gl():
            # Cache miss — fall back to full redraw (first frame after addToWidget)
            self.tree.resyncFromLightweight()
            self.update_joint()
        self.set_state_tools()
    
    def interpolation_slider_released(self):
        """Handle when interpolation slider is released — rebuild geometry."""
        self._resync_after_lightweight()
        self.log_version()
    
    def set_configuration(self, config_index):
        """Set the tree to a saved configuration"""
        if config_index >= len(self.saved_configurations):
            return
        
        config = self.saved_configurations[config_index]
        
        # Apply each state to the corresponding joint
        for i, value in enumerate(config):
            if i >= len(self.config_joint_indices):
                break
            
            joint_index = self.config_joint_indices[i]
            if joint_index >= len(self.tree.Joints):
                continue
            
            joint = self.tree.Joints[joint_index]
            if isinstance(joint, Prismatic):
                actualState = value
            elif isinstance(joint, Revolute):
                actualState = math.radians(value)
            else:
                continue
            
            self.tree.setJointState(joint_index, actualState)
        
        # Update the interpolation slider to match this configuration
        if hasattr(self, 'config_interp_slider') and \
            not sip.isdeleted(self.config_interp_slider):
            # Calculate slider position for this config index
            # Since slider is inverted: max_value at top (config 0), 0 at bottom (last config)
            max_value = self.config_interp_slider.maximum()
            slider_value = max_value - (config_index * 100)
            self.config_interp_slider.blockSignals(True)
            self.config_interp_slider.setValue(slider_value)
            self.config_interp_slider.blockSignals(False)
        
        self.update_joint()
        self.set_state_tools()
        self.log_version()
    
    def delete_configuration(self, config_index):
        """Delete a saved configuration"""
        if config_index < len(self.saved_configurations):
            self.saved_configurations.pop(config_index)
            self.display_saved_configurations()
    
    def move_configuration_up(self, config_index):
        """Move a configuration up in the list"""
        if config_index > 0 and config_index < len(self.saved_configurations):
            # Swap with the one above
            self.saved_configurations[config_index], self.saved_configurations[config_index - 1] = \
                self.saved_configurations[config_index - 1], self.saved_configurations[config_index]
            self.display_saved_configurations()
    
    def move_configuration_down(self, config_index):
        """Move a configuration down in the list"""
        if config_index >= 0 and config_index < len(self.saved_configurations) - 1:
            # Swap with the one below
            self.saved_configurations[config_index], self.saved_configurations[config_index + 1] = \
                self.saved_configurations[config_index + 1], self.saved_configurations[config_index]
            self.display_saved_configurations()
    
    def update_saved_configs_for_joint_added(self, config_position : int = -1):
        """Insert state 0 into all saved configurations when a joint is added.
        
        config_position: index in the config list where the new joint's state
        should go.  -1 (default) appends at the end.
        """
        for config in self.saved_configurations:
            if config_position < 0 or config_position >= len(config):
                config.append(0.0)
            else:
                config.insert(config_position, 0.0)
    
    def update_saved_configs_for_joint_deleted(self, deleted_joint_index):
        """Remove the state for a deleted joint from all saved configurations"""
        if deleted_joint_index not in self.config_joint_indices:
            return  # Joint was a waypoint, doesn't affect saved configs
        
        # Find which position in config corresponds to this joint
        try:
            config_position = self.config_joint_indices.index(deleted_joint_index)
        except ValueError:
            return
        
        # Remove that position from all saved configurations
        for config in self.saved_configurations:
            if config_position < len(config):
                config.pop(config_position)
    
    def update_duration(self, segment_index):
        """Update the duration for a segment from user input"""
        # Find the duration text box and update the stored value
        try:
            # Find the sender widget
            sender = self.sender()
            if sender and isinstance(sender, QLineEdit):
                try:
                    new_duration = float(sender.text())
                    if new_duration > 0:
                        if segment_index < len(self.config_durations):
                            self.config_durations[segment_index] = new_duration
                    else:
                        sender.setText(str(self.config_durations[segment_index]))
                except ValueError:
                    sender.setText(str(self.config_durations[segment_index]))
        except Exception as e:
            print(f"Error updating duration: {e}")
    
    def toggle_animation(self):
        """Toggle between play and pause states"""
        if self.is_animating:
            self.pause_animation()
        else:
            self.start_animation()
    
    def start_animation(self):
        """Start the animation through configurations"""
        if len(self.saved_configurations) < 2:
            return
        
        self.is_animating = True
        
        # Hide selection sphere and gizmo during animation
        self._pre_anim_selected_joint = self.selected_joint
        self._pre_anim_selected_arrow = self.selected_arrow
        self.selected_joint = -1
        self.selected_arrow = -1
        # Redraw once without gizmo so the GL cache is clean
        self.update_joint()
        
        # Create timer if it doesn't exist
        if self.animation_timer is None:
            self.animation_timer = qc.QTimer()
            self.animation_timer.timeout.connect(self.animation_step)
        
        # Record start position and time
        self.animation_start_value = self.config_interp_slider.value()
        self.animation_start_time = qc.QTime.currentTime().msecsSinceStartOfDay() / 1000.0
        
        # Track which segment we're animating through
        max_value = self.config_interp_slider.maximum()
        inverted_value = max_value - self.animation_start_value
        current_config_value = inverted_value / 100.0
        self.animation_current_segment = int(current_config_value)
        self.animation_segment_start_time = self.animation_start_time
        
        # Update button to pause icon
        self.play_pause_button.setText("⏸")
        self.play_pause_button.setToolTip("Pause")
        
        # Start timer (update every 16ms for ~60fps)
        self.animation_timer.start(16)
    
    def pause_animation(self):
        """Pause the animation"""
        self.is_animating = False
        
        if self.animation_timer is not None:
            self.animation_timer.stop()
        
        # Rebuild full geometry after lightweight animation
        # Restore selection state before resync so gizmo reappears
        if hasattr(self, '_pre_anim_selected_joint'):
            self.selected_joint = self._pre_anim_selected_joint
            self.selected_arrow = self._pre_anim_selected_arrow
        self._resync_after_lightweight()
        
        # Update button to play icon
        self.play_pause_button.setText("▶")
        self.play_pause_button.setToolTip("Play")
    
    def toggle_animation_loop(self, state):
        """Toggle animation loop on/off"""
        self.animation_loop = (state == Qt.Checked)
    
    def animation_step(self):
        """Update animation - called by timer"""
        if not self.is_animating or len(self.saved_configurations) < 2:
            return
        
        current_time = qc.QTime.currentTime().msecsSinceStartOfDay() / 1000.0
        
        # Calculate elapsed time in current segment
        elapsed_in_segment = current_time - self.animation_segment_start_time
        
        # Get duration for current segment
        segment_duration = self.config_durations[self.animation_current_segment] if self.animation_current_segment < len(self.config_durations) else 1.0
        
        # Calculate progress through current segment (0.0 to 1.0)
        if segment_duration > 0:
            segment_progress = elapsed_in_segment / segment_duration
        else:
            segment_progress = 1.0
        
        # Check if we've completed this segment
        if segment_progress >= 1.0:
            # Move to next segment
            self.animation_current_segment += 1
            
            # Check if we've reached the end
            if self.animation_current_segment >= len(self.saved_configurations) - 1:
                # Check if we should loop
                if self.animation_loop:
                    # Loop back to start
                    self.animation_current_segment = 0
                    self.animation_segment_start_time = current_time
                    # Set to first config
                    max_value = self.config_interp_slider.maximum()
                    self.config_interp_slider.blockSignals(True)
                    self.config_interp_slider.setValue(max_value)
                    self.config_interp_slider.blockSignals(False)
                    self.interpolate_configurations(max_value)
                    return
                else:
                    # Stop at the end
                    self.pause_animation()
                    return
            else:
                # Start timing the new segment
                self.animation_segment_start_time = current_time
                segment_progress = 0.0
        
        # Calculate config value (segment_index + progress through segment)
        new_config_value = self.animation_current_segment + segment_progress
        
        # Convert to slider value (inverted)
        max_value = self.config_interp_slider.maximum()
        new_inverted_value = new_config_value * 100
        new_slider_value = max_value - new_inverted_value
        
        # Update slider
        self.config_interp_slider.blockSignals(True)
        self.config_interp_slider.setValue(int(new_slider_value))
        self.config_interp_slider.blockSignals(False)
        
        # Trigger interpolation
        self.interpolate_configurations(int(new_slider_value))

    def initialize_grid(self):
        self.grid = gl.GLGridItem()
        self.grid.setColor(self.grid_color)
        self.grid.setSize(self.grid_size, self.grid_size, self.grid_size)
        self.grid.setSpacing(self.grid_spacing, self.grid_spacing, self.grid_spacing)
        self.log_version()
        self.update_joint()

    def add_to_root_func(self, state):
        self.add_to_root = state == Qt.Checked

    def onUpdateJointState(self, value):
        # regenerate joint
        self.update_joint()

    # def onUpdateRadius(self, value):
    #     # Not used in printed trees
    #     pass
        value = value / 10.0
        self.radius = value
        self.tree.changeRadius(value)
        # self.update_joint() # now called in joint_selection_changed
        self.joint_selection_changed(self.selected_joint, force=True)

    @QtCore.pyqtSlot(float)
    def change_mesh_scale(self, scale):
        self.mesh_scale = scale

    def log_version(self):
        log_capacity = 100 #this is what's stored in runtime for undo/redo
        autosave_frequency = 1 #saving everything for analysis: in other circumstances this should be more like 10

        # clear redo history on new version (include version index)
        self.versions = self.versions[:self.version_index + 1]

        if self.total_version_counter % autosave_frequency == 0 and not self.tree is None:
            #self.save_tree(autosave_id=len(self.versions)//autosave_frequency)
            self.save_tree(autosave_id=time.time())  #autosave with timestamp (seconds from unix epoch start)
        if len(self.versions) < log_capacity:
            self.versions.append(copy.deepcopy(self.tree))
        else:
            self.versions.pop(0)
            self.versions.append(copy.deepcopy(self.tree))

        self.version_index = len(self.versions) - 1
        self.total_version_counter += 1

        #print("LOG     version: " + str(self.total_version_counter) + ", size: " + str(len(self.versions)) + ", index: " + str(self.version_index))

    def undo(self):
        if self.version_index > 0:
            self.version_index -= 1
            self.tree = copy.deepcopy(self.versions[self.version_index])
        else:
            if self.version_index == 0:
                self.version_index = -1
            self.tree = None

        self._mark_collision_dirty()
        self.update_joint()
        #print("UNDO    version: " + str(self.total_version_counter) + ", size: " + str(len(self.versions)) + ", index: " + str(self.version_index))

    def redo(self):
        if self.version_index + 1 < len(self.versions):
            self.version_index += 1
            self.tree = copy.deepcopy(self.versions[self.version_index])
            self._mark_collision_dirty()
            self.update_joint()

        #print("REDO    version: " + str(self.total_version_counter) + ", size: " + str(len(self.versions)) + ", index: " + str(self.version_index))

    def _toggle_collision_highlighting(self, checked):
        self.show_collisions = checked
        self.update_joint()

    def toggle_grid_func(self, checked=None):
        if checked is None:
            # Called from keyboard shortcut — toggle the checkbox, which re-enters this method
            self.toggle_grid.setChecked(not self.toggle_grid.isChecked())
            return
        if checked:
            self.plot_widget.addItem(self.grid)
        else:
            self.plot_widget.removeItem(self.grid)
        self.grid_on = checked

    # Success message method with timer
    def show_success(self, message):
        self.status_label.setText(message)
        self.status_label.setStyleSheet("color: " + successColorDefault)
        self.message_display_dock.setVisible(True)
        self.message_display_dock.raise_()

    # Error message method with timer
    def show_error(self, message):
        self.status_label.setText(message)
        self.status_label.setStyleSheet("color: " + errorColorDefault)
        self.message_display_dock.setVisible(True)
        self.message_display_dock.raise_()

    # Method to clear the message
    def clear_message(self):
        self.status_label.setText('')

    def show_delete_widget(self):
        self.delete_joint_dock.setVisible(True)

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

    def set_joint_as_frame(self):
        self.selected_frame = self.selected_joint
        self.frame_label.setText("Joint selected as frame: " + str(self.selected_frame))
        self.local_orient_slider_checkbox.setDisabled(True)
        self.update_joint()

    def remove_frame_button(self):
        self.selected_frame = -1
        self.frame_label.setText("Joint selected as frame: N/A")
        self.local_orient_slider_checkbox.setDisabled(False)
        self.update_joint()

    def init_key_bar(self):
        # Add instructions as a single centered, word-wrapping label
        instructions = (
            "Middle Mouse / Shift+Drag: Pan Camera  |  "
            "T: Translate  |  R: Rotate  |  Delete: Delete Joint  |  "
            "X: Select X Axis  |  Y: Select Y Axis  |  Z: Select Z Axis"
        )
        # ", Ctrl+Z: Undo",  "Ctrl+Y / Ctrl+Shift+Z: Redo"
        # Removed because it's standard enough to guess
        # And it's different Ctrl vs Cmd on Windows vs Mac
        # And we want to save space

        label = QLabel(instructions)
        label.setAlignment(Qt.AlignCenter)
        label.setWordWrap(False)
        self.key_bar_layout.addWidget(label)

    def local_orient_clicked(self):
        self.is_local = self.local_orient_slider_checkbox.isChecked()
        self.update_joint()

    def change_control_type(self):
        if self.translate_joint_radio_button.isChecked():
            self.control_type = "Translate"
        elif self.rotate_joint_radio_button.isChecked():
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
            self.translate_joint_radio_button.setChecked(True)
        elif key == "Rotate":
            self.control_type = key
            self.rotate_joint_radio_button.setChecked(True)
        elif key == "Delete":
            if self.tree and self.selected_joint != -1:
                self.delete_selected_joint()
        elif key == "Escape":
            if self.measure_select_mode:
                self.measure_select_mode = False
                self.measurement_widget.deactivate()
        elif key == "X":
            self.arrow_selection_changed(0)
        elif key == "Y":
            self.arrow_selection_changed(1)
        elif key == "Z":
            self.arrow_selection_changed(2)
        elif key == "G":
            self.toggle_grid_func()

        self.update_joint()

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

    @QtCore.pyqtSlot(bool)
    def mesh_selected_slot(self, is_selected):
        if is_selected != self.mesh_selected:
            self.mesh_selected = is_selected
            self.update_joint()

    @QtCore.pyqtSlot(int)
    def joint_selection_changed(self, index, force : bool = False):
        if force or index != self.selected_joint:
            self.selected_joint = index
            self.selected_arrow = -1
            self.selected_axis_name = 'N/A'
            self.update_joint()
            self.reset_rotation_tools()
            self.reset_translation_tools()
            self.set_state_tools()
            self.highlight_selected_config_box()
            # Sync the select parent prompt dropdown if visible
            if self.select_parent_prompt.isVisible() and index >= 0:
                self.select_parent_combo.blockSignals(True)
                self.select_parent_combo.setCurrentIndex(index)
                self.select_parent_combo.blockSignals(False)

    @QtCore.pyqtSlot(int)
    def arrow_selection_changed(self, index):
        if index != self.selected_arrow and (self.selected_joint != -1 or self.mesh_selected):
            self.selected_arrow = index
                        
            if (index == 0): self.selected_axis_name = 'X'
            elif (index == 1): self.selected_axis_name = 'Y'
            elif (index == 2): self.selected_axis_name = 'Z'
            else: self.selected_axis_name = 'N/A'

            self.update_joint()

            if (self.selected_joint != -1):
                self.reset_rotation_tools()
                self.reset_translation_tools()

    @QtCore.pyqtSlot(int)
    def link_selection_changed(self, index):
        self.select_link_options.setCurrentIndex(self.selected_link)
        if index != self.selected_link:
            self.selected_link = index
            self.selected_arrow = -1
            self.selected_axis_name = 'N/A'
            self.update_joint()
            self.reset_rotation_tools()
            self.reset_translation_tools()
            self.set_state_tools()
            #self.update_radius_slider()

    @QtCore.pyqtSlot(np.ndarray)
    def drag_translate(self, new_position):
        propogate = self.propogate_slider_checkbox.isChecked()

        if (self.mesh_selected):
            old_position = self.referenceMesh.Pose.t
            trans = new_position - old_position
            transformation = SE3.Trans(trans[0], trans[1], trans[2])
            self.referenceMesh.Pose = transformation * self.referenceMesh.Pose

            transform_matrix = QMatrix4x4()
            matrix = self.referenceMesh.Pose.A
            for row in range(4):
                for col in range(4):
                    transform_matrix[row, col] = matrix[row, col]

            transform_matrix[0, 0] = self.mesh_scale
            transform_matrix[1, 1] = self.mesh_scale
            transform_matrix[2, 2] = self.mesh_scale
            self.referenceMesh.mesh.setTransform(transform_matrix)

            self.update_joint()
        else:
            old_position = self.tree.Joints[self.selected_joint].Pose.t
            trans = new_position - old_position
            transformation = SE3.Trans(trans[0], trans[1], trans[2])

            if self.tree.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=False):
                self.update_joint()

            self.reset_rotation_tools()
            self.reset_translation_tools()

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

    @QtCore.pyqtSlot(float)
    def drag_rotate(self, new_rotation):
        transformation = SE3()
        propogate = self.propogate_slider_checkbox.isChecked()

        if (self.selected_axis_name == 'X'):
            transformation = SE3.Rx(new_rotation)
        elif (self.selected_axis_name == 'Y'):
            transformation = SE3.Ry(new_rotation)
        elif (self.selected_axis_name == 'Z'):
            transformation = SE3.Rz(new_rotation)
        
        if (self.mesh_selected):
            self.referenceMesh.Pose = self.referenceMesh.Pose * transformation

            transform_matrix = QMatrix4x4()
            matrix = self.referenceMesh.Pose.A
            for row in range(4):
                for col in range(4):
                    transform_matrix[row, col] = matrix[row, col]
            self.referenceMesh.mesh.setTransform(transform_matrix)

            self.update_joint()
        else: 
            if self.tree.transformJoint(self.selected_joint, transformation, propogate=propogate, safe=True, relative=True):
                self.update_joint()
                self.reset_rotation_tools()
                self.reset_translation_tools()

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

    def rotation_angle_from_matrix(self, rotation_matrix, axis):
        rot = R.from_matrix(rotation_matrix)
        euler_angles = rot.as_euler('xyz', degrees=True)
        return euler_angles[axis]
    
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

    def rotation_textbox_return(self):
        try:
            value = float(self.rotation_textbox.text())
        except ValueError:
            value = 0
        self.adjust_rotation(value)
        self.reset_rotation_tools()
        self.log_version()

    def rotation_slider_released(self):
        self.reset_rotation_tools()
        self.log_version()
    
    def reset_rotation_tools(self):
        r = self.tree.r if self.tree else 1
        self.rotation_slider.setMinimum(-360)
        self.rotation_slider.setMaximum(360)
        self.rotation_slider.setValue(0)
        self.old_rot_val = 0
        self.rotation_textbox.setText("0")
        if self.selected_arrow == -1 or self.selected_joint == -1:
            self.rotation_slider.setDisabled(True)
            self.rotation_textbox.setDisabled(True)
        else:
            self.rotation_slider.setDisabled(False)
            self.rotation_textbox.setDisabled(False) 
    
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
    
    def translation_textbox_return(self):
        try:
            value = float(self.translation_textbox.text()) * 10
        except ValueError:
            value = 0
        self.adjust_translation(value)
        self.reset_translation_tools()
        self.log_version()
    
    def translation_slider_released(self):
        self.reset_translation_tools()
        self.log_version()
    
    def reset_translation_tools(self):
        r = self.tree.r if self.tree else 1
        self.translation_slider.setMinimum(int(-100*r))
        self.translation_slider.setMaximum(int(100*r))
        self.translation_slider.setValue(0)
        self.old_trans_val = 0
        self.translation_textbox.setText("0")
        if self.selected_arrow == -1 or self.selected_joint == -1:
            self.translation_slider.setDisabled(True)
            self.translation_textbox.setDisabled(True)
        else:
            self.translation_slider.setDisabled(False)
            self.translation_textbox.setDisabled(False)

    
    def state_slider_moved(self, value):
        if self.tree and self.selected_joint != -1:
            if isinstance(self.tree.Joints[self.selected_joint], Prismatic):
                # Prismatic joints: slider value is scaled by 100*r
                actualState = value / (100 * self.tree.r)
            elif isinstance(self.tree.Joints[self.selected_joint], Revolute):
                # Revolute joints: slider value is in degrees, actual state is radians
                actualState = math.radians(value)
            else: # Waypoint (but it shouldn't let you move the slider in the first place in that case)
                print("Warning: Tried to move state slider on a waypoint, which should not be possible.")
                return
            if self.tree.setJointState(self.selected_joint, actualState, lightweight=True):
                if not self._fast_update_gl():
                    self.tree.resyncFromLightweight()
                    self.update_joint()
                self.set_state_tools()
            else:
                self.state_slider.blockSignals(True)
                self.state_slider.setValue(int(self.old_state_slider_val))
                self.state_slider.blockSignals(False)
    
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

    def state_textbox_return(self):
        if self.tree and self.selected_joint != -1:
            try:
                value = float(self.state_textbox.text())
            except ValueError:
                return #TODO: handle this with QLineEdit class's setValidator method instead
            if isinstance(self.tree.Joints[self.selected_joint], Prismatic):
                actualState = value
            elif isinstance(self.tree.Joints[self.selected_joint], Revolute):
                actualState = math.radians(value)
            else:
                print("Warning: Tried to edit state textbox on a waypoint, which should not be possible.")
                return
            if self.tree.setJointState(self.selected_joint, actualState):
                self._mark_collision_dirty()
                self.update_joint()
                self.set_state_tools()
                self.log_version()

    def eventFilter(self, obj, event):
        """Intercept Tab in config text boxes: apply value and advance to next."""
        if event.type() == QEvent.KeyPress and event.key() == Qt.Key_Tab:
            if hasattr(self, 'config_text_boxes') and obj in self.config_text_boxes:
                idx = self.config_text_boxes.index(obj)
                joint_index = self.config_joint_indices[idx]
                self.config_textbox_return(joint_index)
                # Move focus to next text box if not at the end
                if idx + 1 < len(self.config_text_boxes):
                    self.config_text_boxes[idx + 1].setFocus()
                    self.config_text_boxes[idx + 1].selectAll()
                return True  # consume the event
        return super().eventFilter(obj, event)

    def config_textbox_return(self, joint_index):
        """Handle configuration textbox input for a specific joint"""
        if self.tree and 0 <= joint_index < len(self.tree.Joints):
            # Find the text box index for this joint
            try:
                text_box_index = self.config_joint_indices.index(joint_index)
                text_box = self.config_text_boxes[text_box_index]
            except (ValueError, IndexError):
                return
            
            try:
                value = float(text_box.text())
            except ValueError:
                return
            
            joint = self.tree.Joints[joint_index]
            if isinstance(joint, Prismatic):
                actualState = value
            elif isinstance(joint, Revolute):
                actualState = math.radians(value)
            else:
                # This shouldn't happen since we skip waypoints
                return
            
            if self.tree.setJointState(joint_index, actualState):
                self._mark_collision_dirty()
                self.update_joint()
                self.set_state_tools()
                self.log_version()

    def config_slider_moved(self, joint_index, value):
        """Handle configuration slider movement for a specific joint"""
        if self.tree and 0 <= joint_index < len(self.tree.Joints):
            joint = self.tree.Joints[joint_index]
            
            if isinstance(joint, Prismatic):
                # Prismatic joints: slider value is scaled by 100*r
                actualState = value / (100 * self.tree.r)
                display_state = actualState
            elif isinstance(joint, Revolute):
                # Revolute joints: slider value is in degrees, actual state is radians
                actualState = math.radians(value)
                display_state = value
            else:
                return
            
            if self.tree.setJointState(joint_index, actualState, lightweight=True):
                # Fast GL update during slider drag
                if not self._fast_update_gl():
                    self.tree.resyncFromLightweight()
                    self.update_joint()
                # Update the corresponding text box
                try:
                    text_box_index = self.config_joint_indices.index(joint_index)
                    self.config_text_boxes[text_box_index].blockSignals(True)
                    self.config_text_boxes[text_box_index].setText(str(round(display_state, 2)))
                    self.config_text_boxes[text_box_index].blockSignals(False)
                except (ValueError, IndexError):
                    pass
                self.set_state_tools()
    
    def config_slider_released(self, joint_index):
        """Handle configuration slider release for a specific joint"""
        self._resync_after_lightweight()
        self.set_state_tools()
        self.log_version()

    def state_slider_released(self):
        self._resync_after_lightweight()
        self.set_state_tools()
        self.log_version()
    
    # Returns tuple of tuples: each inner tuple is (min, max, and state) values
    # of a particular version of joint state information: actual, slider, and textbox in that order. 
    # Revolute joints are actually in radians, but used in degrees in the slider and text box.
    # Prismatic joints are actually in distance units and displayed as such in the text box, 
    # but the slider scales it by 100*r because it needs integer values.
    # Returns None if the joint is a waypoint, no joint is selected, or no tree exists.
    # The state used is the current joint state unless textboxStateInput is provided
    # (in which case it is interpreted as degrees for revolute joints and distance units for prismatic joints).
    def scaled_state_info(self, state=None):
        if self.tree and self.selected_joint != -1:
            joint = self.tree.Joints[self.selected_joint]
            stateRange = self.tree.Joints[self.selected_joint].stateRange()
            if isinstance(joint, Prismatic):
                stateActual = self.tree.Joints[self.selected_joint].state if state is None else max(stateRange[0], min(stateRange[1], state))
                actual = (stateRange[0], stateRange[1], stateActual)
                scale = 100 * self.tree.r
                scaled = (int(stateRange[0] * scale), int(stateRange[1] * scale), int(stateActual * scale))
                actual, slider, textbox = actual, scaled, actual
                return (actual, slider, textbox)
            elif isinstance(joint, Revolute):
                if state is None:
                    stateRadians = self.tree.Joints[self.selected_joint].state
                    stateDegrees = math.degrees(stateRadians)
                else:
                    stateDegrees = state
                    stateRadians = math.radians(stateDegrees)
                radians = (stateRange[0], stateRange[1], stateRadians)
                degrees = (round(math.degrees(stateRange[0])), round(math.degrees(stateRange[1])), round(stateDegrees))
                actual, slider, textbox = radians, degrees, degrees
                return (actual, slider, textbox)
            elif isinstance(joint, Waypoint):
                return None
        else:
            return None
        
    
    def set_state_tools(self):
        scaledInfo = self.scaled_state_info()
        if scaledInfo is None:
            self.state_slider.setDisabled(True)
            self.state_textbox.setDisabled(True)
            self.state_slider.setMinimum(0)
            self.state_slider.setMaximum(0)
            self.state_slider.setValue(0)
            self.state_textbox.setText("0")
        else:
            self.state_slider.setDisabled(False)
            self.state_textbox.setDisabled(False)
            actual, slider, textbox = scaledInfo
            minText, maxText, currentText = textbox
            minSlider, maxSlider, currentSlider = slider
            self.state_slider.setMinimum(minSlider)
            self.state_slider.setMaximum(maxSlider)
            self.state_slider.setValue(currentSlider)
            self.old_state_slider_val = currentSlider
            decimals = 2 if isinstance(self.tree.Joints[self.selected_joint], Prismatic) else 0
            self.state_textbox.setText(str(np.round(currentText, decimals)))
            self.current_state_label.setText(f"Min State: {np.round(minText, decimals)} ≤ Current State: {np.round(currentText, decimals)} ≤ Max State: {np.round(maxText, decimals)}")

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

    def rotate_joint(self, angle, axis):
        # Update the position of the spheres
        joint = self.selected_joint

        propogate = self.propogate_slider_checkbox.isChecked()
        localOrient = self.local_orient_slider_checkbox.isChecked()

        transformation = SE3.AngleAxis(angle, [axis[0], axis[1], axis[2]], unit='deg')
        # transformation = SE3.Trans(0, 1, 0)

        # print(transformation)

        self.tree.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=True, safe=True, localOrient=localOrient)

        # self.selected_joint.translate(-cnt[0], -cnt[1], -cnt[2])
        # self.selected_joint.rotate(angle, axis[0], axis[1], axis[2], local=False)
        # self.selected_joint.translate(cnt[0], cnt[1], cnt[2])

        # for i, a in enumerate(self.axes):
        #     a.translate(-cnt[0], -cnt[1], -cnt[2])
        #     a.rotate(angle, axis[0], axis[1], axis[2], local=False)
        #     a.translate(cnt[0], cnt[1], cnt[2])

    def _fast_update_gl(self) -> bool:
        """Try to update GL items via model-matrix transforms (fast path).
        Returns True if the fast path succeeded, False if a full redraw is needed.
        Skips collision detection, sidebar updates, arrows, etc."""
        if self.tree is None:
            return False
        if not self.tree.updateGLTransforms():
            return False
        # On the first lightweight frame, clear collision highlighting.
        # Link geometry is stale during lightweight mode so collision
        # results would be incorrect; colors are restored on resync.
        if not self._collision_colors_cleared:
            if getattr(self, 'show_collisions', True):
                from style import linkColorDefault
                for link in self.tree.Links:
                    for item in getattr(link, '_gl_items', []):
                        if hasattr(item, 'setColor'):
                            item.setColor(linkColorDefault)
            self._collision_colors_cleared = True
        self._lightweight_dirty = True
        self.plot_widget.update()
        return True

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

    def _resync_after_lightweight(self):
        """Rebuild link geometry after lightweight animation, then do a full
        update_joint so collision capsules, bounding balls, and sidebar
        widgets are all consistent."""
        if self.tree is not None:
            self.tree.resyncFromLightweight()
        self._lightweight_dirty = False
        self._collision_colors_cleared = False
        self._mark_collision_dirty()
        self.update_joint()

    def _on_measurement_dock_visibility(self, visible: bool):
        if not visible and self.measure_select_mode:
            self.measure_select_mode = False
            self.measurement_widget.deactivate()

    def _toggle_measure_select(self, enabled: bool):
        self.measure_select_mode = enabled
        if enabled:
            self.selected_joint = -1
        self.update_joint()

    def _set_measure_active_endpoint(self, ep: int):
        self.measure_active_endpoint = ep
        self.measurement_widget.set_active_endpoint(ep)

    def _measure_joint_selected_slot(self, joint_index: int):
        if joint_index == -1 or self.tree is None:
            return
        ep = self.measure_active_endpoint
        frame_type = self.measurement_widget.frame_types[ep]
        self.measure_points[ep] = (joint_index, frame_type)
        pos = self._get_measure_position(joint_index, frame_type)
        joint = self.tree.Joints[joint_index]
        self.measurement_widget.set_point_info(ep, joint_index, type(joint).__name__, pos)
        next_ep = 1 - ep
        self.measure_active_endpoint = next_ep
        self.measurement_widget.set_active_endpoint(next_ep)
        self.update_joint()

    def _measure_frame_type_changed(self):
        if self.tree is None:
            return
        for ep in (0, 1):
            if self.measure_points[ep] is not None:
                joint_idx = self.measure_points[ep][0]
                new_frame = self.measurement_widget.frame_types[ep]
                self.measure_points[ep] = (joint_idx, new_frame)
                pos = self._get_measure_position(joint_idx, new_frame)
                joint = self.tree.Joints[joint_idx]
                self.measurement_widget.set_point_info(ep, joint_idx, type(joint).__name__, pos)
        self.update_joint()

    def _clear_measurement(self):
        self.measure_points = [None, None]
        self.measure_select_mode = False
        self.measure_active_endpoint = 0
        self.measurement_widget.reset()
        self.update_joint()

    def _get_measure_position(self, joint_index: int, frame_type: str):
        if self.tree is None:
            return None
        joints = self.tree.Joints
        if not (0 <= joint_index < len(joints)):
            return None
        joint = joints[joint_index]
        if frame_type == 'center':
            return np.array(joint.Pose.t)
        elif frame_type == 'proximal':
            return np.array(joint.proximalPosition())
        elif frame_type == 'distal':
            return np.array(joint.distalPosition())
        return None

    def _refresh_measurement_overlay(self):
        if getattr(self, '_measurement_line_item', None) is not None:
            try:
                self.plot_widget.removeItem(self._measurement_line_item)
            except Exception:
                pass
            self._measurement_line_item = None
        _show = (
            self.measure_points[0] is not None
            and self.measure_points[1] is not None
        )
        if _show:
            pos_a = self._get_measure_position(*self.measure_points[0])
            pos_b = self._get_measure_position(*self.measure_points[1])
            if pos_a is not None and pos_b is not None:
                pts = np.array([pos_a, pos_b], dtype=float)
                line = OverlayLine(pos=pts, color=(1.0, 1.0, 0.0, 1.0), width=3, antialias=True)
                self.plot_widget.addItem(line)
                self._measurement_line_item = line
                displacement = pos_b - pos_a
                dist = float(np.linalg.norm(displacement))
                ax = self.measurement_widget.axis_frame
                if ax == 'local_a':
                    R_mat = self.tree.Joints[self.measure_points[0][0]].Pose.R
                    components = R_mat.T @ displacement
                elif ax == 'local_b':
                    R_mat = self.tree.Joints[self.measure_points[1][0]].Pose.R
                    components = R_mat.T @ displacement
                else:
                    components = displacement
                for ep, mp in enumerate(self.measure_points):
                    pos = self._get_measure_position(*mp)
                    joint = self.tree.Joints[mp[0]]
                    self.measurement_widget.set_point_info(ep, mp[0], type(joint).__name__, pos)
                self.measurement_widget.set_result(dist, components)
            else:
                self.measurement_widget.set_result(None, None)
        else:
            self.measurement_widget.set_result(None, None)

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