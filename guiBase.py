import sys, os
import math
import copy
import time
import numpy as np
import pyqtgraph.opengl as gl
import pyqtgraph as pg
from PyQt5 import QtCore as qc, sip
from PyQt5.QtWidgets import (
    QWidget, QPushButton, QDockWidget, QHBoxLayout, QVBoxLayout,
    QLabel, QLineEdit, QCheckBox, QMessageBox, QRadioButton, QSlider, QFileDialog,
    QApplication, QMainWindow, QComboBox, QSizePolicy, QShortcut, QGridLayout, QSpinBox
)
from PyQt5.QtCore import Qt, pyqtSignal, QEvent
from PyQt5.QtGui import QPixmap, QIcon, QSurfaceFormat, QVector3D, QKeyEvent, QKeySequence, QMatrix4x4
from pyqtgraph.Qt import QtCore
from spatialmath import SE3
from scipy.spatial.transform import Rotation as R
from ReferenceMesh import *
from IntersectionHelper import *
from Joint import Prismatic, Revolute, Waypoint


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
        except ValueError:
            self.show_error("Please enter a valid file path and scale factor.")

    def on_clear_clicked(self):
        self.window().referenceMesh = None
        self.window().update_joint()

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

        apply_button = QPushButton('Apply Changes', self)
        apply_button.clicked.connect(self.on_apply_clicked)
        button_layout.addWidget(apply_button)

        self.cancel_button = QPushButton('Cancel', self)
        self.cancel_button.clicked.connect(self.on_cancel_clicked)
        button_layout.addWidget(self.cancel_button)

        layout.addLayout(button_layout)

        self.setLayout(layout)

    def on_apply_clicked(self):
        try:
            spacing = int(self.spacing_input.text())
            line_amt = int(self.line_amt_input.text())

            self.window().grid_size = line_amt * spacing
            self.window().grid_spacing = spacing
            self.window().initialize_grid()

        except ValueError:
            self.show_error("Please enter valid integers.")

    def on_cancel_clicked(self):
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
        import importlib.util
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


class BaseClickableGLViewWidget(gl.GLViewWidget):
    def __init__(self, parent_window, parent=None):
        super().__init__(parent)
        fmt = QSurfaceFormat()
        fmt.setDepthBufferSize(24)
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

    # --- Abstract / hook methods ---

    def get_world_coordinates(self, event):
        raise NotImplementedError

    def _should_pick_items(self) -> bool:
        """Return True if itemsAt should be called on this click."""
        return True

    def _handle_picked_links(self, links) -> None:
        """Called with the list of picked 'Link' items after itemsAt."""
        pass

    def _mpe_profiling_context(self):
        """Return an opaque timing context (or None) for mousePressEvent profiling."""
        return None

    def _mpe_profiling_log(self, ctx, picked_items, skipped: bool) -> None:
        """Log profiling info after itemsAt (or when skipped)."""
        pass

    def _on_drag_release(self) -> None:
        """Called at the start of mouseReleaseEvent when a drag ends."""
        pass

    # --- Identical methods ---

    def toggle_lock(self):
        self.locked = not self.locked
        self.lock_status_changed.emit(self.locked)
        print("Screen lock toggled:", "Locked" if self.locked else "Unlocked")

    def get_normalized_plane_vectors(self, event):
        if self.is_dragging and self.selected_torus:
            origin, dir = self.get_world_coordinates(event)
            selected_joint = self.parent_window.model.Joints[self.parent_window.selected_joint]
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
            selected_joint = self.parent_window.model.Joints[self.parent_window.selected_joint]
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

            if self.mesh_selected:
                selected_joint = self.parent_window.referenceMesh
            else:
                selected_joint = self.parent_window.model.Joints[self.parent_window.selected_joint]

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
        _ctx = self._mpe_profiling_context()

        lpos = event.position() if hasattr(event, 'position') else event.localPos()
        region = [lpos.x()-5, lpos.y()-5, 10, 10]
        dpr = self.devicePixelRatioF()
        region = tuple([x * dpr for x in region])

        links = []
        mesh = []

        if not self.is_dragging and self._should_pick_items():
            import io
            _old_stdout = sys.stdout
            sys.stdout = io.StringIO()
            try:
                picked_items = list(self.itemsAt(region))
            finally:
                sys.stdout = _old_stdout

            self._mpe_profiling_log(_ctx, picked_items, skipped=False)

            for item in picked_items:
                if item.objectName() == "Link":
                    links.append(item)
                if item.objectName() == "Mesh":
                    mesh.append(item)

            self._handle_picked_links(links)
        else:
            self._mpe_profiling_log(_ctx, None, skipped=True)

        if not self.is_dragging:
            self.mesh_selected = len(mesh) > 0

        self.click_signal_link.emit(self.selected_link_index)
        self.click_signal_mesh.emit(self.mesh_selected)

        self.last_drag_pos = event.pos()

        if event.buttons() and event.button() == Qt.MouseButton.MiddleButton:
            self.drag_start_pos = event.pos()

        if event.buttons() and Qt.LeftButton and event.buttons() != QtCore.Qt.MouseButton.MiddleButton:
            self.is_dragging = False
            self.drag_start_pos = event.pos()

            if self.parent_window.measure_select_mode:
                origin, direction = self.get_world_coordinates(event)
                closest_joint_dist = float('inf')
                closest_joint = None
                if self.parent_window.model:
                    for joint in self.parent_window.model.Joints:
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

            if self.mesh_selected and self.parent_window.selected_joint == -1:
                self.selected_axis = None
                self.selected_torus = None
                mesh_obj = self.window().referenceMesh

                mesh_center = mesh_obj.Pose.t
                mesh_center = QVector3D(mesh_center[0], mesh_center[1], mesh_center[2])
                self.compute_joint_axes(mesh_obj)

                if self.parent_window.control_type == "Translate":
                    self.calculate_translation_isect(origin, direction, mesh_center, mesh_obj, event)
                elif self.parent_window.control_type == "Rotate":
                    self.calculate_rotation_isect(origin, direction, mesh_center, mesh_obj)

            elif self.parent_window.selected_joint != -1:
                self.selected_axis = None
                self.selected_torus = None
                selected_joint = self.parent_window.model.Joints[self.parent_window.selected_joint]
                joint_center = selected_joint.Pose.t
                joint_center = QVector3D(joint_center[0], joint_center[1], joint_center[2])

                self.compute_joint_axes(selected_joint)

                if self.parent_window.control_type == "Translate":
                    self.calculate_translation_isect(origin, direction, joint_center, selected_joint, event)
                elif self.parent_window.control_type == "Rotate":
                    self.calculate_rotation_isect(origin, direction, joint_center, selected_joint)

            if self.selected_torus:
                self.is_dragging = True
                self.drag_prev_vector, _ = self.get_normalized_plane_vectors(event)

            closest_joint_dist = float('inf')
            closest_joint = None

            if self.parent_window.model:
                for joint in self.parent_window.model.Joints:
                    center = joint.Pose.t
                    radius = joint.r
                    hit_location = compute_sphere_intersection(origin, direction, center, 1.1 * radius)
                    if hit_location < closest_joint_dist:
                        closest_joint_dist = hit_location
                        closest_joint = joint

            self.selected_joint_temp = closest_joint

    def mouseReleaseEvent(self, event):
        if self.parent_window.measure_select_mode:
            was_dragging = self.is_dragging
            self.is_dragging = False
            if not was_dragging and self.selected_joint_temp is not None:
                joint_index = -1
                for i, j in enumerate(self.parent_window.model.Joints):
                    if j is self.selected_joint_temp:
                        joint_index = i
                        break
                self.measure_joint_selected.emit(joint_index)
            self.selected_joint_temp = None
            return

        if self.is_dragging and (self.selected_axis or self.selected_torus):
            self._drag_backup = None
            self._on_drag_release()
            self.done_transforming.emit(True)

        if self.selected_joint_temp is not None:
            joint_index = -1
            for i, joint in enumerate(self.parent_window.model.Joints):
                if joint is self.selected_joint_temp:
                    joint_index = i
                    break
            self.click_signal.emit(joint_index)
        elif not self.is_dragging:
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

    def _handle_camera_drag(self, event):
        """Shared camera orbit/pan handler for the no-joint-selected drag case.
        Updates last_drag_pos, then dispatches: middle-button or Shift+left → pan;
        left → orbit or pan depending on camera_type."""
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


class BaseWindowKinegamiGUI(QMainWindow):
    """Base class for both GUI windows. Subclasses provide model property and
    fabrication-specific methods; this class holds all shared window logic."""

    LOG_CAPACITY = 20
    AUTOSAVE_FREQUENCY = 10

    # ── Abstract model property ───────────────────────────────────────────────

    @property
    def model(self):
        raise NotImplementedError

    @model.setter
    def model(self, v):
        raise NotImplementedError

    @property
    def model_created(self) -> bool:
        return self.model is not None

    # ── Hook methods (no-op defaults; subclasses override as needed) ──────────

    def _mark_collision_dirty(self) -> None:
        pass

    def _post_undo_hook(self, prev_model) -> None:
        pass

    def _show_message_dock(self) -> None:
        self.message_display_dock.setVisible(True)
        self.message_display_dock.raise_()

    def _on_joint_selection_changed_extra(self, index: int) -> None:
        pass

    def _on_fast_update_gl_extra(self) -> None:
        pass

    def _save_model(self, autosave_id=None) -> None:
        raise NotImplementedError

    def _should_skip_joint_in_config(self, joint) -> bool:
        return isinstance(joint, Waypoint)

    # ── Status messages ───────────────────────────────────────────────────────

    def show_success(self, message):
        self.status_label.setText(message)
        self.status_label.setStyleSheet("color: " + successColorDefault)
        self._show_message_dock()

    def show_error(self, message):
        self.status_label.setText(message)
        self.status_label.setStyleSheet("color: " + errorColorDefault)
        self._show_message_dock()

    def clear_message(self):
        self.status_label.setText('')

    # ── Grid ──────────────────────────────────────────────────────────────────

    def edit_grid_func(self):
        visibility = self.edit_grid_dock.isVisible()
        self.edit_grid_dock.setVisible(not visibility)

    def initialize_grid(self):
        self.grid = gl.GLGridItem()
        self.grid.setColor(self.grid_color)
        self.grid.setSize(self.grid_size, self.grid_size, self.grid_size)
        self.grid.setSpacing(self.grid_spacing, self.grid_spacing, self.grid_spacing)
        self.log_version()
        self.update_joint()

    def toggle_grid_func(self, checked=None):
        if checked is None:
            self.toggle_grid.setChecked(not self.toggle_grid.isChecked())
            return
        if checked:
            self.plot_widget.addItem(self.grid)
        else:
            self.plot_widget.removeItem(self.grid)
        self.grid_on = checked

    # ── UI widget callbacks ───────────────────────────────────────────────────

    def add_to_root_func(self, state):
        self.add_to_root = state == Qt.Checked

    def onUpdateJointState(self, value):
        self.update_joint()

    @QtCore.pyqtSlot(float)
    def change_mesh_scale(self, scale):
        self.mesh_scale = scale

    def init_key_bar(self):
        instructions = (
            "Middle Mouse / Shift+Drag: Pan Camera  |  "
            "T: Translate  |  R: Rotate  |  Delete: Delete Joint  |  "
            "X: Select X Axis  |  Y: Select Y Axis  |  Z: Select Z Axis"
        )
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

    def show_delete_widget(self):
        self.delete_joint_dock.setVisible(True)

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

    # ── Keyboard / key bar ────────────────────────────────────────────────────

    @QtCore.pyqtSlot(str)
    def key_pressed(self, key):
        if key == "Translate":
            self.control_type = key
            self.translate_joint_radio_button.setChecked(True)
        elif key == "Rotate":
            self.control_type = key
            self.rotate_joint_radio_button.setChecked(True)
        elif key == "Delete":
            if self.model and self.selected_joint != -1:
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

    # ── Selection changed ─────────────────────────────────────────────────────

    @QtCore.pyqtSlot(int)
    def joint_selection_changed(self, index, force: bool = False):
        if force or index != self.selected_joint:
            self.selected_joint = index
            self.selected_arrow = -1
            self.selected_axis_name = 'N/A'
            self.update_joint()
            self.reset_rotation_tools()
            self.reset_translation_tools()
            self.set_state_tools()
            self.highlight_selected_config_box()
            self._on_joint_selection_changed_extra(index)

    @QtCore.pyqtSlot(int)
    def arrow_selection_changed(self, index):
        if index != self.selected_arrow and (self.selected_joint != -1 or self.mesh_selected):
            self.selected_arrow = index
            if index == 0: self.selected_axis_name = 'X'
            elif index == 1: self.selected_axis_name = 'Y'
            elif index == 2: self.selected_axis_name = 'Z'
            else: self.selected_axis_name = 'N/A'
            self.update_joint()
            if self.selected_joint != -1:
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

    @QtCore.pyqtSlot(bool)
    def mesh_selected_slot(self, is_selected):
        if is_selected != self.mesh_selected:
            self.mesh_selected = is_selected
            self.update_joint()

    # ── Drag / transform ──────────────────────────────────────────────────────

    @QtCore.pyqtSlot(np.ndarray)
    def drag_translate(self, new_position):
        propogate = self.propogate_slider_checkbox.isChecked()
        if self.mesh_selected:
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
            old_position = self.model.Joints[self.selected_joint].Pose.t
            trans = new_position - old_position
            transformation = SE3.Trans(trans[0], trans[1], trans[2])
            if self.model.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=False):
                self.update_joint()
            self.reset_rotation_tools()
            self.reset_translation_tools()

    @QtCore.pyqtSlot(float)
    def drag_rotate(self, new_rotation):
        transformation = SE3()
        propogate = self.propogate_slider_checkbox.isChecked()
        if self.selected_axis_name == 'X':
            transformation = SE3.Rx(new_rotation)
        elif self.selected_axis_name == 'Y':
            transformation = SE3.Ry(new_rotation)
        elif self.selected_axis_name == 'Z':
            transformation = SE3.Rz(new_rotation)
        if self.mesh_selected:
            self.referenceMesh.Pose = self.referenceMesh.Pose * transformation
            transform_matrix = QMatrix4x4()
            matrix = self.referenceMesh.Pose.A
            for row in range(4):
                for col in range(4):
                    transform_matrix[row, col] = matrix[row, col]
            self.referenceMesh.mesh.setTransform(transform_matrix)
            self.update_joint()
        else:
            if self.model.transformJoint(self.selected_joint, transformation, propogate=propogate, safe=True, relative=True):
                self.update_joint()
                self.reset_rotation_tools()
                self.reset_translation_tools()

    def rotate_joint(self, angle, axis):
        propogate = self.propogate_slider_checkbox.isChecked()
        localOrient = self.local_orient_slider_checkbox.isChecked()
        transformation = SE3.AngleAxis(angle, [axis[0], axis[1], axis[2]], unit='deg')
        self.model.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=True, safe=True, localOrient=localOrient)

    # ── Rotation / translation tools ──────────────────────────────────────────

    def rotation_angle_from_matrix(self, rotation_matrix, axis):
        rot = R.from_matrix(rotation_matrix)
        euler_angles = rot.as_euler('xyz', degrees=True)
        return euler_angles[axis]

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
        r = self.model.r if self.model else 1
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
        r = self.model.r if self.model else 1
        self.translation_slider.setMinimum(int(-100 * r))
        self.translation_slider.setMaximum(int(100 * r))
        self.translation_slider.setValue(0)
        self.old_trans_val = 0
        self.translation_textbox.setText("0")
        if self.selected_arrow == -1 or self.selected_joint == -1:
            self.translation_slider.setDisabled(True)
            self.translation_textbox.setDisabled(True)
        else:
            self.translation_slider.setDisabled(False)
            self.translation_textbox.setDisabled(False)

    def eventFilter(self, obj, event):
        """Intercept Tab in config text boxes: apply value and advance to next."""
        if event.type() == QEvent.KeyPress and event.key() == Qt.Key_Tab:
            if hasattr(self, 'config_text_boxes') and obj in self.config_text_boxes:
                idx = self.config_text_boxes.index(obj)
                joint_index = self.config_joint_indices[idx]
                self.config_textbox_return(joint_index)
                if idx + 1 < len(self.config_text_boxes):
                    self.config_text_boxes[idx + 1].setFocus()
                    self.config_text_boxes[idx + 1].selectAll()
                return True
        return super().eventFilter(obj, event)

    # ── State slider / textbox ────────────────────────────────────────────────

    def state_slider_moved(self, value):
        if self.model and self.selected_joint != -1:
            if isinstance(self.model.Joints[self.selected_joint], Prismatic):
                actualState = value / (100 * self.model.r)
            elif isinstance(self.model.Joints[self.selected_joint], Revolute):
                actualState = math.radians(value)
            else:
                print("Warning: Tried to move state slider on a waypoint, which should not be possible.")
                return
            if self.model.setJointState(self.selected_joint, actualState, lightweight=True):
                if not self._fast_update_gl():
                    self.model.resyncFromLightweight()
                    self.update_joint()
                self.set_state_tools()
            else:
                self.state_slider.blockSignals(True)
                self.state_slider.setValue(int(self.old_state_slider_val))
                self.state_slider.blockSignals(False)

    def state_textbox_return(self):
        if self.model and self.selected_joint != -1:
            try:
                value = float(self.state_textbox.text())
            except ValueError:
                return
            if isinstance(self.model.Joints[self.selected_joint], Prismatic):
                actualState = value
            elif isinstance(self.model.Joints[self.selected_joint], Revolute):
                actualState = math.radians(value)
            else:
                print("Warning: Tried to edit state textbox on a waypoint, which should not be possible.")
                return
            if self.model.setJointState(self.selected_joint, actualState):
                self._mark_collision_dirty()
                self.update_joint()
                self.set_state_tools()
                self.log_version()

    def state_slider_released(self):
        self._resync_after_lightweight()
        self.set_state_tools()
        self.log_version()

    def scaled_state_info(self, state=None):
        """Returns tuple of tuples (actual, slider, textbox) for current joint state,
        or None if no joint selected / joint is a Waypoint."""
        if self.model and self.selected_joint != -1:
            joint = self.model.Joints[self.selected_joint]
            stateRange = joint.stateRange()
            if isinstance(joint, Prismatic):
                stateActual = joint.state if state is None else max(stateRange[0], min(stateRange[1], state))
                actual = (stateRange[0], stateRange[1], stateActual)
                scale = 100 * self.model.r
                scaled = (int(stateRange[0] * scale), int(stateRange[1] * scale), int(stateActual * scale))
                return (actual, scaled, actual)
            elif isinstance(joint, Revolute):
                if state is None:
                    stateRadians = joint.state
                    stateDegrees = math.degrees(stateRadians)
                else:
                    stateDegrees = state
                    stateRadians = math.radians(stateDegrees)
                radians = (stateRange[0], stateRange[1], stateRadians)
                degrees = (round(math.degrees(stateRange[0])), round(math.degrees(stateRange[1])), round(stateDegrees))
                return (radians, degrees, degrees)
            elif isinstance(joint, Waypoint):
                return None
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
            decimals = 2 if isinstance(self.model.Joints[self.selected_joint], Prismatic) else 0
            self.state_textbox.setText(str(np.round(currentText, decimals)))
            self.current_state_label.setText(
                f"Min State: {np.round(minText, decimals)} ≤ "
                f"Current State: {np.round(currentText, decimals)} ≤ "
                f"Max State: {np.round(maxText, decimals)}"
            )

    # ── Configuration widgets ─────────────────────────────────────────────────

    def update_configurations(self):
        """Rebuild the configuration text-boxes/sliders for all non-Waypoint joints."""
        while self.configurations_layout.count():
            item = self.configurations_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
            elif item.layout():
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

        if self.model is None or not self.model_created:
            return

        self.config_joint_widgets = []
        for actual_joint_index, joint in enumerate(self.model.Joints):
            if self._should_skip_joint_in_config(joint):
                continue

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
            text_box.returnPressed.connect(lambda idx=actual_joint_index: self.config_textbox_return(idx))
            text_box.installEventFilter(self)

            slider = QSlider(Qt.Horizontal)
            slider.setFixedWidth(50)

            if isinstance(joint, Prismatic):
                display_state = joint.state
                state_range = joint.stateRange()
                scale = 100 * self.model.r
                slider.setMinimum(int(state_range[0] * scale))
                slider.setMaximum(int(state_range[1] * scale))
                slider.setValue(int(display_state * scale))
            elif isinstance(joint, Revolute):
                display_state = math.degrees(joint.state)
                state_range = joint.stateRange()
                slider.setMinimum(int(math.degrees(state_range[0])))
                slider.setMaximum(int(math.degrees(state_range[1])))
                slider.setValue(int(display_state))
            else:
                display_state = 0.0

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

        self.config_save_button = QPushButton("Define Configuration")
        self.config_save_button.clicked.connect(self.save_current_configuration)
        self.config_save_button.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        self.saved_configs_container = QHBoxLayout()
        self.configurations_layout.addLayout(self.saved_configs_container)

        self.display_saved_configurations()
        self.highlight_selected_config_box()

    def update_config_values(self):
        """Update config text-box/slider values without recreating widgets."""
        if not hasattr(self, 'config_text_boxes') or not self.config_text_boxes:
            return
        if self.model is None or not self.model_created:
            return
        for i, joint_index in enumerate(self.config_joint_indices):
            if joint_index >= len(self.model.Joints):
                continue
            joint = self.model.Joints[joint_index]
            text_box = self.config_text_boxes[i]
            slider = self.config_sliders[i]
            text_box.blockSignals(True)
            slider.blockSignals(True)
            if isinstance(joint, Prismatic):
                display_state = joint.state
                text_box.setText(str(round(display_state, 2)))
                slider.setValue(int(display_state * 100 * self.model.r))
            elif isinstance(joint, Revolute):
                display_state = math.degrees(joint.state)
                text_box.setText(str(round(display_state, 2)))
                slider.setValue(int(display_state))
            text_box.blockSignals(False)
            slider.blockSignals(False)

    def highlight_selected_config_box(self):
        """Highlight the config text-box for the currently selected joint."""
        if not hasattr(self, 'config_text_boxes') or not self.config_text_boxes:
            return
        for text_box in self.config_text_boxes:
            text_box.setStyleSheet("")
        if self.selected_joint != -1 and self.selected_joint in self.config_joint_indices:
            try:
                text_box_index = self.config_joint_indices.index(self.selected_joint)
                self.config_text_boxes[text_box_index].setStyleSheet(
                    "background-color: #FFD700; border: 2px solid #FFA500;"
                )
            except (ValueError, IndexError):
                pass

    def config_textbox_return(self, joint_index):
        """Apply a typed value in a configuration text-box."""
        if self.model and 0 <= joint_index < len(self.model.Joints):
            try:
                text_box_index = self.config_joint_indices.index(joint_index)
                text_box = self.config_text_boxes[text_box_index]
            except (ValueError, IndexError):
                return
            try:
                value = float(text_box.text())
            except ValueError:
                return
            joint = self.model.Joints[joint_index]
            if isinstance(joint, Prismatic):
                actualState = value
            elif isinstance(joint, Revolute):
                actualState = math.radians(value)
            else:
                return
            if self.model.setJointState(joint_index, actualState):
                self._mark_collision_dirty()
                self.update_joint()
                self.set_state_tools()
                self.log_version()

    def config_slider_moved(self, joint_index, value):
        """Handle configuration slider movement."""
        if self.model and 0 <= joint_index < len(self.model.Joints):
            joint = self.model.Joints[joint_index]
            if isinstance(joint, Prismatic):
                actualState = value / (100 * self.model.r)
                display_state = actualState
            elif isinstance(joint, Revolute):
                actualState = math.radians(value)
                display_state = value
            else:
                return
            if self.model.setJointState(joint_index, actualState, lightweight=True):
                if not self._fast_update_gl():
                    self.model.resyncFromLightweight()
                    self.update_joint()
                try:
                    text_box_index = self.config_joint_indices.index(joint_index)
                    self.config_text_boxes[text_box_index].blockSignals(True)
                    self.config_text_boxes[text_box_index].setText(str(round(display_state, 2)))
                    self.config_text_boxes[text_box_index].blockSignals(False)
                except (ValueError, IndexError):
                    pass
                self.set_state_tools()

    def config_slider_released(self, joint_index):
        """Handle configuration slider release."""
        self._resync_after_lightweight()
        self.set_state_tools()
        self.log_version()

    # ── Saved configurations ──────────────────────────────────────────────────

    def save_current_configuration(self):
        """Save current joint states as a named configuration."""
        if not self.config_text_boxes:
            return
        current_config = []
        for text_box in self.config_text_boxes:
            try:
                current_config.append(float(text_box.text()))
            except ValueError:
                current_config.append(0.0)
        insert_index = len(self.saved_configurations)
        if hasattr(self, 'config_interp_slider') and \
                len(self.saved_configurations) > 1 and \
                not sip.isdeleted(self.config_interp_slider):
            slider_value = self.config_interp_slider.value()
            max_value = self.config_interp_slider.maximum()
            inverted_value = max_value - slider_value
            config_value = inverted_value / 100.0
            nearest_config = round(config_value)
            if abs(config_value - nearest_config) > 0.01:
                insert_index = int(config_value) + 1
                insert_index = min(insert_index, len(self.saved_configurations))
        self.saved_configurations.insert(insert_index, current_config)
        self.display_saved_configurations()
        if hasattr(self, 'config_interp_slider') and \
                not sip.isdeleted(self.config_interp_slider):
            max_value = self.config_interp_slider.maximum()
            slider_value = max_value - (insert_index * 100)
            self.config_interp_slider.blockSignals(True)
            self.config_interp_slider.setValue(slider_value)
            self.config_interp_slider.blockSignals(False)

    def display_saved_configurations(self):
        """Render all saved configurations as a grid below the current config."""
        saved_slider_value = None
        if hasattr(self, 'config_interp_slider') and \
                len(self.saved_configurations) > 1 and \
                not sip.isdeleted(self.config_interp_slider):
            saved_slider_value = self.config_interp_slider.value()

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
            self.saved_configs_container.addStretch()
            combined_widget = QWidget()
            grid = QGridLayout(combined_widget)
            grid.setContentsMargins(0, 0, 0, 0)
            grid.setSpacing(0)
            for i, jw in enumerate(self.config_joint_widgets):
                grid.addWidget(jw, 0, i + 2)
            grid.addWidget(self.config_save_button, 0, num_joints + 2, 1, 2)
            combined_widget.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
            self.saved_configs_container.addWidget(combined_widget)
            self.saved_configs_container.addStretch()
            return

        combined_widget = QWidget()
        grid = QGridLayout(combined_widget)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setSpacing(0)

        for i, jw in enumerate(self.config_joint_widgets):
            grid.addWidget(jw, 0, i + 2)
        grid.addWidget(self.config_save_button, 0, num_joints + 2, 1, 2)

        for config_index, config in enumerate(self.saved_configurations):
            row = config_index + 1
            up_button = QPushButton("▲")
            up_button.setFixedSize(30, 22)
            up_button.setEnabled(config_index > 0)
            up_button.clicked.connect(lambda checked, idx=config_index: self.move_configuration_up(idx))
            grid.addWidget(up_button, row, 0)
            down_button = QPushButton("▼")
            down_button.setFixedSize(30, 22)
            down_button.setEnabled(config_index < len(self.saved_configurations) - 1)
            down_button.clicked.connect(lambda checked, idx=config_index: self.move_configuration_down(idx))
            grid.addWidget(down_button, row, 1)
            for j, value in enumerate(config):
                value_label = QLabel(str(round(value, 2)))
                value_label.setAlignment(Qt.AlignCenter)
                value_label.setFixedWidth(55)
                grid.addWidget(value_label, row, j + 2)
            set_button = QPushButton("Set")
            set_button.setFixedSize(60, 22)
            set_button.clicked.connect(lambda checked, idx=config_index: self.set_configuration(idx))
            grid.addWidget(set_button, row, num_joints + 2)
            delete_button = QPushButton("✗")
            delete_button.setFixedSize(30, 22)
            delete_button.setStyleSheet("background-color: #FF4444; color: white; font-weight: bold;")
            delete_button.clicked.connect(lambda checked, idx=config_index: self.delete_configuration(idx))
            grid.addWidget(delete_button, row, num_joints + 3)

        combined_widget.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
        self.saved_configs_container.addStretch()
        self.saved_configs_container.addWidget(combined_widget)

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
            self.config_interp_slider.setMaximum((len(self.saved_configurations) - 1) * 100)
            if saved_slider_value is not None:
                saved_slider_value = min(saved_slider_value, self.config_interp_slider.maximum())
                self.config_interp_slider.setValue(saved_slider_value)
            else:
                self.config_interp_slider.setValue(0)
            self.config_interp_slider.setTickPosition(QSlider.TicksBothSides)
            self.config_interp_slider.setTickInterval(100)
            self.config_interp_slider.sliderMoved.connect(self.interpolate_configurations)
            self.config_interp_slider.sliderReleased.connect(self.interpolation_slider_released)
            slider_layout.addWidget(self.config_interp_slider, 0, Qt.AlignHCenter)
            self.saved_configs_container.addWidget(slider_widget)

            animation_widget = QWidget()
            animation_layout = QVBoxLayout(animation_widget)
            animation_layout.setContentsMargins(10, 0, 10, 0)
            animation_layout.setAlignment(Qt.AlignTop | Qt.AlignHCenter)
            animation_label = QLabel("Animate (s)")
            animation_label.setAlignment(Qt.AlignCenter)
            animation_label.setMinimumWidth(90)
            animation_layout.addWidget(animation_label, 0, Qt.AlignHCenter)
            play_loop_layout = QHBoxLayout()
            play_loop_layout.setSpacing(5)
            if self.is_animating:
                self.play_pause_button = QPushButton("⏸")
                self.play_pause_button.setToolTip("Pause")
            else:
                self.play_pause_button = QPushButton("▶")
                self.play_pause_button.setToolTip("Play")
            self.play_pause_button.setFixedSize(30, 24)
            self.play_pause_button.clicked.connect(self.toggle_animation)
            play_loop_layout.addWidget(self.play_pause_button)
            self.animation_loop_checkbox = QCheckBox("Loop")
            self.animation_loop_checkbox.setMinimumWidth(70)
            self.animation_loop_checkbox.setChecked(self.animation_loop)
            self.animation_loop_checkbox.stateChanged.connect(self.toggle_animation_loop)
            play_loop_layout.addWidget(self.animation_loop_checkbox)
            animation_layout.addLayout(play_loop_layout)
            num_segments = len(self.saved_configurations) - 1
            while len(self.config_durations) < num_segments:
                self.config_durations.append(1.0)
            while len(self.config_durations) > num_segments:
                self.config_durations.pop()
            for i in range(num_segments):
                if i == 0:
                    animation_layout.addSpacing(0)
                else:
                    animation_layout.addSpacing(2)
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
        """Interpolate between saved configurations based on slider value."""
        if len(self.saved_configurations) < 2:
            return
        max_value = self.config_interp_slider.maximum()
        inverted_value = max_value - slider_value
        config_value = inverted_value / 100.0
        config_index = int(config_value)
        if config_index >= len(self.saved_configurations) - 1:
            config_index = len(self.saved_configurations) - 2
        t = config_value - config_index
        config_a = self.saved_configurations[config_index]
        config_b = self.saved_configurations[config_index + 1]
        for i in range(min(len(config_a), len(config_b))):
            if i >= len(self.config_joint_indices):
                break
            joint_index = self.config_joint_indices[i]
            if joint_index >= len(self.model.Joints):
                continue
            interpolated_value = config_a[i] * (1 - t) + config_b[i] * t
            joint = self.model.Joints[joint_index]
            if isinstance(joint, Prismatic):
                actualState = interpolated_value
            elif isinstance(joint, Revolute):
                actualState = math.radians(interpolated_value)
            else:
                continue
            self.model.setJointState(joint_index, actualState, lightweight=True)
        if not self._fast_update_gl():
            self.model.resyncFromLightweight()
            self.update_joint()
        self.set_state_tools()

    def interpolation_slider_released(self):
        """Rebuild geometry after interpolation slider is released."""
        self._resync_after_lightweight()
        self.log_version()

    def set_configuration(self, config_index):
        """Set the model to a saved configuration."""
        if config_index >= len(self.saved_configurations):
            return
        config = self.saved_configurations[config_index]
        for i, value in enumerate(config):
            if i >= len(self.config_joint_indices):
                break
            joint_index = self.config_joint_indices[i]
            if joint_index >= len(self.model.Joints):
                continue
            joint = self.model.Joints[joint_index]
            if isinstance(joint, Prismatic):
                actualState = value
            elif isinstance(joint, Revolute):
                actualState = math.radians(value)
            else:
                continue
            self.model.setJointState(joint_index, actualState)
        if hasattr(self, 'config_interp_slider') and \
                not sip.isdeleted(self.config_interp_slider):
            max_value = self.config_interp_slider.maximum()
            slider_value = max_value - (config_index * 100)
            self.config_interp_slider.blockSignals(True)
            self.config_interp_slider.setValue(slider_value)
            self.config_interp_slider.blockSignals(False)
        self.update_joint()
        self.set_state_tools()
        self.log_version()

    def delete_configuration(self, config_index):
        if config_index < len(self.saved_configurations):
            self.saved_configurations.pop(config_index)
            self.display_saved_configurations()

    def move_configuration_up(self, config_index):
        if config_index > 0 and config_index < len(self.saved_configurations):
            self.saved_configurations[config_index], self.saved_configurations[config_index - 1] = \
                self.saved_configurations[config_index - 1], self.saved_configurations[config_index]
            self.display_saved_configurations()

    def move_configuration_down(self, config_index):
        if config_index >= 0 and config_index < len(self.saved_configurations) - 1:
            self.saved_configurations[config_index], self.saved_configurations[config_index + 1] = \
                self.saved_configurations[config_index + 1], self.saved_configurations[config_index]
            self.display_saved_configurations()

    def update_saved_configs_for_joint_added(self, config_position: int = -1):
        """Insert state 0 into all saved configs when a joint is added."""
        for config in self.saved_configurations:
            if config_position < 0 or config_position >= len(config):
                config.append(0.0)
            else:
                config.insert(config_position, 0.0)

    def update_saved_configs_for_joint_deleted(self, deleted_joint_index):
        """Remove the state for a deleted joint from all saved configs."""
        if deleted_joint_index not in self.config_joint_indices:
            return
        try:
            config_position = self.config_joint_indices.index(deleted_joint_index)
        except ValueError:
            return
        for config in self.saved_configurations:
            if config_position < len(config):
                config.pop(config_position)

    # ── Animation ─────────────────────────────────────────────────────────────

    def update_duration(self, segment_index):
        try:
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
        if self.is_animating:
            self.pause_animation()
        else:
            self.start_animation()

    def start_animation(self):
        if len(self.saved_configurations) < 2:
            return
        self.is_animating = True
        self._pre_anim_selected_joint = self.selected_joint
        self._pre_anim_selected_arrow = self.selected_arrow
        self.selected_joint = -1
        self.selected_arrow = -1
        self.update_joint()
        if self.animation_timer is None:
            self.animation_timer = qc.QTimer()
            self.animation_timer.timeout.connect(self.animation_step)
        self.animation_start_value = self.config_interp_slider.value()
        self.animation_start_time = qc.QTime.currentTime().msecsSinceStartOfDay() / 1000.0
        max_value = self.config_interp_slider.maximum()
        inverted_value = max_value - self.animation_start_value
        current_config_value = inverted_value / 100.0
        self.animation_current_segment = int(current_config_value)
        self.animation_segment_start_time = self.animation_start_time
        self.play_pause_button.setText("⏸")
        self.play_pause_button.setToolTip("Pause")
        self.animation_timer.start(16)

    def pause_animation(self):
        self.is_animating = False
        if self.animation_timer is not None:
            self.animation_timer.stop()
        if hasattr(self, '_pre_anim_selected_joint'):
            self.selected_joint = self._pre_anim_selected_joint
            self.selected_arrow = self._pre_anim_selected_arrow
        self._resync_after_lightweight()
        self.play_pause_button.setText("▶")
        self.play_pause_button.setToolTip("Play")

    def toggle_animation_loop(self, state):
        self.animation_loop = (state == Qt.Checked)

    def animation_step(self):
        if not self.is_animating or len(self.saved_configurations) < 2:
            return
        current_time = qc.QTime.currentTime().msecsSinceStartOfDay() / 1000.0
        elapsed_in_segment = current_time - self.animation_segment_start_time
        segment_duration = self.config_durations[self.animation_current_segment] \
            if self.animation_current_segment < len(self.config_durations) else 1.0
        segment_progress = elapsed_in_segment / segment_duration if segment_duration > 0 else 1.0
        if segment_progress >= 1.0:
            self.animation_current_segment += 1
            if self.animation_current_segment >= len(self.saved_configurations) - 1:
                if self.animation_loop:
                    self.animation_current_segment = 0
                    self.animation_segment_start_time = current_time
                    max_value = self.config_interp_slider.maximum()
                    self.config_interp_slider.blockSignals(True)
                    self.config_interp_slider.setValue(max_value)
                    self.config_interp_slider.blockSignals(False)
                    self.interpolate_configurations(max_value)
                    return
                else:
                    self.pause_animation()
                    return
            else:
                self.animation_segment_start_time = current_time
                segment_progress = 0.0
        new_config_value = self.animation_current_segment + segment_progress
        max_value = self.config_interp_slider.maximum()
        new_slider_value = max_value - new_config_value * 100
        self.config_interp_slider.blockSignals(True)
        self.config_interp_slider.setValue(int(new_slider_value))
        self.config_interp_slider.blockSignals(False)
        self.interpolate_configurations(int(new_slider_value))

    # ── Version history (undo / redo / autosave) ──────────────────────────────

    def log_version(self):
        self.versions = self.versions[:self.version_index + 1]
        if self.total_version_counter % self.AUTOSAVE_FREQUENCY == 0 and self.model is not None:
            self._save_model(autosave_id=time.time())
        if len(self.versions) < self.LOG_CAPACITY:
            self.versions.append(copy.deepcopy(self.model))
        else:
            self.versions.pop(0)
            self.versions.append(copy.deepcopy(self.model))
        self.version_index = len(self.versions) - 1
        self.total_version_counter += 1

    def undo(self):
        prev_model = self.model
        if self.version_index > 0:
            self.version_index -= 1
            self.model = copy.deepcopy(self.versions[self.version_index])
        else:
            if self.version_index == 0:
                self.version_index = -1
            self.model = None
        self._mark_collision_dirty()
        self._post_undo_hook(prev_model)
        self.update_joint()

    def redo(self):
        if self.version_index + 1 < len(self.versions):
            prev_model = self.model
            self.version_index += 1
            self.model = copy.deepcopy(self.versions[self.version_index])
            self._mark_collision_dirty()
            self._post_undo_hook(prev_model)
            self.update_joint()

    # ── Fast / lightweight GL update ──────────────────────────────────────────

    def _fast_update_gl(self) -> bool:
        """Fast path: update GL transforms without full redraw.
        Returns True if the fast path succeeded."""
        if self.model is None:
            return False
        if not self.model.updateGLTransforms():
            return False
        self._on_fast_update_gl_extra()
        self._lightweight_dirty = True
        self.plot_widget.update()
        return True

    def _resync_after_lightweight(self):
        """Rebuild link geometry after lightweight animation, then do a full update_joint."""
        if self.model is not None:
            self.model.resyncFromLightweight()
        self._lightweight_dirty = False
        self._mark_collision_dirty()
        self.update_joint()

    # ── Measurement overlay ───────────────────────────────────────────────────

    def _toggle_measure_select(self, enabled: bool):
        self.measure_select_mode = enabled
        if enabled:
            self.selected_joint = -1
        self.update_joint()

    def _set_measure_active_endpoint(self, ep: int):
        self.measure_active_endpoint = ep
        self.measurement_widget.set_active_endpoint(ep)

    def _on_measurement_dock_visibility(self, visible: bool):
        if not visible and self.measure_select_mode:
            self.measure_select_mode = False
            self.measurement_widget.deactivate()

    def _clear_measurement(self):
        self.measure_points = [None, None]
        self.measure_select_mode = False
        self.measure_active_endpoint = 0
        self.measurement_widget.reset()
        self.update_joint()

    def _measure_joint_selected_slot(self, joint_index: int):
        if joint_index == -1 or self.model is None:
            return
        ep = self.measure_active_endpoint
        frame_type = self.measurement_widget.frame_types[ep]
        self.measure_points[ep] = (joint_index, frame_type)
        pos = self._get_measure_position(joint_index, frame_type)
        joint = self.model.Joints[joint_index]
        self.measurement_widget.set_point_info(ep, joint_index, type(joint).__name__, pos)
        next_ep = 1 - ep
        self.measure_active_endpoint = next_ep
        self.measurement_widget.set_active_endpoint(next_ep)
        self.update_joint()

    def _measure_frame_type_changed(self):
        if self.model is None:
            return
        for ep in (0, 1):
            if self.measure_points[ep] is not None:
                joint_idx = self.measure_points[ep][0]
                new_frame = self.measurement_widget.frame_types[ep]
                self.measure_points[ep] = (joint_idx, new_frame)
                pos = self._get_measure_position(joint_idx, new_frame)
                joint = self.model.Joints[joint_idx]
                self.measurement_widget.set_point_info(ep, joint_idx, type(joint).__name__, pos)
        self.update_joint()

    def _get_measure_position(self, joint_index: int, frame_type: str):
        if self.model is None:
            return None
        joints = self.model.Joints
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
                    R_mat = self.model.Joints[self.measure_points[0][0]].Pose.R
                    components = R_mat.T @ displacement
                elif ax == 'local_b':
                    R_mat = self.model.Joints[self.measure_points[1][0]].Pose.R
                    components = R_mat.T @ displacement
                else:
                    components = displacement
                for ep, mp in enumerate(self.measure_points):
                    pos = self._get_measure_position(*mp)
                    joint = self.model.Joints[mp[0]]
                    self.measurement_widget.set_point_info(ep, mp[0], type(joint).__name__, pos)
                self.measurement_widget.set_result(dist, components)
            else:
                self.measurement_widget.set_result(None, None)
        else:
            self.measurement_widget.set_result(None, None)
