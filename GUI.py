"""
@author: Raymond Feng, Andy Wang, Daniel Feshbach
"""

import sys, os

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
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton, QDockWidget, QComboBox, QHBoxLayout, QLabel, QDialog, QLineEdit, QCheckBox, QMessageBox, QButtonGroup, QRadioButton, QSlider, QSizePolicy, QFileDialog
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QPixmap, QSurfaceFormat, QKeyEvent, QPixmap, QIcon, QMatrix4x4, QVector3D, QMatrix3x3
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
from style import *
from ReferenceMesh import *
from Dialog import *
from JointWidget import *
from IntersectionHelper import *

import warnings
warnings.filterwarnings("ignore")

if hasattr(QtCore.Qt, 'AA_ENnableHighDpiScaling'):
    PyQt5.QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)

if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
    PyQt5.QtWidgets.QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, True)

class DeleteWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Confirm Delete')

        layout = QVBoxLayout()
        layout.addWidget(QLabel('Are you sure you want to delete the joint?'))

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
            self, "Import STL", os.path.join(base_path, "referenceMeshes/"), "STL Files (*.stl)", options=options
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

class AddChainWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Create New Chain')

        layout = QVBoxLayout()

        # Input for the number of sides
        numSides_layout = QHBoxLayout()
        num_sides_label = QLabel("Number of Sides:")
        self.num_sides_input = QLineEdit()
        self.num_sides_input.setPlaceholderText("Enter number of sides")
        numSides_layout.addWidget(num_sides_label)
        numSides_layout.addWidget(self.num_sides_input)
        layout.addLayout(numSides_layout)

        # Input for the radius
        radius_layout = QHBoxLayout()
        radius_label = QLabel("Radius:")
        self.radius_input = QLineEdit()
        self.radius_input.setPlaceholderText("Enter radius")
        radius_layout.addWidget(radius_label)
        radius_layout.addWidget(self.radius_input)
        layout.addLayout(radius_layout)

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
            self.window().joint_range_slider.setEnabled(True)
        except ValueError:
            self.show_error("Please enter valid integers.")

    def on_cancel_clicked(self):
        # Hide the widget if the user cancels
        self.window().add_chain_popup_dock.setVisible(False)
        self.window().add_chain_dock.setVisible(True)

    def show_error(self, message):
        QMessageBox.warning(self, "Invalid Input", message)

class EditGridWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Edit Grid')

        layout = QVBoxLayout()

        # Input for line spacing
        spacing_layout = QHBoxLayout()
        spacing_label = QLabel("Grid line spacing: ")
        self.unit_label = QLabel("(cm)")
        self.spacing_input = QLineEdit()
        self.spacing_input.setPlaceholderText("Enter spacing")
        spacing_layout.addWidget(spacing_label)
        spacing_layout.addWidget(self.spacing_input)
        spacing_layout.addWidget(self.unit_label)
        layout.addLayout(spacing_layout)

        # Input for the amount of lines
        line_amt_layout = QHBoxLayout()
        line_amt_label = QLabel("Size:")
        self.line_amt_input = QLineEdit()
        grid_lines_label = QLabel("grid lines")
        self.line_amt_input.setPlaceholderText("Enter line amount")
        line_amt_layout.addWidget(line_amt_label)
        line_amt_layout.addWidget(self.line_amt_input)
        line_amt_layout.addWidget(grid_lines_label)
        layout.addLayout(line_amt_layout)

        button_layout = QHBoxLayout()

        # Apply button to create the chain
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
            self.window().update_joint()
        except ValueError:
            self.show_error("Please enter valid integers.")

    def on_cancel_clicked(self):
        # Hide the widget if the user cancels
        self.window().edit_grid_dock.setVisible(False)

    def show_error(self, message):
        QMessageBox.warning(self, "Invalid Input", message)

class EditDimensionsWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Edit Dimensions')

        self.layout = QVBoxLayout()

        self.radio_layout = QHBoxLayout()

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
    
    def get_normalized_plane_vectors(self, event):
        if self.is_dragging and self.selected_torus:
            origin, dir = self.get_world_coordinates(event)
            selected_joint = self.parent_window.chain.Joints[self.parent_window.selected_joint]
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
            selected_joint = self.parent_window.chain.Joints[self.parent_window.selected_joint]
            center = selected_joint.Pose.t
            qcenter = QVector3D(center[0], center[1], center[2])
            qaxis = self.selected_torus
            
            plane_vector, normal = self.get_normalized_plane_vectors(event)

            prev_vector = self.drag_prev_vector

            d = QVector3D.dotProduct(plane_vector, prev_vector)
            angle = math.acos(d / (plane_vector.length() * prev_vector.length()))

            cross_product = QVector3D.crossProduct(prev_vector, plane_vector)
            
            if QVector3D.dotProduct(cross_product, normal) < 0:
                angle *= -1

            self.drag_prev_vector = plane_vector

            return math.degrees(angle), normal
    
    def get_closest_point(self, event):
        if self.is_dragging and self.selected_axis:
            origin, dir = self.get_world_coordinates(event)
            selected_joint = self.parent_window.chain.Joints[self.parent_window.selected_joint]
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
        shortest_location = 1000
        shortest_idx = -1
        hit_cylinder = False

        for i in range(3):
            if (self.is_local):
                hit_location = compute_cylinder_intersection(origin, direction, joint_center, self.selected_joint_axes[i], 0.2, selected_joint.boundingBall().r+1)
            else:
                hit_location = compute_cylinder_intersection(origin, direction, joint_center, self.original_axes[i], 0.2, selected_joint.boundingBall().r+1)
            
            if (hit_location < shortest_location):
                shortest_location = hit_location
                shortest_idx = i
                hit_cylinder = True

        if (hit_cylinder):
            self.click_signal_arrow.emit(shortest_idx)
            
            if (self.is_local):
                self.selected_axis = self.selected_joint_axes[shortest_idx]
            else:
                self.selected_axis = self.original_axes[shortest_idx]

            self.hit_cylinder = True
            self.is_dragging = True
        else:
            self.is_dragging = False
            self.selected_axis = None

    def calculate_rotation_isect(self, origin, direction, joint_center, selected_joint):
        shortest_location = 1000
        shortest_idx = -1
        hit_torus = False

        for i in range(3):
            center = joint_center
            tor_rad = selected_joint.r

            if (self.is_local):
                cyl_axis = self.selected_joint_axes[i]
            else:
                cyl_axis = self.original_axes[i]

            hit_location = compute_torus_intersection(origin, direction, center, cyl_axis, major_radius=tor_rad, minor_radius=0.4)
            if (hit_location < shortest_location):
                shortest_location = hit_location
                shortest_idx = i
                hit_torus = True

        if (hit_torus):
            if (self.is_local):
                self.selected_torus = self.selected_joint_axes[shortest_idx]
            else:
                self.selected_torus = self.original_axes[shortest_idx]

            self.selected_axis_orig = self.original_axes[shortest_idx]

            dot = QVector3D.dotProduct(direction, self.selected_torus)
            self.facing_same_dir = dot > 0
                    
            self.click_signal_arrow.emit(shortest_idx)
        else:
            self.is_dragging = False
            self.selected_torus = None

    def mousePressEvent(self, event):
        self.last_drag_pos = event.pos()

        if (event.buttons() and event.button() == Qt.MouseButton.MiddleButton):
            self.drag_start_pos = event.pos()

        if (event.buttons() and Qt.LeftButton and event.buttons() != QtCore.Qt.MouseButton.MiddleButton):
            self.is_dragging = False
            self.drag_start_pos = event.pos()

            origin, direction = self.get_world_coordinates(event)

            self.is_local = self.parent_window.is_local

            # raycasting for widgets
            if (self.parent_window.selected_joint != -1):
                self.selected_axis = None
                self.selected_torus = None
                selected_joint = self.parent_window.chain.Joints[self.parent_window.selected_joint]
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
            closest_joint_dist = 1000
            closest_joint = None

            if (self.parent_window.chain):
                for joint in self.parent_window.chain.Joints:
                    center = joint.Pose.t
                    radius = self.parent_window.radius
                    hit_location = compute_sphere_intersection(origin, direction, center, radius + .2)
                    if (hit_location < closest_joint_dist):
                        closest_joint_dist = hit_location
                        closest_joint = joint
            
            self.selected_joint_temp = closest_joint

    def mouseMoveEvent(self, event):
        if (event.buttons() and (Qt.LeftButton or Qt.MiddleButton)) and (event.pos() - self.drag_start_pos).manhattanLength() >= QApplication.startDragDistance():
            self.is_dragging = True

        if (self.is_dragging):
            if (self.selected_axis):
                new_pos_3D = self.get_closest_point(event)

                selected_joint = self.parent_window.chain.Joints[self.parent_window.selected_joint]
                joint_center = selected_joint.Pose.t
                qsphere_start = QVector3D(joint_center[0], joint_center[1], joint_center[2])

                trans = new_pos_3D - qsphere_start
                transformation = SE3.Trans(trans[0], trans[1], trans[2])

                propogate = self.parent_window.propogate_slider_checkbox.isChecked()

                self.parent_window.chain.transformJoint(self.parent_window.selected_joint, transformation, propogate=propogate, relative=False)
                self.parent_window.update_joint()
            elif (self.selected_torus):
                da, normal = self.get_axis_angle_delta(event)
                if (not self.facing_same_dir):
                    da = -da
                self.parent_window.rotate_joint(da, self.selected_axis_orig)
                self.parent_window.update_joint()
            else:
                curr_pos = event.position() if hasattr(event, 'position') else event.localPos()

                diff = curr_pos - self.last_drag_pos
                self.last_drag_pos = curr_pos

                if event.buttons() == QtCore.Qt.MouseButton.MiddleButton:
                    self.pan(diff.x(), diff.y(), 0, relative='view')
                elif event.buttons() == QtCore.Qt.MouseButton.LeftButton:
                    if (self.camera_type == "Rotate"):
                        self.orbit(-diff.x()*self.orbit_speed, diff.y()*self.orbit_speed)
                    elif (self.camera_type == "Pan"):
                        self.pan(diff.x(), diff.y(), 0, relative='view')

    def mouseReleaseEvent(self, event):
        if (self.selected_joint_temp != None):
            self.click_signal.emit(self.parent_window.chain.Joints.index(self.selected_joint_temp))
        elif (not self.is_dragging): 
            self.click_signal.emit(-1)
            self.selected_axis = None
            self.selected_torus = None

        self.is_dragging = False
        self.drag_prev_vector = None
        
        if self.is_dragging and (self.selected_axis or self.selected_torus):
            self.done_transforming.emit(True)

        else:
            # check to see if link or mesh is selected 
            lpos = event.position() if hasattr(event, 'position') else event.localPos()
            region = [lpos.x()-5, lpos.y()-5, 10, 10]
            dpr = self.devicePixelRatioF()
            region = tuple([x * dpr for x in region])

            links = []
            mesh = []

            for item in self.itemsAt(region):
                if (item.objectName() == "Link"):
                    links.append(item)

                if (item.objectName() == "Mesh"):
                    mesh.append(item)

            if (len(mesh) == 0):
                self.mesh_selected = False
            else:
                self.mesh_selected = True
            
            self.click_signal_link.emit(self.selected_link_index)
            self.click_signal_mesh.emit(self.mesh_selected)

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
 
class PointEditorWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Point Editor")
        self.setGeometry(100, 100, 800, 600)

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
        self.version_index = 0
        self.chain_created = False
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

        self.num_sides = 4
        self.radius = 1.0

        self.plot_widget.click_signal.connect(self.joint_selection_changed)
        self.plot_widget.click_signal_arrow.connect(self.arrow_selection_changed)
        self.plot_widget.click_signal_link.connect(self.link_selection_changed)
        self.plot_widget.click_signal_mesh.connect(self.mesh_selected_slot)
        self.plot_widget.drag_change_position.connect(self.drag_translate)
        self.plot_widget.drag_change_rotation.connect(self.drag_rotate)
        self.plot_widget.done_transforming.connect(self.done_transforming)
        self.plot_widget.key_pressed.connect(self.key_pressed)
        
        # //////////////////////////////////    Keyboard Options    ///////////////////////////////////
        top_dock_widget = QDockWidget("Keyboard Controls", self)
        top_dock_widget.setAllowedAreas(Qt.TopDockWidgetArea)

        self.key_bar = QWidget()
        self.key_bar_layout = QHBoxLayout(self.key_bar)  # Layout is initialized and set to the widget here
        self.key_bar.setFixedHeight(40)
        self.init_key_bar()

        top_dock_widget.setWidget(self.key_bar)
        
        self.add_mesh_widget = AddMeshWidget(self)
        self.add_mesh_dock = QDockWidget("Import Mesh", self)
        self.add_mesh_dock.setWidget(self.add_mesh_widget)
        self.add_mesh_dock.setVisible(True)
        self.mesh_scale = 1.0
        self.add_mesh_widget.change_scale.connect(self.change_mesh_scale)

        # //////////////////////////////////    MESSAGE DISPLAY    ///////////////////////////////////
        message_display_widget = QDockWidget("Messages", self)
        #message_display_widget.setAllowedAreas(Qt.BottomDockWidgetArea)

        # Create a layout specifically for success/error messages
        #self.message_layout = QHBoxLayout()

        # Create a label to display success/error messages
        self.status_label = QLabel('')
        message_display_widget.setWidget(self.status_label)

        # Add the message layout to the main layout
        #message_display_widget.setLayout(self.message_layout)

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

        radius_layout = QHBoxLayout()
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
        add_joints_layout.addWidget(self.add_prismatic)
        add_joints_layout.addWidget(self.add_prismatic_menu)
        add_joints_layout.addWidget(self.add_revolute)
        add_joints_layout.addWidget(self.add_revolute_menu)
        add_joints_layout.addLayout(add_waypoints_layout)
        add_joints_layout.addWidget(self.add_tip)
        add_joints_layout.addWidget(self.add_tip_menu)
        add_joints_layout.addWidget(self.edit_dimension_button)
        add_joints_layout.addWidget(self.edit_dimension_menu)

        add_chain_layout.addWidget(self.create_new_chain)

        self.add_prismatic.clicked.connect(self.add_prismatic_toggle)
        self.add_revolute.clicked.connect(self.add_revolute_toggle)
        self.add_waypoint.clicked.connect(self.add_waypoint_func)
        self.add_tip.clicked.connect(self.add_tip_toggle)
        self.create_new_chain.clicked.connect(self.create_new_chain_func)
        self.edit_dimension_button.clicked.connect(self.edit_joint_dimension)

        self.add_chain_dock = QDockWidget("New Chain", self)
        self.add_chain_button_widget = QWidget()
        self.add_chain_button_widget.setLayout(add_chain_layout)
        self.add_chain_dock.setWidget(self.add_chain_button_widget)

        add_joints_dock = QDockWidget("Add Joints", self)
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

        self.joint_add_root_button = QCheckBox("Add to Root")
        self.joint_add_root_button.setChecked(False)
        self.joint_add_root_button.stateChanged.connect(self.add_to_root_func)
        self.select_joint_options = QComboBox()
        self.select_link_options = QComboBox()
        self.delete_joint_button = QPushButton("Delete Joint")
        self.current_state_label = QLabel('Min State ≤ Current State ≤ Max State')

        #joint_layout = QVBoxLayout()
        self.joint_editing_layout.addWidget(self.joint_add_root_button)
        self.joint_editing_layout.addWidget(self.select_joint_options)
        self.joint_editing_layout.addWidget(self.delete_joint_button)
        self.joint_editing_layout.addWidget(self.select_link_options)
        edit_joints_dock = QDockWidget("Edit Joints", self)
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
        self.rotation_slider.valueChanged.connect(self.adjust_rotation)

        # self. = QLabel('Translate N/A Axis: 0', self)
        self.translate_slider = QSlider(Qt.Horizontal, self)
        self.translate_slider.setMinimum(-100)
        self.translate_slider.setMaximum(100)
        self.translate_slider.setValue(0)
        self.translate_slider.valueChanged.connect(self.adjust_translation)

        # self.state_label = QLabel('Edit Joint N/A State: 0', self)
        self.state_slider = QSlider(Qt.Horizontal, self)
        self.state_slider.setMinimum(-100)
        self.state_slider.setMaximum(100)
        self.state_slider.setValue(0)
        self.state_slider.valueChanged.connect(self.adjust_state)

        translation_layout = QVBoxLayout()
        translation_header_layout = QHBoxLayout()
        translation_header_layout.addWidget(self.translate_joint_radio_button)
        translation_header_layout.addWidget(self.translate_joint_label)
        translation_layout.addLayout(translation_header_layout)
        translation_slider_layout = QHBoxLayout()
        self.translation_input = QLineEdit(self)
        self.translation_input.setPlaceholderText("Enter distance")
        self.translation_input.textChanged.connect(self.adjust_translation)
        translation_slider_layout.addWidget(self.translate_slider)
        translation_slider_layout.addWidget(self.translation_input)
        translation_layout.addLayout(translation_slider_layout)

        rotation_layout = QVBoxLayout()
        rotation_header_layout = QHBoxLayout()
        rotation_header_layout.addWidget(self.rotate_joint_radio_button)
        rotation_header_layout.addWidget(self.rotate_joint_label)
        rotation_layout.addLayout(rotation_header_layout)
        rotation_slider_layout = QHBoxLayout()
        self.rotation_input = QLineEdit(self)
        self.rotation_input.setPlaceholderText("Enter angle in degrees")
        self.rotation_input.textChanged.connect(self.adjust_rotation)
        rotation_slider_layout.addWidget(self.rotation_slider)
        rotation_slider_layout.addWidget(self.rotation_input)
        rotation_layout.addLayout(rotation_slider_layout)

        state_layout = QVBoxLayout()
        state_layout.addWidget(self.current_state_label)
        state_slider_layout = QHBoxLayout()
        self.state_input = QLineEdit(self)
        self.state_input.setPlaceholderText("Enter state")
        self.state_input.textChanged.connect(self.adjust_state)
        state_slider_layout.addWidget(self.state_slider)
        state_slider_layout.addWidget(self.state_input)
        state_layout.addLayout(state_slider_layout)

        checkbox_layout = QHBoxLayout() 
        self.propogate_slider_checkbox = QCheckBox("Propagate")
        self.local_orient_slider_checkbox = QCheckBox("Local Orientation")
        self.propogate_slider_checkbox.setChecked(True)
        self.local_orient_slider_checkbox.setChecked(True)
        self.local_orient_slider_checkbox.stateChanged.connect(self.local_orient_clicked)
        checkbox_layout.addWidget(self.propogate_slider_checkbox)
        checkbox_layout.addWidget(self.local_orient_slider_checkbox)

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

        self.joint_editing_layout.addLayout(state_layout)
        self.joint_editing_layout.addLayout(checkbox_layout)
        self.joint_editing_layout.addLayout(rotation_layout)
        self.joint_editing_layout.addLayout(translation_layout)
        self.joint_editing_layout.addLayout(radius_layout)
        self.joint_editing_layout.addLayout(joint_range_layout)

        self.old_rot_val = 0
        self.old_trans_val = 0
        self.old_state_val = 0
        self.old_radius_val = 1
        self.rotation_slider.setDisabled(True)
        self.translate_slider.setDisabled(True)
        self.state_slider.setDisabled(True)
        self.rotation_input.setDisabled(True)
        self.translation_input.setDisabled(True)
        self.state_input.setDisabled(True)

        # ////////////////////////////////    OPTIONS    ///////////////////////////////////
        self.options_dock = QDockWidget("Options", self)
        self.options_widget = QWidget()
        self.options_layout = QVBoxLayout()

        #self.debug_btn = QPushButton("Debug")
        #self.debug_btn.clicked.connect(self.debug)
        #self.options_layout.addWidget(self.debug_btn)

        self.undo_button = QPushButton("Undo")
        self.undo_button.clicked.connect(self.undo)
        self.options_layout.addWidget(self.undo_button)

        self.redo_button = QPushButton("Redo")
        self.redo_button.clicked.connect(self.redo)
        self.options_layout.addWidget(self.redo_button)
        
        self.toggle_grid = QPushButton("Hide Grid")
        self.toggle_grid.clicked.connect(self.toggle_grid_func)
        self.options_layout.addWidget(self.toggle_grid)

        self.options_widget.setLayout(self.options_layout)
        self.options_dock.setWidget(self.options_widget)
        #self.options_dock.setMaximumSize(300, 150)

        # ////////////////////////////////    CAMERA CONTROLS DOCK    ///////////////////////////////////
        self.camera_controls_dock = QDockWidget("Camera Controls", self)
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
        file_dock = QDockWidget("File", self)
        #file_dock.setAllowedAreas(Qt.RightDockWidgetArea)

        file_dock_widget = QWidget()
        file_dock_layout = QVBoxLayout(file_dock_widget)

        self.save_chain_button = QPushButton('Save Chain')
        self.save_chain_button.clicked.connect(self.save_chain)
        file_dock_layout.addWidget(self.save_chain_button)

        self.load_chain_button = QPushButton('Load Chain')
        self.load_chain_button.clicked.connect(self.load_chain)
        file_dock_layout.addWidget(self.load_chain_button)

        self.save_crease_pattern_button = QPushButton('Export Crease Pattern')
        self.save_crease_pattern_button.clicked.connect(self.save_crease_pattern)  
        file_dock_layout.addWidget(self.save_crease_pattern_button) 

        self.units_layout = QHBoxLayout()
        self.units_label = QLabel(f"Current units: {self.units}")
        file_dock_layout.addWidget(self.units_label)

        self.edit_dims_button = QPushButton("Edit Dimensions")
        self.edit_dims_button.clicked.connect(self.edit_dims_func)
        file_dock_layout.addWidget(self.edit_dims_button)

        file_dock_widget.setLayout(file_dock_layout)
        file_dock.setWidget(file_dock_widget)

        # ////////////////////////////////    STL CONVERSION    ///////////////////////////////////
        # self.random_btn_dock = QDockWidget("Export Options", self)
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
        self.delete_joint_dock = QDockWidget("Confirm Delete", self)
        self.delete_joint_dock.setWidget(self.delete_joint_widget)
        self.delete_joint_dock.setVisible(False)
        
        self.add_chain_button_widget = AddChainWidget(self)
        self.add_chain_popup_dock = QDockWidget("Create New Chain", self)
        self.add_chain_popup_dock.setWidget(self.add_chain_button_widget)
        self.add_chain_popup_dock.setVisible(False)  # Initially hidden

        self.edit_grid_widget = EditGridWidget(self)
        self.edit_grid_dock = QDockWidget("Edit Grid", self)
        self.edit_grid_dock.setWidget(self.edit_grid_widget)
        self.edit_grid_dock.setVisible(False)

        self.edit_dims_widget = EditDimensionsWidget(self)
        self.edit_dims_dock = QDockWidget("Edit Dimensions", self)
        self.edit_dims_dock.setWidget(self.edit_dims_widget)
        self.edit_dims_dock.setVisible(False)
        self.edit_dims_widget.target_units.connect(self.change_units)
        
        self.addDockWidget(Qt.TopDockWidgetArea, top_dock_widget)
        self.addDockWidget(Qt.TopDockWidgetArea, message_display_widget)

        self.addDockWidget(Qt.LeftDockWidgetArea, file_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.edit_dims_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.options_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.add_mesh_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.add_chain_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.add_chain_popup_dock)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.camera_controls_dock)
        
        self.addDockWidget(Qt.RightDockWidgetArea, self.add_chain_dock)
        self.addDockWidget(Qt.RightDockWidgetArea, self.add_chain_popup_dock)
        self.addDockWidget(Qt.RightDockWidgetArea, add_joints_dock)
        self.addDockWidget(Qt.RightDockWidgetArea, edit_joints_dock)
        self.addDockWidget(Qt.RightDockWidgetArea, self.delete_joint_dock)

        self.factor = {
            'Centimeter (cm)': 1.0,
            'Meter (m)': 100.0,
            'Milimeter (mm)': 10.0,
            'Inch (in)': 2.54,
            'Feet (ft)': 30.48
        }

    @QtCore.pyqtSlot(str)
    def change_units(self, key):
        self.units = key
        self.units_label.setText(f"Current units: {self.units}")
        print(key)

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
                    new_chain = KinematicChain(self.chain.Joints[0])
                else:
                    new_chain.append(joint, relative=False, fixedPosition=True,
                                            fixedOrientation=True, safe=False)
                
            self.chain = new_chain

            self.plot_widget.opts['distance'] *= factor

            self.update_joint()

    def edit_dims_func(self):
        visibility = self.edit_dims_dock.isVisible()
        self.edit_dims_dock.setVisible(not visibility)

    def initialize_grid(self):
        self.grid = gl.GLGridItem()
        self.grid.setColor(self.grid_color)
        self.grid.setSize(self.grid_size, self.grid_size, self.grid_size)
        self.grid.setSpacing(self.grid_spacing, self.grid_spacing, self.grid_spacing)
        # self.grid.setSpacing(self.grid_spacing)

    def add_to_root_func(self, state):
        self.add_to_root = state == Qt.Checked

    def onUpdateJointState(self, value):
        # regenerate joint
        self.update_joint()

    def onUpdateRadius(self, value):
        value = value / 10.0
        self.radius = value
        self.chain.changeRadius(value)
        self.update_joint()

    @QtCore.pyqtSlot(float)
    def change_mesh_scale(self, scale):
        self.mesh_scale = scale

    def log_version(self):
        #print("logging version")
        log_capacity = 100
        autosave_frequency = 10

        # clear redo history on new version (include version index)
        self.versions = self.versions[:self.version_index + 1]

        if len(self.versions) % autosave_frequency == 0 and not self.chain is None:
            self.save_chain(autosave_id=len(self.versions)//autosave_frequency)
        if len(self.versions) < log_capacity:
            
            self.versions.append(copy.deepcopy(self.chain))
        else:
            self.versions.pop(0)
            self.versions.append(copy.deepcopy(self.chain))

        self.version_index = len(self.versions) - 1

    def undo(self):
        #print("UNDO, current log length: " + str(len(self.versions)))
        if self.version_index > 0:
            self.version_index -= 1
            self.chain = self.versions[self.version_index]
        else:
            self.chain = None
        # self.reload_IDs()
        self.update_joint()

    def redo(self):
        if self.version_index + 1 < len(self.versions):
            self.version_index += 1
            self.chain = self.versions[self.version_index]
            self.update_joint()
    
    def toggle_grid_func(self):
        if self.grid_on:
            self.plot_widget.removeItem(self.grid)
            self.toggle_grid.setText("Show Grid")
        else:
            self.plot_widget.addItem(self.grid)
            self.toggle_grid.setText("Hide Grid")
        self.grid_on = not self.grid_on

    # Success message method with timer
    def show_success(self, message):
        self.status_label.setText(message)
        self.status_label.setStyleSheet("color: " + successColorDefault)

    # Error message method with timer
    def show_error(self, message):
        self.status_label.setText(message)
        self.status_label.setStyleSheet("color: " + errorColorDefault)

    # Method to clear the message
    def clear_message(self):
        self.status_label.setText('')

    def show_delete_widget(self):
        self.delete_joint_dock.setVisible(True)

    def delete_selected_joint(self):
        self.selected_joint = self.select_joint_options.currentIndex() 
        temp_last = len(self.chain.Joints) - 1
        if len(self.chain.Joints) == 1:
            self.chain = None
            self.chain_created = False
            self.plot_widget.clear()
            self.setCentralWidget(self.plot_widget)
            self.show_success('Joint successfully deleted!')
            self.last_joint = -1
        else:
            backup = copy.deepcopy(self.chain)
            try:
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
        # Create a spacer that will expand to push the content to the center
        left_spacer = QWidget()
        left_spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.key_bar_layout.addWidget(left_spacer)

        # Add instructions centered
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

        # Create a spacer that will expand to push the content to the center from the right side
        right_spacer = QWidget()
        right_spacer.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Preferred)
        self.key_bar_layout.addWidget(right_spacer)

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
            if self.chain and self.selected_joint != -1:
                self.delete_selected_joint()
        elif key == "X":
            self.arrow_selection_changed(0)
        elif key == "Y":
            self.arrow_selection_changed(1)
        elif key == "Z":
            self.arrow_selection_changed(2)
        elif key == "G":
            self.toggle_grid_func()

        self.update_joint()

    def save_crease_pattern(self):
        if self.chain:
            options = QFileDialog.Options()
            try:
                base_path = sys._MEIPASS 
            except AttributeError:
                base_path = os.path.abspath(".")
            file_path, _ = QFileDialog.getSaveFileName(
                self, "Save File", os.path.join(base_path, "save/"), "DXF Files (*.dxf)", options=options
            )
            crease_pattern = self.chain.creasePattern()
            if file_path:
                crease_pattern.show(dxfName=file_path)
            else:
                crease_pattern.show()

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
                    self, "Save File", os.path.join(base_path, "save/"), "Chain Files (*.chain)", options=options
                )
            else:
                try:
                    base_path = sys._MEIPASS
                except AttributeError:
                    base_path = os.path.abspath(".")
                file_path = os.path.join(base_path, f"save/autosave/autosave_{autosave_id}.chain")
                # file_path = f"save/autosave/autosave_{autosave_id}.chain"

            if file_path:
                self.chain.save(file_path)
        
    def load_chain(self):
        options = QFileDialog.Options()
        try:
            base_path = sys._MEIPASS
        except AttributeError:
            base_path = os.path.abspath(".")
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open File", os.path.join(base_path, "save/"), "Chain Files (*.chain)", options=options
        )
        if file_path:
            self.chain = loadKinematicChain(file_path)
            self.radius = self.chain.r
            self.num_sides = self.chain.numSides
            self.chain_created = True
            self.update_joint()
            self.log_version()

    @QtCore.pyqtSlot(bool)
    def mesh_selected_slot(self, is_selected):
        self.mesh_selected = is_selected
        self.update_joint()

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
            #self.update_radius_slider()

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
                self.update_rotation_slider()
                self.update_translate_slider()

    @QtCore.pyqtSlot(int)
    def link_selection_changed(self, index):
        self.select_link_options.setCurrentIndex(self.selected_link)
        if index != self.selected_link:
            self.selected_link = index
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
            old_position = self.chain.Joints[self.selected_joint].Pose.t
            trans = new_position - old_position
            transformation = SE3.Trans(trans[0], trans[1], trans[2])

            if self.chain.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=False):
                self.update_joint()

            self.update_rotation_slider()
            self.update_translate_slider()

    def done_transforming(self, done):
        if done:
            pass
            #self.log_version()

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
            if self.chain.transformJoint(self.selected_joint, transformation, propogate=propogate, safe=False, relative=True):
                self.update_joint()
                self.update_rotation_slider()
                self.update_translate_slider()

    def update_slider(self, slider_type):
        if (slider_type == "rotation"):
            slider = self.rotation_slider
            input = self.rotation_input

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
            slider = self.translate_slider
            input = self.translation_input

            amount = self.chain.Joints[self.selected_joint].Pose.t[self.selected_arrow]
            slider_value = int(amount * 10)
            slider_text = str(slider_value)

            set_slider = self.selected_arrow != -1

            if (set_slider):
                self.old_trans_val = amount

        elif (slider_type == "state"):
            slider = self.state_slider
            input = self.state_input

            min = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[0])
            max = math.degrees(self.chain.Joints[self.selected_joint].stateRange()[1])

            current = int(math.degrees(self.chain.Joints[self.selected_joint].state))
            slider_value = current
            slider_text = str(current)

            set_slider = self.selected_joint != -1 and min != 0 and max != 0

            if (set_slider):
                self.state_slider.setMinimum(int(min))
                self.state_slider.setMaximum(int(max))
                self.old_state_val = current

        if not set_slider:
            slider_value = 0
            slider_text = ""

        slider.blockSignals(True)
        slider.setValue(slider_value)
        slider.setDisabled(False)
        slider.blockSignals(False)
        input.blockSignals(True)
        input.setText(slider_text)
        input.blockSignals(False)
        input.setDisabled(False)

    def update_rotation_slider(self):
        self.update_slider("rotation")

    def update_translate_slider(self):
        self.update_slider("translation")

    def update_state_slider(self):
        self.update_slider("state")

    def rotation_angle_from_matrix(self, rotation_matrix, axis):
        rot = R.from_matrix(rotation_matrix)
        euler_angles = rot.as_euler('xyz', degrees=True)
        return euler_angles[axis]
    
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

        if target_joint.__class__.__name__ == "PrismaticJoint":
            self.edit_dimension_menu.updatePrismatic()
        elif target_joint.__class__.__name__ == "RevoluteJoint":
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
                        new_chain = KinematicChain(new_joint)
                    else:
                        new_chain.append(new_joint, relative=False,
                                        fixedPosition=True, fixedOrientation=True, safe=False)
                else:
                    if new_chain is None:
                        new_chain = KinematicChain(joint)
                    else:
                        new_chain.append(joint, relative=False,
                                        fixedPosition=True, fixedOrientation=True, safe=False)
            
            for idx, joint in enumerate(new_chain.Joints):
                new_chain.setJointState(idx, self._saved_states[idx])
            
            self.chain = new_chain
            self.selected_joint = self._selected_joint
            self.show_success("Chain updated successfully!")
        except Exception as e:
            self.chain = self._backup_chain
            self.show_error("Error rebuilding chain: " + str(e))
        
        self.update_joint()

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
                self.update_joint()
                self.old_rot_val = int(value)
                self.update_rotation_slider()
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
                self.update_joint()
                self.old_trans_val = actualVal
                self.update_translate_slider()
            else:
                self.translate_slider.blockSignals(True)
                self.translate_slider.setValue(int(self.old_trans_val * 10))

    def adjust_state(self, value):
        if not isinstance(value, float) and not isinstance(value, int):
            value = value.strip()
        value = float(value) if value else 0
        actualVal = math.radians(value)
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
                self.state_slider.setValue(int(self.old_state_val))

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

    def rotate_joint(self, angle, axis):
        # Update the position of the spheres
        joint = self.selected_joint

        propogate = self.propogate_slider_checkbox.isChecked()
        localOrient = self.local_orient_slider_checkbox.isChecked()

        transformation = SE3.AngleAxis(angle, [axis[0], axis[1], axis[2]], unit='deg')
        # transformation = SE3.Trans(0, 1, 0)

        # print(transformation)

        self.chain.transformJoint(self.selected_joint, transformation, propogate=propogate, relative=True, safe=False, localOrient=localOrient)

        # self.selected_joint.translate(-cnt[0], -cnt[1], -cnt[2])
        # self.selected_joint.rotate(angle, axis[0], axis[1], axis[2], local=False)
        # self.selected_joint.translate(cnt[0], cnt[1], cnt[2])

        # for i, a in enumerate(self.axes):
        #     a.translate(-cnt[0], -cnt[1], -cnt[2])
        #     a.rotate(angle, axis[0], axis[1], axis[2], local=False)
        #     a.translate(cnt[0], cnt[1], cnt[2])

    def update_joint(self):
        self.edit_dimension_menu.setVisible(False)
        self.edit_dimension_button.setVisible(True)
        self.select_joint_options.blockSignals(True)
        self.select_link_options.blockSignals(True)

        if (not self.stl_generated):
            self.plot_widget.clear()
            self.select_joint_options.clear()
            self.select_link_options.clear()
            self.setCentralWidget(self.plot_widget)

            self.grid = gl.GLGridItem()
            self.grid.setColor(gridColorDefault)

            if self.grid_on:
                self.plot_widget.addItem(self.grid)
        
            if self.mesh_visible and self.referenceMesh is not None:
                self.plot_widget.addItem(self.referenceMesh.mesh)
            
            if self.mesh_selected:
                if (self.control_type == "Translate"):
                    self.referenceMesh.addTranslateArrows(self, selectedArrow=self.selected_arrow, local=self.is_local)
                elif (self.control_type == "Rotate"):
                    self.referenceMesh.addRotateArrows(self, selectedArrow=self.selected_arrow, local=self.is_local)
            
            if self.chain is not None:

                self.select_joint_options.blockSignals(False)
                self.select_link_options.blockSignals(False)

                if (self.control_type == "Translate"):
                    for i, joint in enumerate(self.chain.Joints):
                        if i == self.selected_joint:
                            if (self.selected_frame == -1):
                                joint.addTranslateArrows(self, selectedArrow=self.selected_arrow, local=self.is_local)
                            else:
                                frame_joint = self.chain.Joints[self.selected_frame]
                                joint.addTranslateArrows(self, selectedArrow=self.selected_arrow, local=self.is_local, frame=frame_joint.Pose)
                elif (self.control_type == "Rotate"):
                    for i, joint in enumerate(self.chain.Joints):
                        if i == self.selected_joint:
                            if (self.selected_frame == -1):
                                joint.addRotateArrows(self, selectedArrow=self.selected_arrow, local=self.is_local)
                            else: 
                                frame_joint = self.chain.Joints[self.selected_frame]
                                joint.addRotateArrows(self, selectedArrow=self.selected_arrow, local=self.is_local, frame=frame_joint.Pose)
                
                self.chain.addToWidget(self, selectedJoint=self.selected_joint, selectedLink=self.selected_link, lastJoint = self.last_joint)

            if self.selected_arrow != -1:
                self.rotation_slider.setDisabled(False)
                self.translate_slider.setDisabled(False)
            else:
                self.rotation_slider.setDisabled(True)
                self.translate_slider.setDisabled(True)
                
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
        if not self.add_to_root:
            if (self.chain == None or len(self.chain.Joints) == 0) :
                self.chain = KinematicChain(joint)
            else :
                self.chain.append(joint, relative=True, fixedPosition=True, fixedOrientation=False, safe=False)
                #self.chain.addJoint(self.selected_joint, joint, relative=True, fixedPosition=True, fixedOrientation=False, safe=False)
            self.selected_joint = len(self.chain.Joints) - 1
        else:

            if (self.chain == None or len(self.chain.Joints) == 0) :
                self.chain = KinematicChain(joint)
            else:
                old_root = self.chain.Joints[0]
                joint.Pose = old_root.Pose @ joint.Pose

                new_chain = KinematicChain(joint)

                for i, jt in enumerate(self.chain.Joints):
                    if i == 0:
                        cachedLink = None
                    else:
                        cachedLink = self.chain.Links[i]
                    new_chain.append(jt, relative=False, fixedPosition=True, fixedOrientation=True, safe=False, cachedLink=cachedLink)
                self.chain = new_chain

            self.selected_joint = 0

        self.update_joint()
        self.log_version()
    
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

            waypoint = Waypoint(numSides, self.radius, newPos * newRot)
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
            if (self.chain and len(self.chain.Joints) > 0):
                prevJoint = self.chain.Joints[len(self.chain.Joints)-1]

                if (prevJoint is None):
                    pose = SE3()
                else:
                    prevJoint = self.chain.Joints[0] if self.add_to_root else self.chain.Joints[-1]
                    distance = 4 * self.radius + norm(prevJoint.distalPosition()-prevJoint.Pose.t)
                    if self.add_to_root:
                        distance *= -1
                    pose = SE3(0,0,distance)
                    if prevJoint.pathIndex() == 0:
                        pose = SE3.Ry(np.pi/2) @ pose

                waypoint = Waypoint(numSides, self.radius, pose)
            
            if (self.chain == None):
                waypoint = Waypoint(numSides, self.radius, SE3())
                waypoint_index = 0
            else:
                waypoint_index = len(self.chain.Joints)
            
            if (self.chain == None) or len(self.chain.Joints) == 0:
                waypoint = Waypoint(numSides, self.radius, SE3())
                self.chain = KinematicChain(waypoint)
            elif waypoint_index != 0:
                if self.add_to_root:
                    waypoint.Pose = self.chain.Joints[0].Pose @ waypoint.Pose
                    new_chain = KinematicChain(waypoint)
                    for jt in self.chain.Joints:
                        new_chain.append(jt, relative=False, fixedPosition=True, fixedOrientation=True, safe=False)
                    self.chain = new_chain
                else:
                    self.chain.append(newJoint = waypoint, 
                                    relative=True, fixedPosition=True, fixedOrientation=False, safe=False)
            else:
                self.chain.append(newJoint = waypoint, 
                                    relative=True, fixedPosition=False, fixedOrientation=False, safe=False)

            self.update_joint()
            self.log_version()
            if self.add_to_root:
                self.select_joint_options.setCurrentIndex(0)
            else:
                self.select_joint_options.setCurrentIndex(len(self.chain.Joints) - 1)

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
    window = PointEditorWindow()
    window.show()
    sys.exit(app.exec_())