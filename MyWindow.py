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
from PathCSC import *
from KinematicChain import KinematicChain
from KinematicTree import KinematicTree

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

        # //////////////////////////////////////////////////////////////////////////////////////
        # //////////////////////////////////    AXIS KEY    ////////////////////////////////////
        # //////////////////////////////////////////////////////////////////////////////////////
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

        # //////////////////////////////////////////////////////////////////////////////////////
        # ////////////////////////////////    EDIT JOINTS    ///////////////////////////////////
        # //////////////////////////////////////////////////////////////////////////////////////
        self.select_joint_options = QComboBox()
        self.rotate_joint_button = QPushButton("Rotate Joint Along Axis of Motion")
        self.translate_joint_button = QPushButton("Translate Joint Along Axis of Motion")
        self.delete_joint_button = QPushButton("Delete Joint")

        joint_layout = QVBoxLayout()
        joint_layout.addWidget(self.select_joint_options)
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
        
        if self.chain:
            self.chain.addToWidget(self, selectedJoint=self.selected_joint)

        self.select_joint_options.blockSignals(False)

    def add_point(self):
                
        #initialize pose
        pose = SE3()
        point = gl.GLScatterPlotItem(pos=pose.t, color=(1, 1, 1, 1), size=10)

        self.points.append(point)
        self.poses.append(pose)

        self.select_point_options.addItem(str(len(self.points)))

        self.update_plot()

    def move_point_x_pos(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Tx(.2)
        self.points[self.current_point].pos = self.poses[self.current_point].t 
        self.update_plot()

    def move_point_x_neg(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Tx(-.2)
        self.points[self.current_point].pos = self.poses[self.current_point].t 
        self.update_plot()

    def move_point_y_pos(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Ty(.2)
        self.points[self.current_point].pos = self.poses[self.current_point].t 
        self.update_plot()

    def move_point_y_neg(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Ty(-.2)
        self.points[self.current_point].pos = self.poses[self.current_point].t 
        self.update_plot()

    def move_point_z_pos(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Tz(.2)
        self.points[self.current_point].pos = self.poses[self.current_point].t 
        self.update_plot()

    def move_point_z_neg(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Tz(-.2)
        self.points[self.current_point].pos = self.poses[self.current_point].t 
        self.update_plot()

    def rotate_point_x(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Rx(math.pi/12)
        self.update_plot()

    def rotate_point_y(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Ry(math.pi/12)
        self.update_plot()

    def rotate_point_z(self):
        self.poses[self.current_point] = self.poses[self.current_point] * SE3.Rz(math.pi/12)
        self.update_plot()

    def index_changed(self, index):
        self.points[self.current_point].color = (1,1,1,1)
        self.current_point = index
        self.points[self.current_point].color = (1,0,0,1)
        self.update_plot()

    def update_plot(self):
        self.plot_widget.clear()
        self.setCentralWidget(self.plot_widget)

        grid = gl.GLGridItem()
        self.plot_widget.addItem(grid)

        index = 0
        for point in self.points:
            md = gl.MeshData.sphere(rows=10, cols=20, radius=0.9)
            center = point.pos
            
            #plot sphere
            m1 = gl.GLMeshItem(
                meshdata=md,
                smooth=True,
                color=(0.4, 0.4, 0.4, 0.1),
                shader="balloon",
                glOptions="additive",
            )
            m1.translate(*center)
            self.plot_widget.addItem(m1)

            #plot point
            self.plot_widget.addItem(point)

            # -------- plot all the poses --------
            pose = self.poses[index]
            center = np.append(point.pos, 1)

            #endpoints
            p_x = np.append(point.pos + [1,0,0], 1)
            p_y = np.append(point.pos + [0,1,0], 1)
            p_z = np.append(point.pos + [0,0,1], 1)

            #apply transformation to a centered point
            rotated_point_x = np.dot(pose, p_x - center)
            rotated_point_y = np.dot(pose, p_y - center)
            rotated_point_z = np.dot(pose, p_z - center)
    
            #translate back from center
            rotated_point_x += center
            rotated_point_y += center
            rotated_point_z += center

            xLine = gl.GLLinePlotItem(pos=np.array([point.pos, rotated_point_x[:3]]), color=(1, 0, 0, 1)) 
            yLine = gl.GLLinePlotItem(pos=np.array([point.pos, rotated_point_y[:3]]), color=(0, 1, 0, 1)) 
            zLine = gl.GLLinePlotItem(pos=np.array([point.pos, rotated_point_z[:3]]), color=(0, 0, 1, 1)) 

            self.plot_widget.addItem(xLine)
            self.plot_widget.addItem(yLine)
            self.plot_widget.addItem(zLine)

            index = index + 1

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

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = PointEditorWindow()
    window.show()
    sys.exit(app.exec_())