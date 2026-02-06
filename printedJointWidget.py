import numpy as np
from PyQt5.QtWidgets import QVBoxLayout, QPushButton, QHBoxLayout, QLabel, QDockWidget, QWidget, QLineEdit
from OpenGL.GL import *
from OpenGL.GLU import *
from spatialmath import SE3
import math
from PathCSC import *
from PrintedTube import *
from scipy.spatial.transform import Rotation as R
from style import *

class AddJointMenu(QWidget):
   jointToAdd = None
   

   def __init__(self, parent=None):
       super().__init__(parent)
       self.prevJoint = None
       self.add_to_root = False


   def getJoint(self):
       return self.jointToAdd
      
   def update(self):
       if (self.window().tree and len(self.window().tree.Joints) > 0):
           if not self.window().add_to_root:
               self.prevJoint = self.window().tree.Joints[-1]
               self.add_to_root = False
           else:
               self.prevJoint = self.window().tree.Joints[0]
               self.add_to_root = True
       else:
           self.prevJoint = None
           self.add_to_root = False
  
class AddTransverseRevoluteMenu(AddJointMenu):
   def __init__(self, parent=None):
       super().__init__(parent)
       self.add_joint = self.window().add_joint
       self.initUI()       

   def initUI(self):
       layout = QVBoxLayout()
      
       info_label = QLabel("Add Transverse Revolute Joint (270°)")
       layout.addWidget(info_label)

       apply_button = QPushButton('Add Transverse Revolute Joint')
       apply_button.setAutoDefault(True)
       apply_button.clicked.connect(self.onApplyClicked)
       layout.addWidget(apply_button)

       cancel_button = QPushButton('Cancel')
       cancel_button.clicked.connect(self.window().add_transverse_revolute_toggle)
       layout.addWidget(cancel_button)

       self.setLayout(layout)
       self.update()

   def onApplyClicked(self):
       try:
           self.update()
           
           # Create TransverseRDS3225 joint with 270 degrees, radius 30
           if self.prevJoint is None:
               pose = SE3.Ry(-math.pi/2)
           else:
               # Calculate pose relative to distal Dubins frame of previous joint
               # Transverse revolute faces x direction
               distance = 120  # Approximate spacing for printed joints
               if self.add_to_root: distance *= -1
               pose = SE3.Rt(SE3().R, np.array([distance, 0, 0]))

           self.jointToAdd = TransverseRDS3225(pose, version=270)
           self.add_joint(self.jointToAdd)
           self.window().add_transverse_revolute_toggle()
       except Exception as e:
           self.window().show_error(str(e))


class AddCoaxialRevoluteMenu(AddJointMenu):
   def __init__(self, parent=None):
       super().__init__(parent)
       self.add_joint = self.window().add_joint
       self.initUI()  

   def initUI(self):
       layout = QVBoxLayout()
      
       info_label = QLabel("Add Coaxial Revolute Joint (270°)")
       layout.addWidget(info_label)

       apply_button = QPushButton('Add Coaxial Revolute Joint')
       apply_button.clicked.connect(self.onApplyClicked)
       apply_button.setAutoDefault(True)
       layout.addWidget(apply_button)

       cancel_button = QPushButton('Cancel')
       cancel_button.clicked.connect(self.window().add_coaxial_revolute_toggle)
       layout.addWidget(cancel_button)

       self.setLayout(layout)
       self.update()
      
   def onApplyClicked(self):
       try:
           self.update()
           
           # Create CoaxialRDS3225 joint with 270 degrees, radius 30
           if self.prevJoint is None:
               pose = SE3()
           else:
               # Calculate pose relative to distal Dubins frame of previous joint
               # Coaxial revolute faces x direction
               distance = 120  # Approximate spacing for printed joints
               if self.add_to_root: distance *= -1
               pose = SE3.Rt(SO3.Ry(np.pi/2), np.array([distance, 0, 0]))
           
           self.jointToAdd = CoaxialRDS3225(pose, version=270)
           self.add_joint(self.jointToAdd)
           self.window().add_coaxial_revolute_toggle()
       except Exception as e:
           self.window().show_error(str(e))


class AddTipMenu(AddJointMenu):
   def __init__(self, parent=None):
       super().__init__(parent)
       self.add_joint = self.window().add_joint
       self.initUI()

   def initUI(self):
       layout = QVBoxLayout()

       info_label = QLabel("Add Tip (PrintedHemisphere)")
       layout.addWidget(info_label)
          
       apply_button = QPushButton('Add Tip')
       apply_button.clicked.connect(self.onApplyClicked)
       layout.addWidget(apply_button)
       apply_button.setAutoDefault(True)

       cancel_button = QPushButton('Cancel')
       cancel_button.clicked.connect(self.window().add_tip_toggle)
       layout.addWidget(cancel_button)

       self.setLayout(layout)

   def onApplyClicked(self):
       try:
           self.update()

           # Create PrintedHemisphere with radius 30
           if self.prevJoint is None:
               pose = SE3()
           else:
               # Calculate pose relative to distal Dubins frame of previous joint
               # Tips face z direction, so rotate about y by pi/2
               distance = self.prevJoint.r * 4 + self.window().default_radius/2  # Approximate spacing
               if self.add_to_root: distance *= -1
               pose = SE3.Rt(SE3.Ry(np.pi/2).R, np.array([distance, 0, 0]))

           if self.add_to_root:
               self.jointToAdd = PrintedStartHemisphere(r=30, Pose=pose)
           else:
               self.jointToAdd = PrintedEndHemisphere(r=30, Pose=pose)

           self.add_joint(self.jointToAdd)
           self.window().add_tip_toggle()

       except Exception as e:
           self.window().show_error(str(e))
