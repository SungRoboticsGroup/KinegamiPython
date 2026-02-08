import numpy as np
from PyQt5.QtWidgets import QVBoxLayout, QPushButton, QHBoxLayout, QLabel, QDockWidget, QWidget, QLineEdit
from OpenGL.GL import *
from OpenGL.GLU import *
from spatialmath import SE3
import math
from PathCSC import *
from KinematicChain import *
from OrigamiTube import *
from scipy.spatial.transform import Rotation as R
from style import *

class AddJointMenu(QWidget):
   jointToAdd = None
   

   def __init__(self, parent=None):
       super().__init__(parent)
       self.numSides = 0
       self.r = 0
       self.prevJoint = None
       self.add_to_root = False


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
      
   def update(self):
       self.numSides = self.window().num_sides
       self.r = self.window().radius


       if (self.window().chain and len(self.window().chain.Joints) > 0):
           if not self.window().add_to_root:
               self.prevJoint = self.window().chain.Joints[-1]
               self.add_to_root = False
           else:
               self.prevJoint = self.window().chain.Joints[0]
               self.add_to_root = True
       else:
           self.prevJoint = None
           self.add_to_root = False
  
class AddPrismaticMenu(AddJointMenu):
   def __init__(self, parent=None):
       super().__init__(parent)
       self.add_joint = self.window().add_joint


       self.initUI()       


   def initUI(self):
       layout = QVBoxLayout()
      
       length_layout = QHBoxLayout()
       length_label = QLabel("Neutral Length (default: 3r):")
       self.length_input = QLineEdit()
       self.length_input.returnPressed.connect(self.onApplyClicked)
       length_layout.addWidget(length_label)
       length_layout.addWidget(self.length_input)
       layout.addLayout(length_layout)


       numLayers_layout = QHBoxLayout()
       numLayers_label = QLabel("Number of Layers (default: 3):")
       self.numLayers_input = QLineEdit()
       self.numLayers_input.returnPressed.connect(self.onApplyClicked)
       numLayers_layout.addWidget(numLayers_label)
       numLayers_layout.addWidget(self.numLayers_input)
       layout.addLayout(numLayers_layout)


       angle_layout = QHBoxLayout()
       angle_label = QLabel("Cone Angle (degrees, default: 60):")
       self.angle_input = QLineEdit()
       self.angle_input.returnPressed.connect(self.onApplyClicked)
       angle_layout.addWidget(angle_label)
       angle_layout.addWidget(self.angle_input)
       layout.addLayout(angle_layout)


       apply_button = QPushButton('Add Prismatic Joint')
       apply_button.setAutoDefault(True)
       apply_button.clicked.connect(self.onApplyClicked)
       layout.addWidget(apply_button)


       cancel_button = QPushButton('Cancel')
       cancel_button.clicked.connect(self.window().add_prismatic_toggle)
       layout.addWidget(cancel_button)


       self.setLayout(layout)


       self.update()


   def onApplyClicked(self):
       try:
           self.update()


           neutralLength = 3*self.r if self.length_input.text()=="" else float(self.length_input.text())
           numLayers = 3 if self.numLayers_input.text()=="" else int(self.numLayers_input.text())
           coneAngleText = 60 if self.angle_input.text()=="" else float(self.angle_input.text())


           if (self.prevJoint is None):
               pose = SE3()
           elif self.add_to_root:
               # Compute pose in global coordinates behind the old root
               root = self.window().chain.Joints[0]
               root_proximal_dubins = root.ProximalDubinsFrame()
               r = root.r
               distance = 4 * r + neutralLength / 2
               # Prismatic has pathIndex=2, so rotate by Ry(pi/2) so z-hat aligns with dubins x-hat
               pose = root_proximal_dubins @ SE3.Rt(SE3.Ry(np.pi/2).R, np.array([-distance, 0, 0]))
           else:
               # Calculate pose relative to distal Dubins frame of previous joint
               distance = 4 * self.r + norm(self.prevJoint.distalPosition()-self.prevJoint.Pose.t) + neutralLength/2
               # Prismatic joints face z direction, so rotate about y by pi/2
               pose = SE3.Rt(SE3.Ry(np.pi/2).R, np.array([distance, 0, 0]))


           self.jointToAdd = OrigamiPrismatic(self.numSides, self.r, neutralLength, numLayers, math.radians(coneAngleText), pose)
           self.add_joint(self.jointToAdd)
           self.window().add_prismatic_toggle()
       except ValueError:
           self.window().show_error("Please enter valid integers.")
           # error_dialog = ErrorDialog('Please enter valid integers.')
           # error_dialog.exec_()




class AddRevoluteMenu(AddJointMenu):
   def __init__(self, parent=None):
       super().__init__(parent)
       self.add_joint = self.window().add_joint


       self.initUI()  


   def initUI(self):
       layout = QVBoxLayout()
      
       angle_layout = QHBoxLayout()
       angle_label = QLabel("Total Bending Angle (degrees, default: 180):")
       self.angle_input = QLineEdit()
       self.angle_input.returnPressed.connect(self.onApplyClicked)
       angle_layout.addWidget(angle_label)
       angle_layout.addWidget(self.angle_input)
       layout.addLayout(angle_layout)


       apply_button = QPushButton('Add Revolute Joint')
       apply_button.clicked.connect(self.onApplyClicked)
       apply_button.setAutoDefault(True)
       layout.addWidget(apply_button)


       cancel_button = QPushButton('Cancel')
       cancel_button.clicked.connect(self.window().add_revolute_toggle)
       layout.addWidget(cancel_button)


       self.setLayout(layout)


       self.update()
      
   def onApplyClicked(self):
       self.update()

       bendingAngleText = 180 if self.angle_input.text()=="" else float(self.angle_input.text())
       
       # Validate bending angle
       if bendingAngleText <= 0 or bendingAngleText >= 360:
           from PyQt5.QtWidgets import QMessageBox
           QMessageBox.warning(self, "Invalid Angle", 
                             "Bending angle must be between 0 and 360 degrees (exclusive).\n"
                             "Please enter a value like 180 or 270.")
           return
      
       # Calculate pose
       if self.prevJoint is None:
           pose = SE3()
       elif self.add_to_root:
           # Compute pose in global coordinates behind the old root
           root = self.window().chain.Joints[0]
           root_proximal_dubins = root.ProximalDubinsFrame()
           r = root.r
           # Compute neutralLength for OrigamiRevolute (outerLength=0)
           polygonInnerAngle = np.pi * (self.numSides - 2) / (2 * self.numSides)
           neutralLength = 2 * r * np.sin(polygonInnerAngle) * np.tan(math.radians(bendingAngleText) / 4)
           distance = 4 * r + neutralLength / 2
           # Revolute has pathIndex=0, same orientation as root's proximal dubins frame
           pose = root_proximal_dubins @ SE3.Trans(-distance, 0, 0)
       else:
           # Calculate pose relative to distal Dubins frame of previous joint
           distance = 4 * self.r + norm(self.prevJoint.distalPosition()-self.prevJoint.Pose.t)
           # Revolute joints face x direction, so no rotation needed
           pose = SE3.Rt(SE3().R, np.array([distance, 0, 0]))
       
       self.jointToAdd = OrigamiRevolute(self.numSides, self.r, math.radians(bendingAngleText), pose)
  
       self.add_joint(self.jointToAdd)
       self.window().add_revolute_toggle()


class AddTipMenu(AddJointMenu):
   def __init__(self, parent=None):
       super().__init__(parent)
       self.add_joint = self.window().add_joint


       self.isStart = True


       self.initUI()


   def initUI(self):
       layout = QVBoxLayout()


       length_layout = QHBoxLayout()
       length_label = QLabel("Length:")
       self.length_input = QLineEdit()
       self.length_input.returnPressed.connect(self.onApplyClicked)
       length_layout.addWidget(length_label)
       length_layout.addWidget(self.length_input)
       layout.addLayout(length_layout)
          
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


           length = float(self.length_input.text())
           if self.prevJoint is None:
               self.jointToAdd = OrigamiStartTip(self.numSides, self.r, SE3(), length=length)
           elif self.add_to_root:
               # Compute pose in global coordinates behind the old root
               root = self.window().chain.Joints[0]
               root_proximal_dubins = root.ProximalDubinsFrame()
               r = root.r
               distance = 4 * r + length / 2
               # Tip has pathIndex=2, so rotate by Ry(pi/2) so z-hat aligns with dubins x-hat
               pose = root_proximal_dubins @ SE3.Rt(SE3.Ry(np.pi/2).R, np.array([-distance, 0, 0]))
               self.jointToAdd = OrigamiStartTip(self.numSides, self.r, pose, length=length)
           else:
               # Calculate pose relative to distal Dubins frame of previous joint
               distance = 4*self.r + length/2
               # Tips face z direction, so rotate about y by pi/2
               pose = SE3.Rt(SE3.Ry(np.pi/2).R, np.array([distance, 0, 0]))
               self.jointToAdd = OrigamiEndTip(self.numSides, self.r, pose, length=length)


           self.add_joint(self.jointToAdd)
           self.window().add_tip_toggle()


       except ValueError:
           self.window().show_error("Please enter valid integers.")
           # error_dialog = ErrorDialog('Please enter valid integers.')
           # error_dialog.exec_()


class EditDimensionMenu(AddJointMenu):


   def __init__(self, parent=None):
       super().__init__(parent)
       self.mode = None
       self.main_layout = QVBoxLayout()
       self.setLayout(self.main_layout)


   def clearLayout(self, layout):
       while layout.count():
           item = layout.takeAt(0)
           widget = item.widget()
           if widget is not None:
               widget.setParent(None)
               widget.deleteLater()
           elif item.layout() is not None:
               self.clearLayout(item.layout())


   def rebuildLayout(self, new_layout):
       self.clearLayout(self.main_layout)
       while new_layout.count():
           item = new_layout.takeAt(0)
           if item.widget():
               self.main_layout.addWidget(item.widget())
           elif item.layout():
               self.main_layout.addLayout(item.layout())
       new_layout.deleteLater()


   def updatePrismatic(self):
       self.mode = "prismatic"
       self.numSides = self.window().num_sides
       self.r = self.window().radius


       new_layout = QVBoxLayout()


       length_layout = QHBoxLayout()
       length_label = QLabel("Neutral Length (default: 3r):")
       self.length_input = QLineEdit()
       length_layout.addWidget(length_label)
       length_layout.addWidget(self.length_input)
       new_layout.addLayout(length_layout)


       numLayers_layout = QHBoxLayout()
       numLayers_label = QLabel("Number of Layers (default: 3):")
       self.numLayers_input = QLineEdit()
       numLayers_layout.addWidget(numLayers_label)
       numLayers_layout.addWidget(self.numLayers_input)
       new_layout.addLayout(numLayers_layout)


       angle_layout = QHBoxLayout()
       angle_label = QLabel("Cone Angle (degrees, default: 60):")
       self.angle_input = QLineEdit()
       angle_layout.addWidget(angle_label)
       angle_layout.addWidget(self.angle_input)
       new_layout.addLayout(angle_layout)


       apply_button = QPushButton("Edit Prismatic Joint")
       apply_button.clicked.connect(self.onApplyClicked)
       new_layout.addWidget(apply_button)


       cancel_button = QPushButton("Cancel")
       cancel_button.clicked.connect(self.window().edit_dimension_toggle)
       new_layout.addWidget(cancel_button)


       self.rebuildLayout(new_layout)


   def updateRevolute(self):
       self.mode = "revolute"
       self.numSides = self.window().num_sides
       self.r = self.window().radius
      
       new_layout = QVBoxLayout()


       angle_layout = QHBoxLayout()
       angle_label = QLabel("Total Bending Angle (degrees, default: 180):")
       self.angle_input = QLineEdit()
       angle_layout.addWidget(angle_label)
       angle_layout.addWidget(self.angle_input)
       new_layout.addLayout(angle_layout)


       apply_button = QPushButton("Edit Revolute Joint")
       apply_button.clicked.connect(self.onApplyClicked)
       new_layout.addWidget(apply_button)


       cancel_button = QPushButton("Cancel")
       cancel_button.clicked.connect(self.window().edit_dimension_toggle)
       new_layout.addWidget(cancel_button)


       self.rebuildLayout(new_layout)


   def updateTip(self):
       self.mode = "tip"
       self.numSides = self.window().num_sides
       self.r = self.window().radius
      
       new_layout = QVBoxLayout()


       length_layout = QHBoxLayout()
       length_label = QLabel("Length:")
       self.length_input = QLineEdit()
       length_layout.addWidget(length_label)
       length_layout.addWidget(self.length_input)
       new_layout.addLayout(length_layout)


       apply_button = QPushButton("Edit Tip")
       apply_button.clicked.connect(self.onApplyClicked)
       new_layout.addWidget(apply_button)


       cancel_button = QPushButton("Cancel")
       cancel_button.clicked.connect(self.window().edit_dimension_toggle)
       new_layout.addWidget(cancel_button)


       self.rebuildLayout(new_layout)


   def onApplyClicked(self):
       try:
           self.numSides = self.window().num_sides
           self.r = self.window().radius
           self.prev_joint = self.window().prev_joint


           if self.mode == "prismatic":
               neutral_length = float(self.length_input.text()) if self.length_input.text() != "" else 3 * self.r
               num_layers = int(self.numLayers_input.text()) if self.numLayers_input.text() != "" else 3
               cone_angle = float(self.angle_input.text()) if self.angle_input.text() != "" else 60.0


               if self.prev_joint is None:
                   pose = SE3()
               else:
                   # Get the distal Dubins frame to account for joint state
                   distal_dubins_frame = self.prev_joint.DistalDubinsFrame()
                   distance = 4 * self.r + neutral_length / 2
                   if self.add_to_root:
                       distance *= -1
                   # Translate distance along the distal frame's x-axis (column 0)
                   # Then rotate so prismatic's z-axis (pathIndex 2) aligns with Dubins x-axis
                   new_position = distal_dubins_frame.t + distance * distal_dubins_frame.R[:,0]
                   new_orientation = distal_dubins_frame.R @ SE3.Ry(-np.pi/2).R
                   pose = SE3.Rt(new_orientation, new_position)


               self.editJoint = OrigamiPrismatic(self.numSides, self.r, neutral_length,
                                                 num_layers, math.radians(cone_angle), pose)


           elif self.mode == "revolute":
               bending_angle = float(self.angle_input.text()) if self.angle_input.text() != "" else 180.0


               self.editJoint = OrigamiRevolute(self.numSides, self.r, math.radians(bending_angle), SE3())
               if self.prev_joint is not None:
                   # Get the distal Dubins frame to account for joint state
                   distal_dubins_frame = self.prev_joint.DistalDubinsFrame()
                   distance = 4 * self.r + self.editJoint.neutralLength / 2
                   if self.add_to_root:
                       distance *= -1
                   # Translate distance along the distal frame's x-axis (column 0)
                   new_position = distal_dubins_frame.t + distance * distal_dubins_frame.R[:,0]
                   pose = SE3.Rt(distal_dubins_frame.R, new_position)
                   self.editJoint.Pose = pose


           elif self.mode == "tip":
               length = float(self.length_input.text())
               if self.prev_joint is None:
                   self.editJoint = StartTip(self.numSides, self.r, SE3(), length=length)
               else:
                   # Get the distal Dubins frame to account for joint state
                   distal_dubins_frame = self.prev_joint.DistalDubinsFrame()
                   distance = 4 * self.r + length / 2
                   if self.add_to_root:
                       distance *= -1
                   # Translate distance along the distal frame's x-axis (column 0)
                   new_position = distal_dubins_frame.t + distance * distal_dubins_frame.R[:,0]
                   pose = SE3.Rt(distal_dubins_frame.R, new_position)
                   self.editJoint = EndTip(self.numSides, self.r, pose, length=length)

           self.window().log_version()
           self.window().edit_dimension_toggle()
           self.window().finish_joint_edit(self.editJoint)
       except Exception as e:
           self.window().show_error(str(e))