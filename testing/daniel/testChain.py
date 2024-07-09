import os
import sys
child_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(child_dir, '../..'))
sys.path.append(parent_dir)

from geometryHelpers import *
from MyWindow import *
from Joint import *
from KinematicChain import *
app = QApplication(sys.argv)
window = PointEditorWindow()

numSides = 6
r = 1
"""
chain = KinematicChain(StartTip(numSides, r, SE3.Trans(0,0,5), length=0.5))
chain.append(RevoluteJoint(numSides, r, totalBendingAngle=np.pi, Pose=SE3.Rx(np.pi/3)))
"""

# Example 1A
from KinematicChain import *
r = 1
numSides = 4
revolute = RevoluteJoint(numSides, r, np.pi, SE3())
chain = KinematicChain(revolute)
revolute2 = revolute = RevoluteJoint(numSides, r, np.pi, SE3())
chain.append(revolute2)
revolute3 = RevoluteJoint(numSides, r, np.pi, SE3())
chain.addJoint(0, revolute3)
print(chain.Parents)
print(chain.Children)
chain.delete(0)
print(chain.Parents)
print(chain.Children)
chain.transformJoint(0, SE3(0, 8, 0))
# TODO: switch to 3 layers, replace all figures accordingly
#chain.show(showGlobalFrame=True)


# chain.transformJoint(prismaticIndex, SE3.Trans([0,0,-4]), relative=True, propogate=True, relative=False)

# chain.transformJoint(revoluteIndex, SE3.Trans([0,0,-6]), relative=False, propogate=True, relative=False)


# Example 1B
# chain.translateJointAlongAxisOfMotion(revoluteIndex, -7)
#chain.show(showGlobalFrame=True)
# chain.translateJointAlongAxisOfMotion(endIndex, -10, applyToPreviousWaypoint=True)
#chain.show(showGlobalFrame=True)
# chain.rotateJointAboutAxisOfMotion(revoluteIndex, -np.pi/3)
#chain.show(showGlobalFrame=True)
# chain.translateJointAlongAxisOfMotion(prismaticIndex, -2, propogate=False)
#chain.show(showGlobalFrame=True)
#chain.show(showGlobalFrame=False, showJointPoses=False, showLinkPath=False)
# pattern = chain.creasePattern()


# Example 1C
# minPrismaticState, maxPrismaticState = chain.Joints[prismaticIndex].stateRange()
# chain.setJointState(prismaticIndex, maxPrismaticState)
# chain.setJointState(revoluteIndex, np.pi/2)
#chain.show(showGlobalFrame=True)


chain.addToWidget(window)
window.show()
sys.exit(app.exec_())