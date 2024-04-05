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
chain = KinematicChain(StartTip(numSides, r, Pose=SE3.Trans(-3,-3,0), length=0.5)) 
# TODO: switch to 3 layers, replace all figures accordingly
prismatic = PrismaticJoint(numSides, r, neutralLength=3, numLayers=3, coneAngle=np.pi/4, Pose= SE3.Trans([5,5,0]))
prismaticIndex = chain.append(prismatic)
revolute = RevoluteJoint(numSides, r, np.pi, SE3.Ry(np.pi/4))
revoluteIndex = chain.append(revolute)
end = EndTip(numSides, r, Pose=SE3.Ry(np.pi/2), length=0.5)
endIndex = chain.append(end)
#chain.show(showGlobalFrame=True)


chain.transformJoint(prismaticIndex, SE3.Trans([0,0,-4]), relative=True, propogate=True)

chain.transformJoint(revoluteIndex, SE3.Trans([0,0,-6]), relative=False, propogate=True)

"""
# Example 1B
chain.translateJointAlongAxisOfMotion(revoluteIndex, -7)
#chain.show(showGlobalFrame=True)
chain.translateJointAlongAxisOfMotion(endIndex, -10, applyToPreviousWaypoint=True)
#chain.show(showGlobalFrame=True)
chain.rotateJointAboutAxisOfMotion(revoluteIndex, -np.pi/3)
#chain.show(showGlobalFrame=True)
chain.translateJointAlongAxisOfMotion(prismaticIndex, -2, propogate=False)
#chain.show(showGlobalFrame=True)
#chain.show(showGlobalFrame=False, showJointPoses=False, showLinkPath=False)
pattern = chain.creasePattern()
pattern.save(dxfName="examplePatterns/example1.dxf")

# Example 1C
minPrismaticState, maxPrismaticState = chain.Joints[prismaticIndex].stateRange()
chain.setJointState(prismaticIndex, maxPrismaticState)
chain.setJointState(revoluteIndex, np.pi/2)
#chain.show(showGlobalFrame=True)
"""

chain.addToWidget(window)
window.show()
sys.exit(app.exec_())