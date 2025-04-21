import os
import sys
this_dir = os.path.dirname(__file__)
main_dir = os.path.abspath(os.path.join(this_dir, '../..'))
sys.path.append(main_dir)

# Example 1A
from KinematicTree import *
from testqtgraph import *
from makeKinematicTree import *
from optimizationFunctions import *
from randomTree import *
from multiObjectiveOptimizationFunctions import *

r = 1
numSides = 4
jointLength = 0.45/0.245
unextendedRevoluteJointLength = RevoluteJoint(numSides, r, np.pi, SE3()).neutralLength
extensionLength = (jointLength - unextendedRevoluteJointLength)/2

# spec = JointSpecificationTree(RevoluteJoint(numSides, r, np.pi, SE3()))
# i = spec.addJoint(0, PrismaticJoint(numSides,r,3,3,np.pi/5,SE3(4,0,0)@SE3.Ry(np.pi/2)), relative=True)
# i = spec.addJoint(i, RevoluteJoint(numSides,r,np.pi,SE3.Trans(0,0,4)@SE3.Ry(2*np.pi/3)), relative=True)
# i = spec.addJoint(i, RevoluteJoint(numSides,r,np.pi,SE3.Trans(1,0,0)@SE3.Rx(np.pi/4)), relative=True)
# tree = makeTubularKinematicTree(spec, plotSteps=False, orientUp=True)

spec = JointSpecificationTree(Waypoint(numSides, r, Pose=SE3()))
palmJoint = ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, SE3.Trans(4,0,3)@SE3.Ry(-np.pi/4)@SE3.Rx(np.pi/4))
palm = spec.addJoint(0, palmJoint, relative=True)


thumb0Joint = ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength, SE3.Trans(jointLength,0,0)@SE3.Rx(np.pi/2))
thumb0 = spec.addJoint(palm, thumb0Joint, relative=True)

thumb1 = spec.addJoint(thumb0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
                                           SE3.Trans(4,0,0)), relative=True)

thumbEnd = spec.addJoint(thumb1, EndTip(numSides, r, 
        SE3.Trans(3,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), 1), relative=True)


pointer1 = spec.addJoint(0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
            SE3.Trans(0,0,8)@SE3.Ry(-np.pi/2)@SE3.Rx(np.pi)), relative=True)

pointer2 = spec.addJoint(pointer1, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
                        SE3.Trans(4,0,0)), relative=True)


pointerEnd = spec.addJoint(pointer2, EndTip(numSides, r, 
        SE3.Trans(4,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), 1), relative=True)

        

middle0 = spec.addJoint(0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
        SE3.Trans(-2.25,0,4)@SE3.Ry(-np.pi/2)@SE3.Rx(np.pi/2)), relative=True)

middle1 = spec.addJoint(middle0, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
                        SE3.Trans(4,0,0)@SE3.Rx(np.pi/2)), relative=True)

middle2 = spec.addJoint(middle1, ExtendedRevoluteJoint(numSides, r, np.pi, extensionLength,
                        SE3.Trans(4,0,0)), relative=True)

middleEnd = spec.addJoint(middle2, EndTip(numSides, r, 
        SE3.Trans(4,0,0)@SE3.Ry(np.pi/2)@SE3.Rz(np.pi/2), 1), relative=True)

        
tree = makeTubularKinematicTree(spec, plotSteps=False)

print(len(tree.Joints))

states = [[0]* len(tree.Joints),[np.pi/2,np.pi/4,np.pi/2,np.pi/3,np.pi/3,np.pi/3,np.pi/3,np.pi/3,np.pi/2,np.pi/4,np.pi/2,np.pi/3,np.pi/3,np.pi/3,np.pi/3,np.pi/3]]
for t in method1(tree, states, parallelize=True, guarantee=True):
    t.show()
