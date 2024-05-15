from KinematicTree import *

"""
A tree of Joint objects.
Unlike KinematicTree, this does not have any link geometry connecting the joints.
The point of this class is to represent the types, sizes, and kinematic axes 
specified for the joints: the specific position along and orientation about
this axis will not be considered final in the tree construction algorithm.
"""
class JointSpecificationTree:
    def __init__(self, root: Joint):
        self.Joints = [root]
        self.Children = [[]]
        self.Parents = [-1] # root has no parent

    def addJoint(self, parentIndex : int, newJoint : Joint, relative : bool = True):
        parent = self.Joints[parentIndex]
        if relative:
            newJoint.transformPoseIntoFrame(parent.Pose)
        self.Joints.append(newJoint)
        self.Children.append([])
        self.Parents.append(parentIndex)
        newJointIndex = len(self.Joints)-1
        self.Children[parentIndex].append(newJointIndex)
        return newJointIndex
    
    def depthFirstOrder(self, index : int = 0):
        order = [index]
        for childIndex in self.Children[index]:
            order += self.depthFirstOrder(childIndex)
        return order

    def zHats(self):
        return np.array([joint.Pose.R[:,2] for joint in self.Joints])

"""
Input: 
    jointSpecs, a tree of Joint objects interpreted as specifying their
        types/sizes/etc and kinematic axes (but NOT their specific position 
        along or orientation about those axes) in global coordinates
    r, the tubular radius

Output: a KinematicTree object that instantiates the desired kinematics
        as a tubular tree with radius r
"""
def makeTubularKinematicTree(jointSpecs : JointSpecificationTree, plotSteps : bool = False) -> KinematicTree:
    planeNormal = directionNotOrthogonalToAnyOf(jointSpecs.zHats())
    rootJoint = jointSpecs.Joints[0]

    # set up bounding cylinder
    rootBB = jointSpecs.Joints[0].boundingBall()
    boundingCylinder = Cylinder(rootBB.r, rootBB.c - rootBB.r*planeNormal, 
                                planeNormal, 2*rootBB.r)
    
    cylinders = [copy.deepcopy(boundingCylinder)]

    KT = KinematicTree(rootJoint)

    # map from joint indices in jointSpecs to indices in KT 
    jointIndexSpecsToKT = {0:0} 

    for j in jointSpecs.depthFirstOrder():
        if j == 0:
            continue # we already handled the root above

        newJointSpec = jointSpecs.Joints[j]
        parentIndexInKT = jointIndexSpecsToKT[jointSpecs.Parents[j]]

        if len(KT.Children[parentIndexInKT]) > 0: 
            # if parent already has other children, 
            # add waypoints routing out of the cylinder
            parentJoint = KT.Joints[parentIndexInKT]
            parentJointBB = parentJoint.boundingBall()
            parentJointBBend = parentJointBB.c + parentJointBB.r*planeNormal
            StartPlane = boundingCylinder.startPlane()
            axisToParent = StartPlane.projectionOfPoint(parentJointBBend) - \
                            boundingCylinder.start
            if np.linalg.norm(axisToParent) < 1e-6:
                axisToParent = boundingCylinder.uhat
            axisToParentDir = axisToParent/np.linalg.norm(axisToParent)

            r = jointSpecs.Joints[j].r
            o1 = boundingCylinder.start + (r + boundingCylinder.r)*axisToParentDir + \
                (4*r + StartPlane.signedDistanceToPoint(parentJointBBend))*planeNormal
            Waypoint1 = Waypoint(parentJoint.numSides, r, 
                        SE3.Rt(boundingCylinder.orientation(), o1), pathIndex=2)
            
            parentIndexInKT = KT.addJoint(parentIndexInKT, Waypoint1, 
                                    relative=False, fixedPosition=True, 
                                    fixedOrientation=True, safe=False)
            boundingCylinder.expandToIncludeBall(Waypoint1.boundingBall())

            nDist = boundingCylinder.endPlane().signedDistanceToPoint(o1)
            assert(nDist < 0)
            o2 = o1 - nDist*planeNormal
            Waypoint2 = Waypoint(parentJoint.numSides, r, 
                        SE3.Rt(boundingCylinder.orientation(), o2), pathIndex=2)
            parentIndexInKT = KT.addJoint(parentIndexInKT, Waypoint2, 
                                    relative=False, fixedPosition=True, 
                                    fixedOrientation=True, safe=False)
            boundingCylinder.expandToIncludeBall(Waypoint2.boundingBall())
        
        jInKT = KT.addJoint(parentIndexInKT, newJointSpec,
                    relative=False, fixedPosition=False, fixedOrientation=False, 
                    safe=False, endPlane=boundingCylinder.endPlane())
        jointIndexSpecsToKT[j] = jInKT
        
        if plotSteps:
            ax = plt.figure().add_subplot(projection='3d')
            KT.addToPlot(ax)
            boundingCylinder.addToPlot(ax)
            boundingCylinder.endPlane().addToPlot(ax)
            ax.set_aspect('equal')
            plt.show(block=False)
        
            
        boundingCylinder.expandToIncludeBall(KT.Joints[jInKT].boundingBall())
        boundingCylinder.expandToIncludeBall(KT.Links[jInKT].elbow1BoundingBall)
        boundingCylinder.expandToIncludeBall(KT.Links[jInKT].elbow2BoundingBall)
    
    return KT