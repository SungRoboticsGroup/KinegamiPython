# -*- coding: utf-8 -*-
from KinematicChain import *
r = 1
numSides = 6

chain = KinematicChain(RevoluteJoint(numSides, r, np.pi, SE3()))
prismaticIndex = chain.append(PrismaticJoint(numSides, r, neutralLength=3, 
                                    numLayers=3, coneAngle=np.pi/4, Pose=SE3()))
chain.show()
chain.translateJointAlongKinematicAxis(prismaticIndex, -2)
chain.show()
