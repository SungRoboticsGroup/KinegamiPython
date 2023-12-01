# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 13:29:46 2023

@author: Daniel Feshbach
Tests for tubularOrigami.py
"""
import tubularOrigami
from tubularOrigami import *

import os
import errno
# https://stackoverflow.com/questions/32123394/workflow-to-create-a-folder-if-it-doesnt-exist-already
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

r = 1
numSides = 6
folder = "./examplePatterns/"
make_sure_path_exists(folder)
subfolder = folder+"r"+str(r)+"n"+str(numSides)+"/"
make_sure_path_exists(subfolder)

composed = TubularPattern(numSides, r)
tube = TubeFittingPattern(numSides, r, 2)
#tube.plotRawGraph(saveas=subfolder+"tube_raw_graph", directed=True)
tube.makeDXF(saveas=subfolder+"tube", show=False)

#tube.plotRawGraph(directed=True)
composed.append(TubeFittingPattern(numSides, r, 3))
twist = TwistFittingPattern(numSides, r, 0.45*np.pi, 1)
#twist.plotRawGraph(saveas=subfolder+"twist_raw_graph", directed=True)
twist.makeDXF(saveas=subfolder+"twist", show=False)
composed.append(TwistFittingPattern(numSides, r, 0.9*np.pi, 1))
composed.append(TubeFittingPattern(numSides, r, 3))
composed.makeDXF(saveas=subfolder+"tube_twist_tube", show=False)
#composed.plotRawGraph(saveas=subfolder+"tube_twist_tube_raw_graph", directed=True)

doubleTwist = TubularPattern(numSides, r)
doubleTwist.append(TwistFittingPattern(numSides, r, 0.9*np.pi, 1))
doubleTwist.append(TwistFittingPattern(numSides, r, 0.3*np.pi, 1))
doubleTwist.makeDXF(saveas=subfolder+"twist_twist", show=False)
#doubleTwist.plotRawGraph(saveas=subfolder+"twist_twist_raw_graph", directed=True)

elbow = ElbowFittingPattern(numSides, r, -np.pi/4, np.pi/3)
#elbow.plotRawGraph(saveas=subfolder+"elbow_raw_graph", directed=True)
elbow.makeDXF(saveas=subfolder+"elbow", show=False)
composed.append(elbow)
composed.makeDXF(saveas=subfolder+"tube_twist_tube_elbow", show=False)

prismatic = PrismaticJointPattern(numSides, r, 1, 3, np.pi/3)
#prismatic.plotRawGraph(saveas=subfolder+"prismatic_raw_graph", directed=True)
prismatic.makeDXF(saveas=subfolder+"prismatic", show=False)
composed.append(prismatic)
composed.makeDXF(saveas=subfolder+"tube_twist_tube_elbow_prismatic", show=False)
composed.append(TubeFittingPattern(numSides, r, 1))
composed.makeDXF(saveas=subfolder+"tube_twist_tube_elbow_prismatic_tube", show=False)

revolute = RevoluteJointPattern(numSides, r, totalBendingAngle=np.pi, numSinkLayers=1)
#revolute.plotRawGraph(saveas=subfolder+"revolute_raw_graph", directed=True)
revolute.makeDXF(saveas=subfolder+"revolute_no_sink", show=False)


revolute = RevoluteJointPattern(numSides, r, totalBendingAngle=np.pi, numSinkLayers=3)
#revolute.plotRawGraph(saveas=subfolder+"revolute_raw_graph", directed=True)
revolute.makeDXF(saveas=subfolder+"revolute", show=False)
composed.append(revolute)
composed.append(TubeFittingPattern(numSides, r, 0.5))
composed.makeDXF(saveas=subfolder+"tube_twist_tube_elbow_prismatic_tube_revolute_tube",
                 show=True)
