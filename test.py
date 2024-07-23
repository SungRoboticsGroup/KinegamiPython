from KinematicTree import *
from testqtgraph import *
import time

tree = loadKinematicTree("printedManualHexapod")
plotPrintedTree(tree, "test")
#print(tree.detectCollisions(plot=True, debug=True))