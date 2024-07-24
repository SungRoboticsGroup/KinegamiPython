from KinematicTree import *
from testqtgraph import *
import time

tree = loadKinematicTree("optimizedHandPrinted")
plotPrintedTree(tree, "test")
#print(tree.detectCollisions(plot=True, debug=True))