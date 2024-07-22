from KinematicTree import *
from testqtgraph import *
import time

tree = loadKinematicTree("optimizedHexapod")
tree.show()
plotPrintedTree(origamiToPrinted(tree, 1.5/30), "test")
#print(tree.detectCollisions(plot=True, debug=True))