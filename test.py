from KinematicTree import *
from testqtgraph import *

tree = loadKinematicTree("algorithmHand")

# tree = tree.globalOptimize()

tree = tree.postOptimize()

print(tree.detectCollisions(plot=True))

tree.show()
plotPrintedTree(origamiToPrinted(tree, 0.05), "test3")