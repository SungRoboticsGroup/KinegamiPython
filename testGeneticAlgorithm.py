from KinematicTree import *

tree = loadKinematicTree("test")

tree = tree.globalOptimize()

tree.show()