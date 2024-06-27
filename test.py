from KinematicTree import *
from testqtgraph import *
import time

tree = loadKinematicTree("algorithmHandPostOptimized")

start = time.time()
print(tree.detectCollisions(plot=True))
print(time.time() - start)