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

tree = loadKinematicTree("old/10 Joint Chains/0/Linear (L)/4.545870780944824_1.tree")
plotColoredTrees("old/10 Joint Chains/0/Linear (L)")