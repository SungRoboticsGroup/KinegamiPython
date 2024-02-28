# -*- coding: utf-8 -*-
# Setting up path to be able to import from the parent directory
# https://tutorpython.com/tutorial/python-import-from-parent-directory

import os
import sys
child_dir = os.path.dirname(__file__)
parent_dir = os.path.abspath(os.path.join(child_dir, '..'))
sys.path.append(parent_dir)

from PathCSC import *

path = shortestCSC(1, np.array([0,0,0]), np.array([0,0,1]), np.array([3,0,0]), np.array([0,1,0]))
print(path)
path.show(showCircles=False)