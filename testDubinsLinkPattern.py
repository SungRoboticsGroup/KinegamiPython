# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 14:49:05 2023

@author: Daniel Feshbach
"""
import dubinsLinkPattern
from dubinsLinkPattern import *

import os
import errno
# https://stackoverflow.com/questions/32123394/workflow-to-create-a-folder-if-it-doesnt-exist-already
def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

def shortestCSCFromFrames(r : float, ProximalDubinsFrame : SE3, DistalDubinsFrame : SE3) -> PathCSC:
    return shortestCSC(r, ProximalDubinsFrame.t, ProximalDubinsFrame.R[:,0], DistalDubinsFrame.t, DistalDubinsFrame.R[:,0])


r = 1
numSides = 6
folder = "./examplePatterns/"
make_sure_path_exists(folder)
subfolder = folder+"r"+str(r)+"n"+str(numSides)+"/"
make_sure_path_exists(subfolder)


ProximalDubinsFrame = SE3()
# DistalDubinsFrame = SE3.Rz(-3*np.pi/4) @ SE3.Tx(3*r) @ SE3.Ry(np.pi/2) @ ProximalDubinsFrame 
DistalDubinsFrame = SE3.Rz(-3*np.pi/4) @ SE3.Tx(-3*r) @ SE3.Ry(np.pi/2) @ ProximalDubinsFrame
path = shortestCSCFromFrames(r, ProximalDubinsFrame, DistalDubinsFrame)
path.plot()
assert(np.all(abs(path.error) < 0.0001))

pattern = dubinsLinkPattern(numSides, r, ProximalDubinsFrame, DistalDubinsFrame, path, splitLongElbowsInto=2)
pattern.makeDXF(saveas=subfolder+"dubinsLink", show=True)
