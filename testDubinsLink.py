# -*- coding: utf-8 -*-
"""
Created on Thu Dec 14 18:40:04 2023

@author: Daniel Feshbach
"""
from dubinsLink import *


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

ProximalDubinsFrame = SE3(-2,-3,1) @ SE3.Ry(np.pi/5)
# DistalDubinsFrame = SE3.Rz(-3*np.pi/4) @ SE3.Tx(3*r) @ SE3.Ry(np.pi/2) @ ProximalDubinsFrame 
DistalDubinsFrame = SE3.Rz(-3*np.pi/4) @ SE3.Tx(-3*r) @ SE3.Ry(np.pi/3) @ ProximalDubinsFrame


link = LinkCSC(r, ProximalDubinsFrame, DistalDubinsFrame)
pattern = link.creasePattern(numSides)
pattern.makeDXF(saveas=subfolder+"dubinsLink", show=True)
link.plot(showPath=True, showFrames=False)
