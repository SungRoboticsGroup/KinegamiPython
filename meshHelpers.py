"""
Created on Tue Jun 4 2024

@author: fengray
"""

import pyqtgraph.opengl as gl
from Joint import Joint

from style import *

class LinkMesh(gl.GLMeshItem):
    def __init__(self, id : int = -1, lastJoint : Joint = None, nextJoint : Joint = None, **kwds): 
        self.opts = {
            'meshdata': None,
            'color': meshLinkDefaultColor,
            'drawEdges': False,
            'drawFaces': True,
            'edgeColor': meshLinkEdgeColor,
            'shader': None,
            'smooth': True,
            'computeNormals': True,
        }
        
        super().__init__(**kwds)        
        self.id = id

        self.lastJoint = lastJoint
        self.nextJoint = nextJoint