"""
Created on Tue Jun 4 2024

@author: fengray
"""

import pyqtgraph.opengl as gl
from Joint import Joint

class LinkMesh(gl.GLMeshItem):
    def __init__(self, id : int = -1, lastJoint : Joint = None, nextJoint : Joint = None, **kwds): 
        self.opts = {
            'meshdata': None,
            'color': (1., 1., 1., 1.),
            'drawEdges': False,
            'drawFaces': True,
            'edgeColor': (0.5, 0.5, 0.5, 1.0),
            'shader': None,
            'smooth': True,
            'computeNormals': True,
        }
        
        super().__init__(**kwds)        
        self.id = id

        self.lastJoint = lastJoint
        self.nextJoint = nextJoint