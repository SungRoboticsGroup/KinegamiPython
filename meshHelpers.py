"""
Created on Tue Jun 4 2024

@author: fengray
"""

import pyqtgraph.opengl as gl

class MeshItemWithID(gl.GLMeshItem):
    def __init__(self, id : int = -1, **kwds):
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

class LineItemWithID(gl.GLLinePlotItem):
    def __init__(self, id : int = -1, **kwds):
        """All keyword arguments are passed to setData()"""
        super().__init__(**kwds)
        self.id = id

class LineSphere(gl.GLMeshItem):
    def __init__(self, position = [], rotation=0.0, **kwds):
        super().__init__(**kwds)
        self.position = position
        self.rotation = rotation