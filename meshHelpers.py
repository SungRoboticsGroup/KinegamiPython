"""
Created on Tue Jun 4 2024

@author: fengray
"""

import pyqtgraph.opengl as gl

class MeshItemWithID(gl.GLMeshItem):
    def __init__(self, parentItem=None, id : int = -1, **kwds):
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
        
        super().__init__(parentItem=parentItem)
        glopts = kwds.pop('glOptions', 'opaque')
        self.setGLOptions(glopts)
        shader = kwds.pop('shader', None)
        self.setShader(shader)
        
        self.setMeshData(**kwds)
        
        ## storage for data compiled from MeshData object
        self.vertexes = None
        self.normals = None
        self.colors = None
        self.faces = None
        self.id = id

class LineItemWithID(gl.GLLinePlotItem):
    def __init__(self, parentItem=None, id : int = -1, **kwds):
        """All keyword arguments are passed to setData()"""
        super().__init__(parentItem=parentItem)
        glopts = kwds.pop('glOptions', 'additive')
        self.setGLOptions(glopts)
        self.pos = None
        self.mode = 'line_strip'
        self.width = 1.
        self.color = (1.0,1.0,1.0,1.0)
        self.setData(**kwds)
        self.id = id

class LineSphere(gl.GLMeshItem):
    def __init__(self, parentItem=None, position = [], **kwds):
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
        
        super().__init__(parentItem=parentItem)
        glopts = kwds.pop('glOptions', 'opaque')
        self.setGLOptions(glopts)
        shader = kwds.pop('shader', None)
        self.setShader(shader)
        
        self.setMeshData(**kwds)
        
        ## storage for data compiled from MeshData object
        self.vertexes = None
        self.normals = None
        self.colors = None
        self.faces = None
        self.position = position