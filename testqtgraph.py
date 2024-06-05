import numpy as np
from stl import mesh
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from PyQt5.QtWidgets import QApplication

# Create a QApplication instance
app = QApplication([])

# Load the STL file
stl_mesh = mesh.Mesh.from_file('test.stl')

# Extract vertices and faces from the loaded STL
vertices = stl_mesh.vectors.reshape(-1, 3)
faces = np.arange(vertices.shape[0]).reshape(-1, 3)

# Create a PyQtGraph 3D view
view = pg.Qt.QtWidgets.QMainWindow()
widget = pg.Qt.QtWidgets.QWidget()
view.setCentralWidget(widget)
layout = pg.Qt.QtWidgets.QVBoxLayout()
widget.setLayout(layout)

# Create a 3D graphics view
view3d = gl.GLViewWidget()
layout.addWidget(view3d)

# Create a mesh item and add it to the view
mesh = gl.GLMeshItem(vertexes=vertices, faces=faces, smooth=False, drawEdges=True, edgeColor=(0, 0, 0, 1))
view3d.addItem(mesh)

# Center the view on the mesh
mesh_center = vertices.mean(axis=0)
mesh_size = vertices.max(axis=0) - vertices.min(axis=0)
view3d.setCameraPosition(distance=mesh_size.max() * 2, elevation=30, azimuth=45)
view3d.pan(*mesh_center)

# Show the plot
view.show()
app.exec_()