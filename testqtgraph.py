import numpy as np
from stl.mesh import Mesh
import pyqtgraph as pg
import pyqtgraph.opengl as gl
from PyQt5.QtWidgets import QApplication
from PyQt5.QtGui import QColor
from KinematicTree import *
from spatialmath import SE3

def plotPrintedTree(tree : KinematicTree[PrintedJoint], folder : str):
    assert(tree != None)
    
    if (folder != ""):
        os.makedirs(f"scad_output/{folder}", exist_ok=True)
        os.makedirs(f"scad_output/{folder}/poses", exist_ok=True)
        os.makedirs(f"3d_output/{folder}", exist_ok=True)
        os.makedirs(f"3d_output/{folder}/poses", exist_ok=True)
    # Create a QApplication instance
    app = QApplication([])

    # Create a PyQtGraph 3D view
    view = pg.Qt.QtWidgets.QMainWindow()
    widget = pg.Qt.QtWidgets.QWidget()
    view.setCentralWidget(widget)
    layout = pg.Qt.QtWidgets.QVBoxLayout()
    widget.setLayout(layout)

    # Create a 3D graphics view
    view3d = gl.GLViewWidget()
    layout.addWidget(view3d)

    # Set the background color to white
    view3d.setBackgroundColor(QColor('white'))

    # Enable antialiasing for smoother rendering
    view3d.opts['antialiasing'] = True

    # Center the view on the mesh
    # mesh_center = vertices.mean(axis=0)
    # mesh_size = vertices.max(axis=0) - vertices.min(axis=0)
    # view3d.setCameraPosition(distance=mesh_size.max() * 2, elevation=30, azimuth=45)
    # view3d.pan(*mesh_center)

    real_start = time.time()

    #export all the links and plot them
    for i in range(0,len(tree.Children)):
        start = time.time()
        if len(tree.Children[i]) > 0:
            filepath = tree.exportLink3DFile(i,folder + "/poses", pose=True)
            if filepath:
                plotSTL(view3d, filepath, tree.Joints[i].DistalDubinsFrame() @ SE3.Ry(np.pi/2) @SE3.Rz(-np.pi/2))
            print(f"plotted links from {i}, Time: {time.time() - start}s")
    #export and plot all the joints
    for i in range(0,len(tree.Joints)):
        start = time.time()
        if not isinstance(tree.Joints[i],PrintedWaypoint):
            file1, rot1, file2, rot2 = tree.Joints[i].renderPose(folder)
            plotSTL(view3d, file1, tree.Joints[i].ProximalDubinsFrame() @ rot1, color=(0,0,0.5,0))
            if (file2 != None):
                plotSTL(view3d, file2, tree.Joints[i].DistalDubinsFrame() @ rot2, color=(0,0,0.5,0))
        print(f"plotted joint {i}, Time: {time.time() - start}s")
    
    print(f"TOTAL GENERATION TIME: {time.time() - real_start}s")
    #show axes
    grid = gl.GLGridItem()
    grid.setColor((0,0,0,255))
    view3d.addItem(grid)

    print("try1")

    # Show the plot
    view.show()
    app.exec_()

    print("try2")

def plotSTL(view, filepath : str, pose, color = (0.5,0.5,0.5,1)):
    # Load the STL file
    stl_mesh = Mesh.from_file(filepath)

    # Extract vertices and faces from the loaded STL
    vertices = stl_mesh.vectors.reshape(-1, 3)
    faces = np.arange(vertices.shape[0]).reshape(-1, 3)

    matrix = pose.A#SE3.Trans(pose.t).A
    #apply pose

    homogenous_vertices = np.hstack((vertices, np.ones((vertices.shape[0], 1))))
    #transformed_vertices = homogenous_vertices.dot(matrix.T)[:, :3]
    transformed_vertices = (pose * vertices.T).T

    # Create a mesh item and add it to the view
    meshdata = gl.MeshData(vertexes = transformed_vertices, faces = faces)
    mesh = gl.GLMeshItem(meshdata = meshdata, smooth=True, shader="viewNormalColor", color=color, drawEdges=False)
    mesh.setGLOptions('opaque')
    view.addItem(mesh)