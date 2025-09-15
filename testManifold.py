import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from IPython.display import display
from LinkCSC import *

def plotManifold(manifold, block=True, globalFrame=False):
  # Get mesh representation
  mesh = manifold.to_mesh()
  vertices = mesh.vert_properties[:, :3]
  triangles = mesh.tri_verts

  # Matplotlib 3D plot
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  # Create a list of triangle vertex coordinates
  faces = [vertices[tri] for tri in triangles]
  mesh_collection = Poly3DCollection(faces, alpha=0.7, edgecolor='k')
  ax.add_collection3d(mesh_collection)

  if globalFrame:
     addPosesToPlot(np.array([SE3()]), ax, axisLength=1)

  # Auto scale to the mesh size
  scale = vertices.flatten()
  ax.set_aspect('equal')
  # hide axes
  ax.axis('off')

  plt.show(block=block)

def analyzeProperties(manifold):
  # Get mesh representation
  mesh = manifold.to_mesh()
  vertices = mesh.vert_properties[:, :3]
  triangles = mesh.tri_verts

  volume = manifold.volume()
  surface_area = manifold.surface_area()
  genus = manifold.genus()
  print(f"Volume: {volume}, Surface Area: {surface_area}, Genus: {genus}")
  print(f"Vertices: {vertices.shape}, Triangles: {triangles.shape}")



def save(manifold : m3d.Manifold, filename : str):
    """
    Export a manifold to a 3MF file for 3D printing or CAD applications.
    """
    mesh_data = manifold.to_mesh()
    vertices = mesh_data.vert_properties[:, :3]  # Get XYZ coordinates
    faces = mesh_data.tri_verts
    tri_mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
    
    tri_mesh.export(filename)

"""
# Create a sphere with radius 1.0
sphere = m3d.Manifold.sphere(1.0, 50)

# Create a cube with size 2.0
cube = m3d.Manifold.cube((2.0, 2.0, 2.0))

# Translate the cube
translated_cube = cube.translate((2.0, 0.0, 0.0))

# Perform a boolean union operation
union = sphere + translated_cube

#analyzeProperties(union)
#plotManifold(union)
"""

#link = LinkCSC(r=1, StartDubinsPose=SE3(), EndDubinsPose=SE3.Ry(2*np.pi/3)@SE3.Tx(4.0), maxAnglePerElbow=np.pi/8)
#link.show()

"""
cylinder = m3d.Manifold.cylinder(height=2, radius_low=0.5, radius_high=0.25, circular_segments=6)
plotManifold(cylinder)
"""



  
def sweepRepeatedTransform(points : np.ndarray, Transformation : SE3, 
                           numSteps : int, scale: float = 1.0, StartPose : SE3 = SE3()) -> np.ndarray:
  """
  Repeatedly apply a given SE3 transformation to a set of 2D or 3D points,
  generating a series of shapes along the way.
  Optionally, scale the polygon (with linear interpolation) along the path so the end size is a given
  ratio of the start size.

  Args:
      points (np.ndarray): An Nx2 or Nx3 array of points representing the shape to be swept.
      Transformation (SE3): The SE3 transformation to apply at each step.
      numSteps (int): The number of steps (shapes) to generate. 1 is just the start shape, 2 applies the transformation once, etc.
      scale (float, optional): The ratio of the end size to the start size. Defaults to 1.0.
  """
  if points.shape[1] == 2:
      points = np.hstack((points, np.zeros((points.shape[0], 1))))  # Add Z=0 to make it 3D

  scales = np.linspace(1.0, scale, numSteps)
  sweptPoints = np.zeros((numSteps, points.shape[0], 3))
  sweptPoints[0,:,:] = points
  for i in range(1, numSteps):
      sweptPoints[i,:,:] = (Transformation * sweptPoints[i-1,:,:].T).T
  origin = StartPose.t
  for i in range(1, numSteps):
      origin = Transformation * origin
      sweptPoints[i,:,:] = scales[i] * (sweptPoints[i,:,:] - origin.T) + origin.T
  return sweptPoints

def sweepTransforms(points : np.ndarray, Transformations : list[SE3], 
                    scale: float = 1.0, startOrigin : np.ndarray = None) -> np.ndarray:
  """
  Apply a series of (relative) SE3 transformations to a set of 2D or 3D points,
  generating a series of shapes along the way.
  Optionally, scale the polygon (with linear interpolation) along the path so the end size is a given
  ratio of the start size.

  Args:
      points (np.ndarray): An Nx2 or Nx3 array of points representing the shape to be swept.
      Transformations (list): A list of SE3 transformations to apply.
      scale (float, optional): The ratio of the end size to the start size. Defaults to 1.0.
      relative (bool, optional): If True, each transformation is applied relative to the previous one. Defaults to False.
  """
  if points.shape[1] == 2:
      points = np.hstack((points, np.zeros((points.shape[0], 1))))  # Add Z=0 to make it 3D
  
  numSteps = len(Transformations) + 1
  scales = np.linspace(1.0, scale, numSteps)
  sweptPoints = np.zeros((numSteps, points.shape[0], 3))
  sweptPoints[0,:,:] = points
  for i in range(1, numSteps):
      sweptPoints[i,:,:] = (Transformations[i-1] * sweptPoints[i-1,:,:].T).T
  origin = np.zeros(3) if startOrigin is None else startOrigin
  for i in range(1, numSteps):
      origin = Transformations[i-1] * origin
      sweptPoints[i,:,:] = scales[i] * (sweptPoints[i,:,:] - origin.T) + origin.T
  return sweptPoints


""" 
# Example: Sweep a circle along a bent path
circle = Circle3D(radius=1, center=np.array([0,0,0]), normal=np.array([0,0,1]))
bendingAngle = 3 * np.pi / 4
numSections = 10
anglePerSection = bendingAngle / numSections
distancePerSection = circle.r * np.tan(anglePerSection / 2)
Forward = SE3.Tz(distancePerSection)
Rotate = SE3.AngleAxis(anglePerSection, np.array([0,1,0]))
T = Forward @ Rotate @ Forward
sweep = sweepRepeatedTransform(circle.interpolate(19), T, numSections, scale=1)
print(sweep)

# Plot the swept shape
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(sweep.shape[0]):
    ax.plot(sweep[i,:,0], sweep[i,:,1], sweep[i,:,2], marker='o')

# Plot the origin and the effect of the transformations
Poses = [SE3()]
for i in range(1, numSections):
    Poses.append(Poses[-1] @ T)
Poses = np.array(Poses)
addPosesToPlot(Poses, ax, axisLength=0.2)

ax.set_aspect('equal')
ax.axis('off')
plt.show()
"""
    

"""
# Example usage
elbow = SmoothElbow(arcRadius=1, StartFrame=SE3.Rz(np.pi/4)@SE3(1,2,3), bendingAngle=np.pi/3, rotationalAxisAngle=np.pi/4,
                    startRadius=1, endRadius=0.5, numSides=20, maxSectionAngle=np.pi/10)


# plot the SmoothElbow circles
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
elbow.plotCircles(ax)
ax.set_aspect('equal')
#ax.axis('off')
plt.show()


def plot_trimesh(mesh):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Get vertices and faces
    vertices = mesh.vertices
    faces = mesh.faces

    # Create a list of triangle vertex coordinates
    mesh_faces = [vertices[face] for face in faces]
    mesh_collection = Poly3DCollection(mesh_faces, alpha=0.7, edgecolor='k')
    ax.add_collection3d(mesh_collection)

    # Auto scale to the mesh size
    scale = vertices.flatten()
    ax.auto_scale_xyz(scale, scale, scale)

    ax.set_aspect('equal')
    ax.axis('off')
    plt.show()

# Example usage:
hull = elbow.manifold()
#plotManifold(hull)
"""

"""
bend = Bend(arcRadius=1, StartFrame=SE3.Rz(np.pi/4)@SE3(1,2,3), 
                                bendingAngle=2*np.pi/3, rotationalAxisAngle=np.pi/4,
                                startRadius=1, endRadius=0.5, numSides=20,
                                maxSectionAngle=np.pi/10)
hollowHull = bend.manifold(extendBackward=0.2, extendForward=0.2, wallThickness=0.2, hull=True)
plotManifold(hollowHull)
"""



link = LinkCSC(r=1, StartDubinsPose=SE3.Rz(np.pi/3), 
               EndDubinsPose=SE3.Tx(4.0)@SE3.Ry(3*np.pi/4)@SE3.Rz(3*np.pi/4), 
               maxAnglePerElbow=np.pi/10)
#link.show(showManifold=True, startRadius=1, endRadius=0.75, hullBends=False, wallThickness=0.1, extendBackward=0.2, extendForward=0.2, numSides=20)
module = link.connectableModule(wallThickness=0.1, holeDiameter=0.05, numHoles=4, startRadius=1, endRadius=0.25)
#module = link.manifold(startRadius=1, endRadius=0.75, numSides=20, wallThickness=0.1)
analyzeProperties(module)
plotManifold(module)
link.saveModule("linkModule.obj", wallThickness=0.05, holeDiameter=0.1, numHoles=4, startRadius=1, endRadius=0.5)


"""

holeDiameter = 0.1
wallThickness = 0.2
startRadius = 1.0
endRadius = 0.75
numSides = 20
connectionLength = 2 * holeDiameter
EPS = 0.0001

outset = m3d.Manifold.cylinder(height=2*connectionLength, 
                                       radius_low=endRadius-wallThickness+EPS, 
                                       radius_high=endRadius-wallThickness+EPS,
                                       circular_segments=numSides)
outset -= m3d.Manifold.cylinder(height=2*connectionLength, 
                                       radius_low=endRadius-2*wallThickness, 
                                       radius_high=endRadius-2*wallThickness,
                                       circular_segments=numSides)


holeSlicer = m3d.Manifold()
holeAnglesDegrees = np.linspace(0, 360, 4, endpoint=False)
for angle in holeAnglesDegrees:
    hole = m3d.Manifold.cylinder(height=2*endRadius, radius_low=holeDiameter/2, radius_high=holeDiameter/2, circular_segments=20)
    hole = hole.rotate((0,90,0)).rotate((0,0,angle))
    holeSlicer += hole

outset -= holeSlicer.translate((0,0,1.5*holeDiameter))

plotManifold(outset, globalFrame=True)
"""