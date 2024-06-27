# How the additional code in this branch is structured

## From KinematicTree.py

### Conversion functions
- `origamiToPrinted(tree : KinematicTree[OrigamiJoint], screwRadius: float)`
    - converts an origami tree to a 3d print tree (if not possible, this will throw an error and you may need to adjust spacing)
- `printedToOrigami(tree : KinematicTree[PrintedJoint], numSides: int, numLayers : int = 1)`
    - converts a 3d print tree to an origami tree (if not possible, this will throw and error and you may need to adjust spacing)

### Optimization functions
- `def postOptimize(self)`
    - optimizes joint placement of a tree to make it as small as possible, finding a collision-free local minimum (assuming initial tree is collision-free)
- `def globalOptimize(self, iterations=10)`
    - global optimization to be applied after finding local minimum with postOptimize (**NOTE: STILL IN PROGRESS**)

### Collision detection
- `def detectCollisions(self, specificJointIndices = [], plot=False)`
    - returns the number of collisions detected in the kinematic tree. `specificJointIndices` is a list of specific joint indices to specifically test for collision, and if empty all joints will be evaluated. if plot is true, the tree will be plotted and any bodies that collide will be highlighted.

### Saving functions
- `def save(self, filename : str)`
    - saves a kinematic tree (only works for origami so far) as &lt;filename&gt;.tree in the "**save**" folder
- `def loadKinematicTree(filename : str)`
    - loads a kinematic tree (only works for origami so far) from &lt;filename&gt;.tree in the "**save**" folder

### Exporting
- `def export3DKinematicTree(self, folder = "", fileFormat = "stl")`
    - exports all of the modules required to construct the 3D print tree into the "**3d_output/&lt;folder&gt;**" folder. each joint file is labeled as joint_&lt;index&gt;, and each link file is labeled as link_from_&lt;index1&gt;_to_&lt;index2&gt;
    - the all of the corresponding scad files used to export the final 3d files using openscad are stored in the "**scad_output/&lt;folder&gt;**" folder

### From testqtgraph.py
**NOTE: requires numpy-stl and pyqtgraph**
- `def plotPrintedTree(tree : KinematicTree[PrintedJoint], folder : str)`
    - uses pyqtgraph to plot a 3d print tree. 
    - all stl files used to generate meshes will be in the "**3d_output/&lt;folder&gt;/poses**" folder.
    - all of the corresponding scad files used to export the stl files are stored in the "**scad_output/&lt;folder&gt;/poses**" folder
    - no files are regenerated as long as the joint/link dimensions do not change
