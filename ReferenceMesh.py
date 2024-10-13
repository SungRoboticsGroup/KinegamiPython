import pyqtgraph.opengl as gl
from spatialmath import SE3
import numpy as np
from style import *
from Joint import LineItemWithID, LineSphere

class ReferenceMesh():

    def __init__(self, mesh : gl.GLMeshItem):
        self.mesh = mesh
        self.scale = 1
        self.Pose = SE3()

    def updateScale(self, scale):
        self.mesh.scale(scale, scale, scale)
        self.scale = scale

    def generate_extended_line_points(self, point1, point2, gap):
        point1 = np.array(point1)
        point2 = np.array(point2)
        
        direction = point2 - point1
        distance = np.linalg.norm(direction)
        direction = direction / distance

        point1 = point2 - direction * distance * 4
        
        extended_length = 8 * distance
        num_points = int(extended_length / gap) + 1
        start_point = point1 - direction * distance

        points = [start_point + i * gap * direction for i in range(num_points)]
    
        return np.array(points)

    def generate_extended_axis(self, point1, point2, length):
        direction_vector = point2 - point1
    
        midpoint = (point1 + point2) / 2

        p1 = midpoint - direction_vector * length
        p2 = midpoint + direction_vector * length

        return np.array([p1, p2])

    def rotation_matrix(self, axis, theta):
        # rodrigues rotation formula
        axis = np.asarray(axis)
        axis = axis / np.linalg.norm(axis)
        a = np.cos(theta / 2.0)
        b, c, d = -axis * np.sin(theta / 2.0)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                        [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                        [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])
    
    def generate_circle_points(self, axis, center, rad=1.0, num_points=10, rotation=0.0):
        #angles where the points are placed
        angles = self.generate_angles(num_points)

        axis = np.array(axis)
        axis = axis / np.linalg.norm(axis)
        
        #vector not parallel to the axis
        if (axis == [1, 0, 0]).all() or (axis == [-1, 0, 0]).all():
            not_parallel = np.array([0, 1, 0])
        else:
            not_parallel = np.array([1, 0, 0])

        v1 = np.cross(axis, not_parallel)
        v1 = v1 / np.linalg.norm(v1)
        v1 = v1 * rad

        points = []

        for angle in angles:
            R = self.rotation_matrix(axis, angle + rotation)

            point = center + v1

            line_point = np.array(point.tolist()) - center

            rotated_point = np.dot(R, line_point) + center
            points.append(rotated_point)
        
        return points

    def generate_angles(self, num_points=10):
        angles = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
        return angles
    
    def rotate_vector(self, vector, axis='x'):
        if axis == 'x':
            rotation_matrix = np.array([
                [1, 0, 0],
                [0, 0, -1],
                [0, 1, 0]
            ])
        elif axis == 'y':
            rotation_matrix = np.array([
                [0, 0, 1],
                [0, 1, 0],
                [-1, 0, 0]
            ])

        rotated_vector = np.dot(rotation_matrix, vector)
        return rotated_vector
        
    def addTranslateArrows(self, widget, selectedArrow=-1, local=True, frame : SE3 = None):
        rad = 1
        colors = [(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)]

        if selectedArrow != -1:
            colors[selectedArrow] = (1, 1, 1, 1)

        extended_axis_color = [(1, 0, 0, 0.2), (0, 1, 0, 0.2), (0, 0, 1, 0.2)]
        point1 = self.Pose.t

        if local:
            axes = [self.Pose.R[:, i] for i in range(3)]
        else:
            axes = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]

        if frame:
            axes = [frame.R[:, i] for i in range(3)]

        if selectedArrow != -1:
            point2 = point1 + rad * 1 * axes[selectedArrow]
            extended_line_points = self.generate_extended_line_points(point1, point2, 0.1)
            extended_axis = self.generate_extended_axis(point1, point2, 20)
            extended_axis_line = gl.GLLinePlotItem(pos=extended_axis, color=extended_axis_color[selectedArrow], width=5, antialias=True)
            widget.plot_widget.addItem(extended_axis_line)

            for line in extended_line_points:
                md = gl.MeshData.sphere(rows=3, cols=2)
                sphere = LineSphere(meshdata=md, color=lineSphereColor, shader='shaded', smooth=True, position=line)
                sphere.setObjectName("line_sphere")
                sphere.setGLOptions('translucent')
                sphere.scale(2.0, 2.0, 0.1)
                sphere.translate(*line)
                widget.plot_widget.addItem(sphere)

        for i, axis in enumerate(axes):
            pos = np.array([point1, point1 + rad * 1 * axis])
            line = LineItemWithID(pos=pos, color=colors[i], width=10, antialias=True, id=i)
            line.setObjectName("Arrow")
            widget.plot_widget.addItem(line)

    def addRotateArrows(self, widget, selectedArrow=-1, local=True, frame : SE3 = None):
        if local:
            axes = [self.Pose.R[:, i] for i in range(3)]
        else: 
            axes = [np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])]

        if frame:
            axes = [frame.R[:, i] for i in range(3)]

        rad = 1
        colors = [(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)] 

        if selectedArrow != -1:
            colors[selectedArrow] = selectedArrowColor

        center = self.Pose.t
        extended_circle_color = extendedCircleColors
        num_points = 20

        R = self.Pose.R
        theta = [np.arctan2(R[2, 1], R[2, 2]), np.arctan2(R[0, 2], R[0, 0]), np.arctan2(R[1, 0], R[1, 1])]

        # generate the marker points around the arrows
        if selectedArrow != -1:
            points = self.generate_circle_points(axes[selectedArrow], center, rad, num_points, theta[selectedArrow])
            angles = self.generate_angles(num_points)

            for i, point in enumerate(points[:num_points]):
                md = gl.MeshData.sphere(rows=3, cols=3)
                sphere = LineSphere(meshdata=md, color=lineSphereColor, shader='shaded', smooth=True, position=point, rotation=angles[i])
                sphere.setObjectName("rotate_sphere")
                sphere.setGLOptions('translucent')
                sphere.scale(0.1, 0.1, 0.1)
                sphere.translate(*point)
                widget.plot_widget.addItem(sphere)

        #swap the axes so that the selected axis gets rendered last -> appears above the other axes
        numbered_axes = [[0, axes[0]], [1, axes[1]], [2, axes[2]]]

        if selectedArrow != -1:
            numbered_axes[selectedArrow], numbered_axes[2] = numbered_axes[2], numbered_axes[selectedArrow]

        for tul in numbered_axes:    
            i = tul[0]
            axis = tul[1]

            points = self.generate_circle_points(axis, center, rad, num_points, theta[selectedArrow])
            points.append(points[0])

            #generate the lighter circles
            circle = LineItemWithID(pos=points, color=extended_circle_color[i], width=5, antialias=True, id=i)
            circle.setObjectName("circle")
            widget.plot_widget.addItem(circle)

            arrow = np.array(points[-num_points//4:])
            arrow2 = np.array(points[:num_points//4])
            arrows = np.append(arrow, arrow2, axis=0)

            # generate the arrows
            line = LineItemWithID(pos=arrows, color=colors[i], width=15, antialias=True, id=i)
            line.setObjectName("Arrow")
            widget.plot_widget.addItem(line)