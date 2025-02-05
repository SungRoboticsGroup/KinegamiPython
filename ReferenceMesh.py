import pyqtgraph.opengl as gl
from spatialmath import SE3
import numpy as np
from style import *
from pyqtgraph import Transform3D

class ReferenceMesh():

    def __init__(self, mesh : gl.GLMeshItem):
        self.mesh = mesh
        self.Pose = SE3()

    def updateScale(self, scale):
        prev = self.mesh.transform().matrix()
        prev[0][0] = scale
        prev[1][1] = scale
        prev[2][2] = scale
        new_mat = Transform3D(prev)
        self.mesh.setTransform(new_mat)

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
        pass

    def addRotateArrows(self, widget, selectedArrow=-1, local=True, frame : SE3 = None):
        pass