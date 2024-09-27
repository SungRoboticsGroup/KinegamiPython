import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
from spatialmath import SE3
import time

class CollisionBox:
    def __init__(self, points=None, startDubinsFrame=None, endDubinsFrame=None,r=None):
        if points is None:
            if startDubinsFrame is None or endDubinsFrame is None or r is None:
                self.points = np.zeros((8, 3))
            else:
                self.width = r*2
                self.axis = endDubinsFrame.t - startDubinsFrame.t
                self.height = np.linalg.norm(self.axis)
                self.center = startDubinsFrame.t + self.axis/2

                self.points = np.array([
                    startDubinsFrame.t + r * startDubinsFrame.n + r * startDubinsFrame.o,
                    startDubinsFrame.t - r * startDubinsFrame.n - r * startDubinsFrame.o,
                    startDubinsFrame.t - r * startDubinsFrame.n + r * startDubinsFrame.o,
                    startDubinsFrame.t + r * startDubinsFrame.n - r * startDubinsFrame.o,
                    endDubinsFrame.t + r * endDubinsFrame.n + r * endDubinsFrame.o,
                    endDubinsFrame.t - r * endDubinsFrame.n - r * endDubinsFrame.o,
                    endDubinsFrame.t - r * endDubinsFrame.n + r * endDubinsFrame.o,
                    endDubinsFrame.t + r * endDubinsFrame.n - r * endDubinsFrame.o
                ])
        else:
            self.points = np.array(points)
    

    def rotate(self, angle):
        if self.height == 0:
            return self.points

        newPoints = rotate_box_arbitrary_axis(self.points, angle, self.center - self.axis/2, self.axis)
        return newPoints
    
    def addToPlot(self, ax, color='r', alpha=0.4):
        try:
            hull = ConvexHull(self.points)#, qhull_options='Qs QJ')

            vertices = hull.points[hull.vertices]
            
            faces = hull.simplices
            
            # create the 3D polygon collection
            poly3d = [[vertices[vidx] for vidx in face] for face in faces]
            box_faces = Poly3DCollection(poly3d, facecolors=color, linewidths=1, edgecolors=None, alpha=alpha)
            
            # add box to plot
            ax.add_collection3d(box_faces)
            for plotPoint in self.points:
                ax.scatter(plotPoint[0], plotPoint[1], plotPoint[2], color='black', s=25)
        except Exception as e:
            pass

def separatingAxisTheorem(b1, b2):
    """
    Check if two arbitrarily oriented boxes collide using the Separating Axis Theorem.
    
    :param box1: List of 8 points (x, y, z) representing the corners of the first box
    :param box2: List of 8 points (x, y, z) representing the corners of the second box
    :return: Boolean indicating whether the boxes collide
    """

    if b1.height < 1e-4 or b2.height < 1e-4:
        return False

    if np.linalg.norm(b1.center - b2.center) > np.linalg.norm([b1.height,b1.width,b1.width]) / 2 + np.linalg.norm([b2.height,b2.width,b2.width]) / 2:
        return False

    box1 = np.array(b1.points)
    box2 = np.array(b2.points)

    def get_box_axes(box):
        # Compute three edges from a corner
        edges = [
            box[1] - box[0],
            box[3] - box[0],
            box[4] - box[0]
        ]
        # Normalize the edges, handling zero-length edges
        return [normalize(edge) for edge in edges if not np.allclose(edge, 0)]

    def normalize(v):
        norm = np.linalg.norm(v)
        if norm == 0:
            return v
        return v / norm
    
    def project_box_onto_axis(box, axis):
        projections = np.dot(box, axis)
        return np.min(projections), np.max(projections)

    axes1 = get_box_axes(box1)
    axes2 = get_box_axes(box2)

    # Check all 15 potential separating axes
    for axis in axes1 + axes2 + [np.cross(a, b) for a in axes1 for b in axes2]:
        if np.all(axis == 0):
            continue  # Skip zero vectors resulting from parallel edges

        projection1 = project_box_onto_axis(box1, axis)
        projection2 = project_box_onto_axis(box2, axis)

        if projection1[1] < projection2[0] or projection2[1] < projection1[0]:
            return False  # Separating axis found, no collision

    return True  # No separating axis found, boxes collide

def rotate_box_arbitrary_axis(box, angle_deg, axis_point, axis_direction):
    """
    Rotate a box around an arbitrary axis positioned at a specific point in space.
    
    :param box: Array of shape (8, 3) representing the box corners
    :param angle_deg: Rotation angle in degrees
    :param axis_point: 3D point that the axis passes through
    :param axis_direction: 3D vector representing the direction of the axis
    :return: Rotated box
    """
    angle_rad = np.radians(angle_deg)
    axis_point = np.array(axis_point)
    axis_direction = np.array(axis_direction)
    axis_direction = axis_direction / np.linalg.norm(axis_direction)  # Normalize the axis direction
    
    # Translate box so that rotation axis passes through origin
    translated_box = box - axis_point
    
    # Compute rotation matrix (Rodrigues' rotation formula)
    K = np.array([
        [0, -axis_direction[2], axis_direction[1]],
        [axis_direction[2], 0, -axis_direction[0]],
        [-axis_direction[1], axis_direction[0], 0]
    ])
    rotation_matrix = (
        np.eye(3) + 
        np.sin(angle_rad) * K + 
        (1 - np.cos(angle_rad)) * np.dot(K, K)
    )
    
    # Apply rotation and translate back
    rotated_box = np.dot(translated_box, rotation_matrix.T) + axis_point
    
    return rotated_box

def vectorized_box_collision(boxes1, boxes2):
    """
    Vectorized function to check collisions between multiple pairs of boxes.
    
    :param boxes1: Array of shape (N, 8, 3) representing N boxes
    :param boxes2: Array of shape (N, 8, 3) representing N boxes
    :return: Boolean array of shape (N,) indicating collisions
    """
    def get_box_axes(boxes):
        edges = boxes[:, [1, 3, 4]] - boxes[:, 0][:, np.newaxis]
        norms = np.linalg.norm(edges, axis=2, keepdims=True)
        norms = np.where(norms == 0, 1, norms)  # Avoid division by zero
        return edges / norms

    def project_boxes_onto_axis(boxes, axis):
        projections = np.dot(boxes, axis)
        return np.min(projections, axis=1), np.max(projections, axis=1)

    axes1 = get_box_axes(boxes1)
    axes2 = get_box_axes(boxes2)

    all_axes = np.concatenate([
        axes1.reshape(-1, 3),
        axes2.reshape(-1, 3),
        np.cross(axes1[:, :, np.newaxis], axes2[:, np.newaxis, :]).reshape(-1, 3)
    ])

    collisions = np.ones(len(boxes1), dtype=bool)

    for axis in all_axes:
        if np.allclose(axis, 0):
            continue

        proj1_min, proj1_max = project_boxes_onto_axis(boxes1, axis)
        proj2_min, proj2_max = project_boxes_onto_axis(boxes2, axis)

        collisions &= (proj1_max >= proj2_min) & (proj2_max >= proj1_min)

        if not np.any(collisions):
            break

    return collisions

class CollisionCapsule:
    def __init__(self, base : SE3(), radius : float, height : float):
        self.radius = radius
        self.height = height
        self.base = base @ SE3.Ry(np.pi/2)
        self.otherBase = self.base @ SE3.Trans([0,0,self.height])
        self.start = (self.base @ SE3.Trans([0,0,-self.radius])).t
        self.center = (self.base @ SE3.Trans([0,0,self.height/2])).t
        self.end = (self.base @ SE3.Trans([0,0,self.height + self.radius])).t
        self.halfDist = np.abs(self.height/2) + self.radius

        self.box = CollisionBox(startDubinsFrame=self.base, endDubinsFrame=self.otherBase, r=self.radius)

    
    def addToPlot(self, ax, color='r', alpha=0.25):
        def createCirclePoints(center, radius, normal, num_points=50):
            """Generate points on a circle in 3D."""
            angles = np.linspace(0, 2 * np.pi, num_points)
            circlePoints = []
            for angle in angles:
                point = center + radius * (np.cos(angle) * orthogonalVector(normal) + np.sin(angle) * np.cross(normal, orthogonalVector(normal)))
                circlePoints.append(point)
            return np.array(circlePoints)
        
        def orthogonalVector(v):
            """Find an orthogonal vector to v."""
            if np.isclose(v[0], 0) and np.isclose(v[1], 0):
                if np.isclose(v[2], 0):
                    raise ValueError("Zero vector")
                return np.array([0, 1, 0])
            return np.array([-v[1], v[0], 0])

        try:
            bottomCenter = self.base.t
            topCenter = bottomCenter + self.height * self.base.R[:,2]

            bottomCirclePoints = createCirclePoints(bottomCenter, self.radius, self.base.R[:,2])
            topCirclePoints = createCirclePoints(topCenter, self.radius, self.base.R[:,2])

            theta, phi = np.mgrid[0:np.pi:10j, 0:2*np.pi:10j]
            x = self.radius * np.sin(theta) * np.cos(phi)
            y = self.radius * np.sin(theta) * np.sin(phi)
            z = self.radius * np.cos(theta)

            bottomCap = np.vstack((x.flatten(), y.flatten(), z.flatten())).T
            bottomCap = bottomCap @ self.base.R.T + bottomCenter

            topCap = np.vstack((x.flatten(), y.flatten(), -z.flatten())).T
            topCap = topCap @ self.base.R.T + topCenter

            allPoints = np.vstack((bottomCirclePoints, topCirclePoints, bottomCap, topCap))
            hull = ConvexHull(allPoints)

            # for simplex in hull.simplices:
            #     ax.plot(allPoints[simplex, 0], allPoints[simplex, 1], allPoints[simplex, 2], 'k-')
            ax.add_collection3d(Poly3DCollection(allPoints[hull.simplices], facecolors=color, linewidths=1, edgecolors=None, alpha=alpha))

            # Plot bottom and top circles
            # ax.plot(bottomCirclePoints[:, 0], bottomCirclePoints[:, 1], bottomCirclePoints[:, 2], 'b-')
            # ax.plot(topCirclePoints[:, 0], topCirclePoints[:, 1], topCirclePoints[:, 2], 'b-')
        except:
            pass

    def collidesWith(self, other, includeEnds=False):
        #if too far away to collide, return false
        if np.linalg.norm(self.center - other.center) > self.halfDist + other.halfDist:
            return False, None
        
        result, point = capsuleCollision(self, other)

        return result, point

def pointInCapsuleHemisphere(point, capsule):
    start, end, radius = capsule
    axis = end - start
    axis_length = np.linalg.norm(axis)

    # Check distance from point to start and end
    distToStart = np.linalg.norm(point - start)
    distToEnd = np.linalg.norm(point - end)

    # Check if point is closer to an end than the cylinder part
    if distToStart < distToEnd and distToStart < radius:
        return True  # In start hemisphere
    elif distToEnd < distToStart and distToEnd < radius:
        return True  # In end hemisphere
    
    return False  # Not in either hemisphere

def closest_point_on_line_segment(a, b, p):
    ab = b - a
    t = np.dot(p - a, ab) / np.dot(ab, ab)
    return a + np.clip(t, 0, 1) * ab

    
def capsuleCollision(capsule1, capsule2):
    # find closest points between the line segments of both capsules
    p1, p2 = closestPointsBetweenLineSegments(capsule1.start, capsule1.end, capsule2.start, capsule2.end)
    
    vector = p2 - p1
    distance = np.linalg.norm(vector)
    
    sumRadii = capsule1.radius + capsule2.radius

    if distance < sumRadii:
        collisionPoint = p1 + vector * 0.5
        
        if distance > sumRadii/10000 and (pointInCapsuleHemisphere(collisionPoint, (capsule1.start, capsule1.end, capsule1.radius)) or pointInCapsuleHemisphere(collisionPoint, (capsule2.start, capsule2.end, capsule2.radius))):
            if not separatingAxisTheorem(CollisionBox(startDubinsFrame=capsule1.base, endDubinsFrame=capsule1.otherBase, r = capsule1.radius/(2**0.5) - capsule1.radius/1000), CollisionBox(startDubinsFrame=capsule2.base, endDubinsFrame=capsule2.otherBase, r = capsule2.radius/(2**0.5) - capsule2.radius/1000)):
                return False, None

        return True, collisionPoint
    else:
        return False, None

def closestPointsBetweenLineSegments(a1, a2, b1, b2):
    """
    Find the closest points between two line segments, preferring solutions closer to midpoints.
    
    :param a1: Start point of the first line segment
    :param a2: End point of the first line segment
    :param b1: Start point of the second line segment
    :param b2: End point of the second line segment
    :return: Tuple of closest points (point_on_segment_a, point_on_segment_b)
    """
    # Convert inputs to numpy arrays
    a1, a2, b1, b2 = map(np.array, (a1, a2, b1, b2))
    
    # Vectors
    d1 = a2 - a1
    d2 = b2 - b1
    r = a1 - b1
    
    # Midpoints
    mid_a = (a1 + a2) / 2
    mid_b = (b1 + b2) / 2
    
    a = np.dot(d1, d1)
    e = np.dot(d2, d2)
    f = np.dot(d2, r)
    
    # check if segments degenerate into points
    if a <= 1e-8 and e <= 1e-8:
        # both segments degenerate into points
        return a1, b1
    
    if a <= 1e-8:
        # first segment degenerates into a point
        s = 0
        t = np.clip(f / e, 0, 1)
    elif e <= 1e-8:
        # second segment degenerates into a point
        t = 0
        s = np.clip(-np.dot(d1, r) / a, 0, 1)
    else:
        # non-degenerate case
        c = np.dot(d1, d2)
        b = np.dot(d1, r)
        denom = a * e - c * c
        
        if abs(denom) > 1e-8:
            s = np.clip((b * e - c * f) / denom, 0, 1)
        else:
            #lines are parallel
            s_numer = np.dot(b1 - a1, a2 - a1)
            s_denom = np.dot(a2 - a1, a2 - a1)
            s = np.clip(s_numer / s_denom, 0, 1)
        
        t = np.clip((c * s + f) / e, 0, 1)
        
        # check multiple solutions
        if abs(np.dot(d1, d2)) > 1 - 1e-8:  # nearly parallel
            s_start = max(0, min(1, s, (c - f) / e))
            s_end = min(1, max(0, s, (c + f) / e))
            
            # choose s closest to midpoint
            if abs(s_start - 0.5) < abs(s_end - 0.5):
                s = s_start
            else:
                s = s_end
            
            t = np.clip((c * s + f) / e, 0, 1)
    
    point_on_segment_a = a1 + s * d1
    point_on_segment_b = b1 + t * d2
    
    return point_on_segment_a, point_on_segment_b