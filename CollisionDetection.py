import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
from spatialmath import SE3

class CollisionBox:
    def __init__(self, points=None, startDubinsFrame=None, endDubinsFrame=None,r=None):
        if points is None:
            if startDubinsFrame is None or endDubinsFrame is None or r is None:
                points = np.zeros((8, 3))
            else:
                points = np.array([
                    startDubinsFrame.t + r*startDubinsFrame.o + r*startDubinsFrame.a,
                    startDubinsFrame.t - r*startDubinsFrame.o + r*startDubinsFrame.a,
                    startDubinsFrame.t - r*startDubinsFrame.o - r*startDubinsFrame.a,
                    startDubinsFrame.t + r*startDubinsFrame.o - r*startDubinsFrame.a,
                    endDubinsFrame.t + r*endDubinsFrame.o - r*endDubinsFrame.a,
                    endDubinsFrame.t - r*endDubinsFrame.o - r*endDubinsFrame.a,
                    endDubinsFrame.t - r*endDubinsFrame.o + r*endDubinsFrame.a,
                    endDubinsFrame.t + r*endDubinsFrame.o + r*endDubinsFrame.a,
                ])
    
                # Compute the covariance matrix
                covariance_matrix = np.cov(points, rowvar=False)
                
                # Perform eigen decomposition to get principal components
                eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)
                
                # Sort the eigenvectors by eigenvalues in descending order
                idx = np.argsort(eigenvalues)[::-1]
                eigenvectors = eigenvectors[:, idx]
                
                # Project the points onto the principal components
                projected_points = points @ eigenvectors
                
                # Get the min and max points in the new basis
                min_point = np.min(projected_points, axis=0)
                max_point = np.max(projected_points, axis=0)
                
                # Define the 8 corners of the bounding box in the new basis
                box_corners = np.array([
                    [min_point[0], min_point[1], min_point[2]],
                    [min_point[0], min_point[1], max_point[2]],
                    [min_point[0], max_point[1], min_point[2]],
                    [min_point[0], max_point[1], max_point[2]],
                    [max_point[0], min_point[1], min_point[2]],
                    [max_point[0], min_point[1], max_point[2]],
                    [max_point[0], max_point[1], min_point[2]],
                    [max_point[0], max_point[1], max_point[2]]
                ])
                
                # Transform the box corners back to the original coordinate system
                box_corners = box_corners @ eigenvectors.T

                self.points= box_corners
        else:
            self.points = np.array(points)
    
    def _get_bounds(self):
        min_bound = self.points.min(axis=0)
        max_bound = self.points.max(axis=0)
        return min_bound, max_bound

    def expand_to_envelop(self, other):
        """Expand this box to envelop another box."""
        min_bound_self, max_bound_self = self._get_bounds()
        min_bound_other, max_bound_other = other._get_bounds()

        new_min_bound = np.minimum(min_bound_self, min_bound_other)
        new_max_bound = np.maximum(max_bound_self, max_bound_other)

        # Create new points based on the new bounds
        self.points = np.array([
            [new_min_bound[0], new_min_bound[1], new_min_bound[2]],  # Bottom-front-left
            [new_max_bound[0], new_min_bound[1], new_min_bound[2]],  # Bottom-front-right
            [new_max_bound[0], new_max_bound[1], new_min_bound[2]],  # Bottom-back-right
            [new_min_bound[0], new_max_bound[1], new_min_bound[2]],  # Bottom-back-left
            [new_min_bound[0], new_min_bound[1], new_max_bound[2]],  # Top-front-left
            [new_max_bound[0], new_min_bound[1], new_max_bound[2]],  # Top-front-right
            [new_max_bound[0], new_max_bound[1], new_max_bound[2]],  # Top-back-right
            [new_min_bound[0], new_max_bound[1], new_max_bound[2]]   # Top-back-left
        ])
    
    def addToPlot(self, ax, color='r', alpha=0.4):
        try:
            hull = ConvexHull(self.points)

            # Get the vertices of the convex hull
            vertices = hull.points[hull.vertices]
            
            # Collect the faces from the convex hull
            faces = hull.simplices
            
            # Create the 3D polygon collection
            poly3d = [[vertices[vidx] for vidx in face] for face in faces]
            box_faces = Poly3DCollection(poly3d, facecolors=color, linewidths=1, edgecolors=None, alpha=alpha)
            
            # Add the box to the plot
            ax.add_collection3d(box_faces)
        except:
            pass

def separating_axis_theorem(box1, box2):
    box1 = box1.points
    box2 = box2.points
    # Function to project a box onto an axis and return the min and max projections
    def project_box(box, axis):
        projections = [np.dot(vertex, axis) for vertex in box]
        return min(projections), max(projections)
    
    # Function to check if projections on a given axis overlap
    def overlap(proj1, proj2):
        return proj1[0] <= proj2[1] and proj2[0] <= proj1[1]
    
    # Extract edges and face normals
    def get_edges_and_normals(box):
        edges = []
        normals = []
        for i in range(4):
            edge1 = box[i+4] - box[i]
            edge2 = box[(i+1)%4 + 4] - box[i+4]
            edge3 = box[i] - box[(i+1)%4]
            edges.append(edge1)
            edges.append(edge2)
            edges.append(edge3)
            normals.append(np.cross(edge1, edge2))
        return edges, normals
    
    # Get edges and normals for both boxes
    edges1, normals1 = get_edges_and_normals(box1)
    edges2, normals2 = get_edges_and_normals(box2)
    
    # List of potential separating axes
    axes = normals1 + normals2 + [np.cross(e1, e2) for e1 in edges1 for e2 in edges2]
    
    # Check for a separating axis
    for axis in axes:
        if np.linalg.norm(axis) == 0:
            continue
        axis = axis / np.linalg.norm(axis)
        projection1 = project_box(box1, axis)
        projection2 = project_box(box2, axis)
        if not overlap(projection1, projection2):
            return False  # Found a separating axis, hence boxes do not intersect
    
    return True  # No separating axis found, hence boxes intersect


class CollisionCapsule:
    def __init__(self, base : SE3(), radius : float, height : float):
        self.radius = radius
        self.height = height
        self.base = base @ SE3.Ry(np.pi/2)
        self.otherBase = self.base @ SE3.Trans(self.height * self.base.n)
        self.center = self.base.t + self.base.n * (self.height/2 + np.sign(self.height)*self.radius)
        self.halfDist = np.abs(self.height/2) + self.radius

    
    def addToPlot(self, ax, color='r', alpha=0.25):
        def create_circle_points(center, radius, normal, num_points=50):
            """Generate points on a circle in 3D."""
            angles = np.linspace(0, 2 * np.pi, num_points)
            circle_points = []
            for angle in angles:
                point = center + radius * (np.cos(angle) * orthogonal_vector(normal) + np.sin(angle) * np.cross(normal, orthogonal_vector(normal)))
                circle_points.append(point)
            return np.array(circle_points)
        
        def orthogonal_vector(v):
            """Find an orthogonal vector to v."""
            if np.isclose(v[0], 0) and np.isclose(v[1], 0):
                if np.isclose(v[2], 0):
                    raise ValueError("Zero vector")
                return np.array([0, 1, 0])
            return np.array([-v[1], v[0], 0])

        try:
            # Bottom and top circles centers
            bottom_center = self.base.t
            top_center = bottom_center + self.height * self.base.R[:,2]

            # Generate points on the bottom and top circles
            bottom_circle_points = create_circle_points(bottom_center, self.radius, self.base.R[:,2])
            top_circle_points = create_circle_points(top_center, self.radius, self.base.R[:,2])

            # Create the spherical end caps
            theta, phi = np.mgrid[0:np.pi:10j, 0:2*np.pi:10j]
            x = self.radius * np.sin(theta) * np.cos(phi)
            y = self.radius * np.sin(theta) * np.sin(phi)
            z = self.radius * np.cos(theta)

             # Bottom spherical cap
            bottom_cap = np.vstack((x.flatten(), y.flatten(), z.flatten())).T
            bottom_cap = bottom_cap @ self.base.R.T + bottom_center

            # Top spherical cap
            top_cap = np.vstack((x.flatten(), y.flatten(), -z.flatten())).T
            top_cap = top_cap @ self.base.R.T + top_center

            # Combine all points
            all_points = np.vstack((bottom_circle_points, top_circle_points, bottom_cap, top_cap))
            hull = ConvexHull(all_points)

            # Plot the capsule
            # for simplex in hull.simplices:
            #     ax.plot(all_points[simplex, 0], all_points[simplex, 1], all_points[simplex, 2], 'k-')
            ax.add_collection3d(Poly3DCollection(all_points[hull.simplices], facecolors=color, linewidths=1, edgecolors=None, alpha=alpha))

            # Plot bottom and top circles
            # ax.plot(bottom_circle_points[:, 0], bottom_circle_points[:, 1], bottom_circle_points[:, 2], 'b-')
            # ax.plot(top_circle_points[:, 0], top_circle_points[:, 1], top_circle_points[:, 2], 'b-')
        except:
            pass
    
    def decomposeCapsule(self):
        #segment: Array with two points representing the segment.
        #end_centers: Array with two points representing the centers of the spherical ends.

        bottom_center = self.base.t
        top_center = bottom_center + self.height * self.base.R[:, 2]
        return np.array([bottom_center, top_center]), np.array([bottom_center, top_center])

    def frameOverlap(self, frame, radius):
        return False, None
        #return self.collidesWith(CollisionCapsule(frame, radius, 0))
    
    def collidesWith(self, other):
        #if too far away to collide, return false
        if np.linalg.norm(self.center - other.center) > self.halfDist + other.halfDist:
            return False, None

        segment1, ends1 = self.decomposeCapsule()
        segment2, ends2 = other.decomposeCapsule()

        minDistance = float('inf')
        collisionPoint = None
        
        # Check distance between cylindrical segments
        for point in segment1:
            distance, closestPoint = pointToSegmentDistance(point, segment2)

            pt = 0.5 * (point + closestPoint)
            #if point is on the same side of both circular bases of the capsule, then contact point is not in cylinder
            if distance <= minDistance and signedDistanceToFrame(pt, self.base) * signedDistanceToFrame(pt, self.otherBase) <= 0 and signedDistanceToFrame(pt, other.base) * signedDistanceToFrame(pt, other.otherBase) <= 0:
                minDistance = distance
                collisionPoint = pt

        for point in segment2:
            distance, closestPoint = pointToSegmentDistance(point, segment1)

            pt = 0.5 * (point + closestPoint)
            #if point is on the same side of both circular bases of the capsule, then contact point is not in cylinder
            if distance <= minDistance and signedDistanceToFrame(pt, self.base) * signedDistanceToFrame(pt, self.otherBase) <= 0 and signedDistanceToFrame(pt, other.base) * signedDistanceToFrame(pt, other.otherBase) <= 0:
                minDistance = distance
                collisionPoint = pt
                

        # Check distance between spherical ends
        # for p1 in ends1:
        #     for p2 in ends2:
        #         distance = np.linalg.norm(p1 - p2)
        #         if distance <= minDistance:
        #             minDistance = distance
        #             collisionPoint = 0.5 * (p1 + p2)

        if minDistance <= self.radius + other.radius:
            return True, collisionPoint
        else:
            return False, None

def pointToSegmentDistance(point, segment):
    p = point
    v = segment[0]
    w = segment[1]
    
    l2 = np.sum((v - w) ** 2)
    if l2 == 0:
        return np.linalg.norm(p - v), point
    
    t = max(0, min(1, np.dot(p - v, w - v) / l2))
    projection = v + t * (w - v)

    return np.linalg.norm(p - projection), projection

def signedDistanceToFrame(point, frame : SE3()):
    vector = point - frame.t
    dist = np.dot(vector, frame.n)

    return dist