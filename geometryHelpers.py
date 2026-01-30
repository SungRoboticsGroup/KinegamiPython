"""
Assorted geometry-related helper functions and classes
"""
from __future__ import annotations
import numpy as np
from numpy import cross, dot, arctan2
import scipy
from scipy.spatial.transform import Rotation
from scipy.linalg import null_space
from scipy.optimize import minimize, NonlinearConstraint
from numpy.linalg import norm
from spatialmath import SO3, SE3
import matplotlib.pyplot as plt
import math
from math import remainder
from style import *
from matplotlib import cm
from typing import Optional, Union
from numpy.typing import ArrayLike, NDArray
from types import ModuleType
import manifold3d as m3d
import trimesh
from pytetwild.pytetwild import tetrahedralize
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d import Axes3D

def unit(v):
    return v / np.linalg.norm(v)


"""
Output n points spaced evenly from start to end.
"""
def lerpN(start, end, n):
    assert(start.shape == end.shape)
    d = start.shape[0]
    if n==0:
        return np.empty((0,d))
    elif n==1:
        return np.array([start + (end-start)/2])
    else:
        ts = (np.arange(n)/(n-1)).reshape((-1,1))
        return start + ts*(end-start)


"""
Given a numpy array whose rows are points, and an array of the same points
cycled by a roll of k, find k. If this fails (i.e. if the points don't
actually match via a roll) return None.
"""
def rollToMatch(points, pointsCycled, EPSILON=0.00001):
    n = points.shape[0]
    assert(points.shape == pointsCycled.shape)
    for k in range(n):
        if np.all(np.abs(np.roll(points, k, axis=0) - pointsCycled) < EPSILON):
            return k
    return None


""" 
Given 3D np vectors u and v, return a unit vector orthogonal to both, 
obeying the right-hand-rule if applicable 
(i.e., if neither is 0 and they aren't colinear).
"""
def unitNormalToBoth(u, v):
    # ensure correct dimensions
    u = u.flatten()
    v = v.flatten()
    assert(u.shape==(3,))
    assert(v.shape==(3,))
    
    cp = cross(u,v)
    # if u and v are not colinear and both nonzero
    if norm(cp) > 0:
        return cp / norm(cp)
    else: 
        # return any vector orthogonal to both
        nullSpaceBasis = null_space(np.array([u, v]))
        # any nullspace basis vector will work, 
        # we'll use the first one since it may be the only one
        return nullSpaceBasis[:,0].flatten()


"""
For A and B matrices (of the same shape) storing vectors as rows,
return the unsigned angle between corresponding rows
"""
def unsignedAngles(A,B):
    assert(A.shape == B.shape)
    AnormalizedRows = (A.T * (1/norm(A, axis=1))).T
    BnormalizedRows = (B.T * (1/norm(B, axis=1))).T
    dotProductsOfRows = np.sum(AnormalizedRows*BnormalizedRows, axis=1)
    return np.arccos(np.clip(dotProductsOfRows, -1.0, 1.0))

"""
Unsigned angle between two vectors, based on
https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
"""
def unsignedAngle(a,b):
    return np.arccos(np.clip(np.dot(a/norm(a), b/norm(b)), -1.0, 1.0))
    
    
"""
Signed angle (radians) from vector a to vector b, around the normal vector n.
All inputs should be numpy arrays of shape (3,)
"""
def signedAngle(a, b, n):
    a = a / norm(a)
    b = b / norm(b)
    if norm(n)>0:
        n = n / norm(n)
        
    return arctan2(dot(cross(a,b),n), dot(a,b))

"""
For A and B matrices (of the same shape) storing 2D vectors as rows,
return the signed angle between corresponding rows
"""
def signedAngles2D(A, B):
    assert(A.shape == B.shape)
    k = A.shape[0]
    A3 = np.hstack((A, np.zeros((k,1))))
    B3 = np.hstack((B, np.zeros((k,1))))
    return np.array([signedAngle(A3[i,:], B3[i,:], [0,0,1])
                     for i in range(k)])

def signedAngle2D(a, b):
    return signedAngle(np.hstack((a,[0])), np.hstack((b,[0])), [0,0,1])

"""
For A and B matrices (of the same shape) storing 3D vectors as rows,
return the signed angle about n between corresponding rows
"""
def signedAngles3D(A, B, n):
    assert(A.shape == B.shape)
    k = A.shape[0]
    return np.array([signedAngle(A[i], B[i], n) for i in range(k)]).flatten()

# wrap angles to [0,2pi)
def wrapAngle(angle, EPSILON=0.00001):
    modResult = angle % (2*np.pi)
    if abs(modResult-2*np.pi) < EPSILON:
        return 0
    return modResult


class Line:
    def __init__(self, point, direction, EPSILON=1e-8):
        self.p = np.array(point)
        direction = np.array(direction)
        self.EPSILON = EPSILON
        assert(norm(direction) > self.EPSILON)
        self.dhat = direction / norm(direction)

    def distanceToPoint(self, point):
        point = np.array(point)
        return norm(cross(point - self.p, self.dhat))
    
    def contains(self, point):
        return self.distanceToPoint(point) < self.EPSILON
    
    def projectionOfPoint(self, point):
        return self.p + dot(point - self.p, self.dhat) * self.dhat
    
class Ray:
    def __init__(self, startPoint, direction, EPSILON=1e-8):
        self.startPoint = np.array(startPoint)
        direction = np.array(direction)
        self.EPSILON = EPSILON
        assert(norm(direction) > self.EPSILON)
        self.dhat = direction / norm(direction)

    def distanceToPoint(self, point): 
        if Plane(self.startPoint, self.dhat).sideOfPoint(point) < 0:
            return norm(point - self.startPoint)
        else:
            return Line(self.startPoint, self.dhat).distanceToPoint(point)
    
    def contains(self, point):
        return self.distanceToPoint(point) < self.EPSILON



def unitSphereParameterization(theta, phi):
    return np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])

def directionNotOrthogonalToAnyOf(otherDirections: np.ndarray) -> np.ndarray:
    assert(otherDirections.shape[1] == 3)
    assert(np.min(norm(otherDirections,axis=1)) > 0)
    otherDirections = otherDirections / (norm(otherDirections,axis=1).reshape(-1,1)) #normalize

    # We want a direction not orthogonal to any of the given directions,
    # i.e., with nonzero dot product with all of them.
    # Such a direction must exist since there are finitely many given directions
    # so let's take the one as far from orthogonal as possible, i.e.,
    # maximizing the minimum absolute value of the dot product with 
    # the given directions.
    def negativeMinAbsDotProduct(v):
        return -np.min(np.abs(otherDirections @ v))
    
    if otherDirections.shape[0] == 0:
        return np.array([1,0,0])
    else:
        C = NonlinearConstraint(lambda v : norm(v), 1, 1)
        result = minimize(negativeMinAbsDotProduct, np.random.rand(3), constraints=(C,))
        return result.x
        

class Plane:
    def __init__(self, point, normal, EPSILON=0.00000001):
        self.p = point
        self.nhat = normal / norm(normal)
        self.EPSILON = EPSILON
        
    # sign is + if the point is on the +nhat side, 0 if on plane, - otherwise
    def signedDistanceToPoint(self, point):
        projectionLength = dot(point - self.p, self.nhat)
        if abs(projectionLength) < self.EPSILON: #for numerical stability
            return 0
        else:
            return projectionLength
    
    def projectionOfPoint(self, point):
        return point - self.signedDistanceToPoint(point)*self.nhat

    # return 0 if point is on the plane, +1 if on the +nhat side, -1 otherwise
    def sideOfPoint(self, point):
        return np.sign(self.signedDistanceToPoint(point))
    
    def containsPoint(self, point):
        # plane contains point iff self.p - point is orthogonal to self.nhat
        return self.sideOfPoint(point) == 0
    
    def intersectionWithLine(self, line):
        if abs(dot(self.nhat, line.dhat)) < self.EPSILON: 
            # line is either on the plane or parallel to it
            if self.containsPoint(line.p): 
                # line is on the plane
                return line
            else: 
                # line is parallel to the plane
                return None
        else:
            # line intersects the plane at a single point
            # so the linear system constructed per wikipedia has point solution
            # https://en.wikipedia.org/wiki/Line%E2%80%93plane_intersection
            # what they call cross(p01, p02) is just our nhat
            t = dot(self.nhat, line.p - self.p) / dot(-line.dhat, self.nhat)
            intersect = line.p + t * line.dhat
            assert(self.containsPoint(intersect))
            return intersect
    
    def intersectionsWithParallelLines(self, Points : np.ndarray, 
                                       direction : np.ndarray) -> np.ndarray:
        PtoPoints = Points - self.p.flatten()
        alignment = dot(-direction, self.nhat.flatten())
        ts =  PtoPoints @ self.nhat.reshape(3,1) / alignment
        return Points + ts*direction
    
    """ 
    Return an array containing the sides of the plane that the given line
    intersects with: 0 for on the plane, 1 for +nhat side, -1 for -nhat side
    """
    def sidesOfLine(self, line):
        intersection = self.intersectionWithLine(line)
        if intersection is None:
            # line is parallel to plane, need to check on which side
            return [self.sideOfPoint(line.p)]
        elif type(intersection) == Line:
            # whole line is on plane
            return [0]
        else:
            # line crosses the plane at a point
            return [-1,0,1]
    
    def parallelPlane(self, distance):
        return Plane(self.p + distance*self.nhat, self.nhat, self.EPSILON)

    def grid(self, scale=20, numPoints=9):
        range = np.linspace(-scale, scale, numPoints)
        basis = null_space([self.nhat])
        uhat = basis[:,0]
        vhat = basis[:,1]
        grid = np.zeros((numPoints, numPoints, 3))
        for u in np.arange(numPoints):
            for v in np.arange(numPoints):
                grid[u,v] = self.p + range[u]*uhat + range[v]*vhat
        return grid

    def addToPlot(self, ax, color='red', alpha=0.5, scale=20):
        grid = self.grid(scale, numPoints=9)
        X = grid[:,:,0]
        Y = grid[:,:,1]
        Z = grid[:,:,2]
        ax.plot_surface(X, Y, Z, color=color, alpha=alpha)

def planeFromThreePoints(p1, p2, p3):
    normal = unit(cross(p2-p1, p3-p1))
    return Plane(p1, normal)

class Circle3D:
    def __init__(self, radius, center, normal, radialVector=None):
        assert(norm(normal)>0)
        self.r = radius
        self.c = center
        self.n = normal / norm(normal)
        if radialVector is None:
            self.radialVector = unitNormalToBoth(self.n, self.n).reshape(1,3)
        else:
            self.radialVector = radialVector / norm(radialVector)

    def interpolate(self, count=50):
        angle = np.linspace(0, 2*np.pi, count).reshape(-1,1)
        u = self.r * np.cos(angle)
        v = self.r * np.sin(angle)
        
        # construct basis for circle plane
        uhat = self.radialVector.reshape(1,3)
        vhat = cross(self.n, uhat).reshape(1,3)
        
        # 3d circle points
        return self.c + u @ uhat + v @ vhat


def line_sphere_intersection(line_point: np.ndarray, line_dir: np.ndarray, 
                             sphere_center: np.ndarray, sphere_radius: float,
                             epsilon: float = 1e-8) -> Optional[tuple]:
    """
    Find intersection of a line with a sphere.
    
    Args:
        line_point: A point on the line
        line_dir: Direction vector of the line (should be normalized)
        sphere_center: Center of the sphere
        sphere_radius: Radius of the sphere
        epsilon: Numerical tolerance
        
    Returns:
        Tuple (t1, t2) of parameter values where line intersects sphere,
        where line is parameterized as p(t) = line_point + t * line_dir.
        Returns None if no intersection.
        For tangent intersection, t1 == t2.
        
    Reference:
        https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection
    """
    oc = line_point - sphere_center
    
    a = np.dot(line_dir, line_dir)  # Should be 1 if normalized
    b = 2.0 * np.dot(oc, line_dir)
    c = np.dot(oc, oc) - sphere_radius**2
    
    discriminant = b**2 - 4*a*c
    
    if discriminant < -epsilon:
        return None  # No intersection
    elif abs(discriminant) < epsilon:
        # One intersection (tangent)
        t = -b / (2*a)
        return (t, t)
    else:
        # Two intersections
        sqrt_disc = np.sqrt(discriminant)
        t1 = (-b - sqrt_disc) / (2*a)
        t2 = (-b + sqrt_disc) / (2*a)
        return (min(t1, t2), max(t1, t2))


def discs_cross(disc1: Circle3D, disc2: Circle3D, epsilon: float = 1e-8) -> bool:
    """
    Check if two 3D discs (filled circles) cross through each other
    (i.e., intersect along a line segment of nonzero length).
    
    Args:
        disc1: First disc
        disc2: Second disc
        epsilon: Numerical tolerance
        
    Returns:
        True if discs cross (penetrate), False if separated or just tangent
    """
    # Step 1: Find intersection of the two planes
    plane1 = Plane(disc1.c, disc1.n, EPSILON=epsilon)
    plane2 = Plane(disc2.c, disc2.n, EPSILON=epsilon)
    
    # Check if planes are parallel (or coincident)
    if abs(abs(np.dot(disc1.n, disc2.n)) - 1.0) < epsilon:
        # Planes are parallel or anti-parallel
        return False # If they are coplanar, the discs can press together but not cross
        """
        if plane1.containsPoint(disc2.c):
            # Coplanar: check distance between centers (strictly less than sum of radii to cross)
            center_dist = norm(disc1.c - disc2.c)
            return center_dist < disc1.r + disc2.r - epsilon
        else:
            # Parallel but not coplanar: no intersection
            return False
        """
    
    # Planes intersect in a line
    # Find line direction (perpendicular to both normals)
    line_dir = cross(disc1.n, disc2.n)
    line_dir = line_dir / norm(line_dir)
    
    # Find a point on the line of intersection
    # We solve for a point that lies on both planes
    # Use the method from Plane.intersectionWithLine but in reverse
    # Pick a point on plane1, project to find point on both planes
    # Actually, we can use a more direct approach:
    
    # The line of intersection lies on both planes, so we can find it by
    # solving the system of equations. A simple approach:
    # Find a point on both planes by choosing a convenient coordinate
    n1, n2 = disc1.n, disc2.n
    c1, c2 = disc1.c, disc2.c
    
    # Find largest component of line_dir to avoid division by small numbers
    max_idx = np.argmax(np.abs(line_dir))
    if max_idx == 0:
        # Fix x=0, solve for y,z
        # n1┬╖(p-c1) = 0 and n2┬╖(p-c2) = 0 with p_x = 0
        # This gives us two equations in two unknowns (y, z)
        A = np.array([[n1[1], n1[2]], [n2[1], n2[2]]])
        b = np.array([np.dot(n1, c1), np.dot(n2, c2)])
        if abs(np.linalg.det(A)) > epsilon:
            yz = np.linalg.solve(A, b)
            line_point = np.array([0, yz[0], yz[1]])
        else:
            return False
    elif max_idx == 1:
        # Fix y=0
        A = np.array([[n1[0], n1[2]], [n2[0], n2[2]]])
        b = np.array([np.dot(n1, c1), np.dot(n2, c2)])
        if abs(np.linalg.det(A)) > epsilon:
            xz = np.linalg.solve(A, b)
            line_point = np.array([xz[0], 0, xz[1]])
        else:
            return False
    else:
        # Fix z=0
        A = np.array([[n1[0], n1[1]], [n2[0], n2[1]]])
        b = np.array([np.dot(n1, c1), np.dot(n2, c2)])
        if abs(np.linalg.det(A)) > epsilon:
            xy = np.linalg.solve(A, b)
            line_point = np.array([xy[0], xy[1], 0])
        else:
            return False
    
    # Step 2: Find intersection of line with each sphere
    seg1 = line_sphere_intersection(line_point, line_dir, disc1.c, disc1.r, epsilon)
    seg2 = line_sphere_intersection(line_point, line_dir, disc2.c, disc2.r, epsilon)
    
    if seg1 is None or seg2 is None:
        return False
    
    # Step 3: Check if line segments overlap (strictly, not just tangent)
    # Segments are [seg1[0], seg1[1]] and [seg2[0], seg2[1]]
    # They cross if: max(seg1[0], seg2[0]) < min(seg1[1], seg2[1])
    overlap_start = max(seg1[0], seg2[0])
    overlap_end = min(seg1[1], seg2[1])
    
    return overlap_start < overlap_end - epsilon


class Ball:
    # closed ball centered at self.c of radius self.r
    def __init__(self, center, radius):
        self.c = center
        self.r = radius

    def __repr__(self):
        return (
            "Ball("
            f"center={repr(self.c)},"
            f"radius={repr(self.r)})"
        )
    
    def containsPoint(self, point):
        return norm(self.c - point) <= self.r
    
    def addToPlot(self, ax, color='black', alpha=0.1, frame=False):
        #https://www.tutorialspoint.com/plotting-a-3d-cube-a-sphere-and-a-vector-in-matplotlib
        u, v = np.mgrid[0:2*np.pi:40j, 0:np.pi:20j]
        x = self.c[0] + self.r*np.cos(u)*np.sin(v)
        y = self.c[1] + self.r*np.sin(u)*np.sin(v)
        z = self.c[2] + self.r*np.cos(v)
        if frame:
            return ax.plot_wireframe(x, y, z, color=color, alpha=alpha)
        else:
            return ax.plot_surface(x, y, z, color=color, alpha=alpha)
    
    def addToWidget(self, widget, color=ballDefaultColor, is_waypoint=False):
        import pyqtgraph.opengl as gl
        md = gl.MeshData.sphere(rows=10, cols=10)
        sphere = gl.GLMeshItem(meshdata=md, color=tuple(color), shader='shaded', smooth=True)
        sphere.setGLOptions('translucent')
        sphere.scale(self.r, self.r, self.r)
        sphere.translate(*self.c)
        if (is_waypoint):
            sphere.setObjectName("Waypoint")
            sphere.setGLOptions('opaque')
            sphere.scale(self.r, self.r, self.r)
        widget.plot_widget.addItem(sphere)
    
    def show(self, color='black', alpha=1, frame=False, block=blockDefault):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandles = self.addToPlot(ax, color, alpha, frame)
        ax.set_aspect('equal')
        plt.show(block=block)

    def projectionOntoPlane(self, plane : Plane) -> Circle3D:
        return Circle3D(self.r, plane.projectionOfPoint(self.c), plane.nhat)
    
    def translationToCenterOnPlane(self, plane : Plane):
        return Ball(plane.projectionOfPoint(self.c), self.r)
    
    def expandToCenterOnLine(self, line : Line):
        self.c, self.r = line.projectionOfPoint(self.c), self.r + line.distanceToPoint(self.c)

    # returns whether the balls are tangent (externally, internally, or coincident)
    def isTangentToBall(self, otherBall : Ball, epsilon=1e-8) -> bool:
        return abs(norm(self.c - otherBall.c) - abs(self.r - otherBall.r)) < epsilon

    def newBallTransformedBy(self, T : SE3) -> Ball:
        return Ball(T * self.c, self.r)
    
    def containsBall(self, otherBall : Ball, epsilon=1e-8) -> bool:
        return norm(self.c - otherBall.c) + otherBall.r <= self.r + epsilon

        

"""
Note: we don't need to guarantee minimality of our bounding balls, so we build
bounding balls of bounding balls in a greedy fashion based on the below 
function which takes the minimum bounding ball of 2 balls. The greedy approach
seems to give bounding balls reasonably close to what would be minimal in our 
application anyway. 

Taking the minimum bounding ball of n balls is highly nontrivial, see:
Fischer, Kaspar, and Bernd Gartner. "The smallest enclosing ball of balls: 
combinatorial structure and algorithms." Proceedings of the nineteenth annual 
symposium on Computational geometry. 2003.
"""
# minimum bounding ball of 2 balls
def minBoundingBall(ball1, ball2):
    if ball1.r <= ball2.r:
        smaller, larger = ball1, ball2
    else:
        smaller, larger = ball2, ball1
    
    v = larger.c - smaller.c
    if norm(v) + smaller.r <= larger.r: # if larger contains smaller
        return larger
    else:
        vhat = v / norm(v)
        p = smaller.c - smaller.r*vhat
        q = larger.c + larger.r*vhat
        pq = q-p
        return Ball(p + pq/2, norm(pq)/2)

def distanceBetweenBalls(ball1, ball2):
    centerDistance = norm(ball2.c - ball1.c)
    return max(0, centerDistance - ball1.r - ball2.r)

class Cylinder:
    def __init__(self, radius : float, start : np.ndarray, 
                 direction : np.ndarray, length : float, 
                 uhat : np.ndarray = None):
        assert(length > 0)
        assert(norm(direction)>0)
        self.start = start
        self.direction = direction / norm(direction)
        self.length = length
        self.r = radius
        if uhat is None:
            uhat = null_space([self.direction])[:,0]
        else:
            uhat = uhat.reshape((3))    
        self.uhat = uhat / norm(uhat)
    
    def orientation(self) -> SO3:
        return SO3(np.vstack((self.uhat, 
                              cross(self.direction, self.uhat), 
                              self.direction)).T)
    
    def end(self):
        return self.start + self.length * self.direction

    def startPlane(self) -> Plane:
        return Plane(self.start, self.direction)
    
    def endPlane(self) -> Plane:
        return Plane(self.end(), self.direction)

    # Keeping self.direction constant, expand cylinder to include ball
    def expandToIncludeBall(self, ball : Ball):
        ballStart = ball.c - ball.r * self.direction
        ballEnd = ball.c + ball.r * self.direction
        distanceInBack = self.startPlane().signedDistanceToPoint(ballStart)
        distanceInFront = self.endPlane().signedDistanceToPoint(ballEnd)

        # Update cylinder forward/backward
        if distanceInFront > 0:
            self.length += distanceInFront
        if distanceInBack < 0:
            self.start = self.start + (distanceInBack * self.direction)
            self.length -= distanceInBack
            
        # Update cylinder circular cross-section
        ballTranslatedToStartPlane = ball.translationToCenterOnPlane(self.startPlane())
        cylinderStartBall = Ball(self.start, self.r)
        expandedStartBall = minBoundingBall(cylinderStartBall, ballTranslatedToStartPlane)
        self.start = expandedStartBall.c
        self.r = expandedStartBall.r
    
    def interpolateCircles(self, numPointsPerCircle=32, numCircles=2):
        radialCount = numPointsPerCircle+1 #the +1 is because the first equals the last
        angle = np.linspace(0, 2*np.pi, radialCount) 
        u = self.r * np.cos(angle)
        v = self.r * np.sin(angle)
        """
        circlePlaneBasis = null_space([self.direction])
        uhat = circlePlaneBasis[:,0]
        vhat = circlePlaneBasis[:,1]
        """
        uhat = self.uhat
        vhat = cross(self.direction, uhat)
        circle = u.reshape(-1,1) @ uhat.reshape(1,3) + v.reshape(-1,1) @ vhat.reshape(1,3)
        
        segment = np.linspace(self.start, self.end(), numCircles)
        circlePoints = np.tile(circle, (numCircles,1)) + np.repeat(segment, radialCount, axis=0)
        return circlePoints.reshape((numCircles, radialCount, 3))
    
    def interpolateQtCircles(self, numPointsPerCircle=32, numCircles=10):
        angle = np.linspace(0, 2 * np.pi, numPointsPerCircle, endpoint=False)
        u = self.r * np.cos(angle)
        v = self.r * np.sin(angle)
        
        uhat = self.uhat
        vhat = cross(self.direction, uhat)
        circlePoints = np.outer(u, uhat) + np.outer(v, vhat)

        vertices = []
        for i in range(numCircles):
            t = i / float(numCircles - 1)
            p = self.start + t * self.direction * self.length
            vertices.append(circlePoints + p)

        vertices = np.vstack(vertices)

        indices = []
        for i in range(numCircles - 1):
            for j in range(numPointsPerCircle):
                next_j = (j + 1) % numPointsPerCircle
                indices.extend([
                    [i * numPointsPerCircle + j, i * numPointsPerCircle + next_j, (i + 1) * numPointsPerCircle + j],
                    [(i + 1) * numPointsPerCircle + j, i * numPointsPerCircle + next_j, (i + 1) * numPointsPerCircle + next_j]
                ])

        return vertices, np.array(indices)
    
    def addToPlot(self, ax, numPointsPerCircle=32, color='black', alpha=0.5, frame=False, numCircles=2, edgeColor=None):
        circles = self.interpolateCircles(numPointsPerCircle, numCircles)
        X = circles[:,:,0]
        Y = circles[:,:,1]
        Z = circles[:,:,2]
        if frame:
            return ax.plot_wireframe(X, Y, Z, color=edgeColor, alpha=alpha)
        elif edgeColor is None:
            return ax.plot_surface(X, Y, Z, color=color, alpha=alpha)
        else:
            return ax.plot_surface(X, Y, Z, color=color, alpha=alpha, edgecolor=edgeColor)
    
    def addToWidget(self, widget, numPointsPerCircle=32, numCircles=10, color_list=cylinderColorList, is_joint=False):
        import pyqtgraph.opengl as gl
        vertices, indices = self.interpolateQtCircles(numPointsPerCircle, numCircles)
        meshdata = gl.MeshData(vertexes=vertices, faces=indices)
        meshitem = gl.GLMeshItem(meshdata=meshdata, color=tuple(color_list), shader='shaded', smooth=True)
        meshitem.setGLOptions('translucent')
        if (is_joint):
            meshitem.setObjectName("Joint")
        else:
            meshitem.setObjectName("Link")
        
        widget.plot_widget.addItem(meshitem)

        return vertices, indices
    
    def show(self, numPointsPerCircle=32, color='black', alpha=0.5, frame=False, numCircles=2, block=blockDefault, edgeColor='black'):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandles = self.addToPlot(ax, numPointsPerCircle, color, alpha, frame, numCircles)
        ax.set_aspect('equal')
        plt.show(block=block)

def RotationAboutLine(rotAxisDir : np.ndarray,
                      rotAxisPoint : np.ndarray,
                      angle : float) -> SE3:
    R = SO3.AngleAxis(angle, rotAxisDir)
    t = (np.eye(3) - R) @ rotAxisPoint.reshape(3,1)
    return SE3.Rt(R,t)

class Elbow:
    def __init__(self, radius : float, StartFrame : SE3, bendingAngle : float, 
                 rotationalAxisAngle : float, EPSILON : float = 0.0001):
        rotationalAxisAngle = np.mod(rotationalAxisAngle, 2*np.pi)
        bendingAngle = math.remainder(bendingAngle, 2*np.pi) #wrap to [-pi,pi]
        assert(abs(bendingAngle) < np.pi)
        assert(abs(bendingAngle) > EPSILON)
        if bendingAngle < 0:
            bendingAngle = abs(bendingAngle)
            rotationalAxisAngle = np.mod(rotationalAxisAngle+np.pi, 2*np.pi)
        
        self.StartFrame = StartFrame
        self.r = radius
        self.rotAxisDirLocal = SO3.Rx(rotationalAxisAngle) * np.array([0,1,0])
        self.bendingAngle = bendingAngle
        self.dw = self.r * np.tan(self.bendingAngle / 2)
        
        if bendingAngle < EPSILON:
            self.Forward = SE3()
            self.Rotate = SE3()
            self.Transformation = SE3()
        else:
            self.Forward = SE3.Tx(self.dw)
            self.Rotate = SE3.AngleAxis(self.bendingAngle, self.rotAxisDirLocal) 
            self.Transformation = self.Forward @ self.Rotate @ self.Forward
        
        self.EndFrame = self.StartFrame @ self.Transformation
        
        self.HalfRotate = SE3.AngleAxis(self.bendingAngle / 2, self.rotAxisDirLocal)
        self.midPlaneNormal = (self.StartFrame * self.HalfRotate).R[:,0]
        self.midPoint = (self.StartFrame * self.Forward).t
        
        # self.midPoint is self.dw past the center of StartCircle,
        # so it's sqrt(self.dw**2 + self.r**2) from each point on StartCircle.
        # Meanwhile the elbow tip is 2*self.dw forward from its corresponding 
        # point on StartCircle, so it's self.dw forward and self.r outwards 
        # from self.midPoint, i.e., it's also sqrt(self.dw**2 + self.r**2) away
        self.boundingRadius = np.sqrt(self.dw**2 + self.r**2)
    
    def boundingBall(self) -> Ball:
        return Ball(self.midPoint, self.boundingRadius)
    
    def circleEllipseCircle(self, numSides : int = 32) -> tuple:
        count = numSides+1 # the +1 is because the first and last point will be identical
        angle = np.linspace(0, 2*np.pi, count) 
        ahat = self.StartFrame.R[:,0]
        bhat = self.StartFrame.R[:,1]
        chat = self.StartFrame.R[:,2]
        u = self.r * np.cos(angle).reshape(-1,1)
        v = self.r * np.sin(angle).reshape(-1,1)
        StartCircle = self.StartFrame.t + u @ bhat.reshape(1,3) + v @ chat.reshape(1,3)
        
        midPlane = Plane(self.midPoint, self.midPlaneNormal)
        MidEllipse = midPlane.intersectionsWithParallelLines(StartCircle, ahat)
        
        endBhat = self.EndFrame.R[:,1]
        endChat = self.EndFrame.R[:,2]
        endPoint = self.EndFrame.t
        EndCircle = endPoint + u @ endBhat.reshape(1,3) + v @ endChat.reshape(1,3)
        
        return StartCircle, MidEllipse, EndCircle
    
    def circleEllipseCircleQT(self, numSides : int = 32):
        StartCircle, MidEllipse, EndCircle = self.circleEllipseCircle(numSides)
        vertices = np.vstack((StartCircle[:numSides], MidEllipse[:numSides], EndCircle[:numSides]))

        faces = []
        for i in range(2):
            for j in range(numSides): 
                next_index = (j + 1) % numSides
                faces.append([i * numSides + j, i * numSides + next_index, (i + 1) * numSides + j])
                faces.append([(i + 1) * numSides + j, i * numSides + next_index, (i + 1) * numSides + next_index])

        faces = np.array(faces)

        return vertices, faces
    
    def addToPlot(self, ax, numSides : int = 32, color : str = 'black', 
                  alpha : float = 0.5, wireFrame : bool = False, 
                  showFrames : bool = False):
        
        StartCircle, MidEllipse, EndCircle = self.circleEllipseCircle(numSides)
        ellipses = np.array([StartCircle, MidEllipse, EndCircle])
        X = ellipses[:,:,0]
        Y = ellipses[:,:,1]
        Z = ellipses[:,:,2]
        
        if wireFrame:
            surfaceHandle = ax.plot_wireframe(X, Y, Z, color=color, alpha=alpha)
        else:
            surfaceHandle = ax.plot_surface(X, Y, Z, color=color, alpha=alpha)
        
        frameHandles = []
        if showFrames:
            Fwd = self.StartFrame @ self.Forward 
            FwdRot = Fwd @ self.Rotate
            FwdRotFwd = FwdRot @ self.Forward
            Poses = np.array([self.StartFrame, Fwd, FwdRot, FwdRotFwd])
            aHats, bHats, cHats, origins = addPosesToPlot(Poses, ax, 
                                        axisLength=self.r, xColor='darkred', 
                                        yColor='darkblue', zColor='darkgreen')
            frameHandles = [aHats, bHats, cHats, origins]
            x,y,z = Fwd.t
            u,v,w = np.cross(Fwd.R[:,0], FwdRot.R[:,0])
            ax.quiver(x,y,z,u,v,w,length=self.r,normalize=True)
        
        return frameHandles
    
    def addToWidget(self, widget, numSides : int = 16, color_list=elbowColorList, 
                  alpha : float = 1.0, wireFrame : bool = True, 
                  showFrames : bool = False, debug : bool = False):
        import pyqtgraph.opengl as gl
        
        vertices, faces = self.circleEllipseCircleQT(numSides)

        meshdata = gl.MeshData(vertexes=vertices, faces=faces)
        meshitem = gl.GLMeshItem(meshdata=meshdata, color=tuple(color_list), drawEdges=wireFrame, shader='shaded', smooth=True)
        meshitem.setObjectName("Link")
        
        if showFrames:
            Fwd = self.StartFrame @ self.Forward 
            FwdRot = Fwd @ self.Rotate
            FwdRotFwd = FwdRot @ self.Forward
            Poses = np.array([self.StartFrame, Fwd, FwdRot, FwdRotFwd])

            ax = plt.figure().add_subplot(projection='3d')
            addPosesToPlotQT(Poses, ax, widget.plot_widget,
                                        axisLength=1, xColor=xPoseColor, 
                                        yColor=yPoseColor, zColor=zPoseColor)

            x,y,z = Fwd.t
            u,v,w = np.cross(Fwd.R[:,0], FwdRot.R[:,0])
            line = gl.GLLinePlotItem(pos=np.array([[x,y,z], [0.5*u,0.5*v,0.5*w]]), color=lineColor, width=5) 
            widget.plot_widget.addItem(line)

        meshitem.setGLOptions('translucent')
        widget.plot_widget.addItem(meshitem)
    
    def show(self, numSides : int = 32, color : str = 'black', 
             alpha : float = 0.5, wireFrame : bool = False, 
             showFrames : bool = False, block : bool = False):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandles = self.addToPlot(ax, numSides, color, alpha, wireFrame, 
                                     showFrames)
        ax.set_aspect('equal')
        plt.show(block=block)

class CompoundElbow:
    def __init__(self, radius : float, StartFrame : SE3, bendingAngle : float, 
                 rotationalAxisAngle : float, maxAnglePerElbow : float = np.pi/2, 
                 EPSILON : float = 0.0001):
        rotationalAxisAngle = np.mod(rotationalAxisAngle, 2*np.pi)
        bendingAngle = math.remainder(bendingAngle, 2*np.pi) #wrap to [-pi,pi]
        if abs(bendingAngle) > np.pi+EPSILON:
            raise ValueError("Trying to construct CompoundElbow with abs(bendingAngle) > pi")
        if abs(bendingAngle) < EPSILON:
            raise ValueError("Trying to construct CompoundElbow with abs(bendingAngle) < EPSILON")
        if bendingAngle < 0:
            bendingAngle = abs(bendingAngle)
            rotationalAxisAngle = np.mod(rotationalAxisAngle+np.pi, 2*np.pi)
        
        assert(maxAnglePerElbow >= 0 and maxAnglePerElbow <= np.pi)
        numElbows = (int)(np.ceil(bendingAngle / maxAnglePerElbow)) 
        anglePerElbow = bendingAngle / numElbows
        self.elbows = []
        CurrentFrame = StartFrame
        for i in range(numElbows):
            nextElbow = Elbow(radius, CurrentFrame, anglePerElbow, 
                                     rotationalAxisAngle, EPSILON)
            self.elbows.append(nextElbow)
            CurrentFrame = nextElbow.EndFrame
        
        self.StartFrame = StartFrame
        self.EndFrame = self.elbows[-1].EndFrame
    
    def addToPlot(self, ax, numSides : int = 32, color : str = 'black', 
                  alpha : float = 0.5, wireFrame : bool = False, 
                  showFrames : bool = True, showBoundingBall : bool = False):
        allHandleSets = []
        for elbow in self.elbows:
            handleSet = elbow.addToPlot(ax, numSides, color, alpha, wireFrame, showFrames)
            if showFrames:
                allHandleSets.append(handleSet)
        
        if showBoundingBall:
            self.boundingBall().addToPlot(ax, color=color, alpha = 0.25*alpha)
        
        return allHandleSets
    
    def circleEllipseCircleQT(self, numSides : int = 32):
        allHandleSets = []

        vertices = []
        faces = []
        vertex_offset = 0

        for elbow in self.elbows:
            v, f = elbow.circleEllipseCircleQT(numSides)

            vertices.extend(v.tolist())
            f2 = [[item + vertex_offset for item in sublist] for sublist in f.tolist()]
            faces.extend(f2)
            vertex_offset += len(v)

        return vertices, faces
    
    def addToWidget(self, widget, numSides : int = 32, color_list=compoundElbowColorList, 
                  alpha : float = 1.0, wireFrame : bool = False, 
                  showFrames : bool = True, showBoundingBall : bool = False, debug : bool = False):
        allHandleSets = []
        
        for elbow in self.elbows:
            elbow.addToWidget(widget, numSides, color_list, alpha, wireFrame, showFrames, debug)
        
        if showBoundingBall:
            self.boundingBall().addToWidget(widget, color=tuple(color_list), alpha = 0.25*alpha)
        
        return allHandleSets
    
    def boundingBall(self):
        ball = self.elbows[0].boundingBall()
        for elbow in self.elbows[1:]:
            ball = minBoundingBall(ball, elbow.boundingBall())
        return ball
    
    def show(self, numSides : int = 32, color : str = 'black', 
             alpha : float = 0.5, wireFrame : bool = False, 
             showFrames : bool = True, showBoundingBall : bool = False,
             block : bool = True):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandles = self.addToPlot(ax, numSides, color, alpha, wireFrame, 
                                     showFrames, showBoundingBall)
        ax.set_aspect('equal')
        plt.show(block=block)
        
    
class Arc3D:
    def __init__(self, circleCenter, startPoint, startDir, theta):
        assert(norm(startPoint-circleCenter)>0)
        # verify orthogonality up to numerical stability
        assert(abs(np.dot(startDir, startPoint-circleCenter)) < 0.000001)
        
        self.circleCenter = circleCenter
        self.startPoint = startPoint
        self.startTangent = startDir / norm(startDir)
        self.theta = theta
        
        self.centerToStart = self.startPoint - self.circleCenter
        self.r = norm(self.centerToStart)
        self.startNormal = - self.centerToStart / self.r
                
        self.binormal = cross(self.startTangent, self.startNormal)
        self.rot = Rotation.from_rotvec(self.theta * self.binormal) #TODO: replace this with SO3 to be consistent?
        self.centerToEnd = self.rot.apply(self.centerToStart)
        self.endPoint = self.circleCenter + self.centerToEnd
        self.endNormal = - self.centerToEnd / self.r
        self.endTangent = cross(self.endNormal, self.binormal)
    
    def interpolate(self, count=50, xp : ModuleType = np) -> ArrayLike:
        angle = xp.linspace(0, self.theta, count).reshape(-1,1)
        u = self.r * xp.cos(angle)
        v = self.r * xp.sin(angle)
        
        # construct basis for circle plane
        uhat = xp.asarray(-self.startNormal).reshape(1,3)
        vhat = xp.cross(self.binormal, uhat).reshape(1,3)
        
        # 3d circle points
        return xp.asarray(self.circleCenter) + u @ uhat + v @ vhat
    
    def interpolate_vectorized(self, t_array: np.ndarray) -> np.ndarray:
        """
        Vectorized interpolation: compute 3D positions for multiple t values at once.
        
        Parameters:
        -----------
        t_array : np.ndarray
            Array of parameter values in [0, 1] (shape (n,))
            
        Returns:
        --------
        np.ndarray
            Array of 3D positions (shape (n, 3))
        """
        # Clamp t values to [0, 1]
        t_array = np.clip(t_array, 0.0, 1.0)
        
        # Convert t to angle values: angle = t * theta
        angle = (t_array * self.theta).reshape(-1, 1)
        u = self.r * np.cos(angle)
        v = self.r * np.sin(angle)
        
        # Construct basis for circle plane
        uhat = -self.startNormal.reshape(1, 3)
        vhat = cross(self.binormal, uhat).reshape(1, 3)
        
        # 3D circle points
        return self.circleCenter + u @ uhat + v @ vhat
    
    def interpolateAt(self, t: float) -> np.ndarray:
        """
        Return 3D position at parameter t in [0, 1] along the arc
        """
        t = max(0.0, min(t, 1.0))
        return self.interpolate_vectorized(np.array([t]))[0]
    
    def localOrientation(self) -> np.ndarray:
        """
        Return the precomputed transformation matrix from world coordinates
        to local arc coordinates.
        """
        if not hasattr(self, '_worldToLocalRotation'):
            self._computeLocalFrame()
        return self._worldToLocalRotation

    def sdfSinCos(self) -> np.ndarray:
        """
        Return the precomputed (sin(halfTheta), cos(halfTheta)) for the arc SDF.
        """
        if not hasattr(self, '_sdfSinCos'):
            self._computeLocalFrame()
        return self._sdfSinCos
    
    
    def _computeLocalFrame(self):
        """
        Precompute the transformation matrix from world coordinates to local arc coordinates.
        
        The capped torus SDF formula (from Inigo Quilez) assumes:
        - The torus lies in the XY plane, centered at the origin
        - The arc is SYMMETRIC about the X-axis, spanning angles [-theta, +theta]
        - The formula uses abs(x) to exploit this symmetry
        - The arc "caps" (endpoints) are at angles ┬▒theta from the +X axis
        
        For our Arc3D, we have:
        - startPoint at angle 0 (beginning of arc)
        - endPoint at angle theta (end of arc)
        
        To use the symmetric SDF, we align the LOCAL X-axis with the MIDPOINT
        of the arc (at angle theta/2). This way:
        - The arc spans from -theta/2 to +theta/2 in local coordinates
        - Both endpoints are equidistant from the X-axis
        - The abs(x) symmetry is correctly utilized
        
        The local coordinate system is:
        - Origin: at the arc's circle center
        - Y-axis: points radially outward at the arc's midpoint (IQ's formula expects arc centered on +Y)
        - X-axis: tangent direction at midpoint (perpendicular to Y in the arc plane)
        - Z-axis: binormal (perpendicular to arc plane)
        """
        # Compute the midpoint direction: rotate -startNormal by theta/2 around binormal
        # startNormal points INWARD (toward center), so -startNormal points outward at start
        halfTheta = self.theta / 2
        halfAngleRot = Rotation.from_rotvec(halfTheta * self.binormal)
        centerToMid = halfAngleRot.apply(-self.startNormal)  # radial outward at midpoint
        
        # Y-axis: radial outward at arc midpoint (IQ's formula has arc centered on +Y)
        localY = centerToMid / norm(centerToMid)
        
        # Z-axis: binormal (perpendicular to arc plane)
        localZ = self.binormal
        
        # X-axis: completes right-handed frame, tangent at midpoint
        # cross(Y, Z) gives the tangent direction at the midpoint
        localX = cross(localY, localZ)
        
        # Build rotation matrix: columns are local basis vectors expressed in world coords
        # To transform world -> local, we use the transpose (inverse for orthonormal basis)
        self._worldToLocalRotation = np.column_stack([localX, localY, localZ]).T
        
        # Precompute sin and cos of HALF the arc angle for the symmetric SDF
        # The SDF expects the arc to span [-halfTheta, +halfTheta]
        self._sdfSinCos = np.array([np.sin(halfTheta), np.cos(halfTheta)])
    
    def _sdFlatEndedTorus(self, p: np.ndarray, sc: np.ndarray, ra: float, rb: float) -> float:
        """
        Signed Distance Function for a flat-ended torus in local coordinates.
        
        Modified from Inigo Quilez's capped torus SDF to have flat disc ends
        instead of spherical caps.
        
        The torus is centered at the origin in the XY plane. The arc is SYMMETRIC
        about the Y-axis, spanning angles from (90┬░-halfTheta) to (90┬░+halfTheta).
        With abs(p.x), it handles both sides.
        
        Parameters:
        -----------
        p : np.ndarray
            3D point in local coordinate system (shape (3,))
            - x: tangent direction at arc midpoint
            - y: radial direction at arc midpoint (outward from torus center)
            - z: axial direction (perpendicular to torus plane)
        sc : np.ndarray
            (sin(halfTheta), cos(halfTheta)) where halfTheta = arcAngle/2
        ra : float
            Major radius (distance from torus center to tube center)
        rb : float
            Minor radius (tube radius)
        
        Returns:
        --------
        float
            Signed distance (negative inside, positive outside)
        """
        # Use abs(x) for symmetry, work in XY plane
        p_xy = np.array([abs(p[0]), p[1]])
        pz = p[2]
        
        # Endpoint center on the torus ring (at angle halfTheta from +Y axis)
        # sc = (sin, cos), so endpoint is at ra * sc
        endCenter = ra * sc
        
        # Tangent at endpoint: rotate sc by -90┬░ ΓåÆ (cos, -sin)
        tangent = np.array([sc[1], -sc[0]])
        
        # Vector from endpoint to query point
        toPoint = p_xy - endCenter
        
        # How far past the arc endpoint are we?
        pastEnd = np.dot(toPoint, tangent)
        
        if pastEnd <= 0.0:
            # Inside arc span - standard torus formula
            p_len = np.linalg.norm(p_xy)
            if sc[1] * p_xy[0] > sc[0] * p_xy[1]:
                k = np.dot(sc, p_xy)
            else:
                k = p_len
            return np.sqrt(p_len*p_len + pz*pz + ra*ra - 2.0*ra*k) - rb
        else:
            # Past arc endpoint - distance to flat disc
            # Radial distance from tube axis (project onto sc which points radially)
            radialInPlane = np.dot(toPoint, sc)
            discDist = np.sqrt(radialInPlane*radialInPlane + pz*pz)
            
            # 2D SDF to disc edge
            outsideDisc = max(discDist - rb, 0.0)
            return np.sqrt(pastEnd*pastEnd + outsideDisc*outsideDisc)
    
    def sdf(self, point: np.ndarray, radius: float) -> float:
        """
        Compute the signed distance from a 3D point to this arc's tubular volume.
        
        The arc is treated as a torus section (tube bent along the arc) with 
        flat disc ends instead of spherical caps.
        
        Parameters:
        -----------
        point : np.ndarray
            3D point in world coordinates
        radius : float
            Tube radius around the arc centerline
        
        Returns:
        --------
        float
            Signed distance (negative inside the tube, positive outside)
        """
        # Lazy initialization of the precomputed transformation matrix
        if not hasattr(self, '_worldToLocalRotation'):
            self._computeLocalFrame()
        
        # Transform point from world coordinates to local arc coordinates:
        # 1. Translate so circle center is at origin
        # 2. Rotate so arc midpoint is on +Y axis (for IQ's formula)
        localP = self._worldToLocalRotation @ (point - self.circleCenter)
        
        # Apply the flat-ended torus SDF in local coordinates
        return self._sdFlatEndedTorus(localP, self._sdfSinCos, self.r, radius)
    
    def addToPlot(self, ax, color='black', alpha=1, showDirections=False):
        X,Y,Z = self.interpolate().T
        if showDirections:
            ax.quiver(*self.startPoint, *self.startTangent, length=self.r, color='green')
            ax.quiver(*self.endPoint, *self.endTangent, length=self.r, color='blue')
        return ax.plot(X, Y, Z, color=color, alpha=alpha)

    def show(self, color='black', alpha=1, block=blockDefault, showDirections=False):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandle = self.addToPlot(ax, color, alpha, showDirections)
        ax.set_aspect('equal')
        plt.show(block=block)


def trussManifold(vertices, edges, diameter : float) -> m3d.Manifold:
    T = m3d.Manifold()
    nodes = [m3d.Manifold.sphere(radius=diameter/2, circular_segments=3).translate(vertex[:3]) for vertex in vertices]
    for edge in edges:
        T += (nodes[edge[0]] + nodes[edge[1]]).hull()
    return T


def connectOuterToInner(outerVertices : np.ndarray, outerEdges : np.ndarray, 
                 innerVertices : np.ndarray, innerEdges : np.ndarray, 
                 nearestCount=1) -> tuple[np.ndarray, np.ndarray]:
    """
    Given outer and inner truss specifications (vertices and edges),
    concatenate them into a single truss specification,
    with every outer vertex connected to its nearestCount closest inner vertices.
    
    :param outerVertices: m x 3 numpy array of outer truss vertex positions
    :param outerEdges: p x 2 numpy array of outer truss edges by vertex index
    :param innerVertices: n x 3 numpy array of inner truss vertex positions
    :param innerEdges: q x 2 numpy array of inner truss edges by vertex index
    :param nearestCount: number of nearest inner vertices to connect to each outer vertex
    :return: tuple (vertices, edges) representing the combined truss
    :rtype: (numpy.ndarray, numpy.ndarray)
    """
    if nearestCount <= 0:
        raise ValueError("nearestCount must be positive")
    elif nearestCount > innerVertices.shape[0]:
        raise ValueError("nearestCount cannot exceed number of inner vertices")

    combinedVertices = np.vstack((outerVertices, innerVertices))
    newEdgesList = []
    for outerVertexIndex, outerVertex in enumerate(outerVertices):
        distances = norm(innerVertices - outerVertex.reshape(1,3), axis=1)
        nearestInnerIndices = np.argsort(distances)[:nearestCount]
        for innerIndex in nearestInnerIndices:
            newEdgesList.append([outerVertexIndex, innerIndex + outerVertices.shape[0]])
            #newEdgesList.append([np.where((outerVertices == outerVertex).all(axis=1))[0][0], 
            #                     innerIndex + outerVertices.shape[0]])

    combinedEdges = np.vstack((outerEdges, 
                               innerEdges + outerVertices.shape[0],
                               np.array(newEdgesList)))

    return combinedVertices, combinedEdges

def facesToEdges(faces: np.ndarray) -> np.ndarray:
    """
    Given an array of faces as triangles, quadrilaterals, etc (by vertex index),
    Return an array of edges without duplication
    """
    n = faces.shape[0]  # number of faces
    d = faces.shape[1]  # vertices per face
    
    # Create all edges by pairing each vertex with the next (wrapping around)
    edges = np.stack([faces, np.roll(faces, -1, axis=1)], axis=-1)  # shape: (n, d, 2)
    edges = edges.reshape(-1, 2)  # flatten to (n*d, 2)
    
    # Sort each edge so (a,b) and (b,a) are treated the same
    edges = np.sort(edges, axis=1)
    
    # Remove duplicates
    edges = np.unique(edges, axis=0)
    
    return edges

def tetrahedronsToEdges(tetrahedrons : np.ndarray) -> np.ndarray:
    """
    Compute the edges from a list of tetrahedrons specified by vertex index
    
    :param tetrahedrons: n x 4 numpy array of integers
    :return: array of unique edges as pairs of vertex indices
    :rtype: _ x 2 numpy array of integers
    """
    n = tetrahedrons.shape[0]  # number of tetrahedrons
    
    # Each tetrahedron has 6 edges: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
    edge_pairs = np.array([[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]])
    
    # Create all edges by indexing into tetrahedrons
    # Shape: (n, 6, 2) where n is number of tetrahedrons
    edges = tetrahedrons[:, edge_pairs]
    
    # Flatten to (n*6, 2)
    edges = edges.reshape(-1, 2)
    
    # Sort each edge so (a,b) and (b,a) are treated the same
    edges = np.sort(edges, axis=1)
    
    # Remove duplicates
    edges = np.unique(edges, axis=0)
    
    return edges

def manifoldToGraph(manifold : m3d.Manifold) -> tuple[np.ndarray, np.ndarray]:
    mesh = manifold.to_mesh()
    vertices = mesh.vert_properties
    faces = mesh.tri_verts
    return vertices, facesToEdges(faces)

def manifoldToTruss(manifold : m3d.Manifold, diameter : float, 
                    infill : bool = False) -> m3d.Manifold:
    mesh = manifold.to_mesh()
    vertices = mesh.vert_properties
    faces = mesh.tri_verts
    if infill:
        tetVerts, tets = tetrahedralize(vertices, faces,
                                        edge_length_fac=1,
                                        optimize=True)
        return trussManifold(tetVerts, tetrahedronsToEdges(tets), diameter)
    else:
        edges = facesToEdges(mesh.tri_verts)
        return trussManifold(vertices, edges, diameter)

class Bend:
    def __init__(self, arcRadius : float, StartFrame : SE3, bendingAngle : float, 
                 rotationalAxisAngle : float, startRadius : float = None,
                 endRadius : float = None, numSides : int = 20, 
                 maxSectionAngle : float = np.pi/8, EPSILON : float = 0.0001):     
        self.rotationalAxisAngle = np.mod(rotationalAxisAngle, 2*np.pi)
        bendingAngle = math.remainder(bendingAngle, 2*np.pi) #wrap to [-pi,pi]
        assert(abs(bendingAngle) <= np.pi)
        assert(abs(bendingAngle) > EPSILON)
        if bendingAngle < 0:
            bendingAngle = abs(bendingAngle)
            self.rotationalAxisAngle = np.mod(self.rotationalAxisAngle+np.pi, 2*np.pi)

        self.StartFrame = StartFrame
        self.arcRadius = arcRadius
        self.startRadius = self.arcRadius if startRadius is None else startRadius
        self.endRadius = self.arcRadius if endRadius is None else endRadius
        self.rotAxisDirLocal = SO3.Rx(self.rotationalAxisAngle) * np.array([0,1,0])
        self.bendingAngle = bendingAngle
        self.numSides = numSides
        self.EPSILON = EPSILON
        self.DISTANCE_EPSILON = self.arcRadius * self.EPSILON
        self.maxSectionAngle = maxSectionAngle

        self.numSections = 2 * math.ceil(abs(bendingAngle) / self.maxSectionAngle)
        self.numCircles = self.numSections + 1
        self.anglePerSection = bendingAngle / self.numSections
        self.dwPerSection = self.arcRadius * np.tan(self.anglePerSection / 2)
        
        Forward = SE3.Tx(self.dwPerSection)
        Rotate = SE3.AngleAxis(self.anglePerSection, self.rotAxisDirLocal)
        self.TransformPerSection = Forward @ Rotate @ Forward

        self.poses = [self.StartFrame]
        for i in range(1, self.numCircles):
            self.poses.append(self.poses[-1] @ self.TransformPerSection)
        circles = []
        self.radii = np.linspace(self.startRadius, self.endRadius, self.numCircles)
        for i in range(self.numCircles):
            pose = self.poses[i]
            circles.append(Circle3D(radius=self.radii[i], center=pose.t, normal=pose.R[:,0], radialVector=pose.R[:,1]))
        self.circles = np.array([c.interpolate(self.numSides+1) for c in circles])

    def plotCircles(self, ax):
        for i in range(self.numCircles):
            ax.plot(self.circles[i,:,0], self.circles[i,:,1], self.circles[i,:,2], marker='o')
        addPosesToPlot(np.array(self.poses), ax, axisLength=0.2)
    
    def trimesh(self):
        # Create a mesh by connecting the circles
        vertices = self.circles.reshape((-1, 3))
        faces = []
        for i in range(self.numCircles-1):
            for j in range(self.numSides):
                p0 = i * (self.numSides + 1) + j
                p1 = p0 + 1
                p2 = p0 + (self.numSides + 1)
                p3 = p2 + 1
                faces.append([p0, p2, p1])
                faces.append([p1, p2, p3])
        faces = np.array(faces)
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
        return mesh
    
    def solid_trimesh(self):
        # Create a solid mesh by capping the ends
        mesh = self.trimesh()
        startCircle = self.circles[0]
        endCircle = self.circles[-1]
        startCap = trimesh.Trimesh(vertices=startCircle, 
                                   faces=[[0,i+1,(i+1)%self.numSides+1] for i in range(self.numSides)])
        endCap = trimesh.Trimesh(vertices=endCircle, 
                                 faces=[[0,i+1,(i+1)%self.numSides+1] for i in range(self.numSides)])
        #endCap.apply_translation(endCircle.c - endCap.vertices[0])
        mesh = trimesh.util.concatenate([mesh, startCap, endCap])
        return mesh
    
    def manifold(self, hull : bool = False, extendForward : float = 0, extendBackward : float = 0, thickness : float = None, truss : bool = False,
                 trussNumSides : int = 6, trussMaxSectionAngle : float = np.pi/4) -> m3d.Manifold:       
        if truss:
            if thickness is None:
                raise ValueError("Must specify wall thickness for truss")
            inner = Bend(self.arcRadius, self.StartFrame, self.bendingAngle, self.rotationalAxisAngle,
                                    self.startRadius - thickness/2, self.endRadius - thickness/2,
                                    trussNumSides, trussMaxSectionAngle, self.EPSILON)
            innerSolid = inner.manifold(hull)
            shape = manifoldToTruss(innerSolid, thickness)

        elif hull:
            shape = m3d.Manifold.hull_points(self.circles.reshape((-1,3)))
        else:
            shape = m3d.Manifold.hull_points(self.circles[[0,1]].reshape((-1,3)))
            for i in range(1, self.numCircles-1):
                shape += m3d.Manifold.hull_points(self.circles[[i-1, i, i+1]].reshape((-1,3)))
        
        if extendBackward != 0:
            startCircle = self.circles[0]
            circleBackward = startCircle - self.poses[0].R[:,0]*extendBackward
            circleStack = np.vstack((circleBackward, startCircle, self.circles[1])) if (not truss and self.circles.shape[0] > 1) else np.vstack((circleBackward, startCircle))
            startCap = m3d.Manifold.hull_points(circleStack)
            shape += startCap
        if extendForward != 0:
            endCircle = self.circles[-1]
            circleForward = endCircle + self.poses[-1].R[:,0]*extendForward
            circleStack = np.vstack((self.circles[-2], endCircle, circleForward)) if (not truss and self.circles.shape[0] > 1) else np.vstack((endCircle, circleForward))
            endCap = m3d.Manifold.hull_points(circleStack)
            shape += endCap

        if thickness is not None:
            if not (thickness > 0 and thickness < min(self.startRadius, self.endRadius)):
                raise ValueError("Invalid wall thickness for Bend manifold")
            inner = Bend(self.arcRadius, self.StartFrame, self.bendingAngle, self.rotationalAxisAngle,
                                      self.startRadius - thickness, self.endRadius - thickness,
                                      self.numSides, self.maxSectionAngle, self.EPSILON)
            shape -= inner.manifold(hull, extendForward=extendForward+self.DISTANCE_EPSILON,
                                    extendBackward=extendBackward+self.DISTANCE_EPSILON)
        
        return shape
        
    def rediscretize(self, numSides : int = 20, maxSectionAngle : float = np.pi/8) -> Bend:
        return Bend(self.arcRadius, self.StartFrame, self.bendingAngle, self.rotationalAxisAngle,
                    self.startRadius, self.endRadius, numSides, maxSectionAngle, self.EPSILON)
        
    def truss(self, thickness : float, hull : bool = False, extendForward : float = 0, extendBackward : float = 0) -> m3d.Manifold:
        inner = Bend(self.arcRadius, self.StartFrame, self.bendingAngle, self.rotationalAxisAngle,
                                    self.startRadius - thickness/2, self.endRadius - thickness/2,
                                    self.numSides, self.maxSectionAngle, self.EPSILON)
        innerSolid = inner.manifold(hull)
        truss = manifoldToTruss(innerSolid, thickness)
        



class HollowBend(Bend):
    def __init__(self, arcRadius : float, StartFrame : SE3, bendingAngle : float, 
                 rotationalAxisAngle : float, wallThickness : float,
                 startRadius : float = None, endRadius : float = None, 
                 numSides : int = 20, maxSectionAngle : float = np.pi/8, EPSILON : float = 0.0001):
        super().__init__(arcRadius, StartFrame, bendingAngle, rotationalAxisAngle,
                         startRadius, endRadius, numSides, maxSectionAngle, EPSILON)
        assert(wallThickness > 0 and wallThickness < min(self.startRadius, self.endRadius))
        self.wallThickness = wallThickness
        
        self.inner = Bend(arcRadius, StartFrame, bendingAngle, rotationalAxisAngle,
                                      self.startRadius - wallThickness,
                                      self.endRadius - wallThickness,
                                      numSides, maxSectionAngle, EPSILON)

    def manifold(self, hull : bool = False, extendForward : float = 0, extendBackward : float = 0) -> m3d.Manifold:
        outerHull = super().manifold(hull, extendForward=extendForward, extendBackward=extendBackward)
        innerHull = self.inner.manifold(hull, extendForward=extendForward+self.DISTANCE_EPSILON,
                                        extendBackward=extendBackward+self.DISTANCE_EPSILON)
        return outerHull - innerHull



# arc from a given starting point+direction to a given ending direction 
# (which cannot be parallel to the starting direction)
def arcToDirection(startPoint, startDir, endDir, r) -> Arc3D:
    startDir /= norm(startDir)
    endDir /= norm(endDir)
    if norm(startDir - endDir) < 1e-8:
        # find any direction orthogonal to startDir to use as inward
        # in the nullspace of something
        inward = unitNormalToBoth(startDir, endDir)
        center = startPoint + r*inward
        return Arc3D(center, startPoint, startDir, 0)

    normal = np.cross(startDir, endDir)
    if norm(normal) == 0:
        raise ValueError("startDir and endDir cannot be parallel")
    normal = normal / norm(normal)
    inward = np.cross(normal, startDir)
    inward = inward / norm(inward)
    center = startPoint + r*inward
    angle = np.arccos(np.dot(startDir, endDir))
    return Arc3D(center, startPoint, startDir, angle)
    

# add given reference frames to matplotlib figure ax with a 3d subplot
# pose is a matrix of SE3() objects
# returns the plot handles for the xHats, yHats, zHats, origins
def addPosesToPlot(Poses, ax, axisLength, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, oColors='black', makeAxisLimitsIncludeTips=True):
    if Poses.shape == (4,4): # so it can plot a single frame
        Poses = np.array([Poses])
    
    ux, vx, wx = Poses[:,0:3,0].T # frame xhat coordinates
    uy, vy, wy = Poses[:,0:3,1].T # frame yhat coordinates
    uz, vz, wz = Poses[:,0:3,2].T # frame zhat coordinates
    ox, oy, oz = Poses[:,0:3,3].T # frame origin coordinates
    
    # https://matplotlib.org/stable/gallery/mplot3d/quiver3d.html
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.quiver.html#matplotlib.axes.Axes.quiver
    xHats = ax.quiver(ox, oy, oz, ux, vx, wx, length=axisLength, color=xColor) #plot xhat vectors
    yHats = ax.quiver(ox, oy, oz, uy, vy, wy, length=axisLength, color=yColor) #plot yhat vectors
    zHats = ax.quiver(ox, oy, oz, uz, vz, wz, length=axisLength, color=zColor) #plot zhat vectors
    origins = ax.scatter(ox, oy, oz, c=oColors)

    if makeAxisLimitsIncludeTips:
        oPlusXHats = Poses[:,0:3,3] + axisLength*Poses[:,0:3,0]
        ax.scatter(oPlusXHats[:,0], oPlusXHats[:,1], oPlusXHats[:,2], marker="")
        oPlusYHats = Poses[:,0:3,3] + axisLength*Poses[:,0:3,1]
        ax.scatter(oPlusYHats[:,0], oPlusYHats[:,1], oPlusYHats[:,2], marker="")
        oPlusZHats = Poses[:,0:3,3] + axisLength*Poses[:,0:3,2]
        ax.scatter(oPlusZHats[:,0], oPlusZHats[:,1], oPlusZHats[:,2], marker="")
    
    return (xHats, yHats, zHats, origins)

def showPoses(Poses, axisLength=1, xColor=xColorDefault, yColor=yColorDefault, zColor=zColorDefault, oColors='black', block=blockDefault, makeAxisLimitsIncludeTips=True):
    if type(Poses) == list:
        Poses = np.array(Poses)
    if Poses.shape == (4,4): # so it can plot a single frame
        Poses = np.array([Poses])
    ax = plt.figure().add_subplot(projection='3d')
    handles = addPosesToPlot(Poses, ax, axisLength, xColor, yColor, zColor, oColors, makeAxisLimitsIncludeTips)
    ax.set_xticks(np.arange(3))
    ax.set_yticks(np.arange(3))
    ax.set_zticks(np.arange(3))
    ax.set_aspect('equal')
    plt.show(block=block)



"""
Returns the common normal from line 1 to line 2, input in point-direction form.
If the lines intersect, return the given value representing undefined.
"""
def commonNormal(point1, direction1, point2, direction2, undefined=None):
    direction1 = direction1 / norm(direction1)
    direction2 = direction2 / norm(direction2)
    
    v12 = point2 - point1
    
    cp = cross(direction1, direction2)
    if norm(cp) > 0:
        nhat = cp / norm(cp)
    else: # z axes are parallel
        normalInPlane = cross(direction1, cross(v12, direction1))
        if norm(normalInPlane) == 0: #axes are coincident
            return undefined
        else:
            nhat = normalInPlane / norm(normalInPlane)
            
    projectionAlongNormal = dot(v12, nhat)
    if projectionAlongNormal > 0: # nhat is oriented correctly
        return nhat
    elif projectionAlongNormal < 0: # nhat is opposite orientation
        return -nhat
    else: # axes instersect
        return undefined


def shortestDistanceBetweenLines(point1, direction1, point2, direction2):
    direction1 = direction1 / norm(direction1)
    direction2 = direction2 / norm(direction2)
    cp = cross(direction1, direction2)
    norm_cp = norm(cp)
    if norm_cp == 0: #lines are parallel
        return norm(cross(direction1, point2 - point1)/norm(direction1))
    else:
        return abs(dot(point2 - point1, cp) / norm_cp)


class Torus:
    def __init__(self, majorRadius, minorRadius, center, axisDirection):
        self.R = majorRadius
        self.r = minorRadius
        self.c = center
        self.dhat = axisDirection / norm(axisDirection)
        self.uhat, self.vhat = null_space([self.dhat]).T

    def point(self, turnAngle, azumith):
        u = (self.R + self.r*np.cos(turnAngle))*np.cos(azumith)*self.uhat
        v = (self.R + self.r*np.cos(turnAngle))*np.sin(azumith)*self.vhat
        d = self.r*np.sin(turnAngle)*self.dhat
        return self.c + u + v + d
    
    def interpolate(self, numTurnTicks=32, numAzumithTicks=32):
        numTurnTicks += 1
        numAzumithTicks += 1
        turnAngles = np.linspace(0, 2*np.pi, numTurnTicks)
        azumiths = np.linspace(0, 2*np.pi, numAzumithTicks)
        grid = np.zeros((numTurnTicks, numAzumithTicks, 3))
        for i, turnAngle in enumerate(turnAngles):
            for j, azumith in enumerate(azumiths):
                grid[i,j] = self.point(turnAngle, azumith)
        return grid
    
    def addToPlot(self, ax, numTurnTicks=32, numAzumithTicks=32, colorMap=cm.Blues, alpha=0.5):
        grid = self.interpolate(numTurnTicks, numAzumithTicks)
        X = grid[:,:,0]
        Y = grid[:,:,1]
        Z = grid[:,:,2]
        return ax.plot_surface(X, Y, Z, cmap=colorMap, alpha=alpha)

    def show(self, numTurnTicks=32, numAzumithTicks=32, colorMap=cm.Blues, alpha=0.5, block=blockDefault):
        ax = plt.figure().add_subplot(projection='3d')
        plotHandle = self.addToPlot(ax, numTurnTicks, numAzumithTicks, colorMap, alpha)
        ax.set_aspect('equal')
        plt.show(block=block)

class HornTorus(Torus):
    def __init__(self, radius, center, axisDirection):
        super().__init__(radius, radius, center, axisDirection)


def sdf_aabb(
    min1: ArrayLike,
    max1: ArrayLike,
    min2: ArrayLike,
    max2: ArrayLike,
    xp: ModuleType = np,
) -> Union[float, ArrayLike]:
    """
    Signed distance between axis-aligned bounding boxes.
    
    Args:
        min1: Bottom-left corner(s), shape (3,) or (N, 3)
        max1: Top-right corner(s), shape (3,) or (N, 3)
        min2: Bottom-left corner(s), shape (3,) or (N, 3)
        max2: Top-right corner(s), shape (3,) or (N, 3)
        xp: Array module (np for numpy or cp for cupy)
    
    Returns:
        Signed distance(s). Negative if overlapping, positive if separated.
        Scalar for single box pair, (N,) array for N box pairs.
    """
    min1 = xp.atleast_2d(xp.asarray(min1))
    max1 = xp.atleast_2d(xp.asarray(max1))
    min2 = xp.atleast_2d(xp.asarray(min2))
    max2 = xp.atleast_2d(xp.asarray(max2))
    
    # For each axis: gap = max(0, max(min1, min2) - min(max1, max2))
    lower = xp.maximum(min1, min2)  # (N, 3)
    upper = xp.minimum(max1, max2)  # (N, 3)
    gaps = xp.maximum(0.0, lower - upper)  # (N, 3)
    
    # If separated: distance is L2 norm of gaps
    separated_dist = xp.sqrt(xp.sum(gaps * gaps, axis=-1))
    
    # If overlapping: penetration is smallest overlap
    overlaps = upper - lower  # (N, 3), positive when overlapping
    min_overlap = xp.min(overlaps, axis=-1)
    
    # Return separation distance if separated, else negative penetration
    is_separated = xp.any(gaps > 0, axis=-1)
    result = xp.where(is_separated, separated_dist, -min_overlap)
    
    return float(result[0]) if result.shape[0] == 1 else result


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

def analyzeManifoldProperties(manifold : m3d.Manifold):
  # Get mesh representation
  mesh = manifold.to_mesh()
  vertices = mesh.vert_properties[:, :3]
  triangles = mesh.tri_verts

  volume = manifold.volume()
  surface_area = manifold.surface_area()
  genus = manifold.genus()
  print(f"Volume: {volume}, Surface Area: {surface_area}, Genus: {genus}")
  print(f"Vertices: {vertices.shape}, Triangles: {triangles.shape}")

def saveManifold(manifold : m3d.Manifold, filename : str):
    """
    Export a manifold to a 3MF file for 3D printing or CAD applications.
    """
    mesh_data = manifold.to_mesh()
    vertices = mesh_data.vert_properties[:, :3]  # Get XYZ coordinates
    faces = mesh_data.tri_verts
    tri_mesh = trimesh.Trimesh(vertices=vertices, faces=faces)
    
    tri_mesh.export(filename)
