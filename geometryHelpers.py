# -*- coding: utf-8 -*-
"""
Assorted geometry-related helper functions and classes 
"""

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
import pyqtgraph.opengl as gl
from matplotlib import cm

from style import *

class MeshItemWithID(gl.GLMeshItem):
    def __init__(self, id : int = -1, **kwds):
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
        
        super().__init__(**kwds)        
        self.id = id

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

def planeFromThreePoints(p1, p2, p3):
    normal = unit(cross(p2-p1, p3-p1))
    return Plane(p1, normal)

class Circle3D:
    def __init__(self, radius, center, normal):
        assert(norm(normal)>0)
        self.r = radius
        self.c = center
        self.n = normal / norm(normal)
    
    def interpolate(self, count=50):
        angle = np.linspace(0, 2*np.pi, count).reshape(-1,1)
        u = self.r * np.cos(angle)
        v = self.r * np.sin(angle)
        
        # construct basis for circle plane
        uhat = unitNormalToBoth(self.n, self.n).reshape(1,3)
        vhat = cross(self.n, uhat).reshape(1,3)
        
        # 3d circle points
        return self.c + u @ uhat + v @ vhat

class Ball:
    # closed ball centered at self.c of radius self.r
    def __init__(self, center, radius):
        self.c = center
        self.r = radius
    
    def containsPoint(self, point):
        return norm(self.c - point) <= self.r
        
    def addToWidget(self, widget, color=(0,0,0,0.5), is_waypoint=False, id=-1):
        md = gl.MeshData.sphere(rows=20, cols=20)
        sphere = MeshItemWithID(meshdata=md, color=tuple(color), shader='shaded', smooth=True, id=id)
        sphere.setGLOptions('translucent')
        sphere.scale(self.r, self.r, self.r)
        sphere.translate(*self.c)
        if (is_waypoint):
            sphere.setObjectName("Waypoint")
            sphere.setGLOptions('opaque')
            sphere.scale(self.r * 0.9, self.r * 0.9, self.r * 0.9)
        widget.plot_widget.addItem(sphere)

    def projectionOntoPlane(self, plane : Plane) -> Circle3D:
        return Circle3D(self.r, plane.projectionOfPoint(self.c), plane.nhat)
    
    def translationToCenterOnPlane(self, plane : Plane):
        return Ball(plane.projectionOfPoint(self.c), self.r)

class Cylinder:
    def __init__(self, radius : float, start : np.ndarray, 
                 direction : np.ndarray, length : float, 
                 uhat : np.ndarray = None):
        assert(length > 0)
        assert(norm(direction)>0)
        self.start = start
        self.direction = direction / norm(direction)
        self.length = length
        self.end = start + self.length * self.direction
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
        uhat = self.uhat
        vhat = cross(self.direction, uhat)
        circle = u.reshape(-1,1) @ uhat.reshape(1,3) + v.reshape(-1,1) @ vhat.reshape(1,3)
        
        segment = np.linspace(self.start, self.end, numCircles)
        circlePoints = np.tile(circle, (numCircles,1)) + np.repeat(segment, radialCount, axis=0)
        return circlePoints.reshape((numCircles, radialCount, 3))
        
    def interpolateQtCircles(self, numPointsPerCircle=32, numCircles=10):

        angle = np.linspace(0, 2 * np.pi, numPointsPerCircle, endpoint=False)
        u = self.r * np.cos(angle)
        v = self.r * np.sin(angle)
        
        """
        n = null_space([self.direction])
        uhat = n[:, 0]
        vhat = n[:, 1]
        """
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

    def addToWidget(self, widget, numPointsPerCircle=32, numCircles=10, color_list=(1, 0, 0, 1), is_joint=False, id=-1):
        vertices, indices = self.interpolateQtCircles(numPointsPerCircle, numCircles)
        meshdata = gl.MeshData(vertexes=vertices, faces=indices)
        meshitem = MeshItemWithID(meshdata=meshdata, color=tuple(color_list), shader='shaded', smooth=True, id=id)
        meshitem.setGLOptions('translucent')
        if (is_joint):
            meshitem.setObjectName("Joint")
        else:
            meshitem.setObjectName("Link")
        
        widget.plot_widget.addItem(meshitem)

        return vertices, indices
    
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
    
    def addToWidget(self, widget, numSides : int = 16, color_list=(1, 1, 1, 1), 
                  alpha : float = 1.0, wireFrame : bool = True, 
                  showFrames : bool = False, debug : bool = False):
        
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
                                        axisLength=1, xColor='darkred', 
                                        yColor='darkblue', zColor='darkgreen')

            x,y,z = Fwd.t
            u,v,w = np.cross(Fwd.R[:,0], FwdRot.R[:,0])
            line = gl.GLLinePlotItem(pos=np.array([[x,y,z], [0.5*u,0.5*v,0.5*w]]), color=(1, 0, 0, 1), width=5) 
            widget.plot_widget.addItem(line)

        meshitem.setGLOptions('translucent')
        widget.plot_widget.addItem(meshitem)

def addPosesToPlotQT(Poses, ax, widget, axisLength, xColor=xColorDefault, yColor=yColorDefault, 
                     zColor=zColorDefault, oColors='black', makeAxisLimitsIncludeTips=True):
    if Poses.shape == (4,4): # so it can plot a single frame
        Poses = np.array([Poses])

    ux, vx, wx = Poses[:,0:3,0].T # frame xhat coordinates
    uy, vy, wy = Poses[:,0:3,1].T # frame yhat coordinates
    uz, vz, wz = Poses[:,0:3,2].T # frame zhat coordinates
    ox, oy, oz = Poses[:,0:3,3].T # frame origin coordinates

    xPoints = [[ux[i], vx[i], wx[i]] for i in range(len(ux))]
    yPoints = [[uy[i], vy[i], wy[i]] for i in range(len(ux))]
    zPoints = [[uz[i], vz[i], wz[i]] for i in range(len(ux))]
    origins = [[ox[i], oy[i], oz[i]] for i in range(len(ox))]

    #plot origin points
    originPlot = gl.GLScatterPlotItem(pos=origins, color=(1,1,1,1), size=10)
    widget.addItem(originPlot)

    xPoints = np.column_stack((ox + ux, oy + vx, oz + wx))
    yPoints = np.column_stack((ox + uy, oy + vy, oz + wy))
    zPoints = np.column_stack((ox + uz, oy + vz, oz + wz))
    origins = np.column_stack((ox, oy, oz))

    for i in range(len(ux)):
        xHats = gl.GLLinePlotItem(pos=np.array([origins[i], xPoints[i]]), color=xColor, width=5)
        yHats = gl.GLLinePlotItem(pos=np.array([origins[i], yPoints[i]]), color=yColor, width=5)
        zHats = gl.GLLinePlotItem(pos=np.array([origins[i], zPoints[i]]), color=zColor, width=5)
        widget.addItem(xHats)
        widget.addItem(yHats)
        widget.addItem(zHats)

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
    
    def addToWidget(self, widget, numSides : int = 32, color_list=(1, 1, 1, 0.5), 
                  alpha : float = 1.0, wireFrame : bool = False, 
                  showFrames : bool = True, showBoundingBall : bool = False, debug : bool = False):
        allHandleSets = []
        
        for elbow in self.elbows:
            elbow.addToWidget(widget, numSides, color_list, alpha, wireFrame, showFrames, debug)
        
        if showBoundingBall:
            self.boundingBall().addToWidget(widget, color=tuple(color_list), alpha = 0.25*alpha)
        
        return allHandleSets
    
    def generateMesh(self, numSides : int = 32):
        vertices = []
        faces = []

        for elbow in self.elbows:
            v, f = elbow.circleEllipseCircleQT(numSides)

    def boundingBall(self):
        ball = self.elbows[0].boundingBall()
        for elbow in self.elbows[1:]:
            ball = minBoundingBall(ball, elbow.boundingBall())
        return ball
    
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
    
    def interpolate(self, count=50):
        angle = np.linspace(0, self.theta, count).reshape(-1,1)
        u = self.r * np.cos(angle)
        v = self.r * np.sin(angle)
        
        # construct basis for circle plane
        uhat = -self.startNormal.reshape(1,3)
        vhat = cross(self.binormal, uhat).reshape(1,3)
        
        # 3d circle points
        return self.circleCenter + u @ uhat + v @ vhat
        
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
