"""
Make a figure illustrating the vector fields that a CS path induces on a plane
"""
import numpy as np
from geometryHelpers import *
from PathCS import *
from PathCSC import *
from matplotlib import cm

startColor = 'lightskyblue'
endColor = 'seagreen'

nhat = np.array([0,0,1]) / np.sqrt(2)
midPlane = Plane(np.array([0,0,0]), nhat)
numPointsPerAxis = 16
planeScale=6
GridPoints = midPlane.grid(scale=planeScale, numPoints=numPointsPerAxis).reshape(-1,3)
r=1

sum3D = plt.figure().add_subplot(projection='3d')
dirs2D = plt.figure().add_subplot()
sum2D = plt.figure().add_subplot()
dirs3D = plt.figure().add_subplot(projection='3d')
dense2D = plt.figure().add_subplot()
dense3D = plt.figure().add_subplot(projection='3d')


def plotPathsAndDirs(ax2D, ax3D, Pose, Destinations, color, pathAlpha=0.2):
    Paths = [PathCS(Pose, point, r) for point in Destinations]
    for path in Paths:
        path.addToPlot(ax3D, showWholeCircle=False, startColor=color, 
                       showStart=False, pathColor=color, alpha=pathAlpha, 
                       cscBoundaryMarker='')
    Directions = np.array([path.segmentDir for path in Paths])
    ax2D.quiver(Destinations[:,0], Destinations[:,1], 
                Directions[:,0], Directions[:,1], color=color)
    ax3D.quiver(Destinations[:,0], Destinations[:,1], Destinations[:,2], 
                Directions[:,0], Directions[:,1], Directions[:,2], color=color)
    return Directions

def plotTorusAndCircle(ax2D, ax3D, Pose, color, colorMap):
    HornTorus(r, Pose.t, Pose.R[:,0]).addToPlot(ax3D, colorMap=cm.Greys, alpha=0.3)
    startProjection = np.array([0,0,0])
    startProjection[:-1] = Pose.t[:-1]
    TorusProjectionCircle = Circle3D(2*r, startProjection, np.array([0,0,1]))
    count=60
    CirclePoints = TorusProjectionCircle.interpolate(count=count)
    ax3D.plot(CirclePoints[:,0], 
            CirclePoints[:,1],
            CirclePoints[:,2],
            color=color)
    ax2D.plot(CirclePoints[:,0], CirclePoints[:,1], color=color)
    return TorusProjectionCircle


StartPose = SE3.Trans(1,0,-3)@SE3.Ry(-np.pi/3)
#plotPathsAndDirs(StartPose, grid.reshape(-1,3), startColor)
startCircle = plotTorusAndCircle(dirs2D, dirs3D, StartPose, startColor, cm.Reds)
StartCirclePoints = startCircle.interpolate(60)
dense2D.plot(StartCirclePoints[:,0], StartCirclePoints[:,1], color=startColor)
dense3D.plot(StartCirclePoints[:,0], StartCirclePoints[:,1], StartCirclePoints[:,2], color=startColor)

EndPose = SE3.Trans(-2,-3,2)@SE3.Rz(-5*np.pi/6)
#plotPathsAndDirs(EndPose, grid.reshape(-1,3), endColor)
endCircle = plotTorusAndCircle(dirs2D, dirs3D, EndPose, endColor, cm.Blues)
EndCirclePoints = endCircle.interpolate(60)
dense2D.plot(EndCirclePoints[:,0], EndCirclePoints[:,1], color=endColor)
dense3D.plot(EndCirclePoints[:,0], EndCirclePoints[:,1], EndCirclePoints[:,2], color=endColor)

count = 11
StartCircleSamplePoints = startCircle.interpolate(count)
EndCircleSamplePoints = endCircle.interpolate(count)
plotPathsAndDirs(dirs2D, dirs3D, StartPose, StartCircleSamplePoints, startColor)
plotPathsAndDirs(dirs2D, dirs3D, EndPose, EndCircleSamplePoints, endColor)
#plotPathsAndDirs(StartPose, EndCircleSamplePoints, startColor)
#plotPathsAndDirs(EndPose, StartCircleSamplePoints, endColor)


boundingBall = minBoundingBall(Ball(startCircle.c, startCircle.r),
                               Ball(endCircle.c, endCircle.r))
boundingCircle = Circle3D(boundingBall.r+0.1, boundingBall.c, np.array([0,0,1]))

BoundingCirclePoints = boundingCircle.interpolate(60)
dirs3D.plot(BoundingCirclePoints[:,0],
          BoundingCirclePoints[:,1],
          BoundingCirclePoints[:,2], color='black')
sum3D.plot(BoundingCirclePoints[:,0],
          BoundingCirclePoints[:,1],
          BoundingCirclePoints[:,2], color='black')
dirs2D.plot(BoundingCirclePoints[:,0], BoundingCirclePoints[:,1], color='black')
sum2D.plot(BoundingCirclePoints[:,0], BoundingCirclePoints[:,1], color='black')
dense2D.plot(BoundingCirclePoints[:,0], BoundingCirclePoints[:,1], color='black')
dense3D.plot(BoundingCirclePoints[:,0], BoundingCirclePoints[:,1], BoundingCirclePoints[:,2], color='black')

largerCircle = Circle3D(boundingCircle.r+1, boundingCircle.c, np.array([0,0,1]))
LargerCirclePoints = largerCircle.interpolate(60)
dirs2D.plot(LargerCirclePoints[:,0], LargerCirclePoints[:,1], color='black', alpha=0)
sum2D.plot(LargerCirclePoints[:,0], LargerCirclePoints[:,1], color='black', alpha=0)

"""
numInnerCircles = 4
radii = np.linspace(0, boundingCircle.r, numInnerCircles+2)[1:-1]
print(radii)
SmallerCirclePoints = np.zeros((numInnerCircles,count,3))
for i, smallerRadius in enumerate(radii):
    smallerCircle = Circle3D(smallerRadius, boundingCircle.c, np.array([0,0,1]))
    SmallerCirclePoints[i,:,:] = smallerCircle.interpolate(count)
SmallerCirclePoints = SmallerCirclePoints.reshape(-1,3)
StartToInnerDirs = plotPathsAndDirs(dense2D, dense3D, StartPose, SmallerCirclePoints, startColor)
EndToInnerDirs = plotPathsAndDirs(dense2D, dense3D, EndPose, SmallerCirclePoints, endColor)
InnerDirectionSum = StartToInnerDirs + EndToInnerDirs
dense3D.quiver(SmallerCirclePoints[:,0], SmallerCirclePoints[:,1], SmallerCirclePoints[:,2], InnerDirectionSum[:,0], InnerDirectionSum[:,1], InnerDirectionSum[:,2], color='black')
sum3D.quiver(SmallerCirclePoints[:,0], SmallerCirclePoints[:,1], SmallerCirclePoints[:,2], InnerDirectionSum[:,0], InnerDirectionSum[:,1], InnerDirectionSum[:,2], color='black')
"""
StartToGridDirs = plotPathsAndDirs(dense2D, dense3D, StartPose, GridPoints, startColor)
EndToGridDirs = plotPathsAndDirs(dense2D, dense3D, EndPose, GridPoints, endColor)
GridDirectionSum = StartToGridDirs + EndToGridDirs
dense3D.quiver(GridPoints[:,0], GridPoints[:,1], GridPoints[:,2], GridDirectionSum[:,0], GridDirectionSum[:,1], GridDirectionSum[:,2], color='black')
dense2D.quiver(GridPoints[:,0], GridPoints[:,1], GridDirectionSum[:,0], GridDirectionSum[:,1], color='black')
#sum3D.quiver(GridPoints[:,0], GridPoints[:,1], GridPoints[:,2], GridDirectionSum[:,0], GridDirectionSum[:,1], GridDirectionSum[:,2], color='black')
sum2D.quiver(GridPoints[:,0], GridPoints[:,1], GridDirectionSum[:,0], GridDirectionSum[:,1], color='black')


BoundingSample = boundingCircle.interpolate(count)
StartToBoundingDirs = plotPathsAndDirs(dirs2D, dirs3D, StartPose, BoundingSample, startColor)
EndToBoundingDirs = plotPathsAndDirs(dirs2D, dirs3D, EndPose, BoundingSample, endColor)
DirectionSum = StartToBoundingDirs + EndToBoundingDirs
#dirs3D.quiver(BoundingSample[:,0], BoundingSample[:,1], BoundingSample[:,2], DirectionSum[:,0], DirectionSum[:,1], DirectionSum[:,2], color='black')
sum3D.quiver(BoundingSample[:,0], BoundingSample[:,1], BoundingSample[:,2], DirectionSum[:,0], DirectionSum[:,1], DirectionSum[:,2], color='black')
sum2D.quiver(BoundingSample[:,0], BoundingSample[:,1],
            DirectionSum[:,0], DirectionSum[:,1], color='black')

solution = shortestCSC(r, StartPose.t, StartPose.R[:,0], 
        EndPose.t, -EndPose.R[:,0])
solution.addToPlot(dirs3D, pathColor='black', showCircles=True, 
                startColor=startColor, endColor=endColor, cscBoundaryMarker='')
solution.addToPlot(sum3D, pathColor='black', showCircles=False, 
                startColor=startColor, endColor=endColor, cscBoundaryMarker='')
solutionLine = Line(solution.turn1end, solution.tUnit)
zeroOfDirectionSum = midPlane.intersectionWithLine(solutionLine)
plotPathsAndDirs(dirs2D, dirs3D, StartPose, zeroOfDirectionSum.reshape(1,3), startColor)
plotPathsAndDirs(dirs2D, dirs3D, EndPose, zeroOfDirectionSum.reshape(1,3), endColor)
#plotPathsAndDirs(sum2D, sum3D, StartPose, zeroOfDirectionSum.reshape(1,3), startColor)
#plotPathsAndDirs(sum2D, sum3D, EndPose, zeroOfDirectionSum.reshape(1,3), endColor)
plotPathsAndDirs(dense2D, dense3D, StartPose, zeroOfDirectionSum.reshape(1,3), startColor)
plotPathsAndDirs(dense2D, dense3D, EndPose, zeroOfDirectionSum.reshape(1,3), endColor)

sum3D.plot([zeroOfDirectionSum[0]], [zeroOfDirectionSum[1]], [zeroOfDirectionSum[2]], marker='*', color='black')
sum2D.plot([zeroOfDirectionSum[0]], [zeroOfDirectionSum[1]], marker='*', color='black')
dirs3D.plot([zeroOfDirectionSum[0]], [zeroOfDirectionSum[1]], [zeroOfDirectionSum[2]], marker='*', color='black')
dirs2D.plot([zeroOfDirectionSum[0]], [zeroOfDirectionSum[1]], marker='*', color='black')


Plane(boundingCircle.c, nhat).addToPlot(dirs3D, scale=planeScale, color='black', alpha=0.05)
midPlane.addToPlot(sum3D, scale=planeScale, color='black', alpha=0.05)

dirs3D.set_aspect('equal')
dirs3D.axis('off')
dirs2D.set_aspect('equal')
dirs2D.axis('off')
sum3D.set_aspect('equal')
sum3D.axis('off')
sum2D.set_aspect('equal')
sum2D.axis('off')
dense2D.set_aspect('equal')
dense2D.axis('off')
dense3D.set_aspect('equal')
dense3D.axis('off')
plt.show()