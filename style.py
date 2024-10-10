import ezdxf

# Chosen to be mutually distinguishable in grayscale
xColorDefault = (1.0, 0.27, 0.0, 1.0)  # orangered
yColorDefault = (0.56, 0.93, 0.56, 1.0)  # lightgreen
zColorDefault = (0.0, 0.0, 0.55, 1.0)  # darkblue
proximalColorDefault = (0.0, 1.0, 1.0, 1.0)  # cyan
centerColorDefault = (1.0, 0.0, 1.0, 1.0)  # magenta
distalColorDefault = (1.0, 1.0, 0.0, 1.0)  # yellow 

pathColorDefault = (0.0, 0.0, 0.0, 1.0)  # black
sphereColorDefault = (0.0, 0.0, 0.0, 0.05)  # black, 0.05 opacity

surfaceOpacityDefault = 1.0
linkOpacityDefault = 0.5
linkColorDefault = (0.5, 0.5, 0.5, linkOpacityDefault)  # gray
jointColorDefault = (1.0, 0.0, 1.0, surfaceOpacityDefault)  # magenta 
jointEdgeColorDefault = (0.0, 0.0, 0.55, 1.0)  # darkblue 
groundPlaneColorDefault = (0.0, 0.0, 0.0, 1.0)  # black 

jointAxisScaleDefault = 10
globalAxisScaleDefault = 10
groundPlaneScaleDefault = 10

blockDefault = True

gridColorDefault = (0, 0, 0, 255)  # black
backgroundColorDefault = (255, 255, 255, 255)  # white
successColorDefault = "green"  # green
errorColorDefault = "red"  # red

meshItemDefaultColor = (1.0, 1.0, 1.0, 1.0)  # white
meshItemEdgeColor = (0.5, 0.5, 0.5, 1.0)  # gray
meshLinkDefaultColor = (1.0, 1.0, 1.0, 1.0)  # white
meshLinkEdgeColor = (0.5, 0.5, 0.5, 1.0)  # gray 

ballDefaultColor = (0, 0, 0, 0.5)  # black
cylinderColorList = (1, 0, 0, 1)  # red
elbowColorList = (1, 1, 1, 1)  # white
compoundElbowColorList = (1, 1, 1, 0.5)  # white, 0.5 opacity

xPoseColor = 'darkred'  # darkred
yPoseColor = 'darkblue'  # darkblue
zPoseColor = 'darkgreen'  # darkgreen
oColor = 'black'  # black 

lineColor = (1, 0, 0, 1)  # red
originPointsColor = (1, 1, 1, 1)  # white 
showAxisColor = (0.75, 0.75, 0.75, 1)  # light gray
lineSphereColor = [0, 0, 0, 0]  # black, fully transparent
rotateArrowColors = [(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1)]  # red, green, blue
selectedArrowColor = (1, 1, 1, 1)  # white 

extendedCircleColors = [(1, 0, 0, 0.5), (0, 1, 0, 0.5), (0, 0, 1, 0.5)]  # red, green, blue (0.5 opacity)

selectedLinkColor = (1.0, 1.0, 0.0, 0.5)  # yellow (0.5 opacity)
selectedJointColor = (1.0, 1.0, 0.5, 0.5)  # light yellow (0.5 opacity)

pathStartColor = 'r'  # red
pathEndColor = 'b'  # blue
pathWidgetColor = 'g'  # green
pathPointColor = (1, 1, 1, 1)  # white 
tUnitColor = (1, 1, 1, 1)  # white

# In AutoCAD Color Index 
# https://ezdxf.readthedocs.io/en/stable/concepts/aci.html#aci
mountainColorDefault = 5  # blue
valleyColorDefault = 1  # red 
mountainCColorDefault = 4  # cyan
valleyCColorDefault = 6  # magenta 
boundaryColorDefault = 2  # yellow
cutColorDefault = 3  # green
rawPropertiesColor = "#eaeaea" # light grey
mapPropertiesColor = "#ffffff" # white