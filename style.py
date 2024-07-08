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

surfaceOpacityDefault=1.0
linkOpacityDefault=0.7
linkColorDefault = (0.5, 0.5, 0.5, linkOpacityDefault)  # black
jointColorDefault = (1.0, 0.0, 1.0, surfaceOpacityDefault)  # magenta
jointEdgeColorDefault= (0.0, 0.0, 0.55, 1.0)  # darkblue
groundPlaneColorDefault= (0.0, 0.0, 0.0, 1.0)  # black 

jointAxisScaleDefault=10
globalAxisScaleDefault=10
groundPlaneScaleDefault=10

blockDefault=True

# in AutoCAD Color Index 
# https://ezdxf.readthedocs.io/en/stable/concepts/aci.html#aci
mountainColorDefault=5 #174  #dark blue
valleyColorDefault=1 #red
boundaryColorDefault=2 #yellow
cutColorDefault=3 #101 #green 3