from PathCSC import *

"""
# Far away example that seemed to be violating the 4r theorem but was just numerical issues
startPosition=np.array([ -9.08609321, -30.90958521,  10.29452237])
startDir=np.array([ 0.64074783, -0.44087435, -0.62854756])
endPosition=np.array([ -5.95660638, -19.23582546,   6.42155759])
endDir=np.array([-0.28272321, -0.89664979,  0.34071504])
tDirMag=np.array([ 0.06987554,  0.98920704, -0.12878994, 11.62502318])
"""

# Closer example that we don't know if is solvable with bounded turn angles
tDirMag=np.array([ 0.78061932,  0.41322181, -0.46891492,  3.86815786]) 
startPosition=np.array([ 22.58430219, -48.95899821, -57.77015721])
startDir=np.array([ 0.302911  , -0.61927601, -0.72439088])
endPosition=np.array([ 25.60386095, -47.36059101, -59.58399416]) 
endDir=np.array([ 0.302911  , -0.61927601, -0.72439088])

exampleFromOpt = PathCSC(tDirMag=tDirMag, r=1,
                startPosition=startPosition, startDir=startDir,
                endPosition=endPosition, endDir=endDir,
                circle1sign=1, circle2sign=1)

exampleFromOpt.show()


shortestValid = shortestCSC(r=1, startPosition=startPosition, startDir=startDir,
                            endPosition=endPosition, endDir=endDir, turnAngleLimit=np.pi)

shortestValid.show()