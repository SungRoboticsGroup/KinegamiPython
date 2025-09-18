from PathCSC import *

"""
# Far away example that seemed to be violating the 4r theorem but was just numerical issues
startPosition=np.array([ -9.08609321, -30.90958521,  10.29452237])
startDir=np.array([ 0.64074783, -0.44087435, -0.62854756])
endPosition=np.array([ -5.95660638, -19.23582546,   6.42155759])
endDir=np.array([-0.28272321, -0.89664979,  0.34071504])
tDirMag=np.array([ 0.06987554,  0.98920704, -0.12878994, 11.62502318])
"""

"""
# Closer example that we don't know if is solvable with bounded turn angles
tDirMag=np.array([ 0.78061932,  0.41322181, -0.46891492,  3.86815786]) 
startPosition=np.array([ 22.58430219, -48.95899821, -57.77015721])
startDir=np.array([ 0.302911  , -0.61927601, -0.72439088])
endPosition=np.array([ 25.60386095, -47.36059101, -59.58399416]) 
endDir=np.array([ 0.302911  , -0.61927601, -0.72439088])
"""

"""
# Example where it should output the empty path
tDirMag=np.array([ 2.34688900e-01,  9.12391973e-01, -3.35353556e-01,  2.27336326e-06])
r=1
startPosition=np.array([0.73138697, 1.09109731, 3.79261917])
startDir=np.array([ 0.24498644,  0.9524254 , -0.18129398])
endPosition=np.array([0.73138697, 1.09109731, 3.79261917])
endDir=np.array([ 0.24498644,  0.9524254 , -0.18129398])
"""


# An example where it should give only a segment
tDirMag=np.array([-0.17798786,  0.39034924, -0.90329828,  4.30277564])
r=1 
startPosition=np.array([-0.58780875,  2.97608382, -2.29387138]) 
startDir=np.array([-0.17798786,  0.39034924, -0.90329828])
endPosition=np.array([-1.3536506 ,  4.65566904, -6.18056122])
endDir=np.array([-0.17798786,  0.39034924, -0.90329828])

"""
exampleFromOpt = PathCSC(tDirMag=tDirMag, r=1,
                startPosition=startPosition, startDir=startDir,
                endPosition=endPosition, endDir=endDir,
                circle1sign=1, circle2sign=1)

exampleFromOpt.show()
"""

paths = solveCSC(r=1, startPosition=startPosition, startDir=startDir,
                          endPosition=endPosition, endDir=endDir)

for path in paths:
    print("theta1", path.theta1, "theta2", path.theta2, "error", path.error, "length", path.length)
    path.show(block=False)


shortestValid = shortestCSC(r=1, startPosition=startPosition, startDir=startDir,
                            endPosition=endPosition, endDir=endDir, turnAngleLimit=np.pi)

print("theta1", shortestValid.theta1, "theta2", shortestValid.theta2, "error", shortestValid.error, "length", shortestValid.length)
shortestValid.show()






