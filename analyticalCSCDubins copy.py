import numpy as np
from sympy import simplify, sin, cos, tan, symbols

phi1, theta1, d, phi2, theta2 = sp.symbols('phi1 theta1 d phi2 theta2')
theta1, theta2 = sp.symbols('theta1 theta2')
S1, S2, S3, S4, S5, S6 = sp.symbols('S1 S2 S3 S4 S5 S6')
C1, C2, C3, C4, C5, C6 = sp.symbols('C1 C2 C3 C4 C5 C6')

##FROM FUNCTIONS FOR MATH

def convertDHToDubins(phi1, theta1, d, phi2, theta2):
    dh_param = [phi1, np.pi - theta1, d, phi2 - np.pi, np.pi - theta2]
    #adjust to be nonnegative
    dh_param = [(x + 2 * np.pi) if x < 0 else x for x in dh_param]
    return dh_param

# Convert trigonometric functions to variables
def convertToTrig(func):
    replacements = {
        C1: sp.cos(theta1), C2: sp.cos(theta2),
        C3: sp.cos(theta3), C4: sp.cos(theta4),
        C5: sp.cos(theta5), C6: sp.cos(theta6),
        S1: sp.sin(theta1), S2: sp.sin(theta2),
        S3: sp.sin(theta3), S4: sp.sin(theta4),
        S5: sp.sin(theta5), S6: sp.sin(theta6),
    }
    return func.subs(replacements)

# Convert variables to trigonometric functions
def convertTrigToVar(func):
    replacements = {
        sp.cos(theta1): C1, sp.cos(theta2): C2,
        sp.cos(theta3): C3, sp.cos(theta4): C4,
        sp.cos(theta5): C5, sp.cos(theta6): C6,
        sp.sin(theta1): S1, sp.sin(theta2): S2,
        sp.sin(theta3): S3, sp.sin(theta4): S4,
        sp.sin(theta5): S5, sp.sin(theta6): S6,
    }
    return func.subs(replacements)

# Simplify function by converting to and from trigonometric functions
def trigAndBack(func):
    return convertTrigToVar(simplify(convertToTrig(func)))

# Remove solutions with negative d3 values
def removeD3s(solutionSets):
    return [sol for sol in solutionSets if sol[2] > 0]

##FROM MATH

##FROM FUNCTIONS TO SOLVE DUBINS 3D

def solveDubins3D(x2, v2, PMat, QMat):
    

