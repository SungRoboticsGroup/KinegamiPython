import sympy as sp

# Define symbols
phi1, theta1, d, phi2, theta2 = sp.symbols('phi1 theta1 d phi2 theta2')
theta1, theta2 = sp.symbols('theta1 theta2')
S1, S2, S3, S4, S5, S6 = sp.symbols('S1 S2 S3 S4 S5 S6')
C1, C2, C3, C4, C5, C6 = sp.symbols('C1 C2 C3 C4 C5 C6')

# Convert DH parameters to Dubins parameters
def convert_dh_to_dubins(params):
    phi1, theta1, d, phi2, theta2 = params
    dh_params = [phi1, sp.pi - theta1, d, phi2 - sp.pi, sp.pi - theta2]
    # Adjust to be non-negative
    dh_params = [p + 2 * sp.pi if p < 0 else p for p in dh_params]
    return dh_params

# Convert trigonometric functions to variables
def convert_to_trig(func):
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
def convert_trig_to_var(func):
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
def trig_and_back(func):
    return convert_trig_to_var(sp.simplify(convert_to_trig(func)))

# Remove solutions with negative d3 values
def remove_d3s(solution_sets_in):
    return [sol for sol in solution_sets_in if sol[2] > 0]

# Get 12PMat terms suppressed by C4
def get_12p_mat_terms_sup_c4(func):
    const = func.subs({d3cubxS5: 0, d3cubxC5: 0, d3cub: 0, d3sqxS5: 0, d3sqxC5: 0, d3sq: 0, d3xS5: 0, d3xC5: 0, d3: 0,
