from PyQt5.QtGui import QVector3D, QMatrix3x3
import math
import numpy as np

def compute_sphere_intersection(org, dir, cen, rad):
    # Ensure direction is normalized for numerical stability
    dir_len_sq = dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]
    if dir_len_sq == 0:
        return float('inf')
    
    # Normalize if not already normalized
    if abs(dir_len_sq - 1.0) > 0.001:
        dir_len = math.sqrt(dir_len_sq)
        dir = [dir[0]/dir_len, dir[1]/dir_len, dir[2]/dir_len]
    
    # Ray-sphere intersection with normalized direction (a = 1)
    oc = [org[0] - cen[0], org[1] - cen[1], org[2] - cen[2]]
    
    # For normalized direction: a = 1
    # b = 2 * dot(dir, oc)
    # c = dot(oc, oc) - r^2
    b = 2 * (dir[0] * oc[0] + dir[1] * oc[1] + dir[2] * oc[2])
    c = oc[0] * oc[0] + oc[1] * oc[1] + oc[2] * oc[2] - rad * rad

    discrim = b*b - 4*c  # Since a = 1, this is b*b - 4*1*c

    if discrim < 0:
        return float('inf')
    
    sqrt_discrim = math.sqrt(discrim)
    # For a = 1: t = (-b ± sqrt(discrim)) / 2
    t0 = (-b - sqrt_discrim) / 2
    
    if t0 > 0.001:  # Add small epsilon to avoid self-intersection
        return t0
    else:
        t1 = (-b + sqrt_discrim) / 2
        if t1 > 0.001:
            return t1
        else:
            return float('inf')
    
def compute_cylinder_intersection(org: QVector3D, dir: QVector3D, start: QVector3D, axis: QVector3D, rad, len): 
    # not bothering with computing endcaps because they're too small to matter anyway

    n = dir.normalized()
    a = axis.normalized()
    b = start - org

    n_cross_a = QVector3D.crossProduct(n, a)

    discrim = QVector3D.dotProduct(n_cross_a, n_cross_a) * rad * rad - QVector3D.dotProduct(a, a) * (QVector3D.dotProduct(b, n_cross_a) ** 2)

    if (discrim < -0.00001):
        return float('inf')
    else:
        d = (QVector3D.dotProduct(n_cross_a, QVector3D.crossProduct(b, a)) - math.sqrt(discrim)) / QVector3D.dotProduct(n_cross_a, n_cross_a)
        
        if (d < 0):
            d = (QVector3D.dotProduct(n_cross_a, QVector3D.crossProduct(b, a)) + math.sqrt(discrim)) / QVector3D.dotProduct(n_cross_a, n_cross_a)

        t = QVector3D.dotProduct(a, (n * d - b))
        if (t > 0 and t < len):
            return d
        else:
            return float('inf')
            
def compute_closest_point_on_axis(org, dir, center, axis):
    threshold = 0.001
    # Normalize direction vectors
    direction1 = dir.normalized()
    direction2 = axis.normalized()
    
    # Compute the cross product of the direction vectors
    cross_directions = QVector3D.crossProduct(direction1, direction2)
    cross_directions_norm = cross_directions.length()
    
    # If the cross product is zero, the lines are parallel
    if cross_directions_norm < threshold:
        raise ValueError("The lines are parallel and do not intersect.")
    
    # Compute the vector between the origins
    origin_diff = center - org
    
    # Compute the determinants
    det1 = QVector3D.dotProduct(origin_diff, QVector3D.crossProduct(direction2, cross_directions))
    det2 = QVector3D.dotProduct(origin_diff, QVector3D.crossProduct(direction1, cross_directions))
    
    # Compute the parameters for the points of closest approach
    t1 = det1 / cross_directions_norm**2
    t2 = det2 / cross_directions_norm**2
    
    # Compute the points of closest approach
    point1 = org + t1 * direction1
    point2 = center + t2 * direction2
    
    return point2
    
def compute_plane_intersection(org, dir, normal, point):
    threshold = 0.001
    direction = dir.normalized()
    n = normal.normalized()
    a = QVector3D.dotProduct(n, direction)

    if (abs(a) < threshold):
        return None
    
    b = point - org
    t = QVector3D.dotProduct(n, b) / a

    return org + t * direction

def compute_torus_intersection(org: QVector3D, dir: QVector3D, center: QVector3D, normal: QVector3D, major_radius: float, minor_radius: float):
    n = dir.normalized()
    o = org - center
    n_dot_n = QVector3D.dotProduct(n, n)
    o_dot_n = QVector3D.dotProduct(o, n)
    o_dot_o = QVector3D.dotProduct(o, o)
    n_dot_o = QVector3D.dotProduct(n, o)
    n_dot_normal = QVector3D.dotProduct(n, normal)
    o_dot_normal = QVector3D.dotProduct(o, normal)
    
    R = major_radius
    r = minor_radius
    
    A = n_dot_n * n_dot_n
    B = 4 * n_dot_n * n_dot_o
    C = 2 * n_dot_n * (o_dot_o - R * R - r * r) + 4 * n_dot_o * n_dot_o + 4 * R * R * n_dot_normal * n_dot_normal
    D = 4 * (o_dot_o - R * R - r * r) * n_dot_o + 8 * R * R * o_dot_normal * n_dot_normal
    E = (o_dot_o - R * R - r * r) * (o_dot_o - R * R - r * r) - 4 * R * R * (r * r - o_dot_normal * o_dot_normal)
    
    coeffs = [A, B, C, D, E]
    roots = np.roots(coeffs)
    
    real_roots = [root.real for root in roots if np.isreal(root) and root.real > 0]
    
    if not real_roots:
        return float('inf')
    
    return min(real_roots)