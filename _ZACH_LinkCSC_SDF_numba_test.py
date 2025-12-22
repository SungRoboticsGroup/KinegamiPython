#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numba-Accelerated Ray Marching for LinkCSC SDF

Uses Numba's parallel JIT compilation to accelerate ray marching on the CPU.
This provides significant speedup while keeping the code readable and debuggable
(compared to raw CUDA kernels).

The SDF uses the same capped torus formula from Inigo Quilez, properly aligned
with the Arc3D midpoint for symmetry.

Creates 4 shaded views of a random Dubins path alongside matplotlib ground truth.

Requirements:
    pip install numba

@author: Zach (AI Assistant)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from spatialmath import SE3, SO3
from LinkCSC import LinkCSC
import time

try:
    from numba import njit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    print("WARNING: Numba not available. Install with: pip install numba")
    NUMBA_AVAILABLE = False


# =============================================================================
# Numba-JIT SDF Functions
# =============================================================================

@njit(fastmath=True)
def sd_capped_torus(px: float, py: float, pz: float, 
                    sin_half: float, cos_half: float,
                    major_r: float, minor_r: float) -> float:
    """
    Signed Distance Function for a capped torus in local coordinates.
    
    Reference: Inigo Quilez's capped torus SDF
    https://iquilezles.org/articles/distfunctions/
    
    The torus is centered at origin in XY plane. Arc is SYMMETRIC about X-axis,
    spanning from -halfAngle to +halfAngle. Uses abs(x) for symmetry.
    
    Parameters:
    -----------
    px, py, pz : float
        Point in local coordinates (x=radial at midpoint, y=tangent, z=binormal)
    sin_half, cos_half : float
        sin and cos of half the arc angle
    major_r : float
        Major radius (arc turning radius)
    minor_r : float
        Minor radius (tube radius)
    
    Returns:
    --------
    float
        Signed distance (negative inside, positive outside)
    """
    # Use abs(x) for symmetry about X-axis
    abs_px = abs(px)
    
    # Check if point's angle from +X axis is within arc span
    if cos_half * abs_px > sin_half * py:
        # Within angular span - project onto arc centerline
        k = sin_half * abs_px + cos_half * py
    else:
        # Outside angular span - closest point is at arc cap
        k = np.sqrt(abs_px * abs_px + py * py)
    
    # Distance to torus surface
    return np.sqrt(abs_px*abs_px + py*py + pz*pz + major_r*major_r - 2.0*major_r*k) - minor_r


@njit(fastmath=True)
def sd_capsule(px: float, py: float, pz: float,
               ax: float, ay: float, az: float,
               bx: float, by: float, bz: float,
               r: float) -> float:
    """
    Signed distance to a capsule (line segment with spherical caps).
    
    Parameters:
    -----------
    px, py, pz : float
        Query point
    ax, ay, az : float
        Capsule start point
    bx, by, bz : float
        Capsule end point
    r : float
        Capsule radius
    
    Returns:
    --------
    float
        Signed distance
    """
    # Vector from A to B
    bax = bx - ax
    bay = by - ay
    baz = bz - az
    
    # Vector from A to P
    pax = px - ax
    pay = py - ay
    paz = pz - az
    
    # Project P onto line AB, clamped to [0, 1]
    ba_dot_ba = bax*bax + bay*bay + baz*baz
    if ba_dot_ba < 1e-12:
        # Degenerate capsule - just a sphere at A
        return np.sqrt(pax*pax + pay*pay + paz*paz) - r
    
    t = (pax*bax + pay*bay + paz*baz) / ba_dot_ba
    t = max(0.0, min(1.0, t))
    
    # Closest point on segment
    cx = ax + t * bax
    cy = ay + t * bay
    cz = az + t * baz
    
    # Distance to closest point
    dx = px - cx
    dy = py - cy
    dz = pz - cz
    
    return np.sqrt(dx*dx + dy*dy + dz*dz) - r


@njit(fastmath=True)
def transform_to_local(px: float, py: float, pz: float,
                       cx: float, cy: float, cz: float,
                       r00: float, r01: float, r02: float,
                       r10: float, r11: float, r12: float,
                       r20: float, r21: float, r22: float) -> tuple:
    """
    Transform a world point to local arc coordinates.
    
    Parameters:
    -----------
    px, py, pz : float
        World point
    cx, cy, cz : float
        Arc circle center (local origin)
    r00...r22 : float
        World-to-local rotation matrix elements (row-major)
    
    Returns:
    --------
    tuple (lx, ly, lz)
        Point in local coordinates
    """
    # Translate to center origin
    dx = px - cx
    dy = py - cy
    dz = pz - cz
    
    # Apply rotation
    lx = r00*dx + r01*dy + r02*dz
    ly = r10*dx + r11*dy + r12*dz
    lz = r20*dx + r21*dy + r22*dz
    
    return lx, ly, lz


@njit(fastmath=True)
def sdf_arc(px: float, py: float, pz: float,
            center: np.ndarray, frame: np.ndarray,
            sin_half: float, cos_half: float,
            major_r: float, minor_r: float) -> float:
    """
    Compute SDF for an arc segment.
    
    Parameters:
    -----------
    px, py, pz : float
        World point
    center : np.ndarray
        Arc circle center (3,)
    frame : np.ndarray
        World-to-local rotation matrix (3, 3) row-major
    sin_half, cos_half : float
        sin/cos of half arc angle
    major_r : float
        Arc turning radius
    minor_r : float
        Tube radius
    
    Returns:
    --------
    float
        Signed distance
    """
    # Transform to local coordinates
    lx, ly, lz = transform_to_local(
        px, py, pz,
        center[0], center[1], center[2],
        frame[0, 0], frame[0, 1], frame[0, 2],
        frame[1, 0], frame[1, 1], frame[1, 2],
        frame[2, 0], frame[2, 1], frame[2, 2]
    )
    
    return sd_capped_torus(lx, ly, lz, sin_half, cos_half, major_r, minor_r)


@njit(fastmath=True)
def sdf_link(px: float, py: float, pz: float,
             tube_r: float,
             has_arc1: bool, arc1_center: np.ndarray, arc1_frame: np.ndarray,
             arc1_sin_half: float, arc1_cos_half: float, arc1_major_r: float,
             has_arc2: bool, arc2_center: np.ndarray, arc2_frame: np.ndarray,
             arc2_sin_half: float, arc2_cos_half: float, arc2_major_r: float,
             has_straight: bool, straight_a: np.ndarray, straight_b: np.ndarray) -> float:
    """
    Compute the full LinkCSC SDF (union of arc1, straight, arc2).
    
    Returns the minimum distance to any segment.
    """
    min_dist = 1e10  # Large value
    
    # Arc 1
    if has_arc1:
        d = sdf_arc(px, py, pz, arc1_center, arc1_frame,
                    arc1_sin_half, arc1_cos_half, arc1_major_r, tube_r)
        min_dist = min(min_dist, d)
    
    # Straight capsule
    if has_straight:
        d = sd_capsule(px, py, pz,
                       straight_a[0], straight_a[1], straight_a[2],
                       straight_b[0], straight_b[1], straight_b[2],
                       tube_r)
        min_dist = min(min_dist, d)
    
    # Arc 2
    if has_arc2:
        d = sdf_arc(px, py, pz, arc2_center, arc2_frame,
                    arc2_sin_half, arc2_cos_half, arc2_major_r, tube_r)
        min_dist = min(min_dist, d)
    
    return min_dist


# =============================================================================
# Numba-Parallel Ray Marching
# =============================================================================

@njit(parallel=True, fastmath=True)
def ray_march_parallel(origins: np.ndarray, directions: np.ndarray,
                       tube_r: float,
                       has_arc1: bool, arc1_center: np.ndarray, arc1_frame: np.ndarray,
                       arc1_sin_half: float, arc1_cos_half: float, arc1_major_r: float,
                       has_arc2: bool, arc2_center: np.ndarray, arc2_frame: np.ndarray,
                       arc2_sin_half: float, arc2_cos_half: float, arc2_major_r: float,
                       has_straight: bool, straight_a: np.ndarray, straight_b: np.ndarray,
                       max_steps: int, max_dist: float, epsilon: float) -> np.ndarray:
    """
    Parallel ray marching using Numba's prange.
    
    Parameters:
    -----------
    origins : np.ndarray
        Ray origins (N, 3)
    directions : np.ndarray
        Ray directions (N, 3), normalized
    ... SDF parameters ...
    max_steps : int
        Maximum ray marching iterations
    max_dist : float
        Maximum ray travel distance
    epsilon : float
        Surface intersection threshold
    
    Returns:
    --------
    np.ndarray
        Hit distances (N,), negative if no hit
    """
    num_rays = origins.shape[0]
    hit_distances = np.empty(num_rays, dtype=np.float64)
    
    for i in prange(num_rays):
        ox, oy, oz = origins[i, 0], origins[i, 1], origins[i, 2]
        dx, dy, dz = directions[i, 0], directions[i, 1], directions[i, 2]
        
        t = 0.0
        hit = False
        
        for _ in range(max_steps):
            # Current position
            px = ox + t * dx
            py = oy + t * dy
            pz = oz + t * dz
            
            # Evaluate SDF
            d = sdf_link(px, py, pz, tube_r,
                         has_arc1, arc1_center, arc1_frame,
                         arc1_sin_half, arc1_cos_half, arc1_major_r,
                         has_arc2, arc2_center, arc2_frame,
                         arc2_sin_half, arc2_cos_half, arc2_major_r,
                         has_straight, straight_a, straight_b)
            
            if d < epsilon:
                hit = True
                break
            
            # Step forward (with minimum step size for stability)
            t += max(d, 0.001)
            
            if t > max_dist:
                break
        
        hit_distances[i] = t if hit else -1.0
    
    return hit_distances


@njit(parallel=True, fastmath=True)
def compute_normals_parallel(hit_points: np.ndarray, hit_mask: np.ndarray,
                             tube_r: float,
                             has_arc1: bool, arc1_center: np.ndarray, arc1_frame: np.ndarray,
                             arc1_sin_half: float, arc1_cos_half: float, arc1_major_r: float,
                             has_arc2: bool, arc2_center: np.ndarray, arc2_frame: np.ndarray,
                             arc2_sin_half: float, arc2_cos_half: float, arc2_major_r: float,
                             has_straight: bool, straight_a: np.ndarray, straight_b: np.ndarray,
                             delta: float) -> np.ndarray:
    """
    Compute surface normals via central differences gradient.
    """
    num_points = hit_points.shape[0]
    normals = np.zeros((num_points, 3), dtype=np.float64)
    
    for i in prange(num_points):
        if not hit_mask[i]:
            continue
        
        px, py, pz = hit_points[i, 0], hit_points[i, 1], hit_points[i, 2]
        
        # Central differences for gradient
        dx_pos = sdf_link(px + delta, py, pz, tube_r,
                          has_arc1, arc1_center, arc1_frame,
                          arc1_sin_half, arc1_cos_half, arc1_major_r,
                          has_arc2, arc2_center, arc2_frame,
                          arc2_sin_half, arc2_cos_half, arc2_major_r,
                          has_straight, straight_a, straight_b)
        dx_neg = sdf_link(px - delta, py, pz, tube_r,
                          has_arc1, arc1_center, arc1_frame,
                          arc1_sin_half, arc1_cos_half, arc1_major_r,
                          has_arc2, arc2_center, arc2_frame,
                          arc2_sin_half, arc2_cos_half, arc2_major_r,
                          has_straight, straight_a, straight_b)
        
        dy_pos = sdf_link(px, py + delta, pz, tube_r,
                          has_arc1, arc1_center, arc1_frame,
                          arc1_sin_half, arc1_cos_half, arc1_major_r,
                          has_arc2, arc2_center, arc2_frame,
                          arc2_sin_half, arc2_cos_half, arc2_major_r,
                          has_straight, straight_a, straight_b)
        dy_neg = sdf_link(px, py - delta, pz, tube_r,
                          has_arc1, arc1_center, arc1_frame,
                          arc1_sin_half, arc1_cos_half, arc1_major_r,
                          has_arc2, arc2_center, arc2_frame,
                          arc2_sin_half, arc2_cos_half, arc2_major_r,
                          has_straight, straight_a, straight_b)
        
        dz_pos = sdf_link(px, py, pz + delta, tube_r,
                          has_arc1, arc1_center, arc1_frame,
                          arc1_sin_half, arc1_cos_half, arc1_major_r,
                          has_arc2, arc2_center, arc2_frame,
                          arc2_sin_half, arc2_cos_half, arc2_major_r,
                          has_straight, straight_a, straight_b)
        dz_neg = sdf_link(px, py, pz - delta, tube_r,
                          has_arc1, arc1_center, arc1_frame,
                          arc1_sin_half, arc1_cos_half, arc1_major_r,
                          has_arc2, arc2_center, arc2_frame,
                          arc2_sin_half, arc2_cos_half, arc2_major_r,
                          has_straight, straight_a, straight_b)
        
        nx = dx_pos - dx_neg
        ny = dy_pos - dy_neg
        nz = dz_pos - dz_neg
        
        # Normalize
        length = np.sqrt(nx*nx + ny*ny + nz*nz)
        if length > 1e-10:
            normals[i, 0] = nx / length
            normals[i, 1] = ny / length
            normals[i, 2] = nz / length
    
    return normals


# =============================================================================
# Utility Functions
# =============================================================================

def random_unit_vector():
    """Generate a random unit vector in 3D"""
    vec = np.random.randn(3)
    return vec / np.linalg.norm(vec)


def random_SE3_pose():
    """Generate a random SE3 pose with random position and orientation"""
    position = np.random.uniform(-5, 5, 3)
    axis = random_unit_vector()
    angle = np.random.uniform(0, 2*np.pi)
    rotation = SO3.AngleAxis(angle, axis)
    return SE3.Rt(rotation, position)


def create_random_link(r=1.0, min_separation=3.0, max_attempts=10):
    """Create a random LinkCSC with reasonable parameters"""
    for attempt in range(max_attempts):
        try:
            start_pose = random_SE3_pose()
            separation = np.random.uniform(min_separation, min_separation + 3)
            end_position = start_pose.t + random_unit_vector() * separation
            end_direction = random_unit_vector()
            
            v1 = end_direction
            v2_temp = random_unit_vector()
            v2 = v2_temp - np.dot(v2_temp, v1) * v1
            v2 = v2 / np.linalg.norm(v2)
            v3 = np.cross(v1, v2)
            
            R_end = np.column_stack([v1, v2, v3])
            end_pose = SE3.Rt(SO3(R_end), end_position)
            
            link = LinkCSC(r, start_pose, end_pose, maxAnglePerElbow=np.pi/3)
            return link
        except Exception as e:
            if attempt == max_attempts - 1:
                raise RuntimeError(f"Failed to create valid link: {e}")
            continue
    raise RuntimeError("Failed to create valid link")


def get_link_bounds(link, padding=3.0):
    """Get bounding box of the link"""
    points = np.array([
        link.StartDubinsPose.t,
        link.EndDubinsPose.t,
        link.path.turn1end,
        link.path.turn1end + link.path.tMag * link.path.tUnit
    ])
    
    if link.arc1:
        points = np.vstack([points, link.path.circleCenter1])
    if link.arc2:
        points = np.vstack([points, link.path.circleCenter2])
    
    center = points.mean(axis=0)
    extent = np.max(np.linalg.norm(points - center, axis=1)) + padding
    
    return center, extent


def prepare_link_data(link):
    """
    Prepare link geometry data for Numba functions.
    
    The arc local frame has:
    - X-axis: radial outward at arc MIDPOINT (for abs(x) symmetry in SDF)
    - Y-axis: tangent at midpoint
    - Z-axis: binormal
    
    Returns dict with numpy arrays for all parameters.
    """
    from scipy.spatial.transform import Rotation
    
    data = {
        'tube_r': link.r,
        'has_arc1': link.arc1 is not None and link.path.theta1 > link.EPSILON,
        'has_arc2': link.arc2 is not None and link.path.theta2 > link.EPSILON,
        'has_straight': link.path.tMag > link.DISTANCE_EPSILON,
    }
    
    # Default placeholders (will be ignored if has_arc is False)
    data['arc1_center'] = np.zeros(3)
    data['arc1_frame'] = np.eye(3)
    data['arc1_sin_half'] = 0.0
    data['arc1_cos_half'] = 1.0
    data['arc1_major_r'] = 1.0
    
    data['arc2_center'] = np.zeros(3)
    data['arc2_frame'] = np.eye(3)
    data['arc2_sin_half'] = 0.0
    data['arc2_cos_half'] = 1.0
    data['arc2_major_r'] = 1.0
    
    data['straight_a'] = np.zeros(3)
    data['straight_b'] = np.zeros(3)
    
    # Arc 1 data
    if data['has_arc1']:
        data['arc1_center'] = link.path.circleCenter1.copy()
        data['arc1_major_r'] = link.arc1.r
        
        half_theta = link.arc1.theta / 2.0
        data['arc1_sin_half'] = np.sin(half_theta)
        data['arc1_cos_half'] = np.cos(half_theta)
        
        # Build local frame: X at midpoint, Z is binormal
        halfAngleRot = Rotation.from_rotvec(half_theta * link.arc1.binormal)
        centerToMid = halfAngleRot.apply(-link.arc1.startNormal)
        arcX = centerToMid / np.linalg.norm(centerToMid)
        arcZ = link.arc1.binormal
        arcY = np.cross(arcZ, arcX)
        
        # World-to-local rotation (rows are local axes in world coords)
        data['arc1_frame'] = np.array([arcX, arcY, arcZ])
    
    # Arc 2 data
    if data['has_arc2']:
        data['arc2_center'] = link.path.circleCenter2.copy()
        data['arc2_major_r'] = link.arc2.r
        
        half_theta = link.arc2.theta / 2.0
        data['arc2_sin_half'] = np.sin(half_theta)
        data['arc2_cos_half'] = np.cos(half_theta)
        
        # Build local frame: X at midpoint, Z is binormal
        halfAngleRot = Rotation.from_rotvec(half_theta * link.arc2.binormal)
        centerToMid = halfAngleRot.apply(-link.arc2.startNormal)
        arcX = centerToMid / np.linalg.norm(centerToMid)
        arcZ = link.arc2.binormal
        arcY = np.cross(arcZ, arcX)
        
        # World-to-local rotation (rows are local axes in world coords)
        data['arc2_frame'] = np.array([arcX, arcY, arcZ])
    
    # Straight segment data
    if data['has_straight']:
        data['straight_a'] = link.path.turn1end.copy()
        data['straight_b'] = (link.path.turn1end + link.path.tMag * link.path.tUnit).copy()
    
    return data


def create_camera(center, distance, azimuth, elevation, resolution=(400, 300)):
    """
    Create camera rays for ray marching.
    
    Parameters:
    -----------
    center : np.ndarray
        Point the camera is looking at
    distance : float
        Distance from center
    azimuth : float
        Horizontal angle (radians)
    elevation : float
        Vertical angle (radians)
    resolution : tuple
        Image resolution (width, height)
        
    Returns:
    --------
    origins : np.ndarray (N, 3)
    directions : np.ndarray (N, 3)
    shape : tuple (width, height)
    """
    width, height = resolution
    
    # Camera position
    cam_x = distance * np.cos(elevation) * np.cos(azimuth)
    cam_y = distance * np.cos(elevation) * np.sin(azimuth)
    cam_z = distance * np.sin(elevation)
    cam_pos = center + np.array([cam_x, cam_y, cam_z])
    
    # Camera orientation
    forward = center - cam_pos
    forward = forward / np.linalg.norm(forward)
    
    # World up vector
    world_up = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(forward, world_up)) > 0.99:
        world_up = np.array([1.0, 0.0, 0.0])
    
    right = np.cross(forward, world_up)
    right = right / np.linalg.norm(right)
    
    up = np.cross(right, forward)
    up = up / np.linalg.norm(up)
    
    # Field of view
    fov = 50  # degrees
    aspect = width / height
    
    # Create ray directions
    fov_rad = np.deg2rad(fov)
    tan_fov = np.tan(fov_rad / 2)
    
    # Pixel coordinates (centered at 0)
    px = np.linspace(-1, 1, width)
    py = np.linspace(-1, 1, height)
    px, py = np.meshgrid(px, py)
    
    # Ray directions in camera space
    px = px.flatten() * tan_fov * aspect
    py = py.flatten() * tan_fov
    
    # Transform to world space
    directions = (forward[np.newaxis, :] + 
                 px[:, np.newaxis] * right[np.newaxis, :] + 
                 py[:, np.newaxis] * up[np.newaxis, :])
    
    # Normalize directions
    directions = directions / np.linalg.norm(directions, axis=1, keepdims=True)
    
    # All rays start from camera position
    origins = np.tile(cam_pos, (len(directions), 1))
    
    return origins, directions, (width, height)


def render_view(link_data, center, extent, azimuth, elevation,
                resolution=(400, 300), max_steps=150, epsilon=0.002):
    """
    Render a single view using Numba ray marching.
    
    Returns:
    --------
    image : np.ndarray (H, W, 3) RGB image
    """
    # Create camera
    distance = extent * 2.5
    origins, directions, (width, height) = create_camera(
        center, distance, azimuth, elevation, resolution
    )
    
    # Ray march
    hit_distances = ray_march_parallel(
        origins, directions,
        link_data['tube_r'],
        link_data['has_arc1'], link_data['arc1_center'], link_data['arc1_frame'],
        link_data['arc1_sin_half'], link_data['arc1_cos_half'], link_data['arc1_major_r'],
        link_data['has_arc2'], link_data['arc2_center'], link_data['arc2_frame'],
        link_data['arc2_sin_half'], link_data['arc2_cos_half'], link_data['arc2_major_r'],
        link_data['has_straight'], link_data['straight_a'], link_data['straight_b'],
        max_steps, distance * 2, epsilon
    )
    
    # Compute hit points
    hit_mask = hit_distances > 0
    hit_points = origins + hit_distances[:, np.newaxis] * directions
    
    # Compute normals for shading
    normals = compute_normals_parallel(
        hit_points, hit_mask,
        link_data['tube_r'],
        link_data['has_arc1'], link_data['arc1_center'], link_data['arc1_frame'],
        link_data['arc1_sin_half'], link_data['arc1_cos_half'], link_data['arc1_major_r'],
        link_data['has_arc2'], link_data['arc2_center'], link_data['arc2_frame'],
        link_data['arc2_sin_half'], link_data['arc2_cos_half'], link_data['arc2_major_r'],
        link_data['has_straight'], link_data['straight_a'], link_data['straight_b'],
        delta=0.001
    )
    
    # Lambert shading
    # Light direction (from camera roughly)
    cam_x = extent * 2.5 * np.cos(elevation) * np.cos(azimuth)
    cam_y = extent * 2.5 * np.cos(elevation) * np.sin(azimuth)
    cam_z = extent * 2.5 * np.sin(elevation)
    light_dir = np.array([cam_x, cam_y, cam_z])
    light_dir = light_dir / np.linalg.norm(light_dir)
    
    # Compute brightness
    brightness = np.maximum(0, np.sum(normals * light_dir, axis=1))
    
    # Add ambient
    ambient = 0.3
    brightness = ambient + (1 - ambient) * brightness
    
    # Apply to hit pixels
    image = np.zeros((height * width, 3))
    image[hit_mask] = brightness[hit_mask, np.newaxis] * np.array([0.8, 0.6, 0.4])  # Warm color
    
    # Background color for misses
    image[~hit_mask] = [0.1, 0.1, 0.15]  # Dark blue-gray
    
    # Reshape to image
    image = image.reshape(height, width, 3)
    
    # Flip vertically (camera Y points up)
    image = image[::-1, :, :]
    
    return image


def visualize_with_ground_truth(link, resolution=(400, 300)):
    """
    Create a figure with 4 ray-marched views and matplotlib 3D ground truth.
    """
    print("Preparing link data for Numba...")
    link_data = prepare_link_data(link)
    center, extent = get_link_bounds(link)
    
    # View angles: front, side, top, isometric
    views = [
        (0, 0, "Front"),
        (np.pi/2, 0, "Side"),
        (0, np.pi/2 - 0.01, "Top"),  # Slightly off vertical to avoid gimbal
        (np.pi/4, np.pi/6, "Iso"),
    ]
    
    # Create figure with GridSpec
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1, 1, 1.2])
    
    # Render ray-marched views
    print("Rendering ray-marched views (first call includes JIT compilation)...")
    
    for idx, (azimuth, elevation, title) in enumerate(views):
        row = idx // 2
        col = idx % 2
        
        t0 = time.time()
        image = render_view(link_data, center, extent, azimuth, elevation, resolution)
        t1 = time.time()
        print(f"  {title}: {t1-t0:.2f}s")
        
        ax = fig.add_subplot(gs[row, col])
        ax.imshow(image)
        ax.set_title(f"SDF Ray March - {title}")
        ax.axis('off')
    
    # Ground truth matplotlib 3D view
    print("Creating matplotlib ground truth...")
    ax3d = fig.add_subplot(gs[:, 2], projection='3d')
    link.addToPlot(ax3d, numSides=32, color='coral', alpha=0.7, wireFrame=False,
                   showPath=True, pathColor='blue')
    ax3d.set_title("Ground Truth (Matplotlib 3D)")
    ax3d.set_aspect('equal')
    
    # Set 3D view to match isometric
    ax3d.view_init(elev=30, azim=45)
    ax3d.axis('off')
    
    plt.tight_layout()
    plt.savefig('numba_sdf_render.png', dpi=150, bbox_inches='tight')
    print("Saved to numba_sdf_render.png")
    plt.show()


def test_sdf_values(link):
    """Quick test to verify SDF values at key points."""
    print("\n=== SDF Value Test ===")
    link_data = prepare_link_data(link)
    
    # Test at centerline points (should be -radius)
    test_points = {
        'start': link.StartDubinsPose.t,
        'end': link.EndDubinsPose.t,
        'turn1end': link.path.turn1end,
        'turn2start': link.path.turn1end + link.path.tMag * link.path.tUnit,
    }
    
    if link.arc1:
        mid1 = link.arc1.interpolateAt(0.5)
        test_points['arc1_mid'] = mid1
    
    tube_r = link.r
    
    print(f"Tube radius: {tube_r}")
    print(f"Expected SDF at centerline: -{tube_r}")
    print()
    
    for name, pt in test_points.items():
        # Use the Numba SDF
        sdf_val = sdf_link(pt[0], pt[1], pt[2],
                          link_data['tube_r'],
                          link_data['has_arc1'], link_data['arc1_center'], link_data['arc1_frame'],
                          link_data['arc1_sin_half'], link_data['arc1_cos_half'], link_data['arc1_major_r'],
                          link_data['has_arc2'], link_data['arc2_center'], link_data['arc2_frame'],
                          link_data['arc2_sin_half'], link_data['arc2_cos_half'], link_data['arc2_major_r'],
                          link_data['has_straight'], link_data['straight_a'], link_data['straight_b'])
        
        # Compare with CPU class method
        cpu_val = link.sdf(pt)
        
        diff = abs(sdf_val - cpu_val)
        status = "✓" if diff < 0.001 else "✗"
        print(f"  {name:15s}: Numba={sdf_val:8.4f}, CPU={cpu_val:8.4f}, diff={diff:.6f} {status}")
    
    print()


def main():
    """Main entry point."""
    if not NUMBA_AVAILABLE:
        print("Numba is required. Install with: pip install numba")
        return
    
    print("=" * 60)
    print("Numba-Accelerated LinkCSC SDF Ray Marching")
    print("=" * 60)
    
    # Create random link
    np.random.seed(42)  # For reproducibility
    print("\nCreating random Dubins path link...")
    link = create_random_link(r=1.0, min_separation=4.0)
    
    print(f"  Start: {link.StartDubinsPose.t}")
    print(f"  End:   {link.EndDubinsPose.t}")
    print(f"  Arc1 theta: {np.rad2deg(link.path.theta1):.1f}°" if link.arc1 else "  Arc1: None")
    print(f"  Arc2 theta: {np.rad2deg(link.path.theta2):.1f}°" if link.arc2 else "  Arc2: None")
    print(f"  Straight length: {link.path.tMag:.2f}")
    
    # Test SDF values
    test_sdf_values(link)
    
    # Visualize with ground truth comparison
    print("\n" + "=" * 60)
    visualize_with_ground_truth(link, resolution=(400, 300))


if __name__ == "__main__":
    main()
