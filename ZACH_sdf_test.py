#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CPU Ray Marching Test for LinkCSC SDF

Uses Numba for parallelization while calling the exact same SDF logic
from the Arc3D and LinkCSC classes. This allows visualization to verify
the SDF is correct.

@author: Zach (AI Assistant)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from spatialmath import SE3, SO3
from LinkCSC import LinkCSC
from scipy.spatial.transform import Rotation
import time

try:
    from numba import njit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    print("Numba not available, using pure Python (slow)")
    NUMBA_AVAILABLE = False
    def njit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator
    prange = range


# =============================================================================
# Numba-JIT SDF Functions (matching Arc3D logic exactly)
# =============================================================================

@njit(fastmath=True)
def sd_flat_ended_torus(px: float, py: float, pz: float, 
                        sin_half: float, cos_half: float,
                        major_r: float, minor_r: float) -> float:
    """
    Flat-ended torus SDF in local coordinates.
    Arc is symmetric about Y-axis, with flat disc ends instead of spherical caps.
    """
    abs_px = abs(px)
    
    # Endpoint center on the torus ring
    endCenter_x = major_r * sin_half
    endCenter_y = major_r * cos_half
    
    # Tangent at endpoint (points past the arc end)
    tangent_x = cos_half
    tangent_y = -sin_half
    
    # How far past the endpoint are we?
    toPoint_x = abs_px - endCenter_x
    toPoint_y = py - endCenter_y
    pastEnd = toPoint_x * tangent_x + toPoint_y * tangent_y
    
    if pastEnd <= 0.0:
        # Inside or at arc span - standard torus formula
        if cos_half * abs_px > sin_half * py:
            k = sin_half * abs_px + cos_half * py
        else:
            k = np.sqrt(abs_px * abs_px + py * py)
        return np.sqrt(abs_px*abs_px + py*py + pz*pz + major_r*major_r - 2.0*major_r*k) - minor_r
    else:
        # Past arc endpoint - distance to flat disc
        # Radial distance from tube axis (in the plane of the disc)
        radialInPlane = toPoint_x * sin_half + toPoint_y * cos_half
        discDist = np.sqrt(radialInPlane * radialInPlane + pz * pz)
        
        # 2D SDF to a disc
        outsideDisc = max(discDist - minor_r, 0.0)
        return np.sqrt(pastEnd * pastEnd + outsideDisc * outsideDisc)


@njit(fastmath=True)
def sd_capsule(px: float, py: float, pz: float,
               ax: float, ay: float, az: float,
               bx: float, by: float, bz: float,
               r: float) -> float:
    """Capsule SDF (line segment with spherical caps)."""
    bax, bay, baz = bx - ax, by - ay, bz - az
    pax, pay, paz = px - ax, py - ay, pz - az
    
    ba_dot_ba = bax*bax + bay*bay + baz*baz
    if ba_dot_ba < 1e-12:
        return np.sqrt(pax*pax + pay*pay + paz*paz) - r
    
    t = (pax*bax + pay*bay + paz*baz) / ba_dot_ba
    t = max(0.0, min(1.0, t))
    
    dx = pax - t * bax
    dy = pay - t * bay
    dz = paz - t * baz
    
    return np.sqrt(dx*dx + dy*dy + dz*dz) - r


@njit(fastmath=True)
def sdf_arc(px: float, py: float, pz: float,
            cx: float, cy: float, cz: float,
            frame: np.ndarray,
            sin_half: float, cos_half: float,
            major_r: float, minor_r: float) -> float:
    """Arc SDF: transform to local coords, then flat-ended torus."""
    # Translate to center
    dx, dy, dz = px - cx, py - cy, pz - cz
    
    # Rotate to local (frame is 3x3, row-major)
    lx = frame[0, 0]*dx + frame[0, 1]*dy + frame[0, 2]*dz
    ly = frame[1, 0]*dx + frame[1, 1]*dy + frame[1, 2]*dz
    lz = frame[2, 0]*dx + frame[2, 1]*dy + frame[2, 2]*dz
    
    return sd_flat_ended_torus(lx, ly, lz, sin_half, cos_half, major_r, minor_r)


@njit(fastmath=True)
def sdf_link(px: float, py: float, pz: float,
             tube_r: float,
             has_arc1: bool, arc1_cx: float, arc1_cy: float, arc1_cz: float,
             arc1_frame: np.ndarray, arc1_sin: float, arc1_cos: float, arc1_r: float,
             has_arc2: bool, arc2_cx: float, arc2_cy: float, arc2_cz: float,
             arc2_frame: np.ndarray, arc2_sin: float, arc2_cos: float, arc2_r: float,
             has_straight: bool,
             sa_x: float, sa_y: float, sa_z: float,
             sb_x: float, sb_y: float, sb_z: float,
             start_x: float, start_y: float, start_z: float,
             end_x: float, end_y: float, end_z: float) -> float:
    """Full LinkCSC SDF (union of arc1, straight, arc2)."""
    min_dist = 1e10
    
    if has_arc1:
        d = sdf_arc(px, py, pz, arc1_cx, arc1_cy, arc1_cz, arc1_frame,
                    arc1_sin, arc1_cos, arc1_r, tube_r)
        min_dist = min(min_dist, d)
    else:
        # Sphere at start
        d = np.sqrt((px-start_x)**2 + (py-start_y)**2 + (pz-start_z)**2) - tube_r
        min_dist = min(min_dist, d)
    
    if has_straight:
        d = sd_capsule(px, py, pz, sa_x, sa_y, sa_z, sb_x, sb_y, sb_z, tube_r)
        min_dist = min(min_dist, d)
    
    if has_arc2:
        d = sdf_arc(px, py, pz, arc2_cx, arc2_cy, arc2_cz, arc2_frame,
                    arc2_sin, arc2_cos, arc2_r, tube_r)
        min_dist = min(min_dist, d)
    else:
        # Sphere at end
        d = np.sqrt((px-end_x)**2 + (py-end_y)**2 + (pz-end_z)**2) - tube_r
        min_dist = min(min_dist, d)
    
    return min_dist


# =============================================================================
# Ray Marching
# =============================================================================

@njit(parallel=True, fastmath=True)
def ray_march_parallel(origins: np.ndarray, directions: np.ndarray,
                       tube_r: float,
                       has_arc1: bool, arc1_cx: float, arc1_cy: float, arc1_cz: float,
                       arc1_frame: np.ndarray, arc1_sin: float, arc1_cos: float, arc1_r: float,
                       has_arc2: bool, arc2_cx: float, arc2_cy: float, arc2_cz: float,
                       arc2_frame: np.ndarray, arc2_sin: float, arc2_cos: float, arc2_r: float,
                       has_straight: bool,
                       sa_x: float, sa_y: float, sa_z: float,
                       sb_x: float, sb_y: float, sb_z: float,
                       start_x: float, start_y: float, start_z: float,
                       end_x: float, end_y: float, end_z: float,
                       max_steps: int, max_dist: float, epsilon: float) -> np.ndarray:
    """Parallel ray marching."""
    num_rays = origins.shape[0]
    hit_distances = np.empty(num_rays, dtype=np.float64)
    
    for i in prange(num_rays):
        ox, oy, oz = origins[i, 0], origins[i, 1], origins[i, 2]
        dx, dy, dz = directions[i, 0], directions[i, 1], directions[i, 2]
        
        t = 0.0
        hit = False
        
        for _ in range(max_steps):
            px = ox + t * dx
            py = oy + t * dy
            pz = oz + t * dz
            
            d = sdf_link(px, py, pz, tube_r,
                         has_arc1, arc1_cx, arc1_cy, arc1_cz,
                         arc1_frame, arc1_sin, arc1_cos, arc1_r,
                         has_arc2, arc2_cx, arc2_cy, arc2_cz,
                         arc2_frame, arc2_sin, arc2_cos, arc2_r,
                         has_straight, sa_x, sa_y, sa_z, sb_x, sb_y, sb_z,
                         start_x, start_y, start_z, end_x, end_y, end_z)
            
            if d < epsilon:
                hit = True
                break
            
            t += max(d, 0.001)
            if t > max_dist:
                break
        
        hit_distances[i] = t if hit else -1.0
    
    return hit_distances


@njit(parallel=True, fastmath=True)
def compute_normals_parallel(hit_points: np.ndarray, hit_mask: np.ndarray,
                             tube_r: float,
                             has_arc1: bool, arc1_cx: float, arc1_cy: float, arc1_cz: float,
                             arc1_frame: np.ndarray, arc1_sin: float, arc1_cos: float, arc1_r: float,
                             has_arc2: bool, arc2_cx: float, arc2_cy: float, arc2_cz: float,
                             arc2_frame: np.ndarray, arc2_sin: float, arc2_cos: float, arc2_r: float,
                             has_straight: bool,
                             sa_x: float, sa_y: float, sa_z: float,
                             sb_x: float, sb_y: float, sb_z: float,
                             start_x: float, start_y: float, start_z: float,
                             end_x: float, end_y: float, end_z: float,
                             delta: float) -> np.ndarray:
    """Compute normals via central differences."""
    num_points = hit_points.shape[0]
    normals = np.zeros((num_points, 3), dtype=np.float64)
    
    for i in prange(num_points):
        if not hit_mask[i]:
            continue
        
        px, py, pz = hit_points[i, 0], hit_points[i, 1], hit_points[i, 2]
        
        def sdf_at(x, y, z):
            return sdf_link(x, y, z, tube_r,
                           has_arc1, arc1_cx, arc1_cy, arc1_cz,
                           arc1_frame, arc1_sin, arc1_cos, arc1_r,
                           has_arc2, arc2_cx, arc2_cy, arc2_cz,
                           arc2_frame, arc2_sin, arc2_cos, arc2_r,
                           has_straight, sa_x, sa_y, sa_z, sb_x, sb_y, sb_z,
                           start_x, start_y, start_z, end_x, end_y, end_z)
        
        nx = sdf_at(px + delta, py, pz) - sdf_at(px - delta, py, pz)
        ny = sdf_at(px, py + delta, pz) - sdf_at(px, py - delta, pz)
        nz = sdf_at(px, py, pz + delta) - sdf_at(px, py, pz - delta)
        
        length = np.sqrt(nx*nx + ny*ny + nz*nz)
        if length > 1e-10:
            normals[i, 0] = nx / length
            normals[i, 1] = ny / length
            normals[i, 2] = nz / length
    
    return normals


# =============================================================================
# Data Extraction (using YOUR Arc3D logic)
# =============================================================================

def extract_link_data(link):
    """
    Extract link geometry directly from LinkCSC/Arc3D objects.
    Uses the exact same attributes the classes compute.
    """
    data = {
        'tube_r': link.r,
        'has_arc1': link.arc1 is not None and link.path.theta1 > link.EPSILON,
        'has_arc2': link.arc2 is not None and link.path.theta2 > link.EPSILON,
        'has_straight': link.path.tMag > link.DISTANCE_EPSILON,
        'start_pos': link.StartDubinsPose.t.copy(),
        'end_pos': link.EndDubinsPose.t.copy(),
    }
    
    # Arc 1 - pull data directly from the Arc3D object
    if data['has_arc1']:
        arc = link.arc1
        data['arc1_center'] = arc.circleCenter.copy()
        data['arc1_r'] = arc.r
        
        # Force computation of local frame if not already done
        if not hasattr(arc, '_worldToLocalRotation'):
            arc._computeLocalFrame()
        
        # Use the exact precomputed values from the Arc3D object
        data['arc1_frame'] = arc._worldToLocalRotation.copy()
        data['arc1_sin'] = arc._sdfSinCos[0]  # sin(theta/2)
        data['arc1_cos'] = arc._sdfSinCos[1]  # cos(theta/2)
    else:
        data['arc1_center'] = np.zeros(3)
        data['arc1_frame'] = np.eye(3)
        data['arc1_sin'] = 0.0
        data['arc1_cos'] = 1.0
        data['arc1_r'] = 1.0
    
    # Arc 2 - same approach
    if data['has_arc2']:
        arc = link.arc2
        data['arc2_center'] = arc.circleCenter.copy()
        data['arc2_r'] = arc.r
        
        if not hasattr(arc, '_worldToLocalRotation'):
            arc._computeLocalFrame()
        
        data['arc2_frame'] = arc._worldToLocalRotation.copy()
        data['arc2_sin'] = arc._sdfSinCos[0]
        data['arc2_cos'] = arc._sdfSinCos[1]
    else:
        data['arc2_center'] = np.zeros(3)
        data['arc2_frame'] = np.eye(3)
        data['arc2_sin'] = 0.0
        data['arc2_cos'] = 1.0
        data['arc2_r'] = 1.0
    
    # Straight segment - use PathCSC values directly
    if data['has_straight']:
        data['straight_a'] = link.path.turn1end.copy()
        data['straight_b'] = (link.path.turn1end + link.path.tMag * link.path.tUnit).copy()
    else:
        data['straight_a'] = np.zeros(3)
        data['straight_b'] = np.zeros(3)
    
    return data


# =============================================================================
# Camera and Rendering
# =============================================================================

def create_camera(center, distance, azimuth, elevation, resolution=(400, 300)):
    """Create camera rays for ray marching."""
    width, height = resolution
    
    cam_x = distance * np.cos(elevation) * np.cos(azimuth)
    cam_y = distance * np.cos(elevation) * np.sin(azimuth)
    cam_z = distance * np.sin(elevation)
    cam_pos = center + np.array([cam_x, cam_y, cam_z])
    
    forward = center - cam_pos
    forward = forward / np.linalg.norm(forward)
    
    world_up = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(forward, world_up)) > 0.99:
        world_up = np.array([1.0, 0.0, 0.0])
    
    right = np.cross(forward, world_up)
    right = right / np.linalg.norm(right)
    up = np.cross(right, forward)
    
    fov_rad = np.deg2rad(50)
    tan_fov = np.tan(fov_rad / 2)
    aspect = width / height
    
    px = np.linspace(-1, 1, width)
    py = np.linspace(-1, 1, height)
    px, py = np.meshgrid(px, py)
    px = px.flatten() * tan_fov * aspect
    py = py.flatten() * tan_fov
    
    directions = (forward[np.newaxis, :] + 
                 px[:, np.newaxis] * right[np.newaxis, :] + 
                 py[:, np.newaxis] * up[np.newaxis, :])
    directions = directions / np.linalg.norm(directions, axis=1, keepdims=True)
    
    origins = np.tile(cam_pos, (len(directions), 1))
    
    return origins, directions, (width, height)


def render_view(link_data, center, extent, azimuth, elevation, resolution=(400, 300)):
    """Render a single view using ray marching."""
    distance = extent * 2.5
    origins, directions, (width, height) = create_camera(
        center, distance, azimuth, elevation, resolution
    )
    
    d = link_data
    
    hit_distances = ray_march_parallel(
        origins, directions,
        d['tube_r'],
        d['has_arc1'], d['arc1_center'][0], d['arc1_center'][1], d['arc1_center'][2],
        d['arc1_frame'], d['arc1_sin'], d['arc1_cos'], d['arc1_r'],
        d['has_arc2'], d['arc2_center'][0], d['arc2_center'][1], d['arc2_center'][2],
        d['arc2_frame'], d['arc2_sin'], d['arc2_cos'], d['arc2_r'],
        d['has_straight'],
        d['straight_a'][0], d['straight_a'][1], d['straight_a'][2],
        d['straight_b'][0], d['straight_b'][1], d['straight_b'][2],
        d['start_pos'][0], d['start_pos'][1], d['start_pos'][2],
        d['end_pos'][0], d['end_pos'][1], d['end_pos'][2],
        max_steps=150, max_dist=distance*2, epsilon=0.002
    )
    
    hit_mask = hit_distances > 0
    hit_points = origins + hit_distances[:, np.newaxis] * directions
    
    normals = compute_normals_parallel(
        hit_points, hit_mask,
        d['tube_r'],
        d['has_arc1'], d['arc1_center'][0], d['arc1_center'][1], d['arc1_center'][2],
        d['arc1_frame'], d['arc1_sin'], d['arc1_cos'], d['arc1_r'],
        d['has_arc2'], d['arc2_center'][0], d['arc2_center'][1], d['arc2_center'][2],
        d['arc2_frame'], d['arc2_sin'], d['arc2_cos'], d['arc2_r'],
        d['has_straight'],
        d['straight_a'][0], d['straight_a'][1], d['straight_a'][2],
        d['straight_b'][0], d['straight_b'][1], d['straight_b'][2],
        d['start_pos'][0], d['start_pos'][1], d['start_pos'][2],
        d['end_pos'][0], d['end_pos'][1], d['end_pos'][2],
        delta=0.001
    )
    
    # Light from camera
    cam_dir = np.array([
        np.cos(elevation) * np.cos(azimuth),
        np.cos(elevation) * np.sin(azimuth),
        np.sin(elevation)
    ])
    
    brightness = np.maximum(0, np.sum(normals * cam_dir, axis=1))
    brightness = 0.3 + 0.7 * brightness
    
    image = np.zeros((height * width, 3))
    image[hit_mask] = brightness[hit_mask, np.newaxis] * np.array([0.8, 0.6, 0.4])
    image[~hit_mask] = [0.1, 0.1, 0.15]
    
    image = image.reshape(height, width, 3)[::-1, :, :]
    
    return image


def get_link_bounds(link, padding=3.0):
    """Get bounding box of the link."""
    points = [link.StartDubinsPose.t, link.EndDubinsPose.t, link.path.turn1end]
    if link.path.tMag > 0:
        points.append(link.path.turn1end + link.path.tMag * link.path.tUnit)
    if link.arc1:
        points.append(link.path.circleCenter1)
    if link.arc2:
        points.append(link.path.circleCenter2)
    
    points = np.array(points)
    center = points.mean(axis=0)
    extent = np.max(np.linalg.norm(points - center, axis=1)) + padding
    
    return center, extent


# =============================================================================
# Verification: Compare Numba SDF with CPU class SDF
# =============================================================================

def verify_sdf(link, link_data, num_tests=100):
    """Verify Numba SDF matches CPU class SDF."""
    print("\n=== SDF Verification ===")
    
    center, extent = get_link_bounds(link)
    
    # Random test points
    np.random.seed(123)
    test_points = center + (np.random.rand(num_tests, 3) - 0.5) * extent * 2
    
    max_diff = 0.0
    for pt in test_points:
        # CPU class method
        cpu_val = link.sdf(pt)
        
        # Numba function
        d = link_data
        numba_val = sdf_link(
            pt[0], pt[1], pt[2],
            d['tube_r'],
            d['has_arc1'], d['arc1_center'][0], d['arc1_center'][1], d['arc1_center'][2],
            d['arc1_frame'], d['arc1_sin'], d['arc1_cos'], d['arc1_r'],
            d['has_arc2'], d['arc2_center'][0], d['arc2_center'][1], d['arc2_center'][2],
            d['arc2_frame'], d['arc2_sin'], d['arc2_cos'], d['arc2_r'],
            d['has_straight'],
            d['straight_a'][0], d['straight_a'][1], d['straight_a'][2],
            d['straight_b'][0], d['straight_b'][1], d['straight_b'][2],
            d['start_pos'][0], d['start_pos'][1], d['start_pos'][2],
            d['end_pos'][0], d['end_pos'][1], d['end_pos'][2]
        )
        
        diff = abs(cpu_val - numba_val)
        max_diff = max(max_diff, diff)
    
    print(f"  Tested {num_tests} random points")
    print(f"  Max difference: {max_diff:.2e}")
    print(f"  Status: {'✓ PASS' if max_diff < 1e-6 else '✗ FAIL'}")
    
    return max_diff < 1e-6


# =============================================================================
# Main
# =============================================================================

def random_unit_vector():
    v = np.random.randn(3)
    return v / np.linalg.norm(v)


def create_random_link(r=None, maxAngle=None, seed=None):
    """
    Create a random LinkCSC using the standard shortestCSC path finder.
    
    This relies on LinkCSC's internal PathCSC computation which uses
    scipy.optimize.fsolve to find a valid CSC Dubins path.
    
    Parameters:
    -----------
    r : float, optional
        Turn radius. If None, randomly chosen from [0.5, 2.0]
    maxAngle : float, optional
        Max angle per elbow (radians). If None, randomly chosen from [π/6, π/2]
    seed : int, optional
        Random seed for reproducibility
    """
    if seed is not None:
        np.random.seed(seed)
    
    # Randomize turn radius if not specified
    if r is None:
        r = np.random.uniform(0.5, 2.0)
    
    # Randomize max elbow angle if not specified
    if maxAngle is None:
        maxAngle = np.random.uniform(np.pi/6, np.pi/2)
    
    for _ in range(20):
        try:
            # Random start pose
            pos1 = np.random.uniform(-3, 3, 3)
            axis1 = random_unit_vector()
            R1 = SO3.AngleAxis(np.random.uniform(0, 2*np.pi), axis1)
            start = SE3.Rt(R1, pos1)
            
            # Random end pose - separation should scale with turn radius
            min_sep = 4 * r
            max_sep = 8 * r
            end_pos = start.t + random_unit_vector() * np.random.uniform(min_sep, max_sep)
            end_dir = random_unit_vector()
            
            # Build orthonormal frame for end pose
            v1 = end_dir
            v2 = random_unit_vector()
            v2 = v2 - np.dot(v2, v1) * v1
            v2 = v2 / np.linalg.norm(v2)
            v3 = np.cross(v1, v2)
            
            end = SE3.Rt(SO3(np.column_stack([v1, v2, v3])), end_pos)
            
            # Let LinkCSC use shortestCSC to find a valid path
            link = LinkCSC(r, start, end, maxAnglePerElbow=maxAngle)
            return link
        except:
            continue
    raise RuntimeError("Failed to create link after 20 attempts")


def main():
    print("=" * 60)
    print("CPU Ray Marching Test - Using YOUR Arc3D Logic")
    print("=" * 60)
    
    # Create a random link with variable parameters
    print("\nCreating random Dubins path...")
    # Pass None to randomize, or specify values to test specific cases
    link = create_random_link(r=.5, maxAngle=None, seed=np.random.seed(494))
    
    print(f"  Turn radius: {link.r:.2f}")
    print(f"  Max elbow angle: {np.rad2deg(link.maxAnglePerElbow):.1f}°")
    print(f"  Arc1 theta: {np.rad2deg(link.path.theta1):.1f}°" if link.arc1 else "  Arc1: None")
    print(f"  Arc2 theta: {np.rad2deg(link.path.theta2):.1f}°" if link.arc2 else "  Arc2: None")
    print(f"  Straight length: {link.path.tMag:.2f}")
    
    # Extract data using YOUR logic
    print("\nExtracting geometry (using YOUR Arc3D._computeLocalFrame logic)...")
    link_data = extract_link_data(link)
    
    # Verify SDF matches
    verify_sdf(link, link_data)
    
    # Get bounds
    center, extent = get_link_bounds(link)
    
    # Views
    views = [
        (0, 0, "Front"),
        (np.pi/2, 0, "Side"),
        (0, np.pi/2 - 0.01, "Top"),
        (np.pi/4, np.pi/6, "Iso"),
    ]
    
    # Create figure
    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(2, 3, figure=fig, width_ratios=[1, 1, 1.2])
    
    # Render ray-marched views
    print("\nRendering ray-marched views...")
    for idx, (azimuth, elevation, title) in enumerate(views):
        row, col = idx // 2, idx % 2
        
        t0 = time.time()
        image = render_view(link_data, center, extent, azimuth, elevation, resolution=(400, 300))
        t1 = time.time()
        print(f"  {title}: {t1-t0:.2f}s")
        
        ax = fig.add_subplot(gs[row, col])
        ax.imshow(image)
        ax.set_title(f"SDF Ray March - {title}")
        ax.axis('off')
    
    # Ground truth
    print("Creating matplotlib ground truth...")
    ax3d = fig.add_subplot(gs[:, 2], projection='3d')
    link.addToPlot(ax3d, numSides=32, color='coral', alpha=0.7, showPath=True, pathColor='blue')
    ax3d.set_title("Ground Truth (Matplotlib 3D)")
    ax3d.set_aspect('equal')
    ax3d.view_init(elev=30, azim=45)
    ax3d.axis('off')
    
    plt.tight_layout()
    plt.savefig('sdf_test_render.png', dpi=150, bbox_inches='tight')
    print("\nSaved to sdf_test_render.png")
    plt.show()


if __name__ == "__main__":
    main()
