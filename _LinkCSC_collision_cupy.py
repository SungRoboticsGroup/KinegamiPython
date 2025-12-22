#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPU-Accelerated Collision Detection for LinkCSC paths using CuPy

CuPy is a drop-in replacement for NumPy that runs on GPU.
All array operations are automatically parallelized on the GPU.

Usage:
    from LinkCSC_collision_cupy import compute_arrangement_error, pairwise_distances_gpu
    
    error = compute_arrangement_error(links, num_samples=50, margin=0.5)
    D = pairwise_distances_gpu(links, num_samples=50)

@author: Zach (AI Assistant)
"""

import numpy as np
from typing import List, Tuple, Optional
from scipy.spatial.transform import Rotation

# Try to import CuPy, fall back to NumPy if not available
try:
    import cupy as cp
    GPU_AVAILABLE = True
    print("CuPy GPU acceleration available")
except ImportError:
    import numpy as cp  # Fall back to NumPy
    GPU_AVAILABLE = False
    print("CuPy not available, using NumPy (CPU) instead")
    print("Install with: pip install cupy-cuda11x  (or cupy-cuda12x)")


# =============================================================================
# GPU-Accelerated SDF Functions
# =============================================================================

def sdf_capsule_batch_gpu(points, a, b, r):
    """
    Signed distance from many points to a capsule (GPU).
    
    Parameters:
    -----------
    points : cp.ndarray (N, 3)
        Query points (on GPU)
    a : cp.ndarray (3,)
        Capsule start (on GPU)
    b : cp.ndarray (3,)
        Capsule end (on GPU)
    r : float
        Capsule radius
    
    Returns:
    --------
    cp.ndarray (N,)
        Signed distances (on GPU)
    """
    pa = points - a  # (N, 3)
    ba = b - a       # (3,)
    
    ba_dot_ba = cp.dot(ba, ba)
    if float(ba_dot_ba) < 1e-12:
        # Degenerate: just a sphere
        return cp.linalg.norm(pa, axis=1) - r
    
    # Project onto line segment, clamp to [0, 1]
    h = cp.clip(pa @ ba / ba_dot_ba, 0.0, 1.0)  # (N,)
    
    # Distance to closest point on segment
    closest = pa - cp.outer(h, ba)  # (N, 3)
    return cp.linalg.norm(closest, axis=1) - r


def sdf_capped_torus_batch_gpu(points, center, frame, sin_half, cos_half, major_r, minor_r):
    """
    Signed distance from many points to a capped torus / arc tube (GPU).
    
    Parameters:
    -----------
    points : cp.ndarray (N, 3)
        Query points in world coordinates (on GPU)
    center : cp.ndarray (3,)
        Arc circle center (on GPU)
    frame : cp.ndarray (3, 3)
        World-to-local rotation matrix (on GPU)
    sin_half, cos_half : float
        sin/cos of half the arc angle
    major_r : float
        Arc turning radius
    minor_r : float
        Tube radius
    
    Returns:
    --------
    cp.ndarray (N,)
        Signed distances (on GPU)
    """
    # Transform to local coordinates
    local = (points - center) @ frame.T  # (N, 3)
    
    px = cp.abs(local[:, 0])
    py = local[:, 1]
    pz = local[:, 2]
    
    # Condition: is point within arc angular span?
    within_arc = cos_half * px > sin_half * py
    
    # k = projection factor
    k = cp.where(
        within_arc,
        sin_half * px + cos_half * py,  # Project onto arc
        cp.sqrt(px**2 + py**2)          # Distance to cap
    )
    
    # SDF formula
    return cp.sqrt(px**2 + py**2 + pz**2 + major_r**2 - 2*major_r*k) - minor_r


def sdf_sphere_batch_gpu(points, center, r):
    """Signed distance from many points to a sphere (GPU)."""
    return cp.linalg.norm(points - center, axis=1) - r


# =============================================================================
# LinkCSC Data Container (prepares data for GPU)
# =============================================================================

class LinkCSCDataGPU:
    """
    GPU-ready data container for LinkCSC geometry.
    Extracts geometry and uploads to GPU memory.
    """
    
    def __init__(self, link):
        """Extract geometry from a LinkCSC object and upload to GPU."""
        self.tube_r = float(link.r)
        self.EPSILON = float(link.EPSILON)
        
        # Arc 1 data
        self.has_arc1 = link.arc1 is not None and link.path.theta1 > link.EPSILON
        if self.has_arc1:
            arc = link.arc1
            
            # Build local frame (X at midpoint for symmetry)
            half_theta = arc.theta / 2
            half_rot = Rotation.from_rotvec(half_theta * arc.binormal)
            mid_outward = half_rot.apply(-arc.startNormal)
            
            arcX = mid_outward / np.linalg.norm(mid_outward)
            arcZ = arc.binormal
            arcY = np.cross(arcZ, arcX)
            
            # Upload to GPU
            self.arc1_center = cp.array(arc.circleCenter, dtype=cp.float32)
            self.arc1_frame = cp.array(np.array([arcX, arcY, arcZ]), dtype=cp.float32)
            self.arc1_sin_half = float(np.sin(half_theta))
            self.arc1_cos_half = float(np.cos(half_theta))
            self.arc1_major_r = float(arc.r)
        
        # Arc 2 data
        self.has_arc2 = link.arc2 is not None and link.path.theta2 > link.EPSILON
        if self.has_arc2:
            arc = link.arc2
            
            half_theta = arc.theta / 2
            half_rot = Rotation.from_rotvec(half_theta * arc.binormal)
            mid_outward = half_rot.apply(-arc.startNormal)
            
            arcX = mid_outward / np.linalg.norm(mid_outward)
            arcZ = arc.binormal
            arcY = np.cross(arcZ, arcX)
            
            # Upload to GPU
            self.arc2_center = cp.array(arc.circleCenter, dtype=cp.float32)
            self.arc2_frame = cp.array(np.array([arcX, arcY, arcZ]), dtype=cp.float32)
            self.arc2_sin_half = float(np.sin(half_theta))
            self.arc2_cos_half = float(np.cos(half_theta))
            self.arc2_major_r = float(arc.r)
        
        # Straight segment data
        self.has_straight = link.path.tMag > link.DISTANCE_EPSILON
        if self.has_straight:
            self.straight_a = cp.array(link.path.turn1end, dtype=cp.float32)
            self.straight_b = cp.array(
                link.path.turn1end + link.path.tMag * link.path.tUnit, 
                dtype=cp.float32
            )
        
        # Endpoint positions (for fallback spheres)
        self.start_pos = cp.array(link.StartDubinsPose.t, dtype=cp.float32)
        self.end_pos = cp.array(link.EndDubinsPose.t, dtype=cp.float32)
        
        # Keep reference for sampling (done on CPU)
        self._link = link
    
    def sdf_batch_gpu(self, points_gpu):
        """
        Compute SDF for many points (all on GPU).
        
        Parameters:
        -----------
        points_gpu : cp.ndarray (N, 3)
            Query points already on GPU
        
        Returns:
        --------
        cp.ndarray (N,)
            Signed distances (on GPU)
        """
        N = points_gpu.shape[0]
        
        # Start with large distances
        distances = cp.full(N, cp.inf, dtype=cp.float32)
        
        # Arc 1
        if self.has_arc1:
            d = sdf_capped_torus_batch_gpu(
                points_gpu, self.arc1_center, self.arc1_frame,
                self.arc1_sin_half, self.arc1_cos_half,
                self.arc1_major_r, self.tube_r
            )
            distances = cp.minimum(distances, d)
        else:
            # Sphere at start
            d = sdf_sphere_batch_gpu(points_gpu, self.start_pos, self.tube_r)
            distances = cp.minimum(distances, d)
        
        # Straight section
        if self.has_straight:
            d = sdf_capsule_batch_gpu(
                points_gpu, self.straight_a, self.straight_b, self.tube_r
            )
            distances = cp.minimum(distances, d)
        
        # Arc 2
        if self.has_arc2:
            d = sdf_capped_torus_batch_gpu(
                points_gpu, self.arc2_center, self.arc2_frame,
                self.arc2_sin_half, self.arc2_cos_half,
                self.arc2_major_r, self.tube_r
            )
            distances = cp.minimum(distances, d)
        else:
            # Sphere at end
            d = sdf_sphere_batch_gpu(points_gpu, self.end_pos, self.tube_r)
            distances = cp.minimum(distances, d)
        
        return distances
    
    def sample_points(self, num_samples: int = 50) -> np.ndarray:
        """Sample points along the path centerline (CPU, returns NumPy array)."""
        # Use slightly inset t values to avoid floating point edge cases
        t_values = np.linspace(0.001, 0.999, num_samples)
        return np.array([self._link.interpolateAt(t) for t in t_values])


# =============================================================================
# Main Functions: Error Computation and Pairwise Distances
# =============================================================================

def compute_arrangement_error(links: List, num_samples: int = 50, 
                               margin: float = 0.0) -> float:
    """
    Compute total collision/proximity error for an arrangement of LinkCSC paths.
    
    This function measures how much the paths overlap or come too close.
    Use this as an objective function for optimization.
    
    Parameters:
    -----------
    links : List[LinkCSC]
        The paths to check for collisions
    num_samples : int
        Number of points to sample along each path
    margin : float
        Minimum desired clearance between paths.
        - margin=0: only penalize actual collisions (penetration)
        - margin>0: also penalize paths that are too close
    
    Returns:
    --------
    float
        Total error:
        - 0 = no collisions and all paths have clearance >= margin
        - >0 = some paths overlap or are too close
    """
    N = len(links)
    if N < 2:
        return 0.0
    
    # Step 1: Extract path data and upload to GPU
    path_data = [LinkCSCDataGPU(link) for link in links]
    
    # Step 2: Sample points from all paths (CPU)
    all_points_list = []
    for pd in path_data:
        points = pd.sample_points(num_samples)
        all_points_list.append(points)
    
    all_points_np = np.vstack(all_points_list)  # (N * num_samples, 3)
    
    # Step 3: Upload all points to GPU
    all_points_gpu = cp.array(all_points_np, dtype=cp.float32)
    
    # Step 4: Track which path each point belongs to
    point_sources = np.repeat(np.arange(N), num_samples)  # [0,0,...,1,1,...,2,2,...]
    
    # Step 5: Compute SDF of all points against each path
    total_error_gpu = cp.float32(0.0)
    
    for j in range(N):
        # SDF of ALL points against path j
        distances = path_data[j].sdf_batch_gpu(all_points_gpu)  # (N * num_samples,)
        
        # Mask out self-comparisons (points from path j against path j)
        self_mask = cp.array(point_sources == j)
        distances = cp.where(self_mask, cp.inf, distances)
        
        # Compute error: penalize distances below margin
        # violation = max(margin - distance, 0)
        violations = cp.maximum(margin - distances, 0)
        
        # Accumulate squared penalty (smooth gradient)
        total_error_gpu += cp.sum(violations ** 2)
    
    # Step 6: Get result back to CPU
    if GPU_AVAILABLE:
        return float(total_error_gpu.get())
    else:
        return float(total_error_gpu)


def compute_arrangement_error_linear(links: List, num_samples: int = 50, 
                                      margin: float = 0.0) -> float:
    """
    Same as compute_arrangement_error but with linear penalty (not squared).
    
    Returns sum of penetration depths and margin violations.
    """
    N = len(links)
    if N < 2:
        return 0.0
    
    path_data = [LinkCSCDataGPU(link) for link in links]
    
    all_points_list = [pd.sample_points(num_samples) for pd in path_data]
    all_points_np = np.vstack(all_points_list)
    all_points_gpu = cp.array(all_points_np, dtype=cp.float32)
    
    point_sources = np.repeat(np.arange(N), num_samples)
    
    total_error_gpu = cp.float32(0.0)
    
    for j in range(N):
        distances = path_data[j].sdf_batch_gpu(all_points_gpu)
        self_mask = cp.array(point_sources == j)
        distances = cp.where(self_mask, cp.inf, distances)
        violations = cp.maximum(margin - distances, 0)
        total_error_gpu += cp.sum(violations)
    
    if GPU_AVAILABLE:
        return float(total_error_gpu.get())
    else:
        return float(total_error_gpu)


def pairwise_distances_gpu(links: List, num_samples: int = 50) -> np.ndarray:
    """
    Compute pairwise minimum distances between all links using GPU.
    
    Parameters:
    -----------
    links : List[LinkCSC]
        List of paths
    num_samples : int
        Points to sample per path
    
    Returns:
    --------
    np.ndarray (N, N)
        Distance matrix where D[i,j] = minimum distance between link i and j.
        Diagonal is inf. Negative values indicate collision.
    """
    N = len(links)
    
    # Extract path data
    path_data = [LinkCSCDataGPU(link) for link in links]
    
    # Sample all points and upload to GPU
    all_points_list = [pd.sample_points(num_samples) for pd in path_data]
    all_points_np = np.vstack(all_points_list)
    all_points_gpu = cp.array(all_points_np, dtype=cp.float32)
    
    # Also keep individual point arrays on GPU
    points_per_path_gpu = [
        cp.array(pts, dtype=cp.float32) for pts in all_points_list
    ]
    
    # Distance matrix
    D = np.full((N, N), np.inf, dtype=np.float32)
    
    for i in range(N):
        for j in range(i + 1, N):
            # Points from i, SDF against j
            d_i_to_j = path_data[j].sdf_batch_gpu(points_per_path_gpu[i])
            min_i_to_j = float(cp.min(d_i_to_j).get() if GPU_AVAILABLE else cp.min(d_i_to_j))
            
            # Points from j, SDF against i
            d_j_to_i = path_data[i].sdf_batch_gpu(points_per_path_gpu[j])
            min_j_to_i = float(cp.min(d_j_to_i).get() if GPU_AVAILABLE else cp.min(d_j_to_i))
            
            # Symmetric
            D[i, j] = min(min_i_to_j, min_j_to_i)
            D[j, i] = D[i, j]
    
    return D


def find_collisions(links: List, num_samples: int = 50) -> List[Tuple[int, int, float]]:
    """
    Find all colliding pairs of links.
    
    Returns:
    --------
    List of (i, j, penetration_depth) tuples where links[i] and links[j] collide.
    """
    D = pairwise_distances_gpu(links, num_samples)
    collisions = []
    N = len(links)
    
    for i in range(N):
        for j in range(i + 1, N):
            if D[i, j] < 0:
                collisions.append((i, j, -D[i, j]))  # penetration as positive
    
    return collisions


# =============================================================================
# Demo / Test
# =============================================================================

def demo():
    """Demonstrate GPU collision detection."""
    import time
    from spatialmath import SE3, SO3
    from LinkCSC import LinkCSC
    
    def random_unit_vector():
        v = np.random.randn(3)
        return v / np.linalg.norm(v)
    
    def random_pose(center=np.zeros(3), spread=5.0):
        pos = center + np.random.uniform(-spread, spread, 3)
        axis = random_unit_vector()
        angle = np.random.uniform(0, 2*np.pi)
        R = SO3.AngleAxis(angle, axis)
        return SE3.Rt(R, pos)
    
    def create_random_link(r=1.0):
        for _ in range(20):
            try:
                start = random_pose()
                end_pos = start.t + random_unit_vector() * np.random.uniform(4, 8)
                end_dir = random_unit_vector()
                
                v1 = end_dir
                v2 = random_unit_vector()
                v2 = v2 - np.dot(v2, v1) * v1
                v2 = v2 / np.linalg.norm(v2)
                v3 = np.cross(v1, v2)
                
                end = SE3.Rt(SO3(np.column_stack([v1, v2, v3])), end_pos)
                return LinkCSC(r, start, end, maxAnglePerElbow=np.pi/3)
            except:
                continue
        raise RuntimeError("Failed to create link")
    
    print("=" * 60)
    print("GPU-Accelerated LinkCSC Collision Detection")
    print(f"GPU Available: {GPU_AVAILABLE}")
    print("=" * 60)
    
    # Create random links
    print("\nCreating random links...")
    np.random.seed(42)
    num_links = 10
    links = [create_random_link(r=1.0) for _ in range(num_links)]
    print(f"Created {num_links} links")
    
    # Test arrangement error
    print("\n--- Arrangement Error ---")
    for margin in [0.0, 0.5, 1.0]:
        t0 = time.time()
        error = compute_arrangement_error(links, num_samples=50, margin=margin)
        t1 = time.time()
        print(f"  margin={margin}: error={error:.4f} ({t1-t0:.3f}s)")
    
    # Test pairwise distances
    print("\n--- Pairwise Distances ---")
    t0 = time.time()
    D = pairwise_distances_gpu(links, num_samples=50)
    t1 = time.time()
    print(f"  Computed {num_links}x{num_links} distance matrix in {t1-t0:.3f}s")
    
    # Find collisions
    collisions = find_collisions(links, num_samples=50)
    print(f"\n  Collisions found: {len(collisions)}")
    for i, j, depth in collisions[:5]:  # Show first 5
        print(f"    Links {i} and {j}: penetration = {depth:.3f}")
    
    # Show closest pair
    D_no_collision = D.copy()
    D_no_collision[D_no_collision < 0] = np.inf
    np.fill_diagonal(D_no_collision, np.inf)
    
    if np.any(np.isfinite(D_no_collision)):
        i, j = np.unravel_index(np.argmin(D_no_collision), D.shape)
        print(f"\n  Closest non-colliding pair: {i}, {j} at distance {D[i,j]:.3f}")
    
    # Benchmark scaling
    print("\n--- Scaling Benchmark ---")
    for ns in [20, 50, 100, 200]:
        t0 = time.time()
        _ = compute_arrangement_error(links, num_samples=ns, margin=0.0)
        t1 = time.time()
        print(f"  {ns} samples: {t1-t0:.3f}s")
    
    print("\n" + "=" * 60)
    print("Done!")


if __name__ == "__main__":
    demo()
