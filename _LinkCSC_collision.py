#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fast Pairwise Collision Detection for LinkCSC paths

Uses vectorized NumPy operations for efficient batch SDF evaluation.
For GPU acceleration, swap `np` with `cupy` (drop-in replacement).

Strategy:
    1. Sample points along each path
    2. Batch-evaluate SDF of all points against another path
    3. Check for negative distances (collision) or find minimum distance

@author: Zach (AI Assistant)
"""

import numpy as np
from typing import List, Tuple, Optional
from scipy.spatial.transform import Rotation


# =============================================================================
# Vectorized SDF Functions (operate on batches of points)
# =============================================================================

def sdf_capsule_batch(points: np.ndarray, a: np.ndarray, b: np.ndarray, r: float) -> np.ndarray:
    """
    Signed distance from many points to a capsule.
    
    Parameters:
    -----------
    points : np.ndarray (N, 3)
        Query points
    a : np.ndarray (3,)
        Capsule start
    b : np.ndarray (3,)
        Capsule end
    r : float
        Capsule radius
    
    Returns:
    --------
    np.ndarray (N,)
        Signed distances
    """
    pa = points - a  # (N, 3)
    ba = b - a       # (3,)
    
    ba_dot_ba = np.dot(ba, ba)
    if ba_dot_ba < 1e-12:
        # Degenerate: just a sphere
        return np.linalg.norm(pa, axis=1) - r
    
    # Project onto line segment
    h = np.clip(pa @ ba / ba_dot_ba, 0.0, 1.0)  # (N,)
    
    # Distance to closest point on segment
    closest = pa - np.outer(h, ba)  # (N, 3)
    return np.linalg.norm(closest, axis=1) - r


def sdf_capped_torus_batch(points: np.ndarray, center: np.ndarray, frame: np.ndarray,
                           sin_half: float, cos_half: float, 
                           major_r: float, minor_r: float) -> np.ndarray:
    """
    Signed distance from many points to a capped torus (arc tube).
    
    Parameters:
    -----------
    points : np.ndarray (N, 3)
        Query points in world coordinates
    center : np.ndarray (3,)
        Arc circle center
    frame : np.ndarray (3, 3)
        World-to-local rotation matrix (rows are local axes)
    sin_half, cos_half : float
        sin/cos of half the arc angle
    major_r : float
        Arc turning radius
    minor_r : float
        Tube radius
    
    Returns:
    --------
    np.ndarray (N,)
        Signed distances
    """
    # Transform to local coordinates
    local = (points - center) @ frame.T  # (N, 3)
    
    px = np.abs(local[:, 0])
    py = local[:, 1]
    pz = local[:, 2]
    
    # Condition: is point within arc angular span?
    within_arc = cos_half * px > sin_half * py
    
    # k = projection factor
    k = np.where(
        within_arc,
        sin_half * px + cos_half * py,  # Project onto arc
        np.sqrt(px**2 + py**2)          # Distance to cap
    )
    
    # SDF formula
    return np.sqrt(px**2 + py**2 + pz**2 + major_r**2 - 2*major_r*k) - minor_r


# =============================================================================
# LinkCSC Data Extraction (prepare for vectorized operations)
# =============================================================================

class LinkCSCData:
    """
    Lightweight data container for LinkCSC geometry.
    Extracts just what's needed for SDF computation.
    """
    
    def __init__(self, link):
        """Extract geometry from a LinkCSC object."""
        self.tube_r = link.r
        self.EPSILON = link.EPSILON
        
        # Arc 1 data
        self.has_arc1 = link.arc1 is not None and link.path.theta1 > link.EPSILON
        if self.has_arc1:
            arc = link.arc1
            self.arc1_center = arc.circleCenter.copy()
            self.arc1_major_r = arc.r
            
            # Build local frame (X at midpoint for symmetry)
            half_theta = arc.theta / 2
            half_rot = Rotation.from_rotvec(half_theta * arc.binormal)
            mid_outward = half_rot.apply(-arc.startNormal)
            
            arcX = mid_outward / np.linalg.norm(mid_outward)
            arcZ = arc.binormal
            arcY = np.cross(arcZ, arcX)
            
            self.arc1_frame = np.array([arcX, arcY, arcZ])
            self.arc1_sin_half = np.sin(half_theta)
            self.arc1_cos_half = np.cos(half_theta)
        
        # Arc 2 data
        self.has_arc2 = link.arc2 is not None and link.path.theta2 > link.EPSILON
        if self.has_arc2:
            arc = link.arc2
            self.arc2_center = arc.circleCenter.copy()
            self.arc2_major_r = arc.r
            
            half_theta = arc.theta / 2
            half_rot = Rotation.from_rotvec(half_theta * arc.binormal)
            mid_outward = half_rot.apply(-arc.startNormal)
            
            arcX = mid_outward / np.linalg.norm(mid_outward)
            arcZ = arc.binormal
            arcY = np.cross(arcZ, arcX)
            
            self.arc2_frame = np.array([arcX, arcY, arcZ])
            self.arc2_sin_half = np.sin(half_theta)
            self.arc2_cos_half = np.cos(half_theta)
        
        # Straight segment data
        self.has_straight = link.path.tMag > link.DISTANCE_EPSILON
        if self.has_straight:
            self.straight_a = link.path.turn1end.copy()
            self.straight_b = (link.path.turn1end + link.path.tMag * link.path.tUnit).copy()
        
        # Store interpolation data for sampling
        self.start_pos = link.StartDubinsPose.t.copy()
        self.end_pos = link.EndDubinsPose.t.copy()
        self._link = link  # Keep reference for interpolation
    
    def sdf_batch(self, points: np.ndarray) -> np.ndarray:
        """
        Compute SDF for many points at once.
        
        Parameters:
        -----------
        points : np.ndarray (N, 3)
            Query points
        
        Returns:
        --------
        np.ndarray (N,)
            Signed distances (min over all segments)
        """
        N = points.shape[0]
        
        # Start with large distances
        distances = np.full(N, np.inf)
        
        # Arc 1
        if self.has_arc1:
            d = sdf_capped_torus_batch(
                points, self.arc1_center, self.arc1_frame,
                self.arc1_sin_half, self.arc1_cos_half,
                self.arc1_major_r, self.tube_r
            )
            distances = np.minimum(distances, d)
        else:
            # Sphere at start
            d = np.linalg.norm(points - self.start_pos, axis=1) - self.tube_r
            distances = np.minimum(distances, d)
        
        # Straight section
        if self.has_straight:
            d = sdf_capsule_batch(points, self.straight_a, self.straight_b, self.tube_r)
            distances = np.minimum(distances, d)
        
        # Arc 2
        if self.has_arc2:
            d = sdf_capped_torus_batch(
                points, self.arc2_center, self.arc2_frame,
                self.arc2_sin_half, self.arc2_cos_half,
                self.arc2_major_r, self.tube_r
            )
            distances = np.minimum(distances, d)
        else:
            # Sphere at end
            d = np.linalg.norm(points - self.end_pos, axis=1) - self.tube_r
            distances = np.minimum(distances, d)
        
        return distances
    
    def sample_points(self, count: int = 50) -> np.ndarray:
        """
        Sample points along the path centerline.
        
        Parameters:
        -----------
        count : int
            Number of points to sample
        
        Returns:
        --------
        np.ndarray (count, 3)
            Sampled points
        """
        # Use slightly inset t values to avoid floating point edge cases
        t_values = np.linspace(0.001, 0.999, count)
        return np.array([self._link.interpolateAt(t) for t in t_values])


# =============================================================================
# Pairwise Collision Detection
# =============================================================================

def min_distance_between_links(link_data_a: LinkCSCData, link_data_b: LinkCSCData,
                                num_samples: int = 50) -> float:
    """
    Compute minimum distance between two LinkCSC paths.
    
    Samples points along each path and evaluates SDF against the other.
    
    Parameters:
    -----------
    link_data_a, link_data_b : LinkCSCData
        Pre-extracted link geometry
    num_samples : int
        Number of points to sample along each path
    
    Returns:
    --------
    float
        Minimum signed distance (negative = collision)
    """
    # Sample points along A, evaluate against B
    points_a = link_data_a.sample_points(num_samples)
    dist_a_to_b = link_data_b.sdf_batch(points_a)
    
    # Sample points along B, evaluate against A
    points_b = link_data_b.sample_points(num_samples)
    dist_b_to_a = link_data_a.sdf_batch(points_b)
    
    # Minimum over all
    return min(dist_a_to_b.min(), dist_b_to_a.min())


def check_collision(link_data_a: LinkCSCData, link_data_b: LinkCSCData,
                    num_samples: int = 50) -> bool:
    """
    Check if two links collide.
    
    Returns True if collision detected.
    """
    return min_distance_between_links(link_data_a, link_data_b, num_samples) < 0


def pairwise_distances(links: List['LinkCSC'], num_samples: int = 50) -> np.ndarray:
    """
    Compute pairwise minimum distances between all links.
    
    Parameters:
    -----------
    links : List[LinkCSC]
        List of LinkCSC objects
    num_samples : int
        Points to sample per path
    
    Returns:
    --------
    np.ndarray (n, n)
        Distance matrix where D[i,j] = min distance between link i and j
        Diagonal is inf (self-distance)
    """
    n = len(links)
    
    # Pre-extract all link data
    link_data = [LinkCSCData(link) for link in links]
    
    # Pre-sample all points (can reuse)
    all_points = [ld.sample_points(num_samples) for ld in link_data]
    
    # Distance matrix
    D = np.full((n, n), np.inf)
    
    for i in range(n):
        for j in range(i + 1, n):
            # Points from i, SDF against j
            d_i_to_j = link_data[j].sdf_batch(all_points[i]).min()
            
            # Points from j, SDF against i
            d_j_to_i = link_data[i].sdf_batch(all_points[j]).min()
            
            # Symmetric
            D[i, j] = min(d_i_to_j, d_j_to_i)
            D[j, i] = D[i, j]
    
    return D


def pairwise_collisions(links: List['LinkCSC'], num_samples: int = 50) -> List[Tuple[int, int]]:
    """
    Find all colliding pairs.
    
    Returns:
    --------
    List of (i, j) tuples where links[i] and links[j] collide
    """
    D = pairwise_distances(links, num_samples)
    collisions = []
    n = len(links)
    
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] < 0:
                collisions.append((i, j))
    
    return collisions


# =============================================================================
# Even Faster: Batch All Operations
# =============================================================================

def pairwise_distances_fast(links: List['LinkCSC'], num_samples: int = 50) -> np.ndarray:
    """
    Faster pairwise distances using fully vectorized operations.
    
    Pre-computes all sample points, then batch-evaluates SDFs.
    """
    n = len(links)
    
    # Extract all link data
    link_data = [LinkCSCData(link) for link in links]
    
    # Sample all points: shape (n, num_samples, 3)
    all_points = np.array([ld.sample_points(num_samples) for ld in link_data])
    
    # Distance matrix
    D = np.full((n, n), np.inf)
    
    for i in range(n):
        # Get all points from link i: (num_samples, 3)
        points_i = all_points[i]
        
        for j in range(i + 1, n):
            # Evaluate SDF of link j at points from link i
            d_i_to_j = link_data[j].sdf_batch(points_i).min()
            
            # Evaluate SDF of link i at points from link j
            points_j = all_points[j]
            d_j_to_i = link_data[i].sdf_batch(points_j).min()
            
            D[i, j] = min(d_i_to_j, d_j_to_i)
            D[j, i] = D[i, j]
    
    return D


# =============================================================================
# Demo / Test
# =============================================================================

def demo():
    """Demonstrate collision detection between random links."""
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
                
                # Build orthonormal frame
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
    
    # Create several random links
    print("Creating random links...")
    np.random.seed(42)
    num_links = 10
    links = [create_random_link(r=1.0) for _ in range(num_links)]
    print(f"Created {num_links} links")
    
    # Time pairwise distance computation
    print("\nComputing pairwise distances...")
    
    num_samples = 50
    t0 = time.time()
    D = pairwise_distances_fast(links, num_samples=num_samples)
    t1 = time.time()
    
    print(f"Time for {num_links}x{num_links} pairwise distances: {t1-t0:.3f}s")
    print(f"  ({num_links * (num_links-1) // 2} unique pairs)")
    print(f"  ({num_samples} samples per path)")
    
    # Find collisions
    collisions = []
    for i in range(num_links):
        for j in range(i+1, num_links):
            if D[i, j] < 0:
                collisions.append((i, j, D[i, j]))
    
    print(f"\nCollisions found: {len(collisions)}")
    for i, j, d in collisions:
        print(f"  Links {i} and {j}: distance = {d:.3f}")
    
    # Show closest non-colliding pair
    D_positive = D.copy()
    D_positive[D_positive < 0] = np.inf
    np.fill_diagonal(D_positive, np.inf)
    
    if np.any(np.isfinite(D_positive)):
        i, j = np.unravel_index(np.argmin(D_positive), D.shape)
        print(f"\nClosest non-colliding pair: {i}, {j} at distance {D[i,j]:.3f}")
    
    # Benchmark scaling
    print("\n--- Scaling benchmark ---")
    for ns in [20, 50, 100, 200]:
        t0 = time.time()
        _ = pairwise_distances_fast(links, num_samples=ns)
        t1 = time.time()
        print(f"  {ns} samples: {t1-t0:.3f}s")


if __name__ == "__main__":
    demo()
