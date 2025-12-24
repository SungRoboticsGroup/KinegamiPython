"""cpu_serial_link_sdf.py

CPU serial version using LinkCSC methods directly with AABB pruning.
"""

import numpy as np
from typing import List, Tuple
import signal
from contextlib import contextmanager


class TimeoutException(Exception):
    pass


@contextmanager
def time_limit(seconds):
    """Context manager for timeout"""
    def signal_handler(signum, frame):
        raise TimeoutException("Timed out!")
    
    signal.signal(signal.SIGALRM, signal_handler)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)


class AABB:
    """Axis-Aligned Bounding Box"""
    def __init__(self, min_corner: np.ndarray, max_corner: np.ndarray):
        self.min = min_corner
        self.max = max_corner
        self.center = (self.min + self.max) / 2
        self.half_size = (self.max - self.min) / 2
    
    def distance_to_aabb(self, other: 'AABB') -> float:
        """Compute minimum distance between two AABBs"""
        # Distance in each dimension
        dx = max(0, max(self.min[0] - other.max[0], other.min[0] - self.max[0]))
        dy = max(0, max(self.min[1] - other.max[1], other.min[1] - self.max[1]))
        dz = max(0, max(self.min[2] - other.max[2], other.min[2] - self.max[2]))
        return np.sqrt(dx*dx + dy*dy + dz*dz)
    
    @staticmethod
    def from_points(points: np.ndarray) -> 'AABB':
        """Create AABB from point cloud"""
        return AABB(np.min(points, axis=0), np.max(points, axis=0))
    
    @staticmethod
    def from_link(link) -> 'AABB':
        """Create AABB from LinkCSC by sampling points"""
        # Sample points densely for accurate bounding box
        points = link.interpolate(density=20.0)
        
        # Expand by radius to account for tube thickness
        min_corner = np.min(points, axis=0) - link.r
        max_corner = np.max(points, axis=0) + link.r
        
        return AABB(min_corner, max_corner)


def pairwise_link_distances_serial(
    links: List,
    points: np.ndarray,
    point_link_ids: np.ndarray,
    timeout_seconds: int = 300,
    use_aabb_pruning: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Pairwise minimum distances between all link pairs (CPU serial).
    Uses LinkCSC.sdf() method directly with AABB pruning.
    
    Parameters:
    -----------
    links : List[LinkCSC]
        List of LinkCSC objects
    points : np.ndarray
        Sampled points (N, 3)
    point_link_ids : np.ndarray
        Link ID for each point (N,)
    timeout_seconds : int
        Timeout in seconds
    use_aabb_pruning : bool
        Whether to use AABB pruning to skip distant pairs
    
    Returns:
    --------
    pairwise_dist : np.ndarray
        (L, L) pairwise distances
    point_idx : np.ndarray
        (L, L) point indices achieving minimum
    """
    
    L = len(links)
    N = len(points)
    
    # Initialize output arrays
    pairwise_dist = np.full((L, L), np.inf, dtype=np.float32)
    point_idx = np.full((L, L), -1, dtype=np.int32)
    
    # Precompute AABBs for all links
    aabbs = None
    aabb_distances = None
    if use_aabb_pruning:
        print("  Computing AABBs for pruning...")
        aabbs = [AABB.from_link(link) for link in links]
        
        # Precompute AABB distances
        aabb_distances = np.zeros((L, L), dtype=np.float32)
        for i in range(L):
            for j in range(i+1, L):
                dist = aabbs[i].distance_to_aabb(aabbs[j])
                aabb_distances[i, j] = dist
                aabb_distances[j, i] = dist
    
    try:
        with time_limit(timeout_seconds):
            # Group points by link
            points_by_link = [[] for _ in range(L)]
            point_indices_by_link = [[] for _ in range(L)]
            
            for point_i in range(N):
                link_i = point_link_ids[point_i]
                points_by_link[link_i].append(points[point_i])
                point_indices_by_link[link_i].append(point_i)
            
            # Convert to arrays
            for link_i in range(L):
                if points_by_link[link_i]:
                    points_by_link[link_i] = np.array(points_by_link[link_i])
                    point_indices_by_link[link_i] = np.array(point_indices_by_link[link_i])
                else:
                    points_by_link[link_i] = np.zeros((0, 3))
                    point_indices_by_link[link_i] = np.array([], dtype=np.int32)
            
            # Compute pairwise distances
            total_pairs = L * (L - 1) // 2
            processed_pairs = 0
            skipped_pairs = 0
            
            for link_i in range(L):
                link_i_obj = links[link_i]
                link_i_points = points_by_link[link_i]
                link_i_indices = point_indices_by_link[link_i]
                
                if len(link_i_points) == 0:
                    continue
                
                for link_j in range(link_i + 1, L):
                    # AABB pruning: skip if AABBs are too far apart
                    if use_aabb_pruning:
                        # Conservative threshold: if AABBs are separated by more than 
                        # sum of radii, the links can't be closer
                        aabb_dist = aabb_distances[link_i, link_j]
                        max_possible_proximity = link_i_obj.r + links[link_j].r
                        
                        if aabb_dist > max_possible_proximity:
                            skipped_pairs += 1
                            continue
                    
                    link_j_obj = links[link_j]
                    
                    # Compute distance from link_i points to link_j
                    for idx, point in zip(link_i_indices, link_i_points):
                        dist = link_j_obj.sdf(point)
                        
                        if dist < pairwise_dist[link_i, link_j]:
                            pairwise_dist[link_i, link_j] = dist
                            point_idx[link_i, link_j] = idx
                    
                    # Compute distance from link_j points to link_i
                    link_j_points = points_by_link[link_j]
                    link_j_indices = point_indices_by_link[link_j]
                    
                    for idx, point in zip(link_j_indices, link_j_points):
                        dist = link_i_obj.sdf(point)
                        
                        if dist < pairwise_dist[link_j, link_i]:
                            pairwise_dist[link_j, link_i] = dist
                            point_idx[link_j, link_i] = idx
                    
                    processed_pairs += 1
            
            if use_aabb_pruning:
                print(f"  AABB pruning: skipped {skipped_pairs}/{total_pairs} pairs "
                      f"({100*skipped_pairs/total_pairs:.1f}%)")
        
        return pairwise_dist, point_idx
        
    except TimeoutException:
        print(f"  WARNING: Serial method timed out after {timeout_seconds}s")
        raise
