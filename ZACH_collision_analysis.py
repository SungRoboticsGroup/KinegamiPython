"""collision_analysis.py

Tools for analyzing collisions/intersections between links.
"""

import numpy as np
from typing import List, Tuple, Dict
from dataclasses import dataclass


@dataclass
class CollisionPair:
    """Information about a collision between two links"""
    link_i: int
    link_j: int
    distance: float
    point_idx: int
    point_location: np.ndarray
    
    def __repr__(self):
        return f"Collision(links=({self.link_i},{self.link_j}), dist={self.distance:.4f})"


@dataclass
class CollisionStatistics:
    """Statistics about collisions in the link set"""
    total_pairs: int
    collision_pairs: int
    collision_rate: float
    min_distance: float
    max_distance: float
    mean_distance: float
    median_distance: float
    std_distance: float
    worst_collisions: List[CollisionPair]  # Top 10 worst
    
    def __repr__(self):
        return (
            f"CollisionStatistics(\n"
            f"  total_pairs={self.total_pairs},\n"
            f"  collision_pairs={self.collision_pairs},\n"
            f"  collision_rate={self.collision_rate:.2%},\n"
            f"  min_distance={self.min_distance:.4f},\n"
            f"  max_distance={self.max_distance:.4f},\n"
            f"  mean_distance={self.mean_distance:.4f},\n"
            f"  median_distance={self.median_distance:.4f},\n"
            f"  std_distance={self.std_distance:.4f}\n"
            f")"
        )


def identify_collisions(
    pairwise_distances: np.ndarray,
    point_indices: np.ndarray,
    points: np.ndarray,
    links: List,
    threshold_multiplier: float = 1.0
) -> List[CollisionPair]:
    """
    Identify all collision pairs where distance < threshold.
    
    Parameters:
    -----------
    pairwise_distances : np.ndarray
        (L, L) matrix of pairwise distances
    point_indices : np.ndarray
        (L, L) matrix of point indices
    points : np.ndarray
        (N, 3) array of sampled points
    links : List[LinkCSC]
        List of links
    threshold_multiplier : float
        Multiply link radius by this to get collision threshold
    
    Returns:
    --------
    collisions : List[CollisionPair]
        List of collision pairs
    """
    L = pairwise_distances.shape[0]
    collisions = []
    
    for i in range(L):
        for j in range(i + 1, L):  # Only upper triangle
            dist = pairwise_distances[i, j]
            
            # Check if this is a collision
            # Distance should be less than the radius of the link
            # (negative distance means overlap)
            threshold = threshold_multiplier * links[i].r
            
            if dist < threshold and not np.isinf(dist):
                point_idx = point_indices[i, j]
                
                if point_idx >= 0 and point_idx < len(points):
                    point_loc = points[point_idx]
                else:
                    point_loc = np.array([np.nan, np.nan, np.nan])
                
                collision = CollisionPair(
                    link_i=i,
                    link_j=j,
                    distance=float(dist),
                    point_idx=int(point_idx),
                    point_location=point_loc
                )
                collisions.append(collision)
    
    return collisions


def compute_collision_statistics(
    collisions: List[CollisionPair],
    total_pairs: int
) -> CollisionStatistics:
    """
    Compute statistics about collisions.
    
    Parameters:
    -----------
    collisions : List[CollisionPair]
        List of collision pairs
    total_pairs : int
        Total number of link pairs
    
    Returns:
    --------
    stats : CollisionStatistics
        Statistics about collisions
    """
    if len(collisions) == 0:
        return CollisionStatistics(
            total_pairs=total_pairs,
            collision_pairs=0,
            collision_rate=0.0,
            min_distance=np.inf,
            max_distance=-np.inf,
            mean_distance=np.nan,
            median_distance=np.nan,
            std_distance=np.nan,
            worst_collisions=[]
        )
    
    distances = np.array([c.distance for c in collisions])
    
    # Sort by distance (most negative = worst collision)
    sorted_collisions = sorted(collisions, key=lambda c: c.distance)
    worst_10 = sorted_collisions[:min(10, len(sorted_collisions))]
    
    return CollisionStatistics(
        total_pairs=total_pairs,
        collision_pairs=len(collisions),
        collision_rate=len(collisions) / total_pairs if total_pairs > 0 else 0.0,
        min_distance=float(np.min(distances)),
        max_distance=float(np.max(distances)),
        mean_distance=float(np.mean(distances)),
        median_distance=float(np.median(distances)),
        std_distance=float(np.std(distances)),
        worst_collisions=worst_10
    )


def compare_collision_detections(
    collisions_a: List[CollisionPair],
    collisions_b: List[CollisionPair],
    rtol: float = 1e-4,
    atol: float = 1e-6
) -> Tuple[bool, str]:
    """
    Compare collision detections from two methods.
    
    Parameters:
    -----------
    collisions_a : List[CollisionPair]
        Collisions from method A
    collisions_b : List[CollisionPair]
        Collisions from method B
    rtol : float
        Relative tolerance for distance comparison
    atol : float
        Absolute tolerance for distance comparison
    
    Returns:
    --------
    match : bool
        Whether collision detections match
    message : str
        Description of match/mismatch
    """
    # Create sets of (link_i, link_j) pairs for quick lookup
    pairs_a = {(c.link_i, c.link_j) for c in collisions_a}
    pairs_b = {(c.link_i, c.link_j) for c in collisions_b}
    
    # Check if sets match
    if pairs_a != pairs_b:
        only_in_a = pairs_a - pairs_b
        only_in_b = pairs_b - pairs_a
        
        msg = f"Collision pair mismatch: "
        if only_in_a:
            msg += f"{len(only_in_a)} only in A, "
        if only_in_b:
            msg += f"{len(only_in_b)} only in B"
        
        return False, msg
    
    # Check if distances match for common pairs
    collision_dict_a = {(c.link_i, c.link_j): c for c in collisions_a}
    collision_dict_b = {(c.link_i, c.link_j): c for c in collisions_b}
    
    for pair in pairs_a:
        dist_a = collision_dict_a[pair].distance
        dist_b = collision_dict_b[pair].distance
        
        if not np.isclose(dist_a, dist_b, rtol=rtol, atol=atol):
            return False, f"Distance mismatch for pair {pair}: {dist_a:.6f} vs {dist_b:.6f}"
    
    return True, f"All {len(collisions_a)} collisions match"


def filter_intersecting_pairs(
    pairwise_distances: np.ndarray,
    point_indices: np.ndarray,
    links: List,
    threshold_multiplier: float = 1.0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Filter pairwise results to only include intersecting pairs.
    
    Parameters:
    -----------
    pairwise_distances : np.ndarray
        (L, L) matrix of pairwise distances
    point_indices : np.ndarray
        (L, L) matrix of point indices
    links : List[LinkCSC]
        List of links
    threshold_multiplier : float
        Multiply link radius by this to get intersection threshold
    
    Returns:
    --------
    filtered_distances : np.ndarray
        (L, L) matrix with only intersecting pairs (others set to inf)
    filtered_indices : np.ndarray
        (L, L) matrix with only intersecting pairs (others set to -1)
    intersection_mask : np.ndarray
        (L, L) boolean mask of intersecting pairs
    """
    L = pairwise_distances.shape[0]
    intersection_mask = np.zeros((L, L), dtype=bool)
    
    for i in range(L):
        for j in range(L):
            if i == j:
                continue
            
            dist = pairwise_distances[i, j]
            threshold = threshold_multiplier * links[i].r
            
            if dist < threshold and not np.isinf(dist):
                intersection_mask[i, j] = True
    
    filtered_distances = np.where(intersection_mask, pairwise_distances, np.inf)
    filtered_indices = np.where(intersection_mask, point_indices, -1)
    
    return filtered_distances, filtered_indices, intersection_mask


def print_collision_report(
    method_name: str,
    stats: CollisionStatistics,
    show_worst: int = 5
):
    """Print a formatted collision report"""
    print(f"\n{method_name} Collision Report:")
    print(f"  Total link pairs: {stats.total_pairs}")
    print(f"  Collision pairs: {stats.collision_pairs} ({stats.collision_rate:.2%})")
    
    if stats.collision_pairs > 0:
        print(f"  Distance range: [{stats.min_distance:.4f}, {stats.max_distance:.4f}]")
        print(f"  Mean distance: {stats.mean_distance:.4f}")
        print(f"  Median distance: {stats.median_distance:.4f}")
        print(f"  Std deviation: {stats.std_distance:.4f}")
        
        if show_worst > 0 and stats.worst_collisions:
            print(f"\n  Top {min(show_worst, len(stats.worst_collisions))} worst collisions:")
            for idx, collision in enumerate(stats.worst_collisions[:show_worst]):
                print(f"    {idx+1}. Links ({collision.link_i}, {collision.link_j}): "
                      f"dist = {collision.distance:.4f}")
