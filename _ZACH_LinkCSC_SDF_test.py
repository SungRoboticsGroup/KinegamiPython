#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test file to visualize LinkCSC SDF (Signed Distance Function)

Creates random CSC paths and visualizes the SDF values in 3D space
using a grid of sample points colored by distance.

@author: Zach (AI Assistant)
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from spatialmath import SE3, SO3
from LinkCSC import LinkCSC

def random_unit_vector():
    """Generate a random unit vector in 3D"""
    vec = np.random.randn(3)
    return vec / np.linalg.norm(vec)

def random_SE3_pose():
    """Generate a random SE3 pose with random position and orientation"""
    # Random position in a reasonable range
    position = np.random.uniform(-5, 5, 3)
    
    # Random orientation using random rotation
    axis = random_unit_vector()
    angle = np.random.uniform(0, 2*np.pi)
    rotation = SO3.AngleAxis(angle, axis)
    
    return SE3.Rt(rotation, position)

def create_random_link(r=1.0, min_separation=3.0, max_attempts=10):
    """Create a random LinkCSC with reasonable parameters"""
    for attempt in range(max_attempts):
        try:
            start_pose = random_SE3_pose()
            
            # Create end pose with some minimum separation to ensure valid path
            separation = np.random.uniform(min_separation, min_separation + 3)
            end_position = start_pose.t + random_unit_vector() * separation
            end_direction = random_unit_vector()
            
            # Create rotation matrix with end_direction as first column
            v1 = end_direction
            v2_temp = random_unit_vector()
            v2 = v2_temp - np.dot(v2_temp, v1) * v1
            v2 = v2 / np.linalg.norm(v2)
            v3 = np.cross(v1, v2)
            
            R_end = np.column_stack([v1, v2, v3])
            end_pose = SE3.Rt(SO3(R_end), end_position)
            
            # Use smaller maxAnglePerElbow to avoid theta >= pi issues
            link = LinkCSC(r, start_pose, end_pose, maxAnglePerElbow=np.pi/3)
            
            return link
        except Exception as e:
            if attempt == max_attempts - 1:
                raise RuntimeError(f"Failed to create valid link after {max_attempts} attempts: {e}")
            continue
    
    raise RuntimeError("Failed to create valid link")

def create_sdf_grid(link, grid_resolution=20, padding=2.0):
    """
    Create a 3D grid of points around the link and evaluate SDF at each point.
    
    Parameters:
    -----------
    link : LinkCSC
        The link to evaluate
    grid_resolution : int
        Number of points along each axis
    padding : float
        Extra space around the link bounding box
        
    Returns:
    --------
    grid_points : np.ndarray
        Array of 3D points (N, 3)
    sdf_values : np.ndarray
        SDF value at each point (N,)
    bounds : tuple
        ((xmin, xmax), (ymin, ymax), (zmin, zmax))
    """
    # Get bounding box from link endpoints and path
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
    
    # Compute bounds with padding
    mins = points.min(axis=0) - padding
    maxs = points.max(axis=0) + padding
    
    # Create grid
    x = np.linspace(mins[0], maxs[0], grid_resolution)
    y = np.linspace(mins[1], maxs[1], grid_resolution)
    z = np.linspace(mins[2], maxs[2], grid_resolution)
    
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
    grid_points = np.column_stack([xx.ravel(), yy.ravel(), zz.ravel()])
    
    # Evaluate SDF at each point
    print(f"Evaluating SDF at {len(grid_points)} points...")
    sdf_values = np.array([link.sdf(point) for point in grid_points])
    print(f"SDF evaluation complete. Min: {sdf_values.min():.3f}, Max: {sdf_values.max():.3f}")
    
    bounds = ((mins[0], maxs[0]), (mins[1], maxs[1]), (mins[2], maxs[2]))
    return grid_points, sdf_values, bounds

def plot_sdf_slice(link, ax, grid_points, sdf_values, slice_axis='z', slice_value=None):
    """
    Plot a 2D slice of the SDF field.
    
    Parameters:
    -----------
    link : LinkCSC
        The link being visualized
    ax : matplotlib axis
        2D axis to plot on
    grid_points : np.ndarray
        Grid points (N, 3)
    sdf_values : np.ndarray
        SDF values (N,)
    slice_axis : str
        'x', 'y', or 'z' - which axis to slice along
    slice_value : float, optional
        Value along slice_axis to slice at. If None, uses middle of range
    """
    axis_map = {'x': 0, 'y': 1, 'z': 2}
    slice_idx = axis_map[slice_axis]
    
    # Other two axes
    other_axes = [i for i in range(3) if i != slice_idx]
    
    if slice_value is None:
        slice_value = (grid_points[:, slice_idx].min() + grid_points[:, slice_idx].max()) / 2
    
    # Find points near the slice
    tolerance = (grid_points[:, slice_idx].max() - grid_points[:, slice_idx].min()) / 50
    slice_mask = np.abs(grid_points[:, slice_idx] - slice_value) < tolerance
    
    slice_points = grid_points[slice_mask]
    slice_sdf = sdf_values[slice_mask]
    
    # Plot
    scatter = ax.scatter(slice_points[:, other_axes[0]], 
                        slice_points[:, other_axes[1]], 
                        c=slice_sdf, 
                        cmap='RdBu', 
                        vmin=-link.r*2, vmax=link.r*2,
                        s=50, alpha=0.6)
    
    # Plot link path on this slice
    path_points = link.interpolate(count=100)
    ax.plot(path_points[:, other_axes[0]], 
           path_points[:, other_axes[1]], 
           'k-', linewidth=2, label='Link Path')
    
    axis_names = ['X', 'Y', 'Z']
    ax.set_xlabel(axis_names[other_axes[0]])
    ax.set_ylabel(axis_names[other_axes[1]])
    ax.set_title(f'SDF Slice at {axis_names[slice_idx]}={slice_value:.2f}')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    
    return scatter

def plot_sdf_3d(link, ax, grid_points, sdf_values, threshold=None):
    """
    Plot 3D visualization of SDF values.
    
    Parameters:
    -----------
    link : LinkCSC
        The link being visualized
    ax : matplotlib 3D axis
        3D axis to plot on
    grid_points : np.ndarray
        Grid points (N, 3)
    sdf_values : np.ndarray
        SDF values (N,)
    threshold : float, optional
        Only show points within this distance of surface
    """
    if threshold is None:
        threshold = link.r * 1.5
    
    # Filter points near surface
    near_surface = np.abs(sdf_values) < threshold
    plot_points = grid_points[near_surface]
    plot_sdf = sdf_values[near_surface]
    
    # Color by SDF value
    scatter = ax.scatter(plot_points[:, 0], 
                        plot_points[:, 1], 
                        plot_points[:, 2],
                        c=plot_sdf, 
                        cmap='RdBu',
                        vmin=-link.r, vmax=link.r,
                        s=20, alpha=0.4)
    
    # Plot the actual link
    link.addToPlot(ax, numSides=16, color='gray', alpha=0.7, 
                  wireFrame=False, showPath=True, showBoundary=True)
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'SDF 3D (|dist| < {threshold:.2f})')
    
    return scatter

def test_sdf_properties(link):
    """
    Test key SDF properties:
    1. Points on the path should have distance ≈ -r
    2. SDF should be continuous
    3. Zero-crossing should be near r from centerline
    """
    print("\n" + "="*60)
    print("Testing SDF Properties")
    print("="*60)
    
    # Test 1: Points on centerline should have distance ≈ -r
    print("\nTest 1: Centerline distances")
    centerline_points = link.interpolate(count=20)
    centerline_distances = np.array([link.sdf(p) for p in centerline_points])
    expected = -link.r
    print(f"  Expected distance at centerline: {expected:.3f}")
    print(f"  Actual distances - Mean: {centerline_distances.mean():.3f}, Std: {centerline_distances.std():.3f}")
    print(f"  Min: {centerline_distances.min():.3f}, Max: {centerline_distances.max():.3f}")
    
    # Test 2: Points exactly at radius should be near zero
    print("\nTest 2: Surface distances (should be near 0)")
    if link.arc1:
        # Sample point on arc1 surface
        t_sample = 0.5
        center_point = link.arc1.interpolateAt(t_sample)
        # Move radially outward by r
        radial_dir = (center_point - link.arc1.circleCenter) / np.linalg.norm(center_point - link.arc1.circleCenter)
        surface_point = link.arc1.circleCenter + (link.arc1.r + link.r) * radial_dir
        surface_dist = link.sdf(surface_point)
        print(f"  Arc1 surface point distance: {surface_dist:.6f}")
    
    # Test on straight section
    if link.path.tMag > link.DISTANCE_EPSILON:
        mid_point = link.path.turn1end + 0.5 * link.path.tMag * link.path.tUnit
        # Move perpendicular by r
        perp_dir = np.array([0, 0, 1]) if abs(link.path.tUnit[2]) < 0.9 else np.array([1, 0, 0])
        perp_dir = perp_dir - np.dot(perp_dir, link.path.tUnit) * link.path.tUnit
        perp_dir = perp_dir / np.linalg.norm(perp_dir)
        surface_point = mid_point + link.r * perp_dir
        surface_dist = link.sdf(surface_point)
        print(f"  Straight section surface point distance: {surface_dist:.6f}")
    
    # Test 3: Symmetry - opposite sides should have same absolute distance
    print("\nTest 3: Symmetry test")
    if link.path.tMag > link.DISTANCE_EPSILON:
        mid_point = link.path.turn1end + 0.5 * link.path.tMag * link.path.tUnit
        perp_dir = np.array([0, 0, 1]) if abs(link.path.tUnit[2]) < 0.9 else np.array([1, 0, 0])
        perp_dir = perp_dir - np.dot(perp_dir, link.path.tUnit) * link.path.tUnit
        perp_dir = perp_dir / np.linalg.norm(perp_dir)
        
        offset = 2.0 * link.r
        point_plus = mid_point + offset * perp_dir
        point_minus = mid_point - offset * perp_dir
        dist_plus = link.sdf(point_plus)
        dist_minus = link.sdf(point_minus)
        print(f"  Distance at +offset: {dist_plus:.6f}")
        print(f"  Distance at -offset: {dist_minus:.6f}")
        print(f"  Difference: {abs(dist_plus - dist_minus):.6f} (should be near 0)")

def main():
    """Create and visualize SDF for random LinkCSC paths"""
    np.random.seed(28)  # For reproducibility
    
    num_links = 3
    
    for link_idx in range(num_links):
        print(f"\n{'='*60}")
        print(f"Creating and testing Link {link_idx + 1}")
        print(f"{'='*60}")
        
        try:
            # Create random link
            link = create_random_link(r=1.0, min_separation=3.0)
            print(f"Link created successfully!")
            print(f"  r = {link.r:.2f}")
            print(f"  theta1 = {np.rad2deg(link.path.theta1):.1f}°")
            print(f"  theta2 = {np.rad2deg(link.path.theta2):.1f}°")
            print(f"  straight length = {link.path.tMag:.2f}")
            
            # Test SDF properties
            test_sdf_properties(link)
            
            # Create SDF grid
            print("\nGenerating SDF grid...")
            grid_points, sdf_values, bounds = create_sdf_grid(link, grid_resolution=15, padding=2.0)
            
            # Create visualization
            fig = plt.figure(figsize=(16, 10))
            fig.suptitle(f'Link {link_idx + 1} - SDF Visualization', fontsize=16, fontweight='bold')
            
            # 3D plot
            ax1 = fig.add_subplot(2, 2, 1, projection='3d')
            scatter1 = plot_sdf_3d(link, ax1, grid_points, sdf_values, threshold=link.r*1.5)
            plt.colorbar(scatter1, ax=ax1, label='SDF Value', shrink=0.5)
            
            # 2D slices
            ax2 = fig.add_subplot(2, 2, 2)
            scatter2 = plot_sdf_slice(link, ax2, grid_points, sdf_values, slice_axis='z')
            plt.colorbar(scatter2, ax=ax2, label='SDF Value')
            
            ax3 = fig.add_subplot(2, 2, 3)
            scatter3 = plot_sdf_slice(link, ax3, grid_points, sdf_values, slice_axis='y')
            plt.colorbar(scatter3, ax=ax3, label='SDF Value')
            
            ax4 = fig.add_subplot(2, 2, 4)
            scatter4 = plot_sdf_slice(link, ax4, grid_points, sdf_values, slice_axis='x')
            plt.colorbar(scatter4, ax=ax4, label='SDF Value')
            
            plt.tight_layout()
            plt.show()
            
            print(f"\nLink {link_idx + 1} visualization complete!")
            
        except Exception as e:
            print(f"Error with link {link_idx + 1}: {e}")
            import traceback
            traceback.print_exc()
    
    print("\n" + "="*60)
    print("All tests complete!")
    print("="*60)

if __name__ == "__main__":
    main()
