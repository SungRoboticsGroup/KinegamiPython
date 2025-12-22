#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Sphere Ray Marching visualization for LinkCSC SDF

Uses ray marching (sphere tracing) to render the SDF surface efficiently.
Creates multiple viewpoints to inspect the 3D surface quality.

@author: Zach (AI Assistant)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from spatialmath import SE3, SO3
from LinkCSC import LinkCSC
import time

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
    origins : np.ndarray
        Ray origins (width*height, 3)
    directions : np.ndarray
        Ray directions (width*height, 3)
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
    world_up = np.array([0, 0, 1])
    if abs(np.dot(forward, world_up)) > 0.99:
        world_up = np.array([1, 0, 0])
    
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

def ray_march(link, origins, directions, max_steps=50, epsilon=0.002, max_distance=100.0):
    """
    Perform sphere ray marching for all rays - OPTIMIZED for speed.
    
    Parameters:
    -----------
    link : LinkCSC
        The link to render
    origins : np.ndarray
        Ray origins (N, 3)
    directions : np.ndarray
        Ray directions (N, 3)
    max_steps : int
        Maximum marching steps per ray (reduced from 100 to 50)
    epsilon : float
        Surface hit threshold (increased from 0.001 to 0.002 for looser tolerance)
    max_distance : float
        Maximum ray distance
        
    Returns:
    --------
    hits : np.ndarray
        Boolean array indicating which rays hit (N,)
    depths : np.ndarray
        Distance along ray to hit point (N,)
    iterations : np.ndarray
        Number of steps taken (N,)
    """
    num_rays = len(origins)
    hits = np.zeros(num_rays, dtype=bool)
    depths = np.full(num_rays, max_distance)
    iterations = np.zeros(num_rays, dtype=int)
    
    # Current position for each ray
    positions = origins.copy()
    distances = np.zeros(num_rays)
    active_mask = np.ones(num_rays, dtype=bool)  # Track which rays are still marching
    
    for step in range(max_steps):
        # Skip iteration if all rays are done
        if not np.any(active_mask):
            break
        
        # Get indices of still-active rays
        active_indices = np.where(active_mask)[0]
        
        # Evaluate SDF for all active rays at once (vectorized would be ideal, but SDF needs individual calls)
        for i in active_indices:
            sdf_value = link.sdf(positions[i])
            
            if abs(sdf_value) < epsilon:
                # Hit surface
                hits[i] = True
                depths[i] = distances[i]
                iterations[i] = step
                active_mask[i] = False  # Mark ray as done
            elif distances[i] >= max_distance:
                # Too far - deactivate
                iterations[i] = step
                active_mask[i] = False
            else:
                # March forward with min step to avoid artifacts
                step_size = max(abs(sdf_value), 0.05)  # Clamp minimum step size
                positions[i] += directions[i] * step_size
                distances[i] += step_size
    
    return hits, depths, iterations

def compute_normals(link, hit_points, epsilon=0.005):
    """
    Compute surface normals at hit points using gradient - OPTIMIZED.
    
    Parameters:
    -----------
    link : LinkCSC
        The link
    hit_points : np.ndarray
        Points on surface (N, 3)
    epsilon : float
        Finite difference epsilon (increased for speed)
        
    Returns:
    --------
    normals : np.ndarray
        Surface normals (N, 3)
    """
    normals = np.zeros_like(hit_points)
    
    # Use faster forward differences instead of central differences
    for i, point in enumerate(hit_points):
        # Forward differences (faster, less SDF calls)
        f0 = link.sdf(point)
        fx = link.sdf(point + np.array([epsilon, 0, 0])) - f0
        fy = link.sdf(point + np.array([0, epsilon, 0])) - f0
        fz = link.sdf(point + np.array([0, 0, epsilon])) - f0
        
        gradient = np.array([fx, fy, fz]) / epsilon
        norm = np.linalg.norm(gradient)
        if norm > 1e-6:
            normals[i] = gradient / norm
        else:
            normals[i] = np.array([0, 0, 1])
    
    return normals

def render_sdf(link, center, distance, azimuth, elevation, resolution=(400, 300), 
               shading='depth', show_progress=True):
    """
    Render the SDF surface using ray marching.
    
    Parameters:
    -----------
    link : LinkCSC
        The link to render
    center : np.ndarray
        Camera look-at point
    distance : float
        Camera distance from center
    azimuth : float
        Horizontal viewing angle (radians)
    elevation : float
        Vertical viewing angle (radians)
    resolution : tuple
        Image resolution (width, height)
    shading : str
        'depth', 'normal', 'iterations', or 'ambient_occlusion'
    show_progress : bool
        Print progress information
        
    Returns:
    --------
    image : np.ndarray
        Rendered image (height, width, 3) or (height, width)
    """
    if show_progress:
        print(f"  Rendering view: azimuth={np.rad2deg(azimuth):.0f}°, elevation={np.rad2deg(elevation):.0f}°")
    
    # Create camera rays
    start_time = time.time()
    origins, directions, (width, height) = create_camera(center, distance, azimuth, elevation, resolution)
    
    if show_progress:
        print(f"    Created {len(origins)} rays")
    
    # Ray march
    hits, depths, iterations = ray_march(link, origins, directions)
    
    ray_time = time.time() - start_time
    if show_progress:
        print(f"    Ray marching: {ray_time:.2f}s, {np.sum(hits)} hits ({100*np.sum(hits)/len(hits):.1f}%)")
    
    # Create image based on shading mode
    image = np.zeros((height, width))
    
    if shading == 'depth':
        # Depth map
        image = depths.reshape(height, width)
        # Normalize and invert (closer = brighter)
        max_depth = depths[hits].max() if np.any(hits) else 1.0
        min_depth = depths[hits].min() if np.any(hits) else 0.0
        image = np.where(hits.reshape(height, width), 
                        1.0 - (image - min_depth) / (max_depth - min_depth + 1e-6),
                        0.0)
    
    elif shading == 'iterations':
        # Number of steps taken
        image = iterations.reshape(height, width)
        max_iter = iterations.max()
        image = image / max_iter
    
    elif shading == 'normal':
        # Compute normals for lighting
        hit_points = origins[hits] + directions[hits] * depths[hits, np.newaxis]
        normals = compute_normals(link, hit_points)
        
        # Simple directional lighting
        light_dir = np.array([1, 1, 1])
        light_dir = light_dir / np.linalg.norm(light_dir)
        
        lighting = np.maximum(0, np.dot(normals, light_dir))
        
        image_flat = np.zeros(len(hits))
        image_flat[hits] = lighting
        image = image_flat.reshape(height, width)
    
    elif shading == 'ambient_occlusion':
        # Approximate ambient occlusion using iteration count
        ao = 1.0 - iterations / iterations.max()
        ao = np.where(hits, ao, 0.0)
        image = ao.reshape(height, width)
    
    return image

def visualize_link_multiview(link, resolution=(400, 300)):
    """
    Create a multi-view visualization of the link using ray marching.
    
    Parameters:
    -----------
    link : LinkCSC
        The link to visualize
    resolution : tuple
        Image resolution for each view
    """
    center, extent = get_link_bounds(link)
    distance = extent * 2.0
    
    # Define multiple viewpoints
    views = [
        (0, 0, 'Front'),              # Front view
        (np.pi/2, 0, 'Side'),         # Side view
        (np.pi, 0, 'Back'),           # Back view
        (np.pi/4, np.pi/6, 'Angle 1'), # Angled view 1
        (3*np.pi/4, np.pi/6, 'Angle 2'), # Angled view 2
        (np.pi/2, np.pi/4, 'Top Side'), # Top-side view
    ]
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle(f'LinkCSC Ray Marching Visualization (r={link.r:.2f})', 
                 fontsize=16, fontweight='bold')
    
    axes = axes.flatten()
    
    for idx, (azimuth, elevation, name) in enumerate(views):
        print(f"\nRendering view {idx+1}/{len(views)}: {name}")
        
        # Render with normal shading
        image = render_sdf(link, center, distance, azimuth, elevation, 
                          resolution, shading='normal')
        
        ax = axes[idx]
        im = ax.imshow(image, cmap='gray', interpolation='bilinear')
        ax.set_title(f'{name}\nAz={np.rad2deg(azimuth):.0f}° El={np.rad2deg(elevation):.0f}°')
        ax.axis('off')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
    plt.tight_layout()
    return fig

def visualize_link_interactive(link, resolution=(300, 300)):
    """
    Create a fast interactive visualization with 4 different shading modes - OPTIMIZED.
    
    Parameters:
    -----------
    link : LinkCSC
        The link to visualize
    resolution : tuple
        Image resolution (reduced from 600 to 300 for 4x speedup)
    """
    center, extent = get_link_bounds(link)
    distance = extent * 2.0
    
    # Single view from nice angle
    azimuth = np.pi / 4
    elevation = np.pi / 6
    
    # Different shading modes - only 4
    shading_modes = ['normal', 'depth', 'iterations', 'ambient_occlusion']
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    fig.suptitle(f'LinkCSC SDF Surface - Shading Modes', 
                 fontsize=16, fontweight='bold')
    
    axes = axes.flatten()
    
    for idx, shading in enumerate(shading_modes):
        print(f"Rendering {shading} shading...")
        
        image = render_sdf(link, center, distance, azimuth, elevation, 
                          resolution, shading=shading, show_progress=False)
        
        ax = axes[idx]
        cmap = 'gray' if shading in ['normal', 'ambient_occlusion'] else 'viridis'
        im = ax.imshow(image, cmap=cmap, interpolation='bilinear')
        ax.set_title(f'{shading.replace("_", " ").title()}')
        ax.axis('off')
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
    plt.tight_layout()
    return fig

def main():
    """Main test function"""
    print("="*60)
    print("LinkCSC SDF Ray Marching Visualization")
    print("="*60)
    
    np.random.seed(42)
    
    # Create random link
    print("\nCreating random link...")
    link = create_random_link(r=1.0, min_separation=4.0)
    
    print(f"\nLink properties:")
    print(f"  r = {link.r:.2f}")
    print(f"  theta1 = {np.rad2deg(link.path.theta1):.1f}°")
    print(f"  theta2 = {np.rad2deg(link.path.theta2):.1f}°")
    print(f"  straight length = {link.path.tMag:.2f}")
    
    # Multi-view visualization
    print("\n" + "="*60)
    print("Creating multi-view visualization...")
    print("="*60)
    fig1 = visualize_link_multiview(link, resolution=(300, 300))
    
    # Interactive shading modes
    # print("\n" + "="*60)
    # print("Creating shading comparison...")
    # print("="*60)
    # fig2 = visualize_link_interactive(link, resolution=(50, 50))
    
    plt.show()
    
    print("\n" + "="*60)
    print("Visualization complete!")
    print("="*60)

if __name__ == "__main__":
    main()
