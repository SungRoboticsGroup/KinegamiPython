#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
GPU-Accelerated Sphere Ray Marching visualization for LinkCSC SDF using PyCUDA

Uses CUDA to parallelize ray marching for massive speedup (100x+ over CPU).
Creates 4 shaded views of the 3D surface for quick inspection.

Requirements:
    pip install pycuda

@author: Zach (AI Assistant)
"""

import numpy as np
import matplotlib.pyplot as plt
from spatialmath import SE3, SO3
from LinkCSC import LinkCSC
import time

try:
    import pycuda.autoinit
    import pycuda.driver as cuda
    from pycuda.compiler import SourceModule
    CUDA_AVAILABLE = True
except ImportError:
    print("WARNING: PyCUDA not available. Install with: pip install pycuda")
    CUDA_AVAILABLE = False

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

def create_camera(center, distance, azimuth, elevation, resolution=(800, 600)):
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
        Ray origins (width*height, 3) as float32
    directions : np.ndarray
        Ray directions (width*height, 3) as float32
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
    
    # Convert to float32 for GPU
    return origins.astype(np.float32), directions.astype(np.float32), (width, height)

def prepare_link_data(link):
    """
    Prepare link geometry data for GPU transfer.
    
    The arc's theta is the FULL arc angle (from start to end point), matching
    the Arc3D class definition. The GPU kernel internally handles the symmetric
    transformation by using theta/2 for the capped torus SDF.
    
    The local frame has X-axis pointing to the arc's midpoint for symmetry.
    
    Returns:
    --------
    dict with link geometry parameters as float32 arrays
    """
    from scipy.spatial.transform import Rotation
    
    data = {
        'r': np.float32(link.r),
        'has_arc1': np.int32(1 if link.arc1 is not None and link.path.theta1 > link.EPSILON else 0),
        'has_arc2': np.int32(1 if link.arc2 is not None and link.path.theta2 > link.EPSILON else 0),
        'has_straight': np.int32(1 if link.path.tMag > link.DISTANCE_EPSILON else 0),
    }
    
    # Arc 1 data
    if data['has_arc1']:
        data['arc1_center'] = link.path.circleCenter1.astype(np.float32)
        data['arc1_radius'] = np.float32(link.arc1.r)
        # Pass FULL theta - kernel will halve it internally
        data['arc1_theta'] = np.float32(link.arc1.theta)
        
        # Local frame: X-axis at arc MIDPOINT for symmetry with abs(x) in SDF
        halfAngleRot = Rotation.from_rotvec(link.arc1.theta / 2 * link.arc1.binormal)
        centerToMid = halfAngleRot.apply(-link.arc1.startNormal)
        arcX = centerToMid / np.linalg.norm(centerToMid)
        arcZ = link.arc1.binormal
        arcY = np.cross(arcZ, arcX)
        data['arc1_frame'] = np.column_stack([arcX, arcY, arcZ]).astype(np.float32).flatten()
    else:
        data['arc1_center'] = link.StartDubinsPose.t.astype(np.float32)
        data['arc1_radius'] = np.float32(link.r)
        data['arc1_theta'] = np.float32(0)
        data['arc1_frame'] = np.eye(3, dtype=np.float32).flatten()
    
    # Straight section data
    if data['has_straight']:
        data['straight_start'] = link.path.turn1end.astype(np.float32)
        data['straight_end'] = (link.path.turn1end + link.path.tMag * link.path.tUnit).astype(np.float32)
    else:
        data['straight_start'] = link.path.turn1end.astype(np.float32)
        data['straight_end'] = link.path.turn1end.astype(np.float32)
    
    # Arc 2 data
    if data['has_arc2']:
        data['arc2_center'] = link.path.circleCenter2.astype(np.float32)
        data['arc2_radius'] = np.float32(link.arc2.r)
        # Pass FULL theta - kernel will halve it internally
        data['arc2_theta'] = np.float32(link.arc2.theta)
        
        # Local frame: X-axis at arc MIDPOINT for symmetry with abs(x) in SDF
        halfAngleRot = Rotation.from_rotvec(link.arc2.theta / 2 * link.arc2.binormal)
        centerToMid = halfAngleRot.apply(-link.arc2.startNormal)
        arcX = centerToMid / np.linalg.norm(centerToMid)
        arcZ = link.arc2.binormal
        arcY = np.cross(arcZ, arcX)
        data['arc2_frame'] = np.column_stack([arcX, arcY, arcZ]).astype(np.float32).flatten()
    else:
        data['arc2_center'] = link.EndDubinsPose.t.astype(np.float32)
        data['arc2_radius'] = np.float32(link.r)
        data['arc2_theta'] = np.float32(0)
        data['arc2_frame'] = np.eye(3, dtype=np.float32).flatten()
    
    return data

# CUDA kernel for ray marching
CUDA_KERNEL = """
__device__ float length3(float x, float y, float z) {
    return sqrtf(x*x + y*y + z*z);
}

__device__ float dot3(float ax, float ay, float az, float bx, float by, float bz) {
    return ax*bx + ay*by + az*bz;
}

__device__ float sdf_capsule(float px, float py, float pz,
                             float ax, float ay, float az,
                             float bx, float by, float bz,
                             float r) {
    float pax = px - ax;
    float pay = py - ay;
    float paz = pz - az;
    
    float bax = bx - ax;
    float bay = by - ay;
    float baz = bz - az;
    
    float ba_dot_ba = dot3(bax, bay, baz, bax, bay, baz);
    float h;
    if (ba_dot_ba < 1e-10f) {
        // Degenerate case: start == end, just a sphere
        h = 0.0f;
    } else {
        h = dot3(pax, pay, paz, bax, bay, baz) / ba_dot_ba;
        h = fmaxf(0.0f, fminf(h, 1.0f));
    }
    
    float cx = pax - bax * h;
    float cy = pay - bay * h;
    float cz = paz - baz * h;
    
    return length3(cx, cy, cz) - r;
}

// Capped torus SDF in local coordinates
// The arc is symmetric about the X-axis, spanning [-halfTheta, +halfTheta]
// X-axis points radially outward at the arc's midpoint
__device__ float sdf_capped_torus(float px, float py, float pz,
                                   float sin_half_theta, float cos_half_theta,
                                   float ra, float rb) {
    // abs(x) exploits the symmetry about the X-axis
    float abs_px = fabsf(px);
    float k;
    
    // Check if point's angle from +X is within halfTheta
    if (cos_half_theta * abs_px > sin_half_theta * py) {
        // Inside arc span - project onto arc centerline
        k = sin_half_theta * abs_px + cos_half_theta * py;
    } else {
        // Outside arc span - closest point is at cap (endpoint)
        k = sqrtf(abs_px * abs_px + py * py);
    }
    
    return sqrtf(px*px + py*py + pz*pz + ra*ra - 2.0f*ra*k) - rb;
}

// Transform point to local arc coordinates and evaluate capped torus SDF
// frame is the world-to-local rotation matrix (row-major)
// arc_theta is the FULL arc angle - we halve it internally for the symmetric SDF
__device__ float sdf_arc(float px, float py, float pz,
                        float cx, float cy, float cz,
                        float* frame, float arc_r, float arc_theta, float tube_r) {
    // Translate to arc center
    float lx = px - cx;
    float ly = py - cy;
    float lz = pz - cz;
    
    // Apply frame transformation (frame is row-major 3x3)
    // Each row of frame is a local basis vector
    float lpx = frame[0]*lx + frame[1]*ly + frame[2]*lz;
    float lpy = frame[3]*lx + frame[4]*ly + frame[5]*lz;
    float lpz = frame[6]*lx + frame[7]*ly + frame[8]*lz;
    
    // Halve theta internally for the symmetric capped torus SDF
    float half_theta = arc_theta * 0.5f;
    return sdf_capped_torus(lpx, lpy, lpz, sinf(half_theta), cosf(half_theta), arc_r, tube_r);
}

__device__ float sdf_sphere(float px, float py, float pz,
                           float cx, float cy, float cz, float r) {
    return length3(px - cx, py - cy, pz - cz) - r;
}

__device__ float link_sdf(float px, float py, float pz,
                         int has_arc1, float* arc1_center, float* arc1_frame,
                         float arc1_r, float arc1_theta,
                         int has_straight, float* straight_start, float* straight_end,
                         int has_arc2, float* arc2_center, float* arc2_frame,
                         float arc2_r, float arc2_theta,
                         float tube_r) {
    float dist = 1e10f;
    
    // Arc 1
    if (has_arc1) {
        float d1 = sdf_arc(px, py, pz, arc1_center[0], arc1_center[1], arc1_center[2],
                          arc1_frame, arc1_r, arc1_theta, tube_r);
        dist = fminf(dist, d1);
    } else {
        float d1 = sdf_sphere(px, py, pz, arc1_center[0], arc1_center[1], arc1_center[2], tube_r);
        dist = fminf(dist, d1);
    }
    
    // Straight section
    if (has_straight) {
        float d2 = sdf_capsule(px, py, pz,
                              straight_start[0], straight_start[1], straight_start[2],
                              straight_end[0], straight_end[1], straight_end[2],
                              tube_r);
        dist = fminf(dist, d2);
    }
    
    // Arc 2
    if (has_arc2) {
        float d3 = sdf_arc(px, py, pz, arc2_center[0], arc2_center[1], arc2_center[2],
                          arc2_frame, arc2_r, arc2_theta, tube_r);
        dist = fminf(dist, d3);
    } else {
        float d3 = sdf_sphere(px, py, pz, arc2_center[0], arc2_center[1], arc2_center[2], tube_r);
        dist = fminf(dist, d3);
    }
    
    return dist;
}

__global__ void ray_march_kernel(
    float* origins, float* directions,
    int has_arc1, float* arc1_center, float* arc1_frame, float arc1_r, float arc1_theta,
    int has_straight, float* straight_start, float* straight_end,
    int has_arc2, float* arc2_center, float* arc2_frame, float arc2_r, float arc2_theta,
    float tube_r,
    int* hits, float* depths, int* iterations,
    int num_rays, int max_steps, float epsilon, float max_distance)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= num_rays) return;
    
    // Ray data
    float ox = origins[idx * 3 + 0];
    float oy = origins[idx * 3 + 1];
    float oz = origins[idx * 3 + 2];
    
    float dx = directions[idx * 3 + 0];
    float dy = directions[idx * 3 + 1];
    float dz = directions[idx * 3 + 2];
    
    // Ray march
    float px = ox, py = oy, pz = oz;
    float distance = 0.0f;
    int hit = 0;
    int iter = 0;
    
    for (int step = 0; step < max_steps; step++) {
        iter = step;
        
        float sdf_val = link_sdf(px, py, pz,
                                has_arc1, arc1_center, arc1_frame, arc1_r, arc1_theta,
                                has_straight, straight_start, straight_end,
                                has_arc2, arc2_center, arc2_frame, arc2_r, arc2_theta,
                                tube_r);
        
        if (fabsf(sdf_val) < epsilon) {
            hit = 1;
            break;
        }
        
        if (distance >= max_distance) {
            break;
        }
        
        // March forward - use smaller minimum step for accuracy near surface
        float step_size = fmaxf(sdf_val, 0.001f);
        px += dx * step_size;
        py += dy * step_size;
        pz += dz * step_size;
        distance += step_size;
    }
    
    hits[idx] = hit;
    depths[idx] = distance;
    iterations[idx] = iter;
}
"""

def ray_march_gpu(link, origins, directions, max_steps=200, epsilon=0.001, max_distance=100.0):
    """
    Perform GPU-accelerated ray marching.
    
    Returns:
    --------
    hits : np.ndarray
        Boolean array indicating which rays hit (N,)
    depths : np.ndarray
        Distance along ray to hit point (N,)
    iterations : np.ndarray
        Number of steps taken (N,)
    """
    if not CUDA_AVAILABLE:
        raise RuntimeError("CUDA not available")
    
    num_rays = len(origins)
    
    # Prepare link data
    link_data = prepare_link_data(link)
    
    # Allocate GPU memory
    origins_gpu = cuda.mem_alloc(origins.nbytes)
    directions_gpu = cuda.mem_alloc(directions.nbytes)
    
    # Link geometry on GPU
    arc1_center_gpu = cuda.mem_alloc(link_data['arc1_center'].nbytes)
    arc1_frame_gpu = cuda.mem_alloc(link_data['arc1_frame'].nbytes)
    straight_start_gpu = cuda.mem_alloc(link_data['straight_start'].nbytes)
    straight_end_gpu = cuda.mem_alloc(link_data['straight_end'].nbytes)
    arc2_center_gpu = cuda.mem_alloc(link_data['arc2_center'].nbytes)
    arc2_frame_gpu = cuda.mem_alloc(link_data['arc2_frame'].nbytes)
    
    # Output arrays
    hits = np.zeros(num_rays, dtype=np.int32)
    depths = np.zeros(num_rays, dtype=np.float32)
    iterations_out = np.zeros(num_rays, dtype=np.int32)
    
    hits_gpu = cuda.mem_alloc(hits.nbytes)
    depths_gpu = cuda.mem_alloc(depths.nbytes)
    iterations_gpu = cuda.mem_alloc(iterations_out.nbytes)
    
    # Copy data to GPU
    cuda.memcpy_htod(origins_gpu, origins)
    cuda.memcpy_htod(directions_gpu, directions)
    cuda.memcpy_htod(arc1_center_gpu, link_data['arc1_center'])
    cuda.memcpy_htod(arc1_frame_gpu, link_data['arc1_frame'])
    cuda.memcpy_htod(straight_start_gpu, link_data['straight_start'])
    cuda.memcpy_htod(straight_end_gpu, link_data['straight_end'])
    cuda.memcpy_htod(arc2_center_gpu, link_data['arc2_center'])
    cuda.memcpy_htod(arc2_frame_gpu, link_data['arc2_frame'])
    
    # Compile kernel
    mod = SourceModule(CUDA_KERNEL)
    kernel = mod.get_function("ray_march_kernel")
    
    # Launch kernel
    block_size = 256
    grid_size = (num_rays + block_size - 1) // block_size
    
    kernel(
        origins_gpu, directions_gpu,
        link_data['has_arc1'], arc1_center_gpu, arc1_frame_gpu,
        link_data['arc1_radius'], link_data['arc1_theta'],
        link_data['has_straight'], straight_start_gpu, straight_end_gpu,
        link_data['has_arc2'], arc2_center_gpu, arc2_frame_gpu,
        link_data['arc2_radius'], link_data['arc2_theta'],
        link_data['r'],
        hits_gpu, depths_gpu, iterations_gpu,
        np.int32(num_rays), np.int32(max_steps), np.float32(epsilon), np.float32(max_distance),
        block=(block_size, 1, 1), grid=(grid_size, 1)
    )
    
    # Copy results back
    cuda.memcpy_dtoh(hits, hits_gpu)
    cuda.memcpy_dtoh(depths, depths_gpu)
    cuda.memcpy_dtoh(iterations_out, iterations_gpu)
    
    return hits.astype(bool), depths, iterations_out

def compute_normals_cpu(link, hit_points, epsilon=0.001):
    """
    Compute surface normals at hit points using central differences for better accuracy.
    """
    normals = np.zeros_like(hit_points)
    
    for i, point in enumerate(hit_points):
        # Use central differences for better accuracy
        fx = link.sdf(point + np.array([epsilon, 0, 0])) - link.sdf(point - np.array([epsilon, 0, 0]))
        fy = link.sdf(point + np.array([0, epsilon, 0])) - link.sdf(point - np.array([0, epsilon, 0]))
        fz = link.sdf(point + np.array([0, 0, epsilon])) - link.sdf(point - np.array([0, 0, epsilon]))
        
        gradient = np.array([fx, fy, fz]) / (2 * epsilon)
        norm = np.linalg.norm(gradient)
        if norm > 1e-6:
            normals[i] = gradient / norm
        else:
            normals[i] = np.array([0, 0, 1])
    
    return normals

def render_sdf_gpu(link, center, distance, azimuth, elevation, resolution=(800, 600), 
                   shading='normal', light_dir=None, show_progress=True):
    """
    Render the SDF surface using GPU-accelerated ray marching.
    """
    if show_progress:
        print(f"  Rendering: azimuth={np.rad2deg(azimuth):.0f}°, elevation={np.rad2deg(elevation):.0f}°")
    
    # Create camera rays
    start_time = time.time()
    origins, directions, (width, height) = create_camera(center, distance, azimuth, elevation, resolution)
    
    if show_progress:
        print(f"    Created {len(origins)} rays")
    
    # GPU ray march
    hits, depths, iterations = ray_march_gpu(link, origins, directions)
    
    ray_time = time.time() - start_time
    if show_progress:
        print(f"    GPU ray marching: {ray_time:.3f}s, {np.sum(hits)} hits ({100*np.sum(hits)/len(hits):.1f}%)")
    
    # Create image based on shading mode
    image = np.zeros((height, width, 3))  # RGB image
    
    if shading == 'depth':
        if np.any(hits):
            depth_img = depths.reshape(height, width)
            max_depth = depths[hits].max()
            min_depth = depths[hits].min()
            normalized = np.where(hits.reshape(height, width), 
                            1.0 - (depth_img - min_depth) / (max_depth - min_depth + 1e-6),
                            0.0)
            image[:,:,0] = image[:,:,1] = image[:,:,2] = normalized
    
    elif shading == 'normal':
        # Compute normals for Lambert shading
        hit_indices = np.where(hits)[0]
        if len(hit_indices) > 0:
            hit_origins = origins[hit_indices]
            hit_directions = directions[hit_indices]
            hit_depths = depths[hit_indices]
            hit_points = hit_origins + hit_directions * hit_depths[:, np.newaxis]
            
            normals = compute_normals_cpu(link, hit_points)
            
            # Default light from camera direction + up-right
            if light_dir is None:
                # Camera position
                cam_x = distance * np.cos(elevation) * np.cos(azimuth)
                cam_y = distance * np.cos(elevation) * np.sin(azimuth)
                cam_z = distance * np.sin(elevation)
                cam_dir = np.array([cam_x, cam_y, cam_z])
                cam_dir = cam_dir / np.linalg.norm(cam_dir)
                # Add some light from above-right
                light_dir = cam_dir + np.array([0.3, 0.3, 0.5])
                light_dir = light_dir / np.linalg.norm(light_dir)
            
            # Lambert shading with ambient
            ambient = 0.2
            diffuse = np.maximum(0, np.dot(normals, light_dir))
            lighting = ambient + (1 - ambient) * diffuse
            
            # Create grayscale image
            image_flat = np.zeros((len(hits), 3))
            image_flat[hit_indices, 0] = lighting
            image_flat[hit_indices, 1] = lighting
            image_flat[hit_indices, 2] = lighting
            image = image_flat.reshape(height, width, 3)
    
    return image, hits.reshape(height, width)

def visualize_link_gpu(link, resolution=(800, 600)):
    """
    Create GPU-accelerated visualization with Lambert shading from 4 different views.
    """
    center, extent = get_link_bounds(link)
    distance = extent * 2.5
    
    # Four different viewing angles for comprehensive visualization
    views = [
        {'azimuth': np.pi / 4, 'elevation': np.pi / 6, 'name': 'Front-Right'},
        {'azimuth': 3 * np.pi / 4, 'elevation': np.pi / 6, 'name': 'Back-Right'},
        {'azimuth': -np.pi / 4, 'elevation': np.pi / 6, 'name': 'Front-Left'},
        {'azimuth': 0, 'elevation': np.pi / 3, 'name': 'Top'},
    ]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 14))
    fig.suptitle(f'LinkCSC SDF Surface - Lambert Shading (4 Views)', 
                 fontsize=16, fontweight='bold')
    
    axes = axes.flatten()
    
    total_start = time.time()
    
    for idx, view in enumerate(views):
        print(f"\n[{idx+1}/4] Rendering {view['name']} view...")
        
        image, hits = render_sdf_gpu(link, center, distance, 
                              view['azimuth'], view['elevation'], 
                              resolution, shading='normal', show_progress=True)
        
        ax = axes[idx]
        ax.imshow(image, interpolation='bilinear')
        ax.set_title(f"{view['name']} (az={np.rad2deg(view['azimuth']):.0f}°, el={np.rad2deg(view['elevation']):.0f}°)")
        ax.axis('off')
    
    total_time = time.time() - total_start
    print(f"\n{'='*60}")
    print(f"Total GPU rendering time: {total_time:.3f}s")
    print(f"Average per view: {total_time/4:.3f}s")
    print(f"{'='*60}")
    
    plt.tight_layout()
    return fig

def main():
    """Main test function"""
    print("="*60)
    print("LinkCSC SDF GPU Ray Marching - Lambert Shading Multi-View")
    print("="*60)
    
    if not CUDA_AVAILABLE:
        print("\nERROR: PyCUDA not available!")
        print("Install with: pip install pycuda")
        return
    
    np.random.seed(42)
    
    # Create random link
    print("\nCreating random link...")
    link = create_random_link(r=1.0, min_separation=4.0)
    
    print(f"\nLink properties:")
    print(f"  r = {link.r:.2f}")
    print(f"  theta1 = {np.rad2deg(link.path.theta1):.1f}° (full arc angle)")
    print(f"  theta2 = {np.rad2deg(link.path.theta2):.1f}° (full arc angle)")
    print(f"  straight length = {link.path.tMag:.2f}")
    
    # Debug: Show the link using matplotlib 3D for comparison
    print("\nShowing matplotlib 3D view for comparison...")
    link.show(showPath=True, showFrames=False, showBoundary=True, block=False)
    
    # GPU visualization - Lambert shading from 4 views
    print("\n" + "="*60)
    print("GPU Rendering - Lambert Shading from 4 Views...")
    print("="*60)
    fig = visualize_link_gpu(link, resolution=(800, 600))
    
    plt.show()
    
    print("\n" + "="*60)
    print("Complete!")
    print("="*60)

if __name__ == "__main__":
    main()
