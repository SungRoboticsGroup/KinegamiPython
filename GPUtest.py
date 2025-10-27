import numpy as np
import time
from numba import cuda
import math

# CPU implementation for comparison
def check_collisions_cpu(positions, radii):
    """Check collisions between all pairs of balls on CPU."""
    n = positions.shape[0]
    collisions = []
    
    for i in range(n):
        for j in range(i+1, n):  # Only check each pair once
            # Calculate squared distance between centers
            diff = positions[i] - positions[j]
            dist_squared = np.sum(diff * diff)
            # Check collision (squared distance < squared sum of radii)
            sum_radii = radii[i] + radii[j]
            if dist_squared < (sum_radii * sum_radii):
                collisions.append((i, j))
    
    return collisions

@cuda.jit
def collision_kernel(positions, radii, collisions, collision_count):
    """CUDA kernel: one thread per potential ball pair"""
    # Get thread index
    thread_idx = cuda.grid(1)
    
    # Total number of balls
    n = positions.shape[0]
    
    # Total number of pairs to check
    total_pairs = n * (n - 1) // 2
    
    # Make sure we're within bounds
    if thread_idx < total_pairs:
        # Convert linear thread index to ball pair indices (i,j)
        # This math maps a linear index to upper triangular matrix indices
        i = int(n - 2 - math.floor(math.sqrt(-8*thread_idx + 4*n*(n-1)-7)/2.0 - 0.5))
        j = int(thread_idx + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2)
        
        # Calculate squared distance between centers
        squared_dist = 0.0
        for dim in range(3):  # 3D space
            diff = positions[i, dim] - positions[j, dim]
            squared_dist += diff * diff
        
        # Check collision
        sum_radii = radii[i] + radii[j]
        if squared_dist < (sum_radii * sum_radii):
            # Atomically increment collision counter and get index
            idx = cuda.atomic.add(collision_count, 0, 1)
            if idx < collisions.shape[0]:  # Prevent buffer overflow
                collisions[idx, 0] = i
                collisions[idx, 1] = j

def check_collisions_gpu(positions, radii, max_collisions=100000):
    """Check collisions using CUDA parallelism - one thread per potential pair."""
    n = positions.shape[0]
    
    # Calculate total number of potential pairs
    total_pairs = n * (n - 1) // 2
    
    # Move data to GPU
    positions_gpu = cuda.to_device(positions)
    radii_gpu = cuda.to_device(radii)
    
    # Prepare arrays for results
    collisions_gpu = cuda.device_array((max_collisions, 2), dtype=np.int32)
    collision_count_gpu = cuda.to_device(np.zeros(1, dtype=np.int32))
    
    # Configure CUDA grid
    threads_per_block = 256
    blocks_per_grid = (total_pairs + threads_per_block - 1) // threads_per_block
    
    # Launch kernel
    collision_kernel[blocks_per_grid, threads_per_block](
        positions_gpu, radii_gpu, collisions_gpu, collision_count_gpu
    )
    
    # Get collision count
    collision_count = collision_count_gpu.copy_to_host()[0]
    actual_collisions = min(collision_count, max_collisions)
    
    # Copy collision results back from GPU
    collisions = collisions_gpu[:actual_collisions].copy_to_host()
    
    return collisions

def main():
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Number of balls to simulate
    n_balls = 100
    
    # Generate random balls with positions and radii
    positions = np.random.uniform(-100, 100, (n_balls, 3)).astype(np.float32)
    radii = np.random.uniform(0.5, 5.0, n_balls).astype(np.float32)
    
    # Run CPU version
    print(f"Testing collision detection with {n_balls} balls...")
    print("Running CPU version...")
    start_time = time.time()
    cpu_collisions = check_collisions_cpu(positions, radii)
    cpu_time = time.time() - start_time
    print(f"CPU found {len(cpu_collisions)} collisions in {cpu_time:.4f} seconds")
    
    # Run GPU version
    print("Running GPU version...")
    start_time = time.time()
    gpu_collisions = check_collisions_gpu(positions, radii)
    gpu_time = time.time() - start_time
    print(f"GPU found {len(gpu_collisions)} collisions in {gpu_time:.4f} seconds")
    
    # Calculate speedup
    if cpu_time > 0:
        print(f"GPU speedup: {cpu_time/gpu_time:.2f}x")
    
    # Verify results match
    cpu_set = set((min(i,j), max(i,j)) for i,j in cpu_collisions)
    gpu_set = set((min(i,j), max(i,j)) for i,j in gpu_collisions)
    
    if cpu_set == gpu_set:
        print("Results match! ✓")
    else:
        print("Warning: CPU and GPU results differ!")
        print(f"Only in CPU: {len(cpu_set - gpu_set)}")
        print(f"Only in GPU: {len(gpu_set - cpu_set)}")

if __name__ == "__main__":
    main()
