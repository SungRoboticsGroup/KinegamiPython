"""
Simple profiling benchmark comparing cupyx.profiler.benchmark vs system time.
Tests both timing methods on the same functions to see how they compare.
"""

import numpy as np
import yaml
import time
from pathlib import Path
from typing import Dict, Tuple
import matplotlib.pyplot as plt
from spatialmath import SE3, SO3

try:
    import cupy as cp
    from cupyx.profiler import benchmark as cupyx_benchmark
    HAS_CUPY = True
    HAS_CUPYX_PROFILER = True
except ImportError:
    HAS_CUPY = False
    HAS_CUPYX_PROFILER = False
    print("Warning: CuPy/cupyx.profiler not available")

# Import the vectorized methods
from ZACH_vectorized_link_sdf import (
    pack_links, 
    sample_points_for_links, 
    pairwise_link_distances
)

# Import the serial method
from ZACH_cpu_serial_link_sdf import (
    pairwise_link_distances_serial,
    TimeoutException
)

# Import link generation
from LinkCSC import LinkCSC


def load_config():
    """Load benchmark configuration"""
    config_path = Path(__file__).parent / "ZACH_benchmark_config.yaml"
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def generate_test_links(num_links: int) -> list:
    """Generate random test links"""
    links = []
    for i in range(num_links):
        # Random start and end poses
        start_pos = np.random.uniform(-1, 1, 3)
        end_pos = start_pos + np.random.uniform(-0.5, 0.5, 3)
        
        start_rot = SO3.Rand()
        end_rot = SO3.Rand()
        
        start_pose = SE3(start_pos) * SE3(start_rot)
        end_pose = SE3(end_pos) * SE3(end_rot)
        
        # Random radius
        radius = np.random.uniform(0.01, 0.05)
        
        try:
            link = LinkCSC(
                r=radius,
                StartDubinsPose=start_pose,
                EndDubinsPose=end_pose,
                maxAnglePerElbow=np.pi/3,
                EPSILON=radius * 0.01
            )
            links.append(link)
        except (ValueError, AssertionError):
            # Skip invalid links
            continue
    
    return links


def test_timing_methods():
    """Compare cupyx benchmark vs system time on all three computation methods"""
    
    config = load_config()
    
    # Test parameters
    num_links = 8
    epsilon = 1e-3
    min_radius = 0.02
    timeout = 5
    
    print("="*80)
    print("SIMPLE TIMING COMPARISON: CUPYX vs SYSTEM TIME")
    print("="*80)
    print(f"Test setup: {num_links} links, epsilon={epsilon}")
    print()
    
    # Generate test data
    links = generate_test_links(num_links)
    density = 1.0 / (2.0 * epsilon * min_radius)
    points, point_ids = sample_points_for_links(links, density)
    packed = pack_links(links)
    
    print(f"Generated {len(points)} points for testing")
    print()
    
    # Test 1: GPU Method
    if HAS_CUPY:
        print("1. GPU METHOD (CuPy)")
        print("-" * 40)
        
        def gpu_kernel():
            cp_points = cp.array(points, dtype=cp.float32)
            cp_point_ids = cp.array(point_ids, dtype=cp.int32)
            pairwise_dist, point_idx = pairwise_link_distances(
                cp, cp_points, cp_point_ids, packed, chunk_points=2048, dtype="float32"
            )
            cp.cuda.Stream.null.synchronize()
            return pairwise_dist, point_idx
        
        # CuPy benchmark
        if HAS_CUPYX_PROFILER:
            cupyx_result = cupyx_benchmark(gpu_kernel, (), n_repeat=3)
            print(f"  CuPy benchmark result: {cupyx_result}")
            if hasattr(cupyx_result, 'gpu_times'):
                print(f"  GPU times: {cupyx_result.gpu_times}")
            if hasattr(cupyx_result, 'cpu_times'):
                print(f"  CPU times: {cupyx_result.cpu_times}")
            if hasattr(cupyx_result, 'times'):
                print(f"  Total times: {cupyx_result.times}")
        
        # System time
        system_times = []
        for i in range(3):
            start = time.perf_counter()
            gpu_kernel()
            end = time.perf_counter()
            system_times.append(end - start)
        
        print(f"  System times: {system_times}")
        print(f"  System time mean: {np.mean(system_times):.4f}s")
        print()
    
    # Test 2: CPU Vectorized Method
    print("2. CPU VECTORIZED METHOD (NumPy)")
    print("-" * 40)
    
    def cpu_kernel():
        pairwise_dist, point_idx = pairwise_link_distances(
            np, points, point_ids, packed, chunk_points=2048, dtype="float32"
        )
        return pairwise_dist, point_idx
    
    # CuPy benchmark (works on CPU functions too)
    if HAS_CUPYX_PROFILER:
        cupyx_result = cupyx_benchmark(cpu_kernel, (), n_repeat=3)
        print(f"  CuPy benchmark result: {cupyx_result}")
        if hasattr(cupyx_result, 'gpu_times'):
            print(f"  GPU times: {cupyx_result.gpu_times}")
        if hasattr(cupyx_result, 'cpu_times'):
            print(f"  CPU times: {cupyx_result.cpu_times}")
        if hasattr(cupyx_result, 'times'):
            print(f"  Total times: {cupyx_result.times}")
    
    # System time
    system_times = []
    for i in range(3):
        start = time.perf_counter()
        cpu_kernel()
        end = time.perf_counter()
        system_times.append(end - start)
    
    print(f"  System times: {system_times}")
    print(f"  System time mean: {np.mean(system_times):.4f}s")
    print()
    
    # Test 3: Serial Method
    print("3. SERIAL METHOD")
    print("-" * 40)
    
    def serial_kernel():
        try:
            pairwise_dist, point_idx = pairwise_link_distances_serial(
                links, points, point_ids, timeout_seconds=timeout
            )
            return pairwise_dist, point_idx
        except TimeoutException:
            print("  Serial method timed out")
            return None, None
    
    # CuPy benchmark
    if HAS_CUPYX_PROFILER:
        try:
            cupyx_result = cupyx_benchmark(serial_kernel, (), n_repeat=1)  # Only 1 repeat for serial
            print(f"  CuPy benchmark result: {cupyx_result}")
            if hasattr(cupyx_result, 'gpu_times'):
                print(f"  GPU times: {cupyx_result.gpu_times}")
            if hasattr(cupyx_result, 'cpu_times'):
                print(f"  CPU times: {cupyx_result.cpu_times}")
            if hasattr(cupyx_result, 'times'):
                print(f"  Total times: {cupyx_result.times}")
        except Exception as e:
            print(f"  CuPy benchmark failed: {e}")
    
    # System time
    start = time.perf_counter()
    result = serial_kernel()
    end = time.perf_counter()
    
    if result[0] is not None:
        print(f"  System time: {end - start:.4f}s")
    print()


if __name__ == "__main__":
    test_timing_methods()
