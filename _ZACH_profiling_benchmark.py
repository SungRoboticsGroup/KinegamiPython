"""
Comprehensive profiling benchmark using cupyx.profiler.benchmark.
Performs sanity checks on system time, CPU time, and GPU time measurements.

Sanity checks:
- If cpu_time + gpu_time > system_time: Measurement is weird (possible timing overlap)
- If cpu_time + gpu_time << system_time: Something else is going on (context switches, I/O)
- If cpu_time + gpu_time ≈ system_time (with small difference): System noise is acceptable
"""

import numpy as np
import yaml
import time
import psutil
import os
from pathlib import Path
from typing import Dict, Tuple, Optional
from dataclasses import dataclass
import matplotlib.pyplot as plt
from spatialmath import SE3, SO3

try:
    from cupyx.profiler import benchmark as cupyx_benchmark
    HAS_CUPYX_PROFILER = True
except ImportError:
    HAS_CUPYX_PROFILER = False
    print("Warning: cupyx.profiler not available")

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

# Import collision analysis
from ZACH_collision_analysis import (
    identify_collisions,
    compute_collision_statistics,
)

# Import link generation
from LinkCSC import LinkCSC


@dataclass
class ProfilingResult:
    """Profiling results with multiple timing measurements"""
    method: str
    num_links: int
    num_points: int
    epsilon: float
    
    # Timing measurements
    system_time: float  # Wall clock time
    cpu_time: float     # CPU time (system + user)
    gpu_time: Optional[float]  # GPU time (if available)
    
    # Sanity check results
    cpu_plus_gpu: Optional[float]  # Sum of CPU and GPU times
    time_diff: Optional[float]  # Difference between system time and cpu+gpu
    is_valid: bool  # Whether measurements make sense
    validity_reason: str  # Explanation of validity


def get_cpu_time_diff() -> float:
    """Get the difference in CPU time since last call"""
    process = psutil.Process(os.getpid())
    times = process.cpu_times()
    return times.user + times.system


def profile_gpu_method(links, epsilon: float, min_radius: float, 
                       collision_config: Dict, num_repeat: int = 3) -> ProfilingResult:
    """Profile GPU method using cupyx.profiler.benchmark"""
    
    try:
        import cupy as cp
        
        # Prepare data
        density = 1.0 / (2.0 * epsilon * min_radius)
        packed = pack_links(links)
        points, point_ids = sample_points_for_links(links, density)
        
        # Define the kernel
        def gpu_kernel():
            pairwise_dist, point_idx = pairwise_link_distances(
                cp, points, point_ids, packed, chunk_points=2048, dtype="float32"
            )
            cp.cuda.Stream.null.synchronize()
            return pairwise_dist, point_idx
        
        # Measure system time
        system_start = time.perf_counter()
        
        # Get CPU time before
        cpu_before = get_cpu_time_diff()
        
        # Profile with cupyx
        if HAS_CUPYX_PROFILER:
            gpu_times = cupyx_benchmark(gpu_kernel, (), n_repeat=num_repeat)
            # cupyx_benchmark returns a _PerfCaseResult object
            # Extract times by converting to array
            if hasattr(gpu_times, 'gpu_times'):
                gpu_time = float(np.mean(gpu_times.gpu_times))
            elif hasattr(gpu_times, 'times'):
                gpu_time = float(np.mean(gpu_times.times))
            else:
                # Try to extract values from the object
                gpu_time = float(gpu_times)
        else:
            # Fallback to manual timing
            times = []
            for _ in range(num_repeat):
                t0 = time.perf_counter()
                gpu_kernel()
                times.append(time.perf_counter() - t0)
            gpu_time = np.mean(times)
        
        # Measure system time end
        system_end = time.perf_counter()
        
        # Get CPU time after
        cpu_after = get_cpu_time_diff()
        
        # Calculate measurements
        system_time = system_end - system_start
        cpu_time = cpu_after - cpu_before
        # gpu_time already calculated above
        
        # Sanity checks
        cpu_plus_gpu = cpu_time + gpu_time
        time_diff = abs(system_time - cpu_plus_gpu)
        
        # Determine validity
        is_valid = True
        validity_reason = "OK"
        
        if cpu_plus_gpu > system_time * 1.1:  # 10% margin
            is_valid = False
            validity_reason = "cpu_time + gpu_time > system_time (timing overlap detected)"
        elif cpu_plus_gpu < system_time * 0.5:
            is_valid = False
            validity_reason = f"cpu_time + gpu_time << system_time ({cpu_plus_gpu:.3f} << {system_time:.3f})"
        elif time_diff < system_time * 0.1:  # Less than 10% difference = system noise
            validity_reason = "System noise (acceptable)"
        
        return ProfilingResult(
            method="GPU (CuPy)",
            num_links=len(links),
            num_points=len(points),
            epsilon=epsilon,
            system_time=system_time,
            cpu_time=cpu_time,
            gpu_time=gpu_time,
            cpu_plus_gpu=cpu_plus_gpu,
            time_diff=time_diff,
            is_valid=is_valid,
            validity_reason=validity_reason
        )
        
    except Exception as e:
        print(f"GPU profiling failed: {e}")
        return None


def profile_cpu_vectorized_method(links, epsilon: float, min_radius: float, 
                                  collision_config: Dict, num_repeat: int = 3) -> ProfilingResult:
    """Profile CPU vectorized method using cupyx.profiler.benchmark"""
    
    try:
        # Prepare data
        density = 1.0 / (2.0 * epsilon * min_radius)
        packed = pack_links(links)
        points, point_ids = sample_points_for_links(links, density)
        
        # Define the kernel
        def cpu_kernel():
            pairwise_dist, point_idx = pairwise_link_distances(
                np, points, point_ids, packed, chunk_points=2048, dtype="float32"
            )
            return pairwise_dist, point_idx
        
        # Measure system time
        system_start = time.perf_counter()
        
        # Get CPU time before
        cpu_before = get_cpu_time_diff()
        
        # Profile with cupyx
        if HAS_CUPYX_PROFILER:
            times = cupyx_benchmark(cpu_kernel, (), n_repeat=num_repeat)
            # cupyx_benchmark returns a _PerfCaseResult object
            # Extract times by checking for different attributes
            if hasattr(times, 'gpu_times'):
                cpu_time_measured = float(np.mean(times.gpu_times))
            elif hasattr(times, 'times'):
                cpu_time_measured = float(np.mean(times.times))
            else:
                # Try to extract values from the object
                cpu_time_measured = float(times)
        else:
            # Fallback to manual timing
            times = []
            for _ in range(num_repeat):
                t0 = time.perf_counter()
                cpu_kernel()
                times.append(time.perf_counter() - t0)
            cpu_time_measured = float(np.mean(times))
        
        # Measure system time end
        system_end = time.perf_counter()
        
        # Get CPU time after
        cpu_after = get_cpu_time_diff()
        
        # Calculate measurements
        system_time = system_end - system_start
        cpu_time = cpu_after - cpu_before  # Actual CPU time consumed
        
        # Sanity checks
        cpu_plus_gpu = cpu_time  # GPU is 0 for CPU method
        time_diff = abs(system_time - cpu_plus_gpu)
        
        # Determine validity
        is_valid = True
        validity_reason = "OK"
        
        if cpu_plus_gpu > system_time * 1.1:  # 10% margin
            is_valid = False
            validity_reason = "cpu_time > system_time (timing overlap detected)"
        elif cpu_plus_gpu < system_time * 0.5:
            is_valid = False
            validity_reason = f"cpu_time << system_time ({cpu_plus_gpu:.3f} << {system_time:.3f})"
        elif time_diff < system_time * 0.1:  # Less than 10% difference = system noise
            validity_reason = "System noise (acceptable)"
        
        return ProfilingResult(
            method="CPU Vectorized (NumPy)",
            num_links=len(links),
            num_points=len(points),
            epsilon=epsilon,
            system_time=system_time,
            cpu_time=cpu_time,
            gpu_time=None,  # CPU only
            cpu_plus_gpu=cpu_plus_gpu,
            time_diff=time_diff,
            is_valid=is_valid,
            validity_reason=validity_reason
        )
        
    except Exception as e:
        print(f"CPU vectorized profiling failed: {e}")
        return None


def profile_cpu_serial_method(links, epsilon: float, min_radius: float, 
                             collision_config: Dict, timeout: int = 10,
                             num_repeat: int = 1) -> ProfilingResult:
    """Profile CPU serial method using system time"""
    
    try:
        # Prepare data
        density = 1.0 / (2.0 * epsilon * min_radius)
        points, point_ids = sample_points_for_links(links, density)
        
        # Define the kernel
        def serial_kernel():
            pairwise_dist, point_idx = pairwise_link_distances_serial(
                links, points, point_ids, timeout_seconds=timeout
            )
            return pairwise_dist, point_idx
        
        # For serial, we run once (not repeated as it's slow)
        system_start = time.perf_counter()
        cpu_before = get_cpu_time_diff()
        
        serial_kernel()
        
        system_end = time.perf_counter()
        cpu_after = get_cpu_time_diff()
        
        # Calculate measurements
        system_time = system_end - system_start
        cpu_time = cpu_after - cpu_before
        
        # Sanity checks
        cpu_plus_gpu = cpu_time  # GPU is 0 for serial CPU method
        time_diff = abs(system_time - cpu_plus_gpu)
        
        # Determine validity
        is_valid = True
        validity_reason = "OK"
        
        if cpu_plus_gpu > system_time * 1.1:  # 10% margin
            is_valid = False
            validity_reason = "cpu_time > system_time (timing overlap detected)"
        elif cpu_plus_gpu < system_time * 0.5:
            is_valid = False
            validity_reason = f"cpu_time << system_time ({cpu_plus_gpu:.3f} << {system_time:.3f})"
        elif time_diff < system_time * 0.1:  # Less than 10% difference = system noise
            validity_reason = "System noise (acceptable)"
        
        return ProfilingResult(
            method="CPU Serial",
            num_links=len(links),
            num_points=len(points),
            epsilon=epsilon,
            system_time=system_time,
            cpu_time=cpu_time,
            gpu_time=None,  # CPU only
            cpu_plus_gpu=cpu_plus_gpu,
            time_diff=time_diff,
            is_valid=is_valid,
            validity_reason=validity_reason
        )
        
    except TimeoutException:
        print("Serial method timed out")
        return None
    except Exception as e:
        print(f"Serial profiling failed: {e}")
        return None


def print_profiling_report(results: list):
    """Print a detailed profiling report"""
    print("\n" + "="*80)
    print("PROFILING RESULTS AND SANITY CHECKS")
    print("="*80)
    
    for result in results:
        if result is None:
            continue
            
        print(f"\n{result.method}")
        print(f"  Links: {result.num_links}, Points: {result.num_points}, Epsilon: {result.epsilon:.2e}")
        print(f"\n  Timing Measurements:")
        print(f"    System Time:      {result.system_time:.6f}s")
        print(f"    CPU Time:         {result.cpu_time:.6f}s")
        if result.gpu_time is not None:
            print(f"    GPU Time:         {result.gpu_time:.6f}s")
            print(f"    CPU + GPU Time:   {result.cpu_plus_gpu:.6f}s")
        else:
            print(f"    CPU + GPU Time:   {result.cpu_plus_gpu:.6f}s (CPU only)")
        
        print(f"\n  Sanity Check:")
        print(f"    Time Difference:  {result.time_diff:.6f}s")
        if result.is_valid:
            print(f"    Status:           ✓ VALID ({result.validity_reason})")
        else:
            print(f"    Status:           ✗ INVALID ({result.validity_reason})")
        
        # Additional analysis
        if result.gpu_time is not None:
            gpu_ratio = (result.gpu_time / result.system_time) * 100
            print(f"    GPU Utilization:  {gpu_ratio:.1f}% of system time")


def generate_random_links(num_links: int, config: Dict, seed: int):
    """Generate random LinkCSC objects"""
    rng = np.random.RandomState(seed)
    links = []
    
    r_min, r_max = config['radius_range']
    p_min, p_max = config['position_range']
    max_angle = config['max_angle_per_elbow']
    
    for i in range(num_links):
        r = rng.uniform(r_min, r_max)
        
        # Random poses
        def random_pose():
            pos = rng.uniform(p_min, p_max, 3)
            axis = rng.randn(3)
            axis = axis / np.linalg.norm(axis)
            angle = rng.uniform(0, 2*np.pi)
            R = SO3.AngleAxis(angle, axis)
            return SE3.Rt(R, pos)
        
        start_pose = random_pose()
        end_pose = random_pose()
        
        # Ensure poses are not too close
        max_attempts = 10
        attempts = 0
        while np.linalg.norm(end_pose.t - start_pose.t) < 2 * r and attempts < max_attempts:
            end_pose = random_pose()
            attempts += 1
        
        if attempts >= max_attempts:
            continue
        
        try:
            link = LinkCSC(
                r=r,
                StartDubinsPose=start_pose,
                EndDubinsPose=end_pose,
                maxAnglePerElbow=max_angle,
                EPSILON=r_min * 0.01
            )
            links.append(link)
        except (ValueError, AssertionError):
            continue
    
    return links, r_min


def main():
    """Main profiling benchmark"""
    # Load configuration
    config_path = Path('ZACH_benchmark_config.yaml')
    if not config_path.exists():
        print(f"Error: Configuration file {config_path} not found")
        return
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    link_config = config['link_generation']
    collision_config = config['collision_analysis']
    
    print("="*80)
    print("CUPYX PROFILING BENCHMARK WITH SANITY CHECKS")
    print("="*80)
    
    # Test parameters
    test_cases = [
        {"name": "Small (5 links, 1e-2 epsilon)", "num_links": 5, "epsilon": 1e-2, "seed": 42},
        {"name": "Medium (10 links, 1e-3 epsilon)", "num_links": 10, "epsilon": 1e-3, "seed": 67},
    ]
    
    all_results = []
    
    for test_case in test_cases:
        print(f"\n\nTest: {test_case['name']}")
        print("-" * 80)
        
        # Generate links
        links, min_radius = generate_random_links(
            test_case['num_links'], link_config, test_case['seed']
        )
        
        if len(links) == 0:
            print("Failed to generate links")
            continue
        
        epsilon = float(test_case['epsilon']) * min_radius
        
        # Profile each method
        results = []
        
        print("\nProfiling GPU method...")
        gpu_result = profile_gpu_method(links, test_case['epsilon'], min_radius, collision_config)
        if gpu_result:
            results.append(gpu_result)
        
        print("Profiling CPU vectorized method...")
        cpu_vec_result = profile_cpu_vectorized_method(links, test_case['epsilon'], min_radius, collision_config)
        if cpu_vec_result:
            results.append(cpu_vec_result)
        
        print("Profiling CPU serial method...")
        cpu_serial_result = profile_cpu_serial_method(links, test_case['epsilon'], min_radius, collision_config)
        if cpu_serial_result:
            results.append(cpu_serial_result)
        
        # Print report for this test
        print_profiling_report(results)
        all_results.extend(results)
    
    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Total measurements: {len(all_results)}")
    valid_count = sum(1 for r in all_results if r.is_valid)
    print(f"Valid measurements: {valid_count}/{len(all_results)}")
    
    if valid_count == len(all_results):
        print("\n✓ All measurements are valid!")
    else:
        print("\n✗ Some measurements show issues:")
        for result in all_results:
            if not result.is_valid:
                print(f"  - {result.method}: {result.validity_reason}")


if __name__ == "__main__":
    main()
