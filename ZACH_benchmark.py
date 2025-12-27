"""
Benchmark suite for comparing three methods of pairwise link distance computation:
1. GPU vectorized (CuPy)
2. CPU vectorized (NumPy)
3. CPU serial (for loop with AABB pruning)
"""

import numpy as np
import time
import yaml
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
import matplotlib.pyplot as plt
from spatialmath import SE3, SO3
from collections import defaultdict

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
    compare_collision_detections,
    filter_intersecting_pairs,
    print_collision_report
)

# Import link generation
from LinkCSC import LinkCSC


@dataclass
class BenchmarkResult:
    """Results from a single benchmark run"""
    name: str
    method: str
    num_links: int
    num_points: int
    epsilon: float  # Changed from point_density to epsilon
    execution_time: float
    pairwise_distances: np.ndarray
    pairwise_point_indices: np.ndarray
    collisions: List
    collision_stats: object
    points: np.ndarray = None  # Points array for witness verification
    success: bool = True
    error_message: str = ""
    timed_out: bool = False
    chunk_size: Optional[int] = None  # For chunk size scaling tests


def calculate_point_density(epsilon: float, min_radius: float) -> float:
    """Calculate point density from epsilon and minimum radius"""
    return 1.0 / (2.0 * epsilon * min_radius)


def generate_random_SE3(position_range: Tuple[float, float], seed: Optional[int] = None) -> SE3:
    """Generate a random SE3 pose"""
    if seed is not None:
        rng = np.random.RandomState(seed)
    else:
        rng = np.random
    
    # Random position
    pos = rng.uniform(position_range[0], position_range[1], 3)
    
    # Random orientation (using random rotation)
    axis = rng.randn(3)
    axis = axis / np.linalg.norm(axis)
    angle = rng.uniform(0, 2*np.pi)
    R = SO3.AngleAxis(angle, axis)
    
    return SE3.Rt(R, pos)


def generate_random_links(num_links: int, config: Dict, seed: int) -> Tuple[List[LinkCSC], float]:
    """Generate random LinkCSC objects with total length bounded by max_link_length.
    
    Returns the list of links and the minimum radius from config (for epsilon calculations).
    """
    rng = np.random.RandomState(seed)
    links = []
    
    r_min, r_max = config['radius_range']
    p_min, p_max = config['position_range']
    max_angle = config['max_angle_per_elbow']
    max_link_length = config.get('max_link_length', 10.0)
    
    link_attempt_count = 0
    max_total_attempts = num_links * 100  # Prevent infinite loops
    
    while len(links) < num_links and link_attempt_count < max_total_attempts:
        link_attempt_count += 1
        
        # Random radius
        r = rng.uniform(r_min, r_max)
        
        # Random start and end poses
        start_seed = seed + len(links) * 2 + link_attempt_count
        end_seed = seed + len(links) * 2 + 1 + link_attempt_count
        start_pose = generate_random_SE3((p_min, p_max), start_seed)
        end_pose = generate_random_SE3((p_min, p_max), end_seed)
        
        # Compute the distance between start and end poses
        link_length = np.linalg.norm(end_pose.t - start_pose.t)
        
        if link_length > max_link_length:
            continue
        
        # Ensure poses are not too close
        if link_length < 2 * r:
            continue
        
        try:
            link = LinkCSC(
                r=r,
                StartDubinsPose=start_pose,
                EndDubinsPose=end_pose,
                maxAnglePerElbow=max_angle,
                EPSILON=r_min * 0.01  # Use a small epsilon relative to minimum radius for link generation
            )
            links.append(link)
        except (ValueError, AssertionError) as e:
            continue
    
    return links, r_min


def run_gpu_method(links: List[LinkCSC], epsilon: float, min_radius: float, 
                  collision_config: Dict, chunk_size: Optional[int] = None) -> BenchmarkResult:
    """Run GPU vectorized method"""
    density = calculate_point_density(epsilon, min_radius)
    try:
        import cupy as cp
        xp = cp
        
        # Pack links
        packed = pack_links(links)
        
        # Sample points
        points, point_ids = sample_points_for_links(links, density)
        
        # Time the pairwise computation
        start_time = time.perf_counter()
        pairwise_dist, point_idx = pairwise_link_distances(
            xp, points, point_ids, packed, chunk_points=2048, dtype="float32"
        )
        
        # Convert to CPU
        pairwise_dist_cpu = cp.asnumpy(pairwise_dist)
        point_idx_cpu = cp.asnumpy(point_idx)

        cp.cuda.Stream.null.synchronize()
        
        execution_time = time.perf_counter() - start_time
        
        # Collision analysis
        collisions = identify_collisions(
            pairwise_dist_cpu, point_idx_cpu, points, links,
            collision_config['intersection_threshold_multiplier']
        )
        L = len(links)
        total_pairs = L * (L - 1) // 2
        collision_stats = compute_collision_statistics(collisions, total_pairs)
        
        return BenchmarkResult(
            name="gpu_test",
            method="GPU (CuPy)",
            num_links=len(links),
            num_points=len(points),
            epsilon=epsilon,
            execution_time=execution_time,
            pairwise_distances=pairwise_dist_cpu,
            pairwise_point_indices=point_idx_cpu,
            collisions=collisions,
            collision_stats=collision_stats,
            points=points,
            chunk_size=chunk_size,
            success=True
        )
    except Exception as e:
        return BenchmarkResult(
            name="gpu_test",
            method="GPU (CuPy)",
            num_links=len(links),
            num_points=0,
            epsilon=epsilon,
            execution_time=0.0,
            pairwise_distances=np.array([]),
            pairwise_point_indices=np.array([]),
            collisions=[],
            collision_stats=None,
            success=False,
            error_message=str(e),
            chunk_size=chunk_size
        )


def run_cpu_vectorized_method(links: List[LinkCSC], epsilon: float, min_radius: float, 
                             collision_config: Dict, chunk_size: Optional[int] = None) -> BenchmarkResult:
    """Run CPU vectorized method"""
    density = calculate_point_density(epsilon, min_radius)
    try:
        xp = np
        
        # Pack links
        packed = pack_links(links)
        
        # Sample points
        points, point_ids = sample_points_for_links(links, density)
        
        # Time the pairwise computation
        start_time = time.perf_counter()
        pairwise_dist, point_idx = pairwise_link_distances(
            xp, points, point_ids, packed, chunk_points=2048, dtype="float32"
        )
        execution_time = time.perf_counter() - start_time
        
        # Collision analysis
        collisions = identify_collisions(
            pairwise_dist, point_idx, points, links,
            collision_config['intersection_threshold_multiplier']
        )
        L = len(links)
        total_pairs = L * (L - 1) // 2
        collision_stats = compute_collision_statistics(collisions, total_pairs)
        
        return BenchmarkResult(
            name="cpu_vectorized_test",
            method="CPU Vectorized (NumPy)",
            num_links=len(links),
            num_points=len(points),
            epsilon=epsilon,
            execution_time=execution_time,
            pairwise_distances=pairwise_dist,
            pairwise_point_indices=point_idx,
            collisions=collisions,
            collision_stats=collision_stats,
            points=points,
            chunk_size=chunk_size,
            success=True
        )
    except Exception as e:
        return BenchmarkResult(
            name="cpu_vectorized_test",
            method="CPU Vectorized (NumPy)",
            num_links=len(links),
            num_points=0,
            epsilon=epsilon,
            execution_time=0.0,
            pairwise_distances=np.array([]),
            pairwise_point_indices=np.array([]),
            collisions=[],
            collision_stats=None,
            success=False,
            error_message=str(e),
            chunk_size=chunk_size
        )


def run_cpu_serial_method(links: List[LinkCSC], epsilon: float, min_radius: float,
                          collision_config: Dict, timeout: int, 
                          skip_if_previous_timeout: bool = False) -> BenchmarkResult:
    """Run CPU serial method with timeout tracking"""
    density = calculate_point_density(epsilon, min_radius)
    
    # Skip if previous timeout occurred
    if skip_if_previous_timeout:
        return BenchmarkResult(
            name="cpu_serial_test",
            method="CPU Serial",
            num_links=len(links),
            num_points=0,
            epsilon=epsilon,
            execution_time=0.0,
            pairwise_distances=np.array([]),
            pairwise_point_indices=np.array([]),
            collisions=[],
            collision_stats=None,
            success=False,
            error_message="Skipped due to previous timeout",
            timed_out=True
        )
    try:
        # Sample points (using same method as vectorized for consistency)
        points, point_ids = sample_points_for_links(links, density)
        
        # Time the pairwise computation with timeout
        start_time = time.perf_counter()
        pairwise_dist, point_idx = pairwise_link_distances_serial(
            links, points, point_ids, timeout_seconds=timeout
        )
        execution_time = time.perf_counter() - start_time
        
        # Collision analysis
        collisions = identify_collisions(
            pairwise_dist, point_idx, points, links,
            collision_config['intersection_threshold_multiplier']
        )
        L = len(links)
        total_pairs = L * (L - 1) // 2
        collision_stats = compute_collision_statistics(collisions, total_pairs)
        
        return BenchmarkResult(
            name="cpu_serial_test",
            method="CPU Serial",
            num_links=len(links),
            num_points=len(points),
            epsilon=epsilon,
            execution_time=execution_time,
            pairwise_distances=pairwise_dist,
            pairwise_point_indices=point_idx,
            collisions=collisions,
            collision_stats=collision_stats,
            points=points,
            success=True
        )
    except TimeoutException:
        return BenchmarkResult(
            name="cpu_serial_test",
            method="CPU Serial",
            num_links=len(links),
            num_points=len(points) if 'points' in locals() else 0,
            epsilon=epsilon,
            execution_time=timeout,
            pairwise_distances=np.array([]),
            pairwise_point_indices=np.array([]),
            collisions=[],
            collision_stats=None,
            success=False,
            error_message=f"Timed out after {timeout}s",
            timed_out=True
        )
    except Exception as e:
        return BenchmarkResult(
            name="cpu_serial_test",
            method="CPU Serial",
            num_links=len(links),
            num_points=0,
            epsilon=epsilon,
            execution_time=0.0,
            pairwise_distances=np.array([]),
            pairwise_point_indices=np.array([]),
            collisions=[],
            collision_stats=None,
            success=False,
            error_message=str(e)
        )


def verify_results(results: List[BenchmarkResult], config: Dict, min_radius: float) -> Tuple[bool, str]:
    """Verify that all methods produce similar results.
    
    Uses rtol/atol for comparing method results (should be very close).
    Uses epsilon*min_radius for collision detection tolerance only.
    """
    if len(results) < 2:
        return True, "Not enough results to compare"
    
    successful_results = [r for r in results if r.success]
    if len(successful_results) < 2:
        return False, "Not enough successful results to compare"
    
    link_config = config.get('link_generation', config)
    rtol = float(link_config['distance_rtol'])
    atol = float(link_config['distance_atol'])
    allow_point_mismatch = link_config.get('allow_point_mismatch', False)
    
    # Epsilon tolerance for collision detection only
    epsilon_tol = successful_results[0].epsilon * min_radius
    
    base_result = successful_results[0]
    
    # Print minimum distances for all methods
    print(f"  Minimum Distances:")
    for r in successful_results:
        valid_mask = ~(np.isinf(r.pairwise_distances))
        if np.any(valid_mask):
            min_dist = np.min(r.pairwise_distances[valid_mask])
            print(f"    {r.method}: min_dist={min_dist:.6f}")
        else:
            print(f"    {r.method}: All distances are Inf")
    
    for result in successful_results[1:]:
        # Check distances match (ignore NaN/Inf locations)
        # Use rtol/atol for comparing method results (should be close)
        base_valid = ~(np.isinf(base_result.pairwise_distances))
        result_valid = ~(np.isinf(result.pairwise_distances))
        
        # Compare only valid entries
        common_valid = base_valid & result_valid
        if np.any(common_valid):
            base_valid_vals = base_result.pairwise_distances[common_valid]
            result_valid_vals = result.pairwise_distances[common_valid]
            if not np.allclose(base_valid_vals, result_valid_vals, rtol=rtol, atol=atol, equal_nan=True):
                max_diff = np.max(np.abs(base_valid_vals - result_valid_vals))
                return False, f"Distance mismatch between {base_result.method} and {result.method}: max diff = {max_diff} (rtol={rtol}, atol={atol})"
        else:
            print(f"    No common valid entries to compare!")
        
        # Check point indices match (if required)
        if not allow_point_mismatch:
            if not np.array_equal(base_result.pairwise_point_indices, 
                                 result.pairwise_point_indices):
                return False, f"Point index mismatch between {base_result.method} and {result.method}"
        
        # Check collisions match using epsilon tolerance (collision-specific)
        match, msg = compare_collision_detections(
            base_result.collisions, result.collisions, epsilon_tol, epsilon_tol
        )
        if not match:
            return False, f"Collision mismatch between {base_result.method} and {result.method}: {msg}"
    
    return True, "All methods agree"


def run_single_test(test_config: Dict, link_config: Dict, 
                   collision_config: Dict, timeout: int) -> Tuple[List[BenchmarkResult], bool, str]:
    """Run all three methods on a single test configuration"""
    print(f"\nRunning test: {test_config['name']}")
    print(f"  Links: {test_config['num_links']}, Epsilon: {test_config['epsilon']}")
    
    # Generate links
    links, min_radius = generate_random_links(
        test_config['num_links'],
        link_config,
        test_config['seed']
    )
    
    if len(links) == 0:
        return [], False, "Failed to generate any links"
    
    print(f"  Generated {len(links)} links")
    
    results = []
    serial_timeout_occurred = False
    results = []
    
    # Get default chunk size from config
    chunk_size = link_config.get('chunk_size', 1024)
    
    # GPU method
    print("  Running GPU method...")
    gpu_result = run_gpu_method(links, float(test_config['epsilon']), min_radius, collision_config, chunk_size=chunk_size)
    results.append(gpu_result)
    if gpu_result.success:
        print(f"    Time: {gpu_result.execution_time:.4f}s, Points: {gpu_result.num_points}")
    else:
        print(f"    Failed: {gpu_result.error_message}")
    
    # CPU vectorized method
    print("  Running CPU vectorized method...")
    cpu_vec_result = run_cpu_vectorized_method(links, float(test_config['epsilon']), min_radius, collision_config, chunk_size=chunk_size)
    results.append(cpu_vec_result)
    if cpu_vec_result.success:
        print(f"    Time: {cpu_vec_result.execution_time:.4f}s, Points: {cpu_vec_result.num_points}")
    else:
        print(f"    Failed: {cpu_vec_result.error_message}")
    
    # CPU serial method
    print("  Running CPU serial method...")
    cpu_serial_result = run_cpu_serial_method(links, float(test_config['epsilon']), min_radius,
                                              collision_config, timeout, 
                                              skip_if_previous_timeout=serial_timeout_occurred)
    results.append(cpu_serial_result)
    if cpu_serial_result.success:
        print(f"    Time: {cpu_serial_result.execution_time:.4f}s, Points: {cpu_serial_result.num_points}")
    elif cpu_serial_result.timed_out:
        print(f"    Timed out after {timeout}s")
        serial_timeout_occurred = True
    else:
        print(f"    Failed: {cpu_serial_result.error_message}")
    
    return results, *verify_results(results, link_config, min_radius)


def run_scaling_tests(config: Dict) -> Dict[str, List[BenchmarkResult]]:
    """Run comprehensive scaling tests with epsilon, link count, and chunk size variations"""
    if not config['scaling_tests'].get('enabled', True):
        print("\n Scaling tests disabled")
        return {}
    
    scaling_config = config['scaling_tests']
    link_config = config['link_generation']
    collision_config = config['collision_analysis']
    timeout = config['timeouts']['serial_method_timeout']
    
    # Get default chunk size from config for non-chunk-scaling tests
    default_chunk_size = link_config.get('chunk_size', 1024)
    
    results = {
        'epsilon_scaling': [],
        'link_scaling': [],
        'chunk_size_scaling': []
    }
    
    if scaling_config.get('epsilon_scaling', {}).get('enabled', False):
        # Epsilon scaling tests (3 graphs with different link counts)
        print("\n" + "="*60)
        print("EPSILON SCALING TESTS")
        print("="*60)
        
        epsilon_config = scaling_config['epsilon_scaling']
        for num_links in epsilon_config['num_links']:
            print(f"\nEpsilon scaling with {num_links} links:")

            serial_timeout_occurred = False
            
            links, min_radius = generate_random_links(num_links, link_config, epsilon_config['seed'])
            if len(links) == 0:
                continue
            
            for epsilon_base in epsilon_config['epsilons']:
                epsilon_base_float = float(epsilon_base)
                epsilon = epsilon_base_float * min_radius  # Multiply by minimum radius
                print(f"  Epsilon: {epsilon:.2e} (base={epsilon_base_float:.2e})")
                
                # GPU method
                gpu_result = run_gpu_method(links, epsilon, min_radius, collision_config, chunk_size=default_chunk_size)
                if gpu_result.success:
                    gpu_result.name = f"epsilon_scaling_{num_links}links"
                    results['epsilon_scaling'].append(gpu_result)
                
                # CPU vectorized method
                cpu_vec_result = run_cpu_vectorized_method(links, epsilon, min_radius, collision_config, chunk_size=default_chunk_size)
                if cpu_vec_result.success:
                    cpu_vec_result.name = f"epsilon_scaling_{num_links}links"
                    results['epsilon_scaling'].append(cpu_vec_result)
                
                # CPU serial method (with timeout tracking)
                cpu_serial_result = run_cpu_serial_method(links, epsilon, min_radius, collision_config, 
                                                        timeout, skip_if_previous_timeout=serial_timeout_occurred)
                if cpu_serial_result.success:
                    cpu_serial_result.name = f"epsilon_scaling_{num_links}links"
                    results['epsilon_scaling'].append(cpu_serial_result)
                elif cpu_serial_result.timed_out and not serial_timeout_occurred:
                    serial_timeout_occurred = True
                    print(f"    WARNING: Serial method timed out after {timeout}s")

    if scaling_config.get('link_scaling', {}).get('enabled', False):
        # Link scaling tests (3 graphs with different epsilons)
        print("\n" + "="*60)
        print("LINK SCALING TESTS")
        print("="*60)
        
        link_scaling_config = scaling_config['link_scaling']
        for epsilon_base in link_scaling_config['epsilons']:
            serial_timeout_occurred = False
            epsilon_base_float = float(epsilon_base)
            print(f"\nLink scaling with epsilon={epsilon_base_float:.2e}:")
            
            for num_links in link_scaling_config['num_links_list']:
                print(f"  Num links: {num_links}")
                
                links, min_radius = generate_random_links(num_links, link_config, link_scaling_config['seed'])
                if len(links) == 0:
                    continue
                    
                epsilon = epsilon_base_float * min_radius
                
                # GPU method
                gpu_result = run_gpu_method(links, epsilon, min_radius, collision_config, chunk_size=default_chunk_size)
                if gpu_result.success:
                    gpu_result.name = f"link_scaling_eps{epsilon_base_float:.0e}"
                    results['link_scaling'].append(gpu_result)
                
                # CPU vectorized method
                cpu_vec_result = run_cpu_vectorized_method(links, epsilon, min_radius, collision_config, chunk_size=default_chunk_size)
                if cpu_vec_result.success:
                    cpu_vec_result.name = f"link_scaling_eps{epsilon_base_float:.0e}"
                    results['link_scaling'].append(cpu_vec_result)
                
                # CPU serial method (with timeout tracking)
                cpu_serial_result = run_cpu_serial_method(links, epsilon, min_radius, collision_config, 
                                                        timeout, skip_if_previous_timeout=serial_timeout_occurred)
                if cpu_serial_result.success:
                    cpu_serial_result.name = f"link_scaling_eps{epsilon_base_float:.0e}"
                    results['link_scaling'].append(cpu_serial_result)
                elif cpu_serial_result.timed_out and not serial_timeout_occurred:
                    serial_timeout_occurred = True
                    print(f"    WARNING: Serial method timed out after {timeout}s")
    
    # Chunk size scaling tests
    if scaling_config.get('chunk_size_scaling', {}).get('enabled', False):
        print("\n" + "="*60)
        print("CHUNK SIZE SCALING TESTS")
        print("="*60)
        
        chunk_config = scaling_config['chunk_size_scaling']
        
        # 2 graphs with different epsilons
        for epsilon_base in chunk_config['epsilons']:
            serial_timeout_occurred = False
            epsilon_base_float = float(epsilon_base)
            print(f"\nChunk scaling with epsilon={epsilon_base_float:.2e}:")
            
            links, min_radius = generate_random_links(10, link_config, chunk_config['seed'])  # Fixed link count for chunk testing
            if len(links) == 0:
                continue
                
            epsilon = epsilon_base_float * min_radius
            
            for chunk_size in chunk_config['chunk_sizes']:
                print(f"  Chunk size: {chunk_size}")
                
                # GPU method with chunk size
                gpu_result = run_gpu_method(links, epsilon, min_radius, collision_config, chunk_size=chunk_size)
                if gpu_result.success:
                    gpu_result.name = f"chunk_scaling_eps{epsilon_base_float:.0e}"
                    results['chunk_size_scaling'].append(gpu_result)
                
                # CPU vectorized method with chunk size
                cpu_vec_result = run_cpu_vectorized_method(links, epsilon, min_radius, collision_config, chunk_size=chunk_size)
                if cpu_vec_result.success:
                    cpu_vec_result.name = f"chunk_scaling_eps{epsilon_base_float:.0e}"
                    results['chunk_size_scaling'].append(cpu_vec_result)
        
        # 2 graphs with different link counts
        for num_links in chunk_config['num_links']:
            serial_timeout_occurred = False
            print(f"\nChunk scaling with {num_links} links:")
            
            links, min_radius = generate_random_links(num_links, link_config, chunk_config['seed'])
            if len(links) == 0:
                continue
            
            epsilon = 1e-3 * min_radius  # Fixed epsilon for link count testing
            
            for chunk_size in chunk_config['chunk_sizes']:
                print(f"  Chunk size: {chunk_size}")
                
                # GPU method with chunk size
                gpu_result = run_gpu_method(links, epsilon, min_radius, collision_config, chunk_size=chunk_size)
                if gpu_result.success:
                    gpu_result.name = f"chunk_scaling_{num_links}links"
                    results['chunk_size_scaling'].append(gpu_result)
                
                # CPU vectorized method with chunk size
                cpu_vec_result = run_cpu_vectorized_method(links, epsilon, min_radius, collision_config, chunk_size=chunk_size)
                if cpu_vec_result.success:
                    cpu_vec_result.name = f"chunk_scaling_{num_links}links"
                    results['chunk_size_scaling'].append(cpu_vec_result)
    
    return results


def plot_results(scaling_results: Dict[str, List[BenchmarkResult]], output_dir: Path, config: Dict):
    """Generate comprehensive multi-subplot plots from scaling test results"""
    from collections import defaultdict
    output_dir.mkdir(exist_ok=True)

    if not scaling_results:
        print("\nNo scaling results to plot")
        return

    scaling_config = config.get('scaling_tests', {})

    # Plot epsilon scaling: 3 subplots (one for each link count)
    if scaling_results.get('epsilon_scaling') and scaling_config.get('epsilon_scaling', {}).get('enabled', True):
        print("\nGenerating epsilon_scaling plot...")
        epsilon_config = scaling_config.get('epsilon_scaling', {})
        num_links_list = epsilon_config.get('num_links', [5, 10, 20])
        
        fig, axes = plt.subplots(1, len(num_links_list), figsize=(6*len(num_links_list), 5))
        if len(num_links_list) == 1:
            axes = [axes]
        
        for idx, num_links in enumerate(num_links_list):
            ax = axes[idx]
            methods_data = defaultdict(lambda: {'epsilons': [], 'times': []})
            
            # Filter results for this num_links
            for result in scaling_results['epsilon_scaling']:
                if result.num_links == num_links:
                    methods_data[result.method]['epsilons'].append(result.epsilon)
                    methods_data[result.method]['times'].append(result.execution_time)
            
            # Plot each method
            for method in sorted(methods_data.keys()):
                data = methods_data[method]
                if data['epsilons']:
                    sorted_data = sorted(zip(data['epsilons'], data['times']))
                    epsilons, times = zip(*sorted_data)
                    ax.plot(epsilons, times, marker='o', label=method, linewidth=2, markersize=6)
            
            ax.set_xlabel('Epsilon', fontsize=11)
            ax.set_ylabel('Total Execution Time (s)', fontsize=11)
            ax.set_title(f'Epsilon Scaling (N_links={num_links})', fontsize=12, fontweight='bold')
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)
            ax.set_yscale('log')
            ax.set_xscale('log')
        
        plt.tight_layout()
        filename = output_dir / 'epsilon_scaling.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved plot: {filename}")

    # Plot link scaling: subplots for each epsilon
    if scaling_results.get('link_scaling') and scaling_config.get('link_scaling', {}).get('enabled', True):
        print("\nGenerating link_scaling plot...")
        link_config = scaling_config.get('link_scaling', {})
        epsilons_list = link_config.get('epsilons', [1e-2])
        
        fig, axes = plt.subplots(1, len(epsilons_list), figsize=(6*len(epsilons_list), 5))
        if len(epsilons_list) == 1:
            axes = [axes]
        
        for idx, epsilon_base in enumerate(epsilons_list):
            ax = axes[idx]
            methods_data = defaultdict(lambda: {'num_links': [], 'times': []})
            epsilon_base_float = float(epsilon_base)
            
            # Collect all data for this epsilon
            for result in scaling_results['link_scaling']:
                if np.isclose(result.epsilon, epsilon_base_float):
                    methods_data[result.method]['num_links'].append(result.num_links)
                    methods_data[result.method]['times'].append(result.execution_time)
            
            # Plot each method
            for method in sorted(methods_data.keys()):
                data = methods_data[method]
                if data['num_links']:
                    sorted_data = sorted(zip(data['num_links'], data['times']))
                    num_links_vals, times = zip(*sorted_data)
                    ax.plot(num_links_vals, times, marker='s', label=method, linewidth=2, markersize=6)
            
            ax.set_xlabel('Number of Links', fontsize=11)
            ax.set_ylabel('Total Execution Time (s)', fontsize=11)
            ax.set_title(f'Link Scaling (ε_base={epsilon_base_float:.0e})', fontsize=12, fontweight='bold')
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)
            ax.set_yscale('log')
            ax.set_xscale('log')
        
        plt.tight_layout()
        filename = output_dir / 'link_scaling.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved plot: {filename}")

    # Plot chunk size scaling: subplots for each epsilon
    if scaling_results.get('chunk_size_scaling') and scaling_config.get('chunk_size_scaling', {}).get('enabled', False):
        print("\nGenerating chunk_size_scaling plot...")
        chunk_config = scaling_config.get('chunk_size_scaling', {})
        epsilons_list = chunk_config.get('epsilons', [1e-2])
        
        fig, axes = plt.subplots(1, len(epsilons_list), figsize=(6*len(epsilons_list), 5))
        if len(epsilons_list) == 1:
            axes = [axes]
        
        for idx, epsilon_base in enumerate(epsilons_list):
            ax = axes[idx]
            methods_data = defaultdict(lambda: {'chunk_sizes': [], 'times': []})
            epsilon_base_float = float(epsilon_base)
            
            # Filter results for this epsilon
            for result in scaling_results['chunk_size_scaling']:
                if f'eps{epsilon_base_float:.0e}' in result.name:
                    methods_data[result.method]['chunk_sizes'].append(result.chunk_size or 2048)
                    methods_data[result.method]['times'].append(result.execution_time)
            
            # Plot each method
            for method in sorted(methods_data.keys()):
                data = methods_data[method]
                if data['chunk_sizes']:
                    sorted_data = sorted(zip(data['chunk_sizes'], data['times']))
                    chunk_sizes, times = zip(*sorted_data)
                    ax.plot(chunk_sizes, times, marker='^', label=method, linewidth=2, markersize=6)
            
            ax.set_xlabel('Chunk Size', fontsize=11)
            ax.set_ylabel('Total Execution Time (s)', fontsize=11)
            ax.set_title(f'Chunk Size Scaling (ε_base={epsilon_base_float:.0e})', fontsize=12, fontweight='bold')
            ax.legend(fontsize=9)
            ax.grid(True, alpha=0.3)
            ax.set_yscale('log')
            ax.set_xscale('log')
        
        plt.tight_layout()
        filename = output_dir / 'chunk_size_scaling.png'
        plt.savefig(filename, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved plot: {filename}")

def main():
    """Main benchmark runner"""
    # Load configuration
    config_path = Path('ZACH_benchmark_config.yaml')
    if not config_path.exists():
        print(f"Error: Configuration file {config_path} not found")
        return
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)

    import cupy as cp
    # print(cp.show_config())
    
    print("="*60)
    print("LINK DISTANCE BENCHMARK")
    print("="*60)

    # Dummy run to initialize GPU and avoid first-run overheads
    print("Initializing GPU with dummy call to run_gpu_method, to avoid first-run overheads...")
    dummy_links, _ = generate_random_links(2, config['link_generation'], seed=0)
    _ = run_gpu_method(dummy_links, 1e-3, 0.1, config['collision_analysis'])
    
    if config['verification_tests'].get('enabled', False):
        # Run verification tests
        print("\n" + "="*60)
        print("VERIFICATION TESTS")
        print("="*60)
        
        all_passed = True
        for test_config in config['verification_tests']['tests']:
            results, passed, message = run_single_test(
                test_config, 
                config['link_generation'],
                config['collision_analysis'],
                config['timeouts']['serial_method_timeout']
            )
            
            if passed:
                print(f"  ✓ PASSED: {message}")
            else:
                print(f"  ✗ FAILED: {message}")
                all_passed = False
        
        if all_passed:
            print("\n✓ All verification tests passed!")
        else:
            print("\n✗ Some verification tests failed")
    
    # Run scaling tests
    scaling_results = run_scaling_tests(config)
    
    # Generate plots
    output_dir = Path('benchmark_results')
    plot_results(scaling_results, output_dir, config)
    
    print("\n" + "="*60)
    print("BENCHMARK COMPLETE")
    print("="*60)


if __name__ == '__main__':
    main()
