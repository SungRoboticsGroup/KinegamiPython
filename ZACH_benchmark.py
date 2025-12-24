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
    point_density: float
    execution_time: float
    pairwise_distances: np.ndarray
    pairwise_point_indices: np.ndarray
    collisions: List
    collision_stats: object
    success: bool = True
    error_message: str = ""
    timed_out: bool = False


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


def generate_random_links(num_links: int, config: Dict, seed: int) -> List[LinkCSC]:
    """Generate random LinkCSC objects"""
    rng = np.random.RandomState(seed)
    links = []
    
    r_min, r_max = config['radius_range']
    p_min, p_max = config['position_range']
    max_angle = config['max_angle_per_elbow']
    epsilon = config['epsilon']
    
    for i in range(num_links):
        # Random radius
        r = rng.uniform(r_min, r_max)
        
        # Random start and end poses
        start_seed = seed + i * 2
        end_seed = seed + i * 2 + 1
        start_pose = generate_random_SE3((p_min, p_max), start_seed)
        end_pose = generate_random_SE3((p_min, p_max), end_seed)
        
        # Ensure poses are not too close
        max_attempts = 10
        attempts = 0
        while np.linalg.norm(end_pose.t - start_pose.t) < 2 * r and attempts < max_attempts:
            end_seed += num_links * 2
            end_pose = generate_random_SE3((p_min, p_max), end_seed)
            attempts += 1
        
        if attempts >= max_attempts:
            continue
        
        try:
            link = LinkCSC(
                r=r,
                StartDubinsPose=start_pose,
                EndDubinsPose=end_pose,
                maxAnglePerElbow=max_angle,
                EPSILON=epsilon
            )
            links.append(link)
        except (ValueError, AssertionError) as e:
            continue
    
    return links


def run_gpu_method(links: List[LinkCSC], density: float, collision_config: Dict) -> BenchmarkResult:
    """Run GPU vectorized method"""
    try:
        import cupy as cp
        xp = cp
        
        # Pack links
        packed = pack_links(links)
        
        # Sample points
        points, point_ids = sample_points_for_links(links, density)
        
        # Time the pairwise computation
        start_time = time.time()
        pairwise_dist, point_idx = pairwise_link_distances(
            xp, points, point_ids, packed, chunk_points=2048, dtype="float32"
        )
        
        # Convert to CPU
        pairwise_dist_cpu = cp.asnumpy(pairwise_dist)
        point_idx_cpu = cp.asnumpy(point_idx)
        
        execution_time = time.time() - start_time
        
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
            point_density=density,
            execution_time=execution_time,
            pairwise_distances=pairwise_dist_cpu,
            pairwise_point_indices=point_idx_cpu,
            collisions=collisions,
            collision_stats=collision_stats,
            success=True
        )
    except Exception as e:
        return BenchmarkResult(
            name="gpu_test",
            method="GPU (CuPy)",
            num_links=len(links),
            num_points=0,
            point_density=density,
            execution_time=0.0,
            pairwise_distances=np.array([]),
            pairwise_point_indices=np.array([]),
            collisions=[],
            collision_stats=None,
            success=False,
            error_message=str(e)
        )


def run_cpu_vectorized_method(links: List[LinkCSC], density: float, collision_config: Dict) -> BenchmarkResult:
    """Run CPU vectorized method"""
    try:
        xp = np
        
        # Pack links
        packed = pack_links(links)
        
        # Sample points
        points, point_ids = sample_points_for_links(links, density)
        
        # Time the pairwise computation
        start_time = time.time()
        pairwise_dist, point_idx = pairwise_link_distances(
            xp, points, point_ids, packed, chunk_points=2048, dtype="float32"
        )
        execution_time = time.time() - start_time
        
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
            point_density=density,
            execution_time=execution_time,
            pairwise_distances=pairwise_dist,
            pairwise_point_indices=point_idx,
            collisions=collisions,
            collision_stats=collision_stats,
            success=True
        )
    except Exception as e:
        return BenchmarkResult(
            name="cpu_vectorized_test",
            method="CPU Vectorized (NumPy)",
            num_links=len(links),
            num_points=0,
            point_density=density,
            execution_time=0.0,
            pairwise_distances=np.array([]),
            pairwise_point_indices=np.array([]),
            collisions=[],
            collision_stats=None,
            success=False,
            error_message=str(e)
        )


def run_cpu_serial_method(links: List[LinkCSC], density: float, 
                          collision_config: Dict, timeout: int) -> BenchmarkResult:
    """Run CPU serial method with timeout"""
    try:
        # Sample points (using same method as vectorized for consistency)
        points, point_ids = sample_points_for_links(links, density)
        
        # Time the pairwise computation with timeout
        start_time = time.time()
        pairwise_dist, point_idx = pairwise_link_distances_serial(
            links, points, point_ids, timeout_seconds=timeout
        )
        execution_time = time.time() - start_time
        
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
            point_density=density,
            execution_time=execution_time,
            pairwise_distances=pairwise_dist,
            pairwise_point_indices=point_idx,
            collisions=collisions,
            collision_stats=collision_stats,
            success=True
        )
    except TimeoutException:
        return BenchmarkResult(
            name="cpu_serial_test",
            method="CPU Serial",
            num_links=len(links),
            num_points=len(points) if 'points' in locals() else 0,
            point_density=density,
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
            point_density=density,
            execution_time=0.0,
            pairwise_distances=np.array([]),
            pairwise_point_indices=np.array([]),
            collisions=[],
            collision_stats=None,
            success=False,
            error_message=str(e)
        )


def verify_results(results: List[BenchmarkResult], config: Dict) -> Tuple[bool, str]:
    """Verify that all methods produce the same results"""
    if len(results) < 2:
        return True, "Not enough results to compare"
    
    successful_results = [r for r in results if r.success]
    if len(successful_results) < 2:
        return False, "Not enough successful results to compare"
    
    rtol = float(config['distance_rtol'])
    atol = float(config['distance_atol'])
    allow_point_mismatch = config['allow_point_mismatch']
    
    base_result = successful_results[0]
    
    # Debug: Print distance statistics for all methods
    print(f"  Distance Statistics:")
    for r in successful_results:
        nan_count = np.sum(np.isnan(r.pairwise_distances))
        inf_count = np.sum(np.isinf(r.pairwise_distances))
        valid_mask = ~(np.isnan(r.pairwise_distances) | np.isinf(r.pairwise_distances))
        if np.any(valid_mask):
            valid_dists = r.pairwise_distances[valid_mask]
            print(f"    {r.method}: min={np.min(valid_dists):.6f}, max={np.max(valid_dists):.6f}, " +
                  f"mean={np.mean(valid_dists):.6f}, NaN={nan_count}, Inf={inf_count}")
        else:
            print(f"    {r.method}: All NaN or Inf! (NaN={nan_count}, Inf={inf_count})")
    
    for result in successful_results[1:]:
        # Check distances match (ignore NaN/Inf locations)
        base_valid = ~(np.isnan(base_result.pairwise_distances) | np.isinf(base_result.pairwise_distances))
        result_valid = ~(np.isnan(result.pairwise_distances) | np.isinf(result.pairwise_distances))
        
        # Both should have valid values in same locations
        if not np.array_equal(base_valid, result_valid):
            print(f"    Warning: NaN/Inf locations differ between {base_result.method} and {result.method}")
        
        # Compare only valid entries
        common_valid = base_valid & result_valid
        if np.any(common_valid):
            base_valid_vals = base_result.pairwise_distances[common_valid]
            result_valid_vals = result.pairwise_distances[common_valid]
            if not np.allclose(base_valid_vals, result_valid_vals, rtol=rtol, atol=atol, equal_nan=True):
                max_diff = np.max(np.abs(base_valid_vals - result_valid_vals))
                print(f"    Distance mismatch in valid entries: max diff = {max_diff}")
                return False, f"Distance mismatch between {base_result.method} and {result.method}: max diff = {max_diff}"
        else:
            print(f"    No common valid entries to compare!")
        
        # Check point indices match (if required)
        if not allow_point_mismatch:
            if not np.array_equal(base_result.pairwise_point_indices, 
                                 result.pairwise_point_indices):
                return False, f"Point index mismatch between {base_result.method} and {result.method}"
        
        # Check collisions match
        match, msg = compare_collision_detections(
            base_result.collisions, result.collisions, rtol, atol
        )
        if not match:
            return False, f"Collision mismatch between {base_result.method} and {result.method}: {msg}"
    
    return True, "All methods agree"


def run_single_test(test_config: Dict, link_config: Dict, 
                   collision_config: Dict, timeout: int) -> Tuple[List[BenchmarkResult], bool, str]:
    """Run all three methods on a single test configuration"""
    print(f"\nRunning test: {test_config['name']}")
    print(f"  Links: {test_config['num_links']}, Density: {test_config['point_density']}")
    
    # Generate links
    links = generate_random_links(
        test_config['num_links'],
        link_config,
        test_config['seed']
    )
    
    if len(links) == 0:
        return [], False, "Failed to generate any links"
    
    print(f"  Generated {len(links)} links")
    
    # Run all three methods
    results = []
    
    # GPU method
    print("  Running GPU method...")
    gpu_result = run_gpu_method(links, test_config['point_density'], collision_config)
    results.append(gpu_result)
    if gpu_result.success:
        print(f"    Time: {gpu_result.execution_time:.4f}s, Points: {gpu_result.num_points}")
        if collision_config['report_statistics']:
            print(f"    Collisions: {gpu_result.collision_stats.collision_pairs}")
    else:
        print(f"    Failed: {gpu_result.error_message}")
    
    # CPU vectorized method
    print("  Running CPU vectorized method...")
    cpu_vec_result = run_cpu_vectorized_method(links, test_config['point_density'], collision_config)
    results.append(cpu_vec_result)
    if cpu_vec_result.success:
        print(f"    Time: {cpu_vec_result.execution_time:.4f}s, Points: {cpu_vec_result.num_points}")
        if collision_config['report_statistics']:
            print(f"    Collisions: {cpu_vec_result.collision_stats.collision_pairs}")
    else:
        print(f"    Failed: {cpu_vec_result.error_message}")
    
    # CPU serial method
    print("  Running CPU serial method...")
    cpu_serial_result = run_cpu_serial_method(links, test_config['point_density'], 
                                              collision_config, timeout)
    results.append(cpu_serial_result)
    if cpu_serial_result.success:
        print(f"    Time: {cpu_serial_result.execution_time:.4f}s, Points: {cpu_serial_result.num_points}")
        if collision_config['report_statistics']:
            print(f"    Collisions: {cpu_serial_result.collision_stats.collision_pairs}")
    elif cpu_serial_result.timed_out:
        print(f"    Timed out after {timeout}s")
    else:
        print(f"    Failed: {cpu_serial_result.error_message}")
    
    return results, *verify_results(results, link_config)


def run_scaling_tests(config: Dict) -> Dict[str, List[BenchmarkResult]]:
    """Run scaling tests with varying parameters"""
    scaling_config = config['scaling_tests']
    link_config = config['link_generation']
    collision_config = config['collision_analysis']
    timeout = config['timeouts']['serial_method_timeout']
    
    results = {
        'density_scaling': [],
        'link_scaling': []
    }
    
    # Density scaling
    print("\n" + "="*60)
    print("DENSITY SCALING TESTS")
    print("="*60)
    
    base_config = scaling_config['density_scaling']
    for density in base_config['densities']:
        print(f"\nDensity: {density}")
        
        links = generate_random_links(
            base_config['num_links'],
            link_config,
            base_config['seed']
        )
        
        if len(links) == 0:
            continue
        
        # Run all three methods
        gpu_result = run_gpu_method(links, density, collision_config)
        cpu_vec_result = run_cpu_vectorized_method(links, density, collision_config)
        cpu_serial_result = run_cpu_serial_method(links, density, collision_config, timeout)
        
        if gpu_result.success:
            results['density_scaling'].append(gpu_result)
        if cpu_vec_result.success:
            results['density_scaling'].append(cpu_vec_result)
        if cpu_serial_result.success:
            results['density_scaling'].append(cpu_serial_result)
    
    # Link scaling
    print("\n" + "="*60)
    print("LINK SCALING TESTS")
    print("="*60)
    
    base_config = scaling_config['link_scaling']
    for num_links in base_config['num_links_list']:
        print(f"\nNum links: {num_links}")
        
        links = generate_random_links(
            num_links,
            link_config,
            base_config['seed']
        )
        
        if len(links) == 0:
            continue
        
        # Run all three methods
        gpu_result = run_gpu_method(links, base_config['point_density'], collision_config)
        cpu_vec_result = run_cpu_vectorized_method(links, base_config['point_density'], collision_config)
        cpu_serial_result = run_cpu_serial_method(links, base_config['point_density'], 
                                                  collision_config, timeout)
        
        if gpu_result.success:
            results['link_scaling'].append(gpu_result)
        if cpu_vec_result.success:
            results['link_scaling'].append(cpu_vec_result)
        if cpu_serial_result.success:
            results['link_scaling'].append(cpu_serial_result)
    
    return results


def plot_results(scaling_results: Dict[str, List[BenchmarkResult]], output_dir: Path):
    """Generate plots from scaling test results"""
    output_dir.mkdir(exist_ok=True)
    
    # Plot density scaling
    fig, ax = plt.subplots(figsize=(10, 6))
    
    density_results = scaling_results['density_scaling']
    methods = {}
    for result in density_results:
        if result.method not in methods:
            methods[result.method] = {'densities': [], 'times': []}
        methods[result.method]['densities'].append(result.point_density)
        methods[result.method]['times'].append(result.execution_time)
    
    for method, data in methods.items():
        ax.plot(data['densities'], data['times'], marker='o', label=method, linewidth=2)
    
    ax.set_xlabel('Point Density (points per unit length)', fontsize=12)
    ax.set_ylabel('Execution Time (seconds)', fontsize=12)
    ax.set_title('Execution Time vs Point Density', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'density_scaling.png', dpi=150)
    plt.close()
    
    # Plot link scaling
    fig, ax = plt.subplots(figsize=(10, 6))
    
    link_results = scaling_results['link_scaling']
    methods = {}
    for result in link_results:
        if result.method not in methods:
            methods[result.method] = {'num_links': [], 'times': []}
        methods[result.method]['num_links'].append(result.num_links)
        methods[result.method]['times'].append(result.execution_time)
    
    for method, data in methods.items():
        ax.plot(data['num_links'], data['times'], marker='o', label=method, linewidth=2)
    
    ax.set_xlabel('Number of Links', fontsize=12)
    ax.set_ylabel('Execution Time (seconds)', fontsize=12)
    ax.set_title('Execution Time vs Number of Links', fontsize=14, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    ax.set_xscale('log')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'link_scaling.png', dpi=150)
    plt.close()
    
    print(f"\nPlots saved to {output_dir}")


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
    
    # Run verification tests
    print("\n" + "="*60)
    print("VERIFICATION TESTS")
    print("="*60)
    
    all_passed = True
    for test_config in config['tests']:
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
        
        # Print collision statistics if enabled
        if config['collision_analysis']['report_statistics'] and results:
            for result in results:
                if result.success and result.collision_stats:
                    print_collision_report(result.method, result.collision_stats, show_worst=3)
    
    if all_passed:
        print("\n✓ All verification tests passed!")
    else:
        print("\n✗ Some verification tests failed")
    
    # Run scaling tests
    scaling_results = run_scaling_tests(config)
    
    # Generate plots
    output_dir = Path('benchmark_results')
    plot_results(scaling_results, output_dir)
    
    print("\n" + "="*60)
    print("BENCHMARK COMPLETE")
    print("="*60)


if __name__ == '__main__':
    main()
