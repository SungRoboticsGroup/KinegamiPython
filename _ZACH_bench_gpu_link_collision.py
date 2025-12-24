"""bench_gpu_link_collision.py

Quick speed test for GPU distance checking.
"""

import os
import time

import numpy as np
from spatialmath import SE3, SO3

from LinkCSC import LinkCSC
from ZACH_vectorized_link_sdf import (
    collision_error_from_min_dist,
    min_distance_to_other_links_cupy,
    pack_links,
    sample_points_for_links,
)


def _random_rotation(rng):
    axis = rng.normal(size=3)
    axis /= np.linalg.norm(axis) + 1e-12
    angle = rng.uniform(0, 2 * np.pi)
    return SO3.AngleAxis(angle, axis)


def _random_transform(rng, pos_scale=5.0):
    t = rng.uniform(-pos_scale, pos_scale, size=3)
    R = _random_rotation(rng)
    return SE3.Rt(R, t)


def build_base_link():
    r = 0.5
    start = SE3.Rt(SO3(), np.array([0.0, 0.0, 0.0]))
    end = SE3.Rt(SO3(), np.array([2.0, 0.0, 0.0]))
    return LinkCSC(r, start, end, maxAnglePerElbow=np.deg2rad(45))


def main():
    N = int(os.getenv("N_LINKS", "1000"))
    density = float(os.getenv("DENSITY", "5"))
    chunk_points = int(os.getenv("CHUNK_POINTS", "2048"))

    rng = np.random.default_rng(0)

    base = build_base_link()
    links = [base.newLinkTransformedBy(_random_transform(rng)) for _ in range(N)]

    t0 = time.perf_counter()
    packed = pack_links(links)
    t1 = time.perf_counter()

    pts, ids = sample_points_for_links(links, density=density)
    t2 = time.perf_counter()

    # GPU
    min_other, _ = min_distance_to_other_links_cupy(
        pts, ids, packed, chunk_points=chunk_points, dtype="float32"
    )
    # synchronize by pulling 1 value
    _ = float(min_other[0].get()) if min_other.size > 0 else 0.0
    t3 = time.perf_counter()

    penalty = collision_error_from_min_dist(min_other, margin=0.0, power=2.0)
    total_penalty = float(penalty.sum().get()) if penalty.size > 0 else 0.0
    t4 = time.perf_counter()

    print("=== GPU Link Collision Benchmark ===")
    print(f"Links: {N}")
    print(f"Sample density: {density} pts/unit")
    print(f"Total points: {pts.shape[0]}")
    print(f"Chunk points: {chunk_points}")
    print("")
    print(f"Pack links:      {t1 - t0:.3f}s")
    print(f"Sample points:   {t2 - t1:.3f}s (CPU)")
    print(f"GPU min-other:   {t3 - t2:.3f}s")
    print(f"Penalty reduce:  {t4 - t3:.3f}s")
    print(f"Total penalty:   {total_penalty:.6e}")


if __name__ == "__main__":
    main()
