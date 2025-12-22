import time
import numpy as np
import cupy as cp

def cupy_collisions(positions, threshold):
    # positions: (n,3)
    xp = cp
    pos = xp.asarray(positions, dtype=xp.float32)
    # pairwise squared distances (n,n) -> memory O(n^2)
    D2 = xp.sum((pos[:, None, :] - pos[None, :, :])**2, axis=2)
    thr2 = float(threshold * threshold)
    # mask upper triangle to avoid duplicates and self-comparisons
    n = pos.shape[0]
    mask = xp.triu(xp.ones((n, n), dtype=bool), k=1)
    coll = xp.where((D2 < thr2) & mask)
    return xp.stack(coll, axis=1).get()  # return numpy array of pairs

def cpu_collisions(positions, threshold):
    pos = np.asarray(positions, dtype=np.float32)
    D2 = np.sum((pos[:, None, :] - pos[None, :, :])**2, axis=2)
    thr2 = float(threshold * threshold)
    n = pos.shape[0]
    mask = np.triu(np.ones((n, n), dtype=bool), k=1)
    i, j = np.where((D2 < thr2) & mask)
    return np.stack((i, j), axis=1)

def main(n_points=10000, threshold=1.0):
    np.random.seed(0)
    positions = np.random.uniform(-100, 100, size=(n_points, 3)).astype(np.float32)

    print(f"Checking collisions among {n_points} random 3D points (threshold={threshold})")

    t0 = time.time()
    cpu_pairs = cpu_collisions(positions, threshold)
    t_cpu = time.time() - t0
    print(f"CPU: found {len(cpu_pairs)} pairs in {t_cpu:.4f} s")

    # Warm up GPU and measure
    cp.cuda.Device(0).synchronize()
    t0 = time.time()
    gpu_pairs = cupy_collisions(positions, threshold)
    cp.cuda.Device(0).synchronize()
    t_gpu = time.time() - t0
    print(f"GPU: found {len(gpu_pairs)} pairs in {t_gpu:.4f} s")

    # Compare results (unordered)
    cpu_set = set((int(a), int(b)) for a,b in cpu_pairs)
    gpu_set = set((int(a), int(b)) for a,b in gpu_pairs)

    if cpu_set == gpu_set:
        print("Results match ✓")
    else:
        only_cpu = cpu_set - gpu_set
        only_gpu = gpu_set - cpu_set
        print(f"Mismatch: only in CPU={len(only_cpu)}, only in GPU={len(only_gpu)}")

if __name__ == "__main__":
    main()
