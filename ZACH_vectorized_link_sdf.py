"""vectorized_link_sdf.py

Vectorized SDF for many LinkCSC using array library abstraction.
Can use either CuPy (GPU) or NumPy (CPU) by setting xp parameter.

Idea:
- Pack each LinkCSC into arrays (2 arcs + 1 straight + endpoint spheres).
- Evaluate point->link signed distance on GPU/CPU in chunks.
"""

import numpy as np


class PackedLinks:
    def __init__(
        self,
        arc1_enabled,
        arc1_center,
        arc1_R_w2l,
        arc1_sc,
        arc1_ra,
        arc1_rb,
        arc2_enabled,
        arc2_center,
        arc2_R_w2l,
        arc2_sc,
        arc2_ra,
        arc2_rb,
        seg_enabled,
        seg_a,
        seg_b,
        seg_r,
    ):
        self.arc1_enabled = arc1_enabled
        self.arc1_center = arc1_center
        self.arc1_R_w2l = arc1_R_w2l
        self.arc1_sc = arc1_sc
        self.arc1_ra = arc1_ra
        self.arc1_rb = arc1_rb

        self.arc2_enabled = arc2_enabled
        self.arc2_center = arc2_center
        self.arc2_R_w2l = arc2_R_w2l
        self.arc2_sc = arc2_sc
        self.arc2_ra = arc2_ra
        self.arc2_rb = arc2_rb

        self.seg_enabled = seg_enabled
        self.seg_a = seg_a
        self.seg_b = seg_b
        self.seg_r = seg_r


def pack_links(links):
    """Pack LinkCSC data into numpy arrays."""

    links = list(links)
    L = len(links)

    arc1_enabled = np.zeros((L,), dtype=bool)
    arc1_center = np.zeros((L, 3), dtype=float)
    arc1_R_w2l = np.zeros((L, 3, 3), dtype=float)
    arc1_sc = np.zeros((L, 2), dtype=float)
    arc1_ra = np.zeros((L,), dtype=float)
    arc1_rb = np.zeros((L,), dtype=float)

    arc2_enabled = np.zeros((L,), dtype=bool)
    arc2_center = np.zeros((L, 3), dtype=float)
    arc2_R_w2l = np.zeros((L, 3, 3), dtype=float)
    arc2_sc = np.zeros((L, 2), dtype=float)
    arc2_ra = np.zeros((L,), dtype=float)
    arc2_rb = np.zeros((L,), dtype=float)

    seg_enabled = np.zeros((L,), dtype=bool)
    seg_a = np.zeros((L, 3), dtype=float)
    seg_b = np.zeros((L, 3), dtype=float)
    seg_r = np.zeros((L,), dtype=float)

    for i, link in enumerate(links):
        r = float(link.r)

        # Arc1
        if getattr(link, "arc1", None) is not None and link.path.theta1 > link.EPSILON:
            arc1 = link.arc1
            if not hasattr(arc1, "_worldToLocalRotation"):
                arc1._computeLocalFrame()
            arc1_enabled[i] = True
            arc1_center[i] = np.asarray(arc1.circleCenter, dtype=float)
            arc1_R_w2l[i] = np.asarray(arc1._worldToLocalRotation, dtype=float)
            arc1_sc[i] = np.asarray(arc1._sdfSinCos, dtype=float)
            arc1_ra[i] = float(arc1.r)
            arc1_rb[i] = r

        # Straight segment (capsule). If length ~0, this becomes a sphere at a.
        seg_enabled[i] = True
        a = np.asarray(link.path.turn1end, dtype=float)
        if link.path.tMag > link.DISTANCE_EPSILON:
            b = np.asarray(link.path.turn1end + link.path.tMag * link.path.tUnit, dtype=float)
        else:
            b = a
        seg_a[i] = a
        seg_b[i] = b
        seg_r[i] = r

        # Arc2
        if getattr(link, "arc2", None) is not None and link.path.theta2 > link.EPSILON:
            arc2 = link.arc2
            if not hasattr(arc2, "_worldToLocalRotation"):
                arc2._computeLocalFrame()
            arc2_enabled[i] = True
            arc2_center[i] = np.asarray(arc2.circleCenter, dtype=float)
            arc2_R_w2l[i] = np.asarray(arc2._worldToLocalRotation, dtype=float)
            arc2_sc[i] = np.asarray(arc2._sdfSinCos, dtype=float)
            arc2_ra[i] = float(arc2.r)
            arc2_rb[i] = r

    return PackedLinks(
        arc1_enabled,
        arc1_center,
        arc1_R_w2l,
        arc1_sc,
        arc1_ra,
        arc1_rb,
        arc2_enabled,
        arc2_center,
        arc2_R_w2l,
        arc2_sc,
        arc2_ra,
        arc2_rb,
        seg_enabled,
        seg_a,
        seg_b,
        seg_r,
    )


def sample_points_for_links(links, density):
    """Sample points on each link (CPU)."""

    links = list(links)
    all_points = []
    all_ids = []
    for i, link in enumerate(links):
        pts = link.interpolate(density=density)
        all_points.append(pts)
        all_ids.append(np.full((pts.shape[0],), i, dtype=np.int32))

    points = np.concatenate(all_points, axis=0) if all_points else np.zeros((0, 3), dtype=float)
    ids = np.concatenate(all_ids, axis=0) if all_ids else np.zeros((0,), dtype=np.int32)
    return points, ids


def _sd_capsule(xp, p, a, b, r):
    # p: (N,3)
    # a,b: (L,3)
    # r: (L,)
    # return: (N,L)

    pa = p[:, None, :] - a[None, :, :] # (N,L,3)
    ba = b - a # (L,3)
    ba2 = xp.sum(ba * ba, axis=1) # (L,)

    dot_pa_ba = xp.sum(pa * ba[None, :, :], axis=2)  # (N,L)
    # prevent division by zero
    h = xp.where(ba2[None, :] > 0, dot_pa_ba / ba2[None, :], 0.0)  # (N,L)
    h = xp.clip(h, 0.0, 1.0)                      # (N,L)

    closest = pa - ba[None, :, :] * h[:, :, None] # (N,L,3)
    return xp.sqrt(xp.sum(closest * closest, axis=2)) - r[None, :]  # (N,L)


def _sd_flat_ended_torus(xp, localP, sc, ra, rb):
    # localP: (N,L,3)
    # sc: (L,2)
    # ra, rb: (L,)
    # return: (N,L)

    px = xp.abs(localP[..., 0]) # (N,L)
    py = localP[..., 1] # (N,L)
    pz = localP[..., 2] # (N,L)

    scx = sc[None, :, 0] # (1,L)
    scy = sc[None, :, 1] # (1,L)
    ra_ = ra[None, :] # (1,L)
    rb_ = rb[None, :] # (1,L)

    endCenter_x = ra_ * scx # (1,L)
    endCenter_y = ra_ * scy # (1,L)
    tangent_x = scy # (1,L)
    tangent_y = -scx # (1,L)

    toPoint_x = px - endCenter_x # (N,L)
    toPoint_y = py - endCenter_y # (N,L)
    pastEnd = toPoint_x * tangent_x + toPoint_y * tangent_y  # (N,L)

    p_len = xp.sqrt(px * px + py * py) # (N,L)
    k = xp.where(scy * px > scx * py, scx * px + scy * py, p_len)  # (N,L)
    # Clamp the sqrt argument to prevent NaN from negative values due to numerical precision
    sqrt_arg = p_len * p_len + pz * pz + ra_ * ra_ - 2.0 * ra_ * k  # (N,L)
    sqrt_arg = xp.maximum(sqrt_arg, 0.0)  # Clamp negative values to 0
    base = xp.sqrt(sqrt_arg) - rb_  # (N,L)

    radialInPlane = toPoint_x * scx + toPoint_y * scy  # (N,L)
    discDist = xp.sqrt(radialInPlane * radialInPlane + pz * pz)  # (N,L)
    outsideDisc = xp.maximum(discDist - rb_, 0.0)      # (N,L)
    disc = xp.sqrt(pastEnd * pastEnd + outsideDisc * outsideDisc)  # (N,L)

    return xp.where(pastEnd <= 0.0, base, disc)


def _sd_arc(xp, p, center, R_w2l, sc, ra, rb):
    # localP[n,l,:] = R[l] @ (p[n] - center[l])
    # p: (N,3)
    # center: (L,3)
    # R_w2l: (L,3,3)
    dp = p[:, None, :] - center[None, :, :]          # (N,L,3)
    # matrix multiplication for each N, L pair
    localP = xp.einsum("lij,nlj->nli", R_w2l, dp)   # (N,L,3)
    return _sd_flat_ended_torus(xp, localP, sc, ra, rb)


def distances_points_to_links(xp, points, packed, chunk_points=2048, dtype="float32"):
    """Returns (N,L) signed distances"""
    p = xp.asarray(points, dtype=dtype)
    L = packed.arc1_enabled.shape[0]

    has_arc1 = bool(np.any(packed.arc1_enabled))
    has_arc2 = bool(np.any(packed.arc2_enabled))
    has_seg = True

    # Upload packed arrays
    arc1_enabled = xp.asarray(packed.arc1_enabled)
    arc1_center = xp.asarray(packed.arc1_center, dtype=dtype)
    arc1_R_w2l = xp.asarray(packed.arc1_R_w2l, dtype=dtype)
    arc1_sc = xp.asarray(packed.arc1_sc, dtype=dtype)
    arc1_ra = xp.asarray(packed.arc1_ra, dtype=dtype)
    arc1_rb = xp.asarray(packed.arc1_rb, dtype=dtype)

    arc2_enabled = xp.asarray(packed.arc2_enabled)
    arc2_center = xp.asarray(packed.arc2_center, dtype=dtype)
    arc2_R_w2l = xp.asarray(packed.arc2_R_w2l, dtype=dtype)
    arc2_sc = xp.asarray(packed.arc2_sc, dtype=dtype)
    arc2_ra = xp.asarray(packed.arc2_ra, dtype=dtype)
    arc2_rb = xp.asarray(packed.arc2_rb, dtype=dtype)

    seg_enabled = xp.asarray(packed.seg_enabled)
    seg_a = xp.asarray(packed.seg_a, dtype=dtype)
    seg_b = xp.asarray(packed.seg_b, dtype=dtype)
    seg_r = xp.asarray(packed.seg_r, dtype=dtype)

    # num points (N)
    N = p.shape[0]
    # out is our signed distance from point to sdf dubins (N, L)
    out = xp.empty((N, L), dtype=dtype)

    for st in range(0, N, chunk_points):
        ed = min(N, st + chunk_points) # to handle last case when P is not multiple of 2048
        pp = p[st:ed]

        dist = xp.full((ed - st, L), xp.inf, dtype=dtype)

        if has_arc1:
            d1 = _sd_arc(xp, pp, arc1_center, arc1_R_w2l, arc1_sc, arc1_ra, arc1_rb)
            dist = xp.minimum(dist, xp.where(arc1_enabled[None, :], d1, xp.inf))

        if has_seg:
            d2 = _sd_capsule(xp, pp, seg_a, seg_b, seg_r)
            dist = xp.minimum(dist, xp.where(seg_enabled[None, :], d2, xp.inf))

        if has_arc2:
            d3 = _sd_arc(xp, pp, arc2_center, arc2_R_w2l, arc2_sc, arc2_ra, arc2_rb)
            dist = xp.minimum(dist, xp.where(arc2_enabled[None, :], d3, xp.inf))

        out[st:ed] = dist

    return out


def min_distance_to_other_links(xp, points, point_link_ids, packed, chunk_points=2048, dtype="float32"):
    """For each point: min signed distance to any other link."""
    """Returns (N) minimum signed distances and link id"""
    p = xp.asarray(points, dtype=dtype)
    ids = xp.asarray(point_link_ids, dtype=xp.int32)

    L = packed.arc1_enabled.shape[0] # num links

    has_arc1 = bool(np.any(packed.arc1_enabled))
    has_arc2 = bool(np.any(packed.arc2_enabled))
    has_seg = True

    # Upload packed arrays (same as distances_points_to_links)
    arc1_enabled = xp.asarray(packed.arc1_enabled)
    arc1_center = xp.asarray(packed.arc1_center, dtype=dtype)
    arc1_R_w2l = xp.asarray(packed.arc1_R_w2l, dtype=dtype)
    arc1_sc = xp.asarray(packed.arc1_sc, dtype=dtype)
    arc1_ra = xp.asarray(packed.arc1_ra, dtype=dtype)
    arc1_rb = xp.asarray(packed.arc1_rb, dtype=dtype)

    arc2_enabled = xp.asarray(packed.arc2_enabled)
    arc2_center = xp.asarray(packed.arc2_center, dtype=dtype)
    arc2_R_w2l = xp.asarray(packed.arc2_R_w2l, dtype=dtype)
    arc2_sc = xp.asarray(packed.arc2_sc, dtype=dtype)
    arc2_ra = xp.asarray(packed.arc2_ra, dtype=dtype)
    arc2_rb = xp.asarray(packed.arc2_rb, dtype=dtype)

    seg_enabled = xp.asarray(packed.seg_enabled)
    seg_a = xp.asarray(packed.seg_a, dtype=dtype)
    seg_b = xp.asarray(packed.seg_b, dtype=dtype)
    seg_r = xp.asarray(packed.seg_r, dtype=dtype)

    N = p.shape[0]
    min_other = xp.empty((N,), dtype=dtype)
    argmin_other = xp.empty((N,), dtype=xp.int32)

    for s in range(0, N, chunk_points):
        e = min(N, s + chunk_points)
        pp = p[s:e]
        own = ids[s:e]

        dist = xp.full((e - s, L), xp.inf, dtype=dtype)

        if has_arc1:
            d1 = _sd_arc(xp, pp, arc1_center, arc1_R_w2l, arc1_sc, arc1_ra, arc1_rb)
            dist = xp.minimum(dist, xp.where(arc1_enabled[None, :], d1, xp.inf))

        if has_seg:
            d2 = _sd_capsule(xp, pp, seg_a, seg_b, seg_r)
            dist = xp.minimum(dist, xp.where(seg_enabled[None, :], d2, xp.inf))

        if has_arc2:
            d3 = _sd_arc(xp, pp, arc2_center, arc2_R_w2l, arc2_sc, arc2_ra, arc2_rb)
            dist = xp.minimum(dist, xp.where(arc2_enabled[None, :], d3, xp.inf))

        # Mask out same-link distances
        dist[xp.arange(e - s), own] = xp.inf

        arg = xp.argmin(dist, axis=1) # id of closest link for N points
        dmin = dist[xp.arange(e - s), arg] # the actual signed distance for the N points

        min_other[s:e] = dmin
        argmin_other[s:e] = arg.astype(xp.int32)

    return min_other, argmin_other

def pairwise_link_distances(xp, points, point_link_ids, packed, chunk_points=2048, dtype="float32"):
    """Pairwise minimum distances between all link pairs."""
    """Returns (L,L) distances, (L,L) point indices that achieve them"""
    p = xp.asarray(points, dtype=dtype)
    ids = xp.asarray(point_link_ids, dtype=xp.int32)

    L = packed.arc1_enabled.shape[0]

    has_arc1 = bool(np.any(packed.arc1_enabled))
    has_arc2 = bool(np.any(packed.arc2_enabled))
    has_seg = True

    arc1_enabled = xp.asarray(packed.arc1_enabled)
    arc1_center = xp.asarray(packed.arc1_center, dtype=dtype)
    arc1_R_w2l = xp.asarray(packed.arc1_R_w2l, dtype=dtype)
    arc1_sc = xp.asarray(packed.arc1_sc, dtype=dtype)
    arc1_ra = xp.asarray(packed.arc1_ra, dtype=dtype)
    arc1_rb = xp.asarray(packed.arc1_rb, dtype=dtype)

    arc2_enabled = xp.asarray(packed.arc2_enabled)
    arc2_center = xp.asarray(packed.arc2_center, dtype=dtype)
    arc2_R_w2l = xp.asarray(packed.arc2_R_w2l, dtype=dtype)
    arc2_sc = xp.asarray(packed.arc2_sc, dtype=dtype)
    arc2_ra = xp.asarray(packed.arc2_ra, dtype=dtype)
    arc2_rb = xp.asarray(packed.arc2_rb, dtype=dtype)

    seg_enabled = xp.asarray(packed.seg_enabled)
    seg_a = xp.asarray(packed.seg_a, dtype=dtype)
    seg_b = xp.asarray(packed.seg_b, dtype=dtype)
    seg_r = xp.asarray(packed.seg_r, dtype=dtype)

    # (L,L) output: pairwise_dist[i,j] = min distance from link i to link j
    pairwise_dist = xp.full((L, L), xp.inf, dtype=dtype)
    # (L,L) output: point_idx[i,j] = global point index on link i that achieves min to link j
    point_idx = xp.full((L, L), -1, dtype=xp.int32)

    # For GPU execution, keep results on CPU to avoid repeated transfers
    if xp != np:
        pairwise_dist_cpu = np.full((L, L), np.inf, dtype=dtype)
        point_idx_cpu = np.full((L, L), -1, dtype=np.int32)
    else:
        pairwise_dist_cpu = pairwise_dist
        point_idx_cpu = point_idx

    N = p.shape[0]

    for s in range(0, N, chunk_points):
        e = min(N, s + chunk_points)
        pp = p[s:e]
        own = ids[s:e]

        dist = xp.full((e - s, L), xp.inf, dtype=dtype)

        if has_arc1:
            d1 = _sd_arc(xp, pp, arc1_center, arc1_R_w2l, arc1_sc, arc1_ra, arc1_rb)
            dist = xp.minimum(dist, xp.where(arc1_enabled[None, :], d1, xp.inf))

        if has_seg:
            d2 = _sd_capsule(xp, pp, seg_a, seg_b, seg_r)
            dist = xp.minimum(dist, xp.where(seg_enabled[None, :], d2, xp.inf))

        if has_arc2:
            d3 = _sd_arc(xp, pp, arc2_center, arc2_R_w2l, arc2_sc, arc2_ra, arc2_rb)
            dist = xp.minimum(dist, xp.where(arc2_enabled[None, :], d3, xp.inf))

        # Convert to CPU arrays for loop processing (avoid GPU-CPU sync per iteration)
        # If xp is CuPy, use .get() for explicit conversion
        # If xp is NumPy, arrays are already on CPU
        
        if xp != np:
            # CuPy path: use .get() for explicit GPU->CPU transfer
            dist_cpu = dist.get()
            own_cpu = own.get()
        else:
            # NumPy path: arrays already on CPU
            dist_cpu = dist
            own_cpu = own
        
        # VECTORIZED update: eliminate nested Python loops
        # For each unique link in this chunk, find which points originate from it
        unique_links = np.unique(own_cpu)
        for link_i in unique_links:
            link_i = int(link_i)
            # Find all points from this link
            point_indices_in_chunk = np.where(own_cpu == link_i)[0]
            global_indices = s + point_indices_in_chunk
            
            if len(point_indices_in_chunk) == 0:
                continue
            
            # Get distances for all these points to all links
            dist_subset = dist_cpu[point_indices_in_chunk, :]  # (num_points, L)
            
            # For each target link j (j != link_i), find best point
            for link_j in range(L):
                if link_j == link_i:
                    continue
                
                # Get distances from link_i points to link_j
                dists_to_j = dist_subset[:, link_j]  # (num_points,)
                
                # Find minimum distance, ignoring NaN/Inf values
                valid_mask = ~(np.isnan(dists_to_j) | np.isinf(dists_to_j))
                if np.any(valid_mask):
                    valid_indices = np.where(valid_mask)[0]
                    valid_dists = dists_to_j[valid_indices]
                    best_valid_idx = np.argmin(valid_dists)
                    best_idx_in_subset = valid_indices[best_valid_idx]
                    best_dist = valid_dists[best_valid_idx]
                else:
                    # All values are NaN/Inf, skip
                    continue
                
                # Update if better
                if best_dist < pairwise_dist_cpu[link_i, link_j]:
                    pairwise_dist_cpu[link_i, link_j] = best_dist
                    point_idx_cpu[link_i, link_j] = global_indices[best_idx_in_subset]

    # Return CPU arrays
    return pairwise_dist_cpu, point_idx_cpu
