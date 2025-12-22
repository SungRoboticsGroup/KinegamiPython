"""gpu_link_sdf.py

CuPy vectorized SDF for many LinkCSC.

Idea:
- Pack each LinkCSC into arrays (2 arcs + 1 straight + endpoint spheres).
- Evaluate point->link signed distance on GPU in chunks.
"""

import numpy as np
import cupy as cp


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


def _sd_capsule_cupy(p, a, b, r):
    # p: (N,3)
    # a,b: (L,3)
    # r: (L,)
    # return: (N,L)

    pa = p[:, None, :] - a[None, :, :] # (N,L,3)
    ba = b - a # (L,3)
    ba2 = cp.sum(ba * ba, axis=1) # (L,)

    dot_pa_ba = cp.sum(pa * ba[None, :, :], axis=2)  # (N,L)
    # prevent division by zero
    h = cp.where(ba2[None, :] > 0, dot_pa_ba / ba2[None, :], 0.0)  # (N,L)
    h = cp.clip(h, 0.0, 1.0)                      # (N,L)

    closest = pa - ba[None, :, :] * h[:, :, None] # (N,L,3)
    return cp.sqrt(cp.sum(closest * closest, axis=2)) - r[None, :]  # (N,L)


def _sd_flat_ended_torus_cupy(localP, sc, ra, rb):
    # localP: (N,L,3)
    # sc: (L,2)
    # ra, rb: (L,)
    # return: (N,L)

    px = cp.abs(localP[..., 0]) # (N,L)
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

    p_len = cp.sqrt(px * px + py * py) # (N,L)
    k = cp.where(scy * px > scx * py, scx * px + scy * py, p_len)  # (N,L)
    base = cp.sqrt(p_len * p_len + pz * pz + ra_ * ra_ - 2.0 * ra_ * k) - rb_  # (N,L)

    radialInPlane = toPoint_x * scx + toPoint_y * scy  # (N,L)
    discDist = cp.sqrt(radialInPlane * radialInPlane + pz * pz)  # (N,L)
    outsideDisc = cp.maximum(discDist - rb_, 0.0)      # (N,L)
    disc = cp.sqrt(pastEnd * pastEnd + outsideDisc * outsideDisc)  # (N,L)

    return cp.where(pastEnd <= 0.0, base, disc)


def _sd_arc_cupy(p, center, R_w2l, sc, ra, rb):
    # localP[n,l,:] = R[l] @ (p[n] - center[l])
    # p: (N,3)
    # center: (L,3)
    # R_w2l: (L,3,3)
    dp = p[:, None, :] - center[None, :, :]          # (N,L,3)
    # matrix multiplication for each N, L pair
    localP = cp.einsum("lij,nlj->nli", R_w2l, dp)   # (N,L,3)
    return _sd_flat_ended_torus_cupy(localP, sc, ra, rb)


def distances_points_to_links_cupy(points, packed, chunk_points=2048, dtype="float32"):
    """Returns (P,L) signed distances (can be huge)."""
    p = cp.asarray(points, dtype=dtype)
    L = packed.arc1_enabled.shape[0]

    has_arc1 = bool(np.any(packed.arc1_enabled))
    has_arc2 = bool(np.any(packed.arc2_enabled))
    has_seg = True

    # Upload packed arrays
    arc1_enabled = cp.asarray(packed.arc1_enabled)
    arc1_center = cp.asarray(packed.arc1_center, dtype=dtype)
    arc1_R_w2l = cp.asarray(packed.arc1_R_w2l, dtype=dtype)
    arc1_sc = cp.asarray(packed.arc1_sc, dtype=dtype)
    arc1_ra = cp.asarray(packed.arc1_ra, dtype=dtype)
    arc1_rb = cp.asarray(packed.arc1_rb, dtype=dtype)

    arc2_enabled = cp.asarray(packed.arc2_enabled)
    arc2_center = cp.asarray(packed.arc2_center, dtype=dtype)
    arc2_R_w2l = cp.asarray(packed.arc2_R_w2l, dtype=dtype)
    arc2_sc = cp.asarray(packed.arc2_sc, dtype=dtype)
    arc2_ra = cp.asarray(packed.arc2_ra, dtype=dtype)
    arc2_rb = cp.asarray(packed.arc2_rb, dtype=dtype)

    seg_enabled = cp.asarray(packed.seg_enabled)
    seg_a = cp.asarray(packed.seg_a, dtype=dtype)
    seg_b = cp.asarray(packed.seg_b, dtype=dtype)
    seg_r = cp.asarray(packed.seg_r, dtype=dtype)

    P = p.shape[0]
    out = cp.empty((P, L), dtype=dtype)

    for s in range(0, P, chunk_points):
        e = min(P, s + chunk_points)
        pp = p[s:e]

        dist = cp.full((e - s, L), cp.inf, dtype=dtype)

        if has_arc1:
            d1 = _sd_arc_cupy(pp, arc1_center, arc1_R_w2l, arc1_sc, arc1_ra, arc1_rb)
            dist = cp.minimum(dist, cp.where(arc1_enabled[None, :], d1, cp.inf))

        if has_seg:
            d2 = _sd_capsule_cupy(pp, seg_a, seg_b, seg_r)
            dist = cp.minimum(dist, cp.where(seg_enabled[None, :], d2, cp.inf))

        if has_arc2:
            d3 = _sd_arc_cupy(pp, arc2_center, arc2_R_w2l, arc2_sc, arc2_ra, arc2_rb)
            dist = cp.minimum(dist, cp.where(arc2_enabled[None, :], d3, cp.inf))

        out[s:e] = dist

    return out


def min_distance_to_other_links_cupy(points, point_link_ids, packed, chunk_points=2048, dtype="float32"):
    """For each point: min signed distance to any other link."""
    p = cp.asarray(points, dtype=dtype)
    ids = cp.asarray(point_link_ids, dtype=cp.int32)

    L = packed.arc1_enabled.shape[0]

    has_arc1 = bool(np.any(packed.arc1_enabled))
    has_arc2 = bool(np.any(packed.arc2_enabled))
    has_seg = True

    # Upload packed arrays (same as distances_points_to_links_cupy)
    arc1_enabled = cp.asarray(packed.arc1_enabled)
    arc1_center = cp.asarray(packed.arc1_center, dtype=dtype)
    arc1_R_w2l = cp.asarray(packed.arc1_R_w2l, dtype=dtype)
    arc1_sc = cp.asarray(packed.arc1_sc, dtype=dtype)
    arc1_ra = cp.asarray(packed.arc1_ra, dtype=dtype)
    arc1_rb = cp.asarray(packed.arc1_rb, dtype=dtype)

    arc2_enabled = cp.asarray(packed.arc2_enabled)
    arc2_center = cp.asarray(packed.arc2_center, dtype=dtype)
    arc2_R_w2l = cp.asarray(packed.arc2_R_w2l, dtype=dtype)
    arc2_sc = cp.asarray(packed.arc2_sc, dtype=dtype)
    arc2_ra = cp.asarray(packed.arc2_ra, dtype=dtype)
    arc2_rb = cp.asarray(packed.arc2_rb, dtype=dtype)

    seg_enabled = cp.asarray(packed.seg_enabled)
    seg_a = cp.asarray(packed.seg_a, dtype=dtype)
    seg_b = cp.asarray(packed.seg_b, dtype=dtype)
    seg_r = cp.asarray(packed.seg_r, dtype=dtype)

    P = p.shape[0]
    min_other = cp.empty((P,), dtype=dtype)
    argmin_other = cp.empty((P,), dtype=cp.int32)

    for s in range(0, P, chunk_points):
        e = min(P, s + chunk_points)
        pp = p[s:e]
        own = ids[s:e]

        dist = cp.full((e - s, L), cp.inf, dtype=dtype)

        if has_arc1:
            d1 = _sd_arc_cupy(pp, arc1_center, arc1_R_w2l, arc1_sc, arc1_ra, arc1_rb)
            dist = cp.minimum(dist, cp.where(arc1_enabled[None, :], d1, cp.inf))

        if has_seg:
            d2 = _sd_capsule_cupy(pp, seg_a, seg_b, seg_r)
            dist = cp.minimum(dist, cp.where(seg_enabled[None, :], d2, cp.inf))

        if has_arc2:
            d3 = _sd_arc_cupy(pp, arc2_center, arc2_R_w2l, arc2_sc, arc2_ra, arc2_rb)
            dist = cp.minimum(dist, cp.where(arc2_enabled[None, :], d3, cp.inf))

        # Mask out same-link distances
        dist[cp.arange(e - s), own] = cp.inf

        arg = cp.argmin(dist, axis=1)
        dmin = dist[cp.arange(e - s), arg]

        min_other[s:e] = dmin
        argmin_other[s:e] = arg.astype(cp.int32)

    return min_other, argmin_other


def collision_error_from_min_dist(min_other, margin=0.0, power=2.0):
    """Signed distance -> penalty. Negative means penetration."""
    # If margin=0: penalty = relu(-d)^power
    x = cp.maximum(margin - min_other, 0.0)
    if power == 1.0:
        return x
    return x**power
