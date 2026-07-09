"""Synthetic 3D void scenarios with known ground truth, for demos and verification.

A *lumen* / enclosed cavity is an **H2 void**: sample points on the boundary of a
cavity and its alpha complex has a 2-dimensional homology class whose **death
value equals the cavity radius** (verified: a shell of radius R gives H2 death R).
So "how big is the void" and "how many voids" are read directly off H2.

This module builds controllable scenarios where the void history is known by
construction, so a topological pipeline can be checked against ground truth:

  * :func:`coalescence` -- cavities fuse in PAIRS (H2 count 4 -> 2), each death a
    merge.  Fusion is the boundary of the UNION of two balls (a sphere's points
    are kept only where they are outside the other ball), so the two interiors
    genuinely merge into one cavity when the balls overlap.
  * :func:`ripening` -- Ostwald ripening, mass-conserving: the small cavities
    shrink to nothing (their H2 points migrate to the diagonal) while the large
    one absorbs their volume.  H2 count 4 -> 2 -- the SAME drop as coalescence,
    so the two are told apart only by the merge history, not by counting.
  * :func:`null_model` -- the control: cavities of fixed size (H2 count flat).

Each returns a :class:`Scenario` (frames of 3D points + per-frame ground-truth
voids), and :func:`h2_voids` reads the significant voids out of a frame with
GUDHI (the trusted oracle).  Requires gudhi (already a dependency).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Geometry: sampling the boundary of a union of balls
# ---------------------------------------------------------------------------

def _fib_sphere(n: int) -> np.ndarray:
    """``n`` roughly-even points on the unit sphere (Fibonacci spiral)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi),
                            np.sin(theta) * np.sin(phi),
                            np.cos(phi)])


def union_boundary(centers, radii, n_each: int, rng, jitter: float = 0.015):
    """Points on the boundary of the UNION of the given balls: the points of each
    sphere that lie outside every *other* ball.  Overlapping balls therefore share
    a single interior cavity (their voids merge).  Returns ``(points, labels)``
    where ``labels[i]`` is the index of the ball point ``i`` was sampled from."""
    base = _fib_sphere(n_each)
    pts, labels = [], []
    centers = [np.asarray(c, float) for c in centers]
    for k, (c, r) in enumerate(zip(centers, radii)):
        p = base * r + c + rng.normal(0, jitter * max(r, 1e-3), (n_each, 3))
        keep = np.ones(len(p), bool)
        for j, (c2, r2) in enumerate(zip(centers, radii)):
            if j != k:
                keep &= np.linalg.norm(p - c2, axis=1) > r2
        pts.append(p[keep])
        labels.append(np.full(int(keep.sum()), k))
    return np.vstack(pts), np.concatenate(labels)


def _sample_voids(voids, n_each, rng):
    """Sample the boundary of the ALIVE voids in a frame, keeping each point's
    ORIGINAL void index as its label (so a cavity keeps its colour when others
    die)."""
    idx = [i for i, v in enumerate(voids) if v.alive and v.radius > 0.02]
    if not idx:
        return np.zeros((0, 3)), np.zeros(0, int)
    pts, local = union_boundary([voids[i].center for i in idx],
                                [voids[i].radius for i in idx], n_each, rng)
    return pts, np.array([idx[int(l)] for l in local])


def _project_inside(p, c_other, r_other):
    """Move any point of ``p`` that lies inside the other ball onto that ball's
    surface.  Two overlapping shells then share a union boundary with NO points in
    the merged interior (interior points would fill the fused cavity) -- while
    every index survives, so the result is a FIXED-cardinality moving point cloud
    the exact vineyard can follow.  (Plain full spheres would interpenetrate and
    the fused cavity would read as gone, not merged.)"""
    p = np.array(p, float)
    d = np.linalg.norm(p - c_other, axis=1)
    inside = d < r_other
    if inside.any():
        p[inside] = c_other + r_other * (p[inside] - c_other) / d[inside, None]
    return p


# ---------------------------------------------------------------------------
# Ground truth + scenario container
# ---------------------------------------------------------------------------

@dataclass
class Void:
    """A ground-truth cavity in one frame."""
    center: np.ndarray
    radius: float
    alive: bool = True


@dataclass
class Scenario:
    """A synthetic void movie with known ground truth."""
    name: str
    times: np.ndarray                              # (n_frames,)
    frames: List[Tuple[np.ndarray, np.ndarray]]    # per frame: (points, labels)
    truth: List[List[Void]]                        # per frame: ground-truth voids
    caption: str = ""
    merges: List[Tuple[float, int, int]] = field(default_factory=list)
    # ground-truth merge events (time, absorbed_label, partner_label): a FUSION,
    # where one cavity is absorbed by another.  Empty for resorption/deaths.

    def true_count(self) -> List[int]:
        return [sum(v.alive for v in voids) for voids in self.truth]


# ---------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------

def coalescence(n_frames: int = 60, n_each: int = 240, *, radius: float = 1.0,
                seed: int = 0) -> Scenario:
    """Two PAIRS of cavities that each fuse into one (H2 count 4 -> 2) -- the SAME
    drop as mass-conserving :func:`ripening`, so counting cannot tell them apart;
    only the merge history can.  Fusion is the boundary of the union of two balls
    (interiors genuinely merge).  ``merges`` records the two fusion events.

    The merge threshold (``1.50*radius``) is calibrated to where the alpha complex
    of the sampled shells actually resolves one cavity -- past tangency, since two
    just-touching spheres are not one cavity yet -- so the ground-truth count
    matches the computed count except within ~2 frames of the fusion (the
    resolution scale: "shells touch" vs "the complex resolves one cavity")."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, 1.0, n_frames)
    d = np.interp(t, [0, 1], [2.25 * radius, 1.0 * radius])  # separation per pair
    merge_d = 1.50 * radius                                  # -> both pairs fuse at t~0.6
    Y = 3.0 * radius                                         # the two pairs, far apart
    frames, truth = [], []
    merge_time = None
    for tk, dk in zip(t, d):
        separate = dk >= merge_d
        if not separate and merge_time is None:
            merge_time = float(tk)
        C = [(-dk / 2, Y, 0), (dk / 2, Y, 0), (-dk / 2, -Y, 0), (dk / 2, -Y, 0)]
        # Always sample all four shells: once a pair overlaps, its two shells
        # bound ONE cavity (a dumbbell), so the H2 count drops 4 -> 2.  The
        # absorbed shell's points stay (they are part of the fused boundary);
        # only its status as a separate void ends.
        frames.append(union_boundary(C, [radius] * 4, n_each, rng))
        truth.append([Void(np.array(C[0]), radius, True),      # survivor 0
                      Void(np.array(C[1]), radius, separate),  # absorbed 1 -> 0
                      Void(np.array(C[2]), radius, True),      # survivor 2
                      Void(np.array(C[3]), radius, separate)]) # absorbed 3 -> 2
    mt = merge_time if merge_time is not None else 1.0
    merges = [(mt, 1, 0), (mt, 3, 2)]
    cap = "Coalescence: cavities fuse in pairs (count 4 -> 2) -- each death is a merge."
    return Scenario("coalescence", t, frames, truth, cap, merges=merges)


def ripening(n_frames: int = 60, n_each: int = 240, *, seed: int = 0,
             r_die: float = 0.34) -> Scenario:
    """Ostwald ripening, **mass-conserving**: the two small cavities shrink
    continuously to nothing (so their H2 points migrate to the diagonal and die
    at ~zero persistence) while the largest absorbs exactly the volume they lose.
    Ground truth: a void is alive while its radius > ``r_die`` (near the sampling
    resolution, i.e. it dies on the diagonal, not deleted at a finite size).  The
    count goes 4 -> 2 -- the SAME drop as coalescence, so the two are told apart
    only by the merge history, not by counting."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, 1.0, n_frames)
    centers = np.array([(-3.2, -1.6, 0.3), (3.0, -1.1, -0.2),
                        (-2.6, 3.0, -0.3), (3.1, 2.6, 0.2)], float)
    r0 = np.array([1.3, 1.0, 0.85, 0.85])                # two equal small voids
    frames, truth = [], []
    for tk in t:
        r2, r3 = r0[2] * (1 - tk), r0[3] * (1 - tk)      # -> 0; die (r<r_die) at t~0.6
        v_lost = (r0[2] ** 3 - r2 ** 3) + (r0[3] ** 3 - r3 ** 3)
        r_big = (r0[0] ** 3 + v_lost) ** (1 / 3)         # absorbs the lost volume
        radii = np.array([r_big, r0[1], r2, r3])         # void1 holds steady
        alive = radii > r_die
        voids = [Void(c, float(r), bool(a))
                 for c, r, a in zip(centers, radii, alive)]
        truth.append(voids)
        frames.append(_sample_voids(voids, n_each, rng))
    cap = "Ostwald ripening (mass-conserving): small cavities shrink to nothing, the big one grows."
    return Scenario("ripening", t, frames, truth, cap)


def null_model(n_frames: int = 60, n_each: int = 240, *, seed: int = 0) -> Scenario:
    """Control: the same cavities holding a fixed size (H2 count stays flat) -- the
    contrast that lets TDA distinguish real ripening from nothing happening."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, 1.0, n_frames)
    centers = np.array([(-3.2, -1.6, 0.3), (3.0, -1.1, -0.2),
                        (-2.6, 3.0, -0.3), (3.1, 2.6, 0.2)], float)
    r0 = np.array([1.3, 1.0, 0.8, 0.7])
    frames, truth = [], []
    for tk in t:
        radii = r0 + 0.04 * np.sin(2 * np.pi * tk)           # tiny breathing only
        voids = [Void(c, float(r), True) for c, r in zip(centers, radii)]
        truth.append(voids)
        frames.append(_sample_voids(voids, n_each, rng))
    cap = "Null model: cavities hold their size -- the H2 count stays flat."
    return Scenario("null", t, frames, truth, cap)


# ---------------------------------------------------------------------------
# Fixed-cardinality moving point clouds, for the EXACT engine (moving_vineyard)
#
# The scenarios above resample and cull points every frame (dense boundaries for
# a legible picture), so their point count varies and there is no point identity
# across frames -- fine for per-frame GUDHI, impossible for a vineyard.  The
# exact engine needs a moving point cloud: a FIXED set of points, each carried by
# index from frame to frame.  These builders provide that -- coarse shells (the
# exact engine needs no dense sampling) whose points move continuously -- so the
# cavity history comes from the maintained pairing, not from matching.
# ---------------------------------------------------------------------------

def coalescence_exact(n_each: int = 30, n_frames: int = 8, *, radius: float = 1.0,
                      seed: int = 0) -> Scenario:
    """Two cavities that FUSE into one (count 2 -> 1) as a fixed-cardinality moving
    point cloud.  As the two shells overlap, the absorbed cavity's vine dies
    OFF-diagonal (an abrupt fusion), the signature a fixed-complex vineyard cannot
    produce.  Interior points are projected to the union boundary so the merged
    cavity survives (see :func:`_project_inside`)."""
    base = _fib_sphere(n_each)
    t = np.linspace(0.0, 1.0, n_frames)
    a_sep = np.interp(t, [0, 1], [1.5 * radius, 0.35 * radius])   # half-separation
    a_merge = 0.55 * radius                                       # alpha resolves one cavity
    frames, truth, merge_time = [], [], None
    for tk, a in zip(t, a_sep):
        c0, c1 = np.array([-a, 0.0, 0.0]), np.array([a, 0.0, 0.0])
        p0 = _project_inside(base * radius + c0, c1, radius)
        p1 = _project_inside(base * radius + c1, c0, radius)
        pts = np.vstack([p0, p1])
        labels = np.concatenate([np.zeros(n_each, int), np.ones(n_each, int)])
        frames.append((pts, labels))
        merged = a < a_merge
        if merged and merge_time is None:
            merge_time = float(tk)
        truth.append([Void(c0, radius, True),           # survivor 0
                      Void(c1, radius, not merged)])     # absorbed 1 -> 0
    mt = merge_time if merge_time is not None else 1.0
    cap = "Coalescence (exact engine): two cavities fuse into one (2 -> 1)."
    return Scenario("coalescence_exact", t, frames, truth, cap,
                    merges=[(mt, 1, 0)])


def ripening_exact(n_each: int = 30, n_frames: int = 7, *, radius: float = 1.0,
                   seed: int = 0, r_die: float = 0.25) -> Scenario:
    """Two cavities, one holding its size and one shrinking to nothing (Ostwald
    resorption), as a fixed-cardinality moving point cloud.  The shrinking cavity's
    vine slides continuously to the DIAGONAL (persistence -> 0) -- the resorption
    signature, contrasting with coalescence's off-diagonal fusion death.  Count
    2 -> 1."""
    base = _fib_sphere(n_each)
    t = np.linspace(0.0, 1.0, n_frames)
    c0, c1 = np.array([-1.9 * radius, 0.0, 0.0]), np.array([1.9 * radius, 0.0, 0.0])
    frames, truth = [], []
    for tk in t:
        r1 = radius * (1.0 - 0.94 * tk)                  # shell 1 shrinks toward 0
        pts = np.vstack([base * radius + c0, base * r1 + c1])
        labels = np.concatenate([np.zeros(n_each, int), np.ones(n_each, int)])
        frames.append((pts, labels))
        truth.append([Void(c0, radius, True),
                      Void(c1, float(r1), bool(r1 > r_die))])
    cap = "Ostwald ripening (exact engine): a cavity shrinks to nothing (2 -> 1)."
    return Scenario("ripening_exact", t, frames, truth, cap)


def null_exact(n_each: int = 30, n_frames: int = 7, *, radius: float = 1.0,
               seed: int = 0) -> Scenario:
    """Two cavities holding their size (count stays 2), as a fixed-cardinality
    moving point cloud -- the control for the exact engine."""
    base = _fib_sphere(n_each)
    t = np.linspace(0.0, 1.0, n_frames)
    c0, c1 = np.array([-1.9 * radius, 0.0, 0.0]), np.array([1.9 * radius, 0.0, 0.0])
    frames, truth = [], []
    for tk in t:
        r0 = radius * (1.0 + 0.03 * np.sin(2 * np.pi * tk))   # tiny breathing only
        r1 = 0.85 * radius * (1.0 + 0.03 * np.cos(2 * np.pi * tk))
        pts = np.vstack([base * r0 + c0, base * r1 + c1])
        labels = np.concatenate([np.zeros(n_each, int), np.ones(n_each, int)])
        frames.append((pts, labels))
        truth.append([Void(c0, float(r0), True), Void(c1, float(r1), True)])
    return Scenario("null_exact", t, frames, truth,
                    "Null (exact engine): two cavities hold their size (stays 2).")


def two_pair_coalescence(n_each: int = 16, n_frames: int = 5, *, radius: float = 1.0,
                         seed: int = 0) -> Scenario:
    """The ADVERSARIAL multi-cavity case for the merge partner (>2 cavities, where
    the partner is a genuine CHOICE, not forced).  Two large survivors, far apart;
    a small cavity fuses into EACH -- ground truth: absorbed 1 -> survivor 0, and
    absorbed 3 -> survivor 2, with the two absorbed cavities of DIFFERENT size
    (0.7 vs 0.95) so a correct partner rule must discriminate.

    This scenario currently DEFEATS the genealogy: the H2 survivor vines fragment
    at each fusion flip, so :func:`~chromatic_cells.genealogy.cavity_genealogy`
    returns spurious/mislabelled records (not a clean 4 cavities / 2 fusions) and
    the volume-bookkeeping partner has no signal (a survivor's death radius does
    not grow through a merge).  See the xfail test in ``tests/test_genealogy.py``:
    a clean multi-lumen genealogy is the open research this pins."""
    base = _fib_sphere(n_each)
    t = np.linspace(0.0, 1.0, n_frames)
    cS1, cS2 = np.array([-5.0 * radius, 0, 0]), np.array([5.0 * radius, 0, 0])
    rS, rA, rB = 1.3 * radius, 0.7 * radius, 0.95 * radius
    frames, truth, merge_time = [], [], None
    for tk in t:
        cA = np.array([np.interp(tk, [0, 1], [-3.2 * radius, -4.3 * radius]), 0, 0])
        cB = np.array([np.interp(tk, [0, 1], [3.2 * radius, 4.3 * radius]), 0, 0])
        fused = tk > 0.6
        if fused and merge_time is None:
            merge_time = float(tk)
        pts = np.vstack([
            _project_inside(base * rS + cS1, cA, rA),      # survivor 0 (S1)
            _project_inside(base * rA + cA, cS1, rS),      # absorbed 1 (A) -> 0
            _project_inside(base * rS + cS2, cB, rB),      # survivor 2 (S2)
            _project_inside(base * rB + cB, cS2, rS)])     # absorbed 3 (B) -> 2
        labels = np.concatenate([np.full(n_each, k, int) for k in range(4)])
        frames.append((pts, labels))
        truth.append([Void(cS1, rS, True), Void(cA, rA, not fused),
                      Void(cS2, rS, True), Void(cB, rB, not fused)])
    mt = merge_time if merge_time is not None else 1.0
    return Scenario("two_pair_coalescence", t, frames, truth,
                    "Two small cavities fuse into two different survivors (2 fusions).",
                    merges=[(mt, 1, 0), (mt, 3, 2)])   # ground-truth partners


# ---------------------------------------------------------------------------
# Reading the voids out with GUDHI (the trusted oracle)
# ---------------------------------------------------------------------------

def h2_voids(points, min_persistence: float = 0.35):
    """The significant H2 voids of a 3D point set, as ``[(birth, death), ...]``
    sorted by decreasing persistence.  Alpha-filtration values are square-rooted so
    ``death`` is the cavity radius.  ``min_persistence`` drops diagonal noise."""
    import gudhi
    if len(points) < 4:
        return []
    st = gudhi.AlphaComplex(points=np.asarray(points, float)).create_simplex_tree()
    st.compute_persistence(persistence_dim_max=True)
    out = []
    for b, d in st.persistence_intervals_in_dimension(2):
        if d < np.inf:
            b, d = float(np.sqrt(b)), float(np.sqrt(d))
            if d - b > min_persistence:
                out.append((b, d))
    return sorted(out, key=lambda bd: -(bd[1] - bd[0]))


def void_series(scenario: Scenario, min_persistence: float = 0.28):
    """Per-frame significant H2 voids for a scenario: ``(list_of_diagrams,
    computed_counts, true_counts)`` -- the ground-truth check in one call."""
    diagrams = [np.array(h2_voids(pts, min_persistence)).reshape(-1, 2)
                for pts, _ in scenario.frames]
    computed = [len(d) for d in diagrams]
    return diagrams, computed, scenario.true_count()


def void_tracks(scenario: Scenario, *, significance: float = 0.15):
    """The cavity vines of a scenario, computed by the EXACT engine
    (:func:`~vineyards.vineyard.moving_vineyard`) -- no ``track_features``, no
    frame-to-frame matching, and no duration / persistence *stitching* filter.
    Vine identity, including THROUGH the fusion flip, comes from the maintained
    persistence pairing, so the transient voids at a fusing neck never need
    suppressing: they simply have near-zero persistence.

    Requires a fixed-cardinality moving point cloud, i.e. one of the ``*_exact``
    scenarios (the dense resample-and-cull scenarios have no point identity across
    frames and cannot be fed to the exact engine).

    H2 alpha values are square-rooted so ``death`` is the cavity radius.  A vine is
    a real cavity if its persistence EVER exceeds ``significance`` -- a single
    wide-margin significance floor, not a tuned filter: the alpha complex of any
    finite sample also carries near-zero-persistence tetrahedral voids, and here
    the real cavities sit ~50x above them, so any cut in a broad band works.
    Returns the cavity vines (radius-valued), most persistent first."""
    from vineyards.vineyard import moving_vineyard, Vine

    frames = [np.asarray(p, float) for p, _ in scenario.frames]
    n = {len(f) for f in frames}
    if len(n) != 1:
        raise ValueError(
            "void_tracks needs a fixed-cardinality moving point cloud (same number "
            f"of points every frame); got frame sizes {sorted(n)}.  Use a *_exact "
            "scenario (coalescence_exact / ripening_exact / null_exact).")
    npts = next(iter(n))
    if npts > 300:                              # exact engine is ~O(n^3) per frame
        raise ValueError(
            f"void_tracks runs the EXACT vineyard (~O(n^3) per frame); {npts} "
            "points/frame is past its feasible range (> 300).  Use a coarser "
            "*_exact scenario, or track features at scale with track_features.")

    vy, _events = moving_vineyard(frames)
    cavities = []
    for v in vy.by_dimension(2):
        rv = Vine(dim=2, essential=v.essential)
        for p in v.points:
            if p.death < 1e18:                       # skip essential/at-infinity
                rv.append(float(np.sqrt(max(p.birth, 0.0))),
                          float(np.sqrt(max(p.death, 0.0))), p.t)
        if rv.points and max(p.death - p.birth for p in rv.points) > significance:
            cavities.append(rv)
    cavities.sort(key=lambda v: -max(p.death - p.birth for p in v.points))
    return cavities
