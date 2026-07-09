"""Imaging ingest: labeled segmentation masks -> the lumen-genealogy pipeline.

For real light-sheet embryo data -- e.g. BlastoSPIM (Nunley, Posfai, Shvartsman,
Brown et al., blastospim.flatironinstitute.org): 3D nucleus LABEL masks per
timepoint (each voxel an integer cell id, 0 = background).  We extract cell
CENTROIDS (per-label center of mass) and RADII (voxel count -> equivalent-sphere
radius); the radii give the weighted (regular) alpha its weights, so the H2 voids
are the interstitial cavities (candidate lumens).  The CLI is examples/blastospim.py.

THE ONE REAL GOTCHA -- correspondence.  ``moving_vineyard`` is a MOVING point
cloud: row i at frame t must be the SAME cell at t+1 (it re-triangulates
consecutive frames and diffs by vertex index).  It supplies VINE (feature)
identity, NOT cell correspondence -- that is an INPUT.  So consistent cell
identity across frames is required:

  * pass ``tracks`` -- per-frame global cell ids (e.g. from BlastoSPIM's lineage
    tracking); the same set of ids must appear every frame; or
  * ``tracks=None`` -- nearest-centroid correspondence, valid only over a window
    with NO cell division/death (constant count).  Cell division changes the count
    and breaks the fixed-cardinality requirement of the frame-bound engine;
    :func:`correspond_frames` raises with guidance in that case (restrict to a
    division-free window, or supply tracks).

Honest expectations: early embryo (8-32 cell) may show FEW/NO persistent H2
cavities -- coarsening is later (~32-64+ cell), so use late consecutive series.
This validates the ENGINE on real embryo geometry; the pumping-rate biology needs
lumen-resolved, ion-perturbed imaging (Turlier / Maitre), not nuclei masks.
Requires scipy (``pip install .[imaging]``).
"""
from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

import numpy as np

from chromatic_cells.genealogy import cavity_genealogy, coalescence_fraction


def mask_centroids_radii(labels, *, voxel_size=(1.0, 1.0, 1.0), min_voxels: int = 1):
    """A 3D integer label volume (0 = background) -> ``(ids, centroids, radii)``.

    ``centroids`` (n, 3) are per-label centres of mass scaled by ``voxel_size``;
    ``radii`` (n,) are equivalent-sphere radii of each label's volume (voxel count
    x voxel volume).  Labels with fewer than ``min_voxels`` voxels are dropped."""
    from scipy import ndimage
    labels = np.asarray(labels)
    counts_all = np.bincount(labels.ravel())
    ids = np.array([i for i in range(1, len(counts_all)) if counts_all[i] >= min_voxels],
                   dtype=int)
    if len(ids) == 0:
        return ids, np.zeros((0, 3)), np.zeros(0)
    vs = np.asarray(voxel_size, float)
    coms = ndimage.center_of_mass(np.ones_like(labels, dtype=np.uint8), labels, list(ids))
    centroids = np.asarray(coms, float) * vs
    counts = counts_all[ids].astype(float)
    radii = (3.0 / (4.0 * np.pi) * counts * float(np.prod(vs))) ** (1.0 / 3.0)
    return ids, centroids, radii


def estimate_anisotropy(labels):
    """Estimate the z:xy voxel-spacing ratio of a 3D label volume from cell shapes.

    Light-sheet microscopy (e.g. BlastoSPIM) has anisotropic voxels -- z coarser
    than xy -- and the wrong voxel_size distorts the H2 cavities badly.  When the
    imaging metadata is unavailable, this recovers the ratio geometrically: cells
    are roughly round, so a cell's (xy-extent / z-extent) in voxels approximates
    how much larger the z spacing is.  Returns ``r`` for ``voxel_size=(r, 1, 1)``
    (the median over cells with a resolvable z-extent)."""
    labels = np.asarray(labels)
    ratios = []
    for i in np.unique(labels):
        if i == 0:
            continue
        zz, yy, xx = np.where(labels == i)
        zext = zz.max() - zz.min() + 1
        if zext > 2:
            xyext = 0.5 * ((yy.max() - yy.min() + 1) + (xx.max() - xx.min() + 1))
            ratios.append(xyext / zext)
    return float(np.median(ratios)) if ratios else 1.0


def _reorder(prev_c, cur_c):
    """Hungarian assignment of current centroids to previous ones -> a permutation
    ``perm`` with ``cur_c[perm][k]`` matched to ``prev_c[k]``."""
    from scipy.optimize import linear_sum_assignment
    d = np.linalg.norm(prev_c[:, None, :] - cur_c[None, :, :], axis=2)
    _rows, cols = linear_sum_assignment(d)          # rows are 0..n-1 in order
    return cols


def correspond_frames(per_frame, *, tracks: Optional[Sequence] = None):
    """Turn per-frame ``(ids, centroids, radii)`` into fixed-cardinality corresponded
    ``(centroids_frames, radii_frames)`` -- row i is the SAME cell in every frame,
    which is what :func:`~vineyards.vineyard.moving_vineyard` requires.

    ``tracks[t][k]`` gives the global id of the k-th cell of ``per_frame[t]``; the
    same set of global ids must appear in every frame.  With ``tracks=None`` a
    nearest-centroid correspondence is used, valid only when the cell count is
    constant (no division/death) -- otherwise a ValueError explains the fix."""
    counts = [len(c) for _ids, c, _r in per_frame]
    if tracks is not None:
        common = set(np.asarray(tracks[0]).tolist())
        for gt in tracks[1:]:
            if set(np.asarray(gt).tolist()) != common:
                raise ValueError(
                    "tracks: the set of global cell ids changes between frames "
                    "(cell division / death / appearance).  The frame-bound engine "
                    "needs a fixed cell set -- restrict to a division-free window.")
        order = sorted(common)
        cf, rf = [], []
        for (_ids, c, r), gt in zip(per_frame, tracks):
            pos = {int(g): k for k, g in enumerate(np.asarray(gt))}
            idx = [pos[g] for g in order]
            cf.append(c[idx]); rf.append(r[idx])
        return cf, rf

    if len(set(counts)) != 1:
        raise ValueError(
            f"cell count varies across frames {counts} (cell division / death), so "
            "nearest-centroid correspondence is undefined.  Supply tracks (global "
            "cell ids per frame), or restrict to a division-free window.")
    cf = [per_frame[0][1]]
    rf = [per_frame[0][2]]
    for t in range(1, len(per_frame)):
        perm = _reorder(cf[-1], per_frame[t][1])
        cf.append(per_frame[t][1][perm]); rf.append(per_frame[t][2][perm])
    return cf, rf


def lumen_genealogy(centroids_frames, radii_frames=None, *, significance: float = 0.15,
                    check_hidden: bool = True):
    """Lumen genealogy of a moving, ALREADY-CORRESPONDED cell population.

    ``centroids_frames``: list of ``(n, 3)`` cell-centre arrays with cell identity
    by row.  ``radii_frames`` (optional): matching ``(n,)`` cell radii per frame;
    if given, the weighted (regular) alpha is used so the H2 voids are the
    interstitial lumens.  Cell radii vary per frame; the frame-bound engine takes a
    single (mean) weight per cell -- adequate while radii vary slowly, flagged
    honestly.  Returns ``(records, coalescence_fraction)``."""
    from vineyards.regular import assert_no_hidden
    weights = None
    if radii_frames is not None:
        r = np.asarray([np.asarray(rr, float) for rr in radii_frames])   # (T, n)
        weights = (r.mean(axis=0)) ** 2                                  # squared radii
        if check_hidden:
            assert_no_hidden(centroids_frames[0], weights)               # scope guard
    records = cavity_genealogy(centroids_frames, weights=weights, significance=significance)
    return records, coalescence_fraction(records)


def masks_to_genealogy(volumes, *, voxel_size=(1.0, 1.0, 1.0), tracks=None,
                       min_voxels: int = 1, significance: float = 0.15,
                       normalize: bool = True):
    """End to end: a time series of 3D label volumes -> ``(records, coalescence
    fraction)``.  Extracts centroids + radii per frame, corresponds cells across
    frames (see :func:`correspond_frames`), then runs the weighted lumen genealogy.

    With ``normalize`` (default), centroids + radii are rescaled to units of the
    mean cell radius before the genealogy, so the persistence ``significance`` floor
    is scale-appropriate -- raw voxel/micron radii would otherwise sit far from the
    genealogy's calibration.  (Resorption is judged by a relative fraction of each
    cavity's own peak, so it is scale-invariant regardless.)  ``significance`` may
    still want per-dataset tuning on real data."""
    per_frame = [mask_centroids_radii(v, voxel_size=voxel_size, min_voxels=min_voxels)
                 for v in volumes]
    centroids_frames, radii_frames = correspond_frames(per_frame, tracks=tracks)
    if normalize:
        scale = float(np.mean([r.mean() for r in radii_frames if len(r)]))
        if scale > 0:
            centroids_frames = [c / scale for c in centroids_frames]
            radii_frames = [r / scale for r in radii_frames]
    return lumen_genealogy(centroids_frames, radii_frames, significance=significance)
