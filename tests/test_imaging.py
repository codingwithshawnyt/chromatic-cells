"""Imaging ingest: label masks -> centroids + radii -> corresponded frames -> genealogy.

The pipeline that lets the engine run on real light-sheet embryo masks (BlastoSPIM).
These check the mask extraction is correct, that correspondence handles relabeling
and refuses a changing cell count (the frame-bound engine needs a fixed cell set),
and that a full masks -> genealogy run reads a shrinking shell as a resorption.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy")

from chromatic_cells.imaging import (
    mask_centroids_radii, correspond_frames, masks_to_genealogy,
)


def test_mask_centroids_radii_matches_known_geometry():
    vol = np.zeros((30, 30, 30), int)
    vol[5:9, 5:9, 5:9] = 1            # 64-voxel cube, centre (6.5, 6.5, 6.5)
    vol[20:26, 20:26, 20:26] = 2      # 216-voxel cube, centre (22.5, 22.5, 22.5)
    ids, cen, rad = mask_centroids_radii(vol)
    assert list(ids) == [1, 2]
    assert np.allclose(cen, [[6.5] * 3, [22.5] * 3])
    # equivalent-sphere radius of 64 and 216 voxels
    assert np.allclose(rad, [(3 / (4 * np.pi) * v) ** (1 / 3) for v in (64, 216)])


def test_correspond_frames_nearest_recovers_identity_through_relabel():
    # two cells, relabeled and moved slightly between frames -> nearest-centroid
    # correspondence must restore frame 0's row order.
    f0 = (np.array([1, 2]), np.array([[0, 0, 0.], [5, 0, 0]]), np.array([1., 1]))
    f1 = (np.array([9, 7]), np.array([[5.2, 0, 0.], [0.1, 0, 0]]), np.array([1., 1]))
    cf, _rf = correspond_frames([f0, f1])
    assert cf[1][0, 0] < 1.0 and cf[1][1, 0] > 4.0     # reordered to match frame 0


def test_correspond_frames_rejects_changing_cell_count():
    f0 = (np.array([1, 2]), np.zeros((2, 3)), np.ones(2))
    f1 = (np.array([1, 2, 3]), np.zeros((3, 3)), np.ones(3))   # a cell divided/appeared
    with pytest.raises(ValueError, match="cell count varies"):
        correspond_frames([f0, f1])


def test_correspond_frames_with_tracks_uses_global_ids():
    # same two cells, different local label order, tracks give global identity
    f0 = (np.array([10, 20]), np.array([[0, 0, 0.], [5, 0, 0]]), np.array([1., 2]))
    f1 = (np.array([20, 10]), np.array([[6, 0, 0.], [0.5, 0, 0]]), np.array([2., 1]))
    tracks = [np.array([100, 200]), np.array([200, 100])]      # global ids per local cell
    cf, rf = correspond_frames([f0, f1], tracks=tracks)
    assert cf[0][0, 0] == 0 and cf[1][0, 0] == 0.5             # global 100 first, both frames
    assert rf[0][0] == 1 and rf[1][0] == 1


@pytest.mark.slow
def test_masks_to_genealogy_reads_a_shrinking_shell_as_resorption():
    # cells on a shrinking spherical shell, rasterised as labeled balls (the
    # BlastoSPIM mask format); the extracted centroids enclose one H2 cavity that
    # resorbs -> coalescence fraction 0.0.
    def fib(n):
        i = np.arange(n) + 0.5
        phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
        return np.column_stack([np.cos(th) * np.sin(phi),
                                np.sin(th) * np.sin(phi), np.cos(phi)])
    G, base, centre = 48, fib(22), np.array([24.] * 3)
    zz, yy, xx = np.mgrid[0:G, 0:G, 0:G]
    vols = []
    for t in np.linspace(0, 1, 5):
        rv = np.interp(t, [0, 1], [13.0, 2.0])       # shell shrinks -> cavity resorbs
        vol = np.zeros((G, G, G), np.int32)
        for k, d in enumerate(base):
            c = centre + rv * d
            vol[(zz - c[0]) ** 2 + (yy - c[1]) ** 2 + (xx - c[2]) ** 2 <= 1.6 ** 2] = k + 1
        vols.append(vol)
    records, fraction = masks_to_genealogy(vols)
    assert records                                             # found a cavity
    assert fraction == 0.0                                     # its death is a resorption
    assert any(r.fate == "resorption" for r in records)
