"""The resolution floor of centroid TDA (the gating fact for the biology).

A lumen registers as a distinct persistent H2 void of the CELL-CENTROID alpha
complex only once it is large relative to the cell spacing -- below ~1.25-1.5x it
drowns in the packing's interstitial voids.  So nuclear-only segmentation (e.g.
BlastoSPIM) can see the blastocoel but NOT sub-cell-spacing microlumens; measuring
microlumen coarsening this way needs them grown past the floor (or lumen-resolved
data).  These characterization tests keep that honest limit visible.  Uses gudhi
directly (a static geometric fact -- no vineyard needed).
"""
from __future__ import annotations

import numpy as np
import pytest

gudhi = pytest.importorskip("gudhi")

S = 2.0                      # cell spacing
RCELL = S / 2.0             # touching cells


def _packing(lumens, halo=4.0, jitter=0.2, seed=0):
    rng = np.random.default_rng(seed)
    reach = max(np.linalg.norm(c) + r for c, r in lumens)
    R_out = reach + halo * S
    n = int(np.ceil(R_out / S)) + 1
    pts = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                p = np.array([i, j, k], float) * S + rng.normal(0, jitter * S, 3)
                if np.linalg.norm(p) >= R_out:
                    continue
                if any(np.linalg.norm(p - c) <= r for c, r in lumens):
                    continue
                pts.append(p)
    return np.array(pts)


def _h2_voids(pts):
    ac = gudhi.AlphaComplex(points=pts, precision="exact")
    st = ac.create_simplex_tree()
    st.compute_persistence(persistence_dim_max=True)
    out = []
    for bs, ds in st.persistence_pairs():
        if len(bs) == 3 and len(ds) == 4 and st.filtration(ds) < np.inf:
            pers = np.sqrt(max(st.filtration(ds), 0.0)) - np.sqrt(max(st.filtration(bs), 0.0))
            cen = np.mean([ac.get_point(v) for v in ds], axis=0)
            out.append((pers, cen))
    return out


def _resolved(centre, radius, voids, sites, factor=2.0):
    at = max((v[0] for v in voids if np.linalg.norm(v[1] - centre) < radius), default=0.0)
    noise = max((v[0] for v in voids
                 if all(np.linalg.norm(v[1] - c) > 1.5 * r for c, r in sites)), default=0.0)
    return at > factor * noise and at > 0.3


def test_subcell_microlumen_is_not_resolved():
    # a lumen the size of the cell spacing (interface-scale microlumen) drowns
    R = 0.75 * S
    site = [(np.zeros(3), R)]
    assert not _resolved(np.zeros(3), R, _h2_voids(_packing(site)), site)


def test_lumen_above_threshold_is_resolved():
    # a lumen ~2x the cell spacing (blastocoel / grown microlumen) is clear
    R = 2.0 * S
    site = [(np.zeros(3), R)]
    assert _resolved(np.zeros(3), R, _h2_voids(_packing(site)), site)


def test_multiple_supra_threshold_lumens_are_separated():
    # three well-separated 2x-spacing lumens -> three distinct H2 voids
    R = 2.0 * S
    d = 3.5 * R
    sites = [(np.array([d, 0, 0.]), R),
             (np.array([-0.5 * d, 0.87 * d, 0]), R),
             (np.array([-0.5 * d, -0.87 * d, 0]), R)]
    voids = _h2_voids(_packing(sites))
    assert sum(_resolved(c, r, voids, sites) for c, r in sites) == 3
