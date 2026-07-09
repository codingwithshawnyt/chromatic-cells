"""Lumen genealogy from segmented cells -- the measurement pipeline (plumbing).

The target measurement (Le Verge-Serandour & Turlier, PLoS Comp Biol 2021): in a
coarsening hydro-osmotic embryo, lumens disappear by two routes -- coalescence
(merge) and discharge/resorption -- and their RATIO constrains the active
ion-pumping rate.  Counting lumens cannot separate the routes; the merge history
can.  Pipeline:

    cell centroids + radii per timepoint
        -> weighted (regular / Laguerre) moving_vineyard   [H2 voids = interstitial lumens]
        -> cavity_genealogy                                 [fusion vs resorption]
        -> coalescence fraction                             [merges / (merges + resorptions)]

REAL DATA plugs into :func:`lumen_genealogy` unchanged: segmented cell centroids +
radii per timepoint, from the Maitre lab (light-sheet), or from Turlier's public
hydro-osmotic simulator (VirtualEmbryo / hydroosmotic_chain, on GitHub / Zenodo)
whose own event labels give merge ground truth to validate the genealogy before
microscopy.  That data is not shipped here.

HONEST SCOPE.  The coalescence-fraction readout is VALIDATED on the controlled
two-lumen regimes below (a fusion-dominated regime reads 1.0, a resorption-
dominated regime reads 0.0 -- :func:`two_regime_readout`).  On many lumens at once
the per-cavity genealogy is NOT yet clean: the H2 survivor's vine fragments at
each fusion and the merge partner is only provisional (see chromatic_cells.genealogy),
so the multi-lumen :func:`synthetic_embryo` below is an ILLUSTRATION of the target
measurement, not a validated result.  Getting it clean at scale (and validating on
Turlier's labeled events) is the open work this pipeline is built for.

Run:  python examples/blastocyst.py
"""
from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np

from chromatic_cells.synthetic import _fib_sphere, _project_inside, coalescence_exact, ripening_exact
from chromatic_cells.genealogy import cavity_genealogy, coalescence_fraction
from vineyards.regular import assert_no_hidden


def lumen_genealogy(centroids_frames: List[np.ndarray],
                    radii_frames: Optional[List[np.ndarray]] = None,
                    *, significance: float = 0.15, check_hidden: bool = True):
    """Lumen genealogy of a moving cell population.

    ``centroids_frames``: list of ``(n, 3)`` cell-centre arrays, one per timepoint,
    with cell identity by row (a fixed cell set moving over time).  ``radii_frames``
    (optional): matching ``(n,)`` cell radii per timepoint; if given, the weighted
    (regular) alpha complex is used, so the H2 voids are the interstitial lumens.

    Cell radii generally vary per frame; the frame-bound engine here takes a single
    (constant) weight per cell, so we use each cell's mean radius -- adequate while
    radii vary slowly, and flagged honestly (per-frame weights are the exact
    extension).  Returns ``(records, coalescence_fraction)``."""
    weights = None
    if radii_frames is not None:
        r = np.asarray([np.asarray(rr, float) for rr in radii_frames])   # (T, n)
        weights = (r.mean(axis=0)) ** 2                                  # squared radii
        if check_hidden:
            assert_no_hidden(centroids_frames[0], weights)               # scope guard
    records = cavity_genealogy(centroids_frames, weights=weights,
                               significance=significance)
    return records, coalescence_fraction(records)


def two_regime_readout() -> dict:
    """The VALIDATED readout: a fusion-dominated regime reads coalescence fraction
    1.0, a resorption-dominated regime reads 0.0 -- the clean two-lumen ground
    truth that anchors the measurement (cf. Turlier's pumping-rate regimes)."""
    out = {}
    for name, scen in [("fusion-dominated", coalescence_exact()),
                       ("resorption-dominated", ripening_exact())]:
        recs, frac = lumen_genealogy([p for p, _ in scen.frames])
        out[name] = frac
    return out


def synthetic_embryo(n_each: int = 26, n_frames: int = 8, *, R: float = 1.0
                     ) -> Tuple[List[np.ndarray], List[np.ndarray]]:
    """A synthetic coarsening 'embryo' (ILLUSTRATION only -- see module scope): a
    fusing lumen pair AND a resorbing lumen, as a fixed-cardinality moving cell
    cloud with per-cell radii.  Ground truth by construction: one fusion, one
    resorption.  Returns ``(centroids_frames, radii_frames)``."""
    base = _fib_sphere(n_each)
    cell_r = np.full(3 * n_each, 0.12 * R)                # small, comparable cells
    centroids, radii = [], []
    for t in np.linspace(0, 1, n_frames):
        a = np.interp(t, [0, 1], [1.5 * R, 0.35 * R])    # A, B approach and fuse
        cA, cB = np.array([-a, 0.0, 0.0]), np.array([a, 0.0, 0.0])
        pA = _project_inside(base * R + cA, cB, R)
        pB = _project_inside(base * R + cB, cA, R)
        rC = R * (1.0 - 0.94 * t)                         # C shrinks / resorbs
        pC = base * rC + np.array([0.0, 4.5 * R, 0.0])
        centroids.append(np.vstack([pA, pB, pC]))
        radii.append(cell_r.copy())
    return centroids, radii


def main():
    print("VALIDATED two-lumen readout (coalescence fraction):")
    for regime, frac in two_regime_readout().items():
        print(f"  {regime:22s} -> {frac}")

    print("\nSynthetic mixed 'embryo' (ILLUSTRATION -- 1 fusion + 1 resorption built in):")
    centroids, radii = synthetic_embryo()
    recs, frac = lumen_genealogy(centroids, radii)
    fusions = sum(r.fate == "fusion" for r in recs)
    resorptions = sum(r.fate == "resorption" for r in recs)
    print(f"  {len(recs)} cavity records; {fusions} fusion, {resorptions} resorption")
    print(f"  aggregate coalescence fraction = {frac}")
    print("  NOTE: aggregate fraction is indicative; the per-cavity genealogy at "
          ">2 lumens is\n  not yet clean (H2 fragmentation + provisional partner "
          "-- see chromatic_cells.genealogy).")


if __name__ == "__main__":
    main()
