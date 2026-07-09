"""Cavity genealogy from the exact vineyard: birth, death, and how a cavity dies.

Consumes the vineyard of a moving point cloud (:func:`moving_vineyard`) and emits,
per significant H2 cavity, when it was born and died, its death radius, and -- the
point -- HOW it died:

    * ``resorption`` -- the cavity's vine slides to the DIAGONAL (persistence -> 0):
      it shrank away.
    * ``fusion``     -- the vine terminates OFF the diagonal, at substantial
      persistence: it was absorbed at a combinatorial (flip) event.
    * ``boundary``   -- still alive (substantial) at the last frame.

Why this needs the exact engine (the theorem-shaped fact).
    On a FIXED complex a vine can only die at the diagonal (a birth/death pair
    collapses).  An OFF-diagonal vine death is possible only through a change of
    the complex -- a flip.  So fusion is not merely blurred but INVISIBLE to every
    fixed-complex vineyard and to every frame-to-frame matching heuristic, by
    construction.  Reading ``fate == 'fusion'`` off the pairing is exactly the
    measurement a matching pipeline cannot make.

What is exact here, and what is a stated definition (the honest boundary).
    Exact, straight from the maintained pairing: each vine's birth, death, death
    radius, and whether the death is off-diagonal (fusion) or at the diagonal
    (resorption).  NOT canonical for H2: cavity IDENTITY THROUGH a fusion.  For H0
    "which merged into which" is the elder rule; for H2 there is no canonical
    "which vine continues" -- at a fusion flip the pairing reorganises and the
    survivor's vine fragments (an artifact this module must resolve, not hide).
    So two things are DEFINITIONS, provided and to be justified, not theorems:

      1. Cavity continuation through a fusion -- :func:`stitch_fragments` re-links a
         vine that ends mid-motion to one that begins at the next frame with a
         near-equal death radius (the same cavity, split by the flip).  A narrow,
         geometric re-link across a single-frame gap -- NOT the broad diagram
         matching the project indicts.
      2. The merge PARTNER (who absorbed whom) -- :func:`assign_partners_by_volume`,
         the surviving cavity whose volume gains ~ the dying cavity's volume
         (mass bookkeeping).  With two cavities this is FORCED (one survivor) and
         correct.  With three or more it is a hypothesis, and it does NOT yet hold:
         on the adversarial ``two_pair_coalescence`` scenario (two small cavities
         fusing into two different survivors) the genealogy does not even return
         clean records -- the H2 survivor vines fragment at each fusion flip, so
         the count and fates are wrong and the volume partner has no signal (a
         survivor's death radius does not grow through a merge).  Pinned by the
         xfail test ``test_partner_choice_survives_the_adversarial_multi_cavity_case``.
         A clean multi-lumen genealogy -- resolving the fragmentation from the
         pairing (not a radius heuristic) and a partner observable that carries
         signal -- is the open research; do NOT report multi-cavity as settled.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np


@dataclass
class CavityRecord:
    """One cavity's life, read off the exact vineyard (radii, not squared)."""
    born_t: float
    died_t: float
    born_radius: float
    died_radius: float
    max_persistence: float
    fate: str                       # 'boundary' | 'resorption' | 'fusion'
    partner: Optional[int] = None   # provisional merge partner index (fusion only)


def _radius_cavity_vines(vineyard, significance):
    """Significant H2 vines with alpha values square-rooted to radii."""
    from vineyards.vineyard import Vine
    out = []
    for v in vineyard.by_dimension(2):
        rv = Vine(dim=2, essential=v.essential)
        for p in v.points:
            if p.death < 1e18:
                rv.append(float(np.sqrt(max(p.birth, 0.0))),
                          float(np.sqrt(max(p.death, 0.0))), p.t)
        if rv.points and max(p.death - p.birth for p in rv.points) > significance:
            out.append(rv)
    return out


def stitch_fragments(vines, *, radius_tol: float = 0.08, gap: float = 1.0):
    """Re-link cavity vines split by the pairing reorganisation at a fusion flip.

    DEFINITION (not a theorem): if one vine ends at frame ``t`` (before the last
    frame) and another begins at ``t + gap`` with a near-equal death radius
    (within ``radius_tol``), they are taken to be the same cavity, fragmented by
    the flip -- the exact pairing offers no canonical continuation for H2, so this
    supplies one.  It is a narrow single-frame geometric re-link, NOT diagram
    matching across the whole trajectory.  Returns a new list of merged vines."""
    from vineyards.vineyard import Vine
    if not vines:                                    # no significant cavity to stitch
        return []
    vines = sorted(vines, key=lambda v: v.points[0].t)
    last_t = max(p.t for v in vines for p in v.points)
    merged, consumed = [], set()
    for i, v in enumerate(vines):
        if i in consumed:
            continue
        chain = list(v.points)
        end = chain[-1]
        # extend the chain across flip gaps while a near-equal-radius vine follows
        progressed = True
        while progressed and end.t < last_t - 1e-9:
            progressed = False
            for j, w in enumerate(vines):
                if j in consumed or j == i or not w.points:
                    continue
                start = w.points[0]
                if (abs(start.t - (end.t + gap)) < 1e-9
                        and abs(start.death - end.death) < radius_tol):
                    chain.extend(w.points)
                    consumed.add(j)
                    end = chain[-1]
                    progressed = True
                    break
        nv = Vine(dim=2, essential=v.essential)
        for p in chain:
            nv.append(p.birth, p.death, p.t)
        merged.append(nv)
    return merged


def _classify(vine, last_t, diagonal_tol):
    pers = [p.death - p.birth for p in vine.points]
    if pers[-1] < diagonal_tol:
        return "resorption"                       # slid to the diagonal
    if vine.points[-1].t >= last_t - 1e-9:
        return "boundary"                         # still alive at the end
    return "fusion"                               # off-diagonal death at a flip


def assign_partners_by_volume(records: List[CavityRecord]) -> None:
    """PROVISIONAL merge-partner assignment by mass bookkeeping (candidate (c)).

    For each ``fusion`` cavity, the partner is the cavity alive at its death whose
    volume (radius**3) most plausibly absorbs the dying cavity's volume.  Forced
    when only one survivor exists (two-cavity coalescence); a hypothesis, not a
    result, with three or more cavities -- do not report as settled."""
    for i, r in enumerate(records):
        if r.fate != "fusion":
            continue
        cands = [j for j, s in enumerate(records)
                 if j != i and s.born_t <= r.died_t <= s.died_t + 1e-9
                 and s.fate in ("boundary", "fusion")]
        if not cands:
            r.partner = None
        elif len(cands) == 1:
            r.partner = cands[0]                  # forced
        else:
            donated = r.died_radius ** 3
            r.partner = min(
                cands,
                key=lambda j: abs((records[j].died_radius ** 3
                                   - records[j].born_radius ** 3) - donated))


# The exact vineyard is ~O(n^3) per frame; past a few hundred points it is
# infeasible (and OOMs).  cavity_genealogy is the EXACT-engine path (blastocyst
# scale, ~30-100 cells) -- for hundreds+ of points use track_features instead.
MAX_EXACT_POINTS = 300


def cavity_genealogy(points_frames, *, weights=None, significance: float = 0.15,
                     diagonal_tol: float = 0.12, stitch: bool = True,
                     radius_tol: float = 0.08,
                     max_points: int = MAX_EXACT_POINTS) -> List[CavityRecord]:
    """The cavity genealogy of a fixed-cardinality moving point cloud.

    Runs the exact vineyard, extracts the significant H2 cavities, classifies each
    death (resorption / fusion / boundary), and -- provisionally -- assigns fusion
    partners by volume.  With ``stitch=True`` the survivor's fragments (see
    :func:`stitch_fragments`) are re-linked first, so a fused-into survivor reads
    as one ``boundary`` cavity spanning the motion, not several pieces.

    ``points_frames`` is a list of ``(n, d)`` arrays (a moving point cloud;
    e.g. ``[p for p, _ in scenario.frames]`` for a ``*_exact`` scenario, or cell
    centroids per timepoint).  ``weights`` (squared radii) switches to the weighted
    engine -- check :func:`vineyards.regular.assert_no_hidden` first."""
    from vineyards.vineyard import moving_vineyard
    frames = [np.asarray(f, float) for f in points_frames]
    # Guard the exact engine before it runs (mirrors void_tracks): it needs a
    # fixed-cardinality moving cloud, and it is ~O(n^3) per frame so a large cloud
    # OOMs.  Fail loudly with guidance rather than crash.
    sizes = {len(f) for f in frames}
    if len(sizes) != 1:
        raise ValueError(
            "cavity_genealogy needs a fixed-cardinality moving point cloud (same "
            f"number of points every frame); got frame sizes {sorted(sizes)}.  The "
            "dense void scenarios cull points per frame -- use a *_exact scenario, "
            "or cell centroids with stable identity.")
    n = sizes.pop()
    if n > max_points:
        raise ValueError(
            f"cavity_genealogy runs the EXACT vineyard (~O(n^3) per frame); {n} "
            f"points/frame is past its feasible range (> {max_points}; blastocyst "
            "scale is ~30-100 cells).  Coarsen the sampling, or track prominent "
            "features at scale with vineyards.track_features.")
    vineyard, _events = moving_vineyard(frames, weights=weights)
    vines = _radius_cavity_vines(vineyard, significance)
    if stitch:
        vines = stitch_fragments(vines, radius_tol=radius_tol)
    if not vines:
        return []
    last_t = max(p.t for v in vines for p in v.points)
    records = []
    for v in vines:
        pers = [p.death - p.birth for p in v.points]
        records.append(CavityRecord(
            born_t=v.points[0].t, died_t=v.points[-1].t,
            born_radius=v.points[0].death, died_radius=v.points[-1].death,
            max_persistence=max(pers), fate=_classify(v, last_t, diagonal_tol)))
    assign_partners_by_volume(records)
    return records


def coalescence_fraction(records: List[CavityRecord]) -> Optional[float]:
    """Fraction of cavity DEATHS that are fusions rather than resorptions -- the
    readout the biology needs (merges vs discharges).  ``None`` if nothing died."""
    deaths = [r for r in records if r.fate in ("fusion", "resorption")]
    if not deaths:
        return None
    return sum(r.fate == "fusion" for r in deaths) / len(deaths)
