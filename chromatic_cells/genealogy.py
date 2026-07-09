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
    The advance -- the pairing's own geometry carries the partner signal.  Each H2
    vine's DESTROYER is a tetrahedron (the simplex that fills the cavity);
    moving_vineyard exposes it per frame (``Vine.destroyer_frames``), and its
    centroid LOCALISES the cavity in space -- the information the (birth, death)
    diagram discards.  :func:`assign_partners_by_location` assigns an absorbed
    cavity to the surviving cavity whose centre is NEAREST to where it died: the
    signal ``died_radius**3`` lacked (a survivor's radius does not grow through a
    merge, but its LOCATION is exactly where the absorbed cavity vanished).  On the
    adversarial ``two_pair_coalescence`` (two coalescence pairs, all four cavities
    the same size so volume cannot choose) the partners come out correct -- each
    fusion pairs WITHIN its own pair, never across.  The partner SIGNAL is solved.

    What is NOT yet robust -- the fragmentation cleanup.  The survivor's vine still
    fragments at the flip, so getting exactly-clean records needs
    :func:`_stitch_by_location` (re-link fragments at the same place and size) plus
    :func:`_drop_survivor_fragments` (a 'resorption' co-located with a survivor is a
    fragment, since a real resorption leaves no cavity).  These are HEURISTIC and
    parameter-sensitive: a sweep of two_pair_coalescence is clean at n_each>=22 and
    <=6 frames but over- or under-counts at coarser sampling / more frames
    (e.g. n_each=18 nf=7 misses the fusion).  So the location PARTNER is the solved
    part; a fragmentation cleanup robust across all sampling remains open.
    (``stitch_fragments`` / ``assign_partners_by_volume`` remain as the earlier
    radius/volume versions, superseded by the location ones above.)
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
    partner: Optional[int] = None   # merge partner index (fusion only)
    died_location: Optional[np.ndarray] = None      # cavity centre at its death
    location_frames: dict = None                     # {t: cavity centre} over its life

    def location_at(self, t):
        """Cavity centre nearest in time to ``t`` (destroyer-tet centroid)."""
        if not self.location_frames:
            return None
        return self.location_frames[min(self.location_frames, key=lambda s: abs(s - t))]


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


def _classify(vine, last_t, diagonal_frac):
    # "resorbed" = decayed to a small FRACTION of its own peak persistence -- a
    # scale-invariant test, so it works on unit synthetic data and on voxel/micron
    # real data alike (an absolute tolerance would never fire on real-scale radii).
    pers = [p.death - p.birth for p in vine.points]
    if pers[-1] < diagonal_frac * max(pers):
        return "resorption"                       # slid (relatively) to the diagonal
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


def _cavity_items(vineyard, frames, significance):
    """(radius-vine, {t: destroyer-tet centroid}) per significant H2 cavity.

    The centroid of the DESTROYER simplex (the tetrahedron that fills the cavity,
    exposed per frame by moving_vineyard) localises the cavity in space -- the
    pairing geometry the (birth, death) diagram throws away, and the signal that
    resolves the merge partner where volume bookkeeping fails."""
    from vineyards.vineyard import Vine
    items = []
    for v in vineyard.by_dimension(2):
        rv = Vine(dim=2, essential=v.essential)
        loc = {}
        for p in v.points:
            if p.death < 1e18:
                rv.append(float(np.sqrt(max(p.birth, 0.0))),
                          float(np.sqrt(max(p.death, 0.0))), p.t)
                dsimp = v.destroyer_frames.get(p.t)
                if dsimp is not None:
                    loc[p.t] = np.mean([frames[int(round(p.t))][i] for i in dsimp], axis=0)
        if rv.points and max(pt.death - pt.birth for pt in rv.points) > significance:
            items.append((rv, loc))
    return items


def _stitch_by_location(items, loc_tol, radius_tol):
    """Merge cavity fragments split at a fusion flip, matched by BOTH destroyer
    location and death radius across a single-frame gap.  A survivor stays put and
    keeps its size, so its fragments re-link; the absorbed cavity is elsewhere (or
    a different size), so it stays a separate death.  Location resolves what radius
    alone could not -- two equal-radius survivors far apart.  Not a diagram
    heuristic: it uses the pairing's own destroyer geometry."""
    from vineyards.vineyard import Vine
    items = sorted(items, key=lambda it: it[0].points[0].t)
    last_t = max(p.t for rv, _ in items for p in rv.points)
    merged, consumed = [], set()
    for i, (rv, loc) in enumerate(items):
        if i in consumed:
            continue
        pts, L, tail = list(rv.points), dict(loc), rv.points[-1]
        progressed = True
        while progressed and tail.t < last_t - 1e-9:
            progressed = False
            for j, (rv2, loc2) in enumerate(items):
                if j in consumed or j == i or not rv2.points:
                    continue
                head = rv2.points[0]
                if (abs(head.t - (tail.t + 1)) < 1e-9
                        and tail.t in L and head.t in loc2
                        and np.linalg.norm(L[tail.t] - loc2[head.t]) < loc_tol
                        and abs(head.death - tail.death) < radius_tol):
                    pts.extend(rv2.points); L.update(loc2)
                    consumed.add(j); tail = pts[-1]; progressed = True
                    break
        nv = Vine(dim=2, essential=rv.essential)
        for p in pts:
            nv.append(p.birth, p.death, p.t)
        merged.append((nv, L))
    return merged


def assign_partners_by_location(records: List[CavityRecord]) -> None:
    """Merge-partner assignment from the pairing GEOMETRY: an absorbed cavity's
    partner is the surviving cavity whose centre is nearest to where the absorbed
    cavity died (destroyer-tet centroid).  This is the signal ``died_radius**3``
    lacked -- a survivor's radius does not grow through a merge, but its LOCATION
    is exactly where the absorbed cavity vanished.  Falls back to ``None`` (or the
    forced single survivor) when no location is available."""
    for i, r in enumerate(records):
        if r.fate != "fusion":
            continue
        cands = [j for j, s in enumerate(records)
                 if j != i and s.born_t <= r.died_t <= s.died_t + 1e-9
                 and s.fate in ("boundary", "fusion")]
        if not cands:
            r.partner = None
        elif r.died_location is None:
            r.partner = cands[0] if len(cands) == 1 else None
        else:
            def gap(j):
                loc = records[j].location_at(r.died_t)
                return np.inf if loc is None else float(np.linalg.norm(r.died_location - loc))
            r.partner = min(cands, key=gap)


# The exact vineyard is ~O(n^3) per frame; past a few hundred points it is
# infeasible (and OOMs).  cavity_genealogy is the EXACT-engine path (blastocyst
# scale, ~30-100 cells) -- for hundreds+ of points use track_features instead.
MAX_EXACT_POINTS = 300


def cavity_genealogy(points_frames, *, weights=None, significance: float = 0.15,
                     diagonal_frac: float = 0.2, stitch: bool = True,
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
    items = _cavity_items(vineyard, frames, significance)
    if not items:
        return []
    med_r = float(np.median([p.death for rv, _ in items for p in rv.points])) or 1.0
    loc_tol = 2.5 * med_r
    if stitch:
        items = _stitch_by_location(items, loc_tol=loc_tol, radius_tol=0.25 * med_r)
    last_t = max(p.t for rv, _ in items for p in rv.points)
    records = []
    for rv, loc in items:
        pers = [p.death - p.birth for p in rv.points]
        died_t = rv.points[-1].t
        records.append(CavityRecord(
            born_t=rv.points[0].t, died_t=died_t,
            born_radius=rv.points[0].death, died_radius=rv.points[-1].death,
            max_persistence=max(pers), fate=_classify(rv, last_t, diagonal_frac),
            died_location=loc.get(died_t), location_frames=dict(loc)))
    if stitch:
        records = _drop_survivor_fragments(records, loc_tol)
    assign_partners_by_location(records)
    return records


def _drop_survivor_fragments(records, loc_tol):
    """Remove 'resorption' records that are really survivor fragments: a genuine
    resorption leaves NO cavity behind, so one co-located with a surviving cavity
    at its death is a fragmentation artifact of the pairing, not a real death."""
    survivors = [s for s in records if s.fate == "boundary" and s.location_frames]
    kept = []
    for r in records:
        if r.fate == "resorption" and r.died_location is not None:
            # a fragment only if a MORE-persistent surviving cavity sits there: a
            # real resorption of the dominant cavity keeps its (higher) persistence,
            # so we must not drop it just because a smaller boundary fragment is near.
            near = any(s.location_at(r.died_t) is not None
                       and s.max_persistence >= r.max_persistence
                       and np.linalg.norm(r.died_location - s.location_at(r.died_t)) < loc_tol
                       for s in survivors)
            if near:
                continue                              # a fragment of a surviving cavity
        kept.append(r)
    return kept


def coalescence_fraction(records: List[CavityRecord]) -> Optional[float]:
    """Fraction of cavity DEATHS that are fusions rather than resorptions -- the
    readout the biology needs (merges vs discharges).  ``None`` if nothing died."""
    deaths = [r for r in records if r.fate in ("fusion", "resorption")]
    if not deaths:
        return None
    return sum(r.fate == "fusion" for r in deaths) / len(deaths)
