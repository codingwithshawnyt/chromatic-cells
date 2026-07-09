"""Synthetic void scenarios: the engine recovers known ground truth.

The whole point of synthetic data is that we control the answer, so these tests
are the honest, public check that H2 persistence reports the cavities we built:
a shell's H2 death is its radius; coalescence and ripening BOTH drop the count
4 -> 2 (so counting alone cannot tell fusion from resorption); ripening is
mass-conserving; the null model holds the count flat.  Requires gudhi.
"""

from __future__ import annotations

import numpy as np
import pytest

from chromatic_cells.synthetic import (
    union_boundary, h2_voids, coalescence, ripening, null_model,
    void_series, void_tracks,
    coalescence_exact, ripening_exact, null_exact,
)


# ------------------------------------------------------------------
# The core mapping: a shell is an H2 void whose death is its radius
# ------------------------------------------------------------------

@pytest.mark.parametrize("R", [1.0, 1.5, 2.0])
def test_shell_is_one_void_with_death_equal_radius(R):
    rng = np.random.default_rng(0)
    pts, _ = union_boundary([(0, 0, 0)], [R], 260, rng)
    voids = h2_voids(pts, min_persistence=0.35)
    assert len(voids) == 1                       # exactly one significant cavity
    birth, death = voids[0]
    assert abs(death - R) < 0.05 * R             # death value == cavity radius
    assert birth < death


def test_two_separated_shells_are_two_voids():
    rng = np.random.default_rng(0)
    pts, _ = union_boundary([(-2.5, 0, 0), (2.5, 0, 0)], [1.0, 1.0], 240, rng)
    assert len(h2_voids(pts, 0.35)) == 2


# ------------------------------------------------------------------
# Scenario ground truth is recovered by the computed H2 count
# ------------------------------------------------------------------

def test_null_count_matches_ground_truth_every_frame():
    scen = null_model(24, 300)
    _, computed, true = void_series(scen)
    assert true == [4] * len(true)               # built as four fixed cavities
    assert computed == true                       # ... and recovered exactly


def test_coalescence_count_drops_4_to_2_matching_ground_truth():
    scen = coalescence(40, 320)
    _, computed, true = void_series(scen)
    assert true[0] == 4 and true[-1] == 2         # two pairs each fuse: 4 -> 2
    assert computed[0] == 4 and computed[-1] == 2
    # matches ground truth except within ~2 frames of the fusion (the resolution
    # scale: "shells touch" vs "the complex resolves one cavity")
    assert sum(int(c != t) for c, t in zip(computed, true)) <= 2
    assert len(scen.merges) == 2                   # two fusion events recorded


def test_coalescence_and_ripening_have_the_same_count_but_differ_in_merges():
    """The necessity of the merge history: both scenes drop the H2 count 4 -> 2,
    so counting cannot tell them apart -- but coalescence's deaths are MERGES
    (recorded partners) and ripening's are resorptions (no merge)."""
    coal = coalescence(40, 320)
    rip = ripening(40, 320)
    assert coal.true_count()[-1] == rip.true_count()[-1] == 2
    assert coal.true_count()[0] == rip.true_count()[0] == 4
    assert len(coal.merges) == 2 and len(rip.merges) == 0


def test_ripening_kills_small_voids_null_does_not():
    rip = ripening(24, 320)
    nul = null_model(24, 320)
    _, comp_r, true_r = void_series(rip, 0.42)
    _, comp_n, true_n = void_series(nul, 0.42)
    assert true_r[0] == 4 and true_r[-1] < 4      # ripening loses voids
    assert comp_r[0] == 4 and comp_r[-1] < 4
    assert comp_n[0] == comp_n[-1] == 4           # null keeps all four


# ------------------------------------------------------------------
# Tracking from the EXACT engine (moving_vineyard) -- no matching heuristic
#
# void_tracks reads cavities off the exact vineyard, not track_features, and
# applies only a wide-margin persistence significance floor (no duration /
# stitching filter).  These check it recovers the cavity history, and -- the
# scientific point -- that a fusion is an OFF-diagonal vine death while a
# resorption slides to the DIAGONAL, the distinction a fixed-complex vineyard or
# a matching heuristic cannot make.
# ------------------------------------------------------------------

def _sig_count_per_frame(tracks, n_frames, thr=0.15):
    """Number of cavity vines significant (persistence > thr) at each frame."""
    return [sum(1 for v in tracks for p in v.points
                if abs(p.t - float(k)) < 1e-6 and (p.death - p.birth) > thr)
            for k in range(n_frames)]


def test_exact_engine_recovers_per_frame_cavity_count():
    # fusion and resorption both drop 2 -> 1, null stays 2 -- matched exactly by
    # the exact engine with no matching and no stitching filter.
    for build in (coalescence_exact, ripening_exact, null_exact):
        scen = build()
        tracks = void_tracks(scen)
        assert _sig_count_per_frame(tracks, len(scen.times)) == scen.true_count()


def test_exact_fusion_is_off_diagonal_resorption_is_at_diagonal():
    coal = void_tracks(coalescence_exact())
    rip = void_tracks(ripening_exact())
    last = max(p.t for v in coal for p in v.points)
    # coalescence: a cavity vine terminates BEFORE the final frame while still at
    # substantial persistence -- an off-diagonal fusion death.
    off_diagonal = [v for v in coal
                    if v.points[-1].t < last - 1e-6
                    and (v.points[-1].death - v.points[-1].birth) > 0.2]
    assert off_diagonal, "expected an off-diagonal (fusion) vine death"
    # ripening: the resorbing cavity slides to the diagonal (persistence -> ~0),
    # and no cavity terminates early (both vines span every frame).
    assert min(v.points[-1].death - v.points[-1].birth for v in rip) < 0.1
    assert all(abs(v.points[-1].t - max(p.t for p in v.points)) < 1e-6 for v in rip)


def test_significance_floor_is_wide_margin_not_tuned():
    # the real cavities sit ~50x above the alpha-complex noise, so any cut in a
    # broad band gives the same cavities -- the floor is not a delicate filter.
    scen = null_exact()
    assert (len(void_tracks(scen, significance=0.05))
            == len(void_tracks(scen, significance=0.15))
            == len(void_tracks(scen, significance=0.30)) == 2)


def test_ripening_is_mass_conserving():
    """Total cavity volume is (very nearly) constant -- the physical content of
    Ostwald ripening, and the check that a shrinking void is resorbed, not lost."""
    scen = ripening(40, 320)
    vol = [sum((4 / 3) * np.pi * v.radius ** 3 for v in voids if v.alive)
           for voids in scen.truth]
    assert abs(vol[-1] - vol[0]) < 0.02 * vol[0]
