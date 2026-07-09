"""Cavity genealogy from the exact vineyard: the fusion-vs-resorption readout.

The genealogy extractor reads, off the maintained pairing, HOW each cavity dies:
a fusion is an off-diagonal vine death (only possible through a flip -- invisible
to a fixed-complex vineyard or a matching heuristic), a resorption slides to the
diagonal.  These check the robust part on controlled ground truth: coalescence is
all fusion, ripening all resorption, null neither.  The merge PARTNER is validated
only where it is forced (two cavities); multi-cavity disambiguation is provisional
(see chromatic_cells.genealogy) and deliberately not asserted here.
"""

from __future__ import annotations

import numpy as np
import pytest

from chromatic_cells.synthetic import (
    coalescence, coalescence_exact, ripening_exact, null_exact, two_pair_coalescence,
)
from chromatic_cells.genealogy import cavity_genealogy, coalescence_fraction, stitch_fragments


def _genealogy(scenario):
    return cavity_genealogy([p for p, _ in scenario.frames])


# ------------------------------------------------------------------
# Robustness: the exact-engine path must fail gracefully, not crash/OOM
# ------------------------------------------------------------------

def test_stitch_fragments_handles_empty():
    # no significant cavity -> nothing to stitch (was: max() of empty -> ValueError)
    assert stitch_fragments([]) == []


def test_cavity_genealogy_empty_when_no_significant_cavity():
    # a shell too sparse to enclose an H2 void yields zero cavities, cleanly
    # (this is the path that reached the empty stitch_fragments crash)
    assert cavity_genealogy([p for p, _ in coalescence_exact(n_each=6, n_frames=4).frames]) == []


def test_cavity_genealogy_rejects_dense_scenario():
    # the dense demo scenarios cull points per frame (varying cardinality); the
    # exact engine cannot take them -- fail loudly BEFORE running, not OOM.
    with pytest.raises(ValueError, match="fixed-cardinality"):
        cavity_genealogy([p for p, _ in coalescence(4, 60).frames])


def test_cavity_genealogy_rejects_oversized_cloud():
    # a fixed-cardinality but too-large cloud would OOM the O(n^3) engine
    big = [np.random.default_rng(0).random((400, 3)) for _ in range(2)]
    with pytest.raises(ValueError, match="feasible range"):
        cavity_genealogy(big)


def test_fusion_and_resorption_are_distinguished():
    coal = _genealogy(coalescence_exact())
    rip = _genealogy(ripening_exact())
    nul = _genealogy(null_exact())
    # coalescence: the one death is a fusion; ripening: the one death is a
    # resorption; null: nothing dies.
    assert coalescence_fraction(coal) == 1.0
    assert coalescence_fraction(rip) == 0.0
    assert coalescence_fraction(nul) is None


def test_fusion_death_is_off_diagonal():
    # the theorem-shaped fact: a fusion cavity terminates BEFORE the final frame
    # while still at substantial persistence (an off-diagonal death), which a
    # fixed-complex vineyard cannot produce.
    recs = _genealogy(coalescence_exact())
    last = max(r.died_t for r in recs)
    fusions = [r for r in recs if r.fate == "fusion"]
    assert fusions
    assert all(r.died_t < last - 1e-9 and r.max_persistence > 0.2 for r in fusions)


def test_resorption_reaches_the_diagonal():
    recs = _genealogy(ripening_exact())
    resorbed = [r for r in recs if r.fate == "resorption"]
    assert resorbed
    # a resorbed cavity's death radius is near zero (it shrank away)
    assert all(r.died_radius < 0.2 for r in resorbed)


def test_merge_partner_is_forced_for_two_cavities():
    # with a single survivor the partner is not a guess -- the absorbed cavity's
    # partner must be the surviving (boundary) cavity.
    recs = _genealogy(coalescence_exact())
    survivor = next(i for i, r in enumerate(recs) if r.fate == "boundary")
    for r in recs:
        if r.fate == "fusion":
            assert r.partner == survivor


@pytest.mark.slow
@pytest.mark.xfail(strict=True, reason=(
    "OPEN: the multi-cavity genealogy is not clean.  With >2 cavities the H2 "
    "survivor vines fragment at each fusion flip; cavity_genealogy returns "
    "spurious/mislabelled records (not 4 cavities / 2 fusions) and the "
    "volume-bookkeeping partner has no signal (a survivor's death radius does not "
    "grow through a merge).  This is the SoCG-relevant research still open; the "
    "assertion below is the target a real fix must meet."))
def test_partner_choice_survives_the_adversarial_multi_cavity_case():
    # two small cavities fuse into two DIFFERENT large survivors (ground-truth
    # partners: absorbed 1 -> survivor 0, absorbed 3 -> survivor 2).  A correct
    # genealogy would report exactly four cavities, two of them fusions, each
    # partnered to a distinct large survivor.  It does not (yet).
    recs = cavity_genealogy([p for p, _ in two_pair_coalescence().frames])
    fusions = [r for r in recs if r.fate == "fusion"]
    survivors = [i for i, r in enumerate(recs) if r.fate == "boundary"]
    assert len(recs) == 4                                  # not fragmented
    assert len(fusions) == 2                               # both small cavities fuse
    assert len({r.partner for r in fusions}) == 2          # to two DIFFERENT survivors
    assert all(r.partner in survivors for r in fusions)    # each to a real survivor
