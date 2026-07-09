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

from chromatic_cells.synthetic import coalescence_exact, ripening_exact, null_exact
from chromatic_cells.genealogy import cavity_genealogy, coalescence_fraction


def _genealogy(scenario):
    return cavity_genealogy([p for p, _ in scenario.frames])


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
