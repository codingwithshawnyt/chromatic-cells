# chromatic-cells

Application- and paper-specific work for the **kinetic chromatic vineyard** —
the dynamic chromatic 6-pack with vine identity, toward cell segregation
(zebrafish germ-layer sorting; the *Medusa* line × the chromatic 6-pack).

Broadly-applicable engine code lives in the main `vineyards` repo
(`vineyards/chromatic.py`, `vineyards/chromatic_vineyard.py`, tests). This repo
holds what is specific to the paper/application: the theory (theorems, proofs,
complexity), and (to come) the cell-data interface and biological analysis.

Remote: `git@github.com:codingwithshawnyt/chromatic-vineyards.git`.

## theory/

The paper's theoretical spine (drafts — proofs argued, careful steps flagged for
formal completion; the geometry/empirics are verified in the engine's test suite):

- **`flip-handoff.md`** — the chromatic flip-handoff theorem: at a generic
  bistellar flip in the lifted chromatic Delaunay, dying and arriving simplices
  share the empty-stack radius. Generic case proved; degenerate case = SoS
  corollary. (Empirically certified in `vineyards/tests/test_chromatic_handoff.py`.)
- **`transposition-update.md`** — the transposition update rules for the 6-pack
  (the algorithmic contribution): how an adjacent transposition propagates through
  the six interlinked reductions; the locality lemma (`O(1)` dirty columns) and the
  CEM06-style case structure for image/kernel/cokernel. Framework derived; the
  coupled-pack case enumeration is the remaining algebra.
- **`complexity.md`** — the per-transposition cost bound: target `O(1)` amortised
  (≤ 6× CEM06), reduced to the locality lemma; current placeholder is `O(suffix)`;
  empirical protocol via column-addition counts.

## Status (2026-06-28)

SoS-free engine pieces are built and gated in the main repo: the static 6-pack +
simplex-level pairing (bit-exact vs `chromatic_tda`), the dynamic 6-pack
maintained through transpositions (all six packs, gated vs recompute), the
flip-handoff certified for bistellar flips, and **the `O(1)`-amortised
per-transposition update (#1) implemented and gated** (`IncrementalChromaticSixPack`;
~65× fewer column ops than re-reduce, flat in `n`). Open paper work: the
complexity proof write-up (#2, empirically confirmed), the flip-handoff
formalisation (#3, drafted), and 6-pack vine-identity tracking (#4). SoS
(degenerate flips) is the deferred Edelsbrunner item.
