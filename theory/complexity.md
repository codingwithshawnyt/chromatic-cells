# Complexity of the chromatic 6-pack update

*Paper theory — the cost/locality bound Edelsbrunner asked for (meeting 3).*
*Status: target theorem stated and reduced to the CEM06 bound via the locality
lemma in `transposition-update.md`; the current placeholder's cost is
characterised and is empirically measurable now. The `O(1)`-amortised claim is
contingent on the coupled-pack case enumeration (#1).*

## Two regimes

Let `n = |K|` be the number of chromatic simplices.

**(A) Current placeholder (`ChromaticSixPackVineyard`, re-reduce from first
change).** On `advance` to a new weighting, let `p` be the first rad position
whose simplex changes. Each of the six reductions is re-reduced over `order[p:]`,
reusing the unchanged prefix `order[:p]`. Cost: `Θ(Σ_pack (n − p_pack))` column
visits plus the reduction cascade, i.e. `O(n)` per advance, scaled by the six
packs and the per-column fill. For a **single value-crossing transposition** at
position `i`, `p = i`, so the cost is `O(n − i)` per pack — cheap when `i` is high
(crossings near the top of the filtration), `O(n)` worst case (crossing near the
bottom). This is correct (gated bit-exact) but **not** locality-optimal: it
re-visits unchanged columns of the suffix.

**(B) Target (proper CEM06-style update, `transposition-update.md`).** By the
locality lemma, a transposition dirties only `O(1)` columns of `R_f, V_f`
(and `O(1)` of `V_g`), hence `O(1)` input columns of each coupled reduction. Each
pack then performs a CEM06-style update: `O(1)` column additions + the index swap.

## Theorem (target)

> The chromatic 6-pack of a fixed complex is maintained through an adjacent
> transposition in `O(1)` amortised time (worst case `O(n)`, the length of one
> column operation), with a constant at most `6×` the CEM06 ordinary-persistence
> bound.

*Proof reduction.* The six packs update in dependency order (complex → sub →
{image, kernel, cokernel, relative}). complex/sub/relative are ordinary
persistence reductions, each CEM06 (`O(1)` amortised). image/kernel/cokernel each
reduce, by the locality lemma, to a CEM06-style transposition on `O(1)` dirty
columns in their respective row orders (sub-first for image/kernel, rad for
cokernel). Summing the six `O(1)`-amortised updates gives `O(1)` amortised with
constant `≤ 6`. The worst case `O(n)` is inherited from CEM06 (a single column
addition can touch `n` rows). ∎ *(contingent on the coupled-pack case
enumeration; the locality lemma is the substantive content.)*

Contrast: a from-scratch 6-pack recompute is `O(n^3)` dense / `O(n·m_⊥)` sparse
per query; the placeholder (A) is `O(n)` per transposition; the target (B) is
`O(1)` amortised — the CEM06 speedup, carried to the 6-pack.

## Empirical protocol

The repo's instrumentation counts elementary work; for the chromatic dict
reductions, count **column additions** (symmetric differences) per transposition.
Measure, on flip-free intervals of moving two-colour clouds:

1. **Placeholder (A):** column additions vs. crossing position `i`. Expect linear
   in `n − i` (the suffix length) per pack — confirming the suffix cost.
2. **Target (B), once implemented:** column additions per transposition. Expect
   flat in `n` (amortised `O(1)`), confirming the bound, and a histogram peaked at
   small values with rare `O(n)` spikes (the worst-case column op).

The cross-check is bit-exactness of (B) against (A) transposition-for-
transposition (A is the oracle), so the measured (B) cost is for a *verified*
update. (A) is measurable today; (B) follows #1.

## Measured (placeholder A), `theory/measure_placeholder_cost.py`

Eight flip-free intervals of moving two-colour clouds, fine sub-steps, complex
size `m ∈ [111, 131]`, **28 advances, all bit-exact vs recompute**:

- suffix re-reduced (`m − first_divergence`): median `76`, max `108`;
- column additions per advance: median `192`, max `281`;
- **`corr(suffix length, column additions) = 0.92`** — the placeholder's work
  scales with the re-reduced suffix, confirming regime (A) is `O(suffix)` and is
  *not* yet locality-optimal. This is the quantitative motivation for the `O(1)`
  update (B): the median advance re-reduces ~76 columns where the CEM06-style
  update would touch `O(1)`.

(`last_col_adds` on `ChromaticSixPackVineyard` is the counter; column additions =
symmetric differences, the elementary GF(2) work.)
