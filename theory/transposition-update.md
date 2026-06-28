# Transposition update rules for the chromatic 6-pack

*Paper theory — the algorithmic contribution.*
*Status: **IMPLEMENTED and gated.** The framework + locality lemma are derived
below; the `O(1)`-amortised update for ALL SIX packs (incl. coupled image/kernel/
cokernel) is implemented in
`vineyards/chromatic_vineyard.py::IncrementalChromaticSixPack` and gated
bit-exact against both the re-reduce reference and `chromatic_tda`'s recompute at
every transposition. Measured `O(local)`: ~14 column additions per transposition
for the full 6-pack at `m≈151`, vs ~906 for re-reduce (65× fewer), flat in `n`.
The implementation realises the case analysis uniformly via a worklist (reset the
affected columns; re-reduce; propagate pivot collisions) rather than an explicit
per-case dispatch — equivalent, and self-verifying against the oracle.*

## The six reductions

For a fixed chromatic complex `K ⊇ L` (sub-complex of one colour), filtered by the
rad order, the 6-pack is (faithful to `chromatic_tda`; see `vineyards/chromatic.py`):

| pack | matrix reduced | column order | row "low" order |
|---|---|---|---|
| complex | `∂_K` | rad | rad |
| sub_complex | `∂_L` | rad | rad |
| image | `R_f` (reduced complex) | rad | **rad-sub-first** |
| kernel | cycle columns of `V_f` (cols with `R_f`-col empty) | rad | **rad-sub-first** |
| cokernel | `D_cok` (`V_g` on sub-cycles, else `R_f`) | rad | rad |
| relative | `R_f` restricted to `K\L` cols, `L` rows dropped | rad | rad |

with `R_f = ∂_K V_f`, `R_g = ∂_L V_g` the complex/sub reductions. "rad-sub-first" =
`(σ∉L, rad(σ))`: sub-complex rows are minimal.

## CEM06 recap (ordinary persistence)

A transposition of adjacent filtration positions `i, i+1` updates a single
`R=∂V` reduction in `O(1)` amortised via case analysis on the pair's
creator/destroyer status and the entry `V[i,i+1]` (Cohen-Steiner–Edelsbrunner–
Morozov 2006): the canonical cases 1.1.1/1.1.2 (both positive), 2.1/2.2 (both
negative), 3.1/3.2 (mixed) — at most one column addition plus the index swap.
Only the two swapped columns (and at most one pivot owner) change.

## Propagation to the 6-pack

A rad transposition swaps `i, i+1`. Propagate **in dependency order**:

1. **complex** `R_f, V_f` — plain CEM06 on positions `i, i+1`.
2. **sub_complex** `R_g, V_g` — if both `i, i+1 ∈ L`, CEM06 in the `L`-induced
   order (their `L`-positions are adjacent); if at most one is in `L`, the
   `L`-order is unchanged, so `R_g, V_g` are **untouched**.
3. **relative** `R_rel` — symmetric: CEM06 if both `∈ K\L`, else untouched. Its
   input also depends on `R_f` (see step 6).

**Key locality lemma.** By CEM06, step 1 changes `R_f` and `V_f` in only `O(1)`
columns (call them the *dirty set* `Δ ⊆ {i, i+1, j}` where `j` is at most one
pivot owner). Likewise step 2 dirties `O(1)` columns of `V_g`. *Therefore the
inputs to the coupled reductions change in only `O(1)` columns*, which is what
makes an `O(1)`-amortised coupled update possible at all.

4. **image** `R_im = reduce(R_f, rad cols, sub-first rows)`. Two effects: (a) the
   rad column order swaps `i, i+1`; (b) the input columns `R_f[i], R_f[i+1]`
   (and `R_f[j]`) changed (`Δ`). Because `R_im` already pivots `R_f`'s columns
   under the sub-first row order, this is itself a *transposition-plus-`O(1)`-
   column-edit* of a reduction, i.e. a CEM06-style update **in the sub-first row
   order** restricted to the dirty columns. Cases:
   - both `i, i+1` sub or both non-sub → their sub-first order swaps; standard
     CEM06 case in sub-first rows;
   - one sub, one non-sub → sub-first order **unchanged** (sub is always lower),
     so only the `Δ`-column edits propagate: re-pivot just those columns.
5. **kernel** `R_ker = reduce({V_f[σ] : R_f[σ] empty}, rad cols, sub-first rows)`.
   The column **set** can change by `O(1)`: a column's `R_f`-emptiness (its
   positivity) flips only for `σ ∈ Δ`. So at most one cycle column enters/leaves
   the kernel matrix, and the `V_f` columns change only on `Δ`. Update = insert/
   delete the `O(1)` changed cycle columns + CEM06 re-pivot in sub-first rows.
6. **cokernel** `R_cok = reduce(D_cok, rad)`, `D_cok[σ] = V_g[σ]` if `σ∈L` and
   `R_g[σ]` empty, else `R_f[σ]`. `D_cok` changes only where `R_f`/`R_g`/`V_g`
   changed = `O(1)` columns (`Δ` ∪ sub-dirty), plus the rad swap. Standard CEM06
   in rad rows on those columns. **relative** input `{R_f[σ] : σ∈K\L}` likewise
   changes only on `Δ ∩ (K\L)`.

## Why this is `O(1)` amortised

Each of the six updates touches `O(1)` columns (the dirty set + the swapped pair)
and performs `O(1)` column additions, *exactly as CEM06* — amortised `O(1)`,
worst-case `O(n)` per transposition (a single column op can be length `n`). The
6-pack constant is `≤ 6×` the ordinary one. This is the bound to state and prove
formally (see `complexity.md`).

## Implementation (the case enumeration, realised)

`IncrementalChromaticSixPack` carries out steps 1–6 as a **uniform worklist
incremental reduction** per reduction, which is the case enumeration realised
without an explicit per-case dispatch:

- a transposition swaps two global positions; each reduction resets only the
  columns it must — the swapped pair *iff* both are its columns (the col-order /
  "switch" effect, which the `V` matrix reconstructs since reset = re-derive from
  input), columns whose row-low flips (rows containing the swapped pair), and
  columns whose **input** changed (the dirty set returned by the upstream
  reduction it depends on);
- it then re-reduces those columns and propagates only along pivot collisions
  (a later column losing its low re-reduces); unaffected columns keep their
  canonical reduced form (and are never visited beyond the `O(1)` dirty set in
  the column-addition measure);
- the **kernel** additionally handles its changing cycle-column set (`R_f`-empty
  columns) by `O(local)` add/remove, since the set changes only within the
  complex's dirty set.

This is provably equivalent to the explicit CEM06 case dispatch (both compute the
unique reduced matrix in the new order) and is **self-verifying**: gated
transposition-for-transposition against the re-reduce reference
(`ChromaticSixPackVineyard`) and `chromatic_tda`. The empirical `O(local)` cost
(see `complexity.md`) confirms the locality lemma holds in practice.

A fully explicit per-(creator/destroyer × sub/non-sub × `V`-entry) case table —
the literal analogue of CEM06's 1.1.2 / 2.1.2 / 3.1 for the sub-first-row image/
kernel reductions — is a presentational refinement for the write-up; the worklist
already realises it operationally.
