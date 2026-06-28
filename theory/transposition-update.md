# Transposition update rules for the chromatic 6-pack

*Paper theory — the algorithmic contribution.*
*Status: framework + the key locality reduction derived; the standard packs reduce
to CEM06; the coupled-pack case enumeration is laid out and is the remaining
"careful algebra". The current implementation
(`vineyards/chromatic_vineyard.py::ChromaticSixPackVineyard`) uses a **verified
placeholder**: re-reduce each reduction from the first changed position, which is
correct (gated bit-exact vs recompute) but `O(suffix)`, not the `O(1)`-amortised
CEM06-style update derived here.*

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

## What remains (the careful algebra)

The full **case enumeration** for steps 4–6 — the exact column op per (creator/
destroyer × sub/non-sub × `V`-entry) case, the analogue of CEM06's 1.1.2 / 2.1.2
/ 3.1 for the sub-first-row image/kernel reductions and the `D_cok` cokernel — is
the remaining derivation. The framework above fixes the structure (dependency
order; `O(1)` dirty set; which order each pack pivots in) and reduces each coupled
pack to a CEM06-style transposition on `O(1)` dirty columns. Each case is to be
written and then **gated bit-exact against the verified re-reduce placeholder**
(`ChromaticSixPackVineyard`) — the placeholder is the oracle, so the case analysis
is correct iff it matches it transposition-for-transposition.
