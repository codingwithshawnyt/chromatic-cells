"""Empirical cost of the placeholder 6-pack update (complexity.md, regime A).

Drives ChromaticSixPackVineyard through fine sub-steps of a flip-free interval of
a moving two-colour cloud, recording the elementary work (column additions,
`last_col_adds`) per advance against the suffix length re-reduced (n - first
divergence position).  Also re-verifies bit-exactness vs recompute at every step
(so the measured cost is for a *correct* update).

Run:  python chromatic-cells/theory/measure_placeholder_cost.py
"""

from __future__ import annotations

import numpy as np

from vineyards.chromatic import PACKS, chromatic_filtration, six_pack
from vineyards.chromatic_vineyard import ChromaticSixPackVineyard
from vineyards.comparison import diagrams_equal


def _arr(bars):
    return np.array(bars, float) if bars else np.zeros((0, 3))


def main():
    rows = []
    mismatches = 0
    for seed in range(8):
        rng = np.random.default_rng(seed)
        n = 12
        p0 = rng.random((n, 2))
        vel = rng.standard_normal((n, 2)) * 0.4
        labels = list(rng.integers(0, 2, n))
        labels[0], labels[1] = 0, 1

        base = 0.05 * (seed % 6)
        ts = [base + 0.002 * k for k in range(8)]
        cfgs = [chromatic_filtration(p0 + t * vel, labels) for t in ts]
        if len({frozenset(w) for w, _ in cfgs}) != 1:        # not flip-free
            continue
        m = len(cfgs[0][0])
        w0, sub = cfgs[0]
        vy = ChromaticSixPackVineyard(w0, sub)
        prev_order = list(vy.order)
        for w, _ in cfgs[1:]:
            new_order = sorted(w, key=lambda s: (w[s], len(s), s))
            gi = next((i for i, (a, b) in enumerate(zip(prev_order, new_order))
                       if a != b), len(prev_order))
            vy.advance(w)
            # correctness: maintained == recompute
            got, ref = vy.diagrams(), six_pack(w, sub)
            if not all(diagrams_equal(_arr(got[p]), _arr(ref[p])) for p in PACKS):
                mismatches += 1
            rows.append((m, m - gi, vy.last_col_adds))
            prev_order = new_order

    a = np.array(rows, dtype=float)
    print(f"advances measured: {len(rows)}   bit-exact vs recompute: "
          f"{len(rows) - mismatches}/{len(rows)}")
    print(f"complex size m: {int(a[:,0].min())}..{int(a[:,0].max())}")
    print(f"suffix re-reduced (m - first_divergence): "
          f"median {np.median(a[:,1]):.0f}, max {int(a[:,1].max())}")
    print(f"column additions / advance: median {np.median(a[:,2]):.0f}, "
          f"max {int(a[:,2].max())}")
    # correlation of work with suffix length (confirms O(suffix) regime)
    if a[:, 1].std() > 0:
        corr = np.corrcoef(a[:, 1], a[:, 2])[0, 1]
        print(f"corr(suffix length, column additions) = {corr:.2f}  "
              f"(positive => work scales with the re-reduced suffix, as derived)")


if __name__ == "__main__":
    main()
