"""The resolution floor of centroid TDA: how big must a lumen be, relative to cell
spacing, to register as a distinct persistent H2 void?

This is the question that gates whether nuclear-only (centroid) segmentation --
like BlastoSPIM -- can see microlumen coarsening AT ALL, or only the final
blastocoel.  A microlumen at a cell-cell interface sits at sub-cell scale, on an
edge of the centroid complex; it may fill in at the scale where the two cells
connect and never become a persistent enclosed void.  Only a cavity large relative
to cell spacing, wrapped by a shell of cells, registers.

Setup: cells on a jittered lattice (spacing s, radius s/2 = touching) fill a ball,
with spherical lumen(s) carved out (cells inside a lumen removed).  For each we
find the H2 void whose DESTROYER-tetrahedron centroid is nearest the lumen (that is
the lumen void -- the same pairing geometry the genealogy uses), and compare its
persistence to the packing's interstitial voids (the noise floor).

Finding (both weighted and unweighted alpha, robust across packing disorder):
    lumen radius <~ 1.0 * cell spacing  ->  NOT resolved (drowns in the packing)
    lumen radius >~ 1.5 * cell spacing  ->  resolved, and multiple are separated.

Implication: the blastocoel (large) is visible; interface-scale microlumens are
not; a microlumen becomes visible only once it has grown past ~1.5 cell spacings.
Whether the coarsening microlumens reach that scale -- or need the lumen-resolved
(membrane/fluid) channel instead of nuclei -- is the empirical question this frames.
"""
import numpy as np
import gudhi

S = 2.0                      # cell spacing
RCELL = S / 2.0              # touching cells


def packing_with_lumens(lumens, *, halo_cells=4.0, jitter=0.2, seed=0):
    """Jittered-lattice cell centroids filling a ball, with each ``(centre, radius)``
    lumen carved out (cells inside removed).  ``lumens`` may be one central lumen or
    several at different sites."""
    rng = np.random.default_rng(seed)
    reach = max(np.linalg.norm(c) + r for c, r in lumens) if lumens else 0.0
    R_out = reach + halo_cells * S
    n = int(np.ceil(R_out / S)) + 1
    pts = []
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            for k in range(-n, n + 1):
                p = np.array([i, j, k], float) * S + rng.normal(0, jitter * S, 3)
                if np.linalg.norm(p) >= R_out:
                    continue
                if any(np.linalg.norm(p - c) <= r for c, r in lumens):
                    continue
                pts.append(p)
    return np.array(pts)


def h2_voids(pts, *, weighted=False):
    """H2 voids as ``(persistence_radius, centre, death_radius)`` -- centre is the
    destroyer-tetrahedron centroid (localises the void)."""
    if weighted:
        ac = gudhi.AlphaComplex(points=pts, weights=[RCELL ** 2] * len(pts),
                                precision="exact")
    else:
        ac = gudhi.AlphaComplex(points=pts, precision="exact")
    st = ac.create_simplex_tree()
    st.compute_persistence(persistence_dim_max=True)
    out = []
    for bs, ds in st.persistence_pairs():
        if len(bs) == 3 and len(ds) == 4 and st.filtration(ds) < np.inf:
            br = np.sqrt(max(st.filtration(bs), 0.0))
            dr = np.sqrt(max(st.filtration(ds), 0.0))
            cen = np.mean([ac.get_point(v) for v in ds], axis=0)
            out.append((dr - br, cen, dr))
    return out


def void_at(centre, radius, voids):
    """Persistence of the H2 void located at ``centre`` (within ``radius``), or 0."""
    near = [v[0] for v in voids if np.linalg.norm(v[1] - centre) < radius]
    return max(near, default=0.0)


def packing_noise(voids, sites, *, away=1.5):
    """Max persistence of interstitial voids -- those far from every lumen site."""
    far = [v[0] for v in voids
           if all(np.linalg.norm(v[1] - c) > away * r for c, r in sites)]
    return max(far, default=0.0)


def resolved(centre, radius, voids, sites, *, factor=2.0):
    """A lumen is resolved if the void at its site clears 2x the packing noise."""
    noise = packing_noise(voids, sites)
    p = void_at(centre, radius, voids)
    return p > factor * noise and p > 0.3, p, noise


def main():
    print(f"cell spacing s={S}, cell radius {RCELL} (touching); origin lumen sweep:\n")
    print(f"{'R_lumen/s':>9} | {'lumen persist':>13} | {'noise':>7} | resolved")
    for ratio in [0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]:
        R = ratio * S
        site = [(np.zeros(3), R)]
        voids = h2_voids(packing_with_lumens(site))
        ok, pers, noise = resolved(np.zeros(3), R, voids, site)
        print(f"{ratio:9.2f} | {pers:13.3f} | {noise:7.3f} | {'YES' if ok else 'no'}")

    print("\nthree lumens at separate sites (jitter 0.2), below vs above threshold:")
    for ratio in [0.75, 2.0]:
        R = ratio * S
        d = 3.5 * R                                  # well-separated triangle, no overlap
        sites = [(np.array([d, 0, 0.]), R),
                 (np.array([-0.5 * d, 0.87 * d, 0]), R),
                 (np.array([-0.5 * d, -0.87 * d, 0]), R)]
        voids = h2_voids(packing_with_lumens(sites, halo_cells=4.0))
        nres = sum(resolved(c, r, voids, sites)[0] for c, r in sites)
        print(f"  R_lumen/s={ratio}: {nres}/3 lumens resolved as distinct H2 voids")

    print("\n=> centroid TDA resolves lumens only above ~1.5x cell spacing; "
          "sub-cell-spacing\n   microlumens are invisible to it (need lumen-resolved data).")


if __name__ == "__main__":
    main()
