"""Run the lumen-genealogy pipeline on BlastoSPIM (or any labeled-mask) embryo data.

BlastoSPIM (Nunley, Posfai, Shvartsman, Brown et al.;
blastospim.flatironinstitute.org) is a public light-sheet dataset of preimplantation
mouse embryos with ground-truth 3D nucleus segmentation -- 573 annotated images
across 31 embryos, one series of 89 consecutive timepoints, 8-cell to >100-cell.
The ground truth is voxel LABEL masks (each cell an integer region), which is
exactly what chromatic_cells.imaging consumes: centroids + radii come straight from
the masks, no GPU segmentation needed.

    python examples/blastospim.py /path/to/series   # a folder of per-frame label volumes
    python examples/blastospim.py                   # no data: a synthetic-mask self-test

READ THIS BEFORE POINTING IT AT REAL DATA (see chromatic_cells.imaging for the full
version):
  * Correspondence.  moving_vineyard needs the SAME cell in row i every frame.  The
    vineyard supplies vine identity, not cell correspondence.  Use BlastoSPIM's
    lineage tracking (pass as ``tracks``), or -- over a division-free window -- the
    built-in nearest-centroid fallback.  Cell DIVISION changes the count and is
    rejected with a clear message (restrict the window).
  * Pick a LATE, CONSECUTIVE series.  Early embryo (8-32 cell) may have no
    persistent H2 cavity; coarsening is ~32-64+ cell.  The download mixes single
    snapshots and series -- you want the consecutive one (motion needs adjacent
    frames), i.e. the 89-timepoint embryo.
  * What this shows.  The ENGINE running on real embryo geometry (real noise, real
    cell counts, real rearrangements) -- NOT yet the pumping-rate result, which
    needs lumen-resolved, ion-perturbed imaging (Turlier / Maitre).  Nuclei masks
    give the moving point cloud; the cavities are inferred from it, which is the
    whole method.

Requires scipy (`pip install .[imaging]`); TIFF input additionally needs tifffile.
"""
import sys
from pathlib import Path

import numpy as np

from chromatic_cells.imaging import masks_to_genealogy


def load_volumes(folder):
    """Load a time series of 3D label volumes from ``folder`` (sorted by name).
    Supports ``.npy`` natively; ``.tif``/``.tiff`` via tifffile (optional)."""
    folder = Path(folder)
    paths = sorted(p for p in folder.iterdir()
                   if p.suffix.lower() in (".npy", ".tif", ".tiff"))
    if not paths:
        raise FileNotFoundError(f"no .npy/.tif label volumes in {folder}")
    vols = []
    for p in paths:
        if p.suffix.lower() == ".npy":
            vols.append(np.load(p))
        else:
            try:
                import tifffile
            except ImportError as e:
                raise ImportError("reading .tif masks needs tifffile "
                                  "(`pip install tifffile`)") from e
            vols.append(tifffile.imread(p))
    print(f"loaded {len(vols)} frames from {folder} (shape {vols[0].shape})")
    return vols


def _synthetic_shell_masks(n_cells=22, n_frames=5, grid=48, cell_r=1.6):
    """Self-test data: cells on a shrinking spherical shell, rasterised as labeled
    balls in a voxel grid -- the same LABEL-MASK format as BlastoSPIM.  The
    extracted centroids enclose one H2 cavity that resorbs, so the pipeline should
    read coalescence fraction 0.0.  (Grid/cell sizes kept generous so the cells
    stay separate voxel regions as the shell shrinks -- crowding would merge labels
    and change the cell count, which the correspondence step rightly rejects.)"""
    i = np.arange(n_cells) + 0.5
    phi = np.arccos(1 - 2 * i / n_cells); th = np.pi * (1 + 5 ** 0.5) * i
    base = np.column_stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)])
    centre = np.array([grid / 2.0] * 3)
    zz, yy, xx = np.mgrid[0:grid, 0:grid, 0:grid]
    vols = []
    for t in np.linspace(0, 1, n_frames):
        radius_vox = np.interp(t, [0, 1], [13.0, 2.0])       # shell shrinks -> resorbs
        vol = np.zeros((grid, grid, grid), np.int32)
        for k, d in enumerate(base):
            c = centre + radius_vox * d
            ball = (zz - c[0]) ** 2 + (yy - c[1]) ** 2 + (xx - c[2]) ** 2 <= cell_r ** 2
            vol[ball] = k + 1                                 # consistent label per cell
        vols.append(vol)
    return vols


def run(volumes, *, voxel_size=(1.0, 1.0, 1.0), tracks=None):
    records, fraction = masks_to_genealogy(volumes, voxel_size=voxel_size, tracks=tracks)
    print(f"\n{len(records)} cavity record(s); coalescence fraction = {fraction}")
    for i, r in enumerate(records):
        print(f"  cavity {i}: fate={r.fate:11s} born_r={r.born_radius:.2f} "
              f"died_r={r.died_radius:.2f} t=[{r.born_t:.0f}..{r.died_t:.0f}]")
    print("\nReminder: this is the ENGINE on real geometry, not the pumping-rate "
          "result (which needs lumen-resolved imaging).  Multi-cavity partner "
          "genealogy is still open -- see chromatic_cells.genealogy.")
    return records, fraction


def main():
    if len(sys.argv) > 1:
        run(load_volumes(sys.argv[1]))
    else:
        print("no data path given -- running the synthetic label-mask self-test "
              "(cells on a shrinking shell; expect resorption, fraction 0.0).")
        run(_synthetic_shell_masks())


if __name__ == "__main__":
    main()
