"""chromatic-cells: application- and paper-specific work on top of the general
``vineyards`` engine.

The broadly-applicable engine (moving/kinetic vineyards, weighted/regular alpha,
the chromatic 6-pack) lives in the ``vineyards`` package.  This package holds the
CELL / VOID application built on it: synthetic void scenarios with known ground
truth (coalescence / Ostwald ripening / null), and the cavity genealogy that reads
fusion vs resorption off the exact pairing (toward lumen coarsening / pumping-rate
inference).
"""

from chromatic_cells.synthetic import (
    Scenario,
    Void,
    coalescence,
    ripening,
    null_model,
    coalescence_exact,
    ripening_exact,
    null_exact,
    two_pair_coalescence,
    h2_voids,
    void_series,
    void_tracks,
)
from chromatic_cells.genealogy import (
    CavityRecord,
    cavity_genealogy,
    coalescence_fraction,
    stitch_fragments,
    assign_partners_by_volume,
    assign_partners_by_location,
)
from chromatic_cells.imaging import (
    mask_centroids_radii,
    correspond_frames,
    lumen_genealogy,
    masks_to_genealogy,
)

__all__ = [
    "Scenario",
    "Void",
    "coalescence",
    "ripening",
    "null_model",
    "coalescence_exact",
    "ripening_exact",
    "null_exact",
    "two_pair_coalescence",
    "h2_voids",
    "void_series",
    "void_tracks",
    "CavityRecord",
    "cavity_genealogy",
    "coalescence_fraction",
    "stitch_fragments",
    "assign_partners_by_volume",
    "assign_partners_by_location",
    "mask_centroids_radii",
    "correspond_frames",
    "lumen_genealogy",
    "masks_to_genealogy",
]
