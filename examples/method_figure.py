"""The methods figure: recompute+match vs the exact vineyard, on one fusion.

Both engines run on the SAME coalescence data (two H2 cavities fusing into one).

  * Left  -- recompute + match (``track_features``): per-frame GUDHI diagrams
    matched across frames.  At the fusing neck it throws off short spurious
    cavity tracks, so getting "two clean cavities" needs a hand-tuned duration /
    persistence filter -- i.e. the heuristic suppresses exactly the merge-event
    artifacts.
  * Right -- the exact vineyard (``moving_vineyard``): vine identity comes from
    the maintained pairing, so there is nothing to suppress -- two cavities, and
    the absorbed one dies OFF the diagonal (a fusion), the death a fixed-complex
    vineyard / matching cannot produce.

This is the argument, not an illustration of it.  It is a METHODS figure (for
reviewers), separate from the plain-language scene demos in void_demo.py.

Run:  python examples/method_figure.py       # writes examples/media/exact_vs_matching.png
"""
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm
from matplotlib.lines import Line2D

from chromatic_cells.synthetic import coalescence_exact, h2_voids, void_tracks
from vineyards.feature_tracking import track_features
from chromatic_cells.genealogy import stitch_fragments

_ARIAL = next((p for p in ("/System/Library/Fonts/Supplemental/Arial.ttf",
                           "/Library/Fonts/Arial.ttf",
                           "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf")
               if os.path.exists(p)), None)
if _ARIAL:
    _fm.fontManager.addfont(_ARIAL)
plt.rcParams.update({"font.family": "sans-serif",
                     "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.size": 12})
FG, GRID, SPINE = "#1b2733", "#e4e9f0", "#c2c9d2"
BLUE, RED, MUTE = "#4477AA", "#EE6677", "#b8c0cc"


def render_exact_vs_matching(filename, *, n_each=30, n_frames=8, min_frames=3):
    scen = coalescence_exact(n_each, n_frames)
    frames = [p for p, _ in scen.frames]
    idx = list(range(n_frames))                       # frame-index time axis for both

    # --- LEFT: recompute + match, no hand-tuned duration filter ---
    diagrams = [np.array([[2.0, b, d] for b, d in h2_voids(p, 0.15)]).reshape(-1, 3)
                for p in frames]
    tracks = track_features(frames, times=idx, dims=(2,),
                            min_persistence=0.15, diagrams=diagrams)

    # --- RIGHT: the exact vineyard (fragments stitched to cavity identity) ---
    vines = stitch_fragments(void_tracks(scen))        # clean cavity vines
    last_t = max(p.t for v in vines for p in v.points)

    def pers(pt):
        return pt.death - pt.birth

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12.5, 5.2), dpi=100, sharey=True)
    fig.patch.set_facecolor("white")

    # LEFT panel: persistence over time; short tracks (the neck artifacts) in red
    n_spur = 0
    for tr in tracks:
        ts = [p[0] for p in tr.points]
        ps = [p[2] - p[1] for p in tr.points]
        spurious = len(tr.points) < min_frames
        n_spur += spurious
        axL.plot(ts, ps, "--" if spurious else "-",
                 color=RED if spurious else BLUE,
                 lw=1.8 if spurious else 3.0, alpha=0.95,
                 marker="o", ms=5 if spurious else 6, zorder=4 if spurious else 2)
    axL.set_title("Recompute + match (track_features)\n"
                  f"{len(tracks)} tracks: {n_spur} spurious at the fusing neck "
                  "— needs a hand-tuned filter", fontsize=12, color=FG)
    axL.legend(handles=[
        Line2D([0], [0], color=BLUE, lw=3, label="cavity track"),
        Line2D([0], [0], color=RED, lw=1.8, ls="--", marker="o",
               label=f"spurious (< {min_frames} frames)")],
        loc="upper right", frameon=False, fontsize=10)

    # RIGHT panel: two clean vines; the absorbed one ends OFF the diagonal
    for v in vines:
        ts = [p.t for p in v.points]
        ps = [pers(p) for p in v.points]
        axR.plot(ts, ps, "-", color=BLUE, lw=3.0, marker="o", ms=6, zorder=2)
        ends_early = v.points[-1].t < last_t - 1e-9
        if ends_early and pers(v.points[-1]) > 0.15:            # off-diagonal death
            axR.scatter([v.points[-1].t], [pers(v.points[-1])], marker="X",
                        color=RED, s=150, zorder=5, edgecolor="white", lw=1.2)
            axR.annotate("fusion: off-diagonal death\n(absorbed into the survivor)",
                         (v.points[-1].t, pers(v.points[-1])),
                         textcoords="offset points", xytext=(6, 12), fontsize=9.5,
                         color=FG, arrowprops=dict(arrowstyle="->", color=FG, lw=1))
    axR.set_title("Exact vineyard (moving_vineyard)\n"
                  f"{len(vines)} cavities, no filter — identity from the pairing",
                  fontsize=12, color=FG)
    axR.legend(handles=[
        Line2D([0], [0], color=BLUE, lw=3, label="cavity vine (exact)"),
        Line2D([0], [0], color=RED, lw=0, marker="X", ms=10,
               label="fusion death")],
        loc="upper right", frameon=False, fontsize=10)

    for ax in (axL, axR):
        ax.set_xlabel("frame", color=FG)
        ax.set_xlim(-0.3, n_frames - 0.7)
        ax.set_ylim(0, 0.62)
        ax.axhline(0, color=MUTE, lw=1)                 # the diagonal (persistence 0)
        for sp in ax.spines.values():
            sp.set_color(SPINE)
        ax.tick_params(colors=FG)
        ax.grid(True, color=GRID, lw=0.8)
    axL.set_ylabel("cavity persistence  (death − birth, radius)", color=FG)
    fig.suptitle("Same fusion, two engines: matching needs a filter to hide the "
                 "merge-event artifacts; the exact engine has none",
                 fontsize=13.5, fontweight="bold", y=1.02, color=FG)

    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(filename, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    return filename


if __name__ == "__main__":
    out = Path(__file__).parent / "media"
    print("wrote", render_exact_vs_matching(out / "exact_vs_matching.png"))
