"""Watch the cavities change -- and watch the count change with them.

The whole point, in one picture per scene: on the LEFT you see the cavities
(lumens) as translucent blobs that move, merge or shrink; on the RIGHT the method's
output -- the number of enclosed voids it finds (H2) -- draws itself as a curve and
a big running number that steps down at the exact frame a blob merges or vanishes.

Three scenes, each an animated GIF, all with the SAME layout so the outputs are
visibly different:

  * coalescence -- cavities fuse in pairs        (voids 4 -> 2)
  * ripening    -- Ostwald ripening: small ones shrink to nothing, the big one grows (4 -> 2)
  * null        -- the control: cavities hold their size                (stays 4)

coalescence / ripening / null is the contrast Ishihara asked for: you can SEE the
data change, and you can SEE the topology respond.  Because the data is synthetic
we know the answer, so the count on the right is the method's output plotted
against a ground truth we built (they agree; see tests/test_synthetic.py).

Run:  python examples/void_demo.py            # writes examples/media/*.gif + genus.png
      python examples/void_demo.py ripening   # just one scene
"""
import os
import sys
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm
from PIL import Image
import pyvista as pv

from chromatic_cells.synthetic import coalescence, ripening, null_model, void_series

# Publication figure font: Arial (the Nature/Cell figure standard), with clean
# sans-serif fallbacks so the demo still renders elsewhere.
_ARIAL_CANDIDATES = [
    "/System/Library/Fonts/Supplemental/Arial.ttf",   # macOS
    "/Library/Fonts/Arial.ttf",
    "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",  # Linux
]
ARIAL_TTF = next((p for p in _ARIAL_CANDIDATES if os.path.exists(p)), None)
if ARIAL_TTF:
    _fm.fontManager.addfont(ARIAL_TTF)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 13,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "mathtext.default": "regular",
})


def _vtk_arial(actor):
    """Render a PyVista text actor in Arial (from the .ttf, to match matplotlib)."""
    if actor is None or ARIAL_TTF is None:
        return
    try:
        tp = actor.GetTextProperty()
        tp.SetFontFamily(4)          # VTK_FONT_FILE
        tp.SetFontFile(ARIAL_TTF)
    except Exception:
        pass


# Colour-blind-safe (Tol) palette, indexed by ground-truth void.
PALETTE = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"]
BG = "white"
FG = "#1b2733"        # dark text/axes on white
ACCENT = "#EE6677"    # the count line
GRID = "#e4e9f0"
SPINE = "#c2c9d2"


# ---------------------------------------------------------------------------
# Left panel: the cavities you can see (translucent blobs + their sampled skin)
# ---------------------------------------------------------------------------

def _render_3d(pl, scenario, k, dynamic):
    """Draw the cavities of frame ``k`` as translucent solid blobs (the lumen you
    can see) with a faint dusting of the boundary points the method actually gets.
    Swaps only these actors; the title actor (with its font) persists across
    frames."""
    for a in dynamic:
        pl.remove_actor(a, render=False)
    dynamic.clear()
    pts, labels = scenario.frames[k]
    present = set(int(l) for l in np.unique(labels))
    # the solid, visible cavities -- the star of the panel
    for lab, v in enumerate(scenario.truth[k]):
        if lab in present and v.radius > 0.25:
            dynamic.append(pl.add_mesh(
                pv.Sphere(radius=v.radius, center=v.center,
                          theta_resolution=48, phi_resolution=48),
                color=PALETTE[lab % len(PALETTE)], opacity=0.42,
                smooth_shading=True, specular=0.3, diffuse=0.9))
    # a faint dusting of the boundary sample (what the method is actually handed)
    for lab in np.unique(labels):
        sel = pts[labels == lab]
        dynamic.append(pl.add_mesh(pv.PolyData(sel),
                                   color=PALETTE[int(lab) % len(PALETTE)],
                                   render_points_as_spheres=True, point_size=3.5,
                                   opacity=0.55))
    pl.camera.azimuth = 20 + 34 * scenario.times[k]    # slow orbit for depth
    pl.camera.elevation = 15
    return np.asarray(pl.screenshot(return_img=True))[..., :3]


# ---------------------------------------------------------------------------
# Right panel: the method's output -- the number of voids, drawn over time
# ---------------------------------------------------------------------------

def _render_count(fig, ax, big, scenario, counts, k, ymax):
    """The output: a big running count of voids, and the count curve drawing itself
    up to the current frame (so it steps down exactly when a blob merges/vanishes)."""
    times = np.asarray(scenario.times)
    now = times[k]
    n_now = counts[k]

    # the big number
    big.clear()
    big.axis("off")
    big.text(0.5, 0.5, str(n_now), ha="center", va="center", fontsize=64,
             fontweight="bold", color=ACCENT)
    big.text(0.5, 0.03, "enclosed voids found now  (H2)", ha="center", va="bottom",
             fontsize=12.5, color=FG)

    # the curve so far
    ax.clear()
    ax.plot(times[: k + 1], counts[: k + 1], drawstyle="steps-post",
            color=ACCENT, lw=3.4, solid_capstyle="round")
    ax.scatter([now], [n_now], color=ACCENT, s=70, zorder=5)
    ax.axvline(now, color=FG, lw=1, alpha=0.35)
    ax.set_xlim(times[0], times[-1])
    ax.set_ylim(0, ymax)
    ax.set_yticks(range(0, ymax + 1))
    ax.set_xlabel("time", color=FG)
    ax.set_ylabel("number of voids  (H2)", color=FG)
    ax.set_title("What the method outputs", color=FG, fontsize=13, fontweight="bold")
    ax.text(0.5, -0.29, scenario.caption, transform=ax.transAxes, ha="center",
            color=FG, fontsize=11, style="italic")
    for spine in ax.spines.values():
        spine.set_color(SPINE)
    ax.tick_params(colors=FG)
    ax.grid(True, color=GRID, lw=0.8)
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()


# ---------------------------------------------------------------------------
# Compose the movie: [ cavities you can see | count that responds ]
# ---------------------------------------------------------------------------

def render_void_movie(scenario, filename, *, side=700, fps=11):
    _, counts, _ = void_series(scenario)
    ymax = max(max(counts), max(scenario.true_count())) + 1

    pl = pv.Plotter(off_screen=True, window_size=(side, side))
    pl.set_background("white")
    allpts = np.vstack([p for p, _ in scenario.frames])
    pl.camera_position = "iso"
    pl.reset_camera(bounds=[allpts[:, 0].min(), allpts[:, 0].max(),
                            allpts[:, 1].min(), allpts[:, 1].max(),
                            allpts[:, 2].min(), allpts[:, 2].max()])
    pl.camera.zoom(1.1)
    # add the 3D title ONCE and fix its font once (re-adding it every frame made
    # VTK's font cache apply the .ttf on only one frame).
    _vtk_arial(pl.add_text("The cavities  (blob = the void, dots = the sample)",
                           position="upper_edge", font_size=11, color=FG, font="arial"))

    # right panel: a big number on top, the count curve below
    fig = plt.figure(figsize=(side / 100, side / 100), dpi=100)
    fig.patch.set_facecolor(BG)
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.6],
                          left=0.13, right=0.95, top=0.95, bottom=0.20, hspace=0.28)
    big = fig.add_subplot(gs[0]); big.set_facecolor(BG)
    ax = fig.add_subplot(gs[1]); ax.set_facecolor(BG)

    frames = []
    dynamic = []
    for k in range(len(scenario.frames)):
        left = _render_3d(pl, scenario, k, dynamic)
        right = _render_count(fig, ax, big, scenario, counts, k, ymax)
        h = min(left.shape[0], right.shape[0])
        frames.append(Image.fromarray(np.hstack([left[:h], right[:h]])))
    pl.close()
    plt.close(fig)

    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    pal = [f.convert("P", palette=Image.ADAPTIVE, colors=128) for f in frames]
    pal[0].save(filename, save_all=True, append_images=pal[1:],
                duration=int(1000 / fps), loop=0, optimize=True, disposal=2)
    return filename


# ---------------------------------------------------------------------------
# Genus: the OTHER observable (H1), Ishihara's second sanity check
# ---------------------------------------------------------------------------

def render_genus_figure(filename, *, side=560, n=340):
    """A sphere and a torus both enclose ONE void, so the void count (H2) cannot
    tell them apart -- but the torus has a handle (genus 1), which H1 sees and H2
    does not.  This is Ishihara's second criterion: in one condition the lumen
    count goes DOWN while the genus goes UP, and the two are different letters
    (H2 and H1)."""
    from chromatic_cells.synthetic import _fib_sphere

    def counts(pts, mp=0.3):
        import gudhi
        st = gudhi.AlphaComplex(points=pts).create_simplex_tree()
        st.compute_persistence(persistence_dim_max=True)
        return {k: sum(1 for b, d in
                       [(np.sqrt(x), np.sqrt(y))
                        for x, y in st.persistence_intervals_in_dimension(k) if y < np.inf]
                       if d - b > mp) for k in (1, 2)}

    rng = np.random.default_rng(0)
    sphere = _fib_sphere(n) * 1.3
    u, v = rng.uniform(0, 2 * np.pi, 2 * n), rng.uniform(0, 2 * np.pi, 2 * n)
    Rt, rt = 1.6, 0.9
    torus = np.column_stack([(Rt + rt * np.cos(v)) * np.cos(u),
                             (Rt + rt * np.cos(v)) * np.sin(u), rt * np.sin(v)])
    cs, ct = counts(sphere), counts(torus)

    pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(2 * side, side))
    pl.set_background("white")
    for col, (pts, color, name, c) in enumerate(
            [(sphere, "#4477AA", "Sphere", cs), (torus, "#EE6677", "Torus", ct)]):
        pl.subplot(0, col)
        pl.add_mesh(pv.PolyData(pts), color=color, render_points_as_spheres=True,
                    point_size=6)
        genus = c[1] // 2
        _vtk_arial(pl.add_text(f"{name}   (genus {genus})\nvoids (H2): {c[2]}     "
                               f"handles (H1): {c[1]}", position="upper_left",
                               font_size=12, color="#1b2733", font="arial"))
        pl.camera_position = [(1.7, -1.7, 4.2), (0, 0, 0), (0, 1, 0)]  # near top-down: torus hole shows
    top = np.asarray(pl.screenshot(return_img=True))[..., :3]
    pl.close()

    W = top.shape[1]
    fig, ax = plt.subplots(figsize=(W / 100, 1.5), dpi=100)
    fig.patch.set_facecolor("white"); ax.axis("off")
    ax.text(0.5, 0.72, "One void either way -- but a different genus", ha="center",
            fontsize=15, fontweight="bold", color="#1b2733")
    ax.text(0.5, 0.28,
            "Both shapes enclose ONE void, so counting voids (H2 = 1) cannot tell them apart.  "
            "The torus has a handle (genus 1):\nH1 = 2, while the sphere has H1 = 0.  So "
            "\"lumens down, genus up\" is two different measurements -- H2 and H1 -- and the method reports both.",
            ha="center", va="center", fontsize=10.5, color="#1b2733")
    fig.canvas.draw()
    banner = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    plt.close(fig)
    h = min(top.shape[1], banner.shape[1])
    out = np.vstack([top[:, :h], banner[:, :h]])
    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(out).save(filename)
    return filename


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    builders = {"coalescence": coalescence, "ripening": ripening, "null": null_model}
    out = Path(__file__).parent / "media"
    for name, build in builders.items():
        if which not in ("all", name):
            continue
        scen = build(40, 320)
        _, computed, true = void_series(scen)
        agree = sum(int(c == t) for c, t in zip(computed, true))
        print(f"[{name}] void count (computed) == ground truth on "
              f"{agree}/{len(true)} frames")
        path = render_void_movie(scen, out / f"{name}.gif")
        print(f"        wrote {path}")
    if which in ("all", "genus"):
        print("[genus] wrote", render_genus_figure(out / "genus.png"))


if __name__ == "__main__":
    main()
