"""How do you get from a point cloud to a topology reading?  Grow balls.

The step everyone misses -- and the reason the output looks abstract until you've
seen it.  Two short animations, each ONE frozen frame of data, sweeping a single
knob (the ball radius r):

  * explain_2d_loop   -- a loop (H1) in 2D, where the hole is unambiguous.
  * explain_3d_sphere -- the same idea on a hollow sphere: a void (H2), a lumen.

Grow a ball on every point.  Nothing, then the balls meet and a hole / void
appears (it is BORN), then the balls fill it in and it closes (it DIES).  That
birth-to-death lifespan is one bar -- one topological feature.  Counting the long
bars is the void count that moves over time in the scene demos (void_demo.py).

Left: the balls growing on the points.  Right: the bar filling in, from birth to
death.  No persistence-diagram scatter, no vines -- just the one idea.

Run:  python examples/pipeline_explainer.py            # writes examples/media/*.gif
"""
import os
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as _fm
from matplotlib.patches import Circle
from PIL import Image
import gudhi

# Arial (Nature/Cell figure font), matching the void demo.
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
BLUE, ORANGE, MUTE = "#4477AA", "#EE6677", "#b8c0cc"


def _fig_to_rgb(fig):
    fig.canvas.draw()
    return np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()


# ---------------------------------------------------------------------------
# The right panel, shared: a single barcode bar (birth -> death of one feature)
# ---------------------------------------------------------------------------

def _barcode(ax, r, b0, d0, kind, xmax):
    """One bar on a 'ball radius' axis: it starts at the birth radius and grows to
    the current radius while the feature is open, then locks at the death radius."""
    ax.clear()
    end = min(r, d0)
    born, filled = r >= b0, r >= d0
    # the bar (the feature's lifespan so far)
    if born:
        ax.plot([b0, end], [0.5, 0.5], lw=16, solid_capstyle="butt",
                color=(ORANGE if not filled else "#8fa1b3"),
                alpha=(1.0 if not filled else 0.9))
    # birth / death guide ticks
    ax.plot([b0, b0], [0.30, 0.70], color=MUTE, lw=1.5)
    ax.text(b0, 0.20, f"birth\nr={b0:.2f}", ha="center", va="top", fontsize=9, color=FG)
    if filled:
        ax.plot([d0, d0], [0.30, 0.70], color=MUTE, lw=1.5)
        ax.text(d0, 0.20, f"death\nr={d0:.2f}", ha="center", va="top", fontsize=9, color=FG)
    # the sweep marker (current radius)
    ax.axvline(r, color=FG, lw=1.4, alpha=0.55)
    ax.text(r, 0.80, f"r = {r:.2f}", ha="center", va="bottom", fontsize=10, color=FG)
    # state label
    if not born:
        msg, col = f"no {kind.split()[0]} yet  —  keep growing", MUTE
    elif not filled:
        msg, col = f"the {kind} is OPEN  (born)", ORANGE
    else:
        msg, col = f"the {kind} has FILLED  (died)", "#5b6b7b"
    ax.text(xmax / 2, 0.96, msg, ha="center", va="top", fontsize=13,
            fontweight="bold", color=col)
    ax.set_xlim(0, xmax); ax.set_ylim(0, 1.05)
    ax.set_yticks([])
    ax.set_xlabel("ball radius  r")
    ax.set_title(f"the barcode:  one bar = one {kind}", fontsize=12, color=FG)
    for sp in ax.spines.values():
        sp.set_color(SPINE)
    ax.spines["left"].set_visible(False)
    ax.tick_params(colors=FG)
    ax.grid(True, axis="x", color=GRID, lw=0.7)


# ---------------------------------------------------------------------------
# 2D: a loop (H1)
# ---------------------------------------------------------------------------

def _ring(n, R=1.0, center=(0, 0)):
    a = np.linspace(0, 2 * np.pi, n, endpoint=False)
    return np.column_stack([R * np.cos(a), R * np.sin(a)]) + np.array(center)


def _h1(points):
    st = gudhi.AlphaComplex(points=points).create_simplex_tree()
    st.compute_persistence()
    ivs = [(np.sqrt(b), np.sqrt(d))
           for b, d in st.persistence_intervals_in_dimension(1) if d < np.inf]
    return max(ivs, key=lambda bd: bd[1] - bd[0])


def explain_2d_loop(filename, *, n=22, side=520, fps=11):
    pts = _ring(n)
    b0, d0 = _h1(pts)
    xmax = d0 * 1.35
    grow = np.concatenate([np.linspace(0.0, d0 * 1.18, 30), [d0 * 1.18] * 7])
    frames = []
    for r in grow:
        fig, (axL, axR) = plt.subplots(1, 2, figsize=(2 * side / 100, side / 100), dpi=100)
        fig.patch.set_facecolor("white")
        born, filled = r >= b0, r >= d0
        for p in pts:                                          # growing balls
            axL.add_patch(Circle(p, r, color=BLUE,
                                 alpha=0.16 if not filled else 0.22, lw=0))
        axL.scatter(pts[:, 0], pts[:, 1], s=24, color=BLUE, zorder=5)
        if born and not filled:                                # name the hole
            axL.text(0, 0, "hole", ha="center", va="center", color=ORANGE,
                     fontsize=13, fontweight="bold")
        axL.set_xlim(-2.1, 2.1); axL.set_ylim(-2.1, 2.1); axL.set_aspect("equal")
        axL.set_xticks([]); axL.set_yticks([])
        axL.set_title("grow a ball on each point", fontsize=12, color=FG)
        for sp in axL.spines.values():
            sp.set_color(SPINE)
        _barcode(axR, float(r), b0, d0, "loop (H1)", xmax)
        fig.tight_layout()
        frames.append(_fig_to_rgb(fig)); plt.close(fig)
    _save(frames, filename, fps)
    return filename


# ---------------------------------------------------------------------------
# 3D: a void (H2) -- a hollow sphere, the lumen setting
# ---------------------------------------------------------------------------

def _fib_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)])


def _h2(points):
    st = gudhi.AlphaComplex(points=points).create_simplex_tree()
    st.compute_persistence(persistence_dim_max=True)
    ivs = [(np.sqrt(b), np.sqrt(d))
           for b, d in st.persistence_intervals_in_dimension(2) if d < np.inf]
    return max(ivs, key=lambda bd: bd[1] - bd[0])


def explain_3d_sphere(filename, *, n=120, side=520, fps=11):
    """Growing balls on a hollow sphere close a VOID (H2) then fill it.  The top
    half is cut away and we look straight down into the bowl, so the cavity in the
    middle (and its filling in) is unmistakable."""
    import pyvista as pv
    shell = _fib_sphere(n)
    b0, d0 = _h2(shell)
    xmax = d0 * 1.35
    grow = np.concatenate([np.linspace(0.06, d0 * 1.18, 26), [d0 * 1.18] * 7])

    def panel3d(r):
        pl = pv.Plotter(off_screen=True, window_size=(side, side))
        pl.set_background("white")
        balls = pv.PolyData(shell).glyph(
            scale=False, orient=False,
            geom=pv.Sphere(radius=r, theta_resolution=18, phi_resolution=18))
        # cut the top half off and look down into the bowl: opaque so the hollow
        # (and its filling in) reads clearly, not a translucent mush.
        pl.add_mesh(balls.clip(normal="z", origin=(0, 0, 0.05)), color=BLUE,
                    opacity=1.0, smooth_shading=True, specular=0.2)
        pl.camera_position = [(1.2, -1.1, 5.6), (0, 0, -0.15), (0, 1, 0)]
        return np.asarray(pl.screenshot(return_img=True))[..., :3]

    frames = []
    for r in grow:                                     # left panel via matplotlib
        left = panel3d(float(r))                       # so the title never sits on the bowl
        fig, (axL, axR) = plt.subplots(1, 2, figsize=(2 * side / 100, side / 100), dpi=100)
        fig.patch.set_facecolor("white")
        axL.imshow(left); axL.set_xticks([]); axL.set_yticks([])
        axL.set_title("grow a ball on each point  (top cut away)", fontsize=12, color=FG)
        for sp in axL.spines.values():
            sp.set_color(SPINE)
        _barcode(axR, float(r), b0, d0, "void (H2)", xmax)
        fig.tight_layout()
        frames.append(_fig_to_rgb(fig)); plt.close(fig)
    _save(frames, filename, fps)
    return filename


def _save(frames, filename, fps):
    Path(filename).parent.mkdir(parents=True, exist_ok=True)
    imgs = [Image.fromarray(f).convert("P", palette=Image.ADAPTIVE, colors=128) for f in frames]
    imgs[0].save(filename, save_all=True, append_images=imgs[1:],
                 duration=int(1000 / fps), loop=0, optimize=True, disposal=2)


if __name__ == "__main__":
    out = Path(__file__).parent / "media"
    print("2D loop ->", explain_2d_loop(out / "explain_loop.gif"))
    print("3D void ->", explain_3d_sphere(out / "explain_void.gif"))
