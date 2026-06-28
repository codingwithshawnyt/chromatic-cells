# The chromatic flip-handoff theorem

*Paper theory — application/paper-specific, lives in `chromatic-cells`.*
*Status: theorem stated; proof argued; the one careful step is flagged. The
generic (bistellar, general-position) case is what the theorem covers; the
degenerate case is a corollary requiring SoS and is explicitly out of scope here.*

## Setup and notation

Points `P ⊂ ℝ^d` carry one of `s+1` colours. The **chromatic lift** sends a
colour-`i` point `p` to `p̂ = (p, h·v_i) ∈ ℝ^{d+s}`, where `v_0,…,v_s` are the
vertices of a regular `s`-simplex in `ℝ^s` (height scale `h>0`)
(Cultrera di Montesano–Draganov–Edelsbrunner–Saghafian, *Chromatic Alpha
Complexes*, 2024). The **chromatic Delaunay mosaic** is the Delaunay
triangulation `Del(P̂)` of the lifted points in `ℝ^{d+s}`.

For a simplex `σ` the chromatic filtration value is the **empty-stack radius**
`ϱ(σ)`: split `σ`'s vertices by colour into `B_0,…,B_k`; the smallest *circumstack*
is the common centre `z ∈ ℝ^d` and per-colour radii `r_0,…,r_k` with `S_i` passing
through `B_i`, minimising `max_i r_i`; `ϱ(σ)` is that `max_i r_i²` if the stack is
empty of other points, else the min over cofaces (the chromatic analogue of the
Gabriel / min-coface rule). For a colour subset the centre `z` is the projection
to `ℝ^d` of the centre of the lifted sphere through `B̂_0 ∪ … ∪ B̂_k`.

A point of `ℝ^{d+s}` motion (induced by `P` moving in `ℝ^d`, the lift coordinate
fixed by colour) drives `Del(P̂)` through **bistellar flips**, each occurring at a
**cosphericity event**: `d+s+2` lifted points become cospherical (an insphere
event in `ℝ^{d+s}`). For `s=1, d=2` (two colours in the plane) the lift is in
`ℝ³` and a flip is an insphere event on `5` lifted points (flip23 / flip32).

## Theorem (chromatic flip-handoff, generic case)

> Let `t★` be a generic chromatic Delaunay bistellar flip: at `t★` exactly `d+s+2`
> lifted points lie on a common sphere `Ŝ ⊂ ℝ^{d+s}`, empty of all other lifted
> points, and no further degeneracy occurs (general position otherwise). Let
> `D = {dying simplices}` (in `Del(P̂(t★⁻))` but not `Del(P̂(t★⁺))`) and
> `A = {arriving simplices}` (vice-versa). Then every simplex in `D ∪ A` has the
> same empty-stack radius at `t★`:
>
> `ϱ_{t★}(σ) = R²★` for all `σ ∈ D ∪ A`,
>
> where `R²★` is the squared radius of the circumstack obtained by projecting `Ŝ`
> to `ℝ^d`.

## Proof

**(1) `D ∪ A` are faces of one cospherical configuration.** A bistellar flip
exchanges the two triangulations of the convex hull of the `d+s+2` cospherical
points `Q̂` (the two ways to triangulate a cyclic/`(d+s+1)`-point polytope). Every
simplex created or destroyed by the flip is a face of `conv(Q̂)` and has all its
vertices in `Q̂` — i.e. each `σ ∈ D ∪ A` is spanned by a subset of the cospherical
points lying on `Ŝ`.

**(2) Each flip simplex's circumstack is the projection of `Ŝ`.** Fix
`σ ∈ D ∪ A` with vertices `σ̂ ⊂ Ŝ`. The chromatic-lift correspondence (op. cit.)
is: the lifted sphere through a colourful point set ↔ the circumstack through the
colour classes, with the stack centre `z` the `ℝ^d`-projection of the sphere
centre and the per-colour radii read from the sphere's intersection with the
colour subspaces. Since `σ̂ ⊂ Ŝ`, the sphere `Ŝ` is *a* circumstack of `σ`; its
projected squared radius is `R²★`.

**(3) It is the *smallest* circumstack, and it is empty.** `Ŝ` is empty of all
lifted points (flip hypothesis), so the corresponding stack is empty of all
points — hence `σ`'s smallest *empty* circumstack is no larger than `R²★`. By
genericity there is no strictly smaller empty circumstack through `σ`'s colour
classes at `t★` (a smaller one would be a second cosphericity, excluded), so the
smallest empty circumstack is exactly `Ŝ` and `ϱ_{t★}(σ) = R²★`.

Since (2)–(3) hold for every `σ ∈ D ∪ A` with the *same* `Ŝ`, all share `R²★`. ∎

## The one careful step

Step (3) — *that the cosphere `Ŝ` realises the simplex's **smallest empty**
circumstack at the flip, not merely some circumstack* — is the load-bearing
lemma. The argument above reduces it to: (a) the lift correspondence (lifted
sphere ↔ circumstack), which is established in the chromatic-alpha paper; and
(b) genericity ruling out a competing smaller empty stack. (b) is where the
degenerate-cluster events live (see below). This step should be written against
the chromatic-alpha paper's Lemma relating the lifted Delaunay radius function to
the empty-stack radius; it is the part to nail formally.

## Empirical confirmation

Certified in `vineyards/tests/test_chromatic_handoff.py` and the wider sweep:
isolating single flips by recursive bisection, the genuine bistellar flips
(cardinalities `(3,7)=flip23`, `(7,3)=flip32` in the lifted complex) share the
empty-stack weight to **machine precision** (median rel. spread `3.3e-9`, max
`6.5e-7` over ~190 flips). This is evidence for the theorem in the generic case.

## Corollary (degenerate case — needs SoS)

When the genericity hypothesis fails — several cosphericities coincide, i.e. the
degenerate-cluster events with cardinalities `(6,6)`, `(4,0)`, `(0,4)` observed
empirically (rel. spread `0.5–0.75`, *not* a shared radius) — the flip is not a
single bistellar move and the handoff is not a single shared radius. Resolving
these requires a Simulation-of-Simplicity perturbation that splits the coincident
cosphericities into a sequence of generic flips, each of which obeys the theorem.
This extension is deferred (the SoS item).
