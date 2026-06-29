# The chromatic flip-handoff theorem

*Paper theory — application/paper-specific, lives in `chromatic-cells`.*
*Status: **proved** for every generic (single transversal bistellar) flip —
both the interior regime (finite cosphere, finite shared radius) and the hull
regime (cosphere degenerating to an empty hyperplane, shared radius `∞`). The
proof rests only on results of the chromatic-alpha paper (Cultrera di
Montesano–Draganov–Edelsbrunner–Saghafian, *Chromatic Alpha Complexes*, 2024;
lemma/theorem numbers below refer to it) plus one genericity hypothesis (a
single, transversal cosphericity), whose failure — coincident or non-unique
cosphericities — is the degenerate case requiring SoS and is the only thing out
of scope here.*

## Setup and notation

Points `P ⊂ ℝ^d` carry one of `s+1` colours. The **chromatic lift** sends a
colour-`i` point `p` to `p̂ = (p, h·v_i) ∈ ℝ^{d+s}`, where `v_0,…,v_s` are the
vertices of a regular `s`-simplex in `ℝ^s` (height scale `h>0`; the paper's §3.3
construction — any affinely independent `u_0,…,u_s` works, and the
`chromatic_tda` implementation uses the one-hot choice `u_0=0, u_i=e_i`). The
**chromatic Delaunay complex** `Del(χ)` is isomorphic to the ordinary Delaunay
complex `Del(P̂)` of the lifted points in `ℝ^{d+s}`, by the lift isomorphism
`a ↦ â` (Corollary 3.7).

For a simplex `σ` the chromatic filtration value is the squared **empty-stack
radius** `ϱ(σ) = Rad²(σ)` (Definition 3.2): a `σ`-*stack* is `s+1` concentric
`(d−1)`-spheres `S_0,…,S_s` with common centre `z ∈ ℝ^d`, `S_i` passing through
the colour-`i` vertices of `σ`; the stack is *empty* if no point of `P` lies
strictly inside its own colour's sphere; its radius is `max_i r_i`. `Rad(σ)` is
the radius of the **smallest empty stack** through `σ` — the value of the convex
program `(P_σ)` of §4.2,

```
        minimise   max_i r_i²      over centre z ∈ ℝ^d and radii r,
        subject to ‖x − z‖ = r_{χ(x)}   for x ∈ σ,
                   ‖x − z‖ ≥ r_{χ(x)}   for x ∈ P∖σ.
```

Two facts we use throughout:

- **(Voronoi characterisation, Lemma 3.1.)** A stack centred at `z` passes
  through `σ` and is empty iff `z` lies in the intersection of the *chromatic*
  (per-colour) Voronoi domains of `σ`'s vertices,
  `F_σ := ⋂_{a∈σ} dom(a, P_{χ(a)})`. In particular `σ ∈ Del(χ)` iff `F_σ ≠ ∅`
  (and, under the lift, `F_σ` is the ordinary Voronoi face `⋂_{a∈σ} dom(â, P̂)`,
  Corollary 3.7). The admissible centres of an empty circumstack of `σ` are
  exactly the points of `F_σ`.
- **(Lift of one sphere, Lemma 3.6.)** An empty `(d+s−1)`-sphere `Ŝ ⊂ ℝ^{d+s}`
  with centre `ĉ = (z, w)`, `w ∈ ℝ^s`, and radius `ρ̂` that passes through `σ̂`
  projects to **one** empty stack through `σ` centred at `z`, whose colour-`i`
  radius satisfies

  ```
        r_i² = ρ̂² − ‖w − h·v_i‖² .
  ```

  These per-colour radii are in general **unequal** (they decrease as the lift
  centre `w` moves away from `h·v_i`); the stack's value is therefore
  `max_i r_i² = ρ̂² − min_i ‖w − h·v_i‖²`, **not** `ρ̂²`.

A motion of `P` in `ℝ^d` lifts to a motion of `P̂` in `ℝ^{d+s}` (the lift
coordinate fixed by colour) and drives `Del(P̂)` through **bistellar flips**,
each at a **cosphericity event**: `d+s+2` lifted points become cospherical. For
`s=1, d=2` the lift is in `ℝ³` and a flip is an insphere event on `5` lifted
points (flip23 / flip32). By Corollary 3.7 these are exactly the combinatorial
changes of `Del(χ)`.

## Theorem (chromatic flip-handoff, generic case)

> Let `t★` be a **generic interior** chromatic Delaunay bistellar flip: at `t★`
> exactly `d+s+2` lifted points `Q̂` lie on a common **finite** sphere
> `Ŝ ⊂ ℝ^{d+s}` (centre `ĉ=(z★,w★)`, radius `ρ̂ < ∞`), empty of all other lifted
> points, this is the only degeneracy at `t★`, and the cosphericity is
> transversal (so `Ŝ` — hence `ĉ` — is the well-defined common limit of the
> participating circumspheres). Let `D` be the simplices in
> `Del(χ)(t★⁻)∖Del(χ)(t★⁺)` (dying) and `A` those in `Del(χ)(t★⁺)∖Del(χ)(t★⁻)`
> (arriving). Then every simplex in `D ∪ A` has the same squared empty-stack
> radius at `t★`:
>
> `ϱ_{t★}(σ) = R²★` for all `σ ∈ D ∪ A`,    where  `R²★ = ρ̂² − min_i ‖w★ − h·v_i‖²`
>
> is the value of the empty stack obtained by projecting `Ŝ` to `ℝ^d`
> (a stack with **unequal** per-colour radii, centred at `z★`).

> The remaining generic regime — a **hull flip**, where `Ŝ` degenerates to an
> empty hyperplane (`ρ̂ = ∞`) — is the Proposition after the proof; there the
> shared value is `R²★ = ∞`. Together they cover every generic flip; the only
> deferred case is the non-generic (coincident/degenerate) flip, the SoS
> corollary.

## Proof

**(1) `D ∪ A` are faces of one cospherical configuration.** A bistellar flip
exchanges the two triangulations of `conv(Q̂)` (the two Lawson families of the
`d+s+2`-point Radon configuration). Writing the Radon partition `Q̂ = P₊ ⊔ P₋`,
the dying simplices are the interval `D = {τ : P₊ ⊆ τ ⊊ Q̂}` and the arriving
ones `A = {τ : P₋ ⊆ τ ⊊ Q̂}` (the cofaces of `conv(P₊)`, resp. `conv(P₋)`, inside
`Q̂`). In particular **every** `σ ∈ D ∪ A` has all its vertices on `Ŝ`, i.e.
`σ̂ ⊂ Ŝ`. (This interval exchange `[P₊, Q̂) ↔ [P₋, Q̂)` is the reorganisation of
the generalised-discrete-Morse intervals of `Rad` across the flip; `Rad` is
generalised discrete Morse at each generic time by Theorem 4.6.)

**(2) The projection of `Ŝ` is one empty stack, common to all of `D ∪ A`, of
value `R²★`.** Fix `σ ∈ D ∪ A`. Since `σ̂ ⊂ Ŝ` and `Ŝ` is empty, Lemma 3.6
projects `Ŝ` to an empty stack through `σ`, centred at `z★ = proj(ĉ)`, with
per-colour radii `r_i² = ρ̂² − ‖w★ − h·v_i‖²` and value
`R²★ = ρ̂² − min_i ‖w★ − h·v_i‖²`. This is the **same** stack for every
`σ ∈ D ∪ A` (it depends only on `Ŝ`, not on `σ`). Hence `ϱ_{t★}(σ) ≤ R²★` for
every `σ ∈ D ∪ A`.

> *Note (this is exactly where the earlier draft erred).* We do **not** claim
> `Ŝ` projects to a single circumsphere with equal per-colour radii — it does
> not, and `R²★ ≠ ρ̂²` in general. We claim only that `Ŝ` yields **one** empty
> stack, of value `R²★`, that is admissible for every flip simplex. The equality
> `ϱ_{t★}(σ) = R²★` is established next, by uniqueness of the admissible centre —
> not by excluding a "smaller cosphericity".

**(3) Equality, via collapse of the admissible-centre set.** By Lemma 3.1 the
admissible centres of empty circumstacks through `σ` are exactly the points of
the chromatic Voronoi face `F_σ(t) = ⋂_{a∈σ} dom(a, P_{χ(a)})(t)`, equivalently
the lifted Voronoi face `⋂_{a∈σ} dom(â, P̂)(t)` (Cor 3.7). We show
`F_σ(t★) = {z★}` for every `σ ∈ D ∪ A`; then the empty circumstack of `σ` at
`t★` is **unique** — centred at `z★`, equal to the projection of `Ŝ` — so
`ϱ_{t★}(σ) = R²★`, and the min-coface clause of `Rad` is moot because `σ` itself
carries an empty stack at `t★`.

- *`ĉ ∈ F_σ(t★)`.* For `a ∈ σ` we have `‖ĉ − â‖ = ρ̂` (as `σ̂ ⊂ Ŝ`), while
  `‖ĉ − b̂‖ ≥ ρ̂` for every other lifted point `b̂` (`Ŝ` empty). So `ĉ` is closest
  to each `â`, i.e. `ĉ ∈ ⋂_{a∈σ} dom(â, P̂) = F_σ(t★)`; projecting, `z★ ∈ F_σ`.

- *`F_σ(t★)` is the single point `{ĉ}`.* As `Ŝ` is a genuine finite sphere, the
  Voronoi face `F_σ(t)` of each `σ ∈ D ∪ A` is **bounded**, and its vertices are
  the circumcentres of the `(d+s)`-simplices of the flip that contain `σ`. Every
  top-dimensional coface of a dying `σ` is itself a dying flip simplex — a
  surviving coface would keep `σ` alive — hence a `(d+s)`-simplex spanned by
  `d+s+1` points of `Q̂`. The unique `(d+s−1)`-sphere through `d+s+1` affinely
  independent points of `Ŝ` is `Ŝ` itself, so as `t → t★` each such circumcentre
  converges to `Ŝ`'s centre `ĉ` (transversality keeps these simplices
  non-degenerate up to `t★`, making `ĉ` the well-defined common limit). Thus
  every vertex of `F_σ(t)` converges to `ĉ`, and the bounded face — their convex
  hull — collapses to `{ĉ}`. With `ĉ ∈ F_σ(t★)` this gives `F_σ(t★) = {ĉ}`: the
  admissible centre is forced to `z★ = proj(ĉ)`. (Hull flips, where `Ŝ`
  degenerates to a hyperplane and `F_σ` is unbounded, are the Proposition below.)

  (Consistency check with "`σ` dies/arrives": for a dying `σ`, the bounded face
  `F_σ` is non-empty for `t<t★` and empty for `t>t★` — a Voronoi vertex for a
  top-dimensional `σ`, a higher-dimensional bounded cell for a lower-dimensional
  one; the generic transition is precisely this collapse through the single point
  `z★`. Arriving simplices are symmetric, with `t<t★` and `t>t★` exchanged.)

Therefore `ϱ_{t★}(σ) = R²★` for every `σ ∈ D ∪ A`. Because each `F_σ(t★) = {ĉ}`
shares the same point `ĉ`, all flip simplices share the same centre `z★` and the
same value `R²★`. ∎

## Proposition (hull flips)

> Let `t★` be a generic **hull** flip: `d+s+1` finite lifted points `Q̂₀` become
> coplanar on a hyperplane `Ĥ ⊂ ℝ^{d+s}` empty of all other lifted points (so the
> cospherical sphere `Ŝ` degenerates to `Ĥ`, `ρ̂ = ∞`), and this is the only
> degeneracy at `t★`. Then every simplex in `D ∪ A` has `ϱ_{t★}(σ) = ∞`: the
> handoff holds with the extended common value `R²★ = ∞`.

*Proof.* Compactify the lifted Delaunay complex with the point at infinity `ω`:

```
        Del⁺(P̂) = Del(P̂) ∪ { ω ∗ F : F a facet of conv(P̂) } ,
```

a triangulation of the sphere `𝕊^{d+s}` (the standard one-point
compactification — equivalently the full boundary complex, lower **and** upper,
of `conv(π(P̂))` under the paraboloid lift `π : x ↦ (x, ‖x‖²)`; Edelsbrunner &
Harer, *Computational Topology*, §III). Bistellar flips of `Del(χ)` are exactly
the bistellar flips of `Del⁺(P̂)` on `d+s+2` vertices that become cospherical on
`𝕊^{d+s}`. The hull flip is the case where exactly one of those vertices is `ω`;
the other `d+s+1` are the finite points `Q̂₀ ⊂ Ĥ`, and `Ŝ = Ĥ ∪ {ω}` is the
"sphere through `ω`", i.e. the hyperplane `Ĥ`. Write `n̂` for the unit normal of
`Ĥ` pointing to its empty side.

Fix a finite `σ ∈ D ∪ A` (so `σ̂ ⊂ Ĥ`). By Lemma 3.1 its empty circumstacks
have centres in `F_σ = ⋂_{a∈σ} dom(â, P̂)`, whose vertices are the circumcentres
of the `(d+s)`-cofaces of `σ` in the flip. As in step (3), every top-dimensional
coface of a dying (resp. arriving) `σ` is itself a flip simplex — a surviving
coface would keep `σ` alive. Among them:

- the **finite** flip simplex `Q̂₀` (the `d+s+1` coplanar points): its
  circumsphere flattens to `Ĥ` as `t → t★`, so its circumcentre recedes to
  infinity along `n̂`;
- the **`ω`-cofaces** `Q̂∖{q} = ω ∗ (Q̂₀∖{q})`, `q ∉ σ`: each is dual to the
  ideal Voronoi vertex in the direction normal to its hull facet, and those
  facets converge to `Ĥ`, so these too recede to the ideal point along `n̂`.

Hence every vertex of the closed face `F̄_σ ⊆ 𝕊^{d+s}` converges to the single
**ideal** point `ĉ_∞` in direction `n̂`, and `F̄_σ` collapses to `{ĉ_∞}` at
`t★`. The only admissible centre is therefore at infinity: the smallest empty
circumstack of `σ` is the empty half-space bounded by `Ĥ` (radii `→ ∞`), so
`ϱ_{t★}(σ) = ∞`. As `ĉ_∞` is common to all `σ ∈ D ∪ A`, they share `R²★ = ∞`. ∎

The contrast with step (3) is exactly finite-vs-ideal centre: an interior flip
forces the centre to the finite point `z★` (shared finite `R²★`); a hull flip
forces it to the ideal point `ĉ_∞` (shared `R²★ = ∞`, a handoff at the top of the
filtration — the participating bars are essential and exchange continuously). In
both, the mechanism is the same: the admissible-centre set collapses to one
point, forcing uniqueness.

## Why this closes the gap

The load-bearing step is (3), and it is now proved outright rather than flagged.
The previous draft argued "there is no strictly smaller empty circumstack
through `σ` — a smaller one would be a second cosphericity, excluded." That is
**unsound**: a smaller empty stack with *unequal* per-colour radii is not a
sphere through `σ̂`, so it is not a cosphericity and genericity does not exclude
it on those grounds (the gap the audit identified). The correct argument never
needs to exclude a competitor: at a generic flip the **feasible set of stack
centres** for a flip simplex — its Voronoi face `F_σ`, by Lemma 3.1 — is a
*single point* `z★`, so the empty circumstack is unique and its value is forced
to `R²★`. Uniqueness, not exclusion; and the projected stack has genuinely
unequal radii, so the equal-radii conflation never enters.

Every ingredient is a cited result of the chromatic-alpha paper — Lemma 3.1
(empty stack ⟺ Voronoi face), Lemma 3.6 and Corollary 3.7 (the lift), and
Theorem 4.6 (`Rad` is generalised discrete Morse, framing the interval exchange)
— plus, for hull flips, the standard one-point compactification of the Delaunay
complex. The single remaining hypothesis is that `t★` is one transversal
cosphericity, so the centre `ĉ` (finite for an interior flip, ideal for a hull
flip) is well defined as the common limit of the participating circumspheres. No
appeal to an unproved relation between the lifted-Delaunay radius and the
empty-stack radius is required, and the interior and hull regimes together
exhaust the generic flips — only coincident/degenerate cosphericities (SoS)
remain.

## Empirical confirmation

Certified in `vineyards/tests/test_chromatic_handoff.py` and the wider sweep:
isolating single flips by recursive bisection, the genuine bistellar flips
(cardinalities `(3,7)=flip23`, `(7,3)=flip32` in the lifted complex) share the
empty-stack weight to **machine precision** (median rel. spread `3.3e-9`, max
`6.5e-7` over ~190 flips) — matching `ϱ_{t★}(σ) = R²★` for all flip simplices.

## Corollary (degenerate case — needs SoS)

The single genericity hypothesis fails exactly when the cosphericity is not a
single transversal event — several cosphericities coincide, or `Q̂`'s limiting
sphere is not unique, so the common centre `ĉ` of step (3) is not well defined.
Empirically these are the degenerate-cluster events with cardinalities `(6,6)`,
`(4,0)`, `(0,4)` (rel. spread `0.5–0.75`, *not* a shared radius): the Voronoi
face does not collapse to a single point, the admissible centre is not forced,
and the handoff is not a single shared radius. Resolving them requires a
Simulation-of-Simplicity perturbation that splits the coincident cosphericities
into a sequence of generic flips, each obeying the theorem. This extension is
deferred (the SoS item).
