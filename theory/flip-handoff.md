# The chromatic flip-handoff theorem

*Paper theory — application/paper-specific, lives in `chromatic-cells`.*
*Status: **proved** for the generic (single transversal bistellar) flip — the
case the theorem covers. The proof rests only on results of the chromatic-alpha
paper (Cultrera di Montesano–Draganov–Edelsbrunner–Saghafian, *Chromatic Alpha
Complexes*, 2024; lemma/theorem numbers below refer to it) plus one genericity
hypothesis (a single, transversal cosphericity), whose failure is the
degenerate case requiring SoS and is explicitly out of scope here.*

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

> Let `t★` be a **generic** chromatic Delaunay bistellar flip: at `t★` exactly
> `d+s+2` lifted points `Q̂` lie on a common sphere `Ŝ ⊂ ℝ^{d+s}` (centre
> `ĉ=(z★,w★)`, radius `ρ̂`), empty of all other lifted points, this is the only
> degeneracy at `t★`, and the cosphericity is transversal (so `Ŝ` — hence `ĉ` —
> is the well-defined common limit of the participating circumspheres). Let
> `D` be the simplices in `Del(χ)(t★⁻)∖Del(χ)(t★⁺)` (dying) and `A` those in
> `Del(χ)(t★⁺)∖Del(χ)(t★⁻)` (arriving). Then every simplex in `D ∪ A` has the
> same squared empty-stack radius at `t★`:
>
> `ϱ_{t★}(σ) = R²★` for all `σ ∈ D ∪ A`,    where  `R²★ = ρ̂² − min_i ‖w★ − h·v_i‖²`
>
> is the value of the empty stack obtained by projecting `Ŝ` to `ℝ^d`
> (a stack with **unequal** per-colour radii, centred at `z★`).

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

- *`F_σ(t★)` is the single point `{ĉ}`.* The flip is *interior*: `Ŝ` is a
  genuine finite sphere (`ρ̂ < ∞`), so every Voronoi face `F_σ(t)` for
  `σ ∈ D ∪ A` is **bounded** (a hull flip, where `Ŝ` degenerates to a
  hyperplane and `F_σ` is unbounded, is reduced to the interior case by the
  standard vertex-at-infinity compactification — the ghost-simplex device the
  kinetic-Delaunay engine already uses — under which `Ŝ` passes through the ghost
  vertex). The vertices of `F_σ(t)` are the circumcentres of the
  top-dimensional simplices of the flip that contain `σ`: every
  `(d+s)`-dimensional coface of a dying `σ` is itself a dying flip simplex (a
  surviving coface would keep `σ` alive), so each is a `(d+s)`-simplex spanned by
  `d+s+1` points of `Q̂`. As `t → t★` such a simplex's `d+s+1` vertices lie on
  the limiting sphere `Ŝ`, and the unique `(d+s−1)`-sphere through `d+s+1`
  affinely independent points of `Ŝ` is `Ŝ` itself, so its circumcentre converges
  to `Ŝ`'s centre `ĉ` (transversality keeps the simplex non-degenerate up to
  `t★`, making `ĉ` the well-defined common limit). Hence **all** vertices of
  `F_σ(t)` converge to the single point `ĉ`, and the bounded face `F_σ` — their
  convex hull — collapses to `{ĉ}` at `t★`. With `ĉ ∈ F_σ(t★)` this gives
  `F_σ(t★) = {ĉ}`, i.e. the admissible centre is forced to `z★ = proj(ĉ)`.

  (Consistency check with "`σ` dies/arrives": for a dying `σ`, `F_σ` is a
  positive-dimensional polytope for `t<t★` and empty for `t>t★`; the generic
  transition is precisely this collapse through the single point `z★`. Arriving
  simplices are symmetric, with `t<t★` and `t>t★` exchanged.)

Therefore `ϱ_{t★}(σ) = R²★` for every `σ ∈ D ∪ A`. Because each `F_σ(t★) = {ĉ}`
shares the same point `ĉ`, all flip simplices share the same centre `z★` and the
same value `R²★`. ∎

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
— together with one genericity hypothesis: that `t★` is a single transversal
cosphericity, so `Ŝ` and its centre `ĉ` are well defined as the common limit of
the participating circumspheres. No appeal to an unproved relation between the
lifted-Delaunay radius and the empty-stack radius is required.

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
