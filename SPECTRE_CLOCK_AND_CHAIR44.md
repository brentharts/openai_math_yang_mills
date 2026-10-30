# Symmetry, Doubling and "Clocks": Spectre, Chair44 and the Preprints

A deeper look at the connections between the Spectre and Chair44 monotiles and the repository's preprints that [HIERARCHIES.md](HIERARCHIES.md) listed, with exact computations and figures. It is a fork-level note, not an official OpenAI document. The word "clock" comes from [HR] §10, which tests whether the Spectre's substitution hierarchy "can serve as a clock" (abstract).

**How to read the evidence.** Each claim below is tagged by where it comes from:
- **[source]**: stated in a paper I read. For [T] and [HR] I read fetch-tool summaries of the PDFs, so section numbers and short quotes are as the summaries gave them. For the repository's preprints (O3, E229, E304, E155) I read the extracted text directly.
- **[computed]**: reproduced by [tools/hierarchy_figures.py](tools/hierarchy_figures.py) (`make figures`), which builds on the transcribed Chair44 atlas of [CHAIR44_MODEL.md](CHAIR44_MODEL.md). These are exact checks on a transcription of [T]'s Figure 7 and Table 1, not a proof about [T]'s tile.
- **[mine]**: an analogy or construction of this fork.
- *(background)*: general knowledge, not from a source here.

Citation keys follow [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md#citation-key): [T] Tsiokos, [HR] Hartshorn–Romanov (Research Square, unrefereed), [F] Flicker, [E155] entry 155. New keys: **[O3]** the [two-dimensional O(3) paper](preprints/The-canonical-massive-continuum-limit-of-the-two-dimensional-O3-model-October-4-2026/) (entry 215), **[E229]** the [three-state reconstruction paper](preprints/The-exact-reconstruction-threshold-for-the-three-state-symmetric-channel-September-25-2026/), **[E304]** the [Hilbert–Smith paper](preprints/The-Hilbert-Smith-conjecture-in-every-finite-dimension-September-23-2026/).

## Summary

| Connection | What the closer look found |
|---|---|
| Symmetry: "only mildly aperiodic" | The term is defined in [T]; every tiling has a finite symmetry group, and some have order 8 (D4). I reproduce the order-8 case computationally on substitution fixed points. |
| Block renormalization | Both sides have an exact "doubling defect". In [O3] it obeys a **quadratic** recursion, so it decays doubly exponentially. The Chair44 hierarchy has an exact contact count whose defect decays only **geometrically** (2⁻ⁿ). The analogy is in form, not in rate. |
| Hilbert–Smith (304) | Confirmed to be vocabulary only; the theorem is about manifolds, a tiling hull is not one. |
| Reconstruction on trees (229) and the "clock" tests | Different questions. A Chair44 version of the "records" test gives a radius growing about 2.3 to 2.8 times per level; a noisy-patch version shows smooth degradation, with no threshold. |
| Anderson (261), Parisi (221), triangular lattice (090, 232, 237) | Still only thematic or background. |

## 1. How aperiodic is Chair44?

### What the sources say

- **Definitions [source].** [T] §1.1 (citing Coulbois, Gajardo, Guillon and Lutfalla): "a tile all of whose tilings have finite symmetry group is mildly aperiodic". Mozes's sense, quoted in [T] §1.5 from the hat paper: "strongly aperiodic if it admits tilings but none with any infinite cyclic symmetry". The stricter sense of Coulbois et al. would need every tiling to have a **trivial** symmetry group.
- **[T] results.** Theorem 7.2: Q is mildly aperiodic. Corollary 7.3: "No tiling by Q has a symmetry of infinite order", so Q is strongly aperiodic in Mozes's sense. §7 also says Q does not meet the stricter sense: "the substitution fixed points of Proposition 8.2 include tilings whose symmetry group has order 8". I could not see Proposition 8.2 or its proof.
- **[HR] §7.4, Theorem 5 and Remark 1.** Sym(T) is conjugate to one of 1, C₂ (a half-turn about a face-diagonal axis), C₄ or D₄, and "each occurs". Remark 1: "|Sym(T)| ≤ 8, and the bound is attained", and no atlas tiling has "a 3-fold axis, a Klein four-group of symmetries, or exactly a 2-fold coordinate axis". The proof (paraphrased in the summary) first excludes elements of order three, then concludes Sym(T) "is a 2-group, contained in a copy of D4". For tilings by the solid, the theorem also depends on an input (R) quoted from [T] and audited but not refereed.

So "only mildly aperiodic" means: **no translation periods and no infinite symmetry, but a finite symmetry group, possibly nontrivial, is allowed.** Order 8 is the largest that occurs.

### Subgroup structure [computed]

The 24 proper cubic rotations form a group isomorphic to S₄. The script enumerates all its subgroups: 30 subgroups in 11 conjugacy classes. Seven classes have no element of order 3, and every one of those lies in a copy of D4 (there are three copies, the Sylow 2-subgroups). So **once order-3 elements are excluded, the bound 8 is automatic**; the order-3 exclusion is the content of [HR]'s argument and I have not reproduced it.

![Subgroup classes of the cubic rotations, and the fixed-point seed sets of the substitution](figures/clock_fig1_symmetry.png)

*Left: the 11 classes. Green are the four listed in [HR] Theorem 5; yellow are the other 2-groups, which [HR] Remark 1 excludes; red contain an element of order 3.*

### Symmetric tilings from the substitution [computed]

A pose (G, t) is a **fixed seed** if it is one of its own children under (G, t) ↦ {(GH, 2t + Gu)}; there are 48 (G, 0) and (G, −G(1,1,1)) for the 24 rotations G. Refining a set of fixed seeds gives a self-similar patch. The script takes all 219 maximal sets of pairwise compatible fixed seeds (no overlap, any contact an atlas contact), refines each three times, and keeps those that overlap nothing, create no illegal contact (E = 0) and fill a cube around the origin. It then measures the symmetry group on that window (rotations about half-integer centres near the origin).

| Seeds | Sets | Fill the window | Symmetry group |
|---|---|---|---|
| 2 | 24 | all 24 | trivial |
| 4 | 144 | none | — |
| 5 | 48 | none | — |
| 8 | 3 | all 3 | order 8: the identity, five half-turns (three about coordinate axes, two about face diagonals) and two quarter-turns, i.e. **D4** |

- One of the 24 two-seed sets is the pose at the origin together with atlas pose 9, (0,1,2)(+,+,+) offset (−1,−1,−1), which is a fixed point of the central child. Its level-n patch has bounding box [−2ⁿ, 2ⁿ⁺¹)³ (computed for n = 1 to 4) and fills the window checked. This is consistent with the fixed point that [T] §6 is summarized as describing, whose two near-origin tiles give supertile unions filling an octant and the closure of its complement; I have not matched my example to [T]'s.
- The result persists at level 4 on a window twice as wide. No proper subset of an 8-seed set fills the window (738 subsets tested).
- The **three** 8-seed sets match the **three** copies of D4 in the rotation group. That is a pattern I noticed in the output, not something I proved.
- **What this does not reproduce.** The C₂ and C₄ cases of [HR] Theorem 5 did not appear. This search covers only substitution fixed points built from seeds that are themselves fixed poses; seeds of higher period, and tilings that are not self-similar, were not searched. The three D4 tilings are an instance of [T]'s "order 8" statement, found independently on my transcription; they are not the tilings of Proposition 8.2.

## 2. Block renormalization and repeated doubling

The earlier comparison table is in [HIERARCHIES.md](HIERARCHIES.md#1-block-renormalization-and-substitution-hierarchies). Here is the part that can be made quantitative.

### The doubling defect in [O3] [source]

For spatial circumference w, [O3] §8.2 defines a positive transfer operator K_w and Z(n, w) = Tr K_wⁿ. The one-direction defects are p(n, w) = 2 log Z(n, w) − log Z(2n, w) and q(n, w) = 2 log Z(n, w) − log Z(n, 2w), and

> Δ(n, w) = 2p(n, w) + q(2n, w) = 4 log Z(n, w) − log Z(2n, 2w) ≥ 0,

so Δ is exactly the **extensivity defect** under doubling both sides (the identity is my algebra from (8.5) and (8.7)). Lemma 8.2 gives p(2n, w) ≤ 4p(n, w)², from the concentration of the transfer trace. Proposition 8.3: if Δ(n, w) ≤ 2⁻¹⁰ (n, w ≥ 4) then for every k ≥ 0,

> Δ(2ᵏn, 2ᵏw) ≤ 64⁻¹ · 2^(−2ᵏ).

The proof shows Δ(2n, 2w) ≤ 64 Δ(n, w)², iterates to 64⁻¹(64·2⁻¹⁰)^(2ᵏ), and weakens it to the stated bound. This is a **quadratic recursion**: the defect is squared (up to a constant) at each doubling. By the text, the proof "uses only" positivity of the transfer operator and the rectangular test.

### The Chair44 analogue [computed]

Count the face contacts C_n of the level-n patch of the one-seed hierarchy (8ⁿ chairs, the corner chair refined n times). The pose-based count matches the cell-based count at levels 1 to 4 (16, 183, 1663, 14079; E = 0 there), and levels 5 to 7 are counted from poses (115 711; 937 983; 7 553 023). All seven values satisfy, exactly,

> C_n = (29/8) · 8ⁿ − 3 · 4ⁿ − 1  (fit on levels 1 to 4, which leaves the 2ⁿ coefficient equal to zero; levels 5 to 7 then predicted).

The volume term is 29/8 contacts per chair. The −3·4ⁿ term scales with the **area** of the supertile, and the constant −1 is the corner. The extensivity defect of this count is exact:

> 8 C_n − C_{n+1} = −(12 · 4ⁿ + 7), so |8C_n − C_{n+1}| / C_{n+1} ≈ (12/29) · 2⁻ⁿ.

![Contact counts and the extensivity defect](figures/clock_fig2_doubling.png)

*Left: computed contacts and the formula. Right: the relative defect of the chair count (blue, geometric) against the bound from the O(3) recursion Δ_{k+1} ≤ 64Δ_k² starting at Δ₀ = 2⁻¹⁰ (red). The horizontal axis is the level for the chairs and the doubling number k (shifted by one) for O(3).*

**Reading it.** [mine]
- Both are "volume term plus boundary term" statements, and both are exact in a sense: a self-similar structure has a defect that can be written down. That is the parallel.
- The rates differ completely: 2⁻ⁿ against a doubly exponential bound. The O(3) defect is a difference of free energies of a positive transfer operator, and the recursion squares it; the contact count is a deterministic count with a boundary term of fixed relative size.
- The sign differs too: Δ ≥ 0 in [O3] (a torus, so no boundary), and the chair defect is negative because the free-boundary count has a deficit.
- So "the same shape of problem" holds for the *structure* (an exact coarse-graining step, an idealized unique-parent version, a boundary-interaction difficulty), and does **not** extend to the quantitative self-improvement. Nothing here transfers a theorem.

### The solvable cousin and non-exhaustion [source, computed]

- [O3] §1.1 contrasts the nearest-neighbour model (blocks interact across boundaries, "whose spatial decay must be controlled") with the hierarchical model (Gawędzki–Kupiainen), where "the hierarchical covariance reduces the flow to a single-spin recursion".
- [T] Theorem 6.4: every tiling partitions uniquely into the eight-chair parent clusters. Theorem 6.5: coarsening gives a tiling by Q again. Theorem 6.6: an intrinsic nested family of supertiles of every level. §6 adds that this "does not say that the supertiles containing one fixed tile grow to cover space" and that they "need not"; "fault lines" are not excluded. Whether every tiling by Q lies in the hull is left open (§10, as the summary reports).
- [computed] The self-similar fixed points of section 1 are the idealized case: their supertiles do fill the window around the origin at every level checked. Generic tilings need not.

### Entry 155's descent [source]

[E155] Proposition 2.4: no array satisfying the line rule, with every column nonconstant, has a positive vertical period. The proof takes the smallest such period M (Theorem 2.3 gives p | M), normalizes the array by a shear, and passes to W₁(n, m) = W(n, pm), which again satisfies the line rule (a line of slope d and intercept e becomes the old line of slope pd and intercept pe) and has period M/p. Minimality is contradicted. This is a rescaling that *shrinks a period*, a third kind of scale step (see the table below).

| | Scale step | What improves |
|---|---|---|
| Chair44 hierarchy | factor 2 (volume ×8) | the coarse object is again a tiling by Q |
| [O3] | doubling (area ×4) | Δ ↦ ≤ 64Δ² |
| [E155] | m ↦ pm | period M ↦ M/p |

## 3. Hilbert–Smith (entry 304)

- [E304] Theorem 1.2 (p-adic exclusion): every jointly continuous action of ℤ_p = lim ℤ/pᵏℤ on a connected Hausdorff second-countable topological n-manifold without boundary has nonzero kernel, so it factors through a finite quotient. Theorem 1.1 then gives that any locally compact second-countable group acting faithfully on such a manifold is a Lie group.
- [HR] §7.4 Corollary 2 puts the internal space of the Chair44 hull at ℤ₂³ (as recorded in [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md); I did not re-fetch that corollary for this page).
- **Why there is no interaction [mine, with background].** Theorem 1.2 restricts actions on *manifolds*. A tiling hull is not a manifold: *(background)* it is locally a Euclidean factor times a totally disconnected one. A 2-adic odometer acting on the hull is therefore consistent with Theorem 1.2, and so is a faithful ℤ₂ action on the Cantor-like part. The shared word is "2-adic".
- If someone wanted a real link, the question would be whether the hull's ℤ₂-action extends equivariantly to some manifold containing it. Theorem 1.2 says such an extension cannot be faithful. I do not know of a source asking this, and it has no bearing on [T]'s results.

## 4. Reconstruction on trees (229) and the clock tests

### What [HR] §10 does [source]

Section 10 is labelled "Computation, not Fact"; the results below are numerical, on finite patches, from the fetch summary of the PDF. The fetch gave no sample sizes, and defines distance only "in tile scales ℓ" without defining ℓ. There are four tests.

| Test | Result as summarized | Conclusion quoted |
|---|---|---|
| (a) Clustering (the Barbour–Koslowski–Mercati scale-invariant measure on tile centroids) | 0.360, 0.539, 0.604, 0.626 for levels 1 to 4; excess over a triangular-lattice control 0.028 at level 1, at most 0.001 after | "Positional clustering carries no arrow." |
| (b) Records (leave-one-out prediction of ancestors from a radius-r neighbourhood) | the level-1 ancestor is "determined at 3ℓ (error 0.000, against 0.151 at 2ℓ)"; the level-2 ancestor is "still wrong in 0.446 of cases at 4ℓ and 0.072 at 8ℓ"; the radius "grows by about a factor of three per level" | finer levels are free, coarser ones "cost a patch three times larger"; an arrow only "if an observer's patch grows in one direction, which the tiling does not supply" |
| (c) Look-back (h = level of the lowest common ancestor of two tiles at distance d) | if depth were look-back time, ⟨h⟩ would grow as log_λ d; fitted slope 0.92 (0.89 for a quadtree with λ = 2 through the same code) | "The dictionary holds on average and is shared by any self-similar hierarchy." |
| (d) Parity (reflection-odd statistics against level) | alternation lives in the shapes of supertiles, not in tiles; a level-k outline is nearer its mirror image than a rotation in "24 of 24 comparisons" | |

Overall (quoted): "None supplies a rule for what propagates, and none supplies a physical clock."

### A Chair44 version of tests (b) and (c) [computed, mine]

I ran analogues on the 8-seed D4 fixed point (a tiling of space, built to level 6, 2 097 152 chairs). Definitions here are mine and differ from [HR]'s, which I do not know in detail:
- A chair's **signature at radius r** is the set of chairs whose centres lie within Euclidean distance r, as relative poses in the chair's own frame. ℓ = the chair's edge, 2 cube units.
- The **label** to predict is the chair's slot in its level-L supertile (8, 64 or 512 possibilities).
- The **error** is the leave-one-out error of majority vote among the other chairs with the same signature, on the 4688 chairs with centres in a cube of radius 16 around the origin. I call a radius "sufficient" when the error is zero and no chair has a signature that no other chair shares.

![Recognizability radius and the look-back statistic for Chair44](figures/clock_fig3_clock_tests.png)

| Level L | Sufficient radius (cube units) | In units of ℓ | Ratio to previous level |
|---|---|---|---|
| 1 | √3 ≈ 1.73 | 0.87 | — |
| 2 | √24 ≈ 4.90 | 2.45 | 2.83 |
| 3 | √131 ≈ 11.45 | 5.72 | 2.34 |

The growth is **2.3 to 2.8 per level**, between the inflation factor 2 and [HR]'s "about three" for the Spectre, whose inflation factor is λ ≈ 2.806. Three levels is too few to separate these. The absolute radii are not comparable with [HR]'s 3ℓ, since the label and ℓ differ.

**Look-back.** For pairs of chairs in four cubes of radius 12 away from the origin (so that the origin's special structure does not matter), inside one level-6 supertile, the fitted slope of ⟨h⟩ against log₂(d/ℓ) over d from 2ℓ to 8ℓ is **0.62** for the Chair44 hierarchy and **0.67** for an octree control (dyadic boxes, same code). Both are below 1 and nearly equal, and the curves flatten as d approaches the window size, so the slope is a finite-window artefact here. This agrees with [HR]'s remark that the statistic "is shared by any self-similar hierarchy"; it does not reproduce their 0.92 and 0.89, which come from a different tile, window and fit range.

### Noisy patches [computed, mine]

The tests in [HR] §10 are noiseless as summarized, and in section 4 of [HIERARCHIES.md](HIERARCHIES.md) I said a noisy version of the clock tests would connect them to reconstruction on trees. Here is a toy version: each chair in the radius-r patch is deleted independently with probability ε, and the label is guessed by the Bayes rule against the empirical set of signatures (in-sample: the dictionary includes the test chairs).

![Kesten–Stigum line for entry 229, and the noisy-patch toy](figures/clock_fig4_noisy_and_trees.png)

| ε | 0.0 | 0.3 | 0.5 | 0.7 | 0.9 |
|---|---|---|---|---|---|
| level-1 slot, radius ℓ (chance 0.125) | 1.00 | 0.95 | 0.85 | 0.68 | 0.39 |
| level-2 slot, radius 2.5ℓ (chance 0.016) | 1.00 | 0.97 | 0.93 | 0.83 | 0.52 |

Accuracy falls smoothly. There is **no threshold** in this toy, because the radius is fixed and the signature is redundant, whereas a reconstruction threshold concerns the depth going to infinity.

### What [E229] proves, and how far the analogy goes [source]

[E229] Theorem 1.1: for the symmetric three-state channel K_λ(i, i) = (1 + 2λ)/3, K_λ(i, j) = (1 − λ)/3 (i ≠ j; −1/2 ≤ λ ≤ 1) on the regular b-ary tree (d = b) or a Poisson Galton–Watson tree of mean d > 1, reconstruction occurs **if and only if dλ² > 1**, and there is non-reconstruction at equality. So the positive threshold is λ = d^(−1/2). On the negative side reconstruction occurs precisely when d > 4 and −1/2 ≤ λ < −d^(−1/2). The abstract's community-detection corollary is "(a − b)² > 3(a + 2b)"; this is the same condition with d = (a + 2b)/3 and λ = (a − b)/(a + 2b), which I checked algebraically.

- **Common ground [mine].** Both ask how much local data far from a position in a tree says about that position. The chair substitution tree has branching 8.
- **Differences.** [E229]'s data are noisy copies of a spin at the leaves and the question is the limit of infinite depth. The clock tests have no noise and ask how large a finite patch must be. The deterministic hierarchy always "reconstructs" with enough radius; the noisy toy has no sharp transition. I know of no source that connects the two beyond this analogy.

## 5. The other connections

- **Anderson model (261) and Parisi/ultrametric (221).** The Anderson paper compares quantities "along a hierarchy of increasingly sparse ports" (introduction; §7); the perceptron papers use Panchenko's ultrametricity and finite-depth hierarchical trial laws (entry 221). A substitution hierarchy does define an ultrametric on tiles (distance = level of the smallest common supertile) [mine]. Anderson's pure-point spectrum concerns a random operator, whereas [HR]'s is pure-point *diffraction*. Detail as in [HIERARCHIES.md](HIERARCHIES.md#5-multiscale-and-ultrametric-arguments).
- **Triangular lattice (090, 232, 237).** [HR] §6 puts Spectre vertices in ℤ[ζ₁₂] (Fact 15) and *(background)* the triangular lattice sits inside it. The repository's papers on it do not mention monotile tilings. No link beyond that.

## Evidence status

| Statement | Status |
|---|---|
| Definitions of mild and strong aperiodicity; Theorem 7.2, Corollary 7.3 of [T]; [HR] Theorem 5 and Remark 1 | [source], via fetch summaries; proofs not checked |
| The 11 subgroup classes, 7 two-group classes, three D4's | [computed], exact |
| 24 + 3 covering fixed-point seed sets; D4 symmetry of the three; no symmetric C₂ or C₄ tilings found | [computed], on a transcribed atlas, within the family searched |
| Contact formula (29/8)8ⁿ − 3·4ⁿ − 1 and defect −(12·4ⁿ + 7) | [computed]; levels 1 to 4 cell-checked, 5 to 7 pose-counted assuming E = 0 |
| [O3] Proposition 8.3 and its quadratic recursion; [E229] Theorem 1.1; [E304] Theorems 1.1, 1.2 | [source], read in the preprint text |
| [HR] §10 numbers | [source], via the fetch summary; not recomputed |
| Chair44 recognizability radii, look-back slopes, noisy toy | [computed], my definitions, not a replication of [HR] |
| The strength of every analogy | [mine] |

## Caveats

- The atlas, the child poses and the pose convention come from a text extraction of [T] and from a fit (see [CHAIR44_MODEL.md](CHAIR44_MODEL.md)). The checks passed, but they test consistency, not correctness against the paper's own files. [T]'s legality test was not re-derived.
- [HR] is unrefereed and its Theorem 5 depends on input (R) from [T]. Its §10 text was available only through a fetch summary, with no definition of ℓ and no sample sizes.
- Pose-based contact counts for levels 5 to 7 assume that every face contact in those patches is an atlas contact, as verified for levels 1 to 4.
- The fixed-point search covers seeds that are fixed poses only. The clock analogues use the D4 fixed point, which is a special tiling. The look-back and noisy computations use small samples (one core of 4688 chairs; four windows of at most 1500 chairs; 1500 test chairs per ε) and have no error bars.
- Chair44 and the Spectre differ; the Spectre was not computed here.

## Open questions

1. Obtain [T] Proposition 8.2 and compare its order-8 tilings with the three D4 fixed points. Search seeds of higher period for the C₂ and C₄ cases of [HR] Theorem 5.
2. Read [T] §7's proof that Q is mildly aperiodic (Theorem 7.2) alongside [HR]'s exclusion of order-3 elements.
3. Repeat the recognizability test with more levels (level 4 and beyond) to separate a factor of 2 from a factor of 3, and on a second tile if Spectre data can be obtained.
4. The 14 atlas poses never realised in hierarchical tilings ([HR] Corollary 1; [CHAIR44_MODEL.md](CHAIR44_MODEL.md)): do they correspond to places where the recognizability test fails?

## Reproducing this page

```
make figures        # about 90 seconds; writes figures/clock_fig*.png and build/data/hierarchy_figures.json
make figures FULL=1 # a larger core for the recognizability test
```

The script asserts its exact results (subgroup counts, the contact formula, the symmetry orders, E = 0 on the fixed points) and stops if one fails.

## Related pages

- [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md), [CHAIR44_MODEL.md](CHAIR44_MODEL.md), [HIERARCHIES.md](HIERARCHIES.md), [MONOTILES.md](MONOTILES.md), [PERIODIC_TILING.md](PERIODIC_TILING.md)
- [tools/README.md](tools/README.md), for how the scripts are organised.
