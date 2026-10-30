# The Spectre, Chair44 and Entry 155

A fork-level, citation-level comparison of three single-tile results: the planar **Spectre**, the three-dimensional **Chair44**, and this repository's **entry 155** (a finite tile in ℤ³ with no fully periodic tiling). It is not an official OpenAI document. It expands the four connections summarised in [MONOTILES.md](MONOTILES.md); for entry 155's proof see [PERIODIC_TILING.md](PERIODIC_TILING.md).

**How it was built.** I fetched each source below and asked for its section and theorem numbers. The fetch tool summarises pages, so the section numbers and the short quotes come from those summaries and I have not re-checked them against the full texts. Where a number was not in the text I saw, I say so. I did not check any proof. The Chair44 results are new and mostly unrefereed (see [Evidence status](#evidence-status)).

## Citation key

| Key | Source |
|---|---|
| **[E155]** | OpenAI, [A translational tile with no fully periodic tiling in dimension three](preprints/A-translational-tile-with-no-fully-periodic-tiling-in-dimension-three-September-23-2026/paper.pdf), 23 Sep 2026 (entry 155). Section numbers as in the paper. |
| **[T]** | I. Tsiokos, [A Strongly Aperiodic Monotile in Three Dimensions](https://arxiv.org/abs/2609.19214), arXiv:2609.19214, revision 2 of 16 Sep 2026. |
| **[GS]** | C. Goodman-Strauss, [Notes on a strongly aperiodic monotile in E³](https://arxiv.org/abs/2609.24779), arXiv:2609.24779. Unnumbered theorem; I cite its sections 1 (Background) and 2 (The tile). |
| **[F]** | F. Flicker, [Matching Rules for a Three-Dimensional Strongly Aperiodic Monotile](https://arxiv.org/abs/2609.23783), arXiv:2609.23783v2, 30 Sep 2026. Unnumbered sections; I cite by heading and figure. |
| **[HR]** | B. Hartshorn and T. Romanov, [The Internal Space and Boundary Growth of the Spectre and Tsiokos's Chair44 Tilings](https://www.researchsquare.com/article/rs-11516041/v1), Research Square preprint, DOI 10.21203/rs.3.rs-11516041/v1, 28 Oct 2026. Not peer reviewed. |
| **[SM]** | Smith, Myers, Kaplan and Goodman-Strauss, [A chiral aperiodic monotile](https://arxiv.org/abs/2305.17743v2), arXiv:2305.17743 (the Spectre). The authors are as cited in [HR] refs [1, 2]; I read only the abstract. |
| **[GT22]** | Greenfeld and Tao, [A Counterexample to the Periodic Tiling Conjecture (Announcement)](https://arxiv.org/abs/2209.08451), arXiv:2209.08451, Sep 2022. |
| **[Tao]** | T. Tao, [blog post, 29 Nov 2022](https://terrytao.wordpress.com/2022/11/29/a-counterexample-to-the-periodic-tiling-conjecture-2/), on the full paper (arXiv:2211.15847, Annals 2024, which I could not fetch). |

## The objects

| | Spectre | Chair44 | Entry 155 |
|---|---|---|---|
| Ambient space | Plane | ℝ³ | ℤ³ (and ℝ³ by thickening) |
| Allowed motions | Translations and rotations; reflections neither needed nor permitted ([SM] abstract) | Isometries including reflections; registered poses use the 24 proper cubic rotations ([T] §1.2, §3.1) | Translations only ([E155] Theorem 1.1) |
| The tile | A modification of the edges of the equilateral member of the hat's continuum ([SM] abstract) | A 2×2×2 block of cubes with one corner removed (seven cubes, volume 7), with 8 small square pyramids on each of the 24 exposed unit panels, 192 in all ([T] §2.1, §2.3, §2.4) | A finite subset of ℤ³, huge and not drawn ([PERIODIC_TILING.md](PERIODIC_TILING.md), "Size of the tile") |
| Inflation | λ² = 4 + √15 in ℤ[√15] ([HR] §6, Fact 16) | Factor 2, matrix 2I₃, determinant 8 ([T] §3.1; [HR] §6) | None; a rescaling m ↦ pm is used inside the proof ([E155] Prop. 2.4) |
| Twist per level | θ_S = 11.565° ([HR] §8.1, Fact 6) | 0 ([HR] §8.1, Fact 19) | Not applicable |
| Internal space | Euclidean plane inside a 4D vertex module ℤ[ζ₁₂] ([HR] §6, Facts 15 and 17) | 2-adic, ℤ₂³; no Euclidean internal space ([HR] §6, §7.3, §7.4 Corollary 2) | Not discussed |

## 1. They are different results

**What each theorem says.**
- [E155] Theorem 1.1: a finite T ⊂ ℤ³ tiles ℤ³ by translations, but **no tiling set is invariant under a subgroup of finite index**. The thickening Ω = T + [0,1]³ has the same property in ℝ³ even with arbitrary real translations. The paper states that "neither connectedness of the tile nor absence of every individual nonzero period is asserted" (end of its introduction).
- [T] Theorem 1.1: "The solid Q ... is a strongly aperiodic monotile of ℝ³: Q admits tilings of ℝ³ by congruent copies, and no such tiling has a symmetry of infinite order." Theorem 1.2: a tiling exists; every tiling T has Per(T) = {0} and |Sym(T)| ≤ 24; every tiling is homochiral; each tile lies in a unique infinite hierarchy of 2ⁿ-scaled chair supertiles in registered poses. Isometries include reflections (§1.2); Per and Sym are defined in §1.5.
- "Strongly aperiodic" is used in Mozes's sense, quoted from the hat paper: "a set of tiles is strongly aperiodic if it admits tilings but none with any infinite cyclic symmetry" ([T] §1.5; [GS] footnote 2: "no tessellation has any period whatsoever").

**Why Chair44 does not bear on the periodic tiling conjecture.** [T] §1.1 discusses Greenfeld and Tao, says their counterexample is for a dimension d that "will be extremely large", and records that whether ℤ³ admits such a tile is left open. It also notes Greenfeld and Kolountzakis made the tile connected. The paper "does not claim any result on translation-only tiles." Chair44 is a shape that may be rotated and reflected, so it answers a different question. Entry 155, dated a week after [T]'s revision 2, answers the ℤ³ question that [T] §1.1 leaves open. Neither paper cites the other.

**The strengths are not comparable by inclusion.**
- [E155] excludes only *full* periodicity, for translations only.
- [T] excludes any period at all, but allows rotations.
- [T] §1.1 and §1.5 add that under the stricter taxonomy of Coulbois, Gajardo, Guillon and Lutfalla, a tile all of whose tilings have finite symmetry group is "mildly aperiodic", and "strongly aperiodic" is reserved for trivial stabilizers. Chair44 is then only mildly aperiodic, because an explicit tiling has a symmetry group of order eight (§7). [HR] §7.4 Theorem 5 gives the same picture: Sym(T) is conjugate to 1, C₂ about a face diagonal, C₄ or D₄, so |Sym(T)| ≤ 8 (Remark 1).

## 2. Three ways to force a global property from a local rule

| Method | Source | Rule, with citation |
|---|---|---|
| **Shape only** | Chair44, [T] | "no matching rule of any kind is imposed" ([T] §1.1, last paragraph). Hypotheses: "No face-to-face, lattice, common-orientation, connected-contact-graph or local-finiteness hypothesis is imposed" ([T] Theorem 1.2, §1.2); local finiteness is proved, not assumed (Lemma 4.1). Of 2,388 candidate face-neighbour poses, "Exactly 44 ... are legal", a Lean theorem proved with `native_decide` ([T] §3.1); [HR] §7.1 lists 1194 of the 2388 as mirror images. Every tiling is registered ([T] §5, Theorem 5.6). [T] §1.5 defines an aperiodic monotile as one with no "non-geometric matching rule constraining tile adjacencies". |
| **Explicit decoration** | [GS], [F] | [GS] §2: the L-tile (a large cube with a half-sized cube deleted from a corner, §1) carries three markings, red, green and blue, in a fixed cyclic order at each of its seven corners and at the front socket. "Red may only meet green, and green may only meet red"; blue arrows must match. [F]: three arrow colours (black and white along square edges, blue diagonal), and "A square may meet another square only if their arrows align" (Fig. 2a); or three face colours, with black meeting only white and blue only blue (Chair611, 611 contacts, 310 involving a mirrored Chair). |
| **Algebraic encoding** | Entry 155, [E155] | The tile carries no decoration. The rule lives in the equations on the unknown translation set A: the graph tile (§3.2), dependence tiles (§3.3, Lemma 3.1), cycle tests (§3.4, Lemma 3.2), activation tiles (§4.2), then stacking into one tile (§6, Proposition 6.2). The word rule they encode is the p-adic Sudoku of §2. |

**The first two are interconvertible in principle.** [F] Fig. 4b–d give *geometric encodings* of the face colours, and Fig. 4c and 4d were proposed by J. Walton. Chair44's pyramids are such a geometric encoding of Chair44's arrow rules. [F]'s reading is that three face colours "alone suffice to force the aperiodic tilings", with mirror images excluded at the shell-of-shells level. The counts there: 33 one-layer shells with arrows against 49 with face colours (the 16 new ones each contain a mirrored Chair), and only the same 15 shells appear in infinite tilings. [F]'s Discussion conjectures that simpler matching rules do not exist and gives some supporting counts: letting blue meet any colour gives more than 20,000 one-layer clusters, and merging black with white while keeping arrow directions gives 6400. Second-level clusters were not enumerated for those variants, so the conjecture is open.

**The third is different in kind.** In entry 155 the rule is not on the boundary of the piece. Tile and rule are both encoded in a finite set F ⊂ ℤ² × (finite cyclic group), and the "matching" happens in how translates of F cover a group. That is why the tile has no picture.

The [GS] argument is the classical one: the markings force a hierarchy because "supertiles can only meet fully face-to-face, and they may only meet in the same way that the tiles do", and the corner markings of a supertile equal those of the inner tile ([GS] §2). He describes it as a "reasonably simple Berger-style proof". [F]'s arrow rules force Chairs into Superchairs of doubled lengths, with the Superchair's arrows uniquely determined by the Chair's (Fig. 2c, 2d), and "slides" are forced to be multiples of 2, then 4, and so on, leaving only the zero slide.

## 3. The proof structures parallel each other

This is my observation. None of the sources states a link.

| | Entry 155 | Chair44 |
|---|---|---|
| Renormalisation map | Given an array W obeying the line rule, W₁(n, m) = W(n, pm) also obeys it ([E155] Prop. 2.4) | Halving parent origins and scale gives another tiling by Q, with D(T + v) = D(T) + v/2 ([T] Theorem 6.5, "coarsening") |
| Scale factor | A prime p > 200 ([E155] §2) | 2 ([T] §3.1, Table 1) |
| Structure above one level | A vertical period M must be divisible by p ([E155] Lemma 2.3) | Each tiling partitions uniquely into eight-chair parent clusters ([T] Theorem 6.4), with a nested supertile family at all levels ([T] Theorem 6.6) |
| How it ends | M/p is a smaller vertical period of W₁, contradicting minimality of M | Per(T) = {0} ([T] Theorem 1.2) |
| Rotation during rescaling | None: the normalisation step shifts rows by un + v, a shear | None: "the twist is 0"; inflation is a scaling by 2 with no rotation or mirror, though child rotations are 0°, 90°, 120° and 180° ([HR] §8.1, Fact 19) |
| Arithmetic | The line rule uses the last nonzero p-adic digit, f(t) = t / p^ord_p(t) mod p ([E155] Lemma 2.5) | Internal space ℤ₂³; hull uniquely ergodic with pure-point spectrum ([HR] §7.4 Corollary 2); modular coincidence at depth 3 for 336 of 512 addresses ([HR] §7.3) |

**Where the 2-adic thread appears.**
- [Tao]: the 2022 Greenfeld–Tao construction switched from an odd-prime p-adic Sudoku to a **2-adic** Sudoku, with a large power of 2, which the authors believe to be 2¹⁰, because the 2-group-valued functions and the ℤ/pℤ "clock" functions could not be made to interact.
- [E155]'s introduction says its central construction instead keeps the finite factor *cyclic*, and its Sudoku uses an odd prime p > 200 (§2, §3.1: all finite factors have pairwise coprime orders).
- [HR] §7.3 and §7.4 put the Chair44 hull's internal space at ℤ₂³.

**What the analogy does not say.** Three different things carry the "2-adic" or "p-adic" label. In [E155] the arithmetic is in the *rule* that the tile must encode. In [Tao] it is in the construction of a tile in ℤ² times a finite group. In [HR] it is an *invariant* of an existing substitution hull, an infinite compact group, not a finite one. Hierarchy and descent are the same shape of argument, a structure that maps solutions to rescaled solutions cannot be periodic. Nothing here shows they are the same mechanism.

## 4. Allowed motions change the minimal dimension

| Setting | Result | Source |
|---|---|---|
| Translations, finite tile in ℤ | Every tiling is periodic | [E155] §1 (Newman) |
| Translations, finite tile in ℤ² | A periodic tiling always exists | [E155] §1 (Bhattacharya); Greenfeld–Tao give a quantitative version and a structure theorem |
| Translations, finite tile in ℤ³ | No periodic tiling for a specific tile | [E155] Theorem 1.1 |
| Translations, large ℤ^d | False, via ℤ² × G₀ with G₀ a finite abelian group | [GT22] abstract; [Tao] (2-group) |
| Plane, rotations and translations | The Spectre: strictly chiral, tilings use rotations and translations without reflected tiles, hierarchical substitution | [SM] abstract; [HR] §2, §6 |
| Plane, with reflections | The hat mixes unreflected and reflected tiles in all its tilings | [SM] abstract |
| ℝ³, isometries | Chair44, claimed strongly aperiodic | [T] Theorem 1.1 |

**Reading the table.**
- With translations only, a finite lattice tile first forces non-periodicity in dimension three. Allowing rotations lets a single planar tile (the Spectre) force aperiodicity already in two dimensions.
- The Spectre's tilings are non-periodic outright. The ℤ³ tile of entry 155 is only shown to avoid *full* periodicity.
- **The two settings are not the same question.** Entry 155 is about finite subsets of a lattice, or their unit-cube thickenings. The Spectre is a planar shape under rotations and translations. [E155] §1 also says the positive planar results "do not establish the conjecture for arbitrary bounded measurable Euclidean tiles", so the Euclidean plane is not settled for translations alone.

**The Spectre and Chair44 side by side ([HR]).**
- Boundary: Chair44's supertile has 24·4ᴺ boundary faces and boundary dimension D_b = 2 = d − 1, so no fractal excess ([HR] §7.2). The Spectre's supertile boundary is fractal, with box dimension about 1.3993 ([HR] abstract), and its window has dimension log(5 + 2√6)/log(4 + √15) ≈ 1.111 ([HR] Theorem 3, §6).
- Perimeter: the Spectre's supertile perimeter is 2F₃ₖ₊₅ + 4, growing by the unit 2 + √5 of ℚ(√5) ([HR] Theorem 1; I did not see its section number).
- Diffraction: the Chair44 centre set has an exact pure-point diffraction on dyadic wave vectors, with Bragg masses summing to the density 1/7 ([HR] §9, Fact 22). The Spectre's pure-point diffraction is quoted from the literature ([HR] refs [8, 9]).
- Chirality: the undecorated chair is achiral and decorated Chair44 is chiral, with no mirror contact legal ([HR] §7.6, Fact 18). [T] Theorem 7.4 and Corollary 7.5 give Chair44 as a strictly chiral aperiodic monotile in the Spectre paper's sense.

## Evidence status

| Result | What is machine-checked | Review status |
|---|---|---|
| Entry 155 | A Lean proof tree in this repository, 122 files; my search found no `sorry`, `admit` or `axiom`. I did not build it. | Release by OpenAI; no outside review noted here. |
| Chair44 [T] | "Finite checks and the logical assembly are Lean 4 theorems, kernel-checked modulo the named compiler hooks" (abstract, §1.2); `native_decide` adds a trust hook per theorem, 21 named hooks for the main theorem, "No theorem in the development depends on sorryAx" ([T] §8.1). | Presented as "a proof submission" until refereed (§1.2). Reviewed adversarially twice; "the reviewer could not run Lean"; a third round was requested on 2026-09-10 and its verdict was pending (§8.1). Geometric arguments "still under review". |
| Notes [GS] | No Lean. | The author finds Tsiokos's paper hard to follow and says it has "hallmarks of an LLM-driven research paper" (abstract). He says Tsiokos "discovered the Chair44 monotile through prompts to an LLM" (§1) and that "No LLMs were used in preparation of this paper" (Acknowledgements). |
| Matching rules [F] | Computer enumeration of contacts and shells. | Preprint. Minimality is a conjecture. |
| [HR] | Fixed seeds; finite certificates. Proved there: Lemma 1 (Kamenetsky's unique-hierarchy criterion), Theorem 4 (every atlas tiling lies in the hull), Theorem 5, Corollaries 1 and 2. Computed: 360,189 admissible blocks in 15,214 orbits (Lemma 2). | Not peer reviewed. Registration premise (R) is quoted from [T] Theorem 5.6 and Lemma 5.3, audited by the second author "but not refereed". The einstein property remains Tsiokos's claim. |

[HR] Corollary 1 notes that only 30 of the 44 legal contacts occur in any tiling, a point [F] also makes (footnote 15 there).

## Ideas and questions

Each is mine. None is claimed by a source.

1. **Poses form a group.** [T] §3.1 defines a registered pose as one of the 24 proper cubic rotations plus an integer translation. Registered tilings are therefore subsets of a *non-abelian* group, the integer translations extended by the cubic rotation group. Whether Chair44 can be treated as a single tile in that group is a question I cannot answer. Greenfeld and Tao have a paper whose title mentions "one nonabelian tile" ([E155] reference [4]); I have not read it and do not claim it addresses this.
2. **Does a hierarchy yield a Sudoku?** Entry 155's rule is a finite constraint on words along lines, and Chair44's is a finite atlas of contacts. Whether Chair44's atlas can be rewritten as a line rule with a descent like Proposition 2.4 is unknown to me.
3. **A hands-on comparison.** A script that prints the last-nonzero-digit array for small p and checks the line rule, next to a count of Chair44 legal contacts or inflation addresses from [T] §3.1 and [HR] §7.3. Not built.
4. **Read before citing.** Kamenetsky's criterion and Tsiokos's registration theorem carry much of [HR]'s Chair44 claims. They should be read directly.

## Caveats

- Section and theorem numbers come from fetch-tool summaries of the PDFs. Numbers I did not see are marked as such, for example the statements of [T] Theorems 7.1 and 7.2, which I did not see.
- I could not read the last ~2,600 characters of [T]. The part I did read has no reference list or acknowledgements.
- I could not open the Greenfeld–Tao Annals paper (arXiv:2211.15847), Kamenetsky's criterion, or [F]'s "Note added". The Spectre and hat facts come from [SM]'s abstract and [HR]'s citations.
- The hat paper's arXiv number (2303.10798) and Khare's paper (arXiv:2212.02429) came from comments under [Tao] and are not used for claims here.
- Section 3 and the Ideas section are interpretation.

## Related pages

- [MONOTILES.md](MONOTILES.md), the shorter comparison this page expands.
- [PERIODIC_TILING.md](PERIODIC_TILING.md), the proof of entry 155 and its Lean tree.
- [TILING.md](TILING.md), the wider tiling, dimer and packing context.
