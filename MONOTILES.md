# Aperiodic Monotiles and Entry 155

A fork-level note connecting entry 155 (a finite tile in ℤ³ with no fully periodic tiling) to the 2026 work on three-dimensional aperiodic monotiles, in particular Tsiokos's Chair44 and the Hartshorn–Romanov preprint. It is not an official OpenAI document. For the proof of entry 155 see [PERIODIC_TILING.md](PERIODIC_TILING.md); for the wider tiling picture see [TILING.md](TILING.md).

**How it was built.** I fetched the [Research Square preprint](https://www.researchsquare.com/article/rs-11516041/v1) (page and PDF) and the arXiv abstract pages of three Chair44 papers, then compared them with entry 155's paper, which I had read earlier. The web fetches were summarised by the fetch tool, so quotes below are short phrases from those summaries, not verified against the full texts. I did not read the Chair44 papers' proofs, and I have not checked any claim. The Chair44 results are new, mostly unrefereed, and some are explicitly "under review" (see Caveats).

## The sources

| Date (2026) | Source | What it is |
|---|---|---|
| 16 Sep (rev. 2) | Tsiokos, [A Strongly Aperiodic Monotile in Three Dimensions](https://arxiv.org/abs/2609.19214), arXiv:2609.19214 | Introduces Chair44. |
| 23 Sep | Entry 155, [A translational tile with no fully periodic tiling in dimension three](preprints/A-translational-tile-with-no-fully-periodic-tiling-in-dimension-three-September-23-2026/) | The ℤ³ counterexample to the periodic tiling conjecture. |
| 30 Sep | Flicker, [Matching Rules for a Three-Dimensional Strongly Aperiodic Monotile](https://arxiv.org/abs/2609.23783), arXiv:2609.23783 | Simpler decorations that force the same structure. |
| undated on its page | Goodman-Strauss, [Notes on a strongly aperiodic monotile in E³](https://arxiv.org/abs/2609.24779), arXiv:2609.24779 | A cleaner "Berger-style" presentation. |
| 28 Oct | Hartshorn and Romanov, [The Internal Space and Boundary Growth of the Spectre and Tsiokos's Chair44 Tilings](https://www.researchsquare.com/article/rs-11516041/v1) (DOI 10.21203/rs.3.rs-11516041/v1) | Exact computation on the Spectre and Chair44 hierarchies. Not peer reviewed. The page body lists Brent Hartshorn and Tristan Romanov; the citation metadata lists only Romanov. |

Entry 155's paper does not cite Chair44, and the Hartshorn–Romanov text I saw does not discuss the periodic tiling conjecture, Greenfeld–Tao, or tiling equations. Every connection below is mine.

## Two different statements

Both results give a single tile in three dimensions with no periodic tilings, but they are not the same kind of result.

| | Entry 155 | Chair44 |
|---|---|---|
| Allowed moves | **Translations only** | Isometries; the paper's registered poses use the 24 proper cubic rotations |
| Tile | Finite subset of ℤ³, enormous and not drawn | A seven-cube chair with 24 exposed faces carrying small pyramid bumps and dents; 2138 mesh vertices, volume 7 |
| What is excluded | Every tiling set with a **finite-index** period subgroup | Every tiling with an **infinite cyclic symmetry** (strongly aperiodic in Mozes's sense, as Tsiokos's abstract uses it) |
| Question answered | The periodic tiling conjecture (Lagarias–Wang, Grünbaum–Shephard) | Existence of a strongly aperiodic monotile in E³ |
| How aperiodicity is forced | Algebra: a p-adic "Sudoku" encoded as tiling equations | Geometry and hierarchy: every tile sits in a unique infinite hierarchy of nested supertiles |
| Tilings exist | Yes (explicit, via last-nonzero digit) | Yes, claimed: Q tiles ℝ³ |
| Status | In this repository, with a Lean proof tree | Tsiokos: Lean 4 for finite checks and assembly "modulo named compiler hooks", geometric arguments "still under review" |

Tsiokos's §1.1 discusses Greenfeld–Tao, notes their dimension is "extremely large" and that ℤ³ was left open, and the paper "does not claim any result on translation-only tiles", so Chair44 does not bear on the periodic tiling conjecture, which is about translational tiles. Section-level citations for this and the other comparisons are in [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md). The two results sit side by side: one shows what translations alone can do, the other what a rotation-allowed shape can do.

A note on strength. The two statements are not comparable by inclusion. Entry 155 excludes only *full* periodicity, and the paper asserts nothing about nonzero periods. The Chair44 claim, if correct, excludes any translational period but allows rotations. Tsiokos's paper also notes that under a stricter taxonomy (Coulbois, Gajardo, Guillon and Lutfalla) Chair44 is only "mildly aperiodic", because some tilings have a symmetry group of order eight.

## Three ways to enforce a local rule

Your "matching rules in the loose sense" point sharpens into a three-way comparison. All three make a finite local condition force a global property.

| Method | Example | What carries the rule |
|---|---|---|
| **Shape only** | Chair44 | Bumps and dents. Tsiokos states that "no matching rule of any kind is imposed"; the 44 legal contacts are a theorem about what the shape permits, not a rule added to it. |
| **Explicit decoration** | Goodman-Strauss's marked L-tile; Flicker's arrow and face-colour rules | Three colours that must match in a fixed pattern force a hierarchy. Goodman-Strauss calls the resulting proof "Berger-style". Flicker gives arrow rules equivalent to the originals and face-colour rules (black meets only white, blue meets only blue) with 611 meeting configurations instead of 44; mirror-image tilings are excluded at the shell-of-shells level. |
| **Algebraic encoding** | Entry 155 | Several tiling equations with one unknown translation set (Sections 3–4 of the paper), encoding a p-adic word rule. No geometry is added to the tile. |

Flicker's paper asks for the minimal decoration needed "to ease physical realisation" and argues, without proof, that simpler colour rules likely do not exist. Entry 155 sits at the opposite end: nothing about it is physical or small.

## Parallels worth noting

These are my observations, not claims of either source.

1. **A scaling-by-a-prime argument on both sides.** Entry 155 proves non-periodicity by descent: rescaling the rows m ↦ pm turns a vertical period M into M/p, contradicting minimality (Proposition 2.4). Chair44's hierarchy is an inflation: the Hartshorn–Romanov preprint says it inflates by 2I₃ on ℤ³ (determinant 8) with no twist, and that its internal space is **2-adic** (ℤ₂³) rather than Euclidean. In both, a prime scale and p-adic digits control the structure. Tao's blog post (see [PERIODIC_TILING.md](PERIODIC_TILING.md)) adds that the 2022 Greenfeld–Tao construction also moved to a 2-adic Sudoku, with a large power of 2. The three 2-adic or p-adic features are different objects, and I am not claiming they are related.
2. **Hierarchy versus descent.** A supertile hierarchy and an infinite descent are the same shape of argument: a structure that reproduces itself at a larger scale cannot be periodic. Chair44's induction (Superchairs obey the same rules as Chairs) and entry 155's minimal-period argument both end in "a smaller or larger copy of the same thing exists".
3. **Dimension three and the role of isometries.** Entry 155 is minimal at dimension three because translational tiles of ℤ and ℤ² always tile periodically. Allowing rotations and reflections changes the picture in the plane: the Hartshorn–Romanov preprint treats the planar **Spectre** as an aperiodic monotile and compares it with the three-dimensional Chair44 in a section titled "From Plane to Bulk". So the dimension at which a single tile first forces aperiodicity depends on which motions are allowed *(background for the planar case)*.
4. **Verification culture.** Entry 155 ships with a Lean proof tree ([PERIODIC_TILING.md](PERIODIC_TILING.md)). Tsiokos's paper puts its finite checks and logical assembly in Lean 4, modulo compiler hooks, and Goodman-Strauss's notes say the discovery was reported "with the help of LLM tools". Hartshorn and Romanov use fixed seeds, finite certificates and an audit of a registration theorem that "has not been refereed". The same questions arise for each: what is machine-checked, and what is cited.
5. **Computation as evidence.** The Hartshorn–Romanov paper separates *proved* results (Kamenetsky's unique-hierarchy criterion by finite certificate, every atlas tiling lies in the hull, symmetry groups conjugate to 1, C₂, C₄ or D₄), *computed* ones (360,189 admissible blocks in 15,214 orbits, 2-adic counts, the Spectre perimeter recurrence) and *cited* ones (Tsiokos's registration theorem). That is the same hygiene this fork uses elsewhere, for example the checks in [tools/](tools/README.md).

## Where physics enters

Entry 155 has no physical content. The Chair44 literature does: Flicker's paper suggests physical systems that might realise these tilings; Hartshorn and Romanov compute a pure-point diffraction measure for the Chair44 centre set (Bragg masses summing to the density 1/7) and run four finite-patch tests of whether the Spectre's hierarchy could act as a "clock". They report no arrow of time in positional clustering. See [PHYSICS.md](PHYSICS.md) for the physics-related papers that are in this repository; none of them is about aperiodic order.

## What does not connect

- Chair44 does not refute or imply the periodic tiling conjecture, and entry 155 does not bear on the existence of a strongly aperiodic monotile.
- Entry 155 does not use hierarchical substitution, and the Hartshorn–Romanov preprint does not use tiling equations. I found no common technique beyond the loose parallels above.
- I did not confirm what the paper's "lattice-registered" means: that phrase did not appear in the text I saw. The closest terms are "registered" and "atlas tiling", defined there as tilings by chairs in poses with integer translation parts in which every face-adjacent pair has a relative pose among the 44 legal contacts.

## Possible next steps

- A `tools/` script that generates the last-nonzero-digit array of entry 155 for small p and checks the line rule, next to a script that counts Chair44 legal contacts or inflation addresses, as a hands-on comparison. This is not built.
- Read the Chair44 papers' proofs, especially Kamenetsky's criterion and the registration theorem, before citing any of their statements as established.
- Ask whether a shape-only or colour-matching tile can be rephrased as a system of tiling equations, which would link the two traditions. I know of no source that does this.

## Caveats

- Everything about the Chair44 papers is from abstract pages and fetch-tool summaries. I could not open Tsiokos's last ~2,600 characters of text, and the Hartshorn–Romanov PDF summary shows no reference list with arXiv numbers.
- The Hartshorn–Romanov preprint says Tsiokos's einstein property and aperiodicity claims "remain Tsiokos's claim", and calls its own registration premise "audited but not refereed".
- The search that found these papers also returned a Quanta Magazine update and several news re-posts about Chair44. I did not open them.
- The parallels section is interpretation. Neither source states any of it.

## Related pages

- [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md) for the section-by-section citations.
- [PERIODIC_TILING.md](PERIODIC_TILING.md), [TILING.md](TILING.md), [PHYSICS.md](PHYSICS.md).
