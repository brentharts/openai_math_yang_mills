# The Periodic Tiling Conjecture and Local Rules

A fork-level guide to entry 155, the finite tile in ℤ³ that tiles space but admits no periodic tiling, and to the idea that makes its proof work: *a finite set of local rules that forces a global property*. It is not an official OpenAI document. For the wider tiling, dimer and packing picture see [TILING.md](TILING.md).

**How it was built.** I extracted the text of the [paper](preprints/A-translational-tile-with-no-fully-periodic-tiling-in-dimension-three-September-23-2026/paper.pdf) and read the introduction and Sections 2, 3, 4 and 6 in full, plus Sections 5.1, 5.4, 7.3 and 7.4. I skimmed the rest (the proofs of Sections 5.2 to 5.3, 7.1 and 7.2). I did not check proofs line by line. I also read the Lean statement and listed the Lean proof tree. I did not build the Lean code. Statements marked *(background)* are general knowledge, not from this repository.

## The question

Take a finite set F in an abelian group G, such as ℤ³. A set A of translations *tiles* with F when every element of G is written uniquely as a + f with a in A and f in F, written A ⊕ F = G. A tiling set A is *fully periodic* when it is invariant under a subgroup of finite index, so it repeats along a full lattice of translations. The **periodic tiling conjecture** asks whether a finite tile that tiles must also tile periodically. (The paper cites Lagarias and Wang for an explicit Euclidean version.)

| Dimension | Status, as stated in the paper |
|---|---|
| ℤ | Every tiling by a finite tile is periodic (Newman, finite-state argument). |
| ℤ² | Every finite tile admits a periodic tiling (Bhattacharya). Greenfeld and Tao added a quantitative version and a structure result. |
| ℤ³ | **Entry 155: false.** A finite tile exists that tiles, but with no fully periodic tiling. |
| large ℤ^d | False, by Greenfeld and Tao, via a tile in ℤ² × G₀ with G₀ a finite abelian 2-group. Announced in [arXiv:2209.08451](https://arxiv.org/abs/2209.08451) (September 2022); full paper [arXiv:2211.15847](https://arxiv.org/abs/2211.15847) (Annals of Mathematics, 2024). |

Entry 155 therefore moves the counterexample to the smallest possible lattice dimension.

## The theorem

**Theorem 1.1.** There is a finite nonempty T ⊂ ℤ³ such that:
1. T tiles ℤ³ by translations, but no tiling set A is invariant under a finite-index subgroup;
2. the thickening Ω = T + [0,1]³ tiles ℝ³ by translations up to null sets, but no such tiling has a translation set invariant under a full-rank lattice, **even with arbitrary real translation vectors**.

Because (i) rules out every periodic tiling set, *every* tiling by T is non-periodic in this sense. The paper asserts neither that T is connected nor that no tiling of T has a nonzero period.

## Local rules that force a global property

The idea behind a matching rule is that checking only small neighbourhoods can force a property of the whole tiling. This proof does the same thing, with algebra in place of coloured edges.

**Stage 1: a Sudoku-type word rule (Section 2).** Fix a prime p > 200, let N = p², and let Σ = F_p^× be the nonzero residues. An array W(n, m) with n in {0, …, N−1} and m in ℤ obeys the *line rule* if, along every line of integer slope d and intercept e, the word n ↦ W(n, dn + e) is an *allowed word*. A word w is allowed if there are integers a, b, not both divisible by p, with

- w(n) = an + b (mod p) when p ∤ an + b,
- w(n) = (an + b)/p (mod p) when p ∣ an + b but p² ∤ an + b,
- no constraint when p² ∣ an + b.

So a line looks like a linear function mod p, except where that function vanishes, where the *next p-adic digit* takes over. This is a finite constraint on words over a finite alphabet, which is why it can be encoded in a tile.

*Why it forbids periodicity.* Lemma 2.2 shows any array obeying the line rule agrees with a single affine form An + Bm + C wherever that form is nonzero. If every column is nonconstant then B ≠ 0, and then any vertical period M must be divisible by p. Now rescale the rows, m ↦ pm. The new array still obeys the line rule and still has nonconstant columns, but has the smaller period M/p. Choosing M minimal gives a contradiction. This is a descent argument.

*Why the rule is satisfiable.* The model array is W(n, m) = f(m), where f(t) is the *last nonzero p-adic digit* of t, that is t/p^ord(t) mod p. Lemma 2.5 shows it obeys the line rule and has nonconstant columns.

A small instance I computed (not from the paper): for p = 3, f(1), f(2), … = 1, 2, 1, 1, 2, 2, 1, 2, 1, 1, 2, 1, … Any candidate period fails, for instance f(3) = 1 but f(6) = 2, and f(9) = 1 but f(18) = 2. For p = 3, 5 and 7 (with N = p²) I also checked that the model satisfies the allowed-word condition on every line with |d| ≤ 2p² and a window of intercepts, and that no M below 3000 is a period of f over a window of ±3000. That is a sanity check of Lemma 2.5 on a finite window, not a proof.

**Stage 2: turn the rule into tiling equations (Sections 3 and 4).** The unknown is one translation set A in G = ℤ² × V, where V is a finite *cyclic* group. Several tiles F_ν are imposed at once, each as an equation A ⊕ F_ν = G with the *same* A. Each tile tests one local property of A:

| Tile | Section | What it enforces |
|---|---|---|
| Graph tile H₀ | 3.2 | A contains exactly one point above each base point (x, k), so A is the graph of "output" functions. |
| Dependence tiles | 3.3 | A channel's output c_i does not depend on a given input: A ⊕ F_{i,σ} = G exactly when c_i(b − σ) = c_i(b) for all b. |
| Cycle tests | 3.4 | A list of channels with chosen labels cannot all be "active" at the same point. Used to ban every word not allowed by the line rule. |
| Activation tiles | 4.2 | Every ordinary channel is active at every point, and each of two *seed* channels is active somewhere. |

A symbol is encoded by *which input digit an output depends on*. Each of the N ordinary channels carries one symbol of the word at that position, and a "seed" forces two columns to take two different values, so no column of the extracted array is constant.

The forcing step is a counting argument. The activation tile uses a multiplicity function μ on a product of two prime-order groups with values 2, 0 and 1. Because μ is not constant, exact coverage cannot hold if a channel has no active digit (Lemma 4.2). For the seeds, exact coverage reduces to a convolution identity on a finite group. That identity forces a prime a_i, which is 2 or 3 for the two seeds, to divide p⁴, which is impossible (Lemma 4.3).

**Stage 3: decode a periodic tiling back into the rule (Proposition 4.4).** From any common solution A, read off an array W(n, m) by taking the smallest symbol active at (0, m). The cycle tests show W obeys the line rule, and the seeds show its columns are nonconstant. If A were invariant under a finite-index subgroup, some multiple of the vertical step (0, 1) would be a period of A, and so W would have a vertical period. That contradicts Proposition 2.4. Hence **no common solution is fully periodic**.

**Stage 4: a common solution exists (Section 5).** The last-nonzero-digit array supplies an explicit A. Block labellings of the finite residue digits and a Chinese-remainder choice of offsets let every bijection in the activation criterion be checked directly. The exclusion tests hold because the one active ordinary word at each x is f(L_n(x)), which is allowed.

**Stage 5: one tile instead of many (Section 6).** Take a prime q and a partition ℤ/qℤ = E₁ ⊔ … ⊔ E_s with E_ν − E_ν = ℤ/qℤ for every ν (found by random colouring). The stacked tile F^st = ∪ F_ν × E_ν tiles exactly when all the equations have a common solution, and a periodic tiling by F^st gives a periodic common solution and conversely. With q a fresh prime, the finite factor stays cyclic, so the result is one tile in ℤ² × ℤ/Qℤ.

**Stage 6: from ℤ² × ℤ/Q to ℤ³ and ℝ³ (Section 7).** The group ℤ² × ℤ/Qℤ is a quotient of ℤ³ via π(a₁, a₂, a₃) = (a₁, a₂, a₃ mod Q), with kernel generated by w = (0, 0, Q). The paper lifts the tile to ℤ³ and then builds two shapes with a large common part and one moved marked cube. For real translations it shows: the common part forces all piece centres onto one translated grid (Lemma 7.4); the marked cube then forces the tiling to be invariant under w; and the tiling descends to a periodic tiling of ℤ² × ℤ/Qℤ, contradicting Stage 5. This handles arbitrary real translation vectors, which integer vertices alone would not.

## Compared with Wang tiles and Penrose tilings

| | Wang tiles, Penrose *(background)* | Entry 155 |
|---|---|---|
| Pieces | Several prototiles, with edge or decoration rules | **One** finite tile |
| Allowed moves | Often rotations, reflections or tile-edge matching | **Translations only** |
| Rule enforced by | Matching edges, so local compatibility | A system of tiling equations on one translation set |
| Conclusion | Tilings exist, none is periodic | A tiling exists, **none is fully periodic** |
| Dimension | Plane | Dimension three (and ℤ² is known to be periodic) |
| Explicit picture | Small, drawable | None practical; see below |

The paper does not mention Wang tiles, Penrose tilings, or matching rules by name, and I found none in the summaries. The comparison above is mine. A reasonable one-line description is that the tile of entry 155 is an *algebraically encoded* local rule, where Wang and Penrose tilings use *geometrically drawn* ones. The paper lists the Greenfeld–Tao papers on undecidability of translational tilings among its references. I have not read them, and entry 155 does not claim an undecidability result.

## Size of the tile

The tile is not small. In Section 6, |F| = q · |H₀|, where |H₀| = D · ∏ r_i |P_i| and D = p⁴. Since p > 200 is prime, p ≥ 211 and D ≥ 211⁴ = 1,982,119,441. There are N + 2 = p² + 2 channels, at least 44,523, and each contributes a factor r_i |P_i|. Remark 7.6 says the data can be chosen by terminating finite searches and that "no practical bound on the size of the finite construction is required". So the tile exists and is finite, but no picture of it is given, and I would not expect one.

## The Lean formalization

- Statement: [PeriodicTilingThree.lean](lean/ComparatorChallenges/PeriodicTilingThree.lean) (it ends in `sorry` because it is the challenge to be solved) with its [Comparator config](lean/ComparatorChallenges/PeriodicTilingThree.json). The config points at the solution module `OAI.Geometry.PeriodicTiling.Main`.
- Scope note: [lean/docs/155.md](lean/docs/155.md).
- Proof tree: [lean/OAI/Geometry/PeriodicTiling/](lean/OAI/Geometry/PeriodicTiling/), 122 files. A search finds no `sorry`, `admit` or `axiom` in it. I did not build it.
- The statement includes `IsLeast {d | HasTileWithoutPeriodicComplement d} 3`, so minimality of dimension three is part of what is formalized.
- The proof tree contains its own proof of the lower-dimensional cases (`OneDimensional.lean`, `PlanarPeriodicity.lean`, `WeakPeriodicity.lean`, `LowerDimensions.lean`), while the paper cites the literature for them.

File names that appear to mirror the paper, matched by name only:

| Paper | Lean file |
|---|---|
| §2 word rule, Lemma 2.5 model | [WordObstruction.lean](lean/OAI/Geometry/PeriodicTiling/WordObstruction.lean), [LastDigit.lean](lean/OAI/Geometry/PeriodicTiling/LastDigit.lean) |
| §3 cyclic group, dependence, cycles | [CyclicCoordinates.lean](lean/OAI/Geometry/PeriodicTiling/CyclicCoordinates.lean), [Dependence.lean](lean/OAI/Geometry/PeriodicTiling/Dependence.lean), [CycleTests.lean](lean/OAI/Geometry/PeriodicTiling/CycleTests.lean) |
| §4 activation | [ActivationTiles.lean](lean/OAI/Geometry/PeriodicTiling/ActivationTiles.lean), [Activation.lean](lean/OAI/Geometry/PeriodicTiling/Activation.lean), [ShiftMultiplicity.lean](lean/OAI/Geometry/PeriodicTiling/ShiftMultiplicity.lean) |
| §5 common solution | [BlockActivationInverse.lean](lean/OAI/Geometry/PeriodicTiling/BlockActivationInverse.lean), [CommonSolution.lean](lean/OAI/Geometry/PeriodicTiling/CommonSolution.lean) |
| §6 stacking | [CyclicStack.lean](lean/OAI/Geometry/PeriodicTiling/CyclicStack.lean), [FullDifferencePartition.lean](lean/OAI/Geometry/PeriodicTiling/FullDifferencePartition.lean) |
| §7 residues, rigidity, marked tile | [ResidueTransversal.lean](lean/OAI/Geometry/PeriodicTiling/ResidueTransversal.lean), [GridRigidity.lean](lean/OAI/Geometry/PeriodicTiling/GridRigidity.lean), [MarkedTile.lean](lean/OAI/Geometry/PeriodicTiling/MarkedTile.lean), [ConstituentAssembly.lean](lean/OAI/Geometry/PeriodicTiling/ConstituentAssembly.lean) |
| Top-level theorem | [Counterexample.lean](lean/OAI/Geometry/PeriodicTiling/Counterexample.lean) |

## References the paper builds on

Bhattacharya (periodicity of tilings of ℤ²); Greenfeld and Tao (structure of translational tilings in ℤ^d; a counterexample to the periodic tiling conjecture, Annals 2024; undecidability results); Meyerovitch, Sanadhya and Solomon (reduction of tiling problems); de Dios Pont, Grebík, Greenfeld and Madrid (rational polygonal sets); Newman; Lagarias and Wang; Wijshoff and van Leeuwen. See the paper's reference list for full citations.

## Further reading

- **Greenfeld and Tao, [A Counterexample to the Periodic Tiling Conjecture (Announcement)](https://arxiv.org/abs/2209.08451), arXiv:2209.08451, 18 September 2022.** Announces the disproof of the periodic tiling conjecture for sufficiently large dimension d. It also disproves the analogous continuous conjecture for Euclidean space and gives a counterexample in a group ℤ² × G₀ with G₀ a finite abelian group. Its method encodes a class of "p-adically structured functions" through functional equations (the arXiv abstract page, as fetched; I did not read the full text). This is the starting point that entry 155 lowers to dimension three.
- **Terence Tao, [A counterexample to the periodic tiling conjecture](https://terrytao.wordpress.com/2022/11/29/a-counterexample-to-the-periodic-tiling-conjecture-2/), blog post, 29 November 2022.** Announces the full version, arXiv:2211.15847. The post says the original plan was to encode a "p-adic Sudoku puzzle", but the authors could not make the 2-group-valued functions interact with the Z/pZ "clock" functions, so they switched to a **2-adic** Sudoku with a large power of 2 (which they believe to be 2^10). This matters for reading entry 155. That paper's introduction says its central construction keeps the finite factor *cyclic*, and its Sudoku uses an odd prime p > 200 (Section 2 above). I have not compared the two constructions in detail. The post does not state the dimension of the counterexample or mention Wang tiles. Its comments link arXiv:2212.02429 (Khare, on a proof and quantitative refinement) and arXiv:2303.10798 (Smith et al., an aperiodic monotile); I only saw those through the post.
- **IAS Members' Colloquium, [Einstein, P-adic Sudoku, Domino, and Decidability](https://www.ias.edu/video/einstein-p-adic-sudoku-domino-and-decidability)**, 4 December 2023 (video). Speaker: [Rachel Greenfeld](https://www.ias.edu/scholars/rachel-greenfeld) (supplied by you; the abstract pages I could read do not name the speaker). The abstract describes translational tiling, asks whether one can decide if a set is a translational tile and whether every translational tile admits a periodic tiling, says these questions are closely related, and presents new developments "joint with Terence Tao" answering both. The abstract page is [ias.edu/node/274291](https://ias.edu/node/274291). This is the decidability side of the story that the paper's references to Greenfeld–Tao undecidability results point to.

## Related pages

- [TILING.md](TILING.md) for the other tiling, dimer, matching and packing papers.
- [COMPLEXITY.md](COMPLEXITY.md) for the complexity side. Nothing in the repository connects entry 155 to it.

## Caveats

- Section-level statements come from my reading of the paper; I have not verified any proof, and the Lean tree is unbuilt here.
- The small numerical check above is mine and covers a finite window only.
- The comparison with Wang and Penrose tilings is background, not a claim of the paper.
- The paper concerns *full* periodicity. It does not rule out tilings with some nonzero period.

## To do

- [ ] Read Sections 5.2 to 5.3, 7.1 and 7.2 closely.
- [ ] Build the Lean tree with the Comparator and confirm the verdict.
- [ ] Read the Greenfeld–Tao Annals paper to see how much of Section 2 is theirs.
- [ ] Optional: a `tools/` script that prints the last-nonzero-digit array and checks the line rule for small p, with a figure.
