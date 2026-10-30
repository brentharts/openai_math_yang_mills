# Tilings, Matching Rules and Dimers

A fork-level guide to the results in this repository about tilings, packings and perfect matchings, and to how they connect. It is not an official OpenAI document.

**How it was built.** I searched the numbered summaries in [CONTENTS.md](CONTENTS.md) and the preprint titles for tiling, matching, dimer, packing and lattice-model vocabulary. I read the opening pages (abstract, introduction, and in the first case the proof strategy) of the papers marked "read" below. Every other row describes its CONTENTS summary only. Statements marked *(background)* are general knowledge and are not from this repository.

Entry numbers are those of CONTENTS.md. "Lean" links go to the formalization scope notes where one exists.

## The one tiling paper

A dedicated page, [PERIODIC_TILING.md](PERIODIC_TILING.md), walks through this paper's proof and its Lean tree in detail. It also lists outside reading: the Greenfeld–Tao [announcement](https://arxiv.org/abs/2209.08451), [Tao's blog post](https://terrytao.wordpress.com/2022/11/29/a-counterexample-to-the-periodic-tiling-conjecture-2/) and an [IAS colloquium](https://www.ias.edu/video/einstein-p-adic-sudoku-domino-and-decidability) on the same problem.

| Entry | Result | Status |
|---|---|---|
| **155** | [A translational tile with no fully periodic tiling in dimension three](preprints/A-translational-tile-with-no-fully-periodic-tiling-in-dimension-three-September-23-2026/) ([Lean](lean/docs/155.md), [statement](lean/ComparatorChallenges/PeriodicTilingThree.lean)) | read |

**The theorem.** There is a finite set T ⊂ ℤ³ that tiles ℤ³ by translations, but no translation set A with A ⊕ T = ℤ³ is invariant under a finite-index subgroup. The unit-cube thickening T + [0,1]³ tiles ℝ³ up to null sets, and no such tiling has a translation set invariant under a full-rank lattice, even allowing arbitrary real translation vectors. This disproves the periodic tiling conjecture in dimension three.

**Why dimension three is the end of the line (as stated in the paper).**
- ℤ: every tiling by a finite tile is periodic (Newman's finite-state argument).
- ℤ²: every finite tile admits a periodic tiling (Bhattacharya), with a quantitative version by Greenfeld and Tao.
- ℤ³: this paper. Greenfeld and Tao had earlier disproved the conjecture in sufficiently large dimension; the new result lowers it to the smallest possible lattice dimension.

**Where "matching rules" appear.** The proof does not use Wang tiles or Penrose-style edge labels. Its engine is a *local rule that forces a global property*, which is the same idea behind matching rules:
- A "Sudoku" array W(n, m) over a finite alphabet, where a finite list of allowed words constrains the samples along every integer-slope line. The paper's word rule has two features: an array exists with every column nonconstant, but no such array has a vertical period.
- Those word constraints are turned into several tiling equations sharing one unknown translation set, then combined into a single tile in ℤ² × ℤ/Qℤ, and finally passed to ℤ³ (ℤ² times a finite cyclic group is a quotient of ℤ³, and the finite factors are chosen with coprime orders so the product stays cyclic).
- A fully periodic tiling is "decoded" back into a periodic Sudoku array, which the word rule forbids.

The paper says it follows the p-adic Sudoku strategy of Greenfeld and Tao and does not claim connectedness of the tile, or absence of every individual nonzero period.

*(background)* Aperiodic tile sets such as Wang tiles and Penrose tilings force non-periodicity of *every* tiling of a given tile set. Entry 155 is a different statement: the tile has tilings, but none is periodic. The paper does not discuss Wang tiles or Penrose tilings, and I did not find them in any summary.

## Perfect matchings and dimers

A perfect matching covers every vertex by exactly one edge. On a planar lattice such as the square grid, perfect matchings are the same objects as dimer coverings, which are dominoes tiling a region *(background; the double-dimer paper itself states "A dimer covering is a perfect matching")*.

| Entry | Result | Status |
|---|---|---|
| **113** | [FPRAS for counting perfect matchings](preprints/A-Fully-Polynomial-Randomized-Approximation-Scheme-for-Perfect-Matchings-in-General-Graphs-September-23-2026/) ([Lean](lean/docs/113.md)) | summary |
| **113** | [Entropy and face dimension of the perfect-matching polytope](preprints/Entropy-and-Face-Dimension-of-the-Perfect-Matching-Polytope-September-23-2026/) | read |
| **120** | [Almost-linear-time maximum-cardinality matching](preprints/Almost-Linear-Time-Maximum-Cardinality-Matching-in-Sparse-General-Graphs-September-24-2026/) | read (abstract) |
| **126** | [Exponential PSD rank of matching matrices](preprints/Exponential-PSD-rank-of-positively-shifted-matching-matrices-October-5-2026/) ([Lean](lean/docs/126.md)) | summary |
| **226** | [Curve scaling limit of half-plane double dimers](preprints/The-curve-scaling-limit-of-half-plane-double-dimers-September-23-2026/) | read |
| **378** | [The ETH is false](preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/) ([Lean](lean/docs/378.md)) | read earlier; see [ETH.md](ETH.md) |

- **Entropy (113).** The paper proves the Anari–Oveis Gharan–Vinzant conjecture: for edge marginals x in the matching polytope, F(x) − (2 − 2/m)B(x) ≤ H(x) ≤ F(x), where H is the maximum entropy of a matching law with those marginals, and F, B are explicit sums of one-coordinate functions. It holds on the boundary too, and the paper also proves |supp x| − dim F_x ≤ 3m − 2. It states that the result yields lower counts from uniform marginals.
- **Algorithms (113, 120).** Counting perfect matchings approximately is polynomial (113); finding a maximum matching takes (n+m)^(1+o(1)) randomized time (120). Both sit on the tractable side.
- **Exact counting (378).** [ETH.md](ETH.md) records that Corollary 1.3 of the ETH paper gives *exact* counting of perfect matchings in 2^(o(v+e)) time and polynomial space. Approximate counting (113) and exact counting (378) are different guarantees on the same problem.
- **Lower bound (126).** Every exact semidefinite lift of the perfect-matching polytope has exponential size, even for the PSD rank of its odd-cut slack matrix after a fixed positive shift. So matching is easy for algorithms but not for small semidefinite programs. The summary does not say anything about dimers.
- **Scaling limit (226).** Two independent dimer coverings of the Temperleyan square lattice in the upper half-plane are superposed. Almost surely the superposition is doubled edges plus finite simple even cycles. The paper proves the complete loop ensemble converges to nested CLE₄, matching every macroscopic loop as an unparametrized curve. It uses the standard infinite-volume law, equivalently the inverse Kasteleyn operator vanishing at infinity. This is a statement about uniform random tilings at large scale. No Lean doc is listed.

## Height functions and lattice models

These results describe random lattice configurations through a height function, the same device used for random domino tilings *(background)*. The CONTENTS summaries state the height-field limits below; none of them says it is about tilings.

| Entry | Result | Status |
|---|---|---|
| **225** | [GFF limit of the balanced six-vertex model](preprints/The-Gaussian-free-field-limit-of-the-balanced-six-vertex-model-with-variance-multiplier-1-over-arcsin-c-over-2-September-23-2026/) | summary |
| **232** | [Triangular-lattice Lipschitz heights](preprints/Uniform-real-Lipschitz-surfaces-on-the-triangular-lattice-October-6-2026/) | summary |
| **237** | [Honeycomb self-avoiding walk, exponent 3/4](preprints/Critical-strip-crossing-mass-on-the-honeycomb-lattice-September-26-2026/) | summary |
| **223** | Random-cluster interfaces and SLE | summary |

- 225: with a = b = 1 and 0 < c ≤ 2 the six-vertex height field converges to a Gaussian free field, with exact variance multiplier 1/arcsin(c/2), including the endpoint c = 2.
- 232: Gaussian free field limits for triangular-lattice height models with Lipschitz constraints; at a tuned boundary amplitude the interface converges to chordal SLE₄.
- 226 and 225/232 share the same limit objects (CLE₄, Gaussian free field, SLE₄), which is the clearest thread connecting the dimer paper to the height-model papers. I have not read the papers to see whether they cite each other.
- For the physics side, see [PHYSICS.md](PHYSICS.md) and the BKT paper (entry 216) listed there.

## Packings and coverings

Packing is the dual problem to covering: how densely can copies of a shape fit without overlap?

| Entry | Result | Status |
|---|---|---|
| **090** | [Triangular-lattice optimality](preprints/Universal-optimality-of-the-triangular-lattice-September-23-2026/) ([Lean](lean/docs/090.md)) | summary |
| **090** | [Sharp Fourier certificate for planar circle packing](preprints/A-sharp-Fourier-certificate-for-planar-circle-packing-September-23-2026/) | read (abstract) |
| **343** | [Symplectic ball packing in higher dimensions](preprints/Symplectic-Ball-Packings-in-Higher-Dimensions-September-23-2026/) ([Lean](lean/docs/343.md)) | summary |
| **092** | Optimal order Θ(n log n) of convex-body covering density | summary |
| **100** | Cylinder coverings of the tetrahedron below half the projection area | summary |
| **156** | Borsuk's conjecture fails in dimension nine | summary |
| **118** | Bin packing and configuration-LP gaps | summary |

- The circle-packing paper resolves the planar Cohn–Elkies sharpness conjecture: the two-point Fourier bound attains the optimal density π/(2√3), via a radial Schwartz certificate whose signs are proved with rigorous interval arithmetic. The certificate also recovers uniqueness of the triangular packing among periodic equality cases. It ships a [verification folder](preprints/A-sharp-Fourier-certificate-for-planar-circle-packing-September-23-2026/verification/) that you can run.
- 090 shows the triangular lattice minimizes energy for a large class of potentials (completely monotone in squared distance), which is the physics-flavoured counterpart of the packing result.
- 343 gives an exact criterion for packing symplectic balls: they embed disjointly when Σ R_i^n < R^n and R_i + R_j < R for all pairs.
- A tiling is a packing that leaves no gaps; a covering allows overlaps. The summaries do not relate 155 to 090 or 092, and I do not claim a technical link.

## Matchings, Hamiltonian cycles and Tutte states

| Entry | Result | Status |
|---|---|---|
| **180** | [Barnette's conjecture](preprints/Paired-states-and-Hamiltonian-cycles-in-cubic-bipartite-planar-graphs-September-24-2026/) ([Lean](lean/docs/180.md), [statement](lean/ComparatorChallenges/BarnetteHamiltonian.lean)) | read |
| **111, 114, 185** | Matroid results (prophet inequalities, counting common bases, infinite intersection) | summary |

Barnette's conjecture is proved: every finite simple cubic bipartite planar 3-vertex-connected graph has a Hamiltonian cycle. Matchings enter the proof directly. In the dual triangulation, the paper assigns each remaining black triangle one of its vertices so that every non-root vertex is used exactly once. These face–vertex incidence matchings are Tutte states. It then looks for a *pair* of states making distinct choices at every triangle, whose edge set is acyclic. A pair of matchings whose union is a family of cycles is also what the double-dimer paper (226) studies in the plane. The Barnette paper says Hamiltonicity corresponds to splitting the triangulation's vertices into two induced trees.

Corollary 7.2 of the paper also covers cubic 3-vertex-connected bipartite *Pfaffian* graphs: every three-edge path lies in a Hamiltonian cycle. *(background)* Pfaffian graphs are the graphs for which Kasteleyn's method counts perfect matchings, which is the setting of the inverse Kasteleyn operator in entry 226. The paper itself cites a published equivalence rather than the dimer connection.

## Connections at a glance

| Thread | Entries | What is stated |
|---|---|---|
| Local rules force a global property | 155 | Sudoku word rule: nonconstant columns forbid a vertical period. |
| Matchings as tilings | 226, 113 | The dimer paper: a dimer covering is a perfect matching. |
| Pairs of matchings give loops | 226, 180 | Double dimers (superposition of two matchings) and paired Tutte states. |
| Approximate vs exact counting | 113, 378 | FPRAS vs exact 2^(o(v+e)) algorithm for the same count. |
| Matchings easy, lifts hard | 113, 120, 126 | Polynomial counting and near-linear search, but exponential SDP lifts. |
| Common scaling limits | 223, 225, 226, 232 | GFF, CLE₄ and SLE limits for lattice height and loop models. |
| Optimal packings | 090, 343 | Triangular lattice and ball packing criteria. |

## What the repository does not contain

A search of summaries and titles found no papers on Wang tiles, Penrose or other aperiodic tile sets with matching rules, quasicrystals, domino-tiling enumeration (such as the arctic circle) or polyomino tilings. Several of the papers above were not opened, so a relevant paragraph inside one of them is possible.

## Related pages

- [PERIODIC_TILING.md](PERIODIC_TILING.md) for entry 155 in depth.
- [COMPLEXITY.md](COMPLEXITY.md) and [ETH.md](ETH.md) for the algorithmic side of perfect matchings.
- [PHYSICS.md](PHYSICS.md) for the lattice-model and height-function papers.

## Caveats

- Rows marked "summary" are described from CONTENTS only; "read" means I read the abstract and introduction, not the proofs.
- I have not checked any proof, and the Lean notes say only which statements are formalized.
- Connections marked *(background)* come from general knowledge. The papers may or may not mention them.
- Entry 155's cyclic-group remark is my paraphrase of the paper's sentence that Z² times a finite cyclic group is a quotient of Z³.

## To do

- [ ] Read the double-dimer paper (226) for citations to the six-vertex and height-function papers (225, 232).
- [ ] Read entry 155's Sections 3 to 5 (the dependence tests) and Section 7 (the passage to ℤ³).
- [ ] Check whether the Lean statement for 155 (`PeriodicTilingThree.lean`) matches Theorem 1.1(i).
