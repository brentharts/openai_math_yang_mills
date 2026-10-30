# Scale Hierarchies: Tilings, Renormalization and Yang–Mills

A second scan of the preprints for links to the Spectre and Chair44 aperiodic monotiles (see [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md)) and to Yang–Mills physics. It is a fork-level note, not an official OpenAI document. The thread that survived the scan is **scale hierarchies**: structures that map a solution at one scale to a solution at another. Entry 155 uses one in a proof, Chair44 is built from one, and the 2D O(N) papers use one to construct a field theory.

**How it was done.**
- I extracted the text of all 751 preprint PDFs (every one gave text) and searched it for about 60 terms: tilings and aperiodic order, hierarchy and p-adic vocabulary, and gauge theory. I read the matching passages in context, and the opening pages of the most promising papers.
- My first pass had a regular-expression bug (an en-dash inside a bracket), which silently missed "Yang–Mills" and "Bałaban". I caught it because I had already seen both in paper 365, and reran those terms with a corrected pattern.
- The scan finds words, not mathematics. Content that appears only in formulas or figures is invisible to it.

Entry numbers are those of [CONTENTS.md](CONTENTS.md). Section numbers inside papers come from their tables of contents or text. Anything about Chair44 or the Spectre is cited by the keys in [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md#citation-key): [T] Tsiokos, [HR] Hartshorn–Romanov, [F] Flicker, [GS] Goodman-Strauss.

## First, what the scan did not find

| Term searched | Result |
|---|---|
| Chair44, Spectre (as a tile), supertile, quadtree, quasicrystal, diffraction, Wang tile, Delone set, Meyer set, Ammann | **No hits.** The three "spectre" hits are the French word *spectre* in reference titles. |
| Penrose | 31 papers, but all are other Penroses: the spacetime Penrose inequality (entry 260), Mathew Penrose's *Lectures on the Poisson Process*, the Lebowitz–Penrose limit, Penrose's tree-partition identity. **None is about Penrose tilings.** |
| "tiling" | 43 papers, mostly in the sense of plane subdivisions: the honeycomb and triangular tilings behind the honeycomb-walk and lattice papers, and Voronoi cells. [Entry 155](PERIODIC_TILING.md) is the only paper on translational tilings. The perfect-matching paper (113) uses "tilings" of gadgets by dimers. |
| Wilson loop, quark, gluon, lattice QCD | **No hits** (the 8 case-insensitive "QCD" matches were not real instances). |
| "lattice gauge" | 3 hits, all in one paper (the O(3) continuum limit, entry 215), as a description of Bałaban's earlier work. |

So the preprints contain no aperiodic-order content and never mention Chair44 or the Spectre (Tsiokos's revision 2 is dated 16 September, a week before most of the batch), and none treats a gauge theory. Any connection below is a connection of method or vocabulary, not a citation.

## Spectre and Chair44: connections found

| # | Connection | Strength | Where |
|---|---|---|---|
| 1 | Block renormalization in the O(N) papers | Structural analogy | below |
| 2 | p-adic exclusion on manifolds (304) | Vocabulary, with a caution | below |
| 3 | Pisot numbers and self-similar sets (148, 153) | Checked; no direct application | below |
| 4 | Reconstruction on trees (229) | Thematic | below |
| 5 | Multiscale and ultrametric arguments (261, 221, perceptron papers) | Thematic | below |
| 6 | Triangular-lattice results (090, 232, 237) | Weak; background fact only | below |
| 7 | Flicker's 24-state lattice model | Actionable idea | below |

### 1. Block renormalization and substitution hierarchies

This is the strongest link. The sources say nothing about each other, so the comparison is my own.

| | Chair44 / Spectre | 2D O(3) and O(4) papers (entry 215) | Entry 155 |
|---|---|---|---|
| Scale step | Factor 2, matrix 2I₃, determinant 8 ([T] §3.1; [HR] §6); Spectre linear factor λ = √(4+√15) ≈ 2.806 ([HR] §6, Fact 16) | "Positive transfer and repeated doubling" ([O3 paper](preprints/The-canonical-massive-continuum-limit-of-the-two-dimensional-O3-model-October-4-2026/) §8.2); rectangles are doubled again and again | Rescaling m ↦ pm ([E155] Prop. 2.4) |
| Exactness | The substitution is exact: eight chairs partition the doubled chair ([T] §3.1, Table 1) | "An exact renormalization step" (O3 §5); "Exact renormalization of the spin sources" (§9) | The rescaled array obeys the same line rule exactly |
| Self-improvement | Each level is again a tiling by Q ([T] Theorem 6.5) | Proposition 8.3: if the doubling defect Δ ≤ 2⁻¹⁰ then Δ(2ᵏn, 2ᵏw) ≤ 64⁻¹·2^(−2ᵏ), and λ₂/λ₁ ≤ e^(−(log 2)/n) at circumferences 2ᵏn | A period M becomes M/p, so a minimal period gives a contradiction |
| Extensivity | Volume grows by 8 per level | Proposition 2.1(ii), eq. (2.3): 0 ≤ 4 log Z(n,w) − log Z(2n,2w) ≤ C(1+b+n+w)ᵖ e^(−c·min(n,w)/X(b)); area grows by 4 per doubling | — |
| The hard part | Supertiles need not cover space; there is no exhaustion clause and "fault lines" run along supertile boundaries ([T] Theorem 6.6, §6) | The nearest-neighbour model "generates interactions between blocks whose spatial decay must be controlled" (O3 §1.1) | Stacking and decoding ([E155] §§4–6) |
| Solvable cousin | Unique parent cluster at every level ([T] Theorem 6.4) | Gawędzki–Kupiainen's sigma model with a *hierarchical kinetic term*, where "the hierarchical covariance reduces the flow to a single-spin recursion" (O3 §1.1) | — |
| Outcome | Per(T) = {0}: no translation period ([T] Theorem 1.2) | Positive mass gap, Osterwalder–Schrader reconstruction (O3 §11) | No fully periodic tiling |

**Reading it.** In both tiling hierarchies and block renormalization, a coarse object is built exactly from finer ones and has the same form as the original. In both, the idealized version (unique parent, hierarchical covariance) is clean, and the difficulty lies in what crosses block boundaries. The O(3) paper places itself in the tradition of Kadanoff, Wilson and Bałaban (§1.1). What differs: block renormalization seeks a *limit* of a flow of couplings and controls errors, whereas a substitution is a fixed combinatorial rule and the question is what it forces. Nothing here suggests a theorem transfers.

### 2. p-adic exclusion on manifolds (entry 304)

- [The Hilbert–Smith paper](preprints/The-Hilbert-Smith-conjecture-in-every-finite-dimension-September-23-2026/) proves (Theorem 1.2, "p-adic exclusion") that every jointly continuous action of ℤ_p = lim ℤ/pᵏℤ on a connected manifold has nonzero kernel and factors through a finite quotient. Theorem 1.1 then gives that any locally compact group acting faithfully on such a manifold is a Lie group.
- [HR] §7.4, Corollary 2 puts the internal space of the Chair44 hull at ℤ₂³.
- **Caution.** These do not interact. A tiling hull is not a manifold: it is locally a Euclidean factor times a Cantor-like one. A 2-adic odometer living on the hull is consistent with Theorem 1.2, which restricts actions on manifolds only. I list it because 2-adic groups appear in both, not because one constrains the other.

### 3. Pisot numbers and self-similar sets (entries 148, 153)

I checked whether either paper applies to the Spectre window.
- [Entry 153](preprints/Arithmetic-classification-and-non-Pisot-singularity-for-Bernoulli-convolutions-October-3-2026/) (§1) says: for λ < 1/2 a direct interval cover proves singularity; at λ = 1/2 the Bernoulli convolution is uniform on [−2, 2]; the interesting range is (1/2, 1); a Pisot number is a real algebraic integer > 1 whose other conjugates have modulus < 1; Erdős proved singularity when λ⁻¹ is Pisot.
- Chair44's factor 2 is the λ = 1/2 case. The Spectre window is the attractor of contracting similarities of ratio 1/λ ([HR] Theorem 3) with 1/λ ≈ 0.356 < 1/2, the trivially singular range. So the range 153 is about does not apply to either tile.
- By my computation from the polynomial x⁴ − 8x² + 1 in [HR] Fact 16: λ² = 4 + √15 is a Pisot unit (conjugate 4 − √15 ≈ 0.127), but λ itself is neither Pisot nor Salem, because its conjugate −λ has modulus > 1.
- [Entry 148](preprints/The-entropy-rate-dimension-formula-for-self-similar-measures-on-the-line-September-24-2026/) is for self-similar measures *on the line*; the Spectre window is planar.

*(background)* In substitution tiling theory, Pisot inflation factors are tied to pure-point spectrum. The repository says nothing on this.

### 4. Reconstruction on trees (entry 229)

- The three-state paper proves that the symmetric broadcast process on a b-ary tree reconstructs exactly when dλ² > 1, with non-reconstruction at equality: information about the root survives to depth n or it does not.
- [HR] §10 runs "clock" tests on the Spectre hierarchy: the radius needed to name a tile's level-1 ancestor is 3ℓ and grows about threefold per level; the look-back slope is 0.92 for the Spectre and 0.89 for a quadtree.
- **Analogy.** Both ask how much finite local data reveals about position in a tree. **Difference.** A substitution hierarchy is deterministic; reconstruction is about a noisy channel. A noisy-patch version of the clock tests would connect them, and I know of no source that does this.

### 5. Multiscale and ultrametric arguments

- The Anderson-model paper (entry 261, [pure-point spectrum in 2D](preprints/Pure-Point-Spectrum-for-the-Two-Dimensional-Anderson-Model-at-Every-Positive-Disorder-September-23-2026/)) compares quantities "along a hierarchy of increasingly sparse ports" (introduction; §7, "A hierarchy of probability splits and energy comparisons"), in the multiscale tradition it cites from Fröhlich and Spencer.
- The perceptron free-energy papers (22 ultrametric hits in the Ising version) use "Panchenko's ultrametricity mechanism" and Ruelle probability cascades. Entry 221 states the Mézard–Parisi formula as an infimum over "finite-depth hierarchical trial laws".
- A substitution hierarchy also defines an ultrametric on tiles (distance = level of the smallest common supertile). This is my observation. The spectral notions differ: Anderson's pure-point spectrum is for a random operator, and [HR] Fact 22 is pure-point *diffraction*. The repository has no paper on Schrödinger operators on aperiodic tilings.

### 6. Triangular-lattice results

[HR] §6 puts Spectre vertices in ℤ[ζ₁₂] (Fact 15). *(background)* The triangular lattice sits inside it. The repository's triangular-lattice papers (090 optimality, 232 Lipschitz heights, the honeycomb-walk papers) do not mention tilings by monotiles, and I found no stated link.

### 7. An actionable idea: Flicker's 24-state model

[F]'s Discussion proposes a 24-state lattice model on a cubic lattice, with nearest-neighbour penalties for illegal contacts, whose ground states should match Chair tilings. Twenty-four is the order of the proper cubic rotation group ([T] §3.1), a finite subgroup of SO(3). The repository already studies lattice spin models with the same ingredients (O(3) with continuous SO(3)-symmetry, entry 215; four-state Ashkin–Teller and Potts, 233 and 229; XY, 216), and the fork's [tools/](tools/README.md) has a Monte Carlo for O(N). A first model built from the 44-contact atlas of [T] is in [CHAIR44_MODEL.md](CHAIR44_MODEL.md); it uses the atlas as the legality rule and does not reproduce [F]'s own penalty table.

## Yang–Mills and gauge theory

The full-text scan finds Yang–Mills in 14 papers (25 hits), and none treats a gauge theory.

| Group | Entry | What the text says |
|---|---|---|
| 2D O(N) models | 215 (five papers) | All five cite Polyakov's 1975 paper "Interaction of Goldstone particles in two dimensions. Applications to ferromagnets and massive Yang–Mills fields". The O(3) paper also cites Bałaban's lattice gauge theory work (refs [3], [4]) and N-vector work ([5], [6]) as methodological ancestors, and also cites Dybalski–Stottmeister–Tanimoto and Goswami on Bałaban's variational problem. |
| BFSS matrix model | 270 (two papers) | "supersymmetric matrix quantum mechanics obtained by reducing ten-dimensional supersymmetric Yang–Mills theory to one time dimension" ([positive eigenvalues, SU(2)](preprints/Positive-eigenvalues-of-the-relative-SU-2-BFSS-Hamiltonian-October-5-2026/), §1). That paper cites Lüscher's estimates for bosonic Yang–Mills Hamiltonians ([7], Appendix D of Lüscher 1983, on the mass spectrum on a torus). The SU(N) paper cites Graf–Hoppe and Hasler–Hoppe on reduced super-Yang–Mills. |
| Inverse problems | 365, 372 | The [metric and connection paper](preprints/Determination-of-a-metric-and-a-unitary-connection-from-one-boundary-patch-October-5-2026/) recovers a smooth unitary connection on a trivial rank-two bundle from one boundary patch, up to a unitary gauge equal to the identity on the patch. It cites Cekić's Calderón problem for Yang–Mills connections and states that it "imposes no Yang–Mills equation". The elasticity paper (372) cites Eskin's inverse scattering with external Yang–Mills potentials. |
| Geometry and topology | 050, 055, 060, 306 | Citations of the classical gauge-theoretic geometry: Demailly's Hermitian–Yang–Mills approach in the Griffiths-positivity counterexample on the quadric surface (050); Donaldson's anti-self-dual connections in large-volume stability on threefolds (055); Donaldson on orientation of Yang–Mills moduli spaces in the class VII surfaces paper (060); Donaldson's Floer homology in Yang–Mills theory and Uhlenbeck's removable singularities in the cosmetic-surgery paper (306). |
| General relativity | 260 | One citation of Einstein–Yang–Mills black holes (Sudarsky–Wald) in the Penrose inequality with charge. |

**Related terms.** "Higgs" appears in 12 papers, which are Higgs bundles and Hecke eigensheaves (P = W, the Kodaira inequality, geometric Langlands, Hikita); the Hikita paper (044) mentions gauge theory 8 times in the Coulomb-branch sense. "Instanton" appears in 6 papers, 76 times in the cosmetic-surgery paper alone.

**Three observations.**
1. **The mass-gap analogy is cited, not used.** The O(N) papers cite Polyakov for the analogy between 2D sigma models and Yang–Mills. [YANGMILLS.md](YANGMILLS.md) already records that the repository does not bridge 2D to 4D.
2. **The hierarchical trick.** Gawędzki–Kupiainen's hierarchical sigma model is the exactly solvable cousin of nearest-neighbour block RG (O(3) §1.1). The tiling parallel in section 1 above is that solvable-versus-real contrast, but for tilings.
3. **Finite gauge groups.** *(background)* Lattice gauge theory with a finite subgroup of SU(2) is a standard approximation to the continuous theory. The 24-element rotation group in Chair44 is such a subgroup of SO(3), so Flicker's 24-state model (section 7) is a lattice model with a discrete rotation group as its variable. The fork's [yang-mills-mass-gap-lab](https://github.com/brentharts/yang-mills-mass-gap-lab) works with SU(2). This is a suggestion, not a result from any source.

## Caveats

- This is a word-level scan. A connection stated only in formulas or in a paper's proofs would be missed.
- Section numbers for [O3 paper] come from its table of contents; I read §1.1, §8.2 and Propositions 2.1 and 8.3 but not the whole 100-page paper.
- All Chair44 and Spectre statements come from the fetch-tool summaries described in [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md#caveats).
- The "strength" ratings and every parallel are mine.

## To do

- [x] Extract the 44-contact atlas and build a 24-state model: done as a pose-level model, see [CHAIR44_MODEL.md](CHAIR44_MODEL.md). Still open: re-derive the atlas from the pyramid geometry, and test finite-temperature behaviour.
- [ ] Read O3 §5 (the exact renormalization step) next to [T] §6 to see whether the "defect shrinks under doubling" structure of Proposition 8.3 has a tiling analogue.
- [ ] Read the Anderson-model hierarchy (§7) and the Parisi cascade construction as examples of tree-structured arguments.
- [ ] Optional `tools/` script: print the last-nonzero-digit array of entry 155 and test the line rule for small p.

## Related pages

- [SPECTRE_CLOCK_AND_CHAIR44.md](SPECTRE_CLOCK_AND_CHAIR44.md), a closer look at connections 1, 2, 4 and 5 with computations
- [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md), [MONOTILES.md](MONOTILES.md), [PERIODIC_TILING.md](PERIODIC_TILING.md), [TILING.md](TILING.md)
- [YANGMILLS.md](YANGMILLS.md), [PHYSICS.md](PHYSICS.md)
