# Complexity Theory Results

A fork-level guide to the computational-complexity and algorithms results in this repository, organized around the ETH refutation in [ETH.md](ETH.md). It is not an official OpenAI document.

**How it was built.** I searched the numbered summaries in [CONTENTS.md](CONTENTS.md) for hardness, counting, algorithm and lower-bound results, read the full summaries of the ones listed here, and read the ETH paper's introduction and technique overview. Apart from the ETH paper, none of these papers was opened, so each row describes its CONTENTS summary only. The "Relation to ETH" column says what is stated in the sources I read, and where I add general background it is marked.

Entry numbers are those of CONTENTS.md. "Lean" links go to the formalization scope notes where one exists.

## What connects directly to ETH

| Entry | Result | Relation to ETH |
|---|---|---|
| **378** | [The ETH is false](preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/) ([Lean](lean/docs/378.md)) | The result itself. See [ETH.md](ETH.md). |
| **102** | [The Unique Games Theorem](preprints/The-Unique-Games-Theorem-September-23-2026/) ([Lean](lean/docs/102.md)) | **Cited by the ETH paper.** Its reference `[Ope26]` is this preprint, used for the NP-hardness bounds on Max-k-SAT in the context of Corollary 1.2. |
| **113** | [FPRAS for perfect matchings](preprints/A-Fully-Polynomial-Randomized-Approximation-Scheme-for-Perfect-Matchings-in-General-Graphs-September-23-2026/) ([Lean](lean/docs/113.md)) | **Same problem, different guarantee.** Entry 113 gives a polynomial-time randomized *approximation* of the number of perfect matchings. The ETH paper's Corollary 1.3 gives *exact* counting in 2^(o(v+e)) time and polynomial space. |

## Hardness and approximation

| Entry | Result | Notes |
|---|---|---|
| **102** | Unique Games Conjecture; optimal approximation thresholds | Proves Khot's conjecture. Independent direct reductions give NP-hardness, on unweighted graphs, of approximating Max-Cut beyond the Goemans–Williamson ratio, Vertex Cover below factor two, and Min-UnCut and directed feedback vertex set within any fixed constant factor. The summary says these direct proofs use established PCP and Label Cover results. See also [Max-Cut hardness](preprints/A-Direct-Proof-of-Optimal-Max-Cut-Hardness-September-23-2026/). |
| **105** | [Perfect completeness for 2-to-1 games](preprints/Perfect-completeness-for-2-to-1-games-September-23-2026/) ([Lean](lean/docs/105.md)) | Khot's 2-to-1 Games Conjecture with perfect completeness: NP-hard to separate satisfiable games from games of value at most δ, for every fixed rational δ in (0,1). |
| **106** | [Hardness of coloring three-colorable graphs](preprints/Hardness-of-finding-large-independent-sets-in-three-colorable-graphs-September-24-2026/) ([Lean](lean/docs/106.md)) | NP-hard to color a three-colorable graph with any fixed number of colors. A deterministic polynomial-time reduction from 3SAT produces graphs that are three-colorable, or else have no independent set of size δn. |
| **117** | [Uniform sparsest cut](preprints/Constant-factor-hardness-of-uniform-sparsest-cut-September-24-2026/) ([Lean](lean/docs/117.md)) | Constant-factor approximation is NP-hard, and the Goemans–Linial relaxation has integrality gaps of order at least √(log n)/(log log n)^3. |
| **136** | [Quasilinear PCP for PPAD](preprints/The-PCP-for-PPAD-conjecture-a-quasilinear-reduction-September-25-2026/) | End-of-Line instances reduce to numerical circuit constraints of near-linear size with robust local verification. |

Several of these are stated as reductions from 3SAT, and the Lean note for entry 102 says the formalized reductions start from binary 3SAT. No connection to ETH is stated in the summaries or the ETH introduction. As general background (not from the repository), NP-hardness is a statement about polynomial-time reductions, whereas ETH is a claim about exponential running times, so results of this kind do not depend on ETH.

## Counting and algorithms

| Entry | Result | Notes |
|---|---|---|
| **378** | Subexponential exact #k-SAT | See [ETH.md](ETH.md). |
| **113** | FPRAS for perfect matchings; perfect-matching entropy | Also proves the Anari–Oveis Gharan–Vinzant entropy conjecture. |
| **114** | Approximate counting of common integer polymatroid bases | FPRAS with oracle access to rank values. |
| **115** | Sampling and counting contingency tables with arbitrary margins | Exact uniform sampling in expected polynomial bit time. |
| **138** | [Subset Sum in O(2^(0.49n))](preprints/Subset-Sum-in-Time-2-power-0-49n-October-4-2026/) | Randomized, with n counting the integers. No Lean doc is listed. |

On entry 138: the ETH paper mentions SETH-based bounds for *pseudopolynomial* algorithms for Subset Sum. Entry 138 concerns running time exponential in the number of integers, which appears to be a different regime. I have not confirmed this against either paper.

## Lower bounds (unconditional)

The ETH paper's Corollary 1.4 derives a circuit lower bound for E^NP from its algorithm, the algorithms-to-lower-bounds route. The entries below use other methods.

| Entry | Result | Notes |
|---|---|---|
| **112** | [Beyond the square-root exponent for depth-three circuits](preprints/Beyond-the-Square-Root-Exponent-for-Depth-Three-Boolean-Circuits-September-23-2026/) ([Lean](lean/docs/112.md)) | A polynomial-time language needing 2^(ω(√n)) gates in unbounded-fan-in OR–AND–OR circuits. |
| **135** | Homogeneous depth-five lower bounds for iterated matrix multiplication | n^(Θ(√n)) gates over characteristic-zero fields. |
| **108** | Cubic permanent–determinant lower bound | Ω(n^3) for the border determinantal complexity of the permanent. |
| **133** | [Weisfeiler–Leman complexity](preprints/Unconditional-time-lower-bounds-for-Weisfeiler-Leman-equivalence-September-25-2026/) ([Lean](lean/docs/133.md)) | Unconditional n^(Ω(k)) time lower bounds in specified models, and EXPTIME-completeness with dimension as input. |
| **129** | Exponential state costs for two-way automata | Resolves the Sakoda–Sipser state-succinctness conjecture over growing finite alphabets. |
| **243** | Separating choiceless counting from polynomial time | Confirms the Blass–Gurevich–Shelah noncapture conjecture. |

## Derandomization and other

| Entry | Result | Notes |
|---|---|---|
| **103** | [L = RL = BPL](preprints/Exact-Derandomization-of-Logarithmic-Space-L-equals-RL-equals-BPL-September-23-2026/) ([Lean](lean/docs/103.md)) | An effective compiler from randomized logspace machines to deterministic logspace deciders. The ETH paper's algorithm is also deterministic and uses polynomial space, but the two results are about different resources and no link is stated. |

## Random SAT

| Entry | Result | Notes |
|---|---|---|
| **235** | [Limiting random SAT thresholds, variance and computability](preprints/Computing-the-Random-3-SAT-Threshold-September-27-2026/) ([Lean](lean/docs/235.md)) | Finite positive limiting thresholds for every fixed k ≥ 3 and hitting-time variance Θ_k(n), plus computability of the 3-SAT threshold. The summary credits Gaia Carenini with priority for the threshold-existence result (ECCC TR26-229). This concerns random formulas, and ETH is a worst-case statement, so no link is implied. |

## What the ETH result does and does not touch here

- No other numbered summary in `CONTENTS.md` states a result conditional on ETH or SETH, so nothing in this repository is stated as depending on a hypothesis the new paper refutes.
- The hardness entries above are stated as NP-hardness results, which per the general background above are not conditional on ETH. I did not read those papers, so a conditional lemma inside one of them cannot be excluded.
- The new paper's own limits (no general subexponential algorithm for NP, other hypotheses such as the Orthogonal Vectors Conjecture still open) are listed in [ETH.md](ETH.md).

## Related pages

- [ETH.md](ETH.md) for the ETH refutation.
- [TILING.md](TILING.md) for the perfect-matching results (entries 113, 120, 126) in the context of dimers and tilings.

## Caveats

- Rows other than 378 and 102 are described from their CONTENTS summaries only.
- "Relation to ETH" states only what the sources I read say. Where nothing is stated, the column says so rather than guessing.
- Folder links go to a preprint's folder (README and PDF). Where an entry has several papers, I linked the one whose title matches the summary best.

## To do

- [ ] Read the Lean scope note for entry 102 in full and check which of the five Unique Games papers the ETH paper's section 8.1 citation points to.
- [ ] Check the ETH paper's reference list for any other repository results it cites.
- [ ] Read entries 113 and 138 against the ETH paper's Corollary 1.3 and its Subset Sum discussion.
