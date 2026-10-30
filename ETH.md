# The Exponential Time Hypothesis Result

A fork-level guide to *The Exponential Time Hypothesis is False* (October 25, 2026) and the news around it. It is not an official OpenAI document.

**What was read.** The paper's abstract, introduction, notation, the full technique overview (section 2), theorem and lemma statements throughout, and Appendix C.1 (circuit lower bounds). The detailed proofs in sections 3 to 6 and Appendices A, B and D were not read, so those are described only through the overview's own account of them. News claims come from secondary sources and are labeled.

## The problem in brief

The **Exponential Time Hypothesis (ETH)** says 3-SAT on n variables needs 2^(Ω(n)) time in the worst case. Variants strengthen it: **SETH** (k-SAT needs 2^((1−o(1))n) as k grows), **Gap-ETH**, **NSETH**, and the counting version **#ETH**. Many conditional hardness results in parameterized and fine-grained complexity are built on them. For the other complexity results in this repository, see [COMPLEXITY.md](COMPLEXITY.md).

## What the repository contains

| Item | Link |
|---|---|
| Paper (67 pages) | [paper.pdf](preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/paper.pdf) |
| Preprint folder and citation | [README](preprints/The-Exponential-Time-Hypothesis-is-False-October-25-2026/README.md) |
| Lean formalization scope | [lean/docs/378.md](lean/docs/378.md) |
| Comparator statement | [SubexponentialSAT.lean](lean/ComparatorChallenges/SubexponentialSAT.lean) |
| Index entry | **378** in [CONTENTS.md](CONTENTS.md) |
| Release note | October 26, 2026 in [history.md](history.md) |

A search of `CONTENTS.md`, `history.md`, `README.md` and the Lean docs finds no other mention of ETH or SETH.

## The main theorem

From the paper (Theorem 1.1): for each fixed k ≥ 2 there is a deterministic algorithm that counts the satisfying assignments of a k-CNF with n ≥ 20 variables and m clauses in time (n+m)^(O_k(1)) · 2^(O_k(n log log n / √(log n))), using (n+m)^(O_k(1)) bits of space. When m ≤ Dn for a fixed D, the time improves to (n+m)^(O_{k,D}(1)) · 2^(O_{k,D}(n/√(log n))). Counting also decides satisfiability, so the 3-SAT case is subexponential.

**The Lean statement is weaker than the manuscript.** According to the [Lean scope note](lean/docs/378.md), the formalized bound is 2^(O_k(n (log log n)^2 / √(log n))), with at most polynomially many tape cells visited, for a canonically encoded k-CNF. The sharper manuscript bounds above are not what Lean proves. The paper's rational-interpolation argument gives the exponent n (log log n)^2 / √(log n), which it calls ρ(n), and Appendix D reaches the sharper bounds by tensor contraction.

CONTENTS entry 378 writes the exponent as n/(log n)^(1/2+o(1)), which is the same quantity as ρ(n) up to how the o(1) is packaged.

## Consequences stated in the paper

- **Corollary 1.2:** for every fixed k and ε > 0, a randomized algorithm finds an assignment satisfying at least OPT − εm clauses of a k-CNF in 2^(O(ρ(n))) time. The paper's discussion cites matching NP-hardness bounds for Max-k-SAT from the Unique Games Theorem, and the citation `[Ope26]` is OpenAI's own [The Unique Games Theorem](preprints/The-Unique-Games-Theorem-September-23-2026/) (September 23, 2026), entry **102** in this repository ([Lean](lean/docs/102.md)). So this corollary's context depends on another result in the release.
- **Corollary 1.3:** counting independent sets, and for each fixed q ≥ 2 proper q-colorings, in 2^(o(v)) time and polynomial space on graphs with v vertices, and counting perfect matchings exactly in 2^(o(v+e)).
- **Corollary 1.4:** a language in E^NP with circuit size greater than c · n · √(log n) for infinitely many n, a superlinear lower bound for circuits with arbitrary fan-in-two gates. Corollary C.1 gives the qualitative version E^NP ⊄ SIZE(O(n)).
- Appendix C also covers inverting Boolean maps with constant output locality, including constant-locality Goldreich-style candidate one-way functions. The paper notes this does not rule out one-wayness against polynomial-time algorithms.

## What it refutes and what it does not

The paper says Theorem 1.1 refutes ETH, SETH, Gap-ETH, NSETH and #ETH, and the bounded-error randomized versions of ETH, SETH and Gap-ETH, since the algorithm is deterministic. It notes that Ryan Williams had refuted the Merlin–Arthur version of SETH in 2016.

It also states several limits:

- **No general speedup for NP.** A 2^(o(n)) running time does not transfer under general polynomial-time reductions (that would need 2^(n^(o(1)))), so the paper says the result does not imply subexponential algorithms for general NP problems.
- **Other hypotheses survive.** Some lower bounds that follow from SETH also follow from weaker hypotheses that may still hold, such as the Orthogonal Vectors Conjecture, and the hypothesis that unrestricted CNF-SAT has no 2^((1−δ)n) algorithm. The theorem is for each fixed clause width k, while that hypothesis lets the width grow with n.
- **Proof-system lower bounds are unaffected.** The algorithm does not conflict with unconditional exponential lower bounds for restricted proof systems such as regular resolution.
- **Limited transfer to other problems.** The paper says a faster algorithm for a target problem needs a reduction from that problem to fixed-width SAT that controls the number of variables, and gives Set Splitting with bounded set size as an example where one exists.

## How the algorithm works

This follows section 2 of the paper. Symbols are as in the paper.

**1. Reduce to a sparse core case.** The core algorithm handles a 3-CNF in which every variable occurs in at most three clauses, so the total clause size M is at most 3n.

**2. A polynomial for the count.** Each clause has a satisfaction polynomial Φ_e = 1 − ∏(1 − literal), equal to 1 on Boolean inputs exactly when the clause holds. Then #F = Σ over x in {0,1}^n of ∏ Φ_e(x), a polynomial of total degree at most M. Evaluating it term by term is exhaustive search.

**3. The easy case.** If the variables split into r groups with no clause crossing groups, the count factors as a product of per-group counts, costing about r · 2^(2n/r) work instead of 2^n. General formulas do not split this way, which is the problem the rest of the construction solves.

**4. Make groups independent by deformation.** Introduce a parameter T and replace every clause input by a linear form whose coefficients are polynomials in T of degree at most R. At T = 0 the forms equal the original variables, so H(0) = #F, where H(T) is the deformed sum. At selected points t = 1, …, Q the forms are built so that each clause depends on variables outside its assigned group only through a few shared integer **linear measurements** of the other groups (an (r, τ, D)-design). Fix the measurements and the sum factors group by group. The cost is the number of measurement tuples, (2D+1)^(2τn), times the group size 2^(2n/r).

**5. Recover the answer at T = 0.** Since deg H ≤ 3nR, the values at Q points are not enough, so the algorithm also computes several Taylor coefficients at each point (about 1 + nR/Q of them) and recovers H(0) by Hermite interpolation. The Taylor coefficients come from sums with some inputs fixed to the all-ones pattern on a subset U, combined by inclusion–exclusion over subsets of at most q = ⌊3nR/Q⌋ variables, giving a third cost of about Σ_j C(n,j)·2^j.

**6. Parameters.** Reading the paper's (2.10): R is about √(log n), Q about log n / log log n, and r about √(log n) / (log log n)^2. These make each of the three costs 2^(O(n (log log n)^2 / √(log n))), which gives the stated bound for this case.

**7. The deformation itself (section 2.7).** Each replacement row is built from a matrix M_0 = V_0 f(W_0), made of commuting swap permutation matrices on coordinates indexed by pairs of group choices across the Q points. The polynomial f interpolates (−1)^j on an interval of R consecutive integers near the mean of a binomial distribution. This makes M_0 fix most basis directions, so only a small fraction of directions (at most Dτ/(3r) rows) need to be measured. Multipliers of the form ∏(1 − T/s) switch off, at each point s, any term that would read a variable outside the assigned group, at a cost of one degree each. Keeping R much smaller than Q is what makes Taylor recovery cheap.

**8. From the core case to general #k-SAT.** Counting-preserving sparsification writes the count of any k-CNF as a sum of at most poly · 2^(σn) counts of formulas with boundedly many occurrences per variable. Each is converted to a 3-CNF with at most three occurrences per variable by replacing long clauses with chains of OR gates and replacing repeated variables with copies tied together by implication cycles; both steps preserve the count exactly. Choosing σ near (log n)^(−1/(2k)) already gives a 2^(o(n)) bound for every fixed k.

**9. The sharper bound.** Sparsification alone leaves a polylogarithmic occurrence bound at the right σ. The paper instead turns each leaf into a system of O(n) sparse low-degree polynomial equations over F_2 (parity tests of clause-violation counts, chosen deterministically using expander-code parity checks and a conditional-expectations argument), and counts common zeros of such systems with the same algorithm (Theorem 2.2). Appendix D recasts the whole count as contracting a Boolean tensor network.

The paper credits related ideas to split-and-list algorithms, Björklund's interpolation methods, and matrix-rigidity constructions that interpolate (−1)^j on a short interval.

## How the circuit lower bound is derived

Appendix C.1 of the paper, read in full:

- **Qualitative bound (Corollary C.1):** Williams proved in 2013 that if ETH is false then E^NP does not have linear-size circuits. Applying that implication to Theorem 1.1 gives E^NP ⊄ SIZE(O(n)). The connection is the general idea of deriving circuit lower bounds from faster satisfiability algorithms. The paper stresses it does not assert a superpolynomial lower bound for NP.
- **Quantitative bound (Corollary 1.4):** the paper reruns Williams's witness-circuit argument, combined with a local reduction of Jahanjou, Miles and Viola and a nondeterministic time hierarchy, using the constant-density case of Theorem 1.1. Only the decision version of the algorithm is used, so the strongest counting guarantee is not needed for these circuit bounds.

## Review and verification status

- The repository's [release note](history.md) says an independent working group of computer scientists led by Russell Impagliazzo, Ramamohan Paturi and Ryan Williams will examine and explain the proof, and that the materials were released immediately at the group's request.
- Lean checks the formal statement above (the weaker bound). As with the other results, it does not by itself show that the formal statement matches the intended claim, which the Comparator file and expert review address.
- I found no report of the working group's findings. [Runtime Wire](https://runtimewire.com/article/openai-ai-eth-disproof-external-review) (Oct 26) says the page does not say whether the review has started, and quotes no outside experts.

## Reported independent claim from Anthropic

Secondary sources report that Anthropic posted its own refutation, a paper titled "A 2^o(n) Time Algorithm for 3-SAT" credited to Claude, with a Lean 4 formalization. I have not read that paper, and neither source gives a link to it. What the pages report:

- **Provenance.** [OfficeChai](https://officechai.com/ai/anthropic-says-claude-too-independently-disproved-exponential-time-hypothesis-on-the-same-day-as-openai/) says an internal Claude research model found the core algorithm on October 16, with no human input in the discovery, using about 14 million output tokens, and that the paper acknowledges OpenAI's similar claim. [Runtime Wire](https://runtimewire.com/article/anthropic-posts-its-own-eth-refutation-hours-after-openai-as-lance-fortnow-says-) says the paper and Lean formalization were added to Anthropic's formal-math GitHub repository late on Monday, October 26.
- **Posting date differs.** [NEOTEO](https://www.neoteo.com/en/anthropics-3-sat-preprint-details-its-ai-claim-after-openais-math-announcement) says the preprint was deposited on October 27. The sources therefore disagree by about a day, and I can't tell which is right.
- **Lean.** OfficeChai reports about 1,050 Lean modules on Mathlib, checked by Lean's kernel and an independent checker, with the Lean code written by AI models and reviewed by other instances of the same models, and no person known to have checked it.
- **Expert reaction.** Runtime Wire's headline reports Lance Fortnow saying the approach "looks promising", while its summary has him unable to vouch for it yet. Neither OfficeChai nor Runtime Wire reports any other expert review.
- **Practical effect.** OfficeChai reports the paper saying the improvement is purely asymptotic and says nothing about solving formulas in practice, and does not advance P versus NP.

**Open discrepancy on the break-even point.** [CellCog](https://cellcog.ai/blog/openai-yang-mills-eth-kakeya/) says the bound beats brute force only past about 10^108 variables, and OfficeChai says roughly 10^10,000 for a typical 3-SAT formula. Runtime Wire and NEOTEO give no figure at all. I could not resolve this without the Anthropic paper, so none of these numbers should be relied on. For the OpenAI paper, the sections I read leave the constants implicit (the bounds are stated with O_k(·)), so a break-even point cannot be computed from the theorem statement alone.

## Connections to other pages

- **[COMPLEXITY.md](COMPLEXITY.md):** the other complexity-theory results in this repository, including entry 102 that this paper cites.
- **[KAKEYA.md](KAKEYA.md):** the ETH paper cites Björklund, Kaski and Williams for using generalized Kakeya sets (a finite-field notion) to recover polynomial values from structured evaluations. That is a different object from the Euclidean Kakeya sets of the harmonic-analysis result, and the repository claims no link between them.
- **[PHYSICS.md](PHYSICS.md):** nothing in the repository connects this result to the physics papers. As general background (not from the repository), counting satisfying assignments and counting perfect matchings are partition-function computations for constraint and dimer models, which is the setting of Corollary 1.3.
- **[TILING.md](TILING.md):** the repository's perfect-matching and dimer papers, including the approximate-counting result (entry 113) next to this paper's exact counting.

## Caveats

- The detailed proofs of sections 3 to 6 and Appendices A, B and D were not read. The algorithm description above is the paper's own overview and is not a check of the proof.
- The parameter values in step 6 were read from the paper's equation (2.10) as extracted from the PDF text, and I confirmed they are consistent with the paper's final bound, but I did not check them against the typeset page.
- The news items are from secondary sources, and Scientific American, Reddit and OpenAI's own announcement could not be retrieved when this page was written.
- The weaker Lean bound and the sharper manuscript bound should not be conflated.

## To do

- [ ] Add OpenAI's primary announcement link once it can be retrieved.
- [ ] Add the working group's write-up when it appears.
- [ ] Read the Anthropic paper directly, settle the break-even figure and the posting date.
- [ ] Read sections 3 to 6 and Appendix D of the OpenAI paper and check steps 4 to 9 against them.
- [ ] Add a section on how the Comparator statement is phrased.
