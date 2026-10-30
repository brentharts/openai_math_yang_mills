# Reconstruction, Decidability and Set Theory: Entries 009, 004, 027 and 297

A reading guide to four results that look unrelated but share a question: *how much does a small amount of data determine, and what does proving it cost?* Entry **009** recovers a function field from mod-ℓ Milnor K-theory or Galois data. Entry **004** proves that no algorithm decides whether an integer polynomial has a rational zero. Entry **027** shows integral points are potentially dense on curve character varieties. Entry **297** builds, in ZFC alone, a counterexample to Naimark's problem. It is a fork-level note, not an official OpenAI document, and it assesses no proof.

**How to read the evidence.** Claims are tagged:
- **[source]**: stated in a paper, read in its extracted text; section, theorem and page numbers are the printed ones. Lean facts come from `lean/docs/*.md` and from looking at the files on disk.
- **[mine]**: a comparison or observation of this fork.
- *(background)*: general knowledge, not from a source here.

**Numbering.** Entry numbers are those of [CONTENTS.md](CONTENTS.md); for 009, 004, 027 and 297 they match the numbers in the request that started this page.

## Summary

| Entry | Result | Real dependencies in the repository |
|---|---|---|
| **009** | Three papers: Bogomolov–Pop reconstruction from pro-ℓ abelian-by-central Galois data; a mod-ℓ Milnor K-theory version; a mod-characteristic version | [MK] adapts the method of [BP]; [MKp] cites no OpenAI paper |
| **004** | No algorithm decides rational solvability of integer polynomials (number of variables part of the input) | Needs entry 010's Fontaine–Mazur theorem and entry 004's own second paper (the pointwise 2-converse); may use entry 006's Goldfeld paper instead |
| **027** | Integral points are Zariski dense on every component of SL_r character varieties of curves, over one number field, with prescribed quasi-unipotent boundary classes | None |
| **297** | A unital, infinite-dimensional, simple C\*-algebra with a faithful tracial state and exactly one irreducible representation up to equivalence, in ZFC | None; credits an independent ZFC construction by Tanaka |

### Citation keys

| Key | Paper | Entry | Date | Pages |
|---|---|---|---|---|
| **[BP]** | [The Bogomolov–Pop reconstruction theorem](preprints/The-Bogomolov-Pop-reconstruction-theorem-September-23-2026/) | 009 | 23 Sep | 26 |
| **[MK]** | [Reconstruction of function fields from mod-ℓ Milnor K-theory](preprints/Reconstruction-of-Function-Fields-from-Mod-ell-Milnor-K-Theory-October-5-2026/) | 009 | 5 Oct | 36 |
| **[MKp]** | [Reconstruction from Milnor K-theory modulo the characteristic](preprints/Reconstruction-from-Milnor-K-theory-modulo-the-characteristic-October-5-2026/) | 009 | 5 Oct | 11 |
| **[H10]** | [Hilbert's tenth problem over the rational numbers](preprints/Hilberts-tenth-problem-over-the-rational-numbers-October-6-2026/) | 004 | 6 Oct | 62 |
| **[2C]** | [A pointwise 2-converse for elliptic curves with rational two-torsion](preprints/A-pointwise-2-converse-for-elliptic-curves-with-rational-two-torsion-October-7-2026/) | 004 | 7 Oct | 91 |
| **[CV]** | [Integral points on character varieties of curves](preprints/Integral-points-on-character-varieties-of-curves-September-25-2026/) | 027 | 25 Sep | 27 |
| **[NAI]** | [A counterexample to Naimark's problem in ZFC](preprints/A-counterexample-to-Naimarks-problem-in-ZFC-September-24-2026/) | 297 | 24 Sep | 23 |

Outside keys used below: [FM] the [Fontaine–Mazur modularity at the prime 2](preprints/Fontaine-Mazur-modularity-at-the-prime-2-October-6-2026/) paper (entry 010; see [LANGLANDS.md](LANGLANDS.md)), [GOLD] the [Goldfeld / 2-converse](preprints/Goldfelds-analytic-density-conjecture-and-the-2-converse-for-elliptic-curves-October-7-2026/) paper (entry 006).

## 1. Entry 009: function fields from Milnor K-theory and Galois data

### The three theorems [source]

All three concern finitely generated extensions K/k and L/l of transcendence degree at least 2 over algebraically closed fields.

| | Datum | Characteristic | Conclusion | Ambiguity |
|---|---|---|---|---|
| **[BP]** Thm 1.1 | the pro-ℓ abelian-by-central quotients Π^a, Π^c with the commutator structure ("bracket-compatible" isomorphisms) | ≠ ℓ | Φ : Isom^i_F(K,L) → Isom^c(Π^a_L, Π^a_K)/ℤ_ℓ^× is a **bijection**, where Isom^i_F(K,L) is the set of isomorphisms of perfect closures carrying k to l, modulo Frobenius powers | ℤ_ℓ^× scalars; Frobenius in characteristic p |
| **[MK]** Thm 1.1 | K^M_1/ℓ and K^M_2/ℓ **with their full product** (the relation space R_F spanned by [x]⊗[1−x]) | ≠ ℓ (ℓ any prime, including 2) | a bijection Isom^i_F(K,L) → Isom^M(V_K, V_L)/𝔽_ℓ^× | one global scalar in 𝔽_ℓ^× (trivial at ℓ = 2); Frobenius in characteristic p |
| **[MKp]** Thm 1.1 | K^M_1/p with the Steinberg relations | = p | a bijection Isom(K,k;L,l) → Isom_M(V_K,V_L)/𝔽_p^×, with the **field itself** recovered, not its perfect closure | one scalar; "There is no further Frobenius ambiguity" (Remark 5.3) |

[BP]'s abstract calls it a resolution of the Bogomolov–Pop reconstruction conjecture for these fields. The papers do not state a combined theorem covering ℓ ≠ char and ℓ = char together; reading them as complementary is my inference.

### How the proofs go [source]

- **[BP]** (§§3–7). Dualize a bracket-compatible isomorphism to a map of completed multiplicative groups; commutators correspond to "alternating characters" (Lemma 3.1, valid for ℓ = 2), which recover quasi-prime divisors; Pop's residual-curve criterion recognizes residual curves whose constants are algebraic over finite fields. The new geometric step is the **uniform pencil bound** (Theorem 5.1): for fixed t ∈ K∖k, the whole pencil Θ[t − a], a ∈ k, has uniformly bounded divisor-support degree. It uses Riemann–Hurwitz on a smooth residual curve of bounded genus with "five fixed members of the pencil"; "The use of five anchors works for ℓ = 2". Then an incidence argument produces curve subfields (Theorem 6.1), and §7 (Prop 7.1) finishes by quoting Pop's two reconstruction theorems as black boxes.
- **[MK]** has the same shape with finite coefficients but "uses no pro-ℓ field-reconstruction theorem as a premise" and never lifts the mod-ℓ isomorphism to a pro-ℓ one. It adapts [BP] §§4–6 (bounded support, incidence) and reconstructs field operations from simultaneous values on a blow-up.
- **[MKp]** is elementary by comparison: projective lines over the subfield F^p of p-th powers are recovered by a calculation with derivations (Prop 3.1 forces a unique exponent n with 1 ≤ n < p), a single scalar propagates over a connected graph (Prop 4.2), and the fundamental theorem of projective geometry (Prop 5.1) turns line-preserving maps into field isomorphisms.

### What it rests on [source]

| Input | Used for | Status in the paper |
|---|---|---|
| Pop, minimized decomposition theory (Bull. Math. Soc. Sci. Math. Roumanie 2015) with Topaz's appendix | residual-curve density criterion ([BP], [MK] Thm 3.1) | Published |
| Pop, *Recovering function fields from their decomposition graphs* (2012); Pop, *On the birational anabelian program initiated by Bogomolov I* (Invent. 2012) | [BP] Prop 7.1 (surjectivity) | Published; used as black boxes |
| Pop, *Reconstruction of function fields from their pro-ℓ abelian divisorial inertia*, 2018 | [BP] Prop 7.1(1) | Preprint, "Version 1" |
| Topaz, *Recovering function fields from their integral ℓ-adic cohomology…* (arXiv 2019, v1), and *A Torelli theorem for higher-dimensional function fields* | alternating-pair valuation theory ([MK] Thm 2.1); true-divisor criterion | Preprints; the Torelli paper is cited as "preprint, 2021, version 3" in [MK] and as JEMS 27 (2025) in [BP] |
| H. Tanaka, *Boundedness of regular del Pezzo surfaces over imperfect fields*, Prop 5.3 | finitely many Hilbert polynomials | Published (a different Tanaka from entry 297's) |

### Lean [source and observation]

`lean/docs/009.md` links **[BP] only** and says the formalization proves *injectivity*: if two isomorphisms of perfect closures preserving constants induce the same Galois-data isomorphism up to an ℤ_ℓ-unit, they agree modulo Frobenius. "Existence of a field isomorphism for every admissible Galois-data isomorphism is outside it." The challenge file `ComparatorChallenges/BogomolovPopInjectivity.lean` states `OAI.BogomolovPop.main_graph_injective` with `sorry`, which is the comparator stub; the solution directory `lean/OAI/Algebra/BogomolovPop/` has 5 files (about 3,300 lines) with no `sorry`, `axiom` or `native_decide` by search. **I did not build anything**; no build log is in the repository.

### Neighbours [mine]

- **Anabelian reconstruction (019, 031).** Entry 019 (the p-adic section conjecture) and entry 031 (Uchida's conjecture) also recover arithmetic objects from Galois groups. Neither cites 009 and 009 does not cite them. Entry 019 cites Pop's *birational* p-adic section conjecture and Pop–Stix localization, and entry 031 follows Neukirch–Uchida for number fields. The shared theme is Galois groups determining fields; the methods differ (valuation theory and Milnor K-theory for 009; Kummer theory and decomposition groups for 031).
- **Five branch points.** See §5.

## 2. Entry 004: Hilbert's tenth problem over ℚ

### What is proved [source]

[H10] Theorem 1.1 (p. 3): there is no algorithm which, given f ∈ ℤ[X₁,…,X_n] with n part of the input, decides whether f has a zero in ℚⁿ. Corollary 7.1 (p. 56): the rational-solvability set has Turing degree 0′, and so do two restricted problems: polynomials of total degree at most four, and "does Σqᵢ² have a zero in ℚᴺ" for lists of quadratics. No bound is placed on variables, summands or coefficient heights.

What it **does not** claim: the proof gives "no single existential definition of ℤ in ℚ" and "does not establish many-one completeness" (proof of Cor 7.1). Mazur's conjecture, which obstructs Diophantine definitions, is not addressed.

### The mechanism [source], §1.2 (pp. 5–6)

The argument is *not* a Diophantine definition of ℤ. For each nonconstant integral polynomial f, an effectively generated sequence of finite tests is built, each answerable by finitely many rational-solvability queries, with

f has an integer zero ⟺ every test succeeds.

An integer zero supplies witnesses for all tests. If a rational-solvability algorithm existed, one could run two searches in parallel (enumerate integer zeros; enumerate tests until one fails), so this is a Turing reduction from H10(ℤ) to H10(ℚ).

- **Compactness (§2.4, Lemma 2.9).** If all finite tests succeed, there is a ring R containing a zero of f inside an elementary extension of ℚ (starred sorts, nonstandard integers). To each a ∈ R the tests attach a point T_a on E : y² = x³ − d²x over F = ℚ(√2), d = 75 − 53√2, a curve of rank one (Prop 2.1), additively: T_{a+a'} = T_a + T_{a'}. This gives an injective additive map η : R → \*ℤ with η(1) = 1 that need not respect multiplication.
- **Height contradiction (§2.7, Prop 2.12).** Put δ = f(η(a₁,…,a_n)). If δ = 0 an ordinary integer zero follows by elementarity. Otherwise δ is a nonstandard integer of polynomial size in B = max(1, |η(aᵢ)|), a uniform height bound gives h(Ã) ≪_f B, and canonical-height comparison makes h(x(η(a_j)P)) grow quadratically in |η(a_j)|, giving B² ≪_f B, a contradiction.
- **Where the height bound comes from (§4).** A finite cover of P¹ ramified of index 2 over **five** points of P¹(ℚ), carrying principally polarized abelian surfaces with quaternionic multiplication by the maximal order of the indefinite quaternion algebra of discriminant 210 (a Shimura curve). Fibers over rational points are isogenous to Galois conjugates; a Galois representation gives a weight-two newform of level bounded by a power of M_S(s), then an isogeny factor of a modular Jacobian over a field of degree O(M_S(s)), then height bounds by a quantitative isogeny theorem (Gaudron–Rémond).
- **Primes (§6).** A prime-pattern theorem (Green–Tao, Green–Tao–Ziegler inverse theorem with the 2024 erratum, Selberg sieve) supplies primes with prescribed behavior for the "ordinary integer model".
- **Parity of denominators (§3).** A separate existential condition excludes poles of odd order. Its passage from two-descent data to rational points uses the pointwise 2-converse and the Cassels–Tate pairing: for curves y² = x(x − ℓ)(x + 3ℓ) with Selmer corank 0 or 1 and finite Sha, the alternating perfect pairing makes dim 𝔽₂ Sha[2] even, so Sha[2] = 0 and the Selmer class comes from a rational point (Lemma 3.4).

### Dependencies [source]

```
[H10] ──► [2C]   pointwise 2-converse, Thm 1.1            (Lemma 3.4, §3)
[H10] ──► [FM]   odd regular Fontaine–Mazur at 2, Thm 1.1 (§4.7: makes a 2-adic representation modular;
                  "whose statement has no hypothesis on the residual image")
[H10] ··► [GOLD] Remark 3.5: the unrestricted 2-converse "can replace the pointwise theorem"
[2C]  ··► [GOLD] cited as the "broader companion"; "the broader converse is not an input to the proof"
```

Solid arrows are used; dotted arrows are alternatives or context. [FM] §1 in turn names [H10] §4 as a consumer of its Theorem 1.1, so the dependence is stated on both ends.

The date order is odd but explained: [H10] is dated 6 October and cites [2C], dated 7 October; the READMEs say these versions "update citations to revised manuscripts in this collection". The repository also holds the 24 September versions of both papers, which CONTENTS does not link.

[2C] (p. 2): for E/ℚ with E(ℚ)[2] ≠ 0 and Selmer corank s₂(E) ∈ {0, 1}, ord_{s=1} L(E,s) = rank E(ℚ) = s₂(E) and Sha(E/ℚ) is finite, for every reduction type including at 2. It "does not assert the full Birch–Swinnerton-Dyer leading-coefficient formula". Its inputs are classical (Kato, Ferrero–Washington, Howard–Kolyvagin, Gross–Zagier, Shimura–Waldspurger, Chevalley–Warning, Chaiken's matrix-tree theorem) plus constructions "within this paper"; §7 says its odd detector "remains conditional until Proposition 9.10", an internal ordering of results, not an external hypothesis.

### The elliptic-curve cluster [source unless marked]

| Entry | Role | Link to 004 |
|---|---|---|
| **004** [H10], [2C] | Undecidability, and the pointwise 2-converse it needs | The hub |
| **010** [FM] | Modularity at 2 | Used in [H10] §4.7; [H10] is [FM]'s named consumer |
| **006** [GOLD] | Goldfeld's conjecture and the unrestricted 2-converse | Optional replacement for [2C] in [H10] Remark 3.5 |
| **002** (full BSD from low Selmer corank; Selmer converse at every prime) | Corank 0/1 converse for every prime and the full BSD formula | Not cited by [H10]. Entry 019 cites the Selmer converse, "a further arithmetic test for an elliptic target over ℚ" |
| **019** (local p-adic section conjecture) | Anabelian | Consumes entry 002's Selmer converse, not 004 |

CONTENTS says that 002 together with 006 gives full BSD for a density-one set of quadratic twists. Entry 242 (single-fold Diophantine representations and an at-most-one-solution promise) is also about Hilbert's tenth problem, over the non-negative integers; [H10] does not cite it [mine: thematic only].

### Lean [source]

There is no `lean/docs/004.md` and no 004 challenge file, so nothing about [H10] or [2C] is formalized in the repository.

## 3. Entry 027: potential integral density on curve character varieties

### What is proved [source]

[CV] Theorem 1.1 (p. 2): let X be a smooth connected complex curve, r ≥ 1, and C₁,…,C_s quasi-unipotent conjugacy classes (all eigenvalues roots of unity) with K ⊂ ℂ a number field containing their eigenvalues. There is a finite extension L/K such that M_r(X)(O_L) ∩ Y(C) is Zariski dense in Y, hence in every irreducible component. The same holds with only characteristic polynomials prescribed, and with no boundary conditions (then K = ℚ). For projective X and for the relative statements one may take F = K(√2, ζ_{pq}) for distinct primes p, q > r and L = F(ζ_r, {u^{1/r} : u ∈ O_F^×}).

"Integrality" here is *ambient*: a point of the affine ℤ-scheme M_r(X) over O_L whose complex semisimple representation has the prescribed boundary classes. It "need not avoid smaller Jordan classes after reduction at every finite prime". The paper says it settles the SL_r (determinant-one) curve case of Litt's integral-density question and the SL_r curve case of the Coccia–Litt conjecture for quasi-unipotent boundary data; "the higher-dimensional, higher-rank question remains outside its scope". The non-projective, no-boundary case is already dense over ℤ (free fundamental group); the substantive cases are projective curves and prescribed boundary monodromy.

### The proof [source], §1.2 (p. 4)

Work with GL_n; a boundary matrix of given Jordan type has an invariant flag with scalar root-of-unity quotients. Build *surface diagrams* of facets of various ranks glued along seams by flags. The diagram density theorem (Theorem 4.1) is by induction on the largest rank, cutting a maximal-rank facet to disks and replacing each disk by a subspace and quotient of lower rank. Closed surfaces enter the same step by cutting along a nonseparating circle (Lemma 7.1); genus 1 by density of commuting pairs (Motzkin–Taussky, Gerstenhaber, Guralnick–Sethuraman); genus ≥ 2 by irreducibility of {Π[A_j, B_j] = cI} of dimension (2g−1)n² + 1 (Prop 6.3), via Liebeck–Shalev character counting and Lang–Weil; Dehn-twist powers enlarge the Zariski closure. In §8 determinants are made units and an r-th-root extension gives determinant one, and the exact Jordan classes are selected *after* semisimplification. Three standard inputs: finite generation of integral invariants (§2), Dirichlet's unit theorem (§3), Lang–Weil (§6).

### Lean [source and observation]

`lean/docs/027.md` is explicit that the formalization proves only "a supporting compatibility result for its marked-surface construction": given a diagram with a chosen seam of maximal adjacent facet rank and a positive-genus parent, the produced marked values satisfy every seam equation. "Existence of integral points and Zariski density are not established. The paper's main theorem … is outside this comparison." The tree `lean/OAI/AlgebraicGeometry/CharacterVarieties/` has about 29,600 lines and no `sorry` by search. `lean/formalization.yaml` lists this paper. I could not map the Lean lemma to a numbered statement of the paper, and I did not build anything.

### Neighbours [mine]

Entry 043 (P = W for fixed-determinant SL_n moduli) and entry 014's tame Hecke eigensheaf papers also work on character varieties and local systems with unipotent monodromy on punctured curves (see [LANGLANDS.md](LANGLANDS.md#4-links-across-the-three-entries)). No paper cites any of the others. [CV] cites Coccia–Litt, Simpson, Esnault–Groechenig and de Jong–Esnault for the integrality side.

## 4. Entry 297: Naimark's problem in ZFC

### What is proved [source]

[NAI] Theorem 1.1 (p. 2): there exists a unital infinite-dimensional simple complex C\*-algebra A with a faithful tracial state such that all its nonzero irreducible representations are unitarily equivalent; in particular A is not isomorphic to K(H) for any Hilbert space H. Naimark's question (1951) asks whether the property "every nonzero irreducible representation is equivalent to the defining one" characterizes K(H). "No separability assumption is placed on A." The proof uses only ZFC: ordinary choice and transfinite recursion with set-sized codes, no diamond, no continuum hypothesis, no forcing or large cardinals.

The CONTENTS entry says "the counterexample is nonseparable". **The paper never says so.** It follows from the cited result of Calderón–Farah [Prop 7.2] (if the unique irreducible representation acts on a separable Hilbert space then A ≅ K(H)) together with Theorem 1.1. I recorded this as an inference; I found no use of the word "nonseparable" in the paper.

### The construction [source], "The proof in outline" (pp. 2–3), §§3–6

- **Bridge (Theorem 4.1).** From a faithfully tracial algebra A with no finite-dimensional irreducible representations and two decreasing projection sequences peaking at inequivalent pure states f₁, f₂ with equal traces tending to 0, form a reduced HNN extension H(A) (Ueda's construction via a reduced amalgamated free product). Every pure state of A then keeps a unique extension, and two extended pure states become equivalent iff the old ones were equivalent or form the pair {f₁, f₂}. Faithful trace, a trace-preserving conditional expectation and simplicity pass to H(A). Tools: pure-state excision (Akemann–Anderson–Pedersen) and cross-excision (Akemann–Weaver Lemma 1).
- **Tower (§5–6).** Start from A₀ = ⊗_{m≥1} M₂ (the Lean files call it `CAR`); repeatedly join a chosen pair of pure-state equivalence classes by a bridge, taking norm-closed unions at limits; the length is κ = 𝔠⁺. *Countable determination* (Lemma 5.6, Corollary 5.7): every pure state of the tower is the unique extension of its restriction to a separable dependency-closed subalgebra on which it is pure, so a tower of length ≤ 𝔠 has at most 𝔠 pure states. This replaces "the state-guessing role of a set-theoretic prediction principle". Kishimoto–Ozawa–Sakai homogeneity (Theorem 6.1) is applied only to separable subalgebras. A fixed injection schedule σ : κ × κ → κ, chosen before building, gives every pure state a stage at which it is bridged to the reference state.
- The paper says its distinctive features are "the selective bridge for all old pure states and its compatibility with dependency-closed subiterations". It makes no claim about the density of the resulting algebra.

### Relation to earlier work [source]

Akemann–Weaver (PNAS 2004) built a counterexample under Jensen's ◇_{ℵ₁} and showed independence of the existence of an ℵ₁-generated counterexample. Calderón–Farah replaced ◇ by a Cohen-forcing version of it together with CH (their Theorem 5.4). Vaccaro studied the trace spaces of counterexamples. The paper credits **Tanaka's earlier ZFC construction** (arXiv:2609.26930, 22 September 2026) as already giving the existence conclusion: "Tanaka gives another ZFC construction with a faithful tracial state, so the existence conclusion of Theorem 1.1 is already present there. His construction starts from the CAR algebra and uses strong shell sums and rank-one defects." Tanaka's algebra also has a faithful representation on a separable Hilbert space, which is not irreducible. [NAI] presents an alternative proof. The repository summary for 297 says the same.

### Lean [source and observation]

The Lean note states the full ZFC theorem. `ComparatorChallenges/Naimark.lean` states `OAI.NaimarkZFC.main` with `sorry` (the comparator stub); the solution tree `lean/OAI/Analysis/Naimark/` has 16 files (about 25,900 lines) with no `sorry`, `axiom` or `native_decide` by search. An intermediate hypothesis `KOS.SimpleSeparableRows` appears in `Recursion.lean`; I did not trace how it is discharged for `main`. `formalization.yaml` has no entry for this paper. No build was run.

### Cross-references inside the paper

Statements are numbered "Lemma n.m" but cited as "Theorem n.m" in the proofs, and local labels in the proof of Theorem 4.1 collide with the global Lemmas 2.5 and 2.6. Cosmetic, but it makes cross-checking harder.

### Neighbours [mine]

- **ZFC versus independence.** Entry 297 turns a ◇-conditional counterexample into a ZFC one. Entry 323 shows the *separable quotient problem* is independent of ZFC (relative to a measurable cardinal; counterexamples under CH); entry 244 builds a model where the Partition Principle holds and Choice fails; entry 240 proves Shelah's eventual categoricity conjecture in ZFC; entry 185 disproves infinite matroid intersection and packing/covering in ZFC. These are thematic: none cites [NAI] and [NAI] cites none of them. A text search of all extracted PDFs finds no other paper mentioning Naimark's problem (the other "Naimark" hits are Gelfand–Naimark–Segal and Stone–Naimark–Ambrose–Godement).
- **Operator algebras.** Entries 285 (Baum–Connes and Kadison–Kaplansky counterexamples), 286, 300 and 302 are also about C\*- and von Neumann algebras. I did not look for citations.

## 5. What the four have in common

| Theme | 009 | 004 | 027 | 297 | Strength |
|---|---|---|---|---|---|
| **A finite or bounded witness replaces a global object** | bounded divisor-support degree for a whole pencil (Thm 5.1) | finite tests plus compactness (§2.4) | diagram density by induction on rank | countable determination of pure states (Cor 5.7) | [mine], structural analogy |
| **Five branch points** | Riemann–Hurwitz with five members of a pencil; c_ℓ = 5(1−1/ℓ) − 2 | cover of P¹ ramified of index 2 over five points of P¹(ℚ) | none | none | [mine], see below |
| **Compactness or limits to remove a finiteness assumption** | no cardinality hypothesis on constants | elementary extension of ℚ | none | norm-closed unions along a tower of length 𝔠⁺ | [mine] |
| **Depends on other repository papers** | within the entry only | yes (010, and its own second paper) | no | no | [source] |
| **Lean coverage** | injectivity only | none | a supporting lemma only | the full statement | [source] |

### The five-branch-point coincidence [mine]

In [BP] the support bound comes from Riemann–Hurwitz applied to five fixed members of the pencil, with constant c_ℓ = 5(1 − 1/ℓ) − 2, which equals ½ for ℓ = 2 ("The use of five anchors works for ℓ = 2"). In [H10] §4 the height bound comes from a cover of P¹ ramified of index 2 over exactly five points. *(Background)* Index-2 ramification over m points gives orbifold Euler characteristic 2 − m/2, which is zero for m = 4 (an elliptic or Euclidean situation) and negative, a hyperbolic one, first at m = 5. That is the only common reason both papers need five points, and each proof uses its own machinery. The papers do not cite each other and share no lemma; this is a numerical coincidence with a standard explanation, not a connection between the arguments.

## Evidence status

| Statement | Status |
|---|---|
| Theorem statements, abstracts, proof outlines and bibliographies of the seven papers | [source], read in extracted text |
| Which results [H10] depends on, and in which sections | [source] (citation contexts checked against the text) |
| Details of [H10] §§3–6 proofs, [2C] §§3–9, [BP] §§4–6, [NAI] §§3–5 | **Not read line by line** |
| Correctness of any proof | **Not assessed** |
| Lean facts | [source] from the notes; file counts and `sorry` search are observations; **no build was run** |
| Nonseparability of the entry-297 algebra | Inference, not stated in the paper |
| Thematic links and the five-branch-point remark | [mine] |

## Caveats

- Several 009 inputs are preprints (Pop 2018 v1, Topaz 2019 v1, Topaz's Torelli paper), and the Topaz Torelli paper is cited with different bibliographic data in [MK] and [BP].
- [H10] relies on a long chain (Shimura curve, modularity at 2, isogeny bounds, prime patterns) and does not say which steps are effective. It also relies on [FM], which in turn relies on an outside 2026 preprint (Thorne's Theorem D); see [LANGLANDS.md](LANGLANDS.md#what-it-rests-on-source).
- [2C] is 91 pages with an ultrafilter-stabilization argument in §9 that I did not examine.
- Page numbers are printed page numbers.

## Open questions

- Is there a single statement covering ℓ ≠ char and ℓ = char in 009, with the Frobenius ambiguities reconciled?
- [H10] proves undecidability by a Turing reduction. Does any entry give a Diophantine definition of ℤ in ℚ? (The papers say no.)
- Can [CV]'s result be pushed to higher-dimensional varieties, or to the non-quasi-unipotent boundary classes covered in rank two by Coccia–Litt?
- How does the density of the [NAI] algebra relate to its length 𝔠⁺? The paper makes no claim.

## Related pages

- [LANGLANDS.md](LANGLANDS.md) for entry 010 (used by 004), 014 and 030
- [COMBINATORICS_REPRESENTATIONS.md](COMBINATORICS_REPRESENTATIONS.md) for the third group of new entries
- [COMPLEXITY.md](COMPLEXITY.md) and [ETH.md](ETH.md) for the other decidability and complexity results
- [CONTENTS.md](CONTENTS.md) for the full list of entries
