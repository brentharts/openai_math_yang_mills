# The Langlands Program in the Preprints: Modularity at 2, Imaginary Quadratic Fields, Restricted Geometric Langlands

A reading guide to the three Langlands-program results in this repository: entry **010** (modularity at the prime 2), entry **030** (modularity of elliptic curves over imaginary quadratic fields) and entry **014** (restricted geometric Langlands, Arthur enhancements and Ramanujan). It is a fork-level note, not an official OpenAI document, and it checks nothing about correctness. It records what the papers state, how they depend on each other and on 2026 outside preprints, and which links to other entries are real citations and which are only thematic.

**How to read the evidence.** Claims are tagged:
- **[source]**: stated in a paper, read in its extracted text. Section, theorem and page numbers are the printed ones. Where an entry has several papers, I read the abstract, introduction, every numbered statement and the proof outline of each. For the 72-page Fontaine–Mazur paper I did not read every proof in §§3.3–3.6, 5, 6 and 7.3–7.4.
- **[mine]**: a comparison or observation of this fork.
- *(background)*: general knowledge, not from a source here.

**Numbering.** Entry numbers are those of [CONTENTS.md](CONTENTS.md). For 010, 014 and 030 they agree with the numbers in the request that started this page. (For the combinatorics entries they do not; see [COMBINATORICS_REPRESENTATIONS.md](COMBINATORICS_REPRESENTATIONS.md#a-numbering-note).)

**Lean.** None of the three entries has a Lean note: `lean/docs/010.md`, `014.md` and `030.md` do not exist.

## Summary

| Entry | Result | Setting | Depends on (inside the repo) |
|---|---|---|---|
| **010** | Fontaine–Mazur modularity at p = 2, with no residual-image hypothesis. A dimension theorem for the 2-adic Hecke algebra. "Pro-modularity" of every odd absolutely irreducible 2-adic representation. | Galois representations of G_ℚ | Three papers; [PM] uses [FM] and [DIM] |
| **030** | Every elliptic curve over every imaginary quadratic field is modular, with matching local parameters. | Number fields | None; uses Caraiani–Newton as an outside input |
| **014** | Restricted geometric Langlands (full support), Arthur filtration rationality, Ramanujan–Arthur decomposition, global Arthur enhancements, generic Ramanujan, tame Hecke eigensheaves. | Function fields and curves over F̄_q | Eight papers in two chains (below) |

Outside the Langlands cluster, **entry 004** (Hilbert's tenth problem over ℚ) uses the 010 Fontaine–Mazur theorem in a height argument. That link is the only citation from a result outside these entries into them, and it is the subject of [RECONSTRUCTION_AND_DECIDABILITY.md](RECONSTRUCTION_AND_DECIDABILITY.md).

### Citation keys

| Key | Paper | Entry | Date |
|---|---|---|---|
| **[FM]** | [Fontaine–Mazur modularity at the prime 2](preprints/Fontaine-Mazur-modularity-at-the-prime-2-October-6-2026/), 72 pp. | 010 | 6 Oct |
| **[DIM]** | [The dimension of the two-adic Hecke algebra at odd level](preprints/The-Dimension-of-the-Two-Adic-Hecke-Algebra-at-Odd-Level-October-5-2026/), 12 pp. | 010 | 5 Oct |
| **[PM]** | [Unrestricted pro-modularity at the prime two](preprints/Unrestricted-pro-modularity-at-the-prime-two-October-6-2026/), 29 pp. | 010 | 6 Oct |
| **[ECIQ]** | [Modularity of elliptic curves over imaginary quadratic fields](preprints/Modularity-of-elliptic-curves-over-imaginary-quadratic-fields-October-4-2026/), 30 pp. | 030 | 4 Oct |
| **[RAM]** | [Ramanujan–Arthur decompositions of cuspidal functions at full finite level](preprints/Ramanujan-Arthur-Decompositions-of-Cuspidal-Functions-at-Full-Finite-Level-September-24-2026/), 64 pp. | 014 | 24 Sep |
| **[ARTH]** | [Global Arthur enhancements of cuspidal excursion parameters](preprints/Global-Arthur-Enhancements-of-Cuspidal-Excursion-Parameters-October-5-2026/), 69 pp. | 014 | 5 Oct |
| **[TEMP]** | [Temperedness at ramified places for globally generic exceptional groups](preprints/Temperedness-at-ramified-places-for-globally-generic-exceptional-groups-October-5-2026/), 19 pp. | 014 | 5 Oct |
| **[RAT]** | [Rationality of the canonical unramified Arthur filtration](preprints/Rationality-of-the-Canonical-Unramified-Arthur-Filtration-September-24-2026/), 43 pp. | 014 | 24 Sep |
| **[RGL]** | [The restricted geometric Langlands equivalence in positive characteristic](preprints/The-Restricted-Geometric-Langlands-Equivalence-in-Positive-Characteristic-September-24-2026/), 39 pp. | 014 | 24 Sep |
| **[CT]** | [Constructible tame Hecke eigensheaves in positive characteristic](preprints/Constructible-tame-Hecke-eigensheaves-in-positive-characteristic-October-5-2026/), 46 pp. | 014 | 5 Oct |
| **[SEV]** | [Tame Hecke eigensheaves with several marked points](preprints/Tame-Hecke-Eigensheaves-with-Several-Marked-Points-October-5-2026/), 40 pp. | 014 | 5 Oct |
| **[FROB]** | [Frobenius structures on tame Hecke eigensheaves](preprints/Frobenius-Structures-on-Tame-Hecke-Eigensheaves-October-5-2026/), 28 pp. | 014 | 5 Oct |

The papers label each other inconsistently: "[Ope26]" means [RAM] in [ARTH] but [RAT] in [RAM]. The keys above are this page's own.

## 1. Entry 010: modularity at the prime 2

### What is proved [source]

- **[FM] Theorem 1.1 (p. 1).** Let r : G_ℚ → GL₂(Q̄₂) be continuous and irreducible, odd, unramified outside finitely many primes, and de Rham at 2 with distinct Hodge–Tate weights. Then r is a Tate twist of the 2-adic representation of a classical cuspidal eigenform. The abstract adds "including all residual representations", and §1 says scalar and reducible residual representations are covered. It is the odd, regular, two-dimensional Fontaine–Mazur conjecture over ℚ at 2.
- **[DIM] Theorem 1.1 (p. 1).** For every odd N, every irreducible component of Spec 𝕋₂(N) has Krull dimension exactly 4, where 𝕋₂(N) is the inverse limit over weights k of ℤ₂ ⊗ (the ℤ-algebra generated by T_ℓ and ℓS_ℓ acting on M₁ ⊕ … ⊕ M_k of level Γ₁(N)), Eisenstein forms included. The paper calls this the p = 2 case of Emerton's dimension conjecture.
- **[PM] Theorem 1.1 (p. 1).** Let r : G_ℚ → GL₂(E) be continuous, absolutely irreducible and odd, unramified outside finitely many primes. Then for some odd N, divisible by every odd prime where r ramifies, there is a continuous ℤ₂-algebra map λ : 𝕋₂(N) → O_E with λ(T_ℓ) = tr r(Frob_ℓ) and λ(ℓS_ℓ) = det r(Frob_ℓ) for ℓ ∤ 2N. The paper calls this *pro-modularity*. There is no de Rham hypothesis, and the abstract says it "does not assert classical modularity". The level N may contain auxiliary tame primes.

The headline of the CONTENTS entry for 010 is a conflation: its second sentence ("also proves classical modularity up to Tate twist …") describes [FM], not [PM].

### How the three papers fit [source]

```
[DIM]  dimension 4 of every component           ─┐
[FM]   regular de Rham ⇒ classical (Thm 1.1)     ├─►  [PM]  every odd abs. irreducible r occurs in 𝕋₂(N)
Newton–Thorne, adjoint Selmer vanishing (Thm 5.4) ┘
```

- [PM] cites [FM] Theorem 1.1 to make one auxiliary point classical (Prop 4.5, p. 19–20), and borrows [FM]'s field-preparation, connectedness and propagation lemmas (Prop 2.3, pp. 6–9). It says that [FM]'s Theorem 7.1 "assumes a regular de Rham target", so its Prop 2.3 re-runs the argument rather than quoting a theorem.
- [PM] cites [DIM] Theorem 1.1 and its Prop 4.1. The last step (§6, pp. 24–27) compares two dimensions: at the auxiliary classical point the global deformation problem has tangent dimension at most 3, while "a dimension-four integral Hecke component has local dimension three" there. Equality makes the local ring regular and kills the Hecke kernel, which then vanishes on the whole irreducible trace family D. So the regular de Rham point is needed only "at the auxiliary point".
- [FM] and [DIM] cite none of the other 010 papers. [FM] cites no OpenAI paper except entry 004's, as a consumer (§1, p. 1).

### The three stages of [FM] (§1.2, p. 3–4) [source]

1. **Potential occurrence for a whole family.** *Localized propagation* (Theorem 5.1) starts from a characteristic-2 curve in a deformation space; if the curve is pro-modular, patching with "finite jets" propagates pro-modularity to every irreducible deformation locus containing it, with the dyadic deformation problem unrestricted. Curves are produced by abelian 2-extensions with controlled ramification; in the scalar and reducible cases a small proportion of dyadic places are *marked* with an invariant line. A Selmer calculation (§6) shows the globally reducible locus is too small to disconnect the charts, and a modular seed (§7.1, using Serre's conjecture and Billerey–Menares for the reducible case) and a curve cut from a component containing the target lie in one connected residual fiber (Theorem 7.1).
2. **Classical vectors at the target** (§8). Locally irreducible target: local Langlands gives locally algebraic vectors. Locally reducible target: take the family of trace plus the larger-exponent character, finite over a one-dimensional weight space. A dense set of ordinary points with Hodge–Tate weights {0, 1} is made modular by **Thorne's theorem [59, Theorem D]**, which is an outside 2026 preprint (arXiv:2608.07186). The paper then builds a genuine eigenvector and shows it is locally algebraic.
3. **Descent.** Prime-degree cyclic base-change descent along a solvable tower, correcting finite-character discrepancies by Hecke-character twists. Descent is applied only after classicality, "no descent statement for pro-modular loci is needed".

### What it rests on [source]

| Input | Where used | Status in the paper |
|---|---|---|
| Thorne, *Towards the Fontaine–Mazur conjecture for GL₂*, arXiv:2608.07186v1 (7 Aug 2026), Theorem D | [FM] Lemma 8.5, p. 63, the locally reducible branch | Outside preprint, version 1 |
| Pan (JAMS 2022), Zhang (arXiv:2412.06812) for odd primes | [FM] p. 1: "Together with the odd-prime results … it gives the stated modularity conclusion at every prime" | The 2-adic theorem itself does not depend on them; the "every prime" sentence does |
| Newton–Thorne, adjoint Selmer groups (JEMS 2023), Theorem 5.4 "specialized to GL₂/ℚ" | [DIM] Lemma 5.1, [PM] Prop 6.1 | Published; the specialization to p = 2 is asserted, not reproved |
| Emerton, *p-adic families of modular forms*, Corollary 2.28 (dimension ≥ 4) | [DIM] Theorem 2.4 | Cited; only the pseudodeformation step is re-justified at p = 2 |
| Paškūnas–Tung, block finiteness incl. scalar blocks at p = 2 | [PM] Theorem 3.2, [FM] Prop 4.2 | Published (Forum Sigma 2021) |

If Emerton's Corollary 2.28 were unavailable at p = 2, [DIM]'s statement would reduce to "dimension at most 4" (my reading of the paper's structure: the lower bound is cited from Emerton, not proved there). Several lemmas are numbered "Theorem" but are restated external results ([DIM] 2.4, [PM] 3.2).

### Notes on versions [source]

The repository also holds earlier versions that CONTENTS does not link: a 23 September Fontaine–Mazur paper (72 pp., same abstract) and a 4 October pro-modularity paper. The 6 October [FM] text says its Theorem 1.1 "supplies the weight-two modularity input in the height argument of Hilbert's tenth problem over the rational numbers [44, Section 4]".

## 2. Entry 030: modularity over imaginary quadratic fields

### What is proved [source]

[ECIQ] Theorem 1.1 (p. 1): for K imaginary quadratic and E/K any elliptic curve there is an isobaric automorphic representation π_E of GL₂(𝔸_K) with trivial central character and weight zero such that for every prime ℓ and every ι : Q̄_ℓ → ℂ, r_ι(π_E) ≅ H¹_ét(E_{K̄}, Q̄_ℓ), with the stated equality of local parameters at every finite place away from ℓ and the stated unitary parameter at the complex place. Hence L(E/K, s) = L(π_E, s − ½). π_E is cuspidal if E has no CM, and in the CM case is automorphically induced from a Hecke character (or an isobaric sum when the CM field is K).

### The "geometric prime switch" (§1.2, pp. 2–3) [source]

Fix a non-CM E/K and a large prime p with image of G_K on E[p] containing SL₂(𝔽_p) (the paper takes p > 100). Put n = 3p and F = K(μ_n). For a parameter t in a nonempty affine open, build a genus-n cyclic cover C_n(t) of a quadratic twist B_t of E, branched at two points. Its character pieces V_θ = H¹(C_n(t))_θ are rank-two compatible systems with Hodge–Tate weights 0, 1. With characters η, ψ of orders 3 and p and χ = ηψ, an integral description of the cohomology gives a chain of residual isomorphisms

V_η ≡ V_χ ≡ V_ψ ≡ V_1 (mod p, 3, p respectively), with V_η = H¹(A) for an elliptic curve A/F and V_1 = H¹(B_t).

Automorphy is lifted left to right in three steps (Prop 9.3), starting from the elliptic curve A, whose modularity uses the prime 5. Two degenerations of the cover supply the needed local structure: one with compact-type limit E^n gives a free integral group-ring description and lets one impose local reduction behavior on all pieces at once; the other has explicit vanishing cycles (opposite unipotents), giving absolute irreducibility and Frobenius elements in residual characteristic 3, 5 or p. Two specialization arguments finish: finite torsion is locally constant in a local family, and Hilbert irreducibility applied to the restriction of scalars of the parameter line from K to ℚ.

### What it rests on [source]

- **Caraiani–Newton** (arXiv:2301.10509v3, 2025), cited as [11]. The paper says (p. 2) that its Theorem 1.1 "removes the remaining restrictions on the curve and the imaginary quadratic field" and that "Caraiani–Newton's lifting theorem and residual-image criterion … remain essential inputs". Theorems 6.2–6.4 are restated, specialized versions of [11, Theorem 5.2, Theorem 6.1, Lemma 6.1.3]. The digest I made flags that the specialization, including the ℓ = 3 case with residual image SL₂(𝔽₃) and the ordinary-of-weight-zero check for quadratic twists of Steinberg, is *asserted* in the paper (p. 15) rather than proved.
- Local-global compatibility for ℓ ≠ p (Varma), Galois representations over CM fields (Harris–Lan–Taylor–Thorne, Scholze), cyclic base change (Langlands, Arthur–Clozel), Serre open image, Hilbert irreducibility with weak approximation.
- The CM case is a short separate argument (§9.4), not the cover construction.

No other repository paper is cited by [ECIQ], and none of the 010 papers cites it.

## 3. Entry 014: restricted geometric Langlands and the Ramanujan–Arthur line

CONTENTS lists eight papers under 014. They fall into two chains.

```
Arthur / Ramanujan line (function fields, automorphic side)

  [RAM]  Thm 1.1 decomposition of cusp functions by nilpotent orbit; Cor 1.2 unramified Ramanujan
    ├──► [ARTH]  "Assume Theorem 1.1": global Arthur enhancement          (formally conditional on [RAM])
    ├──► [TEMP]  Thm 1.1: Ramanujan at every place, types G2 F4 E6 E7 E8  (uses [RAM] Cor 1.2 as its Thm 6.1)
    └──► [RAT]   Cor 7.5 only (cuspidal comparison); [RAM] cites [RAT] back for context only

Geometric line (curves over F̄_q, sheaves)

  Gaitsgory–Raskin (arXiv:2508.02237, outside) ──► [RGL] full support, regimes A and B
                                                    ├──► [RAT] one remark (optional identification of the stack)
                                                    └──► [CT]  methods of §§3–7 adapted ──► [SEV], [FROB]
```
Source: the bibliographies and where each citation is used. An arrow means "is used by". [RAM] also cites [RAT] back for context only. Nobody cites [ARTH], and [SEV] and [FROB] do not cite each other.

### The statements [source]

| Paper | Main statement | Scope |
|---|---|---|
| **[RAM]** Thm 1.1 | For split connected semisimple G over 𝔽_q, any effective divisor D (arbitrary multiplicities) and any ℓ ≠ p, the cusp space C_D is the direct sum of subspaces C_{D,O} over nilpotent orbits O of the dual group, where each summand is governed by one orbit at every unramified place and every complex embedding. At D = 0 this is Gaitsgory–Lafforgue–Raskin Conjectures 3.4.5 and 3.4.6. | No hypothesis on p. Rational (ℚ) and Q̄_ℓ forms. |
| **[RAM]** Cor 1.2 | G split connected adjoint absolutely simple; a cuspidal π with one generic unramified place is tempered at every unramified place. | Proof in §6.5 combines the common orbit with the Ciubotaru–Harris local criterion. |
| **[TEMP]** Thm 1.1 | Globally generic cuspidal automorphic representations of split adjoint G of type G₂, F₄, E₆, E₇, E₈ over a global function field are tempered at *every* place. | Any characteristic, any ramification depth. Local criterion (Thm 5.2): a generic representation whose parameter has a pure weight-zero Weil–Deligne completion is tempered. |
| **[ARTH]** Thm 1.2 | Assuming the decomposition [RAM], each occurring cuspidal excursion character ν has an algebraic SL₂ and a commuting Weil-group map recovering the Lafforgue parameter by diagonal specialization on the whole Weil group, including inertia. | Arbitrary full level; "does not assert ellipticity, Arthur-packet classification, or a multiplicity formula". |
| **[RAT]** Thm 1.2 | The canonical Arthur filtration on finitely supported unramified automorphic functions is defined over ℚ (proves Gaitsgory–Lafforgue–Raskin's rationality conjecture 2.7.3 in this setting). | Needs four Lie-theoretic hypotheses on G_{F̄_q}. Unramified only. |
| **[RGL]** Thm 1.4 | The restricted geometric Langlands functor L^restr_G : Shv_Nilp(Bun_G) → IndCoh_Nilp(LS^restr_Ǧ) is an equivalence (Gaitsgory–Raskin's full-support conjecture). | Regime A: k = F̄_q, four hypotheses. Regime B: any algebraically closed field, four hypotheses, p very good and p ∤ |W_G|. |
| **[CT]** Thm 1.1 | Nonzero locally constructible perverse Hecke eigensheaves with Borel level at one marked point, for a Zariski-dense geometric ℓ-adic parameter with regular unipotent tame monodromy, genus ≥ 2, G simple simply connected over F̄_q. | Four hypotheses; "no Frobenius structure is required". |
| **[SEV]** Thm 1.1 | Same for SL_n, p > n, two or more marked points, any unipotent tame monodromy (Jordan type arbitrary, including zero). | Geometric only. |
| **[FROB]** Thm 1.1 | A Weil-sheaf version for SL_n, p > n, one marked point and regular nilpotent N, after passing to 𝔽_{q^m} for some m ≥ 1. | The integer m is not made explicit. |

### The four Lie-theoretic characteristic hypotheses [source]

They are the same four conditions in [RAT] (Assumption 1.1), [RGL] (Hypothesis 1.1) and [CT] (H1–H4): an invariant nondegenerate form that stays nondegenerate on the center of every Levi's Lie algebra; Chevalley restriction for every Levi; centralizers of semisimple elements are Levi subgroups; every nilpotent element lies in the nilradical of a parabolic defined over its field of definition. The first three come from Arinkin–Gaitsgory–Kazhdan–Raskin–Rozenblyum–Varshavsky §14.4.1, the fourth from their §D.1.1.

| Paper | Extra condition on p |
|---|---|
| [RAT] | None stated beyond the four conditions |
| [RGL] regime A | None ("No additional restriction on the order of the Weyl group") |
| [RGL] regime B | p very good for G and p ∤ \|W_G\|, on an arbitrary algebraically closed field |
| [CT] | "We impose no additional condition that p be prime to the order of a Weyl group" |
| [SEV], [FROB] | p > n; both papers *verify* the four conditions for SL_n (invertibility of n and of block factorials) |
| [RAM], [ARTH], [TEMP] | None; none of the four conditions appears |

### Observations from the cross-reading [source unless marked]

1. **[ARTH] is formally conditional.** Its Theorem 1.2 begins "Assume Theorem 1.1", the finite-level decomposition, attributed to [RAM]. [RAM] states its Theorem 1.1 unconditionally. But [ARTH]'s displayed "Theorem 1.1" is a *reformulation* (the class of σ_ν(Frob_x) is φ(diag(a^{d_x}, a^{-d_x}))·c_x with c_x compact) rather than [RAM]'s Theorem 1.1 (the decomposition of C_D). The bridge is [RAM] Theorem 5.1, Lemma 6.2 and Props 6.5–6.6, which [ARTH] does not cite at that point. The cited location is loose.
2. **[TEMP]'s argument does not appear to use the exceptional hypothesis.** Its Theorem 5.2 (local criterion) and Theorem 6.1 (the unramified input, a restatement of [RAM] Cor 1.2) are stated for every split adjoint absolutely simple group, and the only use of "G₂, F₄, E₆, E₇, E₈" is in Theorem 1.1's statement; classical groups are attributed to Lomelí's earlier work (p. 2). This is a reading of the text, not a claim that the theorem extends.
3. **[CT] seems to need less than its hypothesis says.** Its parameter is required to be *regular* unipotent (1.2), but p. 38 says "The descent argument uses unipotence of boundary monodromy, without a regularity requirement", and Prop 8.3 is stated for unipotent monodromy. [SEV] repeats that "unipotence, rather than regular unipotence, is sufficient". So [SEV]'s genuinely new content over [CT] is several marked points and the explicit SL_n, p > n verification, not arbitrary Jordan type as such. [FROB] again assumes regular N.
4. **No paper combines Frobenius structure with several marked points.** [SEV] is geometric; [FROB] has one marked point.
5. **[RGL] Corollary 9.5 is conditional.** It assumes a rational primed trace comparison over ℚ_ℓ that is "not constructed here", and Theorem 1.5 (rational spectral support) takes the rational spectral action as given data, for k = F̄_q only.
6. **Reading the notation.** Text extraction dropped overlines: "Fq" in the abstracts of [RGL] and [RAT] is F̄_q, and "Qℓ" is Q̄_ℓ. I confirmed this against the PDF pages for [RAT] p. 2 and [RGL] p. 4. Read those as F̄_q; CONTENTS states it correctly.
7. **Nobody names "the Arinkin–Gaitsgory conjecture".** [RGL] p. 2 describes characteristic-zero geometric Langlands as established by the five-part project of Arinkin, Beraldo, Campbell, Chen, Færgeman, Gaitsgory, Lin, Raskin and Rozenblyum, with Gaitsgory–Raskin's construction of the functor and multiplicity-one theorem as its first and last parts. [RAT] cites Arinkin–Gaitsgory only as the origin of nilpotent singular support.

### What 014 rests on outside the repository [source]

| Input | Used by |
|---|---|
| Gaitsgory–Raskin, *Geometric Langlands in positive characteristic from characteristic zero*, arXiv:2508.02237v1 (4 Aug 2025) | [RAT], [RGL], [CT], [SEV], [FROB] |
| Gaitsgory–Lafforgue–Raskin, *Tempered vs generic automorphic functions …*, arXiv:2603.26925v2 (2026) | [RAM], [ARTH], [RAT] (the conjectures being resolved) |
| Raskin, *Conjectures of Arthur and Ramanujan for unramified automorphic forms* (ICM 2026 announcement) | [ARTH], [RAT]. [ARTH] p. 4 notes that a global enhancement was announced there with irreducibility, under standing characteristic hypotheses; [ARTH] differs by working at arbitrary full level and not asserting irreducibility |
| V. Lafforgue (JAMS 2018) and the 2017 arXiv survey; Genestier–Lafforgue | [ARTH] cites the JAMS paper, [RAM] the arXiv survey, [TEMP] Genestier–Lafforgue |
| Ciubotaru–Harris, generalized Ramanujan and Arthur conjectures over function fields (arXiv 2023; Annals 204 (2026) per [TEMP]) | [RAM] Cor 1.2, [TEMP] |
| Nadler–Yun gluing; Dhillon–Taylor (2025 preprint); Hansen–Scholze | [CT], [SEV], [FROB] |

The digest also records that [RAM] cites the 2008 Bezrukavnikov–Finkelberg paper as "corrected arXiv version 4, 9 September 2026" and that PDF creation dates differ from printed dates for [ARTH], [RAT] and [RAM]. I did not verify these.

## 4. Links across the three entries

| Link | Kind | Strength |
|---|---|---|
| 010 [FM] → 004 (Hilbert's tenth problem over ℚ) | Citation: 004 §4.7 uses [FM] Theorem 1.1 "whose statement has no hypothesis on the residual image" to make a Galois representation modular | **Real dependency** [source]. 004 cites [FM]; [FM] §1 names 004 as a consumer |
| 010, 030, 014 with each other | None. Nobody cites across; [ECIQ] cites neither 010 nor 014, and 014 cites neither | No citations |
| 030 and 010: lifting by congruences | [ECIQ] chains congruences V_η ≡ V_χ ≡ V_ψ ≡ V_1 and lifts automorphy along them (a 3–5–p version of Wiles's switch); [PM] and [FM] propagate along connected families of deformations. Both are modularity-lifting in the broad sense, with different engines: congruence chains of a single geometric family versus patching on deformation rings | [mine], thematic |
| Zariski density and "not virtually solvable" | The same non-degeneracy condition appears throughout: [PM] and [FM] reduce to representations that are not virtually solvable; [ECIQ] needs SL₂(𝔽_p) in the residual image; [CT], [SEV], [FROB] require Zariski-dense parameters | [mine], thematic |
| Specialization between characteristics | [FM] propagates from a characteristic-2 curve; [CT]–[FROB] lift a curve and its parameter from characteristic p to mixed characteristic and take nearby cycles; [RAM] Lemma 2.1 changes characteristic by a strictly henselian trait | [mine], thematic. Different mathematics under a shared word |
| Number fields versus function fields | 010 and 030 concern G_ℚ and G_K (Galois side, deformation theory); 014 concerns split groups over curves over finite fields (automorphic and sheaf side, Lafforgue's excursion operators) | The two halves of the Langlands program in the preprints share no citations or lemmas |
| 014 and the Betti side | [CT]–[FROB] use the Betti spectral action and gluing of Nadler–Yun, and parameters that are local systems with unipotent monodromy on punctured curves. Entry 027 (integral points on SL_r character varieties of curves with prescribed quasi-unipotent boundary monodromy) and entry 043 (P = W) live on the same moduli spaces. Neither 027 nor 043 cites 014, nor conversely | [mine], thematic. See [RECONSTRUCTION_AND_DECIDABILITY.md](RECONSTRUCTION_AND_DECIDABILITY.md) for 027 |
| 014 and entry 069 (global quantum geometric Langlands at irrational level) | Both are geometric Langlands statements, 069 in characteristic zero for twisted D-modules. I did not read 069's paper (CONTENTS summary only) | [mine], thematic; no citation check made |

## Evidence status

| Statement | Status |
|---|---|
| Theorem statements, abstracts and proof outlines of the twelve papers | [source], read in extracted text |
| Dependency graphs (010 and 014) | [source], from the citations listed in each paper's bibliography and their use locations |
| Internal proofs of [FM] §§3.3–3.6, 5, 6, 7.3–7.4 | **Not read** |
| Correctness of any proof | **Not assessed** |
| Observations 1–5 in §3 | [source] for the quotes, a reading for the interpretation |
| Cross-entry thematic links in §4 | [mine] |

## Caveats

- Several key outside inputs are 2026 arXiv preprints or revised-in-2026 versions: Thorne's Theorem D (for [FM]), Caraiani–Newton v3 (for [ECIQ]), Gaitsgory–Raskin v1 and Gaitsgory–Lafforgue–Raskin v2 (for most of 014). I recorded their versions as printed in the bibliographies. I did not check them.
- The 010 theorems state results "at 2". The "every prime" claim in [FM] p. 1 relies on odd-prime results of Pan and Zhang that are not part of these papers.
- [ECIQ]'s reliance on Caraiani–Newton is through restated, specialized theorems. Whether the specialization is complete is a question for a specialist.
- Page numbers are printed page numbers; the PDFs begin counting at the title page.

## Open questions

- Does [TEMP]'s argument give the same conclusion for all split adjoint absolutely simple groups, given that its two key statements are not restricted to exceptional types?
- Is the regular-unipotent hypothesis in [CT] removable, as its own p. 38 and [SEV] suggest?
- Can a Frobenius structure be combined with several marked points ([SEV] plus [FROB])? The papers do not do it.
- How explicit is the extension degree m in [FROB]?

## Related pages

- [RECONSTRUCTION_AND_DECIDABILITY.md](RECONSTRUCTION_AND_DECIDABILITY.md) for entry 004, which depends on [FM], and for 009, 027 and 297
- [COMBINATORICS_REPRESENTATIONS.md](COMBINATORICS_REPRESENTATIONS.md) for the third group of new entries
- [YANGMILLS.md](YANGMILLS.md), [PHYSICS.md](PHYSICS.md), [HIERARCHIES.md](HIERARCHIES.md) for the physics-facing pages
- [CONTENTS.md](CONTENTS.md) for the full list of entries
