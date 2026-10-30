# Tiling Beyond Entry 155: A Reference Map

A fork-level reading guide that places entry 155 (a finite tile in ℤ³ with no fully periodic tiling) among twelve recent references on translational tilings, polygonal and connected tiles, undecidability, AI-assisted search, packings, a space-filling polyhedron and an optics experiment on an aperiodic monotile. It is not an official OpenAI document. Each section lists the references it uses and cites them by key. For entry 155's proof see [PERIODIC_TILING.md](PERIODIC_TILING.md); for the 3D aperiodic monotiles see [MONOTILES.md](MONOTILES.md) and [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md); for the wider picture see [TILING.md](TILING.md).

**How it was built, and what I could and could not read.** I fetched the abstract pages (or a mirror) of each reference and, for the Greenfeld survey, its HTML version. The fetch tool summarises pages, so quotes are short phrases from those summaries and I have not checked them against full texts. I did not read any proof. Entry 155 was read earlier from its extracted text. Claims are tagged **[source]** (stated in a reference as I saw it), **[repo]** (stated in a repository paper or page), **[mine]** (my comparison or inference) and *(background)* (general knowledge, not checked here).

## Reference key

| Key | Reference | What I read |
|---|---|---|
| **[Gr25]** | R. Greenfeld, [Translational tilings: structured or wild?](https://arxiv.org/abs/2509.25576), arXiv:2509.25576, 29 Sep 2025 (survey, math.CA) | HTML version: section list, Conjecture 3.8, Questions 5.4, 6.5, 7.4 and statements on ℤ³ and connected tiles, as summarised by the fetch tool |
| **[dGGM26]** | J. de Dios Pont, J. Grebík, R. Greenfeld, J. Madrid, [Translational tilings by axis-parallel polygonal sets](https://arxiv.org/abs/2610.15785), arXiv:2610.15785, 11 Oct 2026 | Abstract page only |
| **[dGGM-rat]** | The same four authors, [Periodicity and decidability of translational tilings by rational polygonal sets](https://arxiv.org/abs/2408.02151), arXiv:2408.02151, *Expositiones Mathematicae* 44(1):125620 (2026) | Not read. Cited as ref. [2] of entry 155; its Theorem 1.1 is quoted there |
| **[GK]** | R. Greenfeld and M. N. Kolountzakis, [Tiling, spectrality and aperiodicity of connected sets](https://arxiv.org/abs/2305.14028), arXiv:2305.14028; published as *Israel J. Math.*, May 2026, [DOI 10.1007/s11856-026-2916-5](https://doi.org/10.1007/s11856-026-2916-5) (journal data as you supplied it) | Abstract (arXiv pdf); the fetch page gave no journal reference |
| **[GT-JEMS]** | R. Greenfeld and T. Tao, [Undecidability of translational monotilings](https://doi.org/10.4171/jems/1673), *J. Eur. Math. Soc.*, online 21 Jun 2025; arXiv:2309.09504 | Abstract on the EMS Press page |
| **[GGTW]** | B. Georgiev, J. Gómez-Serrano, T. Tao, A. Z. Wagner, [Mathematical exploration and discovery at scale](https://arxiv.org/abs/2511.02864), arXiv:2511.02864 (v3, 22 Dec 2025) | Abstract page only |
| **[Tao-SD]** | T. Tao, [Sum-difference exponents for boundedly many slopes, and rational complexity](https://arxiv.org/abs/2511.15135), arXiv:2511.15135, Nov 2025 | **arXiv page rate-limited, not read.** I read Tao's blog post about it (republished on [Goodreads](https://www.goodreads.com/author_blog_posts/26222982-sum-difference-exponents-for-boundedly-many-slopes-and-rational-complex)); authorship is from that post |
| **[PS]** | C. Ma, T. T. Zhaowei, P. Li, M. Liu, H. Chen, Z. Mao, Y. Cheng, Y. Qi, Y. Yang, [Finding Kissing Numbers with Game-theoretic Reinforcement Learning](https://arxiv.org/abs/2511.13391), arXiv:2511.13391 (v5, 22 Sep 2026) | Abstract page only |
| **[Gyro]** | [Gyrobifastigium](https://en.wikipedia.org/wiki/Gyrobifastigium), Wikipedia | Page summary (about the polyhedron J₂₆) |
| **[Alv17]** | S. Alvarez, [The gyrobifastigium, not an uncommon shape in chemistry](https://doi.org/10.1016/j.ccr.2017.03.012), *Coordination Chemistry Reviews* 350 (2017), 3–13 | Abstract, via an index page |
| **[IJMS26]** | [Gyrobifastigium-based TPMS lattices with tuneable stiffness and energy absorption](https://www.sciencedirect.com/science/article/pii/S0020740326000767), ScienceDirect (the PII prefix S0020-7403 is *International Journal of Mechanical Sciences*, which I infer from the code, not from the page) | **Could not read.** Only the title appeared in a search result; the fetch failed |
| **[Mor26]** | Y. Moritake et al., [Chiral Diffraction from Aperiodic Monotile Structure](https://doi.org/10.1038/s41467-026-75023-7), *Nature Communications* 17 (2026) | Two press releases ([U-Tokyo IIS](https://www.iis.u-tokyo.ac.jp/en/news/5114), [ScienceDaily](https://sciencedaily.com/releases/2026/09/260914102432.htm)); not the paper |
| **[E155]** | OpenAI, [A translational tile with no fully periodic tiling in dimension three](preprints/A-translational-tile-with-no-fully-periodic-tiling-in-dimension-three-September-23-2026/paper.pdf), 23 Sep 2026 (entry 155) | Extracted text, read earlier |

## 1. The survey, and the question entry 155 answers

*References: [Gr25], [E155], [GT-JEMS].*

**What the survey covers [source, Gr25].** Sections on the domino problem and multi-tile undecidability (§2), discrete monotilings and the discrete periodic tiling conjecture (§3), the reversible structure theorem (§4), periodicity in dimensions one and two (§5), "wild" monotilings in virtually-ℤ² spaces via p-adic Sudoku and undecidability (§6), and continuous tilings (§7). It cites the hat and spectre only in the introduction, as examples that need rotations and reflections, so they lie outside its translation-only scope.

**Its open questions [source, Gr25].**

| Statement | Text as summarised | Where entry 155 stands |
|---|---|---|
| Conjecture 3.8 | Discrete periodic tiling conjecture: "the tiling equation F ⊕ X = G is not aperiodic". Known in ℤ, virtually-ℤ and ℤ²; false in virtually-ℤ² and, by reduction, in ℤᵈ for large d | Entry 155, Theorem 1.1(i), is a finite T ⊂ ℤ³ whose tilings are all non-periodic [E155] |
| Question 5.4 | "Does Conjecture 3.8 hold in ℤ³?" The survey says the ℤ² argument does not extend (Remark 5.3) | **Answered no, if entry 155's statement is read as the survey reads aperiodicity [mine].** I matched the two by wording; I did not check the survey's definition of "aperiodic" |
| Question 6.5 | Is there d ≥ 3 such that translational monotilings are undecidable in ℤᵈ? Known only when d is part of the input | **Not addressed.** A tile with no periodic tiling does not by itself give undecidability [mine] |
| Question 7.4 | Does the continuous conjecture (7.3) hold in ℝ²? "Still open" | Not addressed; entry 155 is dimension three (§2 below) |

**Timing [source and repo].** The survey is dated 29 September 2025 and entry 155 is dated 23 September 2026. Entry 155 is a repository (OpenAI) preprint and has not been refereed. The survey is therefore evidence for what was open a year earlier, and says nothing about entry 155. Entry 155's own introduction states the planar and lattice background that the survey sets out [E155, §1].

**The link through [GT-JEMS] [repo].** Entry 155 follows "the p-adic Sudoku formulation" of [GT-JEMS, Section 4] (its ref. [6]), and the finite-cyclic-factor reduction of Meyerovitch–Sanadhya–Solomon. The abstract of [GT-JEMS] says monotilings are undecidable in virtually-ℤ² spaces ℤ²×G₀ with G₀ a finite abelian group, and hence in ℤᵈ for d ≥ 3 [source]. The survey's §6.4 and §6.5 cover the same constructions [Gr25]. Entry 155 keeps G₀ cyclic and ends in ℤ³ [repo, PERIODIC_TILING.md].

## 2. Polygonal tiles: the plane is open, the cube-union case in ℝ³ is not

*References: [dGGM26], [dGGM-rat], [Gr25], [E155].*

**What the new paper proves [source, dGGM26 abstract].** For bounded open sets in ℝ² whose boundaries are finitely many axis-parallel segments:
- every *topologically minimal* tiling is weakly periodic;
- any such set that tiles the plane admits a *quasi-periodic* tiling, a finite union of translated lattices;
- as a consequence, the translational tiling problem is decidable for axis-parallel polygonal sets whose vertex coordinates are real algebraic numbers.

It uses geometric decompositions of tilings, rational approximation and discrete tiling theory. It states that the discrete ℤ² case was settled by Bhattacharya and that ℝ² is still open.

**What it extends [source].** [dGGM-rat] proved the existence of periodic tilings for *rational* polygonal sets, including disconnected sets and sets with holes. Entry 155 cites it as the positive planar Euclidean result and notes that it does not settle the conjecture for arbitrary bounded measurable tiles [E155, §1]. The survey states the same result and adds that "discretizing polygonal sets assuming only rational slopes is significantly more subtle" [Gr25, §7.1]. [dGGM26] concerns a different class (axis-parallel with real vertices) and a weaker conclusion (quasi-periodic, not periodic).

**Comparison with entry 155 [mine].**
- The ℝ³ statement in entry 155 is about Ω = T + [0,1]³, a finite union of unit cubes with integer corners. Such a set is the three-dimensional analogue of an "axis-parallel polygonal set" with rational vertices.
- For axis-parallel sets, [dGGM26] gives quasi-periodic tilings in dimension two. For the cube-union Ω in dimension three, entry 155 says no tiling has a translation set invariant under a full-rank lattice [E155, Thm 1.1(ii)]. So "tiles ⇒ has a (quasi-)periodic tiling" is available in the plane for this class and is, for *full* periodicity, false in ℝ³ for a set of this form.
- That is a statement about *fully periodic* tilings. [dGGM26] concerns weakly periodic and quasi-periodic tilings, and entry 155 does not claim those fail for Ω. I do not know whether Ω has a quasi-periodic tiling. Greenfeld–Tao's structure theorem in ℤ² (a tiling is a finite union of singly-periodic sets) is cited in [E155, §1]; the survey's reversible structure theorem is its generalisation to ℤᵈ [Gr25, §4], and I did not read how far it reaches in ℤ³.
- Decidability does not transfer. [dGGM26] decides tiling for algebraic axis-parallel polygons in ℝ². Whether tiling by Ω-type sets in ℝ³ is decidable is not addressed by any reference here; Question 6.5 of [Gr25] is the nearest open problem.

## 3. Connected sets and spectrality

*References: [GK], [E155], [Gr25].*

**What [GK] proves [source, abstract].** Both the periodic tiling conjecture and Fuglede's conjecture (tiling by translation versus an orthogonal basis of exponentials) fail in high dimension, but the known counterexamples are highly disconnected. The authors build connected counterexamples:
- a d-dimensional counterexample to the periodic tiling conjecture yields a connected one in dimension d+2;
- connected counterexamples to the "spectral ⇒ tiling" direction of Fuglede's conjecture in dimension d+2, and counterexamples to "tiling ⇒ spectral" in some dimension;
- the tools are "folded bridge" and "spiral bridge" constructions that join the components of a set while preserving its tiling, spectral and aperiodicity properties.

The survey records the same connected-tile result, "even when restricted to connected tiles" in high dimension, and a connected counterexample showing that the stronger version of the continuous conjecture with tiles required to be open and connected is false [Gr25, §3.4.1, §7.2].

**Comparison with entry 155 [mine].** Entry 155 states that it does not claim connectedness of the tile [E155; PERIODIC_TILING.md]. T is a finite subset of ℤ³, and T + [0,1]³ is a union of unit cubes that need not be connected. If the d → d+2 theorem of [GK] applies to an ℝ³ counterexample of this kind (I read only its abstract, not its hypotheses), it would give a connected counterexample in ℝ⁵. I have not checked that the hypotheses fit, so this is a question to ask, not a result. I also found no statement in any reference here about Fuglede's conjecture in dimension three.

**In the repository [repo].** No CONTENTS entry is about Fuglede's conjecture or spectral sets. There is one citation link: the paper behind entry 266 (exactly three mutually unbiased bases in dimension six; "Exact Fourier certificates for complex Hadamard matrices of order six") cites Tao's *Fuglede's conjecture is false in 5 and higher dimensions* (Math. Res. Lett. 11 (2004), 251–258) and says its cubic representative is "the matrix in Tao's construction [23, p. 254]". So the order-six Hadamard matrix that the entry works with is the one from Tao's counterexample to the spectral ⇒ tiling direction. The entry's paper uses it as a matrix, and does not discuss tilings or spectral sets in the text I searched. Other "Fuglede" hits in the repository are Fuglede–Kadison determinants (entry 197) and Fuglede's modulus (the Koebe circle-domain paper), which are unrelated.

## 4. AI-assisted search: AlphaEvolve and the sum-difference exponents

*References: [GGTW], [Tao-SD], [E155], and the Kakeya page [KAKEYA.md](KAKEYA.md).*

**[GGTW] [source, abstract].** An evolutionary coding agent, AlphaEvolve, pairs LLM-generated code proposals with automated evaluation. It is applied to 67 problems in analysis, combinatorics, geometry and number theory, matches the best known results in most cases and improves several. It is combined with Deep Think and AlphaProof for proof generation. The abstract claims that "LLM-guided evolutionary search can find mathematical constructions that complement human intuition".

**[Tao-SD] [source, via Tao's blog post].** This is the "second spinoff" of the AlphaEvolve project. It studies sum-difference constants, which are tied to Kakeya-type problems and to Shannon-entropy inequalities. For the slopes {0, 1, ∞} and a further rational slope a/b in lowest terms, the constant lies between 2 − C/log(2+|a|+|b|) and 2 − C′/log(2+|a|+|b|). AlphaEvolve's experiments had suggested the logarithmic decay and the paper confirms it. With more projections controlled, the parameter becomes a "rational complexity" of the slope and the constant decays to 2 at a polynomial rate whose exponent was not determined. The post says the results do not advance the arithmetic Kakeya conjecture, which concerns constants near 1. I could not read the arXiv page itself.

**Connections [repo and mine].**
- **A direct citation.** The paper "Asymptotically minimal maxima of real Littlewood polynomials" (under CONTENTS entry 076) cites [GGTW, §6.13, Problem 6.28] as having used AlphaEvolve on the same real-sign minimum up to degree 100, "their computations do not give an asymptotic theorem" [repo, paper text]. The matrix-multiplication papers (entry 107 and neighbours) cite a *different* AlphaEvolve paper, arXiv:2608.16884, so they are related only by tool.
- **Kakeya.** [KAKEYA.md](KAKEYA.md) is about the Euclidean Kakeya set and local smoothing (entries 074, 079). The sum-difference constants of [Tao-SD] are an arithmetic, entropy-flavoured relative, and the blog post explicitly says they do not advance the arithmetic Kakeya conjecture. I see no technical bridge to entry 074 [mine].
- **Tiling.** Neither reference is about tiling in its abstract or blog summary. The thematic link is slopes and projections: entry 155's engine is a word rule along every integer-slope line [E155], and [Tao-SD] works with a finite set of projection slopes and the arithmetic of one extra slope. This is a shared vocabulary, not a shared argument [mine]. The repository's trilinear Hilbert transform (entry 086) and ergodic averages "with distinct integer slopes" use the word "slope" in a third sense.

## 5. Packings: AI-found kissing configurations

*References: [PS], [GGTW]; repository entries 090 and the packings section of [TILING.md](TILING.md).*

**[PS] [source, abstract].** The kissing number problem asks how many non-overlapping spheres can touch a central one. The authors cast it as a cooperative matrix-completion game: one player fills cosine entries, the other corrects weak ones. Their reinforcement-learning system PackingStar improves 15 long-standing bounds on kissing numbers and related generalizations, "several provably optimal under natural inner products", and includes an explicit spherical-code realization of the Fischer group Fi22. The page lists v5 dated 22 September 2026 (first submitted 17 November 2025).

**Connections [mine].**
- Packing and tiling are dual: a tiling is a packing with no gaps ([TILING.md](TILING.md) says the same). The repository's packing results are the planar Cohn–Elkies certificate (entry 090, optimal density π/(2√3)) and symplectic ball packing (entry 343). No repository paper mentions kissing numbers (text search).
- [PS] and [GGTW] share a method, machine-guided search for extremal configurations, where the output is a certificate that can be checked by exact arithmetic. The repository's circle-packing paper also ships a verification folder for its certificate [repo, TILING.md]. This is a parallel in how results are checked, not a mathematical link.
- Neither [PS] nor [GGTW] is cited by entry 155, and I found no tiling statement in either abstract.

## 6. The gyrobifastigium

*References: [Gyro], [Alv17], [IJMS26], [Mor26 context], and the monotile pages [MONOTILES.md](MONOTILES.md), [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md).*

**The solid [source, Gyro].** The gyrobifastigium (Johnson solid J₂₆) is made by joining two triangular prisms along a square face with a quarter turn: 4 equilateral triangles, 4 squares, 14 edges, 8 vertices, symmetry group D₂d of order eight. Copies pack to form the gyrated triangular prismatic honeycomb. The page calls it "the only Johnson solid outside of the uniform polyhedra that can tile three-dimensional space" and one of five convex regular-faced space-filling polyhedra (the others are the cube, truncated octahedron, triangular prism and hexagonal prism).

**A correction to the one-line description.** "A unique 3D convex polyhedron that can tessellate" is too strong as written: cubes and many other convex polyhedra fill space. The safer form is the Wikipedia statement above: the only *Johnson solid* beyond the uniform ones that fills space.

**The nearby aperiodic object [source, Gyro].** The page says the Schmitt–Conway–Danzer (SCD) biprism is combinatorially equivalent to the gyrobifastigium (with parallelogram and irregular triangle faces) and "fills space only aperiodically or with screw symmetry", giving a partial answer to the 3D Einstein problem. It cites Senechal's *Quasicrystals and Geometry* in that section. This is the page's claim and I have not checked it.

**Chemistry [source, Alv17].** The author finds eight-atom groups shaped like a gyrobifastigium in inorganic A₄B₄ molecules, organic C₈ groups, tetranuclear M₄X₄ complexes and the eight-coordinate spheres of metal atoms. The ThSi₂ structure of many MX₂ compounds has M atoms in nearly ideal gyrobifastigia, high-pressure forms of K, Rb and Cs contain the shape, and several structures lie along a minimal-distortion pathway between it and the C₈ fragment of the diamond structure.

**Engineering [IJMS26].** The title is "Gyrobifastigium-based TPMS lattices with tuneable stiffness and energy absorption". I could not read it, so I do not know how it uses the polyhedron (TPMS usually denotes triply periodic minimal surfaces) or what it measures.

**Comparison with the repository's tiles [mine].**
- The gyrobifastigium is a convex tile with an ordinary periodic honeycomb. It illustrates the contrast the repository's pages keep making between a tile that fills space periodically and the three kinds of non-periodic behaviour: entry 155 (a tile whose tilings, under translations only, are never fully periodic), Chair44 (a non-convex strongly aperiodic monotile using isometries) and the SCD biprism (combinatorially a gyrobifastigium).
- The SCD biprism is the one place where a *convex*, gyrobifastigium-like shape appears in the aperiodic-monotile story, per [Gyro]. Chair44 is non-convex (a seven-cube body with pyramids) [repo, SPECTRE_CHAIR44.md]. I did not look for a connection between SCD and Chair44 beyond that.
- Entry 155's tile is a finite subset of ℤ³, "huge and not drawn" [repo, PERIODIC_TILING.md], and its tilings use translations only. The gyrobifastigium honeycomb is a statement about space-filling by congruent copies. I did not check whether translates of one orientation alone fill space, so I do not place it in the translational-tile setting [mine].

## 7. Optics on an aperiodic monotile

*References: [Mor26], [HR] and [SM] as cited in [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md), [MONOTILES.md](MONOTILES.md).*

**What the press releases say [source, Mor26 via U-Tokyo IIS and ScienceDaily].** A team led by Yuto Moritake with Masaya Notomi (University of Tokyo, Institute of Industrial Science) fabricated nanoscale patterns on silicon nitride films by electron-beam lithography, based on the "Smith hat", the aperiodic monotile found in 2023. Laser illumination gave pinwheel-like diffraction patterns that showed the structure's chirality, changing with the direction and polarization of the light, and mirrored structures reversed the optical behaviour. The authors describe the response as different from conventional quasicrystalline materials, "a new type of symmetry-controlled optical behavior that combines quasiperiodic order with chirality", and suggest applications in light and polarization control; the press pieces call these possibilities speculative. They do not mention the spectre. I did not read the paper, and a search also returned an arXiv preprint titled "Chiral Diffraction from Aperiodic Monotile Lattice" (arXiv:2506.07561) which I could not open.

**Connections [repo and mine].**
- [MONOTILES.md](MONOTILES.md) and [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md) already record Hartshorn–Romanov's exact pure-point diffraction for the Chair44 centre set (Bragg masses summing to the density 1/7, Fact 22) and the Spectre's pure-point diffraction as quoted from the literature [repo]. [Mor26] is an *experimental* diffraction result for the hat, a different tile.
- The Chair44 papers propose that physical systems might realise their tilings [repo, MONOTILES.md]. [Mor26] is a physical realisation of a planar aperiodic monotile as an optical pattern. It is a concrete precedent for that kind of proposal, in two dimensions.
- The press releases attribute the chiral diffraction to the pattern lacking mirror symmetry, and report that mirrored patterns reverse the optical behaviour [Mor26]. The hat is the tile that is tiled with reflected copies as well as rotated ones *(background)*; the spectre is a monotile for which reflections are neither needed nor permitted [SPECTRE_CHAIR44.md, citing SM]. I found no experiment on the spectre in the sources here.
- Entry 155 has no physical content [repo, MONOTILES.md], and nothing in [Mor26] relates to translational tilings. The link is only "aperiodic monotile structure".

## Connections at a glance

| Reference | Connects to | Kind of link |
|---|---|---|
| [Gr25] | Entry 155 | Poses the ℤ³ question that entry 155's theorem addresses |
| [GT-JEMS] | Entry 155 | Source of the p-adic Sudoku method that entry 155 follows (cited in the paper) |
| [dGGM-rat] | Entry 155 | Cited by the paper |
| [dGGM26] | Entry 155 | Follow-up to a cited paper; planar (axis-parallel) versus entry 155's cube-union set in ℝ³ [mine] |
| [GK] | Entry 155; entry 266 | Possible route to a connected ℝ⁵ counterexample, hypotheses unchecked [mine]; entry 266's paper cites Tao's Fuglede counterexample for a Hadamard matrix |
| [GGTW] | Entry 076 | The Littlewood paper cites it |
| [Tao-SD] | KAKEYA.md | Thematic (slopes; Kakeya-type), no technical bridge |
| [PS] | TILING.md (packings), entry 090 | Thematic (packing; machine-guided search) |
| [Gyro], [Alv17] | MONOTILES.md | The SCD biprism is combinatorially a gyrobifastigium |
| [IJMS26] | — | Unread |
| [Mor26] | MONOTILES.md, SPECTRE_CHAIR44.md | Experimental diffraction on the hat, versus exact diffraction results on Chair44 and the Spectre |

## Evidence status

| Statement | Status |
|---|---|
| Abstract-level descriptions of [dGGM26], [GK], [GT-JEMS], [GGTW], [PS], [Alv17] | [source], from the arXiv or publisher pages as fetched |
| [Gr25] conjectures and questions | [source], from the HTML version as summarised |
| [Tao-SD] results | Second-hand, from Tao's blog post; the arXiv page was not read |
| [IJMS26] | Title only |
| [Mor26] | Press releases only; the paper and arXiv preprint were not read |
| Comparisons with entries 155, 076, 090 | [mine]; entry 155 and 076 were read from their extracted texts |
| Correctness of any theorem | Not assessed |

## Caveats

- Entry 155 is an unrefereed OpenAI preprint and its Lean tree has not been built here. Statements that it "answers" Question 5.4 depend on that.
- The fetch tool summarises pages. Page and theorem numbers in [Gr25] come from those summaries.
- My suggested comparisons in §2, §3 and §6 are conditional on definitions and hypotheses I did not read.
- The Wikipedia claims in §6 are as fetched; I did not check them against the primary sources.

## Open questions

- Does the cube-union set Ω of entry 155 admit a quasi-periodic or weakly periodic tiling? (§2)
- Do the hypotheses of the [GK] d → d+2 theorem hold for entry 155's ℝ³ counterexample? (§3)
- Does Fuglede's conjecture have a counterexample in dimension three or four? No reference here says; Tao's 2004 counterexample cited by entry 266's paper is in dimension five and higher.
- What does [IJMS26] do with the gyrobifastigium? (§6)

## To do

- [ ] Read [Tao-SD] on arXiv and [Mor26] in full.
- [ ] Read the [GK] and [dGGM26] introductions for exact hypotheses.
- [ ] Add [Gr25] and [dGGM26] to the "Further reading" of [PERIODIC_TILING.md](PERIODIC_TILING.md) if you want them there.
- [ ] Confirm the SCD claim about the gyrobifastigium against a primary source.

## Related pages

- [PERIODIC_TILING.md](PERIODIC_TILING.md) for entry 155's proof.
- [MONOTILES.md](MONOTILES.md) and [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md) for aperiodic monotiles.
- [TILING.md](TILING.md) for dimers, matchings and packings.
- [KAKEYA.md](KAKEYA.md) for the Kakeya and local-smoothing results.
- [COMPLEXITY.md](COMPLEXITY.md) for the complexity side; nothing in the repository ties entry 155 to it.
