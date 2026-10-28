# Kakeya and Local Smoothing

This page collects the Kakeya-related material in this repository and the October 2026 news around it. It is a fork-level guide, not an official OpenAI document. For the authoritative statements, use the linked papers and Lean files.

## The problem in brief

A **Kakeya set** in ℝⁿ is a set containing a unit line segment in every direction. The **Kakeya conjecture** says every such set has Hausdorff dimension n. The stronger **Kakeya maximal conjecture** bounds the maximal average over δ-tubes with an L^n estimate that loses only δ^(-ε) for every ε > 0, and it implies the dimension statement.

Kakeya sits inside a web of harmonic-analysis conjectures (local smoothing, Bochner–Riesz, Fourier restriction), which is why one result can settle several at once.

## Results in this repository

| Date | Result | Paper | Lean |
|---|---|---|---|
| Sep 23, 2026 | Kakeya maximal conjecture in ℝ³ | [paper](preprints/The-Kakeya-maximal-conjecture-in-three-dimensions-September-23-2026/paper.pdf) | see CONTENTS entry 074 |
| Sep 24, 2026 | Every Kakeya set in ℝ⁴ has Hausdorff dimension 4 | [paper](preprints/Every-four-dimensional-Kakeya-set-has-full-Hausdorff-dimension-September-24-2026/paper.pdf) | see CONTENTS entry 074 |
| Sep 24, 2026 | Critical L³ local smoothing for the 3D wave equation | [paper](preprints/Critical-local-smoothing-for-the-three-dimensional-wave-equation-September-24-2026/paper.pdf) | [079](lean/docs/079.md) |
| Oct 25, 2026 | Local smoothing for the Euclidean wave equation, all dimensions n ≥ 2 | [paper](preprints/Local-smoothing-for-the-Euclidean-wave-equation-October-25-2026/paper.pdf) | [079](lean/docs/079.md) |

Index entries: **074** (Kakeya in three and four dimensions) and **079** (Local smoothing for the Euclidean wave equation) in [CONTENTS.md](CONTENTS.md). The release is described in [history.md](history.md) under October 26, 2026.

## What the October 25 paper claims

According to this repo's own summary (CONTENTS entry 079 and the Lean scope note):

- Sogge's local smoothing conjecture for the Euclidean wave equation in every spatial dimension n ≥ 2, with arbitrarily small loss above the sharp Sobolev threshold. The proof uses an auxiliary partial differential equation.
- Consequences: the Bochner–Riesz conjecture, maximal bounds for these means for 2 ≤ p < ∞, and sharp bounds for Stein's square function.
- Fourier restriction for the sphere, compact paraboloid and truncated cone **throughout the interior of the conjectured range** (not the full range, including endpoints), plus Schrödinger local smoothing.
- The Kakeya and Nikodym maximal conjectures in all dimensions. In particular, every Euclidean Kakeya set has full Hausdorff dimension.

The [Lean formalization scope note](lean/docs/079.md) lists the formal statements, and [Comparator links](lean/ComparatorChallenges/LocalSmoothing.lean) point to each one: half-wave local smoothing, Bochner–Riesz bounds, spherical extension, the Kakeya maximal estimate, full Hausdorff dimension of Kakeya sets, and Stein's square function.

## Recent news

Source: [CellCog summary, Oct 26–27, 2026](https://cellcog.ai/blog/openai-yang-mills-eth-kakeya/). This is a secondary source. The primary OpenAI announcement, "Building scientific understanding with academic communities", could not be retrieved when this page was written, so please check these items against it.

- OpenAI's announcement of Oct 26, 2026 reports local smoothing in all dimensions, with Kakeya following, alongside an ETH disproof and a Yang-Mills/QCD existence claim.
- Reported timeline for local smoothing/Kakeya: first AI output Oct 18, Lean check Oct 19, then public release of the paper and Lean proof.
- One unnamed outside expert reportedly advised on open release.
- Background cited there: Sogge proposed local smoothing in 1991, Guth–Wang–Zhang proved the 2+1 dimensional case in 2019, and Wang–Zahl proved Kakeya in three dimensions in 2025.

Note one difference in framing. CellCog's phrasing suggests the Stein restriction conjectures are resolved in all dimensions, while this repository's summaries say restriction is obtained in the interior of the conjectured range. Prefer the repository's wording.

## Caveats

- The manuscripts are released as model output, which OpenAI has described as raw research data rather than finished, refereed science (per CellCog's report).
- Lean verifies the proof of the formal statement. It does not by itself show the formal statement matches the original conjecture, which is why the Comparator files exist and why expert review still matters.
- The 3D result of Sep 23 and the Oct 25 paper overlap in scope. Check each paper's introduction for how they relate.

## To do

- [ ] Add the primary OpenAI announcement link once confirmed.
- [ ] Confirm whether Lean docs exist for entry 074 and link them.
- [ ] Add a short reading guide to the auxiliary-PDE argument.
