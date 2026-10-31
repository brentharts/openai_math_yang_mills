# Machine-checked and rigorous lattice gauge theory: companion programmes (links only)

A fork-level note, not part of the OpenAI release. It points to two outside efforts to formalize quantum field theory or lattice gauge theory in Lean 4, says what each one says about its own scope, and says how each relates to the fork's own pages. The fork has **read the top-level documents of both and run none of their code**.

**Nothing on this page proves, supports or bears on Yang–Mills existence or a mass gap in the continuum.** Both efforts say the same of themselves, as quoted below.

## What is linked, and how

| Source | Licence | How this fork uses it |
|---|---|---|
| [lluiseriksson/THE-ERIKSSON-PROGRAMME](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME) (checked at commit `0745f2f`, 31 Oct 2026) | AGPL-3.0 | links only. No file is copied into this Apache-2.0 repository. |
| Douglas, Hoback, Mei and Nissim, [*Formalization of QFT*, arXiv:2603.15770v1](https://arxiv.org/abs/2603.15770) (16 Mar 2026) | not stated on the abstract page | a reference. Summarised from the abstract page only. |
| [brentharts/4D_Wilson_Loop_RG](https://github.com/brentharts/4D_Wilson_Loop_RG) | CC0-1.0 | described in [WILSON_LOOP_RG.md](WILSON_LOOP_RG.md); listed here for completeness. |

All three are nodes in the fork's knowledge graph (`sympy_models/knowledgebase.py`, pseudo-entry `formal-gauge`, keys `eriksson2026`, `dhmn-2026`, `wilsonrg2026`), linked to that entry by a "treats" edge. The Eriksson Programme is also cited from the two-dimensional string-tension equation of [SU2_LATTICE_GAP.md](SU2_LATTICE_GAP.md), as "compare, not used".

## The Eriksson Programme

A Lean 4 / Mathlib project on SU(N) lattice Yang–Mills: polymer cluster expansions, Wilson-loop area laws, exponential clustering and source-specific renormalization-group estimates, with an active research branch on one open estimate, called hRpoly.

### What it says about itself

Taken from its README and status pages, as of the commit above. These are the project's statements, not checked by the fork.

- **Recorded core results:** a volume-uniform Wilson-loop area law in an explicit strong-coupling window (theorem `normalized_wilson_loop_area_law`), a finite-volume exact-activity area law (`finite_volume_area_law_exp`), and an exponential clustering bound for truncated plaquette correlations at strong coupling (`gibbs_truncated_correlation_bound`). It reports every headline result as depending on no axiom beyond Lean's three standard ones, with no `sorry`.
- **Open, by its own account:** the physical activity bound `hRpoly` (checklist at 20 of 41 items, no constructed `TermSource`, window 15 "compatible; not attained"), physical reflection positivity and transfer, the four-dimensional continuum limit, and reconstruction.
- **Scope statement:** "No four-dimensional continuum mass gap is claimed", and the project's own historical distance-to-Clay shorthand is "~0% (<0.1%)". Its README says the counts are construction milestones, not a percentage of a Millennium Problem. Its "lattice mass gap" lemma is, in its own words, an elementary lemma about two real sequences and not a spectral gap of any SU(N) operator.
- **Honesty of record:** the repository keeps incident reports and "no-go" notes in the open (the `docs/INCIDENT-*` and `*-NOGO-*` files), and labels results audit-pending where they are.

The fork's reading is only that these are the stated claims. It has not built the Lean project, rerun an oracle check, or read the proofs.

### Where to find things

Pinned to commit [`0745f2f`](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/tree/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd). Every path below existed at that commit.

| Topic | File |
|---|---|
| Start here, headline results | [README.md](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/README.md) |
| The Lean source tree | [YangMills/README.md](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/YangMills/README.md) |
| Volume-uniform and exact-activity area law (Lean) | [RestrictedGate.lean](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/YangMills/L1_GibbsMeasure/RestrictedGate.lean), [WilsonLoopMonomial.lean](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/YangMills/ClayCore/WilsonLoopMonomial.lean) |
| The area-law paper | [paper/area-law/paper.pdf](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/paper/area-law/paper.pdf) |
| What assumptions the core still carries | [HYPOTHESIS_FRONTIER.md](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/HYPOTHESIS_FRONTIER.md) |
| The open estimate hRpoly, current status | [docs/HRPOLY-STATUS.md](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/docs/HRPOLY-STATUS.md) |
| Which claims are verified, and how | [docs/VERIFICATION-LEDGER.md](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/docs/VERIFICATION-LEDGER.md), [REPRODUCIBILITY.md](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/REPRODUCIBILITY.md), [project-state.json](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/project-state.json) |
| The Bałaban-style gauge RG campaign | [docs/BALABAN-RG-PLAN.md](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/docs/BALABAN-RG-PLAN.md), [docs/BALABAN-SOURCE-BOUNDS.md](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/docs/BALABAN-SOURCE-BOUNDS.md) |
| Two-dimensional Bessel and sector-gap note | [papers/bessel-amos-fh/README.md](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/papers/bessel-amos-fh/README.md) |
| Satellite repositories (transfer matrix, 2D Yang–Mills, numerics) | [SATELLITES.md](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/SATELLITES.md) (a dated July 2026 checkpoint, not live) |

### How it relates to the fork

| Fork page | Connection | Strength |
|---|---|---|
| [SU2_LATTICE_GAP.md](SU2_LATTICE_GAP.md), Part A | The same two-dimensional lattice quantities: character coefficients as Bessel functions `I_{2j+1}`, and the ring and sector spectra. The Eriksson note [papers/bessel-amos-fh](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME/blob/0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd/papers/bessel-amos-fh/README.md) states that every character-sector mass gap of two-dimensional Wilson theory (U(1) and SU(2)) is strictly decreasing in `β`. The fork's `σ(β)` series decreases in `β`, the same direction, but the fork does not check or use the monotonicity result. | Same objects; a statement of theirs, not rederived here. **The two use different β conventions** (the note quotes `c_j = I_{2j+1}(2β)/β` for its normalisation, the fork `(2/β) I_{2j+1}(β)`), so numbers must not be compared across them without converting. |
| [YM_CONNECTIONS.md](YM_CONNECTIONS.md), rows 1 and 2 | Both are rigorous or exact lattice statements at strong coupling or in two dimensions, with no continuum limit. The Eriksson area law and clustering are stated for general lattice dimension `d` in an explicit strong-coupling window; the fork's row 1 is Romanov's forest-gauge result in 3 and 4 dimensions, and row 2 is exact in 2D or for one or two plaquettes. | Same regime, different statements and different methods. Neither says anything about the continuum. |
| [WILSON_LOOP_RG.md](WILSON_LOOP_RG.md) | Both build on Bałaban's renormalization group for lattice Yang–Mills. The Peev–Dogra programme works with Wilson loops through the RG; the Eriksson Programme works with cluster expansions and a gauge RG campaign, with hRpoly as its open estimate. No cross-reference between the two was found in the markdown files of the Eriksson repository. | Common source (Bałaban), no shared results. |
| [`lean/`](lean/) | The fork's Lean directory holds Lean formalizations of results in the OpenAI release; it does not use the Eriksson Programme's Lean. | None. |

## Douglas, Hoback, Mei and Nissim: *Formalization of QFT*

[arXiv:2603.15770](https://arxiv.org/abs/2603.15770), submitted 16 March 2026 (hep-th, cs.LO, math-ph, math.MP). From the arXiv abstract page as retrieved: the paper formalizes, in Lean 4, the construction of the free bosonic quantum field theory in four-dimensional Euclidean spacetime and a proof that it satisfies the Glimm–Jaffe axioms (a variant of the Osterwalder–Schrader axioms). It presents this as a proof of concept that extended arguments in mathematical physics can be turned into machine-checked proofs with existing AI tools, describes the design decisions and the use of coding assistants, and discusses how AI-assisted formalization may affect theoretical physics. The page also says the first release assumed three results (Minlos' theorem, the nuclear property of Schwartz space, Goursat's theorem) that later releases proved or avoided.

What is not known to the fork: the content of the paper beyond the abstract, the current state of its Lean library, and its licence.

How it relates:

- **It is about a free field.** It concerns the Gaussian theory, which has no interaction, no gauge group and no mass-gap problem. It is a statement about what can be formalized, not about Yang–Mills.
- **The Eriksson Programme cites it.** The Eriksson repository's `docs/M3-BRIDGE-AUDIT-20260713.md` names it as having formalized "the FREE field in Lean" and uses it to scope what a machine-checked interacting-RG estimate would still be new. That audit is a project-internal document and contains its own risk estimates, which the fork does not reproduce or endorse.
- **Fork parallel.** The fork's [`lean/`](lean/) directory and [YANGMILLS.md](YANGMILLS.md) are the nearest comparisons: a formalization is a check of a statement, not a validation of its physical interpretation, a point all three efforts make.

## What would be next

- Read the Douglas et al. paper in full and record what the Lean library proves and what it assumes.
- Check whether the Eriksson Programme's Lean builds as stated (needs Lean `v4.29.0-rc6` and the pinned Mathlib; the fork has neither in this environment).
- Convert the β convention of the Eriksson 2D note to the fork's and compare the sector gaps numerically, using the fork's `su2_lattice_gap.py` as the independent computation.
- Further references worth adding once read: Romanov's Lean checks (already linked in [YM_CONNECTIONS.md](YM_CONNECTIONS.md)) against the Eriksson Lean for the same strong-coupling statements.

## Sources

- L. Eriksson, THE ERIKSSON PROGRAMME, [github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME](https://github.com/lluiseriksson/THE-ERIKSSON-PROGRAMME), commit `0745f2fe64daeb0bb53da43abcf5bfb4425e7dfd` (31 Oct 2026). AGPL-3.0. Documents read: `README.md`, `YangMills/README.md`, `HYPOTHESIS_FRONTIER.md`, `docs/HRPOLY-STATUS.md` (opening), `docs/README.md`, `SATELLITES.md`, `papers/bessel-amos-fh/README.md`, `docs/M3-BRIDGE-AUDIT-20260713.md` (the passage on arXiv:2603.15770).
- M. R. Douglas, S. Hoback, A. Mei and R. Nissim, *Formalization of QFT*, [arXiv:2603.15770v1](https://arxiv.org/abs/2603.15770), 16 March 2026. Abstract page only.
- M. Peev and A. S. Dogra, [4D_Wilson_Loop_RG](https://github.com/brentharts/4D_Wilson_Loop_RG): see [WILSON_LOOP_RG.md](WILSON_LOOP_RG.md).

## See also

- [YANGMILLS.md](YANGMILLS.md), [YM_CONNECTIONS.md](YM_CONNECTIONS.md), [SU2_LATTICE_GAP.md](SU2_LATTICE_GAP.md), [WILSON_LOOP_RG.md](WILSON_LOOP_RG.md).
