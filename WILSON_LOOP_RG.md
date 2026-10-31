# Wilson loops through the renormalisation group (companion programme)

A fork-level note, not part of the OpenAI release. It describes a companion repository, [brentharts/4D_Wilson_Loop_RG](https://github.com/brentharts/4D_Wilson_Loop_RG) (a fork of the programme of Martin Peev and Akshunna S. Dogra), says what it claims and how it labels its own status, and says how it relates to the exact lattice results in [SU2_LATTICE_GAP.md](SU2_LATTICE_GAP.md).

**Nothing on this page proves, supports or bears on Yang–Mills existence or a mass gap.** The fork has not checked any of the programme's mathematics. This page reports the programme's own description of itself, as of its commit `70d9a77` (29 October 2026), and one small comparison with the fork's exact two-dimensional result.

## Who wrote what

The programme's README states that the [research outline](https://github.com/brentharts/4D_Wilson_Loop_RG/blob/main/Research_Outline.pdf) and the project announcement were written in full by Peev and Dogra, and that all other documents were generated in full by Claude and/or ChatGPT, including all potentially novel mathematical work. It also says the repository "should [not] be considered an academically rigorous, 'publishable' mathematical paper", that most of it has not been verified directly, and that the authors' judgment is that the manuscript content is "essentially correct". The fork repeats those qualifications and adds none of its own mathematics here. Everything below is attributed to Peev and Dogra and their repository. The programme is licensed CC0-1.0.

## What the programme sets out to do

Bałaban's renormalisation group gives a gauge-fixed lattice measure for Yang–Mills at every finite ultraviolet cutoff. The programme's stated difficulty is that the gauge transformations introduced at successive steps accumulate and seem to obstruct convergence of the gauge-fixed fields as the cutoff is removed. It therefore works with gauge-invariant observables, Wilson loops, which in dimension above two are ultraviolet singular. The aim, at fixed finite physical volume and for a compact connected simple gauge group, is to:

1. identify within the RG the scalar perimeter and corner factors known from the physics literature,
2. divide them out, and
3. prove ultraviolet convergence of the renormalised Wilson-loop correlation functions.

The programme does not pursue infinite volume, and says so, on the ground that it sees no viable route to a mass gap there. It conjectures that if a route to the mass gap exists near the current literature, it lies in a correlation-function inequality for renormalised Wilson loops.

## The first-step statement

For a square contour `C`, the programme's first RG step has the schematic form

```
W_fine(C) ≈ κ_C · W_coarse(C) + Δ_C,        κ_C = exp(−δm · P(C)) · z_c⁴
```

with `P(C)` the perimeter, `δm` the perimeter counterterm and `z_c` the geometric mean of the four corner factors. The bounds it reports for the first step, up to logarithms in the bare coupling `g_ε`, are

```
|δm| ≲ g_ε² / ε,        |log z_c| ≲ g_ε²,        |E[Δ_C]| ≲ g_ε³ / ε
```

with `ε` the fine lattice spacing. Its own assessment is that the remainder bound is not yet sufficient for ultraviolet convergence, since the factor `1/ε` is an artefact of a per-bond estimate, and that eliminating it to reach `g_ε³` is current work. The loop after one step is a "dressed" loop that retains part of the large-field history, and the programme proposes to keep that dressing through the flow and compare with an ordinary loop only at the final infrared scale. It treats two vertex-disjoint square loops at the first step, with cross terms that decay exponentially in the separation.

## The programme's own status table

The README grades each step. This is the programme's assessment, reproduced in compressed form and not checked here.

| Step | The programme's label |
|---|---|
| First RG step, one loop and two disjoint loops | In manuscript (V5, "Complete") |
| Small-field iteration and large-field insertions | Written; verification outstanding |
| R-steps for one and two loops | Conditional on the R-eval, R-round, R-att-ext and R-member estimates |
| Relative weights of large-field contributions | Written; verification outstanding |
| Cutoff-uniform bounds and subsequential limits | Conditional on the R-step induction |
| Comparison of successive cutoffs | Partial / open |
| Full-sequence continuum limit | Not established |
| Arbitrarily many disjoint loops | Not established |
| Decay of connected two-loop correlations | Open |

Their own 27 October status report adds that item-by-item review rounds (by language models and by a CLI referee, recorded in the triage files of the `drafts/` folder) found fatal or serious defects in several drafts, that repairs were made, and that many repairs have not been re-reviewed. A reader should take the table together with that record. The programme also mentions Lean certificates being generated in parallel; the fork has not looked at them.

## How it relates to the fork's lattice results

The two share only the object, a Wilson loop on a lattice with a coupling `g`. They differ in dimension, in method and in what they establish.

| | This fork ([SU2_LATTICE_GAP.md](SU2_LATTICE_GAP.md)) | Wilson-loop RG programme |
|---|---|---|
| Dimension | 2 (Part A); one and two plaquettes (Parts B, C) | 4 |
| Status of results | exact or machine-precision, each rederived by a check | proof drafts, mostly unverified, status graded by the authors |
| Wilson loop | exact: `⟨W(R×T)⟩ = (I₂(β)/I₁(β))^{RT}`, a pure area law | renormalised loop, perimeter and corner factors extracted |
| Perimeter term | none | `δm ∼ g²/ε` per unit length |
| Corner term | none | `log z_c ∼ g²` |
| Finite-volume gap | computed for 1 and 2 plaquettes | not addressed |
| Continuum limit | none, and none claimed | the target |

Two remarks follow.

1. **The two-dimensional loop has no perimeter or corner factor.** In the fork's Part A the lattice SU(2) Wilson loop is exactly `exp(−σ · area)` with `σ = −ln(I₂/I₁)`, so `κ_C` is trivially 1 and the whole of the programme's extraction step is absent. That is a statement about two dimensions: the perimeter and corner renormalisations the programme isolates are a four-dimensional phenomenon, tied to the ultraviolet singularity of the loop. This is a consistency observation, not a test of the programme's estimates.
2. **A concrete place where the two could meet is the strong-coupling end.** The programme works at small `g_ε` (weak coupling). The fork's exact results are strongest at strong coupling and at one or two plaquettes, and have no scaling window. They do not overlap with the programme's regime.

Neither remark supports or tests the programme's claims.

## What is not here

- No part of the programme's mathematics has been verified by the fork, and no result from it is used in the fork's models.
- The fork has not reproduced any of the bounds `g²/ε`, `g²`, `g³/ε`.
- A check of the programme's first-step constants against a small explicit computation would be a natural next addition (see below), but none exists yet.
- The programme's own repository, not this page, is the authority on its status. It is updated often.

## What would be next

- **A numerical check of the first-step shape.** In a small lattice system, compute the ratio of a loop expectation to its leading perimeter-and-corner factor and compare its scaling in `g` with the programme's first-step bounds. Four-dimensional SU(2) at small volume is too large for the exact methods used for the ladder, so this would be a Monte Carlo check with its own statistical error, and it would test only the leading scaling, not the proof.
- **Cross-reading the Bałaban label dictionary.** The programme's [paper-label dictionary](https://github.com/brentharts/4D_Wilson_Loop_RG/blob/main/drafts/BALABAN_PAPER_LABELS.md) lists Bałaban's papers with citations and links. A cross-reference to the fork's existing Bałaban mentions in [YANGMILLS.md](YANGMILLS.md) would be a cheap documentation step.

## Sources

- [brentharts/4D_Wilson_Loop_RG](https://github.com/brentharts/4D_Wilson_Loop_RG), commit `70d9a77`: `README.md` (progress tables), `Research_Outline.pdf` (27 October 2026), `drafts/WILSON_IT1_STATUS_REPORT_2026_10_27.md`. Authors: Martin Peev and Akshunna S. Dogra; material generated by Claude and/or ChatGPT as stated in the repository; licence CC0-1.0.
- This fork: [SU2_LATTICE_GAP.md](SU2_LATTICE_GAP.md), Part A.

## See also

- [YANGMILLS.md](YANGMILLS.md): the status of the reported claim.
- [YM_CONNECTIONS.md](YM_CONNECTIONS.md): every connection to Yang–Mills the fork has checked, ranked by strength.
- [SU2_LATTICE_GAP.md](SU2_LATTICE_GAP.md): the exact lattice results.
