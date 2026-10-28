# Yang-Mills and QCD

This page tracks the Yang-Mills / QCD claim reported in October 2026. It is a fork-level note, not an official OpenAI document.

## Status in this repository

**No Yang-Mills manuscript, Lean file or reasoning trace is included in this repository.** A search of the markdown files finds no entry for it in [CONTENTS.md](CONTENTS.md), and [history.md](history.md) does not list it among the releases. Everything below comes from secondary reporting and should be checked against OpenAI's own announcement.

## What has been reported

Source: [CellCog summary, Oct 26–27, 2026](https://cellcog.ai/blog/openai-yang-mills-eth-kakeya/). I could not retrieve the primary OpenAI page, "Building scientific understanding with academic communities", so treat these points as unverified.

- OpenAI is reported to claim an **existence proof for quantum chromodynamics (QCD) in four dimensions**, including pure Yang-Mills theory. The pure Yang-Mills part is described as a subproblem of the Clay Millennium Prize problem.
- Reported timeline: first AI output on Sep 9 for Yang-Mills and Sep 18 for QCD, with the Yang-Mills work finished on Oct 16.
- The artifacts are reportedly **held by an expert group** rather than released. The group is reported to be led by Nobel laureate David Gross, with Arthur Jaffe (co-author of the Clay problem statement), Douglas, Klebanov and Nekrasov (surnames only, as given in the source). It is reported to plan its own write-up before any AI materials are published.

## What the claim does and does not cover

- The Clay Millennium problem asks for both **existence** of a quantum Yang-Mills theory on ℝ⁴ and a **mass gap**. The reporting describes an existence claim only and does not mention a mass gap result. See the [Clay problem statement](https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf) for the exact requirements.
- Without the manuscript or a formalization, nothing here can be checked. In particular, there is no Lean file whose statement could be compared with the Clay formulation.

## Related results already in this repository

None of the papers below is about Yang-Mills, and the repository does not link them to it. They are grouped here because they prove the kinds of statements a rigorous Yang-Mills construction would need: constructing a field theory from a lattice model, proving a spectral gap, and handling a gauge-type matrix model. See [PHYSICS.md](PHYSICS.md) for the wider list.

- **Mass gap from a lattice model, in two dimensions.** [The canonical massive continuum limit of the 2D O(3) model](preprints/The-canonical-massive-continuum-limit-of-the-two-dimensional-O3-model-October-4-2026/) builds an interacting relativistic theory from a lattice model. Per its abstract, the fields satisfy the Osterwalder–Schrader axioms, with a unique vacuum and a positive Hamiltonian gap ([Lean](lean/docs/215.md)). The 2D sigma model is a standard toy for Yang-Mills because it is asymptotically free and has a dynamically generated mass. The step from there to 4D gauge theory is the hard part. Details from the paper's introduction:
  - It sits in the block renormalization-group tradition and cites Bałaban's work on lattice gauge theory among its methodological ancestors. That is the same tradition a rigorous lattice Yang-Mills construction would draw on.
  - The proof works by exact blocking: partition the lattice into squares, keep one averaged spin per square, and integrate out the rest, step after step, with the bare coupling growing linearly in the number of steps (the paper's equation 1.4 has β_N ≈ H + (log L / 2π) N).
  - It does not transfer the earlier O(4) mass theorem to O(3). The paper says the different target dimension changes the current identities, the Gaussian normalization and the coefficient of the coupling drift, and it re-establishes those directly.
  - It proves a nonzero connected four-point function, so the limiting theory is interacting rather than free. It does not use the integrable factorized S-matrix as a premise.
  - The result concerns one lattice action (nearest-neighbor, no topological term, no external field).
- **Isolated particle mass.** [An isolated particle pole for the 2D O(3) spin field](preprints/An-Isolated-Particle-Pole-for-the-Two-Dimensional-O3-Spin-Field-October-4-2026/) stresses that a mass gap and an isolated lowest mass are different statements, and proves the second.
- **Asymptotic-freedom scaling.** [Exact mass asymptotics for the 2D O(4) lattice model](preprints/Exact-mass-asymptotics-for-the-two-dimensional-O4-lattice-model-October-5-2026/) gives a transfer gap that scales like √β·e^(−πβ).
- **Four-dimensional QFT structure.** [Scale and conformal symmetry in 4D operational QFT](preprints/Scale-and-conformal-symmetry-in-four-dimensional-operational-quantum-field-theory-September-26-2026/) works in four dimensions, but under hypotheses that define a class of theories (a discrete, bounded-below scaling spectrum with finite multiplicities, finite scaling support, and a bounded local net). It does not construct an interacting gauge theory. Two points from its introduction bear on gauge theories: it states that in a gauge description the construction takes place in the chosen gauge-invariant local algebras, and its conclusion is a local conformal Ward action, with integration to a global conformal action left as a separate question. It also notes the framework needs no Haag duality or split property.
- **Gauge-group matrix model.** [The unique threshold bound state of the SU(N) BFSS model](preprints/The-unique-threshold-bound-state-of-the-SU-N-BFSS-model-September-24-2026/) and [positive eigenvalues for SU(2)](preprints/Positive-eigenvalues-of-the-relative-SU-2-BFSS-Hamiltonian-October-5-2026/) concern SU(N) matrix quantum mechanics. This is quantum mechanics, not a field theory, and it is a useful contrast with a mass gap: the paper notes the continuous spectrum starts at zero, so the zero-energy bound state lies at a threshold rather than below a gapped continuum.
- **Spectral gaps in finite systems.** The [Haldane gap](preprints/The-periodic-spin-one-Haldane-gap-September-24-2026/) and [Laughlin gap](preprints/Uniform-Stability-of-the-Spherical-Laughlin-Gap-October-5-2026/) papers prove uniform gaps as system size grows, a technical cousin of a mass gap.
- **Maxwell fields.** [Global smoothness for relativistic Vlasov–Maxwell](preprints/Global-classical-solutions-of-the-three-dimensional-relativistic-Vlasov-Maxwell-system-September-23-2026/) ([Lean](lean/docs/362.md)) is the classical, abelian analogue: a 3D system coupling matter to electromagnetic fields.

## What the repository does not contain

A keyword search of the CONTENTS summaries, the Lean documentation and the preprint titles found no papers on lattice gauge theory, Wilson loops, confinement or the area law, and none constructing an interacting four-dimensional quantum field theory. This was a search of summaries and titles, not of every paper's full text. The gap between the two-dimensional sigma-model results above and 4D Yang-Mills is therefore the open part, and the repository gives no indication of how it would be bridged.

## Why this page exists

This repository publishes the model's math manuscripts and their Lean formalizations. The Yang-Mills material is deliberately not part of that release, so this page records the gap and gives one place to update when something is published.

## Related pages

- [KAKEYA.md](KAKEYA.md) for the local smoothing and Kakeya results announced alongside this claim.
- [history.md](history.md) for the release log.

## To do

- [ ] Add the primary OpenAI announcement link once confirmed.
- [ ] Replace this page's secondary sourcing with the expert group's write-up when it appears.
- [ ] If a manuscript or Lean formalization is released, link it and compare the formal statement with the Clay problem.
