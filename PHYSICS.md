# Physics-Related Results

A fork-level guide to the preprints in this repository that connect to physics: quantum field theory, statistical mechanics, quantum many-body systems, kinetic theory, and general relativity. It is not an official OpenAI document.

**How it was built.** I read the numbered summaries in [CONTENTS.md](CONTENTS.md) for the whole index and skimmed the opening pages of a few papers (the O(3) continuum limit, the isolated particle pole, and 4D scale/conformal symmetry). Most other entries are described from their CONTENTS summaries only, so check the papers before relying on a detail. Links go to each preprint's folder, which holds its README and PDF.

The numbers in brackets are the CONTENTS entry numbers. "Lean" links go to the formalization scope notes where one exists.

## 1. Quantum field theory and mass gaps

The most direct physics content. These are constructions of massive interacting field theories and spectral-gap proofs.

| Entry | Result | Notes |
|---|---|---|
| [215] | [Canonical massive continuum limit of the 2D O(3) model](preprints/The-canonical-massive-continuum-limit-of-the-two-dimensional-O3-model-October-4-2026/) | Constructs a non-Gaussian relativistic theory. The paper's abstract says the fields satisfy the Osterwalder–Schrader axioms, with a unique vacuum and a positive Hamiltonian gap. [Lean](lean/docs/215.md) |
| [215] | [An isolated particle pole for the 2D O(3) spin field](preprints/An-Isolated-Particle-Pole-for-the-Two-Dimensional-O3-Spin-Field-October-4-2026/) | A mass gap alone does not give an isolated particle. This paper proves the lowest mass is an atom of the spectral measure, separated from the rest of the spectrum. |
| [215] | [Exact mass asymptotics for the 2D O(4) lattice model](preprints/Exact-mass-asymptotics-for-the-two-dimensional-O4-lattice-model-October-5-2026/), [sharp mass bounds](preprints/Sharp-mass-bounds-for-the-two-dimensional-O4-model-September-23-2026/) | The lattice transfer gap behaves like a constant times √β·e^(−πβ), the nonperturbative scaling expected from asymptotic freedom. |
| [215] | [Exponential decay in 2D classical O(n) models](preprints/Exponential-decay-in-two-dimensional-classical-On-models-September-23-2026/) | Exponential spin-correlation decay for n ≥ 3 at every positive temperature. |
| none | [Scale and conformal symmetry in 4D operational QFT](preprints/Scale-and-conformal-symmetry-in-four-dimensional-operational-quantum-field-theory-September-26-2026/) | Scale symmetry implies local conformal symmetry for a stated class of 4D unitary, positive-energy theories. The conclusion is local, and the class is defined by hypotheses (discrete scaling spectrum, finite scaling support). It does not construct a specific interacting theory. |
| [280] | [Unitary vertex operator algebras and conformal nets](preprints/Strongly-rational-unitary-vertex-operator-algebras-and-conformal-nets-September-25-2026/) | Every simple unitary strongly rational VOA generates a completely rational conformal net, linking the algebraic and operator-algebraic formulations of 2D conformal field theory. [Lean](lean/docs/280.md) |
| [270] | [Unique threshold bound state of the SU(N) BFSS model](preprints/The-unique-threshold-bound-state-of-the-SU-N-BFSS-model-September-24-2026/), [positive eigenvalues for SU(2)](preprints/Positive-eigenvalues-of-the-relative-SU-2-BFSS-Hamiltonian-October-5-2026/) | The first paper proves the threshold-bound-state conjecture: exactly one normalizable zero-energy state for each N ≥ 2, after removing the center of mass. The continuous spectrum starts at zero, so this state sits at a threshold rather than below a gapped continuum. The second shows infinitely many positive-energy bound states for SU(2), contradicting an exclusion in the original BFSS paper. |

## 2. Statistical mechanics, critical phenomena and conformal limits

| Entry | Result | Notes |
|---|---|---|
| [216] | [Critical and near-critical XY model, BKT universality](preprints/BKT-universality-for-height-and-planar-spin-fields-September-24-2026/) | Critical correlation decay r^(−1/4)(log r)^(1/8) and the BKT essential singularity for the square-lattice XY model. Some conclusions carry stated assumptions. Companions: [critical center magnetization](preprints/Critical-Center-Magnetization-in-the-Planar-XY-Model-October-5-2026/), [critical spin field](preprints/The-Critical-Spin-Field-of-the-Planar-XY-Model-October-5-2026/). |
| [218] | [Conformal universality for weakly interacting Ising models](preprints/Conformal-universality-of-bulk-Ising-correlations-under-weak-interactions-September-23-2026/) | Critical limits and SLE₃ interfaces survive weak perturbations and weak bond disorder. [Lean](lean/docs/218.md) |
| [223] | [Random-cluster interfaces, SLE and CLE limits](preprints/Thermal-FK-Ising-interfaces-and-massive-SLE-October-5-2026/) | SLE_κ limits for critical square-lattice random-cluster interfaces for 0 < q ≤ 4, plus thermal and disordered variants. |
| [211] | [Planar maps and Liouville quantum gravity](preprints/The-critical-Liouville-quantum-sphere-and-geometric-limits-of-FK-maps-at-q-equals-4-September-24-2026/) | Critical FK planar maps converge to Liouville quantum gravity spheres for 0 < q ≤ 4. [Lean](lean/docs/211.md) |
| [225], [232], [233] | Gaussian free field limits for six-vertex and Lipschitz height models; Ashkin–Teller currents | Height-function models converging to Gaussian fields. |
| [213], [214], [224] | Percolation on quasi-transitive graphs and Voronoi percolation | Critical behavior and Cardy's formula for Poisson–Voronoi percolation. |
| [217], [221], [234] | Sherrington–Kirkpatrick, Mézard–Parisi and orthogonally invariant spin glasses | Free energies and fluctuation laws for mean-field spin glasses. |

## 3. Quantum many-body systems and spectral gaps

Gap proofs are the closest technical cousin to a mass-gap statement, though in finite lattice systems rather than continuum field theory.

| Entry | Result | Notes |
|---|---|---|
| [268] | [Spin-one Haldane gap](preprints/The-periodic-spin-one-Haldane-gap-September-24-2026/) | Uniformly positive gap for the antiferromagnetic spin-1 Heisenberg chain on even periodic rings. Companion: [boundary-field gap](preprints/A-boundary-field-gap-for-the-spin-one-Heisenberg-chain-September-24-2026/). |
| [269] | [Uniform Laughlin gap](preprints/Uniform-Stability-of-the-Spherical-Laughlin-Gap-October-5-2026/) | Fermionic Laughlin gap at filling 1/3 on the sphere, stable under weak bounded disorder. Companion: [Fock-space inequality](preprints/A-Fock-space-inequality-and-the-Laughlin-spectral-gap-September-24-2026/). [Lean](lean/docs/269.md) |
| [265] | [2D area law from a spectral gap](preprints/A-two-dimensional-area-law-from-a-global-spectral-gap-September-24-2026/) | Entanglement area law for gapped 2D ground states, with a [polynomial PEPS approximation](preprints/Polynomial-PEPS-approximation-of-gapped-square-grid-ground-states-September-24-2026/). |
| [267] | [Bose–Einstein condensation, dilute hard-sphere gas](preprints/Bose-Einstein-condensation-at-positive-temperature-in-the-dilute-hard-sphere-gas-October-5-2026/) | Positive-temperature condensation and exact quantum depletion. |
| [271] | [Quantum Heisenberg ferromagnet: Bloch's law](preprints/Spontaneous-magnetization-in-the-quantum-Heisenberg-ferromagnet-September-24-2026/) | Bloch's T^(3/2) law with its exact coefficient. |
| [261] | [Anderson localization](preprints/Pure-Point-Spectrum-for-the-Two-Dimensional-Anderson-Model-at-Every-Positive-Disorder-September-23-2026/) | Pure-point spectrum in 2D at every disorder, and absolutely continuous spectrum for weak disorder in d ≥ 3. [Lean](lean/docs/261.md) |
| [263], [275], [278] | Coulomb atoms and molecules: ionization conjectures, QMA-hardness, failure of Kohn–Sham representation | Electronic-structure questions for the full Coulomb model. |
| [262] | [Sharp Lieb–Thirring inequalities](preprints/Sharp-One-Dimensional-Lieb-Thirring-Constants-September-23-2026/) | Sharp constants and equality cases in one dimension. |

## 4. Classical fields, kinetic theory and fluids

| Entry | Result | Notes |
|---|---|---|
| [362] | [Global smoothness for relativistic Vlasov–Maxwell](preprints/Global-classical-solutions-of-the-three-dimensional-relativistic-Vlasov-Maxwell-system-September-23-2026/) | Large-data global existence and uniqueness for the 3D one-species system. The only Maxwell-equation result in the index. [Lean](lean/docs/362.md) |
| [363], [364] | Hard-sphere Boltzmann equation: [nonuniqueness](preprints/Hard-sphere-fluctuations-on-the-regular-Boltzmann-lifespan-September-23-2026/) and kinetic limits | Entropy-solution nonuniqueness, and derivation of Boltzmann from many-particle systems over its lifespan. |
| [376] | [Universal computation in forced Navier–Stokes flows](preprints/Universal-Computation-with-Eventually-Stationary-Navier-Stokes-Forcing-September-27-2026/) | Fluid flows that carry out computation under forcing. |

## 5. General relativity

| Entry | Result | Notes |
|---|---|---|
| [260] | [Spacetime Penrose inequalities](preprints/Boundary-graph-deformations-for-the-spacetime-Penrose-inequality-September-27-2026/) | Sharp enclosing-area inequality in every dimension n ≥ 3, with charge, rotation and AdS extensions. [Lean](lean/docs/260.md) Related: [Kerr–Newman Penrose inequality](preprints/The-Kerr-Newman-Penrose-Inequality-for-Axisymmetric-Electrovacuum-Exteriors-October-5-2026/). |
| [264] | Strong cosmic censorship near Kerr | Local strong cosmic censorship near two-ended Kerr data. |
| [348] | Einstein four-manifolds | Classification of nonnegatively curved Einstein four-manifolds and an L² topological gap. |

## Reasoning summaries

The repository releases abridged summaries of the model's reasoning for ten results (the table is in the main README). Three are physics-related, all in [reasoning_traces](reasoning_traces/):

| Entry | Trace | What it shows |
|---|---|---|
| [221] | [Mézard–Parisi formula for diluted spin glasses](reasoning_traces/mezard-parisi-formula.pdf) | Also available as [LaTeX source](reasoning_traces/mezard-parisi-formula.tex). |
| [271] | [Spontaneous magnetization in the quantum Heisenberg ferromagnet](reasoning_traces/spontaneous-magnetization-quantum-heisenberg-ferromagnet.pdf) | A 6-page summary that begins from the prompt (zero-field KMS states with nonzero magnetization at low temperature for d ≥ 3) and walks through attempted approaches. |
| [362] | [Relativistic Vlasov–Maxwell](reasoning_traces/relativistic-vlasov-maxwell.pdf) | A 6-page summary starting from the full system with its constraint equations (Gauss's law and no magnetic monopoles). |

No reasoning summary exists for the QFT and mass-gap papers in section 1.

## What is not in the repository

A search of the CONTENTS summaries, the Lean documentation, and the preprint titles found no papers on lattice gauge theory, Wilson loops, confinement, or the area law for gauge fields, and none on four-dimensional interacting quantum field theory construction. (This was a keyword search of summaries and titles, not a full-text search of every paper.) The closest material is the two-dimensional sigma-model work in section 1. The tangential mathematical neighbor is entry [043], the P = W conjecture for Higgs moduli spaces, which concerns the geometry of moduli of Higgs bundles rather than quantum field theory.

## Caveats

- Where a summary or abstract mentions hypotheses or stated assumptions (for example in the XY and 4D conformal-symmetry papers), those matter for what is actually proved.
- Only entries with a Lean link have a formalization in this repository, and Lean checks the formal statement, not that it matches the physical conjecture.
- This page groups results by topic. The repository does not claim these results combine into anything beyond what each paper states.

## Related pages

- [YANGMILLS.md](YANGMILLS.md) for the Yang-Mills claim and the subset of these results closest to it.
- [KAKEYA.md](KAKEYA.md) for the harmonic-analysis results, which have uses in PDE and wave propagation.
