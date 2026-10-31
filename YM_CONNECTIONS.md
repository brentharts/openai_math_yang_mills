# Connections to Yang–Mills across the companion repositories

A fork-level note, not part of the OpenAI release. It collects, in one place, every exact connection the fork has found between its companion repositories and Yang–Mills theory, says how strong each connection is, and says where it stops.

**Nothing on this page proves, supports or bears on Yang–Mills existence or a mass gap.** Every item is classical, perturbative, topological or an analogy. The status of the reported claim is in [YANGMILLS.md](YANGMILLS.md).

## The connections, strongest first

| # | Connection | Strength | Checked in | Not checked |
|---|---|---|---|---|
| 1 | The one-loop coefficient of `SU(3)` is the asymptotic-freedom number: `b₃ = −(11 − 2n_f/3)`, which is `−7` for six flavours; the pure-gauge limits are `−11` and `−22/3` | Same equation | `cicy_unification.py`, `ym_asymptotic_freedom.py` | Higher loops beyond the library's two-loop exponent; thresholds |
| 2 | Tree-level gluon amplitudes: the colour-ordered MHV amplitude, its symmetries, its soft limit and the analytic structure BCFW uses, all as exact rationals | Same theory, weak coupling | `cicy_amplitudes.py` | Loops, beyond-MHV amplitudes, the three-point sign conventions |
| 3 | The BPST instanton: a self-dual solution of the classical `SU(2)` equations with second Chern number exactly 1, whose boundary map is the three-sphere of the Hopf fibration | Same object, classical | `hopf_gauge.py` | Anything quantum: the measure, the moduli integral, the vacuum angle |
| 4 | Gauge-coupling unification: the one-loop prediction of `α₃` depends on the field content but, exactly, not on the number of generations | Consistency of a gauge theory with matter | `cicy_unification.py` | Compactifications, couplings, masses |
| 5 | The chirality repository's anomaly sums and the `E₈` branching that contains `SO(10) × SU(4)` | Consistency conditions of a gauge theory | `chirality_lisi.py`, `chirality_frame.py` | Real forms, the Clifford side, the uniqueness theorem (not reproduced) |
| 6 | The two-dimensional `O(N)` models: a gap constructed from a lattice model, and a shared linear drift of the coupling under blocking | Analogy, kept as a correspondence | `explore/ym_mass_gap.py`, [YANGMILLS.md](YANGMILLS.md) | Any transfer from `O(N)` to `SU(N)` |
| 7 | The Hopf insulator: a band gap and an integer invariant of a non-interacting lattice model | None beyond the three-sphere | `hopf_insulator.py` | It has no gauge field |

Rows 1 and 2 are the two places where the library states something about the Yang–Mills Lagrangian itself, and both are exact only in the weak-coupling expansion. Row 3 is exact but classical. Rows 4 and 5 are about whether a gauge theory with given matter is consistent, not about its dynamics. Rows 6 and 7 are listed so that nobody reads them as stronger than they are.

## What the figures show

The concise paper's Figures 5 and 6 are drawn from these models by [`sympy_models/pictures.py`](sympy_models/pictures.py), and every plotted number is asserted before it is drawn.

- **Figure 5a** is row 7: the Hopf invariant `0, −1, +2, −1, 0` against the mass parameter, with the gap closing at `h = ±1, ±3` where the sum jumps by `−1, +3, −3, +1`.
- **Figure 5b** is the section of the fibration that the BPST boundary map lives on, in stereographic projection.
- **Figure 6a** is rows 1 and 4: the one-loop running of the three couplings, which meet near `2 × 10¹⁶` GeV for the MSSM coefficients and do not for the Standard Model's.
- **Figure 6b** is row 3: the instanton's radial density and its cumulative second Chern number, which is exactly `1/2` at `r = ρ` and tends to `1`.

## What each connection does not give

- The limit `c = lim m_L/(a Λ_L)` and its positivity, which is the open step in [YANGMILLS.md](YANGMILLS.md) and the fork's Section 6, appear in none of the rows.
- The instanton number and the Hopf invariant are different integers of different maps. The page [HOPF_GAUGE.md](HOPF_GAUGE.md) says so, and nothing here relates them.
- One-loop running is a perturbative statement. It is evidence that the coupling gets small at short distance, and says nothing about the spectrum at long distance.
- Tree amplitudes are rational functions with poles at massless kinematics. They contain no gap.

## What would be progress

Items the fork could check next, in increasing distance from what is here:

1. **Loop amplitudes.** The CICY repository computes none. A one-loop `N = 4` check of the cusp anomalous dimension would connect rows 1 and 2.
2. **Beyond MHV.** The six-point NMHV amplitudes are checked by CICY only for cyclicity and reflection, since there is no closed formula to compare against.
3. **A lattice-side exact statement for `SU(2)`** beside the `O(3)` result, in the form of a finite-volume transfer-matrix gap that can be computed exactly for a small lattice and compared with the continuum scaling. This is a numerical bridge, not a proof.

## Sources

- [brentharts/CICY](https://github.com/brentharts/CICY) (GPL-3.0), pinned at `67ce063`: README sections "Twistor theory and scattering amplitudes" and "The number 137, and what unification actually predicts".
- [brentharts/HopfFibration](https://github.com/brentharts/HopfFibration) and [brentharts/quantum-lattice](https://github.com/brentharts/quantum-lattice): see [HOPF_GAUGE.md](HOPF_GAUGE.md).
- [brentharts/chirality-standard-model](https://github.com/brentharts/chirality-standard-model): the author's `lisi/` directory only; the repository's theorem core is by another author and is not reproduced.
- Textbook: Parke and Taylor, Phys. Rev. Lett. 56, 2459 (1986); Belavin, Polyakov, Schwartz and Tyupkin, Phys. Lett. B 59, 85 (1975); Georgi, Quinn and Weinberg, Phys. Rev. Lett. 33, 451 (1974).

The fork copies no code from these repositories. Each model rederives the closed forms symbolically or in exact rational arithmetic, and its docstring lists what is not tested.

## Running

```bash
make sympy-pdf      # all models, both PDFs; the pictures need matplotlib and are skipped without it
```
