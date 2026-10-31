# Connections to Yang–Mills across the companion repositories

A fork-level note, not part of the OpenAI release. It collects, in one place, every exact connection the fork has found between its companion repositories and Yang–Mills theory, says how strong each connection is, and says where it stops.

**Nothing on this page proves, supports or bears on Yang–Mills existence or a mass gap in the continuum.** Every item is classical, perturbative, topological, an analogy, (row 1) a rigorous statement about a lattice at strong coupling in three dimensions, or (row 2) an exact computation on a lattice system too small to have a continuum limit. Neither of the last two reaches the continuum limit. The status of the reported claim is in [YANGMILLS.md](YANGMILLS.md).

## The connections, strongest first

| # | Connection | Strength | Checked in | Not checked |
|---|---|---|---|---|
| 1 | Lattice `SU(N)` gauge theory at strong coupling, in `Z³` and `Z⁴`. In `Z³`: in a period-3 forest gauge in which every free link has 7 plaquette partners (8 in temporal gauge, 12 ungauged): Dobrushin uniqueness and a mass bound `m ≥ −log(7β_W/4)` for `SU(2)` at `β_W < 4/7`, and a spectral-gap range `|β| < 1/(2h)` from a certified Hessian bound `h = 9.599991`. In `Z⁴` (period-2 forest): the same window from `h ≤ 14.420001`, which sharpens Bhattacharjee–Bhattacharya–Lockwood's `29/2` and moves the range from `1/29` to `1/28.84` | The same lattice theory, but only at strong coupling and in three dimensions | `forest_gauge.py` (counts, arithmetic, the exact identity-field lower bound `9.3479`); `forest_gauge_4d.py` (the four-dimensional sequel: counts, the identity-field value `9 + √17`, the windows) | The theorems, the 26-million-square interval certificate (3D), the `14.420001` certificate with its Lean and interval checks (4D) and the floating-point lower bound `13.974` are quoted, not rerun; the continuum limit; any statement at weak coupling; the conjecture that the optimal constant is `9.3479` in 3D |
| 2 | Lattice `SU(2)` on systems small enough to solve to full precision. In two dimensions: the string tension `σ = −ln(I₂(β)/I₁(β))`, the ring spectrum `−L ln(I_{2j+1}/I₁)` and Casimir scaling `→ 8/3` (exact, confining, no finite infinite-volume gap). For one plaquette: the Hamiltonian gap, a Mathieu spectrum `(g²/8)[b₄(8/g⁴) − b₂(8/g⁴)]`, with its exact rational strong-coupling series, the weak-coupling limit `2√2 − 5g²/8` and a minimum `2.1962` at `g² ≈ 0.9875`. For two plaquettes sharing a link (a ladder, built from 6j symbols): weak-coupling limit `√6`, exact strong-coupling splitting `26/35` (even) and `14/15` (odd), and a gap below the one-plaquette gap on every grid point, ratio `→ √3/2` at weak coupling | Same lattice theory, but in two dimensions or with one or two plaquettes; computed, not a theorem about a field theory. The spanning-tree reduction is related in method (not in result) to Romanov's forest gauge of row 1, which it does not use | `su2_lattice_gap.py` (34 checks) and `su2_lattice_ladder.py` (17 checks); see [SU2_LATTICE_GAP.md](SU2_LATTICE_GAP.md) | Three or more plaquettes (the strip formula is untested), any scaling window; the 2+1D and 3+1D theories |
| 3 | The one-loop coefficient of `SU(3)` is the asymptotic-freedom number: `b₃ = −(11 − 2n_f/3)`, which is `−7` for six flavours; the pure-gauge limits are `−11` and `−22/3` | Same equation | `cicy_unification.py`, `ym_asymptotic_freedom.py` | Higher loops beyond the library's two-loop exponent; thresholds |
| 4 | Tree-level gluon amplitudes: the colour-ordered MHV amplitude, its symmetries, its soft limit and the analytic structure BCFW uses, all as exact rationals | Same theory, weak coupling | `cicy_amplitudes.py` | Loops, beyond-MHV amplitudes, the three-point sign conventions |
| 5 | The BPST instanton: a self-dual solution of the classical `SU(2)` equations with second Chern number exactly 1, whose boundary map is the three-sphere of the Hopf fibration | Same object, classical | `hopf_gauge.py` | Anything quantum: the measure, the moduli integral, the vacuum angle |
| 6 | Gauge-coupling unification: the one-loop prediction of `α₃` depends on the field content but, exactly, not on the number of generations | Consistency of a gauge theory with matter | `cicy_unification.py` | Compactifications, couplings, masses |
| 7 | The chirality repository's anomaly sums and the `E₈` branching that contains `SO(10) × SU(4)` | Consistency conditions of a gauge theory | `chirality_lisi.py`, `chirality_frame.py` | Real forms, the Clifford side, the uniqueness theorem (not reproduced) |
| 8 | The two-dimensional `O(N)` models: a gap constructed from a lattice model, and a shared linear drift of the coupling under blocking | Analogy, kept as a correspondence | `explore/ym_mass_gap.py`, [YANGMILLS.md](YANGMILLS.md) | Any transfer from `O(N)` to `SU(N)` |
| 9 | The Hopf insulator: a band gap and an integer invariant of a non-interacting lattice model | None beyond the three-sphere | `hopf_insulator.py` | It has no gauge field |

Rows 1 and 2 are the only items about the spectrum of a lattice gauge theory. Row 1 holds only in a regime that excludes the continuum limit, where a gap is expected anyway; the fork rederives its counts and arithmetic and does not rerun its certificate. Row 2 is computed to full precision, but for systems with no continuum limit: the two-dimensional theory confines and has no finite infinite-volume gap, and the one-plaquette gap is positive for the trivial reason that it is a one-dimensional Schrödinger problem. Rows 3 and 4 are the two places where the library states something about the Yang–Mills Lagrangian itself, and both are exact only in the weak-coupling expansion. Row 5 is exact but classical. Rows 6 and 7 are about whether a gauge theory with given matter is consistent, not about its dynamics. Rows 8 and 9 are listed so that nobody reads them as stronger than they are.

## What the figures show

The concise paper's Figures 5 and 6 are drawn from these models by [`sympy_models/pictures.py`](sympy_models/pictures.py), and every plotted number is asserted before it is drawn.

- **Figure 5a** is row 9: the Hopf invariant `0, −1, +2, −1, 0` against the mass parameter, with the gap closing at `h = ±1, ±3` where the sum jumps by `−1, +3, −3, +1`.
- **Figure 5b** is the section of the fibration that the BPST boundary map lives on, in stereographic projection.
- **Figure 6a** is rows 3 and 6: the one-loop running of the three couplings, which meet near `2 × 10¹⁶` GeV for the MSSM coefficients and do not for the Standard Model's.
- **Figure 6b** is row 5: the instanton's radial density and its cumulative second Chern number, which is exactly `1/2` at `r = ρ` and tends to `1`.

## What each connection does not give

- The limit `c = lim m_L/(a Λ_L)` and its positivity, which is the open step in [YANGMILLS.md](YANGMILLS.md) and the fork's Section 6, appear in none of the rows. Row 1 is the closest: it is a proven gap, but at strong coupling `β → 0`, which is the opposite end from the continuum limit `β → ∞`, and in three dimensions.
- In four dimensions the window gains about 0.55% over the earlier proved bound, at most about 3% more is possible with this method, and the repository records that a gap bounded below in lattice units uniformly in the spacing cannot be the continuum step (BBL's no-go theorem), so this route does not lead to the continuum limit. Row 1's three-dimensional ranges are small: `β_W < 4/7` for `SU(2)` by Dobrushin and `β_W < 0.208` by the Hessian route, and since any valid Hessian constant is at least `9.3479` (the exact lower bound checked here), the Hessian route cannot exceed `0.214` for `SU(2)`. The repository's comparison table shows the same ceiling.
- Row 2 has no scaling window. The one-plaquette gap tends to the constant `2√2` at weak coupling, not to zero, so it cannot follow `m ~ e^{−1/(2 b₀ g²)}`; and the two-dimensional string tension has no running coupling, so there is nothing to set beside the O(N) local slopes of `ln m`. What the row does give is a worked, fully checked example of a lattice gauge spectrum at both ends of the coupling, with the expansions failing at the minimum near `g² ≈ 1`.
- The instanton number and the Hopf invariant are different integers of different maps. The page [HOPF_GAUGE.md](HOPF_GAUGE.md) says so, and nothing here relates them.
- One-loop running is a perturbative statement. It is evidence that the coupling gets small at short distance, and says nothing about the spectrum at long distance.
- Tree amplitudes are rational functions with poles at massless kinematics. They contain no gap.

## What would be progress

Items the fork could check next, in increasing distance from what is here:

1. **Loop amplitudes.** The CICY repository computes none. A one-loop `N = 4` check of the cusp anomalous dimension would connect rows 1 and 3.
2. **Beyond MHV.** The six-point NMHV amplitudes are checked by CICY only for cyclicity and reflection, since there is no closed formula to compare against.
3. **A lattice-side exact statement for `SU(2)`** beside the `O(3)` result, in the form of a finite-volume transfer-matrix gap that can be computed exactly for a small lattice and compared with the continuum scaling. This is a numerical bridge, not a proof. *Started, in the two smallest cases (row 2, [SU2_LATTICE_GAP.md](SU2_LATTICE_GAP.md)): the exact two-dimensional string tension and the one-plaquette Hamiltonian gap. The comparison with continuum scaling is not reached, since neither system has a scaling window. The two-plaquette ladder is now done (gap `√6` at weak coupling, volume-dependent, no scaling window). The next step is three plaquettes, to test the predicted strip formula `2√(2 − cos(π/(N+1)))`. The spanning-tree gauge reduction used there is related in method to Romanov's forest gauge (row 1), but none of his results or code are used.*

## Sources

- [romanov360/yang-mills-4d-forest-certificate](https://github.com/romanov360/yang-mills-4d-forest-certificate) (code MIT, text CC BY 4.0), after Bhattacharjee, Bhattacharya and Lockwood, doi [10.5281/zenodo.23729004](https://doi.org/10.5281/zenodo.23729004): README, `forest4d.py`, `bloch4d.py`, `NEXT.md`. The fork rederives the forest, the counts, the identity-field value and the window arithmetic; the certificate, the Lean checks and the interval cover are quoted.
- [romanov360/yang-mills-3d-forest-gauge](https://github.com/romanov360/yang-mills-3d-forest-gauge) (code MIT, paper CC BY 4.0), v1.0.0, doi [10.5281/zenodo.23783317](https://doi.org/10.5281/zenodo.23783317): README, `structure3d.py`, `torus3d.py`, `identity_exact.py`, `bloch_identity.py`. The fork rederives the forest, the partner counts, the threshold arithmetic and the identity-field lower bound itself; the interval-arithmetic certificate and the Lean file are quoted.
- [brentharts/CICY](https://github.com/brentharts/CICY) (GPL-3.0), pinned at `67ce063`: README sections "Twistor theory and scattering amplitudes" and "The number 137, and what unification actually predicts".
- [brentharts/HopfFibration](https://github.com/brentharts/HopfFibration) and [brentharts/quantum-lattice](https://github.com/brentharts/quantum-lattice): see [HOPF_GAUGE.md](HOPF_GAUGE.md).
- [brentharts/chirality-standard-model](https://github.com/brentharts/chirality-standard-model): the author's `lisi/` directory only; the repository's theorem core is by another author and is not reproduced.
- Textbook: Parke and Taylor, Phys. Rev. Lett. 56, 2459 (1986); Belavin, Polyakov, Schwartz and Tyupkin, Phys. Lett. B 59, 85 (1975); Georgi, Quinn and Weinberg, Phys. Rev. Lett. 33, 451 (1974); Kogut and Susskind, Phys. Rev. D 11, 395 (1975); Drouffe and Zuber, Phys. Rep. 102, 1 (1983); NIST DLMF chapter 28 (Mathieu functions), for row 2.

The fork copies no code from these repositories. Each model rederives the closed forms symbolically or in exact rational arithmetic, and its docstring lists what is not tested.

## Running

```bash
make sympy-pdf      # all models, both PDFs; the pictures need matplotlib and are skipped without it
make su2            # only the row 2 models (51 checks, a few seconds)
```
