# The Hopf fibration, the Hopf insulator and gauge theory

A fork-level note, not part of the OpenAI release. It records what the fork checks about two companion repositories and one textbook solution, and where the connection to Yang-Mills stops.

**Nothing on this page proves, supports or bears on Yang-Mills existence or a mass gap.** The Hopf insulator's gap is a single-particle band gap of a non-interacting lattice model. The BPST instanton is a classical solution. See [YANGMILLS.md](YANGMILLS.md) for the status of the actual claim.

## Why these belong together

| Object | Where it lives | What it is |
|---|---|---|
| Hopf map `h : S³ → S²` | `S³ = SU(2)` = unit quaternions; fibres are `U(1)` orbits | a principal `U(1)` bundle, Hopf invariant 1 |
| Hopf insulator | `d(k) : T³ → S²` with `d = h(z₁, z₂)` | a lattice band structure classified by the Hopf invariant |
| BPST instanton | classical SU(2) field on `ℝ⁴` | pure gauge at infinity, a map `S³ → SU(2) = S³` of degree 1 |

The only shared fact is the first row: all three use the same three-sphere. The insulator's integer is the Hopf invariant of a map from the Brillouin zone; the instanton's integer is the degree of a map from the boundary of spacetime. They are different integers of different maps, and the fork does not equate them.

## Companion repositories

The fork **links** to these and copies none of their code (`HopfFibration` is GPL-2.0, `quantum-lattice` is GPL-3.0, this fork is Apache-2.0). Everything the fork checks is rederived symbolically in `sympy_models/models/`.

- [brentharts/HopfFibration](https://github.com/brentharts/HopfFibration), pinned at `bd099a6`: numpy and Blender scripts for fibres, stereographic projection and linking, with `hopf_figures.py` and `HOPF_FIGURES.md`.
- [brentharts/quantum-lattice](https://github.com/brentharts/quantum-lattice), pinned at `b5e5736`: PySide6 front end over pyqula, with `HOPF_INSULATOR.md`, `HOPF_INSULATOR_FIGURES.md` and `extras/ql_hopf_figures.py`.

The figure scripts and deep-dive notes (`HOPF_FIGURES.md`, `HOPF_INSULATOR_FIGURES.md`) are at those commits. The model docstrings name them as the comparison point; the fork checks the closed forms, and the companion repositories' own test runs cover the numerics.

## Newer companion work (3-4 Nov 2026)

Both companion repositories have moved on since the pins above: `HopfFibration` at `20c0e7c` and [brentharts/Quantum_Gravity_Lab](https://github.com/brentharts/Quantum_Gravity_Lab) (CC0-1.0) at `4c1474e`. The pins are unchanged. The fork reads the newer scripts from the repositories' own notes and rederives the geometric facts; it copies no code.

| Companion item | What it is | Strength | What the fork does |
|---|---|---|---|
| `HopfFibration/hopf_knots.py` | `f_{p,q}`: preimages are `(q,p)` torus knots, linking number `pq` | direct (the Hopf invariant) | `hopf_fibration.py` already proves the invariant `pq`; the companion's Gauss-integral check is not rerun |
| `HopfFibration/hopf_berger.py`, `hopf_pancake.py`; `Quantum_Gravity_Lab/qgl_pancake.py`, `qgl_hopf.py` | the spatial slice of a biaxial Bianchi IX universe is a Berger sphere (Hopf fibre scaled by `s = R1/R2`); the fibre collapses at the pancake while the slice curvature stays finite | direct for the biaxial case only | `su2_berger_top.py` rederives the volume `2 pi^2 R1 R2^2`, Milnor's Ricci curvatures and the full Laplace spectrum ([SU2_BERGER_TOP.md](SU2_BERGER_TOP.md)); the companions' numerical pancake solution is not rerun |
| `HopfFibration/hopf_be_window.py` | one-link Bakry-Emery tensor of lattice `SU(2)`; the `1/3` threshold and its Berger-metric window | indirect (shared `S^3 = SU(2)`) | not rerun; the caveat on Berger metrics below applies |
| `Quantum_Gravity_Lab/tests/test_s3_identifications.py` | Euler-angle `S^3`: fibre period `4 pi`, fibre length `2 pi R1`, base radius `R2/2` | geometry of the same `S^3` | the same facts are in `su2_berger_top.py` (volume, Laplacian in Euler angles) and need the period `4 pi` for half-integer `j` |

The companions report numerical cross-checks of this fork's Hopf statements (23 of 23) and of the geometry behind the pancake modules (18 of 18). Those were scratch scripts and are in neither repository's suite; the fork has not rerun them.

## What is checked

Labels follow the companion notes: **Exact** is symbolic or exhaustive arithmetic done here, **Counted** is an exact integral or finite count, **Quoted** is taken from a paper and not rederived, **Standard** is textbook material.

### `models/hopf_fibration.py` (entry `hopf-fibration`)

| Result | Label |
|---|---|
| `\|h\| = \|z\|²`, and `h` is constant on `e^{it}(z₁, z₂)` | Exact |
| `\|⟨z, w⟩\|² = ½(\|z\|²\|w\|² + h(z)·h(w))`, so fibres are Clifford parallel and their distance is half the base angle | Exact |
| `h*σ = −2 dA` in Hopf coordinates, and `(1/4π²)∫ A∧dA = ±1` | Exact / Counted |
| `f_{p,q} = Hopf(z₁ᵖ, z_2^q)` normalised has Hopf invariant `pq` (integrand `pq S′(η)`) | Exact |
| half-plane section: base point, closed form `ρ = s/(1 + c sin a)`, `z = c cos a/(1 + c sin a)`, circles of centre `(1/s, 0)` and radius `c/s` | Exact |

The sign in `h*σ = −2 dA` is the orientation convention of this `h`; the invariant is `±1` and its modulus is what is used.

### `models/hopf_insulator.py` (entry `hopf-insulator`)

| Result | Label |
|---|---|
| gap closes only at `h ∈ {−3, −1, 1, 3}` (`z₁ = z₂ = 0` at a TRIM with `m = 0`) | Exact |
| invariant `Σ_TRIM χ [m < 0]` equals `0, −1, +2, −1, 0` in the five regimes, with TRIM multiplicities `1, 3, 3, 1` and jumps `−1, +3, −3, +1` | Exact |
| `\|z₂\|² = 1 + (h−2)² + 2(h−2)cos k_z` on `k_x = k_y = π` | Exact |
| the invariant equals the linking number of preimages; the model's invariant values | Quoted (Moore–Ran–Wen 2008; Deng–Wang–Shen–Duan 2013) |

The TRIM formula was found for the companion figures and compared there with a numerical degree integral at 13 values of `h`; the fork evaluates it exactly but does not redo the degree integral.

### `models/hopf_gauge.py` (entry `hopf-fibration`)

| Result | Label |
|---|---|
| BPST field `A = 2η x/(x² + ρ²)`: `F = −4ηρ²/(x² + ρ²)²`, self-dual, `D_μF_{μν} = 0` | Standard (checked symbolically) |
| `F F̃ = 192ρ⁴/(x² + ρ²)⁴` pointwise and `c₂ = (1/32π²)∫ F F̃ d⁴x = 1` for every `ρ` | Standard / Counted |
| the boundary map `x/\|x\|` is the identity of `S³` (volume density `sin²θ₁ sinθ₂`, total `2π²`), degree 1 | Exact |

This item is textbook material (Belavin–Polyakov–Schwartz–Tyupkin 1975), included so the page can state in one place how the integer is obtained and what it is. It is not a new result.

## What this does not say

- It does not construct, or move toward constructing, a quantum Yang-Mills theory on `ℝ⁴`, and it makes no statement about a mass gap. The instanton number is a topological integer of a classical field.
- The Hopf insulator is non-interacting. Its gap is a band gap, and its invariant says nothing about gauge fields.
- The Hopf invariant and the instanton number are not related by anything on this page.
- A Berger metric (Hopf fibre scaled) is left-invariant but not bi-invariant. A lattice link carries two `SU(2)` actions, the gauge transformations at its two ends, and only the round metric respects both. So results for squashed links (`su2_berger_top.py`, and the Berger-metric window of `hopf_be_window.py`) are geometry of `S^3`, not statements about a gauge-invariant theory. See [SU2_BERGER_TOP.md](SU2_BERGER_TOP.md).
- The numerical parts of the companion repositories (Berry-curvature grids, preimage tracing, the degree integral, the figures) are not rerun by the fork. Their own test suites cover them.
- Orientation and sign conventions of the Blender scripts and of pyqula (quantum-lattice's `d` uses `conj(z₁)·z₂`, the reflection `y → −y` of the fibration's `h`) are noted, not unified.

## Running

```bash
make sympy      # runs every model, including hopf_fibration, hopf_insulator and hopf_gauge
make hopf       # only the checks of these three models, no export
make berger     # only the Berger-sphere model (SU2_BERGER_TOP.md)
```

Each check raises on failure, so a model that does not reproduce its numbers stops the export.

## See also

- [YANGMILLS.md](YANGMILLS.md): the status of the reported claim and the companion tooling.
- [SU2_BERGER_TOP.md](SU2_BERGER_TOP.md): the Laplace spectrum of `SU(2)` with a squashed Hopf fibre.
- [YM_CONNECTIONS.md](YM_CONNECTIONS.md): all the connections to Yang–Mills the fork has checked, ranked by strength.
- [sympy_models/README.md](sympy_models/README.md): the model table and the labels.
