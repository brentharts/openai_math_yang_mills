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
- The numerical parts of the companion repositories (Berry-curvature grids, preimage tracing, the degree integral, the figures) are not rerun by the fork. Their own test suites cover them.
- Orientation and sign conventions of the Blender scripts and of pyqula (quantum-lattice's `d` uses `conj(z₁)·z₂`, the reflection `y → −y` of the fibration's `h`) are noted, not unified.

## Running

```bash
make sympy      # runs every model, including hopf_fibration, hopf_insulator and hopf_gauge
make hopf       # only the checks of these three models, no export
```

Each check raises on failure, so a model that does not reproduce its numbers stops the export.

## See also

- [YANGMILLS.md](YANGMILLS.md): the status of the reported claim and the companion tooling.
- [sympy_models/README.md](sympy_models/README.md): the model table and the labels.
