# tools/: physics models and the PDF report

Python scripts that model results from the physics-related preprints in this repository, using numpy/scipy for numerics, sympy for symbolic derivations, and matplotlib for figures. Each script also writes LaTeX, and `make` assembles it into a PDF with `pdflatex`.

These are fork-level illustrations. They are not part of the OpenAI release, and none of them re-proves a repository theorem. See [PHYSICS.md](../PHYSICS.md) for the physics-related papers in the repository and [YANGMILLS.md](../YANGMILLS.md) for the Yang-Mills discussion.

## Quick start

From the repository root:

```sh
make setup     # pip install -r tools/requirements.txt
make           # run every model, then build build/physics_models.pdf
make full      # longer Monte Carlo runs and one more ring size (several minutes)
make help      # list all targets
```

Requirements: Python 3, a TeX Live install with `pdflatex`, and the packages in [requirements.txt](requirements.txt). numba is optional but makes the Monte Carlo roughly 100 times faster.

Output goes to `build/` (ignored by git): `build/figures/` (PDF and PNG), `build/tex/` (LaTeX fragments), `build/data/` (raw numbers as JSON) and `build/physics_models.pdf`.

## The models

| Script | Repo entry | What it computes |
|---|---|---|
| [bloch_law.py](bloch_law.py) | 271 | Ideal-magnon density of the 3D Heisenberg ferromagnet as an exact Bessel-function sum. Sympy derives its large-β expansion and checks that the leading coefficient is Bloch's and the next is the lattice correction `3ζ(5/2)/(128π^{3/2})` quoted in the repository's paper. |
| [on_model_mc.py](on_model_mc.py) | 215 | Wolff cluster Monte Carlo for the 2D O(3) and O(4) lattice models: susceptibility and second-moment correlation length against β, compared with the asymptotic-freedom form `β^{1/(N-2)} e^{-2πβ/(N-2)}`. Sympy handles the scaling forms and the coupling drift. |
| [kt_flow.py](kt_flow.py) | 216 | Kosterlitz RG flow for the XY model. Sympy solves it exactly (conserved quantity, closed-form correlation length, `μ ln ξ → π`) and derives the `r^{-1/4} (ln r)^{1/8}` form of the critical correlation function. Scipy checks the closed forms. |
| [haldane_gap.py](haldane_gap.py) | 268 | Exact diagonalization of periodic spin-1 rings, gap against L, and the gap across the bilinear-biquadratic family. Validated against the exact AKLT energy and known spin-1/2 ring energies. |
| [chair44_model.py](chair44_model.py) | external | Pose-level constraint model from the 44-contact atlas of the Chair44 tile (Tsiokos, arXiv:2609.19214, transcribed). Checks that the eight child poses fill the doubled chair, that substitution patches up to 512 chairs have no illegal contact, and runs zero-temperature growth. Not a repository entry; see [../CHAIR44_MODEL.md](../CHAIR44_MODEL.md). |
| [hierarchy_figures.py](hierarchy_figures.py) | external | Exact computations and figures for [../SPECTRE_CLOCK_AND_CHAIR44.md](../SPECTRE_CLOCK_AND_CHAIR44.md): subgroup classes of the cubic rotations, symmetry groups of substitution fixed points, the exact contact-count formula, clock-style recognizability tests and a noisy-patch toy on Chair44. Writes `figures/*.png` (`make figures`), not part of the PDF. |
| [combinatorics_checks.py](combinatorics_checks.py) | 168, 169, 205 | Exact small-case checks for [../COMBINATORICS_REPRESENTATIONS.md](../COMBINATORICS_REPRESENTATIONS.md): equal poset-isomorphism classes of Bruhat intervals in S₃–S₅ carry one Kazhdan–Lusztig polynomial; elementary expansion of chromatic quasisymmetric functions of all natural unit interval graphs on n ≤ 6 vertices (positivity, graph-inversion bookkeeping, palindromicity); universal tensor squares of S_n for n ≤ 22 (exact to n = 15, modulo two primes beyond). Usage: `python3 tools/combinatorics_checks.py [chromatic\|kl\|saxl\|all] [--figures] [--nmax=N]`; `make combinatorics` writes `figures/comb_*.png`, not part of the PDF. Needs networkx (listed in requirements.txt). |
| [zeta_strip_figure.py](zeta_strip_figure.py) | 003 | Critical-strip figure with the 7/8 and 11/12 lines and the displacement range 0 < δ ≤ 3/8 for [../ZETA_ZEROS.md](../ZETA_ZEROS.md); the first 200 zeros of ζ via mpmath compared with the Riemann–von Mangoldt term M(t) used by the Riemann Research Laboratory; functional-equation orbit check. Not a test of the Riemann hypothesis. `make zeta` writes `figures/zeta_*.png`. |

| [sym2lean.py](sym2lean.py), [texmath.py](texmath.py), [lean_emit.py](lean_emit.py) | tooling | sympy → LaTeX → Lean 4 for the equations in [../sympy_models/](../sympy_models/README.md): each equation is printed as LaTeX, read back, compared with its sympy expression, and only then written as a Lean theorem. See [sympy → LaTeX → Lean](#sympy--latex--lean) below. `make sym2lean`, `make sym2lean-pdf`, `make sym2lean-test`; not part of the physics PDF. |

Support files: [common.py](common.py) (paths, plotting style, LaTeX helpers), [make_summary.py](make_summary.py) (headline-results table), and [templates/physics_models.tex](templates/physics_models.tex) (the report skeleton).

## How the checks work

Every script ends with assertions through `common.check`. If one fails, the script exits nonzero and `make` stops, so a PDF only exists if the numerics agreed with the exact results they are tested against (for example the AKLT ground-state energy, the Bloch correction coefficient, and numerical against closed-form correlation lengths). Monte Carlo seeds are fixed, so results are reproducible on the same software versions.

Each section of the PDF states what its model does not test. In particular the Monte Carlo compares a second-moment mass with the shape of an asymptotic formula, not with the transfer-matrix gap, and the Haldane-gap extrapolation uses an assumed finite-size form.

## sympy → LaTeX → Lean

The aim is to make Lean statements about physics readable: a model is written once in sympy, the LaTeX a reader sees and the Lean a prover checks are both generated from it, and the Lean theorem's doc comment carries that LaTeX and its source. The pipeline is adapted from the author's [RosettaMath](https://github.com/brentharts/RosettaMath) (`rosettalean.py`, `leanexport.py`) and [crust](https://github.com/brentharts/crust) (`LEANOS.md`, `tools/latex2ocaml.py`), rewritten to be self-contained here.

```sh
make sym2lean                                   # default models, both flavors -> build/lean/*.lean
make sym2lean SYM2LEAN_ARGS="all --list"        # census of every registered equation, nothing written
make sym2lean SYM2LEAN_ARGS="bloch_law --prove --check"
make sym2lean-pdf                               # build/sym2lean_document.pdf
python3 tools/sym2lean.py --latex 'E = m c^{2}' # one hand-written statement (unverified, see below)
```

| Module | Does |
|---|---|
| [texmath.py](texmath.py) | Reads the physics subset of LaTeX (fractions, powers, roots, subscripted and Greek quantities, elementary functions, one relation) into a small tree, and back into sympy. Refuses, with a reason, anything with a silently wrong reading: `f(x)` (function or product), `\sin 2x`, `x^23`, a bare `e^x`, `o(...)`, sums, integrals, derivatives, `\approx`. |
| [lean_emit.py](lean_emit.py) | Prints a tree as a Lean 4 theorem. Flavor `mathlib` (`import Mathlib`, the flavor that belongs beside [../lean/](../lean/)) or `core` (no imports, `Real` declared as axioms, so a bare `lean` can check that the statements parse and elaborate). Also runs `lean` and reads its messages back per theorem. |
| [sym2lean.py](sym2lean.py) | The driver. Takes `sympy_models/eqlib.REGISTRY` and writes `build/lean/SymModels_{mathlib,core}.lean`, `build/tex/sym2lean.tex` (a fragment) and `build/data/sym2lean.json`. |

Every equation ends in exactly one status: **verified** (the LaTeX reads back to the model's sympy expression, symbolically or, failing that, numerically; a Lean theorem is written), **mismatch** (it does not: a documentation bug, `--strict` exits nonzero), **refused** (outside the subset, with the reason), or **tex-only** (hand-written LaTeX with no sympy relation to compare to; held back unless `--tex-only`).

What a theorem says, and what it does not:

- It is a *statement*. Its body is `sorry` unless `--prove` found a rational-function identity for `ring1` / `field_simp` to close (mathlib flavor only). Most model equations are laws or definitions, not identities, and stay conjectures. Lean accepting a conjecture says it parses, elaborates and is a `Prop`; the report distinguishes `proved` from `statement` from what Lean was never asked.
- Lean's `x / 0 = 0` makes a statement with a denominator easier to satisfy than the physics is, so sympy assumptions (`positive=True`, `nonzero=True`, `integer=True`, ...) become hypotheses and binder types (`ℝ`, `ℤ`, `ℕ`) of the theorem rather than being dropped. Without assumptions a quantity is just a real number.
- `x ^ (y : ℝ)` is Lean's real power, which is the familiar one only for `x > 0`. Integer exponents are printed as natural or integer powers.
- Without a `lean` on the machine `--check` claims nothing. The `mathlib` flavor additionally needs a built Mathlib in `lean/` (`lake exe cache get`).

The fork names sympy symbols with LaTeX source (`Symbol('\\beta_k')`, `Symbol('a\\Lambda_L')`). A name that is one quantity (`\beta_k`) is read as such; a composite name (`a\Lambda_L`, `S-m_S(\beta)`) is one quantity that LaTeX shows as juxtaposed letters, and is read as that one quantity rather than a product. A name that is really an operator (`\lim \frac{d...}{d...}`) is refused.

To use it for a new model, register `Equation(key, title, source_kind, source, expr=sp.Eq(lhs, rhs), ...)` as usual (give symbols their assumptions, and build relations with `evaluate=False` if sympy would otherwise decide them), leave `latex_str` unset so the printed LaTeX is the model, and run `make sym2lean SYM2LEAN_ARGS="<model>"`. Most of the existing registry is hand-written `latex_str` with no sympy expression to check it against, which is why the verified count is low; `make sym2lean-test` runs the 126 selftests without running any model.

### The `@lean` DSL

[leandsl.py](leandsl.py) lets a model be a decorated Python function instead of a registered equation, which reaches beyond the LaTeX-readable subset: natural-number types with casts, `binom`, `fact`, finite sums, user definitions, and hypotheses. Examples are in [dsl_models/](dsl_models/); `make dsl` writes `build/lean/Dsl_mathlib.lean` and `build/tex/dsl.tex`.

```python
@lean.theorem(source="OAI...Clocks.clock_tail")
def clock_tail(N: Nat, r: Nat, B: Real, delta: Real):
    require(B >= 0, delta >= 0, r > 0, 4*E*B*(N + r)*delta <= r)
    return binom(N + r, r) * (2*B*delta)**r <= Rational(1, 2)**r
```

`@lean.define` makes a definition that later code can call. Unlike `sym2lean.py`, Lean is printed straight from the sympy tree with no LaTeX read-back, so the output is "generated", not "verified"; every body is `sorry` unless `proof=` is given, and nothing has been checked by Lean here. Real exponents print as `Real.rpow` (right only for positive bases); floats are refused.

### Python source via `ast`: loop invariants from asserts

[pyast2lean.py](pyast2lean.py) parses ordinary typed Python with `ast` (no sympy in the path) and writes Lean `Id.run do` code, LaTeX pseudo-code and the verification conditions of a loop. `make pyast` runs [py_models/](py_models/) and writes `build/lean/PyAst_mathlib.lean` and `build/tex/pyast.tex`; `PYAST_ARGS=--prove` tries `omega`/`nlinarith` before `sorry`.

```python
def sum_to(n: Nat) -> int:
    assert n >= 0                          # precondition
    s = 0
    i = 0
    assert 2 * s == i * (i - 1) and 0 <= i and i <= n   # invariant (asserts right before the loop)
    decreases(n - i)                       # optional variant
    while i < n:
        s = s + i
        i = i + 1
    assert 2 * s == n * (n - 1)            # postcondition (asserts after the loop)
    return s
```

This yields `sum_to_init`, `sum_to_preserve`, `sum_to_exit` and `sum_to_variant`, with the loop body executed symbolically (parallel assignment and `if` handled exactly). Supported: `int`, `Nat` (as ℤ with `0 ≤ n`, never ℕ), `float` (as exact ℝ), arithmetic, comparisons, `and`/`or`/`not`, `abs min max sqrt exp log factorial comb`, `if`, `while`, `for i in range(...)`. `//` and `%` print as `Int.fdiv`/`Int.fmod`. Everything else (lists, strings, `break`, nested loops in the VC part) is refused with its line number. Nothing here has been checked by Lean; that the `do` block and the VC transition agree is by construction, not proved.

## Running one model

```sh
make bloch       # or: on, kt, haldane, chair44
python3 tools/haldane_gap.py --full
```

Typical run times on a laptop: Bloch and KT about 3 seconds each, Haldane about 6 seconds, the quick Monte Carlo about one minute.

## Adding a model

1. Add `tools/<name>.py` that writes `build/tex/<name>.tex` (a LaTeX section), figures through `common.savefig`, and a JSON file through `common.save_json`.
2. Add the script name to `SCRIPTS` in the Makefile, and `\input{tex/<name>.tex}` to the template.
3. Add a row to [make_summary.py](make_summary.py) if it has a headline number.

## Used by other repositories

- RosettaMath's `rosettaoai.py` reads `build/data/*.json` (set `OPENAI_MATH_DIR` if the checkouts are not side by side) and checks nine library equations against it. Run `make models` first.
- yang-mills-mass-gap-lab imports `build/data/on_model_mc.json` with `ymlab.scaling_comparison.import_fork()`, recording this repository's commit hash with the data.

Changing a script's JSON keys breaks both readers.
