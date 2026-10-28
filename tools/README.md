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

Support files: [common.py](common.py) (paths, plotting style, LaTeX helpers), [make_summary.py](make_summary.py) (headline-results table), and [templates/physics_models.tex](templates/physics_models.tex) (the report skeleton).

## How the checks work

Every script ends with assertions through `common.check`. If one fails, the script exits nonzero and `make` stops, so a PDF only exists if the numerics agreed with the exact results they are tested against (for example the AKLT ground-state energy, the Bloch correction coefficient, and numerical against closed-form correlation lengths). Monte Carlo seeds are fixed, so results are reproducible on the same software versions.

Each section of the PDF states what its model does not test. In particular the Monte Carlo compares a second-moment mass with the shape of an asymptotic formula, not with the transfer-matrix gap, and the Haldane-gap extrapolation uses an assumed finite-size form.

## Running one model

```sh
make bloch       # or: on, kt, haldane
python3 tools/haldane_gap.py --full
```

Typical run times on a laptop: Bloch and KT about 3 seconds each, Haldane about 6 seconds, the quick Monte Carlo about one minute.

## Adding a model

1. Add `tools/<name>.py` that writes `build/tex/<name>.tex` (a LaTeX section), figures through `common.savefig`, and a JSON file through `common.save_json`.
2. Add the script name to `SCRIPTS` in the Makefile, and `\input{tex/<name>.tex}` to the template.
3. Add a row to [make_summary.py](make_summary.py) if it has a headline number.
