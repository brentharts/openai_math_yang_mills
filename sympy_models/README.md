# sympy_models/

Symbolic models of results from the preprints and other references in this repository, written with sympy. Each model builds its equations as sympy expressions, checks them symbolically, and exports them as LaTeX so one document can later pull all the equations together.

This is fork-level tooling, not part of the OpenAI release. The checks test that formulas are consistent with each other. They do not re-prove any paper.

## Quick start

```sh
pip install -r sympy_models/requirements.txt
python3 sympy_models/export.py              # run every model, write sympy_models/build/
python3 sympy_models/export.py --pdf        # also build build/equations_document.pdf (needs pdflatex)
make sympy-pdf                              # same, from the repo root
python3 sympy_models/export.py bloch_law kt_flow   # run selected models only
```

Output in `sympy_models/build/` (regenerate it, do not edit it):

| File | Contents |
|---|---|
| `equations.tex` | One numbered `equation` per registered equation, labelled `eq:<key>`, with its source as a LaTeX comment. Use `\input{equations.tex}`. |
| `equations.json` | Key, title, LaTeX, source kind, source, caveats, model and the checks each equation passed. |
| `equations.md` | A short index table with CONTENTS entry numbers. |
| `equations_summary.tex`, `equations_by_entry.tex` | Generated fragments for the document below: a summary table and one section per CONTENTS entry. |
| `equations_document.pdf` | The built document (only with `--pdf`). |

## The equations document

[equations_document.tex](equations_document.tex) is the hand-editable wrapper. It holds the preamble, abstract and provenance legend, and `\input`s the two generated fragments. Build it with `python3 sympy_models/export.py --pdf` or `make sympy-pdf`. The PDF has one section per CONTENTS entry (titles are read from [../CONTENTS.md](../CONTENTS.md)), a last section for equations with no entry, and for each equation its source kind, source, caveats and the checks it passed. Each equation keeps its label `eq:<key>`, so it can be cited with `\ref` or `\eqref`.

- Set `\checksfalse` in the wrapper for a compact version without the lists of checks.
- Running a subset of models (`export.py bloch_law --pdf`) builds a partial document.
- The wrapper uses the default OT1 font encoding on purpose. T1 falls back to bitmap fonts unless cm-super or lmodern is installed, and neither is in the standard minimal TeX Live used here.

Long statements use an optional `%%BR%%` marker in `latex_str` to mark where they may break across lines.

The export exits nonzero if any check fails, so a build never contains an equation whose check failed.

## Layout

- [eqlib.py](eqlib.py): the `Equation` record, the registry, `check()` and the LaTeX block writer.
- [export.py](export.py): discovers modules in `models/`, runs them and writes `build/`.
- [models/](models/): one module per model.

## Adding a model

1. Create `models/<name>.py` with `NAME = "<name>"` and a `build_and_verify()` function.
2. Build expressions with sympy. For each equation worth keeping, call `register(NAME, Equation(...))`.
3. Verify with `check(NAME, "what was checked", condition, eq=the_equation)`. Prefer `sp.simplify(a - b) == 0`, limits and series over numeric comparisons.
4. Run `python3 sympy_models/export.py <name>`.

## Provenance

Every equation records where it came from, because the sources carry different weight:

| `source_kind` | Meaning |
|---|---|
| `paper` | Read from a manuscript in `preprints/`. Give the entry, section and equation number. |
| `repo-doc` | Quoted from this fork's markdown notes, not re-read from the paper. Check the paper before citing. |
| `textbook` | Standard background that is not in any manuscript here. |
| `derived` | Obtained by sympy from other registered equations. |

## Models

Equations from the manuscripts below were transcribed from the PDFs, with theorem and equation numbers recorded in each `source`. The checks verify the consistency of the stated formulas. They do not re-prove any theorem, and each module's docstring says what it leaves untested.

| Module | Entry | Manuscripts used | Verified here | Not tested |
|---|---|---|---|---|
| [bloch_law.py](models/bloch_law.py) | 271 | Bloch's Law paper; first lattice correction; spherical magnetization law | D and the dispersion's quadratic part (anisotropic and next-neighbour examples); the Bloch coefficient and the lattice term 3ζ(5/2)/(128π^{3/2}) derived independently from the ideal-magnon integral via the large-argument series of e^{-z}I₀(z); that series numerically against the exact Bessel sum at βS = 30 and 60; sinh(m\|t\|)/(m\|t\|) as the sphere average and the m² second moment | The interacting spin-system limit itself; the Dyson interaction terms (order T⁴) |
| [kt_flow.py](models/kt_flow.py) | 216 | Critical logarithmic correction; essential singularity; BKT height/spin fields; critical spin field; center magnetization | The paper's equations (1.4)–(1.8) against each other: fixed profile (1,1)/j, eigenvector (1,−1/2) with eigenvalue 1+1/j, flow solution k·cot and k·csc, its small-w expansion, the pole T = π/k and A_XY = (log L)π/√(3𝒫); v² = 1/4 and 8πv² = 2π; the exponent 1/(2πK) and K = 2/π giving 1/4 | The constants A_XY, B_XY and 𝒫 (not evaluated in the papers); the center-magnetization result is conditional on companion inputs |
| [haldane_gap.py](models/haldane_gap.py) | 268 | The periodic spin-one Haldane gap | Spin-one commutators and Casimir; f, G and f(1−(1−u)³) = G(u)u²; the purity identities (4.5)–(4.6); the squaring lemma on random exact weights; both parameter rows of (4.18) with hypotheses (4.10) and the resulting bounds (4/105)·log(80/79) and log(20)/784; exact diagonalization for L = 4…10 as a sanity check of H_L | The finite inputs S(60, 105/4) > 7/8 and the like (Propositions 5.4–5.5, Appendices A–C), which are the heavy part of the proof |
| [periodic_tiling.py](models/periodic_tiling.py) | 155 | A translational tile with no fully periodic tiling in dimension three | Section 2 at p = 3 and p = 5 against an independent enumeration of all allowed words: the last-nonzero-digit model satisfies the line rule, has nonconstant columns and no vertical period; the exceptional-index count and the seed-test remark | The tiling construction (Sections 3–7); the paper fixes p > 200 |
| [on_scaling.py](models/on_scaling.py) | 215 | Fork notes only (`repo-doc`) | N = 4 reduction, log-slope, blocking drift | Everything in the papers |
| [ym_asymptotic_freedom.py](models/ym_asymptotic_freedom.py) | none | Textbook | Two-loop lattice scale against the beta function | No manuscript here treats a gauge theory (see [../YANGMILLS.md](../YANGMILLS.md)) |

Notes from reading the sources:
- In the essential-singularity paper, δ = √(b_c − b). The radical is detached from its argument in the PDF text layer, so this was confirmed on the rendered page. Only with this reading do (1.3), (1.6) and (1.8) agree.
- The Bloch coefficient and the 3ζ(5/2)/(128π^{3/2}) correction are in two separate papers, not in the spherical-magnetization paper that CONTENTS lists first.
- A check that passes at p = 3, 5 for the tiling model is not evidence about p > 200. A mutation test (replacing the last nonzero digit by `t mod p`) makes the line-rule check fail on hundreds of words, so the check is not vacuous.

## Caveats

- Display forms are kept symbolic (for example `b_0`, `b_1` instead of their numeric values) so the LaTeX stays readable. Checks run on the substituted form, and a check ties the two together.
- sympy chooses the LaTeX layout. Rearrange in the document if you want a different presentation.
- A `repo-doc` equation is only as reliable as the note it came from.
