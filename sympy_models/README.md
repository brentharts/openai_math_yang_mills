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
| `equations_summary.tex`, `equations_by_entry.tex`, `knowledge_graph.tex`, `bibliography.tex` | Generated fragments for the document below: a summary table, one section per CONTENTS entry, the knowledge-graph part with its appendices, and the references. |
| `equations_document.pdf` | The built document (only with `--pdf`). |

## The equations document

[equations_document.tex](equations_document.tex) is the hand-editable wrapper. It holds the preamble, abstract and provenance legend, and `\input`s the two generated fragments. Build it with `python3 sympy_models/export.py --pdf` or `make sympy-pdf`. The PDF has one section per CONTENTS entry (titles are read from [../CONTENTS.md](../CONTENTS.md)), a last section for equations with no entry, and for each equation its source kind, source, caveats and the checks it passed. Each equation keeps its label `eq:<key>`, so it can be cited with `\ref` or `\eqref`.

- Set `\checksfalse` in the wrapper for a compact version without the lists of checks.
- Running a subset of models (`export.py bloch_law --pdf`) builds a partial document.
- The wrapper uses the default OT1 font encoding on purpose. T1 falls back to bitmap fonts unless cm-super or lmodern is installed, and neither is in the standard minimal TeX Live used here.

Long statements use an optional `%%BR%%` marker in `latex_str` to mark where they may break across lines.

The export exits nonzero if any check fails, so a build never contains an equation whose check failed.

## The knowledge graph

[knowledgebase.py](knowledgebase.py) is a knowledge graph of equations in the manner of [RosettaMath](https://github.com/brentharts/RosettaMath) (`rosettaphys.py`, `rosettapaper.py`), with sympy doing the algebra. The data at the bottom of the file, 155 equations, 150 physical quantities, 125 per-equation "readings" tables and 92 hand-stated links, is copied from RosettaMath's knowledge base (MIT licence, same author); the engine above it is new.

- An equation declares which quantity each of its letters denotes (its *reading*). Two equations are joined only on a quantity both declare, and only if every other letter they share is read the same way. A glyph is not a quantity: the `W` of S = k_B log W counts microstates, the `W` of W = F d is an energy.
- The reader turns each equation's LaTeX into sympy (94 of the 155 are accepted; derivatives, sums, traces, primes and equations without readings stay in the graph but are not rearranged, and `refusals()` says why). `solve_for`, `join`, `all_joins` and `families` then work on sympy expressions; an even power is refused because sympy returns more than one root.
- `check_dimensions` tests each readable equation for dimensional homogeneity from the quantities' dimensions; `collisions()` lists the letters the library uses for two different quantities.

```python
import knowledgebase as kb
kb.graph().census()                          # nodes and edges by kind
j = kb.join("Ideal gas law", "Thermal energy")   # P V/(R n) = E/k_B, pivot: temperature
kb.selftest()                                # structural checks of the tables
```

[models/knowledge.py](models/knowledge.py) uses it in three ways: the selftest; three cross-checks of knowledge-base entries against this repository's models (the Kosterlitz flow against the explicit solution k cot kw, k csc kw of entry 216; the asymptotic-freedom mass law, transmutation slope and blocking of entry 215; the magnon dispersion and Bloch density of entry 271; the Spectre inflation numbers; and the fork's own additions to the library, namely the Chair44 counts, the doubling defect and the reconstruction thresholds, each tied to the model it comes from); and the registration, as `derived` equations, of the joins on the pivots in `REGISTERED_PIVOTS` (statistical mechanics, lattice scaling, running couplings). The library yields more joins than are registered; the rest are counted in the document. [kbdoc.py](kbdoc.py) writes `build/knowledge_graph.tex`, the last part of the document: graph census, readings, dimension report, joins by pivot, collisions, families and an appendix of the readable equations. Every number there is computed when the document is built.

**Fork additions.** `EXTRA_QUANTITIES`, `EXTRA_EQUATIONS`, `EXTRA_READINGS` and `EXTRA_LINKS` in `knowledgebase.py` add 9 equations and 12 quantities from the companion notes to the copied tables, under the same rules (readings, dimensions, collisions). Pseudo-entries such as `spectre-chair44` (`PSEUDO_ENTRIES`) give such results a section of their own in the document.

**References and links.** The graph also holds reference nodes (the manuscripts of entries 155, 215, 216, 229, 268, 271 and 304, read from the `README.md` of each preprint directory, the outside articles the companion notes cite (`EXTERNAL`: [T], [HR], [F], [GS], [SM], [GT22], Tao's blog post, with DOI and review status where given), plus the repositories RosettaMath cites and RosettaMath itself) and one Wikipedia node per article (`slug` of an equation or quantity). Edges `sourced-from`, `treats` and `explained-by` tie equations, entries and quantities to them. In the document, a `paper` equation cites the manuscripts its source string tags (`[Log]`, `[Ess]`, `[B]`, ...), each entry section lists its manuscripts, every equation of the library and every quantity is a labelled row (`K12`, `Q66`, `\label{kbeq:...}`, `\label{kbq:...}`) in the appendices, and bridges and joins link back to those rows and to the paper equations they cross-check with `\ref`/`\eqref`, with names as Wikipedia `\href`s. `build/bibliography.tex` holds the cited references. To cite a new manuscript, add it to `MANUSCRIPTS` in `knowledgebase.py` with the tag used in the model's source strings.

A join is the statement that two quoted equations imply a third by transitivity, with the same letter read the same way taken to be the same variable. Some joins state only a condition (m c^2 = 1/2 m v^2 on energy); that is why registration is by an explicit list of pivots. The checks test the algebra and the dimensions, not whether RosettaMath's entries quote their sources correctly beyond the three cross-checks above.

## The concise document

`export.py --pdf` also builds `build/concise_document.pdf` (hand-editable wrapper `concise_document.tex`,
generated `build/concise_body.tex` and `build/concise_bibliography.tex`): a short paper with no tables,
about physics and geometry only. `concise.py` generates large displays in the manner of Figures 1 and 2
of RosettaMath's `rosettapaper.py`: one quantity per row, overbraced, equal to each library equation
solved for it and underbraced with that equation's name, closed by a side brace with a tally.

- Figure 1 (physics): energy, force, G, k_B, hbar, h and the lattice mass, from `PHYSICS_ROWS`.
- Figure 2 (geometry): the Spectre and hat inflation factors, checked equal as numbers under the
  stated values of the other letters, the window dimension and the Chair44 contact count.
- Figure 3 (transmutation): slopes, coupling and mass, with the O(N) scaling form joined to the gauge-theory
  lattice mass only by arrows taken from the library's links (never by `=`; `NOT_IDENTICAL` in `concise.py` enforces it); the Yang-Mills exponent 1/(2 b_0 g^2) equals r_g*beta is checked.
- Figure 4 (Lie algebras and chirality): half-spinor, generation count, Weyl fields and hypercharge rows, each
  checked equal in one concrete case (the SO(10) 27, three generations, the up quark), E8 counted three ways, and
  arrows for dimension counts that are links and not equations.

Every term is substituted back into its source equation and must give an identity; a failure raises, so
nothing unverified is written. The pivots are named, not computed, because "every family" does not fit
a page and the two inflation factors only look alike. Change `PHYSICS_ROWS`, `GEOMETRY_ROWS`,
`TRANSMUTATION_ROWS` and `SELECTED` in `concise.py` to choose what is shown.

## Exploration scripts

`explore/` holds scripts that test connections the knowledge graph suggests. Each defines `run()`, returns its
prose, a display equation and the checks that passed (a failed check raises), and runs alone or all together:

```sh
python3 sympy_models/explore/su_so.py          # one script
python3 sympy_models/explore/__init__.py        # all eight
```

`concise.py` runs all of them and prints their results in the concise document's section "Explorations".

| Script | Question | Result |
| :--- | :--- | :--- |
| `su_so.py` | where does dim su(N) = dim so(n)? | Pell equation (2n-1)^2 - 8N^2 = -7; only (2,3) and (4,6) are isomorphisms; S(10) = 2 S(6) S(4) = 16 is the Pati-Salam split |
| `spinor_series.py` | which exceptional algebras are so(n) plus spinors? | F4 (n=9), E6 (n=10, with u(1) and two 16s), E7 (n=12, with su(2)), E8 (n=16); n=5 gives 14 = dim G2 by coincidence (checked against G2's closed subsystems) |
| `e8_gradings.py` | which gradings of E8 by one simple root reproduce the repository's Z3 x Z3 table? | scan of all 8 nodes and moduli 2 to 6; the pair (7, 1) exactly, three more pairs in a different arrangement |
| `generation_bound.py` | how many generations fit in an eigenspace? | 4 F N_gen <= 128: two at most; 4 x 16 = 64 = S(14) and 3 x 64 = 192 (an observation about numbers) |
| `weak_mixing.py` | sin^2 = 3/8 at unification, run down | one-loop, from the library's running law: 0.208 (SM), 0.231 (MSSM), measured 0.2312 |
| `hat_golden.py` | why is the hat expansion ratio phi^4? | charpoly (x^2-1)(x^2-7x+1), spectrum {1,-1,phi^4,phi^-4}, tr M^n = 1+(-1)^n+L_4n; N_u/N_r = kappa because the unreflected-hat vector is the H row of M; Spectre unit 4+sqrt15 is not golden; f_H = 1/3 is rational (a CICY docstring says all frequencies are irrational) |
| `ym_mass_gap.py` | where would a Yang-Mills mass gap sit? | two-loop exponent -51/121 for all N; E_gap = c hbar c Lambda_L iff m_L/(a Lambda_L) -> c; blocking drift 2 b0 ln L vs the O(3) ln L/2pi (a correspondence); one-loop Lambda ~ 0.09 GeV. Proves nothing about the gap |
| `kb_paths.py` | what chains join different repositories' results? | shortest chains over mathematical edges (never through where an equation was published, never through pi) |

`explore/_roots.py` builds the root systems of A, D, E6, E7 and E8 from their Dynkin diagrams.

## Layout

- [eqlib.py](eqlib.py): the `Equation` record, the registry, `check()` and the LaTeX block writer.
- [export.py](export.py): discovers modules in `models/`, runs them and writes `build/`.
- [knowledgebase.py](knowledgebase.py), [kbdoc.py](kbdoc.py): the knowledge graph and its document section.
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
| `derived` | Obtained by sympy from other registered equations, including the compositions generated by the knowledge graph. |
| `knowledge-base` | Quoted from the library copied from RosettaMath and cross-checked against a model of the manuscripts. A second transcription, not a reading of the paper. |

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
| [spectre_chair44.py](models/spectre_chair44.py) | companion notes (SC) | [SPECTRE_CHAIR44.md](../SPECTRE_CHAIR44.md), [SPECTRE_CLOCK_AND_CHAIR44.md](../SPECTRE_CLOCK_AND_CHAIR44.md), citing [T], [HR], [SM] | The Chair44 tile's numbers (7 cubes, 24 panels, 192 pyramids, inflation 2I, boundary 24 4^N); the fork's contact formula against its seven counts, its exact defect and the 12/29 2^-n decay; the cube's rotation group is S_4 with 30 subgroups in 11 classes and every subgroup without an element of order 3 inside a D_4; the Spectre's 4 + sqrt 15, window dimension 1.111 and perimeter unit 2 + sqrt 5 | Aperiodicity of either tile; the order-3 exclusion of [HR] Theorem 5; the atlas of 44 contacts; the contact counts themselves (the fork's, from `tools/hierarchy_figures.py`) |
| [chirality_lisi.py](models/chirality_lisi.py) | companion repository (CL) | [brentharts/chirality-standard-model](https://github.com/brentharts/chirality-standard-model): `lisi/THEORY.md`, `lisi/README.md`, citing [DG] Distler-Garibaldi and [Lisi26] | E8 rebuilt from its Cartan matrix alone (240 roots), then read off the roots: the two involutions (120 + 128, 136 + 112), the family grading 86 + 81 + 81, Lisi's grading 92 + 78 + 78 with its five-grading 14 + 64 + 84 + 64 + 14, and the Z3 x Z3 table of THEORY.md section 9 (reproduced exactly); E6 x A2 with the 27 repeating over six classes of 27 and the 45 zero-sum triples of the cubic invariant; D4 x D4 with three blocks of 64 by triality class; so(n) and half-spinor dimensions, so(16) + 128 = 248, 27 = 1 + 16 + 10, spin(11,3) = 6 + 45 + 40, Cl(11,3) = M(128); the magic-square table; the 16's anomaly sums, Q = T3 + Y, the trace ratio 5/3 and sin^2 = 3/8; the dimension count 180 > 128 | Real forms and signatures; the Distler-Garibaldi real-form argument and its factor 2 x 2 (taken from THEORY.md); the triality-algebra nullspace (only the dimensions 0, 2, 9, 28 are used); that Lisi's linearly dependent generations evade the bound |
| [doubling_defect.py](models/doubling_defect.py) | 215 | The O(3) paper, Section 8.2 | The identity Delta = 4 log Z - log Z; the closed form of Proposition 8.3 and its weaker stated bound; that this bound beats the Chair44 relative defect eventually (analogy of form, not of rate) | Lemma 8.2 and Proposition 8.3 |
| [reconstruction.py](models/reconstruction.py) | 229 | The three-state reconstruction paper | The channel is stochastic with eigenvalues 1, lambda, lambda; nonnegative entries on -1/2 <= lambda <= 1; the negative branch exists exactly when d > 4; the community condition is the same line | The theorem |
| [hilbert_smith.py](models/hilbert_smith.py) | 304 | The Hilbert-Smith paper | Nothing: the theorem is recorded as a statement so the notes' discussion can cite it | Everything |
| [knowledge.py](models/knowledge.py) | 215, 216, 229, 271, SC, CL | RosettaMath's library (`knowledge-base`) | Selftest of the tables; the KT flow, asymptotic-freedom scaling and Bloch density against the other models; every registered join substituted back into both equations | RosettaMath's other 150-odd entries beyond structure and dimensions; the physics behind a join |

Notes from reading the sources:
- In the essential-singularity paper, δ = √(b_c − b). The radical is detached from its argument in the PDF text layer, so this was confirmed on the rendered page. Only with this reading do (1.3), (1.6) and (1.8) agree.
- The Bloch coefficient and the 3ζ(5/2)/(128π^{3/2}) correction are in two separate papers, not in the spherical-magnetization paper that CONTENTS lists first.
- A check that passes at p = 3, 5 for the tiling model is not evidence about p > 200. A mutation test (replacing the last nonzero digit by `t mod p`) makes the line-rule check fail on hundreds of words, so the check is not vacuous.

## Caveats

- Display forms are kept symbolic (for example `b_0`, `b_1` instead of their numeric values) so the LaTeX stays readable. Checks run on the substituted form, and a check ties the two together.
- sympy chooses the LaTeX layout. Rearrange in the document if you want a different presentation.
- A `repo-doc` equation is only as reliable as the note it came from.


The concise paper (`concise_document.tex`) now has a dedicated section for `hat_golden` and for `ym_mass_gap`, a conclusion, and omits the knowledge-graph internals (`kb_paths` still runs in the build).
