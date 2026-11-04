# Quantum gravity on the SU(2) ladder (fork note)

Reads the two-plaquette SU(2) lattice model (`sympy_models/models/su2_lattice_ladder.py`) as loop-quantum-gravity kinematics: gauge-invariant states are spin networks on the theta graph, and the link Casimir `j(j+1)` is the squared area in units of `8πγℓ_P²`.

## Contents
- `sympy_models/models/qg_spin_network.py`: 31 checks. Area spectrum, Casimir dictionary, intertwiner counts by Weyl integration (Catalan numbers for spin-1/2), Immirzi parameter `ln2/(π√3)` from spin-1/2 counting, quantum tetrahedron (`±√3/4`), ladder cut area (strong `2/(3g^8)` exact; weak `3/(√10 g^2)` by the virial theorem).
- `sympy_models/knowledgebase.py`: five new equations, four quantities, six references (pseudo-entry `spin-network-qg`).
- `sympy_models/qg_pictures.py`: two figures, every plotted number asserted.
- `sympy_models/qgravity.py` + `qgravity_document.tex`: the third PDF (`build/qgravity_document.pdf`), derived from `concise.py`. Build with `make qg` or `python3 sympy_models/export.py --pdf`.

## Added in the second round
- `models/su2_lattice_strip.py` (15 checks): the three-plaquette strip. Basis, signed 6j plaquette operators (signs fitted and checked against an independent Wigner 3j construction on 5047 transitions), weak-coupling gap `2 sqrt(2 - cos(pi/4)) = 2.2741` (numerics 2.2749 extrapolated), exact strong-coupling splitting (`68/105`, `14/15`, `88/105`), gap ordering strip <= ladder <= plaquette.
- `models/qg_volume.py` (12 checks): the triple-product operator on four-valent intertwiner spaces. A 1D strip has only trivalent vertices; the first four-valent vertex is the corner of an L-shaped three-plaquette cluster. The strong-coupling `<Q^2> = 1/(27 g^16)` there is a leading-order estimate that assumes the b-c two-string amplitude equals the disjoint-string one; the L-shape Hamiltonian is not built.
- `models/qg_einstein4.py` (19 checks): entry 348 (Einstein four-manifolds). so(4) = su(2) + su(2), the Weyl blocks as spin-2 objects (Casimir 6), the curvature blocks of S^4, CP^2, S^2 x S^2 computed from their metrics, chi and tau, the L^2 energy. Four new knowledge-graph equations and a figure of generated rows.
- References: three OpenAI Einstein papers, Gursky-LeBrun, EPRL, Hartshorn-Romanov (Spectre/Chair44, rs-11516041, read), and Romanov's forest-gauge paper (Zenodo 23783316 as supplied, read as the attached file 2610.0069v1.pdf; the paper cites its code record as ...317, not checked whether ...316 is the same deposit). The RosettaMath citation was removed from this paper.

## Third round
- `models/qg_penrose.py` (11 checks, entry 260): Penrose bound saturated by Schwarzschild; charged, 4D and AdS equality cases as algebra; with the LQG area, `M/M_P >= sqrt(N ln2/4 pi)`, 0.2349 at N = 1. A heuristic: the theorem is about classical data. New kb equation `Penrose mass at equality`.
- `models/forest_gauge_3d.py` (23 checks, `forest-gauge`): arithmetic of Romanov's paper only; one rounding remark (0.468 printed for 0.46875). Not part of the quantum gravity paper.

- `models/su2_lattice_lshape.py` (L-shaped cluster with the four-valent corner vertex): phases pinned by a pointwise check against explicit tensor contraction (45 state-plaquette pairs). The strong-coupling estimate `<Q^2> = 1/(27 g^16)` is confirmed (ratio 0.9808, 0.9939, 0.99961, 0.999976 at g^2 = 3, 4, 8, 16). Not reached: the weak-coupling limit; the volume prefactor.

## Fourth round
- `sympy_models/models/forest_gauge_4d_theorems.py` (34 checks, `forest-gauge`): the finite statements of Bhattacharjee-Bhattacharya-Lockwood-Romanov (30 Oct 2026; pages 1-6 and Appendix A read, proofs not read): Table III line field, partner counts 12/14/18, Table I from root systems (17 groups) with Lichnerowicz (equality only for SU(2): the j = 1/2 Casimir), and the numbers of Theorems 1-7. Counted in the paper (new section "Gauge fixing and the curvature of the link"). One printed SU(3) value (0.0177) is 0.00004 above what the stated alpha_3 gives.
- New reference `bbll-forest-4d`. The Zenodo record 10.5281/zenodo.23841532 could not be fetched (host not on the session allowlist), so it is not cited; attach the PDF or paste the abstract to add it.

## Fifth round
- `models/qg_kerr.py` (20 checks, entry 264): Kerr horizon algebra exact; the Penrose gap of Kerr is `M r_-/2`, so entry 260 is strict for Kerr exactly when the inner horizon exists (the condition `r_- > 0` of the near-Kerr papers); `A_+ A_- = (8 pi J)^2`. Heuristic: the Kerr bound `A >= 8 pi |J|` against the spin-1/2 puncture area needs N >= 5 at J = 1/2, N >= 10 at J = 1. New paper section "Strong cosmic censorship near Kerr"; three OpenAI manuscripts added (K2, K1, K0).
- `models/qg_planar_maps.py` (17 checks, entry 211): brute-force rooted planar maps and tree-rooted maps against Tutte and Mullin; KPZ exponents; SU(2)_k invariant counts (closed walks on `A_{k+1}`), loop weight `sqrt q = 2 cos(pi/(k+2))`, q = 1, 2, 3, 4 at k = 1, 2, 4, infinity. New paper section "Liouville quantum gravity and the Catalan counts"; five OpenAI manuscripts added. Liouville and loop quantum gravity are different theories; the link is the Temperley-Lieb loop weight and the Catalan counts only.
- Classical formulas (Tutte, Mullin, KPZ) are quoted by name and checked by brute force or exact algebra; no source was retrieved for them. Quantum-group-deformed LQG is mentioned and not cited.
- `document.py` gains the superscript one (U+00B9) for a manuscript title.

## Sixth round
- `models/su2_berger_top.py` (25 checks) and [SU2_BERGER_TOP.md](SU2_BERGER_TOP.md): the Laplace spectrum of SU(2) with a squashed Hopf fibre (Berger sphere), computed in the Peter-Weyl labels used here. It is a quantum symmetric top, `4[(j(j+1) - m'^2)/R2^2 + m'^2/R1^2]`; the fibre-constant sector is the round S^2 of radius R2/2, and the first nonzero eigenvalue is `min(8, 2 + 1/s^2)/R2^2`. The metric breaks one of a link's two SU(2) symmetries, so this is not a gauge-invariant deformation of the Casimir `j(j+1)` of the round case used above. Not part of the quantum gravity paper.

## Limits
Not a model of quantum gravity: fixed graph, gauge-theory Hamiltonian, no diffeomorphism invariance, volume only for one four-valent vertex. The all-spin Immirzi value (about 0.2375) is quoted from Domagala-Lewandowski, not rederived. References are cited by arXiv ID; the Kogut-Susskind DOI is unchecked. Title and author line are placeholders.

## Extending
Add generated rows to `QG_ROWS` or new figure functions in `qgravity.py`; add results to the model and a section in `write()`.
