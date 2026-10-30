"""The Spectre and Tsiokos's Chair44, from SPECTRE_CHAIR44.md and SPECTRE_CLOCK_AND_CHAIR44.md.

These are the fork's two companion notes, and they cite four sources:
  [T]   I. Tsiokos, A Strongly Aperiodic Monotile in Three Dimensions (arXiv:2609.19214)
  [HR]  B. Hartshorn and T. Romanov, The Internal Space and Boundary Growth of the Spectre and
        Tsiokos's Chair44 Tilings (Research Square, DOI 10.21203/rs.3.rs-11516041/v1; not peer reviewed)
  [F]   F. Flicker, Matching Rules for a Three-Dimensional Strongly Aperiodic Monotile
  [SM]  Smith, Myers, Kaplan and Goodman-Strauss, A chiral aperiodic monotile (the Spectre)

The notes read those papers through fetch-tool summaries, and say so; section numbers here are as
the notes give them. This module checks arithmetic and group theory that the notes state or imply,
exactly, and does not touch the aperiodicity theorems themselves:

  * the Chair44 tile's numbers: 7 cubes, 24 exposed panels, 192 pyramids; inflation 2I, determinant 8;
    the boundary of a supertile has 24 4^n faces, boundary dimension 2;
  * the contact formula C_n = (29/8) 8^n - 3 4^n - 1 against the seven counts the notes computed
    (the counts themselves come from tools/hierarchy_figures.py and are not recomputed here), its
    exact defect, and the 12/29 2^-n relative decay;
  * the rotation group of the cube is S_4, has 30 subgroups in 11 conjugacy classes, and every
    subgroup without an element of order three lies in a copy of D_4 (so |Sym| <= 8 once order three
    is excluded, which is the part of [HR] Theorem 5 that is not group theory);
  * the Spectre's inflation factor 4 + sqrt 15, its window dimension about 1.111, and the unit
    2 + sqrt 5 that the perimeter grows by.

Not tested: that Chair44 is aperiodic, that the Spectre is, [HR]'s exclusion of order-three symmetry,
the atlas of 44 contacts (the notes transcribe it from a text extraction), and anything about the
counts of contacts beyond the seven the notes computed.
"""
import itertools

import sympy as sp
from sympy.combinatorics import Permutation, PermutationGroup

from eqlib import Equation, check, register

NAME = "spectre_chair44"
ENTRY = "spectre-chair44"
T = "[T] Tsiokos"
HR = "[HR] Hartshorn-Romanov"
NOTE_C = "SPECTRE_CLOCK_AND_CHAIR44.md"


def _rotations():
    """The 24 proper cubic rotations as signed permutation matrices of determinant +1."""
    mats = []
    for perm in itertools.permutations(range(3)):
        for signs in itertools.product((1, -1), repeat=3):
            m = sp.zeros(3, 3)
            for row, (col, sg) in enumerate(zip(perm, signs)):
                m[row, col] = sg
            if m.det() == 1:
                mats.append(m)
    return mats


def build_and_verify():
    # ---- the Chair44 tile ------------------------------------------------------------------
    cubes = {c for c in itertools.product((0, 1), repeat=3)} - {(1, 1, 1)}
    exposed = sum(1 for c in cubes for ax in range(3) for d in (-1, 1)
                  if tuple(c[i] + (d if i == ax else 0) for i in range(3)) not in cubes)
    eq_tile = register(NAME, Equation(
        key="chair44_tile_numbers", entry=ENTRY,
        title="The Chair44 tile: cubes, panels and pyramids",
        source_kind="paper", source=f"{T} Sections 2.1, 2.3, 2.4 (as the fork's notes cite them)",
        latex_str=(r"\#\text{cubes}=7,\qquad \#\text{exposed panels}=24,"
                   r"\qquad \#\text{pyramids}=8\cdot 24=192"),
        notes="A 2x2x2 block of unit cubes with one corner removed, volume 7, with eight small square "
              "pyramids on each exposed unit panel. The pyramids' profile is what makes the tile aperiodic "
              "and is not modelled here.",
    ))
    check(NAME, "the corner-removed 2x2x2 block has 7 cubes", len(cubes) == 7, eq=eq_tile)
    check(NAME, "its surface has 24 exposed unit panels (removing a corner swaps 3 outer panels for 3 inner)",
          exposed == 24, eq=eq_tile)
    check(NAME, "8 pyramids on each of 24 panels gives 192", 8 * exposed == 192, eq=eq_tile)

    lam = 2 * sp.eye(3)
    eq_infl = register(NAME, Equation(
        key="chair44_inflation_matrix", entry=ENTRY,
        title="The Chair44 inflation",
        source_kind="paper", source=f"{T} Section 3.1; {HR} Section 6",
        latex_str=r"\Lambda = 2\,I_3,\qquad \det\Lambda = 8,\qquad \mu_C = 2",
        notes="Scaling by 2 with no rotation or mirror: the twist is zero ([HR] Section 8.1, Fact 19). "
              "A level-n supertile is made of 8^n chairs.",
    ))
    check(NAME, "det(2 I_3) = 8, the number of child chairs", lam.det() == 8, eq=eq_infl)
    check(NAME, "the eigenvalues are all 2 (a pure scaling)", set(lam.eigenvals()) == {2}, eq=eq_infl)

    n, N = sp.symbols("n N", positive=True, integer=True)
    boundary = 24 * 4**N
    eq_bd = register(NAME, Equation(
        key="chair44_boundary_faces", entry=ENTRY,
        title="Boundary of a Chair44 supertile",
        source_kind="paper", source=f"{HR} Section 7.2",
        latex_str=r"B_N = 24\cdot 4^{N},\qquad D_b=\lim_{N\to\infty}\frac{\log B_N}{\log 2^{N}}=2=d-1",
        notes="No fractal excess: the boundary dimension equals the ambient dimension minus one. "
              "Contrast the Spectre's supertile boundary, of box dimension about 1.3993 ([HR] abstract).",
    ))
    check(NAME, "log(24 4^N) / log(2^N) tends to 2 as N -> infinity",
          sp.limit(sp.log(boundary) / sp.log(2**N), N, sp.oo) == 2, eq=eq_bd)
    check(NAME, "the surface of a level-N supertile at N = 0 is the 24 panels of one chair",
          boundary.subs(N, 0) == 24, eq=eq_bd)

    # ---- the contact formula of the fork's notes -----------------------------------------------
    counts = {1: 16, 2: 183, 3: 1663, 4: 14079, 5: 115711, 6: 937983, 7: 7553023}
    C = sp.Rational(29, 8) * 8**n - 3 * 4**n - 1
    eq_c = register(NAME, Equation(
        key="chair44_contact_formula", entry=ENTRY,
        title="Face contacts in the Chair44 hierarchy",
        source_kind="repo-doc", source=f"{NOTE_C} Section 2 [computed], from tools/hierarchy_figures.py",
        expr=sp.Eq(sp.Symbol("C_n"), C),
        latex_str=(r"C_n=\frac{29}{8}\,8^{n}-3\cdot 4^{n}-1,%%BR%%"
                   r"C_{n+1}-8C_n=12\cdot 4^{n}+7"),
        notes="Level-n patch of the one-seed hierarchy (8^n chairs). Levels 1 to 4 were counted cell by "
              "cell and levels 5 to 7 from poses, assuming every contact is an atlas contact. The counts "
              "are the fork's, not recomputed here; the formula was fitted on levels 1 to 4. The volume "
              "term is 29/8 contacts per chair, the 4^n term scales with the area of the supertile, and "
              "-1 is the corner.",
    ))
    check(NAME, "the formula reproduces all seven counts of the notes: 16, 183, 1663, 14079, 115711, 937983, 7553023",
          all(C.subs(n, k) == v for k, v in counts.items()), eq=eq_c)
    check(NAME, "it was fitted on levels 1 to 4 and predicts levels 5 to 7 exactly",
          all(C.subs(n, k) == counts[k] for k in (5, 6, 7)), eq=eq_c)
    defect = sp.simplify(C.subs(n, n + 1) - 8 * C)
    check(NAME, "C_{n+1} - 8 C_n = 12 4^n + 7 exactly (the extensivity defect)",
          sp.simplify(defect - (12 * 4**n + 7)) == 0, eq=eq_c)
    rel = sp.limit(sp.simplify(defect / C.subs(n, n + 1)) * 2**n, n, sp.oo)
    check(NAME, "the relative defect |8 C_n - C_{n+1}| / C_{n+1} times 2^n tends to 12/29, so it decays like 2^-n",
          rel == sp.Rational(12, 29), eq=eq_c)
    check(NAME, "contacts per chair tend to 29/8 (the volume term)",
          sp.limit(C / 8**n, n, sp.oo) == sp.Rational(29, 8), eq=eq_c)

    # ---- the atlas of poses ------------------------------------------------------------------
    eq_atlas = register(NAME, Equation(
        key="chair44_atlas_counts", entry=ENTRY,
        title="Legal contacts of Chair44",
        source_kind="paper", source=f"{T} Section 3.1; {HR} Section 7.1 and Corollary 1",
        latex_str=(r"44\ \text{legal of}\ 2388\ \text{candidate poses},\qquad 1194\ \text{are mirror images},"
                   r"\qquad 30\ \text{of the 44 occur in a tiling}"),
        notes="A Lean theorem in [T] (native_decide) gives 44; [HR] Corollary 1 gives 30. The legality test "
              "(pyramid profiles on a shared panel) is not modelled; the fork's tools/chair44_model.py uses "
              "the transcribed atlas.",
    ))
    check(NAME, "1194 + 1194 = 2388: exactly half of the candidate poses are mirror images", 2 * 1194 == 2388, eq=eq_atlas)
    check(NAME, "44 - 30 = 14 legal contacts never occur in a tiling", 44 - 30 == 14, eq=eq_atlas)

    # ---- the rotation group, its subgroups, and the bound |Sym| <= 8 ---------------------------------
    rots = _rotations()
    eq_sym = register(NAME, Equation(
        key="chair44_symmetry_bound", entry=ENTRY,
        title="Symmetry groups of Chair44 tilings",
        source_kind="paper", source=f"{HR} Section 7.4, Theorem 5 and Remark 1",
        latex_str=(r"\mathrm{Sym}(T)\ \sim\ 1,\ C_2,\ C_4\ \text{or}\ D_4,\qquad |\mathrm{Sym}(T)|\le 8,"
                   r"\qquad\text{attained}"),
        notes="Each occurs. Elements of order three are excluded by the paper's argument (not reproduced "
              "here); the bound 8 then follows from the group theory checked below. The result depends on "
              "an input (R) quoted from [T], audited but not refereed. [T] Theorem 7.2 and Corollary 7.3 give "
              "the weaker 'mildly aperiodic' and 'no symmetry of infinite order'.",
    ))
    check(NAME, "there are 24 proper cubic rotations", len(rots) == 24, eq=eq_sym)
    diagonals = [(1, 1, 1), (1, 1, -1), (1, -1, 1), (-1, 1, 1)]

    def diag_perm(m):
        images = []
        for d in diagonals:
            v = tuple(m * sp.Matrix(d))
            for j, e in enumerate(diagonals):
                if v == e or v == tuple(-x for x in e):
                    images.append(j)
        return Permutation(images)

    perms = [diag_perm(m) for m in rots]
    check(NAME, "the action on the 4 body diagonals is faithful, so the rotation group is a subgroup of S_4 of order 24",
          len(set(map(tuple, (p.array_form for p in perms)))) == 24, eq=eq_sym)
    S4 = PermutationGroup(perms)
    check(NAME, "that group is all of S_4 (order 24)", S4.order() == 24, eq=eq_sym)
    elems = list(S4.elements)
    subgroups = set()
    for g, h in itertools.product(elems, repeat=2):
        subgroups.add(frozenset(PermutationGroup([g, h]).elements))
    check(NAME, "S_4 has 30 subgroups (every subgroup of S_4 is generated by at most two elements)",
          len(subgroups) == 30, eq=eq_sym)

    def conj(sub, g):
        return frozenset(g * x * g**-1 for x in sub)
    classes, seen = [], set()
    for sub in subgroups:
        if sub in seen:
            continue
        cl = {conj(sub, g) for g in elems}
        seen |= cl
        classes.append(cl)
    check(NAME, "the 30 subgroups fall into 11 conjugacy classes", len(classes) == 11, eq=eq_sym)
    no3 = [cl for cl in classes if all(x.order() != 3 for x in next(iter(cl)))]
    check(NAME, "7 classes contain no element of order 3 (the 2-groups: 1, C2 twice, C4, V4 twice, D4)",
          len(no3) == 7, eq=eq_sym)
    d4 = [s for s in subgroups if len(s) == 8]
    check(NAME, "there are 3 copies of D_4 (the Sylow 2-subgroups)", len(d4) == 3, eq=eq_sym)
    check(NAME, "every subgroup with no element of order 3 lies in a copy of D_4, so has order at most 8",
          all(any(sub <= big for big in d4) for cl in no3 for sub in cl)
          and max(len(next(iter(cl))) for cl in no3) == 8, eq=eq_sym)
    check(NAME, "the classes that do contain an element of order 3 include groups of order 12 and 24, "
                "so the order-3 exclusion is what gives the bound 8",
          max(len(next(iter(cl))) for cl in classes if cl not in no3) == 24, eq=eq_sym)

    # ---- the Spectre ---------------------------------------------------------------------------------
    mu = 4 + sp.sqrt(15)
    eq_mu = register(NAME, Equation(
        key="spectre_inflation", entry=ENTRY,
        title="The Spectre's area inflation factor",
        source_kind="paper", source=f"{HR} Section 6, Fact 16",
        latex_str=r"\lambda^2=\mu=4+\sqrt{15}\in\mathbb Z[\sqrt{15}],\qquad \mu^2=8\mu-1,\qquad \mu^{-1}=4-\sqrt{15}",
        notes="lambda is the linear factor, (sqrt 6 + sqrt 10)/2. The knowledge base carries the same numbers "
              "('Spectre inflation factor' and its neighbours), checked against this entry in the model knowledge.",
    ))
    x = sp.symbols("x")
    check(NAME, "mu = 4 + sqrt 15 is a root of x^2 - 8x + 1", sp.expand(mu**2 - 8 * mu + 1) == 0, eq=eq_mu)
    check(NAME, "its reciprocal is 4 - sqrt 15, so mu is a unit of norm 1 in Z[sqrt 15]",
          sp.simplify(1 / mu - (4 - sp.sqrt(15))) == 0 and 16 - 15 == 1, eq=eq_mu)
    ell = (sp.sqrt(6) + sp.sqrt(10)) / 2
    check(NAME, "the linear factor (sqrt 6 + sqrt 10)/2 squares to mu", sp.simplify(sp.expand(ell**2) - mu) == 0, eq=eq_mu)
    check(NAME, "the minimal polynomial of mu over Q is x^2 - 8x + 1",
          sp.minimal_polynomial(mu, x) == x**2 - 8 * x + 1, eq=eq_mu)

    window = sp.log(5 + 2 * sp.sqrt(6)) / sp.log(mu)
    eq_win = register(NAME, Equation(
        key="spectre_window_dimension", entry=ENTRY,
        title="Dimension of the Spectre's window",
        source_kind="paper", source=f"{HR} Theorem 3, Section 6",
        expr=sp.Eq(sp.Symbol("D"), window),
        latex_str=r"D=\frac{\log(5+2\sqrt6)}{\log(4+\sqrt{15})}\approx 1.111",
        notes="The window of the cut-and-project description in the 4D vertex module Z[zeta_12]. The "
              "supertile boundary has box dimension about 1.3993 ([HR] abstract), not modelled here.",
    ))
    check(NAME, "D is about 1.111", abs(sp.N(window, 12) - sp.Float("1.1111", 5)) < 1e-3, eq=eq_win)
    check(NAME, "5 + 2 sqrt 6 = (sqrt 2 + sqrt 3)^2, so D = 2 log(sqrt 2 + sqrt 3) / log mu",
          sp.expand((sp.sqrt(2) + sp.sqrt(3))**2 - (5 + 2 * sp.sqrt(6))) == 0, eq=eq_win)
    check(NAME, "1 < D < 2: strictly between a line and the plane", 1 < sp.N(window) < 2, eq=eq_win)

    k = sp.symbols("k", positive=True, integer=True)
    phi = (1 + sp.sqrt(5)) / 2
    fib = lambda m: (phi**m - (1 - phi)**m) / sp.sqrt(5)
    eq_per = register(NAME, Equation(
        key="spectre_perimeter_growth", entry=ENTRY,
        title="Perimeter of a Spectre supertile",
        source_kind="paper", source=f"{HR} Theorem 1 (the notes do not give its section number)",
        latex_str=r"P_k=2F_{3k+5}+4,\qquad \frac{P_{k+1}}{P_k}\longrightarrow 2+\sqrt5=\phi^{3}",
        notes="F_n are the Fibonacci numbers. The perimeter grows by the unit 2 + sqrt 5 of Q(sqrt 5), while the "
              "area grows by mu = 4 + sqrt 15, which is why the boundary is fractal.",
    ))
    check(NAME, "F_{3k+5} from Binet's formula matches the Fibonacci numbers 5, 21, 89 at k = 0, 1, 2",
          [sp.simplify(fib(3 * kk + 5)) for kk in range(3)] == [5, 21, 89], eq=eq_per)
    check(NAME, "phi^3 = 2 + sqrt 5", sp.simplify(sp.expand(phi**3) - (2 + sp.sqrt(5))) == 0, eq=eq_per)
    # Binet: sqrt5 P_k / phi^(3k) = 2 (phi^m - psi^m u) + 4 sqrt5 t with m = 5 (for P_k) or 8 (for P_{k+1}),
    # t = phi^(-3k) and u = (psi/phi)^(3k); both tend to 0 because |psi/phi| = phi^-2 < 1.
    psi = (1 - sp.sqrt(5)) / 2
    t, u = sp.symbols("t u")
    num = 2 * (phi**8 - psi**8 * u) + 4 * sp.sqrt(5) * t
    den = 2 * (phi**5 - psi**5 * u) + 4 * sp.sqrt(5) * t
    check(NAME, "|psi / phi| < 1, so t and u tend to 0 along k -> infinity", abs(sp.N(psi / phi)) < 1, eq=eq_per)
    check(NAME, "hence P_{k+1} / P_k tends to phi^3 = 2 + sqrt 5 (Binet's formula, ratio taken at t = u = 0)",
          sp.simplify(sp.expand((num / den).subs({t: 0, u: 0})) - (2 + sp.sqrt(5))) == 0, eq=eq_per)
    check(NAME, "(2 + sqrt 5)^2 is about 17.9 against mu about 7.87: the perimeter squared outgrows the area, "
                "which a smooth boundary (perimeter^2 ~ area) would not allow",
          sp.N((2 + sp.sqrt(5))**2) > sp.N(mu), eq=eq_per)
