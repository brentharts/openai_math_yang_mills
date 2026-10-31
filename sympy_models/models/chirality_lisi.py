"""E8, the chirality theorem and Lisi's triality, from the fork brentharts/chirality-standard-model.

The repository (lisi/README.md, lisi/THEORY.md, paper/paper.tex) checks finite facts about E8 at the
root level and about the Standard Model content of the 27. This module does not copy its scripts.
It rebuilds the root system of E8 from the Cartan matrix alone (Bourbaki numbering, closure under
the simple reflections), reads every decomposition off the coefficients of the roots, and compares
with the numbers the repository states:

  * 240 roots, and the gradings by the coefficient of one simple root: node 8 and node 1 give the
    two involutions (E7 x A1 and D8) and node 7 the order-three family grading E6 x A2;
  * E8 = E6 x A2 + (27,3) + (27bar,3bar); the 27 repeats over the three weights of the A2; the 45
    zero-sum triples of the E6 cubic invariant;
  * E8 > D4 x D4 with 240 = 24 + 24 + 3.64, and E6 roots 72 = 24 + 3.16;
  * E8 > D7 x U(1), the five-grading 14 + 64 + 84 + 64 + 14 (Lisi's order-three grading);
  * the two order-three gradings together make Z3 x Z3 with the table of THEORY.md section 9;
  * classical dimension formulas: so(n), half-spinors, so(16) + 128 = 248, 27 = 1 + 16 + 10,
    spin(11,3) = 6 + 45 + 40, Cl(11,3) = M(128, R);
  * the 16 of SO(10): anomaly sums, Q = T3 + Y, the hypercharge normalisation 5/3 and sin^2 = 3/8;
  * the Distler-Garibaldi dimension count, with the eigenspace dimensions taken from the roots.

Not tested: that these structures are the right ones for physics, the real forms (signatures) that
the repository says a root computation cannot see, the Distler-Garibaldi real-form argument, the
triality-algebra nullspace computation (only the dimensions 0, 2, 9, 28 are used, and the classical
formulas below recover them), and everything the repository marks [not checked].
"""
import itertools
from collections import Counter
from fractions import Fraction

import sympy as sp

from eqlib import Equation, check, register

NAME = "chirality_lisi"
ENTRY = "chirality-lisi"
REPO = "brentharts/chirality-standard-model"
THEORY = f"{REPO}: lisi/THEORY.md"
DG = "[DG] Distler-Garibaldi"
LISI = "[Lisi26] Lisi"

# Bourbaki labelling: 1-3-4-5-6-7-8 and 2 attached to 4.
EDGES = [(1, 3), (3, 4), (4, 2), (4, 5), (5, 6), (6, 7), (7, 8)]
N = 8


def cartan():
    a = [[2 if i == j else 0 for j in range(N)] for i in range(N)]
    for i, j in EDGES:
        a[i - 1][j - 1] = a[j - 1][i - 1] = -1
    return a


A = cartan()


def ip(c, d):
    """Inner product of two roots given by simple-root coefficients (E8 is simply laced)."""
    return sum(c[i] * A[i][j] * d[j] for i in range(N) for j in range(N))


def roots_of_e8():
    simple = [tuple(1 if j == i else 0 for j in range(N)) for i in range(N)]
    found, frontier = set(simple), list(simple)
    while frontier:
        new = []
        for beta in frontier:
            for i in range(N):
                pairing = sum(beta[j] * A[j][i] for j in range(N))
                r = tuple(b - (pairing if k == i else 0) for k, b in enumerate(beta))
                if r not in found:
                    found.add(r)
                    new.append(r)
        frontier = new
    return found | {tuple(-x for x in r) for r in found}


def classify(roots, key):
    return Counter(key(r) for r in roots)


def build_and_verify():
    roots = roots_of_e8()
    positive = [r for r in roots if all(x >= 0 for x in r)]

    # ---- E8 from its Cartan matrix ---------------------------------------------------------
    eq_roots = register(NAME, Equation(
        key="cl_e8_roots", entry=ENTRY, title="E8 rebuilt from its Cartan matrix",
        source_kind="derived",
        source=f"Bourbaki Cartan matrix of E8, closed under simple reflections; compared with {THEORY} section 0",
        latex_str=r"|\Phi(E_8)|=240=2\cdot 120,\qquad \dim E_8=240+8=248",
        notes="Every later count of this section is read off the coefficients of these 240 roots, not taken "
              "from the repository's coordinate scripts, so the two constructions agree or the build fails.",
    ))
    check(NAME, "the simple-reflection closure of the eight simple roots has 240 elements", len(roots) == 240, eq=eq_roots)
    check(NAME, "120 are positive (all coefficients non-negative) and every root has norm 2",
          len(positive) == 120 and all(ip(r, r) == 2 for r in roots), eq=eq_roots)
    check(NAME, "roots plus rank is the dimension: 240 + 8 = 248", len(roots) + N == 248, eq=eq_roots)
    check(NAME, "the Cartan matrix is symmetric with determinant 1 (E8 is simply laced and unimodular)",
          sp.Matrix(A).is_symmetric() and sp.Matrix(A).det() == 1, eq=eq_roots)

    # ---- the three one-node gradings -------------------------------------------------------
    by7 = classify(roots, lambda r: r[6] % 3)
    by1 = classify(roots, lambda r: r[0] % 3)
    five = classify(roots, lambda r: r[0])
    par1 = classify(roots, lambda r: r[0] % 2)
    par8 = classify(roots, lambda r: r[7] % 2)
    eq_grad = register(NAME, Equation(
        key="cl_e8_gradings", entry=ENTRY, title="Three gradings of E8: two involutions and two order-three gradings",
        source_kind="repo-doc", source=f"{THEORY} sections 7-9; {DG} section 3; {LISI} section 6 (as the repository reads them)",
        latex_str=(r"\begin{aligned}"
                   r"248&=120+128=136+112\\"
                   r"248&=86+81+81\\"
                   r"248&=92+78+78=(14+64)+(84+8)+(64+14)"
                   r"\end{aligned}"),
        notes="The first line is the two inner involutions (fixed algebras so(16) and e7 + su(2)); the second is the "
              "family grading of the chirality theorem; the third is Lisi's order-three grading, whose fixed algebra "
              "is so(14) + u(1) and whose five-grading by the U(1) charge has 14, 64, 84, 64, 14 roots. They are "
              "different gradings of one algebra. Root level only: a real form is not seen.",
    ))
    check(NAME, "node 1 even/odd: 112 + 128 roots, so with the 8 Cartan directions the fixed algebra of the first "
                "involution has dimension 120 = dim so(16) and the other eigenspace has 128",
          par1[0] + 8 == 120 and par1[1] == 128, eq=eq_grad)
    check(NAME, "node 8 even/odd: 128 + 112 roots, fixed algebra 136 = 133 + 3 (e7 + su(2)), other eigenspace 112 "
                "(the dimensions the dimension no-go of Distler-Garibaldi uses)",
          par8[0] + 8 == 136 and par8[1] == 112, eq=eq_grad)
    check(NAME, "node 7 mod 3: 86 + 81 + 81, the family grading with fixed algebra e6 + su(3)",
          (by7[0] + 8, by7[1], by7[2]) == (86, 81, 81), eq=eq_grad)
    check(NAME, "node 1 mod 3: 92 + 78 + 78, Lisi's grading with fixed algebra so(14) + u(1)",
          (by1[0] + 8, by1[1], by1[2]) == (92, 78, 78), eq=eq_grad)
    check(NAME, "node 1 gives the five-grading 14 + 64 + 84 + 64 + 14 (84 = roots of so(14))",
          [five[k] for k in (-2, -1, 0, 1, 2)] == [14, 64, 84, 64, 14], eq=eq_grad)
    check(NAME, "92 + 78 + 78 = 248 and 86 + 81 + 81 = 248 and the two gradings have the same total",
          92 + 78 + 78 == 86 + 81 + 81 == 248, eq=eq_grad)

    # ---- the Z3 x Z3 table -----------------------------------------------------------------
    joint = Counter((r[6] % 3, r[0] % 3) for r in roots)
    joint[(0, 0)] += 8
    table = [[joint[(i, j)] for j in range(3)] for i in range(3)]
    eq_z33 = register(NAME, Equation(
        key="cl_z3_z3_table", entry=ENTRY, title="The two order-three gradings together: a Z3 x Z3 grading of E8",
        source_kind="repo-doc", source=f"{THEORY} section 9 (check_lisi_theta.py)",
        latex_str=(r"\begin{array}{c|ccc} & k=0 & k=1 & k=2\\\hline \chi=0 & 50 & 18 & 18\\ \chi=1 & 21 & 42 & 18\\ "
                   r"\chi=2 & 21 & 18 & 42\end{array}"
                   r"\qquad\text{rows: }86,\,81,\,81;\ \text{columns: }92,\,78,\,78"),
        notes="Rows are the family grading chi (e6 + su(3) fixed), columns Lisi's grading k (so(14) + u(1) fixed). "
              "The rows and columns reproduce the two gradings above. The repository computes the table from "
              "its own root coordinates; here it comes from the coefficients of two simple roots.",
    ))
    check(NAME, "the joint table of node 7 mod 3 against node 1 mod 3 is exactly the repository's",
          table == [[50, 18, 18], [21, 42, 18], [21, 18, 42]], eq=eq_z33)
    check(NAME, "row sums 86, 81, 81 and column sums 92, 78, 78",
          [sum(r) for r in table] == [86, 81, 81] and [sum(c) for c in zip(*table)] == [92, 78, 78], eq=eq_z33)

    # ---- E6 x A2, the 27 over the 3, the cubic invariant -----------------------------------
    a2 = [tuple(1 if j == 2 else 0 for j in range(N)), tuple(1 if j == 3 else 0 for j in range(N))]  # alpha_3, alpha_4
    perp = [r for r in roots if all(ip(r, s) == 0 for s in a2)]          # the centraliser E6
    in_a2 = [r for r in roots if all(r[k] == 0 for k in range(N) if k not in (2, 3))]
    rest = [r for r in roots if r not in set(perp) and r not in set(in_a2)]
    classes = classify(rest, lambda r: (ip(r, a2[0]), ip(r, a2[1])))
    eq_fam = register(NAME, Equation(
        key="cl_family_module", entry=ENTRY, title="E8 over E6 x A2: the 27 repeats over the three weights",
        source_kind="repo-doc", source=f"{THEORY} section 1 and lisi/check_e6_su3_bridge.py",
        latex_str=r"240=72+6+6\cdot 27,\qquad 27\boxtimes 3:\ 81\ \text{weights},\qquad 45\ \text{zero-sum triples}",
        notes="The roots orthogonal to an A2 form E6 (72 roots); the other 162 fall into six classes of 27 by their "
              "pairing with the A2, the three weights of the 3 and of its dual. The 45 triples of weights of the 27 "
              "that sum to zero are the monomials of the E6 cubic invariant, which is why E6 is chiral without "
              "a cubic anomaly.",
    ))
    check(NAME, "the roots orthogonal to an A2 of E8 number 72 (E6) and the A2 has 6, so 240 - 78 = 162 remain",
          len(perp) == 72 and len(in_a2) == 6 and len(rest) == 162, eq=eq_fam)
    check(NAME, "the 162 fall into six classes, each of 27 roots: the 27 repeats over three weights and their negatives",
          len(classes) == 6 and set(classes.values()) == {27}, eq=eq_fam)
    # weights of the 27: the class with pairings (1, 0), projected off the A2
    w27 = [r for r in rest if (ip(r, a2[0]), ip(r, a2[1])) == (1, 0)]
    triples = [t for t in itertools.combinations(range(27), 3) if _zero_sum_in_e6(w27, t)]
    check(NAME, "there are 45 triples of weights of the 27 whose E6 parts sum to zero (the cubic invariant)",
          len(triples) == 45, eq=eq_fam)
    check(NAME, "the six classes are exchanged by negation, so 27bar = -27 over the dual weights",
          all(tuple(-x for x in r) in set(rest) for r in w27), eq=eq_fam)

    # ---- D4 x D4, E6 over D4 ---------------------------------------------------------------
    d4 = [tuple(1 if j == i else 0 for j in range(N)) for i in (1, 2, 3, 4)]       # alpha_2..alpha_5
    d4_roots = [r for r in roots if all(r[k] == 0 for k in range(N) if k not in (1, 2, 3, 4))]
    d4_perp = [r for r in roots if all(ip(r, s) == 0 for s in d4)]
    mixed = [r for r in roots if r not in set(d4_roots) and r not in set(d4_perp)]
    eq_d4 = register(NAME, Equation(
        key="cl_d4_d4", entry=ENTRY, title="E8 over D4 x D4 and E6 over D4: the triality blocks",
        source_kind="repo-doc", source=f"{THEORY} section 1 and lisi/check_triality_d4d4.py",
        latex_str=r"240=24+24+3\cdot 64,\qquad 72=24+3\cdot 16",
        notes="A D4 inside E8 has a D4 as centraliser. The 192 remaining roots are three blocks of 64, one for each "
              "of the three eight-dimensional representations 8v, 8s, 8c that triality permutes. Inside E6 the same "
              "D4 leaves 3.16.",
    ))
    check(NAME, "a D4 spanned by alpha_2..alpha_5 has 24 roots and its centraliser has 24 more",
          len(d4_roots) == 24 and len(d4_perp) == 24, eq=eq_d4)
    def triality_class(r):
        """The class of the first D4's weight in P/Q = Z2 x Z2: 8v, 8s or 8c (the outer nodes 2, 3, 5 of D4)."""
        a2, a3, _a4, a5 = (ip(r, s) for s in d4)
        return ((a2 + a5) % 2, (a3 + a5) % 2)
    blocks = classify(mixed, triality_class)
    check(NAME, "the other 192 roots fall into the three nontrivial classes 8v, 8s, 8c of the first D4's weight "
                "lattice, 64 roots each (3 x 8 x 8)",
          len(mixed) == 192 and sorted(blocks.values()) == [64, 64, 64] and (0, 0) not in blocks, eq=eq_d4)
    in_e6 = set(perp)
    d4_in_e6 = [r for r in d4_perp if r in in_e6]
    check(NAME, "the centraliser D4 lies inside E6, and E6 = 24 + 48 = 24 + 3 x 16",
          len(d4_in_e6) == 24 and len(in_e6) - 24 == 48, eq=eq_d4)

    # ---- classical dimension formulas ------------------------------------------------------
    n = sp.Symbol("n", positive=True)
    D = sp.Symbol("D")
    so_dim = lambda m: sp.Rational(m * (m - 1), 2)
    half = lambda m: 2 ** (sp.Rational(m, 2) - 1)
    eq_so = register(NAME, Equation(
        key="cl_orthogonal_dimension", entry=ENTRY, title="Dimension of so(n)",
        source_kind="textbook", source="the dimension of the antisymmetric n x n matrices",
        expr=sp.Eq(D, n * (n - 1) / 2), latex_str=None,
        notes="With n = 8 this is 28, the dimension of the triality algebra tri(O); with n = 10, 45; with n = 14 "
              "(the signature (11,3)), 91; with n = 16, 120.",
    ))
    check(NAME, "n = 8, 10, 12, 14, 16 give 28, 45, 66, 91, 120",
          [so_dim(m) for m in (8, 10, 12, 14, 16)] == [28, 45, 66, 91, 120], eq=eq_so)
    check(NAME, "so(8) has 24 roots and rank 4: 24 + 4 = 28, so the D4 above is the triality algebra of the octonions",
          so_dim(8) == len(d4_roots) + 4, eq=eq_so)
    check(NAME, "so(14) has 84 roots and rank 7: 84 + 7 = 91, and the 84 are the roots at level 0 of the five-grading",
          so_dim(14) == five[0] + 7, eq=eq_so)

    S = sp.Symbol("S")
    eq_half = register(NAME, Equation(
        key="cl_half_spinor", entry=ENTRY, title="Dimension of a half-spinor of so(n), n even",
        source_kind="textbook", source="Spin(n) for even n has two half-spinor representations of dimension 2^(n/2-1)",
        expr=sp.Eq(S, 2 ** (n / 2 - 1)), latex_str=None,
        notes="n = 8 gives the 8s and 8c of triality, n = 10 the 16 of the Standard Model generation, n = 12 the 32, "
              "n = 14 the 64 that Lisi's model uses and n = 16 the 128.",
    ))
    check(NAME, "n = 8, 10, 12, 14, 16 give 8, 16, 32, 64, 128", [half(m) for m in (8, 10, 12, 14, 16)] == [8, 16, 32, 64, 128],
          eq=eq_half)
    check(NAME, "the 64 of so(14) is the 64 of the five-grading, and 128 + 120 = 248 (so(16) plus a half-spinor)",
          half(14) == five[1] and so_dim(16) + half(16) == 248, eq=eq_half)

    eq_so16 = register(NAME, Equation(
        key="cl_e8_so16", entry=ENTRY, title="E8 as so(16) plus a half-spinor",
        source_kind="derived", source=f"cl_orthogonal_dimension and cl_half_spinor at n = 16; {DG} section 3",
        expr=sp.Eq(D, n * (n - 1) / 2 + 2 ** (n / 2 - 1)), latex_str=None,
        notes="This is a dimension count, not a theorem about all n: it is the Lie algebra E8 only at n = 16. The "
              "fixed algebra of the first involution is so(16) and its other eigenspace is the 128.",
    ))
    values = {m: so_dim(m) + half(m) for m in (8, 10, 12, 14, 16)}
    check(NAME, "the count is 248 at n = 16 and not at n = 8, 10, 12, 14 (36, 61, 98, 155)",
          values == {8: 36, 10: 61, 12: 98, 14: 155, 16: 248}, eq=eq_so16)
    check(NAME, "its split 120 + 128 is the node 1 involution of the roots", (120, 128) == (par1[0] + 8, par1[1]), eq=eq_so16)

    DV = sp.Symbol("D_V")
    eq_27 = register(NAME, Equation(
        key="cl_e6_fundamental", entry=ENTRY, title="The 27 of E6 under SO(10) x U(1)",
        source_kind="repo-doc", source=f"{THEORY} section 2 (sm_charges.py)",
        expr=sp.Eq(DV, 1 + 2 ** (n / 2 - 1) + n), latex_str=None,
        notes="At n = 10: 27 = 1 + 16 + 10, the singlet, the half-spinor and the vector. The 16 is one Standard "
              "Model generation with a right-handed neutrino; the 10 is vector-like.",
    ))
    f27 = lambda m: 1 + half(m) + m
    check(NAME, "n = 10 gives 27 (and n = 8 gives 17, n = 12 gives 45): the fundamental of E6 is the SO(10) case",
          (f27(10), f27(8), f27(12)) == (27, 17, 45), eq=eq_27)
    check(NAME, "3 x 27 = 81 and 2 x 81 + 86 = 248, the adjoint branching of E8 over E6 x SU(3)",
          3 * f27(10) == 81 and 2 * 81 + 86 == 248, eq=eq_27)

    # ---- spin(11,3), Pati-Salam, Clifford ---------------------------------------------------
    eq_spin = register(NAME, Equation(
        key="cl_spin113", entry=ENTRY, title="spin(11,3), Pati-Salam, the Standard Model and the Majorana-Weyl 64",
        source_kind="repo-doc", source=f"{THEORY} section 3; {LISI} section 2 (as the repository reads it)",
        latex_str=(r"91=6+45+40,\quad 21=15+3+3,\quad 12=8+3+1,\quad 2^{14}=128^2,\quad 128=64+64"),
        notes="spin(11,3) has the dimension of so(14); it splits as spin(1,3) + spin(10) + the 40 frame-Higgs "
              "generators (4 times the 10). Cl(11,3) is the full real matrix algebra of size 128; its chirality "
              "element splits the 128 in two halves of 64. The signature, which this arithmetic cannot see, is the "
              "point of the repository's discussion.",
    ))
    check(NAME, "6 + 45 + 40 = 91 = dim so(14), with 6 = dim so(4), 45 = dim so(10) and 40 = 4 x 10",
          6 + 45 + 40 == so_dim(14) == 91 and so_dim(4) == 6 and so_dim(10) == 45, eq=eq_spin)
    su = lambda m: m * m - 1
    check(NAME, "Pati-Salam su(4) + su(2) + su(2) is 15 + 3 + 3 = 21 and the Standard Model su(3) + su(2) + u(1) is 8 + 3 + 1 = 12",
          su(4) + su(2) + su(2) == 21 and su(3) + su(2) + 1 == 12, eq=eq_spin)
    check(NAME, "the Clifford algebra on 14 generators has dimension 2^14 = 16384 = 128^2 and the half-spinor of so(14) is 64 = 128/2",
          2 ** 14 == 128 ** 2 and half(14) == 128 // 2, eq=eq_spin)

    # ---- the magic square -------------------------------------------------------------------
    tri = {"R": 0, "C": 2, "H": 9, "O": 28}
    dim = {"R": 1, "C": 2, "H": 4, "O": 8}
    L = {(a, b): tri[a] + tri[b] + 3 * dim[a] * dim[b] for a in tri for b in tri}
    eq_magic = register(NAME, Equation(
        key="cl_magic_square_table", entry=ENTRY, title="The Freudenthal-Tits magic square",
        source_kind="repo-doc", source=f"{THEORY} section 7 and lisi/check_magic_square.py",
        latex_str=(r"\begin{array}{c|cccc} & \mathbb R & \mathbb C & \mathbb H & \mathbb O\\\hline "
                   r"\mathbb R & 3 & 8 & 21 & 52\\ \mathbb C & 8 & 16 & 35 & 78\\ \mathbb H & 21 & 35 & 66 & 133\\ "
                   r"\mathbb O & 52 & 78 & 133 & 248\end{array}"),
        notes="L(A,B) = tri(A) + tri(B) + 3 dim A dim B with tri = 0, 2, 9, 28. The repository computes tri(A) as an "
              "exact nullspace; this module takes the four dimensions and checks they agree with the classical "
              "algebras (so(8), su(2)^3, u(1)^2). The corner is E8 = 28 + 28 + 3.64.",
    ))
    rows = [[L[(a, b)] for b in "RCHO"] for a in "RCHO"]
    check(NAME, "the formula reproduces the table 3 8 21 52 / 8 16 35 78 / 21 35 66 133 / 52 78 133 248",
          rows == [[3, 8, 21, 52], [8, 16, 35, 78], [21, 35, 66, 133], [52, 78, 133, 248]], eq=eq_magic)
    check(NAME, "tri(O) = so(8) = 28, tri(H) = su(2)^3 = 9, tri(C) = u(1)^2 = 2: the classical dimensions agree",
          (tri["O"], tri["H"], tri["C"]) == (so_dim(8), 3 * su(2), 2), eq=eq_magic)
    check(NAME, "the corner 28 + 28 + 192 = 248 is the root count of the D4 x D4 blocks (24 + 24 + 192) plus the rank 8",
          L[("O", "O")] == 248 == len(d4_roots) + len(d4_perp) + len(mixed) + N, eq=eq_magic)

    # ---- the 16 of SO(10), anomaly sums, Q = T3 + Y -------------------------------------------
    F = Fraction
    sixteen = [   # name, colour dim, isospin dim, hypercharge Y = q/6
        ("Q_L", 3, 2, F(1, 6)), ("u^c", 3, 1, F(-2, 3)), ("d^c", 3, 1, F(1, 3)),
        ("L", 1, 2, F(-1, 2)), ("e^c", 1, 1, F(1)), ("nu^c", 1, 1, F(0))]
    states = []
    for name, c, w, y in sixteen:
        t3s = [F(0)] if w == 1 else [F(1, 2), F(-1, 2)]
        states += [(name, y, t3) for _ in range(c) for t3 in t3s]
    eq_16 = register(NAME, Equation(
        key="cl_sm_sixteen", entry=ENTRY, title="The 16: anomaly sums and the electric charge",
        source_kind="repo-doc", source=f"{THEORY} section 2 (sm_charges.py)",
        latex_str=(r"\sum Y=\sum Y^{3}=\sum Q=\sum Q^{3}=0\ \text{over the }16,\qquad Q=T_{3}+Y"),
        notes="One generation of Q_L, u^c, d^c, L, e^c and nu^c, with Y = q/6. The sums are exact rational arithmetic "
              "over the 16 states. The mixed sums SU(2)^2 Y and SU(3)^2 Y vanish too, and the charges are 0, "
              "+-1/3, +-2/3, +-1.",
    ))
    check(NAME, "the 16 states: 6 + 3 + 3 + 2 + 1 + 1", len(states) == 16, eq=eq_16)
    check(NAME, "the sums of Y and of Y^3 over the 16 vanish (the gravitational and cubic hypercharge anomalies)",
          sum(s[1] for s in states) == 0 and sum(s[1] ** 3 for s in states) == 0, eq=eq_16)
    charges = [s[1] + s[2] for s in states]
    check(NAME, "with Q = T3 + Y the sums of Q and of Q^3 vanish and the charges are 0, +-1/3, +-2/3, +-1",
          sum(charges) == 0 and sum(q ** 3 for q in charges) == 0
          and set(charges) == {F(0), F(1, 3), F(-1, 3), F(2, 3), F(-2, 3), F(1), F(-1)}, eq=eq_16)
    su2_y = sum(c * y for _n, c, w, y in sixteen if w == 2)           # doublets, weighted by colour
    su3_y = sum(w * y for _n, c, w, y in sixteen if c == 3)           # triplets, weighted by isospin
    check(NAME, "the mixed anomalies vanish: sum over doublets of Y = 0 (SU(2)^2 Y) and over triplets of Y = 0 (SU(3)^2 Y)",
          su2_y == 0 and su3_y == 0, eq=eq_16)

    # ---- hypercharge normalisation and the weak mixing angle ----------------------------------
    k_norm = sum(s[1] ** 2 for s in states) / sum(s[2] ** 2 for s in states)
    s2 = sp.Symbol("s_W")
    kk = sp.Symbol("k")
    eq_mix = register(NAME, Equation(
        key="cl_weak_mixing", entry=ENTRY, title="The weak mixing angle at unification",
        source_kind="derived", source="g^2 sum T3^2 = g'^2 sum Y^2 at the unification scale; sum over the 16",
        expr=sp.Eq(s2, 1 / (1 + kk)), latex_str=None,
        notes="s_W = sin^2 theta_W = g'^2 / (g^2 + g'^2). If the two traces carry equal couplings at the unification "
              "scale then g'^2 / g^2 = 1 / k with k the ratio of the traces, and s_W = 1 / (1 + k). For one "
              "generation k = 5/3 and s_W = 3/8, the SU(5) value at the unification scale; it runs to about 0.23 "
              "at low energy, which this module does not compute.",
    ))
    check(NAME, "k = sum Y^2 / sum T3^2 over the 16 is exactly 5/3", k_norm == F(5, 3), eq=eq_mix)
    check(NAME, "s_W = 1 / (1 + 5/3) = 3/8", sp.Rational(1) / (1 + sp.Rational(5, 3)) == sp.Rational(3, 8), eq=eq_mix)
    check(NAME, "the relation inverts to g'^2/g^2 = 1/k: s_W = (1/k) / (1 + 1/k)",
          sp.simplify(1 / (1 + kk) - (1 / kk) / (1 + 1 / kk)) == 0, eq=eq_mix)

    # ---- the Distler-Garibaldi dimension count ---------------------------------------------------
    W, ngen, fpg = sp.symbols("W N_gen F")
    eq_dg = register(NAME, Equation(
        key="cl_dg_dimension_bound", entry=ENTRY, title="The dimension no-go for three generations",
        source_kind="repo-doc", source=f"{DG} section 3, as {THEORY} section 8 states it",
        expr=sp.Eq(W, 4 * fpg * ngen), latex_str=None,
        notes="W counts the linearly independent weight vectors the fermions need: 2 x 2 x 45 = 180 for three "
              "generations of 15 Weyl fields. The factor 2 x 2 is as THEORY.md states it; the paper was not re-read "
              "here. The eigenspaces of the two involutions have dimensions 128 and 112 (checked above from the "
              "roots), so 180 does not fit. Lisi's model says its generations are not linearly independent, which "
              "the repository does not examine.",
    ))
    need = 4 * 15 * 3
    check(NAME, "three generations of 15 Weyl fields need 4 x 45 = 180 weight vectors", need == 180, eq=eq_dg)
    check(NAME, "the larger eigenspace of an inner involution is 128 (from the roots), and 180 > 128",
          max(par1[1], par8[1]) == 128 and need > 128, eq=eq_dg)
    check(NAME, "at most 128 / 60 = 2.13 generations would fit, so the count allows two and not three",
          sp.floor(sp.Rational(128, 60)) == 2, eq=eq_dg)


def _zero_sum_in_e6(w27, t):
    """True if the E6 components of three weights sum to zero, i.e. the three roots sum to a root-lattice
    element of the A2 plane (their total pairing with the A2 is 3 times one weight, and the total has no
    other component)."""
    total = [sum(w27[i][k] for i in t) for k in range(N)]
    # the total must lie in the span of the A2: only coefficients 2 and 3 (alpha_3, alpha_4) may be nonzero
    return all(total[k] == 0 for k in range(N) if k not in (2, 3))
