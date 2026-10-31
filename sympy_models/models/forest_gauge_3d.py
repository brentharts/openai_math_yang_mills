"""Arithmetic of T. Romanov, "Three-dimensional lattice Yang-Mills theory in a forest gauge" (Zenodo 10.5281/zenodo.23783316 as supplied;
attached file 2610.0069v1.pdf). A companion to the su2 lattice models: it checks the finite counts and the arithmetic of that paper,
nothing else.

What is checked (every number is asserted):

  forest3d_counts    Proposition 2.2 and Lemma 2.1, recomputed here independently in Python from the string
                     210102120102120210120210102 (lexicographic over {0,1,2}^3, first coordinate most significant, F = {(x - e_mu(x), mu(x))}):
                     no cycle in a box; per period cell 81 plaquettes (54 with one link of F, 27 with two), 54 free links; every free link
                     in four plaquettes with exactly seven partners (and sum over plaquettes through it of the free links a_p equal to 11,
                     so 11 - 4 = 7); eight partners in the temporal gauge mu = 0 and twelve without gauge fixing; sums over the cell
                     (sum a_p = 4 x 54, sum a_p (a_p - 1) = 7 x 54) and the integrality identity k(k-1) - (4k-6) = (k-2)(k-3) >= 0 of
                     Proposition 4.1, with the bound 16 - 6P/a = 7 at P/a = 3/2, against 20/3 from convexity.
  forest3d_ceilings  The ceilings beta_W < 4/p for p partners: 4/7 (forest), 1/2 (temporal), 1/3 (no gauge fixing), 2/9 in d = 4 with 18
                     partners, and the rule 6(d-1) partners giving 2/(3(d-1)); the gap bound m >= -log(7 beta_W/4) and its value at sample points.
  forest3d_table     Table 1 from beta_W = N^2 beta and |beta| < kappa/(N h'), kappa_SU(N) = N/2, i.e. beta_W < N^2/(2 h'), for
                     h' = 16, 12, 11, 9.599991, 9.3479; the stated 1/19.199982 and 1/18.6959; 9.599991 = 9599991/1000000; the ranges
                     0.2083 and 0.4687; "from 1/3 to 4/7" and "from 1/24 to 1/19.2" of the discussion; SO(N): kappa = (N-2)/4.
  forest3d_bloch     Proposition 1.3 and its proof: the 81 x 54 integer matrix B of the coboundary on free links at the Bloch momentum
                     (pi, pi, pi) of the period-3 cell is rebuilt; the characteristic polynomial of B^T B is factored over Q and compared with the
                     statement (factors (x-6)^2 (x-4)^2, two squared octics, two nonics); the nonic of the statement is one of them; its
                     largest real root is 9.347955586999... and is the largest eigenvalue of B^T B.

Not tested: every theorem. The Dobrushin/Kantorovich estimates (the sharp quarter), the passage from a Dobrushin bound to a state and a
gap, reflection positivity, the certified Hessian bound h = 9.599991 (semidefinite certificate, outward-rounded interval arithmetic,
26,371,504 squares), the Bakry-Emery and log-Sobolev conclusions, the decay rate m_G, the lower bound proof beyond the matrix
computation above, and every statement about existence or a gap. Nothing here bears on Yang-Mills in the continuum.

Discrepancies found: see the notes of forest3d_table (a rounding of the SU(3) entry).
"""
import itertools
import sys
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register

NAME = "forest_gauge_3d"
ENTRY = "forest-gauge"
MU_STRING = "210102120102120210120210102"
PLANES = ((0, 1), (0, 2), (1, 2))


# ---------------------------------------------------------------------------------------------------------------
def e(a):
    return tuple(1 if i == a else 0 for i in range(3))


def add(x, y):
    return tuple(x[i] + y[i] for i in range(3))


MU = {xyz: int(MU_STRING[i]) for i, xyz in enumerate(itertools.product(range(3), repeat=3))}


def mu_forest(x):
    return MU[tuple(c % 3 for c in x)]


def mu_temporal(x):
    return 0


def mu_none(x):
    return None


def is_fixed(mu, y, a):
    """The link (y, a) is in F_mu iff mu(y + e_a) = a (the link into the site y + e_a along the direction mu)."""
    return mu is not mu_none and mu(add(y, e(a))) == a


def plaquette(x, a, b):
    return [(x, a), (add(x, e(a)), b), (add(x, e(b)), a), (x, b)]


def census(mu, window=9):
    """Counts for the links of the central cell, computed on the unfolded lattice (all plaquettes meeting them)."""
    cell = list(itertools.product(range(3), repeat=3))
    sites = list(itertools.product(range(-1, window - 1), repeat=3))
    free_links = [(x, a) for x in cell for a in range(3) if not is_fixed(mu, x, a)]
    through = {l: [] for l in free_links}
    for x in sites:
        for a, b in PLANES:
            p = plaquette(x, a, b)
            for l in p:
                if l in through:
                    through[l].append(p)
    out = {}
    for l, plist in through.items():
        partners = set()
        ap_sum = 0
        for p in plist:
            fr = [q for q in p if not is_fixed(mu, *q)]
            ap_sum += len(fr)
            partners.update(q for q in fr if q != l)
        out[l] = (len(plist), len(partners), ap_sum)
    return free_links, out


def has_cycle(links):
    parent = {}

    def find(u):
        while parent.setdefault(u, u) != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u
    for (y, a) in links:
        u, v = find(y), find(add(y, e(a)))
        if u == v:
            return True
        parent[u] = v
    return False


def bloch_matrix():
    """B: 81 plaquettes x 54 free links of the cell, (d xi)_p = xi_e1 + xi_e2 - xi_e3 - xi_e4 with the Bloch factor (-1)^(sum of the
    cell translation) for a link outside the base cell, at momentum (pi, pi, pi) of the period-3 cell."""
    cell = list(itertools.product(range(3), repeat=3))
    free = [(x, a) for x in cell for a in range(3) if not is_fixed(mu_forest, x, a)]
    col = {l: i for i, l in enumerate(free)}
    rows = []
    for x in cell:
        for a, b in PLANES:
            r = [0] * len(free)
            for sgn, (y, c) in zip((1, 1, -1, -1), plaquette(x, a, b)):
                if is_fixed(mu_forest, y, c):
                    continue
                base = tuple(t % 3 for t in y)
                shift = tuple(t // 3 for t in y)
                r[col[(base, c)]] += sgn * (-1) ** sum(shift)
            rows.append(r)
    return sp.Matrix(rows)


# ---------------------------------------------------------------------------------------------------------------
def build_and_verify():
    # =====================================================================================================
    # Proposition 2.2 and Lemma 2.1
    # =====================================================================================================
    cell = list(itertools.product(range(3), repeat=3))
    n_pl = len(cell) * len(PLANES)
    by_fixed = {}
    for x in cell:
        for a, b in PLANES:
            k = sum(is_fixed(mu_forest, *l) for l in plaquette(x, a, b))
            by_fixed[k] = by_fixed.get(k, 0) + 1
    free, cen = census(mu_forest)
    _, cen_t = census(mu_temporal)
    _, cen_n = census(mu_none)
    box = range(9)
    box_links = [((x, y, z), a) for x in box for y in box for z in box for a in range(3)
                 if max(((x, y, z)[i] + (1 if i == a else 0)) for i in range(3)) <= 8 and is_fixed(mu_forest, (x, y, z), a)]
    mu_ok = len(MU_STRING) == 27 and set(MU_STRING) == {"0", "1", "2"}
    every_site_head = all(sum(is_fixed(mu_forest, add(x, tuple(-t for t in e(a))), a) and mu_forest(x) == a for a in range(3)) == 1 for x in cell)

    eq_cnt = register(NAME, Equation(
        key="forest3d_counts", entry=ENTRY, title="Forest gauge on Z^3: the counts of Proposition 2.2",
        source_kind="derived", source="derived here by an independent enumeration; the statements are Proposition 2.2 and Lemma 2.1 of Romanov's paper",
        latex_str=(r"\text{per period cell:}\ \ 81\ \text{plaquettes}=54\,(\text{one link of }F)+27\,(\text{two}),\qquad 54\ \text{free links}"
                   r"%%BR%%\text{partners of a free link:}\ \ 7\ (F),\quad 8\ (\mu\equiv0),\quad 12\ (\text{no gauge fixing}),\qquad "
                   r"\sum_{p\ni e}a_p=11\ \Rightarrow\ \sum_{p\ni e}(a_p-1)=7"
                   r"%%BR%%k(k-1)-(4k-6)=(k-2)(k-3)\ge0,\qquad 16-6\tfrac{P}{a}\Big|_{P/a=3/2}=7\ \ \bigl(\text{convexity alone: }\tfrac{20}{3}\bigr)"),
        notes="Verifies the finite counts of the paper by an independent enumeration (Python, not the paper's scripts and not its Lean file): the "
              "string is read lexicographically over {0,1,2}^3 with the first coordinate most significant, F = {(x - e_mu(x), mu(x))}, and a link "
              "(y, a) is free iff mu(y + e_a) != a. Partners are free links that share a plaquette, counted on the unfolded lattice. It does not test "
              "the theorems that use the counts (the Dobrushin argument, its ceiling, the Hessian certificate).",
    ))
    check(NAME, "the string 210102120102120210120210102 has 27 entries in {0,1,2}, and every site of the cell is the head of exactly one link of F",
          mu_ok and every_site_head, eq=eq_cnt)
    check(NAME, "Lemma 2.1: the links of F inside a 9 x 9 x 9 box form a forest (union-find finds no cycle)", not has_cycle(box_links), eq=eq_cnt)
    check(NAME, f"a period cell has {n_pl} plaquettes, {by_fixed.get(1, 0)} with one link of F and {by_fixed.get(2, 0)} with two, and none with zero or "
                f"three; {len(free)} free links of 81",
          n_pl == 81 and by_fixed == {1: 54, 2: 27} and len(free) == 54, eq=eq_cnt)
    check(NAME, "every one of the 54 free links of the cell lies in 4 plaquettes and has exactly 7 partners; the sum over those plaquettes of the number a_p of "
                "free links is 11 for each (so the partners are 11 - 4 = 7)",
          all(v == (4, 7, 11) for v in cen.values()) and len(cen) == 54, eq=eq_cnt)
    check(NAME, "temporal gauge mu = 0: every free link (2 of 3 per site, 54 per cell) has 8 partners in 4 plaquettes; without gauge fixing every link "
                "(81 per cell) has 12 partners in 4 plaquettes",
          len(cen_t) == 54 and all(v[:2] == (4, 8) for v in cen_t.values()) and len(cen_n) == 81 and all(v[:2] == (4, 12) for v in cen_n.values()),
          eq=eq_cnt)
    sum_ap = sum(len([q for q in plaquette(x, a, b) if not is_fixed(mu_forest, *q)]) for x in cell for a, b in PLANES)
    sum_ap2 = sum(len([q for q in plaquette(x, a, b) if not is_fixed(mu_forest, *q)]) * (len([q for q in plaquette(x, a, b) if not is_fixed(mu_forest, *q)]) - 1)
                  for x in cell for a, b in PLANES)
    check(NAME, f"over a cell sum a_p = {sum_ap} = 4 x 54 and sum a_p (a_p - 1) = {sum_ap2} = 7 x 54 (a_p = 3 on the 54 plaquettes with one fixed link, 2 on the 27 "
                "with two)", sum_ap == 4 * 54 and sum_ap2 == 7 * 54 and 54 * 3 * 2 + 27 * 2 * 1 == sum_ap2, eq=eq_cnt)
    k = sp.Symbol("k")
    check(NAME, "k(k-1) - (4k-6) = (k-2)(k-3) as a polynomial identity, and it is >= 0 for every integer k (the integrality step of Proposition 4.1)",
          sp.expand(k * (k - 1) - (4 * k - 6) - (k - 2) * (k - 3)) == 0 and all((j - 2) * (j - 3) >= 0 for j in range(-50, 51)), eq=eq_cnt)
    check(NAME, "with sum_p a_p = 4a and P = 3a/2 plaquettes per free link, |beta| (1/a) sum a_p(a_p - 1) >= |beta| (16 - 6P/a) a/a = 7 |beta|; convexity (mean "
                "a_p = 8/3) alone gives (3/2)(8/3)(5/3) = 20/3 < 7 (the paper's remark)",
          16 - 6 * Fr(3, 2) == 7 and Fr(3, 2) * Fr(8, 3) * (Fr(8, 3) - 1) == Fr(20, 3), eq=eq_cnt)

    # =====================================================================================================
    # ceilings and the gap formula
    # =====================================================================================================
    ceil = lambda p: Fr(4, p)
    eq_ceil = register(NAME, Equation(
        key="forest3d_ceilings", entry=ENTRY, title="Dobrushin ceilings 4/p and the gap formula of the forest-gauge paper",
        source_kind="derived", source="derived here from the paper's statements (Theorem 1.1, Remark 4.2, the introduction)",
        latex_str=(r"\beta_W<\frac{4}{p}:\quad p=7:\ \tfrac47,\qquad p=8:\ \tfrac12,\qquad p=12:\ \tfrac13,\qquad p=18\ (d=4):\ \tfrac29,\qquad "
                   r"p=6(d-1)\ \Rightarrow\ \frac{2}{3(d-1)}"
                   r"%%BR%%m\ \ge\ -\log\frac{7\beta_W}{4},\qquad \alpha=7|\beta|=\tfrac74|\beta_W|\quad(\beta_W=4\beta\ \text{for}\ SU(2))"),
        notes="Verifies the arithmetic that turns a partner count p into the range 4/p (alpha = p |beta| < 1 with beta_W = N^2 beta = 4 beta for SU(2)) and the "
              "sample values of the gap bound. It does not test that the Dobrushin condition holds, that the Kantorovich coefficient is |beta| "
              "(the sharp quarter), or that a gap follows.",
    ))
    check(NAME, "p partners give alpha = p|beta| = p beta_W/4 < 1, i.e. beta_W < 4/p: 4/7, 1/2, 1/3 for p = 7, 8, 12 (forest, temporal gauge, no gauge fixing)",
          [ceil(7), ceil(8), ceil(12)] == [Fr(4, 7), Fr(1, 2), Fr(1, 3)], eq=eq_ceil)
    check(NAME, "d = 4 without gauge fixing: 18 partners give 4/18 = 2/9; in d dimensions 6(d-1) partners give 2/(3(d-1)): 1/3 at d = 3 and 2/9 at d = 4",
          ceil(18) == Fr(2, 9) and all(ceil(6 * (d - 1)) == Fr(2, 3 * (d - 1)) for d in range(2, 12)) and Fr(2, 3 * 2) == Fr(1, 3) and 6 * 2 == 12 and 6 * 3 == 18,
          eq=eq_ceil)
    sample = [Fr(1, 10), Fr(1, 4), Fr(5, 14), Fr(1, 2)]
    mvals = [-mp.log(7 * mp.mpf(b.numerator) / b.denominator / 4) for b in sample]
    check(NAME, "the gap bound -log(7 beta_W/4) is positive exactly for beta_W < 4/7: " +
          ", ".join(f"beta_W = {b}: {mp.nstr(v, 6)}" for b, v in zip(sample, mvals)) + "; it is 0 at beta_W = 4/7 and decreases",
          all(v > 0 for v in mvals) and abs(-mp.log(7 * mp.mpf(4) / 7 / 4)) < 1e-15 and mvals[0] > mvals[1] > mvals[2] > mvals[3], eq=eq_ceil)
    check(NAME, "the discussion: the SU(2) range improves from 1/3 to 4/7 (ratio 12/7) and the curvature range from 1/24 to 1/19.2 (1/(2h) with h = 9.6 = 48/5 gives 5/96)",
          Fr(4, 7) / Fr(1, 3) == Fr(12, 7) and Fr(1, 2) / Fr(48, 5) == Fr(5, 96) and abs(float(Fr(5, 96)) - 1 / 19.2) < 1e-15, eq=eq_ceil)

    # =====================================================================================================
    # Table 1
    # =====================================================================================================
    h = Fr(9599991, 1000000)
    hF = mp.mpf("9.347955586999")
    printed = {  # h' -> (SU(2), SU(3)) as printed in Table 1
        "16": (0.125, 0.281), "12": (0.167, 0.375), "11": (0.182, 0.409), "9.599991": (0.208, 0.468), "9.3479": (0.214, 0.481)}
    exact = {}
    for hp in ("16", "12", "11", "9.599991", "9.3479"):
        hv = Fr(hp)
        exact[hp] = tuple(Fr(N * N, 2) / hv for N in (2, 3))                           # beta_W < N^2 kappa/(N h') = N^2/(2 h')
    rounded_ok, truncated = [], []
    for hp, (a2, a3) in exact.items():
        p2, p3 = printed[hp]
        for ex, pr in ((a2, p2), (a3, p3)):
            r3 = round(float(ex), 3)
            if abs(r3 - pr) < 1e-12:
                rounded_ok.append((hp, pr))
            elif abs(float(ex) // 0.001 * 0.001 - pr) < 1e-9:
                truncated.append((hp, pr, float(ex)))
            else:
                raise AssertionError(f"Table 1 entry for h' = {hp}: printed {pr}, computed {float(ex)}")
    eq_tab = register(NAME, Equation(
        key="forest3d_table", entry=ENTRY, title="Table 1 of the forest-gauge paper: strong-coupling ranges from beta_W = N^2 beta",
        source_kind="derived", source="derived here from the paper's Table 1 and Theorem 1.2",
        latex_str=(r"|\beta|<\frac{\kappa_G}{N h'},\ \ \kappa_{SU(N)}=\tfrac N2,\ \ \beta_W=N^{2}\beta\ \Rightarrow\ \beta_W<\frac{N^{2}}{2h'}"
                   r"%%BR%%h'=16,\,12,\,11,\,9.599991,\,9.3479:\quad SU(2):\ 0.125,\,0.167,\,0.182,\,0.2083,\,0.2139;\quad SU(3):\ 0.28125,\,0.375,\,0.4091,\,0.46875,\,0.4814"
                   r"%%BR%%\tfrac{1}{2h}=\tfrac{1}{19.199982},\qquad \tfrac{1}{2h_F}=\tfrac{1}{18.6959},\qquad h=\tfrac{9599991}{1000000},\ \ h_F=9.347955587\ldots"),
        notes="Every Table 1 entry that rests on a curvature constant is recomputed as N^2/(2h') in exact rationals. Nine of the ten printed entries are the "
              "value rounded to three decimals. One is not: for SU(3) at h = 9.599991 the table prints 0.468 and the text of Theorem 1.2 prints "
              "0.4687, while N^2/(2h) = 4.5/9.599991 = 0.468750(4), which rounds to 0.469 (three decimals) and 0.4688 (four): the printed values are truncations "
              "of 0.46875. The SU(3) condition itself is correct and slightly better than printed, so this is a rounding discrepancy, not an error of the "
              "method. The Dobrushin rows (0.333, 0.5, 0.571 = 4/7) are the ceilings above. Not tested: that any of these constants is valid (Theorem 1.2 "
              "and the Bakry-Emery inputs are quoted).",
    ))
    check(NAME, "h = 9599991/1000000 = 9.599991 exactly, and 2h = 19.199982 (the stated |beta| < 1/19.199982)", h == Fr("9.599991") and 2 * h == Fr("19.199982"), eq=eq_tab)
    check(NAME, "Table 1, SU(2) and SU(3), computed as N^2/(2h') from h' = 16, 12, 11, 9.599991, 9.3479: " +
          "; ".join(f"h'={hp}: {float(a):.5f}, {float(b):.5f}" for hp, (a, b) in exact.items()) +
          f". Of the ten printed entries {len(rounded_ok)} equal the three-decimal rounding and {len(truncated)} is a truncation "
          f"(SU(3), h = 9.599991: printed 0.468, exact 0.468750(4))",
          len(rounded_ok) == 9 and truncated == [("9.599991", 0.468, float(exact["9.599991"][1]))] and abs(float(exact["9.599991"][1]) - 0.46875) < 1e-6,
          eq=eq_tab)
    check(NAME, "the text of Theorem 1.2: beta_W < 0.2083 (SU(2)) is the four-decimal rounding of 2/h = 0.208333; the SU(3) value printed 0.4687 is 0.46875(04) truncated "
                "(rounds to 0.4688): discrepancy recorded in the notes",
          round(float(exact["9.599991"][0]), 4) == 0.2083 and round(float(exact["9.599991"][1]), 4) == 0.4688 and int(float(exact["9.599991"][1]) * 1e4) / 1e4 == 0.4687,
          eq=eq_tab)
    check(NAME, "Proposition 1.3 corollary: 2 h_F = 18.69591117 (printed 1/18.6959 after rounding); no curvature argument in this gauge goes beyond |beta| < 1/18.6959; "
                "the optimal constant lies in [9.3479, 9.6000]",
          abs(2 * hF - mp.mpf("18.695911173998")) < 1e-9 and abs(round(float(2 * hF), 4) - 18.6959) < 1e-12 and 9.3479 < hF < 9.6, eq=eq_tab)
    check(NAME, "SO(N): kappa = (N-2)/4 gives |beta| < (N-2)/(4 N h) (Remark 5.3), e.g. N = 3: 1/(12 h) and kappa_SO(3) = kappa_SU(2)/... = 1/4 against kappa_SU(2) = 1",
          Fr(3 - 2, 4) == Fr(1, 4) and Fr(1, 4 * 3) / h == Fr(1, 12) / h, eq=eq_tab)

    # =====================================================================================================
    # Proposition 1.3: the Bloch matrix and the nonic
    # =====================================================================================================
    B = bloch_matrix()
    G = B.T * B
    x = sp.Symbol("x")
    from sympy.polys.matrices import DomainMatrix
    cp = DomainMatrix.from_Matrix(G).to_field().charpoly() if False else DomainMatrix.from_Matrix(G).charpoly()
    P = sp.Poly(sp.Poly(cp, x).as_expr(), x)
    fac = sp.factor_list(P.as_expr())
    shape = sorted((sp.degree(f, x), m) for f, m in fac[1])
    paper = x ** 9 - 35 * x ** 8 + 505 * x ** 7 - 3901 * x ** 6 + 17538 * x ** 5 - 46770 * x ** 4 + 72344 * x ** 3 - 60490 * x ** 2 + 22788 * x - 1944
    in_list = any(sp.expand(f - paper) == 0 for f, m in fac[1])
    roots = sp.Poly(paper, x).nroots(n=30, maxsteps=200)
    real_roots = sorted(float(r.as_real_imag()[0]) for r in roots if abs(sp.im(r)) < 1e-20)
    top = mp.polyroots([int(c) for c in sp.Poly(paper, x).all_coeffs()], maxsteps=500, extraprec=200)
    top_real = max(mp.re(r) for r in top if abs(mp.im(r)) < mp.mpf("1e-25"))
    lam_max = max(np.linalg.eigvalsh(np.array(G.tolist(), dtype=float)))
    eq_blo = register(NAME, Equation(
        key="forest3d_bloch", entry=ENTRY, title="The lower bound 9.3479: the Bloch matrix and its nonic",
        source_kind="derived", source="derived here by rebuilding the matrix; the statement is Proposition 1.3 and its proof (section 6) of Romanov's paper",
        latex_str=(r"B\in\mathbb Z^{81\times54},\quad \chi_{B^{T}B}=(x-6)^{2}(x-4)^{2}\,q_8^{2}\,q_8'^{2}\,q_9\,q_9',\qquad "
                   r"q_9=x^{9}-35x^{8}+505x^{7}-3901x^{6}+17538x^{5}-46770x^{4}+72344x^{3}-60490x^{2}+22788x-1944"
                   r"%%BR%%\lambda_{\max}(B^{T}B)=9.347955586999\ldots=\text{largest root of }q_9,\qquad h_F\ge9.3479"),
        notes="Rebuilds the 81 x 54 integer matrix of (d xi)_p on the free links of the period cell with Bloch factor (-1)^(cell translation) at momentum "
              "(pi, pi, pi), factors the characteristic polynomial of B^T B over Q, and compares with the statement. The reading of 'two squared octics and "
              "two nonics' is: squares of two octics, and two nonics (2 + 2 + 16 + 16 + 9 + 9 = 54 when the squares are counted: (x-6)^2(x-4)^2 has degree 4, "
              "two octics squared 32, two nonics 18; total 54). The claim that the supremum of ||d xi||^2/||xi||^2 over free fields is >= this eigenvalue "
              "(Bloch waves on large boxes) and that both signs of the Hessian reach it (the J-flux field) are not tested here.",
    ))
    check(NAME, "B is an 81 x 54 integer matrix (entries in {-2..2}) and B^T B is 54 x 54 symmetric", B.shape == (81, 54) and all(v == int(v) for v in B)
          and G == G.T and max(abs(v) for v in B) <= 2, eq=eq_blo)
    check(NAME, f"the characteristic polynomial of B^T B has degree 54 and factors over Q with degrees (multiplicities) {shape}; the statement's shape is "
                "(x-6)^2 (x-4)^2, two octics squared and two nonics",
          sp.degree(P.as_expr(), x) == 54 and shape == [(1, 2), (1, 2), (8, 2), (8, 2), (9, 1), (9, 1)], eq=eq_blo)
    check(NAME, "the factors of degree one are x - 6 and x - 4, each squared", sorted([f for f, m in fac[1] if sp.degree(f, x) == 1 and m == 2], key=str) in
          ([x - 6, x - 4], [x - 4, x - 6]), eq=eq_blo)
    check(NAME, "the nonic printed in Proposition 1.3 is one of the two nonic factors of the characteristic polynomial, coefficient by coefficient", in_list, eq=eq_blo)
    check(NAME, f"its largest real root is {mp.nstr(top_real, 13)} (mpmath, 200 extra digits); the paper prints 9.347955586999..., agreement to 1e-12; it is the largest root "
                f"of the polynomial and the largest eigenvalue of B^T B (numpy {lam_max:.12f})",
          abs(top_real - mp.mpf("9.347955586999")) < mp.mpf("1e-12") and abs(top_real - mp.mpf(lam_max)) < mp.mpf("1e-9")
          and top_real == max(mp.re(r) for r in top) and abs(real_roots[-1] - float(top_real)) < 1e-12, eq=eq_blo)
    check(NAME, "the real roots of the nonic are all positive, and 9.3479 < 9.347955587 < 9.6 (so the lower bound lies below the certified h = 9.599991)",
          all(r > 0 for r in real_roots) and mp.mpf("9.3479") < top_real < mp.mpf(h.numerator) / h.denominator, eq=eq_blo)
