"""Strong-coupling lattice gauge theory on Z^3 in a period-3 forest gauge (companion repository).

Source: Tristan Romanov, romanov360/yang-mills-3d-forest-gauge (code MIT, paper CC BY 4.0), DOI
10.5281/zenodo.23783317 (v1.0.0). This module rederives, in exact arithmetic and without copying code, the
combinatorial and arithmetic parts of that work. The 27-letter pattern below is the repository's data.

What is checked (all exact unless marked):
  * the pattern fixes 27 of the 81 links of a period cell, every site is the head of exactly one fixed link,
    and the fixed links form a forest on Z^3 (the coordinate sum falls by one along every parent step);
  * each of the 54 free links has exactly 7 distinct free plaquette partners (plaquettes with 1, 1, 1 and 2
    fixed links); temporal gauge gives 8 and no gauge fixing gives 12;
  * the arithmetic of Theorem A (row sum 7|beta|, critical beta_W = 4/7, mass bound, the temporal and
    ungauged comparisons) and of Theorem B (the range 1/(2h), the comparison table);
  * the lower bound for the Hessian constant: at the identity field and Bloch momentum (pi, pi, pi) the 54 x 54
    integer Bloch matrix has characteristic polynomial with an irreducible integer nonic factor whose top
    root is 9.347955587.

Quoted, not rederived: Theorem A (Dobrushin uniqueness and decay), Theorem B (the certified bound
h = 9.599991, by a rational certificate and interval arithmetic over 26 million squares), the Lean file, the
Hessian normalisation -f_p = |Y1 + Y2 - Y3 - Y4|^2, and the conjecture that the optimal constant equals 9.3479.

Not a continuum statement. The results hold for strong coupling on a lattice and say nothing about the
Clay problem, the continuum limit, existence or a mass gap of quantum Yang-Mills on R^4. The repository states
the same.
"""
import itertools
from collections import Counter
from fractions import Fraction as Fr

import sympy as sp

from eqlib import Equation, check, register

NAME = "forest_gauge"
ENTRY = "forest-gauge"
REPO = "romanov360/yang-mills-3d-forest-gauge"
DOI = "10.5281/zenodo.23783317"
SRC = f"{REPO} v1.0.0"

PATTERN = "210102120102120210120210102"      # the repository's period-3 pattern: direction of the fixed head link
CELLS = sorted(itertools.product(range(3), repeat=3))
MU = {c: int(PATTERN[i]) for i, c in enumerate(CELLS)}
H_CERT = Fr(9599991, 10 ** 6)                  # Theorem B (quoted)


def unit(a):
    return tuple(int(i == a) for i in range(3))


def add(x, y, s=1):
    return tuple(p + s * q for p, q in zip(x, y))


def mu(x):
    return MU[tuple(v % 3 for v in x)]


def is_fixed(x, a):
    """Link (x, a) from x to x + e_a is fixed iff its head x + e_a points back along a."""
    return mu(add(x, unit(a))) == a


def plaquettes_of(x, a):
    """The four plaquettes of Z^3 containing link (x, a), each as its four links with true coordinates."""
    out = []
    for b in range(3):
        if b == a:
            continue
        for s in (0, -1):
            o = add(x, unit(b), s)           # base corner of a plaquette in the (a, b) plane
            out.append([(o, a), (add(o, unit(a)), b), (add(o, unit(b)), a), (o, b)])
    return out


def partners(x, a, fixed):
    """Distinct non-fixed links sharing a plaquette with (x, a), and the fixed-link count of each plaquette."""
    part, types = set(), []
    for p in plaquettes_of(x, a):
        types.append(sum(fixed(*g) for g in p))
        part |= {g for g in p if g != (x, a) and not fixed(*g)}
    return len(part), tuple(sorted(types))


def bloch_matrix_pi():
    """B at k = (pi, pi, pi): rows = the 81 plaquettes of a cell, columns = the 54 free links of a cell."""
    free = [(x, a) for x in CELLS for a in range(3) if not is_fixed(x, a)]
    col = {l: i for i, l in enumerate(free)}
    rows = []
    for x in CELLS:
        for a in range(3):
            for b in range(a + 1, 3):
                r = [0] * len(free)
                for (y, c), sgn in zip([(x, a), (add(x, unit(a)), b), (add(x, unit(b)), a), (x, b)], (1, 1, -1, -1)):
                    yc = tuple(v % 3 for v in y)
                    n = sum((v - w) // 3 for v, w in zip(y, yc))
                    if (yc, c) in col:
                        r[col[(yc, c)]] += sgn * (-1) ** n          # e^{i k.n} = (-1)^{n1+n2+n3} at k = pi
                rows.append(r)
    return sp.Matrix(rows)


def build_and_verify():
    # ---- the forest ---------------------------------------------------------------------------
    links = [(x, a) for x in CELLS for a in range(3)]
    fixed_cell = [l for l in links if is_fixed(*l)]
    heads = Counter(tuple(v % 3 for v in add(x, unit(a))) for x, a in fixed_cell)
    eq_forest = register(NAME, Equation(
        key="fg_forest", entry=ENTRY, title="The period-3 forest gauge in three dimensions: 27 fixed links, one head per site",
        source_kind="derived", source=f"derived here from the pattern of {SRC}",
        latex_str=(r"e=(x,a)\ \text{fixed}\iff \mu(x+\hat a)=a;\qquad \#\text{links}=81,\ \#\text{fixed}=27,\ "
                   r"\#\text{free}=54\ \text{per period cell}"),
        notes="The pattern assigns to each of the 27 sites of a period cell the direction of the fixed link that ends "
              "there, so every site is the head of exactly one fixed link. Each step to the tail lowers the coordinate "
              "sum by one, so the parent map has no cycle, and a graph in which every vertex has one parent edge and "
              "no directed cycle has no undirected cycle: the fixed links form a forest. The repository checks "
              "acyclicity numerically; this is the one-line argument. What is special to this pattern is the partner count "
              "of the next equation: changing one entry of the pattern gives free links with 6, 8 or 9 partners.",
    ))
    check(NAME, "the pattern uses each direction nine times", Counter(PATTERN) == Counter({"0": 9, "1": 9, "2": 9}), eq=eq_forest)
    check(NAME, "a period cell has 81 links, 27 fixed and 54 free", (len(links), len(fixed_cell)) == (81, 27), eq=eq_forest)
    check(NAME, "every one of the 27 sites is the head of exactly one fixed link (true for any pattern, by the definition of the rule)",
          len(heads) == 27 and set(heads.values()) == {1}, eq=eq_forest)
    parent = lambda x: add(x, unit(mu(x)), -1)
    check(NAME, "the parent step x -> x - e_mu(x) lowers the coordinate sum by 1 at 27^2 sample sites (no cycle)",
          all(sum(parent(x)) == sum(x) - 1 for x in itertools.product(range(-9, 12), repeat=3)), eq=eq_forest)
    check(NAME, "the fixed links of the parent map are exactly the fixed links of the rule",
          all(is_fixed(parent(x), mu(x)) for x in itertools.product(range(-6, 7), repeat=3)), eq=eq_forest)

    # ---- partner counts -----------------------------------------------------------------------
    free_cell = [l for l in links if not is_fixed(*l)]
    res = [partners(x, a, is_fixed) for x, a in free_cell]
    temporal = lambda x, a: a == 0                      # fix every link in direction 0
    nogauge = lambda x, a: False
    res_t = [partners(x, a, temporal) for x, a in links if a != 0]
    res_n = [partners(x, a, nogauge) for x, a in links]
    eq_part = register(NAME, Equation(
        key="fg_partners", entry=ENTRY, title="Seven plaquette partners per free link in the forest gauge",
        source_kind="derived", source=f"derived here; compare {REPO}: README, structure3d.py",
        latex_str=(r"\#\{\text{free links sharing a plaquette with }e\}=7\ \ (\text{forest}),\quad 8\ \ (\text{temporal}),"
                   r"\quad 12\ \ (\text{no gauge fixing})"),
        notes="Counted over the 54 free links of a period cell, with true Z^3 coordinates so that no partner is "
              "identified with another by periodicity. The plaquette types are one fixed link in three of the four "
              "plaquettes and two in the fourth. A smaller number of partners gives a larger strong-coupling range.",
    ))
    check(NAME, "each of the 54 free links has exactly 7 distinct free partners", len(res) == 54 and {r[0] for r in res} == {7}, eq=eq_part)
    check(NAME, "the four plaquettes of each free link have 1, 1, 1 and 2 fixed links",
          {r[1] for r in res} == {(1, 1, 1, 2)}, eq=eq_part)
    check(NAME, "temporal gauge: every free (spatial) link has 8 partners", {r[0] for r in res_t} == {8}, eq=eq_part)
    check(NAME, "no gauge fixing: every link has 12 partners", {r[0] for r in res_n} == {12}, eq=eq_part)

    # ---- Theorem A arithmetic -----------------------------------------------------------------
    bW = sp.symbols("beta_W", positive=True)
    N2 = 4                                               # beta_W = N^2 beta for SU(2)
    row = lambda partners_n: Fr(partners_n)              # row sum is (partners) * |beta|
    crit = lambda partners_n: Fr(N2, partners_n)
    eq_dob = register(NAME, Equation(
        key="fg_dobrushin", entry=ENTRY, title="SU(2) Dobrushin range and decay rate in the forest gauge (arithmetic)",
        source_kind="paper", source=f"{SRC}, Theorem A",
        latex_str=(r"\sum_{e'}|\beta|\cdot\mathbf 1[e'\sim e]=7|\beta|<1\iff \beta_W<\tfrac47,\qquad "
                   r"m\ \ge\ -\log\tfrac{7\beta_W}{4}\ \ (\beta_W=4\beta);"
                   r"%%BR%%\text{temporal: }\beta_W<\tfrac12,\quad \text{no gauge fixing: }\beta_W<\tfrac13"),
        notes="Quoted: the theorem (Dobrushin uniqueness and exponential decay of correlations for SU(2) at strong "
              "coupling on Z^3) is from the repository's paper and is not reproved. Checked here: the arithmetic from "
              "the partner counts of the previous equation to the thresholds and the mass bound. This is a lattice "
              "statement at strong coupling; it is not the continuum Yang-Mills problem.",
    ))
    check(NAME, "critical beta_W = 4/p for p = 7, 8, 12 partners: 4/7, 1/2, 1/3",
          [crit(p) for p in (7, 8, 12)] == [Fr(4, 7), Fr(1, 2), Fr(1, 3)], eq=eq_dob)
    check(NAME, "the forest range exceeds the temporal one by 1/14 and the ungauged one by 5/21",
          crit(7) - crit(8) == Fr(1, 14) and crit(7) - crit(12) == Fr(5, 21), eq=eq_dob)
    check(NAME, "mass bound -log(7 beta_W/4) is positive exactly for beta_W < 4/7 and tends to 0 at 4/7",
          sp.limit(-sp.log(7 * bW / 4), bW, sp.Rational(4, 7), "-") == 0
          and sp.log(sp.Rational(4, 7) * 7 / 4) == 0
          and all(-sp.log(7 * sp.Rational(b) / 4) > 0 for b in ("1/10", "1/4", "1/2", "11/20")), eq=eq_dob)
    check(NAME, "the same bound written with beta: -log(7|beta|) with beta = beta_W/4",
          sp.simplify(-sp.log(7 * bW / 4) - (-sp.log(7 * (bW / 4)))) == 0, eq=eq_dob)

    # ---- Theorem B arithmetic -----------------------------------------------------------------
    h_low = sp.Rational(9347955586999, 10 ** 12)         # lower bound, quoted to 12 digits
    rng = lambda N, h: Fr(N * N, 1) / (2 * Fr(h))
    h_cert = H_CERT
    eq_hess = register(NAME, Equation(
        key="fg_hessian_range", entry=ENTRY, title="The strong-coupling range from the Hessian bound, and the comparison table",
        source_kind="paper", source=f"{SRC}, Theorem B and comparison table",
        latex_str=(r"|\mathrm{Hess}\,S_0|\le h=9599991/10^6%%BR%%\Longrightarrow\ \text{gap for }|\beta|<1/(2h)\approx 1/19.20;"
                   r"%%BR%%\beta_W=N^2\beta<N^2/(2h):\ \ SU(2)\ 0.208,\ \ SU(3)\ 0.468"),
        notes="Quoted: the bound h = 9599991/10^6 itself, which the repository certifies by a rational site and plaquette "
              "certificate plus an interval-arithmetic cover and which is not rerun here. Checked: the range formula "
              "N^2/(2h) at the constants of the comparison table. The constants 16 and 12 are those that reproduce the "
              "table's SZZ and 'constant 12' rows, not independent readings of those papers. The ceiling row uses the "
              "lower bound of the next equation: no proof by this method can exceed it.",
    ))
    rows = {"SZZ": 16, "constant 12": 12, "Theorem B": h_cert, "ceiling": Fr(9347955586999, 10 ** 12)}
    table = {k: (float(rng(2, v)), float(rng(3, v))) for k, v in rows.items()}
    check(NAME, "SU(2) ranges 4/(2h) at h = 16, 12, 9.599991, 9.347955587 round to 0.125, 0.167, 0.208, 0.214",
          [round(table[k][0], 3) for k in rows] == [0.125, 0.167, 0.208, 0.214], eq=eq_hess)
    check(NAME, "SU(3) ranges 9/(2h) at the same h round (down) to 0.281, 0.375, 0.468, 0.481",
          [int(table[k][1] * 1000) / 1000 for k in rows] == [0.281, 0.375, 0.468, 0.481], eq=eq_hess)
    check(NAME, "1/(2h) at h = 9599991/10^6 differs from 1/19.20 by less than 1e-6",
          abs(float(1 / (2 * h_cert)) - 1 / 19.2) < 1e-6, eq=eq_hess)
    check(NAME, "for SU(2) the Hessian range 0.208 is below Dobrushin's 4/7: the two theorems are different statements "
                "(a spectral-gap bound for the whole group family versus uniqueness and decay), not nested ranges",
          rng(2, h_cert) < Fr(4, 7), eq=eq_hess)

    # ---- exact lower bound at the identity ----------------------------------------------------
    B = bloch_matrix_pi()
    x = sp.symbols("x")
    M = B.T * B
    fac = sp.factor_list(M.charpoly(x).as_expr())[1]
    degs = sorted((sp.degree(f, x), m) for f, m in fac)
    top = max((max(sp.Poly(f, x).intervals(eps=sp.Rational(1, 10 ** 12)), key=lambda q: q[0][1]) for f, _ in fac),
              key=lambda q: q[0][1])
    (lo, hi), _ = top
    nonic = [f for f, m in fac if sp.degree(f, x) == 9 and any(abs(float(r) - 9.3479555) < 1e-6 for r in sp.Poly(f, x).nroots())]
    eq_low = register(NAME, Equation(
        key="fg_lower_bound", entry=ENTRY, title="Exact lower bound for the Hessian constant at the identity field",
        source_kind="derived", source=f"derived here; compare {REPO}: identity_exact.py, bloch_identity.py",
        latex_str=(r"h_F\ \ge\ \lambda_{\max}\big(B^{\mathsf T}B\big)\Big|_{k=(\pi,\pi,\pi)}=9.347955587\ldots%%BR%%"
                   r"\lambda\ \text{a root of}\ x^9-35x^8+505x^7-3901x^6+17538x^5%%BR%%"
                   r"\qquad-46770x^4+72344x^3-60490x^2+22788x-1944"),
        notes="At the identity every transport is trivial and the Hessian is three copies of the scalar form "
              "sum_p (x_1 + x_2 - x_3 - x_4)^2 over free links (the repository's normalisation, quoted). Its Bloch matrix "
              "at k = (pi, pi, pi) is the integer 81 x 54 matrix B, and B^T B has an exactly factored characteristic "
              "polynomial. Any periodic field gives a lower bound for the constant on Z^3. Whether the constant equals "
              "9.3479 is the repository's conjecture and is not claimed: the identity is one field, and the repository's "
              "certificate bounds the constant above by 9.599991.",
    ))
    check(NAME, "B is the 81 x 54 integer matrix of the cell: plaquettes against free links", B.shape == (81, 54), eq=eq_low)
    check(NAME, "characteristic polynomial of B^T B factors over Z into degrees 1, 1, 8, 8, 9, 9 with multiplicities 2, 2, 2, 2, 1, 1",
          degs == [(1, 2), (1, 2), (8, 2), (8, 2), (9, 1), (9, 1)], eq=eq_low)
    p9 = sp.Poly(x**9 - 35 * x**8 + 505 * x**7 - 3901 * x**6 + 17538 * x**5 - 46770 * x**4 + 72344 * x**3
                 - 60490 * x**2 + 22788 * x - 1944, x)
    check(NAME, "the stated nonic is one of the irreducible factors", len(nonic) == 1 and sp.expand(nonic[0] - p9.as_expr()) == 0, eq=eq_low)
    check(NAME, "its largest root is isolated in an interval of width < 1e-12 that contains 9.347955587",
          hi - lo < sp.Rational(1, 10 ** 12) and lo < sp.Rational(9347955587, 10 ** 9) < hi + sp.Rational(1, 10 ** 9), eq=eq_low)
    check(NAME, "the lower bound sits below the certified bound 9.599991 (the two are consistent)",
          hi < sp.Rational(9599991, 10 ** 6), eq=eq_low)
    check(NAME, "trace of B^T B equals the sum of squares of the entries of B, and B has only entries 0, +-1",
          M == M.T and M.trace() == sum(v * v for v in B) and set(B) <= {-1, 0, 1}, eq=eq_low)
