"""Arithmetic and finite structure of "Explicit Lattice Yang-Mills Mass Gap for Every Compact Simple Group Reducing the Continuum Gap to One Cube"
(D. Bhattacharjee, U. Bhattacharya, J. M. Lockwood and T. Romanov; dated 30 October 2026; attached file eb95ffa3-YangMillsForestGauge.pdf, 98 pages;
Hexagon preprint 2610.00341v4 as the paper gives it). The fork read pages 1-6 and the appendix with Table III; it did not read the proofs.

This is the companion of forest_gauge.py, forest_gauge_4d.py (the repository's period-2 pattern and the bound h_F <= 14.420001) and forest_gauge_3d.py
(the three-dimensional paper). It rebuilds, without copying code, the finite statements of the four-dimensional paper that can be recomputed:

  fg4t_linefield   Table III (the line field sigma on Z^4 of period two): at each of the 16 residue classes pi_x fixes sigma(x) and is a three-cycle on the other
                   three directions (Lemma 8), every plaquette of a periodic lattice contains exactly one link of F, every site has one outgoing and one incoming
                   link of F in the direction sigma(x) (Lemma 9), F meets a box in a forest, a link outside F shares a plaquette with exactly twelve other links
                   outside F, two in each plaquette, four of its plaquettes containing a link of F that touches it (Lemma 9 iii, iv); the count of period-two
                   direction fields with the property (footnote 11: eight), and the dimension count 2(d-1) = d(d-1)/2, which holds only for d = 4.
  fg4t_partners    The partner numbers: 12 in the line field, 14 in the temporal gauge, 18 with no gauge fixing; the rule 4/p and the ceilings 1/3 and 2/9
                   (Theorem 2 and the introduction); the 1/4 of the Kantorovich modulus is quoted.
  fg4t_table1      Table I from Lie-algebra data computed here: the root systems, the dual Coxeter number h^v, the Dynkin index l_rho of the smallest
                   faithful representation (by the Weyl dimension formula and the Casimir), kappa = h^v/(2 l_rho), and the range 2 N kappa/29, for nine families;
                   Lichnerowicz's bound against the Casimir of the defining representation, with equality exactly for SU(2).
  fg4t_theorems    The numbers of Theorems 1, 3, 4, 5, 6 and 7, Table II and footnote 9: gaps at the stated points, the range 4/29 and 9/29, the decay rate
                   m_G(beta) = arccosh(1 + (kappa/(N|beta|) - 29/2)/12) and its limits, the factor 29/28.840002.

Not tested: every theorem and every estimate of the paper. The proof that the Kantorovich modulus of the SU(2) one-link law is one quarter, the
Dobrushin argument, the plaquette-block constants alpha(beta_W) and alpha_3(beta_W) (quoted at the stated points only), the certified Hessian bound
h_F <= 14.420001 (see forest_gauge_4d.py), reflection positivity, the transfer matrix, the three open inputs (I1)-(I3), and every statement about the
continuum. The paper itself states that the continuum limit and the gap in physical units remain open.
"""
import itertools
import math
import sys
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import sympy as sp

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register

NAME = "forest_gauge_4d_theorems"
ENTRY = "forest-gauge"

# Table III of the paper: x modulo two (digits x0 x1 x2 x3), sigma(x), then sigma(x + e_mu) for mu = 0..3
TABLE_III = """
0000 0  0 2 3 1
1000 0  0 3 1 2
0100 2  3 0 2 1
1100 3  2 0 1 3
0010 3  1 2 0 3
1010 1  3 1 0 2
0110 2  1 3 2 0
1110 1  2 1 3 0
0001 1  2 1 3 0
1001 2  1 3 2 0
0101 1  3 1 0 2
1101 3  1 2 0 3
0011 3  2 0 1 3
1011 2  3 0 2 1
0111 0  0 3 1 2
1111 0  0 2 3 1
"""
SIGMA, PI = {}, {}
for line in TABLE_III.strip().splitlines():
    f = line.split()
    x = tuple(int(c) for c in f[0])
    SIGMA[x] = int(f[1])
    PI[x] = tuple(int(c) for c in f[2:])


def e(a):
    return tuple(int(i == a) for i in range(4))


def add(x, y, s=1):
    return tuple(p + s * q for p, q in zip(x, y))


def sig(x):
    return SIGMA[tuple(v % 2 for v in x)]


PLANES = tuple(itertools.combinations(range(4), 2))


def plaquette(x, a, b):
    return [(x, a), (add(x, e(a)), b), (add(x, e(b)), a), (x, b)]


def is_F(link):
    x, a = link
    return sig(x) == a


def torus_links(L):
    return [(x, a) for x in itertools.product(range(L), repeat=4) for a in range(4)]


def wrap(link, L):
    x, a = link
    return (tuple(v % L for v in x), a)


def plaquettes_torus(L):
    return [[wrap(l, L) for l in plaquette(x, a, b)] for x in itertools.product(range(L), repeat=4) for a, b in PLANES]


def partners(fixed, L):
    """For every free link of the L^4 torus (L >= 4, so no two links share two plaquettes): the number of plaquettes through it and its set of free partners."""
    through = {}
    for p in plaquettes_torus(L):
        for l in p:
            through.setdefault(l, []).append(p)
    out = {}
    for l, plist in through.items():
        if fixed(l):
            continue
        part, per = set(), []
        for p in plist:
            fr = [q for q in p if not fixed(q) and q != l]
            per.append(len(fr))
            part.update(fr)
        out[l] = (len(plist), len(part), tuple(per))
    return out


# ---------------------------------------------------------------------------------------------------------------
# the line field
# ---------------------------------------------------------------------------------------------------------------
def three_cycle_on_rest(s, pi_row):
    """pi_x(s) = s, and pi_x permutes the other three directions with no fixed point (a three-cycle)."""
    rest = [m for m in range(4) if m != s]
    if pi_row[s] != s:
        return False
    img = [pi_row[m] for m in rest]
    return sorted(img) == rest and all(pi_row[m] != m for m in rest)


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


def count_period_two_fields():
    """All sigma: {0,1}^4 -> {0,..,3} such that at every site pi_x fixes sigma(x) and is a three-cycle on the rest, pi_x(mu) = sigma(x + e_mu).
    Depth-first search over the 16 sites, each constraint tested when its five sites are assigned."""
    sites = list(itertools.product(range(2), repeat=4))
    idx = {s: i for i, s in enumerate(sites)}
    nb = {s: [tuple(v ^ (m == k) for k, v in enumerate(s)) for m in range(4)] for s in sites}
    need = {s: max(idx[t] for t in [s] + nb[s]) for s in sites}
    by_last = {}
    for s in sites:
        by_last.setdefault(need[s], []).append(s)
    sol, found = {}, []

    def ok(s):
        return three_cycle_on_rest(sol[s], tuple(sol[t] for t in nb[s]))

    def rec(i):
        if i == 16:
            found.append(dict(sol))
            return
        for v in range(4):
            sol[sites[i]] = v
            if all(ok(s) for s in by_last.get(i, [])):
                rec(i + 1)
        del sol[sites[i]]
    rec(0)
    return found


def build_and_verify():
    L = 4
    fields = count_period_two_fields()
    rows_ok = all(three_cycle_on_rest(SIGMA[x], PI[x]) for x in SIGMA)
    table_consistent = all(PI[x][m] == SIGMA[tuple(v ^ (k == m) for k, v in enumerate(x))] for x in SIGMA for m in range(4))
    planes = plaquettes_torus(L)
    per_plaq = [sum(is_F(l) for l in p) for p in planes]
    out_cnt = {x: sum(1 for a in range(4) if is_F((x, a))) for x in itertools.product(range(L), repeat=4)}
    in_cnt = {x: [a for a in range(4) if is_F((add(x, e(a), -1), a))] for x in itertools.product(range(L), repeat=4)}
    box = [(x, a) for x in itertools.product(range(7), repeat=4) for a in range(4) if is_F((x, a)) and all(v + (i == a) <= 6 for i, v in enumerate(x))]
    par_F = partners(is_F, L)
    shape_F = {(v[0], v[1], v[2]) for v in par_F.values()}
    # Lemma 9 (iii): of the six plaquettes of e = (x, mu) not in F, four hold a link of F touching e (in planes {mu, sigma(x)} and {mu, pi_x(mu)}), two hold the opposite one
    touch_ok = True
    for x in itertools.product(range(L), repeat=4):
        for mu in range(4):
            if is_F((x, mu)):
                continue
            s, t = sig(x), PI[tuple(v % 2 for v in x)][mu]
            cnt_touch = 0
            for nu in range(4):
                if nu == mu:
                    continue
                for base in (x, add(x, e(nu), -1)):
                    p = plaquette(base, *sorted((mu, nu)))
                    F_link = [l for l in p if is_F(l)][0]
                    ends = {F_link[0], add(F_link[0], e(F_link[1]))}
                    touching = bool(ends & {x, add(x, e(mu))})
                    if touching:
                        cnt_touch += 1
                        touch_ok &= nu in (s, t)
            touch_ok &= cnt_touch == 4
    d = sp.Symbol("d", positive=True)
    d_sol = sp.solve(sp.Eq(2 * (d - 1), d * (d - 1) / 2), d)

    eq_lf = register(NAME, Equation(
        key="fg4t_linefield", entry=ENTRY, title="The line field of period two in four dimensions: one fixed link in every plaquette and twelve free partners",
        source_kind="derived", source="derived here from Table III and Lemmas 8 and 9 of the four-dimensional forest-gauge paper (dated 30 October 2026)",
        latex_str=(r"\pi_x(\sigma(x))=\sigma(x),\ \pi_x|_{R}\ \text{a 3-cycle},\ \ \pi_x(\mu)=\sigma(x+e_\mu)\ \Rightarrow\ \#\{F\cap\partial p\}=1\ \ \forall p"
                   r"%%BR%%2(d-1)=\frac{d(d-1)}{2}\iff d=4,\qquad \#\text{partners}=6\cdot2=12,\qquad \#\{\text{period-2 fields}\}=" + str(len(fields))),
        notes="Table III is transcribed from the paper's Appendix A and every row is tested: pi_x fixes sigma(x) and is a fixed-point-free permutation of the other "
              "three directions, which is Lemma 8, and pi_x(mu) = sigma(x + e_mu) is consistent across the table. The plaquette count, the outgoing and incoming link "
              "of F at each site, the forest property in a box and the twelve partners are computed on the 4^4 torus (links share at most one plaquette when L >= 4) "
              "and on a box, not read from the paper. The count of period-two fields is a depth-first search over all 4^16 assignments with the Lemma 8 "
              "constraint tested at each site; the paper states eight, equivalent under the symmetries of the lattice. Not tested: that the eight are equivalent "
              "under the lattice symmetries, and the Lean certificate of the line field that the paper mentions.",
    ))
    check(NAME, "Table III has 16 rows, one per residue class of x modulo two, with sigma(x) in {0,1,2,3}", len(SIGMA) == 16 and set(SIGMA.values()) <= {0, 1, 2, 3}, eq=eq_lf)
    check(NAME, "the table is self-consistent: the entry pi_x(mu) equals sigma(x + e_mu) read from the row of x + e_mu, for all 16 x and 4 mu", table_consistent, eq=eq_lf)
    check(NAME, "Lemma 8: in every row pi_x fixes sigma(x) and is a three-cycle (no fixed point) on the other three directions", rows_ok, eq=eq_lf)
    check(NAME, f"every one of the {len(planes)} plaquettes of the 4^4 torus contains exactly one link of F (counted directly from F = {{(x, sigma(x))}})",
          len(planes) == 256 * 6 and set(per_plaq) == {1}, eq=eq_lf)
    check(NAME, "Lemma 9 (i): every site has exactly one outgoing link of F, and exactly one incoming one, both in the direction sigma(x)",
          all(c == 1 for c in out_cnt.values()) and all(in_cnt[x] == [sig(x)] for x in in_cnt), eq=eq_lf)
    check(NAME, f"Lemma 9 (ii): F meets a box of {7 ** 4} sites in {len(box)} links that contain no cycle (a forest; straight segments of the lines)", not has_cycle(box) and len(box) > 0, eq=eq_lf)
    check(NAME, f"Lemma 9 (iv): every link outside F lies in 6 plaquettes, with exactly two free partners in each, so exactly 12 distinct partners; observed (planes, partners, per plaquette) "
                f"= {sorted(shape_F)} over {len(par_F)} free links", shape_F == {(6, 12, (2,) * 6)} and len(par_F) == 3 * 4 ** 4, eq=eq_lf)
    check(NAME, "Lemma 9 (iii): of the six plaquettes of a free link, four contain a link of F touching it, lying in the planes {mu, sigma(x)} and {mu, pi_x(mu)}; "
                "in the other two the link of F is opposite it (all free links of the torus)", touch_ok, eq=eq_lf)
    check(NAME, "a link lies in 2(d-1) plaquettes and a site carries d(d-1)/2 of them, so a field with one fixed link per site and one per plaquette needs "
                "2(d-1) = d(d-1)/2: the only solution with d > 1 is d = 4 (sympy)", sorted(d_sol) == [1, 4] or sorted(d_sol) == [4], eq=eq_lf)
    check(NAME, f"footnote 11: exactly eight direction fields of period two have the property of Lemma 8 (exhaustive search over 4^16 assignments, {len(fields)} found), "
                "and Table III is one of them", len(fields) == 8 and any(all(f[x] == SIGMA[x] for x in SIGMA) for f in fields), eq=eq_lf)

    # ---------------------------------------------------------------------------------------------------------------
    # partner numbers and ceilings
    # ---------------------------------------------------------------------------------------------------------------
    par_T = partners(lambda l: l[1] == 0, L)
    par_N = partners(lambda l: False, L)
    shape_T, shape_N = {(v[0], v[1]) for v in par_T.values()}, {(v[0], v[1]) for v in par_N.values()}
    spatial_T = {v[1] for l, v in par_T.items()}
    ceil = {p: Fr(4, p) for p in (12, 14, 18)}
    eq_pt = register(NAME, Equation(
        key="fg4t_partners", entry=ENTRY, title="Partners of a free link and the Dobrushin ceiling 4/p in four dimensions",
        source_kind="derived", source="derived here from the introduction, footnote 4 and Theorem 2 of the four-dimensional forest-gauge paper",
        latex_str=(r"p_{\text{line field}}=12,\quad p_{\text{temporal}}=14,\quad p_{\text{no fixing}}=3\cdot6=18"
                   r"%%BR%%\beta_W<\frac4p:\quad\frac4{12}=\frac13,\quad\frac4{14}=\frac27,\quad\frac4{18}=\frac29;\qquad p=6(d-1)\ \Rightarrow\ \frac{2}{3(d-1)}"),
        notes="With the one-link Kantorovich modulus one quarter (quoted, not rederived), the weighted Dobrushin sum for a link with p partners is p beta_W/4, "
              "so the condition p beta_W/4 < 1 gives beta_W < 4/p. The partner numbers are counted on the 4^4 torus: 12 in the line field, 14 in the temporal gauge "
              "(all links in direction 0 fixed: 2 partners in each of the two plaquettes of the planes with the fixed direction, 3 in each of the four others, "
              "2 x 1 + 4 x 3 = 14, as footnote 4 says), 18 with no gauge fixing (six plaquettes, three partners each). The ceiling statement of Theorem 2 for the "
              "line field (the weighted Dobrushin condition fails at beta_W >= 1/3 with single-link coefficients) is the same arithmetic. Not tested: the "
              "modulus 1/4, and that no better single-link argument exists.",
    ))
    check(NAME, f"line field: 12 partners (counted); temporal gauge: {sorted(spatial_T)} partners for every spatial free link, in 6 plaquettes; no gauge fixing: {sorted({v[1] for v in par_N.values()})} partners",
          shape_T == {(6, 14)} and shape_N == {(6, 18)}, eq=eq_pt)
    check(NAME, "the ceilings 4/12 = 1/3 (line field), 4/14 = 2/7 (temporal), 4/18 = 2/9 (no gauge fixing: the 'sharp quarter' value of the ungauged measure in Table II)",
          ceil[12] == Fr(1, 3) and ceil[14] == Fr(2, 7) and ceil[18] == Fr(2, 9), eq=eq_pt)
    dd = sp.Symbol("d")
    check(NAME, "p = 6(d-1) partners without gauge fixing gives the ceiling 4/(6(d-1)) = 2/(3(d-1)); at d = 4 this is 2/9, at d = 3 it is 1/3 (as the three-dimensional paper has)",
          sp.simplify(4 / (6 * (dd - 1)) - 2 / (3 * (dd - 1))) == 0 and Fr(2, 3 * 3) == Fr(2, 9) and Fr(2, 3 * 2) == Fr(1, 3), eq=eq_pt)
    check(NAME, "with weight beta_W/4 per partner the row sum is 12 beta_W/4 = 3 beta_W, which is < 1 exactly for beta_W < 1/3 (Theorem 1: m >= -log(3 beta_W))",
          sp.solve(sp.Eq(12 * sp.Symbol("b") / 4, 1), sp.Symbol("b")) == [sp.Rational(1, 3)], eq=eq_pt)
    check(NAME, "the Hessian of the Wilson action on a geodesic through Q_e is -6 beta_W times the metric of the unit S^3 whose Ricci curvature is 2, so the Bakry-Emery tensor "
                "2 - 6 beta_W is negative exactly for beta_W > 1/3 (the paper's remark on why curvature in the gauge-fixed SU(2) measure stops at 1/3)",
          sp.solve(sp.Eq(2 - 6 * sp.Symbol("b"), 0), sp.Symbol("b")) == [sp.Rational(1, 3)], eq=eq_pt)

    # ---------------------------------------------------------------------------------------------------------------
    # Table I from root systems
    # ---------------------------------------------------------------------------------------------------------------
    def gram(kind, r):
        """Symmetric Gram matrix of the simple roots, long roots normalised to (alpha, alpha) = 2, Bourbaki numbering."""
        R = sp.Rational
        G = sp.zeros(r, r)

        def put(i, j, v):
            G[i - 1, j - 1] = G[j - 1, i - 1] = v
        ln = {i: 2 for i in range(1, r + 1)}
        edges = []
        if kind == "A":
            edges = [(i, i + 1, -1) for i in range(1, r)]
        elif kind == "B":
            ln[r] = 1
            edges = [(i, i + 1, -1) for i in range(1, r)]
        elif kind == "C":
            for i in range(1, r):
                ln[i] = 1
            edges = [(i, i + 1, R(-1, 2)) for i in range(1, r - 1)] + [(r - 1, r, -1)]
        elif kind == "D":
            edges = [(i, i + 1, -1) for i in range(1, r - 1)] + [(r - 2, r, -1)]
        elif kind == "G":
            ln = {1: R(2, 3), 2: 2}
            edges = [(1, 2, -1)]
        elif kind == "F":
            ln = {1: 2, 2: 2, 3: 1, 4: 1}
            edges = [(1, 2, -1), (2, 3, -1), (3, 4, R(-1, 2))]
        elif kind == "E":
            edges = [(1, 3, -1), (3, 4, -1), (4, 5, -1), (2, 4, -1)] + [(i, i + 1, -1) for i in range(5, r)]
        for i in range(1, r + 1):
            put(i, i, R(ln[i]))
        for i, j, v in edges:
            put(i, j, R(v))
        return G


    def positive_roots(G):
        """G: list of lists of Fractions."""
        r = len(G)
        simple = [tuple(int(i == j) for j in range(r)) for i in range(r)]

        def ip(a, b):
            return sum(a[i] * b[j] * G[i][j] for i in range(r) for j in range(r) if a[i] and b[j])
        roots, frontier = set(simple), list(simple)
        while frontier:
            new = []
            for b in frontier:
                for i in range(r):
                    c = 2 * ip(b, simple[i]) / G[i][i]
                    s = tuple(int(b[k] - c * simple[i][k]) for k in range(r))
                    if s not in roots:
                        roots.add(s)
                        new.append(s)
            frontier = new
        return sorted(a for a in roots if all(v >= 0 for v in a)), ip


    def lie_data(kind, r, hw):
        """(dim g, h^v, dim V, Dynkin index, Casimir) of the representation with highest weight hw = (m_1..m_r) in fundamental weights; exact Fractions."""
        Gs = gram(kind, r)
        G = [[Fr(int(Gs[i, j].p), int(Gs[i, j].q)) for j in range(r)] for i in range(r)]
        pos, ip = positive_roots(G)
        n_g = r + 2 * len(pos)

        def wip(m, a):          # (lambda, alpha) = sum m_i a_i (alpha_i, alpha_i)/2
            return sum(m[i] * a[i] * G[i][i] / 2 for i in range(r))
        Gi = Gs.inv()           # omega_k = sum_j c_j alpha_j with sum_j c_j (alpha_j, alpha_i) = delta_ik (alpha_k, alpha_k)/2
        om = [[Fr(int((Gi[k, j] * Gs[k, k] / 2).p), int((Gi[k, j] * Gs[k, k] / 2).q)) for j in range(r)] for k in range(r)]
        rho_w = (1,) * r

        def ww(m1, m2):
            a1 = [sum(m1[k] * om[k][j] for k in range(r)) for j in range(r)]
            return wip(m2, a1)
        C2 = ww(hw, hw) + 2 * ww(hw, rho_w)
        lam_rho = tuple(m + 1 for m in hw)
        dim = Fr(1)
        for a in pos:
            dim *= wip(lam_rho, a) / wip(rho_w, a)
        longest = max((a for a in pos if ip(a, a) == 2), key=sum)
        hv = (ip(longest, longest) + 2 * wip(rho_w, longest)) / 2
        return n_g, hv, dim, dim * C2 / n_g, C2


    def fund(kind, r, k=1, mult=1):
        m = [0] * r
        m[k - 1] = mult
        return tuple(m)


    FAMILIES = []
    for N in (2, 3, 4, 5):
        FAMILIES.append((f"SU({N})", "A", N - 1, fund("A", N - 1), N, N, 1, Fr(N, 2), Fr(N * N, 29)))
    for N in (5, 6, 7, 8):
        kind, r = ("B", (N - 1) // 2) if N % 2 else ("D", N // 2)
        FAMILIES.append((f"SO({N})", kind, r, fund(kind, r), N, N - 2, 2, Fr(N - 2, 4), Fr(N * (N - 2), 58)))
    FAMILIES.append(("SO(3)", "A", 1, (2,), 3, 2, 4, Fr(1, 4), Fr(3, 58)))
    for n in (2, 3, 4):
        FAMILIES.append((f"Sp({n})", "C", n, fund("C", n), 2 * n, n + 1, 1, Fr(n + 1, 2), Fr(2 * n * (n + 1), 29)))
    FAMILIES.append(("G2", "G", 2, (1, 0), 7, 4, 2, Fr(1), Fr(14, 29)))
    FAMILIES.append(("F4", "F", 4, (0, 0, 0, 1), 26, 9, 6, Fr(3, 4), Fr(39, 29)))
    FAMILIES.append(("E6", "E", 6, fund("E", 6, 1), 27, 12, 6, Fr(1), Fr(54, 29)))
    FAMILIES.append(("E7", "E", 7, fund("E", 7, 7), 56, 18, 12, Fr(3, 4), Fr(84, 29)))
    FAMILIES.append(("E8", "E", 8, fund("E", 8, 8), 248, 30, 60, Fr(1, 4), Fr(124, 29)))

    rows, dims_ok, hv_ok, ell_ok, kap_ok, rng_ok, lich_ok = [], True, True, True, True, True, True
    equal_cases = []
    for name, kind, r, hw, N, hv_p, ell_p, kap_p, rng_p in FAMILIES:
        n_g, hv, dim, ell, C2 = lie_data(kind, r, hw)
        kap = hv / (2 * ell)
        dims_ok &= dim == N
        hv_ok &= hv == hv_p
        ell_ok &= ell == ell_p
        kap_ok &= kap == kap_p
        rng_ok &= Fr(2 * N, 1) * kap / 29 == rng_p
        # Lichnerowicz: lambda_1 >= n/(n-1) Ric on a compact n-manifold; matrix coefficients of the defining representation have eigenvalue dim G / N in the metric -Tr
        lam = Fr(int(n_g), N)
        bound = Fr(int(n_g), int(n_g) - 1) * kap
        lich_ok &= lam >= bound
        if lam == bound:
            equal_cases.append(name)
        rows.append((name, int(n_g), int(hv), int(ell), str(kap), str(rng_p), lam, bound))

    eq_t1 = register(NAME, Equation(
        key="fg4t_table1", entry=ENTRY, title="Table I of the four-dimensional paper from root-system data, and Lichnerowicz's bound",
        source_kind="derived", source="derived here from the Lie-algebra data of the compact simple groups; Table I and Theorem 5 are the paper's",
        latex_str=(r"\kappa_G=\frac{h^{\vee}}{2\ell_\rho},\quad \ell_\rho=\frac{\dim V_\rho\,C_2(\rho)}{\dim\mathfrak g},\quad C_2(\lambda)=(\lambda,\lambda+2\rho),\quad "
                   r"\beta_W<\frac{2N\kappa_G}{29}"
                   r"%%BR%%SU(N):\ \kappa=\tfrac N2,\ \tfrac{N^2}{29};\quad SO(N):\ \tfrac{N-2}{4},\ \tfrac{N(N-2)}{58};\quad Sp(n):\ \tfrac{n+1}{2},\ \tfrac{2n(n+1)}{29};\quad "
                   r"G_2,F_4,E_6,E_7,E_8:\ \kappa=1,\tfrac34,1,\tfrac34,\tfrac14"
                   r"%%BR%%\lambda_\rho=\frac{\dim G}{N}\ \ge\ \frac{\dim G}{\dim G-1}\,\kappa_G\ \ (\text{Lichnerowicz}),\qquad\text{equality only for }SU(2)"),
        notes="For each group the simple roots are built from the Dynkin diagram (long roots of squared length 2), the positive roots are generated by reflections, "
              "dim g = rank + 2 #positive roots, h^v = (theta, theta + 2 rho)/2 with theta the highest long root, the dimension of the representation by the Weyl "
              "dimension formula, and the Dynkin index by dim V (lambda, lambda + 2 rho)/dim g (which is 1 for the defining representation of SU(N) and Sp(n)). "
              "The representations are the paper's: the defining one for SU(N), SO(N) (N >= 5), Sp(n), and the smallest for the exceptional groups; SO(3) is the spin-one "
              "representation of SU(2). The ranges 2 N kappa/29 are the last column of Table I. The metric of the paper is -Tr(rho(X) rho(Y)): in it "
              "the matrix coefficients of the defining representation are eigenfunctions of the Laplacian with eigenvalue dim G/N, since the trace of the "
              "Casimir is dim G. Lichnerowicz's theorem then gives dim G/N >= (dim G/(dim G - 1)) kappa as a consistency test of the tabulated kappa: "
              "none of the seventeen groups violates it. Not tested: the 29/2 Hessian bound, the Bakry-Emery argument, and the exceptional groups' smallest faithful "
              "representations beyond their dimensions (7, 26, 27, 56, 248 are reproduced).",
    ))
    check(NAME, "the Weyl dimension formula reproduces the dimension of every representation of Table I: " + ", ".join(f"{r[0]}" for r in rows), dims_ok, eq=eq_t1)
    check(NAME, "the dual Coxeter number h^v from the root systems equals the table's for all " + str(len(rows)) + " groups (SU(N): N, SO(N): N-2, Sp(n): n+1, G2: 4, F4: 9, E6: 12, E7: 18, E8: 30, SO(3): 2)", hv_ok, eq=eq_t1)
    check(NAME, "the Dynkin index l_rho computed from the Casimir equals the table's (1, 2, 1, 2, 6, 6, 12, 60 and 4 for SO(3))", ell_ok, eq=eq_t1)
    check(NAME, "kappa_G = h^v/(2 l_rho) equals the table's value for all groups (SU(N): N/2; SO(N): (N-2)/4; SO(3): 1/4; Sp(n): (n+1)/2; G2: 1; F4: 3/4; E6: 1; E7: 3/4; E8: 1/4)", kap_ok, eq=eq_t1)
    check(NAME, "the range 2 N kappa/29 equals the last column: SU(N): N^2/29; SO(N): N(N-2)/58; SO(3): 3/58; Sp(n): 2n(n+1)/29; G2: 14/29; F4: 39/29; E6: 54/29; E7: 84/29; E8: 124/29",
          rng_ok, eq=eq_t1)
    check(NAME, "Lichnerowicz: dim G/N >= dim G/(dim G - 1) kappa_G for every group above, a necessary condition on the tabulated curvatures", lich_ok, eq=eq_t1)
    check(NAME, f"equality in Lichnerowicz holds exactly for {equal_cases}: for SU(2) the j = 1/2 Casimir, 4 j(j+1)/2 = 3/2 in the metric -Tr, equals "
                "(3/2) Ric on S^3 (radius sqrt 2, Ric = 1): the smallest electric excitation of one link saturates the Lichnerowicz bound, as the round sphere must (Obata)",
          equal_cases == ["SU(2)"] and Fr(3, 2) == Fr(3, 2) * 1, eq=eq_t1)
    check(NAME, "for SU(3) the same comparison is 8/3 against 12/7: no saturation (SU(3) is not a round sphere)", Fr(8, 3) > Fr(8, 7) * Fr(3, 2) == Fr(12, 7), eq=eq_t1)

    # ---------------------------------------------------------------------------------------------------------------
    # the numbers of the theorems
    # ---------------------------------------------------------------------------------------------------------------
    b = mp.mpf
    def gap_thm1(bw):
        return -mp.log(3 * bw)
    def gap_thm3(alpha):
        return -mp.log(alpha) / 2
    def gap_thm5(kappa, N, beta):
        return mp.acosh(1 + (kappa / (N * abs(beta)) - b(29) / 2) / 12)
    def gap_thm5_bw(kappa, N, bw):
        return gap_thm5(kappa, N, bw / N ** 2)
    m5_edge = [(N, kap) for N, kap in ((2, b(1)), (3, b(3) / 2))]
    ratio = b(29) / b("28.840002")
    eq_th = register(NAME, Equation(
        key="fg4t_theorems", entry=ENTRY, title="The numbers of Theorems 1-7 and Table II of the four-dimensional paper",
        source_kind="derived", source="derived here from the statements of Theorems 1, 3, 4, 5, 6 and 7, Table II and footnote 9 of the four-dimensional forest-gauge paper",
        latex_str=(r"m\ge-\log(3\beta_W)\ \ (\beta_W<\tfrac13);\qquad m\ge-\tfrac12\log\alpha(\beta_W),\ \alpha(\tfrac13)\le0.9561,\ \alpha(0.345)\le0.9973"
                   r"%%BR%%m\ge-\log(1.33319\,\beta_W)\ \ (SU(3),\ \beta_W\le\tfrac34)"
                   r"%%BR%%m_G(\beta)=\operatorname{arccosh}\!\Big(1+\tfrac1{12}\Big(\frac{\kappa_G}{N|\beta|}-\tfrac{29}2\Big)\Big)\to\log\frac{\kappa_G}{6N|\beta|},\qquad "
                   r"m_G\to0\ \text{at}\ \beta_W=\tfrac{2N\kappa_G}{29}"
                   r"%%BR%%\tfrac{29}{28.840002}=1.0055,\qquad SU(2):\ \tfrac1{12},\tfrac18,\tfrac29,\tfrac4{29},\tfrac13,0.345;\quad SU(3):\ \tfrac3{16},\tfrac9{32},\tfrac9{29},\tfrac34,0.76"),
        notes="Every number the paper prints for the theorems is recomputed from the formulas it states. The quantities alpha(beta_W) and alpha_3(beta_W) of the plaquette "
              "blocks are quoted at the stated points (0.9561 and 0.9973 for SU(2); 0.9182, 0.9653, 0.9950 for SU(3)); their derivation is not tested. "
              "The gap -1/2 log alpha(1/3) is 0.022445, which the paper prints as 0.0224; for SU(3) -1/2 log 0.9653 = 0.017658, which the paper prints as 0.0177: the printed number is the rounded value, 0.00004 above the bound that the stated "
              "alpha_3(3/4) <= 0.9653 gives (the alpha_3 is itself a rounded upper bound, so the true gap may be larger; the fork does not know). At 0.76, "
              "-1/2 log 0.99950 = 2.5e-4 as the paper prints. The decay rate of Theorem 5 is evaluated at SU(2) and SU(3) and compared with the single-link Dobrushin gap -log(3 beta_W). Not tested: the "
              "proofs, the Kantorovich modulus, the certified bounds, and the validity of the transfer-matrix reading.",
    ))
    r_su2_edge = 2 * 2 * b(1) / 29
    check(NAME, "Theorem 1 at beta_W = 1/3 gives gap 0 (the end of the range), at beta_W = 0.1 gives -log 0.3 = 1.204, and is positive exactly for beta_W < 1/3",
          abs(gap_thm1(b(1) / 3)) < 1e-12 and abs(gap_thm1(b("0.1")) - b("1.2039728043")) < 1e-9 and gap_thm1(b("0.33")) > 0 > gap_thm1(b("0.34")), eq=eq_th)
    check(NAME, f"Theorem 3: m >= -1/2 log alpha(1/3) = {mp.nstr(gap_thm3(b('0.9561')), 5)}, which the paper prints as 0.0224 (rounded down: the printed number is a valid lower bound); at beta_W = 0.345, "
                f"-1/2 log 0.9973 = {mp.nstr(gap_thm3(b('0.9973')), 4)}",
          b("0.0224") < gap_thm3(b("0.9561")) < b("0.0225") and b("0.00135") < gap_thm3(b("0.9973")) < b("0.00136"), eq=eq_th)
    g7_1 = gap_thm3(b("0.9653"))
    check(NAME, f"Theorem 7: -1/2 log alpha_3(3/4) = {mp.nstr(g7_1, 5)} with alpha_3(3/4) <= 0.9653; the paper prints 0.0177, which is the three-significant-figure rounding "
                "of that value and 0.2 per cent above it, so as a lower bound the printed 0.0177 slightly exceeds what the stated alpha_3 gives (recorded in the notes); at 0.76, "
                "-1/2 log 0.9995 = 2.5e-4 as printed",
          b("0.01765") < g7_1 < b("0.01766") and round(float(g7_1), 4) == 0.0177 and abs(gap_thm3(b("0.99950")) - b("0.00025")) < 1e-6, eq=eq_th)
    check(NAME, f"Theorem 6 (SU(3)): -log(1.33319 x 3/4) = {mp.nstr(-mp.log(b('1.33319') * b('0.75')), 4)}, printed as 1.1e-4; the range ends where 1.33319 beta_W = 1, at beta_W = {mp.nstr(1 / b('1.33319'), 6)} just above 3/4",
          abs(-mp.log(b("1.33319") * b("0.75")) - b("1.1e-4")) < 5e-6 and b("0.75") < 1 / b("1.33319") < b("0.7501"), eq=eq_th)
    check(NAME, "Theorem 4: the range for SU(2) is 4/29 = 0.1379 and for SU(3) 9/29 = 0.3103 (printed 0.310); Table II's ungauged rows are N^2/(2h) at h = 48 (Shen-Zhu-Zhu) and h = 16 (sharp): "
                "1/12 and 1/8 for SU(2), 3/16 and 9/32 for SU(3)",
          Fr(4, 29) == Fr(4, 1) / 29 and round(9 / 29, 3) == 0.310 and Fr(4, 48) == Fr(1, 12) and Fr(4, 32) == Fr(1, 8)
          and Fr(9, 48) == Fr(3, 16) and Fr(9, 32) == Fr(9, 32), eq=eq_th)
    check(NAME, "Table II, order of the ranges: SU(2): 1/12 < 1/8 < 4/29 < 2/9 < 1/3 < 0.345; SU(3): 3/16 < 9/32 < 9/29 < 3/4 < 0.76 (the forest-gauge curvature range lies above the sharp ungauged "
                "curvature range, and the single-link Dobrushin and block ranges lie above it)",
          Fr(1, 12) < Fr(1, 8) < Fr(4, 29) < Fr(2, 9) < Fr(1, 3) < Fr(345, 1000) and Fr(3, 16) < Fr(9, 32) < Fr(9, 29) < Fr(3, 4) < Fr(76, 100), eq=eq_th)
    mm = lambda k, N, bw: gap_thm5_bw(b(k), N, bw)
    check(NAME, "Theorem 5: m_G -> 0 at the end of the range beta_W = 2 N kappa/29 (SU(2): 4/29; SU(3): 9/29), is positive inside it and increasing as beta_W falls",
          abs(mm(1, 2, b(4) / 29)) < 1e-9 and abs(mm("1.5", 3, b(9) / 29)) < 1e-9 and mm(1, 2, b("0.1")) > mm(1, 2, b("0.13")) > 0 and mm("1.5", 3, b("0.2")) > mm("1.5", 3, b("0.3")) > 0, eq=eq_th)
    check(NAME, f"Theorem 5 as beta_W -> 0: m_G - log(kappa/(6N|beta|)) -> 0 (difference {mp.nstr(mm(1, 2, b('1e-8')) - mp.log(b(1) / (6 * 2 * b('1e-8') / 4)), 3)} at beta_W = 1e-8, SU(2))",
          abs(mm(1, 2, b("1e-8")) - mp.log(b(1) / (6 * 2 * b("1e-8") / 4))) < 1e-6, eq=eq_th)
    check(NAME, f"SU(2) at beta_W = 0.1: Theorem 5 gives {mp.nstr(mm(1, 2, b('0.1')), 4)}, Theorem 1 gives {mp.nstr(gap_thm1(b('0.1')), 4)}: in the common range the single-link gap is the larger, "
                "and Theorem 1 reaches 1/3 where Theorem 5 stops at 4/29 (the paper says the two stand side by side)",
          gap_thm1(b("0.1")) > mm(1, 2, b("0.1")) > 0, eq=eq_th)
    check(NAME, f"footnote 9: the factor 29/28.840002 = {mp.nstr(ratio, 6)}, printed 1.0055 (four decimals), the gain in every range of Table I from the certified bound 14.420001",
          round(float(ratio), 4) == 1.0055, eq=eq_th)
    check(NAME, "the block ranges pass the single-link ceilings: SU(2) 0.345 > 1/3, and SU(3) 19/25 = 0.76 > 3/4 by exactly 0.01",
          Fr(345, 1000) > Fr(1, 3) and Fr(19, 25) == Fr(76, 100) and Fr(19, 25) > Fr(3, 4) and Fr(19, 25) - Fr(3, 4) == Fr(1, 100), eq=eq_th)
