"""The four-dimensional period-2 forest gauge and the bound h_F <= 14.420001 (companion repository).

Source: Tristan Romanov, romanov360/yang-mills-4d-forest-certificate (code MIT, text CC BY 4.0), which
sharpens the Bhattacharjee-Bhattacharya-Lockwood (BBL) bound 29/2 for the gauge-fixed Hessian of the Wilson action
on Z^4 (D. Bhattacharjee, U. Bhattacharya and J. M. Lockwood, Zenodo doi:10.5281/zenodo.23729004). The fork
rederives, in exact arithmetic and without copying code, the combinatorics and arithmetic of that work. The
16-letter pattern below is the repository's data.

What is checked (exact unless marked):
  * the pattern fixes 16 of the 64 links of a period cell, every site is the head of one fixed link, the fixed links
    form a forest, and every one of the 96 plaquettes of the cell contains exactly one fixed link;
  * at the identity field and Bloch momentum 0 the 48 x 48 integer matrix B^T B has the factor x^2 - 18x + 64,
    so its top eigenvalue is 9 + sqrt(17) = 13.1231 (BBL's identity-field value); a coarse grid of other momenta
    does not exceed it (floating point);
  * the arithmetic of the ranges: 1/29 (BBL), 1/28.840002 (this bound), and the ceiling 1/(2 * 13.974).

Quoted, not rederived: the proved bound h_F <= 14420001/10^6 (a rational site-plus-plaquette certificate with
exact LDL^T pivots, and an all-phase sum-of-squares check in the Lean kernel, independently an interval-
arithmetic cover of 68,968,000 squares), the lower bound 13.974 (a floating-point hill climb, found by BBL and
reproduced by the repository), the reduction from SU(N) to U(2) phases, and the Hessian normalisation.

Not a continuum statement. The range |beta| < 1/(2 h_F) is a strong-coupling window at a fixed lattice spacing.
The repository states that it says nothing about the continuum limit, and it records that a gap bounded below
in lattice units uniformly in the spacing is incompatible with a nontrivial continuum limit (BBL's no-go
theorem), so no estimate of this kind can be the continuum step. Nothing here is a proof of Yang-Mills
existence or of a mass gap on R^4.
"""
import itertools
from collections import Counter
from fractions import Fraction as Fr

import numpy as np
import sympy as sp

from eqlib import Equation, check, register

NAME = "forest_gauge_4d"
ENTRY = "forest-gauge"
REPO = "romanov360/yang-mills-4d-forest-certificate"
SRC = f"{REPO} v1.0"
BBL = "Bhattacharjee, Bhattacharya and Lockwood, doi:10.5281/zenodo.23729004"

PATTERN = "0212331001332120"                    # the repository's period-2 pattern (direction of the head link)
CELLS = sorted(itertools.product(range(2), repeat=4))
MU = {c: int(PATTERN[i]) for i, c in enumerate(CELLS)}
H_PROVED = Fr(14420001, 10 ** 6)                 # quoted
H_BBL = Fr(29, 2)                                # quoted
H_LOWER = Fr(13974, 1000)                        # quoted, floating point


def unit(a):
    return tuple(int(i == a) for i in range(4))


def add(x, y, s=1):
    return tuple(p + s * q for p, q in zip(x, y))


def mu(x):
    return MU[tuple(v % 2 for v in x)]


def is_fixed(x, a):
    return mu(add(x, unit(a))) == a


def plaquette_links(x, a, b):
    return [(x, a), (add(x, unit(a)), b), (add(x, unit(b)), a), (x, b)]


def partners(x, a):
    """Distinct free links sharing a plaquette with (x, a) in Z^4, and the fixed-link count of each plaquette."""
    part, types = set(), []
    for b in range(4):
        if b == a:
            continue
        for s in (0, -1):
            o = add(x, unit(b), s)
            p = plaquette_links(o, a, b) if a < b else plaquette_links(o, b, a)
            types.append(sum(is_fixed(*g) for g in p))
            part |= {g for g in p if g != (x, a) and not is_fixed(*g)}
    return len(part), tuple(sorted(types))


def bloch_rows(k):
    """Rows (96 plaquettes of the cell) of the Bloch matrix B(k) over the 48 free links, complex entries."""
    free = [(x, a) for x in CELLS for a in range(4) if not is_fixed(x, a)]
    col = {l: i for i, l in enumerate(free)}
    rows = []
    for x in CELLS:
        for a in range(4):
            for b in range(a + 1, 4):
                r = [0] * len(free)
                for (y, c), sg in zip(plaquette_links(x, a, b), (1, 1, -1, -1)):
                    yc = tuple(v % 2 for v in y)
                    n = [(v - w) // 2 for v, w in zip(y, yc)]
                    if (yc, c) in col:
                        r[col[(yc, c)]] += sg * np.exp(1j * float(np.dot(k, n)))
                rows.append(r)
    return np.array(rows, dtype=complex)


def build_and_verify():
    links = [(x, a) for x in CELLS for a in range(4)]
    fixed_cell = [l for l in links if is_fixed(*l)]
    plaqs = [(x, a, b) for x in CELLS for a in range(4) for b in range(a + 1, 4)]
    heads = Counter(tuple(v % 2 for v in add(x, unit(a))) for x, a in fixed_cell)
    eq_forest = register(NAME, Equation(
        key="fg4_forest", entry=ENTRY, title="The period-2 forest gauge in four dimensions: 16 fixed links, one fixed link per plaquette",
        source_kind="derived", source=f"derived here from the pattern of {SRC}",
        latex_str=(r"\#\text{links}=64,\ \#\text{fixed}=16,\ \#\text{free}=48,\ \#\text{plaquettes}=96;"
                   r"%%BR%%\#\{\text{fixed links in }p\}=1\ \text{for every plaquette }p"),
        notes="Every site is the head of exactly one fixed link by the definition of the rule, and the fixed links form a "
              "forest because each parent step lowers the coordinate sum by one. What is special to this pattern is that "
              "every plaquette contains exactly one fixed link, so every plaquette has three free links; the repository "
              "checks this numerically and Lean-checks it on all of Z^4. Here it is counted over the cell, which suffices "
              "because the pattern is periodic.",
    ))
    check(NAME, "the pattern uses each of the four directions four times", Counter(PATTERN) == Counter({c: 4 for c in "0123"}), eq=eq_forest)
    check(NAME, "a period cell has 64 links, 16 fixed and 48 free, and 96 plaquettes",
          (len(links), len(fixed_cell), len(plaqs)) == (64, 16, 96), eq=eq_forest)
    check(NAME, "every one of the 16 sites is the head of exactly one fixed link", len(heads) == 16 and set(heads.values()) == {1}, eq=eq_forest)
    parent = lambda x: add(x, unit(mu(x)), -1)
    check(NAME, "the parent step lowers the coordinate sum by 1 at 1296 sample sites (no cycle)",
          all(sum(parent(x)) == sum(x) - 1 for x in itertools.product(range(-3, 6), repeat=4)), eq=eq_forest)
    cnt = Counter(sum(is_fixed(*g) for g in plaquette_links(x, a, b)) for x, a, b in plaqs)
    check(NAME, "every one of the 96 plaquettes of the cell contains exactly one fixed link", cnt == Counter({1: 96}), eq=eq_forest)
    check(NAME, "the same holds on a window of Z^4 (plaquettes with base point in [-2, 3]^4)",
          all(sum(is_fixed(*g) for g in plaquette_links(x, a, b)) == 1
              for x in itertools.product(range(-2, 4), repeat=4) for a in range(4) for b in range(a + 1, 4)), eq=eq_forest)

    # ---- partners -----------------------------------------------------------------------------
    free_cell = [l for l in links if not is_fixed(*l)]
    res = [partners(x, a) for x, a in free_cell]
    eq_part = register(NAME, Equation(
        key="fg4_partners", entry=ENTRY, title="Plaquette partners of a free link in the 4D forest gauge",
        source_kind="derived", source=f"derived here; compare {REPO}: forest4d.py, PROOF.md",
        latex_str=(r"\#\{\text{plaquettes through a free link}\}=6,\quad "
                   r"\#\{\text{free partners}\}=" + str(sorted({r[0] for r in res})[0]) + r"\ \text{to}\ " + str(sorted({r[0] for r in res})[-1])),
        notes="Counted over the 48 free links of the cell with true Z^4 coordinates. Each of the 6 plaquettes through a free "
              "link contains exactly one fixed link, so it contributes two free partners; the count of distinct partners "
              "can be smaller when plaquettes share one. Unlike the three-dimensional gauge, the four-dimensional proof does "
              "not rest on this count but on the Hessian certificate; the count is recorded for comparison with the 7 of the "
              "three-dimensional forest.",
    ))
    check(NAME, "each free link lies in 6 plaquettes, every one with exactly one fixed link",
          all(r[1] == (1,) * 6 for r in res) and len(res) == 48, eq=eq_part)
    check(NAME, "the number of distinct free partners is the same 12 for every free link (6 plaquettes x 2 others, all distinct)",
          {r[0] for r in res} == {12}, eq=eq_part)

    # ---- identity-field lower bound -----------------------------------------------------------
    B0 = sp.Matrix(np.rint(bloch_rows(np.zeros(4)).real).astype(int).tolist())
    x = sp.symbols("x")
    M = B0.T * B0
    fac = sp.factor_list(M.charpoly(x).as_expr())[1]
    quad = sp.expand(x ** 2 - 18 * x + 64)
    roots = sp.solve(quad, x)
    eq_id = register(NAME, Equation(
        key="fg4_identity", entry=ENTRY, title="The identity field in d = 4: top eigenvalue 9 + sqrt(17)",
        source_kind="derived", source=f"derived here; compare {REPO}: bloch4d.py, and {BBL}", cites=["bbl-2026"],
        latex_str=(r"\lambda_{\max}\big(B^{\mathsf T}B\big)\Big|_{k=0}=9+\sqrt{17}=13.1231\ldots,\qquad "
                   r"\lambda^2-18\lambda+64=0"),
        notes="At Q = I the Hessian is three copies of the scalar form sum_p (x_1 + x_2 - x_3 - x_4)^2 over free links "
              "(the repository's normalisation, quoted), periodic under 2 Z^4. Its Bloch matrix at k = 0 is the integer "
              "48 x 48 matrix B^T B. The grid check of other momenta is floating point. The repository and BBL both find "
              "that the identity is not the worst field in d = 4: the lower bound 13.974 comes from a hill climb and is "
              "quoted, not reproduced here.",
    ))
    check(NAME, "B is the 96 x 48 integer matrix of the cell", B0.shape == (96, 48), eq=eq_id)
    check(NAME, "the characteristic polynomial of B^T B at k = 0 has x^2 - 18x + 64 as a factor",
          any(sp.expand(f - quad) == 0 for f, _ in fac), eq=eq_id)
    check(NAME, "the largest root of that factor is 9 + sqrt(17), and no factor of B^T B has a larger root",
          max(roots, key=lambda r: float(r)) == 9 + sp.sqrt(17)
          and all(float(max(sp.Poly(f, x).nroots(n=30), key=lambda r: sp.re(r))) <= float(9 + sp.sqrt(17)) + 1e-12 for f, _ in fac), eq=eq_id)
    ks = [np.array(k) * np.pi / 2 for k in itertools.product(range(4), repeat=4)]
    grid = max(np.linalg.eigvalsh((lambda B: B.conj().T @ B)(bloch_rows(k)))[-1] for k in ks)
    check(NAME, "numerical: over the 256 momenta k in (pi/2 Z)^4 the top eigenvalue does not exceed 9 + sqrt(17)",
          grid <= float(9 + sp.sqrt(17)) + 1e-9, eq=eq_id)

    # ---- the ranges ---------------------------------------------------------------------------
    rng = lambda N, h: float(Fr(N * N) / (2 * Fr(h)))
    eq_rng = register(NAME, Equation(
        key="fg4_ranges", entry=ENTRY, title="The strong-coupling window from a Hessian bound in d = 4",
        source_kind="paper", source=f"{SRC}, README table; {BBL}", cites=["bbl-2026"],
        latex_str=(r"|\mathrm{Hess}\,S_0|\le h\ \Longrightarrow\ |\beta|<\tfrac1{2h};%%BR%%"
                   r"h=\tfrac{29}2:\ \tfrac1{29};\quad h=\tfrac{14420001}{10^6}:\ \tfrac1{28.840002};\quad h\ge13.974"),
        notes="Quoted: the three constants (BBL's 29/2, the repository's proved bound, the floating-point lower bound). "
              "Checked: the arithmetic from each constant to the window and to beta_W = N^2 beta. By Bakry-Emery this is a "
              "log-Sobolev inequality, decay of correlations and a gap at a fixed lattice spacing and strong coupling. The "
              "improvement over BBL is about 0.55% in the range, and at most about 3.2% more could be gained before the ceiling set by the lower bound. This is "
              "not the continuum limit.",
    ))
    check(NAME, "1/(2h) at h = 29/2 is 1/29 exactly", 1 / (2 * H_BBL) == Fr(1, 29), eq=eq_rng)
    check(NAME, "1/(2h) at h = 14420001/10^6 equals 1/28.840002 exactly", 1 / (2 * H_PROVED) == Fr(10 ** 6, 28840002), eq=eq_rng)
    check(NAME, "SU(2) ranges 4/(2h): 0.1379 (BBL), 0.1387 (proved), 0.1431 (ceiling)",
          [round(rng(2, h), 4) for h in (H_BBL, H_PROVED, H_LOWER)] == [0.1379, 0.1387, 0.1431], eq=eq_rng)
    check(NAME, "SU(3) ranges 9/(2h): 0.3103 (BBL), 0.3121 (proved), 0.3220 (ceiling)",
          [round(rng(3, h), 4) for h in (H_BBL, H_PROVED, H_LOWER)] == [0.3103, 0.3121, 0.3220], eq=eq_rng)
    check(NAME, "the order is lower bound < identity < proved < BBL: 13.974 > 9 + sqrt 17 and 14.420001 < 14.5",
          float(9 + sp.sqrt(17)) < float(H_LOWER) < float(H_PROVED) < float(H_BBL), eq=eq_rng)
    check(NAME, "the range gain over BBL is 0.55% (29/2 against 14.420001), and at most 3.2% more is possible above 13.974",
          abs(float(H_BBL / H_PROVED) - 1 - 0.005548) < 1e-5 and abs(float(H_PROVED / H_LOWER) - 1 - 0.031916) < 1e-5, eq=eq_rng)
