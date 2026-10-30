"""A translational tile with no fully periodic tiling in dimension three (entry 155).

Statements are transcribed from
  [T]  A translational tile with no fully periodic tiling in dimension three
       preprints/A-translational-tile-with-no-fully-periodic-tiling-in-dimension-three-September-23-2026/

The proof builds a finite tile in Z^3 from a two-digit p-adic Sudoku rule (Section 2), encodes
it by simultaneous tiling equations on a cyclic finite factor (Sections 3-5), stacks (Section 6)
and passes to Z^3 and R^3 (Section 7). This module models the combinatorial core of Section 2:
the line rule, the last-nonzero-digit model of Lemma 2.5, and the facts about it that the proof
uses. The tiling construction itself (Sections 3-7) is NOT modelled.

The paper fixes a prime p > 200. The finite checks below run at p = 3 and p = 5 and test the
mechanism (they are independent enumerations of the definitions, not the paper's proofs), so
they illustrate Lemma 2.5 and the vertical-period obstruction at small p. They say nothing
about the p > 200 construction.
"""
import itertools

import sympy as sp

from eqlib import Equation, check, register

NAME = "periodic_tiling"
ENTRY = "155"
T = "[T] translational tile paper"


def _last_nonzero_digit(t, p):
    """f(t) = t / p^{ord_p t} mod p, with f(0) = 1  (Lemma 2.5)."""
    if t == 0:
        return 1
    return (t // p ** int(sp.multiplicity(p, t))) % p


def _patterns(p):
    """All allowed-word patterns of Definition 2.1 over indices 0..p^2-1; None marks a free symbol."""
    N = p * p
    pats = set()
    for a, b in itertools.product(range(N), repeat=2):
        if a % p == 0 and b % p == 0:
            continue
        word = []
        for n in range(N):
            tv = a * n + b
            if tv % p != 0:
                word.append(tv % p)
            elif tv % (p * p) != 0:
                word.append((tv // p) % p)
            else:
                word.append(None)
        pats.add(tuple(word))
    return list(pats)


def _matches(word, pat):
    return all(q is None or q == w for w, q in zip(word, pat))


def build_and_verify():
    # ---- statements -----------------------------------------------------------------------
    register(NAME, Equation(
        key="tiling_sum_definition", entry=ENTRY,
        title="Translational tile and full periodicity",
        source_kind="paper", source=f"{T}, Section 1",
        latex_str=(r"A\oplus F=G\iff \text{every }g\in G\text{ is uniquely }a+f\ (a\in A,\ f\in F);%%BR%%"
                   r"A\ \text{fully periodic}\iff A+H=A\ \text{for a finite-index subgroup }H\le G"),
    ))
    eq_thm = register(NAME, Equation(
        key="tiling_theorem", entry=ENTRY,
        title="Counterexample to the periodic tiling conjecture in dimension three",
        source_kind="paper", source=f"{T}, Theorem 1.1",
        latex_str=(r"\exists\,\text{finite }T\subset\mathbb{Z}^3:\ \ T\oplus A=\mathbb{Z}^3\text{ for some }A,\ "
                   r"\text{but no such }A\text{ is invariant under a finite-index subgroup;}%%BR%%"
                   r"\Omega=T+[0,1]^3\ \text{tiles }\mathbb{R}^3\text{ with no full-rank-lattice-periodic tiling}"),
        notes="The Euclidean assertion allows arbitrary real translation vectors. Full periodicity is the "
              "excluded property; neither connectedness of the tile nor absence of every nonzero period is asserted.",
    ))

    p_sym = sp.Symbol("p", positive=True)
    eq_rule = register(NAME, Equation(
        key="tiling_line_rule", entry=ENTRY,
        title="Allowed words and the line rule (two-digit p-adic Sudoku rule)",
        source_kind="paper", source=f"{T}, Definition 2.1",
        latex_str=(r"N=p^2,\ \ \Sigma=\mathbb{F}_p^{\times},\ \ I=\{0,\dots,N-1\};%%BR%%"
                   r"w(n)=\begin{cases}an+b \pmod p,& p\nmid an+b,\\ (an+b)/p \pmod p,& p\mid an+b,\ p^2\nmid an+b,\end{cases}"
                   r"%%BR%%\text{free if }p^2\mid an+b,\ \ p\nmid\gcd(a,b);%%BR%%"
                   r"W\ \text{satisfies the line rule if }n\mapsto W(n,dn+e)\text{ is allowed for all }d,e\in\mathbb{Z}"),
        notes="A prime p > 200 is fixed in the paper.",
    ))
    eq_f = register(NAME, Equation(
        key="tiling_last_nonzero_digit", entry=ENTRY,
        title="Last-nonzero-digit model of the line rule",
        source_kind="paper", source=f"{T}, Lemma 2.5",
        latex_str=(r"f(t)=\frac{t}{p^{\operatorname{ord}_p(t)}}\pmod p\ \ (t\neq0),\quad f(0)=1;"
                   r"%%BR%%W(n,m)=f(m)\ \text{satisfies the line rule with every column nonconstant}"),
    ))
    eq_obs = register(NAME, Equation(
        key="tiling_vertical_period_obstruction", entry=ENTRY,
        title="No vertical period when every column is nonconstant",
        source_kind="paper", source=f"{T}, Proposition 2.4",
        latex_str=(r"\nexists\,W:I\times\mathbb{Z}\to\Sigma\ \text{with the line rule, every column nonconstant,}%%BR%%"
                   r"\text{and}\ W(n,m+M)=W(n,m)\ \text{for some }M>0"),
        notes="A fully periodic tiling would decode to such an array, which is the contradiction.",
    ))
    register(NAME, Equation(
        key="tiling_line_parametrization", entry=ENTRY,
        title="Parametrization of the lines",
        source_kind="paper", source=f"{T}, eq. (2.1)",
        latex_str=r"L_n(x)=x_2+nx_1,\qquad u_n=(1,-n)\qquad(n\in I)",
        notes="L_n: Z^2 -> Z is onto with kernel generated by u_n.",
    ))

    # ---- finite checks at small primes -----------------------------------------------------
    # Check that L_n is onto Z with kernel <u_n> is symbolic and prime-free.
    x1, x2, nn = sp.symbols("x1 x2 n", integer=True)
    check(NAME, "L_n(u_n) = 0 for u_n = (1,-n)", sp.expand((x2 + nn * x1).subs({x1: 1, x2: -nn})) == 0, eq=eq_f)

    for p in (3, 5):
        N = p * p
        # f is well defined with values in F_p^x and constant along p-power scalings
        vals = {_last_nonzero_digit(t, p) for t in range(-p**4, p**4 + 1)}
        check(NAME, f"p={p}: f takes values exactly in F_p^x = {{1..{p - 1}}} on [-p^4, p^4]",
              vals == set(range(1, p)), eq=eq_f)
        check(NAME, f"p={p}: f(p^r t) = f(t) for r = 0..3 and |t| <= p^3",
              all(_last_nonzero_digit(p**r * t, p) == _last_nonzero_digit(t, p)
                  for r in range(4) for t in range(-p**3, p**3 + 1)), eq=eq_f)
        check(NAME, f"p={p}: f(1)=1 and f(2)=2, so every column of W(n,m)=f(m) is nonconstant",
              _last_nonzero_digit(1, p) == 1 and _last_nonzero_digit(2, p) == 2, eq=eq_f)

        # Line rule for the model, tested against an independent enumeration of ALL allowed words.
        pats = _patterns(p)
        R = p * p
        bad = []
        for d, e in itertools.product(range(-R, R + 1), repeat=2):
            word = [_last_nonzero_digit(d * n + e, p) for n in range(N)]
            if not any(_matches(word, pt) for pt in pats):
                bad.append((d, e))
        check(NAME, f"p={p}: W(n,m)=f(m) satisfies the line rule for all |d|,|e| <= {R} "
                    f"({len(pats)} allowed patterns enumerated)", not bad, eq=eq_f)

        # Model has no vertical period: for each M up to p^4 find a witness m with f(m+M) != f(m)
        no_period = all(any(_last_nonzero_digit(m + M, p) != _last_nonzero_digit(m, p)
                            for m in range(0, p**5)) for M in range(1, p**4 + 1))
        check(NAME, f"p={p}: no vertical period M <= p^4 for the model array (consistent with Proposition 2.4)",
              no_period, eq=eq_obs)

        # Exceptional indices: a nonzero affine form over F_p has at most p zeros on I = {0..p^2-1}
        worst = 0
        for a, b in itertools.product(range(p), repeat=2):
            if a == 0 and b == 0:
                continue
            zeros = sum(1 for n in range(N) if (a * n + b) % p == 0)
            worst = max(worst, zeros)
        check(NAME, f"p={p}: a nonzero affine form has at most N/p = p exceptional indices on I (max {worst})",
              worst <= p, eq=eq_rule)

        # Remark after (2.1): x1 = 0 mod p and x2 = t mod p (t in {1,2}) forces every word entry to equal t
        ok = True
        for t_ in (1, 2):
            for x1v, x2v in itertools.product(range(-3 * p, 3 * p + 1), repeat=2):
                if x1v % p == 0 and x2v % p == t_ % p:
                    if any(_last_nonzero_digit(x2v + n * x1v, p) != t_ for n in range(N)):
                        ok = False
        check(NAME, f"p={p}: if x1 = 0 and x2 = t (mod p), t in {{1,2}}, then f(L_n(x)) = t for every n",
              ok, eq=eq_f)
