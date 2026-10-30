"""The doubling defect of the O(3) paper (entry 215) and its comparison with Chair44.

  [O3]  The canonical massive continuum limit of the two-dimensional O(3) model
        preprints/The-canonical-massive-continuum-limit-of-the-two-dimensional-O3-model-October-4-2026/

Statements follow SPECTRE_CLOCK_AND_CHAIR44.md Section 2, which read the preprint text directly. For a
spatial circumference w, [O3] Section 8.2 defines Z(n, w) = Tr K_w^n, the one-direction defects
p(n, w) = 2 log Z(n, w) - log Z(2n, w) and q(n, w) = 2 log Z(n, w) - log Z(n, 2w), and the doubling defect
Delta(n, w) = 2 p(n, w) + q(2n, w). Lemma 8.2 bounds p(2n, w) <= 4 p(n, w)^2 and Proposition 8.3 turns the
quadratic recursion Delta(2n, 2w) <= 64 Delta(n, w)^2 into a doubly exponential decay.

Verified here: the identity Delta = 4 log Z(n, w) - log Z(2n, 2w) (the fork's algebra from (8.5) and (8.7));
that the closed form in Proposition 8.3 solves the recursion and implies the weaker stated bound; and that,
against the Chair44 contact defect of the model spectre_chair44, the O(3) bound is eventually smaller than
any geometric rate (so the analogy of the notes is one of form and not of rate).

Not tested: Lemma 8.2 and Proposition 8.3 themselves (they concern a positive transfer operator, not a
number), and the relation of the transfer operator to the O(3) model.
"""
import sympy as sp

from eqlib import Equation, check, register

NAME = "doubling_defect"
ENTRY = "215"
O3 = "[O3] O(3) continuum-limit paper"


def build_and_verify():
    L = sp.Function("L")           # L(n, w) = log Z(n, w)
    n, w = sp.symbols("n w", positive=True)
    p = lambda a, b: 2 * L(a, b) - L(2 * a, b)
    q = lambda a, b: 2 * L(a, b) - L(a, 2 * b)
    Delta = 2 * p(n, w) + q(2 * n, w)

    eq_def = register(NAME, Equation(
        key="o3_doubling_defect", entry=ENTRY,
        title="The doubling defect",
        source_kind="paper", source=f"{O3} Section 8.2, (8.5) and (8.7)",
        latex_str=(r"p(n,w)=2\log Z(n,w)-\log Z(2n,w),\qquad q(n,w)=2\log Z(n,w)-\log Z(n,2w),%%BR%%"
                   r"\Delta(n,w)=2p(n,w)+q(2n,w)=4\log Z(n,w)-\log Z(2n,2w)\ \ge 0"),
        notes="Z(n, w) = Tr K_w^n is the trace of the n-th power of a positive transfer operator at spatial "
              "circumference w. Delta is the extensivity defect under doubling both sides; it is nonnegative "
              "on a torus. The last equality is the fork's algebra, checked here.",
    ))
    check(NAME, "Delta = 2 p(n, w) + q(2n, w) equals 4 log Z(n, w) - log Z(2n, 2w)",
          sp.expand(Delta - (4 * L(n, w) - L(2 * n, 2 * w))) == 0, eq=eq_def)
    check(NAME, "the cross terms log Z(2n, w) cancel between p and q", sp.expand(Delta).coeff(L(2 * n, w)) == 0, eq=eq_def)

    # Proposition 8.3: Delta(n, w) <= 2^-10 => Delta(2^k n, 2^k w) <= 64^-1 2^(-2^k)
    k = sp.symbols("k", nonnegative=True, integer=True)
    a = lambda kk: sp.Rational(1, 64) * (64 * sp.Rational(1, 2**10)) ** (2**kk)
    eq_prop = register(NAME, Equation(
        key="o3_prop_8_3", entry=ENTRY,
        title="Doubly exponential decay of the doubling defect",
        source_kind="paper", source=f"{O3} Proposition 8.3 (with Lemma 8.2)",
        latex_str=(r"\Delta(n,w)\le 2^{-10}\ (n,w\ge4)\ \Longrightarrow\ "
                   r"\Delta(2^kn,2^kw)\le \frac{1}{64}\,2^{-2^k}\quad(k\ge 0),"
                   r"%%BR%%\text{from}\quad\Delta(2n,2w)\le 64\,\Delta(n,w)^{2}"),
        notes="The proof iterates the quadratic recursion to 64^{-1}(64 * 2^{-10})^{2^k} and weakens it to the "
              "stated bound. The theorem and the recursion are the paper's; only the iteration is checked here.",
    ))
    check(NAME, "a_0 = 2^-10 (the hypothesis) and a_{k+1} = 64 a_k^2 for a_k = (64 2^-10)^(2^k) / 64",
          a(0) == sp.Rational(1, 2**10) and sp.simplify(a(k + 1) - 64 * a(k) ** 2) == 0, eq=eq_prop)
    check(NAME, "the iterate is 2^(-6 - 4 * 2^k), below the stated bound 2^(-6 - 2^k) for every k >= 0",
          sp.simplify(a(k) - 2 ** (-6 - 4 * 2**k)) == 0 and sp.simplify(4 * 2**k - 2**k - 3 * 2**k) == 0
          and all(a(kk) <= sp.Rational(1, 64) * sp.Rational(1, 2) ** (2**kk) for kk in range(12)), eq=eq_prop)

    # comparison with the Chair44 contact defect (model spectre_chair44): relative defect (12/29) 2^-n
    m = sp.symbols("m", positive=True)
    chair = sp.Rational(12, 29) * 2 ** (-m)
    bound = sp.Rational(1, 64) * 2 ** (-(2**m))
    eq_cmp = register(NAME, Equation(
        key="o3_vs_chair44_rates", entry=ENTRY,
        title="Rates of the two doubling defects",
        source_kind="derived", source="o3_prop_8_3 against chair44_contact_formula",
        latex_str=(r"\frac{\tfrac1{64}\,2^{-2^{m}}}{\tfrac{12}{29}\,2^{-m}}\ \longrightarrow\ 0"
                   r"\quad(m\to\infty)"),
        notes="Both structures have an exact coarse-graining step and a boundary-interaction difficulty, but the "
              "O(3) defect obeys a quadratic recursion (doubly exponential bound) and the Chair44 contact count "
              "has a defect of relative size 12/29 2^-n. Nothing transfers from one to the other; the analogy is "
              "in form. The sign also differs: the torus defect is >= 0, the chair defect arises from a free boundary.",
    ))
    check(NAME, "the ratio of the O(3) bound at doubling number m to the Chair44 relative defect at level m tends to 0",
          sp.limit(sp.simplify(bound / chair), m, sp.oo) == 0, eq=eq_cmp)
    check(NAME, "at m = 4 the O(3) bound is already smaller than the Chair44 defect by a factor of more than 10^3",
          sp.N(chair.subs(m, 4) / bound.subs(m, 4)) > 1e3, eq=eq_cmp)
