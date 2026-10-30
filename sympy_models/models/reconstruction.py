"""The reconstruction threshold for the symmetric three-state channel (entry 229).

  [E229]  The exact reconstruction threshold for the three-state symmetric channel
          preprints/The-exact-reconstruction-threshold-for-the-three-state-symmetric-channel-September-25-2026/

Statements follow SPECTRE_CLOCK_AND_CHAIR44.md Section 4, which read the preprint text directly. Theorem 1.1:
for the channel K_lambda(i, i) = (1 + 2 lambda)/3, K_lambda(i, j) = (1 - lambda)/3 (i != j, -1/2 <= lambda <= 1)
on the regular b-ary tree (d = b) or a Poisson Galton-Watson tree of mean d > 1, the root is reconstructible
from the leaves if and only if d lambda^2 > 1, with non-reconstruction at equality. On the negative side
reconstruction occurs exactly when d > 4 and -1/2 <= lambda < -d^(-1/2). The community-detection corollary
reads (a - b)^2 > 3 (a + 2 b).

Verified here: that K_lambda is a stochastic matrix with eigenvalues 1, lambda, lambda (so lambda is the
channel's second eigenvalue, which is why the Kesten-Stigum line is d lambda^2 = 1), that its entries are
nonnegative exactly on -1/2 <= lambda <= 1, that the negative branch exists exactly when d > 4, and that the
community-detection condition is the same line under d = (a + 2b)/3, lambda = (a - b)/(a + 2b).

Not tested: the theorem, which is about the tree and its leaves, not about the matrix.
"""
import sympy as sp

from eqlib import Equation, check, register

NAME = "reconstruction"
ENTRY = "229"
E229 = "[E229] three-state reconstruction paper"


def build_and_verify():
    lam, a, b = sp.symbols("lambda a b", real=True)
    d = sp.Symbol("d", positive=True)
    K = sp.Matrix(3, 3, lambda i, j: (1 + 2 * lam) / 3 if i == j else (1 - lam) / 3)

    eq_ch = register(NAME, Equation(
        key="e229_channel", entry=ENTRY,
        title="The symmetric three-state channel",
        source_kind="paper", source=f"{E229} Theorem 1.1",
        latex_str=(r"K_\lambda(i,i)=\frac{1+2\lambda}{3},\qquad K_\lambda(i,j)=\frac{1-\lambda}{3}\ (i\ne j),"
                   r"\qquad -\tfrac12\le\lambda\le 1"),
    ))
    check(NAME, "each row of K_lambda sums to 1", sp.simplify(K * sp.ones(3, 1) - sp.ones(3, 1)) == sp.zeros(3, 1), eq=eq_ch)
    check(NAME, "the eigenvalues of K_lambda are 1, lambda and lambda",
          sp.simplify(K.eigenvals().get(1, 0)) == 1 and K.eigenvals().get(lam, 0) == 2, eq=eq_ch)
    check(NAME, "the entries are nonnegative exactly when -1/2 <= lambda <= 1",
          sp.solve_univariate_inequality((1 + 2 * lam) / 3 >= 0, lam, relational=False) == sp.Interval(-sp.Rational(1, 2), sp.oo)
          and sp.solve_univariate_inequality((1 - lam) / 3 >= 0, lam, relational=False) == sp.Interval(-sp.oo, 1), eq=eq_ch)

    eq_thm = register(NAME, Equation(
        key="e229_threshold", entry=ENTRY,
        title="Exact reconstruction threshold",
        source_kind="paper", source=f"{E229} Theorem 1.1",
        latex_str=(r"\text{reconstruction}\iff d\lambda^{2}>1,\ \text{none at equality};%%BR%%"
                   r"\text{for }\lambda<0:\ \text{reconstruction}\iff d>4\ \text{and}\ -\tfrac12\le\lambda<-d^{-1/2}"),
        notes="d is the branching number of the tree. d lambda^2 = 1 is the Kesten-Stigum line, here shown to be "
              "the exact threshold. The theorem is the paper's; the checks test the arithmetic around it.",
    ))
    check(NAME, "the positive threshold is lambda = d^(-1/2)",
          set(sp.solve(sp.Eq(d * lam**2, 1), lam)) == {-1 / sp.sqrt(d), 1 / sp.sqrt(d)}, eq=eq_thm)
    dp = d
    check(NAME, "-d^(-1/2) lies above -1/2 exactly when d > 4, so the negative branch exists only then",
          sp.solve_univariate_inequality(-1 / sp.sqrt(dp) > -sp.Rational(1, 2), dp, relational=False)
          == sp.Interval.open(4, sp.oo), eq=eq_thm)
    check(NAME, "d lambda^2 = 1 means lambda^2 equals the reciprocal of the branching number, "
                "so reconstruction needs the second eigenvalue to beat d^(-1/2)",
          sp.simplify(sp.solve(sp.Eq(d * lam**2, 1), d)[0] - 1 / lam**2) == 0, eq=eq_thm)

    eq_com = register(NAME, Equation(
        key="e229_community", entry=ENTRY,
        title="Community-detection corollary",
        source_kind="paper", source=f"{E229} abstract",
        latex_str=(r"(a-b)^{2}>3(a+2b)\quad\Longleftrightarrow\quad d\lambda^{2}>1\ \text{ with }\ "
                   r"d=\frac{a+2b}{3},\ \lambda=\frac{a-b}{a+2b}"),
        notes="The same line written for a three-community block model with within-community degree a and "
              "between-community degree b.",
    ))
    dd, ll = (a + 2 * b) / 3, (a - b) / (a + 2 * b)
    check(NAME, "d lambda^2 = (a - b)^2 / (3 (a + 2 b)) under the substitution",
          sp.simplify(dd * ll**2 - (a - b) ** 2 / (3 * (a + 2 * b))) == 0, eq=eq_com)
    check(NAME, "for a + 2b > 0, d lambda^2 > 1 is the same inequality as (a - b)^2 > 3 (a + 2b)",
          sp.simplify(sp.expand((dd * ll**2 - 1) * 3 * (a + 2 * b) - ((a - b) ** 2 - 3 * (a + 2 * b)))) == 0, eq=eq_com)
