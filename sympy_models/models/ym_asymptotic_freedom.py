"""Two-loop asymptotic freedom for pure SU(N) Yang-Mills (textbook).

This is standard background for comparing the 2D O(N) scaling in on_scaling.py
with lattice Yang-Mills (see the yang-mills-mass-gap-lab fork). It is NOT from
any manuscript in this repository, which contains no gauge-theory results.
"""
import sympy as sp

from eqlib import Equation, check, register

NAME = "ym_asymptotic_freedom"


def build_and_verify():
    g, Nc = sp.symbols("g N_c", positive=True)
    a = sp.symbols("a", positive=True)  # lattice spacing
    Lam = sp.symbols(r"\Lambda_L", positive=True)

    B0, B1 = sp.symbols("b_0 b_1", positive=True)  # for display only
    b0 = 11 * Nc / (48 * sp.pi**2)
    b1 = sp.Rational(34, 3) * (Nc / (16 * sp.pi**2)) ** 2

    eq_b0 = register(NAME, Equation(
        key="ym_b0",
        title="One-loop coefficient, pure SU(N)",
        expr=sp.Eq(sp.Symbol("b_0"), b0),
        source_kind="textbook",
        source="standard; beta(g) = -b_0 g^3 - b_1 g^5 + ...",
    ))
    eq_b1 = register(NAME, Equation(
        key="ym_b1",
        title="Two-loop coefficient, pure SU(N)",
        expr=sp.Eq(sp.Symbol("b_1"), b1),
        source_kind="textbook",
        source="standard",
    ))

    # Lattice scale: a * Lambda_L = exp(-1/(2 b0 g^2)) * (b0 g^2)^(-b1/(2 b0^2))
    a_lam = sp.exp(-1 / (2 * b0 * g**2)) * (b0 * g**2) ** (-b1 / (2 * b0**2))
    eq_scale = register(NAME, Equation(
        key="ym_lattice_scale",
        title="Lattice spacing from the bare coupling, two loops",
        # Display form keeps b_0, b_1 symbolic so the LaTeX stays readable;
        # a_lam (coefficients substituted) is what the checks use.
        latex_str=(r"a\Lambda_L=\exp\!\Bigl(-\frac{1}{2b_0g^{2}}\Bigr)\,(b_0g^{2})^{-b_1/(2b_0^{2})}"),
        expr=sp.Eq(sp.Symbol(r"a\Lambda_L"),
                   sp.exp(-1 / (2 * B0 * g**2)) * (B0 * g**2) ** (-B1 / (2 * B0**2))),
        source_kind="textbook",
        source="standard two-loop result",
        notes="Universal asymptotic form only. Says nothing about existence of "
              "the continuum limit or a mass gap.",
    ))

    # Check: a dg/da = b0 g^3 + b1 g^5 (a decreases as g decreases), so
    # d ln(a Lambda)/dg = 1/(b0 g^3 + b1 g^5). Compare as a series in g.
    lhs = sp.diff(sp.log(a_lam), g)
    rhs = 1 / (b0 * g**3 + b1 * g**5)
    shown = eq_scale.expr.rhs.subs({B0: b0, B1: b1})
    check(NAME, "display form with b_0, b_1 substituted equals the checked form",
          sp.simplify(sp.log(shown) - sp.log(a_lam)) == 0, eq=eq_scale)

    # The two-loop formula is exact only through O(1/g); the residual starts at O(g)
    # (three-loop order), so compare the series strictly below g^1.
    diff = sp.series(sp.simplify(lhs - rhs), g, 0, 1).removeO()
    check(NAME, "d ln(a Lambda)/dg matches 1/(b0 g^3 + b1 g^5) through O(1/g)",
          sp.simplify(diff) == 0, eq=eq_scale)

    # Mass in lattice units: m a = c * (a Lambda) with c a constant, i.e. the same
    # exponential law. Express the leading slope for comparison with O(N).
    leading = -1 / (2 * b0)
    eq_slope = register(NAME, Equation(
        key="ym_log_slope",
        title="Leading slope of ln(a Lambda) in 1/g^2",
        expr=sp.Eq(sp.Symbol(r"\lim \frac{d\ln(a\Lambda_L)}{d(1/g^2)}"), leading),
        source_kind="derived",
        source="differentiate ym_lattice_scale",
        notes="Compare with -2 pi/(N-2) for O(N) in beta (on_log_slope).",
    ))
    x = sp.symbols("x", positive=True)  # x = 1/g^2
    slope_x = sp.diff(sp.log(a_lam.subs(g, 1 / sp.sqrt(x))), x)
    check(NAME, "slope in 1/g^2 tends to -1/(2 b0)",
          sp.simplify(sp.limit(slope_x, x, sp.oo) - leading) == 0, eq=eq_slope)
