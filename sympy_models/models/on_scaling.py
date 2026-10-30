"""Asymptotic-freedom scaling of the 2D O(N) lattice model (repo entry 215).

Equations are the ones quoted in this repository's own docs (tools/README.md,
YANGMILLS.md), not re-read from the papers. The checks are symbolic consistency
checks of those formulas against each other; they do not re-prove the papers.
"""
import sympy as sp

from eqlib import Equation, check, register

NAME = "on_scaling"


def build_and_verify():
    beta, N = sp.symbols("beta N", positive=True)
    L, H = sp.symbols("L H", positive=True)
    k = sp.symbols("k", integer=True, nonnegative=True)
    m, xi = sp.symbols("m xi", positive=True)

    # Leading asymptotic-freedom form of the mass (correlation length is 1/m).
    mass = beta ** (1 / (N - 2)) * sp.exp(-2 * sp.pi * beta / (N - 2))
    eq_mass = register(NAME, Equation(
        key="on_mass_af", entry="215",
        title="O(N) mass, leading asymptotic-freedom form",
        expr=sp.Eq(sp.Symbol("m"), mass),
        source_kind="repo-doc",
        source="tools/README.md (on_model_mc.py row), entry 215",
        notes="Shape of the scaling form only; up to an N-dependent constant. "
              "Not the transfer-matrix gap.",
    ))

    # N = 4 specialisation should match the O(4) paper's sqrt(beta) e^{-pi beta}.
    mass_o4 = sp.simplify(mass.subs(N, 4))
    target_o4 = sp.sqrt(beta) * sp.exp(-sp.pi * beta)
    eq_o4 = register(NAME, Equation(
        key="o4_gap_scaling", entry="215",
        title="O(4) transfer-gap scaling",
        expr=sp.Eq(sp.Symbol("m_{O(4)}"), sp.sqrt(beta) * sp.exp(-sp.pi * beta)),
        source_kind="repo-doc",
        source="PHYSICS.md section 1 (exact mass asymptotics, O(4))",
        notes="Paper states the gap behaves like a constant times this.",
    ))
    check(NAME, "N=4 reduces to sqrt(beta) exp(-pi beta)",
          sp.simplify(mass_o4 - target_o4) == 0, eq=eq_o4)

    # Local slope of ln m in beta, the quantity the Monte Carlo compares.
    slope = sp.simplify(sp.diff(sp.log(mass), beta))
    expected = -2 * sp.pi / (N - 2) + 1 / ((N - 2) * beta)
    eq_slope = register(NAME, Equation(
        key="on_log_slope", entry="215",
        title="Local slope of ln m against beta",
        expr=sp.Eq(sp.Derivative(sp.log(sp.Symbol("m")), beta),
                   -2 * sp.pi / (N - 2) + 1 / ((N - 2) * beta)),
        source_kind="derived",
        source="differentiate on_mass_af",
        notes="Approaches -2 pi/(N-2) from above as beta grows.",
    ))
    check(NAME, "d ln m / d beta matches closed form",
          sp.simplify(slope - expected) == 0, eq=eq_slope)
    check(NAME, "slope -> -2 pi/(N-2) as beta -> oo",
          sp.limit(expected, beta, sp.oo) == -2 * sp.pi / (N - 2), eq=eq_slope)

    # Coupling drift under exact blocking: beta_N ~ H + (log L / 2 pi) N
    # (eq. 1.4 of the O(3) paper, as quoted in YANGMILLS.md).
    beta_k = H + sp.log(L) / (2 * sp.pi) * k
    eq_drift = register(NAME, Equation(
        key="o3_coupling_drift", entry="215",
        title="Linear drift of the bare coupling under blocking (O(3))",
        expr=sp.Eq(sp.Symbol(r"\beta_k"), beta_k),
        source_kind="repo-doc",
        source="YANGMILLS.md, quoting O(3) paper eq. (1.4)",
        notes="Written with k for the number of blocking steps. Check the paper "
              "for the exact statement and error terms before citing.",
    ))
    check(NAME, "drift per step is log(L)/(2 pi)",
          sp.simplify(beta_k.subs(k, k + 1) - beta_k - sp.log(L) / (2 * sp.pi)) == 0,
          eq=eq_drift)

    # Consistency: blocking by L multiplies the correlation length by L, so with
    # beta_k growing by log(L)/(2 pi) per step, ln(xi) should drop by log(L) per step
    # in lattice units of the blocked lattice. In the N=3 form, d ln m = -2 pi dbeta
    # at leading order: -2 pi * log(L)/(2 pi) = -log(L).
    leading_dlnm = -2 * sp.pi / (3 - 2) * (sp.log(L) / (2 * sp.pi))
    check(NAME, "leading ln m change per blocking step equals -log L (N=3)",
          sp.simplify(leading_dlnm + sp.log(L)) == 0, eq=eq_drift)
