"""Heisenberg/ClockTails.lean of the OpenAI release, written in the leandsl DSL.

Original Lean:
    lemma factorial_lower (r : ℕ) : ((r:ℝ)/Real.exp 1)^r ≤ (r.factorial:ℝ)
    theorem clock_tail {N r : ℕ} {B δ : ℝ} (hB : 0≤B) (hδ : 0≤δ) (hr : 0<r)
        (hR : 4*Real.exp 1*B*((N+r:ℕ):ℝ)*δ ≤ r) :
        (((N+r).choose r:ℝ)*(2*B*δ)^r) ≤ (1/2:ℝ)^r
"""
from leandsl import lean, Nat, Real, require, binom, fact, E, Rational
import sympy as sp


@lean.theorem(source="OAI.MathematicalPhysics.Heisenberg.Clocks.factorial_lower",
              doc="Stirling-type lower bound: $(r/e)^r \\le r!$.")
def factorial_lower(r: Nat):
    ratio = sp.Mul(r, sp.Pow(E, -1, evaluate=False), evaluate=False)
    return sp.Pow(ratio, r, evaluate=False) <= fact(r)


@lean.theorem(source="OAI.MathematicalPhysics.Heisenberg.Clocks.clock_tail",
              doc="Tail bound: if $4eB(N+r)\\delta \\le r$ then $\\binom{N+r}{r}(2B\\delta)^r \\le 2^{-r}$.")
def clock_tail(N: Nat, r: Nat, B: Real, delta: Real):
    require(B >= 0, delta >= 0, r > 0, 4 * E * B * (N + r) * delta <= r)
    return binom(N + r, r) * (2 * B * delta) ** r <= Rational(1, 2) ** r
