"""Heisenberg/ClockTails.lean of the OpenAI release in the V2 syntax: a category with explicit Lean statements,
connected items, and LaTeX prose built from the items themselves. (heisenberg_clocks.py is the V1 form.)"""
from leandsl import lean, Nat, Real, require, binom, fact, E, Rational, wikipedia, chain, le, eq
import sympy as sp
from oai_tree import family   # the generated class of family 271 (the quantum Heisenberg ferromagnet) of the release

# one reading of r for every item here; connect() refuses items that read a letter differently
R_READING = "the exponent, which is also the factorial's argument and the binomial's lower index"


@lean.category(
    r'''
    \textbf{Heisenberg/ClockTails.lean} of the OpenAI release, written in the leandsl DSL. It is a module behind the Lean
    proof of the family's main theorem (the import analysis finds it in the family's closure).
    ''',
    title="Heisenberg clock tails",
)
class HeisenbergClockTails(family("271")):
    r'''The file declares the lemma \texttt{factorial\_lower}, a Stirling-type lower bound for $r!$ (its Lean proof uses
    \texttt{Stirling.le\_factorial\_stirling}), and the theorem \texttt{clock\_tail}, whose proof applies it.'''


@HeisenbergClockTails.lemmas(
    lean='lemma factorial_lower (r : ℕ) : ((r:ℝ)/Real.exp 1)^r ≤ (r.factorial:ℝ)',
    source="OAI.Heisenberg.Clocks.factorial_lower",
    doc="Stirling-type lower bound: $(r/e)^r \\le r!$.",
    refs=[wikipedia("Stirling's approximation")],
    reads={"r": R_READING},
)
def factorial_lower(r: Nat):
    r'''The bound $(r/e)^r\le r!$, for every natural number $r$.'''
    ratio = sp.Mul(r, sp.Pow(E, -1, evaluate=False), evaluate=False)
    return sp.Pow(ratio, r, evaluate=False) <= fact(r)


@HeisenbergClockTails.theorems(
    lean='''theorem clock_tail {N r : ℕ} {B δ : ℝ} (hB : 0≤B) (hδ : 0≤δ) (hr : 0<r)
        (hR : 4*Real.exp 1*B*((N+r:ℕ):ℝ)*δ ≤ r) :
        (((N+r).choose r:ℝ)*(2*B*δ)^r) ≤ (1/2:ℝ)^r''',
    source="OAI.Heisenberg.Clocks.clock_tail",
    doc="Tail bound: if $4eB(N+r)\\delta \\le r$ then $\\binom{N+r}{r}(2B\\delta)^r \\le 2^{-r}$.",
    refs=[wikipedia("Binomial coefficient")],
    reads={"N": "the binomial's upper index minus its lower index", "r": R_READING},
)
def clock_tail(N: Nat, r: Nat, B: Real, delta: Real):
    r'''The proof in the Lean file is a \texttt{calc} chain; it is reproduced and checked step by step below.'''
    require(B >= 0, delta >= 0, r > 0, 4 * E * B * (N + r) * delta <= r)
    return binom(N + r, r) * (2 * B * delta) ** r <= Rational(1, 2) ** r


# connecting concepts: later used to generate LaTeX that explains how the pieces fit
HeisenbergClockTails.theorems.clock_tail.connect(HeisenbergClockTails.lemmas.factorial_lower)



# The calc of `clock_tail` in lean/OAI/MathematicalPhysics/Heisenberg/ClockTails.lean, step for step, with the Lean lemma
# each step uses. Every step is decided by sympy where it can be and otherwise sampled at points that satisfy the
# hypotheses; the printed document says which, and a false step makes this file fail to load. Whether Lean accepts
# the proof is a separate matter: nothing here has been compiled.
@HeisenbergClockTails.theorems.clock_tail.chain(title="The proof is this chain of inequalities:")
def clock_tail_chain(N, r, B, delta):
    n = N + r
    return chain(
        binom(n, r) * (2 * B * delta) ** r,
        le(n ** r / fact(r) * (2 * B * delta) ** r, r"\texttt{Nat.choose\_le\_pow\_div}: $\binom{n}{r}\le n^r/r!$"),
        eq((2 * B * delta * n) ** r / fact(r), r"\texttt{ring}"),
        le((2 * B * delta * n) ** r / (r / E) ** r, r"$(r/e)^r\le r!$",
           using=HeisenbergClockTails.lemmas.factorial_lower),
        eq((2 * B * delta * n / (r / E)) ** r, r"\texttt{div\_pow}"),
        le(Rational(1, 2) ** r,
           r"$2B\delta n/(r/e)\le 1/2$ from $4eB\delta n\le r$, and \texttt{pow\_le\_pow\_left}"),
    )

# a category that depends on another one: lookups see the inherited items
@lean.category(
    r'''
    \section{ Summary }
    The tail bound %s is built on the lemma %s.
    ''' % (HeisenbergClockTails.theorems.clock_tail.toLatex(), HeisenbergClockTails.lemmas.factorial_lower.toLatex()),
)
class ClockTailsSummary(HeisenbergClockTails):
    pass
