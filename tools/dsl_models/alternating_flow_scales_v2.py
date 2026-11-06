"""MathematicalPhysics/AlternatingFlow/Scales.lean of the OpenAI release in the V2 syntax: a category with explicit Lean
statements, connected items, and LaTeX prose built from the items themselves.

The file belongs to the closure of family 376 (universal computation in forced Navier-Stokes flows) through the
alternating-flow construction. It chooses the scales of that construction: stage n has the capacity K = N+2n+2 and the
scale exponent s = (N+1)2^((n+3)^2), so that the exponents leave room for everything stored at the earlier stages, and the
sizes eps_n = b^(-s_n) of the perturbations are summable with a small total.

Everything the DSL can state (statements about natural numbers) is an item with a body; the statements about the real
numbers eps_n = (b^(s N n))^(-1), whose exponent is far too large to evaluate, and the infinite sums are results recorded
with their Lean text. The Lean names use the Greek letter epsilon; the DSL names of the results spell it `eps`.

Not modelled: the other modules of AlternatingFlow."""
from leandsl import lean, Nat, Real, require, Rational, wikipedia, chain, le, lt, eq
import sympy as sp
from oai_tree import family   # the generated class of family 376 of the release

# one reading of each letter for every item here; connect() refuses items that read a letter differently
N_READING = "the offset of the capacities, a natural number: $K_n=N+2n+2$"
I_READING = "the stage of the construction, a natural number"
R_READING = "the order of the derivative that is paid for, a natural number at most the stage"
B_READING = "the base of the perturbation sizes, a natural number: $\\varepsilon_n=b^{-s_n}$"

_MAIN = family("376").results["OAI.AlternatingNS.alternating"]   # the family's Comparator challenge, a generated Item


def _K(N, n):
    return N + 2 * n + 2


def _s(N, n):
    return (N + 1) * 2**((n + 3)**2)


@lean.category(
    r'''
    \textbf{MathematicalPhysics/AlternatingFlow/Scales.lean} of the OpenAI release, written in the leandsl DSL. It is a module
    in the closure of the Lean proof of the family's main theorem %s, through the alternating-flow construction; the family
    itself is [[F376|family 376]]. Stage $n$ has the capacity $K_n=N+2n+2$ and the scale exponent
    $s_n=(N+1)\,2^{(n+3)^2}$, and the perturbations have the sizes $\varepsilon_n=b^{-s_n}$. The scales are chosen so that
    they separate: $s_{n+1}=2^{2n+7}s_n$ ([[AlternatingFlowScales.s_succ]]), a quadratic exponent that is large
    ([[AlternatingFlowScales.pow_lower]], [[AlternatingFlowScales.s_lower]]) and therefore larger than the capacity
    ([[AlternatingFlowScales.capacity_lt]]). From that follow the gap between stages
    ([[AlternatingFlowScales.scale_gap]]) and the gap that pays for $r$ derivatives
    ([[AlternatingFlowScales.derivative_gap]]); the sizes $\varepsilon_n$ then decay at least geometrically
    ([[AlternatingFlowScales.eps_le_geometric]]), their sum is below $2\varepsilon_0<1/4$
    ([[AlternatingFlowScales.total_memory_le]], [[AlternatingFlowScales.initial_memory_small]]), and a derivative costs
    less than the size of the previous stage ([[AlternatingFlowScales.derivative_cost_bound]]).
    ''' % _MAIN.latex_ref(page=True),
    title="Alternating flow scales",
    refs=[wikipedia("Geometric series"), wikipedia("Navier–Stokes equations")],
)
class AlternatingFlowScales(family("376")):
    r'''The file defines $K_n=N+2n+2$, $s_n=(N+1)2^{(n+3)^2}$ and $\varepsilon_n=b^{-s_n}$.'''


@AlternatingFlowScales.lemmas(
    lean='''lemma s_succ (N n : ℕ) : s N (n + 1) = 2 ^ (2 * n + 7) * s N n''',
    source="OAI.AlternatingNS.Scales.s_succ",
    doc="One stage later the scale exponent is $2^{2n+7}$ times larger: $s_{n+1}=2^{2n+7}s_n$.",
    reads={"N": N_READING, "n": I_READING},
)
def s_succ(N: Nat, n: Nat):
    r'''Because $(n+4)^2=(n+3)^2+(2n+7)$.'''
    return sp.Eq(_s(N, n + 1), 2**(2 * n + 7) * _s(N, n), evaluate=False)


@AlternatingFlowScales.lemmas(
    lean='''lemma pow_lower (n : ℕ) : 8 * (n + 1) ≤ 2 ^ ((n + 3) ^ 2)''',
    source="OAI.AlternatingNS.Scales.pow_lower",
    doc="The quadratic exponent is large: $8(n+1)\\le2^{(n+3)^2}$.",
    reads={"n": I_READING},
)
def pow_lower(n: Nat):
    r'''Because $n+1\le2^n$ and $n+3\le(n+3)^2$, so $8(n+1)\le8\cdot2^n=2^{n+3}\le2^{(n+3)^2}$.'''
    return 8 * (n + 1) <= 2**((n + 3)**2)


@AlternatingFlowScales.lemmas(
    lean='''lemma s_lower (N n : ℕ) : 8 * (N + 1) * (n + 1) ≤ s N n''',
    source="OAI.AlternatingNS.Scales.s_lower",
    doc="The scale exponent is at least $8(N+1)(n+1)$.",
    reads={"N": N_READING, "n": I_READING},
)
def s_lower(N: Nat, n: Nat):
    r'''Multiply \texttt{pow\_lower} by $N+1$.'''
    return 8 * (N + 1) * (n + 1) <= _s(N, n)


@AlternatingFlowScales.lemmas(
    lean='''lemma capacity_lt (N n : ℕ) : 1 + 2 * K N n < s N n''',
    source="OAI.AlternatingNS.Scales.capacity_lt",
    doc="The scale exponent exceeds what the capacity needs: $1+2K_n<s_n$.",
    reads={"N": N_READING, "n": I_READING},
)
def capacity_lt(N: Nat, n: Nat):
    r'''$1+2K_n=2N+4n+5$ is less than $8(N+1)(n+1)$, which is at most $s_n$ by \texttt{s\_lower}.'''
    return 1 + 2 * _K(N, n) < _s(N, n)


@AlternatingFlowScales.lemmas(
    lean='''lemma scale_gap (N n : ℕ) : s N n + 1 + 2 * K N n < s N (n + 1)''',
    source="OAI.AlternatingNS.Scales.scale_gap",
    doc="Consecutive scale exponents are separated by more than the capacity: $s_n+1+2K_n<s_{n+1}$.",
    reads={"N": N_READING, "n": I_READING},
)
def scale_gap(N: Nat, n: Nat):
    r'''Since $1+2K_n<s_n$ the left side is below $2s_n$, which is at most $2^{2n+7}s_n=s_{n+1}$.'''
    return _s(N, n) + 1 + 2 * _K(N, n) < _s(N, n + 1)


@AlternatingFlowScales.lemmas(
    lean='''lemma derivative_gap (N r n : ℕ) (hn : r ≤ n) :
    s N n + (r + 1) * (s N n + 1 + K N n) ≤ s N (n + 1)''',
    source="OAI.AlternatingNS.Scales.derivative_gap",
    doc="The gap between stages also pays for $r\\le n$ derivatives: $s_n+(r+1)(s_n+1+K_n)\\le s_{n+1}$.",
    reads={"N": N_READING, "r": R_READING, "n": I_READING},
)
def derivative_gap(N: Nat, r: Nat, n: Nat):
    r'''Since $K_n+1<s_n$ each of the $r+1$ terms is below $2s_n$, so the left side is below $(2r+3)s_n\le2^{2n+7}s_n=s_{n+1}$,
    because $2r+3\le2n+7<2^{2n+7}$.'''
    require(r <= n)
    return _s(N, n) + (r + 1) * (_s(N, n) + 1 + _K(N, n)) <= _s(N, n + 1)


AlternatingFlowScales.results.declare(
    "strictMono_s",
    lean='''lemma strictMono_s (N : ℕ) : StrictMono (s N)''',
    doc=r"The scale exponents increase strictly: $s_{n+1}=2^{2n+7}s_n>s_n$ because $s_n>0$.",
    source="OAI.AlternatingNS.Scales.strictMono_s",
)

AlternatingFlowScales.results.declare(
    "eps_le_geometric",
    lean='''lemma ε_le_geometric (b N n : ℕ) (hb : 2 ≤ b) :
    ε b N n ≤ (1 / 2 : ℝ) ^ n''',
    doc=r"The sizes decay at least geometrically: $b^{-s_n}\le2^{-s_n}\le2^{-n}$, because $n\le s_n$ (\texttt{n\_le\_s}).",
    source="OAI.AlternatingNS.Scales.ε_le_geometric",
)

AlternatingFlowScales.results.declare(
    "total_memory_le",
    lean='''lemma total_memory_le (b N : ℕ) (hb : 2 ≤ b) :
    (∑' n : ℕ, ε b N n) ≤ 2 * ε b N 0''',
    doc=r"The total size is at most twice the first: $\varepsilon_n\le\varepsilon_0 2^{-n}$ (since $s_n\ge s_0+n$) and "
        r"$\sum_n2^{-n}=2$.",
    source="OAI.AlternatingNS.Scales.total_memory_le",
    refs=[wikipedia("Geometric series")],
)

AlternatingFlowScales.results.declare(
    "initial_memory_small",
    lean='''lemma initial_memory_small (b N : ℕ) (hb : 4 ≤ b) : 2 * ε b N 0 < 1 / 4''',
    doc=r"With $b\ge4$ the first size is small: $s_0\ge2$ gives $b^{s_0}\ge16$, so $2\varepsilon_0\le1/8<1/4$.",
    source="OAI.AlternatingNS.Scales.initial_memory_small",
)

AlternatingFlowScales.results.declare(
    "derivative_cost_bound",
    lean='''lemma derivative_cost_bound (b N r n : ℕ) (hb : 2 ≤ b) (hn : r ≤ n) :
    ε b N (n + 1) * (b : ℝ) ^ ((r + 1) * (s N n + 1 + K N n)) ≤ ε b N n''',
    doc=r"Paying for $r\le n$ derivatives at stage $n+1$ costs less than the size of stage $n$: "
        r"\texttt{derivative\_gap} says $(r+1)(s_n+1+K_n)\le s_{n+1}-s_n$.",
    source="OAI.AlternatingNS.Scales.derivative_cost_bound",
)

AlternatingFlowScales.results.declare(
    "weighted_eps_bounded",
    lean='''lemma weighted_ε_bounded (b N J : ℕ) (hb : 2 ≤ b) :
    ∃ C : ℝ, 0 ≤ C ∧ ∀ n : ℕ, (1 + (n : ℝ)) ^ J * ε b N n ≤ C''',
    doc=r"Every polynomial weight is beaten: $(1+n)^J\varepsilon_n$ is summable, so it is bounded (and tends to $0$).",
    source="OAI.AlternatingNS.Scales.weighted_ε_bounded",
)


# connecting concepts: later used to generate LaTeX that explains how the pieces fit
AlternatingFlowScales.lemmas.s_lower.connect(AlternatingFlowScales.lemmas.pow_lower)
AlternatingFlowScales.lemmas.capacity_lt.connect(AlternatingFlowScales.lemmas.s_lower)
AlternatingFlowScales.lemmas.scale_gap.connect(AlternatingFlowScales.lemmas.s_succ, AlternatingFlowScales.lemmas.capacity_lt)
AlternatingFlowScales.lemmas.derivative_gap.connect(AlternatingFlowScales.lemmas.s_succ)
AlternatingFlowScales.results.derivative_cost_bound.connect(AlternatingFlowScales.lemmas.derivative_gap)
AlternatingFlowScales.results.total_memory_le.connect(AlternatingFlowScales.results.eps_le_geometric)
AlternatingFlowScales.results.initial_memory_small.connect(AlternatingFlowScales.lemmas.s_lower)


@AlternatingFlowScales.lemmas.pow_lower.chain(title="The proof is this chain:")
def pow_lower_chain(n):
    return chain(
        8 * (n + 1),
        le(8 * 2**n, r"$n+1\le2^n$ (\texttt{Nat.lt\_two\_pow\_self})"),
        eq(2**(n + 3), r"\texttt{pow\_add}"),
        le(2**((n + 3)**2), r"$n+3\le(n+3)^2$ (\texttt{nlinarith}) and \texttt{Nat.pow\_le\_pow\_right}"),
    )


@AlternatingFlowScales.lemmas.capacity_lt.chain(title="The proof is this chain:")
def capacity_lt_chain(N, n):
    return chain(
        1 + 2 * _K(N, n),
        lt(8 * (N + 1) * (n + 1), r"$2N+4n+5<8Nn+8N+8n+8$ (\texttt{nlinarith})"),
        le(_s(N, n), r"\texttt{s\_lower}", using=AlternatingFlowScales.lemmas.s_lower),
    )


@AlternatingFlowScales.lemmas.scale_gap.chain(title="The proof is this chain:")
def scale_gap_chain(N, n):
    return chain(
        _s(N, n) + 1 + 2 * _K(N, n),
        lt(2 * _s(N, n), r"$1+2K_n<s_n$", using=AlternatingFlowScales.lemmas.capacity_lt),
        le(2**(2 * n + 7) * _s(N, n), r"$2\le2^{2n+7}$ (\texttt{Nat.le\_pow})"),
        eq(_s(N, n + 1), r"\texttt{s\_succ}", using=AlternatingFlowScales.lemmas.s_succ),
    )


@AlternatingFlowScales.lemmas.derivative_gap.chain(title="The proof is this chain:")
def derivative_gap_chain(N, r, n):
    s = _s(N, n)
    return chain(
        s + (r + 1) * (s + 1 + _K(N, n)),
        le(s + (r + 1) * (2 * s), r"$K_n+1<s_n$ (\texttt{K\_add\_one\_lt}) in each of the $r+1$ terms"),
        eq((2 * r + 3) * s, r"\texttt{ring}"),
        le(2**(2 * n + 7) * s, r"$2r+3\le2n+7<2^{2n+7}$ (\texttt{Nat.lt\_two\_pow\_self}, \texttt{omega}), using $r\le n$"),
        eq(_s(N, n + 1), r"\texttt{s\_succ}", using=AlternatingFlowScales.lemmas.s_succ),
    )
