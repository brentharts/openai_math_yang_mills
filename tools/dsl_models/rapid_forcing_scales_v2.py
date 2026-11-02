"""MathematicalPhysics/RapidForcing/Scales.lean of the OpenAI release in the V2 syntax: a category with explicit Lean statements,
connected items, and LaTeX prose built from the items themselves.

The file sits in the closure of family 235 (limiting random SAT thresholds, with the computability result), through the
RapidForcing construction: a fluid-like forcing that is made to run a Turing machine. The file holds the bookkeeping of
scales that construction needs: a record of the machine at stage n is stored at the dyadic scale delta_n = 2^(-s_n), the
number of values it may take is 2^(L_n), and the sizes are chosen so that the total of all stored information is below 1.

The Lean statements are about a structure `d : ScaleData` with natural-number fields `rA`, `rQ`, `B` and derived functions
L, C, s, capacity, D, delta, b. The DSL has no structures, so each item here is stated for the numbers themselves
(C, L, s are letters, or are written out in `rA`, `rQ`, `B` where the statement needs the definition) and the Lean text is
given with `lean=` as it is in the file (a `variable (d : ScaleData)` stands over the section, so `d.s n` is the stage-n
scale). The tsum and the address-index lemmas, which use infinite sums and `Finset` sums over a function k, are
recorded as results with their Lean text and are not given a body.

Not modelled: the rest of the file's 300 lines of analysis (summability of the geometric majorant, reading the newest record
from an address) beyond the results below, and the other forty modules of RapidForcing."""
from leandsl import lean, Nat, Real, require, Rational, wikipedia, chain, le, lt, eq
import sympy as sp
from oai_tree import family   # the generated class of family 235 of the release

# one reading of each letter for every item here; connect() refuses items that read a letter differently
N_READING = "the stage of the construction, a natural number"
C_READING = "the scale constant of the construction: $C=\\max(1,\\max(L_0,1+2r_A))$"
S_READING = "the exponent of the stage-$n$ scale: $s_n=C\\cdot2^{(n+2)^2}$"
L_READING = "the exponent of the stage-$n$ capacity: $L_n=r_Q+n+2+r_A(2n+B+1)$, so the capacity is $2^{L_n}$"
RA_READING = "the binary logarithm of the number of tape symbols, a natural number"
RQ_READING = "the binary logarithm of the number of states, a natural number"
B_READING = "the length of the input word, a natural number"

_MAIN = family("235").results["OAI.FixedClauseThreshold.Computability.main"]   # the family's Comparator challenge, a generated Item


def _L(n, rA, rQ, B):
    """The capacity exponent of stage n."""
    return rQ + n + 2 + rA * (2 * n + B + 1)


def _C(rA, rQ, B):
    return sp.Max(1, sp.Max(_L(0, rA, rQ, B), 1 + 2 * rA))


@lean.category(
    r'''
    \textbf{MathematicalPhysics/RapidForcing/Scales.lean} of the OpenAI release, written in the leandsl DSL. It is a module in
    the closure of the Lean proof of the family's main theorem %s, through the RapidForcing construction (a forcing term that
    makes a fluid-like flow run a Turing machine); the family itself is [[F235|family 235]]. The construction stores the
    record of stage $n$ at the scale $\delta_n=2^{-s_n}$ with $s_n=C\cdot2^{(n+2)^2}$, and the record can take $2^{L_n}$
    values with $L_n$ linear in $n$. Three kinds of fact are needed. \emph{Growth:} $s_n$ grows by the factor $2^{2n+5}$ at each
    stage ([[RapidForcingScales.s_succ]]), faster than any capacity, which is
    [[RapidForcingScales.memory_exponents]] and [[RapidForcingScales.small_motion_exponents]] and follows from the
    lower bound [[RapidForcingScales.s_lower]] on the quadratic exponent
    ([[RapidForcingScales.power_square_lower]]). \emph{Memory:} the product $\delta_n\cdot2^{L_n}$ is
    $2^{-(s_n-L_n)}$ ([[RapidForcingScales.delta_capacity]]) and at most $2^{-15(n+1)}$, so the total over all stages is below
    $1/32767$ ([[RapidForcingScales.total_memory_le]]). \emph{Addresses:} that bound makes the history of the first $n$
    stages a number below 1, which is an integer multiple of $\delta_n$ ([[RapidForcingScales.delta_ratio]],
    [[RapidForcingScales.read_newest_record]]), so the newest record can be read off by a remainder.
    ''' % _MAIN.latex_ref(page=True),
    title="Rapid forcing scales",
    refs=[wikipedia("Geometric series"), wikipedia("Turing machine")],
)
class RapidForcingScales(family("235")):
    r'''The scale data has three natural numbers $r_A,r_Q,B$ and defines $L_n=r_Q+n+2+r_A(2n+B+1)$,
    $C=\max(1,\max(L_0,1+2r_A))$, $s_n=C\,2^{(n+2)^2}$, $\delta_n=2^{-s_n}$ and the capacity $2^{L_n}$.'''


@RapidForcingScales.lemmas(
    lean='''lemma power_square_lower (n : ℕ) : 16 * (n + 1) ≤ 2 ^ ((n + 2) ^ 2)''',
    source="OAI.RapidForcing.ScaleData.power_square_lower",
    doc="The quadratic exponent is large: $16(n+1)\\le2^{(n+2)^2}$.",
    reads={"n": N_READING},
)
def power_square_lower(n: Nat):
    r'''Because $n+1\le2^n$ and $n+4\le(n+2)^2$, so $16(n+1)\le16\cdot2^n=2^{n+4}\le2^{(n+2)^2}$.'''
    return 16 * (n + 1) <= 2**((n + 2)**2)


@RapidForcingScales.lemmas(
    lean='''lemma s_succ (n : ℕ) : d.s (n + 1) = 2 ^ (2 * n + 5) * d.s n''',
    source="OAI.RapidForcing.ScaleData.s_succ",
    doc="One stage later the scale exponent is $2^{2n+5}$ times larger: $s_{n+1}=2^{2n+5}s_n$.",
    reads={"n": N_READING, "C": C_READING},
)
def s_succ(C: Nat, n: Nat):
    r'''Because $(n+3)^2=(2n+5)+(n+2)^2$.'''
    return sp.Eq(C * 2**((n + 3)**2), 2**(2 * n + 5) * (C * 2**((n + 2)**2)), evaluate=False)


@RapidForcingScales.lemmas(
    lean='''lemma s_succ_lower (n : ℕ) : 32 * (n + 1) * d.s n ≤ d.s (n + 1)''',
    source="OAI.RapidForcing.ScaleData.s_succ_lower",
    doc="Each stage multiplies the scale exponent by at least $32(n+1)$: $32(n+1)s_n\\le s_{n+1}$.",
    reads={"n": N_READING, "C": C_READING},
)
def s_succ_lower(C: Nat, n: Nat):
    r'''By \texttt{s\_succ} the factor is $2^{2n+5}=32\cdot2^{2n}$, and $n+1\le2^n\le2^{2n}$.'''
    return 32 * (n + 1) * (C * 2**((n + 2)**2)) <= C * 2**((n + 3)**2)


@RapidForcingScales.lemmas(
    lean='''lemma s_lower (n : ℕ) : 16 * d.C * (n + 1) ≤ d.s n''',
    source="OAI.RapidForcing.ScaleData.s_lower",
    doc="The scale exponent is at least $16C(n+1)$.",
    reads={"n": N_READING, "C": C_READING},
)
def s_lower(C: Nat, n: Nat):
    r'''Multiply \texttt{power\_square\_lower} by $C$.'''
    return 16 * C * (n + 1) <= C * 2**((n + 2)**2)


@RapidForcingScales.lemmas(
    lean='''lemma L_le (n : ℕ) : d.L n ≤ d.C * (n + 1)''',
    source="OAI.RapidForcing.ScaleData.L_le",
    doc="The capacity exponent is at most $C(n+1)$.",
    reads={"n": N_READING, "rA": RA_READING, "rQ": RQ_READING, "B": B_READING},
)
def L_le(n: Nat, rA: Nat, rQ: Nat, B: Nat):
    r'''$L_n=L_0+(1+2r_A)n$ is linear in $n$ with the slope $1+2r_A\le C$ and the intercept $L_0\le C$.'''
    return _L(n, rA, rQ, B) <= _C(rA, rQ, B) * (n + 1)


@RapidForcingScales.lemmas(
    lean='''lemma memory_exponents (n : ℕ) : 15 * (n + 1) + d.L n ≤ d.s n''',
    source="OAI.RapidForcing.ScaleData.memory_exponents",
    doc="The scale exponent exceeds the capacity exponent by at least $15(n+1)$.",
    reads={"n": N_READING, "rA": RA_READING, "rQ": RQ_READING, "B": B_READING},
)
def memory_exponents(n: Nat, rA: Nat, rQ: Nat, B: Nat):
    r'''By \texttt{L\_le} and $C\ge1$ the left side is at most $15(n+1)+C(n+1)\le16C(n+1)$, which is at most $s_n$ by \texttt{s\_lower}.'''
    return 15 * (n + 1) + _L(n, rA, rQ, B) <= _C(rA, rQ, B) * 2**((n + 2)**2)


@RapidForcingScales.lemmas(
    lean='''lemma small_motion_exponents (n : ℕ) :
    (n + 10) * d.s n + d.L (n + 1) ≤ d.s (n + 1)''',
    source="OAI.RapidForcing.ScaleData.small_motion_exponents",
    doc="The next scale exponent exceeds $(n+10)s_n$ by at least the next capacity exponent: $(n+10)s_n+L_{n+1}\\le s_{n+1}$.",
    reads={"n": N_READING, "rA": RA_READING, "rQ": RQ_READING, "B": B_READING},
)
def small_motion_exponents(n: Nat, rA: Nat, rQ: Nat, B: Nat):
    r'''Since $L_{n+1}\le C(n+2)\le(n+2)s_n$, the left side is at most $(2n+12)s_n\le32(n+1)s_n\le s_{n+1}$.'''
    C = _C(rA, rQ, B)
    return (n + 10) * (C * 2**((n + 2)**2)) + _L(n + 1, rA, rQ, B) <= C * 2**((n + 3)**2)


@RapidForcingScales.lemmas(
    lean='''lemma delta_capacity (n : ℕ) :
    d.δ n * (d.capacity n : ℝ) = (1 / 2 : ℝ) ^ (d.s n - d.L n)''',
    source="OAI.RapidForcing.ScaleData.delta_capacity",
    doc="The information stored at stage $n$ is $\\delta_n\\cdot2^{L_n}=2^{-(s_n-L_n)}$.",
    reads={"S": S_READING, "L": L_READING},
)
def delta_capacity(S: Nat, L: Nat):
    r'''Because $L_n\le s_n$ (\texttt{L\_le\_s}), the exponent $s_n-L_n$ is a natural number and $2^{-s_n}\cdot2^{L_n}=2^{-(s_n-L_n)}$.'''
    require(L <= S)
    return sp.Eq(sp.Pow(Rational(1, 2), S, evaluate=False) * 2**L, sp.Pow(Rational(1, 2), S - L, evaluate=False), evaluate=False)


@RapidForcingScales.lemmas(
    lean='''lemma delta_ratio {j n : ℕ} (hjn : j ≤ n) :
    d.δ j / d.δ n = ((2 ^ (d.s n - d.s j) : ℕ) : ℝ)''',
    source="OAI.RapidForcing.ScaleData.delta_ratio",
    doc="The ratio of two scales is a whole number: $\\delta_j/\\delta_n=2^{s_n-s_j}$ for $j\\le n$.",
    reads={"Sj": S_READING, "Sn": S_READING},
)
def delta_ratio(Sj: Nat, Sn: Nat):
    r'''So a multiple of $\delta_j$ is a multiple of $\delta_n$; this is what lets the history be an integer multiple of the newest scale.'''
    require(Sj <= Sn)
    return sp.Eq(sp.Pow(Rational(1, 2), Sj, evaluate=False) / sp.Pow(Rational(1, 2), Sn, evaluate=False), 2**(Sn - Sj), evaluate=False)


@RapidForcingScales.lemmas(
    lean='''lemma memory_majorant_eq (n : ℕ) :
    (1 / 2 : ℝ) ^ (15 * (n + 1)) = ((1 / 2 : ℝ) ^ 15) ^ n * (1 / 2 : ℝ) ^ 15''',
    source="OAI.RapidForcing.memory_majorant_eq",
    doc="The memory majorant is geometric: $2^{-15(n+1)}=(2^{-15})^n\\cdot2^{-15}$.",
    reads={"n": N_READING},
)
def memory_majorant_eq(n: Nat):
    r'''Splitting the exponent $15(n+1)=15n+15$.'''
    return sp.Eq(sp.Pow(Rational(1, 2), 15 * (n + 1), evaluate=False),
                 sp.Pow(Rational(1, 2**15), n, evaluate=False) * Rational(1, 2**15), evaluate=False)


RapidForcingScales.results.declare(
    "memory_majorant_sum",
    lean='''lemma memory_majorant_sum :
    (∑' n : ℕ, (1 / 2 : ℝ) ^ (15 * (n + 1))) = 1 / 32767''',
    doc=r"The majorant sums to $\sum_{n\ge0}2^{-15(n+1)}=2^{-15}/(1-2^{-15})=1/32767$, a geometric series.",
    source="OAI.RapidForcing.memory_majorant_sum",
    refs=[wikipedia("Geometric series")],
)

RapidForcingScales.results.declare(
    "total_memory_le",
    lean='''lemma total_memory_le :
    (∑' n : ℕ, d.δ n * (d.capacity n : ℝ)) ≤ 1 / 32767''',
    doc=r"The information stored over all stages is at most $1/32767$: each term is at most $2^{-15(n+1)}$ "
        r"(\texttt{delta\_capacity\_le}, from \texttt{delta\_capacity} and \texttt{memory\_exponents}) and the majorant sums to $1/32767$.",
    source="OAI.RapidForcing.ScaleData.total_memory_le",
)

RapidForcingScales.results.declare(
    "read_newest_record",
    lean='''theorem read_newest_record (k : ℕ → ℕ) (n : ℕ)
    (hk : ∀ j ≤ n, k j < d.capacity j) :
    ∃ m : ℕ, m < d.D n ∧ 0 ≤ d.history k n ∧ d.history k n < 1 ∧
      d.history k n = (m : ℝ) * d.δ n ∧ m % d.capacity n = k n''',
    doc=r"If each record $k_j$ is below the capacity of its stage, the history $\sum_{j\le n}\delta_jk_j$ lies in $[0,1)$, is an integer "
        r"multiple $m\delta_n$ of the newest scale with $m<2^{s_n}$, and the newest record is recovered as $m\bmod2^{L_n}=k_n$.",
    source="OAI.RapidForcing.ScaleData.read_newest_record",
)


# connecting concepts: later used to generate LaTeX that explains how the pieces fit
RapidForcingScales.lemmas.s_lower.connect(RapidForcingScales.lemmas.power_square_lower)
RapidForcingScales.lemmas.s_succ_lower.connect(RapidForcingScales.lemmas.s_succ)
RapidForcingScales.lemmas.memory_exponents.connect(RapidForcingScales.lemmas.L_le, RapidForcingScales.lemmas.s_lower)
RapidForcingScales.lemmas.small_motion_exponents.connect(RapidForcingScales.lemmas.L_le, RapidForcingScales.lemmas.s_succ_lower)
RapidForcingScales.lemmas.memory_majorant_eq.connect(RapidForcingScales.results.memory_majorant_sum)
RapidForcingScales.results.total_memory_le.connect(
    RapidForcingScales.lemmas.delta_capacity, RapidForcingScales.lemmas.memory_exponents,
    RapidForcingScales.results.memory_majorant_sum)
RapidForcingScales.results.read_newest_record.connect(
    RapidForcingScales.results.total_memory_le, RapidForcingScales.lemmas.delta_ratio)


# The proofs that are one calc each, step for step.
@RapidForcingScales.lemmas.power_square_lower.chain(title="The proof is this chain:")
def power_square_lower_chain(n):
    return chain(
        16 * (n + 1),
        le(16 * 2**n, r"$n+1\le2^n$ (\texttt{Nat.lt\_two\_pow\_self})"),
        eq(2**(n + 4), r"\texttt{pow\_add}"),
        le(2**((n + 2)**2), r"$n+4\le(n+2)^2$ (\texttt{nlinarith}) and \texttt{Nat.pow\_le\_pow\_right}"),
    )


@RapidForcingScales.lemmas.s_lower.chain(title="The proof is this chain:")
def s_lower_chain(C, n):
    return chain(
        16 * C * (n + 1),
        eq(C * (16 * (n + 1)), r"\texttt{ring}"),
        le(C * 2**((n + 2)**2), r"multiply \texttt{power\_square\_lower} by $C$",
           using=RapidForcingScales.lemmas.power_square_lower),
    )


@RapidForcingScales.lemmas.L_le.chain(title="The proof is this chain:")
def L_le_chain(n, rA, rQ, B):
    return chain(
        _L(n, rA, rQ, B),
        eq(_L(0, rA, rQ, B) + (1 + 2 * rA) * n, r"\texttt{L\_eq}: $L_n=L_0+(1+2r_A)n$ (\texttt{ring})"),
        le(_C(rA, rQ, B) + _C(rA, rQ, B) * n,
           r"$L_0\le C$ (\texttt{L\_zero\_le\_C}) and $1+2r_A\le C$ (\texttt{slope\_le\_C}), both from the definition of $C$ as a maximum"),
        eq(_C(rA, rQ, B) * (n + 1), r"\texttt{ring}"),
    )


@RapidForcingScales.lemmas.memory_exponents.chain(title="The proof is this chain:")
def memory_exponents_chain(n, rA, rQ, B):
    C = _C(rA, rQ, B)
    return chain(
        15 * (n + 1) + _L(n, rA, rQ, B),
        le(15 * (n + 1) + C * (n + 1), r"$L_n\le C(n+1)$", using=RapidForcingScales.lemmas.L_le),
        le(16 * C * (n + 1), r"$C\ge1$ (\texttt{one\_le\_C}): $15(n+1)\le15C(n+1)$"),
        le(C * 2**((n + 2)**2), r"$16C(n+1)\le s_n$", using=RapidForcingScales.lemmas.s_lower),
    )


@RapidForcingScales.lemmas.small_motion_exponents.chain(title="The proof is this chain:")
def small_motion_exponents_chain(n, rA, rQ, B):
    C = _C(rA, rQ, B)
    s = C * 2**((n + 2)**2)
    return chain(
        (n + 10) * s + _L(n + 1, rA, rQ, B),
        le((n + 10) * s + (n + 2) * s,
           r"$L_{n+1}\le C(n+2)$ (\texttt{L\_le}) and $C\le s_n$ (\texttt{C\_le\_s})", using=RapidForcingScales.lemmas.L_le),
        eq((2 * n + 12) * s, r"\texttt{ring}"),
        le(32 * (n + 1) * s, r"$2n+12\le32(n+1)$ (\texttt{nlinarith})"),
        le(C * 2**((n + 3)**2), r"\texttt{s\_succ\_lower}", using=RapidForcingScales.lemmas.s_succ_lower),
    )
