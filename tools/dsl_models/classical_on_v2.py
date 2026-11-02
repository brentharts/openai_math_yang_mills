"""ClassicalON/CoarseExponential.lean and ClassicalON/DyadicSums.lean of the OpenAI release in the V2 syntax.

These are two modules behind the Lean proof of family 215 (exponential decay of correlations in the two-dimensional
classical O(n) model, n >= 3), the release's nearest relative of the lattice-gauge-theory constructions behind
Yang-Mills (the repository itself contains no Yang-Mills Lean; see YANGMILLS.md).

Each category is one Lean file: explicit Lean statements, connected items, and LaTeX built from the items themselves.
Where the Lean statement is about a type the DSL cannot express (a sum over the Fintype `DyadicFreq k`, a statement
about every `LatticeGraph`), the Lean text is given as it is and the body is the arithmetic the proof reduces it to,
grouped by shell; the docstrings say which reading is used."""
from leandsl import lean, Nat, Real, require, Rational, wikipedia, chain, le, eq, _sym
import sympy as sp
from oai_tree import family   # the generated class of family 215 (the classical O(n) models) of the release

# one reading of each letter for every item here; connect() refuses items that read a letter differently
K_READING = "the number of dyadic scales: the shells have radius $2^0,\\dots,2^{k-1}$"
I_READING = "the index of a dyadic shell, whose radius is 2 to the power i"
T_READING = "the number of halvings of the connection probability, a natural number"
N_READING = "the scale of the geometric tail: the probability halves for every additional 25N of radius"
D_READING = "the distance between the two sites"


def _idx(name: str):
    """A summation index over the natural numbers (the Lean sums over `Fin k` are written over `range k`)."""
    return _sym(name, Nat)


def _pow2(x):
    return sp.Pow(2, x, evaluate=False)


def _shell_card(i):
    """Exactly how many frequencies the shell of radius 2^i has: the type is `Fin (2^i) × Fin (2*2^i+1)`."""
    return _pow2(i) * (2 * _pow2(i) + 1)


_MAIN = family("215").results["OAI.ClassicalON.main"]   # the family's Comparator challenge, a generated Item


# ---------------------------------------------------------------------------------------------------------------
# CoarseExponential.lean: from "halves every 25N" to an exponential
# ---------------------------------------------------------------------------------------------------------------
@lean.category(
    r'''
    \textbf{ClassicalON/CoarseExponential.lean} of the OpenAI release, written in the leandsl DSL. It is a module behind the
    Lean proof of the family's main theorem %s (the import analysis finds it in the family's closure); the family itself
    is [[F215|family 215]]. The file turns the geometric tail proved earlier (the connection probability is at most
    $(1/2)^t$ once the radius reaches $N(25t+4)$) into an exponential bound $4e^{-m\,d}$ in the distance $d$. The one
    identity it needs is [[ClassicalONCoarseExponential.half_pow_exp]], and the last step of the proof of
    [[ClassicalONCoarseExponential.connection_exponential]] is isolated as the lemma
    [[ClassicalONCoarseExponential.half_pow_le_exp_tail]], whose proof is written out as a chain.
    ''' % _MAIN.latex_ref(page=True),
    title="Classical O(n): geometric tail to exponential",
    refs=[wikipedia("N-vector model"), wikipedia("Exponential decay")],
)
class ClassicalONCoarseExponential(family("215")):
    r'''The file declares the identity \texttt{half\_pow\_exp} and the theorem \texttt{connection\_exponential}, whose proof
    ends in a \texttt{calc} that uses it. The theorem quantifies over every lattice graph, which the DSL cannot state, so it
    is recorded as a result with its Lean text; its closing inequality is stated and checked as an item of its own.'''


@ClassicalONCoarseExponential.theorems(
    lean='theorem half_pow_exp (t : ℕ) : (1/2:ℝ)^t=Real.exp (-Real.log 2*(t:ℝ))',
    source="OAI.ClassicalON.half_pow_exp",
    doc="Halving $t$ times is an exponential: $(1/2)^t = e^{-t\\log 2}$.",
    refs=[wikipedia("Exponential decay")],
    reads={"t": T_READING},
)
def half_pow_exp(t: Nat):
    r'''The identity $(1/2)^t=e^{-t\log 2}$ for every natural number $t$; it lets the geometric tail be compared with an exponential.'''
    return sp.Eq(sp.Pow(Rational(1, 2), t, evaluate=False), sp.exp(-sp.log(2) * t, evaluate=False), evaluate=False)


ClassicalONCoarseExponential.results.declare(
    "connection_exponential",
    lean='''theorem connection_exponential (β : ℝ) (hβ : 0 ≤ β) :
    ∃ m : ℝ,0 < m  ∧  ∀ (G : LatticeGraph) (b : G.edges → ℝ),(∀ e,0 ≤ b e  ∧  b e ≤ β) →
      ∀ x y : G.vertices,G.connectionProbability b x y ≤ 4*Real.exp (-m*siteDistance x.val y.val)''',
    doc=r"For every interaction bound $\beta\ge0$ there is $m>0$ such that the connection probability of two sites is at most "
        r"$4e^{-m\,d}$, $d$ the distance between them, on every lattice graph with edge strengths in $[0,\beta]$.",
    source="OAI.ClassicalON.LatticeGraph.connection_exponential",
)


@ClassicalONCoarseExponential.theorems(
    source="",
    doc="The closing step of the proof of \\texttt{connection\\_exponential}: $(1/2)^t \\le 4e^{-m d}$ with $m=\\log 2/(50N)$ "
        "when $d/(50N) < t+2$.",
    reads={"t": T_READING, "N": N_READING, "d": D_READING},
)
def half_pow_le_exp_tail(t: Nat, N: Nat, d: Real):
    r'''This statement is \emph{not} a declaration of the Lean file: it is the arithmetic that the last \texttt{calc} of
    \texttt{connection\_exponential} carries out, with $f=d/(50N)$ and the hypothesis $f<t+2$ that the proof derives from
    $t=\lfloor f\rfloor-1$. The chain below follows that \texttt{calc} step for step.'''
    require(N > 0, d >= 0, d / (50 * N) < t + 2)
    return sp.Pow(Rational(1, 2), t, evaluate=False) <= 4 * sp.exp(-(sp.log(2) / (50 * N)) * d)


# connecting concepts: later used to generate LaTeX that explains how the pieces fit
ClassicalONCoarseExponential.theorems.half_pow_le_exp_tail.connect(ClassicalONCoarseExponential.theorems.half_pow_exp)
ClassicalONCoarseExponential.results.connection_exponential.connect(
    ClassicalONCoarseExponential.theorems.half_pow_exp,
    note="the closing calc of its proof rewrites the geometric tail with it")


# The calc at the end of `connection_exponential` in lean/OAI/Probability/ClassicalON/CoarseExponential.lean.
@ClassicalONCoarseExponential.theorems.half_pow_le_exp_tail.chain(title="The proof is this chain, the final calc of the Lean proof:")
def half_pow_le_exp_tail_chain(t, N, d):
    f = d / (50 * N)
    return chain(
        sp.Pow(Rational(1, 2), t, evaluate=False),
        eq(sp.exp(-sp.log(2) * t), r"$(1/2)^t=e^{-t\log 2}$",
           using=ClassicalONCoarseExponential.theorems.half_pow_exp),
        le(sp.exp(sp.log(2) * (2 - f)),
           r"\texttt{Real.exp\_le\_exp}: $-t\le 2-f$ because $f<t+2$, and $\log 2>0$"),
        eq(4 * sp.exp(-(sp.log(2) / (50 * N)) * d),
           r"$e^{2\log 2}=4$ and $-m\,d=-f\log 2$ for $m=\log 2/(50N)$ (\texttt{Real.exp\_add}, \texttt{ring})"),
    )


# ---------------------------------------------------------------------------------------------------------------
# DyadicSums.lean: the sums over dyadic frequencies that bound the mode costs
# ---------------------------------------------------------------------------------------------------------------
@lean.category(
    r'''
    \textbf{ClassicalON/DyadicSums.lean} of the OpenAI release, written in the leandsl DSL. It is a module behind the Lean
    proof of the family's main theorem %s (the import analysis finds it in the family's closure); the family itself is
    [[F215|family 215]]. The Fourier modes of the lattice are grouped into dyadic shells: shell $i$ has radius $2^i$ and
    contains $2^i(2\cdot2^i+1)$ frequencies, between $2(2^i)^2$ and $3(2^i)^2$ of them. The file sums the costs of all
    modes and shows that they grow only linearly in the number $k$ of shells. The statements about $\texttt{DyadicFreq}\ k$
    are read shell by shell, which is how the Lean proofs treat them. The two chains
    [[ClassicalONDyadicSums.dyadic_inverse_square_sum]] and [[ClassicalONDyadicSums.dyadic_min_cost_sum]] are the
    \texttt{calc} blocks of those proofs.
    ''' % _MAIN.latex_ref(page=True),
    title="Classical O(n): dyadic sums",
    refs=[wikipedia("N-vector model"), wikipedia("Geometric series")],
)
class ClassicalONDyadicSums(family("215")):
    r'''The file proves four bounds, all linear in $k$: a geometric sum $\sum_j[j\le i]\,2^j/2^i\le2$, its double sum
    $\sum_{i,j}\min(2^i/2^j,2^j/2^i)\le4k$, and from these the two sums over all dyadic frequencies,
    $\sum_p r_p^{-2}\le3k$ and $\sum_{p,q}\min(r_p^{-3}r_q^{-1},r_q^{-3}r_p^{-1})\le36k$, where $r_p=2^i$ is the radius of $p$.'''


ClassicalONDyadicSums.results.declare(
    "dyadic_shell_card_bounds",
    lean='''theorem dyadic_shell_card_bounds (j : ℕ) :
    2*(2^j)^2 ≤ Fintype.card (Fin (2^j) × Fin (2*2^j+1)) ∧
      Fintype.card (Fin (2^j) × Fin (2*2^j+1)) ≤ 3*(2^j)^2''',
    doc=r"The shell of radius $2^j$ has between $2(2^j)^2$ and $3(2^j)^2$ frequencies (it has exactly $2^j(2\cdot2^j+1)$).",
    source="OAI.ClassicalON.dyadic_shell_card_bounds",
)


@ClassicalONDyadicSums.theorems(
    lean='''theorem sum_dyadic_ratio_lower (k i : ℕ) :
    (∑ j ∈ Finset.range k, if j ≤ i then (2:ℝ)^j/2^i else 0) ≤ 2''',
    source="OAI.ClassicalON.sum_dyadic_ratio_lower",
    doc="Only the shells at or below $i$ count, and they sum to less than 2: $\\sum_{j<k}[j\\le i]\\,2^j/2^i\\le 2$.",
    refs=[wikipedia("Geometric series")],
    reads={"k": K_READING, "i": I_READING},
)
def sum_dyadic_ratio_lower(k: Nat, i: Nat):
    r'''A geometric series: $\sum_{j=0}^{i}2^j=2^{i+1}-1<2\cdot2^i$.'''
    j = _idx("j")
    return sp.Sum(sp.Piecewise((_pow2(j) / _pow2(i), j <= i), (0, True)), (j, 0, k - 1)) <= 2


@ClassicalONDyadicSums.theorems(
    lean='''theorem weighted_dyadic_min_pair (a b : ℝ) (ha : 0<a) (hb : 0<b) :
    a^2*b^2*min ((1/a)^3*(1/b)) ((1/b)^3*(1/a)) = min (b/a) (a/b)''',
    source="OAI.ClassicalON.weighted_dyadic_min_pair",
    doc="Weighting the pair of costs by $a^2b^2$ turns it into the ratio $\\min(b/a,a/b)$.",
    reads={"a": "a radius, in the application a power of two $2^i$", "b": "a radius, in the application a power of two $2^j$"},
)
def weighted_dyadic_min_pair(a: Real, b: Real):
    r'''Two radii $a,b>0$; applied with $a=2^i$, $b=2^j$, where $a^2b^2$ is the size of the weight $4^{i+j}$ of two shells.'''
    require(a > 0, b > 0)
    return sp.Eq(a**2 * b**2 * sp.Min((1 / a)**3 * (1 / b), (1 / b)**3 * (1 / a)), sp.Min(b / a, a / b), evaluate=False)


@ClassicalONDyadicSums.theorems(
    lean='''theorem dyadic_min_ratio_sum (k : ℕ) :
    (∑ i : Fin k, ∑ j : Fin k, min ((2:ℝ)^i.val/2^j.val) (2^j.val/2^i.val)) ≤ 4*(k:ℝ)''',
    source="OAI.ClassicalON.dyadic_min_ratio_sum",
    doc="The ratios of two shells decay away from the diagonal, so their sum is linear: "
        "$\\sum_{i,j<k}\\min(2^i/2^j,2^j/2^i)\\le 4k$.",
    reads={"k": K_READING},
)
def dyadic_min_ratio_sum(k: Nat):
    r'''The minimum is at most the part below the diagonal plus the part above it; by symmetry the two double sums agree, and each
    row of one of them is at most 2.'''
    i, j = _idx("i"), _idx("j")
    return sp.Sum(sp.Sum(sp.Min(_pow2(i) / _pow2(j), _pow2(j) / _pow2(i)), (j, 0, k - 1)), (i, 0, k - 1)) <= 4 * k


@ClassicalONDyadicSums.theorems(
    lean='''theorem dyadic_inverse_square_sum (k : ℕ) :
    (∑ p : DyadicFreq k, (1/(dyadicRadius p:ℝ))^2) ≤ 3*(k:ℝ)''',
    source="OAI.ClassicalON.dyadic_inverse_square_sum",
    doc="Summed over all dyadic frequencies, the inverse squared radii are at most $3k$.",
    reads={"k": K_READING},
)
def dyadic_inverse_square_sum(k: Nat):
    r'''Read shell by shell: the shell $j$ has $2^j(2\cdot2^j+1)$ frequencies, each of radius $2^j$, so it contributes
    $2^j(2\cdot2^j+1)\,(1/2^j)^2$, which is at most 3.'''
    j = _idx("j")
    return sp.Sum(_shell_card(j) * (1 / _pow2(j))**2, (j, 0, k - 1)) <= 3 * k


@ClassicalONDyadicSums.theorems(
    lean='''theorem dyadic_min_cost_sum (k : ℕ) :
    (∑ p : DyadicFreq k, ∑ q : DyadicFreq k,
      min ((1/(dyadicRadius p:ℝ))^3*(1/(dyadicRadius q:ℝ)))
          ((1/(dyadicRadius q:ℝ))^3*(1/(dyadicRadius p:ℝ)))) ≤ 36*(k:ℝ)''',
    source="OAI.ClassicalON.dyadic_min_cost_sum",
    doc="The pairwise mode costs summed over all pairs of dyadic frequencies are at most $36k$.",
    reads={"k": K_READING},
)
def dyadic_min_cost_sum(k: Nat):
    r'''Read shell by shell, as in the Lean proof (its step \texttt{he}): the pair of shells $(i,j)$ has $2^i(2\cdot2^i+1)$ times
    $2^j(2\cdot2^j+1)$ pairs of frequencies, all with the same cost.'''
    i, j = _idx("i"), _idx("j")
    cost = sp.Min((1 / _pow2(i))**3 * (1 / _pow2(j)), (1 / _pow2(j))**3 * (1 / _pow2(i)))
    return sp.Sum(sp.Sum(_shell_card(i) * _shell_card(j) * cost, (j, 0, k - 1)), (i, 0, k - 1)) <= 36 * k


# connecting concepts: later used to generate LaTeX that explains how the pieces fit
ClassicalONDyadicSums.theorems.dyadic_min_ratio_sum.connect(ClassicalONDyadicSums.theorems.sum_dyadic_ratio_lower)
ClassicalONDyadicSums.theorems.dyadic_min_cost_sum.connect(
    ClassicalONDyadicSums.theorems.weighted_dyadic_min_pair, ClassicalONDyadicSums.theorems.dyadic_min_ratio_sum)


# The calc of `sum_dyadic_ratio_lower`, step for step.
@ClassicalONDyadicSums.theorems.sum_dyadic_ratio_lower.chain(title="The proof is this chain:")
def sum_dyadic_ratio_lower_chain(k, i):
    j = _idx("j")
    return chain(
        sp.Sum(sp.Piecewise((_pow2(j) / _pow2(i), j <= i), (0, True)), (j, 0, k - 1)),
        le(sp.Sum(_pow2(j) / _pow2(i), (j, 0, i)),
           r"\texttt{Finset.sum\_le\_sum\_of\_subset\_of\_nonneg}: the shells $j\le i$ are among $0,\dots,i$, and every term is $\ge0$"),
        eq((_pow2(i + 1) - 1) / _pow2(i), r"\texttt{geom\_sum\_mul}: $\sum_{j=0}^{i}2^j=2^{i+1}-1$"),
        le(2, r"$(2^{i+1}-1)/2^i\le2$ (\texttt{div\_le\_iff}, \texttt{linarith})"),
    )


# The calc of `weighted_dyadic_min_pair`: the weight is pulled into the minimum, then each side is simplified.
@ClassicalONDyadicSums.theorems.weighted_dyadic_min_pair.chain(title="The proof is this chain:")
def weighted_dyadic_min_pair_chain(a, b):
    return chain(
        a**2 * b**2 * sp.Min((1 / a)**3 * (1 / b), (1 / b)**3 * (1 / a)),
        eq(sp.Min(sp.Mul(a**2 * b**2, (1 / a)**3 * (1 / b), evaluate=False),
                  sp.Mul(a**2 * b**2, (1 / b)**3 * (1 / a), evaluate=False), evaluate=False),
           r"\texttt{mul\_min\_of\_nonneg}, because $a^2b^2\ge0$"),
        eq(sp.Min(b / a, a / b), r"\texttt{field\_simp}, one argument of the minimum at a time"),
    )


# The calc of `dyadic_min_ratio_sum`: bound the minimum by the two triangles, fold them, bound each row.
@ClassicalONDyadicSums.theorems.dyadic_min_ratio_sum.chain(title="The proof is this chain:")
def dyadic_min_ratio_sum_chain(k):
    i, j = _idx("i"), _idx("j")
    below = sp.Piecewise((_pow2(j) / _pow2(i), j <= i), (0, True))
    above = sp.Piecewise((_pow2(i) / _pow2(j), i <= j), (0, True))
    return chain(
        sp.Sum(sp.Sum(sp.Min(_pow2(i) / _pow2(j), _pow2(j) / _pow2(i)), (j, 0, k - 1)), (i, 0, k - 1)),
        le(sp.Sum(sp.Sum(sp.Add(below, above, evaluate=False), (j, 0, k - 1)), (i, 0, k - 1)),
           r"\texttt{point}: the minimum is at most the part with $j\le i$ plus the part with $i\le j$"),
        eq(2 * sp.Sum(sp.Sum(below, (j, 0, k - 1)), (i, 0, k - 1)),
           r"\texttt{Finset.sum\_add\_distrib} and \texttt{Finset.sum\_comm}: the two triangles are mirror images"),
        le(2 * sp.Sum(2, (i, 0, k - 1)),
           r"each row $\sum_j[j\le i]2^j/2^i$ is at most 2",
           using=ClassicalONDyadicSums.theorems.sum_dyadic_ratio_lower),
        eq(4 * k, r"\texttt{simp}, \texttt{ring}: $\sum_{i<k}2=2k$"),
    )


# The calc of `dyadic_inverse_square_sum`, read shell by shell.
@ClassicalONDyadicSums.theorems.dyadic_inverse_square_sum.chain(title="The proof is this chain:")
def dyadic_inverse_square_sum_chain(k):
    j = _idx("j")
    return chain(
        sp.Sum(_shell_card(j) * (1 / _pow2(j))**2, (j, 0, k - 1)),
        le(sp.Sum(3, (j, 0, k - 1)),
           r"\texttt{Finset.sum\_le\_sum}: the shell has at most $3(2^j)^2$ frequencies, and $(2^j)^2(1/2^j)^2=1$",
           using=ClassicalONDyadicSums.results.dyadic_shell_card_bounds),
        eq(3 * k, r"\texttt{simp}, \texttt{ring}: $\sum_{j<k}3=3k$"),
    )


# The calc of `dyadic_min_cost_sum`, read shell by shell.
@ClassicalONDyadicSums.theorems.dyadic_min_cost_sum.chain(title="The proof is this chain:")
def dyadic_min_cost_sum_chain(k):
    i, j = _idx("i"), _idx("j")
    cost = sp.Min((1 / _pow2(i))**3 * (1 / _pow2(j)), (1 / _pow2(j))**3 * (1 / _pow2(i)))
    ratio = sp.Min(_pow2(i) / _pow2(j), _pow2(j) / _pow2(i))
    return chain(
        sp.Sum(sp.Sum(_shell_card(i) * _shell_card(j) * cost, (j, 0, k - 1)), (i, 0, k - 1)),
        le(sp.Sum(sp.Sum(3 * _pow2(i)**2 * (3 * _pow2(j)**2) * cost, (j, 0, k - 1)), (i, 0, k - 1)),
           r"\texttt{gcongr}: the two shells have at most $3(2^i)^2$, resp.\ $3(2^j)^2$, frequencies",
           using=ClassicalONDyadicSums.results.dyadic_shell_card_bounds),
        eq(sp.Sum(sp.Sum(9 * ratio, (j, 0, k - 1)), (i, 0, k - 1)),
           r"with $a=2^i$, $b=2^j$ (the minimum of the two ratios is symmetric)",
           using=ClassicalONDyadicSums.theorems.weighted_dyadic_min_pair),
        le(36 * k,
           r"$9\sum_{i,j}\min(2^i/2^j,2^j/2^i)\le 9\cdot4k$",
           using=ClassicalONDyadicSums.theorems.dyadic_min_ratio_sum),
    )
