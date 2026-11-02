"""PeriodicTiling/ResidueCountBound.lean of the OpenAI release in the V2 syntax: a category with explicit Lean statements,
connected items, and LaTeX prose built from the items themselves.

The file belongs to family 155 (a translational tile with no fully periodic tiling in dimension three). It holds the
numerical estimate behind a probabilistic-method step of that proof: a union bound over at most 129^3 events, each of which
fails with probability at most (15598/15625)^32000.

The Lean statements are over the rationals; the DSL has the types Nat, Int and Real only, so the body of an item (from which
the LaTeX is made) reads its rationals as reals, and the Lean text is given with `lean=` as it is in the file. Where the DSL
cannot state a declaration (a natural-number division) it is recorded as a result with its Lean text."""
from leandsl import lean, Nat, Real, require, Rational, wikipedia, chain, le, lt, eq
import sympy as sp
from oai_tree import family   # the generated class of family 155 of the release

# one reading of each letter for every item here; connect() refuses items that read a letter differently
C_READING = "the probability that one pair of colours is good: in the application $27/125^2=27/15625$"
R_READING = "the number of independent pairs, an exponent"
P_READING = "an upper bound for the probability that one pair is bad: in the application $15598/15625=1-27/15625$"
X_READING = "the probability that a block of 1000 pairs are all bad, which is $(1-c)^{1000}$ in the application"
K_READING = "the number of independent pairs of one event, an exponent"

_P = Rational(15598, 15625)          # 1 - 27/15625, the probability that one pair is bad
_P32000 = lean.raw("((15598 : ℝ) / 15625) ^ 32000", tex=r"\left(\frac{15598}{15625}\right)^{32000}")   # see residue_union_bound_32000


def _pow(base, exp):
    """A power that sympy does not evaluate (15598/15625 to the 32000th power has 134000 digits)."""
    return sp.Pow(base, exp, evaluate=False)


_MAIN = family("155").results["OAI.PeriodicTilingThree.periodic_tiling_counterexample_and_minimality"]   # a generated Item


@lean.category(
    r'''
    \textbf{Geometry/PeriodicTiling/ResidueCountBound.lean} of the OpenAI release, written in the leandsl DSL. It is a module
    behind the Lean proof of the family's main theorem %s (the import analysis finds it in the family's closure); the family
    itself is [[F155|family 155]]. The proof needs a colouring of the residues for which every one of at most $129^3$ requested
    differences is realised by one of $32000$ independent pairs. With $125$ colours, a pair is good for at least $27$ of the
    $125^2=15625$ pairs of colours, so a uniformly random colouring leaves one pair bad with probability at most
    $15598/15625=1-c$ for $c=27/15625$, and all $32000$ pairs of an event bad with probability at most
    $(15598/15625)^{32000}$. The file proves that the union bound over the events is below 1
    ([[PeriodicTilingResidueBounds.residue_union_bound_32000]]), so a good colouring exists
    (the module \texttt{ResidueTransversal.lean} applies it). The estimate is built from Bernoulli's inequality
    [[PeriodicTilingResidueBounds.one_sub_pow_mul_linear_le_one]] and a halving step, in four stages that
    are connected below.
    ''' % _MAIN.latex_ref(page=True),
    title="Residue count bound",
    refs=[wikipedia("Probabilistic method"), wikipedia("Boole's inequality")],
)
class PeriodicTilingResidueBounds(family("155")):
    r'''The file declares nine statements, from Bernoulli's inequality $(1-c)^r(1+rc)\le1$ down to the numbers $64^3=262144$ and
    $32000\le64^3/3-1$. All Lean statements here are over $\mathbb Q$, the exponent of the bad-pair probability is a natural number.'''


@PeriodicTilingResidueBounds.theorems(
    lean='''theorem one_sub_pow_mul_linear_le_one (c : ℚ) (_hc0 : 0 ≤ c) (hc1 : c ≤ 1)
    (r : ℕ) : (1 - c) ^ r * (1 + (r : ℚ) * c) ≤ 1''',
    source="OAI.PeriodicTilingThree.one_sub_pow_mul_linear_le_one",
    doc="Bernoulli's inequality in the form $(1-c)^r(1+rc)\\le 1$ for $0\\le c\\le1$.",
    refs=[wikipedia("Bernoulli's inequality")],
    reads={"c": C_READING, "r": R_READING},
)
def one_sub_pow_mul_linear_le_one(c: Real, r: Nat):
    r'''The probability that $r$ independent pairs are all bad, $(1-c)^r$, times $1+rc$ is at most 1. The Lean proof is an induction on $r$;
    its inductive step is the item \texttt{one\_sub\_pow\_induction\_step} below.'''
    require(c >= 0, c <= 1)
    return (1 - c)**r * (1 + r * c) <= 1


@PeriodicTilingResidueBounds.theorems(
    source="",
    doc="The inductive step of Bernoulli's inequality: if $(1-c)^r(1+rc)\\le1$ then $(1-c)^{r+1}(1+(r+1)c)\\le1$.",
    reads={"c": C_READING, "r": R_READING},
)
def one_sub_pow_induction_step(c: Real, r: Nat):
    r'''This statement is \emph{not} a declaration of the Lean file: it is the \texttt{succ} case of the induction in the proof of
    \texttt{one\_sub\_pow\_mul\_linear\_le\_one}, with the induction hypothesis as an assumption. The chain below is the
    \texttt{calc} of that case, step for step.'''
    require(c >= 0, c <= 1, (1 - c)**r * (1 + r * c) <= 1)
    return (1 - c)**(r + 1) * (1 + (r + 1) * c) <= 1


@PeriodicTilingResidueBounds.theorems(
    lean='''private theorem scalar_half_bound (x c : ℚ) (hx : 0 ≤ x)
    (hc : (1 : ℚ) / 1000 ≤ c) (h : x * (1 + 1000 * c) ≤ 1) : x ≤ 1 / 2''',
    source="OAI.PeriodicTilingThree.scalar_half_bound",
    doc="If $x(1+1000c)\\le1$ and $c\\ge1/1000$ then $x\\le1/2$.",
    reads={"x": X_READING, "c": C_READING},
)
def scalar_half_bound(x: Real, c: Real):
    r'''Because $1+1000c\ge2$ when $c\ge1/1000$: the factor $1+1000c$ at least doubles $x\ge0$ and the product is still at most 1.'''
    require(x >= 0, c >= Rational(1, 1000), x * (1 + 1000 * c) <= 1)
    return x <= Rational(1, 2)


@PeriodicTilingResidueBounds.theorems(
    lean='''theorem residue_block_bound (c : ℚ) (hc0 : 0 ≤ c) (hc1 : c ≤ 1)
    (hc : (1 : ℚ) / 1000 ≤ c) : (1 - c) ^ 1000 ≤ 1 / 2''',
    source="OAI.PeriodicTilingThree.residue_block_bound",
    doc="A block of 1000 independent pairs is all bad with probability at most $1/2$ when each is good with probability $c\\ge1/1000$.",
    reads={"c": C_READING},
)
def residue_block_bound(c: Real):
    r'''The Lean proof applies \texttt{scalar\_half\_bound} to $x=(1-c)^{1000}$; its hypothesis $x(1+1000c)\le1$ is Bernoulli's
    inequality with $r=1000$.'''
    require(c >= 0, c <= 1, c >= Rational(1, 1000))
    return (1 - c)**1000 <= Rational(1, 2)


@PeriodicTilingResidueBounds.theorems(
    lean='''theorem residue_failure_block_bound : ((15598 : ℚ) / 15625) ^ 1000 ≤ 1 / 2''',
    source="OAI.PeriodicTilingThree.residue_failure_block_bound",
    doc="For the application, $(15598/15625)^{1000}\\le1/2$.",
)
def residue_failure_block_bound():
    r'''The block bound at $c=27/15625$: a block of $1000$ pairs is all bad with probability at most $1/2$ (exactly, it is about $0.177$).'''
    return sp.Le(_pow(_P, 1000), Rational(1, 2), evaluate=False)


@PeriodicTilingResidueBounds.theorems(
    lean='''private theorem union_bound_from_block (p : ℚ) (hp : 0 ≤ p)
    (hblock : p ^ 1000 ≤ 1 / 2) : (129 : ℚ) ^ 3 * p ^ 32000 < 1''',
    source="OAI.PeriodicTilingThree.union_bound_from_block",
    doc="If a block of 1000 pairs is all bad with probability at most $1/2$, then $129^3p^{32000}<1$.",
    reads={"p": P_READING},
)
def union_bound_from_block(p: Real):
    r'''The $32000$ pairs are $32$ blocks of $1000$, so all are bad with probability at most $2^{-32}$, and $129^3<2^{32}$; see
    \texttt{thirty\_two\_blocks\_bound}.'''
    require(p >= 0, p**1000 <= Rational(1, 2))
    return sp.Mul(_pow(129, 3), _pow(p, 32000), evaluate=False) < 1


@PeriodicTilingResidueBounds.theorems(
    source="",
    doc="The arithmetic of \\texttt{union\\_bound\\_from\\_block} on the block probability $y=p^{1000}$: if $0\\le y\\le1/2$ then $129^3y^{32}<1$.",
    reads={"y": "the probability that a block of 1000 pairs are all bad, which is $p^{1000}$"},
)
def thirty_two_blocks_bound(y: Real):
    r'''This statement is \emph{not} a declaration of the Lean file: it is what \texttt{union\_bound\_from\_block} comes to when
    $p^{32000}$ is written as $(p^{1000})^{32}$ (\texttt{pow\_mul}), which is how its proof proceeds. The chain below is the
    \texttt{calc} of that proof.'''
    require(y >= 0, y <= Rational(1, 2))
    return sp.Mul(_pow(129, 3), _pow(y, 32), evaluate=False) < 1


@PeriodicTilingResidueBounds.theorems(
    lean='''theorem residue_union_bound_32000 :
    (129 : ℚ) ^ 3 * ((15598 : ℚ) / 15625) ^ 32000 < 1''',
    source="OAI.PeriodicTilingThree.residue_union_bound_32000",
    doc="The union bound for the application: $129^3(15598/15625)^{32000}<1$.",
)
def residue_union_bound_32000():
    r'''At most $129^3$ requested differences, each with $32000$ pairs; the module \texttt{ResidueTransversal.lean} uses exactly this
    inequality to get a colouring that realises every requested difference. Evaluated exactly, the left side is about
    $1.98\times10^{-18}$.'''
    # the power has 134000 digits, which sympy would fold into one rational when printing a product, so it is a literal term
    return sp.Lt(sp.Mul(_pow(129, 3), _P32000, evaluate=False), 1, evaluate=False)


@PeriodicTilingResidueBounds.theorems(
    lean='''theorem residue_union_bound_of_large_exponent {k : ℕ} (hk : 32000 ≤ k) :
    (129 : ℚ) ^ 3 * ((15598 : ℚ) / 15625) ^ k < 1''',
    source="OAI.PeriodicTilingThree.residue_union_bound_of_large_exponent",
    doc="More pairs only help: $129^3(15598/15625)^k<1$ for every $k\\ge32000$.",
    reads={"k": K_READING},
)
def residue_union_bound_of_large_exponent(k: Nat):
    r'''The bad probability $15598/15625$ is below 1, so its powers decrease with the exponent (\texttt{pow\_le\_pow\_of\_le\_one}):
    $129^3p^k\le129^3p^{32000}$ for $k\ge32000$, and the right side is below 1 by \texttt{residue\_union\_bound\_32000}.'''
    require(k >= 32000)
    return sp.Mul(_pow(129, 3), _pow(_P, k), evaluate=False) < 1


@PeriodicTilingResidueBounds.theorems(
    lean='''theorem residue_cube_64 : (64 : ℕ) ^ 3 = 262144''',
    source="OAI.PeriodicTilingThree.residue_cube_64",
    doc="$64^3=262144$.",
)
def residue_cube_64():
    r'''The number of residues modulo $64$ in dimension three.'''
    return sp.Eq(_pow(64, 3), 262144, evaluate=False)


PeriodicTilingResidueBounds.results.declare(
    "residue_matching_budget_64",
    lean='''theorem residue_matching_budget_64 : 32000 ≤ (64 : ℕ) ^ 3 / 3 - 1''',
    doc=r"$32000\le\lfloor64^3/3\rfloor-1$ (natural-number division and subtraction): there is room for $32000$ disjoint pairs of residues.",
    source="OAI.PeriodicTilingThree.residue_matching_budget_64",
)


# connecting concepts: later used to generate LaTeX that explains how the pieces fit
PeriodicTilingResidueBounds.theorems.one_sub_pow_mul_linear_le_one.connect(
    PeriodicTilingResidueBounds.theorems.one_sub_pow_induction_step, note="its inductive step, written out")
PeriodicTilingResidueBounds.theorems.residue_block_bound.connect(
    PeriodicTilingResidueBounds.theorems.scalar_half_bound, PeriodicTilingResidueBounds.theorems.one_sub_pow_mul_linear_le_one)
PeriodicTilingResidueBounds.theorems.union_bound_from_block.connect(
    PeriodicTilingResidueBounds.theorems.thirty_two_blocks_bound, note="its closing arithmetic, written on the block probability")
PeriodicTilingResidueBounds.theorems.residue_failure_block_bound.connect(PeriodicTilingResidueBounds.theorems.residue_block_bound)
PeriodicTilingResidueBounds.theorems.residue_union_bound_32000.connect(
    PeriodicTilingResidueBounds.theorems.union_bound_from_block, PeriodicTilingResidueBounds.theorems.residue_failure_block_bound)
PeriodicTilingResidueBounds.theorems.residue_union_bound_of_large_exponent.connect(
    PeriodicTilingResidueBounds.theorems.residue_union_bound_32000)


# The `succ` case of the induction in `one_sub_pow_mul_linear_le_one`, step for step.
@PeriodicTilingResidueBounds.theorems.one_sub_pow_induction_step.chain(title="The proof is this chain:")
def one_sub_pow_induction_step_chain(c, r):
    return chain(
        (1 - c)**(r + 1) * (1 + (r + 1) * c),
        eq(sp.Mul(_pow(1 - c, r), sp.Mul(1 - c, 1 + (r + 1) * c, evaluate=False), evaluate=False),
           r"\texttt{pow\_succ}, \texttt{ring}: peel one factor $1-c$ off the power"),
        le((1 - c)**r * (1 + r * c),
           r"\texttt{mul\_le\_mul\_of\_nonneg\_left}, because $(1-c)^r\ge0$ and "
           r"$(1-c)(1+(r+1)c)=1+rc-(r+1)c^2\le1+rc$ (\texttt{hstep})"),
        le(1, r"the induction hypothesis"),
    )


# The `h` / `htwo` argument of `scalar_half_bound`.
@PeriodicTilingResidueBounds.theorems.scalar_half_bound.chain(title="The proof is this chain:")
def scalar_half_bound_chain(x, c):
    return chain(
        x,
        eq(2 * x / 2, r"\texttt{ring}"),
        le(x * (1 + 1000 * c) / 2,
           r"\texttt{mul\_le\_mul\_of\_nonneg\_left}: $2\le1+1000c$ (\texttt{htwo}, from $c\ge1/1000$) and $x\ge0$"),
        le(Rational(1, 2), r"the hypothesis $x(1+1000c)\le1$"),
    )


# The block bound at c = 27/15625.
@PeriodicTilingResidueBounds.theorems.residue_failure_block_bound.chain(title="The proof is this chain:")
def residue_failure_block_bound_chain():
    return chain(
        _pow(_P, 1000),
        eq(_pow(sp.Add(1, -Rational(27, 15625), evaluate=False), 1000), r"\texttt{heq}: $1-27/15625=15598/15625$ (\texttt{norm\_num})"),
        le(Rational(1, 2), r"the block bound at $c=27/15625$; its hypotheses are \texttt{norm\_num}",
           using=PeriodicTilingResidueBounds.theorems.residue_block_bound),
    )


# The calc of `union_bound_from_block` runs on p^32000 = (p^1000)^32 <= (1/2)^32. It is checked here on y = p^1000, the
# probability that a block is all bad: sympy's `factor` cannot be run on a polynomial of degree 32000 in p (it exhausts
# memory), and nothing else in the checker would decide the step symbolically.
@PeriodicTilingResidueBounds.theorems.thirty_two_blocks_bound.chain(title="The proof is this chain, the final calc of the Lean proof read on $y=p^{1000}$:")
def thirty_two_blocks_bound_chain(y):
    return chain(
        sp.Mul(_pow(129, 3), _pow(y, 32), evaluate=False),
        le(sp.Mul(_pow(129, 3), _pow(Rational(1, 2), 32), evaluate=False),
           r"\texttt{pow\_le\_pow\_left}: $0\le y\le1/2$ (the step \texttt{hblocks} of the Lean proof)"),
        lt(1, r"\texttt{norm\_num}: $129^3=2146689<2^{32}=4294967296$"),
    )


# a category that depends on another one: lookups see the inherited items, and the links are made from the labels
@lean.category(
    r'''
    The union bound %s rests on the block bound %s, which is Bernoulli's inequality %s followed by a halving step.
    All three are declared in [[PeriodicTilingResidueBounds]], which this class extends; here they are reached as attributes of the
    subclass: \texttt{PeriodicTilingResidueSummary.residue\_union\_bound\_32000} is the item of that class.
    ''' % (PeriodicTilingResidueBounds.residue_union_bound_32000.latex_ref(page=True),
           PeriodicTilingResidueBounds.residue_block_bound.latex_ref(page=True),
           PeriodicTilingResidueBounds.one_sub_pow_mul_linear_le_one.latex_ref(page=True)),
)
class PeriodicTilingResidueSummary(PeriodicTilingResidueBounds):
    pass
