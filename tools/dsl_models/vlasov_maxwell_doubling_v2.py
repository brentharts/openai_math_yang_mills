"""Maxwell/ExtendedField.lean (the doubling-time lemmas) of the OpenAI release in the V2 syntax: a category with explicit Lean
statements, connected items, and LaTeX prose built from the items themselves.

The file belongs to family 362 (global smoothness for the relativistic Vlasov-Maxwell system: classical matter coupled to the
Maxwell field, the abelian analogue of a gauge theory). Its first lemmas are an argument about times: if a quantity doubles from
level 2^(N+n) to level 2^(N+n+1) in a time of length Delta_n that obeys a certain inequality, then Delta_n is at least c/(n+1),
and gaps of that size cannot all fit below a finite time, because the harmonic series diverges. The two arithmetic lemmas are
items with bodies; the two lemmas that conclude `False` (the DSL states relations) are results with their Lean text.

Not modelled: the rest of the file (the extended vector field, its trace, and the volume-preserving flow)."""
from leandsl import lean, Nat, Real, require, Rational, wikipedia, chain, le, eq
import sympy as sp
from oai_tree import family   # the generated class of family 362 of the release

# one reading of each letter for every item here; connect() refuses items that read a letter differently
N_READING = "the offset of the dyadic levels: the levels are $2^{N+n}$"
N_INDEX_READING = "the index of the doubling step, a natural number (the step from level $2^{N+n}$ to the next)"
M_READING = "a positive constant: the coefficient of $\\sqrt{\\Delta}$ in the one-step inequality"
A_READING = "a positive constant: the coefficient of $\\Delta$ times the logarithm in the one-step inequality"
C_READING = "a positive constant for the logarithm: in the application $C=\\log1026+N\\log2+\\log2$"
D_READING = "the length of one doubling step, a time interval of length at most $t_{n+1}-t_n$"

_MAIN = family("362").results["OAI.RVM.global_classical_solution"]   # the family's Comparator challenge, a generated Item


@lean.category(
    r'''
    \textbf{Analysis/VlasovMaxwell/Maxwell/ExtendedField.lean} of the OpenAI release, written in the leandsl DSL. It is a module
    behind the Lean proof of the family's main theorem %s (the import analysis finds it in the family's closure); the family
    itself is [[F362|family 362]]. The Vlasov--Maxwell system couples charged matter to the Maxwell field, the classical abelian
    case of a gauge theory. The first part of the file, modelled here, is the time argument that rules out a blow-up in
    finite time: a quantity whose level doubles at step $n$ needs a time of at least $c/(n+1)$ for that step
    ([[VlasovMaxwellDoublingTimes.doubling_interval_lower_bound]], fed by
    [[VlasovMaxwellDoublingTimes.doubling_log_bound]]), and gaps of that size add up to infinity
    ([[VlasovMaxwellDoublingTimes.no_bounded_harmonic_gaps]]), so the doubling times cannot all stay below a finite time
    ([[VlasovMaxwellDoublingTimes.no_bounded_doubling_times]]). The module \texttt{UniformPopulation.lean} applies the same
    lemmas to the dyadic levels $R_n=2^{N+n}$. The rest of the file (the extended vector field, its trace and the
    volume-preserving flow) is not modelled.
    ''' % _MAIN.latex_ref(page=True),
    title="Vlasov--Maxwell doubling times",
    refs=[wikipedia("Vlasov equation"), wikipedia("Harmonic series (mathematics)")],
)
class VlasovMaxwellDoublingTimes(family("362")):
    r'''The file proves two inequalities, $\log(2+1024\cdot2^{N+n})\le C(n+1)$ and $c/(n+1)\le\Delta$, and two statements that end in
    \texttt{False}, one about gaps $c/(n+1)$ and one about doubling times; the second is built from the first and the two inequalities.'''


@VlasovMaxwellDoublingTimes.lemmas(
    lean='''lemma doubling_log_bound (N n : ℕ) :
    Real.log (2 + 1024 * (2 : ℝ)^(N+n)) ≤
      (Real.log 1026 + N * Real.log 2 + Real.log 2) * (n+1)''',
    source="OAI.RVM.doubling_log_bound",
    doc="The logarithm grows linearly in the step: $\\log(2+1024\\cdot2^{N+n})\\le(\\log1026+N\\log2+\\log2)(n+1)$.",
    refs=[wikipedia("Logarithm")],
    reads={"N": N_READING, "n": N_INDEX_READING},
)
def doubling_log_bound(N: Nat, n: Nat):
    r'''The level $2^{N+n}$ has a logarithm that is linear in $n$, so $\log(2+1024\cdot2^{N+n})$ is at most a constant times $n+1$.'''
    return sp.log(2 + 1024 * 2**(N + n)) <= (sp.log(1026) + N * sp.log(2) + sp.log(2)) * (n + 1)


@VlasovMaxwellDoublingTimes.lemmas(
    lean='''lemma doubling_interval_lower_bound {M A C Δ : ℝ} {n : ℕ}
    (hM : 0 < M) (hA : 0 < A) (hC : 0 < C) (hΔ : 0 ≤ Δ)
    (h : (1 : ℝ)/2 ≤ 1024*M*Real.sqrt Δ + 1024^2*A*Δ*(C*(n+1))) :
    min ((4*1024*M)^2)⁻¹ (4*1024^2*A*C)⁻¹ / (n+1) ≤ Δ''',
    source="OAI.RVM.doubling_interval_lower_bound",
    doc="A step that satisfies the one-step inequality takes a time of at least $c/(n+1)$, "
        "$c=\\min\\bigl((4\\cdot1024M)^{-2},(4\\cdot1024^2AC)^{-1}\\bigr)$.",
    reads={"M": M_READING, "A": A_READING, "C": C_READING, "Delta": D_READING, "n": N_INDEX_READING},
)
def doubling_interval_lower_bound(M: Real, A: Real, C: Real, Delta: Real, n: Nat):
    r'''The Lean proof splits on whether $1024M\sqrt\Delta\ge1/4$. If so, $(4\cdot1024M)^2\Delta\ge1$, which is the first argument of the
    minimum. If not, the inequality gives $1\le4\cdot1024^2AC\,\Delta(n+1)$, which is the second one divided by $n+1$.'''
    require(M > 0, A > 0, C > 0, Delta >= 0,
            Rational(1, 2) <= 1024 * M * sp.sqrt(Delta) + 1024**2 * A * Delta * (C * (n + 1)))
    return sp.Min(1 / (4 * 1024 * M)**2, 1 / (4 * 1024**2 * A * C)) / (n + 1) <= Delta


VlasovMaxwellDoublingTimes.results.declare(
    "no_bounded_harmonic_gaps",
    lean='''lemma no_bounded_harmonic_gaps {t : ℕ → ℝ} {c T : ℝ}
    (hc : 0 < c) (ht : ∀ n, t n ≤ T)
    (hg : ∀ n, c/(n+1) ≤ t (n+1)-t n) : False''',
    doc=r"Times $t_n\le T$ cannot gain $c/(n+1)$ at every step: the gaps would sum to $c\sum_n1/(n+1)=\infty$ "
        r"while $\sum_{i<n}(t_{i+1}-t_i)=t_n-t_0\le T-t_0$.",
    source="OAI.RVM.no_bounded_harmonic_gaps",
    refs=[wikipedia("Harmonic series (mathematics)")],
)

VlasovMaxwellDoublingTimes.results.declare(
    "no_bounded_doubling_times",
    lean='''lemma no_bounded_doubling_times {t Δ : ℕ → ℝ} {M A T : ℝ} (N : ℕ)
    (hM : 0 < M) (hA : 0 < A) (ht : ∀ n, t n ≤ T)
    (hΔ : ∀ n, 0 ≤ Δ n) (hgap : ∀ n, Δ n ≤ t (n+1)-t n)
    (hinc : ∀ n, (1 : ℝ)/2 ≤ 1024*M*Real.sqrt (Δ n) +
      1024^2*A*Δ n*Real.log (2+1024*(2 : ℝ)^(N+n))) : False''',
    doc=r"Doubling times cannot stay below a finite time $T$: each step of length $\Delta_n\le t_{n+1}-t_n$ obeys the one-step "
        r"inequality, so $\Delta_n\ge c/(n+1)$ and the gaps are harmonic.",
    source="OAI.RVM.no_bounded_doubling_times",
)


@VlasovMaxwellDoublingTimes.theorems(
    source="",
    doc="The inequality inside the proof of \\texttt{no\\_bounded\\_doubling\\_times}: every doubling step satisfies "
        "$c/(n+1)\\le\\Delta$ with $c=\\min\\bigl((4\\cdot1024M)^{-2},(4\\cdot1024^2AC)^{-1}\\bigr)$, $C=\\log1026+N\\log2+\\log2$.",
    reads={"M": M_READING, "A": A_READING, "Delta": D_READING, "N": N_READING, "n": N_INDEX_READING},
)
def doubling_gap_lower_bound(M: Real, A: Real, Delta: Real, N: Nat, n: Nat):
    r'''This statement is \emph{not} a declaration of the Lean file: it is the step that the proof of \texttt{no\_bounded\_doubling\_times}
    carries out for every $n$, with the constant $C=\log1026+N\log2+\log2$ of \texttt{doubling\_log\_bound}. Its hypothesis is the
    one-step inequality \texttt{hinc}. The chain below states how the two lemmas combine.'''
    C = sp.log(1026) + N * sp.log(2) + sp.log(2)
    require(M > 0, A > 0, Delta >= 0,
            Rational(1, 2) <= 1024 * M * sp.sqrt(Delta) + 1024**2 * A * Delta * sp.log(2 + 1024 * 2**(N + n)))
    return sp.Min(1 / (4 * 1024 * M)**2, 1 / (4 * 1024**2 * A * C)) / (n + 1) <= Delta


# connecting concepts: later used to generate LaTeX that explains how the pieces fit
VlasovMaxwellDoublingTimes.results.no_bounded_doubling_times.connect(
    VlasovMaxwellDoublingTimes.results.no_bounded_harmonic_gaps,
    VlasovMaxwellDoublingTimes.lemmas.doubling_interval_lower_bound,
    VlasovMaxwellDoublingTimes.lemmas.doubling_log_bound)
VlasovMaxwellDoublingTimes.theorems.doubling_gap_lower_bound.connect(
    VlasovMaxwellDoublingTimes.lemmas.doubling_log_bound, note="its logarithm bound turns the one-step inequality into the hypothesis of the next lemma")


# The first calc of `doubling_log_bound`: bound the sum by a product, take logarithms, then bound by a multiple of n+1.
@VlasovMaxwellDoublingTimes.lemmas.doubling_log_bound.chain(title="The proof is this chain:")
def doubling_log_bound_chain(N, n):
    return chain(
        sp.log(2 + 1024 * 2**(N + n)),
        le(sp.log(1026 * 2**(N + n)),
           r"\texttt{Real.log\_le\_log}: $2+1024\cdot2^{N+n}\le1026\cdot2^{N+n}$ because $1\le2^{N+n}$ (\texttt{one\_le\_pow})"),
        eq(sp.log(1026) + (N + n) * sp.log(2), r"\texttt{Real.log\_mul} and \texttt{Real.log\_pow}"),
        le((sp.log(1026) + N * sp.log(2) + sp.log(2)) * (n + 1),
           r"\texttt{nlinarith}: the difference is $n\log1026+Nn\log2+\log2\ge0$, as $\log1026,\log2\ge0$"),
    )


# How the two lemmas combine in each step of `no_bounded_doubling_times`.
@VlasovMaxwellDoublingTimes.theorems.doubling_gap_lower_bound.chain(title="The step is this chain:")
def doubling_gap_lower_bound_chain(M, A, Delta, N, n):
    C = sp.log(1026) + N * sp.log(2) + sp.log(2)
    return chain(
        sp.Min(1 / (4 * 1024 * M)**2, 1 / (4 * 1024**2 * A * C)) / (n + 1),
        le(Delta,
           r"applied with this $C$: its hypothesis is \texttt{hinc} with $\log(2+1024\cdot2^{N+n})$ "
           r"replaced by the larger $C(n+1)$, by \texttt{doubling\_log\_bound}",
           using=VlasovMaxwellDoublingTimes.lemmas.doubling_interval_lower_bound),
    )
