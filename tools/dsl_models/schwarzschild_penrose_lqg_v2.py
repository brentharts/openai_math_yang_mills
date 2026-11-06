"""Relativity/CKS/Schwarzschild{Area,Charges}.lean of the OpenAI release in the V2 syntax, and its arithmetic read against the
loop-quantum-gravity area spectrum.

The release part. Family 260 proves spacetime Penrose inequalities: the ADM or Bondi mass is at least sqrt(A_min/16 pi), where
A_min is the least area of a surface that encloses the boundary. The modules here show that Schwarzschild data of mass m
is an equality case: the least enclosing area is 16 pi m^2 and the Bondi charge (m, 0) has invariant mass m. The Lean statements are
about Riemannian volumes (values in ENNReal) and a Bondi charge, which the DSL has no types for, so they are results with their
Lean text, except `numerical_equality`, whose real-number content (the square roots agree because the area is 16 pi m^2) is an item.

The quantum-gravity part (NOT in the release). Spectre/Chair44 and the SU(2) ladder paper (Hartshorn and Romanov, Section I and J)
reads the Penrose bound at the area of N spin-1/2 punctures, A_N = 4 sqrt(3) pi gamma N l_P^2 with the Immirzi parameter
gamma = ln 2/(pi sqrt 3), so A_N = 4 N ln 2 l_P^2. Those items are marked derived: they are algebra on the classical
bound and the quantised area, and the paper itself calls the composition a heuristic (a classical inequality about initial
data against an eigenvalue of an area operator). Nothing in the Lean files says anything about spin networks."""
from leandsl import lean, Nat, Real, require, Rational, wikipedia, chain, le, lt, eq
import sympy as sp
from oai_tree import family   # the generated class of family 260 of the release

M_READING = "the mass, in units with $G=c=1$: the Bondi (or ADM) mass of the Schwarzschild data, a positive real number"
N_READING = "the number of spin-$\\tfrac12$ punctures that pierce the surface, a natural number"
J_READING = "the angular momentum, in units with $\\hbar=1$, a positive real number"
W_READING = "$\\sqrt{M^2-a^2}$, so that the horizons of Kerr are at $r_\\pm=M\\pm w$"

_GAMMA = sp.log(2) / (sp.pi * sp.sqrt(3))   # the Immirzi parameter of the spin-1/2 puncture count
_MAIN = family("260").results["OAI.CKSMain.schwarzschild_equality_examples"]   # the family's Comparator challenge, a generated Item


@lean.category(
    r'''
    \textbf{Geometry/Relativity/CKS/SchwarzschildArea.lean} and \textbf{SchwarzschildCharges.lean} of the OpenAI release, written in
    the leandsl DSL, and the arithmetic that [[F260|family 260]] (spacetime Penrose inequalities) and the loop-quantum-gravity
    area spectrum have in common. The Lean modules behind %s show that Schwarzschild data is an equality case of the
    Penrose inequality $m\ge\sqrt{A_{\min}/16\pi}$: the least enclosing area is $16\pi m^2$
    ([[SchwarzschildPenroseLQG.full_minimum_enclosing_area]], from the area of the scaled sphere
    [[SchwarzschildPenroseLQG.scaled_sphere_area]]) and the Bondi charge $(m,0)$ has invariant mass $m$, so both sides
    are the same square root ([[SchwarzschildPenroseLQG.numerical_equality]]). The second half of the page is \emph{not} in the
    release. In loop quantum gravity the area of a surface pierced by $N$ spin-$\tfrac12$ links is a quantised eigenvalue
    $A_N=4\sqrt3\pi\gamma N\ell_P^2$; with the Immirzi parameter $\gamma=\ln2/(\pi\sqrt3)$ that is $4N\ln2\,\ell_P^2$, and the
    Penrose bound at that area is $M/M_P\ge\sqrt{N\ln2/4\pi}$ ([[SchwarzschildPenroseLQG.puncture_penrose_mass]]).
    The Kerr area bound $A_+\ge8\pi|J|$ then asks for $N\ge2\pi|J|/\ln2$ punctures
    ([[SchwarzschildPenroseLQG.puncture_count_threshold]]), and the exact gap
    $M^2-A_+/16\pi=Mr_-/2$ of Kerr ([[SchwarzschildPenroseLQG.kerr_penrose_gap]]) is the algebra behind the statement that the
    Penrose inequality is strict for every rotating black hole. These three are marked as derived: they combine a classical
    inequality about initial data with an eigenvalue of an operator, which is a heuristic.
    ''' % _MAIN.latex_ref(page=True),
    title="Schwarzschild, Penrose and the area spectrum",
    refs=[wikipedia("Penrose inequality"), wikipedia("Loop quantum gravity"), wikipedia("Immirzi parameter")],
)
class SchwarzschildPenroseLQG(family("260")):
    r'''Two results are about the area of the Schwarzschild sphere, one real identity is the equality case, and three derived items put the
    quantised area into the bound.'''


SchwarzschildPenroseLQG.results.declare(
    "scaled_sphere_area",
    lean='''lemma scaled_sphere_area {m : ℝ} (hm : 0 < m) :
    riemannianVolume (boundaryScaledMetric hm).toContinuousRiemannianMetric univ =
      ENNReal.ofReal (16 * Real.pi * m^2)''',
    doc=r"The Schwarzschild horizon, a round sphere of radius $2m$, has area $4\pi(2m)^2=16\pi m^2$.",
    source="OAI.CKSSchwarzschild.scaled_sphere_area",
    refs=[wikipedia("Schwarzschild radius")],
)

SchwarzschildPenroseLQG.results.declare(
    "full_minimum_enclosing_area",
    lean='''theorem full_minimum_enclosing_area {m : ℝ} (hm : 0 < m) :
    CKSFullCutArea.minEnclosingArea (smoothMetric hm) = ENNReal.ofReal (16 * Real.pi * m^2)''',
    doc=r"No surface that encloses the horizon has less area than the horizon itself: the least enclosing area is $16\pi m^2$.",
    source="OAI.CKSSchwarzschild.full_minimum_enclosing_area",
)


@SchwarzschildPenroseLQG.lemmas(
    lean='''lemma numerical_equality {m : ℝ} (hm : 0 < m) :
    Real.sqrt ((bondiCharge (cksData hm).massAspect).1^2 -
      ‖(bondiCharge (cksData hm).massAspect).2‖^2) =
    Real.sqrt (CKSSourceExterior.minimumEnclosingArea (smoothMetric hm) / (16*Real.pi))''',
    source="OAI.CKSSchwarzschild.numerical_equality",
    doc="The equality case of the Penrose inequality: for the Bondi charge $(m,0)$ and the least enclosing area $16\\pi m^2$, "
        "$\\sqrt{m^2-0}=\\sqrt{16\\pi m^2/16\\pi}$.",
    reads={"m": M_READING},
)
def numerical_equality(m: Real):
    r'''The Lean statement is about the Bondi charge and the least enclosing area of the data; the proof rewrites them with
    \texttt{bondi\_charge} $=(m,0)$ and \texttt{source\_minimum\_area} $=16\pi m^2$, after which only this identity of real numbers is left.'''
    require(m > 0)
    return sp.Eq(sp.sqrt(m**2 - 0**2), sp.sqrt(16 * sp.pi * m**2 / (16 * sp.pi)), evaluate=False)


@SchwarzschildPenroseLQG.theorems(
    source="",
    doc="Derived, not in the release: the Penrose bound $\\sqrt{A_N/16\\pi}$ at the area $A_N=4\\sqrt3\\pi\\gamma N$ of $N$ spin-$\\tfrac12$ "
        "punctures (Planck units, $\\gamma=\\ln2/(\\pi\\sqrt3)$) is $\\sqrt{N\\ln2/4\\pi}$.",
    reads={"N": N_READING},
)
def puncture_penrose_mass(N: Nat):
    r'''\emph{Not a declaration of the Lean files, and a heuristic composition:} the left side is the right-hand side of the Penrose
    inequality evaluated at an eigenvalue of the area operator of loop quantum gravity, for the value $\gamma=\ln2/(\pi\sqrt3)$ of the
    Immirzi parameter that one black-hole entropy count gives. For one puncture the bound is $0.2349$ Planck masses.'''
    return sp.Eq(sp.sqrt(4 * sp.sqrt(3) * sp.pi * _GAMMA * N / (16 * sp.pi)), sp.sqrt(N * sp.log(2) / (4 * sp.pi)), evaluate=False)


@SchwarzschildPenroseLQG.theorems(
    source="",
    doc="Derived, not in the release: if the quantised area $A_N=4N\\ln2$ is at least the Kerr bound $8\\pi J$, then $N\\ge2\\pi J/\\ln2$.",
    reads={"N": N_READING, "J": J_READING},
)
def puncture_count_threshold(N: Nat, J: Real):
    r'''\emph{Not a declaration of the Lean files, and a heuristic composition:} the classical inequality $A_+\ge8\pi|J|$ for a Kerr horizon,
    asked of the area of $N$ punctures. It gives $N\ge4.53$ at $J=\tfrac12$ and $N\ge9.06$ at $J=1$, so five and ten punctures.'''
    require(J > 0, 4 * N * sp.log(2) >= 8 * sp.pi * J)
    return N >= 2 * sp.pi * J / sp.log(2)


@SchwarzschildPenroseLQG.theorems(
    source="",
    doc="Derived, not in the release: for Kerr, with $r_\\pm=M\\pm w$ and $A_+=8\\pi Mr_+$, the Penrose gap is $M^2-A_+/16\\pi=Mr_-/2$.",
    reads={"M": M_READING, "w": W_READING},
)
def kerr_penrose_gap(M: Real, w: Real):
    r'''\emph{Not a declaration of the Lean files} (family 260 treats Schwarzschild data, and the Kerr papers of entry 264 have no Lean). The
    horizon area of Kerr is $4\pi(r_+^2+a^2)=8\pi Mr_+$, so the gap vanishes exactly when the inner horizon $r_-=M-w$ does.'''
    return sp.Eq(M**2 - 8 * sp.pi * M * (M + w) / (16 * sp.pi), M * (M - w) / 2, evaluate=False)


# connecting concepts: later used to generate LaTeX that explains how the pieces fit
SchwarzschildPenroseLQG.results.full_minimum_enclosing_area.connect(SchwarzschildPenroseLQG.results.scaled_sphere_area)
SchwarzschildPenroseLQG.lemmas.numerical_equality.connect(SchwarzschildPenroseLQG.results.full_minimum_enclosing_area)
SchwarzschildPenroseLQG.theorems.puncture_penrose_mass.connect(
    SchwarzschildPenroseLQG.lemmas.numerical_equality,
    note="the same square root $\\sqrt{A/16\\pi}$, evaluated at a quantised area instead of $16\\pi m^2$")
SchwarzschildPenroseLQG.theorems.puncture_count_threshold.connect(SchwarzschildPenroseLQG.theorems.puncture_penrose_mass)
SchwarzschildPenroseLQG.theorems.kerr_penrose_gap.connect(
    SchwarzschildPenroseLQG.lemmas.numerical_equality, note="the Schwarzschild equality case is the gap with $w=M$, $r_-=0$")


@SchwarzschildPenroseLQG.lemmas.numerical_equality.chain(title="The proof is this chain:")
def numerical_equality_chain(m):
    return chain(
        sp.sqrt(m**2 - 0**2),
        eq(sp.sqrt(m**2), r"\texttt{simp}: the second component of the Bondi charge $(m,0)$ is $0$"),
        eq(sp.sqrt(16 * sp.pi * m**2 / (16 * sp.pi)), r"\texttt{field\_simp}: $16\pi m^2/16\pi=m^2$"),
    )


@SchwarzschildPenroseLQG.theorems.puncture_penrose_mass.chain(title="The composition is this chain:")
def puncture_penrose_mass_chain(N):
    return chain(
        sp.sqrt(4 * sp.sqrt(3) * sp.pi * _GAMMA * N / (16 * sp.pi)),
        eq(sp.sqrt(4 * N * sp.log(2) / (16 * sp.pi)), r"$\gamma=\ln2/(\pi\sqrt3)$, so $4\sqrt3\pi\gamma=4\ln2$"),
        eq(sp.sqrt(N * sp.log(2) / (4 * sp.pi)), r"cancel the factor $4$"),
    )


@SchwarzschildPenroseLQG.theorems.puncture_count_threshold.chain(title="The count is this chain:")
def puncture_count_threshold_chain(N, J):
    return chain(
        2 * sp.pi * J / sp.log(2),
        eq(8 * sp.pi * J / (4 * sp.log(2)), r"multiply numerator and denominator by $4$"),
        le(4 * N * sp.log(2) / (4 * sp.log(2)), r"the hypothesis $8\pi J\le4N\ln2$, divided by $4\ln2>0$"),
        eq(N, r"cancel"),
    )


@SchwarzschildPenroseLQG.theorems.kerr_penrose_gap.chain(title="The gap is this chain:")
def kerr_penrose_gap_chain(M, w):
    return chain(
        M**2 - 8 * sp.pi * M * (M + w) / (16 * sp.pi),
        eq(M**2 - M * (M + w) / 2, r"$8\pi/16\pi=1/2$"),
        eq((M**2 - M * w) / 2, r"expand: $M^2-M(M+w)/2=(M^2-Mw)/2$"),
        eq(M * (M - w) / 2, r"factor $M$"),
    )
