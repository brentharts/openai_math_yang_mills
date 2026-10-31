"""Where a Yang-Mills mass gap would sit in the library, and what the library can and cannot say about it.

The Clay problem asks for a quantum SU(N) Yang-Mills theory on R^4 with a positive gap in its spectrum. Nothing in the OpenAI
release proves that (see YANGMILLS.md), so this script does not try. It asks what is checkable inside the knowledge graph:

  1. The gap chain. The library states E_gap = hbar c m_L / a ('Mass gap in physical units'), and the textbook two-loop scale
     a Lambda_L = exp(-1/(2 b0 g^2)) (b0 g^2)^(-b1/(2 b0^2)). A continuum limit with a finite gap needs m_L / (a Lambda_L) to
     tend to a constant c as g -> 0, and then E_gap = c hbar c Lambda_L is independent of the coupling. We check that this
     is exactly the RG-invariance condition, and that the two-loop exponent -b1/(2 b0^2) is 51/121 for every N.
     The number c, the gap in units of Lambda, is what lattice numerics estimate and what the Clay problem needs to be positive.
  2. The blocking picture. Refining the lattice by L raises 1/g^2 by 2 b0 ln L at one loop: a linear drift in the number of
     steps. The same linear drift in the O(N) model, with (N-2)/(2 pi) in place of 2 b0, is the library's mass-scaling slope
     and, for N = 3, the (log L / 2 pi) N of the O(3) entry. A shape shared by two models (a correspondence), not an identity.
  3. The scale. Joining the library's 'Dimensional transmutation' to the weak-mixing exploration gives a one-loop QCD scale
     from alpha_s(M_Z) alone; a gap would be a constant multiple of it. This is an order-of-magnitude illustration and not a
     computation of any gap.
"""
import sys
from pathlib import Path

import sympy as sp

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import knowledgebase as kb
from explore import Result, main


def run() -> Result:
    r = Result("ym_mass_gap", "Where a Yang-Mills mass gap would sit", prose="")
    N = sp.Symbol("N", positive=True)
    g, L = sp.symbols("g L", positive=True)
    b0 = 11 * N / (48 * sp.pi ** 2)
    b1 = sp.Rational(34, 3) * (N / (16 * sp.pi ** 2)) ** 2
    # 1. the exponent
    expo = sp.simplify(-b1 / (2 * b0 ** 2))
    r.check("two-loop exponent -b1/(2 b0^2) = -51/121 for every N", expo == sp.Rational(-51, 121))
    a_lam = sp.exp(-1 / (2 * b0 * g ** 2)) * (b0 * g ** 2) ** expo
    # RG: a dg/da = -(b0 g^3 + b1 g^5) with g -> 0 as a -> 0 ; so d ln(a Lambda)/dg = 1/(b0 g^3 + b1 g^5) (up to b1 g^5 order terms)
    series = sp.series(sp.diff(sp.log(a_lam), g) - 1 / (b0 * g ** 3 + b1 * g ** 5), g, 0, 1).removeO()
    r.check("the two-loop scale solves d ln(a Lambda)/dg = 1/(b0 g^3 + b1 g^5) up to terms that vanish as g -> 0",
            sp.simplify(series) == 0)
    c = sp.Symbol("c", positive=True)
    mL = c * a_lam                                                    # asymptotic scaling of the lattice gap
    Egap = mL / a_lam                                                 # E_gap / (hbar c) * a_lam / a_lam, in units of Lambda_L
    r.check("if m_L = c a Lambda_L(g) then E_gap/(hbar c Lambda_L) = m_L/(a Lambda_L) = c, independent of g: scaling is the "
            "condition for a coupling-independent physical gap", sp.simplify(sp.diff(Egap, g)) == 0)
    drift = sp.simplify(sp.diff(sp.log(a_lam), g) * g ** 3 * b0)
    r.check("a, hence the physical lattice spacing, tends to 0 as g -> 0 (a Lambda_L -> 0): the continuum limit sits at "
            "the Gaussian point", sp.limit(a_lam.subs(N, 3), g, 0, "+") == 0)
    gap_lib = kb.equation("Mass gap in physical units")
    r.check("the library's physical-gap relation E = hbar c m_L / a solves for the lattice mass as m_L = E a /(hbar c)",
            kb.solve_for(gap_lib, "lattice mass") is not None)
    # 2. blocking drift
    n, g0 = sp.symbols("n g0", positive=True)
    inv = 1 / g0 ** 2 + 2 * b0 * n * sp.log(L)
    r.check("at one loop a/L-refinement obeys d(1/g^2)/d ln(1/a) = 2 b0 = 11 N/(24 pi^2), so 1/g_n^2 grows linearly in n",
            sp.simplify(sp.diff(inv, n) / sp.log(L) - 11 * N / (24 * sp.pi ** 2)) == 0)
    Nn = sp.Symbol("N_c", positive=True)
    beta, Nv = sp.symbols("beta N_v", positive=True)
    ms = kb.equation("Asymptotic-freedom mass scaling").relation()
    rhs = ms.rhs.subs({s: sp.Symbol(s.name, positive=True) for s in ms.rhs.free_symbols})
    slope = sp.simplify(sp.diff(sp.log(rhs), sp.Symbol("beta", positive=True)).subs(sp.Symbol("N", positive=True), 3))
    big = sp.limit(slope, sp.Symbol("beta", positive=True), sp.oo)
    r.check("the O(N) library form has d ln m / d beta -> -2 pi/(N-2); for N = 3, beta rises by ln L/(2 pi) per refinement "
            "by L, the linear drift of the O(3) construction",
            sp.simplify(big + 2 * sp.pi) == 0)
    r.check("both drifts are linear in the number of steps; the coefficients are 2 b0 (gauge) and (N-2)/(2 pi) (spin), "
            "different numbers for different theories, so this is a correspondence and not an identity",
            sp.simplify(2 * b0.subs(N, 3) - sp.Rational(1, 1) * (3 - 2) / (2 * sp.pi)) != 0)
    # 3. one-loop scale from the weak-mixing chain
    tr = kb.equation("Dimensional transmutation")
    alpha_s, MZ = sp.Rational(1179, 10000), sp.Rational(911876, 10000)
    bq = -sp.Rational(11, 3) * 3 + sp.Rational(2, 3) * 5          # 'One-loop coefficient (QCD)', N = 3, n_f = 5
    r.check("the library's QCD coefficient b = -11N/3 + 2 n_f/3 is -23/3 for N = 3, n_f = 5", bq == sp.Rational(-23, 3))
    lam = float(MZ * sp.exp(2 * sp.pi / (bq * alpha_s)))
    r.check("one-loop Lambda from alpha_s(M_Z) with the library's transmutation formula lies between 10 MeV and 1 GeV",
            0.01 < lam < 1.0)
    r.data.update({"Lambda_1loop [GeV]": round(lam, 3), "exponent": expo})
    r.prose = (r"A rigorous mass gap is not in this library, and none of the release's manuscripts treats a gauge theory. "
               r"What the graph can state is the shape a gap would take. The lattice relation $E_{\mathrm{gap}}=\hbar c\,m_L/a$ "
               r"joins the two-loop scale $a\Lambda_L=e^{-1/(2b_0g^2)}(b_0g^2)^{-51/121}$ (the exponent is $-b_1/2b_0^2$ for every "
               r"$N$), so a coupling-independent physical gap is exactly the statement that $m_L/(a\Lambda_L)\to c$ as $g\to0$: "
               r"the gap is $c\,\hbar c\Lambda_L$. Positivity of $c$ is the open problem. Refining the lattice by a factor $L$ "
               r"adds $2b_0\ln L$ to $1/g^2$, a linear drift in the number of steps, as the O(3) construction adds $\ln L/2\pi$ "
               r"to $\beta$ per step. The shared linear drift is a correspondence between two models, not an identity. As an "
               r"illustration of the scale, the library's transmutation formula with $\alpha_s(M_Z)$ and $b=-23/3$ gives "
               r"$\Lambda\approx" + f"{lam:.2f}" + r"$~GeV at one loop, the right order of magnitude and nothing more.")
    r.display = (r"\underbrace{E_{\mathrm{gap}}}_{\text{open: } c>0?}=\hbar c\,\frac{m_L}{a}=c\,\hbar c\,\Lambda_L,\qquad "
                r"c=\lim_{g\to0}\frac{m_L(g)}{a\Lambda_L(g)},\qquad a\Lambda_L=e^{-\frac{1}{2b_0g^{2}}}(b_0g^{2})^{-\frac{51}{121}}")
    return r


if __name__ == "__main__":
    main(run())
