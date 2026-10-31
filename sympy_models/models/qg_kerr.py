"""Entry 264 (strong cosmic censorship near two-ended Kerr data): the Kerr horizon quantities the near-Kerr papers use, and their place beside the area spectrum.

Source: OpenAI, "Generic Future Inextendibility with Square-Integrable Connection Near a Fixed Kerr Spacetime" (September 23 2026; section 2.1 and the
discussion of the fixed rate kappa), with the C^1 and C^2 companions of the same date. The fork read the abstract, introduction, section 1 and section 2.1-2.2 of the
first; it did not read the proofs. What these papers use of Kerr is: mass M > 0 and rotation 0 < |a| < M (reduced to M = 1, 0 < a < 1 by scaling), inner
horizon r_- > 0, and the positive absolute inner surface gravity kappa, which fixes the decay rate of the lapse (a ~ e^{-kappa T}) and the packet amplitude
A_T = e^{-kappa T/4}, larger than the square root of the lapse. This module computes, in exact algebra (G = c = 1, Boyer-Lindquist r_+ > r_- > 0, so M = (r_+ + r_-)/2
and a^2 = r_+ r_-), the quantities around that:

  qg_kerr_horizons   r_+- = M +- sqrt(M^2 - a^2); area A_+- = 4 pi (r_+-^2 + a^2); surface gravities kappa_+- = (r_+- - r_-+)/(2 (r_+-^2 + a^2)); angular velocity
                     Omega = a/(r^2 + a^2); the Smarr formula and the first law at BOTH horizons (the first law as an identity of one-forms in (r_+, r_-));
                     the scaling that reduces a fixed centre to M = 1; the limits a -> 0 (|kappa_-| -> infinity, a^2 |kappa_-| -> M) and a -> M (both -> 0), which are why the
                     papers' constants are not uniform toward zero rotation or extremality.
  qg_kerr_penrose    the Penrose gap of Kerr: M^2 - A_+/(16 pi) = M r_-/2, so the Penrose inequality of entry 260 is strict for Kerr exactly when the inner horizon exists
                     (the condition the near-Kerr papers need); the Christodoulou-Ruffini form M^2 = M_irr^2 + J^2/(4 M_irr^2); the area product A_+ A_- = (8 pi J)^2; the
                     bound A_+ >= 8 pi |J|, with equality only at extremality.
  qg_kerr_quantum    A HEURISTIC composition with the area spectrum of qg_spin_network.py: N spin-1/2 punctures of area 4 sqrt3 pi gamma N l_P^2 against the Kerr bound
                     A >= 8 pi |J| (hbar = 1) gives N >= 2 pi |J| / ln 2 at gamma = ln2/(pi sqrt3): five punctures for J = 1/2, ten for J = 1, and 2 pi/ln 2 = 9.06 for J = 1 is
                     the same number as the Compton threshold of qg_penrose.

Checked: every identity above exactly (sympy), the three limits, and the scalings quoted from the paper (lapse e^{-kappa T}, amplitude e^{-kappa T/4} > sqrt(lapse)).
Not tested: any theorem of entries 264 or 260 (genericity, inextendibility, the double-null comparison, the averaged holonomy argument), the claim that a
rotating quantum horizon is a Kerr horizon, and the identification of the Kerr bound with a bound on spin-network data: the area A >= 8 pi |J| is a statement about
a classical solution, an LQG area eigenvalue is an eigenvalue of an operator, and nothing here connects the two beyond using one symbol A. The exterior
and interior of the papers are perturbations of Kerr, and the quantities here are those of Kerr itself.
"""
import sys
from pathlib import Path

import sympy as sp

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register

NAME = "qg_kerr"
ENTRY = "264"
PAPER = "[K2] OpenAI, Generic Future Inextendibility with Square-Integrable Connection Near a Fixed Kerr Spacetime (preprints, September 23 2026), sections 1 and 2"


def build_and_verify():
    rp, rm, lam, T, v, p = sp.symbols("r_p r_m lambda T v p", positive=True)
    # Boyer-Lindquist: r_+ > r_- > 0, M = (r_+ + r_-)/2, a^2 = r_+ r_-
    M = (rp + rm) / 2
    a = sp.sqrt(rp * rm)
    J = M * a
    Ap, Am = 4 * sp.pi * (rp ** 2 + a ** 2), 4 * sp.pi * (rm ** 2 + a ** 2)
    kp, km = (rp - rm) / (2 * (rp ** 2 + a ** 2)), (rm - rp) / (2 * (rm ** 2 + a ** 2))      # signed
    Op, Om = a / (rp ** 2 + a ** 2), a / (rm ** 2 + a ** 2)
    Ms, as_ = sp.symbols("M a", positive=True)
    root = sp.sqrt(Ms ** 2 - as_ ** 2)

    eq_h = register(NAME, Equation(
        key="qg_kerr_horizons", entry=ENTRY, title="Kerr horizons: radii, areas, surface gravities, Smarr formula and the first law at both horizons",
        source_kind="derived", source=f"derived here from the Kerr family used in {PAPER}; the genericity theorems are the papers'",
        latex_str=(r"r_\pm=M\pm\sqrt{M^{2}-a^{2}},\quad r_+r_-=a^{2},\quad A_\pm=4\pi(r_\pm^{2}+a^{2})=8\pi Mr_\pm,\quad "
                   r"\kappa_\pm=\frac{r_\pm-r_\mp}{2(r_\pm^{2}+a^{2})}=\pm\frac{\sqrt{M^{2}-a^{2}}}{2Mr_\pm},\quad \Omega_\pm=\frac{a}{r_\pm^{2}+a^{2}}"
                   r"%%BR%%M=2\Omega_\pm J+\frac{\kappa_\pm A_\pm}{4\pi},\qquad dM=\frac{\kappa_\pm}{8\pi}\,dA_\pm+\Omega_\pm\,dJ,\qquad J=Ma"
                   r"%%BR%%\frac{|\kappa_-|}{\kappa_+}=\frac{r_+}{r_-}>1,\qquad \kappa_+|\kappa_-|=\frac{M^{2}-a^{2}}{4M^{2}a^{2}},\qquad a^{2}|\kappa_-|\xrightarrow{a\to0}M,\qquad \kappa_\pm\xrightarrow{a\to M}0"),
        notes="Exact in the variables r_+ > r_- > 0 (so that M = (r_+ + r_-)/2 and a = sqrt(r_+ r_-)), with G = c = 1. The first law is checked as an identity of "
              "one-forms: the partial derivatives with respect to r_+ and r_- of both sides agree, at the outer horizon and at the inner one (with the signed kappa_-). "
              "The scaling M -> lambda M, a -> lambda a sends A to lambda^2 A, kappa to kappa/lambda, J to lambda^2 J, so M kappa depends on a/M alone: this is the "
              "reduction of a fixed centre to M = 1 in section 2.1 of the paper. The limits show why the constants of the papers cannot be uniform: |kappa_-| grows like M/a^2 as a -> 0 "
              "(the inner horizon merges with the singularity) and both surface gravities vanish as a -> M. Not tested: that kappa so defined is the rate "
              "the paper's exact profile rescales to (the paper takes it as the center's absolute inner surface gravity), and anything about perturbed data.",
    ))
    check(NAME, "the radii: r_+ + r_- = 2M and r_+ r_- = a^2 for r_+- = M +- sqrt(M^2 - a^2); r_- > 0 exactly when a != 0 (r_- = a^2/r_+)",
          sp.simplify((Ms + root) + (Ms - root) - 2 * Ms) == 0 and sp.simplify((Ms + root) * (Ms - root) - as_ ** 2) == 0, eq=eq_h)
    check(NAME, "the area: 4 pi (r_+ ^2 + a^2) = 8 pi M r_+ and 4 pi (r_-^2 + a^2) = 8 pi M r_- (exact)",
          sp.simplify(Ap - 8 * sp.pi * M * rp) == 0 and sp.simplify(Am - 8 * sp.pi * M * rm) == 0, eq=eq_h)
    check(NAME, "the surface gravities: kappa_+ = sqrt(M^2-a^2)/(2 M r_+) and kappa_- = -sqrt(M^2-a^2)/(2 M r_-) (exact; sqrt(M^2-a^2) = (r_+ - r_-)/2)",
          sp.simplify(kp - (rp - rm) / 2 / (2 * M * rp)) == 0 and sp.simplify(km + (rp - rm) / 2 / (2 * M * rm)) == 0, eq=eq_h)
    check(NAME, "Smarr at both horizons: M = 2 Omega J + kappa A/(4 pi) with the signed kappa (exact)",
          sp.simplify(M - (2 * Op * J + kp * Ap / (4 * sp.pi))) == 0 and sp.simplify(M - (2 * Om * J + km * Am / (4 * sp.pi))) == 0, eq=eq_h)
    first_p = [sp.simplify(sp.diff(M, x) - (kp / (8 * sp.pi) * sp.diff(Ap, x) + Op * sp.diff(J, x))) for x in (rp, rm)]
    first_m = [sp.simplify(sp.diff(M, x) - (km / (8 * sp.pi) * sp.diff(Am, x) + Om * sp.diff(J, x))) for x in (rp, rm)]
    check(NAME, "the first law dM = (kappa/8 pi) dA + Omega dJ holds as an identity of one-forms in (r_+, r_-), at the outer horizon and at the inner horizon (signed kappa_-)",
          first_p == [0, 0] and first_m == [0, 0], eq=eq_h)
    check(NAME, "|kappa_-|/kappa_+ = r_+/r_- > 1: the inner surface gravity exceeds the outer one for every rotating subextremal Kerr; and kappa_+ |kappa_-| = (M^2 - a^2)/(4 M^2 a^2)",
          sp.simplify(-km / kp - rp / rm) == 0 and sp.simplify(-kp * km - ((M ** 2 - a ** 2) / (4 * M ** 2 * a ** 2))) == 0, eq=eq_h)
    scaled = [sp.simplify(f.subs({rp: lam * rp, rm: lam * rm}, simultaneous=True) / f - w) for f, w in ((Ap, lam ** 2), (kp, 1 / lam), (J, lam ** 2), (M, lam))]
    check(NAME, "scaling (r_+, r_-) -> lambda (r_+, r_-): A -> lambda^2 A, kappa -> kappa/lambda, J -> lambda^2 J, M -> lambda M, so M kappa depends only on a/M: "
                "the reduction to M = 1 of section 2.1 loses nothing at fixed a/M", all(s == 0 for s in scaled), eq=eq_h)
    x = sp.Symbol("x", positive=True)
    # a -> 0 : r_- -> 0 at fixed r_+ = 2M (so M = (2M + r_-)/2 -> M); use r_+ = 2 - x, r_- = x -> M = 1, a^2 = x (2 - x)
    a2 = x * (2 - x)
    km_x = -(2 - 2 * x) / (2 * x * 2)       # |kappa_-| = (r_+ - r_-)/(2 r_- (r_+ + r_-)) with r_+ + r_- = 2
    lim0 = sp.limit(a2 * (-km_x), x, 0)
    check(NAME, f"as a -> 0 (r_- = x -> 0, M = 1): a^2 |kappa_-| -> {lim0} = M and |kappa_-| ~ M/a^2 -> infinity; as a -> M (r_+ - r_- -> 0) both kappa_+- -> 0 (exact limits): the papers' "
                "constants depend on the fixed centre and are not uniform in either direction",
          lim0 == 1 and sp.limit(kp.subs({rp: 1 + x, rm: 1 - x}), x, 0) == 0 and sp.limit(km.subs({rp: 1 + x, rm: 1 - x}), x, 0) == 0, eq=eq_h)
    # the scalings of the paper: lapse a ~ e^{-kappa T}, packet amplitude A_T = e^{-kappa T/4} > sqrt(lapse) = e^{-kappa T/2}
    kk = sp.Symbol("kappa", positive=True)
    AT, lapse = sp.exp(-kk * T / 4), sp.exp(-kk * T)
    check(NAME, "the scalings of the paper: with lapse a ~ e^{-kappa T} and amplitude A_T = e^{-kappa T/4}, A_T / sqrt(lapse) = e^{kappa T/4} -> infinity: the amplitude exceeds the square root of the lapse "
                "(the paper: 'This is larger than sqrt(a)'), and A_T -> 0",
          sp.simplify(AT / sp.sqrt(lapse) - sp.exp(kk * T / 4)) == 0 and sp.limit(sp.exp(kk * T / 4), T, sp.oo) == sp.oo and sp.limit(AT, T, sp.oo) == 0, eq=eq_h)
    check(NAME, "blueshift against decay (the standard heuristic, not the paper's argument): a perturbation decaying as v^{-p} meets a blueshift e^{|kappa_-| v}, and v^{-p} e^{|kappa_-| v} -> infinity for every p > 0",
          sp.limit(v ** (-p) * sp.exp(kk * v), v, sp.oo) == sp.oo, eq=eq_h)

    # ---------------------------------------------------------------------------------------------------------
    eq_p = register(NAME, Equation(
        key="qg_kerr_penrose", entry=ENTRY, title="The Penrose gap of Kerr, the Christodoulou-Ruffini form and the area product",
        source_kind="derived", source=f"derived here; the Penrose inequality is the theorem of entry 260 [Pen]; Kerr is the background of {PAPER}",
        latex_str=(r"M^{2}-\frac{A_+}{16\pi}=\frac{Mr_-}{2}>0\ \ (a\neq0),\qquad M^{2}=M_{\rm irr}^{2}+\frac{J^{2}}{4M_{\rm irr}^{2}},\ \ M_{\rm irr}=\sqrt{\tfrac{A_+}{16\pi}}"
                   r"%%BR%%A_+A_-=(8\pi J)^{2},\qquad A_+-8\pi|J|=8\pi M(r_+-a)\ \ge0,\ \ =0\iff r_+=r_-\ (a=M)"),
        notes="The Penrose inequality m >= sqrt(A/16 pi) (entry 260) is an equality for Schwarzschild and strict for every rotating Kerr, with the exact gap M r_-/2 = "
              "M^2 - M_irr^2: the gap vanishes exactly when the inner horizon r_- does, which is the condition (r_- > 0) the near-Kerr papers need. The Christodoulou-Ruffini "
              "formula is the same identity written with J. The area product of the two horizons is (8 pi J)^2. All exact. Not tested: the Penrose inequality itself "
              "(that is the theorem of entry 260), and any relation of this gap to the censorship theorems beyond the shared condition r_- > 0: that coincidence is "
              "stated, not explained.",
    ))
    Mirr2 = Ap / (16 * sp.pi)
    check(NAME, "Penrose gap: M^2 - A_+/(16 pi) = M r_-/2 (exact), positive iff r_- > 0 iff a != 0; zero at a = 0 (Schwarzschild: r_- = 0, equality in the Penrose bound)",
          sp.simplify(M ** 2 - Mirr2 - M * rm / 2) == 0 and sp.limit((M * rm / 2), rm, 0) == 0, eq=eq_p)
    check(NAME, "Christodoulou-Ruffini: M^2 = M_irr^2 + J^2/(4 M_irr^2) with M_irr^2 = A_+/(16 pi) (exact)", sp.simplify(M ** 2 - (Mirr2 + J ** 2 / (4 * Mirr2))) == 0, eq=eq_p)
    check(NAME, "the area product A_+ A_- = (8 pi J)^2 = 64 pi^2 M^2 a^2 (exact)", sp.simplify(Ap * Am - (8 * sp.pi * J) ** 2) == 0, eq=eq_p)
    check(NAME, "A_+ - 8 pi |J| = 8 pi M (r_+ - a) >= 0, with equality exactly at r_+ = r_- (extremal Kerr): r_+ - a = sqrt(r_+)(sqrt(r_+) - sqrt(r_-)) (exact, positive for r_+ > r_- > 0)",
          sp.simplify(Ap - 8 * sp.pi * J - 8 * sp.pi * M * (rp - a)) == 0 and sp.simplify(rp - a - sp.sqrt(rp) * (sp.sqrt(rp) - sp.sqrt(rm))) == 0
          and sp.limit((Ap - 8 * sp.pi * J).subs({rp: 1 + x, rm: 1 - x}), x, 0) == 0, eq=eq_p)
    check(NAME, "subextremality M >= sqrt|J| (a <= M, J = M a <= M^2): at fixed J the mass is at least sqrt|J| M_P, while the Penrose bound gives only M >= sqrt(A/(16 pi)) >= sqrt(|J|/2): "
                "the Kerr bound is the stronger one at fixed J",
          sp.simplify(M ** 2 - J - M * (sp.sqrt(rp) - sp.sqrt(rm)) ** 2 / 2) == 0
          and sp.simplify(Mirr2 - J / 2 - M * sp.sqrt(rp) * (sp.sqrt(rp) - sp.sqrt(rm)) / 2) == 0, eq=eq_p)

    # ---------------------------------------------------------------------------------------------------------
    gamma_val = sp.log(2) / (sp.pi * sp.sqrt(3))
    N, Jv = sp.symbols("N J", positive=True)
    eq_q = register(NAME, Equation(
        key="qg_kerr_quantum", entry=ENTRY, title="The Kerr area bound at the area of N spin-1/2 punctures (a heuristic composition)",
        source_kind="derived", source="derived here; composes the Kerr bound of the previous display with the area spectrum and gamma = ln 2/(pi sqrt 3) of qg_spin_network.py",
        latex_str=(r"A_N=4\sqrt3\pi\gamma_I N\ell_P^{2}\ \ge\ 8\pi|J|\ (\hbar=1)\ \Longleftrightarrow\ N\ \ge\ \frac{2|J|}{\sqrt3\gamma_I}\xrightarrow{\gamma_I=\ln2/(\pi\sqrt3)}\frac{2\pi|J|}{\ln2}"
                   r"%%BR%%|J|=\tfrac12:\ N\ge\frac{\pi}{\ln2}=4.53\ (N\ge5),\qquad |J|=1:\ N\ge\frac{2\pi}{\ln2}=9.06\ (N\ge10)"),
        notes="HEURISTIC, with the caveats of qg_penrose: the Kerr bound A >= 8 pi |J| is a statement about a classical Kerr solution; an LQG area eigenvalue is an "
              "eigenvalue of an operator on spin networks; the value of gamma is that of counting spin-1/2 punctures only; and J is taken in units of hbar with the ordinary spin "
              "quantisation of angular momentum, which is an assumption. The closed forms are exact consequences. The number 2 pi/ln 2 for |J| = hbar is the same number as the "
              "threshold at which the Schwarzschild radius of the Penrose-saturating mass exceeds its reduced Compton wavelength (qg_penrose): both are the equation "
              "4 N ln 2 = 8 pi at different readings of the left-hand side, so the coincidence is arithmetic.",
    ))
    Nmin = lambda j: float(2 * sp.pi * j / sp.log(2))
    check(NAME, "with gamma = ln2/(pi sqrt3), 4 sqrt3 pi gamma N = 4 N ln2 and the bound 4 N ln2 >= 8 pi |J| reads N >= 2 pi |J|/ln2 (exact)",
          sp.simplify(4 * sp.sqrt(3) * sp.pi * gamma_val * N - 4 * N * sp.log(2)) == 0 and sp.simplify(sp.solve(sp.Eq(4 * N * sp.log(2), 8 * sp.pi * Jv), N)[0] - 2 * sp.pi * Jv / sp.log(2)) == 0, eq=eq_q)
    check(NAME, f"smallest spin J = 1/2: N >= pi/ln2 = {Nmin(0.5):.3f}, so five spin-1/2 punctures; J = 1: N >= {Nmin(1):.3f}, so ten; the area of five punctures is {float(4*5*sp.log(2)):.3f} l_P^2 against the Kerr minimum 4 pi = {float(4*sp.pi):.3f} l_P^2",
          4.53 < Nmin(0.5) < 4.54 and 9.06 < Nmin(1) < 9.07 and 4 * 5 * float(sp.log(2)) > float(4 * sp.pi) > 4 * 4 * float(sp.log(2)), eq=eq_q)
    check(NAME, "the J = hbar threshold 2 pi/ln2 is the same number as the Compton threshold of qg_penrose (N >= 2 pi/ln2 = 9.06, first met at N = 10)",
          sp.simplify(2 * sp.pi * 1 / sp.log(2) - 2 * sp.pi / sp.log(2)) == 0 and 9.06 < float(2 * sp.pi / sp.log(2)) < 9.07, eq=eq_q)
    nn = float(4 * sp.pi ** 2 * sp.Rational(1, 4) / sp.log(2) ** 2)
    check(NAME, f"the area product A_+ A_- = (8 pi J)^2 at spin-1/2 punctures with gamma = ln2/(pi sqrt3): if A_+ = 4 N ln2 and A_- = 4 N' ln2 (l_P = 1) then N N' = 4 pi^2 J^2/ln^2 2 = {nn:.2f} for J = 1/2, "
                f"not an integer (nearest {round(nn)}); no claim for other spins or for a general puncture set, and the inner horizon is not a surface of the spin-network states in any case",
          sp.simplify(Ap * Am / (16 * sp.log(2) ** 2) - 4 * sp.pi ** 2 * J ** 2 / sp.log(2) ** 2) == 0
          and abs(nn - 20.54) < 0.01 and abs(nn - round(nn)) > 0.4, eq=eq_q)
    # mass: the Kerr subextremal mass bound against the Penrose bound for the same punctures, J = 1/2, N = 5
    n5 = 5
    penrose_N = float(sp.sqrt(n5 * sp.log(2) / (4 * sp.pi)))
    kerr_J = float(sp.sqrt(sp.Rational(1, 2)))
    check(NAME, f"for J = 1/2 and N = 5: the Penrose mass bound is M/M_P >= {penrose_N:.4f} and the subextremal Kerr bound M >= sqrt|J| M_P = {kerr_J:.4f}: at the smallest spin the Kerr bound is the larger",
          penrose_N < kerr_J and abs(penrose_N - 0.5253) < 5e-4 and abs(kerr_J - 0.7071) < 5e-4, eq=eq_q)
