"""The Penrose inequality (entry 260 of the OpenAI release) set beside the LQG area spectrum (companion to qg_spin_network.py).

Entry 260 ("Spacetime Penrose inequalities: enclosing area, charge, rotation, and anti-de Sitter extensions"; the paper of
5 Oct 2026, preprints/Spacetime-Penrose-inequalities-enclosing-area-charge-and-rigidity-October-5-2026) proves, for asymptotically
flat initial data satisfying the dominant energy condition, that the invariant ADM mass is bounded below by the minimum enclosing
area of a weakly trapped boundary,

    m >= sqrt(A_min / 16 pi)   (three spatial dimensions, G = c = 1),

with the charged, four-dimensional and anti-de Sitter variants quoted in its section 1.1. Loop quantum gravity gives the area of a
surface a discrete spectrum, with a smallest nonzero value 4 sqrt3 pi gamma l_P^2 (qg_spin_network.py). This module asks what the
bound says if the two are put side by side, and checks every piece of arithmetic that statement needs.

Registered, with their checks:

  qg_penrose_equality  the saturation of the inequality by Schwarzschild, from the library's Schwarzschild radius and horizon area
                       equations; the bound in Planck units, M / M_P >= sqrt(A / (16 pi l_P^2)); and the equality cases of the charged
                       (Reissner-Nordstrom), four-dimensional (Schwarzschild-Tangherlini) and anti-de Sitter (Schwarzschild-AdS) forms,
                       which the paper's section 1.1 states: each is checked as an algebraic identity, not as a theorem.
  qg_penrose_quantum   the bound at a quantised area: for N spin-1/2 punctures with gamma = ln2/(pi sqrt3), A = 4 N ln2 l_P^2 and
                       M / M_P >= sqrt(N ln2 / 4 pi), which is sqrt(S / (4 pi k_B)) in terms of the Bekenstein-Hawking entropy; the
                       smallest value is 0.2349 M_P; and the crude criterion that the Schwarzschild radius of the saturating mass exceeds its
                       reduced Compton wavelength holds first at N = 10.

Checked: the identities above, exactly (sympy), and the numbers quoted.
Not tested: the Penrose inequality itself (a theorem of the paper, not re-proved), anything about data that are not saturated, the charged
bound beyond its algebraic form and the branch r >= Q, any claim that a quantum area eigenvalue is the "minimum enclosing area" of
classical initial data. That last identification is the heuristic of this module: the theorem is a statement about classical Riemannian
initial data, an area eigenvalue is an eigenvalue of an operator on a Hilbert space of spin networks, and nothing here connects the two
beyond putting the same symbol A in both. For N = 1 the area radius is 0.47 l_P, below the Planck length, where a classical horizon is not a
meaningful object: the Compton criterion above is a crude way of saying so, and fails there.
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import knowledgebase as kb
from eqlib import Equation, check, register

NAME = "qg_penrose"
ENTRY = "260"
PAPER = "OpenAI, Spacetime Penrose inequalities (preprints/Spacetime-Penrose-inequalities-enclosing-area-charge-and-rigidity-October-5-2026), section 1.1"


def build_and_verify():
    G, c, hbar, A, M, r, ell, Q, m, N, gam = sp.symbols("G c hbar A M r ell_P Q m N gamma_I", positive=True)
    Mp = sp.sqrt(hbar * c / G)
    lp = sp.sqrt(hbar * G / c ** 3)

    # ---------------------------------------------------------------------------------------------
    # saturation and the variants
    # ---------------------------------------------------------------------------------------------
    eq_eq = register(NAME, Equation(
        key="qg_penrose_equality", entry=ENTRY, title="The Penrose bound: saturation by Schwarzschild, Planck units and the equality cases of its variants",
        source_kind="derived", source=f"derived here from the statements of {PAPER}; the inequalities themselves are the paper's theorems",
        latex_str=(r"m\ \ge\ \sqrt{\frac{A_{\min}}{16\pi}}\ \ (G=c=1),\qquad M\ \ge\ \frac{c^{2}}{G}\sqrt{\frac{A_{\min}}{16\pi}},\qquad "
                   r"\frac{M}{M_P}\ \ge\ \sqrt{\frac{A_{\min}}{16\pi\,\ell_P^{2}}}"
                   r"%%BR%%A=4\pi r_s^{2},\ r_s=\tfrac{2GM}{c^{2}}\ \Longrightarrow\ M=\tfrac{c^{2}}{G}\sqrt{\tfrac{A}{16\pi}}\ \ \text{(Schwarzschild saturates)}"
                   r"%%BR%%\text{charged:}\ \ m\ge\tfrac12\bigl(r+\tfrac{Q^{2}}{r}\bigr)\ \ (r\ge Q),\quad r_+=m+\sqrt{m^{2}-Q^{2}}\ \text{saturates};\qquad "
                   r"4\text{D:}\ \ m\ge\tfrac12\bigl(\tfrac{A_e}{2\pi^{2}}\bigr)^{2/3},\ \ f=1-\tfrac{2m}{r^{2}}\ \text{saturates}"
                   r"%%BR%%\text{AdS}\ (\Lambda=-3):\ \ m\ge\sqrt{\tfrac{A}{16\pi}}\bigl(1+\tfrac{A}{4\pi}\bigr),\ \ f=1+r^{2}-\tfrac{2m}{r}\ \text{saturates}"),
        notes="The three-dimensional inequality bounds the invariant ADM mass below by the minimum enclosing area of the weakly trapped boundary; the "
              "charged, four-dimensional and anti-de Sitter forms are those quoted in the introduction of the paper (with the cosmological constant -3, "
              "the AdS inequality is a local statement near Schwarzschild-AdS data, not a global one). Each is checked at the model solution that "
              "the paper identifies as equality: Schwarzschild, Reissner-Nordstrom with the larger horizon radius, Schwarzschild-Tangherlini, "
              "Schwarzschild-AdS. That is an algebraic consistency check of the displayed formulas, not a test of the theorems. The charged "
              "bound in the form m >= (r + Q^2/r)/2 is equivalent to r <= m + sqrt(m^2 - Q^2) only for r >= Q, as the paper says; for r < Q "
              "the right-hand side is not a bound.",
    ))
    # (a) Schwarzschild saturates, from the library's own equations
    sch = kb.solve_for(kb.equation("Schwarzschild radius"), "length").expr           # 2 G M / c^2
    harea = kb.solve_for(kb.equation("Horizon area"), "area").expr                    # 4 pi r_s^2
    pen = kb.solve_for(kb.equation("Penrose mass at equality"), "mass").expr          # c^2/G sqrt(A/16 pi)
    A_sch = harea.subs(sp.Symbol("r_s"), sch.subs({sp.Symbol("G"): G, sp.Symbol("M"): M, sp.Symbol("c"): c}))
    pen_s = pen.subs({sp.Symbol("c"): c, sp.Symbol("G"): G, sp.Symbol("A"): A_sch})
    check(NAME, "the library's Schwarzschild radius r_s = 2GM/c^2 and horizon area A = 4 pi r_s^2, put into the Penrose mass (c^2/G) sqrt(A/16 pi), return M "
                "exactly: Schwarzschild data saturate the inequality", sp.simplify(pen_s - M) == 0, eq=eq_eq)
    # (b) Planck units
    ratio = sp.simplify((c ** 2 / G * sp.sqrt(A / (16 * sp.pi))) / Mp - sp.sqrt(A / (16 * sp.pi * lp ** 2)))
    check(NAME, "in Planck units M_P = sqrt(hbar c/G), l_P = sqrt(hbar G/c^3) the bound reads M/M_P >= sqrt(A/(16 pi l_P^2)), an identity in G, c, hbar, A",
          ratio == 0, eq=eq_eq)
    # (c) charged: equivalence on r >= Q, equality at r_+
    rp = m + sp.sqrt(m ** 2 - Q ** 2)
    check(NAME, "Reissner-Nordstrom: with r_+ = m + sqrt(m^2 - Q^2) the bound m >= (r + Q^2/r)/2 is an equality at r = r_+ (exact), and r_+ >= Q",
          sp.simplify((rp + Q ** 2 / rp) / 2 - m) == 0 and sp.simplify(rp ** 2 - Q ** 2 - (rp - Q) * (rp + Q)) == 0, eq=eq_eq)
    rng = np.random.default_rng(260)
    bad_equiv, bad_branch = 0, 0
    for _ in range(4000):
        mv = rng.uniform(0.5, 3.0)
        qv = rng.uniform(0.0, 0.999) * mv
        rv = rng.uniform(0.05, 6.0)
        lhs = mv >= 0.5 * (rv + qv ** 2 / rv)                    # the mass form
        rhs = rv <= mv + np.sqrt(mv ** 2 - qv ** 2)              # the area form
        if rv >= qv:
            bad_equiv += lhs != rhs
        elif lhs != rhs:
            bad_branch += 1
    check(NAME, f"the mass form and the area form of the charged bound agree at every one of the sampled points with r >= Q (seed 260; "
                f"{bad_equiv} disagreements) and disagree for some r < Q ({bad_branch} sampled points), which is why the paper restricts the mass "
                f"form to r >= Q", bad_equiv == 0 and bad_branch > 0, eq=eq_eq)
    # (d) four dimensions: Tangherlini, omega_3 = 2 pi^2
    rh, mm = sp.symbols("r_h m_t", positive=True)
    omega3 = 2 * sp.pi ** 2
    m_h = sp.solve(sp.Eq(1 - 2 * mm / rh ** 2, 0), mm)[0]                          # horizon where f = 1 - 2m/r^2 vanishes
    Ae = omega3 * rh ** 3
    check(NAME, "Schwarzschild-Tangherlini f = 1 - 2m/r^2: the horizon r_h^2 = 2m has enclosing volume A_e = 2 pi^2 r_h^3 and m = (1/2)(A_e/2 pi^2)^(2/3) "
                "exactly, so the four-dimensional bound is saturated", sp.simplify(m_h - sp.Rational(1, 2) * (Ae / omega3) ** sp.Rational(2, 3)) == 0, eq=eq_eq)
    # (e) AdS
    m_ads = sp.solve(sp.Eq(1 + rh ** 2 - 2 * mm / rh, 0), mm)[0]
    A_ads = 4 * sp.pi * rh ** 2
    check(NAME, "Schwarzschild-AdS (Lambda = -3) f = 1 + r^2 - 2m/r: the horizon has m = (r_h + r_h^3)/2 = sqrt(A/16 pi)(1 + A/4 pi) with A = 4 pi r_h^2 "
                "(exact), so the AdS form is saturated; at small area it reduces to the asymptotically flat bound",
          sp.simplify(m_ads - sp.sqrt(A_ads / (16 * sp.pi)) * (1 + A_ads / (4 * sp.pi))) == 0
          and sp.limit(sp.sqrt(A / (16 * sp.pi)) * (1 + A / (4 * sp.pi)) / sp.sqrt(A / (16 * sp.pi)), A, 0) == 1, eq=eq_eq)

    # ---------------------------------------------------------------------------------------------
    # the bound at a quantised area
    # ---------------------------------------------------------------------------------------------
    gamma_val = sp.log(2) / (sp.pi * sp.sqrt(3))
    eq_q = register(NAME, Equation(
        key="qg_penrose_quantum", entry=ENTRY, title="The Penrose bound at the area of N spin-1/2 punctures (a heuristic composition)",
        source_kind="derived", source="derived here; composes the equality-form Penrose mass of entry 260 with the LQG area spectrum and gamma = ln 2/(pi sqrt 3) of qg_spin_network.py",
        latex_str=(r"A_N=4\sqrt3\pi\gamma_I N\ell_P^{2}\Big|_{\gamma_I=\ln2/(\pi\sqrt3)}=4N\ln2\,\ell_P^{2},\qquad "
                   r"\frac{M}{M_P}\ \ge\ \sqrt{\frac{N\ln2}{4\pi}}=\sqrt{\frac{S}{4\pi k_B}},\qquad \frac{M_{\min}}{M_P}\Big|_{N=1}=\sqrt{\frac{\ln2}{4\pi}}\approx0.2349"
                   r"%%BR%%\text{general punctures:}\ \ \frac{M}{M_P}\ge\sqrt{\tfrac{\gamma_I}{2}\sum_i\sqrt{j_i(j_i+1)}},\qquad "
                   r"r_s\ge\frac{\hbar}{Mc}\ \Longleftrightarrow\ \frac{M}{M_P}\ge\frac{1}{\sqrt2}\ \Longleftrightarrow\ N\ge\frac{2\pi}{\ln2}=9.06"),
        notes="HEURISTIC. The Penrose inequality is a theorem about classical initial data; an area eigenvalue of loop quantum gravity is an eigenvalue "
              "of an operator on a Hilbert space of spin networks. Putting the eigenvalue in the place of A_min is a heuristic with no derivation here, "
              "and the value of gamma used is the one from counting spin-1/2 punctures only (a different count gives a different gamma). The "
              "closed forms are exact consequences. For N = 1 the area radius sqrt(A/4pi) is 0.47 l_P, smaller than the Planck length, and the reduced "
              "Compton wavelength of the saturating mass is 4.3 l_P, larger than the horizon; the classical statement has no content there. The "
              "criterion r_s >= hbar/(Mc) is a crude marker, not a theorem: it first holds at N = 10.",
    ))
    A_N = 4 * sp.sqrt(3) * sp.pi * gam * N * lp ** 2
    ratio_N = sp.simplify(sp.sqrt(A_N.subs(gam, gamma_val) / (16 * sp.pi * lp ** 2)) - sp.sqrt(N * sp.log(2) / (4 * sp.pi)))
    check(NAME, "with gamma = ln2/(pi sqrt3) the area of N spin-1/2 punctures is 4 N ln2 l_P^2 and the Penrose bound in Planck units is sqrt(N ln2/(4 pi)) "
                "(exact)", sp.simplify(A_N.subs(gam, gamma_val) - 4 * N * sp.log(2) * lp ** 2) == 0 and ratio_N == 0, eq=eq_q)
    bh = kb.equation("Bekenstein-Hawking entropy").relation()
    S_over_k = sp.simplify((bh.rhs / sp.Symbol("k_B")).subs({sp.Symbol("A"): A, sp.Symbol("G"): G, sp.Symbol("c"): c, sp.Symbol("hbar"): hbar}))
    lhs = sp.sqrt(A / (16 * sp.pi * lp ** 2))
    rhs = sp.sqrt(S_over_k / (4 * sp.pi))
    check(NAME, "with the library's Bekenstein-Hawking entropy S/k_B = A/(4 l_P^2) the same bound reads M/M_P >= sqrt(S/(4 pi k_B)) for every A (exact)",
          sp.simplify(lhs - rhs) == 0, eq=eq_q)
    j = sp.symbols("j", positive=True)
    gen = sp.simplify(sp.sqrt(8 * sp.pi * gam * sp.sqrt(j * (j + 1)) * lp ** 2 / (16 * sp.pi * lp ** 2)) - sp.sqrt(gam / 2 * sp.sqrt(j * (j + 1))))
    check(NAME, "for punctures of spins j_i the bound is M/M_P >= sqrt((gamma/2) sum_i sqrt(j_i(j_i+1))) (exact; shown for one puncture of spin j)", gen == 0, eq=eq_q)
    mmin = float(sp.sqrt(sp.log(2) / (4 * sp.pi)))
    check(NAME, f"M_min/M_P at N = 1 is sqrt(ln2/(4 pi)) = {mmin:.4f} (0.2349 to four digits); its area radius sqrt(A/(4 pi)) is {float(sp.sqrt(sp.log(2)/sp.pi)):.4f} l_P "
                f"and its reduced Compton wavelength {1/mmin:.2f} l_P",
          abs(mmin - 0.2349) < 5e-5 and abs(float(sp.sqrt(sp.log(2) / sp.pi)) - 0.4697) < 5e-5 and abs(1 / mmin - 4.257) < 5e-3, eq=eq_q)
    # Compton criterion: r_s = sqrt(A/4 pi) at saturation, hbar/(M c) = l_P/(M/M_P)
    fN = lambda n: float(sp.sqrt(n * sp.log(2) / (4 * sp.pi)))
    r_s_over_lp = lambda n: 2 * fN(n)                                              # r_s = 2 G M/c^2 = 2 (M/M_P) l_P
    compton_over_lp = lambda n: 1 / fN(n)
    first = next(n for n in range(1, 60) if r_s_over_lp(n) >= compton_over_lp(n))
    ok_cross = all(abs(r_s_over_lp(n) - float(sp.sqrt(4 * n * sp.log(2) / (4 * sp.pi)))) < 1e-12 for n in range(1, 30))   # r_s = sqrt(A/4 pi) = l_P sqrt(N ln2/pi)
    check(NAME, f"at saturation r_s = sqrt(A/(4 pi)) = l_P sqrt(N ln2/pi) (checked for N = 1..29), and r_s >= hbar/(Mc) first holds at N = {first} "
                f"(N = 9 gives {r_s_over_lp(9):.4f} against {compton_over_lp(9):.4f} l_P, N = 10 gives {r_s_over_lp(10):.4f} against {compton_over_lp(10):.4f}); "
                f"the threshold 2 pi/ln2 = {float(2*sp.pi/sp.log(2)):.3f}", first == 10 and ok_cross and 9.06 < float(2 * sp.pi / sp.log(2)) < 9.07, eq=eq_q)
