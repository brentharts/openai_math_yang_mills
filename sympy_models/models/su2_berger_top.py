"""The SU(2) link with a squashed Hopf fibre: the Berger-sphere Laplacian is a symmetric top.

Reads SU(2) = S^3 with the left-invariant metric g = (1/4)[R1^2 (dpsi + cos(th) dphi)^2 + R2^2 (dth^2 + sin^2(th) dphi^2)]
(Euler angles; the fibre is the psi-orbit, R1 = R2 is the unit round S^3) and asks for the spectrum of its Laplace-Beltrami
operator in the same Wigner functions that give the Casimir j(j+1) of models/su2_lattice_gap.py and models/qg_spin_network.py.

What is checked:
  * The Laplacian of g, computed from the metric in coordinates, acts on the Wigner function e^{-i(m phi + m' psi)} d^j_{m m'}(theta) with
    the eigenvalue 4[(j(j+1) - m'^2)/R2^2 + m'^2/R1^2]. At R1 = R2 this is 4 j(j+1) = k(k+2), k = 2j, with degeneracy (2j+1)^2.
  * That is the spectrum of a quantum symmetric top (moments I_perp ~ R2^2, I_par ~ R1^2; m' is the body-axis projection).
  * The m' = 0 sector is the Laplacian of a round S^2 of radius R2/2 (the base of the Hopf fibration), with multiplicity 2l+1.
  * Weyl's law with the volume 2 pi^2 R1 R2^2, by exact counting of eigenvalues.
  * The first nonzero eigenvalue is min(8, 2 + 1/s^2)/R2^2 for s = R1/R2; Lichnerowicz's bound (3/2) min Ric holds for s < sqrt(2), with
    equality only at the round sphere, and is far from sharp as the fibre collapses.
  * Adiabatic decoupling: for s < 1/sqrt(L - 2) every eigenvalue below L belongs to the m' = 0 sector, i.e. the S^2 spectrum.

What is not claimed:
  * The Berger metric is left-invariant but not right-invariant: of the two SU(2) actions on a link (the gauge transformations at its two
    ends) it keeps one and breaks the other to U(1). The top Hamiltonian is therefore NOT a gauge-invariant deformation of the
    Kogut-Susskind electric term j(j+1). This model computes a Riemannian spectrum. It is not a lattice gauge theory.
  * Nothing here concerns Yang-Mills existence or a mass gap; the same limits as HOPF_GAUGE.md and YM_CONNECTIONS.md apply.
  * The pancake reading (s = R1/R2 = a0 t / b0 near the pancake in Quantum_Gravity_Lab) is a statement about the spatial slice's
    Laplacian at each time, not about a field solving any equation on that background.

Sources: none retrieved. The Berger-sphere Laplace spectrum is a classical fact and the symmetric-top levels are in every quantum
mechanics text; both are rederived here from the metric, so no check rests on a citation. models/su2_lattice_gap.py's `k^2 - 1` Laplacian of
S^3 on class functions is the same round-sphere spectrum (k = 2j + 1). Written for this fork; no code is copied from another repository.
"""
from fractions import Fraction as Fr
from itertools import product

import sympy as sp
from sympy.physics.quantum.spin import Rotation

from eqlib import Equation, check, register

NAME = "su2_berger_top"
ENTRY = "hopf-fibration"
SRC = "classical (no source retrieved); rederived here from the metric"


def half_integers(jmax2):
    """j = k/2 for k = 0..jmax2 as sympy rationals."""
    return [sp.Rational(k, 2) for k in range(jmax2 + 1)]


def weights(j):
    """m = -j, -j+1, ..., j."""
    k = int(2 * j)
    return [sp.Rational(a, 2) for a in range(-k, k + 1, 2)]


def eigenvalue(j, mp, s2):
    """Eigenvalue in units of 1/R2^2 with s2 = (R1/R2)^2, exact rationals (j, mp may be Fractions)."""
    return 4 * (j * (j + 1) - mp * mp) + 4 * mp * mp / s2


def count_below(s2, lam_max, max_k=3000):
    """Number of eigenvalues <= lam_max (R2 = 1), with multiplicity 2j+1, by exact counting."""
    total = 0
    for k in range(max_k):
        j = Fr(k, 2)
        # lam = 4 j(j+1) + 4 m'^2 (1/s2 - 1) is linear in m'^2, so for every j' >= j it is at least this monotone bound
        # (the minimum over m' is at m' = 0 if s2 <= 1 and at |m'| = j' if s2 > 1):
        low = 4 * j * (j + 1) if s2 <= 1 else 4 * j + 4 * j * j / s2
        if low > lam_max:
            break
        for a in range(-k, k + 1, 2):
            if eigenvalue(j, Fr(a, 2), s2) <= lam_max:
                total += k + 1
    else:
        raise RuntimeError("max_k too small")
    return total


def count_sphere2(lam_max):
    """Eigenvalues <= lam_max of the round S^2 of radius 1/2, with multiplicity 2l+1: 4 l (l+1) <= lam_max."""
    total, l = 0, 0
    while 4 * l * (l + 1) <= lam_max:
        total += 2 * l + 1
        l += 1
    return total


def first_nonzero(s2, k_max=24, integer_j_only=False):
    """Smallest nonzero eigenvalue by brute force over j <= k_max/2 (integer j only: the quotient SO(3) = S^3/{+-1})."""
    best = None
    for k in range(k_max + 1):
        if integer_j_only and k % 2:
            continue
        j = Fr(k, 2)
        for a in range(-k, k + 1, 2):
            if k == 0:
                continue
            v = eigenvalue(j, Fr(a, 2), s2)
            if best is None or v < best:
                best = v
    return best


def build_and_verify():
    th = sp.symbols("theta", positive=True)
    ph, ps = sp.symbols("phi psi", real=True)
    R1, R2 = sp.symbols("R1 R2", positive=True)
    X = [th, ph, ps]

    # ---- the metric, its volume and the Laplacian ------------------------------------------------------------------------------------
    g = sp.Rational(1, 4) * sp.Matrix([[R2**2, 0, 0],
                                       [0, R2**2 * sp.sin(th) ** 2 + R1**2 * sp.cos(th) ** 2, R1**2 * sp.cos(th)],
                                       [0, R1**2 * sp.cos(th), R1**2]])
    dens = R1 * R2**2 * sp.sin(th) / 8                         # sqrt(det g), to be checked
    ginv = sp.simplify(g.inv())

    def laplacian(f):
        return sum(sp.diff(dens * ginv[i, j] * sp.diff(f, X[j]), X[i]) for i in range(3) for j in range(3)) / dens

    eq_lap = register(NAME, Equation(
        key="bt_spectrum", entry=ENTRY, title="Laplacian spectrum of the Berger sphere (symmetric top)",
        source_kind="derived", source=SRC,
        latex_str=(r"-\Delta\,e^{-i(m\varphi+m'\psi)}d^{\,j}_{mm'}(\theta)=4\Big[\frac{j(j+1)-m'^2}{R_2^2}+\frac{m'^2}{R_1^2}\Big]"
                   r"e^{-i(m\varphi+m'\psi)}d^{\,j}_{mm'}(\theta)"),
        notes="Classical, rederived. Metric (1/4)[R1^2 (dpsi + cos th dphi)^2 + R2^2 (dth^2 + sin^2 th dphi^2)]; psi has period 4 pi (the standard S^3). "
              "Degeneracy 2j+1 for each (j, m') in the Peter-Weyl sum. Left-invariant, not right-invariant: this is not a gauge-invariant "
              "electric term. Silent on Yang-Mills existence and any gap.",
    ))
    check(NAME, "sqrt(det g) = R1 R2^2 sin(theta)/8", sp.simplify(g.det() - dens**2) == 0, eq=eq_lap)
    vol = sp.integrate(dens, (th, 0, sp.pi), (ph, 0, 2 * sp.pi), (ps, 0, 4 * sp.pi))
    check(NAME, "volume = 2 pi^2 R1 R2^2 (the unit S^3 has 2 pi^2)", sp.simplify(vol - 2 * sp.pi**2 * R1 * R2**2) == 0, eq=eq_lap)

    ok_exact, n_exact, ok_num, n_num = True, 0, True, 0
    points = [{th: sp.Rational(7, 10), ph: sp.Rational(3, 10), ps: sp.Rational(1, 5), R1: sp.Rational(4, 5), R2: sp.Rational(13, 10)},
              {th: sp.Rational(23, 10), ph: sp.Rational(-11, 10), ps: sp.Rational(37, 10), R1: sp.Rational(3, 2), R2: sp.Rational(1, 2)},
              {th: sp.Rational(1, 3), ph: sp.Rational(5, 2), ps: sp.Rational(-9, 10), R1: sp.Rational(1, 20), R2: sp.Rational(3, 1)}]
    for j in (sp.Rational(1, 2), sp.Integer(1), sp.Rational(3, 2)):
        for m, mp in product(weights(j), weights(j)):
            f = sp.exp(-sp.I * (m * ph + mp * ps)) * Rotation.d(j, m, mp, th).doit()
            lam = 4 * ((j * (j + 1) - mp**2) / R2**2 + mp**2 / R1**2)
            resid = -laplacian(f) - lam * f
            if j <= 1:                                             # exact: sympy reduces these to 0
                ok_exact &= sp.simplify(resid) == 0
                n_exact += 1
            else:                                                  # j = 3/2: simplify does not reach 0 reliably; 40 digits at 3 points
                ok_num &= all(abs(sp.N(resid.subs(pt), 40)) < sp.Float("1e-30") * (1 + abs(sp.N(f.subs(pt), 40))) for pt in points)
                n_num += 1
    check(NAME, f"-Laplacian of the Wigner function = lambda * function, symbolically, for all {n_exact} triples (j, m, m') with j = 1/2, 1",
          ok_exact, eq=eq_lap)
    check(NAME, f"the same for all {n_num} triples with j = 3/2, to 40 digits at 3 points (R1, R2, theta, phi, psi) (sympy does not "
                "simplify these to 0)", ok_num, eq=eq_lap)
    j_, m_ = sp.symbols("j m_", positive=True)
    round_ev = sp.simplify(4 * ((j_ * (j_ + 1) - m_**2) + m_**2) - 4 * j_ * (j_ + 1))
    check(NAME, "round limit R1 = R2 = 1: lambda = 4 j(j+1) = k(k+2) with k = 2j, independent of m'", round_ev == 0, eq=eq_lap)
    check(NAME, "round limit: degeneracy sum_{m'} (2j+1) = (2j+1)^2 = (k+1)^2",
          all(sum(2 * k_j + 1 for _ in weights(k_j)) == (2 * k_j + 1) ** 2 for k_j in half_integers(6)), eq=eq_lap)

    # ---- the symmetric top --------------------------------------------------------------------------------------------------------------
    Ip, Ia = sp.symbols("I_perp I_par", positive=True)
    J, K = sp.symbols("J K", positive=True)
    top = (J * (J + 1) - K**2) / (2 * Ip) + K**2 / (2 * Ia)
    sub = {Ip: R2**2 / 8, Ia: R1**2 / 8}
    eq_top = register(NAME, Equation(
        key="bt_top", entry=ENTRY, title="The Berger-sphere Laplacian is the quantum symmetric top",
        source_kind="derived", source="derived here from bt_spectrum; the top's levels are textbook",
        latex_str=(r"E_{JK}=\frac{J(J+1)-K^2}{2I_\perp}+\frac{K^2}{2I_\parallel}=\lambda_{j,m'}\ \ (J=j,\ K=m'),\qquad "
                   r"I_\perp=\frac{R_2^2}{8},\ \ I_\parallel=\frac{R_1^2}{8}\ \ (\hbar=1)"),
        notes="With hbar = 1 and the identification J = j, K = m', the top energy E = (J(J+1) - K^2)/(2 I_perp) + K^2/(2 I_par) equals the "
              "Laplacian eigenvalue exactly for I_perp = R2^2/8 and I_par = R1^2/8. J half-integer is allowed because SU(2) double-covers the "
              "top's rotation group. An identification of spectra, not of Hamiltonians' symmetry: the top is invariant under the rotation "
              "group it lives on and a body-axis U(1), the link's Kogut-Susskind term under two SU(2)s.",
    ))
    lam_sym = 4 * ((J * (J + 1) - K**2) / R2**2 + K**2 / R1**2)
    check(NAME, "E_{JK} with I_perp = R2^2/8, I_par = R1^2/8 equals the Laplacian eigenvalue 4[(J(J+1)-K^2)/R2^2 + K^2/R1^2]",
          sp.simplify(top.subs(sub) - lam_sym) == 0, eq=eq_top)

    # ---- the m' = 0 sector is the base S^2 ----------------------------------------------------------------------------------------------
    l = sp.symbols("l", nonnegative=True, integer=True)
    eq_base = register(NAME, Equation(
        key="bt_base", entry=ENTRY, title="The fibre-constant sector is the Laplacian of the base S^2 of radius R2/2",
        source_kind="derived", source="derived here from bt_spectrum",
        latex_str=(r"m'=0\ (j=l\in\mathbb Z):\quad \lambda=\frac{4\,l(l+1)}{R_2^2}=\frac{l(l+1)}{(R_2/2)^2},\quad \text{multiplicity }2l+1,"
                   r"\ \text{for every }R_1"),
        notes="The functions constant along the Hopf fibre are functions on S^2. The base radius R2/2 is the one used by Quantum_Gravity_Lab's "
              "horizon code. A half-integer j has only half-integer m', so m' = 0 forces j to be an integer.",
    ))
    check(NAME, "m' = 0 eigenvalue 4 l(l+1)/R2^2 equals l(l+1)/(R2/2)^2, independent of R1",
          sp.simplify(lam_sym.subs({J: l, K: 0}) - l * (l + 1) / (R2 / 2) ** 2) == 0, eq=eq_base)
    check(NAME, "m' = 0 occurs only for integer j (half-integer j has half-integer weights only)",
          all((0 in weights(jj)) == (jj.q == 1) for jj in half_integers(12)), eq=eq_base)
    # the base S^2 of radius R2/2 directly: metric (R2^2/4)(dtheta^2 + sin^2 theta dphi^2), Laplacian of P_l(cos theta)
    u = sp.cos(th)
    lap_s2 = lambda f: (sp.diff(sp.sin(th) * sp.diff(f, th), th) / sp.sin(th)) / (R2**2 / 4)
    check(NAME, "the round S^2 of radius R2/2 has -Laplacian P_l(cos theta) = l(l+1)/(R2/2)^2 P_l, l = 0..4, with the m' = 0 value of the Berger sphere",
          all(sp.simplify(-lap_s2(sp.legendre(l_, u)) - 4 * l_ * (l_ + 1) / R2**2 * sp.legendre(l_, u)) == 0 for l_ in range(5)), eq=eq_base)

    # ---- Weyl's law ---------------------------------------------------------------------------------------------------------------------
    eq_weyl = register(NAME, Equation(
        key="bt_weyl", entry=ENTRY, title="Weyl's law on the Berger sphere",
        source_kind="derived", source="derived here: Weyl (1911) with the volume of bt_spectrum, tested by exact counting",
        latex_str=(r"N(\Lambda)\sim\frac{\mathrm{Vol}}{6\pi^2}\Lambda^{3/2}=\frac{R_1R_2^2}{3}\Lambda^{3/2},\qquad \mathrm{Vol}=2\pi^2R_1R_2^2"),
        notes="Tested at Lambda = 20000 for s = R1/R2 with s^2 = 1/16, 1/4, 1, 4, 9 (R2 = 1): the counts agree with the law to within 0.2%. A "
              "consistency check of the spectrum, the multiplicities and the volume 2 pi^2 R1 R2^2 against each other, not a proof of Weyl's law.",
    ))
    errs = []
    for s2 in (Fr(1, 16), Fr(1, 4), Fr(1), Fr(4), Fr(9)):
        n = count_below(s2, 20000)
        w = float(sp.sqrt(sp.Rational(s2.numerator, s2.denominator))) * 20000**1.5 / 3
        errs.append(abs(n / w - 1))
    check(NAME, f"N(L) / (s L^(3/2)/3) is within 1% of 1 at L = 20000 for s^2 = 1/16, 1/4, 1, 4, 9 (worst {max(errs):.4f})",
          max(errs) < 0.01, eq=eq_weyl)

    # ---- the first nonzero eigenvalue and Lichnerowicz ------------------------------------------------------------------------------------
    eq_gap = register(NAME, Equation(
        key="bt_gap", entry=ENTRY, title="The first nonzero Laplacian eigenvalue of the Berger sphere",
        source_kind="derived", source="derived here from bt_spectrum; compared with Lichnerowicz (1958)",
        latex_str=(r"\lambda_1=\frac{1}{R_2^2}\min\!\Big(8,\ 2+\frac{R_2^2}{R_1^2}\Big),\qquad \text{crossing at }R_1^2/R_2^2=\tfrac16;"
                   r"\qquad \lambda_1\ge\tfrac32\min\mathrm{Ric}\ \ (0<s<\sqrt2),\ \text{equality iff }s=1"),
        notes="For s^2 >= 1/6 the lowest mode is j = 1/2, m' = +-1/2 (a fibre-momentum mode); below, it is the base mode l = 1, which does "
              "not feel the fibre: the gap saturates at 8/R2^2 as the fibre collapses, while the Ricci curvature along the fibre goes to 0 and "
              "Lichnerowicz's bound goes to 0 with it. For s >= sqrt(2) the minimum Ricci curvature (4 - 2 s^2)/R2^2 is not positive and the "
              "bound says nothing.",
    ))
    grid = [Fr(1, 100), Fr(1, 36), Fr(1, 10), Fr(1, 7), Fr(1, 6), Fr(1, 5), Fr(1, 4), Fr(1, 2), Fr(1), Fr(2), Fr(3), Fr(10), Fr(100)]
    check(NAME, "brute-force minimum over j <= 12 equals min(8, 2 + 1/s^2) at 13 rational values of s^2 from 1/100 to 100 (R2 = 1)",
          all(first_nonzero(s2) == min(Fr(8), 2 + 1 / s2) for s2 in grid), eq=eq_gap)
    s = sp.symbols("s", positive=True)
    check(NAME, "the two candidates cross at s^2 = 1/6: 2 + 1/s^2 = 8", sp.solve(sp.Eq(2 + 1 / s**2, 8), s) == [1 / sp.sqrt(6)], eq=eq_gap)
    check(NAME, "round sphere: lambda_1 = 3 (the unit S^3, k = 1)", min(8, 2 + 1) == 3, eq=eq_gap)

    # Ricci of the biaxial metric by Milnor's formula (frame structure constants l_i = 2 R_i/(R_j R_k)) in terms of s = R1/R2 (R2 = 1)
    Rs = (s, sp.Integer(1), sp.Integer(1))
    ell = [2 * Rs[i] / (Rs[(i + 1) % 3] * Rs[(i + 2) % 3]) for i in range(3)]
    mu = [sum(ell) / 2 - e for e in ell]
    ric = [sp.simplify(2 * mu[(i + 1) % 3] * mu[(i + 2) % 3]) for i in range(3)]
    check(NAME, "Milnor: Ricci of the Berger sphere (fibre, base, base) = (2 s^2, 4 - 2 s^2, 4 - 2 s^2)",
          [sp.simplify(a - b) for a, b in zip(ric, (2 * s**2, 4 - 2 * s**2, 4 - 2 * s**2))] == [0, 0, 0], eq=eq_gap)
    # Lichnerowicz: 3/2 min Ric <= lambda_1, as two polynomial inequalities with explicit factorisations
    below = sp.factor(2 + 1 / s**2 - 3 * s**2)                    # s <= 1: min Ric = 2 s^2
    above = sp.factor(2 + 1 / s**2 - (6 - 3 * s**2))              # 1 <= s < sqrt 2: min Ric = 4 - 2 s^2
    check(NAME, "s <= 1: 2 + 1/s^2 - 3 s^2 = (1 - s^2)(1 + 3 s^2)/s^2 >= 0, zero only at s = 1",
          sp.simplify(below - (1 - s**2) * (1 + 3 * s**2) / s**2) == 0, eq=eq_gap)
    check(NAME, "1 <= s: 2 + 1/s^2 - (6 - 3 s^2) = (3 s^2 - 1)(s^2 - 1)/s^2 >= 0, zero only at s = 1",
          sp.simplify(above - (3 * s**2 - 1) * (s**2 - 1) / s**2) == 0, eq=eq_gap)
    check(NAME, "for s <= 1 the other candidate also clears the bound: 8 - 3 s^2 = 5 + 3 (1 - s^2) > 0",
          sp.simplify(8 - 3 * s**2 - (5 + 3 * (1 - s**2))) == 0, eq=eq_gap)
    check(NAME, "for s >= 1 the minimum is the first candidate: 2 + 1/s^2 <= 8 exactly when s^2 >= 1/6 (so for all s >= 1)",
          sp.simplify((8 - (2 + 1 / s**2)) - (6 * s**2 - 1) / s**2) == 0, eq=eq_gap)
    check(NAME, "as s -> 0 the bound 3 s^2 -> 0 while lambda_1 stays at 8: Lichnerowicz is far from sharp when the fibre collapses",
          sp.limit(3 * s**2, s, 0) == 0 and first_nonzero(Fr(1, 10**6)) == 8, eq=eq_gap)

    # ---- adiabatic decoupling -------------------------------------------------------------------------------------------------------------
    eq_dec = register(NAME, Equation(
        key="bt_decoupling", entry=ENTRY, title="As the fibre collapses the low spectrum is that of the base S^2",
        source_kind="derived", source="derived here from bt_spectrum",
        latex_str=(r"s^2<\frac{1}{L-2}\ \Rightarrow\ \text{every eigenvalue}\le L/R_2^2\text{ has }m'=0,\ \text{so } N_s(L)=N_{S^2(R_2/2)}(L)"),
        notes="Any mode with fibre weight m' != 0 has eigenvalue at least 2 + 1/s^2 (R2 = 1). So below any fixed level the fibre-momentum modes "
              "leave and the S^2 Laplacian remains. A statement about the spatial Laplacian at one fixed s; with s = a0 t / b0 near the "
              "Quantum_Gravity_Lab pancake it says which modes stay light at early times. No dynamics is solved.",
    ))
    ok = True
    for L in (10, 50, 100, 1000):
        s2 = Fr(1, L - 2) * Fr(9, 10)                              # safely below the threshold
        ok &= count_below(s2, L) == count_sphere2(L)
    check(NAME, "exact counts at L = 10, 50, 100, 1000 and s^2 = 0.9/(L - 2): N_s(L) equals the S^2 count",
          ok, eq=eq_dec)
    L = sp.Symbol("L", positive=True)
    check(NAME, "the lowest fibre-momentum eigenvalue 2 + 1/s^2 (j = 1/2, m' = 1/2) equals L at s = 1/sqrt(L - 2), and decreases in s "
                "(derivative -2/s^3 < 0), so it exceeds L for every smaller s",
          sp.simplify((2 + 1 / s**2).subs(s, 1 / sp.sqrt(L - 2)) - L) == 0 and sp.simplify(sp.diff(2 + 1 / s**2, s) + 2 / s**3) == 0,
          eq=eq_dec)


    # ---- the quotient SO(3) --------------------------------------------------------------------------------------------------------------
    eq_so3 = register(NAME, Equation(
        key="bt_so3", entry=ENTRY, title="The same on the quotient SO(3) = S^3/{+-1}: integer j only",
        source_kind="derived", source="derived here from bt_spectrum: the antipodal map acts on the j-th Wigner function by (-1)^(2j)",
        latex_str=(r"\mathrm{SO}(3):\ j\in\mathbb Z,\qquad \lambda_1=\frac{1}{R_2^2}\min\!\Big(8,\ 4+\frac{4R_2^2}{R_1^2}\Big),\qquad "
                   r"\text{round }\mathbb{RP}^3:\ \lambda_1=8"),
        notes="Quantum_Gravity_Lab's horizon code uses the standard S^3 (fibre period 4 pi) and records that whether the pancake paper's slices are S^3 "
              "or the quotient SO(3) (period 2 pi) is not checked. This shows what the question changes in the Laplace spectrum: the half-integer "
              "levels, including the lowest fibre-momentum level 2 + 1/s^2, disappear; the first fibre-momentum level becomes j = 1, m' = +-1 "
              "at 4 + 4/s^2, and the crossing with the base level 8 moves to s = 1. The antipodal map is -1 on S^3 and a rotation by 2 pi of "
              "the fibre in the Euler angles, which is why the fibre period halves.",
    ))
    check(NAME, "antipodal map: psi -> psi + 2 pi multiplies e^{-i m' psi} by e^{-2 pi i m'} = (-1)^(2 m') = -1 for half-integer m', +1 for integer",
          all(sp.simplify(sp.exp(-2 * sp.pi * sp.I * mp_) - (-1) ** int(2 * mp_)) == 0 for jj in half_integers(6) for mp_ in weights(jj)),
          eq=eq_so3)
    check(NAME, "brute-force minimum over integer j <= 12 equals min(8, 4 + 4/s^2) at the same 13 values of s^2",
          all(first_nonzero(s2, integer_j_only=True) == min(Fr(8), 4 + 4 / s2) for s2 in grid), eq=eq_so3)
    check(NAME, "on SO(3) the base level 8 and the first fibre-momentum level 4 + 4/s^2 cross at s = 1; the round RP^3 has lambda_1 = 8",
          sp.solve(sp.Eq(4 + 4 / s**2, 8), s) == [1] and first_nonzero(Fr(1), integer_j_only=True) == 8, eq=eq_so3)


if __name__ == "__main__":
    import eqlib
    build_and_verify()
    print(len(eqlib.checks()), "checks passed in", NAME)
