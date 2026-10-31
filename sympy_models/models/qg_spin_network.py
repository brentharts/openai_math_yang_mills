"""Spin networks and loop quantum gravity on top of the SU(2) ladder (companion to su2_lattice_ladder.py).

This module asks what the Hilbert space of the SU(2) ladder looks like if it is read as the kinematics of loop quantum gravity
(LQG): a spin network on a fixed graph, with a geometric operator attached to each link. It is a *kinematical* model. The
Hamiltonian of the ladder is the Kogut-Susskind Hamiltonian of an SU(2) gauge theory, not the Hamiltonian constraint of
gravity, and the graph is fixed, so nothing here is diffeomorphism invariant. What is shared is the Hilbert space
(spin networks of SU(2) holonomies) and the operator that measures a link's flux.

What is registered, with its checks:

  qg_area_spectrum      A = 8 pi gamma l_P^2 sum_i sqrt(j_i (j_i + 1)); the minimal area, the first five single-link eigenvalues
                        as exact radicals, the unequal spacing and its limit 1/2 (sympy limits), the knowledge-base equation
                        agreeing with this one.
  qg_casimir_dictionary the link Casimir is the squared area in units of 8 pi gamma l_P^2, and the ladder's electric energy is a sum
                        of seven such squares; checked against the ladder's own bookkeeping on every state up to spin 6.
  qg_intertwiners       dimension of the SU(2)-invariant subspace of a tensor product of spins: Weyl's integration formula
                        (exact integers) against numerical Haar quadrature and against the triangle-rule count for four legs;
                        trivalent vertices have a unique intertwiner exactly when the ladder's triangle rule holds; N spin-1/2
                        legs give the Catalan number C_{N/2} (all even N up to 40).
  qg_entropy_immirzi    leading-order count of N spin-1/2 punctures, S = A/4 l_P^2 solved for the Immirzi parameter
                        (gamma = ln 2 / (pi sqrt 3)), with the Catalan asymptotics checked against exact integers.
  qg_quantum_tetrahedron the operator Q = J1 . (J2 x J3) on the two-dimensional invariant space of four spin-1/2: built from explicit
                        16 x 16 matrices, Hermitian, commuting with the total spin, spectrum +-sqrt(3)/4.
  qg_ladder_geometry    the mean-square area of the surface separating the two plaquettes, <c2(J)> in the ground state of the
                        ladder: exact strong-coupling coefficient 2/(3 g^8) (first-order perturbation theory, exact rationals)
                        and weak-coupling limit g^2 <c2(J)> -> 3/sqrt(10) (virial theorem for the harmonic limit), both against
                        the numerical ground state, plus truncation convergence and monotonic decrease in g^2.

Not tested: the Hamiltonian constraint of LQG or any dynamics of gravity; diffeomorphism invariance; the volume operator beyond the
four-valent vertex of four spin-1/2 (and there only the spectrum of Q, not the prefactor that turns it into a volume); the entropy
count over all spins (Domagala-Lewandowski, quoted in the notes, not rederived); the Immirzi parameter's physical value;
three or more plaquettes; spin foams; a semiclassical or continuum limit. Nothing here bears on existence or a mass gap of
Yang-Mills on R^4, and nothing here is a model of quantum gravity: it is the kinematics of one, evaluated on one small graph.
"""
import sys
from fractions import Fraction as Fr
from functools import lru_cache
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy.sparse.linalg as sla
import sympy as sp
from scipy.integrate import quad

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import knowledgebase as kb
from eqlib import Equation, check, register
from models import su2_lattice_ladder as ladder

NAME = "qg_spin_network"
ENTRY = "spin-network-qg"


# ---------------------------------------------------------------------------------------------------------------
# spins are carried as k = 2j, so every count below is an integer, as in su2_lattice_ladder.py
# ---------------------------------------------------------------------------------------------------------------
def c2(k):
    """Casimir j(j+1) of spin j = k/2, as an exact Fraction."""
    return Fr(k * (k + 2), 4)


def char(k):
    """Character of spin k/2 as a Laurent polynomial in z = e^{i psi}: z^k + z^(k-2) + ... + z^-k, as {exponent: coefficient}."""
    return {k - 2 * i: 1 for i in range(k + 1)}


def mul(p, q):
    out = {}
    for a, x in p.items():
        for b, y in q.items():
            out[a + b] = out.get(a + b, 0) + x * y
    return out


def invariants(ks):
    """Dimension of the SU(2)-invariant subspace of V_{k1/2} x V_{k2/2} x ...: the Haar integral of the product of characters.

    Weyl's integration formula for SU(2) is (1/2) CT[f(z) (1 - z^2)(1 - z^-2)] for a class function f symmetric under z -> 1/z, and
    CT[f (1 - z^2)(1 - z^-2)] / 2 = CT[f (1 - z^2)] by the symmetry, so the dimension is the coefficient of z^0 minus that of z^-2."""
    p = {0: 1}
    for k in ks:
        p = mul(p, char(k))
    return p.get(0, 0) - p.get(-2, 0)


def invariants_quad(ks):
    """The same integral by numerical quadrature: (2/pi) int_0^pi sin^2(psi) prod_i sin((k_i + 1) psi)/sin(psi) dpsi."""
    def f(psi):
        s = np.sin(psi)
        if abs(s) < 1e-12:
            return 0.0
        return (2 / np.pi) * s * s * np.prod([np.sin((k + 1) * psi) / s for k in ks])
    return quad(f, 0, np.pi, limit=200)[0]


def n_triangle(k1, k2, k3, k4):
    """Number of intermediate spins K in the coupling (1 2)(3 4): the dimension of the four-valent intertwiner space."""
    top = min(k1 + k2, k3 + k4)
    return sum(1 for K in range(abs(k1 - k2), top + 1, 2) if ladder.allowed(k1, k2, K) and ladder.allowed(K, k3, k4))


def entropy_per_puncture(N):
    """(1/N) ln C_{N/2}: the exact entropy per spin-1/2 puncture of the SU(2)-invariant count (N even)."""
    mp.mp.dps = 30
    return float(mp.log(mp.mpf(int(sp.catalan(N // 2)))) / N)


# ---------------------------------------------------------------------------------------------------------------
# the ladder ground state and the geometry of the cut between the two plaquettes
# ---------------------------------------------------------------------------------------------------------------
@lru_cache(maxsize=None)
def geometry(g2, kmax=None):
    """(E0, <c2(J)>, <c2(j_a)>) in the ground state of the ladder at g^2 = g2: the mean-square area of the surface cutting the
    shared link, and of the surface cutting one of a plaquette's own links, in units of (8 pi gamma l_P^2)^2."""
    H, states, _ = ladder.build(kmax or ladder.kmax_for(g2), g2)
    n = H.shape[0]
    w, v = sla.eigsh(H, k=2, sigma=0.0, which="LM", v0=np.ones(n) / np.sqrt(n))      # fixed start vector: reproducible
    i = int(np.argmin(w))
    p = v[:, i] ** 2
    p = p / p.sum()
    cJ = np.array([float(c2(K)) for (_, _, K) in states])
    cA = np.array([float(c2(ka)) for (ka, _, _) in states])
    return float(w[i]), float(p @ cJ), float(p @ cA)


def strong_coefficient():
    """Exact coefficient b in <c2(J)> = b / g^8 + O(g^-16): first-order perturbation theory around the electric vacuum (0, 0, 0).

    The vacuum couples through W_a and W_b to the two strings (1/2, 0, 1/2) and (0, 1/2, 1/2); an admixture a_m = <m|H|0>/(E_0 - E_m)
    with <m|H|0> = -v/g^2 and E_m - E_0 = (g^2/2) e_el(m) puts probability a_m^2 on a state m whose shared link has Casimir c2(J_m)."""
    t = sp.symbols("t", positive=True)       # g^2
    total = 0
    for m, v in ladder.neighbours_exact((0, 0, 0)):
        e = sp.Rational(ladder.e_el(m).numerator, ladder.e_el(m).denominator)
        amp = (v / t) / ((t / 2) * e)
        total += amp ** 2 * sp.Rational(m[2] * (m[2] + 2), 4)
    return sp.simplify(total * t ** 4)


def weak_limits():
    """Exact limits of g^2 <c2(J)> and g^2 <c2(j_a)> at weak coupling, from the harmonic limit of the ladder.

    The electric energy is (g^2/2) E^T K E with K = [[4, -1], [-1, 4]] in the two loop electric fields E = (E_a, E_b), and the shared
    link carries E_a - E_b, so c2(J) -> (E_a - E_b)^2 and c2(j_a) -> E_a^2 at large spin. Each of the three components of a mode of K
    eigenvalue lambda is an oscillator of frequency omega = sqrt(lambda/2) whose kinetic energy is half of its ground-state energy
    omega/2 (virial theorem): (g^2/2) lambda <E_mode^2> = omega/4 per component, so <E_mode^2> = 3 omega / (2 g^2 lambda) for three
    components. Rotating to the modes s = (E_a + E_b)/sqrt 2 and d = (E_a - E_b)/sqrt 2 (eigenvalues 3 and 5) gives both limits."""
    K = sp.Matrix([[4, -1], [-1, 4]])
    modes = {}
    for lam, _, vecs in K.eigenvects():
        v = vecs[0] / sp.sqrt((vecs[0].T * vecs[0])[0])
        modes[lam] = (v, sp.sqrt(sp.Rational(lam, 2)))
    EE = sp.zeros(2)                               # g^2 <E_i E_j>
    for lam, (v, om) in modes.items():
        EE += (3 * om / (2 * lam)) * (v * v.T)
    cJ = sp.simplify((sp.Matrix([[1, -1]]) * EE * sp.Matrix([1, -1]))[0])
    cA = sp.simplify(EE[0, 0])
    return cJ, cA


# ---------------------------------------------------------------------------------------------------------------
# the quantum tetrahedron: four spin-1/2, explicit matrices
# ---------------------------------------------------------------------------------------------------------------
def tetrahedron():
    """Q = eps_{abc} J1_a J2_b J3_c restricted to the invariant subspace of four spin-1/2, from explicit 16 x 16 matrices.

    Returns (dimension of the invariant subspace, Q restricted to it in an orthonormal basis, Hermitian, commutes with total J)."""
    s = [sp.Matrix([[0, 1], [1, 0]]) / 2, sp.Matrix([[0, -sp.I], [sp.I, 0]]) / 2, sp.Matrix([[1, 0], [0, -1]]) / 2]
    I2 = sp.eye(2)

    def op(leg, a):
        mats = [I2] * 4
        mats[leg] = s[a]
        out = mats[0]
        for m in mats[1:]:
            out = sp.kronecker_product(out, m)
        return out

    J = [sum((op(leg, a) for leg in range(4)), sp.zeros(16)) for a in range(3)]
    J2 = sum((Ja * Ja for Ja in J), sp.zeros(16))
    Q = sp.zeros(16)
    for a in range(3):
        for b in range(3):
            for c in range(3):
                e = sp.LeviCivita(a, b, c)
                if e:
                    Q += e * op(0, a) * op(1, b) * op(2, c)
    Q = Q.applyfunc(sp.simplify)
    null = J2.nullspace()
    basis = sp.Matrix.hstack(*sp.GramSchmidt(null, True))
    Qr = (basis.H * Q * basis).applyfunc(sp.simplify)
    herm = (Q - Q.H).applyfunc(sp.simplify) == sp.zeros(16)
    commutes = all((Q * Ja - Ja * Q).applyfunc(sp.simplify) == sp.zeros(16) for Ja in J)
    return len(null), Qr, herm, commutes


def build_and_verify():
    gamma, ell, j = sp.symbols("gamma_I ell_P j")

    # =====================================================================================================
    # the area spectrum
    # =====================================================================================================
    eq_area = register(NAME, Equation(
        key="qg_area_spectrum", entry=ENTRY, title="The area spectrum of a spin-network surface",
        source_kind="textbook", source="[RS], [AL]; the closed-form checks are done here",
        latex_str=(r"A\bigl(\{j_i\}\bigr)=8\pi\gamma_{I}\,\ell_P^{2}\sum_{i}\sqrt{j_i(j_i+1)},\qquad "
                   r"A_{\min}=4\sqrt3\,\pi\gamma_{I}\ell_P^{2}\ \ (j=\tfrac12)"
                   r"%%BR%%\sqrt{j(j+1)}=\Bigl(j+\tfrac12\Bigr)-\frac{1}{8j}+O(j^{-2}),\qquad "
                   r"\sqrt{(j+\tfrac12)(j+\tfrac32)}-\sqrt{j(j+1)}\ \xrightarrow{j\to\infty}\ \tfrac12"),
        notes="Each link of spin j that pierces the surface adds 8 pi gamma l_P^2 sqrt(j(j+1)); the area is the sum over the punctures "
              "(a surface that meets no vertex). The Immirzi parameter gamma is a free parameter of the area spectrum. The spectrum is "
              "discrete with a smallest nonzero value, and is not equally spaced: the gaps shrink to 1/2 in units of 8 pi gamma l_P^2 "
              "for a step of half a unit in j. Quoted from the standard references and checked here in closed form; the derivation "
              "of the operator from the quantised triad is not done. Surfaces through vertices, and the volume, are not covered.",
    ))
    a_unit = 8 * sp.pi * gamma * ell ** 2
    check(NAME, "j = 1/2 gives the smallest nonzero area 8 pi gamma l_P^2 sqrt(3/4) = 4 sqrt 3 pi gamma l_P^2",
          sp.simplify(a_unit * sp.sqrt(sp.Rational(1, 2) * sp.Rational(3, 2)) - 4 * sp.sqrt(3) * sp.pi * gamma * ell ** 2) == 0, eq=eq_area)
    first = [sp.sqrt(sp.Rational(k, 2) * (sp.Rational(k, 2) + 1)) for k in range(1, 6)]
    check(NAME, "the first five single-link eigenvalues sqrt(j(j+1)), j = 1/2 .. 5/2, are sqrt3/2, sqrt2, sqrt15/2, sqrt6, sqrt35/2 "
                "(in units of 8 pi gamma l_P^2)",
          first == [sp.sqrt(3) / 2, sp.sqrt(2), sp.sqrt(15) / 2, sp.sqrt(6), sp.sqrt(35) / 2], eq=eq_area)
    jp = sp.symbols("j_p", positive=True)
    check(NAME, "(j + 1/2)^2 - j(j+1) = 1/4, so sqrt(j(j+1)) lies strictly below j + 1/2 and the difference is 1/(4 (j + 1/2 + sqrt(j(j+1)))) "
                "= 1/(8j) + O(j^-2)",
          sp.simplify((jp + sp.Rational(1, 2)) ** 2 - jp * (jp + 1) - sp.Rational(1, 4)) == 0
          and sp.limit(8 * jp * ((jp + sp.Rational(1, 2)) - sp.sqrt(jp * (jp + 1))), jp, sp.oo) == 1, eq=eq_area)
    gap = sp.sqrt((jp + sp.Rational(1, 2)) * (jp + sp.Rational(3, 2))) - sp.sqrt(jp * (jp + 1))
    steps = [float(gap.subs(jp, sp.Rational(k, 2))) for k in range(1, 9)]
    check(NAME, "the spectrum is not equally spaced: the steps sqrt((j+1/2)(j+3/2)) - sqrt(j(j+1)) for j = 1/2, 1, ..., 4 are " +
          ", ".join(f"{s:.4f}" for s in steps) + ", all above 1/2 and strictly decreasing",
          abs(steps[0] - 0.5482) < 1e-3 and abs(steps[1] - 0.5223) < 1e-3
          and all(s > 0.5 for s in steps) and all(a > b for a, b in zip(steps, steps[1:])), eq=eq_area)
    check(NAME, "the step for a half-unit change of j tends to exactly 1/2 as j -> infinity (sympy limit)",
          sp.limit(gap, jp, sp.oo) == sp.Rational(1, 2), eq=eq_area)
    kb_area = kb.solve_for(kb.equation("LQG area spectrum"), "area").expr
    kb_c2 = kb.solve_for(kb.equation("Spin Casimir"), "casimir eigenvalue").expr
    mine = a_unit * sp.sqrt(j * (j + 1))
    check(NAME, "the knowledge-base equations 'Spin Casimir' and 'LQG area spectrum', joined on the Casimir, give exactly this "
                "formula for one puncture",
          sp.simplify(kb_area.subs(sp.Symbol("c_2"), kb_c2) - mine) == 0 and kb.check_dimensions("LQG area spectrum") is None, eq=eq_area)

    # =====================================================================================================
    # the Casimir is the squared area
    # =====================================================================================================
    eq_dict = register(NAME, Equation(
        key="qg_casimir_dictionary", entry=ENTRY, title="The link Casimir as a squared area, and the ladder's electric energy",
        source_kind="derived", source="derived here from qg_area_spectrum and su2_ladder_basis; [KS] for the Hamiltonian",
        latex_str=(r"E_l^{2}=c_2(j_l)=\Bigl(\frac{A_l}{8\pi\gamma_{I}\ell_P^{2}}\Bigr)^{2}"
                   r"%%BR%%\sum_{l=1}^{7}E_l^{2}=3c_2(j_a)+3c_2(j_b)+c_2(J)=\sum_{l=1}^{7}\Bigl(\frac{A_l}{8\pi\gamma_{I}\ell_P^{2}}\Bigr)^{2}"),
        notes="On a link of the lattice theory the electric field squared is the Casimir j(j+1); on a link of a spin network the area "
              "through the surface it pierces is 8 pi gamma l_P^2 sqrt(j(j+1)). So the electric term of the ladder's Hamiltonian is, link by "
              "link, the square of an area in units of 8 pi gamma l_P^2. This is a relation between two operators on one Hilbert space, "
              "not a statement that the ladder's dynamics is gravity: the Hamiltonian constraint of LQG is a different operator, and "
              "no multiplication by an area squared appears in it.",
    ))
    check(NAME, "(A/(8 pi gamma l_P^2))^2 = j(j+1) when A is the area spectrum of one puncture, as a sympy identity",
          sp.simplify((mine / a_unit) ** 2 - j * (j + 1)) == 0, eq=eq_dict)
    bad = 0
    nstate = 0
    for ka in range(7):
        for kb_ in range(7):
            for K in range(abs(ka - kb_), ka + kb_ + 1, 2):
                nstate += 1
                links = [c2(ka)] * 3 + [c2(kb_)] * 3 + [c2(K)]
                bad += ladder.e_el((ka, kb_, K)) != sum(links)
    check(NAME, f"the ladder's electric Casimir sum equals the sum over its seven links (three of spin j_a, three of j_b, one of J) "
                f"for all {nstate} gauge-invariant states with j_a, j_b <= 3", bad == 0, eq=eq_dict)
    H, st, _ = ladder.build(6, 0.7)
    ref = np.array([0.5 * 0.7 * float(sum([c2(ka)] * 3 + [c2(kb_)] * 3 + [c2(K)])) + 4 / 0.7 for (ka, kb_, K) in st])
    check(NAME, "the diagonal of the ladder Hamiltonian at g^2 = 0.7 is (g^2/2) times the sum of seven squared areas (units 8 pi gamma l_P^2) "
                "plus 4/g^2, to machine precision",
          np.max(np.abs(H.diagonal() - ref)) < 1e-12, eq=eq_dict)

    # =====================================================================================================
    # intertwiners
    # =====================================================================================================
    eq_int = register(NAME, Equation(
        key="qg_intertwiners", entry=ENTRY, title="Intertwiner spaces: the dimension of the invariants of a vertex",
        source_kind="derived", source="derived here; Weyl's integration formula and the triangle rule are textbook; [Bar] for the four-valent case",
        latex_str=(r"\dim\mathrm{Inv}\Bigl(\bigotimes_{i=1}^{n}V_{j_i}\Bigr)=\int_{SU(2)}\prod_{i=1}^{n}\chi_{j_i}(U)\,dU"
                   r"=\bigl[z^{0}\bigr]\Bigl(1-z^{2}\Bigr)\prod_{i}\chi_{j_i}(z)"
                   r"%%BR%%n=3:\ \ \dim=\begin{cases}1&\triangle(j_1,j_2,j_3)\\0&\text{otherwise}\end{cases}\qquad "
                   r"n=4:\ \ \dim=\#\{J:\ \triangle(j_1,j_2,J),\ \triangle(J,j_3,j_4)\}"
                   r"%%BR%%j_i=\tfrac12,\ n=N:\ \ \dim=C_{N/2}=\frac{1}{N/2+1}\binom{N}{N/2}\ :\ 1,\ 2,\ 5,\ 14,\ 42,\ \dots"),
        notes="The invariant subspace of the tensor product of the spins meeting at a vertex is the space of intertwiners, the "
              "quantum numbers of the vertex. At a trivalent vertex it is at most one-dimensional, so a spin network on a graph with "
              "only trivalent vertices (the theta graph of the ladder) is fixed by its link spins; at a four-valent vertex it is the "
              "number of intermediate spins, and that is the first place a vertex has its own quantum numbers. The count is the Haar "
              "integral of the product of characters; here it is evaluated exactly by Weyl's formula and by numerical quadrature.",
    ))
    worst = 0.0
    cases = [[1, 1, 1, 1], [2, 2, 2, 2], [1, 1, 2, 2], [3, 1, 2, 2], [1, 2, 3], [2, 2, 2], [1, 1, 3], [4, 3, 3, 2, 2]]
    for ks in cases:
        worst = max(worst, abs(invariants(ks) - invariants_quad(ks)))
    check(NAME, f"Weyl's formula for the Haar integral agrees with numerical Haar quadrature on {len(cases)} tensor products "
                f"(spins up to 2; largest difference {worst:.1e})", worst < 1e-8, eq=eq_int)
    bad3 = 0
    n3 = 0
    for k1 in range(13):
        for k2 in range(13):
            for k3 in range(13):
                n3 += 1
                bad3 += invariants([k1, k2, k3]) != (1 if (ladder.allowed(k1, k2, k3)) else 0)
    check(NAME, f"a trivalent vertex has exactly one intertwiner when the triangle rule of the ladder holds and none otherwise "
                f"(all {n3} triples with spins up to 6)", bad3 == 0, eq=eq_int)
    bad4 = 0
    n4 = 0
    for ks in np.ndindex(7, 7, 7, 7):
        n4 += 1
        bad4 += invariants(list(ks)) != n_triangle(*ks)
    check(NAME, f"a four-valent vertex has as many intertwiners as intermediate spins J allowed by both triangle rules "
                f"(all {n4} four-tuples with spins up to 3)", bad4 == 0, eq=eq_int)
    check(NAME, "four spin-1/2: 2 intertwiners; four spin-1: 3; spins (1/2, 1/2, 1, 1): 2",
          (invariants([1, 1, 1, 1]), invariants([2, 2, 2, 2]), invariants([1, 1, 2, 2])) == (2, 3, 2), eq=eq_int)
    cat = all(invariants([1] * N) == sp.catalan(N // 2) for N in range(2, 41, 2))
    check(NAME, "N spin-1/2 legs have exactly the Catalan number C_{N/2} of invariants for every even N from 2 to 40 (1, 2, 5, 14, ..., "
                "6564120420)", cat and invariants([1] * 40) == 6564120420, eq=eq_int)
    check(NAME, "the invariants are zero when the total spin is half-odd (odd N of spin-1/2) for N = 1, 3, ..., 39",
          all(invariants([1] * N) == 0 for N in range(1, 40, 2)), eq=eq_int)

    # =====================================================================================================
    # entropy and the Immirzi parameter
    # =====================================================================================================
    N_, g_ = sp.symbols("N gamma_I", positive=True)
    bh = kb.equation("Bekenstein-Hawking entropy").relation()
    pl = kb.equation("Planck length").relation()
    area_N = kb.solve_for(kb.equation("Area of N spin-1/2 punctures"), "area").expr
    S_over_kB = sp.simplify(((bh.rhs / sp.Symbol("k_B")).subs(sp.Symbol("A"), area_N)).subs(sp.Symbol("ell_P"), pl.rhs))
    gam = sp.solve(sp.Eq(S_over_kB.subs({sp.Symbol("N"): N_, sp.Symbol("gamma_I"): g_}), N_ * sp.log(2)), g_)
    eq_ent = register(NAME, Equation(
        key="qg_entropy_immirzi", entry=ENTRY, title="Counting spin-1/2 punctures and the value of the Immirzi parameter",
        source_kind="derived", source="the idea is [Rov]; the exponent, the Catalan asymptotics and the algebra are derived here; [DL] for the count over all spins",
        latex_str=(r"W=C_{N/2}\sim\frac{2^{N}}{\sqrt\pi\,(N/2)^{3/2}},\qquad A=4\sqrt3\,\pi N\gamma_{I}\ell_P^{2},\qquad "
                   r"\frac{S}{k_B}=\frac{c^{3}A}{4G\hbar}=\sqrt3\,\pi\gamma_{I}N"
                   r"%%BR%%\frac{S}{k_B}=\ln W\simeq N\ln2\ \Longrightarrow\ \gamma_{I}=\frac{\ln 2}{\pi\sqrt3}\approx0.1274"),
        notes="If a horizon is pierced by N links of spin 1/2 only, its area is N times the minimal area and its entropy is the logarithm "
              "of the number of invariants, N ln 2 at leading order (the Catalan number, or 2^N if the SU(2) constraint is dropped: the "
              "same exponential). Requiring S = A/4 l_P^2 then fixes gamma = ln 2/(pi sqrt 3) = 0.1274. This is a statement about the "
              "spin-1/2 count alone. Counting all spins, with the area fixed, gives a different number (Domagala and Lewandowski: about "
              "0.2375, quoted, not rederived here), so the value of gamma is a property of the counting, not a constant the area "
              "spectrum determines. The horizon is a given area; this is not a derivation of black-hole entropy from the dynamics.",
    ))
    check(NAME, "with the knowledge-base entropy S = k_B c^3 A/(4 G hbar), the Planck length l_P^2 = hbar G/c^3 and the area of N spin-1/2 "
                "punctures, S/k_B = sqrt3 pi gamma N exactly (all of G, hbar, c cancel)",
          sp.simplify(S_over_kB - sp.sqrt(3) * sp.pi * sp.Symbol("gamma_I") * sp.Symbol("N")) == 0, eq=eq_ent)
    check(NAME, "S/k_B = N ln 2 is solved by gamma = ln2/(pi sqrt3) = 0.127384 (sympy solve, then numerical value to 1e-6)",
          len(gam) == 1 and sp.simplify(gam[0] - sp.log(2) / (sp.pi * sp.sqrt(3))) == 0 and abs(float(gam[0]) - 0.127384) < 1e-6, eq=eq_ent)
    mp.mp.dps = 30
    ratios = []
    for n in (100, 1000, 10000):
        cn = mp.mpf(int(sp.catalan(n)))
        ratios.append(cn * mp.sqrt(mp.pi) * mp.mpf(n) ** mp.mpf("1.5") / mp.mpf(4) ** n)
    check(NAME, f"Catalan asymptotics C_n sqrt(pi) n^(3/2)/4^n -> 1: {', '.join(mp.nstr(r, 6) for r in ratios)} at n = 100, 1000, 10000 "
                f"(error ~ 9/(8n))", all(abs(r - 1) < 9 / (8 * n) * 1.2 for r, n in zip(ratios, (100, 1000, 10000))), eq=eq_ent)
    per = [entropy_per_puncture(N) for N in (100, 200, 400, 800, 1600)]
    check(NAME, f"the entropy per puncture (1/N) ln C_(N/2) increases towards ln 2 = 0.6931 along N = 100, 200, 400, 800, 1600: "
                f"{', '.join(f'{x:.4f}' for x in per)}; the shortfall is (3/(2N)) ln(N/2) + ln(pi)/(2N), still 0.0066 at N = 1600",
          all(a < b for a, b in zip(per, per[1:])) and per[-1] < np.log(2) and np.log(2) - per[-1] < 0.007, eq=eq_ent)
    lnC = float(mp.log(mp.mpf(int(sp.catalan(500)))))
    pred = 1000 * np.log(2) - 1.5 * np.log(500) - 0.5 * np.log(np.pi) - 9 / (4 * 1000)
    check(NAME, f"ln W = N ln 2 - (3/2) ln(N/2) - (1/2) ln pi - 9/(4N) + O(N^-2) at N = 1000 (exact {lnC:.6f}, formula {pred:.6f}; "
                f"without the 9/(4N) term the two differ by 2.2e-3)", abs(lnC - pred) < 1e-5, eq=eq_ent)

    # =====================================================================================================
    # the quantum tetrahedron
    # =====================================================================================================
    dim, Qr, herm, comm = tetrahedron()
    eigs = Qr.eigenvals()
    eq_tet = register(NAME, Equation(
        key="qg_quantum_tetrahedron", entry=ENTRY, title="The four-valent vertex of four spin-1/2 as a quantum tetrahedron",
        source_kind="derived", source="derived here from explicit matrices; [Bar] for the interpretation",
        latex_str=(r"Q=\epsilon_{abc}\,J^{(1)}_{a}J^{(2)}_{b}J^{(3)}_{c}\Big|_{\mathrm{Inv}(V_{1/2}^{\otimes4})},\qquad "
                   r"\dim\mathrm{Inv}=2,\qquad \mathrm{spec}\,Q=\Bigl\{-\tfrac{\sqrt3}{4},\,+\tfrac{\sqrt3}{4}\Bigr\},\qquad V\ \propto\ \sqrt{|Q|}"),
        notes="Four spin-1/2 legs meeting at a vertex carry a two-dimensional intertwiner space, the smallest in which a vertex has quantum "
              "numbers of its own. The operator Q, the triple product of three of the four spins, is the volume squared up to a constant: its "
              "absolute value has a single nonzero eigenvalue sqrt3/4, so the quantum tetrahedron of four spin-1/2 has one possible volume "
              "and two orientations. The prefactor that turns Q into a volume (a power of 8 pi gamma l_P^2) and the comparison with the "
              "tables of [Bar] are not done here, and the volume operator of LQG on general vertices is not built.",
    ))
    check(NAME, "the invariant subspace of four spin-1/2 inside the 16-dimensional tensor product is 2-dimensional (nullspace of J^2)",
          dim == 2, eq=eq_tet)
    check(NAME, "Q is Hermitian on the full 16-dimensional space and commutes with all three components of the total spin, so it acts "
                "inside the invariant subspace", herm and comm, eq=eq_tet)
    check(NAME, f"the restriction of Q to the invariant subspace has eigenvalues {sorted(eigs, key=float)} each once, exactly",
          eigs == {-sp.sqrt(3) / 4: 1, sp.sqrt(3) / 4: 1}, eq=eq_tet)
    check(NAME, "the spectrum is symmetric about 0, as the exchange of two legs sends Q to -Q: one volume sqrt(sqrt3/4) with two orientations",
          sp.simplify(Qr.trace()) == 0 and sp.simplify(Qr * Qr - sp.Rational(3, 16) * sp.eye(2)) == sp.zeros(2), eq=eq_tet)

    # =====================================================================================================
    # the geometry of the ladder
    # =====================================================================================================
    b_strong = strong_coefficient()
    cJ_weak, cA_weak = weak_limits()
    eq_lad = register(NAME, Equation(
        key="qg_ladder_geometry", entry=ENTRY, title="The mean-square area of the cut between the two plaquettes, in the ground state of the ladder",
        source_kind="derived", source="derived here from su2_ladder_basis by perturbation theory and the harmonic limit; numerical ground states",
        latex_str=(r"\bigl\langle A_{\mathrm{cut}}^{2}\bigr\rangle_{0}=\bigl(8\pi\gamma_{I}\ell_P^{2}\bigr)^{2}\,\bigl\langle c_2(J)\bigr\rangle_{0},\qquad "
                   r"\bigl\langle c_2(J)\bigr\rangle_{0}\ \xrightarrow{g\to\infty}\ \frac{" + sp.latex(b_strong) + r"}{g^{8}},\qquad "
                   r"g^{2}\bigl\langle c_2(J)\bigr\rangle_{0}\ \xrightarrow{g\to0}\ " + sp.latex(cJ_weak)
                   + r"%%BR%%g^{2}\bigl\langle c_2(j_a)\bigr\rangle_{0}\ \xrightarrow{g\to0}\ " + sp.latex(cA_weak)
                   + r"=\tfrac34\Bigl(\tfrac1{\sqrt6}+\tfrac1{\sqrt{10}}\Bigr)"),
        notes="The surface that separates the two plaquettes cuts only the shared link, so its area is 8 pi gamma l_P^2 sqrt(c_2(J)) and its "
              "mean square is the ground-state expectation of c_2(J). At strong coupling the ground state is the electric vacuum, with J = 0, "
              "dressed to first order by the two strings (1/2, 0, 1/2) and (0, 1/2, 1/2); the probability of finding the shared link excited is "
              "8/(9 g^8), and c_2(1/2) = 3/4, giving 2/(3 g^8) exactly. At weak coupling the ground state spreads over large spins and "
              "<c_2(J)> grows like 1/g^2: the virial theorem for the two coupled oscillators of su2_ladder_weak fixes the coefficient 3/sqrt10, "
              "because the shared link carries E_a - E_b. Only the leading terms are derived; the first correction in g^2 is not. The cut "
              "is a surface of the ladder's own graph, not a physical horizon, and the ladder's Hamiltonian is not that of gravity.",
    ))
    check(NAME, "exact first-order coefficient: <c2(J)> = (2/3)/g^8 from the two strings coupled to the electric vacuum "
                "(rational arithmetic with the exact 6j matrix elements of su2_lattice_ladder)", b_strong == sp.Rational(2, 3), eq=eq_lad)
    check(NAME, "exact weak-coupling limits from the virial theorem: g^2 <c2(J)> -> 3/sqrt10 and g^2 <c2(j_a)> -> (3/4)(1/sqrt6 + 1/sqrt10)",
          sp.simplify(cJ_weak - 3 / sp.sqrt(10)) == 0 and sp.simplify(cA_weak - sp.Rational(3, 4) * (1 / sp.sqrt(6) + 1 / sp.sqrt(10))) == 0, eq=eq_lad)
    strong = {g: geometry(g, 8) for g in (6.0, 8.0)}
    rel = {g: strong[g][1] / (float(b_strong) / g ** 4) - 1 for g in strong}
    check(NAME, "numerical ground states at g^2 = 6 and 8 give <c2(J)> equal to 2/(3 g^8) to " +
          ", ".join(f"{abs(rel[g]):.1e} (g^2={g:g})" for g in rel), abs(rel[6.0]) < 5e-3 and abs(rel[8.0]) < 2e-3, eq=eq_lad)
    w05, w10 = geometry(0.05), geometry(0.1)
    ext_J = 2 * (0.05 * w05[1]) - 0.1 * w10[1]          # g^2 <c2> is linear in g^2 near 0: extrapolate to g^2 = 0
    ext_A = 2 * (0.05 * w05[2]) - 0.1 * w10[2]
    check(NAME, f"numerical ground states at g^2 = 0.05 and 0.1, extrapolated linearly to g^2 = 0: g^2 <c2(J)> = {ext_J:.5f} against "
                f"{float(cJ_weak):.5f}, g^2 <c2(j_a)> = {ext_A:.5f} against {float(cA_weak):.5f} (tolerance 1e-3)",
          abs(ext_J - float(cJ_weak)) < 1e-3 and abs(ext_A - float(cA_weak)) < 1e-3, eq=eq_lad)
    ref = geometry(0.3, 30)[1]
    d4, d6, d12 = (abs(geometry(0.3, k)[1] - ref) for k in (8, 12, 24))
    check(NAME, f"truncation convergence at g^2 = 0.3: <c2(J)> at j_max = 4, 6 and 12 differs from the j_max = 15 value by "
                f"{d4:.1e}, {d6:.1e} and {d12:.1e}, falling at each step to the rounding level",
          d4 < 1e-4 and d6 < 1e-9 and d12 < 1e-11 and d4 > d6 > d12 * 0.5, eq=eq_lad)
    grid = (0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
    vals = [geometry(g)[1] for g in grid]
    check(NAME, "<c2(J)> decreases strictly with g^2 over g^2 = 0.1 .. 8: " + ", ".join(f"{v:.3g}" for v in vals),
          all(a > b for a, b in zip(vals, vals[1:])), eq=eq_lad)
    check(NAME, "the shared-link fluctuation is larger than a plaquette link's at weak coupling and the ratio of the limits is "
                f"{float(cJ_weak / cA_weak):.3f} (= (3/sqrt10)/((3/4)(1/sqrt6 + 1/sqrt10)))",
          float(cJ_weak) > float(cA_weak) and abs(float(cJ_weak / cA_weak) - 1.7459) < 1e-3, eq=eq_lad)
    return {g: geometry(g) for g in grid}


if __name__ == "__main__":
    tab = build_and_verify()
    print("g^2      E0          <c2(J)>        <c2(j_a)>")
    for g, (e0, cj, ca) in tab.items():
        print(f"{g:6.3f}  {e0:10.5f}  {cj:12.5e}  {ca:12.5e}")
