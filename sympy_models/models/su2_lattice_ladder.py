"""Two-plaquette SU(2) ladder: the first system with a volume dependence (companion to su2_lattice_gap.py).

Two plaquettes sharing one link, open boundary: 7 links, 6 vertices, the Kogut-Susskind Hamiltonian
H = (g^2/2) sum_l E_l^2 + (1/g^2) sum_p (2 - Tr U_p) in lattice units (see su2_lattice_gap.py for conventions).

Basis. Gauge-invariant states are spin networks on the theta graph with three edges: the three outer links of
plaquette a carry the same spin j_a, those of plaquette b carry j_b, and the shared link carries J, with the triangle
rule and j_a + j_b + J integer. In the two plaquette holonomies U_a, U_b (a spanning-tree gauge: five links fixed, two
loop variables left, the same idea as the forest gauge of Romanov's repositories, but in the Hamiltonian setting and
without any of his results) the states are Psi = Tr[P_J (D^{j_a}(U_a) x D^{j_b}(U_b))], with norm^2 d_J/(d_a d_b).
The electric energy is (g^2/2)[3 j_a(j_a+1) + 3 j_b(j_b+1) + J(J+1)].

The plaquette operator. Multiplication by chi_{1/2}(U_a) moves (j_a, J) to (j_a +- 1/2, J +- 1/2) with matrix element
sqrt(d_a d_a' d_J d_J') {j_a 1/2 j_a'; J' j_b J}^2, and chi_{1/2}(U_b) likewise with a and b exchanged. The square of the
6j symbol, not the symbol, appears because the left and right index sets of the loop wavefunction carry the same
recoupling. The squares have the closed forms of Edmonds' Table 5 for a 6j symbol with an entry 1/2.

Checked: the closed forms against sympy's exact 6j symbols (exact rational equality, every transition up to
j = 4); the matrix element against an independent construction from explicit Clebsch-Gordan contractions (every
transition up to j = 5/2, both plaquettes, to 1e-12); the reduction to the one-plaquette matrix element 1 when j_b = 0
(or j_a = 0); Hermiticity; that the two plaquette operators commute exactly on the whole truncated space; that the
reflection a <-> b maps one onto the other. The weak-coupling limit is derived independently from two coupled
three-dimensional oscillators (kinetic matrix [[4, -1], [-1, 4]] in the loop electric fields, modes with
omega_s = sqrt(3/2), omega_d = sqrt(5/2)) and compared with the numerical spectrum, including the reflection parities.
The strong-coupling splitting of the doubly degenerate first excited level is derived by exact degenerate
second-order perturbation theory in exact arithmetic and compared with the numbers.

Not tested: three or more plaquettes (the strip formula in su2_ladder_weak is a prediction from the same harmonic
argument, verified here only for one and two plaquettes), periodic boundary conditions, any scaling window or
continuum limit, higher-order strong-coupling terms, the first correction in g^2 at weak coupling (only the limit is
checked), the 2+1D and 3+1D theories. Nothing here bears on existence or a mass gap of Yang-Mills on R^4.

Related, not used: the forest-gauge results of Romanov (forest_gauge.py, forest_gauge_4d.py) are rigorous statements
about the infinite periodic lattice at strong coupling. Nothing in this module uses their results or code.
"""
import itertools
import sys
from fractions import Fraction as Fr
from functools import lru_cache
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as sla
import sympy as sp
from sympy.physics.wigner import clebsch_gordan, wigner_6j

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register
from models.su2_lattice_gap import gap_num as plaquette_gap

NAME = "su2_lattice_ladder"
ENTRY = "su2-lattice-gap"
KS = "Kogut-Susskind, Phys. Rev. D 11, 395 (1975)"
EDM = "Edmonds, Angular Momentum in Quantum Mechanics (1957), Table 5 (6j symbols with an entry 1/2)"


# ---------------------------------------------------------------------------------------------------------------
# spins are carried as k = 2j, so every quantity below is an integer
# ---------------------------------------------------------------------------------------------------------------
def allowed(ka, kb, K):
    return abs(ka - kb) <= K <= ka + kb and (ka + kb + K) % 2 == 0


def six_sq(ka, kb, K, kap, Kp):
    """Exact square of {ja 1/2 ja'; J' jb J} (a Fraction) from the closed forms for a 6j symbol with an entry 1/2.
    Reordered by two symmetries to {a b c; 1/2 e f} with a = jb, b = J', c = ja', e = ja, f = J."""
    a, b, c, e, f = kb, Kp, kap, ka, K
    if b == 0:                       # zero entry: the symbol is (+-) 1/sqrt((2c+1)(2f+1)) with f = 1/2
        return Fr(1, 2 * (c + 1))
    if c == 0:                       # likewise with e = 1/2
        return Fr(1, 2 * (b + 1))
    if e == c - 1 and f == b + 1:    # {a b c; 1/2 c-1/2 b+1/2}
        return Fr(((a - b + c) // 2) * ((a + b - c) // 2 + 1), (b + 1) * (b + 2) * c * (c + 1))
    if e == c - 1 and f == b - 1:    # {a b c; 1/2 c-1/2 b-1/2}
        return Fr(((a + b + c) // 2 + 1) * ((-a + b + c) // 2), b * (b + 1) * c * (c + 1))
    if e == c + 1 and f == b - 1:    # swap b <-> c in the first form
        return Fr(((a + b - c) // 2) * ((a - b + c) // 2 + 1), (c + 1) * (c + 2) * b * (b + 1))
    if e == c + 1 and f == b + 1:    # swap upper/lower in columns 2, 3: {a f e; 1/2 b c}, second form with (b, c) -> (f, e)
        return Fr(((a + f + e) // 2 + 1) * ((-a + f + e) // 2), f * (f + 1) * e * (e + 1))
    raise ValueError((ka, kb, K, kap, Kp))


def w_a(ka, kb, K, kap, Kp):
    return float(np.sqrt((ka + 1) * (kap + 1) * (K + 1) * (Kp + 1)) * float(six_sq(ka, kb, K, kap, Kp)))


def w_b(ka, kb, K, kbp, Kp):
    return float(np.sqrt((kb + 1) * (kbp + 1) * (K + 1) * (Kp + 1)) * float(six_sq(kb, ka, K, kbp, Kp)))


def w_exact(d_prod, sq):
    """sqrt(d_prod) * sq as a sympy number."""
    return sp.sqrt(sp.Integer(d_prod)) * sp.Rational(sq.numerator, sq.denominator)


def steps(ka, kb, K, kmax):
    """Transitions out of (ka, kb, K): ('a'|'b', new state)."""
    for da, dJ in itertools.product((1, -1), repeat=2):
        kap, Kp = ka + da, K + dJ
        if kap >= 0 and Kp >= 0 and kap <= kmax and allowed(kap, kb, Kp):
            yield "a", (kap, kb, Kp)
    for db, dJ in itertools.product((1, -1), repeat=2):
        kbp, Kp = kb + db, K + dJ
        if kbp >= 0 and Kp >= 0 and kbp <= kmax and allowed(ka, kbp, Kp):
            yield "b", (ka, kbp, Kp)


def element(which, s, t):
    return w_a(*s[:2], s[2], t[0], t[2]) if which == "a" else w_b(s[0], s[1], s[2], t[1], t[2])


def e_el(s):
    ka, kb, K = s
    return Fr(3 * ka * (ka + 2) + 3 * kb * (kb + 2) + K * (K + 2), 4)


def build(kmax, g2):
    idx, states = {}, []
    for ka in range(kmax + 1):
        for kb in range(kmax + 1):
            for K in range(abs(ka - kb), ka + kb + 1, 2):
                idx[(ka, kb, K)] = len(states)
                states.append((ka, kb, K))
    n = len(states)
    rows, cols, vals = [], [], []
    diag = np.array([0.5 * g2 * float(e_el(s)) + 4 / g2 for s in states])
    for i, s in enumerate(states):
        for which, t in steps(*s, kmax):
            rows.append(i)
            cols.append(idx[t])
            vals.append(-element(which, s, t) / g2)
    H = sps.csr_matrix((vals, (rows, cols)), shape=(n, n)) + sps.diags(diag)
    return H, states, idx


def operators(kmax):
    idx, states = {}, []
    for ka in range(kmax + 1):
        for kb in range(kmax + 1):
            for K in range(abs(ka - kb), ka + kb + 1, 2):
                idx[(ka, kb, K)] = len(states)
                states.append((ka, kb, K))
    n = len(states)
    data = {"a": ([], [], []), "b": ([], [], [])}
    for i, s in enumerate(states):
        for which, t in steps(*s, kmax):
            r, c, v = data[which]
            r.append(i)
            c.append(idx[t])
            v.append(element(which, s, t))
    mats = {k: sps.csr_matrix((v, (r, c)), shape=(n, n)) for k, (r, c, v) in data.items()}
    perm = np.array([idx[(kb, ka, K)] for (ka, kb, K) in states])
    return mats["a"], mats["b"], states, idx, perm


def kmax_for(g2):
    return 36 if g2 < 0.15 else 30 if g2 < 0.4 else 24 if g2 < 0.8 else 18 if g2 < 1.5 else 12


def levels(g2, nev=5, kmax=None):
    """Lowest levels and their reflection parity (a <-> b) of the ladder Hamiltonian."""
    H, states, idx = build(kmax or kmax_for(g2), g2)
    perm = np.array([idx[(kb, ka, K)] for (ka, kb, K) in states])
    w, v = sla.eigsh(H, k=nev, sigma=0.0, which="LM")
    o = np.argsort(w)
    w, v = w[o], v[:, o]
    par = np.array([float(v[:, i] @ v[perm, i]) for i in range(nev)])
    return w, par


# ---------------------------------------------------------------------------------------------------------------
# independent construction from explicit Clebsch-Gordan contractions (used only for checks, small spins)
# ---------------------------------------------------------------------------------------------------------------
@lru_cache(None)
def cg(k1, k2, K):
    j1, j2, J = sp.Rational(k1, 2), sp.Rational(k2, 2), sp.Rational(K, 2)
    C = np.zeros((k1 + 1, k2 + 1, K + 1))
    for a in range(k1 + 1):
        for b in range(k2 + 1):
            m1, m2 = -j1 + a, -j2 + b
            M = m1 + m2
            if abs(M) <= J:
                C[a, b, int(M + J)] = float(clebsch_gordan(j1, j2, J, m1, m2, M))
    return C


def norm2(ka, kb, K):
    return (K + 1) / ((ka + 1) * (kb + 1))


def w_a_cg(ka, kb, K, kap, Kp):
    """<Psi'| chi_{1/2}(U_a) |Psi> from the Clebsch-Gordan series D^j D^{1/2} = sum D^{j'} and Schur orthogonality."""
    C1, C2, C3 = cg(ka, kb, K), cg(ka, 1, kap), cg(kap, kb, Kp)
    T = np.einsum("mpM,nqM,msu,nsv->upvq", C1, C1, C2, C2, optimize=True)
    c = np.einsum("upvq,upN,vqN->", T, C3, C3, optimize=True) / (Kp + 1)
    return c * np.sqrt(norm2(kap, kb, Kp) / norm2(ka, kb, K))


def w_b_cg(ka, kb, K, kbp, Kp):
    C1, C2, C3 = cg(ka, kb, K), cg(kb, 1, kbp), cg(ka, kbp, Kp)
    T = np.einsum("mpM,nqM,psu,qsv->munv", C1, C1, C2, C2, optimize=True)
    c = np.einsum("munv,muN,nvN->", T, C3, C3, optimize=True) / (Kp + 1)
    return c * np.sqrt(norm2(ka, kbp, Kp) / norm2(ka, kb, K))


# ---------------------------------------------------------------------------------------------------------------
# strong coupling: exact degenerate second-order perturbation theory
# ---------------------------------------------------------------------------------------------------------------
def neighbours_exact(s):
    """Neighbours of s with exact matrix elements (sympy numbers), any size."""
    out = []
    ka, kb, K = s
    big = 10**6
    for which, t in steps(ka, kb, K, big):
        if which == "a":
            dp = (ka + 1) * (t[0] + 1) * (K + 1) * (t[2] + 1)
            sq = six_sq(ka, kb, K, t[0], t[2])
        else:
            dp = (kb + 1) * (t[1] + 1) * (K + 1) * (t[2] + 1)
            sq = six_sq(kb, ka, K, t[1], t[2])
        out.append((t, w_exact(dp, sq)))
    return out


def strong_splitting():
    """Second-order gap of the two lowest excited states, in the sectors symmetric and antisymmetric under a <-> b.
    Returns (c_plus, c_minus) with gap_pm = (3/2) g^2 + c_pm / g^6 + O(g^-14)."""
    t = sp.symbols("t", positive=True)           # t = g^2
    A, B, vac = (1, 0, 1), (0, 1, 1), (0, 0, 0)
    ES = sp.Rational(3, 2) * t
    El = lambda s: sp.Rational(e_el(s).numerator, e_el(s).denominator) * t / 2
    assert El(A) == ES and El(B) == ES
    nb = {s: dict(neighbours_exact(s)) for s in (A, B, vac)}

    def second(X, Y, E0):
        tot = 0
        for m, vx in nb[X].items():
            if m in (A, B):          # only the degenerate pair is excluded; the vacuum is a legitimate intermediate state
                continue
            vy = nb[Y].get(m)
            if vy is None:
                continue
            tot += (vx / t) * (vy / t) / (E0 - El(m))
        return sp.simplify(tot)
    hAA, hAB = second(A, A, ES), second(A, B, ES)
    vacshift = sum((v / t) ** 2 / (0 - El(m)) for m, v in nb[vac].items())
    cplus = sp.simplify((hAA + hAB - vacshift) * t**3)
    cminus = sp.simplify((hAA - hAB - vacshift) * t**3)
    return cplus, cminus, hAA * t**3, hAB * t**3, vacshift * t**3


def build_and_verify():
    # =====================================================================================================
    # the matrix element
    # =====================================================================================================
    eq_basis = register(NAME, Equation(
        key="su2_ladder_basis", entry=ENTRY, title="Two-plaquette ladder: basis, energies and plaquette operators",
        source_kind="derived", source=f"derived here; {KS} for H; {EDM} for the 6j closed forms",
        latex_str=(r"|j_a,j_b,J\rangle,\ \ \triangle(j_a,j_b,J),\quad "
                   r"H=\frac{g^{2}}{2}\bigl[3C(j_a)+3C(j_b)+C(J)\bigr]+\frac{1}{g^{2}}\bigl(4-W_a-W_b\bigr),\ \ C(j)=j(j+1)"
                   r"%%BR%%\langle j_a{+}\epsilon,j_b,J{+}\eta|W_a|j_a,j_b,J\rangle=\sqrt{d_{j_a}d_{j_a'}d_Jd_{J'}}\,"
                   r"\begin{Bmatrix}j_a&\tfrac12&j_a'\\J'&j_b&J\end{Bmatrix}^{2},\ \ \epsilon,\eta=\pm\tfrac12"),
        notes="Spin network on the theta graph: three outer links of plaquette a in spin j_a, three of b in j_b, the shared "
              "link in J. In loop variables the states are Tr[P_J (D^{j_a}(U_a) x D^{j_b}(U_b))] with norm^2 d_J/(d_a d_b). "
              "The 6j symbol enters squared because the left and right index sets of the wavefunction carry the same "
              "recoupling. W_b is the same with a and b exchanged. Open boundary: seven links, two plaquettes.",
    ))
    worst = Fr(0)
    n6 = 0
    for ka, kb in itertools.product(range(0, 9), repeat=2):
        for K in range(abs(ka - kb), ka + kb + 1, 2):
            for which, t in steps(ka, kb, K, 9):
                if which != "a":
                    continue
                ref = wigner_6j(sp.Rational(ka, 2), sp.Rational(1, 2), sp.Rational(t[0], 2), sp.Rational(t[2], 2),
                                sp.Rational(kb, 2), sp.Rational(K, 2)) ** 2
                mine = six_sq(ka, kb, K, t[0], t[2])
                worst = max(worst, abs(sp.Rational(mine.numerator, mine.denominator) - ref))
                n6 += 1
    check(NAME, f"closed forms for the square of {{ja 1/2 ja'; J' jb J}} equal sympy's exact 6j symbol squared, as exact "
                f"rationals, for all {n6} transitions with ja, jb up to 4", worst == 0, eq=eq_basis)

    ww = 0.0
    nw = 0
    for ka, kb in itertools.product(range(0, 6), repeat=2):
        for K in range(abs(ka - kb), ka + kb + 1, 2):
            for which, t in steps(ka, kb, K, 5):
                ref = w_a_cg(ka, kb, K, t[0], t[2]) if which == "a" else w_b_cg(ka, kb, K, t[1], t[2])
                ww = max(ww, abs(ref - element(which, (ka, kb, K), t)))
                nw += 1
    check(NAME, f"the matrix elements of W_a and W_b equal an independent construction from explicit Clebsch-Gordan "
                f"contractions for all {nw} transitions with spins up to 5/2 (largest difference {ww:.1e})",
          ww < 1e-12, eq=eq_basis)
    red = all(abs(element("a", (ka, 0, ka), (ka + d, 0, ka + d)) - 1) < 1e-14 and
              abs(element("b", (0, ka, ka), (0, ka + d, ka + d)) - 1) < 1e-14
              for ka in range(0, 20) for d in (1, -1) if ka + d >= 0)
    check(NAME, "j_b = 0 (or j_a = 0): the plaquette matrix element is 1, the one-plaquette tridiagonal chain of su2_plaq_hamiltonian, "
                "and the electric energy 3C(j) + C(j) = 4C(j) is that of four links", red, eq=eq_basis)
    Wa, Wb, states, idx, perm = operators(14)
    comm = abs((Wa @ Wb - Wb @ Wa)).max()
    check(NAME, f"W_a and W_b are symmetric, commute exactly on the whole truncated space (largest entry {comm:.1e}), and the "
                f"reflection a <-> b maps one to the other",
          abs(Wa - Wa.T).max() == 0 and comm < 1e-12 and abs(Wa[perm][:, perm] - Wb).max() == 0, eq=eq_basis)
    g_a, g_b = levels(0.3, 5, 24)[0], levels(0.3, 5, 30)[0]
    check(NAME, f"truncation convergence at g^2 = 0.3: the lowest five levels at j_max = 12 and 15 agree to 1e-7 "
                f"(largest difference {np.max(abs(g_a - g_b)):.1e})", np.max(abs(g_a - g_b)) < 1e-7, eq=eq_basis)

    # =====================================================================================================
    # weak coupling: two coupled oscillators
    # =====================================================================================================
    K2 = sp.Matrix([[4, -1], [-1, 4]])
    w_of = lambda lam: sp.sqrt(sp.Rational(2 * lam, 4))   # omega^2 = 2 c, c = (kinetic eigenvalue)/4
    ws, wd = w_of(3), w_of(5)
    levels_pred = [sp.simplify(2 * ws), sp.simplify(ws + wd), sp.simplify(2 * wd), sp.simplify(4 * ws)]
    eq_weak = register(NAME, Equation(
        key="su2_ladder_weak", entry=ENTRY, title="Weak-coupling limit of the ladder spectrum and the strip of N plaquettes",
        source_kind="derived", source="derived here from su2_ladder_basis (harmonic approximation); numerical check",
        latex_str=(r"\Delta_{g\to0}=2\omega_s=\sqrt6,\ \ \omega_s=\sqrt{\tfrac32},\ \ \omega_d=\sqrt{\tfrac52},\quad "
                   r"\text{next levels }\ \omega_s+\omega_d\ (\text{odd}),\ \ 2\omega_d=\sqrt{10},\ \ 4\omega_s%%BR%%"
                   r"\text{strip of }N\text{ plaquettes (prediction): }\ \Delta_{g\to0}=2\sqrt{2-\cos\tfrac{\pi}{N+1}}\ \xrightarrow{N\to\infty}\ 2"),
        notes="Near U = 1 each plaquette angle is a three-dimensional oscillator. The electric energy is (g^2/2)[3E_a^2 + "
              "3E_b^2 + (E_a - E_b)^2], the kinetic matrix [[4, -1], [-1, 4]], with eigenvalues 3 (symmetric mode) and 5; "
              "omega = sqrt(kinetic eigenvalue/2). Gauge invariance keeps total angular momentum 0 of the pair of modes, "
              "so the lowest excitation is the radial one in the soft mode, 2 omega_s. For N plaquettes in a row the "
              "kinetic matrix is 4 - (adjacency of the path), with lowest eigenvalue 4 - 2cos(pi/(N+1)); only N = 1 "
              "(2 sqrt 2) and N = 2 (sqrt 6) are verified here. The limit is the harmonic approximation of a strip "
              "of width one and is not a continuum statement.",
    ))
    check(NAME, "kinetic matrix [[4,-1],[-1,4]] has eigenvalues 3 and 5, so omega_s = sqrt(3/2) and omega_d = sqrt(5/2); "
                "the one-plaquette value (kinetic 4) gives omega = sqrt 2 as in su2_plaq_weak",
          sorted(K2.eigenvals()) == [3, 5] and sp.simplify(ws - sp.sqrt(sp.Rational(3, 2))) == 0
          and sp.simplify(wd - sp.sqrt(sp.Rational(5, 2))) == 0 and sp.simplify(w_of(4) - sp.sqrt(2)) == 0, eq=eq_weak)
    check(NAME, "predicted excitation energies: 2 omega_s = sqrt 6, omega_s + omega_d = sqrt(3/2) + sqrt(5/2), "
                "2 omega_d = sqrt 10, 4 omega_s = 2 sqrt 6",
          levels_pred[0] == sp.sqrt(6) and levels_pred[2] == sp.sqrt(10) and levels_pred[3] == 2 * sp.sqrt(6), eq=eq_weak)
    N_ = sp.symbols("N", positive=True, integer=True)
    delta_N = lambda n: 2 * sp.sqrt(2 - sp.cos(sp.pi / (n + 1)))
    check(NAME, "strip formula: N = 1 gives 2 sqrt 2 (the one plaquette), N = 2 gives sqrt 6 (this ladder), and N -> infinity gives 2",
          sp.simplify(delta_N(1) - 2 * sp.sqrt(2)) == 0 and sp.simplify(delta_N(2) - sp.sqrt(6)) == 0
          and sp.limit(2 * sp.sqrt(2 - sp.cos(sp.pi / (N_ + 1))), N_, sp.oo) == 2, eq=eq_weak)
    num = {}
    for g in (0.05, 0.1):
        w, par = levels(g, 5, 44 if g < 0.07 else 36)
        num[g] = (w[1:] - w[0], par[1:])
    extrap = 2 * num[0.05][0] - num[0.1][0]          # linear extrapolation to g^2 = 0
    target = np.array([float(x) for x in levels_pred[:3]])
    err = np.abs(extrap[:3] - target)
    check(NAME, f"numerical spectrum, linearly extrapolated from g^2 = 0.05 and 0.1 to g^2 = 0: lowest three excitations "
                f"{np.round(extrap[:3], 4).tolist()} against {np.round(target, 4).tolist()} (differences {np.round(err, 4).tolist()}, "
                f"tolerance 6e-3)", bool(np.all(err < 6e-3)), eq=eq_weak)
    check(NAME, f"reflection parities of the three lowest excitations are (+, -, +), as the symmetric mode (2 omega_s), "
                f"the mixed state and 2 omega_d predict (numerical: {np.round(num[0.05][1][:3], 2).tolist()})",
          np.allclose(num[0.05][1][:3], [1, -1, 1], atol=0.01), eq=eq_weak)
    check(NAME, f"the ladder's lowest excitation at g^2 = 0.05 ({num[0.05][0][0]:.4f}) is below the one-plaquette gap "
                f"({plaquette_gap(0.05):.4f}); their weak-coupling ratio is sqrt6/(2 sqrt2) = sqrt3/2",
          num[0.05][0][0] < plaquette_gap(0.05) and sp.simplify(sp.sqrt(6) / (2 * sp.sqrt(2)) - sp.sqrt(3) / 2) == 0, eq=eq_weak)

    # =====================================================================================================
    # strong coupling
    # =====================================================================================================
    cplus, cminus, hAA, hAB, vac = strong_splitting()
    eq_strong = register(NAME, Equation(
        key="su2_ladder_strong", entry=ENTRY, title="Strong-coupling splitting of the first excited level of the ladder",
        source_kind="derived", source="derived here by exact degenerate second-order perturbation theory from su2_ladder_basis",
        latex_str=(r"\Delta_\pm=\frac32\,g^{2}+\frac{c_\pm}{g^{6}}+O(g^{-14}),\qquad c_+=" + sp.latex(cplus)
                   + r",\ \ c_-=" + sp.latex(cminus) + r"\ \ (\text{one plaquette: }\tfrac{14}{15})"),
        notes="The states (1/2, 0, 1/2) and (0, 1/2, 1/2) are a string around plaquette a and around b, degenerate at "
              "3 g^2/2. They couple only at second order, through (1,1,0)/(1,1,2)-type intermediate states, and split into "
              "combinations even (+) and odd (-) under a <-> b. The vacuum, which lies in the even sector, shifts down by 4/(3 g^6) and "
              "pushes the even string state up by the same amount through the direct coupling W_a + W_b; the odd state does not "
              "see it, and at strong coupling its gap equals the one-plaquette gap through this order. Exact rationals.",
    ))
    check(NAME, f"exact second-order coefficients: diagonal {sp.nsimplify(hAA)}, off-diagonal {sp.nsimplify(hAB)}, vacuum "
                f"{sp.nsimplify(vac)} (units g^-6); c_+ = {cplus}, c_- = {cminus}; all rational",
          all(sp.nsimplify(x).is_Rational for x in (cplus, cminus, hAA, hAB, vac)) and vac == sp.Rational(-4, 3), eq=eq_strong)
    obs = {}
    for g in (6.0, 8.0):
        w, par = levels(g, 5, 8)
        obs[g] = (w[1] - w[0], int(np.sign(par[1])), w[2] - w[0], int(np.sign(par[2])))
    pred_ok = True
    detail = []
    for g, (d1, p1, d2, p2) in obs.items():
        pp = 1.5 * g + float(cplus) / g**3
        pm = 1.5 * g + float(cminus) / g**3
        got = {p1: d1, p2: d2}
        e1 = abs(got.get(1, np.nan) - pp)
        e2 = abs(got.get(-1, np.nan) - pm)
        pred_ok &= bool(e1 < 5e-4 and e2 < 5e-4) and p1 * p2 < 0
        detail.append(f"g^2={g}: even {got.get(1, float('nan')):.6f} vs {pp:.6f}, odd {got.get(-1, float('nan')):.6f} vs {pm:.6f}")
    check(NAME, "the two lowest excitations at g^2 = 6 and 8 have opposite reflection parity and match 3g^2/2 + c_pm/g^6 "
                "to 5e-4: " + "; ".join(detail), pred_ok, eq=eq_strong)
    lower = "even" if cplus < cminus else "odd"
    check(NAME, f"the {lower} combination is the lower of the pair at strong coupling (c_+ = {cplus} against c_- = {cminus})",
          (cplus < cminus) == (obs[8.0][1] == 1), eq=eq_strong)

    # =====================================================================================================
    # the gap against g^2, and against one plaquette
    # =====================================================================================================
    grid = (0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
    rows = {}
    for g in grid:
        w, par = levels(g, 4)
        rows[g] = (w[1] - w[0], par[1], plaquette_gap(g))
    ratio = {g: r[0] / r[2] for g, r in rows.items()}
    eq_gap = register(NAME, Equation(
        key="su2_ladder_gap", entry=ENTRY, title="The ladder gap against the one-plaquette gap",
        source_kind="derived", source="derived here from su2_ladder_basis (numerical)",
        latex_str=(r"\Delta_{\mathrm{ladder}}(g)\le\Delta_{\mathrm{plaquette}}(g),\qquad "
                   r"\frac{\Delta_{\mathrm{ladder}}}{\Delta_{\mathrm{plaquette}}}\ \to\ \frac{\sqrt3}{2}\ (g\to0),\quad \to 1\ (g\to\infty)"),
        notes="Numerical on the grid g^2 = 0.1, 0.25, 0.5, 1, 2, 4, 8 with truncation checked. A larger system has a "
              "lower gap in lattice units at weak coupling (the harmonic argument gives sqrt 3/2 for two plaquettes against "
              "one), and the two agree at strong coupling, where the first excitation is a string around one plaquette "
              "of energy 3 g^2/2 whichever system it lives in. This is a volume dependence for two volumes only.",
    ))
    check(NAME, "ladder gap <= one-plaquette gap on the grid g^2 = 0.1 .. 8; ratios " +
          ", ".join(f"{g}: {ratio[g]:.4f}" for g in grid), all(r <= 1 + 1e-9 for r in ratio.values()), eq=eq_gap)
    check(NAME, f"the ratio tends to sqrt3/2 = 0.8660 at weak coupling (0.8660 + 0.05 at g^2 = 0.1: {ratio[0.1]:.4f}) and to 1 at "
                f"strong coupling ({ratio[8.0]:.5f} at g^2 = 8)",
          abs(ratio[0.1] - np.sqrt(3) / 2) < 0.05 and ratio[8.0] > 0.99, eq=eq_gap)
    gs = np.array([0.5, 0.7, 0.9, 1.1, 1.3, 1.6])
    gv = np.array([levels(g, 3)[0][1] - levels(g, 3)[0][0] for g in gs])
    imin = int(np.argmin(gv))
    check(NAME, f"the ladder gap has an interior minimum near g^2 = 1 (scan 0.5 .. 1.6: smallest {gv[imin]:.4f} at g^2 = {gs[imin]}, "
                f"larger at both ends)", 0 < imin < len(gs) - 1, eq=eq_gap)
    return {g: (rows[g][0], rows[g][1], rows[g][2]) for g in grid}


if __name__ == "__main__":
    tab = build_and_verify()
    print("g^2     ladder gap  parity  plaquette gap  ratio")
    for g, (a, p, b) in tab.items():
        print(f"{g:6.3f}  {a:.6f}   {p:+.0f}     {b:.6f}      {a/b:.4f}")
