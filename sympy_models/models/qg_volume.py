"""The volume (triple product) operator on four-valent intertwiner spaces (companion to qg_spin_network.py and su2_lattice_strip.py).

A trivalent vertex carries a unique intertwiner and no volume: the three spins fix a triangle and nothing else. The first vertex with a
volume of its own is four-valent. For four legs of spins j1..j4 the invariant subspace Inv(V_j1 x V_j2 x V_j3 x V_j4) has one state for
each intermediate spin K allowed by both triangle rules, and on it the triple product

    Q = J1 . (J2 x J3) = eps_abc J1_a J2_b J3_c

is a Hermitian operator that commutes with the total spin. Its absolute value is, up to a constant that this module does not fix,
the volume of the quantum tetrahedron (Barbieri 1998; Rovelli-Smolin): V is proportional to sqrt|Q|, so V^2 ~ |Q|. qg_spin_network.py
computed the one case of four spin-1/2 legs (two states, Q = +-sqrt3/4); this module takes general spins.

Where the four-valent vertex lives. The ladder and the strip of plaquettes have no volume. In a strip every vertex of the lattice has
at most three links, and after the gauge reduction of su2_lattice_ladder.py and su2_lattice_strip.py every vertex of the spin-network graph
is trivalent: the strip has no intertwiner multiplicity at all. A four-valent lattice vertex first appears for three plaquettes in an
L shape, the vertex shared by all three (the same happens at the centre of a 2 x 2 block). This module therefore computes the volume
operator on the intertwiner space of that vertex, and one strong-coupling result for the L shape that needs only the strip's own
amplitude. It does NOT build the Hamiltonian of the L shape: that needs recoupling coefficients between the (ab, bt) and (ac, cr)
channels at the four-valent vertex for all spins, which is the next step, and no phase of it is guessed here.

What is registered, with its checks:

  qg_volume_triple_product  Q on the invariant subspace, built from explicit matrices for eight sets of spins from (1/2)^4 up to (2)^4:
                            Hermitian, commuting with the total spin, preserving the invariant subspace, of dimension equal to Weyl's
                            integration formula and to the triangle count, with a spectrum symmetric under Q -> -Q (time reversal);
                            the operator identity Q = i [J1.J2, J2.J3] on the whole tensor product; Q = 0 when a leg has spin 0.
  qg_volume_spectra         exact eigenvalues: the coupled-basis form Q_KK' = i (x_K - x_K') <K|J2.J3|K'> with squared 6j recouplings
                            reproduces the explicit spectrum to 1e-10 for all the sets, gives exact eigenvalues as radicals, e.g.
                            +-sqrt3/4 for four spin-1/2, {0, +-sqrt3} for four spin-1, +-3 sqrt3/4 and +-3 sqrt35/4 for four spin-3/2,
                            and the traces Tr Q = Tr Q^3 = 0, Tr Q^2 = 3/8, 6, 171/4, 198 for k = 1..4.
  qg_volume_semiclassical   for four equal spins j the largest eigenvalue over (4/(3 sqrt 3)) [j(j+1)]^{3/2}, the triple product of a
                            regular tetrahedron of edge-vector length sqrt(j(j+1)), increases from 0.7937 (j = 3/2) to 0.921 (j = 8)
                            and stays below 1.
  qg_volume_l_shape         the L-shaped three-plaquette cluster: the maximum vertex degree of the strip is 3 and of the L shape 4,
                            attained at one vertex; the smallest set of plaquette strings that lights all four legs of that vertex is
                            {b, c}; on the vector that two strings through the vertex give, Q^2 = 3/16 exactly; the amplitude of a
                            two-string state is 4/(9 g^8) (second-order perturbation theory, checked against the disjoint strings
                            a and c of the strip's ground state at g^2 = 6, 8, 10); hence <Q^2>_0 = 1/(27 g^16) at strong coupling.

Not tested: the L-shape Hamiltonian and its ground state; the prefactor that turns Q into a volume (so only the dimensionless
triple product is computed, and the constant relating V to sqrt|Q| depends on the volume operator one chooses: Rovelli-Smolin,
Ashtekar-Lewandowski and the Barbieri tetrahedron differ); the volume operator on vertices with more than four legs; the expectation of
Q^2 beyond the leading order at strong coupling, and any weak-coupling value; the semiclassical limit beyond the monotone trend
shown (the ratio is 0.94 at j = 12, and the claim that it tends to 1 is a trend, not a proof); any dynamics of gravity. The strong-coupling
result uses only the structure of the state (two loops through a vertex) and the strip's amplitude for two disjoint strings: it
assumes, and does not test in the L shape itself, that the amplitude of two strings that share a vertex but no link is the same as for
two that share nothing. Nothing here bears on existence or a mass gap of Yang-Mills on R^4, and nothing is a model of quantum gravity.
"""
import itertools
import sys
from pathlib import Path

import numpy as np
import scipy.sparse.linalg as sla
import sympy as sp
from sympy.physics.wigner import wigner_6j

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register
from models import qg_spin_network as qg
from models import su2_lattice_strip as strip

NAME = "qg_volume"
ENTRY = "spin-network-qg"
R = sp.Rational


# ---------------------------------------------------------------------------------------------------------------
# explicit matrices
# ---------------------------------------------------------------------------------------------------------------
def _jm(k):
    j = k / 2
    m = -j + np.arange(k + 1)
    Jp = np.zeros((k + 1, k + 1), dtype=complex)
    for i in range(k):
        Jp[i + 1, i] = np.sqrt(j * (j + 1) - m[i] * (m[i] + 1))
    return [(Jp + Jp.conj().T) / 2, (Jp - Jp.conj().T) / 2j, np.diag(m).astype(complex)]


def _embed(ks, leg, op):
    mats = [np.eye(k + 1, dtype=complex) for k in ks]
    mats[leg] = op
    out = mats[0]
    for m in mats[1:]:
        out = np.kron(out, m)
    return out


_EPS = {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}


def explicit(ks):
    """Everything about Q for the four spins ks (k = 2j) from explicit matrices of the tensor product (dimension prod (k+1)).

    Returns a dict with the leg operators J, total spin Jt, Q on the whole space, an orthonormal basis V of the invariant subspace
    (the null space of the total J^2), and Q restricted to it."""
    J = [[_embed(ks, leg, _jm(ks[leg])[a]) for a in range(3)] for leg in range(4)]
    Jt = [sum(J[leg][a] for leg in range(4)) for a in range(3)]
    J2 = sum(x @ x for x in Jt)
    Q = sum(s * J[0][a] @ J[1][b] @ J[2][c] for (a, b, c), s in _EPS.items())
    w, v = np.linalg.eigh(J2)
    V = v[:, np.abs(w) < 1e-9]
    return {"J": J, "Jt": Jt, "Q": Q, "V": V, "Qr": V.conj().T @ Q @ V}


def coupled(ks):
    """Exact Q on the invariant subspace in the basis |(12)K, (34)K; 0>: Q = i [J1.J2, J2.J3], where J1.J2 is diagonal and J2.J3 is
    diagonal in the other channel (23)L, (14)L; the two channels are related by sqrt((2K+1)(2L+1)) {j1 j2 K; j3 j4 L}, exactly.
    Returns (the orthogonal recoupling matrix U, the Hermitian sympy matrix Q)."""
    k1, k2, k3, k4 = ks
    j = [R(k, 2) for k in ks]
    Ks = list(range(max(abs(k1 - k2), abs(k3 - k4)), min(k1 + k2, k3 + k4) + 1, 2))
    Ls = list(range(max(abs(k2 - k3), abs(k1 - k4)), min(k2 + k3, k1 + k4) + 1, 2))
    assert len(Ks) == len(Ls)
    x = [(R(K, 2) * (R(K, 2) + 1) - j[0] * (j[0] + 1) - j[1] * (j[1] + 1)) / 2 for K in Ks]
    y = [(R(L, 2) * (R(L, 2) + 1) - j[1] * (j[1] + 1) - j[2] * (j[2] + 1)) / 2 for L in Ls]
    U = sp.Matrix(len(Ks), len(Ls), lambda a, b: sp.sqrt((Ks[a] + 1) * (Ls[b] + 1)) *
                  wigner_6j(j[0], j[1], R(Ks[a], 2), j[2], j[3], R(Ls[b], 2)))
    Yk = U * sp.diag(*y) * U.T
    Qm = sp.Matrix(len(Ks), len(Ks), lambda a, b: sp.I * (x[a] - x[b]) * Yk[a, b])
    return U, Qm.applyfunc(sp.simplify)


def exact_spectrum(ks):
    """Eigenvalues of Q as exact radicals (sorted), from the characteristic polynomial of the coupled-basis matrix."""
    X = sp.Symbol("x")
    p = sp.factor(sp.expand(coupled(ks)[1].charpoly(X).as_expr()))
    roots = []
    for fac, mult in sp.factor_list(p)[1]:
        roots += [r for r in sp.solve(fac, X)] * mult
    return sorted((sp.nsimplify(sp.simplify(r)) for r in roots), key=lambda e: float(e)), p


def _lattice_degrees(plaquettes):
    """Vertex degrees of the lattice made of unit squares with lower-left corners at the given points."""
    links = set()
    for (x, y) in plaquettes:
        c = [(x, y), (x + 1, y), (x + 1, y + 1), (x, y + 1)]
        for i in range(4):
            links.add(frozenset((c[i], c[(i + 1) % 4])))
    deg = {}
    for l in links:
        for v in l:
            deg[v] = deg.get(v, 0) + 1
    return deg


def build_and_verify():
    # =====================================================================================================
    # Q on the invariant subspace
    # =====================================================================================================
    eq_q = register(NAME, Equation(
        key="qg_volume_triple_product", entry=ENTRY, title="Triple-product (volume) operator on a four-valent intertwiner space",
        source_kind="derived", source="derived here; the quantum tetrahedron is Barbieri, Nucl. Phys. B 518 (1998), in the knowledge base",
        latex_str=(r"Q=\epsilon_{abc}J_{1a}J_{2b}J_{3c}=i\,[\,J_1\!\cdot\!J_2,\ J_2\!\cdot\!J_3\,]\ \ \text{on}\ \ "
                   r"\mathrm{Inv}\bigl(V_{j_1}\!\otimes V_{j_2}\!\otimes V_{j_3}\!\otimes V_{j_4}\bigr),\qquad "
                   r"[Q,J_{\mathrm{tot}}]=0,\ \ Q^\dagger=Q,\ \ \mathrm{spec}\,Q=-\mathrm{spec}\,Q"
                   r"%%BR%%\dim\mathrm{Inv}=\#\{K:\ \triangle(j_1,j_2,K),\ \triangle(K,j_3,j_4)\},\qquad V\propto\sqrt{|Q|}"),
        notes="Four legs, one intertwiner for each intermediate spin K. The operator is real up to a factor i in the coupled basis and "
              "anticommutes with time reversal, so its spectrum is symmetric. It vanishes if any leg has spin 0 (then the vertex is "
              "trivalent). The constant relating the volume to sqrt|Q| is not fixed here: it depends on the volume operator chosen.",
    ))
    sets = [(1, 1, 1, 1), (2, 2, 2, 2), (2, 2, 1, 1), (1, 1, 2, 2), (2, 2, 3, 3), (3, 3, 3, 3), (4, 4, 4, 4), (1, 2, 3, 2)]
    info, dev_h, dev_c, dev_p, dev_s = {}, 0.0, 0.0, 0.0, 0.0
    dims_ok = True
    for ks in sets:
        e = explicit(ks)
        info[ks] = e
        Q, Jt, V, Qr = e["Q"], e["Jt"], e["V"], e["Qr"]
        dev_h = max(dev_h, float(np.abs(Q - Q.conj().T).max()))
        dev_c = max(dev_c, max(float(np.abs(Q @ Jt[a] - Jt[a] @ Q).max()) for a in range(3)))
        dev_p = max(dev_p, float(np.abs(Q @ V - V @ (V.conj().T @ Q @ V)).max()))
        ev = np.linalg.eigvalsh((Qr + Qr.conj().T) / 2)
        dev_s = max(dev_s, float(np.abs(np.sort(ev) + np.sort(ev)[::-1]).max()))
        dims_ok &= V.shape[1] == qg.invariants(ks) == qg.n_triangle(*ks)
    check(NAME, f"for eight sets of four spins, (1/2)^4, (1)^4, (1,1,1/2,1/2), (1/2,1/2,1,1), (1,1,3/2,3/2), (3/2)^4, (2)^4 and the mixed set (1/2,1,3/2,1), the explicit matrix of Q is "
                f"Hermitian (largest entry of Q - Q^dagger {dev_h:.0e}), commutes with the total spin ({dev_c:.0e}), maps the invariant "
                f"subspace to itself ({dev_p:.0e}), has the dimension of Weyl's integration formula and of the triangle count in every case, "
                f"and a spectrum symmetric under Q -> -Q ({dev_s:.0e})",
          dev_h == 0 and dev_c < 1e-10 and dev_p < 1e-10 and dev_s < 1e-10 and dims_ok, eq=eq_q)
    ident = 0.0
    for ks in [(1, 1, 1, 1), (2, 2, 2, 2), (1, 2, 3, 2)]:
        e = info[ks]
        J = e["J"]
        X = sum(J[0][a] @ J[1][a] for a in range(3))
        Y = sum(J[1][a] @ J[2][a] for a in range(3))
        ident = max(ident, float(np.abs(e["Q"] - 1j * (X @ Y - Y @ X)).max()))
    check(NAME, f"the operator identity Q = i [J1.J2, J2.J3] holds on the whole tensor product for (1/2)^4, (1)^4 and (1/2, 1, 3/2, 1) "
                f"(largest entry of the difference {ident:.1e}); it is what makes Q diagonalisable in terms of 6j symbols below",
          ident < 1e-12, eq=eq_q)
    zero = []
    for ks in [(1, 1, 2, 0), (2, 2, 2, 0), (1, 2, 3, 0), (2, 2, 4, 0)]:
        e = explicit(ks)
        zero.append((ks, e["V"].shape[1], float(np.abs(e["Qr"]).max())))
    check(NAME, "Q = 0 on the invariant space when the fourth leg has spin 0 (the vertex is then trivalent), for "
          + ", ".join(f"{ks}: dimension {d}, largest entry {m:.0e}" for ks, d, m in zero),
          all(d == 1 and m < 1e-12 for _, d, m in zero), eq=eq_q)
    # =====================================================================================================
    # exact spectra
    # =====================================================================================================
    eq_sp = register(NAME, Equation(
        key="qg_volume_spectra", entry=ENTRY, title="Exact spectra of the triple product for four equal and mixed spins",
        source_kind="derived", source="derived here from qg_volume_triple_product (squared 6j recoupling); numerical cross-check",
        latex_str=(r"\mathrm{spec}\,Q:\ \ (\tfrac12)^4:\ \pm\tfrac{\sqrt3}{4};\quad (1)^4:\ 0,\pm\sqrt3;\quad (\tfrac32)^4:\ \pm\tfrac{3\sqrt3}{4},\ "
                   r"\pm\tfrac{3\sqrt{35}}{4};\quad (1,1,\tfrac12,\tfrac12):\ \pm\tfrac{1}{\sqrt2}"
                   r"%%BR%%\mathrm{Tr}\,Q=\mathrm{Tr}\,Q^3=0,\qquad \mathrm{Tr}\,Q^2=\tfrac38,\ 6,\ \tfrac{171}{4},\ 198\ \ \bigl(j=\tfrac12,1,\tfrac32,2\bigr)"
                   r"%%BR%%Q_{KK'}=i\,(x_K-x_{K'})\sum_L U_{KL}\,y_L\,U_{K'L},\ \ U_{KL}=\sqrt{(2K+1)(2L+1)}\begin{Bmatrix}j_1&j_2&K\\j_3&j_4&L\end{Bmatrix}"),
        notes="x_K = J1.J2 = [K(K+1) - j1(j1+1) - j2(j2+1)]/2 and y_L = J2.J3 likewise with the channel (23)(14). U is orthogonal "
              "(checked) and no sign convention enters: a phase that depends on K or on L alone does not change the spectrum. "
              "The spectrum for four spin-1 is {0, +-sqrt3}: a state of zero volume exists when the number of states is odd.",
    ))
    worst, orth_ok = 0.0, True
    exact = {}
    for ks in sets:
        U, Qm = coupled(ks)
        orth_ok &= (U * U.T - sp.eye(U.shape[0])).applyfunc(sp.simplify) == sp.zeros(*U.shape)
        ev_exact = np.sort(np.linalg.eigvalsh(np.array(Qm.evalf(30), dtype=complex)))
        ev_num = np.sort(np.linalg.eigvalsh((info[ks]["Qr"] + info[ks]["Qr"].conj().T) / 2))
        worst = max(worst, float(np.abs(ev_exact - ev_num).max()))
        exact[ks] = Qm
    check(NAME, f"the coupled-basis form of Q with squared-6j recoupling (orthogonal U, checked exactly) reproduces the spectrum of the explicit "
                f"matrices for all eight sets of spins, largest difference {worst:.1e}", orth_ok and worst < 1e-10, eq=eq_sp)
    spec = {ks: exact_spectrum(ks)[0] for ks in [(1, 1, 1, 1), (2, 2, 2, 2), (3, 3, 3, 3), (2, 2, 1, 1), (2, 2, 3, 3)]}
    want = {(1, 1, 1, 1): [-sp.sqrt(3) / 4, sp.sqrt(3) / 4], (2, 2, 2, 2): [-sp.sqrt(3), 0, sp.sqrt(3)],
            (3, 3, 3, 3): sorted([-3 * sp.sqrt(35) / 4, -3 * sp.sqrt(3) / 4, 3 * sp.sqrt(3) / 4, 3 * sp.sqrt(35) / 4], key=float),
            (2, 2, 1, 1): [-1 / sp.sqrt(2), 1 / sp.sqrt(2)], (2, 2, 3, 3): [-sp.sqrt(26) / 2, 0, sp.sqrt(26) / 2]}
    ok = all(len(spec[ks]) == len(want[ks]) and all(sp.simplify(a - b) == 0 for a, b in zip(spec[ks], want[ks])) for ks in want)
    check(NAME, "exact eigenvalues from the characteristic polynomials: (1/2)^4: +-sqrt3/4; (1)^4: 0, +-sqrt3; (3/2)^4: +-3 sqrt3/4, "
                "+-3 sqrt35/4; (1,1,1/2,1/2): +-1/sqrt2; (1,1,3/2,3/2): 0, +-sqrt26/2 (the four-spin-1/2 value agrees with qg_quantum_tetrahedron)",
          ok, eq=eq_sp)
    tr = {ks: (sp.simplify(exact[ks].trace()), sp.simplify((exact[ks] ** 3).trace()), sp.simplify((exact[ks] ** 2).trace()))
          for ks in [(1, 1, 1, 1), (2, 2, 2, 2), (3, 3, 3, 3), (4, 4, 4, 4)]}
    tr_num = {ks: float(np.real(np.trace(info[ks]["Qr"] @ info[ks]["Qr"]))) for ks in tr}
    check(NAME, f"Tr Q = Tr Q^3 = 0 and Tr Q^2 = " + ", ".join(f"{tr[ks][2]}" for ks in tr) + " for four equal spins j = 1/2, 1, 3/2, 2, exactly; "
                f"the explicit matrices give " + ", ".join(f"{tr_num[ks]:.4f}" for ks in tr),
          all(t[0] == 0 and t[1] == 0 for t in tr.values()) and [tr[k][2] for k in tr] == [R(3, 8), 6, R(171, 4), 198]
          and all(abs(float(tr[ks][2]) - tr_num[ks]) < 1e-9 for ks in tr), eq=eq_sp)

    # =====================================================================================================
    # semiclassical trend
    # =====================================================================================================
    eq_sc = register(NAME, Equation(
        key="qg_volume_semiclassical", entry=ENTRY, title="Largest volume of four equal spins against the regular tetrahedron",
        source_kind="derived", source="derived here from qg_volume_spectra (numerical trend)",
        latex_str=(r"\frac{\lambda_{\max}(Q)}{\frac{4}{3\sqrt3}\,[j(j+1)]^{3/2}}=0.866,\ 0.796,\ 0.794,\ 0.808,\ 0.840,\ 0.867,\ 0.921"
                   r"\ \ (j=\tfrac12,1,\tfrac32,2,3,4,8)\ \ \text{increasing, below }1"),
        notes="Four unit vectors summing to zero with all mutual angles equal (the regular tetrahedron) have triple product "
              "4/(3 sqrt3) = 0.7698 for three of them; with vectors of length sqrt(j(j+1)) it is 0.7698 [j(j+1)]^{3/2}. "
              "This is a trend that the ratio rises toward 1 at about 1/j, not a proof of the semiclassical limit.",
    ))
    cl = 4 / (3 * np.sqrt(3))
    ratios = {}
    for k in (1, 2, 3, 4, 6, 8, 12, 16):
        A = np.array(coupled((k, k, k, k))[1].evalf(), dtype=complex)
        j = k / 2
        ratios[k] = float(np.linalg.eigvalsh(A).max() / (cl * (j * (j + 1)) ** 1.5))
    seq = [ratios[k] for k in (3, 4, 6, 8, 12, 16)]
    check(NAME, "largest eigenvalue over (4/(3 sqrt3))[j(j+1)]^{3/2} for four equal spins: " +
          ", ".join(f"j={k / 2:g}: {ratios[k]:.4f}" for k in ratios) + "; it increases strictly from j = 3/2 to j = 8 and stays below 1",
          all(a < b for a, b in zip(seq, seq[1:])) and max(ratios.values()) < 1 and ratios[16] > 0.92, eq=eq_sc)

    # =====================================================================================================
    # the L shape
    # =====================================================================================================
    eq_l = register(NAME, Equation(
        key="qg_volume_l_shape", entry=ENTRY, title="The four-valent vertex of the L-shaped cluster at strong coupling",
        source_kind="derived", source="derived here from qg_volume_triple_product and su2_strip_basis (leading order)",
        latex_str=(r"\langle Q^2\rangle_0=\frac{3}{16}\,\bigl|\psi_{bc}\bigr|^2=\frac{3}{16}\Bigl(\frac{4}{9\,g^{8}}\Bigr)^{2}=\frac{1}{27\,g^{16}}\,"
                   r"\bigl(1+O(g^{-8})\bigr)\ \ \text{(L shape, strong coupling)}"
                   r"%%BR%%\psi_{bc}=2\,\frac{(1/g^{2})^{2}}{(3g^{2}/2)(3g^{2})}=\frac{4}{9g^{8}},\quad E_{b}=\tfrac32g^{2},\ E_{bc}=3g^{2}"),
        notes="Plaquettes a, b, c in an L: a at the origin, b to its right, c above. The vertex shared by all three has four links, "
              "ab, bt (top of b), ac and cr (right of c). Only strings around b and c light all four legs, with spin 1/2 on each; "
              "the two loops each contract a pair of legs, an invariant vector on which Q^2 = 3/16. Its amplitude in the ground state "
              "is the second-order value 4/(9 g^8), the same as for two disjoint strings, which the strip's ground state shows for "
              "plaquettes a and c (relative deviation about g^-8). The next contribution to <Q^2> is order g^-24 (three strings). "
              "Not derived from the L-shape Hamiltonian, which is not built here.",
    ))
    strip_deg = _lattice_degrees([(0, 0), (1, 0), (2, 0)])
    l_deg = _lattice_degrees([(0, 0), (1, 0), (0, 1)])
    block_deg = _lattice_degrees([(0, 0), (1, 0), (0, 1), (1, 1)])
    four = [v for v, d in l_deg.items() if d == 4]
    check(NAME, f"vertex degrees of the lattice: the three-plaquette strip has maximum degree {max(strip_deg.values())} (every vertex "
                f"trivalent or bivalent, so no four-valent intertwiner and no volume), the L shape has maximum degree {max(l_deg.values())} "
                f"at exactly {len(four)} vertex {four[0]}, the 2 x 2 block at {[v for v, d in block_deg.items() if d == 4]}",
          max(strip_deg.values()) == 3 and max(l_deg.values()) == 4 and four == [(1, 1)] and
          [v for v, d in block_deg.items() if d == 4] == [(1, 1)], eq=eq_l)
    loops = {"a": [(0, 0), (1, 0), (1, 1), (0, 1)], "b": [(1, 0), (2, 0), (2, 1), (1, 1)], "c": [(0, 1), (1, 1), (1, 2), (0, 2)]}
    corner = (1, 1)
    legs_at = {}
    for name, c in loops.items():
        legs_at[name] = {frozenset((c[i], c[(i + 1) % 4])) for i in range(4) if corner in (c[i], c[(i + 1) % 4])}
    lit = {}
    for r in (1, 2, 3):
        for sub in itertools.combinations("abc", r):
            lit[sub] = set().union(*(legs_at[p] for p in sub))
    minimal = [s for s, l in lit.items() if len(l) == 4 and not any(len(t) < len(s) and set(t) <= set(s) and len(lit[t]) == 4 for t in lit)]
    check(NAME, "legs of the corner vertex lit by plaquette strings: " + "; ".join(f"{''.join(s)}: {len(l)}" for s, l in lit.items()
                                                                                    if len(s) <= 2) + f"; the only minimal set of strings that lights all four is {{{', '.join(''.join(m) for m in minimal)}}}",
          minimal == [("b", "c")] and all(len(lit[(p,)]) == 2 for p in "abc"), eq=eq_l)
    # the vector of two strings through the vertex: singlet (ab, bt) x singlet (ac, cr) in (1/2)^4 with legs (ab, bt, ac, cr)
    e = info[(1, 1, 1, 1)]
    sing = np.array([0, 1, -1, 0]) / np.sqrt(2)
    v = np.kron(sing, sing).astype(complex)
    Q, Jt = e["Q"], e["Jt"]
    inv_dev = max(float(np.abs(Jt[a] @ v).max()) for a in range(3))
    q2 = float((v.conj() @ Q @ Q @ v).real)
    q1 = float((v.conj() @ Q @ v).real)
    q2_vec = float(np.abs(Q @ Q @ v - (3 / 16) * v).max())
    check(NAME, f"the vector that two strings through the vertex give, (singlet of the two legs of b) x (singlet of the two legs of c), is "
                f"invariant (J_tot v: {inv_dev:.0e}), has <Q> = {q1:.0e} and <Q^2> = {q2:.6f} = 3/16, and Q^2 v = (3/16) v exactly "
                f"(largest deviation {q2_vec:.0e}), as it must on the two-dimensional space where Q = +-sqrt3/4",
          inv_dev < 1e-14 and abs(q1) < 1e-14 and abs(q2 - 3 / 16) < 1e-14 and q2_vec < 1e-14, eq=eq_l)
    ratios2 = {}
    for g2 in (6.0, 8.0, 10.0):
        H, S, kk, _ = strip.hamiltonian(40, g2)
        n = H.shape[0]
        w, vv = sla.eigsh(H, k=1, which="SA", v0=np.ones(n) / np.sqrt(n), tol=1e-13)
        gs = vv[:, 0]
        key = lambda s: strip._keys(np.array([s]))[0]
        i0, iac = np.searchsorted(kk, key((0, 0, 0, 0, 0, 0))), np.searchsorted(kk, key((1, 1, 0, 0, 1, 1)))
        ratios2[g2] = float(gs[iac] / gs[i0]) / (4 / (9 * g2 ** 4))
    check(NAME, "ground-state amplitude of two disjoint strings (plaquettes a and c of the strip, state (1/2, 1/2, 0, 0, 1/2, 1/2)) "
                "over the vacuum, divided by 4/(9 g^8): " + ", ".join(f"g^2={g}: {r:.5f}" for g, r in ratios2.items()) +
                "; the deviation from 1 falls like g^-8",
          all(abs(r - 1) < 1e-3 for r in ratios2.values()) and abs(ratios2[8.0] - 1) < abs(ratios2[6.0] - 1)
          and abs(ratios2[10.0] - 1) < abs(ratios2[8.0] - 1), eq=eq_l)
    amp = R(4, 9)
    check(NAME, f"<Q^2>_0 = (3/16) (4/9)^2 g^-16 = {R(3, 16) * amp ** 2} g^-16, exactly 1/27, and the second-order amplitude 2 (1/g^2)^2 / "
                f"((3g^2/2)(3g^2)) equals 4/(9 g^8) as a symbolic identity; three-string states (with a) enter at amplitude g^-12 and "
                f"probability g^-24",
          R(3, 16) * amp ** 2 == R(1, 27) and sp.simplify(2 * sp.Symbol("t", positive=True) ** -2 /
                                                          ((3 * sp.Symbol("t", positive=True) / 2) * (3 * sp.Symbol("t", positive=True)))
                                                          - R(4, 9) / sp.Symbol("t", positive=True) ** 4) == 0, eq=eq_l)
    return {"spectra": spec, "ratios": ratios, "amp_ratio": ratios2}


if __name__ == "__main__":
    out = build_and_verify()
    print("exact spectra:", out["spectra"])
    print("semiclassical ratios:", {k: round(v, 4) for k, v in out["ratios"].items()})
    print("two-string amplitude / (4/(9 g^8)):", out["amp_ratio"])
