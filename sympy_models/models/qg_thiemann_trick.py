"""The identities behind Thiemann's regularised Hamiltonian constraint (companion to qg_volume.py).

THIS MODULE DOES NOT CONSTRUCT THE QUANTUM CONSTRAINT OPERATOR. It checks, exactly, the classical steps by which Thiemann
rewrites the Euclidean constraint as a holonomy-and-volume expression that can be quantised on spin networks (QSD, gr-qc/9606089),
and the one quantum property of the result that follows from the volume operator alone.

With a triad E^a_i, V = sqrt|det E| (one cell), and a Poisson bracket {A_c^k, E^d_l} = c delta^d_c delta^k_l:

 (1) the cofactor identity   eps^{kij} E^a_i E^b_j / sqrt|det E|  =  2 eps^{abc} dV/dE^c_k ,
     so that the Euclidean constraint density eps^{ijk} E^a_i E^b_j F^k_ab / sqrt|det E| equals
     2 eps^{abc} F^k_ab dV/dE^c_k = (2/c) eps^{abc} F^k_ab {A_c^k, V}:  the inverse square root of the determinant, which has no
     quantisation, is traded for a Poisson bracket with the volume, which does;
 (2) the same bracket, computed on the phase space of (A, E), gives {A_c^k, V} = c dV/dE^c_k;
 (3) the regularisation: with h_alpha = exp(eps^2 F_ab^k tau_k) (a small loop) and h_s = exp(eps a^k tau_k) (a short segment),
     Tr( h_alpha h_s {h_s^{-1}, V} ) = (eps^3 / 2) F_ab^k {A_c^k, V} + O(eps^4), for tau_k = -i sigma_k/2, so the trace of a
     small loop holonomy times the bracket of a segment with the volume reproduces the integrand of (1) up to the factor
     eps^3/2 that the regularisation absorbs; the leading coefficient does not depend on the base holonomy a.
 (4) a quantum property: the triple product Q = eps_abc J1_a J2_b J3_c vanishes identically on the invariant subspace of a vertex with
     three legs (J3 = -J1 - J2), for every set of spins tested. The volume is therefore zero at trivalent vertices, and so is
     the commutator with the volume that appears in the regularised constraint: a Thiemann-type constraint cannot act
     at a trivalent vertex (it acts only where V is not zero, first at four-valent vertices). This is a statement about V, in the
     convention of qg_volume.py, not about a constructed constraint.

Not done here, and the reason it matters. The quantum operator replaces {h, V} by (1/i hbar)[h, V] and acts by inserting new
edges and vertices: it changes the graph. Its matrix elements need 6j recoupling at every new vertex, with phases that this
module does not guess (qg_volume.py declines for the same reason). Whether the algebra of such operators closes, the
question raised by Nicolai, Peeters and Zamaklar, cannot be tested on a fixed graph. No claim is made here about anomaly
freedom, the Master Constraint, a semiclassical limit, or the relation of this construction to the ladder Hamiltonian.
"""
import itertools
import sys
from pathlib import Path

import numpy as np
import sympy as sp

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register

NAME = "qg_thiemann_trick"
ENTRY = "spin-network-qg"


def eps3():
    return sp.LeviCivita


def triad_symbols():
    return sp.Matrix(3, 3, lambda a, i: sp.Symbol(f"E{a}{i}"))


def pauli_tau():
    s1 = sp.Matrix([[0, 1], [1, 0]])
    s2 = sp.Matrix([[0, -sp.I], [sp.I, 0]])
    s3 = sp.Matrix([[1, 0], [0, -1]])
    return [-sp.I * s / 2 for s in (s1, s2, s3)]


def mexp_series(M, eps, order):
    out = sp.eye(2)
    term = sp.eye(2)
    for n in range(1, order + 1):
        term = term * M / n
        out = out + term
    return out


def spin_matrices(k):
    j = sp.Rational(k, 2)
    dim = k + 1
    ms = [j - i for i in range(dim)]
    Jz = sp.diag(*ms)
    Jp = sp.zeros(dim, dim)
    for col in range(1, dim):
        m = ms[col]
        Jp[col - 1, col] = sp.sqrt(j * (j + 1) - m * (m + 1))
    Jx = (Jp + Jp.T) / 2
    Jy = (Jp - Jp.T) / (2 * sp.I)
    return Jx, Jy, Jz


def kron_all(mats):
    out = mats[0]
    for m in mats[1:]:
        out = sp.kronecker_product(out, m)
    return out


def triple_trivalent(ks):
    """Q = eps_abc J1_a J2_b J3_c on the invariant subspace of the product of three spins k/2 (numpy)."""
    dims = [k + 1 for k in ks]
    ops = []
    for leg in range(3):
        mats = []
        for a in range(3):
            parts = [sp.eye(d) for d in dims]
            parts[leg] = spin_matrices(ks[leg])[a]
            mats.append(np.array(kron_all(parts).evalf(), dtype=complex))
        ops.append(mats)
    Q = sum(float(sp.LeviCivita(a, b, c)) * ops[0][a] @ ops[1][b] @ ops[2][c]
            for a in range(3) for b in range(3) for c in range(3) if len({a, b, c}) == 3)
    Jt = [ops[0][a] + ops[1][a] + ops[2][a] for a in range(3)]
    C2 = sum(J @ J for J in Jt)
    w, v = np.linalg.eigh((C2 + C2.conj().T) / 2)
    inv = v[:, np.abs(w) < 1e-9]
    return Q, inv


def build_and_verify():
    eq = register(NAME, Equation(
        key="qg_thiemann_trick", entry=ENTRY, title="The classical identities behind Thiemann's regularisation of the Euclidean constraint",
        source_kind="derived", source="derived here; the construction is Thiemann, Quantum Spin Dynamics, gr-qc/9606089, whose quantum operator is not built here",
        latex_str=(r"\frac{\epsilon^{kij}E^{a}_{i}E^{b}_{j}}{\sqrt{|\det E|}}=2\,\epsilon^{abc}\frac{\partial V}{\partial E^{c}_{k}},\qquad V=\sqrt{|\det E|},"
                   r"\qquad \{A^{k}_{c},V\}=c\,\frac{\partial V}{\partial E^{c}_{k}}"
                   r"%%BR%%\frac{\epsilon^{ijk}E^{a}_{i}E^{b}_{j}F^{k}_{ab}}{\sqrt{|\det E|}}=\frac{2}{c}\,\epsilon^{abc}F^{k}_{ab}\{A^{k}_{c},V\}"
                   r"%%BR%%\mathrm{Tr}\bigl(h_{\alpha}h_{s}\{h_{s}^{-1},V\}\bigr)=\tfrac{\varepsilon^{3}}{2}F^{k}_{ab}\{A^{k}_{c},V\}+O(\varepsilon^{4}),\qquad "
                   r"\mathrm{Inv}(V_{j_1}\!\otimes V_{j_2}\!\otimes V_{j_3}):\ \ \epsilon_{abc}J_1^{a}J_2^{b}J_3^{c}=0"),
        notes="Classical steps of Thiemann's construction, checked exactly (cofactor identity, the bracket on the (A, E) phase space, the small-loop expansion in exact "
              "series) plus the vanishing of the triple product at a trivalent vertex. The quantum constraint operator, which changes the graph and needs 6j "
              "recoupling at the new vertices, is NOT constructed; whether the algebra of such operators closes is not tested.",
    ))
    E = triad_symbols()
    detE = E.det()
    V = sp.sqrt(detE)

    # (1) the cofactor identity, symbolically, by cofactors and the square root of the determinant squared
    ok = True
    for k in range(3):
        for a in range(3):
            for b in range(3):
                lhs2 = sum(sp.LeviCivita(k, i, j) * E[a, i] * E[b, j] for i in range(3) for j in range(3))
                rhs2 = sum(sp.LeviCivita(a, b, c) * sp.diff(detE, E[c, k]) for c in range(3))
                ok = ok and sp.expand(lhs2 - rhs2) == 0
    check(NAME, "cofactor identity eps^{kij} E^a_i E^b_j = eps^{abc} d(det E)/dE^c_k for all 27 index triples, as a polynomial identity in the nine entries of E; "
                "dividing by sqrt|det E| and using d sqrt(det E)/dE = (1/(2 sqrt(det E))) d(det E)/dE gives the first display", ok, eq=eq)

    rng = np.random.default_rng(3)
    worst = 0.0
    for _ in range(200):
        M = rng.normal(size=(3, 3))
        d = np.linalg.det(M)
        sgn = np.sign(d)
        cof = np.linalg.inv(M).T * d
        dV = sgn * cof / (2 * np.sqrt(abs(d)))
        lhs = np.einsum("kij,ai,bj->abk", np.array([[[sp.LeviCivita(k, i, j) for j in range(3)] for i in range(3)] for k in range(3)], dtype=float), M, M) / np.sqrt(abs(d))
        rhs = 2 * np.einsum("abc,ck->abk", np.array([[[sp.LeviCivita(a, b, c) for c in range(3)] for b in range(3)] for a in range(3)], dtype=float), dV) * sgn
        worst = max(worst, np.max(np.abs(lhs - rhs)))
    check(NAME, f"the same identity in floating point on 200 random triads of both orientations (sign of the determinant included): largest difference {worst:.1e}", worst < 1e-10, eq=eq)

    # (2) the Poisson bracket on the (A, E) phase space
    c = sp.Symbol("c", positive=True)
    A = sp.Matrix(3, 3, lambda c_, k: sp.Symbol(f"A{c_}{k}"))

    def pb(f, g):
        return sum(c * (sp.diff(f, A[a, k]) * sp.diff(g, E[a, k]) - sp.diff(g, A[a, k]) * sp.diff(f, E[a, k])) for a in range(3) for k in range(3))

    okb = all(sp.simplify(pb(A[cc, k], V) - c * sp.diff(V, E[cc, k])) == 0 for cc in range(3) for k in range(3))
    check(NAME, "with {A^k_c, E^d_l} = c delta^d_c delta^k_l the bracket of A with V = sqrt(det E) is c dV/dE for all nine components", okb, eq=eq)

    # (3) the small loop expansion, exact series
    eps = sp.Symbol("eps")
    tau = pauli_tau()
    rows = []
    allok = True
    for a_vec, F_vec, w_vec in [((1, 0, 0), (0, 1, 0), (0, 0, 1)), ((2, -1, 3), (1, 2, -1), (-1, 1, 2)), ((0, 0, 0), (1, 1, 1), (2, 0, 1)),
                                ((sp.Rational(1, 2), 3, -2), (-1, 0, 3), (1, 1, -1))]:
        s = sp.Symbol("s")
        avec = sp.Matrix(a_vec)
        wvec = sp.Matrix(w_vec)
        Fvec = sp.Matrix(F_vec)
        tau_dot = lambda v: sum((v[k] * tau[k] for k in range(3)), sp.zeros(2, 2))
        hs = mexp_series(eps * tau_dot(avec), eps, 6)
        hs_inv_s = mexp_series(-eps * tau_dot(avec + s * wvec), eps, 6)
        brk = sp.diff(hs_inv_s, s).subs(s, 0)          # {h_s^{-1}, V} with {A, V} = w
        ha = mexp_series(eps**2 * tau_dot(Fvec), eps, 3)
        tr = sp.expand((ha * hs * brk).trace())
        co = sp.series(tr, eps, 0, 5).removeO()
        c3 = sp.simplify(co.coeff(eps, 3))
        c1 = sp.simplify(co.coeff(eps, 1))
        c2 = sp.simplify(co.coeff(eps, 2))
        target = sp.Rational(1, 2) * (Fvec.dot(wvec))
        allok = allok and sp.simplify(c3 - target) == 0 and c1 == 0 and c2 == 0
        rows.append(f"a={tuple(a_vec)}: coefficient {c3}")
    check(NAME, "exact series: Tr(h_alpha h_s {h_s^{-1}, V}) has no eps^1 or eps^2 term and its eps^3 coefficient is (1/2) F . {A, V} for four choices of the base "
                "holonomy a, the loop F and the bracket {A, V}, so the leading term is independent of a (" + "; ".join(rows) + ")", allok, eq=eq)

    # (4) trivalent vertices have no volume
    spins = [(1, 1, 2), (1, 2, 3), (2, 2, 2), (2, 2, 4), (1, 3, 4), (3, 3, 2)]
    res = []
    for ks in spins:
        Q, inv = triple_trivalent(ks)
        dim = inv.shape[1]
        nrm = np.max(np.abs(Q @ inv)) if dim else 0.0
        res.append((ks, dim, nrm))
    check(NAME, "a vertex with three legs has no volume: on the invariant subspace of spins " + ", ".join(f"{[k / 2 for k in ks]}" for ks, _, _ in res)
          + f" (dimension one each, as the triangle rule says) the triple product is zero (largest entry {max(r[2] for r in res):.1e})",
          all(d == 1 for _, d, _ in res) and max(r[2] for r in res) < 1e-9, eq=eq)

    ks = (1, 1, 1)
    Q, inv = triple_trivalent((1, 1, 1))
    check(NAME, "the triangle rule leaves nothing to check for three half-integer spins 1/2, 1/2, 1/2 (the product has no invariant, dimension "
                f"{inv.shape[1]}, because the total spin is half-integer)", inv.shape[1] == 0, eq=eq)
    return eq


if __name__ == "__main__":
    build_and_verify()
    print("ok")
