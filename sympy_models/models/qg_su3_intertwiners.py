"""Intertwiner spaces for SU(3): the counting step of the plan to put more of the standard-model gauge group on the spin networks.

The SU(2) ladder of qg_spin_network.py counts the invariants of a vertex (the intertwiners) by Weyl's integration formula. This
module does the same for SU(3). It is a *kinematical, group-theoretic* model: it counts the SU(3)-invariant vectors in a tensor
product of irreducible representations, which is the number of vertex labels of a spin network (a colour-flux network) with
gauge group SU(3) on a graph with that vertex. It says nothing about dynamics, nothing about gravity, nothing about matter
fields or chirality, and nothing about the standard model beyond the group SU(3) itself: the electroweak SU(2) x U(1) and the
matter content are not touched. In particular SU(3) is not the group of loop quantum gravity, so the area spectrum of the SU(2)
model has no counterpart here.

What is registered, with its checks:

  qg_su3_intertwiners   dim Inv(V_1 x ... x V_n) = Haar integral of the product of characters = CT[prod chi_i * prod_{alpha>0}(1 - e^-alpha)],
                        with characters from Gelfand-Tsetlin patterns. Checked three independent ways: against numerical Haar
                        quadrature on the maximal torus with characters from the Weyl character formula (alternants); against
                        Schur-Weyl duality, which gives the invariants of 3^n x 3bar^n as n! minus the squares of the dimensions of
                        the symmetric-group irreps with more than three rows; and against the hook length formula, which gives the
                        invariants of 3^(3k) as the number of standard Young tableaux of shape (k,k,k) (1, 5, 42, 462, 6006).
                        The structural contrasts with SU(2): a trivalent vertex does not have a unique intertwiner (8 x 8 x 8 has
                        two), the vertex count vanishes unless the triality of the legs adds to zero, and the Casimir of a link is
                        checked against the Dynkin index of the weights.

Not tested: any Hamiltonian, ground state or spectrum; a volume or area operator for SU(3); a coupling to matter; the group
SU(3) x SU(2) x U(1); chirality or anomaly cancellation; gravity. Nothing here is a model of the standard model or of quantum gravity.
"""
import itertools
import sys
from fractions import Fraction as Fr
from functools import lru_cache
from math import factorial
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register

NAME = "qg_su3_intertwiners"
ENTRY = "spin-network-qg"


# ---------------------------------------------------------------------------------------------------------------
# irreps (p, q): Dynkin labels, p boxes of the fundamental and q of the antifundamental, highest weight lambda = (p+q, q, 0).
# a weight (n1, n2, n3) is carried in the reduced coordinates (n1 - n3, n2 - n3), so e1 + e2 + e3 = 0.
# ---------------------------------------------------------------------------------------------------------------
def red(w):
    return (w[0] - w[2], w[1] - w[2])


@lru_cache(maxsize=None)
def char(p, q):
    """Character of the irrep (p, q) as {reduced weight: multiplicity}, from the Gelfand-Tsetlin patterns of shape (p+q, q, 0)."""
    l1, l2, l3 = p + q, q, 0
    out = {}
    for a1 in range(l2, l1 + 1):
        for a2 in range(l3, l2 + 1):
            for b in range(a2, a1 + 1):
                w = red((b, a1 + a2 - b, (l1 + l2 + l3) - (a1 + a2)))
                out[w] = out.get(w, 0) + 1
    return out


def dim_irrep(p, q):
    return (p + 1) * (q + 1) * (p + q + 2) // 2


def mul(a, b):
    out = {}
    for u, x in a.items():
        for v, y in b.items():
            k = (u[0] + v[0], u[1] + v[1])
            out[k] = out.get(k, 0) + x * y
    return out


# the positive roots in reduced coordinates: e1 - e2, e2 - e3, e1 - e3
_ROOTS = [(1, -1), (1, 2), (2, 1)]


def invariants(irreps):
    """Dimension of the SU(3)-invariant subspace of V_{(p1,q1)} x V_{(p2,q2)} x ...: the Haar integral of the product of characters.

    Weyl's integration formula for a Weyl-invariant class function f is  int f = CT[ f * prod_{alpha > 0} (1 - e^{-alpha}) ]
    (the same step as the SU(2) formula CT[f (1 - z^2)]); the constant term is the coefficient of the weight 0."""
    f = {(0, 0): 1}
    for (p, q) in irreps:
        f = mul(f, char(p, q))
    for r in _ROOTS:
        f = mul(f, {(0, 0): 1, (-r[0], -r[1]): -1})
    return f.get((0, 0), 0)


# ---------------------------------------------------------------------------------------------------------------
# independent route 1: numerical Haar quadrature on the maximal torus, characters from the Weyl character formula
# ---------------------------------------------------------------------------------------------------------------
def invariants_quad(irreps, M=48):
    """(1/6) int |Delta|^2 prod chi_i  on the torus by an M x M equispaced grid, offset off the zeros of the denominator.

    An equispaced grid integrates a trigonometric polynomial exactly when M exceeds its degree, so this is exact up to rounding."""
    # different offsets on the two axes, so that no grid point lies on a wall x_i = x_j (the denominator of the characters)
    th1 = 2 * np.pi * (np.arange(M) + 0.637) / M
    th2 = 2 * np.pi * (np.arange(M) + 0.891) / M
    t1, t2 = np.meshgrid(th1, th2, indexing="ij")
    x = [np.exp(1j * t1), np.exp(1j * t2), np.exp(-1j * (t1 + t2))]

    def alt(l):                                   # det[x_i^{l_j + 3 - j}], j = 1..3
        e = [l[0] + 2, l[1] + 1, l[2]]
        d = 0
        for pi in itertools.permutations(range(3)):
            sign = 1 if sum(1 for i in range(3) for j in range(i + 1, 3) if pi[i] > pi[j]) % 2 == 0 else -1
            d = d + sign * x[0] ** e[pi[0]] * x[1] ** e[pi[1]] * x[2] ** e[pi[2]]
        return d

    delta = alt((0, 0, 0))
    val = np.ones_like(delta)
    for (p, q) in irreps:
        val = val * alt((p + q, q, 0)) / delta
    val = val * np.abs(delta) ** 2 / 6
    return float(np.mean(val).real)


# ---------------------------------------------------------------------------------------------------------------
# independent route 2 and 3: Schur-Weyl duality and the hook length formula (symmetric group data, no SU(3) characters)
# ---------------------------------------------------------------------------------------------------------------
def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def f_hook(shape):
    """Number of standard Young tableaux of the shape, by the hook length formula."""
    n = sum(shape)
    conj = [sum(1 for r in shape if r > c) for c in range(shape[0])] if shape else []
    prod = 1
    for i, r in enumerate(shape):
        for j in range(r):
            prod *= (r - j - 1) + (conj[j] - i - 1) + 1
    return factorial(n) // prod


def schur_weyl_count(n):
    """dim of the SU(3)-invariants in V^n x Vbar^n: the sum over partitions of n with at most three rows of f^2."""
    return sum(f_hook(lam) ** 2 for lam in partitions(n) if len(lam) <= 3)


# ---------------------------------------------------------------------------------------------------------------
# Casimir and the Dynkin index from the weights
# ---------------------------------------------------------------------------------------------------------------
def casimir(p, q):
    """C2(p, q) = (p^2 + q^2 + pq)/3 + p + q, with C2(3) = 4/3 and C2(8) = 3."""
    return Fr(p * p + q * q + p * q, 3) + p + q


def dynkin_index(p, q):
    """T(R) = Tr_R(T3^2) with T3 = diag(1, -1, 0)/2, read off the weights: the weight (n1, n2, n3) has T3 eigenvalue (n1 - n2)/2.

    The weights are stored reduced, (n1 - n3, n2 - n3), and n1 - n2 is the difference of the two coordinates."""
    return sum(Fr(m) * Fr(w[0] - w[1], 2) ** 2 for w, m in char(p, q).items())


def build_and_verify():
    eq = register(NAME, Equation(
        key="qg_su3_intertwiners", entry=ENTRY, title="Intertwiner spaces for SU(3): the invariants of a vertex",
        source_kind="derived", source="derived here; Weyl's integration formula, Schur-Weyl duality and the hook length formula are textbook",
        latex_str=(r"\dim\mathrm{Inv}\Bigl(\bigotimes_{i=1}^{n}V_{(p_i,q_i)}\Bigr)=\int_{SU(3)}\prod_{i=1}^{n}\chi_{(p_i,q_i)}(U)\,dU"
                   r"=\bigl[e^{0}\bigr]\prod_{i}\chi_{(p_i,q_i)}\prod_{\alpha>0}\bigl(1-e^{-\alpha}\bigr)"
                   r"%%BR%%\dim\mathrm{Inv}\bigl(V^{\otimes n}\otimes\bar V^{\otimes n}\bigr)=\sum_{\lambda\vdash n,\ \ell(\lambda)\le3}(f^{\lambda})^{2},\qquad "
                   r"\dim\mathrm{Inv}\bigl(V^{\otimes3k}\bigr)=f^{(k,k,k)}:\ 1,\ 5,\ 42,\ 462,\ 6006,\ \dots"
                   r"%%BR%%\dim V_{(p,q)}=\tfrac12(p+1)(q+1)(p+q+2),\qquad C_2(p,q)=\tfrac13(p^2+q^2+pq)+p+q"),
        notes="The SU(3)-invariant subspace of the tensor product of the representations meeting at a vertex is the space of "
              "intertwiners of a spin network with gauge group SU(3), the vertex labels of a colour-flux network. It is a purely "
              "group-theoretic count and is evaluated here by Weyl's integration formula for SU(3), whose characters come from "
              "Gelfand-Tsetlin patterns, and checked by three routes that do not use those characters. Unlike SU(2), a trivalent vertex "
              "does not have a unique intertwiner (the adjoint three times has two), so a network on the theta graph needs vertex labels. "
              "SU(3) is not the group of loop quantum gravity; this model has no area spectrum and does not touch the matter content "
              "or the electroweak factor of the standard model.",
    ))

    # --- the characters
    bad = [(p, q) for p in range(7) for q in range(7) if sum(char(p, q).values()) != dim_irrep(p, q)]
    check(NAME, "the Gelfand-Tsetlin character of every irrep (p, q) with p, q <= 6 has dimension (p+1)(q+1)(p+q+2)/2 "
                "(49 irreps)", not bad, eq=eq)
    conj_ok = all({(-w[0], -w[1]): m for w, m in char(q, p).items()} == char(p, q) for p in range(6) for q in range(6))
    check(NAME, "the character of (p, q) is the character of (q, p) with every weight negated (the dual representation), for p, q <= 5",
          conj_ok, eq=eq)

    # --- small cases with known answers
    fund, anti, adj, sext = (1, 0), (0, 1), (1, 1), (2, 0)
    known = [([fund, anti], 1), ([fund, fund, fund], 1), ([adj, adj], 1), ([adj, adj, adj], 2), ([fund, anti, adj], 1),
             ([fund, fund], 0), ([adj], 0), ([fund, fund, anti, anti], 2), ([sext, (0, 2)], 1), ([fund, fund, sext], 0)]
    check(NAME, "known invariants: 3 x 3bar: 1; 3 x 3 x 3: 1; 8 x 8: 1; 8 x 8 x 8: 2; 3 x 3bar x 8: 1; 3 x 3: 0; 8: 0; "
                "3 x 3 x 3bar x 3bar: 2; 6 x 6bar: 1; 3 x 3 x 6: 0",
          all(invariants(ir) == v for ir, v in known), eq=eq)

    # --- independent route 1: numerical Haar quadrature
    cases = [[fund, anti], [fund, fund, fund], [adj, adj], [adj, adj, adj], [fund, anti, adj], [sext, (0, 2)], [(2, 1), (1, 2)],
             [adj, adj, adj, adj], [fund, fund, anti, anti], [(2, 2), adj, (1, 2)], [fund, fund, fund, adj, adj]]
    worst = max(abs(invariants(c) - invariants_quad(c)) for c in cases)
    check(NAME, f"Weyl's formula agrees with numerical Haar quadrature on the torus (characters from the Weyl character formula, "
                f"a different construction) on {len(cases)} tensor products; largest difference {worst:.1e}", worst < 1e-8, eq=eq)

    # --- independent route 2: Schur-Weyl duality
    sw = [(n, invariants([fund] * n + [anti] * n), schur_weyl_count(n)) for n in range(1, 6)]
    check(NAME, "Schur-Weyl: invariants of 3^n x 3bar^n equal the sum of f^2 over partitions of n with at most three rows, for n = 1..5: "
                + ", ".join(f"{a}" for _, a, _ in sw) + " (n! for n <= 3, then 23 and 103, not 24 and 120)",
          all(a == b for _, a, b in sw) and [a for _, a, _ in sw] == [1, 2, 6, 23, 103], eq=eq)

    # --- independent route 3: the hook length formula
    hk = [(k, invariants([fund] * (3 * k)), f_hook((k, k, k))) for k in range(1, 6)]
    check(NAME, "hook length formula: invariants of 3^(3k) equal the number of standard Young tableaux of shape (k,k,k), for k = 1..5: "
                + ", ".join(str(a) for _, a, _ in hk), all(a == b for _, a, b in hk) and [a for _, a, _ in hk] == [1, 5, 42, 462, 6006],
          eq=eq)
    check(NAME, "an invariant needs a total triality zero: no invariants in 3^n for n not divisible by 3, n = 1..14",
          all(invariants([fund] * n) == 0 for n in range(1, 15) if n % 3), eq=eq)

    # --- the contrast with SU(2): trivalent vertices
    triples = list(itertools.product([(p, q) for p in range(4) for q in range(4)], repeat=3))
    vals = {t: invariants(list(t)) for t in triples}
    mx = max(vals.values())
    first = min((t for t, v in vals.items() if v == 2), key=lambda t: sum(map(sum, t)))
    check(NAME, f"over all {len(triples)} triples of irreps with p, q <= 3 the number of invariants takes the values "
                f"{sorted(set(vals.values()))}, so a trivalent SU(3) vertex is not multiplicity free; the smallest triple with two is {first}",
          mx >= 2 and min(vals.values()) == 0 and first == (adj, adj, adj), eq=eq)
    perm = all(vals[t] == vals[(t[1], t[0], t[2])] == vals[(t[2], t[1], t[0])] for t in triples[::7])
    dual = all(vals[t] == vals[tuple((q, p) for (p, q) in t)] for t in triples)
    tri = all((vals[t] == 0) or (sum(p - q for (p, q) in t) % 3 == 0) for t in triples)
    check(NAME, "the count is symmetric under permuting the legs, under replacing every leg by its dual, and vanishes unless the "
                "triality of the three legs adds to zero mod 3 (all triples with p, q <= 3)", perm and dual and tri, eq=eq)

    # --- the link Casimir against the Dynkin index of the weights
    ok_c = all(8 * dynkin_index(p, q) == dim_irrep(p, q) * casimir(p, q) for p in range(7) for q in range(7))
    check(NAME, "the Casimir C2(p,q) = (p^2+q^2+pq)/3 + p + q satisfies dim(G) T(R) = dim(R) C2(R) with T(R) the Dynkin index read off the "
                "weights, for all 49 irreps with p, q <= 6; C2 = 4/3 for 3 and 3 for 8",
          ok_c and casimir(1, 0) == Fr(4, 3) and casimir(1, 1) == 3 and casimir(1, 0) == casimir(0, 1), eq=eq)
    return {"adjoint^3": invariants([adj, adj, adj]), "hook": [a for _, a, _ in hk]}


if __name__ == "__main__":
    print(build_and_verify())
