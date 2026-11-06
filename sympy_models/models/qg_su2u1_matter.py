"""Matter charged under SU(2) x U(1) on the two-plaquette ladder graph (a toy for the weak-isospin / hypercharge factor).

What this is. The ladder graph of su2_lattice_ladder.py (7 links, 6 vertices) is given the gauge group SU(2)_w x U(1)_Y
of the electroweak sector, a group that is deliberately NOT the Ashtekar-Barbero SU(2) of the gravity side: the two
are different groups acting on different things. Static charged doublets (a Higgs doublet or a lepton doublet,
hypercharge in units of 1/2) are put at vertices. Gauss's law at every vertex then selects the gauge-invariant electric
flux configurations: an integer U(1) flux n_l and an SU(2) spin k_l/2 on every link, with the doublet counted as a source
of spin 1/2 and charge q.

The electric (strong-coupling) energy of a flux configuration is
    E = (g^2/2) sum_l (k_l/2)(k_l/2 + 1) + (g'^2/2) sum_l n_l^2 ,
so the lowest energy of a static pair is a minimum over Gauss-law-satisfying configurations. This is exact counting
on a finite graph: no dynamics, no plaquette term, no continuum limit.

Checked: SU(2) vertex invariants by two independent routes (iterated Clebsch-Gordan coupling and the Weyl integration
formula as a Laurent coefficient); with no sources the Gauss-law states reproduce the (j_a, j_b, J) basis of
su2_lattice_ladder.py with the same electric energies and the cycle space of the graph has dimension 2 for U(1); one source
(or two with the same charge) gives no invariant state, two opposite sources give the shortest-path flux; the minimum
energy and its degeneracy equal (3/4 g^2/2 + g'^2/2 for charge 1) times the graph distance, and the number of
shortest paths, by brute force over all configurations up to a cutoff. One generation of the standard model fields
(SU(3) x SU(2) x U(1) representations) is checked for the anomaly sums (Y, Y^3, SU(2)^2 Y, SU(3)^2 Y), the Witten
parity (an even number of doublets), and the Z6 that acts trivially on every field; the three Yukawa vertices
are counted as SU(3) invariants (from qg_su3_intertwiners) times SU(2) invariants times the hypercharge sum.

Not tested here: any dynamics (plaquette terms, matter hopping, a Higgs potential), fermion doubling or chirality on the
lattice, gravity-matter coupling (the Ashtekar-Barbero SU(2) does not appear), the Hamiltonian constraint, anomaly
cancellation beyond the sums listed, three generations, mixing angles. That SU(3) x SU(2) x U(1) vertex invariants exist
for the Yukawa couplings is group theory and says nothing about whether the coupled theory exists.
"""
import itertools
import sys
from fractions import Fraction as Fr
from functools import lru_cache
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register
from models.qg_su3_intertwiners import invariants as su3_invariants

NAME = "qg_su2u1_matter"
ENTRY = "spin-network-qg"

# the ladder: vertices and oriented links (the orientation only fixes the sign of the U(1) flux)
VERTS = ("s1", "s2", "a1", "a2", "b1", "b2")
LINKS = (("s1", "s2"), ("s1", "a1"), ("a1", "a2"), ("a2", "s2"), ("s1", "b1"), ("b1", "b2"), ("b2", "s2"))


# ---------------------------------------------------------------------------------------------------------------
# SU(2) vertex invariants, spins carried as k = 2j
# ---------------------------------------------------------------------------------------------------------------
def inv_coupling(ks):
    """Number of singlets in the tensor product of spins k/2, by iterated coupling."""
    mult = {0: 1}
    for k in ks:
        new = {}
        for t, m in mult.items():
            for r in range(abs(t - k), t + k + 1, 2):
                new[r] = new.get(r, 0) + m
        mult = new
    return mult.get(0, 0)


def inv_weyl(ks):
    """The same number as the coefficient of z^0 in prod chi_k(z) (1 - z^2), chi_k = z^k + z^(k-2) + ... + z^-k."""
    poly = {0: 1}
    for k in ks:
        new = {}
        for e, c in poly.items():
            for d in range(-k, k + 1, 2):
                new[e + d] = new.get(e + d, 0) + c
        poly = new
    return poly.get(0, 0) - poly.get(2, 0)


# ---------------------------------------------------------------------------------------------------------------
# Gauss-law configurations on the ladder
# ---------------------------------------------------------------------------------------------------------------
def at(v):
    return [i for i, (a, b) in enumerate(LINKS) if v in (a, b)]


def su2_configs(src, K):
    """All spin assignments k_l in 0..K whose vertex invariant counts are nonzero, with the sources (a dict vertex -> number of
    doublets) included; yields (ks, degeneracy)."""
    for ks in itertools.product(range(K + 1), repeat=len(LINKS)):
        d = 1
        for v in VERTS:
            legs = [ks[i] for i in at(v)] + [1] * src.get(v, 0)
            n = inv_coupling(legs)
            if n == 0:
                d = 0
                break
            d *= n
        if d:
            yield ks, d


def u1_configs(src, N):
    """All integer fluxes n_l in -N..N with out-flow minus in-flow equal to the source charge at every vertex."""
    for ns in itertools.product(range(-N, N + 1), repeat=len(LINKS)):
        ok = True
        for v in VERTS:
            div = sum(ns[i] if LINKS[i][0] == v else -ns[i] for i in at(v))
            if div != src.get(v, 0):
                ok = False
                break
        if ok:
            yield ns


def e_su2(ks):
    return sum(Fr(k * (k + 2), 4) for k in ks)


def e_u1(ns):
    return sum(n * n for n in ns)


def distance(u, v):
    dist, frontier = {u: 0}, [u]
    while frontier:
        nxt = []
        for x in frontier:
            for a, b in LINKS:
                for p, q in ((a, b), (b, a)):
                    if p == x and q not in dist:
                        dist[q] = dist[x] + 1
                        nxt.append(q)
        frontier = nxt
    return dist[v]


def count_shortest(u, v):
    """Number of shortest simple paths from u to v."""
    d = distance(u, v)
    def rec(x, seen, left):
        if x == v:
            return 1 if left == 0 else 0
        if left == 0:
            return 0
        tot = 0
        for a, b in LINKS:
            for p, q in ((a, b), (b, a)):
                if p == x and q not in seen:
                    tot += rec(q, seen | {q}, left - 1)
        return tot
    return rec(u, {u}, d)


# ---------------------------------------------------------------------------------------------------------------
# One generation of standard model fields, all left-handed Weyl: (name, SU(3) irrep (p, q), SU(2) k = 2j, hypercharge Y)
# ---------------------------------------------------------------------------------------------------------------
FIELDS = (("Q", (1, 0), 1, Fr(1, 6)), ("u^c", (0, 1), 0, Fr(-2, 3)), ("d^c", (0, 1), 0, Fr(1, 3)),
          ("L", (0, 0), 1, Fr(-1, 2)), ("e^c", (0, 0), 0, Fr(1)))
HIGGS = ("H", (0, 0), 1, Fr(1, 2))
HTILDE = ("Htilde", (0, 0), 1, Fr(-1, 2))     # i sigma_2 H^*: the same SU(2) representation, opposite hypercharge


def dim3(pq):
    p, q = pq
    return (p + 1) * (q + 1) * (p + q + 2) // 2


def yukawa(fields):
    return su3_invariants([f[1] for f in fields]) * inv_coupling([f[2] for f in fields]) * (1 if sum(f[3] for f in fields) == 0 else 0)


def build_and_verify():
    eq = register(NAME, Equation(
        key="qg_su2u1_matter", entry=ENTRY, title="Static SU(2) x U(1) charges on the ladder: Gauss's law and the shortest flux path",
        source_kind="derived", source="derived here; Gauss's law on a graph, Clebsch-Gordan coupling and the standard model charge table are textbook",
        latex_str=(r"\sum_{l\ni v}\pm n_l=q_v,\qquad \dim\mathrm{Inv}_v=\dim\mathrm{Inv}\Bigl(\bigotimes_{l\ni v}V_{k_l/2}\otimes V_{1/2}^{\otimes s_v}\Bigr)\ge1"
                   r"%%BR%%E=\tfrac{g^{2}}{2}\sum_l\tfrac{k_l}{2}\Bigl(\tfrac{k_l}{2}+1\Bigr)+\tfrac{g'^{2}}{2}\sum_l n_l^{2},\qquad "
                   r"E_{\min}=\Bigl(\tfrac{3}{8}g^{2}+\tfrac12g'^{2}\Bigr)\,d(u,v)\ \text{for a doublet of charge 1 at $u$ and $v$}"
                   r"%%BR%%\sum Y=\sum Y^{3}=0,\quad \tfrac{t}{3}+\tfrac{k}{2}+Y\in\mathbb{Z}\ \text{for every field}"),
        notes="The group of this model is SU(2)_w x U(1)_Y, the weak-isospin and hypercharge factors, NOT the Ashtekar-Barbero SU(2) "
              "of the gravity side. Static doublets at vertices are sources for both Gauss laws; in the ladder basis the "
              "electric flux of an opposite pair must run along a path, and the lowest energy is the graph distance times the "
              "energy per link. The numbers are exact finite counts over all configurations below a cutoff. Nothing here "
              "couples matter to the gravitational spin networks, defines a matter Hamiltonian or a Higgs potential, or says "
              "anything about chirality on a lattice.",
    ))

    # --- SU(2) invariants, two routes
    tup = [t for n in range(1, 6) for t in itertools.product(range(5), repeat=n)]
    check(NAME, f"SU(2) vertex invariants: iterated Clebsch-Gordan coupling equals the Weyl coefficient formula on all {len(tup)} "
                f"tuples of 1 to 5 spins k <= 4", all(inv_coupling(t) == inv_weyl(t) for t in tup), eq=eq)
    check(NAME, "known invariants: 1/2 x 1/2: 1; 1/2 x 1/2 x 1: 1; 1 x 1 x 1: 1; 1/2 x 1/2 x 1/2 x 1/2: 2; a single doublet: 0; "
                "1/2 x 1: 0", all(inv_coupling(t) == v for t, v in [((1, 1), 1), ((1, 1, 2), 1), ((2, 2, 2), 1),
                                                                   ((1, 1, 1, 1), 2), ((1,), 0), ((1, 2), 0)]), eq=eq)

    # --- no sources: the ladder basis
    cfgs = list(su2_configs({}, 4))
    nonzero_basis = {(ks[1], ks[4], ks[0]) for ks, _ in cfgs}
    basis = {(ja, jb, J) for ja in range(0, 5) for jb in range(0, 5) for J in range(0, 5)
             if abs(ja - jb) <= J <= ja + jb and (ja + jb + J) % 2 == 0}
    equal_sides = all(ks[1] == ks[2] == ks[3] and ks[4] == ks[5] == ks[6] for ks, _ in cfgs)
    check(NAME, f"no sources: the Gauss-law configurations with spins k <= 4 are exactly the (j_a, j_b, J) states of the ladder "
                f"basis (three equal outer links per plaquette, the shared link J, triangle rule): {len(cfgs)} states, "
                f"all with multiplicity one", equal_sides and {(ks[1], ks[4], ks[0]) for ks, _ in cfgs} == {
                    (a, b, J) for (a, b, J) in basis if all(x <= 4 for x in (a, b, J))} and all(d == 1 for _, d in cfgs), eq=eq)
    en = all(e_su2(ks) == Fr(3 * ks[1] * (ks[1] + 2), 4) + Fr(3 * ks[4] * (ks[4] + 2), 4) + Fr(ks[0] * (ks[0] + 2), 4)
             for ks, _ in cfgs)
    check(NAME, "their electric energies are 3 j_a(j_a+1) + 3 j_b(j_b+1) + J(J+1) in units of g^2/2, as in su2_lattice_ladder.py", en, eq=eq)

    flows = list(u1_configs({}, 2))
    na_nb = {(ns[1], ns[4], ns[0]) for ns in flows}
    free = all(ns[0] == -(ns[1] + ns[4]) and ns[1] == ns[2] == ns[3] and ns[4] == ns[5] == ns[6] for ns in flows)
    check(NAME, f"U(1) with no sources: the divergence-free fluxes with |n| <= 2 are the cycle space of the graph, which has dimension "
                f"7 - 6 + 1 = 2 (two plaquette fluxes; {len(flows)} configurations, the shared link carries minus the sum of the two outer fluxes when both are "
                f"oriented away from s1)", len(flows) == len({(ns[1], ns[4]) for ns in flows}) and free and len(LINKS) - len(VERTS) + 1 == 2, eq=eq)

    # --- selection rules for the sources
    one = {"a1": 1}
    two_same = {"a1": 1, "b2": 1}
    check(NAME, "a single doublet cannot be neutralised: no SU(2) Gauss-law state (the centre of SU(2) acts as -1) and no U(1) flux "
                "state (a net charge on a finite graph), for every flux k <= 3, |n| <= 2",
          not list(su2_configs(one, 3)) and not list(u1_configs({"a1": 1}, 2)), eq=eq)
    check(NAME, "two doublets of the SAME charge: SU(2) alone allows a neutral state, U(1) does not, so it is the U(1) charge that "
                "separates a particle from its antiparticle (the doublet is pseudoreal)",
          bool(list(su2_configs(two_same, 3))) and not list(u1_configs(two_same, 2)), eq=eq)

    # --- opposite pair: shortest path
    rows = []
    allok = True
    for u, v in [("a1", "b2"), ("a1", "a2"), ("s1", "s2"), ("a1", "b1"), ("a2", "b1"), ("s1", "a2")]:
        d = distance(u, v)
        paths = count_shortest(u, v)
        su = list(su2_configs({u: 1, v: 1}, 3))
        e2 = min(e_su2(ks) for ks, _ in su)
        deg2 = sum(dg for ks, dg in su if e_su2(ks) == e2)
        fl = list(u1_configs({u: 1, v: -1}, 2))
        e1 = min(e_u1(ns) for ns in fl)
        deg1 = sum(1 for ns in fl if e_u1(ns) == e1)
        ok = e2 == Fr(3, 4) * d and e1 == d and deg2 == paths and deg1 == paths
        allok = allok and ok
        rows.append(f"{u}-{v}: d={d}, {paths} path(s)")
    check(NAME, "an opposite pair of charge-1 doublets: over all configurations (k <= 3, |n| <= 2) the lowest SU(2) flux energy is 3/4 per link "
                "and the lowest U(1) flux energy 1 per link times the graph distance, each with degeneracy equal to the number of "
                "shortest paths (" + "; ".join(rows) + ")", allok, eq=eq)

    # --- the standard model fields
    allf = FIELDS
    sy = sum(f[3] * dim3(f[1]) * (f[2] + 1) for f in allf)
    sy3 = sum(f[3] ** 3 * dim3(f[1]) * (f[2] + 1) for f in allf)
    s_su2 = sum(f[3] * dim3(f[1]) for f in allf if f[2] == 1)
    s_su3 = sum(f[3] * (f[2] + 1) for f in allf if f[1] != (0, 0))
    check(NAME, f"one generation of left-handed fields (Q, u^c, d^c, L, e^c): the sums of Y and Y^3 over all components vanish "
                f"({sy}, {sy3}); the SU(2)^2 Y sum over doublets is {s_su2}, and the SU(3)^2 Y sum over triplets and antitriplets is {s_su3}",
          sy == 0 and sy3 == 0 and s_su2 == 0 and s_su3 == 0, eq=eq)
    ndoub = sum(dim3(f[1]) for f in allf if f[2] == 1)
    check(NAME, f"Witten's SU(2) anomaly needs an even number of doublets: one generation has {ndoub} (3 colours of Q and one L)",
          ndoub % 2 == 0, eq=eq)
    z6 = all(((Fr(1, 3) * (0 if f[1] == (0, 0) else (1 if f[1] == (1, 0) else 2)) + Fr(f[2], 2) + f[3]).denominator == 1) for f in allf + (HIGGS,))
    check(NAME, "the Z6 generated by (e^{2 pi i/3}, -1, e^{2 pi i Y}) acts trivially on every field and on the Higgs: t/3 + k/2 + Y is an "
                "integer (t = triality, k = 2 isospin), so the group acting faithfully is (SU(3) x SU(2) x U(1))/Z6", z6, eq=eq)

    Q, uc, dc, L, ec = (dict((f[0], f) for f in FIELDS)[n] for n in ("Q", "u^c", "d^c", "L", "e^c"))
    ys = {"Q H u^c": yukawa([Q, HIGGS, uc]), "Q Htilde d^c": yukawa([Q, HTILDE, dc]), "L Htilde e^c": yukawa([L, HTILDE, ec])}
    bad = {"Q H d^c": yukawa([Q, HIGGS, dc]), "Q Htilde u^c": yukawa([Q, HTILDE, uc]), "L H e^c": yukawa([L, HIGGS, ec]),
           "Q Q": yukawa([Q, Q]), "Q H": yukawa([Q, HIGGS])}
    check(NAME, "Yukawa vertices as invariants of SU(3) x SU(2) x U(1): the number of invariants (SU(3) count from qg_su3_intertwiners, "
                "times SU(2) count, times 1 if the hypercharges sum to zero) is 1 for " + ", ".join(ys) + "; it is 0 for "
                + ", ".join(bad), all(v == 1 for v in ys.values()) and all(v == 0 for v in bad.values()), eq=eq)
    return eq


if __name__ == "__main__":
    build_and_verify()
    print("ok")
