"""Entry 211 (random planar maps, critical Fortuin-Kasteleyn maps and Liouville quantum gravity) beside the Catalan counts of the spin-network paper.

Source: the 211 papers of the OpenAI release: "Canonical conformal limits of subcritical FK planar maps" [FKs] (Sep 24 2026, abstract, section 1.1 and section 3: the
parametrisation q = 2 + 2 cos(pi gamma^2/2), kappa = 16/gamma^2, loop weight q^{l/2}), "The critical Liouville quantum sphere and geometric limits of FK maps at q = 4" [FK4]
(Sep 24 2026, abstract and introduction: gamma = 2 at q = 4), and the Ising companions [FKw, FKe] (gamma = sqrt3, the sqrt3-quantum sphere). The fork read the abstracts and
introductions and the page of the statements above; it did not read the proofs. This module checks what a reader can recompute:

  qg_maps_counts     Rooted planar maps with n edges, by brute force over pairs of permutations (sigma, alpha) of 2n half-edges with genus 0 and connectedness, for n = 1..4, against
                     Tutte's M_n = 2 * 3^n (2n)!/(n!(n+2)!) (2, 9, 54, 378); the same enumeration weighted by the number of spanning trees (Kirchhoff) against Mullin's
                     Cat_n Cat_{n+1} (2, 10, 70, 588); growth 12^n n^{-5/2} and 16^n n^{-3} (exact ratios, large-n asymptotics).
  qg_maps_exponents  The dictionary q -> kappa -> gamma -> c: sqrt q = -2 cos(4 pi/kappa), c = (3 kappa - 8)(6 - kappa)/(2 kappa) = 25 - 6 Q^2 (Q = 2/gamma + gamma/2), KPZ
                     gamma_str = (c - 1 - sqrt((1-c)(25-c)))/12 = 1 - 4/gamma^2; at q = 0, 1, 2, 3, 4 these are c = -2, 0, 1/2, 4/5, 1 and gamma_str = -1, -1/2, -1/3, -1/5, 0; rooted
                     counts grow as mu^n n^(gamma_str - 2), which reproduces the exponents -3 (Mullin, q = 0) and -5/2 (Tutte, pure gravity) found by brute force above.
  qg_maps_su2k       The link to the spin networks of qg_spin_network.py: the invariants of N spin-1/2 legs have dimension Catalan(N/2) (undeformed SU(2), loop weight 2 = sqrt q at q = 4) and at level k
                     of the quantum group the dimension of the invariant space is the number of closed walks of length N on the Dynkin diagram A_{k+1}, which equals the Catalan number while N/2 <= k
                     and grows as (2 cos(pi/(k+2)))^N = (sqrt q_k)^N with q_k = 4 cos^2(pi/(k+2)): k = 1, 2, 4, infinity give q = 1, 2, 3, 4 (percolation, Ising, three-state Potts, the critical end).

Not tested: any theorem of entry 211 (convergence to LQG spheres, the q > 4 tree phase, Liouville Brownian motion, spectral convergence), the physical identification of an
SU(2)_k spin network with a loop-weighted map (the numbers agree; no map between the objects is constructed here), and the Immirzi-type constant for k < infinity (stated only as the
growth rate ln sqrt q of a count, not as a physical constant). The planar maps of entry 211 are random surfaces with a loop decoration; the spin networks are states of a Hilbert space.
"""
import itertools
import math
import sys
from fractions import Fraction as Fr
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register

NAME = "qg_planar_maps"
ENTRY = "211"
PAPER = "[FKs] OpenAI, Canonical conformal limits of subcritical FK planar maps (preprints, September 24 2026), section 1.1"


def cat(n):
    return math.comb(2 * n, n) // (n + 1)


def tutte(n):
    return 2 * 3 ** n * math.factorial(2 * n) // (math.factorial(n) * math.factorial(n + 2))


def cycles(perm):
    seen, out = [False] * len(perm), []
    for s in range(len(perm)):
        if not seen[s]:
            c, x = [], s
            while not seen[x]:
                seen[x] = True
                c.append(x)
                x = perm[x]
            out.append(c)
    return out


def rooted_counts(n):
    """(maps, tree-weighted maps) with n edges: sum over sigma in S_{2n} for the fixed alpha (0 1)(2 3)...; genus 0 and connected; the factor turning
    labelled counts into rooted counts is (2n-1)!!/(2n-1)! = 1/(2^(n-1) (n-1)!)."""
    h = 2 * n
    alpha = [i ^ 1 for i in range(h)]
    count, trees = 0, 0
    for sigma in itertools.permutations(range(h)):
        sc = cycles(sigma)
        V = len(sc)
        sa = [sigma[alpha[i]] for i in range(h)]
        F = len(cycles(sa))
        if V - n + F != 2:
            continue
        vert = [0] * h
        for vi, c in enumerate(sc):
            for x in c:
                vert[x] = vi
        # connected? the graph on vertices with an edge per alpha-pair
        parent = list(range(V))

        def find(u):
            while parent[u] != u:
                parent[u] = parent[parent[u]]
                u = parent[u]
            return u
        for i in range(0, h, 2):
            parent[find(vert[i])] = find(vert[i + 1])
        if len({find(u) for u in range(V)}) != 1:
            continue
        count += 1
        L = np.zeros((V, V))
        for i in range(0, h, 2):
            u, w = vert[i], vert[i + 1]
            if u != w:
                L[u, u] += 1
                L[w, w] += 1
                L[u, w] -= 1
                L[w, u] -= 1
        trees += int(round(np.linalg.det(L[1:, 1:]))) if V > 1 else 1
    scale = Fr(1, 2 ** (n - 1) * math.factorial(n - 1))
    return count * scale, trees * scale


def build_and_verify():
    # ------------------------------------------------------------------------------------------------ counts
    brute = {n: rooted_counts(n) for n in (1, 2, 3, 4)}
    eq_c = register(NAME, Equation(
        key="qg_maps_counts", entry=ENTRY, title="Rooted planar maps and tree-rooted maps: brute force against Tutte and Mullin",
        source_kind="derived", source=f"derived here; the objects are the planar maps of entry 211 ({PAPER}); the formulas are Tutte's and Mullin's, rechecked",
        latex_str=(r"M_n=\frac{2\cdot3^{n}(2n)!}{n!\,(n+2)!}=2,9,54,378,\ldots\ \sim\ \frac{2}{\sqrt\pi}\,\frac{12^{n}}{n^{5/2}},\qquad "
                   r"T_n=\mathrm{Cat}_n\mathrm{Cat}_{n+1}=2,10,70,588,\ldots\ \sim\ \frac{4}{\pi}\,\frac{16^{n}}{n^{3}}"
                   r"%%BR%%\mathrm{Cat}_n=\frac1{n+1}\binom{2n}{n}\sim\frac{4^{n}}{\sqrt\pi\,n^{3/2}}\quad(\text{plane trees, and the invariants of }2n\text{ spin-}\tfrac12\text{ legs})"),
        notes="Rooted planar maps with n edges are counted as pairs (sigma, alpha) of permutations of the 2n half-edges, alpha a fixed-point-free involution (the edges) and "
              "sigma the vertex rotation, of genus 0 (V - n + F = 2, with F the number of cycles of sigma alpha) and connected; for fixed alpha the sum over sigma in S_{2n} is turned into "
              "the rooted count by the factor 1/(2^(n-1) (n-1)!). The same enumeration weighted by the number of spanning trees of the map's graph (Kirchhoff's determinant, loops "
              "contribute nothing) gives the tree-rooted maps. Brute force reaches n = 4 (8! permutations); the closed forms are checked against it and their growth is checked "
              "at n = 10^6 with log-gamma. The tree-rooted maps are the q -> 0 (spanning-tree-weighted) maps of entry 211. Not tested: Mullin's bijection with shuffles of two "
              "Dyck words, and any statement about the scaling limit.",
    ))
    check(NAME, f"brute force over all pairs (sigma, alpha) at n = 1..4 gives {[int(brute[n][0]) for n in brute]} rooted planar maps, equal to Tutte's 2 * 3^n (2n)!/(n! (n+2)!) = {[tutte(n) for n in (1, 2, 3, 4)]}",
          all(brute[n][0] == tutte(n) for n in brute) and [tutte(n) for n in (1, 2, 3, 4)] == [2, 9, 54, 378], eq=eq_c)
    check(NAME, f"the same enumeration weighted by the number of spanning trees gives {[int(brute[n][1]) for n in brute]}, equal to Mullin's Cat_n Cat_(n+1) = {[cat(n) * cat(n + 1) for n in (1, 2, 3, 4)]}",
          all(brute[n][1] == cat(n) * cat(n + 1) for n in brute) and [cat(n) * cat(n + 1) for n in (1, 2, 3, 4)] == [2, 10, 70, 588], eq=eq_c)
    mp.mp.dps = 30
    n_big = mp.mpf(10) ** 6
    M_big = 2 * mp.mpf(3) ** n_big * mp.e ** (mp.loggamma(2 * n_big + 1) - mp.loggamma(n_big + 1) - mp.loggamma(n_big + 3))
    T_big_log = (mp.loggamma(2 * n_big + 1) - 2 * mp.loggamma(n_big + 1) - mp.log(n_big + 1)) + (mp.loggamma(2 * n_big + 3) - 2 * mp.loggamma(n_big + 2) - mp.log(n_big + 2))
    # ratio to the asymptotic forms, in logs
    lM = mp.log(2) + n_big * mp.log(3) + mp.loggamma(2 * n_big + 1) - mp.loggamma(n_big + 1) - mp.loggamma(n_big + 3)
    lM_as = mp.log(2 / mp.sqrt(mp.pi)) + n_big * mp.log(12) - mp.mpf(5) / 2 * mp.log(n_big)
    lT_as = mp.log(4 / mp.pi) + n_big * mp.log(16) - 3 * mp.log(n_big)
    check(NAME, f"asymptotics at n = 10^6 (log-gamma): log M_n - log[(2/sqrt pi) 12^n n^(-5/2)] = {mp.nstr(lM - lM_as, 3)} and log T_n - log[(4/pi) 16^n n^(-3)] = {mp.nstr(T_big_log - lT_as, 3)} (both of order 1/n)",
          abs(lM - lM_as) < 1e-5 and abs(T_big_log - lT_as) < 1e-5, eq=eq_c)
    check(NAME, "growth rates: M_(n+1)/M_n -> 12 and T_(n+1)/T_n -> 16 (exact ratios at n = 2000), and Cat_(n+1)/Cat_n -> 4",
          abs(float(Fr(tutte(2001), tutte(2000))) - 12) < 0.02 and abs(float(Fr(cat(2001) * cat(2002), cat(2000) * cat(2001))) - 16) < 0.04 and abs(float(Fr(cat(2001), cat(2000))) - 4) < 0.01, eq=eq_c)

    shuffles = [sum(math.comb(2 * n, 2 * k) * cat(k) * cat(n - k) for k in range(n + 1)) for n in range(0, 41)]
    check(NAME, "Mullin's count as shuffles of two Dyck words: sum_k binom(2n, 2k) Cat_k Cat_(n-k) = Cat_n Cat_(n+1) for n = 0..40 (exact integers): a tree-rooted map is a pair of "
                "plane-tree words (the tree and its dual) interleaved, each counted by the Catalan numbers that count the invariants of 2k and 2(n-k) spin-1/2 legs",
          all(shuffles[n] == cat(n) * cat(n + 1) for n in range(41)), eq=eq_c)

    # ------------------------------------------------------------------------------------------------ exponents
    g, kap = sp.symbols("gamma kappa", positive=True)
    Q = 2 / g + g / 2
    c_of_kappa = lambda k: (3 * k - 8) * (6 - k) / (2 * k)
    kpg = 16 / g ** 2
    gam_str = lambda c: (c - 1 - sp.sqrt((1 - c) * (25 - c))) / 12
    rows = []
    for qv, kv in ((0, sp.Integer(8)), (1, sp.Integer(6)), (2, sp.Rational(16, 3)), (3, sp.Rational(24, 5)), (4, sp.Integer(4))):
        gv = 4 / sp.sqrt(kv)
        rows.append((qv, kv, sp.simplify(gv ** 2), c_of_kappa(kv), sp.simplify(gam_str(c_of_kappa(kv)))))
    q_of_g = lambda gv: sp.simplify(2 + 2 * sp.cos(sp.pi * gv ** 2 / 2))
    eq_e = register(NAME, Equation(
        key="qg_maps_exponents", entry=ENTRY, title="The FK dictionary q, kappa, gamma, c and the string susceptibility",
        source_kind="derived", source=f"derived here from the parametrisation of {PAPER}, q = 2 + 2 cos(pi gamma^2/2), kappa = 16/gamma^2",
        latex_str=(r"q=2+2\cos\frac{\pi\gamma^{2}}2,\ \ \kappa'=\frac{16}{\gamma^{2}},\ \ \sqrt q=-2\cos\frac{4\pi}{\kappa'},\qquad c=\frac{(3\kappa'-8)(6-\kappa')}{2\kappa'}=25-6Q^{2},\ \ Q=\frac2\gamma+\frac\gamma2"
                   r"%%BR%%\gamma_{\rm str}=\frac{c-1-\sqrt{(1-c)(25-c)}}{12}=1-\frac4{\gamma^{2}},\qquad \#\{\text{rooted maps}\}\sim\mu^{n}n^{\gamma_{\rm str}-2}"
                   r"%%BR%%q=0,1,2,3,4:\ \ \gamma^{2}=2,\tfrac83,3,\tfrac{10}3,4;\ \ c=-2,0,\tfrac12,\tfrac45,1;\ \ \gamma_{\rm str}=-1,-\tfrac12,-\tfrac13,-\tfrac15,0"),
        notes="The parametrisation is the paper's (equation 1.1: q in (0,4), gamma in (sqrt2, 2), kappa in (4, 8) with kappa = 16/gamma^2). The central charge of the matter "
              "c = (3 kappa - 8)(6 - kappa)/(2 kappa) and the KPZ string susceptibility are standard and are used here as such; the identities between the forms are exact. The rooted-count "
              "rule mu^n n^(gamma_str - 2) is the standard scaling form: it is checked against the two brute-force cases above (Tutte, c = 0, exponent -5/2, and Mullin, c = -2, exponent -3), "
              "not derived. The q = 4 end is the critical Liouville sphere of [FK4] (gamma = 2, c = 1, gamma_str = 0, with logarithmic corrections); the paper says q > 4 leaves the surface phase "
              "for the Brownian tree, which is outside this table. Not tested: that the limits of entry 211 are LQG surfaces of these parameters (that is its theorem).",
    ))
    check(NAME, "q = 2 + 2 cos(pi gamma^2/2) with gamma^2 = 16/kappa' gives q = 0, 1, 2, 3, 4 at kappa' = 8, 6, 16/3, 24/5, 4, i.e. gamma^2 = 2, 8/3, 3, 10/3, 4; and sqrt q = -2 cos(4 pi/kappa') at the same points (exact)",
          [q_of_g(sp.sqrt(r[2])) for r in rows] == [0, 1, 2, 3, 4] and all(sp.simplify(-2 * sp.cos(4 * sp.pi / r[1]) - sp.sqrt(r[0])) == 0 for r in rows), eq=eq_e)
    check(NAME, "the Ising point of the entry: q = 2 gives gamma = sqrt3 (the 'sqrt3-quantum sphere'), kappa' = 16/3 and c = 1/2; q = 4 gives gamma = 2, c = 1",
          rows[2][2] == 3 and rows[2][3] == sp.Rational(1, 2) and rows[4][2] == 4 and rows[4][3] == 1, eq=eq_e)
    check(NAME, "c = (3 kappa' - 8)(6 - kappa')/(2 kappa') with kappa' = 16/gamma^2 equals 25 - 6 Q^2 = 26 - c_L with c_L = 1 + 6 Q^2, Q = 2/gamma + gamma/2: matter and Liouville central charges sum to 26 (exact, all gamma)",
          sp.simplify(c_of_kappa(kpg) - (25 - 6 * Q ** 2)) == 0 and sp.simplify(c_of_kappa(kpg) + (1 + 6 * Q ** 2) - 26) == 0, eq=eq_e)
    t = sp.Symbol("t", positive=True)       # t = gamma^2/4 in (0, 1)
    cg = 13 - 6 * (t + 1 / t)
    check(NAME, "KPZ: gamma_str = (c - 1 - sqrt((1-c)(25-c)))/12 equals 1 - 4/gamma^2 for gamma < 2 (with t = gamma^2/4 < 1, sqrt((1-c)(25-c)) = 6 (1/t - t); exact)",
          sp.simplify(((1 - cg) * (25 - cg)) - 36 * (1 / t - t) ** 2) == 0 and sp.simplify((cg - 1 - 6 * (1 / t - t)) / 12 - (1 - 1 / t)) == 0, eq=eq_e)
    check(NAME, f"the table: c = {[str(r[3]) for r in rows]} and gamma_str = {[str(r[4]) for r in rows]} at q = 0, 1, 2, 3, 4 "
                "(known values: -2, 0, 1/2, 4/5, 1 and -1, -1/2, -1/3, -1/5, 0)",
          [r[3] for r in rows] == [-2, 0, sp.Rational(1, 2), sp.Rational(4, 5), 1] and [r[4] for r in rows] == [-1, sp.Rational(-1, 2), sp.Rational(-1, 3), sp.Rational(-1, 5), 0], eq=eq_e)
    check(NAME, "the brute-force exponents are those of the rule mu^n n^(gamma_str - 2): Tutte's n^(-5/2) is gamma_str = -1/2 (c = 0, q = 1) and Mullin's n^(-3) is gamma_str = -1 (c = -2, q = 0); "
                "the exponents -5/2 and -3 were found in the previous display from the exact formulas",
          rows[1][4] - 2 == sp.Rational(-5, 2) and rows[0][4] - 2 == -3, eq=eq_e)

    # ------------------------------------------------------------------------------------------------ SU(2)_k
    def walks(k, Nlen):
        """closed walks of length Nlen from node 0 on the path graph with nodes 0..k (the Dynkin diagram A_{k+1})."""
        v = [0] * (k + 1)
        v[0] = 1
        for _ in range(Nlen):
            v = [(v[i - 1] if i > 0 else 0) + (v[i + 1] if i < k else 0) for i in range(k + 1)]
        return v[0]
    qk = lambda k: 4 * sp.cos(sp.pi / (k + 2)) ** 2
    eq_k = register(NAME, Equation(
        key="qg_maps_su2k", entry=ENTRY, title="Spin-1/2 invariants at level k and the FK loop weight",
        source_kind="derived", source="derived here; the spin-1/2 invariant counts of qg_spin_network.py (Catalan numbers) and the loop weight sqrt q of entry 211",
        latex_str=(r"\dim\mathrm{Inv}_k(V_{1/2}^{\otimes N})=\#\{\text{closed walks of length }N\text{ on }A_{k+1}\}=\mathrm{Cat}_{N/2}\ \ (N\le2k),\qquad "
                   r"\#\sim\Big(2\cos\frac{\pi}{k+2}\Big)^{N}=(\sqrt{q_k})^{N}"
                   r"%%BR%%q_k=4\cos^{2}\frac{\pi}{k+2};\qquad k=1,2,4,\infty:\ \ q_k=1,2,3,4,\qquad \tfrac1N\ln\#\ \to\ \ln2\cos\tfrac{\pi}{k+2}=\tfrac12\ln q_k\ \ (\ln2\ \text{at }k=\infty)"),
        notes="At level k of the quantum group U_q(sl_2), q = exp(i pi/(k+2)), the invariants of N spin-1/2 legs are counted by closed walks of length N from the trivial "
              "representation on the Dynkin diagram A_{k+1} (the fusion graph of spin 1/2, truncated at spin k/2); for N <= 2k the truncation does not act and the count is the Catalan number of "
              "the undeformed SU(2) (the count of qg_spin_network.py). The loop weight of the spin-1/2 line is the quantum dimension 2 cos(pi/(k+2)) = sqrt(q_k), which is the weight "
              "sqrt q per loop of the FK model of entry 211 at q = q_k: q = 1, 2, 3 are k = 1, 2, 4 (percolation, Ising, three-state Potts), and the undeformed SU(2), loop weight 2 = dim V_{1/2}, "
              "is q = 4, the critical end of the range 0 < q <= 4 of the entry. The numbers agree because the same Temperley-Lieb algebra with loop weight sqrt q acts in both; no map from a "
              "spin network to a decorated planar map is constructed here, and the growth rate (1/2) ln q is the growth rate of a count, not an Immirzi-type constant (for k = infinity "
              "it is the ln 2 that fixes gamma_I in qg_spin_network.py by a different, physical, argument). The Beraha value k = 3 gives q = ((1 + sqrt5)/2)^2, not an integer.",
    ))
    check(NAME, "q_k = 4 cos^2(pi/(k+2)) is 1, 2, 3 at k = 1, 2, 4 and tends to 4 as k -> infinity; k = 3 gives ((1 + sqrt5)/2)^2 = 2.618 (exact)",
          [sp.simplify(qk(k)) for k in (1, 2, 4)] == [1, 2, 3] and sp.limit(qk(sp.Symbol("k", positive=True)), sp.Symbol("k", positive=True), sp.oo) == 4 and sp.simplify(qk(3) - ((1 + sp.sqrt(5)) / 2) ** 2) == 0, eq=eq_k)
    check(NAME, "for N = 2m <= 2k the count of closed walks on A_(k+1) is the Catalan number Cat_m (k = 6, m = 0..6: " + ", ".join(str(walks(6, 2 * m)) for m in range(7)) + "); "
                "it is smaller once the truncation acts (k = 2, m = 3: " + str(walks(2, 6)) + " against Cat_3 = 5)",
          all(walks(6, 2 * m) == cat(m) for m in range(7)) and walks(2, 6) == 4 < cat(3), eq=eq_k)
    check(NAME, "the largest eigenvalue of the adjacency matrix of A_(k+1) is 2 cos(pi/(k+2)) = sqrt(q_k), for k = 1..8 (numerical eigenvalues against the closed form, to 1e-12)",
          all(abs(max(np.linalg.eigvalsh(np.diag([1.0] * k, 1) + np.diag([1.0] * k, -1))) - 2 * math.cos(math.pi / (k + 2))) < 1e-12 for k in range(1, 9)), eq=eq_k)
    ratios = {k: walks(k, 2 * 61) / walks(k, 2 * 60) for k in (1, 2, 4, 10)}
    check(NAME, f"growth of the invariant counts: W_(N+2)/W_N at N = 120 is {ratios[1]:.6f}, {ratios[2]:.6f}, {ratios[4]:.6f}, {ratios[10]:.4f} for k = 1, 2, 4, 10, against q_k = 1, 2, 3, {float(qk(10)):.4f} "
                "(the ratio of walk counts approaches the square of the top eigenvalue)",
          abs(ratios[1] - 1) < 1e-9 and abs(ratios[2] - 2) < 1e-9 and abs(ratios[4] - 3) < 1e-6 and abs(ratios[10] - float(qk(10))) < 1e-3, eq=eq_k)
    check(NAME, f"the entropy per leg (1/N) ln W_N at N = 400 is {math.log(walks(2, 400)) / 400:.4f}, {math.log(walks(4, 400)) / 400:.4f} for k = 2, 4 against (1/2) ln q_k = {math.log(2) / 2:.4f}, {math.log(3) / 2:.4f} "
                "(corrections of order 1/N); and for k = infinity, (1/N) ln Cat_(N/2) = "+f"{math.log(cat(200)) / 400:.4f}"+" at N = 400 (corrections of order ln N/N) tends to ln 2 = 0.6931, the per-puncture entropy of the section on counting punctures",
          abs(math.log(walks(2, 400)) / 400 - math.log(2) / 2) < 5e-3 and abs(math.log(walks(4, 400)) / 400 - math.log(3) / 2) < 5e-3
          and abs(math.log(cat(200)) / 400 - math.log(2)) < 0.03, eq=eq_k)
    check(NAME, "the Catalan numbers are the numbers of plane trees (and of Dyck paths): Cat_n counts rooted plane trees with n edges, the discrete trees whose scaling limit is the Brownian tree of the q > 4 phase of "
                "entry 211; sum over n of Cat_n x^n satisfies C = 1 + x C^2 (exact), with radius of convergence 1/4",
          sp.simplify(sp.series(sum(cat(i) * sp.Symbol("x") ** i for i in range(8)), sp.Symbol("x"), 0, 8).removeO()
                      - sp.series((1 - sp.sqrt(1 - 4 * sp.Symbol("x"))) / (2 * sp.Symbol("x")), sp.Symbol("x"), 0, 8).removeO()) == 0, eq=eq_k)
