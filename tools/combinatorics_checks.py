#!/usr/bin/env python3
"""Small exact checks that accompany COMBINATORICS_REPRESENTATIONS.md.

Three independent computations, each tied to one statement in the preprints:

  chromatic  Entry 169 (Shareshian-Wachs, elementary positivity).
             For every natural unit interval graph on n <= NMAX vertices,
             expand the chromatic quasisymmetric function in the elementary
             basis by brute force, then test
               (a) every coefficient c_lambda(q) lies in N[q],
               (b) sum_lambda c_lambda(q) = sum_{sigma in D'_G} q^{ginv_G(sigma)}
                   (the bookkeeping behind the paper's identity (1.2)),
               (c) c_lambda(q) = q^M c_lambda(1/q) with M = |E(G)|  (Lemma 2.5),
               (d) unimodality of c_lambda(q) (NOT proved in the paper).
  kl         Entry 168 (combinatorial invariance of Kazhdan-Lusztig polynomials).
             In S_4 and S_5, group all Bruhat intervals by isomorphism class of
             the abstract poset and check that each class carries one polynomial.
  saxl       Entries 205 (Saxl / tensor square conjecture).
             For n <= NMAX, find every partition lambda of n with g(lambda,lambda,nu)>0
             for all nu |- n, using exact Murnaghan-Nakayama characters.

Usage:  python3 tools/combinatorics_checks.py [chromatic|kl|saxl|all] [--figures]
Outputs go to stdout; figures (optional) go to figures/comb_fig*.png.
This reproduces known small cases.  It is not a proof of anything in the papers.
"""
import itertools
import math
import sys
import time
import warnings
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache

import numpy as np

warnings.filterwarnings("ignore", category=UserWarning)

HERE = __import__("os").path.dirname(__import__("os").path.abspath(__file__))
FIGDIR = __import__("os").path.join(HERE, "..", "figures")


# ----------------------------------------------------------------- partitions
def partitions(n, maxpart=None):
    if maxpart is None:
        maxpart = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def conj(lam):
    return tuple(sum(1 for x in lam if x > j) for j in range(lam[0])) if lam else ()


# ============================================================ 1. chromatic
def natural_unit_interval_graphs(n):
    """h weakly increasing with i <= h(i) <= n  (1-indexed); edges {i,j}, i<j<=h(i)."""
    def rec(i, prev, acc):
        if i > n:
            yield tuple(acc)
            return
        for v in range(max(prev, i), n + 1):
            yield from rec(i + 1, v, acc + [v])
    return rec(1, 1, [])


def edges_of(h):
    n = len(h)
    return {(i, j) for i in range(1, n + 1) for j in range(i + 1, h[i - 1] + 1)}


def e_in_m(lam, mu):
    """coefficient of x^mu in e_lam = number of 0-1 matrices, row sums lam, col sums mu."""
    l = len(mu)

    @lru_cache(maxsize=None)
    def go(i, rem):
        if i == len(lam):
            return 1 if all(r == 0 for r in rem) else 0
        tot = 0
        for S in itertools.combinations(range(l), lam[i]):
            r = list(rem)
            ok = True
            for a in S:
                r[a] -= 1
                if r[a] < 0:
                    ok = False
                    break
            if ok:
                tot += go(i + 1, tuple(r))
        return tot
    return go(0, tuple(mu))


def monomial_coefficients(h):
    """chi_G coefficient of x^mu as a list (index = power of q) for each partition mu of n."""
    n = len(h)
    E = edges_of(h)
    M = len(E)
    out = {}
    for mu in partitions(n):
        arr = np.zeros(M + 1, dtype=np.int64)
        seen = set()
        labels = [a for a, m in enumerate(mu) for _ in range(m)]
        for f in set(itertools.permutations(labels)):
            if seen is None:
                pass
            ok = True
            asc = 0
            for (a, b) in E:
                fa, fb = f[a - 1], f[b - 1]
                if fa == fb:
                    ok = False
                    break
                if fa < fb:
                    asc += 1
            if ok:
                arr[asc] += 1
        out[mu] = arr
    return out


def e_expansion(h):
    n = len(h)
    parts = list(partitions(n))
    import sympy
    Mm = sympy.Matrix([[e_in_m(lam, mu) for lam in parts] for mu in parts])  # rows mu, cols lam
    Minv = Mm.inv()
    mono = monomial_coefficients(h)
    Mq = len(next(iter(mono.values()))) - 1
    c = {}
    for j, lam in enumerate(parts):
        coef = np.zeros(Mq + 1, dtype=object)
        for i, mu in enumerate(parts):
            r = Minv[j, i]
            if r != 0:
                coef = coef + np.array([int(r) * int(x) for x in mono[mu]], dtype=object)
        c[lam] = [int(x) for x in coef]
    return c, Mq


def d_prime_ginv(h):
    n = len(h)
    E = edges_of(h)
    cnt = Counter()
    for s in itertools.permutations(range(1, n + 1)):
        if any(s[i] > s[i + 1] and (min(s[i], s[i + 1]), max(s[i], s[i + 1])) not in E for i in range(n - 1)):
            continue
        g = sum(1 for i in range(n) for j in range(i + 1, n)
                if s[i] > s[j] and (s[j], s[i]) in E)
        cnt[g] += 1
    return cnt


def unimodal(a):
    i = 0
    while i + 1 < len(a) and a[i] <= a[i + 1]:
        i += 1
    while i + 1 < len(a) and a[i] >= a[i + 1]:
        i += 1
    return i == len(a) - 1


def run_chromatic(nmax=5):
    print("== chromatic: elementary expansion of chi_G for natural unit interval graphs ==")
    print("n  graphs  all c in N[q]  sum c = sum q^ginv  palindromic  non-unimodal")
    for n in range(1, nmax + 1):
        gs = list(natural_unit_interval_graphs(n))
        ok_pos = ok_sum = ok_pal = 0
        nonuni = 0
        for h in gs:
            c, Mq = e_expansion(h)
            pos = all(all(x >= 0 for x in v) for v in c.values())
            tot = [0] * (Mq + 1)
            for v in c.values():
                for k, x in enumerate(v):
                    tot[k] += x
            gi = d_prime_ginv(h)
            ref = [gi.get(k, 0) for k in range(Mq + 1)]
            ok_pos += pos
            ok_sum += (tot == ref)
            ok_pal += all(v == v[::-1] for v in c.values())
            nonuni += sum(0 if unimodal(v) else 1 for v in c.values())
        print(f"{n}  {len(gs):5d}  {ok_pos:5d}          {ok_sum:5d}               {ok_pal:5d}        {nonuni}")
    # the paper's example (1.3): path 1-2-3
    c, _ = e_expansion((2, 3, 3))
    print("example (1.3), path 1-2-3 (h=(2,3,3)):",
          {lam: v for lam, v in c.items() if any(v)})
    # complete graph: chi = [n]_q! e_n, trivial dot action on the flag variety
    for n in (3, 4):
        c, _ = e_expansion(tuple([n] * n))
        print(f"complete graph K{n}: nonzero c_lambda =",
              {lam: v for lam, v in c.items() if any(v)})


# ================================================================== 2. KL
def perm_length(w):
    return sum(1 for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j])


def bruhat_leq(x, w):
    n = len(x)
    for k in range(1, n):
        a = sorted(x[:k])
        b = sorted(w[:k])
        if any(ai > bi for ai, bi in zip(a, b)):
            return False
    return True


def kl_polys(n):
    """P_{x,w}(q) for S_n as dict[(x,w)] -> tuple of coefficients, via the standard recursion."""
    perms = list(itertools.permutations(range(1, n + 1)))
    L = {w: perm_length(w) for w in perms}
    perms.sort(key=lambda w: L[w])

    def rmul(w, i):  # w * s_i  (swap positions i,i+1)
        w = list(w)
        w[i], w[i + 1] = w[i + 1], w[i]
        return tuple(w)

    def padd(a, b):
        m = max(len(a), len(b))
        return tuple((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(m))

    def shift(a, k):
        return (0,) * k + tuple(a)

    def psub(a, b):
        m = max(len(a), len(b))
        r = [(a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0) for i in range(m)]
        while len(r) > 1 and r[-1] == 0:
            r.pop()
        return tuple(r)

    below = {w: [x for x in perms if L[x] <= L[w] and bruhat_leq(x, w)] for w in perms}
    P = {}
    mu = {}
    for w in perms:
        if L[w] == 0:
            P[(w, w)] = (1,)
            continue
        i = next(i for i in range(n - 1) if w[i] > w[i + 1])
        v = rmul(w, i)
        for x in below[w]:
            xs = rmul(x, i)
            xs_less = L[xs] < L[x]
            # P_{x,w} = q^{1-c} P_{xs,v} + q^c P_{x,v} - sum_z mu(z,v) q^{(l(w)-l(z))/2} P_{x,z}
            c = 1 if xs_less else 0
            t1 = P.get((xs, v), (0,))
            t2 = P.get((x, v), (0,))
            res = padd(shift(t1, 1 - c), shift(t2, c))
            for z in below[v]:
                zs = rmul(z, i)
                if L[zs] < L[z] and bruhat_leq(x, z) and z != v:
                    m = mu.get((z, v), 0)
                    if m:
                        d = (L[w] - L[z]) // 2
                        pz = P.get((x, z), (0,))
                        res = psub(res, tuple(m * t for t in shift(pz, d)))
            while len(res) > 1 and res[-1] == 0:
                res = res[:-1]
            P[(x, w)] = res
        for x in below[w]:
            d = L[w] - L[x]
            if d % 2 == 1:
                p = P[(x, w)]
                k = (d - 1) // 2
                mu[(x, w)] = p[k] if k < len(p) else 0
    return perms, L, P, below


def interval_digraph(x, w, below_w, L):
    import networkx as nx
    elems = [y for y in below_w if bruhat_leq(x, y)]
    G = nx.DiGraph()
    G.add_nodes_from(elems)
    for a in elems:
        for b in elems:
            if L[b] == L[a] + 1 and bruhat_leq(a, b):
                G.add_edge(a, b)
    return G


def run_kl(nmax=5, figures=False):
    import networkx as nx
    from networkx.algorithms.isomorphism import DiGraphMatcher
    print("== kl: poset isomorphism classes of Bruhat intervals carry one KL polynomial ==")
    for n in range(3, nmax + 1):
        t0 = time.time()
        perms, L, P, below = kl_polys(n)
        buckets = defaultdict(list)
        for w in perms:
            for x in below[w]:
                G = interval_digraph(x, w, below[w], L)
                key = (len(G), L[w] - L[x], G.number_of_edges(),
                       nx.weisfeiler_lehman_graph_hash(G))
                buckets[key].append(((x, w), G))
        classes = bad = 0
        nontriv = set()
        ncomp = 0
        for key, lst in buckets.items():
            reps = []
            for (xw, G) in lst:
                ncomp += 1
                for rep in reps:
                    if DiGraphMatcher(G, rep[1]).is_isomorphic():
                        rep[2].append(xw)
                        break
                else:
                    reps.append([xw, G, [xw]])
            for rep in reps:
                classes += 1
                polys = {P[xw] for xw in rep[2]}
                if len(polys) > 1:
                    bad += 1
                nontriv |= polys
        print(f"S_{n}: {ncomp} intervals, {classes} poset-isomorphism classes, "
              f"classes with >1 polynomial: {bad}, distinct polynomials: {len(nontriv)} "
              f"({time.time()-t0:.0f}s)")
        if n == 4:
            ex = [((1, 3, 2, 4), (3, 4, 1, 2)), ((2, 1, 4, 3), (4, 2, 3, 1))]
            for (x, w) in ex:
                G = interval_digraph(x, w, below[w], L)
                print(f"   P_{{{''.join(map(str,x))},{''.join(map(str,w))}}} = {P[(x,w)]}, "
                      f"interval size {len(G)}, rank {L[w]-L[x]}")
            if figures:
                kl_figure(perms, L, P, below, ex)


def kl_figure(perms, L, P, below, ex):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.6))
    for ax, (x, w) in zip(axs, ex):
        G = interval_digraph(x, w, below[w], L)
        layers = defaultdict(list)
        for y in G.nodes:
            layers[L[y] - L[x]].append(y)
        pos = {}
        for r, ys in layers.items():
            ys.sort()
            for k, y in enumerate(ys):
                pos[y] = ((k - (len(ys) - 1) / 2) * 1.6, r)
        for a, b in G.edges:
            ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], color="#8899aa", lw=1, zorder=1)
        for y, (px, py) in pos.items():
            ax.text(px, py, "".join(map(str, y)), ha="center", va="center", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#335", lw=0.8), zorder=2)
        p = P[(x, w)]
        poly = " + ".join((str(c) if k == 0 else (("" if c == 1 else str(c)) + "q" + ("" if k == 1 else f"^{k}")))
                          for k, c in enumerate(p) if c)
        ax.set_title(f"[{''.join(map(str,x))}, {''.join(map(str,w))}] in $S_4$:  P = {poly}\n"
                     f"{len(G)} elements, rank {L[w]-L[x]}", fontsize=10)
        ax.axis("off")
    fig.suptitle("The two Bruhat intervals of S4 with a nontrivial Kazhdan-Lusztig polynomial", fontsize=11)
    fig.tight_layout()
    import os
    os.makedirs(FIGDIR, exist_ok=True)
    fig.savefig(os.path.join(FIGDIR, "comb_fig1_kl_intervals.png"), dpi=160)
    print("   wrote figures/comb_fig1_kl_intervals.png")


# =============================================================== 3. Saxl
def beta_set(lam, size):
    lam = list(lam) + [0] * (size - len(lam))
    return {lam[i] + (size - 1 - i) for i in range(size)}


@lru_cache(maxsize=None)
def chi(lam, rho):
    """Murnaghan-Nakayama character chi^lam(rho); lam, rho tuples (partitions of same n)."""
    if not rho:
        return 1 if not lam else 0
    r = rho[0]
    rest = rho[1:]
    size = len(lam)
    beads = beta_set(lam, size)
    total = 0
    for b in sorted(beads):
        nb = b - r
        if nb < 0 or nb in beads:
            continue
        between = sum(1 for x in beads if nb < x < b)
        new = (beads - {b}) | {nb}
        srt = sorted(new, reverse=True)
        newlam = tuple(x for x in (srt[i] - (size - 1 - i) for i in range(size)) if x > 0)
        total += (-1) ** between * chi(newlam, rest)
    return total


def zee(rho):
    c = Counter(rho)
    z = 1
    for i, m in c.items():
        z *= i ** m * math.factorial(m)
    return z


def universal_squares(n):
    parts = list(partitions(n))
    fact = math.factorial(n)
    csize = [fact // zee(r) for r in parts]
    X = [[chi(lam, r) for r in parts] for lam in parts]  # X[lam][rho]
    univ = []
    for i, lam in enumerate(parts):
        w = [csize[k] * X[i][k] ** 2 for k in range(len(parts))]
        ok = True
        for j in range(len(parts)):
            g = sum(w[k] * X[j][k] for k in range(len(parts)))
            assert g % fact == 0
            if g <= 0:
                ok = False
                break
        if ok:
            univ.append(lam)
    return univ


PRIMES = (2147483629, 2147483587)  # two primes below 2^31


def universal_squares_mod(n):
    """Same test as universal_squares, but g is computed modulo two primes (as in the paper's
    Section 7).  A nonzero residue certifies g > 0; if both residues vanish the test is
    inconclusive and the partition is NOT counted.  Returns (certified, inconclusive)."""
    parts = list(partitions(n))
    fact = math.factorial(n)
    csize = [fact // zee(r) for r in parts]
    X = [[chi(lam, r) for r in parts] for lam in parts]
    res = []
    for p in PRIMES:
        Xp = np.array([[x % p for x in row] for row in X], dtype=np.int64)
        cp = np.array([c % p for c in csize], dtype=np.int64)
        invf = pow(fact % p, p - 2, p)
        res.append((p, Xp, cp, invf))
    cert, inconc = [], []
    for i, lam in enumerate(parts):
        ok = True
        for (p, Xp, cp, invf) in res:
            w = (cp * Xp[i] % p) * Xp[i] % p            # c_rho chi_lam(rho)^2 mod p
            g = ((Xp * w[None, :]) % p).sum(axis=1) % p  # one value per nu
            g = g * invf % p
            bad = np.nonzero(g == 0)[0]
            if len(bad):
                ok = False
                break
        # (a partition passes only if both primes give nonzero residues for every nu)
        if ok:
            cert.append(lam)
    return cert


def run_saxl(nmax=14, figures=False):
    print("== saxl: partitions lambda with a universal tensor square, n <= %d ==" % nmax)
    print("n   #universal  self-conjugate  triangular?  example")
    counts = {}
    for n in range(1, nmax + 1):
        t0 = time.time()
        U = universal_squares(n) if n <= 15 else universal_squares_mod(n)
        counts[n] = len(U)
        tri = any(m * (m + 1) // 2 == n for m in range(1, 10))
        stair = None
        if tri:
            m = next(m for m in range(1, 10) if m * (m + 1) // 2 == n)
            stair = tuple(range(m, 0, -1))
        sc = all(conj(l) == l for l in U)
        ex = U[0] if U else None
        extra = ""
        if stair is not None:
            extra = f"  staircase {stair} universal: {stair in U}"
        print(f"{n:2d}  {len(U):6d}      {str(sc):5s}          {str(tri):5s}       {ex}{extra}   ({time.time()-t0:.0f}s)")
    zero = [n for n, c in counts.items() if c == 0]
    print("degrees with NO universal tensor square:", zero)
    if figures:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import os
        ns = sorted(counts)
        fig, ax = plt.subplots(figsize=(8.2, 4.0))
        tri = lambda n: any(m * (m + 1) // 2 == n for m in range(1, 10))
        cols = ["#c0392b" if counts[n] == 0 else ("#1f5f99" if tri(n) else "#7aa6c2") for n in ns]
        bars = ax.bar(ns, [counts[n] for n in ns], color=cols)
        for n, b in zip(ns, bars):
            if n > 15:
                b.set_hatch("//")
                b.set_edgecolor("white")
            if counts[n] == 0:
                ax.text(n, 0.15, "none", ha="center", fontsize=9, color="#c0392b", fontweight="bold")
        ax.set_xticks(ns)
        ax.set_yticks(range(0, max(counts.values()) + 1))
        ax.set_xlabel("n")
        ax.set_ylabel("number of lambda with a universal square")
        ax.set_title("Partitions lambda of n whose tensor square S^lambda (x) S^lambda contains every irreducible\n"
                     "red: none (n = 2, 4, 9); dark blue: triangular n; hatched: certified modulo two primes", fontsize=9)
        fig.tight_layout()
        os.makedirs(FIGDIR, exist_ok=True)
        fig.savefig(os.path.join(FIGDIR, "comb_fig2_universal_squares.png"), dpi=160)
        print("wrote figures/comb_fig2_universal_squares.png")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    fig = "--figures" in sys.argv
    what = args[0] if args else "all"
    if what in ("chromatic", "all"):
        run_chromatic(6)
    if what in ("kl", "all"):
        run_kl(5, fig)
    if what in ("saxl", "all"):
        nmax = 14
        for a in sys.argv[1:]:
            if a.startswith("--nmax="):
                nmax = int(a.split("=")[1])
        run_saxl(nmax, fig)
