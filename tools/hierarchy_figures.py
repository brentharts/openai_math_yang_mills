#!/usr/bin/env python3
"""Figures and exact computations for SPECTRE_CLOCK_AND_CHAIR44.md (external data; not a repository entry).

Built on the Chair44 pose model of tools/chair44_model.py (atlas and child poses of Tsiokos, arXiv:2609.19214,
transcribed from a text extraction; see that script for the convention and its checks).  Writes four PNG files
to ``figures/`` (committed, so the markdown pages can show them) and numbers to ``build/data/hierarchy_figures.json``.

What each part computes
-----------------------
A. The 11 conjugacy classes of subgroups of the 24 proper cubic rotations, which classes contain an element of
   order 3, and which are 2-groups.  Used to read the symmetry statements of [T] section 7 and [HR] section 7.4.
B. Fixed points of the substitution.  A pose (G,t) is a fixed seed if it is one of its own children.  Maximal sets
   of pairwise compatible fixed seeds are refined 3 and 4 times; those that overlap nothing, create no illegal contact
   and cover a window around the origin are tilings of a growing cube.  Their symmetry groups are measured on the
   window (rotations about half-integer centres near the origin).
C. Contact counts of the one-seed hierarchy to level 7, fitted exactly, and the extensivity defect 8 C_n - C_{n+1}.
D. Clock-style tests on the 8-seed fixed point: the radius at which a chair's neighbourhood fixes its slot in its
   level-k supertile (analogue of [HR] section 10, "records"), and the look-back statistic against an octree control.
E. A noisy-patch toy (each neighbour deleted with probability eps) and the Kesten-Stigum line of entry 229.

Everything in D and E is this fork's own construction and is not a replication of [HR] (different tile, different
definition of "ancestor", no Spectre data).  Parts A to C are exact computations on the transcribed atlas.
"""

from __future__ import annotations

import itertools
import sys
import time
from collections import Counter, defaultdict

import numpy as np
from scipy.spatial import cKDTree

import common as C
import chair44_model as M

FULL = "--full" in sys.argv
FIGDIR = C.ROOT / "figures"
RM = np.array(M.ROTS)  # (24, 3, 3)
ONE = np.array([1, 1, 1])


def idx(m: np.ndarray) -> int:
    return M.ROT_INDEX[m.tobytes()]


COMP = np.array([[idx(RM[i] @ RM[j]) for j in range(24)] for i in range(24)])  # G_i G_j
TR = np.array([idx(RM[i].T) for i in range(24)])
I24 = idx(np.eye(3, dtype=int))
CHILD = [(idx(M.frame(p, s)), np.array(u)) for p, s, u in M.CHILDREN_RAW]
ATLAS = [(idx(G), np.array(t)) for G, t in M.parse_atlas()]
ATLAS_KEYS = {M.pkey(G, t) for G, t in M.parse_atlas()}
CHILD_M = [(M.frame(p, s), np.array(u)) for p, s, u in M.CHILDREN_RAW]


def savefig(fig, name: str) -> None:
    FIGDIR.mkdir(exist_ok=True)
    fig.savefig(FIGDIR / f"{name}.png", dpi=150, bbox_inches="tight")
    import matplotlib.pyplot as plt
    plt.close(fig)


def key(r, t, off=1 << 12):
    t = t + off
    assert t.min() >= 0 and t.max() < (1 << 13)
    return (r.astype(np.int64) << 39) | (t[:, 0] << 26) | (t[:, 1] << 13) | t[:, 2]


def refine_vec(r, t):
    rs, ts = [], []
    for k, u in CHILD:
        rs.append(COMP[r, k])
        ts.append(2 * t + np.einsum("nij,j->ni", RM[r], u))
    return np.stack(rs, 1).reshape(-1), np.stack(ts, 1).reshape(-1, 3)  # index = parent*8 + slot


def centres(r, t):
    return np.einsum("nij,j->ni", RM[r], ONE) + t


# ----------------------------------------------------------------------------------------------
# A. Subgroup classes of the cubic rotation group
# ----------------------------------------------------------------------------------------------
def involution_axis(m: np.ndarray) -> str:
    v = [np.array(v) for v in itertools.product((-1, 0, 1), repeat=3) if any(v) and np.array_equal(m @ np.array(v), np.array(v))][0]
    return "coord" if np.count_nonzero(v) == 1 else "diag"


def elt_order(i: int) -> int:
    k, x = 1, i
    while x != I24:
        x = COMP[x, i]
        k += 1
    return k


def part_a():
    def closure(gens):
        S = {I24} | set(gens)
        while True:
            new = {int(COMP[a, b]) for a in S for b in S} | S
            if new == S:
                return frozenset(S)
            S = new

    subs = {closure(g) for k in (1, 2, 3) for g in itertools.combinations_with_replacement(range(24), k)}
    # conjugacy classes
    def conj(H, g):
        return frozenset(int(COMP[COMP[g, h], TR[g]]) for h in H)
    classes, seen = [], set()
    for H in sorted(subs, key=lambda s: (len(s), sorted(s))):
        if H in seen:
            continue
        orb = {conj(H, g) for g in range(24)}
        seen |= orb
        classes.append(sorted(orb, key=sorted))

    def name(H):
        ords = sorted(elt_order(h) for h in H)
        n = len(H)
        inv = [h for h in H if elt_order(h) == 2]
        axes = sorted(involution_axis(RM[h]) for h in inv)
        if n == 1:
            return "1"
        if n == 2:
            return "C2 (coordinate axis)" if axes == ["coord"] else "C2 (face diagonal)"
        if n == 3:
            return "C3"
        if n == 4:
            return "C4" if 4 in ords else ("V4 (normal)" if axes == ["coord"] * 3 else "V4 (other)")
        if n == 6:
            return "D3"
        if n == 8:
            return "D4"
        if n == 12:
            return "A4"
        return "O (24)"

    info = []
    for orb in classes:
        H = orb[0]
        info.append({"name": name(H), "order": len(H), "size": len(orb),
                     "has3": any(elt_order(h) == 3 for h in H), "members": orb})
    C.check(len(info) == 11, "the cubic rotation group (S_4) has 11 conjugacy classes of subgroups")
    C.check(sum(c["size"] for c in info) == len(subs) == 30, "S_4 has 30 subgroups in total")
    two_groups = [c for c in info if not c["has3"]]
    d4 = next(c for c in info if c["name"] == "D4")
    C.check(d4["size"] == 3, "there are 3 subgroups of order 8 (the Sylow 2-subgroups, all conjugate to D4)")
    C.check(all(any(set(H) <= set(D) for D in d4["members"]) for c in two_groups for H in c["members"]),
            "every subgroup with no element of order 3 lies in a copy of D4, so its order is at most 8")
    C.check(len(two_groups) == 7, "7 of the 11 classes have no element of order 3")
    # cover relations
    rel = defaultdict(set)
    for a, ca in enumerate(info):
        for b, cb in enumerate(info):
            if a != b and ca["order"] < cb["order"] and any(set(ca["members"][0]) <= set(D) for D in cb["members"]):
                rel[a].add(b)
    cover = {(a, b) for a in list(rel) for b in rel[a] if not any(b in rel.get(m, set()) for m in rel[a] if m != b)}
    return info, cover


HR_ALLOWED = {"1", "C2 (face diagonal)", "C4", "D4"}  # [HR] Theorem 5 (as quoted in the fetched text)


def draw_a(ax, info, cover):
    pos = {"1": (0.5, 0), "C2 (coordinate axis)": (0.12, 1), "C2 (face diagonal)": (0.5, 1), "C3": (0.88, 1),
           "C4": (0.04, 2), "V4 (normal)": (0.3, 2), "V4 (other)": (0.56, 2), "D3": (0.88, 2),
           "D4": (0.2, 3), "A4": (0.7, 3), "O (24)": (0.45, 4)}
    for a, b in cover:
        (x1, y1), (x2, y2) = pos[info[a]["name"]], pos[info[b]["name"]]
        ax.plot([x1, x2], [y1, y2], color="0.7", lw=1, zorder=1)
    for c in info:
        x, y = pos[c["name"]]
        if c["name"] in HR_ALLOWED:
            fc, ec = "#b8e0c8", "#1b7a43"
        elif c["has3"]:
            fc, ec = "#f2d0d0", "#a33"
        else:
            fc, ec = "#f6e7b4", "#a77b00"
        ax.text(x, y, f'{c["name"]}\norder {c["order"]}, {c["size"]} conj.', ha="center", va="center", fontsize=6.6,
                bbox=dict(boxstyle="round,pad=0.3", fc=fc, ec=ec), zorder=2)
    ax.set_xlim(-0.12, 1.08); ax.set_ylim(-0.5, 4.5); ax.axis("off")
    ax.set_title("Subgroup classes of the 24 cubic rotations", fontsize=9)
    ax.text(0.0, -0.45, "green: listed in [HR] Thm 5   yellow: 2-group, excluded by [HR] Rem. 1\nred: contains an element of order 3 (excluded)",
            fontsize=6.5, va="bottom")


# ----------------------------------------------------------------------------------------------
# B. Fixed points of the substitution and their symmetry
# ----------------------------------------------------------------------------------------------
def bron_kerbosch(adj, R, P, X, out):
    if not P and not X:
        out.append(sorted(R))
        return
    u = max(P | X, key=lambda v: len(adj[v] & P))
    for v in list(P - adj[u]):
        bron_kerbosch(adj, R | {v}, P & adj[v], X & adj[v], out)
        P = P - {v}
        X = X | {v}


def window_symmetry(r, t, win_radius):
    """Rotations (index) with a half-integer centre c such that x -> R(x-c)+c maps the window poses into the pose set."""
    ks = np.sort(key(r, t))
    cen = centres(r, t)
    w = np.where(np.abs(cen).max(1) <= win_radius)[0]
    out = defaultdict(list)
    for R in range(24):
        for c2 in itertools.product(range(-2, 3), repeat=3):
            c2 = np.array(c2)
            sh = c2 - RM[R] @ c2
            if (sh % 2).any():
                continue
            sh //= 2
            r2 = COMP[R, r[w]]
            t2 = np.einsum("ij,nj->ni", RM[R], t[w]) + sh
            k2 = key(r2, t2)
            pos = np.searchsorted(ks, k2)
            pos[pos == len(ks)] = len(ks) - 1
            if (ks[pos] == k2).all():
                out[tuple(int(x) for x in c2)].append(R)
    return out


def part_b():
    fixed = []
    for R in range(24):
        for t in (np.zeros(3, int), -(RM[R] @ ONE)):
            if any(COMP[R, k] == R and np.array_equal(2 * t + RM[R] @ u, t) for k, u in CHILD):
                fixed.append((R, t))
    C.check(len(fixed) == 48, "48 poses are fixed seeds (one of their own children): (G,0) and (G,-G(1,1,1)) for each of the 24 G")
    cells = [{tuple(c) for c in M.place(RM[R], t)} for R, t in fixed]
    n = len(fixed)
    adj = {i: set() for i in range(n)}
    for a, b in itertools.combinations(range(n), 2):
        if cells[a] & cells[b]:
            continue
        ow = {c: 0 for c in cells[a]} | {c: 1 for c in cells[b]}
        ok = (0, 1) not in M.contacts(ow) or M.legal([(RM[fixed[a][0]], fixed[a][1]), (RM[fixed[b][0]], fixed[b][1])], 0, 1, ATLAS_KEYS)
        if ok:
            adj[a].add(b); adj[b].add(a)
    cl = []
    bron_kerbosch(adj, set(), set(range(n)), set(), cl)
    sizes = Counter(len(c) for c in cl)
    results = []
    for c in cl:
        r = np.array([fixed[i][0] for i in c]); t = np.array([fixed[i][1] for i in c])
        for _ in range(3):
            r, t = refine_vec(r, t)
        poses = [(RM[a], b) for a, b in zip(r, t)]
        st = M.patch_stats(poses, ATLAS_KEYS)
        if st["overlap"] or st["illegal"]:
            results.append({"seeds": len(c), "ok": False})
            continue
        own = st["owner"]
        cov = all((x, y, z) in own for x in range(-4, 4) for y in range(-4, 4) for z in range(-4, 4))
        entry = {"seeds": len(c), "ok": True, "cover": cov, "clique": c}
        if cov:
            g = window_symmetry(r, t, 4)
            best = max(g.values(), key=len) if g else []
            entry["order"] = len(best)
            entry["rots"] = best
        results.append(entry)
    cov = [e for e in results if e.get("cover")]
    C.check(all(e["ok"] for e in results), "every maximal compatible seed set refines to an overlap-free patch with E = 0 (level 3)")
    summary = Counter((e["seeds"], e["order"]) for e in cov)
    # types of the symmetry groups of the covering sets
    d4_sets = [e for e in cov if e["order"] == 8]
    C.check(len(d4_sets) == 3, "exactly 3 covering seed sets have a symmetry group of order 8")
    for e in d4_sets:
        ords = sorted(elt_order(i) for i in e["rots"])
        axes = sorted(involution_axis(RM[i]) for i in e["rots"] if elt_order(i) == 2)
        C.check(ords == [1, 2, 2, 2, 2, 2, 4, 4] and axes.count("coord") == 3 and axes.count("diag") == 2,
                "that group has the element orders and half-turn axes of D4 (identity, 5 half-turns, 2 quarter-turns)")
    # level 4 recheck of the covering sets, and the 'only the full sets cover' test on subsets of the 8-seed sets
    for e in cov:
        r = np.array([fixed[i][0] for i in e["clique"]]); t = np.array([fixed[i][1] for i in e["clique"]])
        for _ in range(4):
            r, t = refine_vec(r, t)
        g = window_symmetry(r, t, 8)
        best = max(g.values(), key=len) if g else []
        C.check(len(best) == e["order"], f"level 4: the symmetry order {e['order']} of the {e['seeds']}-seed set persists on a window of twice the radius")
    sub_cover = Counter()
    for e in d4_sets:
        for k in range(2, 9):
            for sub in itertools.combinations(e["clique"], k):
                r = np.array([fixed[i][0] for i in sub]); t = np.array([fixed[i][1] for i in sub])
                for _ in range(2):
                    r, t = refine_vec(r, t)
                cen = centres(r, t)
                st = M.patch_stats([(RM[a], b) for a, b in zip(r, t)], ATLAS_KEYS)
                own = st["owner"]
                # coverage test on the radius-2 window at level 2 (the window shrinks with the level)
                full = all((x, y, z) in own for x in range(-2, 2) for y in range(-2, 2) for z in range(-2, 2))
                sub_cover[(k, full)] += 1
    C.check(sub_cover[(8, True)] == 3 and all(not v_ or k_ == 8 for (k_, v_) in sub_cover), "on the level-2 window test, only the full 8-seed sets cover; no proper subset does")
    by_size = {k: {"covering": sum(1 for e in results if e["seeds"] == k and e.get("cover")), "total": v} for k, v in sorted(sizes.items())}
    return {"by_size": by_size, "fixed": len(fixed), "cliques": dict(sizes), "covering": {f"{k[0]} seeds, order {k[1]}": v for k, v in summary.items()},
            "total_cliques": len(cl), "subset_cover": {f"{k[0]}:{k[1]}": v for k, v in sub_cover.items()}}


# ----------------------------------------------------------------------------------------------
# C. Contact counts and extensivity defect
# ----------------------------------------------------------------------------------------------
def contacts_count(r, t):
    ks = np.sort(key(r, t))
    tot = 0
    for a, ta in ATLAS:
        r2 = COMP[r, a]
        t2 = np.einsum("nij,j->ni", RM[r], ta) + t
        k2 = key(r2, t2)
        pos = np.searchsorted(ks, k2)
        pos[pos == len(ks)] = len(ks) - 1
        tot += int((ks[pos] == k2).sum())
    return tot // 2  # the atlas is inverse-closed, so each contact is counted from both ends


def part_c(levels=7):
    r, t = np.array([I24]), np.zeros((1, 3), int)
    counts = {}
    for n in range(1, levels + 1):
        r, t = refine_vec(r, t)
        counts[n] = contacts_count(r, t)
    # cell-based cross-check where it is affordable
    poses = [(np.eye(3, dtype=int), (0, 0, 0))]
    for n in range(1, 5):
        poses = M.refine(poses, CHILD_M)
        st = M.patch_stats(poses, ATLAS_KEYS)
        C.check(st["illegal"] == 0 and st["contacts"] == counts[n], f"level {n}: pose-based and cell-based contact counts agree ({counts[n]})")
    import sympy as sp
    n_ = sp.symbols("n")
    a, b, c, d = sp.symbols("a b c d")
    f = a * 8**n_ + b * 4**n_ + c * 2**n_ + d
    sol = sp.solve([f.subs(n_, k) - counts[k] for k in range(1, 5)], [a, b, c, d])
    C.check(sol == {a: sp.Rational(29, 8), b: -3, c: 0, d: -1},
            "levels 1-4 fit C_n = (29/8) 8^n - 3 * 4^n - 1 with no 2^n term")
    for k in range(5, levels + 1):
        C.check(sp.Rational(29, 8) * 8**k - 3 * 4**k - 1 == counts[k], f"level {k}: the formula predicts the contact count {counts[k]}")
    for k in range(1, levels):
        C.check(8 * counts[k] - counts[k + 1] == -(12 * 4**k + 7), f"extensivity defect 8 C_{k} - C_{k+1} = -(12*4^{k} + 7)")
    return counts


def draw_c(axs, counts):
    ns = np.array(sorted(counts))
    Cn = np.array([counts[n] for n in ns], float)
    ax = axs[0]
    ax.semilogy(ns, Cn, "o", label="computed contacts $C_n$")
    ax.semilogy(ns, 29 / 8 * 8.0**ns - 3 * 4.0**ns - 1, "-", lw=1, label=r"$\frac{29}{8}8^n-3\cdot4^n-1$")
    ax.semilogy(ns, 29 / 8 * 8.0**ns, "--", lw=1, color="0.5", label=r"volume term $\frac{29}{8}8^n$")
    ax.set_xlabel("level $n$ (8$^n$ chairs)"); ax.set_ylabel("face contacts"); ax.legend(fontsize=7)
    ax.set_title("exact contact count of the hierarchy", fontsize=9)
    ax = axs[1]
    ks = np.arange(1, max(ns))
    defect = np.array([(12 * 4.0**k + 7) / counts[k + 1] for k in ks])
    ax.semilogy(ks, defect, "o-", label=r"chairs: $|8C_n-C_{n+1}|/C_{n+1}$")
    ax.semilogy(ks, 12 / 29 * 2.0**(-ks), "--", color="0.5", lw=1, label=r"$\frac{12}{29}2^{-n}$")
    kk = np.arange(0, 6)
    ax.semilogy(kk + 1, 2.0**(-4 * 2.0**kk) / 64, "s-", color="C3", ms=4, label=r"O(3): bound from $\Delta_{k+1}\leq64\Delta_k^2$, $\Delta_0=2^{-10}$")
    ax.set_ylim(1e-40, 3)
    ax.set_xlabel("level $n$ (chairs) or doubling number $k$ (O(3), shifted by 1)"); ax.set_ylabel("relative defect")
    ax.legend(fontsize=6.5, loc="lower left")
    ax.set_title("extensivity defect under doubling", fontsize=9)


# ----------------------------------------------------------------------------------------------
# D. Clock-style tests on the 8-seed fixed point
# ----------------------------------------------------------------------------------------------
class Tiling:
    def __init__(self, seeds, levels):
        r = np.array([s[0] for s in seeds]); t = np.array([s[1] for s in seeds])
        for _ in range(levels):
            r, t = refine_vec(r, t)
        self.r, self.t, self.N = r, t, levels
        self.cen = centres(r, t)
        self.tree = cKDTree(self.cen.astype(float))

    def signature(self, i, rad):
        nb = np.array(self.tree.query_ball_point(self.cen[i].astype(float), rad))
        rr = COMP[TR[self.r[i]], self.r[nb]]
        dt = np.einsum("ji,nj->ni", RM[self.r[i]], self.t[nb] - self.t[i])
        return np.sort(key(rr, dt))


def loo(sigs, labels):
    grp = defaultdict(Counter)
    for s, lab in zip(sigs, labels):
        grp[s][lab] += 1
    err = un = n = 0
    for s, lab in zip(sigs, labels):
        cn = grp[s]
        if sum(cn.values()) == 1:
            un += 1
            continue
        cn2 = cn.copy(); cn2[lab] -= 1; n += 1
        err += cn2.most_common(1)[0][0] != lab
    return err / max(n, 1), un / len(labels)


def part_d(seeds):
    core_r = 20 if FULL else 16
    T = Tiling(seeds, 6)
    core = np.where(np.abs(T.cen).max(1) <= core_r)[0]
    levels = {1: 8, 2: 64, 3: 512}
    labels = {L: [int(i % m) for i in core] for L, m in levels.items()}
    cache = {}

    def sigs(rad):
        if rad not in cache:
            cache[rad] = [hash(T.signature(i, rad).tobytes()) for i in core]
        return cache[rad]

    grid = [1, 2, 3, 4, 5, 6, 7, 8, 10, 12]
    curves = {L: [] for L in levels}
    unseen = {L: [] for L in levels}
    for rad in grid:
        sg = sigs(rad + 1e-9)
        for L in levels:
            e, u = loo(sg, labels[L])
            curves[L].append(e); unseen[L].append(u)
    # first radius with zero error and nothing unseen, by bisection over squared centre distances
    r2s = sorted({a * a + b * b + c * c for a in range(13) for b in range(13) for c in range(13) if 0 < a * a + b * b + c * c <= 144})
    thr = {}
    for L in levels:
        lo, hi = 0, len(r2s) - 1
        def good(j):
            e, u = loo(sigs(r2s[j] ** 0.5 + 1e-9), labels[L])
            return e == 0 and u == 0
        if not good(hi):
            thr[L] = None
            continue
        while lo < hi:
            mid = (lo + hi) // 2
            if good(mid):
                hi = mid
            else:
                lo = mid + 1
        thr[L] = r2s[lo] ** 0.5
        C.check(lo == 0 or not good(lo - 1), f"level {L}: the zero-error radius {thr[L]:.3f} is the first on the list")
    # look-back statistic.  Near the fixed point every chair shares the same large supertiles, so the pairs are drawn
    # from cubes of radius 12 centred away from the origin, and only pairs inside one level-N supertile are used.
    rng = np.random.default_rng(7)
    centres_ = [np.array(c) for c in ((40, 8, -12), (-36, 20, 14), (10, -42, 22), (24, 30, -34))]
    d_all, h_all, ho_all = [], [], []
    for cc in centres_:
        reg = np.where(np.abs(T.cen - cc).max(1) <= 12)[0]
        sub = rng.choice(reg, size=min(1500, len(reg)), replace=False)
        P = T.cen[sub].astype(float)
        D = np.sqrt(((P[:, None] - P[None]) ** 2).sum(-1))
        h = np.zeros(D.shape, int)
        A, B = sub[:, None], sub[None, :]
        for k in range(1, T.N + 1):
            h = np.where((h == 0) & ((A // 8**k) == (B // 8**k)), k, h)
        ho = np.zeros(D.shape, int)
        for k in range(1, 9):
            same = np.all(np.floor(P[:, None] / 2 ** (k + 1)) == np.floor(P[None] / 2 ** (k + 1)), axis=-1)
            ho = np.where((ho == 0) & same, k, ho)
        iu = np.triu_indices(len(sub), 1)
        ok = (h[iu] > 0) & (ho[iu] > 0)
        d_all.append(D[iu][ok]); h_all.append(h[iu][ok]); ho_all.append(ho[iu][ok])
    d, hh, hho = map(np.concatenate, (d_all, h_all, ho_all))
    sel = (d >= 4) & (d <= 16)
    x = np.log2(d[sel] / 2)
    slope_c = np.polyfit(x, hh[sel], 1)[0]
    slope_o = np.polyfit(x, hho[sel], 1)[0]
    lb = []
    for lo_, hi_ in [(2, 4), (4, 8), (8, 16), (16, 24)]:
        m = (d >= lo_) & (d < hi_)
        lb.append((lo_, hi_, int(m.sum()), float(hh[m].mean()), float(hho[m].mean())))
    return {"grid": grid, "curves": curves, "unseen": unseen, "thr": thr, "core": int(len(core)), "core_radius": core_r,
            "slope_chair": float(slope_c), "slope_octree": float(slope_o), "lookback": lb}


def draw_d(axs, d):
    colors = {1: "C0", 2: "C1", 3: "C2"}
    ax = axs[0]
    for L in (1, 2, 3):
        ax.plot(np.array(d["grid"]) / 2, d["curves"][str(L)] if isinstance(next(iter(d["curves"])), str) else d["curves"][L], "o-", ms=3, color=colors[L], label=f"level {L} slot")
    ax.set_xlabel(r"patch radius $r/\ell$  ($\ell$ = chair edge = 2 cubes)"); ax.set_ylabel("leave-one-out error")
    ax.set_title("how much of the tiling fixes your place in it", fontsize=9); ax.legend(fontsize=7)
    ax = axs[1]
    Ls = np.array([1, 2, 3])
    th = np.array([d["thr"][L] if L in d["thr"] else d["thr"][str(L)] for L in Ls]) / 2
    ax.semilogy(Ls, th, "o-", label="zero-error radius (Chair44, this fork)")
    ax.semilogy(Ls, th[0] * 2.0 ** (Ls - 1), "--", color="0.5", label=r"$2^{L-1}$ (inflation factor 2)")
    ax.semilogy(Ls, th[0] * 3.0 ** (Ls - 1), ":", color="C3", label=r"$3^{L-1}$ ([HR] report about 3 per level, Spectre)")
    ax.set_xticks(Ls); ax.set_xlabel("level $L$ of the supertile"); ax.set_ylabel(r"radius $/\ell$"); ax.legend(fontsize=6.5)
    ax.set_title("radius needed grows with the inflation factor", fontsize=9)
    ax = axs[2]
    lb = d["lookback"]
    xs = [np.log2(np.sqrt(a * b) / 2) for a, b, *_ in lb]
    ax.plot(xs, [q[3] for q in lb], "o-", label=f"Chair44 hierarchy (slope {d['slope_chair']:.2f})")
    ax.plot(xs, [q[4] for q in lb], "s--", label=f"octree control (slope {d['slope_octree']:.2f})")
    ax.set_xlabel(r"$\log_2(d/\ell)$, bin centres"); ax.set_ylabel(r"mean level $h$ of lowest common supertile")
    ax.set_title("look-back statistic", fontsize=9); ax.legend(fontsize=7)


# ----------------------------------------------------------------------------------------------
# E. Noisy patches and the Kesten-Stigum line
# ----------------------------------------------------------------------------------------------
def part_e(seeds):
    T = Tiling(seeds, 5)
    core = np.where(np.abs(T.cen).max(1) <= 10)[0]
    rng = np.random.default_rng(11)
    out = {}
    eps_grid = np.linspace(0, 0.9, 10)
    for L, m, rad in ((1, 8, 2.0), (2, 64, 5.0)):
        sigs = [T.signature(i, rad + 1e-9) for i in core]
        uni = np.unique(np.concatenate(sigs))
        col = {k: j for j, k in enumerate(uni)}
        def vec(s):
            v = np.zeros(len(uni), bool); v[[col[k] for k in s]] = True
            return v
        keyd = {}
        for s, i in zip(sigs, core):
            keyd.setdefault(s.tobytes(), [s, Counter()])[1][int(i % m)] += 1
        reps = list(keyd.values())
        S = np.array([vec(s) for s, _ in reps])           # distinct signatures
        sizes = S.sum(1)
        cnt = np.array([[c[l] for l in range(m)] for _, c in reps], float)  # signature x label counts
        acc = []
        test = rng.choice(len(core), size=min(1500, len(core)), replace=False)
        for eps in eps_grid:
            ok = 0
            for j in test:
                keep = rng.random(len(sigs[j])) >= eps
                q = vec(sigs[j][keep])
                comp = (S | ~q).all(1)                     # signatures containing the observed neighbours
                w = np.where(comp, eps ** (sizes - q.sum()) if eps > 0 else (sizes == q.sum()).astype(float), 0.0)
                post = (w[:, None] * cnt).sum(0)
                ok += int(np.argmax(post) == int(core[j] % m)) if post.sum() > 0 else 0
            acc.append(ok / len(test))
        out[L] = {"radius": rad, "signatures": len(reps), "accuracy": acc}
    out["eps"] = [float(e) for e in eps_grid]
    return out


def draw_e(axs, e):
    ax = axs[0]
    b = np.arange(2, 13)
    ax.plot(b, b ** -0.5, "-", label=r"$\lambda_c=b^{-1/2}$ ($b\lambda^2=1$)")
    ax.fill_between(b, b ** -0.5, 1, color="C0", alpha=0.12)
    ax.fill_between(b, 0, b ** -0.5, color="C3", alpha=0.10)
    ax.axvline(8, color="0.5", ls=":")
    ax.text(8.1, 0.82, "$b=8$ (chair\nsubstitution tree)", fontsize=7)
    ax.text(4, 0.55, "reconstruction\nof the root", fontsize=7, color="C0")
    ax.text(7, 0.08, "none", fontsize=7, color="C3")
    ax.set_xlabel("branching $b$"); ax.set_ylabel(r"channel parameter $\lambda$"); ax.set_ylim(0, 1)
    ax.set_title("entry 229: three-state channel, tree", fontsize=9); ax.legend(fontsize=7, loc="upper right")
    ax = axs[1]
    for L, c in ((1, "C0"), (2, "C1")):
        ax.plot(e["eps"], e[L]["accuracy"] if L in e else e[str(L)]["accuracy"], "o-", ms=3, color=c,
                label=f"level-{L} slot, radius {e[L]['radius'] / 2 if L in e else e[str(L)]['radius'] / 2:g}$\\ell$ (chance {1 / (8 if L == 1 else 64):.3f})")
    ax.set_xlabel("probability each neighbouring chair is deleted"); ax.set_ylabel("Bayes accuracy")
    ax.set_title("noisy-patch version of the clock test (toy)", fontsize=9); ax.legend(fontsize=6.5)


# ----------------------------------------------------------------------------------------------
def main() -> None:
    C.setup_style()
    import matplotlib.pyplot as plt
    t0 = time.time()
    payload = {}

    info, cover = part_a()
    pb = part_b()
    payload["B"] = pb
    print("B", pb, round(time.time() - t0), "s", flush=True)
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.6), gridspec_kw={"width_ratios": [1.25, 1]})
    draw_a(axs[0], info, cover)
    ax = axs[1]
    bs = pb["by_size"]
    ks = sorted(bs)
    cov = [bs[k]["covering"] for k in ks]; non = [bs[k]["total"] - bs[k]["covering"] for k in ks]
    xs_ = np.arange(len(ks))
    ax.bar(xs_, cov, color="#b8e0c8", edgecolor="k", label="fills the window around the origin")
    ax.bar(xs_, non, bottom=cov, color="#d9d9d9", edgecolor="k", label="leaves holes")
    ax.set_xticks(xs_); ax.set_xticklabels([f"{k} seeds" for k in ks]); ax.set_ylabel("maximal compatible seed sets")
    for x_, k in zip(xs_, ks):
        if bs[k]["covering"]:
            order = 8 if k == 8 else 1
            ax.text(x_, bs[k]["total"] + 2, f"symmetry order {order}", ha="center", fontsize=7)
    ax.legend(fontsize=7, loc="upper right")
    ax.set_title("fixed-point seed sets of the substitution\n(219 maximal sets of the 48 fixed poses)", fontsize=8.5)
    savefig(fig, "clock_fig1_symmetry")

    counts = part_c()
    payload["C"] = counts
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 4.0))
    draw_c(axs, counts)
    fig.tight_layout()
    savefig(fig, "clock_fig2_doubling")
    print("C done", round(time.time() - t0), "s", flush=True)

    # seeds: the first covering 8-seed set
    fixed = []
    for R in range(24):
        for t in (np.zeros(3, int), -(RM[R] @ ONE)):
            if any(COMP[R, k] == R and np.array_equal(2 * t + RM[R] @ u, t) for k, u in CHILD):
                fixed.append((R, t))
    # re-derive the covering 8-seed sets for use as seeds
    cells = [{tuple(c) for c in M.place(RM[R], t)} for R, t in fixed]
    adj = {i: set() for i in range(len(fixed))}
    for a, b in itertools.combinations(range(len(fixed)), 2):
        if cells[a] & cells[b]:
            continue
        ow = {c: 0 for c in cells[a]} | {c: 1 for c in cells[b]}
        if (0, 1) not in M.contacts(ow) or M.legal([(RM[fixed[a][0]], fixed[a][1]), (RM[fixed[b][0]], fixed[b][1])], 0, 1, ATLAS_KEYS):
            adj[a].add(b); adj[b].add(a)
    cl = []
    bron_kerbosch(adj, set(), set(range(len(fixed))), set(), cl)
    eight = [c for c in cl if len(c) == 8]
    seeds = [fixed[i] for i in eight[0]]

    pd = part_d(seeds)
    payload["D"] = pd
    print("D", {k: pd[k] for k in ("thr", "slope_chair", "slope_octree", "core")}, round(time.time() - t0), "s", flush=True)
    fig, axs = plt.subplots(1, 3, figsize=(13.5, 3.9))
    draw_d(axs, pd)
    fig.tight_layout()
    savefig(fig, "clock_fig3_clock_tests")

    pe = part_e(seeds)
    payload["E"] = pe
    fig, axs = plt.subplots(1, 2, figsize=(11, 3.9))
    draw_e(axs, pe)
    fig.tight_layout()
    savefig(fig, "clock_fig4_noisy_and_trees")
    C.save_json("hierarchy_figures", payload)
    print("done", round(time.time() - t0), "s")
    print("thresholds (r/l):", {L: v / 2 for L, v in pd["thr"].items()})
    print("lookback:", pd["lookback"], pd["slope_chair"], pd["slope_octree"])
    print("noisy:", {L: [round(a, 3) for a in pe[L]["accuracy"]] for L in (1, 2)})


if __name__ == "__main__":
    main()
