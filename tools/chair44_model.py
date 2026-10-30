#!/usr/bin/env python3
"""Pose-level constraint model built from the Chair44 contact atlas (not a repository entry).

Source of the data (external, not part of the OpenAI release):
    Tsiokos, arXiv:2609.19214, Figure 7 (the 44 legal contact poses) and Table 1 (the eight child poses of
    the doubled chair).  The numbers below were transcribed from a text extraction of that paper; they are
    not copied from a machine-readable file, so the checks in this script double as a transcription test.

What this script is -- and is not
---------------------------------
Flicker's paper proposes a 24-state cubic-lattice model (24 = order of the proper cubic rotation group)
with nearest-neighbour penalties for illegal contacts.  Its penalty table is NOT reproduced here.  This
script builds the closest model that the Tsiokos atlas supports by itself:

    configuration : a set of chairs, each with a state (G, t), G one of the 24 proper cubic rotations
                    and t an integer translation, acting on the unit chair by x -> G x + t;
    energy        : E = (number of overlapping cells) + (number of face-contacting chair pairs whose
                    relative pose is not in the 44-pose atlas, in either direction);
    ground states : E = 0.

The atlas is taken as the definition of "legal contact".  The pyramid-geometry test that produced the atlas
(Tsiokos section 3, panel labels in Figure 3) is not re-derived: those labels were not in the extracted text.

Pose convention (inferred, then verified)
-----------------------------------------
The extracted text gives frames (p)(s) with G(x)_i = s_i x_{p_i} and offsets, but not the unit chair or the
action explicitly.  We take the chair to be the seven unit cubes a + [0,1]^3, a in {0,1}^3 \\ {(1,1,1)}, and the
pose to act as x -> G x + t about the origin.  Among the conventions tried (unit or doubled cubes, rotation
about the origin or the chair centre, G or G^T) this is the only one under which all 44 atlas poses are
overlap-free and in face contact with the root chair; check (1) below asserts that.

Checks (each stops the build if it fails)
-----------------------------------------
  1. the 44 poses are distinct, have det G = +1, and each is overlap-free and face-touching the root;
  2. the 8 child poses of Table 1 fill the doubled chair exactly (56 cells, no overlap);
  3. the refinement rule (G,t) -> {(GH, 2t+Gu)} applied n times gives overlap-free patches of 8^n chairs
     whose union is the 2^n-times scaled chair, and every face contact is an atlas contact (E = 0).
Results that are only reported, not asserted: how many atlas poses are realised, how many absolute
rotation states occur, and how far zero-temperature random growth gets before it is stuck.
"""

from __future__ import annotations

import itertools
import sys
from collections import Counter, defaultdict

import numpy as np
import sympy as sp

import common as C

# ----------------------------------------------------------------------------------------------
# Data transcribed from Tsiokos arXiv:2609.19214 (Figure 7 and Table 1)
# ----------------------------------------------------------------------------------------------
# (permutation p, signs s, offset): G(x)_i = s_i * x_{p_i}
ATLAS_RAW = """
012 --+ 2,2,-2 | 012 --+ 2,2,2 | 012 --+ 3,3,-1 | 012 --+ 3,3,1 | 012 -+- 2,-2,2 | 012 -+- 2,2,2
012 +-- -2,2,2 | 012 +-- 2,2,2 | 012 +++ -1,-1,-1 | 012 +++ 1,1,1
021 -++ 4,0,0 | 021 +-+ -1,3,-1 | 021 +-+ 0,4,0 | 021 ++- 0,0,4 | 021 ++- 1,1,3
102 -++ 0,0,0 | 102 -++ 4,0,0 | 102 +-+ 0,0,0 | 102 +-+ 0,4,0 | 102 ++- -1,-1,3 | 102 ++- 0,0,0
102 ++- 0,0,4 | 102 ++- 1,1,3
120 --+ 2,2,-2 | 120 --+ 2,2,2 | 120 --+ 3,3,1 | 120 -+- 2,-2,2 | 120 -+- 2,2,2 | 120 -+- 3,-1,3
120 +-- -2,2,2 | 120 +-- 2,2,2
201 --+ 2,2,-2 | 201 --+ 2,2,2 | 201 --+ 3,3,1 | 201 -+- 2,-2,2 | 201 -+- 2,2,2 | 201 +-- -2,2,2
201 +-- -1,3,3 | 201 +-- 2,2,2
210 -++ 3,-1,-1 | 210 -++ 4,0,0 | 210 +-+ 0,4,0 | 210 ++- 0,0,4 | 210 ++- 1,1,3
"""
# Table 1: the eight child poses of the doubled chair (p, s, u)
CHILDREN_RAW = [
    ("012", "+++", (0, 0, 0)), ("102", "++-", (0, 0, 4)), ("021", "+-+", (0, 4, 0)), ("201", "+--", (0, 4, 4)),
    ("210", "-++", (4, 0, 0)), ("120", "-+-", (4, 0, 4)), ("012", "--+", (4, 4, 0)), ("012", "+++", (1, 1, 1)),
]


def frame(p: str, s: str) -> np.ndarray:
    """3x3 integer matrix of G(x)_i = s_i x_{p_i}."""
    M = np.zeros((3, 3), dtype=int)
    for i in range(3):
        M[i, int(p[i])] = 1 if s[i] == "+" else -1
    return M


def parse_atlas() -> list[tuple[np.ndarray, tuple[int, int, int]]]:
    out = []
    for item in ATLAS_RAW.replace("\n", "|").split("|"):
        item = item.strip()
        if not item:
            continue
        p, s, o = item.split()
        out.append((frame(p, s), tuple(int(x) for x in o.split(","))))
    return out


def all_rotations() -> list[np.ndarray]:
    rots = []
    for p in itertools.permutations("012"):
        for s in itertools.product("+-", repeat=3):
            M = frame("".join(p), "".join(s))
            if round(np.linalg.det(M)) == 1:
                rots.append(M)
    return rots


ROTS = all_rotations()
ROT_INDEX = {M.tobytes(): i for i, M in enumerate(ROTS)}


def pkey(G: np.ndarray, t) -> tuple:
    return (G.tobytes(), tuple(int(x) for x in t))


# ----------------------------------------------------------------------------------------------
# Geometry on the integer cell lattice
# ----------------------------------------------------------------------------------------------
def chair_cells(k: int) -> np.ndarray:
    """Cells (unit cubes, lower corners) of the chair built from cubes of side k."""
    cells = [
        tuple(a[i] * k + d[i] for i in range(3))
        for a in itertools.product((0, 1), repeat=3) if a != (1, 1, 1)
        for d in itertools.product(range(k), repeat=3)
    ]
    return np.array(sorted(cells), dtype=int)


UNIT = chair_cells(1)


def place(G: np.ndarray, t) -> np.ndarray:
    """Cells of the unit chair under x -> G x + t.  Work with doubled cell centres so rotations stay integral."""
    Y = (UNIT * 2 + 1) @ G.T + 2 * np.asarray(t)
    return (Y - 1) // 2


def contacts(owner: dict) -> set[tuple[int, int]]:
    """Unordered pairs of chairs that share a face, given a map cell -> chair index."""
    pairs = set()
    for c, i in owner.items():
        for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            j = owner.get((c[0] + d[0], c[1] + d[1], c[2] + d[2]))
            if j is not None and j != i:
                pairs.add((min(i, j), max(i, j)))
    return pairs


def relative(poses, i, j):
    """Pose of chair j in the frame of chair i: (G_i^T G_j, G_i^T (t_j - t_i))."""
    Gi, ti = poses[i]
    Gj, tj = poses[j]
    return pkey(Gi.T @ Gj, Gi.T @ (np.asarray(tj) - np.asarray(ti)))


def legal(poses, i, j, atlas_keys) -> bool:
    return relative(poses, i, j) in atlas_keys or relative(poses, j, i) in atlas_keys


# ----------------------------------------------------------------------------------------------
# Substitution patches
# ----------------------------------------------------------------------------------------------
def refine(poses, children):
    out = []
    for G, t in poses:
        for H, u in children:
            out.append((G @ H, 2 * np.asarray(t) + G @ np.asarray(u)))
    return out


def patch_stats(poses, atlas_keys):
    owner, overlap = {}, 0
    for i, (G, t) in enumerate(poses):
        for c in map(tuple, place(G, t)):
            if c in owner:
                overlap += 1
            owner[c] = i
    pairs = contacts(owner)
    bad = [(i, j) for i, j in pairs if not legal(poses, i, j, atlas_keys)]
    # both orientations of every contact are recorded (the atlas is inverse-closed, see main())
    used = {k for i, j in pairs for k in (relative(poses, i, j), relative(poses, j, i))}
    rots = Counter(ROT_INDEX[G.tobytes()] for G, _ in poses)
    return {"chairs": len(poses), "cells": len(owner), "overlap": overlap, "contacts": len(pairs),
            "illegal": len(bad), "atlas_used": len(used & atlas_keys), "rot_counts": rots, "owner": owner}


# ----------------------------------------------------------------------------------------------
# Zero-temperature growth: add chairs at atlas poses, never creating overlap or an illegal contact
# ----------------------------------------------------------------------------------------------
def grow(atlas, atlas_keys, n_max, rng):
    """Random sequential growth keeping E = 0.  Exact stuck test: a rejected candidate stays rejected
    (overlaps and contacts only accumulate), so a candidate pool that empties means no legal addition exists
    at any atlas pose of any placed chair.  Returns (number placed, stuck?)."""
    poses = [(np.eye(3, dtype=int), (0, 0, 0))]
    owner = {tuple(c): 0 for c in place(*poses[0])}
    seen = {pkey(*poses[0])}
    pool = []

    def push(idx):
        G, t = poses[idx]
        for Ga, ta in atlas:
            g2, t2 = G @ Ga, G @ np.asarray(ta) + np.asarray(t)
            k = pkey(g2, t2)
            if k not in seen:
                seen.add(k)
                pool.append((g2, t2))

    push(0)
    while pool and len(poses) < n_max:
        G, t = pool.pop(int(rng.integers(len(pool))))
        cells = [tuple(c) for c in place(G, t)]
        if any(c in owner for c in cells):
            continue
        # neighbouring chairs touched by this one
        nb = set()
        for c in cells:
            for d in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
                j = owner.get((c[0] + d[0], c[1] + d[1], c[2] + d[2]))
                if j is not None:
                    nb.add(j)
        trial = poses + [(G, t)]
        new = len(poses)
        if not all(legal(trial, j, new, atlas_keys) for j in nb):
            continue
        poses.append((G, t))
        for c in cells:
            owner[c] = new
        push(new)
    return len(poses), (not pool and len(poses) < n_max), poses


# ----------------------------------------------------------------------------------------------
def main() -> None:
    full = "--full" in sys.argv
    C.setup_style()
    atlas = parse_atlas()
    atlas_keys = {pkey(G, t) for G, t in atlas}
    inv_keys = {pkey(G.T, -(G.T @ np.asarray(t))) for G, t in atlas}

    # (1) the atlas itself ------------------------------------------------------------------
    C.check(len(atlas) == 44, "atlas has 44 poses")
    C.check(len(atlas_keys) == 44, "the 44 atlas poses are distinct")
    C.check(len(ROTS) == 24, "there are 24 proper cubic rotations (signed permutation matrices with det +1)")
    dets = {int(round(np.linalg.det(G))) for G, _ in atlas}
    C.check(dets == {1}, "all atlas frames are proper rotations (det +1)")
    root_cells = {tuple(c) for c in UNIT}
    n_free = n_touch = 0
    for G, t in atlas:
        pc = {tuple(c) for c in place(G, t)}
        n_free += not (pc & root_cells)
        owner = {c: 0 for c in root_cells} | {c: 1 for c in pc}
        n_touch += (0, 1) in contacts(owner)
    C.check(n_free == 44, "all 44 atlas poses are overlap-free with the root chair")
    C.check(n_touch == 44, "all 44 atlas poses are in face contact with the root chair")
    rel_rots = {ROT_INDEX[G.tobytes()] for G, _ in atlas}
    inv_closed = len(atlas_keys & inv_keys)

    # (2) Table 1: eight children fill the doubled chair ----------------------------------------
    children = [(frame(p, s), np.array(u)) for p, s, u in CHILDREN_RAW]
    cc = Counter(tuple(c) for G, u in children for c in place(G, u))
    doubled = {tuple(c) for c in chair_cells(2)}
    C.check(set(cc) == doubled and max(cc.values()) == 1,
            "the 8 child chairs fill the doubled chair exactly (56 cells, no overlap)")

    # (3) iterated refinement ---------------------------------------------------------------
    levels = 4 if full else 3
    stats = {}
    poses = [(np.eye(3, dtype=int), (0, 0, 0))]
    for n in range(1, levels + 1):
        poses = refine(poses, children)
        st = patch_stats(poses, atlas_keys)
        stats[n] = st
        C.check(st["chairs"] == 8**n and st["overlap"] == 0, f"level {n}: {8**n} chairs, no overlapping cells")
        C.check(st["illegal"] == 0,
                f"level {n}: all {st['contacts']} face contacts are atlas contacts (E = 0)")
        if n <= 3:
            C.check({c for c in st["owner"]} == {tuple(c) for c in chair_cells(2**n)},
                    f"level {n}: the patch fills the {2**n}x-scaled chair")

    # growth ------------------------------------------------------------------------------
    n_runs, n_max = (60, 600) if full else (24, 300)
    sizes, stuck_flags = [], []
    rot_growth = Counter()
    for seed in range(n_runs):
        rng = np.random.default_rng(1000 + seed)
        n, stuck, ps = grow(atlas, atlas_keys, n_max, rng)
        sizes.append(n)
        stuck_flags.append(stuck)
        if seed == 0:
            rot_growth = Counter(ROT_INDEX[G.tobytes()] for G, _ in ps)
            own = {}
            for i, (G, t) in enumerate(ps):
                for c in map(tuple, place(G, t)):
                    own[c] = i
            first = {}
            for i, j in contacts(own):
                for k in (relative(ps, i, j), relative(ps, j, i)):
                    if k in atlas_keys:
                        first[k] = min(first.get(k, 10**9), max(i, j))
            growth_used = len(first)
            growth_curve = [sum(1 for v in first.values() if v <= m) for m in range(len(ps))]
    sizes = np.array(sizes)

    # symbolic note: group order via sympy permutation group of the 24 rotations acting on the cube's 4 diagonals
    from sympy.combinatorics import Permutation, PermutationGroup
    # rotations about axes: generators acting on the 4 body diagonals
    g1 = Permutation([1, 2, 3, 0])   # 4-fold rotation permutes diagonals as a 4-cycle
    g2 = Permutation([1, 0, 2, 3])   # a 2-fold rotation swaps two diagonals
    order = PermutationGroup([g1, g2]).order()
    C.check(order == 24, "sympy: the group generated on the 4 body diagonals has order 24 (S_4 = cubic rotations)")

    top = stats[levels]
    rot_hist = [top["rot_counts"].get(i, 0) for i in range(24)]
    payload = {
        "atlas_poses": 44,
        "distinct_relative_rotations": len(rel_rots),
        "atlas_inverse_closed": inv_closed,
        "levels": {str(n): {"chairs": s["chairs"], "contacts": s["contacts"], "illegal": s["illegal"],
                            "atlas_used": s["atlas_used"],
                            "rotations_present": sum(1 for v in s["rot_counts"].values() if v)}
                   for n, s in stats.items()},
        "rotation_histogram_top_level": rot_hist,
        "growth": {"runs": n_runs, "n_max": n_max, "sizes": [int(x) for x in sizes],
                   "stuck": [bool(x) for x in stuck_flags],
                   "run0_rotations_present": len(rot_growth), "run0_atlas_used": growth_used, "run0_curve": growth_curve},
    }
    C.save_json("chair44_model", payload)

    # figure --------------------------------------------------------------------------------
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(11, 3.8))
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    grid = np.zeros((4, 4, 4), dtype=bool)
    colors = np.empty((4, 4, 4), dtype=object)
    cmap = plt.get_cmap("tab10")
    for i, (G, u) in enumerate(children):
        for c in place(G, u):
            grid[tuple(c)] = True
            colors[tuple(c)] = cmap(i)
    ax.voxels(grid, facecolors=colors, edgecolor="k", linewidth=0.2)
    ax.set_title("8 children fill the doubled chair", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.bar(range(24), rot_hist, color="C0")
    ax2.set_xlabel("rotation state (of 24)")
    ax2.set_ylabel(f"chairs at level {levels}")
    ax2.set_title("state occupancy of the substitution patch", fontsize=9)
    ax3 = fig.add_subplot(1, 3, 3)
    ax3.plot(growth_curve, color="C1", label="random growth, run 0")
    ax3.axhline(stats[levels]["atlas_used"], color="C0", ls="--", label="substitution patch")
    ax3.axhline(44, color="k", ls=":", label="all 44")
    ax3.set_xlabel("chairs placed")
    ax3.set_ylabel("atlas poses realised")
    ax3.legend(fontsize=7, loc="lower right")
    ax3.set_title(f"zero-temperature growth (cap {n_max})", fontsize=9)
    fig.tight_layout()
    C.savefig(fig, "chair44_model")

    # LaTeX ---------------------------------------------------------------------------------
    rows = [[str(n), str(s["chairs"]), str(s["contacts"]), str(s["illegal"]), str(s["atlas_used"]),
             str(sum(1 for v in s["rot_counts"].values() if v))] for n, s in stats.items()]
    stuck_n = int(sum(stuck_flags))
    tex = r"""
\section{A pose-level model from the Chair44 atlas (external data)}
\label{sec:chair44}

\paragraph{What this is.} Flicker's paper proposes a 24-state cubic-lattice model with nearest-neighbour penalties for illegal contacts whose ground states
should be Chair tilings. His penalty table is not reproduced here. This section builds the closest model that the contact atlas of
Tsiokos (arXiv:2609.19214, Figure~7 and Table~1) supports by itself. A state is a rotation $G$ among the 24 proper cubic rotations together with an integer
translation $t$; the chair is $\bigcup_{a\in\{0,1\}^3\setminus\{(1,1,1)\}}(a+[0,1]^3)$ acted on by $x\mapsto Gx+t$. The energy is
""" + C.equation(r"E=\#\{\text{overlapping cells}\}+\#\{\text{face-contacting pairs whose relative pose is not in the 44-pose atlas}\},", "eq:chair44-energy") + r"""
and the ground states are the configurations with $E=0$. The atlas is taken as the definition of a legal contact; the pyramid-geometry test that produced it is not re-derived.

\paragraph{Data and convention.} The 44 poses and the eight child poses were transcribed from a text extraction of the paper. The action $x\mapsto Gx+t$ on the unit chair was not stated
there; it was chosen because it is the only convention among those tried (unit or doubled cubes, rotation about the origin or the chair centre, $G$ or $G^{T}$) under which all 44 poses are
overlap-free and in face contact with the root chair, which the script asserts. The atlas contains """ + str(len(rel_rots)) + r""" distinct rotations out of 24, and """ + str(inv_closed) + r""" of its 44 poses have their inverse in the atlas.

\paragraph{Substitution check.} The refinement rule $(G,t)\mapsto\{(GH,\,2t+Gu)\}$ over the eight child poses $(H,u)$ was iterated. The script asserts that the eight children fill the doubled chair exactly, that each
level is overlap-free, that level $n$ fills the $2^n$-times scaled chair (levels up to 3), and that every face contact in the patch is an atlas contact. Hence $E=0$ for the substitution patches.
""" + C.table(
        ["level", "chairs", "face contacts", "illegal", "atlas poses used", "rotations present"], rows,
        caption=r"Substitution patches. `Illegal' counts face contacts whose relative pose is not one of the 44 atlas poses (zero at every level).",
        colspec="rrrrrr",
    ) + C.figure(
        "chair44_model",
        rf"Left: the eight children of Table~1 filling the doubled chair. Middle: occupancy of the 24 absolute rotation states in the level-{levels} patch. "
        rf"Right: atlas poses realised against chairs placed in one zero-temperature growth run (cap {n_max}), against the substitution patch ({stats[levels]['atlas_used']}) and all 44.",
        "\\linewidth",
    ) + r"""
\paragraph{Zero-temperature growth.} Starting from one chair, a chair is added at an atlas pose of an existing chair whenever it overlaps nothing and every face it touches is an atlas contact.
A rejected candidate stays rejected, since overlaps and contacts only accumulate, so an empty candidate pool is an exact dead end for this rule. Of """ + str(n_runs) + r""" runs with a cap of """ + str(n_max) + r""" chairs,
""" + str(stuck_n) + r""" stopped early (median size """ + f"{int(np.median(sizes))}" + r"""). In the first run the patch uses """ + str(len(rot_growth)) + r""" of the 24 rotation states and """ + str(growth_used) + r""" of the 44 atlas poses.
The substitution patches realise """ + str(stats[levels]["atlas_used"]) + r""" of the 44 poses at level """ + str(levels) + r""" (both orientations of each contact counted); Hartshorn and Romanov (Research Square, unrefereed) state in their Corollary~1 that only 30 of the 44 contacts occur, which this count matches.
The growth runs use more poses than the hierarchical patches do, so the pairwise atlas rule alone does not confine local growth to the hierarchical structure, at least up to the cap.

\paragraph{Caveat.} This tests the consistency of a transcription and a convention, and shows that the substitution patches are ground states of the energy above. It does not test that the atlas is
the complete set of legal contacts, does not use the face-colour or arrow rules of Flicker, has no temperature and no thermodynamic limit, and says nothing about whether Flicker's 24-state lattice model has Chair tilings
as its ground states. A growth run that stops early only shows a dead end for random local growth; it does not bound what a global search could place.
"""
    C.write_tex("chair44", tex)
    print("wrote build/tex/chair44.tex")
    print("levels:", {n: (s["chairs"], s["contacts"], s["illegal"], s["atlas_used"]) for n, s in stats.items()})
    print("rotations in atlas:", len(rel_rots), "inverse-closed:", inv_closed)
    print("growth sizes:", sorted(sizes.tolist()), "stuck:", stuck_n, "of", n_runs)


if __name__ == "__main__":
    main()
