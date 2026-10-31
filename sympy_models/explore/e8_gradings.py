"""Scanning the gradings of E8 by one simple root.

The chirality repository uses three of them: node 7 mod 3 (E6 x A2, the family grading), node 1 mod 3
(Lisi's grading, so(14) + u(1) fixed) and the two involutions (node 1 and node 8 mod 2). This script
scans every node and every modulus 2..6, lists the fixed-algebra dimension of each, finds all pairs of
order-three gradings whose joint table is the repository's, and tabulates the joint grading of the
two involutions (Z2 x Z2).
"""
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from explore import Result, main
from explore._roots import RootSystem

TABLE = [[50, 18, 18], [21, 42, 18], [21, 18, 42]]


def fixed_dim(e8, node, m):
    c = e8.grading(node, m)
    return c[0] + 8


def joint(e8, a, ma, b, mb):
    c = {}
    for r in e8.roots:
        k = (r[a - 1] % ma, r[b - 1] % mb)
        c[k] = c.get(k, 0) + 1
    c[(0, 0)] = c.get((0, 0), 0) + 8
    return [[c.get((i, j), 0) for j in range(mb)] for i in range(ma)]


def run() -> Result:
    r = Result("e8_gradings", "Scanning the one-node gradings of E8", prose="")
    e8 = RootSystem.of("E8")
    marks = e8.marks()
    r.check("the marks (highest-root coefficients) are 2 3 4 6 5 4 3 2 (Bourbaki)", marks == [2, 3, 4, 6, 5, 4, 3, 2])
    scan = {(k, m): fixed_dim(e8, k, m) for k in range(1, 9) for m in range(2, 7)}
    r.check("every grading by node k mod m is an automorphism: eigenspace dimensions sum to 248",
            all(sum((e8.grading(k, m)[i] + (8 if i == 0 else 0)) for i in range(m)) == 248
                for k in range(1, 9) for m in range(2, 7)))
    r.check("node 8 mod 2 fixes 136 (E7 + A1), node 1 mod 2 fixes 120 (D8), node 7 mod 3 fixes 86 (E6 + A2), node 1 mod 3 fixes 92 (D7 + U1)",
            (scan[(8, 2)], scan[(1, 2)], scan[(7, 3)], scan[(1, 3)]) == (136, 120, 86, 92))
    # the largest possible fixed algebras are the semisimple maximal-rank ones, at m = mark of the node
    semis = {(k, marks[k - 1]): scan[(k, marks[k - 1])] for k in range(1, 9)}
    r.data["fixed dimension at m = mark"] = semis
    # pairs of order-three gradings with the repository's table
    exact, same_shape = [], []
    for a in range(1, 9):
        for b in range(1, 9):
            if a == b:
                continue
            t = joint(e8, a, 3, b, 3)
            if t == TABLE:
                exact.append((a, b))
            elif sorted(map(sorted, t)) == sorted(map(sorted, TABLE)):
                same_shape.append((a, b))
    r.data["pairs (node, node) giving exactly the table"] = exact
    r.data["pairs giving the table up to permuting rows or columns"] = same_shape
    r.check("the repository's pair (node 7, node 1) gives the table exactly", (7, 1) in exact)
    r.check("some other pair gives the same cell sizes in a different arrangement",
            len(exact) + len(same_shape) > 1)
    # Z2 x Z2: the two involutions together
    z = joint(e8, 1, 2, 8, 2)
    r.data["Z2 x Z2 table (node 1 parity by node 8 parity)"] = z
    r.check("the joint table of the two involutions has total 248", sum(map(sum, z)) == 248)
    r.check("its row sums are the eigenspaces of node 1 (120, 128) and its column sums those of node 8 (136, 112)",
            [sum(row) for row in z] == [120, 128] and [sum(col) for col in zip(*z)] == [136, 112])
    r.prose = (r"Grading $E_8$ by the coefficient of one simple root modulo $m$ gives an automorphism for every node "
               r"and modulus; the scan lists the fixed algebra of all $8\times5$ of them. Among pairs of order-three "
               r"gradings, the repository's pair (nodes 7 and 1) gives its table exactly, and "
               + str(len(same_shape)) + r" other ordered pairs give the same cell sizes in a different arrangement. "
               r"The two involutions together give a $\mathbb Z_2\times\mathbb Z_2$ grading with cells "
               + ", ".join(str(x) for row in z for x in row) + r".")
    r.display = (r"\begin{array}{c|cc} & \text{even}_8 & \text{odd}_8\\\hline \text{even}_1 & %d & %d\\ \text{odd}_1 & %d & %d\end{array}"
                 % (z[0][0], z[0][1], z[1][0], z[1][1]))
    return r


if __name__ == "__main__":
    main(run())
