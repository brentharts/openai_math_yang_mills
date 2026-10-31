"""Root systems of the simply laced simple Lie algebras from their Dynkin diagrams.

The closure of the simple roots under the simple reflections, in simple-root coordinates. Used by
the exploration scripts; models/chirality_lisi.py carries its own copy for E8 so that the model
stays self-contained.
"""
from __future__ import annotations

from collections import Counter

DIAGRAMS = {   # name -> (rank, edges) in Bourbaki numbering
    "A": lambda n: (n, [(i, i + 1) for i in range(1, n)]),
    "D": lambda n: (n, [(i, i + 1) for i in range(1, n - 1)] + [(n - 2, n)]),
    "E6": lambda: (6, [(1, 3), (3, 4), (4, 2), (4, 5), (5, 6)]),
    "E7": lambda: (7, [(1, 3), (3, 4), (4, 2), (4, 5), (5, 6), (6, 7)]),
    "E8": lambda: (8, [(1, 3), (3, 4), (4, 2), (4, 5), (5, 6), (6, 7), (7, 8)]),
}


def cartan(rank, edges):
    a = [[2 if i == j else 0 for j in range(rank)] for i in range(rank)]
    for i, j in edges:
        a[i - 1][j - 1] = a[j - 1][i - 1] = -1
    return a


class RootSystem:
    def __init__(self, rank, edges):
        self.rank, self.edges = rank, edges
        self.A = cartan(rank, edges)
        n = rank
        simple = [tuple(1 if j == i else 0 for j in range(n)) for i in range(n)]
        found, frontier = set(simple), list(simple)
        while frontier:
            new = []
            for beta in frontier:
                for i in range(n):
                    pairing = sum(beta[j] * self.A[j][i] for j in range(n))
                    r = tuple(b - (pairing if k == i else 0) for k, b in enumerate(beta))
                    if r not in found:
                        found.add(r)
                        new.append(r)
            frontier = new
        self.roots = found | {tuple(-x for x in r) for r in found}

    @classmethod
    def of(cls, kind, n=None):
        spec = DIAGRAMS[kind]
        return cls(*(spec(n) if n is not None else spec()))

    @property
    def dim(self):
        return len(self.roots) + self.rank

    def ip(self, c, d):
        n = self.rank
        return sum(c[i] * self.A[i][j] * d[j] for i in range(n) for j in range(n))

    def grading(self, node, modulus=None):
        """Counts of roots by the coefficient of one simple root (mod `modulus` if given)."""
        return Counter((r[node - 1] % modulus) if modulus else r[node - 1] for r in self.roots)

    def marks(self):
        return [max(r[k] for r in self.roots) for k in range(self.rank)]
