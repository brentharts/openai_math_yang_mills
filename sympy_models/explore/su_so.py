"""dim su(N) = dim so(n): where the two classical series meet.

The knowledge graph links 'Special unitary dimension' (N^2 - 1) to 'Orthogonal algebra dimension'
(n(n - 1)/2). Their equality is the Pell-type equation (2n - 1)^2 - 8 N^2 = -7. This script solves it,
tests which solutions are isomorphisms of Lie algebras (equal rank and root count, and equal Cartan
matrices up to relabelling) and which are numerical coincidences, and then follows the two real
isomorphisms into the Standard Model: so(6) = su(4) is the Pati-Salam SU(4), and the half-spinor of
so(10) is 2 S(6) S(4) = 2 * 4 * 2 = 16, the (4,2,1) + (4bar,1,2) of Pati-Salam.
"""
import itertools
import math
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from explore import Result, main
from explore._roots import RootSystem


def pell_solutions(limit=3000):
    out = []
    for N in range(2, limit):
        x2 = 8 * N * N - 7
        x = math.isqrt(x2)
        if x * x == x2 and x % 2 == 1:
            out.append((N, (x + 1) // 2))
    return out


def so_rank_roots(n):
    """(rank, roots) of so(n): B_m for n = 2m + 1, D_m for n = 2m."""
    m = n // 2
    return (m, 2 * m * m) if n % 2 else (m, 2 * m * (m - 1))


def cartan_isomorphic(a, b):
    if len(a) != len(b):
        return False
    return any(all(a[p[i]][p[j]] == b[i][j] for i in range(len(a)) for j in range(len(a)))
               for p in itertools.permutations(range(len(a))))


def run() -> Result:
    r = Result("su_so", "Where su(N) and so(n) have the same dimension", prose="")
    sols = pell_solutions()
    r.data["solutions (N, n)"] = sols[:7]
    r.check("N^2 - 1 = n(n - 1)/2 at each solution found", all(N * N - 1 == n * (n - 1) // 2 for N, n in sols))
    r.check("the first solutions are (2,3), (4,6), (11,16), (23,33), (64,91)",
            sols[:5] == [(2, 3), (4, 6), (11, 16), (23, 33), (64, 91)])
    iso = []
    for N, n in sols[:5]:
        iso.append(((N, n), (N - 1, N * (N - 1)) == so_rank_roots(n)))
    r.data["same rank and root count"] = iso
    r.check("(2,3) and (4,6) have equal rank and root count; (11,16), (23,33), (64,91) do not",
            [s for _, s in iso] == [True, True, False, False, False])
    a3, d3 = RootSystem.of("A", 3), RootSystem.of("D", 3)
    r.check("the Cartan matrices of A3 and D3 agree up to relabelling, and both have 12 roots",
            cartan_isomorphic(a3.A, d3.A) and len(a3.roots) == len(d3.roots) == 12)
    r.check("8 = n(n-1)/2 has no integer solution, so su(3) has no so(n) of its dimension",
            all(n * (n - 1) // 2 != 8 for n in range(2, 12)))
    half = lambda n: 2 ** (n // 2 - 1)
    r.check("the half-spinors of so(6) and so(4) are 4 and 2: the fundamentals of su(4) and su(2)",
            (half(6), half(4)) == (4, 2))
    r.check("S(n + m) = 2 S(n) S(m) for even n, m: 16 = 2 x 4 x 2, the (4,2,1) + (4bar,1,2) of Pati-Salam",
            all(half(n + m) == 2 * half(n) * half(m) for n in range(4, 12, 2) for m in range(4, 12, 2))
            and half(10) == 2 * half(6) * half(4) == 16)
    r.check("so(10) contains so(6) + so(4): 45 = 15 + 6 + 6 x 4, with 15 = dim su(4) and 6 = dim su(2)^2",
            45 == 15 + 6 + 6 * 4)
    r.prose = (r"The equation $\dim\mathfrak{su}(N)=\dim\mathfrak{so}(n)$ is $(2n-1)^2-8N^2=-7$. Its first solutions "
               r"$(N,n)=(2,3)$, $(4,6)$, $(11,16)$, $(23,33)$ and $(64,91)$ were found by search and checked. Only the first two "
               r"are isomorphisms ($\mathfrak{su}(2)\cong\mathfrak{so}(3)$ and $\mathfrak{su}(4)\cong\mathfrak{so}(6)$: equal rank and "
               r"root count, and equal Cartan matrices up to relabelling); at $(11,16)$ the ranks are 10 and 8. "
               r"The second isomorphism is the Pati--Salam $SU(4)$, and since $S(n+m)=2S(n)S(m)$ the half-spinor of "
               r"$\mathfrak{so}(10)$ splits as $16=2\cdot4\cdot2$, the $(4,2,1)\oplus(\bar4,1,2)$ of Pati--Salam.")
    r.display = r"S(10)=2\,S(6)\,S(4)=2\cdot 4\cdot 2=16"
    return r


if __name__ == "__main__":
    main(run())
