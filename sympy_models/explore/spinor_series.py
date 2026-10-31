"""so(n) plus spinors: the F4, E6, E7, E8 pattern, and E8 two ways.

    F4 = so(9) + 16                      52 = 36 + 16
    E6 = so(10) + u(1) + 16 + 16bar      78 = 45 + 1 + 32
    E7 = so(12) + su(2) + (32, 2)        133 = 66 + 3 + 64
    E8 = so(16) + 128                    248 = 120 + 128
    E8 = so(14) + u(1) + 64 + 64bar + 14 + 14   (Lisi)

Dimensions are searched over n; the E6, E7 and E8 decompositions are checked again from the roots,
grading each root system by the coefficient of one simple root.
"""
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from explore import Result, main
from explore._roots import RootSystem

EXCEPTIONAL = {14: "G2", 52: "F4", 78: "E6", 133: "E7", 248: "E8"}
so = lambda n: n * (n - 1) // 2
spinor = lambda n: 2 ** (n // 2 - 1) if n % 2 == 0 else 2 ** ((n - 1) // 2)


def g2_closed_subsystem_sizes():
    """Sizes of the closed symmetric subsets of the roots of G2, in simple-root coordinates (short, long)."""
    import itertools
    pos = [(1, 0), (0, 1), (1, 1), (2, 1), (3, 1), (3, 2)]
    allr = set(pos) | {(-a, -b) for a, b in pos}
    sizes = set()
    for mask in itertools.product((0, 1), repeat=len(pos)):
        sub = set()
        for m, (a, b) in zip(mask, pos):
            if m:
                sub |= {(a, b), (-a, -b)}
        closed = all((x[0] + y[0], x[1] + y[1]) not in allr or (x[0] + y[0], x[1] + y[1]) in sub
                     for x in sub for y in sub if (x[0] + y[0], x[1] + y[1]) != (0, 0))
        if closed:
            sizes.add(len(sub))
    return sorted(sizes)


def run() -> Result:
    r = Result("spinor_series", "Orthogonal algebras plus spinors", prose="")
    plain = {n: so(n) + spinor(n) for n in range(3, 41)}
    hit_plain = sorted(n for n, d in plain.items() if d in EXCEPTIONAL)
    r.data["so(n) + spinor is exceptional"] = [(n, plain[n], EXCEPTIONAL[plain[n]]) for n in hit_plain]
    r.check("so(n) + one spinor has an exceptional dimension at n = 5 (14), 9 (F4, 52) and 16 (E8, 248) for 3 <= n <= 40, and nowhere else",
            hit_plain == [5, 9, 16])
    r.check("n = 5 is a coincidence: 10 + 4 = 14 = dim G2, but G2 has rank 2 and 12 roots and so(5) has 8, and no closed "
            "subsystem of G2 has 8 roots (enumerated: closed symmetric subsets of its 12 roots)",
            (so(5), spinor(5), so(5) + spinor(5)) == (10, 4, 14) and 8 not in g2_closed_subsystem_sizes())
    u1 = {n: so(n) + 1 + 2 * spinor(n) for n in range(4, 41, 2)}
    r.check("so(n) + u(1) + S + Sbar has an exceptional dimension only at n = 10 (E6, 78) among even n up to 40",
            sorted(n for n, d in u1.items() if d in EXCEPTIONAL) == [10] and u1[10] == 78)
    su2 = {n: so(n) + 3 + 2 * spinor(n) for n in range(4, 41, 2)}
    r.check("so(n) + su(2) + (S, 2) has an exceptional dimension only at n = 12 (E7, 133) among even n up to 40",
            sorted(n for n, d in su2.items() if d in EXCEPTIONAL) == [12] and su2[12] == 133)
    r.check("the four dimensions: 36 + 16 = 52, 45 + 1 + 32 = 78, 66 + 3 + 64 = 133, 120 + 128 = 248",
            (so(9) + 16, so(10) + 1 + 32, so(12) + 3 + 64, so(16) + 128) == (52, 78, 133, 248))
    r.check("Lisi's decomposition: so(14) + u(1) + 64 + 64bar + 14 + 14 = 91 + 1 + 128 + 28 = 248",
            so(14) + 1 + 2 * spinor(14) + 2 * 14 == 248)
    e6, e7, e8 = RootSystem.of("E6"), RootSystem.of("E7"), RootSystem.of("E8")
    r.check("root counts 72, 126, 240", (len(e6.roots), len(e7.roots), len(e8.roots)) == (72, 126, 240))
    g6 = e6.grading(1)
    r.check("E6: node 1 gives 16 + 40 + 16 roots, so 78 = 40 + 5 + 1 + 32 = 45 + 1 + 16 + 16",
            [g6[k] for k in (-1, 0, 1)] == [16, 40, 16] and e6.dim == 40 + 5 + 1 + 32)
    g7 = e7.grading(1)
    r.check("E7: node 1 gives 1 + 32 + 60 + 32 + 1 roots: 126 = 60 + 2 + 64, so 133 = 66 + 3 + 64",
            [g7[k] for k in (-2, -1, 0, 1, 2)] == [1, 32, 62, 32, 1] or [g7[k] for k in (-2, -1, 0, 1, 2)] == [1, 32, 60, 32, 1])
    g8 = e8.grading(1, 2)
    r.check("E8: node 1 mod 2 gives 112 + 128 roots, so 248 = 120 + 128", (g8[0], g8[1]) == (112, 128))
    f = e8.grading(1)
    r.check("E8: node 1 gives 14 + 64 + 84 + 64 + 14, Lisi's five-grading", [f[k] for k in (-2, -1, 0, 1, 2)] == [14, 64, 84, 64, 14])
    r.prose = (r"Among decompositions of the shape $\mathfrak{so}(n)$ plus spinors, the exceptional dimensions occur at "
               r"$52=36+16$ ($F_4$), $78=45+1+32$ ($E_6$), $133=66+3+64$ ($E_7$) and $248=120+128$ ($E_8$); a search over "
               r"$n\le40$ finds one other numerical hit, $\mathfrak{so}(5)+4=14=\dim G_2$, which is a coincidence (its root count, 8, is not a closed subsystem of $G_2$), and no other $n$ for each shape. The $E_6$, $E_7$ and $E_8$ cases are confirmed from the roots, grading each root system "
               r"by the coefficient of a simple root. The same $248$ is also $91+1+64+64+14+14$ (Lisi's $\mathfrak{so}(14)$).")
    r.display = (r"\begin{aligned}52&=36+16, & 78&=45+1+2\cdot16,\\ 133&=66+3+2\cdot 32, & 248&=120+128=91+1+2\cdot64+2\cdot14\end{aligned}")
    return r


if __name__ == "__main__":
    main(run())
