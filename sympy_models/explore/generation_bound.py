"""The generation count from two sides: the chirality theorem and the dimension bound.

The chirality theorem gives N_gen = chi_V / chi_F = 3. The Distler-Garibaldi dimension count, as THEORY.md
states it, asks for W = 4 F N_gen linearly independent weight vectors, with F Weyl fields per generation,
against an eigenspace of dimension 128 or 112. This script tabulates the largest N_gen that fits, for
F = 15 (no right-handed neutrino) and F = 16 (with), and the eigenspace dimension three generations need.
It also records an arithmetic observation: 4 x 16 = 64 is the half-spinor of so(14), the Majorana-Weyl 64
of Lisi's model, and three of them are 192 = 3 x 64 = the roots of the three triality blocks of D4 x D4.
"""
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from explore import Result, main
from explore._roots import RootSystem


def run() -> Result:
    r = Result("generation_bound", "How many generations fit in an eigenspace of E8", prose="")
    e8 = RootSystem.of("E8")
    g1, g8 = e8.grading(1, 2), e8.grading(8, 2)
    eigen = {"node 1": g1[1], "node 8": g8[1]}                       # the -1 eigenspaces: 128 and 112
    r.check("the -1 eigenspaces of the two involutions have dimensions 128 and 112, read from the roots",
            (eigen["node 1"], eigen["node 8"]) == (128, 112))
    table = {}
    for F in (15, 16):
        for name, dim in eigen.items():
            table[(F, name)] = dim // (4 * F)
    r.data["largest N_gen with 4 F N_gen <= eigenspace"] = table
    r.check("with F = 15 at most 2 generations fit in 128 and 1 in 112; with F = 16, 2 and 1",
            (table[(15, "node 1")], table[(15, "node 8")], table[(16, "node 1")], table[(16, "node 8")]) == (2, 1, 2, 1))
    need = {F: 4 * F * 3 for F in (15, 16)}
    r.check("three generations need 180 weight vectors (F = 15) or 192 (F = 16)", need == {15: 180, 16: 192})
    r.check("both exceed the larger eigenspace (128)", all(v > 128 for v in need.values()))
    # the observation
    S14 = 2 ** (14 // 2 - 1)
    r.check("4 x 16 = 64 is the half-spinor of so(14) (the Majorana-Weyl 64 of Lisi's model)", 4 * 16 == S14 == 64)
    d4 = [tuple(1 if j == i else 0 for j in range(8)) for i in (1, 2, 3, 4)]
    in_d4 = [x for x in e8.roots if all(x[k] == 0 for k in range(8) if k not in (1, 2, 3, 4))]
    perp = [x for x in e8.roots if all(e8.ip(x, s) == 0 for s in d4)]
    r.check("3 x 64 = 192 is the number of roots outside D4 x D4, the three triality blocks",
            3 * S14 == 240 - len(in_d4) - len(perp) == 192)
    r.check("numerically, 4 F at F = 16 equals the size of one triality block, 8 x 8 = 64 (an equality of numbers, not an explanation)",
            4 * 16 == 8 * 8)
    r.prose = (r"The chirality theorem's three generations need $4\cdot 15\cdot 3=180$ weight vectors, or $4\cdot16\cdot3=192$ with "
               r"right-handed neutrinos, and the largest eigenspace of an inner involution of $E_8$ has $128$ (read from the roots; the other has $112$): "
               r"at most two generations fit. This is the dimension count only, with the factor $4$ as THEORY.md states it, and it assumes "
               r"linearly independent generations, which Lisi's model denies. The arithmetic is suggestive: $4\cdot16=64$ is the half-spinor of "
               r"$\mathfrak{so}(14)$, and $192=3\cdot64$ is the number of roots in the three triality blocks of $D_4\times D_4$. "
               r"It is an observation about numbers, not an explanation.")
    r.display = r"4\cdot 16 = 64 = S(14),\qquad 3\cdot 64 = 192 = 240-24-24"
    return r


if __name__ == "__main__":
    main(run())
