"""Why is the hat expansion ratio the fourth power of the golden ratio?

Romanov's observation: the hat inflation factor kappa = (7 + 3 sqrt 5)/2 equals phi^4, and the same number is the ratio of
unreflected to reflected hats. This script takes the four-metatile substitution (H, T, P, F; quoted from CICY
pyCICY.monotile.SUBSTITUTION, columns = parent, rows = children) and derives what can be derived, exactly, in sympy:

  * the characteristic polynomial, the Perron root and its algebra (kappa^2 - 7 kappa + 1 = 0, kappa = 3 phi + 2, Lucas traces);
  * the Perron eigenvector (the metatile frequencies), irrational, hence no period;
  * the chirality ratio N_u / N_r = kappa, and a structural reason for it: the unreflected-hat count of each metatile equals
    the number of H metatiles in its inflation, i.e. the H row of the substitution matrix, so N_u = (M f)_H = kappa f_H while
    N_r = f_H. This is a relation between two quoted tables, checked here; it is not derived from the hat's geometry;
  * the linear factor: the area of the whole tiling scales by kappa, so lengths scale by sqrt(kappa) = phi^2;
  * the contrast with the Spectre, whose unit 4 + sqrt 15 is not a power of phi (trace 8 is not a Lucas number).
It does not say why the substitution has these entries; that is the geometry of the hat and is not attempted.
"""
import sys
from pathlib import Path

import sympy as sp

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import knowledgebase as kb
from explore import Result, main

M = sp.Matrix([[3, 1, 2, 2], [1, 0, 0, 0], [3, 0, 1, 1], [3, 0, 2, 3]])
NAMES = ("H", "T", "P", "F")
HATS = {"H": (4, 1), "T": (1, 0), "P": (2, 0), "F": (2, 0)}       # (total, reflected) per metatile
PHI = (1 + sp.sqrt(5)) / 2


def lucas(n):
    a, b = 2, 1
    for _ in range(n):
        a, b = b, a + b
    return a


def is_phi_power(x, nmax=40):
    """x is phi^n (n >= 1) iff x is a root of t^2 - L_n t + (-1)^n with x > 1."""
    for n in range(1, nmax):
        if sp.simplify(x ** 2 - lucas(n) * x + (-1) ** n) == 0 and x > 1:
            return n
    return None


def run() -> Result:
    r = Result("hat_golden", "The hat expansion ratio is the fourth power of the golden ratio", prose="")
    x = sp.Symbol("x")
    cp = sp.factor(M.charpoly(x).as_expr())
    r.check("charpoly of the hat substitution is (x-1)(x+1)(x^2-7x+1)", sp.expand(cp - (x - 1) * (x + 1) * (x ** 2 - 7 * x + 1)) == 0)
    kappa = (7 + 3 * sp.sqrt(5)) / 2
    r.check("kappa = (7+3 sqrt5)/2 is a root of x^2-7x+1 and is the Perron root of M",
            sp.simplify(kappa ** 2 - 7 * kappa + 1) == 0 and max(M.eigenvals(), key=lambda e: float(e)) == kappa)
    r.check("kappa = phi^4 = 3 phi + 2", sp.simplify(kappa - PHI ** 4) == 0 and sp.simplify(kappa - (3 * PHI + 2)) == 0)
    r.check("kappa + 1/kappa = 7 = L_4, the fourth Lucas number", sp.simplify(kappa + 1 / kappa) == lucas(4) == 7)
    r.check("kappa is phi^n for n = 4 by the Lucas-quadratic test", is_phi_power(kappa) == 4)
    # powers of the matrix: trace and entries are Lucas-type
    for n in range(1, 8):
        tr = (M ** n).trace()
        if tr != 1 + (-1) ** n + lucas(4 * n):
            raise AssertionError(f"tr M^{n}")
    r.checks.append("tr M^n = 1 + (-1)^n + L_{4n} for n = 1..7: the spectrum {1, -1, phi^4, phi^-4} exactly")
    # Fibonacci form of kappa^n
    fib = sp.fibonacci
    r.check("kappa^n = (L_{4n} + F_{4n} sqrt5)/2 for n = 1..6",
            all(sp.simplify(kappa ** n - (lucas(4 * n) + fib(4 * n) * sp.sqrt(5)) / 2) == 0 for n in range(1, 7)))
    # Perron eigenvector (frequencies)
    v = (M - kappa * sp.eye(4)).nullspace()
    r.check("the Perron eigenspace is one-dimensional", len(v) == 1)
    v = sp.simplify(v[0] / sum(v[0]))
    f = dict(zip(NAMES, [sp.nsimplify(sp.radsimp(sp.simplify(c))) for c in v]))
    irr = {k: sp.simplify(c).coeff(sp.sqrt(5)) != 0 for k, c in f.items()}
    r.check("the Perron frequencies are positive; T, P, F have non-zero sqrt5 part, so they are irrational: no period "
            "(H is rational, f_H = 1/3, so 'every frequency is irrational' is too strong, though one irrational suffices)",
            all(float(c) > 0 for c in f.values()) and f["H"] == sp.Rational(1, 3) and
            irr == {"H": False, "T": True, "P": True, "F": True})
    # chirality ratio
    fv = sp.Matrix([f[n] for n in NAMES])
    total = sum(HATS[n][0] * f[n] for n in NAMES)
    refl = sum(HATS[n][1] * f[n] for n in NAMES)
    ratio = sp.simplify((total - refl) / refl)
    r.check("unreflected/reflected hat ratio, weighted by the Perron frequencies, equals kappa = phi^4", sp.simplify(ratio - kappa) == 0)
    r.check("the reflected fraction is 1/(1+phi^4)", sp.simplify(refl / total - 1 / (1 + kappa)) == 0)
    # structural reason
    u = sp.Matrix([[HATS[n][0] - HATS[n][1] for n in NAMES]])
    r.check("structural: the unreflected-hat count per metatile is exactly the H row of the substitution matrix (3,1,2,2)", u == M[0, :])
    r.check("hence N_u = (M f)_H = kappa f_H, while N_r = f_H (the one reflected hat sits in each H): ratio kappa for any positive "
            "Perron vector", sp.simplify((u * fv)[0] - kappa * f["H"]) == 0 and sp.simplify(refl - f["H"]) == 0)
    # the same ratio, from the knowledge graph's library
    chi = kb.equation("Hat chirality ratio")
    cst = kb.equation("Hat constant as a golden power")
    r.check("the library holds both statements: N_u/N_r = kappa and kappa = phi^4", chi is not None and cst is not None)
    # linear factor
    lin = sp.sqrt(kappa)
    r.check("the area of the tiling scales by kappa per inflation, so lengths scale by sqrt(kappa) = phi^2 = (3+sqrt5)/2",
            sp.simplify(lin - PHI ** 2) == 0 and sp.simplify(PHI ** 2 - (3 + sp.sqrt(5)) / 2) == 0)
    # spectre contrast
    mu = 4 + sp.sqrt(15)
    r.check("the Spectre unit mu = 4+sqrt15 has trace 8, not a Lucas number, so it is not a power of phi",
            sp.simplify(mu ** 2 - 8 * mu + 1) == 0 and 8 not in [lucas(n) for n in range(0, 30)] and is_phi_power(mu) is None)
    r.check("the Spectre's perimeter growth 2+sqrt5 is phi^3 (trace 4 = L_3)", sp.simplify(2 + sp.sqrt(5) - PHI ** 3) == 0 and is_phi_power(2 + sp.sqrt(5)) == 3)
    r.check("phi^3 and phi^4 are adjacent Lucas quadratics x^2-4x-1 and x^2-7x+1, so their sharing Q(sqrt5) carries no further information",
            sp.simplify((2 + sp.sqrt(5)) ** 2 - 4 * (2 + sp.sqrt(5)) - 1) == 0)
    r.data.update({"charpoly": cp, "frequencies": {k: round(float(c), 5) for k, c in f.items()},
                   "kappa": float(kappa)})
    r.prose = (r"The factor is not a coincidence of decimals. The four-metatile substitution has characteristic polynomial "
               r"$(x-1)(x+1)(x^2-7x+1)$, so its Perron root solves $\kappa^2-7\kappa+1=0$; and $x^2-L_nx+(-1)^n$ is the "
               r"minimal polynomial of $\varphi^n$, with $L_4=7$, so $\kappa=\varphi^4=3\varphi+2$ and the whole spectrum is "
               r"$\{1,-1,\varphi^4,\varphi^{-4}\}$ (verified through $\operatorname{tr}M^n=1+(-1)^n+L_{4n}$). Because the area of "
               r"the tiling scales by $\kappa$ at each step, lengths scale by $\varphi^2$: the fourth power is the square of a "
               r"linear factor. The same number is the chirality ratio, and there is a structural reason in the data: the number of "
               r"unreflected hats in each metatile equals the number of $H$ metatiles in its inflation (the $H$ row $(3,1,2,2)$ of "
               r"$M$), while the only reflected hat sits in each $H$, so $N_u=(Mf)_H=\kappa f_H$ and $N_r=f_H$. This links two "
               r"quoted tables and is not derived from the hat's geometry. The Spectre is different: its unit $4+\sqrt{15}$ has "
               r"trace $8$, not a Lucas number, so it is no power of $\varphi$; only its perimeter growth $2+\sqrt5=\varphi^3$ "
               r"is golden. That $\varphi^3$ and $\varphi^4$ have adjacent Lucas quadratics is forced by sharing $\mathbb{Q}(\sqrt5)$ "
               r"and is not evidence of a deeper link.")
    r.display = (r"\underbrace{\frac{N_u}{N_r}}_{\text{chirality}}=\frac{(Mf)_H}{f_H}=\kappa=\frac{7+3\sqrt5}{2}=3\varphi+2=\varphi^{4}"
                 r"=\bigl(\varphi^{2}\bigr)^{2},\qquad \kappa+\kappa^{-1}=L_4=7")
    return r


if __name__ == "__main__":
    main(run())
