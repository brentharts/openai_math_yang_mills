#!/usr/bin/env python3
"""Bloch's T^{3/2} law and its first lattice correction (CONTENTS entry 271).

Model
-----
Nearest-neighbour quantum Heisenberg ferromagnet on Z^3, H = -sum_<xy> S_x . S_y.
Around the fully polarised state a single magnon has energy

    S * eps(k),   eps(k) = 2 * sum_a (1 - cos k_a).

The *ideal-magnon density* is

    n(x) = int_{[-pi,pi]^3} d^3k/(2 pi)^3  1 / (exp(x eps(k)) - 1),   x = beta*S.

The repo paper "The first lattice correction to Bloch's law" proves that the
magnetisation deficit agrees with this density up to o(beta^{-5/2}), and that
the first correction beyond Bloch's term is 3 zeta(5/2) (beta S)^{-5/2} / (128 pi^{3/2}).

What this script does
---------------------
* Writes n(x) as an exact lattice sum of Bessel functions,
      n(x) = sum_{m>=1} [ exp(-2 x m) I_0(2 x m) ]^3,
  and evaluates it to machine precision (direct sum + analytic tail).
* Uses sympy to derive the large-x expansion
      n(x) = sum_j A_j x^{-3/2-j}
  from the Bessel asymptotic series, and checks that A_0 is Bloch's coefficient
  and A_1 is the correction stated in the paper.
* Plots the exact density against the two truncations.

It computes the *ideal magnon density* only. It does not simulate the many-body
magnetisation; that this density controls the deficit is the content of the paper.
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from scipy.special import ive, zeta

import common as C

ORDER = 3  # number of correction terms derived symbolically


# ----------------------------------------------------------------------------
# Symbolic part
# ----------------------------------------------------------------------------
def symbolic_expansion():
    """Large-x expansion coefficients A_j of n(x), derived with sympy."""
    u = sp.symbols("u", positive=True)  # u = 1/z
    # e^{-z} I_0(z) ~ (2 pi z)^{-1/2} * sum_k ((2k-1)!!)^2 / (k! 8^k) z^{-k}
    ck = [sp.factorial2(2 * k - 1) ** 2 / (sp.factorial(k) * 8**k) for k in range(ORDER + 1)]
    poly = sum(c * u**k for k, c in enumerate(ck))
    cube = sp.expand(poly**3)
    a = [cube.coeff(u, j) for j in range(ORDER + 1)]  # exact up to u^ORDER
    # m-th term is (4 pi x m)^{-3/2} * sum_j a_j (2 x m)^{-j}; summing over m gives zeta(3/2+j)
    A = [
        a[j] * sp.Rational(1, 2**j) * sp.zeta(sp.Rational(3, 2) + j) / (8 * sp.pi ** sp.Rational(3, 2))
        for j in range(ORDER + 1)
    ]
    return a, A


def verify_symbolic(a, A):
    bloch = sp.zeta(sp.Rational(3, 2)) / (8 * sp.pi ** sp.Rational(3, 2))
    paper_corr = 3 * sp.zeta(sp.Rational(5, 2)) / (128 * sp.pi ** sp.Rational(3, 2))
    C.check(sp.simplify(A[0] - bloch) == 0, "A_0 = zeta(3/2)/(8 pi^{3/2}) (Bloch coefficient)")
    C.check(sp.simplify(A[1] - paper_corr) == 0,
            "A_1 = 3 zeta(5/2)/(128 pi^{3/2}) (correction stated in the repo paper)")


# ----------------------------------------------------------------------------
# Numerical part
# ----------------------------------------------------------------------------
def n_exact(x: float, m_direct: int = 200_000, a=None) -> float:
    """Ideal-magnon density via the Bessel sum, with an analytic tail."""
    m = np.arange(1, m_direct + 1, dtype=float)
    direct = np.sum(ive(0, 2.0 * x * m) ** 3)
    q = m_direct + 1
    tail = 0.0
    for j, aj in enumerate(a):
        tail += aj * 2.0**-j * zeta(1.5 + j, q) * x ** (-1.5 - j)
    return float(direct + tail / (4.0 * np.pi) ** 1.5)


def n_grid(x: float, N: int = 192) -> float:
    """Independent check: midpoint rule on the Brillouin zone (offset grid avoids k=0)."""
    k = 2.0 * np.pi * (np.arange(N) + 0.5) / N - np.pi
    e1 = 2.0 * (1.0 - np.cos(k))
    eps = e1[:, None, None] + e1[None, :, None] + e1[None, None, :]
    return float(np.mean(1.0 / np.expm1(x * eps)))


def main() -> None:
    C.ensure_dirs()
    C.setup_style()
    import matplotlib.pyplot as plt

    a_sym, A_sym = symbolic_expansion()
    verify_symbolic(a_sym, A_sym)
    a_num = [float(v) for v in a_sym]
    A_num = [float(v.evalf(30)) for v in A_sym]

    # --- numerical sanity checks ------------------------------------------------
    x_test = 5.0
    C.check(abs(n_exact(x_test, 200_000, a_num) / n_exact(x_test, 2_000_000, a_num) - 1) < 1e-12,
            "Bessel sum independent of the direct/tail split")
    xg = 0.4
    rel = abs(n_grid(xg) / n_exact(xg, 200_000, a_num) - 1)
    C.check(rel < 1e-2, f"Bessel sum agrees with Brillouin-zone quadrature at x={xg} (rel. diff {rel:.1e})")

    # --- the curve ----------------------------------------------------------------
    xs = np.geomspace(1.0, 400.0, 60)
    ex = np.array([n_exact(x, 200_000, a_num) for x in xs])
    lead = A_num[0] * xs**-1.5
    corr1 = lead + A_num[1] * xs**-2.5
    corr2 = corr1 + A_num[2] * xs**-3.5
    resid = (ex - corr1) * xs**3.5  # should tend to A_2

    big = xs > 20
    C.check(abs(resid[-1] / A_num[2] - 1) < 0.05,
            f"(n - Bloch - correction) x^(7/2) -> A_2 = {A_num[2]:.5f} (got {resid[-1]:.5f} at x={xs[-1]:.0f})")
    C.check(np.all(np.abs(ex - corr1)[big] < np.abs(ex - lead)[big]),
            "adding the paper's correction improves the Bloch approximation for x > 20")

    # --- figures ------------------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.2))
    s = xs**1.5  # scale out Bloch's x^{-3/2} so the corrections are visible
    ax[0].semilogx(xs, ex * s, "k-", lw=2, label="exact ideal-magnon density")
    ax[0].semilogx(xs, lead * s, "--", label=r"Bloch: $A_0$")
    ax[0].semilogx(xs, corr1 * s, ":", lw=2, label=r"+ correction $A_1/x$")
    ax[0].semilogx(xs, corr2 * s, "-.", lw=1.2, label=r"+ $A_2/x^2$")
    ax[0].set_xlabel(r"$x=\beta S$")
    ax[0].set_ylabel(r"$n(x)\,x^{3/2}$")
    ax[0].set_title("Magnon density, Bloch factor removed")
    ax[0].legend()

    ax[1].semilogx(xs, resid, "o-", ms=3, label=r"$(n-A_0x^{-3/2}-A_1x^{-5/2})\,x^{7/2}$")
    ax[1].axhline(A_num[2], color="C3", ls="--", label=rf"$A_2={A_num[2]:.4f}$ (sympy)")
    ax[1].set_xlabel(r"$x=\beta S$")
    ax[1].set_title("Next-order coefficient emerges")
    ax[1].legend()
    C.savefig(fig, "bloch_law")

    # --- data and LaTeX -------------------------------------------------------------
    sample = [2, 5, 10, 20, 50, 100, 200]
    rows = []
    table_data = []
    for x in sample:
        e = n_exact(float(x), 200_000, a_num)
        l0 = A_num[0] * x**-1.5
        l1 = l0 + A_num[1] * x**-2.5
        rows.append([f"{x}", f"{e:.6e}", f"{l0:.6e}", f"{l1:.6e}", f"{(e - l1) * x**3.5:.4f}"])
        table_data.append({"x": x, "exact": e, "bloch": l0, "bloch_plus_correction": l1})
    C.save_json("bloch_law", {"A": A_num, "table": table_data})

    powers = ["x^{-3/2}", "x^{-5/2}", "x^{-7/2}"]
    expansion = " + ".join(
        [sp.latex(A_sym[j]) + r"\," + powers[j] for j in range(3)] + [r"O(x^{-9/2})"]
    )

    tex = r"""
\section{Bloch's law and its first lattice correction}
\label{sec:bloch}

\paragraph{Source in the repository.}
CONTENTS entry 271; preprints \texttt{Blochs-Law-for-Finite-Range-Heisenberg-Ferromagnets-in-Three-Dimensions-October-5-2026},
\texttt{The-first-lattice-correction-to-Blochs-law-October-5-2026}; Lean notes \texttt{lean/docs/271.md}.
The second paper states that the magnetisation deficit of the 3D nearest-neighbour quantum Heisenberg ferromagnet agrees
with the \emph{full ideal-magnon density} up to $o(\beta^{-5/2})$, and gives the correction
$3\zeta(5/2)(\beta S)^{-5/2}/(128\pi^{3/2})$ beyond Bloch's term.

\paragraph{Model computed here.}
Single-magnon energies around the polarised state are $S\,\varepsilon(\mathbf k)$ with
$\varepsilon(\mathbf k)=2\sum_{a=1}^{3}(1-\cos k_a)$. With $x=\beta S$ the ideal-magnon density is
""" + C.equation(
        r"n(x)=\int_{[-\pi,\pi]^3}\frac{d^3k}{(2\pi)^3}\,\frac{1}{e^{x\,\varepsilon(\mathbf k)}-1}"
        r"=\sum_{m\ge 1}\Big[e^{-2xm}I_0(2xm)\Big]^{3}."
    ) + r"""
The second form follows from $1/(e^{a}-1)=\sum_{m\ge1}e^{-ma}$ and
$\int_{-\pi}^{\pi}\frac{dk}{2\pi}e^{2xm\cos k}=I_0(2xm)$, and it is what the code evaluates.

\paragraph{Symbolic expansion (sympy).}
Inserting the Bessel asymptotics $e^{-z}I_0(z)\sim(2\pi z)^{-1/2}\sum_k\frac{((2k-1)!!)^2}{k!\,8^k}z^{-k}$,
cubing, and summing over $m$ with $\sum_m m^{-s}=\zeta(s)$ gives
""" + C.equation("n(x)\\sim " + expansion + ".", "eq:bloch-expansion") + r"""
The script checks symbolically that the leading coefficient is Bloch's
$\zeta(3/2)/(8\pi^{3/2})$ and that the first correction coefficient is exactly the
$3\zeta(5/2)/(128\pi^{3/2})$ quoted in the repository paper. This is a consistency check of the
quoted coefficient against the ideal-magnon expansion, not a re-proof of the theorem.

""" + C.figure("bloch_law",
               r"Left: exact ideal-magnon density with Bloch's term and the first correction. "
               r"Right: the remainder scaled by $x^{7/2}$ approaches the next coefficient $A_2$ derived with sympy.",
               "\\linewidth") + C.table(
        ["$x=\\beta S$", "exact $n(x)$", "Bloch", "Bloch + corr.", "$(n-\\text{approx})\\,x^{7/2}$"],
        rows,
        caption=r"Ideal-magnon density and its truncations. The last column tends to $A_2$ as $x\to\infty$.",
    ) + r"""
\paragraph{Caveat.} Only the ideal-magnon density is computed. The magnetisation of the interacting quantum
system is controlled by it only through the theorem in the repository paper (with the order of limits stated there:
volume limit first, then the zero-field derivative, then low temperature).
"""
    C.write_tex("bloch", tex)
    print("wrote build/tex/bloch.tex")


if __name__ == "__main__":
    main()
