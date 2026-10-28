#!/usr/bin/env python3
"""Two-dimensional O(N) lattice models: mass scaling from Monte Carlo (CONTENTS entry 215).

Model
-----
Nearest-neighbour O(N) model on an L x L periodic square lattice with unit spins
q_x in S^{N-1} and Gibbs weight exp(beta * sum_<xy> q_x . q_y), no field, no topological term.
This is the model of the repo papers on the canonical O(3) continuum limit and on the
exact O(4) mass asymptotics.

What the repo papers state (from the CONTENTS summary of entry 215):
    O(4):  m_lat(beta) ~ 32 exp(pi/4 - 1/2) * sqrt(beta) * exp(-pi beta)      (full transfer gap)
    n >= 3: exponential decay of spin correlations at every positive temperature.

Asymptotic freedom predicts, for general N >= 3,
    m(beta) = C_N * beta^{1/(N-2)} * exp(-2 pi beta / (N-2)) * (1 + O(1/beta)).
(For N=4 this is the form of the O(4) formula above; for N=3 it is beta * exp(-2 pi beta).)

What this script does
---------------------
* Runs a Wolff single-cluster Monte Carlo for N = 3 and N = 4.
* Measures the susceptibility chi and the second-moment correlation length
      xi^2 = (chi/F - 1) / (4 sin^2(pi/L)),
  with F the structure factor at the smallest nonzero momentum.
* Forms R_N(beta) = (1/xi) / (beta^{1/(N-2)} exp(-2 pi beta/(N-2))), which should flatten
  at large beta if asymptotic scaling holds, and compares the O(4) plateau with the
  constant in the repo's O(4) formula.

1/xi_2nd is a second-moment mass, not the transfer-matrix gap; the two differ by an O(1)
factor, so this is a shape comparison and not a test of the constant.

Usage:  python tools/on_model_mc.py [--full]
"""

from __future__ import annotations

import sys
import time

import numpy as np

import common as C

try:
    from numba import njit

    HAVE_NUMBA = True
except ImportError:  # pragma: no cover - fallback path
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return lambda f: f


# ----------------------------------------------------------------------------
# Wolff cluster algorithm for O(N)
# ----------------------------------------------------------------------------
@njit(cache=True)
def _wolff_cluster(s, beta, L, N, mark, stamp, stack, r):
    """Grow and flip one Wolff cluster; return its size."""
    V = L * L
    # random reflection axis r (unit vector)
    nrm = 0.0
    for c in range(N):
        r[c] = np.random.normal()
        nrm += r[c] * r[c]
    nrm = np.sqrt(nrm)
    for c in range(N):
        r[c] /= nrm

    seed = np.random.randint(V)
    # reflect the seed
    p = 0.0
    for c in range(N):
        p += s[seed, c] * r[c]
    for c in range(N):
        s[seed, c] -= 2.0 * p * r[c]
    mark[seed] = stamp
    top = 0
    stack[top] = seed
    top += 1
    size = 1
    while top > 0:
        top -= 1
        i = stack[top]
        xi_ = i // L
        yi_ = i % L
        # projection of the (already reflected) spin i; original projection is minus this
        pi_cur = 0.0
        for c in range(N):
            pi_cur += s[i, c] * r[c]
        for d in range(4):
            if d == 0:
                j = ((xi_ + 1) % L) * L + yi_
            elif d == 1:
                j = ((xi_ - 1) % L) * L + yi_
            elif d == 2:
                j = xi_ * L + (yi_ + 1) % L
            else:
                j = xi_ * L + (yi_ - 1) % L
            if mark[j] == stamp:
                continue
            pj = 0.0
            for c in range(N):
                pj += s[j, c] * r[c]
            prod = -pi_cur * pj  # = p_i(original) * p_j
            if prod > 0.0:
                if np.random.random() < 1.0 - np.exp(-2.0 * beta * prod):
                    for c in range(N):
                        s[j, c] -= 2.0 * pj * r[c]
                    mark[j] = stamp
                    stack[top] = j
                    top += 1
                    size += 1
    return size


@njit(cache=True)
def _run(N, L, beta, ntherm, nmeas, seed):
    np.random.seed(seed)
    V = L * L
    s = np.empty((V, N))
    for i in range(V):
        nrm = 0.0
        for c in range(N):
            s[i, c] = np.random.normal()
            nrm += s[i, c] ** 2
        nrm = np.sqrt(nrm)
        for c in range(N):
            s[i, c] /= nrm
    mark = np.zeros(V, dtype=np.int64)
    stack = np.empty(V, dtype=np.int64)
    r = np.empty(N)
    stamp = 0
    chi = np.empty(nmeas)
    Fk = np.empty(nmeas)
    kk = 2.0 * np.pi / L
    cx = np.empty(L)
    sx = np.empty(L)
    for t in range(L):
        cx[t] = np.cos(kk * t)
        sx[t] = np.sin(kk * t)
    for sweep in range(ntherm + nmeas):
        flipped = 0
        while flipped < V:  # one "sweep" = clusters until ~V sites were flipped
            stamp += 1
            flipped += _wolff_cluster(s, beta, L, N, mark, stamp, stack, r)
        if sweep >= ntherm:
            m = np.zeros(N)
            ax_c = np.zeros(N)
            ax_s = np.zeros(N)
            ay_c = np.zeros(N)
            ay_s = np.zeros(N)
            for i in range(V):
                x = i // L
                y = i % L
                for c in range(N):
                    v = s[i, c]
                    m[c] += v
                    ax_c[c] += v * cx[x]
                    ax_s[c] += v * sx[x]
                    ay_c[c] += v * cx[y]
                    ay_s[c] += v * sx[y]
            msq = 0.0
            fx = 0.0
            fy = 0.0
            for c in range(N):
                msq += m[c] ** 2
                fx += ax_c[c] ** 2 + ax_s[c] ** 2
                fy += ay_c[c] ** 2 + ay_s[c] ** 2
            idx = sweep - ntherm
            chi[idx] = msq / V
            Fk[idx] = 0.5 * (fx + fy) / V
    return chi, Fk


def run_point(N: int, L: int, beta: float, ntherm: int, nmeas: int, seed: int = 12345, nblocks: int = 20):
    """Return (chi, xi, xi_err) for one (N, L, beta)."""
    chi_t, F_t = _run(N, L, beta, ntherm, nmeas, seed)
    # jackknife over blocks, which absorbs autocorrelation if blocks are long enough
    blk = nmeas // nblocks
    chi_b = chi_t[: blk * nblocks].reshape(nblocks, blk).mean(axis=1)
    F_b = F_t[: blk * nblocks].reshape(nblocks, blk).mean(axis=1)

    def xi_of(chi, F):
        return np.sqrt(max(chi / F - 1.0, 1e-12) / (4.0 * np.sin(np.pi / L) ** 2))

    xi = xi_of(chi_b.mean(), F_b.mean())
    jk = np.array([
        xi_of(np.delete(chi_b, b).mean(), np.delete(F_b, b).mean()) for b in range(nblocks)
    ])
    err = np.sqrt((nblocks - 1) / nblocks * np.sum((jk - jk.mean()) ** 2))
    return float(chi_b.mean()), float(xi), float(err)


def high_temperature_check() -> None:
    """chi = 1 + 4 beta/N + O(beta^2) for the O(N) model; verify the simulation reproduces it."""
    N, L, beta = 3, 16, 0.05
    chi, _, _ = run_point(N, L, beta, 200, 20000, seed=7)
    pred = 1 + 4 * beta / N + 12 * (beta / N) ** 2
    C.check(abs(chi / pred - 1) < 0.01,
            f"high-temperature series: chi={chi:.4f} vs 1+4b/N+12(b/N)^2={pred:.4f} (O(3), beta=0.05)")


# ----------------------------------------------------------------------------
# Symbolic part
# ----------------------------------------------------------------------------
def symbolic_part():
    """Asymptotic-freedom form, the O(4) constant, and the coupling drift beta(a)."""
    import sympy as sp

    b, Nn = sp.symbols("beta N", positive=True)
    Cc, mphys, a = sp.symbols("C m_phys a", positive=True)
    m_general = Cc * b ** (1 / (Nn - 2)) * sp.exp(-2 * sp.pi * b / (Nn - 2))
    m_O3 = sp.simplify(m_general.subs(Nn, 3))
    m_O4 = sp.simplify(m_general.subs(Nn, 4))
    C4 = 32 * sp.exp(sp.pi / 4 - sp.Rational(1, 2))  # constant in the repo's O(4) formula
    # Coupling drift for N=3: demand a fixed physical mass, m_lat(beta) = a * m_phys.
    drift = sp.Eq(sp.log(1 / a), 2 * sp.pi * b - sp.log(b) - sp.log(Cc / mphys))
    return {"m_general": m_general, "m_O3": m_O3, "m_O4": m_O4, "C4": C4,
            "C4_num": float(C4.evalf()), "drift": drift}


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------
def choose_L(xi_pilot: float, cap: int) -> int:
    return int(min(cap, max(32, 8 * np.ceil(xi_pilot))))


def main(full: bool = False) -> None:
    import matplotlib.pyplot as plt
    import sympy as sp

    C.ensure_dirs()
    C.setup_style()
    high_temperature_check()

    if full:
        plan = {3: [1.1, 1.2, 1.3, 1.4, 1.5, 1.6], 4: [1.4, 1.6, 1.8, 2.0, 2.2, 2.4]}
        nmeas, ntherm, cap = 60000, 3000, 320
    else:
        plan = {3: [1.2, 1.3, 1.4, 1.5], 4: [1.6, 1.8, 2.0, 2.2]}
        nmeas, ntherm, cap = 12000, 1500, 192

    results = {}
    t0 = time.time()
    for N, betas in plan.items():
        rows = []
        for beta in betas:
            _, xi_p, _ = run_point(N, 64, beta, 500, 3000, seed=99)
            L = choose_L(xi_p, cap)
            chi, xi, err = run_point(N, L, beta, ntherm, nmeas, seed=2024 + int(100 * beta))
            rows.append({"beta": beta, "L": L, "chi": chi, "xi": xi, "xi_err": err,
                         "L_over_xi": L / xi})
            print(f"N={N} beta={beta:.2f} L={L:3d} chi={chi:9.3f} xi={xi:7.3f}+-{err:.3f}"
                  f"  L/xi={L / xi:5.1f}  [{time.time() - t0:.0f}s]", flush=True)
        results[N] = rows

    # --- analysis ----------------------------------------------------------------
    sym = symbolic_part()
    fits = {}
    for N, rows in results.items():
        beta = np.array([r["beta"] for r in rows])
        xi = np.array([r["xi"] for r in rows])
        err = np.array([r["xi_err"] for r in rows])
        for r in rows:
            C.check(r["L_over_xi"] > 4.5,
                    f"finite-size control: L/xi={r['L_over_xi']:.1f} for N={N}, beta={r['beta']}")
        C.check(bool(np.all(np.diff(xi) > 0)), f"xi increases with beta for N={N}")
        y = np.log(xi) + np.log(beta) / (N - 2)  # = 2 pi beta/(N-2) + const under asymptotic scaling
        slope, _ = np.polyfit(beta, y, 1, w=xi / err)
        theory = 2 * np.pi / (N - 2)
        fits[N] = {"slope": float(slope), "theory": float(theory)}
        # Asymptotic scaling is approached slowly (O(1/beta) corrections), so only a loose check.
        C.check(abs(slope / theory - 1) < 0.35,
                f"N={N}: fitted slope of corrected ln(xi) = {slope:.2f} vs asymptotic-freedom {theory:.2f}")
        for r in rows:
            R = 1.0 / (r["xi"] * r["beta"] ** (1 / (N - 2)) * np.exp(-2 * np.pi * r["beta"] / (N - 2)))
            r["R"] = float(R)

    C.save_json("on_model_mc", {"results": {str(k): v for k, v in results.items()},
                                "fits": {str(k): v for k, v in fits.items()},
                                "full": full, "numba": HAVE_NUMBA})

    # --- figures -------------------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.2))
    colors = {3: "C0", 4: "C3"}
    labels = {3: r"$\propto\beta^{-1}e^{2\pi\beta}$ (matched at first point)",
              4: r"$\propto\beta^{-1/2}e^{\pi\beta}$ (matched at first point)"}
    for N, rows in results.items():
        beta = np.array([r["beta"] for r in rows])
        xi = np.array([r["xi"] for r in rows])
        err = np.array([r["xi_err"] for r in rows])
        ax[0].errorbar(beta, xi, yerr=err, fmt="o", color=colors[N], label=f"O({N}) Monte Carlo")
        y0 = np.log(xi[0]) + np.log(beta[0]) / (N - 2)
        bb = np.linspace(beta.min() - 0.05, beta.max() + 0.05, 100)
        pred = np.exp(y0 + 2 * np.pi * (bb - beta[0]) / (N - 2)) * bb ** (-1 / (N - 2))
        ax[0].plot(bb, pred, "-", color=colors[N], alpha=0.5, label=labels[N])
    ax[0].set_yscale("log")
    ax[0].set_xlabel(r"$\beta$")
    ax[0].set_ylabel(r"$\xi_{2\mathrm{nd}}$")
    ax[0].set_title("Correlation length vs coupling")
    ax[0].legend(fontsize=7)

    for N, rows in results.items():
        ax[1].plot([r["beta"] for r in rows], [r["R"] for r in rows], "o-", color=colors[N],
                   label=f"O({N})")
    ax[1].axhline(sym["C4_num"], color="C3", ls=":", label=r"repo O(4) gap constant $32e^{\pi/4-1/2}$")
    ax[1].set_yscale("log")
    ax[1].set_xlabel(r"$\beta$")
    ax[1].set_ylabel(r"$R_N=\xi^{-1}\beta^{-1/(N-2)}e^{2\pi\beta/(N-2)}$")
    ax[1].set_title("Asymptotic-scaling ratio")
    ax[1].legend(fontsize=8)
    C.savefig(fig, "on_model_mc")

    # --- LaTeX -----------------------------------------------------------------------
    rows_tex = []
    for N, rows in results.items():
        for r in rows:
            rows_tex.append([f"O({N})", f"{r['beta']:.2f}", f"{r['L']}", f"{r['chi']:.2f}",
                             f"{r['xi']:.2f} $\\pm$ {r['xi_err']:.2f}", f"{r['R']:.2f}"])
    fit_text = ", ".join(f"O({k}) {v['slope']:.2f} (theory {v['theory']:.2f})" for k, v in fits.items())
    stats = "full" if full else "quick"

    tex = r"""
\section{Two-dimensional \texorpdfstring{O($N$)}{O(N)} models: mass scaling}
\label{sec:on}

\paragraph{Source in the repository.}
CONTENTS entry 215; preprints
\texttt{The-canonical-massive-continuum-limit-of-the-two-dimensional-O3-model-October-4-2026},
\texttt{Exact-mass-asymptotics-for-the-two-dimensional-O4-lattice-model-October-5-2026},
\texttt{An-Isolated-Particle-Pole-for-the-Two-Dimensional-O3-Spin-Field-October-4-2026}; Lean notes \texttt{lean/docs/215.md}.
The O(4) paper determines the full transfer-gap asymptotics (as summarised in CONTENTS)
""" + C.equation(
        r"m_{\mathrm{lat}}(\beta)\sim " + sp.latex(sym["C4"]) + r"\,\sqrt{\beta}\,e^{-\pi\beta},", "eq:o4"
    ) + r"""
and the O(3) paper builds a massive continuum limit along couplings $\beta_N\approx H+\frac{\log L}{2\pi}N$
with lattice spacings $a_N=L^{-N}$.

\paragraph{Model computed here.}
Spins $q_x\in S^{N-1}$ on an $L\times L$ torus with weight $\exp(\beta\sum_{\langle xy\rangle}q_x\cdot q_y)$.
A Wolff single-cluster Monte Carlo (numba-compiled when available) measures the susceptibility $\chi$ and the
second-moment correlation length
""" + C.equation(
        r"\xi^2=\frac{\chi/F-1}{4\sin^2(\pi/L)},\qquad "
        r"F=\frac{1}{V}\Big\langle\Big|\sum_x q_x e^{ik\cdot x}\Big|^2\Big\rangle_{k=(2\pi/L,0)}."
    ) + r"""
\paragraph{Asymptotic freedom (symbolic).}
Perturbative renormalisation predicts, for $N\ge3$,
""" + C.equation("m(\\beta)=" + sp.latex(sym["m_general"]) + r"\,(1+O(1/\beta)),", "eq:af") + r"""
which for $N=3$ is $""" + sp.latex(sym["m_O3"]) + r"""$ and for $N=4$ is $""" + sp.latex(sym["m_O4"]) + r"""$,
the form of the repository's O(4) formula~\eqref{eq:o4} with $C=""" + f"{sym['C4_num']:.2f}" + r"""$.
Demanding a fixed physical mass $m_{\rm lat}=a\,m_{\rm phys}$ at lattice spacing $a$ in the $N=3$ form gives
""" + C.equation(sp.latex(sym["drift"]), "eq:drift") + r"""
so $\beta(a)=\frac{1}{2\pi}\ln\frac1a+O(\ln\ln\frac1a)$: a drift of $\frac{\ln L}{2\pi}$ per unit of $N$ when $a_N=L^{-N}$,
plus a logarithmic correction. This has the same shape as the trajectory in the repository paper; the script does not
reproduce that paper's construction.

""" + C.figure(
        "on_model_mc",
        r"Left: second-moment correlation length from Wolff Monte Carlo with the asymptotic-freedom form matched at the "
        r"first point. Right: the ratio $R_N$ flattens at large $\beta$ if asymptotic scaling holds; the dotted line is the "
        r"constant in the repository's O(4) transfer-gap formula.",
        "\\linewidth",
    ) + C.table(
        ["model", r"$\beta$", "$L$", r"$\chi$", r"$\xi_{2\mathrm{nd}}$", "$R_N$"],
        rows_tex,
        caption=("Monte Carlo results (" + stats + " statistics). Fitted slopes of "
                 r"$\ln\xi+\frac{1}{N-2}\ln\beta$ against $\beta$: " + fit_text + "."),
    ) + r"""
\paragraph{Caveats.} $1/\xi_{2\mathrm{nd}}$ is a second-moment mass, not the transfer-matrix gap, so $R_N$ is not expected
to equal the O(4) gap constant; the comparison tests the \emph{shape} $\beta^{1/(N-2)}e^{-2\pi\beta/(N-2)}$.
Corrections of order $1/\beta$ are large at the couplings reachable here: the fitted slopes fall short of the
asymptotic-freedom values, and $R_3$ is still rising at the largest coupling, so no plateau is claimed for O(3).
Statistical errors are jackknife errors over blocks. This is a numerical illustration, not a proof.
"""
    C.write_tex("on_model", tex)
    print("wrote build/tex/on_model.tex")


if __name__ == "__main__":
    main(full="--full" in sys.argv)
