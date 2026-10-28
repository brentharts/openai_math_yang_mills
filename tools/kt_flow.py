#!/usr/bin/env python3
"""Kosterlitz renormalisation-group flow for the planar XY model (CONTENTS entry 216).

What the repo states (CONTENTS entry 216, square-lattice nearest-neighbour cosine XY model):
    critical axis correlations  C_{beta_c}(r) ~ A r^{-1/4} (log r)^{1/8},
    BKT essential singularity   sqrt(beta_c - beta) * log xi(beta) -> B,   with A, B > 0.
(Some further conclusions there carry stated assumptions; see the paper.)

Model computed here
-------------------
The standard Kosterlitz RG equations in the variables x = pi*K - 2 (K = J/T, so x = 0 is the
universal-jump value K = 2/pi) and y (vortex fugacity, suitably normalised):

        dx/dl = -y^2,        dy/dl = -x y,                  l = ln(length scale).

sympy is used to
  * show x^2 - y^2 is conserved and solve the flow in closed form on the disordered side,
  * obtain ln xi = l*(mu) in closed form, and prove mu * l* -> pi as mu -> 0,
  * integrate the spin scaling dimension along the critical flow x = y = 1/(l + c) exactly,
    which produces the r^{-1/4} (ln r)^{1/8} form.

scipy then integrates the flow numerically and compares with the closed forms.

This is the textbook RG picture. It reproduces the *form* of the two statements above (exponent 1/4,
log-correction exponent 1/8, essential singularity), not the non-universal constants A and B of the
lattice model, and it is not a proof of what the repo paper proves.
"""

from __future__ import annotations

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

import common as C


# ----------------------------------------------------------------------------
# Symbolic part
# ----------------------------------------------------------------------------
def symbolic():
    l, mu, y0, c, th0 = sp.symbols("l mu y_0 c theta_0", positive=True)
    x = sp.Function("x")(l)
    y = sp.Function("y")(l)

    # conserved quantity
    inv = x**2 - y**2
    dinv = sp.simplify(sp.diff(inv, l).subs({sp.Derivative(x, l): -(y**2), sp.Derivative(y, l): -x * y}))
    C.check(dinv == 0, "x^2 - y^2 is conserved by the flow")

    # closed-form solution on the disordered side:  x^2 - y^2 = -mu^2
    xs = -mu * sp.tan(mu * l + th0)
    ys = mu / sp.cos(mu * l + th0)
    # (rewrite in terms of cos so sympy can use sin^2 + cos^2 = 1 reliably)
    res_x = sp.simplify((sp.diff(xs, l) + ys**2).rewrite(sp.cos))
    res_y = sp.simplify((sp.diff(ys, l) + xs * ys).rewrite(sp.cos))
    C.check(res_x == 0, "x(l) = -mu tan(mu l + theta_0) solves dx/dl = -y^2")
    C.check(res_y == 0, "y(l) = mu sec(mu l + theta_0) solves dy/dl = -x y")

    # initial data (x0 > 0, y(0)=y0): cos(theta_0) = mu/y0 with theta_0 < 0.
    # Correlation length: y(l*) = 1  =>  cos(mu l* + theta_0) = mu.
    lstar = (sp.acos(mu) + sp.acos(mu / y0)) / mu
    lim = sp.limit(mu * lstar, mu, 0, "+")
    C.check(sp.simplify(lim - sp.pi) == 0, "mu * l*(mu) -> pi as mu -> 0 (essential singularity)")

    # critical flow x = y = 1/(l + c); spin scaling dimension eta = 1/(2 pi K) = 1/(2 (2 + x))
    u = sp.symbols("u", positive=True)
    eta = 1 / (2 * (2 + 1 / u))
    prim = sp.integrate(sp.simplify(eta), u)
    C.check(sp.simplify(prim - (u / 4 - sp.log(2 * u + 1) / 8)) == 0,
            "int eta dl = u/4 - (1/8) ln(2u+1), u = l + c")
    return {"lstar": lstar, "lim": lim, "prim": prim, "xs": xs, "ys": ys, "mu": mu, "y0": y0}


# ----------------------------------------------------------------------------
# Numerical part
# ----------------------------------------------------------------------------
def flow_lstar(mu: float, y0: float) -> float:
    """Integrate dx/dl=-y^2, dy/dl=-xy from (x0,y0), x0>0 on x^2-y^2=-mu^2, until y=1."""
    x0 = np.sqrt(y0**2 - mu**2)

    def rhs(_l, v):
        return [-(v[1] ** 2), -v[0] * v[1]]

    def hit(_l, v):
        return v[1] - 1.0

    hit.terminal = True
    hit.direction = 1
    sol = solve_ivp(rhs, [0, 10.0 / mu], [x0, y0], events=hit, rtol=1e-11, atol=1e-13, method="DOP853")
    C.check(len(sol.t_events[0]) == 1, f"flow reached y=1 (mu={mu}, y0={y0})")
    return float(sol.t_events[0][0])


def lstar_closed(mu: float, y0: float) -> float:
    return (np.arccos(mu) + np.arccos(mu / y0)) / mu


def main() -> None:
    import matplotlib.pyplot as plt

    C.ensure_dirs()
    C.setup_style()
    sym = symbolic()

    # --- numerics vs closed form ---------------------------------------------------
    mus = np.array([0.4, 0.2, 0.1, 0.05, 0.02, 0.01])
    table = {}
    for y0 in (0.3, 0.1):
        num = np.array([flow_lstar(m, y0) for m in mus if m < y0])
        cl = np.array([lstar_closed(m, y0) for m in mus if m < y0])
        C.check(np.max(np.abs(num / cl - 1)) < 1e-6, f"numerical l* matches closed form for y0={y0}")
        table[y0] = (mus[mus < y0], num)
    m_small = table[0.1][0][-1]
    C.check(abs(table[0.1][0][-1] * table[0.1][1][-1] / np.pi - 1) < 0.1,
            f"mu*l* within 10% of pi at mu={m_small}")

    # --- critical correlation: r^{1/4} C(r) ~ (ln r)^{1/8} -------------------------------
    cc = 1.0
    ls = np.array([10, 1e2, 1e3, 1e4, 1e6, 1e9])
    u = ls + cc
    local_exp = (1 / 8) * (2 * ls) / (2 * u + 1)  # d ln(r^{1/4} C)/d ln l
    C.check(abs(local_exp[-1] - 1 / 8) < 1e-6, "log-correction exponent -> 1/8 as ln r -> infinity")

    # --- figure ---------------------------------------------------------------------------
    fig, ax = plt.subplots(1, 3, figsize=(14.5, 4.2))

    # (a) phase portrait
    xx = np.linspace(-0.35, 0.5, 20)
    yy = np.linspace(0.0, 0.5, 14)
    X, Y = np.meshgrid(xx, yy)
    ax[0].streamplot(X, Y, -(Y**2), -X * Y, color="0.7", density=1.1, linewidth=0.8, arrowsize=0.8)
    for sgn, st in ((1, "-"), (-1, "-")):
        pass
    ax[0].plot([0, 0.5], [0, 0.5], "k--", lw=1.2, label="separatrix $x=y$")
    ax[0].plot([-0.35, 0], [0.35, 0], "k--", lw=1.2)
    for mu_ in (0.03, 0.08):  # disordered-side trajectories x^2 - y^2 = -mu^2
        y0_ = 0.3
        x0_ = np.sqrt(y0_**2 - mu_**2)
        sol = solve_ivp(lambda _l, v: [-(v[1] ** 2), -v[0] * v[1]], [0, 4.0 / mu_], [x0_, y0_],
                        rtol=1e-9, atol=1e-12, dense_output=True)
        t = np.linspace(0, sol.t[-1], 600)
        v = sol.sol(t)
        keep = v[1] < 0.5
        ax[0].plot(v[0][keep], v[1][keep], lw=2, label=rf"$\mu={mu_}$")
    ax[0].set_xlim(-0.35, 0.5)
    ax[0].set_ylim(0, 0.5)
    ax[0].set_xlabel(r"$x=\pi K-2$")
    ax[0].set_ylabel(r"$y$ (vortex fugacity)")
    ax[0].set_title("Kosterlitz flow")
    ax[0].legend(fontsize=8, loc="upper left")

    # (b) mu * l*  -> pi
    mm = np.geomspace(0.003, 0.3, 80)
    for y0, st in ((0.3, "C0"), (0.1, "C1")):
        sel = mm < y0
        ax[1].semilogx(mm[sel], mm[sel] * lstar_closed(mm[sel], y0), "-", color=st,
                       label=rf"closed form, $y_0={y0}$")
        mu_n, l_n = table[y0]
        ax[1].semilogx(mu_n, mu_n * l_n, "o", color=st, ms=5)
    ax[1].axhline(np.pi, color="k", ls="--", lw=1, label=r"$\pi$")
    ax[1].set_xlabel(r"$\mu\propto\sqrt{T-T_c}$")
    ax[1].set_ylabel(r"$\mu\,\ln\xi$")
    ax[1].set_title(r"$\sqrt{t}\,\ln\xi\to$ const (dots: ODE solver)")
    ax[1].legend(fontsize=8)

    # (c) log correction
    lgrid = np.geomspace(5, 1e9, 200)
    ax[2].semilogx(lgrid, (1 / 8) * 2 * lgrid / (2 * (lgrid + cc) + 1), "C2-")
    ax[2].axhline(1 / 8, color="k", ls="--", lw=1)
    ax[2].set_xlabel(r"$\ln r$")
    ax[2].set_ylabel(r"$d\ln(r^{1/4}C)/d\ln\ln r$")
    ax[2].set_title(r"Critical $C(r)\sim r^{-1/4}(\ln r)^{1/8}$")
    ax[2].text(8, 0.1255, "$1/8$", fontsize=9)
    C.savefig(fig, "kt_flow")

    # --- data and LaTeX ------------------------------------------------------------------------
    rows = []
    for m, n in zip(*table[0.1]):
        rows.append([f"{m:.3f}", f"{n:.2f}", f"{lstar_closed(m, 0.1):.2f}", f"{m * n:.3f}"])
    def fmt_l(v: float) -> str:
        return f"{int(v)}" if v < 1e4 else f"$10^{{{int(round(np.log10(v)))}}}$"

    rows2 = [[fmt_l(l_), f"{e:.5f}"] for l_, e in zip(ls, local_exp)]
    C.save_json("kt_flow", {"mu_lstar_to_pi": [(float(m), float(m * n)) for m, n in zip(*table[0.1])],
                            "log_exponent": [(float(a), float(b)) for a, b in zip(ls, local_exp)]})

    tex = r"""
\section{XY model: Kosterlitz flow, BKT scaling and the \texorpdfstring{$(\log r)^{1/8}$}{(log r) exponent 1/8} correction}
\label{sec:kt}

\paragraph{Source in the repository.}
CONTENTS entry 216; preprints \texttt{BKT-universality-for-height-and-planar-spin-fields-September-24-2026},
\texttt{Critical-Center-Magnetization-in-the-Planar-XY-Model-October-5-2026},
\texttt{The-Critical-Spin-Field-of-the-Planar-XY-Model-October-5-2026}. For the square-lattice nearest-neighbour cosine XY model the
summary states critical axis correlations $C_{\beta_c}(r)\sim A\,r^{-1/4}(\log r)^{1/8}$ and the BKT essential singularity
$\sqrt{\beta_c-\beta}\,\log\xi(\beta)\to B$ with $A,B>0$. Some conclusions in that family carry stated assumptions.

\paragraph{Model computed here.} The textbook Kosterlitz RG equations in the variables $x=\pi K-2$ ($K=J/T$) and the vortex fugacity $y$:
""" + C.equation(r"\frac{dx}{d\ell}=-y^{2},\qquad \frac{dy}{d\ell}=-xy,\qquad \ell=\ln(\text{length}).", "eq:kt") + r"""
\paragraph{Exact solution (sympy).} The quantity $x^2-y^2$ is conserved. On the disordered side, $x^2-y^2=-\mu^2$, the solution is
""" + C.equation(
        "x(\\ell)=" + sp.latex(sym["xs"]) + ",\\qquad y(\\ell)=" + sp.latex(sym["ys"]) + ".", "eq:kt-sol"
    ) + r"""
Starting from $y(0)=y_0$ with $x(0)>0$ and stopping when $y=1$ (strong coupling) gives the correlation length
""" + C.equation(r"\ln\xi=\ell^*(\mu)=" + sp.latex(sym["lstar"]) + r",\qquad \mu\,\ell^*\;\longrightarrow\;" + sp.latex(sym["lim"]) + r"\ \ (\mu\to0),", "eq:lstar") + r"""
and since $\mu^2\propto T-T_c$ this is the essential singularity $\ln\xi\sim \pi/\mu\propto 1/\sqrt{T-T_c}$.
On the critical separatrix $x=y=1/(\ell+c)$. The spin scaling dimension is $\eta=1/(2\pi K)=1/(2(2+x))$, so with $u=\ell+c$
""" + C.equation(
        r"\int\eta\,d\ell=" + sp.latex(sym["prim"]) + r"\quad\Longrightarrow\quad C(r)\propto e^{-\int\eta\,d\ell}\propto r^{-1/4}\,(\ln r)^{1/8}.", "eq:corr"
    ) + r"""
The leading $r^{-1/4}$ comes from $\eta=1/4$ at $K=2/\pi$, and the exponent $1/8$ of the logarithm comes from the first correction $\eta\approx\frac14-\frac x8$.
Both match the form quoted from the repository paper; the RG does not predict the lattice constants $A$ and $B$.

""" + C.figure(
        "kt_flow",
        r"Left: Kosterlitz flow with the separatrices $x=\pm y$ and two disordered-side trajectories. Middle: $\mu\,\ln\xi$ from the closed form "
        r"(lines) and from numerical integration (dots) tends to $\pi$. Right: the local exponent of the logarithmic correction tends to $1/8$.",
        "\\linewidth",
    ) + C.table(
        [r"$\mu$", r"$\ell^*$ (ODE)", r"$\ell^*$ (closed form)", r"$\mu\,\ell^*$"], rows,
        caption=r"Numerical integration versus the closed form for $y_0=0.1$; $\mu\,\ell^*\to\pi$ as $\mu\to0$.",
    ) + C.table(
        [r"$\ln r$", r"$d\ln(r^{1/4}C)/d\ln\ln r$"], rows2,
        caption=r"Local exponent of the logarithmic correction on the critical flow ($c=1$); the limit is $1/8$.",
        colspec="rr",
    ) + r"""
\paragraph{Caveat.} This is the standard continuum RG picture, with a conserved quantity and exact formulas from sympy.
It illustrates why the exponents in the repository's statements take the values $1/4$ and $1/8$ and why $\ln\xi\sim1/\sqrt{t}$; it is not
the lattice proof, and the conventions ($x$, $y$ normalisation, $\mu\propto\sqrt{t}$) are those of the standard flow, not of the repository paper.
"""
    C.write_tex("kt_flow", tex)
    print("wrote build/tex/kt_flow.tex")


if __name__ == "__main__":
    main()
