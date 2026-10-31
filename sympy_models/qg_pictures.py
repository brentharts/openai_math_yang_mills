"""Pictures for the quantum gravity paper (Figures 4 and 5 of qgravity_document.tex), drawn from models/qg_spin_network.py.

Nothing is copied from another repository. Every panel is computed from a function or closed form that the model already checks, and
the script asserts the numbers it plots before drawing them (the assertions are returned as checks and counted in the paper).

  Figure 4  Quantum geometry of a spin-network link and of its count
      (a) the single-link area eigenvalues sqrt(j(j+1)) against j, with the line j + 1/2 they approach, and the step between
          neighbouring levels, which falls to 1/2;
      (b) the entropy per puncture (1/N) ln C_{N/2} of N spin-1/2 punctures, exact Catalan values and the asymptotic formula, against N,
          with the limit ln 2 that fixes the Immirzi parameter.
  Figure 5  The cut between the two plaquettes of the ladder
      (a) the mean-square area <c2(J)> of the surface cutting the shared link against g^2, from the ground state of the ladder, with the
          exact strong-coupling coefficient 2/(3 g^8) and the weak-coupling law 3/(sqrt10 g^2);
      (b) g^2 <c2(J)> and g^2 <c2(j_a)> against g^2, approaching their exact limits as g -> 0.

Requires matplotlib; qgravity.py skips the pictures with a message if it is missing.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import sympy as sp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from models import qg_spin_network as qg  # noqa: E402

INK = "#222222"
C = {"blue": "#1b6ca8", "orange": "#d98200", "purple": "#7a3b8f", "grey": "#888888"}

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK,
                     "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.bbox": "tight", "pdf.fonttype": 42})


def panel_label(ax, text):
    ax.text(-0.12, 1.04, text, transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK)


# ---------------------------------------------------------------------------------------------
def figure4(out: Path, checks: list):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.1))

    # (a) area levels and their steps
    js = np.arange(1, 13) / 2                                          # j = 1/2, 1, ..., 6
    lev = np.sqrt(js * (js + 1))
    assert np.all(lev < js + 0.5) and np.all((js + 0.5) - lev > 0)
    assert abs(float(sp.sqrt(sp.Rational(3, 4))) - lev[0]) < 1e-15 and abs(lev[1] - np.sqrt(2)) < 1e-15
    checks.append("Fig. 4a: sqrt(j(j+1)) < j + 1/2 for the twelve plotted spins j = 1/2 .. 6, with the first two levels sqrt3/2 and sqrt2")
    step = np.sqrt((js + 0.5) * (js + 1.5)) - lev
    assert np.all(step > 0.5) and np.all(np.diff(step) < 0) and step[-1] - 0.5 < 0.01
    checks.append("Fig. 4a: the steps between neighbouring levels are above 1/2, strictly decreasing, and within 0.01 of 1/2 by j = 6")
    jj = np.linspace(0, 6.5, 100)
    a.plot(jj, jj + 0.5, color=C["grey"], ls=(0, (3, 3)), lw=1.0, label=r"$j+\frac{1}{2}$")
    a.plot(js, lev, "o", color=C["blue"], ms=4.5, label=r"$\sqrt{j(j+1)}$")
    a.set_xlim(0, 6.5)
    a.set_ylim(0, 7)
    a.set_xlabel(r"spin $j$")
    a.set_ylabel(r"area in units of $8\pi\gamma_I\ell_P^2$")
    a2 = a.twinx()
    a2.spines["right"].set_visible(True)
    a2.plot(js[:-1] + 0.25, step[:-1], "-", color=C["orange"], lw=1.4, label="step (right axis)")
    a2.axhline(0.5, color=C["orange"], ls=":", lw=0.8)
    a2.set_ylim(0.45, 0.6)
    a2.set_ylabel(r"step to the next level", color=C["orange"])
    a2.tick_params(axis="y", colors=C["orange"])
    h1, l1 = a.get_legend_handles_labels()
    h2, l2 = a2.get_legend_handles_labels()
    a.legend(h1 + h2, l1 + l2, frameon=False, loc="upper left", fontsize=7.5)
    panel_label(a, "(a)")

    # (b) entropy per puncture
    ns = np.array([2, 4, 6, 8, 10, 14, 20, 30, 50, 80, 130, 200, 320, 500, 800, 1300])
    ex = np.array([qg.entropy_per_puncture(int(n)) for n in ns])
    asym = np.log(2) - 1.5 * np.log(ns / 2) / ns - 0.5 * np.log(np.pi) / ns - 9 / (4 * ns ** 2)
    assert np.all(np.diff(ex) > 0) and np.all(ex < np.log(2))
    checks.append("Fig. 4b: the exact entropy per puncture (1/N) ln C_{N/2} at sixteen values of N from 2 to 1300 increases and stays below ln 2")
    big = ns >= 50
    assert np.max(np.abs(ex[big] - asym[big])) < 2e-4
    checks.append("Fig. 4b: for N >= 50 the asymptotic formula ln 2 - (3/2N) ln(N/2) - ln(pi)/(2N) - 9/(4N^2) is within 2e-4 of the exact value")
    b.axhline(np.log(2), color=C["grey"], ls=(0, (3, 3)), lw=1.0)
    b.plot(ns[big], asym[big], "-", color=C["orange"], lw=1.4, label="asymptotic formula")
    b.plot(ns, ex, "o", color=C["blue"], ms=3.8, label=r"exact $\frac{1}{N}\ln C_{N/2}$")
    b.set_xscale("log")
    b.set_ylim(0.2, 0.77)
    b.set_xlabel(r"number of spin-$\frac{1}{2}$ punctures $N$")
    b.set_ylabel("entropy per puncture")
    b.text(2.2, np.log(2) + 0.02, r"$\ln 2\ \Rightarrow\ \gamma_I=\ln 2/(\pi\sqrt{3})$", fontsize=8, color=INK)
    b.legend(frameon=False, loc="lower right", fontsize=8)
    panel_label(b, "(b)")
    fig.tight_layout(w_pad=2.4)
    fig.savefig(out / "fig_qg_geometry.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
def figure5(out: Path, checks: list):
    b_strong = float(qg.strong_coefficient())
    cJ_w, cA_w = (float(x) for x in qg.weak_limits())
    grid = np.array([0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0])
    rows = [qg.geometry(float(g)) for g in grid]
    cJ = np.array([r[1] for r in rows])
    cA = np.array([r[2] for r in rows])

    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.1))

    # (a) mean-square area of the cut against g^2, with both asymptotes
    assert np.all(np.diff(cJ) < 0)
    checks.append("Fig. 5a: <c2(J)> decreases strictly along the thirteen plotted values of g^2 from 0.05 to 8")
    weak = cJ_w / grid
    strong = b_strong / grid ** 4
    assert abs(cJ[0] / weak[0] - 1) < 0.03 and abs(cJ[-1] / strong[-1] - 1) < 1e-3 and abs(cJ[-2] / strong[-2] - 1) < 5e-3
    checks.append("Fig. 5a: <c2(J)> is within 3% of 3/(sqrt10 g^2) at g^2 = 0.05 and within 0.1% and 0.5% of 2/(3 g^8) at g^2 = 8 and 6")
    gg = np.geomspace(0.04, 12, 200)
    a.loglog(gg, cJ_w / gg, color=C["orange"], ls=(0, (4, 2)), lw=1.2, label=r"weak: $3/(\sqrt{10}\,g^2)$")
    a.loglog(gg[gg > 0.9], b_strong / gg[gg > 0.9] ** 4, color=C["purple"], ls=":", lw=1.4, label=r"strong: $2/(3g^8)$")
    a.loglog(grid, cJ, "o-", color=C["blue"], ms=4, lw=1.4, label=r"ground state of the ladder")
    a.set_xlim(0.04, 12)
    a.set_ylim(1e-5, 1e3)
    a.set_xlabel(r"$g^2$")
    a.set_ylabel(r"$\langle c_2(J)\rangle_0$, mean-square area of the cut")
    a.legend(frameon=False, loc="lower left", fontsize=7.5)
    panel_label(a, "(a)")

    # (b) g^2 <c2> against g^2: the approach to the harmonic limits
    m = grid <= 1.0
    gJ, gA = (grid * cJ)[m], (grid * cA)[m]
    assert np.all(gJ < cJ_w) and np.all(gA < cA_w) and np.all(np.diff(gJ) < 0) and np.all(np.diff(gA) < 0)
    checks.append("Fig. 5b: g^2 <c2(J)> and g^2 <c2(j_a)> lie below their weak-coupling limits and fall away from them as g^2 grows from 0.05 to 1")
    lin = np.polyfit(grid[:3], (grid * cJ)[:3], 1)
    assert abs(lin[1] - cJ_w) < 1e-3
    checks.append("Fig. 5b: a straight-line fit to g^2 <c2(J)> at g^2 = 0.05, 0.1, 0.2 meets the vertical axis within 1e-3 of 3/sqrt10")
    b.plot(grid[m], gJ, "o-", color=C["blue"], ms=4, lw=1.4, label=r"$g^2\langle c_2(J)\rangle$ (shared link)")
    b.plot(grid[m], gA, "s-", color=C["purple"], ms=3.8, lw=1.4, label=r"$g^2\langle c_2(j_a)\rangle$ (plaquette link)")
    b.plot([0], [cJ_w], marker=">", color=C["orange"], ms=6, clip_on=False, zorder=5)
    b.plot([0], [cA_w], marker=">", color=C["orange"], ms=6, clip_on=False, zorder=5)
    b.text(0.06, cJ_w - 0.012, r"$3/\sqrt{10}$", fontsize=8, color=C["orange"])
    b.text(0.06, cA_w + 0.035, r"$\frac{3}{4}\left(\frac{1}{\sqrt{6}}+\frac{1}{\sqrt{10}}\right)$", fontsize=8, color=C["orange"])
    b.set_xlim(-0.03, 1.05)
    b.set_ylim(0, 1.08)
    b.set_xlabel(r"$g^2$")
    b.set_ylabel(r"$g^2\langle c_2\rangle_0$")
    b.legend(frameon=False, loc="lower left", fontsize=7.5)
    panel_label(b, "(b)")
    fig.tight_layout(w_pad=2.4)
    fig.savefig(out / "fig_qg_ladder.pdf")
    plt.close(fig)


def make(out: Path) -> list[str]:
    checks: list[str] = []
    figure4(out, checks)
    figure5(out, checks)
    return checks
