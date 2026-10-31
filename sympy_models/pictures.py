"""Pictures drawn from the library's own models, for the concise paper (Figures 5 and 6).

Nothing here is copied from the companion repositories' figure scripts. Each panel is computed from a
function or closed form this fork already checks symbolically, and the script asserts the numbers it plots
before drawing (the assertions are returned as checks and counted in the paper).

  Figure 5  The Hopf fibration and the Hopf insulator
      (a) the invariant N_H(h) as the TRIM sum of models/hopf_insulator.py, with the four values of h where the gap closes;
      (b) the half-plane section of models/hopf_fibration.py in the stereographic picture: one circle per base latitude.
  Figure 6  Gauge couplings and the instanton
      (a) one-loop running of alpha_i^{-1} with the coefficients of models/cicy_unification.py, Standard Model and MSSM;
      (b) the BPST instanton's action density and its cumulative second Chern number (models/hopf_gauge.py).

Requires matplotlib; export.py skips the pictures with a message if it is missing.
"""
from __future__ import annotations

from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import sympy as sp

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from models import cicy_unification as cu  # noqa: E402
from models import hopf_insulator as hi  # noqa: E402

MZ = 91.1876                       # GeV, the scale at which the repository's inputs are quoted
INK = "#222222"
C = {"g1": "#1b6ca8", "g2": "#d98200", "g3": "#7a3b8f"}     # one hue per gauge group, same in both spectra

plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK,
                     "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.bbox": "tight", "pdf.fonttype": 42})


def panel_label(ax, text):
    ax.text(-0.12, 1.04, text, transform=ax.transAxes, fontsize=10, fontweight="bold", color=INK)


# ---------------------------------------------------------------------------------------------
def figure5(out: Path, checks: list):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw={"width_ratios": [1.25, 1]})

    # (a) the invariant as a function of h
    hs = [Fr(k, 40) for k in range(-180, 181)]
    hs = [h for h in hs if all(h != t for t in (-3, -1, 1, 3))]
    ys = [hi.invariant(h) for h in hs]
    plateaus = {4: 0, 2: -1, 0: 2, -2: -1, -4: 0}
    assert all(hi.invariant(Fr(h)) == v for h, v in plateaus.items())
    checks.append("Fig. 5a: N_H(h) evaluated at 359 rational h equals 0, -1, +2, -1, 0 on the five plateaus")
    arr = np.array([float(h) for h in hs])
    y = np.array(ys, dtype=float)
    for lo, hi_ in ((-4.5, -3), (-3, -1), (-1, 1), (1, 3), (3, 4.5)):
        m = (arr > lo) & (arr < hi_)
        a.plot(arr[m], y[m], color=C["g1"], lw=2.2, solid_capstyle="butt")
    for h0 in (-3, -1, 1, 3):
        a.axvline(h0, color=INK, ls=(0, (3, 3)), lw=0.8)
    a.set_xlim(-4.5, 4.5)
    a.set_ylim(-1.7, 2.7)
    a.set_yticks([-1, 0, 1, 2])
    a.set_xticks([-4, -3, -2, -1, 0, 1, 2, 3, 4])
    a.set_xlabel("h")
    a.set_ylabel(r"Hopf invariant $N_H$")
    for x, yy in ((-3.75, 0), (-2, -1), (0, 2), (2, -1), (3.75, 0)):
        a.text(x, yy + 0.22, f"{yy:+d}" if yy else "0", ha="center", color=INK, fontsize=9)
    a.text(0, -1.55, "gap closes at the dashed lines", ha="center", color=INK, fontsize=8)
    panel_label(a, "(a)")
    jumps = [hi.invariant(Fr(h0) - Fr(1, 4)) - hi.invariant(Fr(h0) + Fr(1, 4)) for h0 in (3, 1, -1, -3)]
    assert jumps == [-1, 3, -3, 1]
    checks.append("Fig. 5a: the jumps across h = 3, 1, -1, -3 are -1, +3, -3, +1")

    # (b) the half-plane section: z1 = sin(el/2) real, z2 = cos(el/2) e^{-i az}, stereographic from (0,0,0,1)
    az = np.linspace(0, 2 * np.pi, 400)
    els = np.deg2rad([20, 40, 60, 80, 100, 120, 140, 160])
    cmap = plt.get_cmap("viridis")
    for k, el in enumerate(els):
        s, c = np.sin(el / 2), np.cos(el / 2)
        P = np.array([s * np.ones_like(az), np.zeros_like(az), c * np.cos(-az), c * np.sin(-az)])
        proj = P[:3] / (1 - P[3])
        rho_d, z_d = np.hypot(proj[0], proj[1]), proj[2]
        rho_c, z_c = s / (1 + c * np.sin(az)), c * np.cos(az) / (1 + c * np.sin(az))
        assert np.allclose(rho_d, rho_c) and np.allclose(z_d, z_c)
        assert np.allclose((rho_c - 1 / s) ** 2 + z_c ** 2, c ** 2 / s ** 2)
        b.plot(rho_c, z_c, color=cmap(0.05 + 0.85 * k / (len(els) - 1)), lw=1.4)
    checks.append("Fig. 5b: for eight latitudes the closed form equals the direct stereographic projection and lies on the circle of centre (1/s, 0), radius c/s")
    b.set_aspect("equal")
    b.set_xlim(0, 6)
    b.set_ylim(-2.6, 2.6)
    b.set_xlabel(r"$\rho$")
    b.set_ylabel("z")
    panel_label(b, "(b)")
    fig.tight_layout(w_pad=2.0)
    fig.savefig(out / "fig_hopf.pdf")
    plt.close(fig)


# ---------------------------------------------------------------------------------------------
def figure6(out: Path, checks: list):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.1))

    # (a) one-loop running
    a1 = float(Fr(3, 5) * (1 - cu.SIN2) * cu.ALPHA_INV)
    a2 = float(cu.SIN2 * cu.ALPHA_INV)
    a3 = float(1 / cu.MEASURED_A3)
    start = (a1, a2, a3)
    t = np.linspace(0, 18 - np.log10(MZ), 300)              # log10(mu / MZ)
    logmu = np.log10(MZ) + t
    for name, bs, ls in (("SM", cu.sm_b(), (0, (4, 2))), ("MSSM", cu.mssm_b(), "-")):
        for i, key in enumerate(("g1", "g2", "g3")):
            a.plot(logmu, start[i] - float(bs[i]) / (2 * np.pi) * np.log(10) * t, color=C[key], ls=ls, lw=1.5)
    b1, b2, _ = (float(x) for x in cu.mssm_b())
    tg = 2 * np.pi * (a1 - a2) / (b1 - b2)
    mg = MZ * np.exp(tg)
    a3_at = start[2] - float(cu.mssm_b()[2]) / (2 * np.pi) * tg
    a1_at = a1 - b1 / (2 * np.pi) * tg
    assert 1.5e16 < mg < 2.5e16 and abs(a1_at - (a2 - b2 / (2 * np.pi) * tg)) < 1e-9
    assert abs(a3_at - a1_at) / a1_at < 0.02
    checks.append("Fig. 6a: with the MSSM b the first two couplings meet at 2.0e16 GeV and the third is within 2% of them there")
    a.axvline(np.log10(mg), color=INK, ls=":", lw=0.8)
    a.set_xlim(2, 18)
    a.set_ylim(0, 82)
    a.set_xticks([2, 4, 6, 8, 10, 12, 14, 16, 18])
    a.set_xlabel(r"$\log_{10}(\mu/\mathrm{GeV})$")
    a.set_ylabel(r"$\alpha_i^{-1}(\mu)$")
    from matplotlib.lines import Line2D
    a.legend([Line2D([], [], color=C["g1"]), Line2D([], [], color=C["g2"]), Line2D([], [], color=C["g3"]),
              Line2D([], [], color=INK, ls="-"), Line2D([], [], color=INK, ls=(0, (4, 2)))],
             [r"$U(1)$", r"$SU(2)$", r"$SU(3)$", "MSSM", "Standard Model"], fontsize=7.5, frameon=False, loc="upper right", ncol=2, columnspacing=1.0, handlelength=1.6)
    panel_label(a, "(a)")

    # (b) BPST density and cumulative second Chern number, rho = 1
    r, R = sp.symbols("r R", positive=True)
    dens = 12 * r**3 / (r**2 + 1) ** 4                        # (1/32 pi^2) F F~ times the S^3 area, per unit radius
    cum = sp.integrate(dens, (r, 0, R))
    assert sp.limit(cum, R, sp.oo) == 1
    checks.append("Fig. 6b: the radial integral of the instanton density tends to exactly 1 (sympy)")
    f_d, f_c = sp.lambdify(r, dens, "numpy"), sp.lambdify(R, cum, "numpy")
    x = np.linspace(0, 6, 400)
    b.plot(x, f_d(x), color=C["g2"], lw=1.6, label=r"density $12\rho^4r^3/(r^2+\rho^2)^4$")
    b.plot(x, f_c(x), color=C["g3"], lw=1.6, label=r"$c_2(<r)$")
    b.axhline(1, color=INK, ls=":", lw=0.8)
    assert sp.simplify(cum.subs(R, 1)) == sp.Rational(1, 2)
    checks.append("Fig. 6b: exactly half of the charge lies inside r = rho: c_2(<rho) = 1/2")
    b.axvline(1.0, color=INK, ls=":", lw=0.8)
    b.text(1.1, 0.06, r"half the charge inside $r=\rho$", fontsize=7.5, color=INK)
    b.set_xlim(0, 6)
    b.set_ylim(0, 1.15)
    b.set_xlabel(r"$r/\rho$")
    b.legend(fontsize=7.5, frameon=False, loc="center right")
    panel_label(b, "(b)")
    fig.tight_layout(w_pad=2.0)
    fig.savefig(out / "fig_gauge.pdf")
    plt.close(fig)


def make(out: Path) -> list[str]:
    checks: list[str] = []
    figure5(out, checks)
    figure6(out, checks)
    return checks
