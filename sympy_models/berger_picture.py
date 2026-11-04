"""Picture for SU2_BERGER_TOP.md, drawn from models/su2_berger_top.py.

  (a) the Laplacian levels of the Berger sphere against s = R1/R2 (R2 = 1): the base levels (m' = 0, flat, the S^2 spectrum) and the
      fibre-momentum levels (m' != 0, rising like 1/s^2 as the fibre collapses), the first nonzero eigenvalue min(8, 2 + 1/s^2) and
      Lichnerowicz's bound (3/2) min Ric;
  (b) the counting function N(L) at s = 0.1: it equals the S^2 count up to L = 2 + 1/s^2, where the first fibre-momentum level enters.

Every plotted number is asserted before it is drawn. Nothing is copied from another repository. Requires matplotlib.

    cd sympy_models && python3 berger_picture.py ../figures/su2_berger_levels.png
"""
from __future__ import annotations

import sys
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from models import su2_berger_top as bt  # noqa: E402

INK = "#222222"
C = {"blue": "#1b6ca8", "orange": "#d98200", "purple": "#7a3b8f", "grey": "#888888", "green": "#2a8a4a"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK, "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "savefig.bbox": "tight"})


def draw(path):
    s = np.geomspace(0.04, 3.0, 400)
    fig, ax = plt.subplots(1, 2, figsize=(10.6, 4.2), constrained_layout=True)

    a = ax[0]
    # levels with j <= 2 and m' >= 0 (m' and -m' coincide)
    drawn = []
    for k in range(0, 5):
        j = Fr(k, 2)
        for kk in range(k % 2, k + 1, 2):
            mp = Fr(kk, 2)
            lam = np.array([float(bt.eigenvalue(j, mp, Fr(float(x) ** 2).limit_denominator(10**9))) for x in s])
            assert np.allclose(lam, 4 * (float(j * (j + 1) - mp * mp)) + 4 * float(mp * mp) / s**2, rtol=1e-6)
            if mp == 0:
                a.plot(s, lam, color=C["blue"], lw=1.8, zorder=3)
                a.text(3.05, lam[-1], f"j={int(j)}", fontsize=7.5, color=C["blue"], va="center")
            elif j <= 1:
                a.plot(s, lam, color=C["orange"], lw=1.1, alpha=0.9)
                drawn.append((j, mp))
    lam1 = np.minimum(8.0, 2 + 1 / s**2)
    a.plot(s, lam1, color=INK, lw=1.6, ls="--", zorder=4, label=r"$\lambda_1=\min(8,\,2+1/s^2)$")
    ric = np.minimum(2 * s**2, 4 - 2 * s**2)
    lich = 1.5 * ric
    ok = (ric > 0)
    assert np.all(lam1[ok] >= lich[ok] - 1e-12)                 # Lichnerowicz holds wherever Ric_min > 0
    a.plot(s[ok], lich[ok], color=C["purple"], lw=1.2, ls=":", label=r"Lichnerowicz $1.5\,\min\mathrm{Ric}$")
    a.plot([], [], color=C["blue"], lw=1.8, label=r"$m'=0$: base $S^2$ levels $4l(l+1)$")
    a.plot([], [], color=C["orange"], lw=1.1, label=r"$m'\neq0$, $j\leq1$: fibre-momentum levels")
    a.axvline(1.0, color=C["grey"], lw=0.8)
    a.text(1.0, 1.005, "round $S^3$", transform=a.get_xaxis_transform(), fontsize=8, color=C["grey"], ha="center", va="bottom")
    a.axvline(1 / np.sqrt(6), color=C["green"], lw=0.8, ls="-.")
    a.text(1 / np.sqrt(6), 1.005, r"$s^2=1/6$", transform=a.get_xaxis_transform(), fontsize=8, color=C["green"], ha="center", va="bottom")
    a.set(xscale="log", ylim=(0, 36), xlim=(0.04, 3.0), xlabel=r"$s=R_1/R_2$ (fibre / base)", ylabel=r"eigenvalue $\times R_2^2$",
          title="Squashing the fibre lifts the fibre-momentum levels", )
    a.set_title(a.get_title(), pad=16)
    a.legend(frameon=False, fontsize=8, loc="upper right")

    b = ax[1]
    s2 = Fr(1, 100)
    L_thr = 2 + 1 / s2
    Ls = np.linspace(0, 160, 1601)
    n_s = np.array([bt.count_below(s2, Fr(L).limit_denominator(1000)) for L in Ls])
    n_b = np.array([bt.count_sphere2(float(L)) for L in Ls])
    assert np.all(n_s[Ls < float(L_thr) - 1e-9] == n_b[Ls < float(L_thr) - 1e-9])   # decoupled below the first fibre-momentum level
    assert np.all(n_s >= n_b)
    b.step(Ls, n_s, where="post", color=C["orange"], lw=1.8, label=r"Berger sphere, $s=0.1$")
    b.step(Ls, n_b, where="post", color=C["blue"], lw=1.2, ls="--", label=r"round $S^2$ of radius $R_2/2$")
    b.axvline(float(L_thr), color=C["green"], lw=0.8, ls="-.")
    b.text(float(L_thr) + 1.5, 4, r"$2+1/s^2$: first fibre-" + "\nmomentum level", fontsize=8, color=C["green"])
    b.set(xlabel=r"$L$ (in units of $1/R_2^2$)", ylabel=r"number of eigenvalues $\leq L$ (with multiplicity)",
          title="Below that level the spectrum is exactly that of the base")
    b.legend(frameon=False, fontsize=8, loc="upper left")
    for q in ax:
        q.grid(alpha=0.25)
    fig.savefig(path, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "../figures/su2_berger_levels.png")
    out.parent.mkdir(parents=True, exist_ok=True)
    draw(out)
    print("wrote", out)
