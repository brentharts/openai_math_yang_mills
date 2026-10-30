#!/usr/bin/env python3
"""Critical-strip figure and a zero-counting check for ZETA_ZEROS.md (not a repository entry).

Usage:  python3 tools/zeta_strip_figure.py [--nzeros=N] [--figures]

What it does
------------
1. Computes the first N nontrivial zeros of zeta with mpmath (default N = 200).
2. Compares the counting function N(t) with the smooth term used in the Riemann Research Laboratory,
       M(t) = (t/2pi) log(t/2pi) - t/2pi + 7/8,
   and reports E(t) = N(t) - M(t) at the midpoints between consecutive zeros.
3. Checks the functional-equation orbit {rho, 1-rho, conj(rho), 1-conj(rho)} for a synthetic displaced zero
   at delta = 3/8 (the largest displacement not excluded by a zero-free region Re s > 7/8; see the page).
4. With --figures, writes figures/zeta_fig1_strip.png (strip with 1/2, 7/8, 11/12, 1 and the delta range)
   and figures/zeta_fig2_counting.png (N(t), M(t) and E(t)).

What it is not
--------------
mpmath's zetazero locates the n-th zero on the critical line; finding zeros there is not a test of the
Riemann hypothesis, and this script does not claim one.  The counting comparison is a sanity check on the
indexing and on M(t); it proves nothing about E(t) for large t.
"""
from __future__ import annotations

import math
import sys
import warnings

import mpmath as mp

warnings.filterwarnings("ignore")


def M(t: float) -> float:
    return t / (2 * math.pi) * math.log(t / (2 * math.pi)) - t / (2 * math.pi) + 7.0 / 8.0


def zeros(n: int) -> list[float]:
    mp.mp.dps = 20
    return [float(mp.zetazero(k).imag) for k in range(1, n + 1)]


def orbit(beta: float, gamma: float):
    return [complex(beta, gamma), complex(1 - beta, -gamma), complex(beta, -gamma), complex(1 - beta, gamma)]


def run(nz: int, figures: bool) -> None:
    print(f"== zeta zeros: first {nz} zeros, counting function against M(t) ==")
    g = zeros(nz)
    assert all(g[i] < g[i + 1] for i in range(nz - 1)), "zero ordinates must increase"
    assert abs(g[0] - 14.134725) < 1e-5 and abs(g[1] - 21.022040) < 1e-5

    # N(t) at the midpoint between zero k and k+1 equals k.
    mids = [(g[k] + g[k + 1]) / 2 for k in range(nz - 1)]
    E = [(k + 1) - M(mids[k]) for k in range(nz - 1)]
    print(f"gamma_1 = {g[0]:.6f}, gamma_{nz} = {g[-1]:.6f}")
    print(f"E(t) = N(t) - M(t) at midpoints: min {min(E):+.3f}, max {max(E):+.3f}, "
          f"mean {sum(E) / len(E):+.3f}")
    print("(M(t) includes the constant 7/8, as in the laboratory's definition.)")
    # The step function N(t) jumps by one at each zero, so the continuous sup of |E| is larger.
    import numpy as _np
    ts = _np.linspace(g[0] - 1, g[-1] + 1, 4000)
    Nt = _np.searchsorted(_np.array(g), ts, side="right")
    Ec = Nt - _np.array([M(t) for t in ts])
    print(f"E(t) on a fine grid (just before and after jumps): min {Ec.min():+.3f}, max {Ec.max():+.3f}")

    # Orbit check for the largest displacement not excluded by a zero-free region Re s > 7/8.
    d = 3 / 8
    q = orbit(0.5 + d, g[0])
    assert abs(q[1] - (1 - q[0])) < 1e-15 and abs(q[2] - q[0].conjugate()) < 1e-15
    assert abs(q[3] - (1 - q[2])) < 1e-15
    assert max(z.real for z in q) == 0.875 and min(z.real for z in q) == 0.125
    print(f"orbit at delta = 3/8 on gamma_1: real parts {sorted(set(z.real for z in q))}")

    if figures:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        ax.axvspan(0, 1, color="#eef2f7")
        ax.axvspan(7 / 8, 1, color="#d9534f", alpha=0.25, label="excluded if the 7/8 claim holds (and by symmetry Re s < 1/8)")
        ax.axvspan(0, 1 / 8, color="#d9534f", alpha=0.25)
        ax.axvspan(1 / 8, 7 / 8, color="#5bc0de", alpha=0.18, label="not excluded: displacement 0 < delta <= 3/8")
        for x, lab, st, dy in [(0.0, "0", "-", 0), (1 / 8, "1/8", "--", 0), (0.5, "1/2", "-", 0),
                               (7 / 8, "7/8", "--", 0), (11 / 12, "11/12 (weaker stage)", ":", 0.07), (1.0, "1", "-", 0)]:
            ax.axvline(x, color="k", ls=st, lw=1.1)
            ax.text(x, 1.02 + dy, lab, ha="center", va="bottom", fontsize=8.5, transform=ax.get_xaxis_transform())
        ax.plot([0.5] * 30, g[:30], "k.", ms=5, label="first 30 zeros (on the line)")
        ax.set_xlim(-0.05, 1.05)
        ax.set_ylim(0, g[29] + 3)
        ax.set_xlabel("Re s")
        ax.set_ylabel("Im s")
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), fontsize=8, frameon=False)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        fig.savefig("figures/zeta_fig1_strip.png", dpi=150)
        plt.close(fig)

        import numpy as np
        ts = np.linspace(g[0] - 1, g[-1] + 1, 4000)
        Nt = np.searchsorted(np.array(g), ts, side="right")
        Mt = np.array([M(t) for t in ts])
        fig, (a, b) = plt.subplots(2, 1, figsize=(7.2, 5.2), sharex=True)
        a.step(ts, Nt, where="post", label="N(t) (zeros counted)", lw=1.2)
        a.plot(ts, Mt, "--", label="M(t) with 7/8", lw=1.2)
        a.set_ylabel("count")
        a.legend(fontsize=8)
        b.plot(ts, Nt - Mt, lw=0.9)
        b.axhline(0, color="k", lw=0.6)
        b.set_ylabel("E(t) = N(t) - M(t)")
        b.set_xlabel("t")
        fig.tight_layout()
        fig.savefig("figures/zeta_fig2_counting.png", dpi=150)
        plt.close(fig)
        print("wrote figures/zeta_fig1_strip.png, figures/zeta_fig2_counting.png")


if __name__ == "__main__":
    nz = 200
    for a in sys.argv[1:]:
        if a.startswith("--nzeros="):
            nz = int(a.split("=")[1])
    run(nz, "--figures" in sys.argv)
