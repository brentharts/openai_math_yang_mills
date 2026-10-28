#!/usr/bin/env python3
"""Collect headline numbers from build/data/*.json into build/tex/summary.tex."""

from __future__ import annotations

import sys

import common as C


def need(name: str) -> dict:
    data = C.load_json(name)
    if data is None:
        print(f"missing build/data/{name}.json; run the model scripts first (make models)", file=sys.stderr)
        raise SystemExit(1)
    return data


def main() -> None:
    bloch = need("bloch_law")
    mc = need("on_model_mc")
    kt = need("kt_flow")
    hal = need("haldane_gap")

    A = bloch["A"]
    fits = mc["fits"]
    mu_l = kt["mu_lstar_to_pi"][-1]
    L12 = hal["gap"][hal["L"].index(12)]

    rows = [
        ["271", "Bloch law, lattice correction",
         r"ideal-magnon density; sympy series",
         rf"$A_0={A[0]:.5f}$, $A_1={A[1]:.5f}$, $A_2={A[2]:.5f}$; $A_1$ equals the coefficient quoted in the paper"],
        ["215", r"2D O($N$) mass scaling",
         "Wolff Monte Carlo; sympy scaling forms",
         "fitted slopes of $\\ln\\xi$: " + ", ".join(
             f"O({k}) {v['slope']:.2f} (asymptotic freedom {v['theory']:.2f})" for k, v in fits.items())],
        ["216", "XY model, BKT scaling",
         "Kosterlitz flow; sympy exact solution",
         rf"$\mu\ln\xi={mu_l[1]:.3f}\to\pi$ at $\mu={mu_l[0]}$; log exponent $\to1/8$"],
        ["268", "Spin-1 Haldane gap",
         "exact diagonalisation; sympy bond algebra",
         rf"$\Delta(12)={L12:.4f}$; fit $\Delta_\infty={hal['aitken']:.4f}$ (literature $\approx0.4105$)"],
    ]
    tex = C.table(
        ["entry", "topic", "method", "headline result"], rows, colspec="llp{2.9cm}p{6.2cm}",
        caption=("Headline results. Each script asserts its own consistency checks, "
                 "so the report only exists if they passed."),
    )
    C.write_tex("summary", tex)
    print("wrote build/tex/summary.tex")


if __name__ == "__main__":
    main()
