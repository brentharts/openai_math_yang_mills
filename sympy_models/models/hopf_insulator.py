"""The three-dimensional Hopf insulator, from the companion repository brentharts/quantum-lattice.

The repository (HOPF_INSULATOR.md, HOPF_INSULATOR_FIGURES.md, extras/ql_hopf_figures.py) is GPL-3.0 and
wraps a vendored pyqula; this fork is Apache-2.0, so this module copies none of it. It rederives in exact
algebra the closed forms that the repository's figures rely on, for the model H = d . sigma with
d = Hopf(z1, z2), z1 = sin kx + i sin ky, z2 = sin kz + i (cos kx + cos ky + cos kz + h):

  * the gap closes (z1 = z2 = 0) only at h = -3, -1, 1, 3, where some time-reversal-invariant
    momentum (TRIM) has m = cos kx + cos ky + cos kz + h = 0;
  * the Hopf invariant is the TRIM sum  sum_i chi_i [m_i < 0]  with chi = cos kx cos ky cos kz,
    which gives 0, -1, +2, -1, 0 for h > 3, 1 < h < 3, |h| < 1, -3 < h < -1, h < -3;
  * the multiplicities of the TRIM by the value of sum cos are 1, 3, 3, 1;
  * on the line kx = ky = pi, |z2|^2 = 1 + (h - 2)^2 + 2 (h - 2) cos kz.

Quoted, not rederived: that the Hopf invariant of the map d equals the linking number of the preimages
of two points (Moore-Ran-Wen 2008) and that it is the invariant of the insulator class (Deng, Wang,
Shen, Duan 2013). The repository checks the linking numbers numerically (-1 at h = 2, +2 at h = 0.5).

Not tested: the numerical Berry curvature and the preimage tracing (the repository's figures), the
degree integral over the zone, and any statement about interacting systems. The gap of this
insulator is a single-particle band gap of a lattice model. It has no relation to the Yang-Mills
mass gap, and nothing here is a statement about Yang-Mills existence.
"""
import itertools
from collections import Counter

import sympy as sp

from eqlib import Equation, check, register

NAME = "hopf_insulator"
ENTRY = "hopf-insulator"
REPO = "brentharts/quantum-lattice"
DOC = f"{REPO}: HOPF_INSULATOR_FIGURES.md"
MRW = "Moore, Ran and Wen, arXiv:0804.4527"
DWSD = "Deng, Wang, Shen and Duan, Phys. Rev. B 88, 201105 (2013), arXiv:1307.7206"

TRIM = list(itertools.product((0, 1), repeat=3))        # k_i in {0, pi}


def cos_sum(k):
    return sum(1 if ki == 0 else -1 for ki in k)


def chi(k):
    return (-1) ** sum(k)


def invariant(h):
    """Sum over the eight TRIM of chi [m < 0], for exact rational h not at a transition."""
    return sum(chi(k) for k in TRIM if cos_sum(k) + h < 0)


def build_and_verify():
    kx, ky, kz, h = sp.symbols("kx ky kz h", real=True)
    m = sp.cos(kx) + sp.cos(ky) + sp.cos(kz) + h

    # ---- where the gap closes -----------------------------------------------------------------
    eq_gap = register(NAME, Equation(
        key="hi_gap_closing", entry=ENTRY, title="The Hopf insulator is gapped except at h = -3, -1, 1, 3",
        source_kind="derived", source=f"derived here; compare {REPO}: HOPF_INSULATOR.md and {DOC} figure A",
        latex_str=(r"z_1=\sin k_x+i\sin k_y,\ z_2=\sin k_z+i\,m,\ m=\sum_i\cos k_i+h;"
                   r"\qquad z_1=z_2=0\iff k\in\mathrm{TRIM},\ m=0\iff h\in\{-3,-1,1,3\}"),
        notes="H = d . sigma with d the unit Hopf vector, so the single-particle bands are +-1 and the only way for the "
              "gap to close is for d to become undefined, i.e. z1 = z2 = 0. That is a band-gap statement for a "
              "non-interacting lattice model. It is unrelated to a Yang-Mills mass gap.",
    ))
    # z1 = z2 = 0 needs sin kx = sin ky = sin kz = 0, i.e. every k_i in {0, pi} (mod 2 pi), and m = 0 there.
    closing = sorted({-cos_sum(k) for k in TRIM})
    check(NAME, "sin k = 0 at k = 0 and pi, where cos k = +1 and -1: the eight TRIM are the only points with z1 = 0 and sin kz = 0",
          all(sp.sin(a * sp.pi) == 0 and sp.cos(a * sp.pi) == (-1) ** a for a in (0, 1)) and len(TRIM) == 8, eq=eq_gap)
    check(NAME, "m = 0 at a TRIM exactly for h = -(sum of the three cosines) in {-3, -1, 1, 3}",
          closing == [-3, -1, 1, 3], eq=eq_gap)
    check(NAME, "away from those four values every TRIM has m != 0 (examples h = -4, -2, 0, 2, 4, 1/2)",
          all(cos_sum(k) + hv != 0 for k in TRIM for hv in (-4, -2, 0, 2, 4, sp.Rational(1, 2))), eq=eq_gap)

    # ---- the TRIM formula for the invariant ---------------------------------------------------
    mult = Counter(cos_sum(k) for k in TRIM)
    sign_by_sum = {s: {chi(k) for k in TRIM if cos_sum(k) == s} for s in mult}
    regimes = [(sp.Rational(7, 2), 0), (2, -1), (0, 2), (-2, -1), (-sp.Rational(7, 2), 0)]
    eq_trim = register(NAME, Equation(
        key="hi_trim_formula", entry=ENTRY, title="The Hopf invariant as a sum over the eight TRIM",
        source_kind="derived", source=f"derived here; checked against the repository's degree integral in {DOC} figure A",
        latex_str=(r"N_H(h)=\sum_{k\in\mathrm{TRIM}}\chi(k)\,\Theta(-m(k)),\quad \chi=\cos k_x\cos k_y\cos k_z;"
                   r"%%BR%%N_H=0,-1,+2,-1,0\ \text{for}\ h>3,\ 1<h<3,\ |h|<1,\ -3<h<-1,\ h<-3"),
        cites=["moore-ran-wen", "deng-wang-shen-duan"],
        notes="Each TRIM is a point where z1 = 0 and d is the south pole, with orientation chi, so the invariant counts "
              "signed TRIM with m < 0. This closed form was found for the repository's figures and checked against "
              "its numerical degree integral at 13 values of h; here it is evaluated exactly. The identification of "
              "this number with the linking number of preimages is quoted from " + MRW + " and " + DWSD + ".",
    ))
    check(NAME, "the eight TRIM have sum cos k_i = 3, 1, -1, -3 with multiplicities 1, 3, 3, 1",
          [mult[s] for s in (3, 1, -1, -3)] == [1, 3, 3, 1], eq=eq_trim)
    check(NAME, "chi = cos kx cos ky cos kz equals +1, -1, +1, -1 on the four classes (sum cos = 3, 1, -1, -3)",
          [sign_by_sum[s] for s in (3, 1, -1, -3)] == [{1}, {-1}, {1}, {-1}], eq=eq_trim)
    check(NAME, "the TRIM sum gives 0, -1, +2, -1, 0 in the five regimes (h = 7/2, 2, 0, -2, -7/2)",
          [invariant(hv) for hv, _ in regimes] == [v for _, v in regimes], eq=eq_trim)
    check(NAME, "constant within each regime: the same value at six rational h per regime",
          all(len({invariant(hv + d) for d in (sp.Rational(-1, 4), sp.Rational(-1, 8), 0, sp.Rational(1, 8),
                                               sp.Rational(1, 4), sp.Rational(1, 3))}) == 1
              for hv, _ in regimes), eq=eq_trim)
    R = sp.Rational
    jumps = [invariant(h0 - R(1, 4)) - invariant(h0 + R(1, 4)) for h0 in (3, 1, -1, -3)]
    check(NAME, "going down in h, the invariant changes by -1, +3, -3, +1 across h = 3, 1, -1, -3: "
                "chi times the multiplicity of the TRIM class whose m changes sign",
          jumps == [-1, 3, -3, 1]
          and jumps == [sign_by_sum[s_].copy().pop() * mult[s_] for s_ in (-3, -1, 1, 3)], eq=eq_trim)

    # ---- |z2|^2 on the line kx = ky = pi ------------------------------------------------------
    z2sq = sp.sin(kz) ** 2 + m.subs({kx: sp.pi, ky: sp.pi}) ** 2
    eq_line = register(NAME, Equation(
        key="hi_line_z2", entry=ENTRY, title="|z2|^2 along the line kx = ky = pi",
        source_kind="derived", source=f"derived here; the curvature statement is from {DOC} figure C",
        latex_str=r"|z_2|^2\big|_{k_x=k_y=\pi}=1+(h-2)^2+2(h-2)\cos k_z",
        notes="On the lines z1 = 0 the Berry curvature of the repository's figure C is 4/|z2|^2, so the closed form "
              "for |z2|^2 gives the curvature along the line; it equals 4 only at h = 2, and the relative error of the "
              "repository's grid value shrinks with the grid size. The curvature itself is not recomputed here.",
    ))
    check(NAME, "|z2|^2 = 1 + (h - 2)^2 + 2 (h - 2) cos kz on kx = ky = pi",
          sp.simplify(sp.expand(z2sq) - (1 + (h - 2) ** 2 + 2 * (h - 2) * sp.cos(kz))) == 0, eq=eq_line)
    check(NAME, "at h = 2, |z2|^2 = 1 for every kz, so the curvature 4/|z2|^2 is exactly 4",
          sp.simplify(z2sq.subs(h, 2)) == 1, eq=eq_line)
