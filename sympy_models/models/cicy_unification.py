"""One-loop gauge couplings: the Standard Model and MSSM coefficients, and what unification predicts.

Companion repository: brentharts/CICY (GPL-3.0), README section "The number 137, and what unification
actually predicts" (theories/couplings). This module does not copy its code. It derives the one-loop
coefficients b_i from the field content with the textbook formula, eliminates the unified coupling and
scale from the three running equations, and compares with the numbers that section reports:

  * b = (41/10, -19/6, -7) for the Standard Model with three generations and one Higgs doublet, and
    b = (33/5, 1, -3) for the MSSM, derived from group-theory data (Dynkin indices and hypercharges), with
    the SU(3) coefficient the same asymptotic-freedom number as the library's Yang-Mills law (11 minus the
    flavour term);
  * eliminating alpha_G and M_G leaves alpha_3^{-1} = alpha_2^{-1} - (alpha_1^{-1} - alpha_2^{-1})
    (b_2 - b_3)/(b_1 - b_2), exactly rational for rational inputs: 0.1169 for the MSSM and 0.0712 for the
    Standard Model, against the measured 0.1181 (the MSSM value is about 1% low);
  * a complete generation adds 2 to every b_i, so it cancels from every difference: the prediction is
    independent of the number of generations (exact);
  * the prediction is sensitive to incomplete multiplets: a second Higgs pair or one vector-like colour triplet
    pair moves it by a large factor.

Not tested: thresholds, two-loop running, the dilaton, any compactification, and the claim that
heterotic models reach these spectra (the repository's own caveat: "no physical couplings, masses or
predictions"). Inputs are the measured alpha_em^{-1}(M_Z) = 127.95 and sin^2 theta_W = 0.2312 as the
repository's section uses them. This is perturbative running of a gauge coupling: it is about
asymptotic freedom at one loop and says nothing about a mass gap or Yang-Mills existence.
"""
from fractions import Fraction as Fr

import sympy as sp

from eqlib import Equation, check, register

NAME = "cicy_unification"
ENTRY = "cicy-gauge"
REPO = "brentharts/CICY"
SRC = f"{REPO}: README.md section 'The number 137, and what unification actually predicts'"

ALPHA_INV = Fr(12795, 100)      # alpha_em^{-1}(M_Z), as in the repository's section
SIN2 = Fr(2312, 10000)          # sin^2 theta_W(M_Z)
MEASURED_A3 = Fr(1181, 10000)


def sm_b(ng=3, higgs=1):
    """(b1, b2, b3) of the Standard Model: b = -11/3 C2(G) + 2/3 sum_Weyl T + 1/3 sum_complex scalars T, U(1) as GUT-normalised."""
    # per generation, Weyl fermions (colour, isospin, Y): Q (3,2,1/6), u^c (3,1,-2/3), d^c (3,1,1/3), L (1,2,-1/2), e^c (1,1,1)
    fer = [(3, 2, Fr(1, 6)), (3, 1, Fr(-2, 3)), (3, 1, Fr(1, 3)), (1, 2, Fr(-1, 2)), (1, 1, Fr(1))]
    b3 = Fr(-11) + Fr(2, 3) * ng * sum((Fr(1, 2) if c == 3 else 0) * i for c, i, _ in fer)
    b2 = Fr(-22, 3) + Fr(2, 3) * ng * sum((Fr(1, 2) if i == 2 else 0) * c for c, i, _ in fer) + Fr(1, 3) * higgs * Fr(1, 2)
    y2 = sum(c * i * y * y for c, i, y in fer)
    b1 = Fr(3, 5) * (Fr(2, 3) * ng * y2 + Fr(1, 3) * higgs * 2 * Fr(1, 4))
    return (b1, b2, b3)


def mssm_b(ng=3, higgs_pairs=1, triplet_pairs=0):
    """MSSM: b = -3 C2(G) + sum over chiral multiplets of T(R); U(1) GUT-normalised."""
    chir = [(3, 2, Fr(1, 6)), (3, 1, Fr(-2, 3)), (3, 1, Fr(1, 3)), (1, 2, Fr(-1, 2)), (1, 1, Fr(1))]
    b3 = Fr(-9) + ng * sum((Fr(1, 2) if c == 3 else 0) * i for c, i, _ in chir) + triplet_pairs * 2 * Fr(1, 2)
    b2 = Fr(-6) + ng * sum((Fr(1, 2) if i == 2 else 0) * c for c, i, _ in chir) + higgs_pairs * 2 * Fr(1, 2)
    y2 = sum(c * i * y * y for c, i, y in chir)
    b1 = Fr(3, 5) * (ng * y2 + higgs_pairs * 2 * 2 * Fr(1, 4) + triplet_pairs * 2 * 3 * Fr(1, 9))
    return (b1, b2, b3)


def alpha3_inv(b, alpha_inv=ALPHA_INV, sin2=SIN2):
    a1 = Fr(3, 5) * (1 - sin2) * alpha_inv
    a2 = sin2 * alpha_inv
    b1, b2, b3 = b
    return a2 - (a1 - a2) * (b2 - b3) / (b1 - b2)


def build_and_verify():
    eq_b = register(NAME, Equation(
        key="cu_one_loop_b", entry=ENTRY, title="One-loop coefficients of the Standard Model and the MSSM from the field content",
        source_kind="textbook", source=f"{SRC}; the formula is standard (Georgi-Quinn-Weinberg 1974)",
        latex_str=(r"\mu\frac{d\alpha_i^{-1}}{d\mu}=-\frac{b_i}{2\pi},\quad b^{\rm SM}=\bigl(\tfrac{41}{10},-\tfrac{19}{6},-7\bigr),"
                   r"\quad b^{\rm MSSM}=\bigl(\tfrac{33}{5},1,-3\bigr)"),
        notes="Derived from Dynkin indices and hypercharges with three generations and the stated Higgs content. The "
              "SU(3) entry -7 is the library's asymptotic-freedom number: -11 plus 4/3 per generation from six quark "
              "flavours. One loop only; this is perturbative running, not a statement about a gap.",
    ))
    sm, ms = sm_b(), mssm_b()
    check(NAME, "Standard Model b = (41/10, -19/6, -7) from three generations and one Higgs doublet",
          sm == (Fr(41, 10), Fr(-19, 6), Fr(-7)), eq=eq_b)
    check(NAME, "MSSM b = (33/5, 1, -3) from three generations and two Higgs doublets",
          ms == (Fr(33, 5), Fr(1), Fr(-3)), eq=eq_b)
    nf = 6
    check(NAME, "b_3 = -(11 - 2 n_f/3) with n_f = 6 quark flavours: the library's SU(3) asymptotic-freedom number",
          sm[2] == -(11 - Fr(2, 3) * nf), eq=eq_b)
    check(NAME, "the pure-gauge limit (no matter) is b_3 = -11 for SU(3), b_2 = -22/3 for SU(2)",
          sm_b(ng=0, higgs=0)[2] == -11 and sm_b(ng=0, higgs=0)[1] == Fr(-22, 3), eq=eq_b)

    eq_p = register(NAME, Equation(
        key="cu_alpha3", entry=ENTRY, title="What one-loop unification predicts for alpha_3",
        source_kind="derived", source=SRC,
        latex_str=(r"\alpha_3^{-1}(M_Z)=\alpha_2^{-1}-(\alpha_1^{-1}-\alpha_2^{-1})\frac{b_2-b_3}{b_1-b_2},\qquad "
                   r"\alpha_3^{\rm MSSM}=0.117,\ \alpha_3^{\rm SM}=0.071\ \ (\text{measured }0.1181)"),
        notes="Eliminates alpha_G and M_G from the three one-loop equations; the only theory input is the b_i. With "
              "alpha_em^{-1} = 127.95 and sin^2 theta_W = 0.2312 as inputs. One loop, no thresholds.",
    ))
    a_ms, a_sm = 1 / alpha3_inv(ms), 1 / alpha3_inv(sm)
    check(NAME, "MSSM gives alpha_3 = 0.1169 (repository: 0.1170), about 1% (1.02%) below the measured 0.1181",
          abs(a_ms - Fr(1170, 10000)) < Fr(1, 2000) and abs(a_ms - MEASURED_A3) / MEASURED_A3 < Fr(3, 200), eq=eq_p)
    check(NAME, "Standard Model gives alpha_3 = 0.0712, about 40% below the measured value",
          abs(a_sm - Fr(712, 10000)) < Fr(1, 2000) and abs(a_sm - MEASURED_A3) / MEASURED_A3 > Fr(1, 3), eq=eq_p)

    eq_g = register(NAME, Equation(
        key="cu_generation_independence", entry=ENTRY, title="The unification prediction does not depend on the number of generations",
        source_kind="derived", source=SRC,
        latex_str=r"b_i\to b_i+2n_g\ \ \text{for all }i\ \Rightarrow\ b_i-b_j\ \text{unchanged},\quad \alpha_3\ \text{independent of }n_g",
        notes="A complete generation fills SU(5) multiplets and adds the same amount to every b_i, so it cancels from "
              "all differences. So the index the repository computes best, the chiral generation count, is exactly the "
              "quantity this prediction cannot see.",
    ))
    base = {ng: mssm_b(ng) for ng in (2, 3, 4, 5, 6)}
    check(NAME, "each added generation shifts every MSSM b_i by exactly 2",
          all(tuple(base[n + 1][i] - base[n][i] for i in range(3)) == (2, 2, 2) for n in (2, 3, 4, 5)), eq=eq_g)
    check(NAME, "alpha_3 is exactly the same rational number for 2, 3, 4, 5 and 6 generations",
          len({alpha3_inv(b) for b in base.values()}) == 1, eq=eq_g)
    check(NAME, "the same holds for the Standard Model spectrum (each generation adds 4/3 to b_3, 4/3 to b_2, 4/3 to b_1)",
          len({alpha3_inv(sm_b(ng)) for ng in (2, 3, 4, 5)}) == 1, eq=eq_g)

    eq_v = register(NAME, Equation(
        key="cu_vectorlike_sensitivity", entry=ENTRY, title="The prediction is sensitive to incomplete multiplets",
        source_kind="derived", source=SRC,
        latex_str=r"\Delta b=(\tfrac35,1,0)\ (\text{extra Higgs pair})\ \Rightarrow\ \alpha_3\approx0.78;\quad \Delta b=(\tfrac25,0,1)\ (3+\bar3)\ \Rightarrow\ \alpha_3\approx0.067",
        notes="Vector-like matter does not change the chiral index, so a compactification's index cannot see it, yet "
              "it moves alpha_3 by an order of magnitude. The repository quotes 0.786 and 0.067; the first is near a "
              "pole of the formula (alpha_3^{-1} is small) and varies at the third digit with the inputs.",
    ))
    hp = 1 / alpha3_inv(mssm_b(higgs_pairs=2))
    tp = 1 / alpha3_inv(mssm_b(triplet_pairs=1))
    check(NAME, "a second Higgs pair gives alpha_3 near 0.78 (repository 0.786)", Fr(76, 100) < hp < Fr(80, 100), eq=eq_v)
    check(NAME, "one vector-like colour triplet pair gives alpha_3 = 0.067 (repository 0.067)",
          abs(tp - Fr(67, 1000)) < Fr(1, 1000), eq=eq_v)
    check(NAME, "the shifts are Delta b = (3/5, 1, 0) and (2/5, 0, 1)",
          tuple(a - b for a, b in zip(mssm_b(higgs_pairs=2), ms)) == (Fr(3, 5), Fr(1), Fr(0))
          and tuple(a - b for a, b in zip(mssm_b(triplet_pairs=1), ms)) == (Fr(2, 5), Fr(0), Fr(1)), eq=eq_v)
