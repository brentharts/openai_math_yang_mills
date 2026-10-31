"""Tree-level gluon amplitudes in exact rational arithmetic, from the CICY repository's twistor section.

Companion repository: brentharts/CICY (GPL-3.0), README sections "Twistor theory and scattering amplitudes"
and "Amplitudes on surfaces" (pyCICY.twistor). This module does not copy its code. It builds kinematics
from rational spinors, writes the Parke-Taylor formula for maximally helicity-violating (MHV) colour-ordered
gluon amplitudes, and checks, with Fractions and sympy and no tolerance:

  * momentum conservation as a linear condition (the columns of the antiholomorphic spinors lie in the kernel of
    the 2 x n matrix of holomorphic ones), the Schouten identity, cyclicity and reflection (-1)^n;
  * the U(1) decoupling identity and a Kleiss-Kuijf relation;
  * the soft-gluon factorisation, an algebraic identity A_n = S . A_{n-1};
  * the analytic input to BCFW recursion: after a shift of one holomorphic spinor the amplitude falls as z^-2, and
    the residues at its two poles sum to A(0) with no term at infinity;
  * the rank of the matrix of all colour orderings across kinematic points is (n-2)!, the Kleiss-Kuijf count (not
    the BCJ count (n-3)!, which a constant-coefficient rank cannot see);
  * the decorated permutations of [n] number sum_j n!/j!.

Not tested: loop amplitudes (the repository says none are computed), beyond-MHV amplitudes, the three-point
building blocks and their sign conventions (the repository found a (-1)^n sign error there), the
relation of any of this to a quantum Yang-Mills theory. These are classical-limit, perturbative statements
about gluon scattering at tree level; they say nothing about existence or a mass gap.
"""
import itertools
import random
from fractions import Fraction as Fr

import sympy as sp

from eqlib import Equation, check, register

NAME = "cicy_amplitudes"
ENTRY = "cicy-gauge"
REPO = "brentharts/CICY"
SRC = f"{REPO}: README.md section 'Twistor theory and scattering amplitudes' (pyCICY.twistor)"


def br(a, b):
    """The 2 x 2 bracket <ab> = a1 b2 - a2 b1 of two spinors."""
    return a[0] * b[1] - a[1] * b[0]


def kinematics(n, seed):
    """Rational spinors lam, lamt (n of each) with sum_i lam_i lamt_i^T = 0 (complexified momenta)."""
    rnd = random.Random(seed)
    while True:
        lam = [(Fr(rnd.randint(-12, 12)), Fr(rnd.randint(-12, 12))) for _ in range(n)]
        M = sp.Matrix([[sp.Rational(l[a].numerator, l[a].denominator) for l in lam] for a in range(2)])
        ker = M.nullspace()
        if M.rank() == 2 and len(ker) == n - 2 and all(br(lam[i], lam[j]) != 0 for i in range(n) for j in range(i + 1, n)):
            break
    lamt = [[None, None] for _ in range(n)]
    for b in range(2):
        coef = [sp.Rational(rnd.randint(-5, 5)) for _ in ker]
        v = sum((c * k for c, k in zip(coef, ker)), sp.zeros(n, 1))
        for i in range(n):
            lamt[i][b] = Fr(int(v[i].p), int(v[i].q))
    return lam, [tuple(x) for x in lamt]


def parke_taylor(lam, order, neg):
    """MHV colour-ordered amplitude <ab>^4 / prod <i i+1> for the cyclic `order` of leg labels, with negative helicities `neg`."""
    a, b = neg
    num = br(lam[a], lam[b]) ** 4
    den = Fr(1)
    for k in range(len(order)):
        den *= br(lam[order[k]], lam[order[(k + 1) % len(order)]])
    return num / den


def exact_rank(rows):
    """Rank of a matrix of Fractions by Gaussian elimination (exact)."""
    m = [list(r) for r in rows]
    rank, ncol = 0, len(m[0]) if m else 0
    for c in range(ncol):
        piv = next((i for i in range(rank, len(m)) if m[i][c] != 0), None)
        if piv is None:
            continue
        m[rank], m[piv] = m[piv], m[rank]
        for i in range(rank + 1, len(m)):
            if m[i][c] != 0:
                f = m[i][c] / m[rank][c]
                m[i] = [a - f * b for a, b in zip(m[i], m[rank])]
        rank += 1
    return rank


def shuffles(xs, ys):
    if not xs:
        yield tuple(ys)
        return
    if not ys:
        yield tuple(xs)
        return
    for rest in shuffles(xs[1:], ys):
        yield (xs[0],) + rest
    for rest in shuffles(xs, ys[1:]):
        yield (ys[0],) + rest


def decorated_k(sigma):
    """Number of anti-exceedances of a permutation (sigma(i) < i), fixed points counted separately by colour."""
    return sum(1 for i, s in enumerate(sigma) if s < i)


def build_and_verify():
    eq_kin = register(NAME, Equation(
        key="ca_kinematics", entry=ENTRY, title="Momentum conservation as a linear condition; Schouten identity",
        source_kind="repo-doc", source=f"{SRC}, subsection 'Momentum conservation as a linear condition'",
        latex_str=(r"\sum_i\lambda_i\tilde\lambda_i=0\iff\tilde\lambda^{(b)}\in\ker[\lambda_1\cdots\lambda_n],\qquad "
                   r"\langle ij\rangle\langle kl\rangle+\langle ik\rangle\langle lj\rangle+\langle il\rangle\langle jk\rangle=0"),
        notes="With lambda and lambda-tilde independent (complexified or split signature) a kinematic point is constructed, "
              "not solved for. Here with rational spinors so every check below is an exact equality of rationals.",
    ))
    pts = {n: kinematics(n, seed=n) for n in range(4, 8)}
    ok_mom = True
    for n, (lam, lamt) in pts.items():
        for a in range(2):
            for b in range(2):
                ok_mom &= sum(lam[i][a] * lamt[i][b] for i in range(n)) == 0
    check(NAME, "sum_i lambda_i lambda-tilde_i = 0 exactly for n = 4, 5, 6, 7", ok_mom, eq=eq_kin)
    lam4, lamt4 = pts[4]
    ok_sch = all(br(lam4[i], lam4[j]) * br(lam4[k], lam4[l]) + br(lam4[i], lam4[k]) * br(lam4[l], lam4[j])
                 + br(lam4[i], lam4[l]) * br(lam4[j], lam4[k]) == 0 for i, j, k, l in itertools.permutations(range(4)))
    check(NAME, "the Schouten identity holds for all 24 orderings of four spinors", ok_sch, eq=eq_kin)
    ok_ward = all(sum(br(lam4[k], lam4[i]) * br(lamt4[i], lamt4[l]) for i in range(4)) == 0 for k in range(4) for l in range(4))
    check(NAME, "sum_i <k i>[i l] = 0 follows from momentum conservation, without being imposed", ok_ward, eq=eq_kin)

    eq_pt = register(NAME, Equation(
        key="ca_parke_taylor", entry=ENTRY, title="Parke-Taylor MHV amplitude: cyclicity, reflection, U(1) decoupling, Kleiss-Kuijf",
        source_kind="textbook", source=f"{SRC}; Parke and Taylor, Phys. Rev. Lett. 56, 2459 (1986)",
        latex_str=(r"A_n^{\rm MHV}(1,\dots,n)=\frac{\langle ab\rangle^4}{\langle12\rangle\langle23\rangle\cdots\langle n1\rangle},\quad "
                   r"A(n,\dots,1)=(-1)^nA(1,\dots,n),\quad \sum_{\rm cyc}A(1,\dots,n)=0"),
        notes="The colour-ordered tree amplitude of n gluons with two negative helicities a, b. The checks use exact "
              "rationals at random complex kinematics, so a failure is a failure, never a tolerance.",
    ))
    cyc = refl = u1 = kk = True
    for n in (5, 6, 7):
        lam, _ = pts[n]
        neg = (0, 1)
        base = list(range(n))
        A0 = parke_taylor(lam, base, neg)
        cyc &= all(parke_taylor(lam, base[k:] + base[:k], neg) == A0 for k in range(n))
        refl &= parke_taylor(lam, base[::-1], neg) == (-1) ** n * A0
        # U(1) decoupling: insert leg n-1 in every slot among the others, in cyclic order
        others = base[:-1]
        u1 &= sum(parke_taylor(lam, others[:k] + [n - 1] + others[k:], neg) for k in range(len(others))) == 0
    check(NAME, "cyclicity: all n rotations of the colour ordering give the same amplitude (n = 5, 6, 7)", cyc, eq=eq_pt)
    check(NAME, "reflection: A(n, ..., 1) = (-1)^n A(1, ..., n) (n = 5, 6, 7)", refl, eq=eq_pt)
    check(NAME, "U(1) decoupling: the sum over insertions of one leg into the cyclic order of the others is exactly 0", u1, eq=eq_pt)
    # Kleiss-Kuijf: A(1, alpha, n, beta) = (-1)^|beta| sum over alpha shuffle beta^T of A(1, sigma, n)
    for n, split in ((5, 2), (6, 2), (6, 3), (7, 3)):
        lam, _ = pts[n]
        neg = (1, 2)
        last = n - 1
        mid = list(range(1, n - 1))
        alpha, beta = mid[:split - 1], mid[split - 1:]
        lhs = parke_taylor(lam, [0] + alpha + [last] + beta, neg)
        rhs = (-1) ** len(beta) * sum(parke_taylor(lam, [0] + list(s) + [last], neg) for s in shuffles(alpha, beta[::-1]))
        kk &= lhs == rhs
    check(NAME, "Kleiss-Kuijf: A(1, alpha, n, beta) = (-1)^|beta| sum_{alpha shuffle beta^T} A(1, sigma, n) for four splittings", kk, eq=eq_pt)

    eq_soft = register(NAME, Equation(
        key="ca_soft", entry=ENTRY, title="Soft-gluon factorisation of the Parke-Taylor amplitude",
        source_kind="derived", source=f"{SRC}, subsection 'Soft factorisation is linear algebra on a graph Laplacian' (only the algebraic identity is used)",
        latex_str=(r"A_n(\dots,n-1,n,1,\dots)=S(n-1,n,1)\,A_{n-1}(\dots,n-1,1,\dots),\qquad "
                   r"S(a,s,b)=\frac{\langle ab\rangle}{\langle as\rangle\langle sb\rangle}"),
        notes="For the MHV formula this is an exact algebraic identity for any lambda_n. Scaling lambda_n by epsilon gives the "
              "epsilon^-2 behaviour of a soft positive-helicity gluon. The repository's graph-Laplacian formulation of "
              "soft factorisation is not rederived here.",
    ))
    soft = True
    for n in (5, 6, 7):
        lam, _ = pts[n]
        neg = (0, 1)
        An = parke_taylor(lam, list(range(n)), neg)
        Am = parke_taylor(lam, list(range(n - 1)), neg)
        S = br(lam[n - 2], lam[0]) / (br(lam[n - 2], lam[n - 1]) * br(lam[n - 1], lam[0]))
        soft &= An == S * Am
        eps = Fr(1, 7)
        lam2 = lam[:n - 1] + [(lam[n - 1][0] * eps, lam[n - 1][1] * eps)]
        soft &= parke_taylor(lam2, list(range(n)), neg) == An / eps ** 2
    check(NAME, "A_n = S(n-1, n, 1) A_{n-1} exactly, and scaling lambda_n by eps multiplies A_n by eps^-2 (n = 5, 6, 7)", soft, eq=eq_soft)

    eq_bcfw = register(NAME, Equation(
        key="ca_bcfw_analytic", entry=ENTRY, title="The analytic input to BCFW: z^-2 falloff and a two-term residue sum",
        source_kind="derived", source=f"{SRC}, subsection 'BCFW against Parke-Taylor' (only the analytic structure is rederived)",
        latex_str=(r"|\hat s\rangle=|s\rangle-z|r\rangle:\quad A(z)=O(z^{-2}),\qquad "
                   r"A(0)=-\sum_{p}\operatorname*{Res}_{z=z_p}\frac{A(z)}{z}\ \ (\text{two poles})"),
        notes="Shift one holomorphic spinor of a positive-helicity leg. The MHV amplitude then has simple poles only where "
              "the two brackets containing that leg vanish, and falls as z^-2, so A(0) is the sum of two residues with no "
              "boundary term. This is the analytic fact BCFW recursion rests on. The three-point factors that make each "
              "residue a product of lower amplitudes, and their sign convention, are not checked here (the repository "
              "found a (-1)^n sign error exactly there).",
    ))
    z = sp.symbols("z")
    bc = True
    decay = True
    for n in (5, 6, 7):
        lam, _ = pts[n]
        R = lambda f: sp.Rational(f.numerator, f.denominator)
        L = [(R(a), R(b)) for a, b in lam]
        r, s = 1, n - 1                      # shift lambda_s -> lambda_s - z lambda_r; s is positive, r is not adjacent to s
        L[s] = (L[s][0] - z * L[r][0], L[s][1] - z * L[r][1])
        num = br(L[1], L[2]) ** 4            # negative helicities 1, 2, neither of them s
        den = sp.Integer(1)
        for k in range(n):
            den *= br(L[k], L[(k + 1) % n])
        A = sp.cancel(num / den)
        numer, denom = sp.fraction(A)
        decay &= sp.degree(denom, z) - sp.degree(numer, z) == 2
        poles = sp.roots(sp.Poly(denom, z))
        total = sum(sp.residue(A / z, z, p) for p in poles) + sp.residue(A / z, z, 0)
        bc &= len(poles) == 2 and sp.simplify(total) == 0
    check(NAME, "A(z) falls as z^-2 under the shift (n = 5, 6, 7): numerator degree 0, denominator degree 2", decay, eq=eq_bcfw)
    check(NAME, "A(z)/z has exactly two finite poles besides z = 0 and the residues at all three sum to exactly 0", bc, eq=eq_bcfw)

    eq_rank = register(NAME, Equation(
        key="ca_kk_rank", entry=ENTRY, title="The rank of all colour orderings across kinematic points is (n-2)!",
        source_kind="repo-doc", source=f"{SRC}, subsection 'What a rank computation can and cannot see'",
        latex_str=r"\operatorname{rank}\bigl[A(1,\sigma)\bigr]_{\sigma\in S_{n-1},\ \rm kinematic\ points}=(n-2)!\quad(\text{not }(n-3)!)",
        notes="A rank across kinematic points sees only relations with constant coefficients. Kleiss-Kuijf coefficients are "
              "+-1, so it counts (n-2)!. BCJ coefficients are built from Mandelstam invariants, vary from point to point, and "
              "cut the count to (n-3)! only if fed in as input; reporting (n-2)! with that explanation is the honest result.",
    ))
    for n, reps in ((5, 12), (6, 40)):
        orders = [(0,) + p for p in itertools.permutations(range(1, n))]
        rows = []
        for t in range(reps):
            lam, _ = kinematics(n, seed=100 * n + t)
            rows.append([parke_taylor(lam, list(o), (1, 2)) for o in orders])
        rank = exact_rank(rows)
        check(NAME, f"n = {n}: rank of the {len(orders)} orderings over {reps} kinematic points is {sp.factorial(n - 2)} = (n-2)!",
              rank == sp.factorial(n - 2), eq=eq_rank)

    eq_cells = register(NAME, Equation(
        key="ca_decorated", entry=ENTRY, title="Decorated permutations: the cells of the positive Grassmannian",
        source_kind="repo-doc", source=f"{SRC}, subsection 'The positive Grassmannian'; Postnikov, arXiv:math/0609764",
        latex_str=r"\#\{\text{decorated permutations of }[n]\}=\sum_\sigma2^{\mathrm{fix}(\sigma)}=\sum_{j=0}^n\frac{n!}{j!}\quad(16,\,65,\,326,\,1957)",
        notes="Counted by enumeration: permutations of [n] with each fixed point coloured one of two ways. The bijection to "
              "cells of the non-negative Grassmannian is Postnikov's and is quoted. The k-grading is not checked.",
    ))
    totals = {}
    for n in (3, 4, 5, 6):
        totals[n] = sum(2 ** sum(1 for i, s in enumerate(p) if s == i) for p in itertools.permutations(range(n)))
    check(NAME, "the totals for n = 3, 4, 5, 6 are 16, 65, 326, 1957",
          [totals[n] for n in (3, 4, 5, 6)] == [16, 65, 326, 1957], eq=eq_cells)
    check(NAME, "each total equals sum_{j=0}^{n} n!/j!",
          all(totals[n] == sum(sp.factorial(n) // sp.factorial(j) for j in range(n + 1)) for n in totals), eq=eq_cells)
