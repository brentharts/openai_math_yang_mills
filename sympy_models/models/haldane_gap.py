"""The periodic spin-one Haldane gap (entry 268).

Statements and the proof's scalar skeleton are transcribed from
  [H]  The periodic spin-one Haldane gap
       preprints/The-periodic-spin-one-Haldane-gap-September-24-2026/

The paper proves gamma_L > (4/105) log(80/79) for even L >= 60 and gamma_L > log(20)/784 for
even L >= 2304 by a purity bootstrap (Proposition 4.3) started from two finite computer-assisted
inputs (Propositions 5.4 and 5.5, with the thermal and trial bounds of Appendices A-C).

What this module checks, exactly (rational arithmetic):
  * the spin-one matrices satisfy the stated commutator and Casimir normalization;
  * the algebra of the bootstrap: f, G, the identity f(1-(1-u)^3) = G(u) u^2, the coupled
    identities (4.5)-(4.6), the closed form of u_j, and the squaring lemma on random weights;
  * the scalar hypotheses (4.10) for both parameter rows (4.18) and the two resulting bounds.

What it does NOT check: the finite inputs S(60, 105/4) > 7/8 etc. (Propositions 5.4-5.5), which are
the heavy part of the proof and rely on certified matrix-exponential enclosures. The small exact-
diagonalization check at the end (L = 4..10, numpy/scipy if installed) is a sanity check of the
Hamiltonian only; the theorem's explicit bounds start at L = 60.
"""
import random
from fractions import Fraction

import sympy as sp

from eqlib import Equation, check, register

NAME = "haldane_gap"
ENTRY = "268"
H = "[H] periodic spin-one Haldane gap paper"


def _spin1():
    sz = sp.diag(1, 0, -1)
    sp_ = sp.sqrt(2) * sp.Matrix([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
    sm = sp_.T
    sx = (sp_ + sm) / 2
    sy = (sp_ - sm) / (2 * sp.I)
    return sx, sy, sz


def build_and_verify():
    # ---- the model and the theorem ---------------------------------------------------------
    eq_alg = register(NAME, Equation(
        key="haldane_spin_algebra", entry=ENTRY,
        title="Spin-one normalization",
        source_kind="paper", source=f"{H}, Section 1, before eq. (1.1)",
        latex_str=(r"[S^\alpha,S^\beta]=i\sum_\eta\varepsilon_{\alpha\beta\eta}S^\eta,%%BR%%"
                   r"\sum_\alpha(S^\alpha)^2=2I_3"),
    ))
    sx, sy, sz = _spin1()
    comm = lambda a, b: a * b - b * a
    check(NAME, "spin-one matrices: [Sx,Sy] = i Sz",
          sp.simplify(comm(sx, sy) - sp.I * sz) == sp.zeros(3, 3), eq=eq_alg)
    check(NAME, "spin-one matrices: [Sy,Sz] = i Sx",
          sp.simplify(comm(sy, sz) - sp.I * sx) == sp.zeros(3, 3), eq=eq_alg)
    check(NAME, "spin-one matrices: [Sz,Sx] = i Sy",
          sp.simplify(comm(sz, sx) - sp.I * sy) == sp.zeros(3, 3), eq=eq_alg)
    check(NAME, "spin-one Casimir is S(S+1) I = 2 I",
          sp.simplify(sx**2 + sy**2 + sz**2 - 2 * sp.eye(3)) == sp.zeros(3, 3), eq=eq_alg)

    register(NAME, Equation(
        key="haldane_hamiltonian", entry=ENTRY,
        title="Periodic spin-one Heisenberg chain",
        source_kind="paper", source=f"{H}, eq. (1.1)",
        latex_str=(r"H_L=\sum_{j=0}^{L-1}\sum_{\alpha=x,y,z}S^\alpha_jS^\alpha_{j+1},"
                   r"\qquad S^\alpha_L=S^\alpha_0,\qquad L\ge4\ \text{even}"),
        notes="Coupling constant one, pure bilinear (no biquadratic term).",
    ))
    register(NAME, Equation(
        key="haldane_gap_definition", entry=ENTRY,
        title="Gap above the full ground energy",
        source_kind="paper", source=f"{H}, eq. (1.2)",
        latex_str=r"\gamma_L=E_1(L)-E_0(L),\qquad \Delta_1=\liminf_{L\to\infty,\ L\ \mathrm{even}}\gamma_L",
        notes="E_1 is the next DISTINCT energy, so ground-state degeneracy counts against the gap.",
    ))

    b1 = sp.Rational(4, 105) * sp.log(sp.Rational(80, 79))
    b2 = sp.log(20) / 784
    eq_thm = register(NAME, Equation(
        key="haldane_theorem", entry=ENTRY,
        title="Uniform gap with explicit constants",
        source_kind="paper", source=f"{H}, Theorem 1.1, eqs. (1.3)-(1.4)",
        latex_str=(r"\gamma_L>\frac{4}{105}\log\frac{80}{79}\ \ (L\ge60),%%BR%%"
                   r"\gamma_L>\frac{\log 20}{784}\ \ (L\ge2304),%%BR%%"
                   r"\Delta_1\ge\frac{\log20}{784}>0"),
        notes="The ground state is also unique for every even L >= 60. The bound is far below the "
              "numerical bulk gap quoted in the paper (0.41050, White-Huse DMRG); it proves positivity, "
              "not a sharp value.",
    ))

    # ---- bootstrap algebra --------------------------------------------------------------------
    u = sp.symbols("u", positive=True)
    f = u**2 / (2 * (1 - u) ** 2)
    G = sp.Rational(1, 2) * ((1 - u) ** -1 + (1 - u) ** -2 + (1 - u) ** -3) ** 2
    eq_fG = register(NAME, Equation(
        key="haldane_squaring_functions", entry=ENTRY,
        title="Squaring function f and its composite G",
        source_kind="paper", source=f"{H}, eqs. (4.2), (4.8), (4.9)",
        expr=sp.Eq(sp.Function("f")(sp.Symbol("u")), f),
        latex_str=(r"f(u)=\frac{u^{2}}{2(1-u)^{2}},%%BR%%"
                   r"G(u)=\frac12\Bigl((1-u)^{-1}+(1-u)^{-2}+(1-u)^{-3}\Bigr)^{2},%%BR%%"
                   r"f\bigl(1-(1-u)^{3}\bigr)=G(u)\,u^{2}"),
    ))
    check(NAME, "f'(u) = u/(1-u)^3 >= 0 on [0,1)", sp.simplify(sp.diff(f, u) - u / (1 - u) ** 3) == 0, eq=eq_fG)
    check(NAME, "identity (4.9): f(1-(1-u)^3) = G(u) u^2",
          sp.simplify(f.subs(u, 1 - (1 - u) ** 3) - G * u**2) == 0, eq=eq_fG)
    # With a = 1/(1-u) > 1 on (0,1) and a' = a^2:  G = (a + a^2 + a^3)^2 / 2 has
    # G' = (a + a^2 + a^3)(1 + 2a + 3a^2) a^2, a product of positive factors.
    a_ = 1 / (1 - u)
    G_prime_factored = (a_ + a_**2 + a_**3) * (1 + 2 * a_ + 3 * a_**2) * a_**2
    check(NAME, "G'(u) = (a+a^2+a^3)(1+2a+3a^2) a^2 with a = 1/(1-u): positive, so G is increasing on (0,1)",
          sp.simplify(sp.diff(G, u) - G_prime_factored) == 0, eq=eq_fG)

    # Squaring lemma 4.1 on random exact weights: 1 - P' <= (1-P)^2/(2P^2) <= f(1-P)
    rng = random.Random(268)
    worst = Fraction(0)
    for _ in range(400):
        m = rng.randint(2, 7)
        raw = [Fraction(rng.randint(1, 50), 1) for _ in range(m)]
        tot = sum(raw)
        w = [r / tot for r in raw]
        P = sum(x * x for x in w)
        Pp = sum((x * x / P) ** 2 for x in w)
        lhs = 1 - Pp
        rhs = (1 - P) ** 2 / (2 * P**2)
        worst = max(worst, lhs - rhs)
        assert lhs <= rhs, "Lemma 4.1 violated"
    eq_sq = register(NAME, Equation(
        key="haldane_squaring_lemma", entry=ENTRY,
        title="Squaring probabilities",
        source_kind="paper", source=f"{H}, Lemma 4.1, eq. (4.2)",
        latex_str=(r"P=\sum_i w_i^2,\ \ P'=\sum_i\Bigl(\frac{w_i^2}{P}\Bigr)^{2}"
                   r"\ \Longrightarrow\ 1-P'\le f(1-P)"),
    ))
    check(NAME, "Lemma 4.1 inequality holds on 400 random exact rational weight vectors (max slack <= 0)",
          worst <= 0, eq=eq_sq)

    # Identities (4.5)-(4.6) from S(n,b) = Z_{2n}(b)/Z_n(b)^2 and T(n,b) = Z_n(2b)/Z_n(b)^2
    n, b = sp.symbols("n b", positive=True)
    Z = sp.Function("Z")
    S_ = lambda nn, bb: Z(2 * nn, bb) / Z(nn, bb) ** 2
    T_ = lambda nn, bb: Z(nn, 2 * bb) / Z(nn, bb) ** 2
    eq_pur = register(NAME, Equation(
        key="haldane_purities", entry=ENTRY,
        title="Spatial and physical purities and their coupling",
        source_kind="paper", source=f"{H}, eqs. (4.1), (4.5), (4.6)",
        latex_str=(r"S(n,b)=\frac{Z_{2n}(b)}{Z_n(b)^{2}},\quad T(n,b)=\frac{Z_n(2b)}{Z_n(b)^{2}};%%BR%%"
                   r"S(n,2b)=S(n,b)^{2}\,\frac{T(2n,b)}{T(n,b)^{2}},%%BR%%"
                   r"T(4n,b)=T(2n,b)^{2}\,\frac{S(2n,2b)}{S(2n,b)^{2}}"),
    ))
    check(NAME, "identity (4.5): S(n,2b) = S(n,b)^2 T(2n,b) / T(n,b)^2",
          sp.simplify(S_(n, 2 * b) - S_(n, b) ** 2 * T_(2 * n, b) / T_(n, b) ** 2) == 0, eq=eq_pur)
    check(NAME, "identity (4.6): T(4n,b) = T(2n,b)^2 S(2n,2b) / S(2n,b)^2",
          sp.simplify(T_(4 * n, b) - T_(2 * n, b) ** 2 * S_(2 * n, 2 * b) / S_(2 * n, b) ** 2) == 0, eq=eq_pur)

    # Closed form of the defect sequence, and the Bernoulli-type bound used in (4.16)
    C, r, jj = sp.symbols("C r j", positive=True)
    uj = lambda k: r ** (2**k) / C
    eq_seq = register(NAME, Equation(
        key="haldane_defect_sequence", entry=ENTRY,
        title="Defect sequence of the bootstrap",
        source_kind="paper", source=f"{H}, Proposition 4.3, eqs. (4.10)-(4.14)",
        latex_str=(r"n_j=2^jn_0,\ \ b_j=2^jb_0,\ \ u_j=C^{-1}r^{2^j},\ \ r=Cu_0<1,\ \ "
                   r"u_{j+1}=Cu_j^{2};\qquad \gamma_L>\frac{\log(1/r)}{b_0}\ \ (L\ge n_0)"),
        notes="Hypotheses: 0<u_0<1/6, G(u_0)<=C, r<1, C(1-3u_0)>3, and the two initial purity bounds (4.11).",
    ))
    check(NAME, "u_{j+1} = C u_j^2 for u_j = r^(2^j)/C",
          sp.simplify(sp.powsimp(sp.expand_power_base(C * uj(jj) ** 2, force=True), force=True) - uj(jj + 1)) == 0,
          eq=eq_seq)
    check(NAME, "u_0 = r/C equals the starting defect u_0 when r = C u_0",
          sp.simplify(uj(0).subs(r, C * u) - u) == 0, eq=eq_seq)
    check(NAME, "(1-u)^3 - (1-3u) = u^2 (3-u) >= 0 on [0,1], giving (4.16)",
          sp.expand((1 - u) ** 3 - (1 - 3 * u) - u**2 * (3 - u)) == 0, eq=eq_seq)
    w_ = sp.symbols("w", positive=True)
    check(NAME, "(1-w)/w is strictly decreasing, so w_0 > 1-3u_j gives (1-w_0)/w_0 <= 3u_j/(1-3u_j)",
          sp.simplify(sp.diff((1 - w_) / w_, w_) + 1 / w_**2) == 0, eq=eq_seq)
    b0s, k2 = sp.symbols("b_0 k", positive=True)
    check(NAME, "e^{-b_j gamma} < r^(2^j) with b_j = 2^j b_0 gives gamma > log(1/r)/b_0 for every j",
          sp.simplify(sp.expand_log(-sp.log(r ** (2**k2)), force=True) / (2**k2 * b0s)
                      - (-sp.log(r)) / b0s) == 0, eq=eq_seq)

    # ---- the two parameter rows, exactly ----------------------------------------------------
    rows = [
        ("row 1 (n0=60, b0=105/4)", 60, sp.Rational(105, 4), sp.Rational(1, 8), sp.Rational(79, 10),
         sp.Rational(4, 105) * sp.log(sp.Rational(80, 79))),
        ("row 2 (n0=2304, b0=784)", 2304, sp.Integer(784), sp.Rational(1, 100), sp.Integer(5),
         sp.log(20) / 784),
    ]
    for label, n0, b0v, u0, Cv, gap in rows:
        rv = Cv * u0
        Gu0 = G.subs(u, u0)
        check(NAME, f"{label}: u0 < 1/6", u0 < sp.Rational(1, 6), eq=eq_seq)
        check(NAME, f"{label}: G(u0) = {Gu0} <= C = {Cv}", Gu0 <= Cv, eq=eq_seq)
        check(NAME, f"{label}: r = C u0 = {rv} < 1", rv < 1, eq=eq_seq)
        check(NAME, f"{label}: C(1-3u0) = {Cv * (1 - 3 * u0)} > 3", Cv * (1 - 3 * u0) > 3, eq=eq_seq)
        check(NAME, f"{label}: 3/(C(1-3u0)) = {3 / (Cv * (1 - 3 * u0))} < 1", 3 / (Cv * (1 - 3 * u0)) < 1, eq=eq_seq)
        check(NAME, f"{label}: log(1/r)/b0 equals the stated gap bound",
              sp.simplify(sp.log(1 / rv) / b0v - gap) == 0, eq=eq_thm)
    # Values the paper quotes explicitly in Section 4.3
    check(NAME, "G(1/8) = 913952/117649 (paper, Section 4.3)", G.subs(u, sp.Rational(1, 8)) == sp.Rational(913952, 117649), eq=eq_seq)
    check(NAME, "row 2: 9/(2 (1-3u)^2) at u=1/100 equals 45000/9409, and G(1/100) is below it",
          (sp.Rational(9, 2) / (1 - 3 * sp.Rational(1, 100)) ** 2 == sp.Rational(45000, 9409))
          and G.subs(u, sp.Rational(1, 100)) <= sp.Rational(45000, 9409), eq=eq_seq)
    check(NAME, "second bound is the stronger one: log(20)/784 > (4/105) log(80/79)",
          sp.N(b2) > sp.N(b1) > 0, eq=eq_thm)
    check(NAME, "dyadic intervals [2^j n0, 2^(j+1) n0] cover every length L >= n0 (both rows)",
          all(any(2**jx * n0 <= Lx <= 2 ** (jx + 1) * n0 for jx in range(0, 14))
              for n0 in (60, 2304) for Lx in range(n0, 40 * n0, 2)), eq=eq_seq)

    # ---- optional exact-diagonalization sanity check of the Hamiltonian -------------------
    try:
        import numpy as np
        import scipy.sparse as sps
        import scipy.sparse.linalg as spla
    except ImportError:
        return
    Sx, Sy, Sz = (np.array(sp.N(m, 20).tolist(), dtype=complex) for m in (sx, sy, sz))
    ops = [sps.csr_matrix(M) for M in (Sx, Sy, Sz)]
    I3 = sps.identity(3, format="csr")

    def site(op, j, Lc):
        out = sps.identity(1, format="csr")
        for i in range(Lc):
            out = sps.kron(out, op if i == j else I3, format="csr")
        return out

    results = {}
    for Lc in (4, 6, 8, 10):
        Hm = sum(site(o, j, Lc) @ site(o, (j + 1) % Lc, Lc) for j in range(Lc) for o in ops)
        if 3**Lc <= 729:
            ev = np.sort(np.linalg.eigvalsh(Hm.toarray()))
        else:
            ev = np.sort(spla.eigsh(Hm, k=8, which="SA", return_eigenvectors=False))
        tol = 1e-8
        e0 = ev[0]
        distinct = [e for e in ev if e - e0 > tol]
        g = distinct[0] - e0
        n_ground = int(sum(1 for e in ev if abs(e - e0) < tol))
        results[Lc] = (g, n_ground)
        check(NAME, f"exact diagonalization L={Lc}: ground state unique and gap gamma_L = {g:.4f} > log(20)/784",
              n_ground == 1 and g > float(sp.N(b2)), eq=eq_thm)
    eq_thm.notes += (" Sanity values from exact diagonalization (not a test of the theorem, whose explicit "
                     "bounds start at L = 60): " + "; ".join(f"L={Lc}: gap {g:.4f}" for Lc, (g, _) in results.items()) + ".")
