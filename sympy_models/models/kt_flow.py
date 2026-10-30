"""Critical and near-critical planar XY model: the KT flow in the papers (entry 216).

Statements are transcribed from manuscripts in this repository:

  [Log]  The critical logarithmic correction for the planar XY model
  [Ess]  Essential Singularity of the Correlation Length in the Planar XY Model
  [Ht]   BKT universality for height and planar spin fields
  [Ctr]  Critical Center Magnetization in the Planar XY Model   (conditional on companion inputs)
  [Spin] The Critical Spin Field of the Planar XY Model

all under preprints/ with the October 5 or September 24, 2026 suffixes.

The central object is the quadratic recursion (x, y) -> (x - y^2, y - x y) of [Ess] eq. (1.4)
and its rescaled limiting flow, eq. (1.6). The checks verify that the paper's equations
(1.4)-(1.8) are consistent with one another: the critical profile (1,1)/j, the growth
direction (1,-1/2), the small-w expansion of the flow solution, the pole T = pi/k and the
amplitude A_XY = (log L) pi / sqrt(3 P). They do not re-prove the theorems, and they do not
compute the non-universal constants A_XY, B_XY or P.

Reading note: in [Ess] the scaling variable is delta = sqrt(b_c - b); the radical is set
detached from its argument in the PDF text layer, so this was confirmed on the rendered page.
"""
import sympy as sp

from eqlib import Equation, check, register

NAME = "kt_flow"
ENTRY = "216"
LOG = "[Log] critical logarithmic correction paper"
ESS = "[Ess] essential singularity paper"
HT = "[Ht] BKT universality paper"
CTR = "[Ctr] critical center magnetization paper"
SPIN = "[Spin] critical spin field paper"


def build_and_verify():
    # ---- the theorems about the lattice XY model (statements only) ------------------------
    register(NAME, Equation(
        key="xy_definitions", entry=ENTRY,
        title="Free-box correlation, mass and critical inverse temperature",
        source_kind="paper", source=f"{LOG}, eq. (2); {ESS}, eq. (1.2)",
        latex_str=(r"C_b(r)=\lim_{R\to\infty}\mathbb{E}^{\mathrm{free}}_{R,b}\cos(\theta_0-\theta_{re_1}),%%BR%%"
                   r"m(b)=\lim_{r\to\infty}-\tfrac1r\log C_b(r),\quad b_c=\inf\{b>0:m(b)=0\}"),
        notes="Nearest-neighbour cosine XY model on Z^2 with free boundary; the thermodynamic limit "
              "precedes r -> infinity. xi(b) = 1/m(b) on the massive side.",
    ))
    eq_log = register(NAME, Equation(
        key="xy_critical_log_correction", entry=ENTRY,
        title="Critical correlation with logarithmic correction",
        source_kind="paper", source=f"{LOG}, Theorem 1.1",
        latex_str=(r"C_{b_c}(r)=B_{\mathrm{XY}}\,r^{-1/4}(\log r)^{1/8}\bigl(1+o(1)\bigr),"
                   r"\qquad B_{\mathrm{XY}}\in(0,\infty)"),
        notes="B_XY is a non-universal constant that is not evaluated.",
    ))
    eq_ess = register(NAME, Equation(
        key="xy_essential_singularity", entry=ENTRY,
        title="BKT essential singularity of the correlation length",
        source_kind="paper", source=f"{ESS}, Theorem 1.1, eq. (1.3)",
        latex_str=r"\lim_{b\uparrow b_c}\sqrt{b_c-b}\,\log\xi(b)=A_{\mathrm{XY}},\qquad A_{\mathrm{XY}}\in(0,\infty)",
        notes="The constant depends on the normalization of b and is not evaluated in closed form.",
    ))

    # ---- the recursion and its critical profile -------------------------------------------
    xs, ys, j = sp.symbols("x y j", positive=True)
    eq_rec = register(NAME, Equation(
        key="kt_recursion", entry=ENTRY,
        title="Leading recursion of the scalar coordinates",
        source_kind="paper", source=f"{ESS}, eq. (1.4)",
        expr=sp.Eq(sp.Symbol(r"(x,y)\mapsto"), sp.Symbol(r"(x-y^{2},\;y-xy)")),
        latex_str=r"(x,y)\longmapsto\bigl(x-y^{2},\;y-xy\bigr)",
        notes="x is a quadratic gradient perturbation, y the first Fourier harmonic of the coarse height weight.",
    ))
    eq_crit = register(NAME, Equation(
        key="kt_critical_profile", entry=ENTRY,
        title="Critical trajectory",
        source_kind="paper", source=f"{ESS}, text after eq. (1.4)",
        latex_str=r"(x_j,y_j)=\frac{(1,1)}{j}+O\!\left(\frac{\log j}{j^{2}}\right)\qquad(b=b_c)",
    ))
    u = sp.symbols("u", positive=True)  # u = 1/j
    step = u - u**2                       # x - y^2 with x = y = u
    series_gap = sp.series(u / (1 + u) - step, u, 0, 5).removeO()   # (1/(j+1)) - image of 1/j
    check(NAME, "(1,1)/j is a fixed profile of the recursion: its image differs from (1,1)/(j+1) only at O(j^-3)",
          sp.expand(series_gap).coeff(u, 1) == 0 and sp.expand(series_gap).coeff(u, 2) == 0
          and sp.expand(series_gap).coeff(u, 3) == 1, eq=eq_crit)

    # ---- linear response along the critical trajectory -------------------------------------
    jac = sp.Matrix([xs - ys**2, ys - xs * ys]).jacobian([xs, ys])
    at_crit = jac.subs({xs: 1 / j, ys: 1 / j})
    direction = sp.Matrix([1, -sp.Rational(1, 2)])
    eq_resp = register(NAME, Equation(
        key="kt_response", entry=ENTRY,
        title="Temperature derivative of the critical trajectory",
        source_kind="paper", source=f"{ESS}, eq. (1.5)",
        latex_str=r"\frac1j\,\partial_b(x_j,y_j)\Big|_{b=b_c}\longrightarrow\mathcal{P}\,\bigl(1,-\tfrac12\bigr),\qquad 0<\mathcal{P}<\infty",
        notes="Positivity of P is the essential input in the paper (a finite-volume variance derivative bound).",
    ))
    check(NAME, "linearization at (1/j,1/j) has (1,-1/2) as eigenvector with eigenvalue 1+1/j",
          sp.simplify(at_crit * direction - (1 + 1 / j) * direction) == sp.zeros(2, 1), eq=eq_resp)
    check(NAME, "eigenvalue product prod_{i=1}^{j-1}(1+1/i) = j for j = 2..30, so the response grows like j",
          all(sp.prod([1 + sp.Rational(1, i) for i in range(1, jj)]) == jj for jj in range(2, 31)), eq=eq_resp)
    check(NAME, "inductive step j (1+1/j) = j+1", sp.simplify(j * (1 + 1 / j) - (j + 1)) == 0, eq=eq_resp)

    # ---- limiting flow ---------------------------------------------------------------------
    w, k, P, L, delta = sp.symbols("w k P L delta", positive=True)
    X = k * sp.cot(k * w)
    Y = k * sp.csc(k * w)
    eq_flow = register(NAME, Equation(
        key="kt_limiting_flow", entry=ENTRY,
        title="Limiting rescaled flow and its explicit solution",
        source_kind="paper", source=f"{ESS}, eq. (1.6)",
        expr=sp.Eq(sp.Symbol(r"(X,Y)"), sp.Symbol(r"(k\cot kw,\;k\csc kw)")),
        latex_str=(r"X'=-Y^{2},\quad Y'=-XY,\qquad X(w)=k\cot(kw),\quad Y(w)=k\csc(kw),"
                   r"\quad k=\sqrt{3\mathcal{P}}"),
        notes="w = delta j with delta = sqrt(b_c - b); X = x/delta, Y = y/delta.",
    ))
    check(NAME, "X' = -Y^2", sp.simplify(sp.diff(X, w) + Y**2) == 0, eq=eq_flow)
    check(NAME, "Y' = -X Y", sp.simplify(sp.diff(Y, w) + X * Y) == 0, eq=eq_flow)
    check(NAME, "X^2 - Y^2 = -k^2 is conserved", sp.simplify(X**2 - Y**2 + k**2) == 0, eq=eq_flow)

    # Small-w expansion against the critical profile and its first-order response.
    # x_j = 1/j - delta^2 P j,  y_j = 1/j + delta^2 P j / 2   (b - b_c = -delta^2), with j = w/delta.
    x_j = 1 / j - delta**2 * P * j
    y_j = 1 / j + delta**2 * P * j / 2
    Xj = sp.simplify((x_j / delta).subs(j, w / delta))
    Yj = sp.simplify((y_j / delta).subs(j, w / delta))
    kP = sp.sqrt(3 * P)
    X_ser = sp.series(X.subs(k, kP), w, 0, 3).removeO()
    Y_ser = sp.series(Y.subs(k, kP), w, 0, 3).removeO()
    check(NAME, "critical profile + response (1.4)-(1.5) reproduce the small-w expansion of X = k cot kw",
          sp.simplify(X_ser - Xj) == 0, eq=eq_flow)
    check(NAME, "critical profile + response (1.4)-(1.5) reproduce the small-w expansion of Y = k csc kw",
          sp.simplify(Y_ser - Yj) == 0, eq=eq_flow)

    # ---- pole, exit profile and the amplitude --------------------------------------------
    T = sp.pi / k
    eq_pole = register(NAME, Equation(
        key="kt_pole_time", entry=ENTRY,
        title="Outgoing pole of the flow",
        source_kind="paper", source=f"{ESS}, text after eq. (1.6)",
        expr=sp.Eq(sp.Symbol("T"), T),
        latex_str=r"T=\frac{\pi}{k}\quad\text{and}\quad \delta\,j_e\to T",
    ))
    check(NAME, "Y -> +infinity as w -> T from below", sp.limit(Y, w, T, "-") == sp.oo, eq=eq_pole)
    check(NAME, "X -> -infinity as w -> T from below", sp.limit(X, w, T, "-") == -sp.oo, eq=eq_pole)
    check(NAME, "sin(kT) = 0 first at T = pi/k", sp.simplify(sp.sin(k * T)) == 0, eq=eq_pole)

    register(NAME, Equation(
        key="kt_height_profile", entry=ENTRY,
        title="Limiting height-variance coefficient profile",
        source_kind="paper", source=f"{ESS}, eq. (1.7)",
        latex_str=r"\mathcal{A}(w)=8\pi\ \ (0<w<T),\qquad \mathcal{A}(w)=0\ \ (w>T)",
    ))
    P_pos = sp.symbols("P_pos", positive=True)
    A_xy = sp.log(L) * T.subs(k, sp.sqrt(3 * P_pos))
    eq_A = register(NAME, Equation(
        key="kt_amplitude", entry=ENTRY,
        title="The constant A_XY in terms of the block ratio and P",
        source_kind="paper", source=f"{ESS}, eq. (1.8)",
        expr=sp.Eq(sp.Symbol(r"A_{\mathrm{XY}}"), A_xy),
        latex_str=r"A_{\mathrm{XY}}=(\log L)\,\frac{\pi}{\sqrt{3\mathcal{P}}}",
        notes="L is the fixed block ratio of the map; the paper notes the intermediate coordinates "
              "depend on the construction while the left side does not.",
    ))
    # Derivation: xi ~ L^{j_e} with delta j_e -> T and delta = sqrt(b_c - b)
    #   => sqrt(b_c - b) log xi = delta j_e log L -> T log L = (pi/k) log L.
    je = T.subs(k, sp.sqrt(3 * P)) / delta     # delta * j_e -> T, with k = sqrt(3 P)
    check(NAME, "sqrt(b_c-b) log xi = delta * j_e * log L -> (pi/sqrt(3P)) log L  (equation (1.8))",
          sp.simplify(delta * je * sp.log(L) - A_xy.subs(P_pos, P)) == 0, eq=eq_A)

    # ---- the BKT height threshold --------------------------------------------------------
    Jnn = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    v2_nn = sp.Rational(1, 2 * len(Jnn)) * sum(y[0] ** 2 for y in Jnn)
    eq_ht = register(NAME, Equation(
        key="kt_height_threshold", entry=ENTRY,
        title="Height model: variance of the Laplacian and critical effective temperature",
        source_kind="paper", source=f"{HT}, eqs. (1.1) and (1.8), Theorem 1.1",
        latex_str=(r"v_J^2=\frac{1}{2|J|}\sum_{y\in J}y_1^{2},%%BR%%"
                   r"\beta_{\mathrm{eff}}\bigl(J,\beta_c(J)\bigr)=8\pi v_J^{2}"),
        notes="J is a finite square-symmetric interaction set containing the nearest neighbours.",
    ))
    check(NAME, "v^2 for the nearest-neighbour set is 1/4", v2_nn == sp.Rational(1, 4), eq=eq_ht)
    check(NAME, "8 pi v^2 = 2 pi at nearest neighbours (stated in the paper)",
          sp.simplify(8 * sp.pi * v2_nn - 2 * sp.pi) == 0, eq=eq_ht)

    # ---- spin-field exponents ------------------------------------------------------------
    K, r = sp.symbols("K r", positive=True)
    cov = sp.log(1 / r) / (2 * sp.pi)    # local covariance singularity of the Dirichlet GFF, [Ht] after (1.14)
    two_point = sp.exp(cov / K)          # E[V(z) conj V(w)] = exp(Cov/K) from the normalization in (1.15)
    exponent = sp.simplify(-sp.log(two_point) / sp.log(r))
    eq_exp = register(NAME, Equation(
        key="kt_imaginary_chaos_exponent", entry=ENTRY,
        title="Two-point decay of the imaginary exponential of a GFF",
        source_kind="derived", source=f"{HT}, eqs. (1.15), Theorem 1.2; {SPIN}, abstract",
        latex_str=r"\mathbb{E}\bigl[V_K(z)\,\overline{V_K(w)}\bigr]\sim|z-w|^{-1/(2\pi K)}",
        notes="Uses the local covariance (2 pi)^{-1} log|z-w|^{-1} and the exp(E Phi^2 /(2K)) normalization. "
              "Theorem 1.2 has K_M(b)/b -> 1, the low-temperature spin-wave exponent 1/(2 pi b). "
              "[Spin] gives stiffness 2/pi at criticality.",
    ))
    check(NAME, "decay exponent is 1/(2 pi K)", sp.simplify(exponent - 1 / (2 * sp.pi * K)) == 0, eq=eq_exp)
    check(NAME, "stiffness K = 2/pi gives the critical exponent 1/4",
          sp.simplify((1 / (2 * sp.pi * K)).subs(K, 2 / sp.pi) - sp.Rational(1, 4)) == 0, eq=eq_exp)

    eq_ctr = register(NAME, Equation(
        key="xy_center_magnetization", entry=ENTRY,
        title="Center magnetization at criticality (conditional)",
        source_kind="paper", source=f"{CTR}, abstract",
        latex_str=r"M_n=A_{\mathrm{XY}}'\,n^{-1/8}(\log n)^{1/16}\bigl(1+o(1)\bigr)",
        notes="CONDITIONAL: the paper assumes the stated critical-height, local-renormalization and "
              "spin-field inputs from its companions. The paper's own symbol for the constant is A_XY; "
              "primed here to avoid a clash with the constant in the essential-singularity theorem.",
    ))
    n_ = sp.symbols("n", positive=True)
    check(NAME, "square of n^(-1/8)(log n)^(1/16) is n^(-1/4)(log n)^(1/8), the exponents of the critical correlation",
          sp.simplify(sp.powsimp(sp.expand_power_base((n_ ** sp.Rational(-1, 8) * sp.log(n_) ** sp.Rational(1, 16)) ** 2, force=True), force=True)
                      - n_ ** sp.Rational(-1, 4) * sp.log(n_) ** sp.Rational(1, 8)) == 0, eq=eq_ctr)
    eq_ctr.notes += " The squaring check is exponent bookkeeping only; the papers do not derive one from the other this way."
