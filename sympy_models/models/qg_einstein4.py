"""Einstein four-manifolds and the two SU(2)s of a spin network (entry 348 of the OpenAI release; companion to qg_spin_network.py).

Why this entry. A spin network is built from SU(2) holonomies. A Riemannian four-manifold has a second, geometric SU(2): the
group Spin(4) of its frames is SU(2)_+ x SU(2)_-, the two-forms split as Lambda^2 = Lambda^+ (+) Lambda^-, and each half is
the adjoint (spin 1) of one SU(2). Entry 348 (OpenAI, "Positively curved Einstein four-manifolds", 23 Sep 2026, with the L^2
gap paper of 5 Oct 2026 and the zero-plane rigidity paper of 4 Oct 2026) is a statement about the curvature operator in exactly
these terms: for Ric = 3g it is Id + T on Lambda^+ and Id + U on Lambda^-, with T, U symmetric trace-free (the two Weyl
blocks, each a spin-2 object of its SU(2)), and the Einstein equation removes the block that mixes Lambda^+ with Lambda^-.
This module computes those blocks for the three model spaces, and checks every number the paper and the library use about them.

Registered, with their checks:

  qg_so4_casimirs      so(4) = su(2)_+ + su(2)_- from explicit generators: the Casimirs (3/4, 3/4) on the vector, (2, 0) and (0, 2) on
                       Lambda^+ and Lambda^-, and 6 = j(j+1) at j = 2 on the Weyl blocks Sym^2_0 Lambda^+ (5 dimensional); the
                       two copies commute; a simple two-plane has |sigma_+| = |sigma_-|.
  qg_einstein_blocks   the curvature operator of the round S^4, of Fubini-Study CP^2 and of S^2(1/sqrt3) x S^2(1/sqrt3), all rescaled to
                       Ric = 3g, computed from the metric in a chart (not from a formula for the curvature of the model): the block
                       Id + T on Lambda^+, Id + U on Lambda^-, the mixed block (zero exactly when the metric is Einstein, shown by the
                       unequal product S^2(sqrt(1/3)) x S^2(1)), the Hilbert-Schmidt norms, the minimal sectional curvature
                       1 - (x + y)/2 against brute force over planes, and the second moments of the two Weyl blocks.
  qg_einstein_topology the Chern-Gauss-Bonnet and signature integrals from these curvatures and the exact volumes: chi = 2, 3, 4 and
                       tau = 0, 1, 0, written through the second moments E v^2 and E b^2 and the volume, together with the
                       Hitchin-Thorpe numbers and the scale invariance of the L^2 energy of the trace-free Ricci tensor.

Checked: the curvature operator from second derivatives of the chart metric (exact, sympy); Ric = 3g; the block structure; the
topological numbers; the closed forms; the range of sectional curvature over a seeded random sample of planes.
Not tested: any claim of the three papers beyond these three spaces: the classification itself, the existence and value of the L^2
constant eps_0 (the paper does not give a number), half-conformal flatness, and any statement about a metric other than the three
model spaces and the one unbalanced product. Nothing here says that a spin network or the ladder has anything to do with an Einstein
metric: what is checked is that the group that the curvature of an Einstein four-manifold lives in is the group the links of a spin
network carry, and the numbers that go with it.
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register

NAME = "qg_einstein4"
ENTRY = "348"
PAPER = "OpenAI, Positively curved Einstein four-manifolds (preprints/Positively-curved-Einstein-four-manifolds-September-23-2026)"

X = sp.symbols("x0 x1 x2 x3", real=True)
ORIGIN = {x: 0 for x in X}
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
S2 = sp.sqrt(2)


# ---------------------------------------------------------------------------------------------------------------
# metrics in a chart centred at a point, and the curvature operator there
# ---------------------------------------------------------------------------------------------------------------
def metric_sphere4(K=1):
    r2 = sum(x ** 2 for x in X)
    return sp.eye(4) * 4 / (K * (1 + r2) ** 2)


def metric_product(K1, K2):
    a = X[0] ** 2 + X[1] ** 2
    b = X[2] ** 2 + X[3] ** 2
    return sp.diag(*([4 / (K1 * (1 + a) ** 2)] * 2 + [4 / (K2 * (1 + b) ** 2)] * 2))


def metric_cp2(scale=2):
    """Fubini-Study in the affine chart (z1, z2) = (x0 + i x1, x2 + i x3): g = delta/(1+r^2) - (u u + v v)/(1+r^2)^2 with
    u = sum(x dx + y dy), v = sum(x dy - y dx). Holomorphic sectional curvature 4; times `scale` is Ric = 6/scale g."""
    x0, y0, x1, y1 = X
    r2 = sum(x ** 2 for x in X)
    u = sp.Matrix([x0, y0, x1, y1])
    v = sp.Matrix([-y0, x0, -y1, x1])
    return scale * (sp.eye(4) / (1 + r2) - (u * u.T + v * v.T) / (1 + r2) ** 2)


def riemann_origin(g):
    """R_{ijkl} in the coordinate basis at the origin (R_{ijij} = sectional curvature), from the metric's first and second
    derivatives there. R_ijkl = 1/2 (d_j d_k g_il + d_i d_l g_jk - d_i d_k g_jl - d_j d_l g_ik) + g^{mn}(G_{m,jk} G_{n,il} - G_{m,jl} G_{n,ik})."""
    g0 = g.subs(ORIGIN)
    gi = g0.inv()
    d1 = [[[sp.diff(g[i, j], X[k]).subs(ORIGIN) for k in range(4)] for j in range(4)] for i in range(4)]       # d_k g_ij
    d2 = {}
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(k, 4):
                    d2[(i, j, k, l)] = d2[(i, j, l, k)] = sp.diff(g[i, j], X[k], X[l]).subs(ORIGIN)
    G = [[[sp.Rational(1, 2) * (d1[j][k][i] + d1[i][k][j] - d1[i][j][k]) for k in range(4)] for j in range(4)] for i in range(4)]  # G[m][i][j] = G_{m,ij}
    R = {}
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    val = sp.Rational(1, 2) * (d2[(i, l, j, k)] + d2[(j, k, i, l)] - d2[(j, l, i, k)] - d2[(i, k, j, l)])
                    for m in range(4):
                        for n in range(4):
                            if gi[m, n] != 0:
                                val += gi[m, n] * (G[m][j][k] * G[n][i][l] - G[m][j][l] * G[n][i][k])
                    R[(i, j, k, l)] = sp.simplify(val)
    return R, g0


def frame_curvature(g):
    """Curvature operator on the 6 bivectors e_a ^ e_b (a < b) in the orthonormal frame e_i = d_i / sqrt(g_ii), and the Ricci matrix."""
    R, g0 = riemann_origin(g)
    assert all(g0[i, j] == 0 for i in range(4) for j in range(4) if i != j)
    s = [sp.sqrt(g0[i, i]) for i in range(4)]
    op = sp.Matrix(6, 6, lambda A, B: sp.simplify(R[(PAIRS[A][0], PAIRS[A][1], PAIRS[B][0], PAIRS[B][1])]
                                                   / (s[PAIRS[A][0]] * s[PAIRS[A][1]] * s[PAIRS[B][0]] * s[PAIRS[B][1]])))
    ric = sp.Matrix(4, 4, lambda j, k: sp.simplify(sum(R[(i, j, i, k)] / s[i] ** 2 for i in range(4)) / (s[j] * s[k])))
    return op, ric


def hodge_basis():
    """Orthonormal bases of Lambda^+ and Lambda^- in the basis e_a ^ e_b of the 6 bivectors, orientation e0 e1 e2 e3.
    *(01) = (23), *(02) = -(13), *(03) = (12)."""
    def vec(*entries):
        v = sp.zeros(6, 1)
        for pair, c in entries:
            v[PAIRS.index(pair)] = c
        return v / S2
    plus = [vec(((0, 1), 1), ((2, 3), 1)), vec(((0, 2), 1), ((1, 3), -1)), vec(((0, 3), 1), ((1, 2), 1))]
    minus = [vec(((0, 1), 1), ((2, 3), -1)), vec(((0, 2), 1), ((1, 3), 1)), vec(((0, 3), 1), ((1, 2), -1))]
    return plus, minus


def blocks(op):
    P, M = hodge_basis()
    Bp = sp.Matrix(3, 3, lambda a, b: sp.simplify((P[a].T * op * P[b])[0]))
    Bm = sp.Matrix(3, 3, lambda a, b: sp.simplify((M[a].T * op * M[b])[0]))
    C = sp.Matrix(3, 3, lambda a, b: sp.simplify((P[a].T * op * M[b])[0]))
    return Bp, Bm, C


def volume_chart(g):
    """Volume of a space whose chart is radially symmetric about the origin (S^4, CP^2): integral of sqrt(det g) over R^4 in polar
    form, with the density read off on the ray (r, 0, 0, 0)."""
    r = sp.symbols("r", positive=True)
    point = {X[0]: r, X[1]: 0, X[2]: 0, X[3]: 0}
    det = sp.factor(sp.cancel(g.subs(point).det()))
    dens = sp.sqrt(det)
    dens = sp.simplify(dens) if dens.is_Pow or dens.is_Mul else dens
    return sp.simplify(sp.integrate(2 * sp.pi ** 2 * r ** 3 * dens, (r, 0, sp.oo)))


def volume_sphere2(K):
    """Area of the sphere of curvature K from its stereographic chart: integral of 4/(K (1 + r^2)^2) over R^2 = 4 pi / K."""
    r = sp.symbols("r", positive=True)
    return sp.integrate(2 * sp.pi * r * 4 / (K * (1 + r ** 2) ** 2), (r, 0, sp.oo))


def hs(A):
    return sp.simplify(sp.sqrt((A * A.T).trace()))


# ---------------------------------------------------------------------------------------------------------------
# so(4) = su(2)_+ + su(2)_-
# ---------------------------------------------------------------------------------------------------------------
def so4_generators():
    def L(a, b):
        m = sp.zeros(4)
        m[a, b], m[b, a] = 1, -1
        return m
    eps = sp.LeviCivita
    Kp, Km = [], []
    for a in range(3):
        rot = sum((sp.Rational(1, 2) * eps(a, b, c) * L(b, c) for b in range(3) for c in range(3)), sp.zeros(4))
        Kp.append((rot + L(a, 3)) / 2)
        Km.append((rot - L(a, 3)) / 2)
    return Kp, Km


def comm(a, b):
    return a * b - b * a


def casimir(gens, rep):
    """-sum_a rep(K_a)^2 (the generators are real antisymmetric, so this is the positive Casimir j(j+1))."""
    return sum((-(rep(k) * rep(k)) for k in gens), sp.zeros(*rep(gens[0]).shape))


# ---------------------------------------------------------------------------------------------------------------
def build_and_verify():
    # =====================================================================================================
    # so(4) and the Casimirs
    # =====================================================================================================
    Kp, Km = so4_generators()
    eps = sp.LeviCivita
    struct_p = all(comm(Kp[a], Kp[b]) == sum((-eps(a, b, c) * Kp[c] for c in range(3)), sp.zeros(4)) for a in range(3) for b in range(3))
    struct_m = all(comm(Km[a], Km[b]) == sum((-eps(a, b, c) * Km[c] for c in range(3)), sp.zeros(4)) for a in range(3) for b in range(3))
    commute = all(comm(Kp[a], Km[b]) == sp.zeros(4) for a in range(3) for b in range(3))
    vec = lambda k: k
    cp_vec = casimir(Kp, vec)
    cm_vec = casimir(Km, vec)

    # Lambda^2 = antisymmetric 4x4 matrices: the adjoint action B -> [K, B] on the 6-dim space
    def basis_antisym():
        out = []
        for a, b in PAIRS:
            m = sp.zeros(4)
            m[a, b], m[b, a] = 1, -1
            out.append(m)
        return out

    AB = basis_antisym()

    def coords(m):                      # coordinates of an antisymmetric matrix in AB
        return sp.Matrix([m[a, b] for a, b in PAIRS])

    def ad(k):
        return sp.Matrix.hstack(*[coords(comm(k, B)) for B in AB])

    cp_adj = casimir(Kp, ad)
    cm_adj = casimir(Km, ad)

    # Sym^2 Lambda^+: the adjoint of su(2)_+ acts on Lambda^+ (3-dim) through its structure constants; on symmetric 3x3 T by T -> A T + T A^T
    P, _ = hodge_basis()
    Pm = sp.Matrix.hstack(*P)                     # 6 x 3 isometric embedding of Lambda^+
    A_p = [sp.simplify(Pm.T * ad(k) * Pm) for k in Kp]
    sym_basis = []
    for i in range(3):
        for j in range(i, 3):
            m = sp.zeros(3)
            m[i, j] = m[j, i] = 1
            sym_basis.append(m)

    def sym_coords(m):
        return sp.Matrix([m[i, j] for i in range(3) for j in range(i, 3)])

    def act_sym(a):
        return sp.Matrix.hstack(*[sym_coords(a * B + B * a.T) for B in sym_basis])

    cas_sym = sum((-(act_sym(a) * act_sym(a)) for a in A_p), sp.zeros(6))
    ev_sym = {k: v for k, v in sp.Matrix(cas_sym).eigenvals().items()}

    # the trace-free subspace is spin 2: the identity is the singlet, so the eigenvalue 6 has multiplicity 5, the eigenvalue 0 multiplicity 1
    eq_so4 = register(NAME, Equation(
        key="qg_so4_casimirs", entry=ENTRY, title="so(4) = su(2) + su(2): the two SU(2)s and the Casimirs of the curvature blocks",
        source_kind="derived", source="derived here from explicit generators; the decomposition is textbook (Singer-Thorpe 1969)",
        latex_str=(r"\Lambda^{2}=\Lambda^{+}\oplus\Lambda^{-}=(1,0)\oplus(0,1),\qquad (c_2^{+},c_2^{-})\Big|_{V}=\bigl(\tfrac34,\tfrac34\bigr),\quad "
                   r"\Big|_{\Lambda^{+}}=(2,0),\quad \Big|_{\Lambda^{-}}=(0,2)"
                   r"%%BR%%\mathrm{Sym}^{2}_{0}\Lambda^{\pm}=(2,0)\oplus(0,2)\ \ (\dim 5)\ \text{each},\qquad c_2=j(j+1)\big|_{j=2}=6"
                   r"%%BR%%\Sigma=\sigma_{+}+\sigma_{-}\ \text{simple}\ \Longrightarrow\ |\sigma_{+}|=|\sigma_{-}|"),
        notes="The frame group of an oriented Riemannian four-manifold is Spin(4) = SU(2) x SU(2); the two factors act on the self-dual and "
              "anti-self-dual two-forms, each as the adjoint. The Weyl curvature of one chirality is a trace-free symmetric endomorphism of "
              "Lambda^+ (or Lambda^-), the five-dimensional spin-2 representation of its SU(2), Casimir 6. This is the group of the links "
              "of a spin network with one factor dropped: the self-dual (Ashtekar) connection is the SU(2)_+ part of the Spin(4) connection. "
              "The last line is the algebraic form of the simplicity of a two-plane: the self-dual and anti-self-dual parts of a decomposable "
              "bivector have equal norm. Checked on explicit matrices; nothing here is a statement about gravity.",
    ))
    check(NAME, "the six generators K^+_a = (rot_a + L_a4)/2 and K^-_a = (rot_a - L_a4)/2 of so(4) satisfy [K_a, K_b] = -eps_abc K_c within each "
                "set and the two sets commute, as exact 4 x 4 matrices (so(4) = su(2) + su(2))", struct_p and struct_m and commute, eq=eq_so4)
    check(NAME, "on the vector representation the Casimirs are (c_+, c_-) = (3/4, 3/4): both matrices are (3/4) times the identity",
          cp_vec == sp.Rational(3, 4) * sp.eye(4) and cm_vec == sp.Rational(3, 4) * sp.eye(4), eq=eq_so4)
    check(NAME, "on the 6-dimensional space of two-forms the Casimirs have eigenvalues {2, 2, 2, 0, 0, 0} for su(2)_+ and for su(2)_- (the three "
                "that are 2 are exchanged by the Hodge star: Lambda^+ = (1,0), Lambda^- = (0,1)), j(j+1) = 2 at j = 1",
          sorted(cp_adj.eigenvals().items()) == [(0, 3), (2, 3)] and sorted(cm_adj.eigenvals().items()) == [(0, 3), (2, 3)]
          and sp.simplify(Pm.T * cp_adj * Pm - 2 * sp.eye(3)) == sp.zeros(3), eq=eq_so4)
    check(NAME, "the action of su(2)_+ on symmetric 3 x 3 matrices (the curvature blocks, 6-dimensional) has Casimir eigenvalue 6 = 2 (2 + 1) with "
                "multiplicity 5 (the trace-free part, spin 2) and 0 with multiplicity 1 (the trace, spin 0)",
          ev_sym == {6: 5, 0: 1}, eq=eq_so4)
    rng = np.random.default_rng(348)
    worst = 0.0
    Pn = np.array(Pm.evalf(), dtype=float)
    Mm = np.array(sp.Matrix.hstack(*hodge_basis()[1]).evalf(), dtype=float)
    for _ in range(200):
        u, v = np.linalg.qr(rng.normal(size=(4, 2)))[0].T
        biv = np.array([u[a] * v[b] - u[b] * v[a] for a, b in PAIRS])
        sp_ = Pn.T @ biv
        sm_ = Mm.T @ biv
        worst = max(worst, abs(sp_ @ sp_ - sm_ @ sm_), abs(sp_ @ sp_ - 0.5))
    check(NAME, f"for 200 random orthonormal pairs (u, v) in R^4 (seed 348) the bivector u ^ v has |sigma_+|^2 = |sigma_-|^2 = 1/2 (largest deviation "
                f"{worst:.1e}): a decomposable two-plane has equal self-dual and anti-self-dual parts", worst < 1e-12, eq=eq_so4)

    # =====================================================================================================
    # the three model spaces, from their metrics
    # =====================================================================================================
    models = {
        "S4": metric_sphere4(1),                        # sectional curvature 1, Ric = 3g
        "CP2": metric_cp2(2),                           # holomorphic sectional curvature 2, Ric = 3g
        "S2xS2": metric_product(3, 3),                  # two spheres of curvature 3 (radius 1/sqrt3), Ric = 3g
        "S2xS2unbal": metric_product(3, 1),             # not Einstein: Ric = diag(3, 3, 1, 1)
    }
    data = {}
    for key, g in models.items():
        op, ric = frame_curvature(g)
        Bp, Bm, C = blocks(op)
        data[key] = dict(op=op, ric=ric, Bp=Bp, Bm=Bm, C=C)

    three = ("S4", "CP2", "S2xS2")
    ric_ok = all(data[k]["ric"] == 3 * sp.eye(4) for k in three)
    spectra = {k: (sorted(data[k]["Bp"].eigenvals().items()), sorted(data[k]["Bm"].eigenvals().items())) for k in data}
    T = {k: data[k]["Bp"] - sp.eye(3) for k in three}
    U = {k: data[k]["Bm"] - sp.eye(3) for k in three}

    eq_blk = register(NAME, Equation(
        key="qg_einstein_blocks", entry=ENTRY, title="Curvature blocks of the three model Einstein four-manifolds",
        source_kind="derived", source=f"derived here from the metrics in a chart; the block form is the one of {PAPER}, section 1",
        latex_str=(r"\mathcal R=\begin{pmatrix}\mathrm{Id}+T&0\\0&\mathrm{Id}+U\end{pmatrix},\qquad \mathrm{Ric}=3g"
                   r"%%BR%%S^{4}:\ (T,U)=(0,0);\qquad \mathbb{CP}^{2}:\ T\sim\mathrm{diag}(2,-1,-1),\ U=0;"
                   r"\qquad S^{2}\!\times\! S^{2}:\ T\sim U\sim\mathrm{diag}(2,-1,-1)"
                   r"%%BR%%|T|^{2}=0,\ 6,\ 6;\quad |U|^{2}=0,\ 0,\ 6;\qquad K_{\min}=1-\tfrac{x+y}{2}=1,\ \tfrac12,\ 0,\qquad "
                   r"x=-\lambda_{\min}(T),\ y=-\lambda_{\min}(U)"),
        notes="For each space the metric is written in a chart centred at a point, the curvature operator on the six bivectors is computed "
              "from the first and second derivatives of the metric there (exact, sympy), and it is split into its action on the self-dual and "
              "anti-self-dual two-forms. The metrics are rescaled to Ric = 3g, the normalisation of the paper. All three spaces are "
              "homogeneous, so the same blocks hold at every point. The trace-free parts T, U are the Weyl blocks. The unbalanced product "
              "S^2(K=3) x S^2(K=1) has Ric = diag(3,3,1,1), not Einstein, and a nonzero block between Lambda^+ and Lambda^-. The minimal "
              "sectional curvature 1 - (x+y)/2 is the condition 'x + y < 2' of the paper for strict positivity: it is violated, with equality, "
              "by S^2 x S^2, which has planes of zero curvature. Nothing is checked about metrics other than these four.",
    ))
    check(NAME, "the metrics, in their charts, give Ric = 3g for the round S^4, Fubini-Study CP^2 (holomorphic sectional curvature 2) and "
                "S^2(1/sqrt3) x S^2(1/sqrt3), computed from second derivatives of the metric (exact)", ric_ok, eq=eq_blk)
    check(NAME, "the mixed block <Lambda^+, R Lambda^-> vanishes for the three Einstein spaces, and is nonzero for the unbalanced product "
                "S^2(K=3) x S^2(K=1), whose Ricci tensor diag(3,3,1,1) is not a multiple of g (so the vanishing of the mixed block is the "
                "trace-free Ricci equation, and it fails when that does)",
          all(data[k]["C"] == sp.zeros(3) for k in three) and data["S2xS2unbal"]["C"] != sp.zeros(3)
          and data["S2xS2unbal"]["ric"] == sp.diag(3, 3, 1, 1), eq=eq_blk)
    check(NAME, "the blocks have the spectra {1,1,1 | 1,1,1} (S^4), {3,0,0 | 1,1,1} (CP^2, self-dual orientation), {3,0,0 | 3,0,0} (S^2 x S^2): "
                "Id + T and Id + U with T, U trace-free",
          spectra["S4"] == ([(1, 3)], [(1, 3)]) and spectra["CP2"] == ([(0, 2), (3, 1)], [(1, 3)])
          and spectra["S2xS2"] == ([(0, 2), (3, 1)], [(0, 2), (3, 1)]) and all(M.trace() == 0 for M in (*T.values(), *U.values())),
          eq=eq_blk)
    hsT = {k: hs(T[k]) ** 2 for k in three}
    hsU = {k: hs(U[k]) ** 2 for k in three}
    check(NAME, "the squared Hilbert-Schmidt norms |T|^2 and |U|^2 are 0, 6, 6 and 0, 0, 6 for S^4, CP^2, S^2 x S^2; with v = |T|/sqrt6, b = |U|/sqrt6 "
                "this is v^2 = 1 and b^2 = 1 on the nonzero blocks: the gap E v^2 >= 1 of the paper (Gursky-LeBrun) holds with equality on "
                "each nonzero block of the three spaces",
          [hsT[k] for k in three] == [0, 6, 6] and [hsU[k] for k in three] == [0, 0, 6], eq=eq_blk)
    lam = lambda M: min(sp.Matrix(M).eigenvals())
    lam_max = lambda M: max(sp.Matrix(M).eigenvals())
    xy = {k: (-lam(T[k]), -lam(U[k])) for k in three}
    vb = {k: (sp.sqrt(hsT[k] / 6), sp.sqrt(hsU[k] / 6)) for k in three}
    check(NAME, "x = -lambda_min(T), y = -lambda_min(U) are (0,0), (1,0), (1,1); x + y = 0, 1, 2, so x + y < 2 (strict sectional positivity) for S^4 and "
                "CP^2 and not for S^2 x S^2; v <= x <= 2v holds on every nonzero block (v = 1, x = 1)",
          [xy[k] for k in three] == [(0, 0), (1, 0), (1, 1)]
          and all((vb[k][0] == 0 and xy[k][0] == 0) or (vb[k][0] <= xy[k][0] <= 2 * vb[k][0]) for k in three)
          and all((vb[k][1] == 0 and xy[k][1] == 0) or (vb[k][1] <= xy[k][1] <= 2 * vb[k][1]) for k in three), eq=eq_blk)

    # sectional curvature of every plane: K = (1/2) (n_+^T (1+T) n_+ + n_-^T (1+U) n_-) for unit n_pm; min over the two spheres
    kmin_formula = {k: 1 - sp.Rational(1, 2) * (xy[k][0] + xy[k][1]) for k in three}
    kmax_formula = {k: 1 + sp.Rational(1, 2) * (lam_max(T[k]) + lam_max(U[k])) for k in three}
    rng = np.random.default_rng(2026)
    sec = {}
    for k in three:
        op = np.array(data[k]["op"].evalf(), dtype=float)
        vals = []
        for _ in range(4000):
            u, v = np.linalg.qr(rng.normal(size=(4, 2)))[0].T
            biv = np.array([u[a] * v[b] - u[b] * v[a] for a, b in PAIRS])
            vals.append(biv @ op @ biv)
        sec[k] = (min(vals), max(vals))
    ok_sec = all(sec[k][0] >= float(kmin_formula[k]) - 1e-12 and sec[k][1] <= float(kmax_formula[k]) + 1e-12
                 and abs(sec[k][0] - float(kmin_formula[k])) < 3e-2 and abs(sec[k][1] - float(kmax_formula[k])) < 3e-2 for k in three)
    check(NAME, "the sectional curvature of 4000 random two-planes (seed 2026) lies inside [1 - (x+y)/2, 1 + (x'+y')/2] and comes within 3e-2 of both ends "
                "(the ends are reached only on a null set of planes): S^4 [1, 1], CP^2 [1/2, 2], S^2 x S^2 [0, 3] (measured "
                + "; ".join(f"{k} [{sec[k][0]:.3f}, {sec[k][1]:.3f}]" for k in three) + ")", ok_sec, eq=eq_blk)
    Pm_, Mm_ = (sp.Matrix.hstack(*b) for b in hodge_basis())
    exact_ok = True
    for k in three:
        for pick, formula in ((min, kmin_formula[k]), (max, kmax_formula[k])):
            evp = pick(sp.Matrix(data[k]["Bp"]).eigenvects(), key=lambda t: t[0])
            evm = pick(sp.Matrix(data[k]["Bm"]).eigenvects(), key=lambda t: t[0])
            npl = evp[2][0] / evp[2][0].norm()
            nmi = evm[2][0] / evm[2][0].norm()
            sig = (Pm_ * npl + Mm_ * nmi) / S2
            Kval = sp.simplify((sig.T * data[k]["op"] * sig)[0])
            simple = sp.simplify((Pm_.T * sig).dot(Pm_.T * sig) - (Mm_.T * sig).dot(Mm_.T * sig)) == 0
            exact_ok = exact_ok and simple and sp.simplify(Kval - formula) == 0
    check(NAME, "exactly: the plane built from the lowest (highest) eigenvectors n_+ of Id + T and n_- of Id + U, sigma = (n_+ + n_-)/sqrt2, has "
                "|sigma_+| = |sigma_-| and sectional curvature (1/2)(lambda_min(Id+T) + lambda_min(Id+U)) = 1 - (x+y)/2 (and with lambda_max the maximum) "
                "on S^4, CP^2 and S^2 x S^2", exact_ok, eq=eq_blk)

    # =====================================================================================================
    # topology: Chern-Gauss-Bonnet and the signature theorem on the three spaces
    # =====================================================================================================
    vols = {"S4": volume_chart(models["S4"]), "CP2": volume_chart(models["CP2"]), "S2xS2": volume_sphere2(3) ** 2}
    check(NAME, "volumes from the integral of sqrt(det g) over the chart (exact): 8 pi^2/3 (S^4), 2 pi^2 (CP^2, Ric = 3g), 16 pi^2/9 (S^2 x S^2)",
          vols["S4"] == 8 * sp.pi ** 2 / 3 and vols["CP2"] == 2 * sp.pi ** 2 and vols["S2xS2"] == 16 * sp.pi ** 2 / 9, eq=None)
    rm2 = {k: 4 * (data[k]["op"] * data[k]["op"]).trace() for k in three}              # R_ijkl R^ijkl = 4 tr(op^2)
    ric2 = {k: (data[k]["ric"] * data[k]["ric"]).trace() for k in three}
    scal = {k: data[k]["ric"].trace() for k in three}
    chi = {k: sp.simplify(vols[k] * (rm2[k] - 4 * ric2[k] + scal[k] ** 2) / (32 * sp.pi ** 2)) for k in three}
    wp2 = {k: 4 * hsT[k] for k in three}                                               # |W^+|^2 = 4 |T|^2 (tensor norm)
    wm2 = {k: 4 * hsU[k] for k in three}
    tau = {k: sp.simplify(vols[k] * (wp2[k] - wm2[k]) / (48 * sp.pi ** 2)) for k in three}
    s_p = {k: hsT[k] / 6 for k in three}                                               # E v^2 (constant on a homogeneous space)
    s_m = {k: hsU[k] / 6 for k in three}
    eq_top = register(NAME, Equation(
        key="qg_einstein_topology", entry=ENTRY, title="Euler characteristic and signature of the three model spaces from their curvature",
        source_kind="derived", source=f"derived here; the identities in this form are the ones the paper uses to tie the two Weyl second moments to the volume ({PAPER}, sections 3-4)",
        latex_str=(r"\chi=\frac{1}{32\pi^{2}}\int\bigl(|\mathrm{Rm}|^{2}-4|\mathrm{Ric}|^{2}+R^{2}\bigr)\,d\mu,\qquad "
                   r"\tau=\frac{1}{48\pi^{2}}\int\bigl(|W^{+}|^{2}-|W^{-}|^{2}\bigr)\,d\mu"
                   r"%%BR%%\mathrm{Ric}=3g:\quad \chi=\frac{3V}{4\pi^{2}}\bigl(1+s_{+}+s_{-}\bigr),\qquad \tau=\frac{V}{2\pi^{2}}\bigl(s_{+}-s_{-}\bigr),"
                   r"\qquad s_{+}=\mathbb E\,v^{2},\ \ s_{-}=\mathbb E\,b^{2}"
                   r"%%BR%%(\chi,\tau)=(2,0),\ (3,1),\ (4,0),\qquad V=\tfrac{8\pi^{2}}{3},\ 2\pi^{2},\ \tfrac{16\pi^{2}}{9}"),
        notes="With Ric = 3g the Ricci terms of the Gauss-Bonnet integrand cancel against the scalar curvature (R = 12, |Ric|^2 = 36, so "
              "-4 |Ric|^2 + R^2 = 0) and |Rm|^2 = |W|^2 + 24, which gives the first form; the second is the signature theorem. They are "
              "the 'characteristic identities' that relate the two second moments s_+ = E v^2 and s_- = E b^2 to the volume and to chi and tau. "
              "Here they are evaluated on the three homogeneous spaces only, where the integrands are constants and the volume comes from the "
              "chart metric. The orientation is the complex one on CP^2 (tau = +1; the opposite orientation gives -1 and exchanges the blocks). "
              "The Hitchin-Thorpe inequality 2 chi >= 3 |tau| holds with room on all three, and the Gursky-Malchiodi condition 2 chi - 3 |tau| <= 4 "
              "(the paper) holds for S^4 (4) and CP^2 (3) and not for S^2 x S^2 (8). The conclusions of the papers, which apply the identities "
              "to an arbitrary metric, are not tested.",
    ))
    check(NAME, "Chern-Gauss-Bonnet from the computed curvature: chi = 2, 3, 4 for S^4, CP^2, S^2 x S^2 (|Rm|^2 = 24, 48, 72; |Ric|^2 = 36; R = 12; "
                "volumes as above), exact",
          [chi[k] for k in three] == [2, 3, 4] and [rm2[k] for k in three] == [24, 48, 72]
          and all(ric2[k] == 36 and scal[k] == 12 for k in three), eq=eq_top)
    check(NAME, "the signature theorem from the computed Weyl blocks (|W^+|^2 = 4|T|^2, |W^-|^2 = 4|U|^2): tau = 0, 1, 0, exact; the self-dual block of "
                "CP^2 carries the whole signature",
          [tau[k] for k in three] == [0, 1, 0], eq=eq_top)
    s_pp, s_mm, V, chi_s, tau_s = sp.symbols("s_p s_m V chi tau")
    closed_chi = 3 * V / (4 * sp.pi ** 2) * (1 + s_pp + s_mm)
    closed_tau = V / (2 * sp.pi ** 2) * (s_pp - s_mm)
    ok_closed = all(sp.simplify(closed_chi.subs({V: vols[k], s_pp: s_p[k], s_mm: s_m[k]}) - chi[k]) == 0
                    and sp.simplify(closed_tau.subs({V: vols[k], s_pp: s_p[k], s_mm: s_m[k]}) - tau[k]) == 0 for k in three)
    check(NAME, "the closed forms chi = 3V(1 + s_+ + s_-)/(4 pi^2) and tau = V(s_+ - s_-)/(2 pi^2) reproduce chi and tau on the three spaces with "
                "(s_+, s_-) = (0,0), (1,0), (1,1) and the volumes above (exact)", ok_closed
          and [(s_p[k], s_m[k]) for k in three] == [(0, 0), (1, 0), (1, 1)], eq=eq_top)
    # invert: the Weyl integrals from chi, tau and V
    wp_from = 16 * sp.pi ** 2 * chi_s - 12 * V + 24 * sp.pi ** 2 * tau_s
    wm_from = 16 * sp.pi ** 2 * chi_s - 12 * V - 24 * sp.pi ** 2 * tau_s
    ok_inv = all(sp.simplify(wp_from.subs({chi_s: chi[k], tau_s: tau[k], V: vols[k]}) - vols[k] * wp2[k]) == 0
                 and sp.simplify(wm_from.subs({chi_s: chi[k], tau_s: tau[k], V: vols[k]}) - vols[k] * wm2[k]) == 0 for k in three)
    check(NAME, "inverting: the integrals of |W^+|^2 and |W^-|^2 over the space are 16 pi^2 chi - 12 V +- 24 pi^2 tau (exact on the three spaces); they "
                "are non-negative, which is the Hitchin-Thorpe-type restriction chi >= (3/(4 pi^2)) V + (3/2) |tau|", ok_inv, eq=eq_top)
    ht = {k: 2 * chi[k] - 3 * abs(tau[k]) for k in three}
    check(NAME, "2 chi - 3 |tau| = 4, 3, 8 for S^4, CP^2, S^2 x S^2, all >= 0 (Hitchin-Thorpe); only the first two satisfy 2 chi - 3 |tau| <= 4",
          [ht[k] for k in three] == [4, 3, 8], eq=eq_top)

    # L^2 energy of the trace-free Ricci tensor (the quantity of the L^2 gap paper): zero on Einstein metrics, scale invariant in dimension 4
    c, a = sp.symbols("c a", positive=True)
    E2_scaling = sp.simplify((c ** -2) * c ** 2)                                      # |E|^2_{cg} = c^-2 |E|^2_g, d mu_{cg} = c^2 d mu_g
    ric_u = data["S2xS2unbal"]["ric"]
    R_u = ric_u.trace()
    E_u = ric_u - R_u / 4 * sp.eye(4)
    vol_u = volume_sphere2(3) * volume_sphere2(1)
    D_u = sp.simplify((E_u * E_u).trace() * vol_u)
    check(NAME, "the energy D(g) = integral of |Ric - (R/4) g|^2 is scale invariant in dimension 4 (the integrand scales as c^-2 and the volume as c^2), "
                "vanishes on the three Einstein spaces (Ric = 3g gives the trace-free part 0), and is 64 pi^2/3 on S^2(K=3) x S^2(K=1) "
                "(trace-free Ricci diag(1,1,-1,-1), |E|^2 = 4, volume 16 pi^2/3)",
          E2_scaling == 1 and all((data[k]["ric"] - scal[k] / 4 * sp.eye(4)) == sp.zeros(4) for k in three)
          and E_u == sp.diag(1, 1, -1, -1) and D_u == 64 * sp.pi ** 2 / 3, eq=eq_top)
