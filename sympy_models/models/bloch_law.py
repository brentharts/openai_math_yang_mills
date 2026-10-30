"""Bloch's law, its first lattice correction, and the spherical magnetization law (entry 271).

Statements are transcribed from three manuscripts in this repository:

  [B]  Bloch's Law for Finite-Range Heisenberg Ferromagnets in Three Dimensions
       preprints/Blochs-Law-for-Finite-Range-Heisenberg-Ferromagnets-in-Three-Dimensions-October-5-2026/
  [C]  The first lattice correction to Bloch's law
       preprints/The-first-lattice-correction-to-Blochs-law-October-5-2026/
  [Sph] The spherical magnetization law for the three-dimensional quantum Heisenberg ferromagnet
       preprints/The-spherical-magnetization-law-for-the-three-dimensional-quantum-Heisenberg-ferromagnet-October-5-2026/

The checks do not re-prove those theorems (which concern the interacting quantum spin
system). They verify that the stated coefficients follow from the ideal-magnon integral
that [B] and [C] identify as the source of the leading terms, and that the identities
connecting the statements to each other hold.
"""
import mpmath as mp
import sympy as sp

from eqlib import Equation, check, register

NAME = "bloch_law"
ENTRY = "271"
B = "[B] Bloch's Law paper"
C = "[C] first lattice correction paper"
SPH = "[Sph] spherical magnetization paper"


def _D_matrix(support, S):
    """D_{S,J} = (S/2) sum_z J(z) z z^T for a symmetric finite support {z: J(z)}."""
    D = sp.zeros(3, 3)
    for z, w in support.items():
        v = sp.Matrix(z)
        D += sp.Rational(1, 2) * S * w * v * v.T
    return D


def _dispersion(support, S, k):
    """eps_{S,J}(k) = S sum_z J(z) (1 - cos(k . z))."""
    return S * sum(w * (1 - sp.cos(sum(ki * zi for ki, zi in zip(k, z)))) for z, w in support.items())


def build_and_verify():
    S, beta = sp.symbols("S beta", positive=True)
    k = sp.symbols("k1 k2 k3", real=True)
    lam = sp.symbols("lambda", positive=True)
    kvec = sp.Matrix(k)

    # ---- [B] model quantities ------------------------------------------------------
    eq_D = register(NAME, Equation(
        key="bloch_D", entry=ENTRY,
        title="Quadratic dispersion matrix",
        source_kind="paper", source=f"{B}, eq. (1.4)",
        latex_str=r"D_{S,J}=\frac{S}{2}\sum_{z\in\mathbb{Z}^3}J(z)\,zz^{T}",
        notes="J >= 0 finite-range, J(0)=0, J(z)=J(-z), support generates Z^3 (eq. 1.1).",
    ))
    eq_disp = register(NAME, Equation(
        key="bloch_dispersion", entry=ENTRY,
        title="One-magnon dispersion and its quadratic part",
        source_kind="paper", source=f"{B}, eq. (1.5)",
        latex_str=(r"\varepsilon_{S,J}(k)=S\sum_{z}J(z)\bigl(1-\cos(k\cdot z)\bigr)"
                   r"=k^{T}D_{S,J}\,k+O(|k|^{4})"),
    ))

    # Check eq. (1.5) on two supports: anisotropic nearest neighbour, and nn plus a diagonal bond.
    a1, a2, a3, c = sp.symbols("a1 a2 a3 c", positive=True)
    unit = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

    def sym(vec):  # both z and -z
        return [vec, tuple(-x for x in vec)]

    support_aniso = {}
    for w, e in zip((a1, a2, a3), unit):
        for z in sym(e):
            support_aniso[z] = w
    support_nnn = dict(support_aniso)
    for z in sym((1, 1, 0)):
        support_nnn[z] = c

    for label, sup in (("anisotropic nearest-neighbour", support_aniso),
                       ("nearest-neighbour plus (1,1,0) bond", support_nnn)):
        Dm = _D_matrix(sup, S)
        eps = _dispersion(sup, S, [lam * ki for ki in k])
        quad = sp.series(eps, lam, 0, 3).removeO().coeff(lam, 2)
        check(NAME, f"eq (1.5): quadratic part of eps equals k^T D k ({label})",
              sp.simplify(quad - (kvec.T * Dm * kvec)[0, 0]) == 0, eq=eq_disp)

    # Nearest-neighbour specialisation: D = S I, as stated after Theorem 1.1 of [B].
    support_nn = {z: 1 for e in unit for z in sym(e)}
    D_nn = _D_matrix(support_nn, S)
    check(NAME, "J(+-e_j)=1 gives D = S I", D_nn == S * sp.eye(3), eq=eq_D)

    # ---- Bloch coefficient from the ideal-magnon integral -------------------------------
    t_, d1, d2, d3 = sp.symbols("t d1 d2 d3", positive=True)
    kk = sp.symbols("kk", real=True)
    one_d = [sp.integrate(sp.exp(-t_ * d * kk**2) / (2 * sp.pi), (kk, -sp.oo, sp.oo)) for d in (d1, d2, d3)]
    gauss = sp.simplify(one_d[0] * one_d[1] * one_d[2])
    target = (4 * sp.pi * t_) ** sp.Rational(-3, 2) / sp.sqrt(d1 * d2 * d3)
    gauss_ok = sp.simplify(gauss - target) == 0

    zeta32 = sp.zeta(sp.Rational(3, 2))
    detD = sp.Symbol(r"\det D_{S,J}", positive=True)
    coeff_expr = zeta32 / (8 * sp.pi ** sp.Rational(3, 2) * sp.sqrt(detD))
    eq_bloch = register(NAME, Equation(
        key="bloch_theorem", entry=ENTRY,
        title="Bloch's T^{3/2} law with exact coefficient",
        source_kind="paper", source=f"{B}, Theorem 1.1, eq. (1.6)",
        expr=sp.Eq(sp.Symbol(r"c_{\mathrm{Bloch}}"), coeff_expr),
        latex_str=(r"\lim_{\beta\to\infty}\beta^{3/2}\bigl(S-m_{S,J}(\beta)\bigr)"
                   r"=\frac{\zeta(3/2)}{8\pi^{3/2}\sqrt{\det D_{S,J}}}"),
        notes="Order of limits is part of the statement: volume, then field h -> 0+, then beta -> infinity.",
    ))
    check(NAME, "Gaussian integral over R^3 of exp(-t k^T D k)/(2 pi)^3 is (4 pi t)^(-3/2)/sqrt(det D)",
          gauss_ok, eq=eq_bloch)
    # Summing e^{-n beta eps} over n >= 1 turns the Bose factor into zeta(3/2):
    #   sum_n (4 pi n beta)^(-3/2) / sqrt(det D) = zeta(3/2) / (8 pi^(3/2) beta^(3/2) sqrt(det D)).
    n = sp.symbols("n", positive=True)
    summand = (4 * sp.pi * n * beta) ** sp.Rational(-3, 2)
    prefactor = sp.simplify(summand * n ** sp.Rational(3, 2) * beta ** sp.Rational(3, 2))
    check(NAME, "Bose sum prefactor (4 pi)^(-3/2) equals 1/(8 pi^(3/2))",
          sp.simplify(prefactor - 1 / (8 * sp.pi ** sp.Rational(3, 2))) == 0, eq=eq_bloch)

    eq_nn = register(NAME, Equation(
        key="bloch_nearest_neighbour", entry=ENTRY,
        title="Nearest-neighbour Bloch coefficient",
        source_kind="paper", source=f"{B}, remark after Theorem 1.1",
        expr=sp.Eq(sp.Symbol(r"c_{\mathrm{Bloch}}^{\mathrm{nn}}"),
                   zeta32 / (8 * sp.pi ** sp.Rational(3, 2) * S ** sp.Rational(3, 2))),
    ))
    check(NAME, "sqrt(det(S I)) = S^(3/2) reproduces the nearest-neighbour coefficient",
          sp.simplify(sp.sqrt(D_nn.det()) - S ** sp.Rational(3, 2)) == 0, eq=eq_nn)

    # ---- [C] ideal-magnon density and the lattice correction ---------------------------
    eq_nsw = register(NAME, Equation(
        key="bloch_nsw", entry=ENTRY,
        title="Full ideal-magnon density (nearest neighbour)",
        source_kind="paper", source=f"{C}, eq. (1.1)",
        latex_str=(r"n^{\mathrm{sw}}_S(\beta)=\int_{[-\pi,\pi]^3}\frac{d^3k}{(2\pi)^3}"
                   r"\,\frac{1}{e^{\beta\varepsilon_S(k)}-1},%%BR%%"
                   r"\varepsilon_S(k)=2S\sum_{j=1}^{3}(1-\cos k_j)"),
    ))
    eq_cthm = register(NAME, Equation(
        key="bloch_correction_theorem", entry=ENTRY,
        title="Magnetization deficit matches the ideal-magnon density",
        source_kind="paper", source=f"{C}, Theorem 1.1",
        latex_str=r"\lim_{\beta\to\infty}(\beta S)^{5/2}\Bigl(S-m_S(\beta)-n^{\mathrm{sw}}_S(\beta)\Bigr)=0",
    ))

    # Heat-kernel step: int_{-pi}^{pi} dk/(2 pi) e^{-w(1 - cos k)} = e^{-w} I_0(w).
    # Large-w expansion of f(w) = e^{-w} I_0(w), from the Bessel ODE  z I0'' + I0' - z I0 = 0,
    # which for f reads  z f'' + (2z + 1) f' + f = 0.
    z, x = sp.symbols("z x", positive=True)
    a = sp.symbols("a1:4")
    g = sum(([1] + list(a))[j] * z ** (-sp.Rational(1, 2) - j) for j in range(4))
    ode = sp.expand((z * sp.diff(g, z, 2) + (2 * z + 1) * sp.diff(g, z) + g) * sp.sqrt(z))
    poly = sp.Poly(sp.expand(ode.subs(z, 1 / x)), x)
    sol = sp.solve([poly.coeff_monomial(x**j) for j in (1, 2, 3)], a, dict=True)[0]
    a1v, a2v, a3v = sol[a[0]], sol[a[1]], sol[a[2]]
    eq_bessel = register(NAME, Equation(
        key="bloch_bessel_asymptotics", entry=ENTRY,
        title="Large-argument expansion of e^{-z} I_0(z)",
        source_kind="derived", source="solve z f'' + (2z+1) f' + f = 0 by series",
        expr=sp.Eq(sp.Function("f")(z),
                   (2 * sp.pi * z) ** sp.Rational(-1, 2)
                   * (1 + a1v / z + a2v / z**2 + a3v / z**3 + sp.Symbol(r"O(z^{-4})") )),
        latex_str=(r"e^{-z}I_0(z)=\frac{1}{\sqrt{2\pi z}}\Bigl(1+\frac{1}{8z}+\frac{9}{128z^{2}}"
                   r"+\frac{75}{1024z^{3}}+O(z^{-4})\Bigr)"),
        notes="Used to expand the lattice heat kernel; standard asymptotics, derived here from the ODE.",
    ))
    check(NAME, "e^{-z} I0(z) series coefficients are 1/8, 9/128, 75/1024",
          (a1v, a2v, a3v) == (sp.Rational(1, 8), sp.Rational(9, 128), sp.Rational(75, 1024)), eq=eq_bessel)

    # n^sw = sum_{n>=1} [e^{-2 n t} I0(2 n t)]^3 with t = beta S. Cube the series, put z = 2 n t.
    F = 1 + a1v * x + a2v * x**2 + a3v * x**3
    F3 = sp.expand(sp.series(F**3, x, 0, 4).removeO())
    g3 = [F3.coeff(x, j) for j in range(4)]
    pref = (4 * sp.pi) ** sp.Rational(-3, 2)  # (2 pi z)^(-3/2) with z = 2 n t gives (4 pi t n)^(-3/2)
    c_k = [sp.simplify(pref * g3[j] / 2**j) for j in range(4)]  # coefficient of zeta(3/2+j) t^(-3/2-j)

    t = sp.symbols("t", positive=True)  # t = beta S
    lead = c_k[0] * sp.zeta(sp.Rational(3, 2)) * t ** sp.Rational(-3, 2)
    corr = c_k[1] * sp.zeta(sp.Rational(5, 2)) * t ** sp.Rational(-5, 2)
    eq_cor = register(NAME, Equation(
        key="bloch_correction_corollary", entry=ENTRY,
        title="First lattice correction to Bloch's law",
        source_kind="paper", source=f"{C}, Corollary 1.2",
        expr=sp.Eq(sp.Symbol(r"S-m_S(\beta)"), lead + corr),
        latex_str=(r"S-m_S(\beta)=\frac{\zeta(3/2)}{8\pi^{3/2}}(\beta S)^{-3/2}"
                   r"+\frac{3\zeta(5/2)}{128\pi^{3/2}}(\beta S)^{-5/2}+o\bigl((\beta S)^{-5/2}\bigr)"),
        notes="Spin-wave (ideal-magnon) terms only. Dyson's interaction correction enters at order T^4 "
              "(beta^{-4}), beyond this statement.",
    ))
    check(NAME, "leading coefficient from the ideal-magnon sum is 1/(8 pi^(3/2))",
          sp.simplify(c_k[0] - 1 / (8 * sp.pi ** sp.Rational(3, 2))) == 0, eq=eq_cor)
    check(NAME, "first correction coefficient from the ideal-magnon sum is 3/(128 pi^(3/2))",
          sp.simplify(c_k[1] - 3 / (128 * sp.pi ** sp.Rational(3, 2))) == 0, eq=eq_cor)

    eq_next = register(NAME, Equation(
        key="bloch_nsw_next_term", entry=ENTRY,
        title="Next term of the ideal-magnon density (not in the paper)",
        source_kind="derived", source="same heat-kernel expansion as bloch_correction_corollary",
        expr=sp.Eq(sp.Symbol(r"[n^{\mathrm{sw}}_S]_{(\beta S)^{-7/2}}"),
                   c_k[2] * sp.zeta(sp.Rational(7, 2)) * t ** sp.Rational(-7, 2)),
        latex_str=(r"n^{\mathrm{sw}}_S(\beta)=\frac{\zeta(3/2)}{8\pi^{3/2}}(\beta S)^{-3/2}"
                   r"+\frac{3\zeta(5/2)}{128\pi^{3/2}}(\beta S)^{-5/2}"
                   r"+\frac{33\,\zeta(7/2)}{4096\,\pi^{3/2}}(\beta S)^{-7/2}+O\bigl((\beta S)^{-9/2}\bigr)"),
        notes="Coefficient of the ideal-magnon integral only. It is NOT a claim about the magnetization "
              "at this order; the papers prove Theorem 1.1 up to o((beta S)^{-5/2}).",
    ))

    # Numerical test of the expansion of n^sw at moderate t.
    # mp.nsum is inaccurate on this slowly decaying (n^(-3/2)) tail, so sum exactly up to a cutoff
    # and add the tail from the derived series with Hurwitz zeta functions. The series is guarded
    # by comparing it with mpmath's own Bessel function at the cutoff.
    mp.mp.dps = 30
    g3_mp = [mp.mpf(sp.N(gj, 30)) for gj in g3]
    pref_mp = (4 * mp.pi) ** mp.mpf(-1.5)
    CUT = 400

    def term(m, tv):
        return (mp.exp(-2 * m * tv) * mp.besseli(0, 2 * m * tv)) ** 3

    def term_series(m, tv):
        return pref_mp * (m * tv) ** mp.mpf(-1.5) * sum(g3_mp[j] / (2 * m * tv) ** j for j in range(4))

    def nsw(tv):
        tv = mp.mpf(tv)
        direct = mp.fsum(term(m, tv) for m in range(1, CUT + 1))
        tail = pref_mp * tv ** mp.mpf(-1.5) * sum(
            g3_mp[j] * (2 * tv) ** (-j) * mp.zeta(mp.mpf(3) / 2 + j, CUT + 1) for j in range(4))
        return direct + tail

    for tv in (30, 60):
        guard = abs(term_series(CUT, mp.mpf(tv)) / term(CUT, mp.mpf(tv)) - 1)
        check(NAME, f"tail series agrees with the true Bessel term at n = {CUT}, beta S = {tv}",
              guard < mp.mpf("1e-12"), eq=eq_next)
        exact = nsw(tv)
        two_term = (mp.mpf(sp.N(c_k[0], 30)) * mp.zeta(1.5) * mp.mpf(tv) ** (-1.5)
                    + mp.mpf(sp.N(c_k[1], 30)) * mp.zeta(2.5) * mp.mpf(tv) ** (-2.5))
        predicted_residual = mp.mpf(sp.N(c_k[2], 30)) * mp.zeta(3.5) * mp.mpf(tv) ** (-3.5)
        ratio = (exact - two_term) / predicted_residual
        check(NAME, f"numerical n^sw minus two terms matches the derived t^(-7/2) term at beta S = {tv} "
                    f"(ratio {float(ratio):.3f})", abs(ratio - 1) < 0.08, eq=eq_next)

    # ---- [Sph] spherical law -----------------------------------------------------------
    m, u = sp.symbols("m u", positive=True)
    th, ph = sp.symbols("theta phi", real=True)
    mgf = sp.integrate(sp.integrate(sp.exp(u * sp.cos(th)) * sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2 * sp.pi)) / (4 * sp.pi)
    eq_sph = register(NAME, Equation(
        key="bloch_spherical_law", entry=ENTRY,
        title="Spherical magnetization law",
        source_kind="paper", source=f"{SPH}, Theorem 1.1, eq. (1.4)",
        latex_str=(r"\lim_{L\to\infty}\Bigl\langle\exp\Bigl(\tfrac1V\sum_{a=1}^{3}t_aM^a_L\Bigr)\Bigr\rangle_{\beta,L}"
                   r"=\frac{1}{4\pi}\int_{S^2}e^{m\,t\cdot n}\,d\sigma(n)=\frac{\sinh(m|t|)}{m|t|}"),
        notes="Holds for every fixed spin S and beta >= beta_0(S); m = m_{S,beta} is the zero-field right "
              "derivative of the infinite-volume pressure.",
    ))
    check(NAME, "average of exp(u cos theta) over the sphere is sinh(u)/u  (u = m|t|)",
          sp.simplify(mgf - sp.sinh(u) / u) == 0, eq=eq_sph)

    t1, t2, t3 = sp.symbols("t1 t2 t3", real=True)
    r2 = t1**2 + t2**2 + t3**2
    mgf_series = sp.series(sp.sinh(u) / u, u, 0, 4).removeO()  # 1 + u^2/6 + ...
    quad = sp.expand(mgf_series.subs(u, m * sp.sqrt(r2)))
    laplacian = sum(sp.diff(quad, ti, 2) for ti in (t1, t2, t3))
    eq_two = register(NAME, Equation(
        key="bloch_spherical_second_moment", entry=ENTRY,
        title="Mean-square magnetization from the spherical law",
        source_kind="paper", source=f"{SPH}, Theorem 1.1, eq. (1.5)",
        latex_str=(r"\lim_{L\to\infty}\frac{1}{V^{2}}\sum_{x,y\in\Lambda_L}"
                   r"\langle S_x\cdot S_y\rangle_{\beta,L}=m^{2}"),
    ))
    check(NAME, "Laplacian of the t^2 term of sinh(m|t|)/(m|t|) at 0 equals m^2",
          sp.simplify(laplacian.subs({t1: 0, t2: 0, t3: 0}) - m**2) == 0, eq=eq_two)
