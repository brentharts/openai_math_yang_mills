"""Exactly computable SU(2) lattice gauge theories: a 1+1D string tension and a one-plaquette Hamiltonian gap.

This is the numerical bridge asked for in YM_CONNECTIONS.md ("What would be progress", item 3): a lattice-side
statement for SU(2) that can be computed to full precision on a small system, set beside the O(3) results. It is
a finite-volume computation. It has no continuum limit, no scaling window and no weak-coupling plateau of the kind
a four-dimensional gap would need, and it says nothing about existence or a mass gap of Yang-Mills on R^4.

Conventions. Wilson weight exp((beta/2) Tr U) per plaquette with beta = 4/g^2 (Euclidean); Kogut-Susskind
Hamiltonian H = (g^2/2) sum_l E_l^2 + (1/g^2) sum_p (2 - Tr U_p) with E_l^2 = j(j+1) on a link carrying spin j and
Tr in the fundamental, in lattice units. Other texts differ from these by factors of 2 in g^2.

Part A, two dimensions (Euclidean, character expansion). The plaquette weight has character coefficients
(1/d_j) int dU w chi_j = (2/beta) I_{2j+1}(beta), so a ring of L plaquettes has transfer-matrix eigenvalues
(I_{2j+1}/I_1)^L and the string tension is sigma = -ln(I_2/I_1). Checked: the Bessel identity behind it; the Haar
integral by quadrature; the large-beta series sigma = 3/(2b) + 3/(4b^2) + 3/(16b^3) - 9/(16b^4) against the exact
ratio; Casimir scaling (ln of the j = 1 ratio over the j = 1/2 one tends to 8/3); the match with the Hamiltonian
value (g^2/2) j(j+1) per link; and the whole chain (the d_j factors included) against Monte Carlo on the 1 x 1 and
2 x 1 periodic lattices, where Z = sum_j (2 I_{2j+1}/beta)^A.

Part B, one plaquette (Hamiltonian). Gauge-invariant states are class functions of the plaquette holonomy, spanned
by the characters chi_j with all four links in the same spin j. In u = sin(theta) chi_j the problem is
H u = -(g^2/2)(u'' + u) + (2/g^2)(1 - cos theta) u on [0, pi] with u(0) = u(pi) = 0, which is Mathieu's equation, so
E_k = 2/g^2 - g^2/2 + (g^2/8) b_{2k}(8/g^4), k = 2j + 1 = 1, 2, ... Checked: the fusion rule chi_{1/2} chi_j =
chi_{j-1/2} + chi_{j+1/2}, the Casimir and the Haar normalisation, the reduction to Mathieu's equation (symbolically),
that the character-basis matrix IS the Mathieu matrix (symbolically); the gap computed three ways (the tridiagonal
character matrix, a finite-difference solution in theta with Richardson extrapolation, scipy's Mathieu b_{2k}) to
1e-8 over g^2 in [0.05, 8]; truncation convergence; the exact rational strong-coupling series of the gap (computed by
a continued-fraction recursion in exact arithmetic, odd orders vanishing, second order 28/15 matching a hand
calculation) against the numbers; the weak-coupling limit 2 sqrt 2 from the three-dimensional oscillator and its
first correction -5 g^2/8 from first-order perturbation theory in the quartic term (sympy, exact Hermite
integrals), with the next term -41 sqrt 2 g^4/512 from the large-q expansion of Mathieu's equation checked
numerically; and that the gap never closes, with a minimum 2.1962 near g^2 = 0.9875.

Not tested: anything with more than one plaquette, any spatial volume dependence, any link to the continuum limit
or to asymptotic freedom (the two-dimensional gauge theory is super-renormalizable; its string tension has no
dimensional transmutation), the 2+1D and 3+1D theories. Sturm-Liouville simplicity of the spectrum is quoted as a
standard theorem. The next exact step beyond this module is a two-plaquette ladder, which needs Wigner 6j symbols.
"""
import sys
from fractions import Fraction as Fr
from pathlib import Path

import numpy as np
import sympy as sp
import mpmath as mp
from scipy import integrate, special
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import minimize_scalar

if __package__ in (None, ""):          # allow `python3 sympy_models/models/su2_lattice_gap.py`
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register

NAME = "su2_lattice_gap"
ENTRY = "su2-lattice-gap"
KS = "Kogut-Susskind, Phys. Rev. D 11, 395 (1975)"
DZ = "Drouffe-Zuber, Phys. Rep. 102, 1 (1983)"
DLMF = "NIST DLMF ch. 28 (Mathieu functions), sec. 28.8 for large q"


# ---------------------------------------------------------------------------------------------------------
# Part A helpers
# ---------------------------------------------------------------------------------------------------------
def haar_su2(rng, n):
    """n Haar-distributed SU(2) matrices from normalised Gaussian quaternions."""
    q = rng.normal(size=(n, 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    a, b, c, d = q.T
    U = np.empty((n, 2, 2), dtype=complex)
    U[:, 0, 0], U[:, 0, 1] = a + 1j * b, c + 1j * d
    U[:, 1, 0], U[:, 1, 1] = -c + 1j * d, a - 1j * b
    return U


def dag(U):
    return np.conj(np.transpose(U, (0, 2, 1)))


def real_tr(P):
    return np.real(np.trace(P, axis1=1, axis2=2))


# ---------------------------------------------------------------------------------------------------------
# Part B helpers
# ---------------------------------------------------------------------------------------------------------
def trunc_for(g2):
    return 40 if g2 >= 1 else (120 if g2 >= 0.1 else 400)


def levels_jbasis(g2, nlev=4, K=None):
    """Lowest levels of the one-plaquette Hamiltonian in the character basis k = 2j+1 = 1..K (tridiagonal)."""
    K = K or trunc_for(g2)
    k = np.arange(1, K + 1)
    d = 0.5 * g2 * (k**2 - 1) + 2 / g2
    e = np.full(K - 1, -1 / g2)
    return eigh_tridiagonal(d, e, select="i", select_range=(0, nlev - 1), eigvals_only=True)


def _fd(g2, Np1, nlev):
    h = np.pi / Np1
    th = h * np.arange(1, Np1)
    d = 0.5 * g2 * (2 / h**2 - 1) + (2 / g2) * (1 - np.cos(th))
    e = np.full(Np1 - 2, -0.5 * g2 / h**2)
    return eigh_tridiagonal(d, e, select="i", select_range=(0, nlev - 1), eigvals_only=True)


def levels_fd(g2, nlev=4, N1=2000):
    """Finite differences for u(theta) on (0, pi), Richardson-extrapolated in the step h (error O(h^4))."""
    return (4 * _fd(g2, 2 * N1, nlev) - _fd(g2, N1, nlev)) / 3


def levels_mathieu(g2, nlev=4):
    q = 8 / g2**2
    return np.array([2 / g2 - g2 / 2 + (g2 / 8) * special.mathieu_b(2 * k, q) for k in range(1, nlev + 1)])


def gap_num(g2):
    v = levels_jbasis(g2, 2)
    return float(v[1] - v[0])


def strong_series(a, M):
    """Eigenvalue a (k = a+1 in the character labelling, i.e. spin j = a/2) of C2 + eps (2 - chi_{1/2}) as a power series in
    eps through eps^M with exact Fraction coefficients. C2 = j(j+1), the chain is tridiagonal with off-diagonal -eps, and
    the Brillouin-Wigner continued fraction E = d_a + eps^2 (G_up + G_down) is iterated to a fixed point order by order."""
    zero = [Fr(0)] * (M + 1)

    def mul(x, y):
        o = list(zero)
        for i, xi in enumerate(x):
            if xi:
                for j in range(M + 1 - i):
                    o[i + j] += xi * y[j]
        return o

    def inv(x):
        o = [1 / x[0]]
        for n in range(1, M + 1):
            o.append(-sum(x[k] * o[n - k] for k in range(1, n + 1)) / x[0])
        return o

    def sub(x, y):
        return [p - q for p, q in zip(x, y)]

    def add(x, y):
        return [p + q for p, q in zip(x, y)]

    eps = list(zero)
    eps[1] = Fr(1)
    eps2 = mul(eps, eps)

    def d(k):
        c = list(zero)
        c[0], c[1] = Fr(k * (k + 2), 4), Fr(2)
        return c

    E = d(a)
    for _ in range(M + 2):
        G = list(zero)
        for k in range(a + M + 2, a, -1):
            G = inv(sub(sub(E, d(k)), mul(eps2, G)))
        term = mul(eps2, G)
        if a >= 1:
            G = list(zero)
            for k in range(0, a):
                G = inv(sub(sub(E, d(k)), mul(eps2, G)))
            term = add(term, mul(eps2, G))
        E = add(d(a), term)
    return E


def series_gap(M):
    E0, E1 = strong_series(0, M), strong_series(1, M)
    return [b - a for a, b in zip(E0, E1)]


def eval_gap_series(coeffs, g2):
    eps = 1 / (2 * g2**2)
    return 2 * g2 * float(sum(float(c) * eps**n for n, c in enumerate(coeffs)))


def build_and_verify():
    mp.mp.dps = 40
    beta, g2, th = sp.symbols("beta g2 theta", positive=True)

    # =====================================================================================================
    # Part A. Euclidean 2D SU(2): character coefficients, string tension, Casimir scaling
    # =====================================================================================================
    eq_coef = register(NAME, Equation(
        key="su2_char_coeff", entry=ENTRY, title="Character coefficients of the SU(2) Wilson plaquette weight",
        source_kind="textbook", source=f"{DZ} (character expansion); the Bessel form is derived here",
        latex_str=(r"\frac{1}{d_j}\int_{SU(2)}\!dU\,e^{\frac{\beta}{2}\operatorname{Tr}U}\chi_j(U)"
                   r"=\frac{2}{\beta}\,I_{2j+1}(\beta),\qquad d_j=2j+1,\quad \beta=\frac{4}{g^{2}}"),
        notes="Haar measure (2/pi) sin^2(theta) dtheta on the class angle theta in [0, pi], Tr U = 2 cos theta and "
              "chi_j = sin((2j+1) theta)/sin theta. The integral equals I_{2j} - I_{2j+2} and the Bessel recurrence "
              "turns this into (2(2j+1)/beta) I_{2j+1}. Two-dimensional Euclidean lattice, no matter.",
    ))
    pts = [(n, b) for n in range(1, 7) for b in ("1/2", "2", "7", "30")]
    def bessel_gap(n, b):
        x = mp.mpf(sp.Rational(b).p) / mp.mpf(sp.Rational(b).q)
        lhs = mp.besseli(n - 1, x) - mp.besseli(n + 1, x)
        return abs(lhs - 2 * n / x * mp.besseli(n, x)) / abs(lhs)
    rec = all(bessel_gap(n, b) < mp.mpf(10) ** -30 for n, b in pts)
    check(NAME, "Bessel recurrence I_{n-1} - I_{n+1} = (2n/beta) I_n to 30 relative digits for n = 1..6 and beta = 1/2, 2, 7, 30",
          rec, eq=eq_coef)
    quad_ok = True
    for b in (0.7, 2.5, 6.0):
        for n in range(1, 7):
            # (2/pi) sin^2(t) e^{b cos t} chi_j with chi_j = sin(nt)/sin(t): the sin^2/sin leaves sin(t) sin(nt)
            val, _ = integrate.quad(lambda t: (2 / np.pi) * np.sin(t) * np.sin(n * t) * np.exp(b * np.cos(t)),
                                    0, np.pi, epsabs=1e-13, epsrel=1e-13)
            quad_ok &= abs(val / n - 2 * float(mp.besseli(n, b)) / b) < 1e-10 * max(1, abs(val))
    check(NAME, "Haar quadrature of exp(beta cos t) chi_j reproduces (2/beta) I_{2j+1} for 2j+1 = 1..6, beta = 0.7, 2.5, 6",
          quad_ok, eq=eq_coef)

    # asymptotic series of sigma = -ln(I2/I1) from the standard large-argument expansion of I_nu
    w = sp.symbols("w", positive=True)

    def i_asym(nu, K=6):
        s = 0
        for k in range(K + 1):
            a = sp.Integer(1)
            for m in range(1, k + 1):
                a *= 4 * nu**2 - (2 * m - 1) ** 2
            s += (-1) ** k * a / (sp.factorial(k) * 8**k) * w**k
        return s
    sig_series = sp.expand(sp.series(-sp.log(i_asym(2) / i_asym(1)), w, 0, 5).removeO())
    sig_expected = sp.Rational(3, 2) * w + sp.Rational(3, 4) * w**2 + sp.Rational(3, 16) * w**3 - sp.Rational(9, 16) * w**4

    eq_sigma = register(NAME, Equation(
        key="su2_string_tension_2d", entry=ENTRY,
        title="Two-dimensional SU(2): string tension, ring spectrum and Casimir scaling",
        source_kind="derived", source=f"derived here from su2_char_coeff; {DZ}; torus partition function by the heat-kernel/character sum",
        latex_str=(r"\sigma(\beta)=-\ln\frac{I_2(\beta)}{I_1(\beta)}=\frac{3}{2\beta}+\frac{3}{4\beta^{2}}+\frac{3}{16\beta^{3}}"
                   r"-\frac{9}{16\beta^{4}}+O(\beta^{-5}),%%BR%%"
                   r"E_j(L)=-L\ln\frac{I_{2j+1}(\beta)}{I_1(\beta)},\qquad Z_{\mathrm{torus}}(A)=\sum_{j}\Bigl(\frac{2I_{2j+1}(\beta)}{\beta}\Bigr)^{A},"
                   r"%%BR%%\frac{E_1}{E_{1/2}}\longrightarrow\frac{1\cdot2}{\frac12\cdot\frac32}=\frac83\ \ (\beta\to\infty)"),
        notes="Exactly solvable and confining at every beta, with a string tension, an area law and Casimir scaling. "
              "The gap on a ring of L plaquettes grows linearly with L, so there is no infinite-volume gap in the "
              "sense of a finite mass: this is the opposite of the two-dimensional O(3) model, whose gap is finite "
              "and appears by dimensional transmutation. The continuum limit (beta -> infinity with beta a^2 fixed) "
              "exists trivially because the theory is super-renormalizable, with sigma = g^2 a^2 (3/8) in lattice "
              "units. No asymptotic freedom, no running coupling.",
    ))
    check(NAME, "large-argument series of the Bessel ratio gives sigma = 3/(2b) + 3/(4b^2) + 3/(16b^3) - 9/(16b^4) + O(b^-5)",
          sp.simplify(sig_series - sig_expected) == 0, eq=eq_sigma)
    ok = True
    for z, tol in ((50, 3e-8), (200, 1e-10), (1000, 1e-13)):
        exact = -mp.log(mp.besseli(2, z) / mp.besseli(1, z))
        ok &= abs(float(exact) - float(sig_expected.subs(w, sp.Rational(1, z)))) < tol
    check(NAME, "the series agrees with -ln(I_2/I_1) to O(beta^-5) at beta = 50, 200, 1000 (errors 6e-9, 6e-12, 2e-15)",
          ok, eq=eq_sigma)
    casimir_ok = True
    for b in (1e3, 1e4, 1e5):
        r = mp.log(mp.besseli(3, b) / mp.besseli(1, b)) / mp.log(mp.besseli(2, b) / mp.besseli(1, b))
        casimir_ok &= abs(float(r) - 8 / 3) < 3.0 / b
    check(NAME, "Casimir scaling: ln(I_3/I_1)/ln(I_2/I_1) -> 8/3 = C2(1)/C2(1/2) with error below 3/beta",
          casimir_ok, eq=eq_sigma)
    order_ok = all(float(mp.besseli(n, b)) > float(mp.besseli(n + 1, b))
                   for n in range(1, 12) for b in (0.1, 1.0, 5.0, 20.0, 100.0))
    check(NAME, "I_n(beta) decreases in n, so the ring levels E_j(L) increase with j: the lowest excitation is j = 1/2",
          order_ok, eq=eq_sigma)
    gg = sp.symbols("g", positive=True)
    check(NAME, "beta = 4/g^2: sigma -> 3/(2 beta) = 3 g^2/8 = (g^2/2) j(j+1) at j = 1/2, the Hamiltonian electric energy per link",
          sp.simplify(sp.Rational(3, 2) / (4 / gg**2) - (gg**2 / 2) * sp.Rational(1, 2) * sp.Rational(3, 2)) == 0, eq=eq_sigma)

    # Monte Carlo on small periodic lattices tests the whole chain, d_j factors included
    rng = np.random.default_rng(20261031)
    N = 400_000
    bmc = 1.5
    U1, U2, V1, V2 = (haar_su2(rng, N) for _ in range(4))
    w1 = np.exp(bmc / 2 * real_tr(U1 @ V1 @ dag(U1) @ dag(V1)))
    w2 = (np.exp(bmc / 2 * real_tr(U1 @ V2 @ dag(U1) @ dag(V1)))
          * np.exp(bmc / 2 * real_tr(U2 @ V1 @ dag(U2) @ dag(V2))))
    for A, wv in ((1, w1), (2, w2)):
        exact = sum((2 * float(mp.besseli(n, bmc)) / bmc) ** A for n in range(1, 40))
        err = wv.std() / np.sqrt(N)
        check(NAME, f"Monte Carlo on the {'1 x 1' if A == 1 else '2 x 1'} periodic lattice (beta = 1.5, 4e5 samples): "
                    f"Z = sum_j (2 I_{{2j+1}}/beta)^{A} within 4 standard errors",
              abs(float(wv.mean()) - exact) < 4 * err, eq=eq_sigma)

    # =====================================================================================================
    # Part B. One plaquette, Hamiltonian
    # =====================================================================================================
    k = sp.symbols("k", positive=True, integer=True)
    eq_ham = register(NAME, Equation(
        key="su2_plaq_hamiltonian", entry=ENTRY, title="One-plaquette SU(2) Hamiltonian in the character basis",
        source_kind="textbook", source=f"{KS}; the reduction to one plaquette and the tridiagonal form are derived here",
        latex_str=(r"H=\frac{g^{2}}{2}\sum_{l=1}^{4}E_l^{2}+\frac{1}{g^{2}}\bigl(2-\operatorname{Tr}U_p\bigr),%%BR%%"
                   r"H_{kk}=\frac{g^{2}}{2}(k^{2}-1)+\frac{2}{g^{2}},\quad H_{k,k\pm1}=-\frac{1}{g^{2}},\qquad k=2j+1=1,2,\dots"),
        notes="Four links, four vertices, open boundary. Gauge invariance leaves class functions of the plaquette "
              "holonomy, spanned by the characters; every link carries the same spin j, so sum E^2 = 4 j(j+1) and "
              "Tr U_p multiplies chi_j by chi_{1/2}. This is one spatial plaquette of the 2+1D theory in lattice "
              "units, with no dependence on any volume.",
    ))
    t = sp.symbols("t", real=True)
    fusion = all(sp.expand((2 * sp.cos(t) * sp.sin(kk * t) - sp.sin((kk - 1) * t) - sp.sin((kk + 1) * t)).rewrite(sp.exp)) == 0
                 for kk in range(1, 12))
    check(NAME, "chi_{1/2} chi_j = chi_{j-1/2} + chi_{j+1/2}: 2 cos t sin(kt) = sin((k-1)t) + sin((k+1)t) for k = 1..11 "
                "(the k = 1 term has no j = -1/2 partner because sin 0 = 0)", fusion, eq=eq_ham)
    u = sp.Function("u")
    lap = sp.simplify(sp.diff(sp.sin(t) ** 2 * sp.diff(u(t) / sp.sin(t), t), t) / sp.sin(t) ** 2 - (sp.diff(u(t), t, 2) + u(t)) / sp.sin(t))
    check(NAME, "Laplacian of S^3 on a class function f = u/sin t is (u'' + u)/sin t", lap == 0, eq=eq_ham)
    j_sym = sp.symbols("j", positive=True)
    cas = sp.simplify(-(g2 / 2) * (sp.diff(sp.sin(k * t), t, 2) + sp.sin(k * t)) - (g2 / 2) * (k**2 - 1) * sp.sin(k * t))
    check(NAME, "-(g^2/2)(u'' + u) on sin(kt) is (g^2/2)(k^2 - 1) sin(kt) = 2 g^2 j(j+1) sin(kt) with k = 2j+1 (four links, g^2/2 each)",
          cas == 0 and sp.simplify(((2 * j_sym + 1) ** 2 - 1) / 4 - j_sym * (j_sym + 1)) == 0, eq=eq_ham)
    orth = all(sp.integrate(2 / sp.pi * sp.sin(a * t) * sp.sin(b * t), (t, 0, sp.pi)) == (1 if a == b else 0)
               for a in range(1, 7) for b in range(1, 7))
    check(NAME, "the characters are orthonormal under (2/pi) sin^2 t dt: (2/pi) int sin(at) sin(bt) = delta_ab for a, b = 1..6",
          orth, eq=eq_ham)
    check(NAME, "g^2 -> infinity: the first gap is 4 links x (g^2/2) j(j+1) at j = 1/2 = 3 g^2/2, equal to the 2D ring gap L sigma "
                "with L = 4 and sigma = 3 g^2/8",
          sp.simplify(4 * sp.Rational(3, 8) * g2 - sp.Rational(3, 2) * g2) == 0, eq=eq_ham)

    eq_math = register(NAME, Equation(
        key="su2_plaq_mathieu", entry=ENTRY, title="The one-plaquette spectrum is a Mathieu spectrum",
        source_kind="derived", source=f"derived here from su2_plaq_hamiltonian; {DLMF}",
        latex_str=(r"y''+\bigl(a-2q\cos 2x\bigr)y=0,\ \ y(0)=y(\tfrac\pi2)=0,\quad a=\frac{8}{g^{2}}\Bigl(E+\frac{g^{2}}{2}-\frac{2}{g^{2}}\Bigr),\ \ q=-\frac{8}{g^{4}}"
                   r"%%BR%%E_k=\frac{2}{g^{2}}-\frac{g^{2}}{2}+\frac{g^{2}}{8}\,b_{2k}\!\Bigl(\frac{8}{g^{4}}\Bigr),\qquad "
                   r"\Delta(g)=E_2-E_1=\frac{g^{2}}{8}\bigl[b_4(q)-b_2(q)\bigr],\ \ q=\frac{8}{g^{4}}"),
        notes="theta = 2x. The odd, period-pi Mathieu functions se_{2k} vanish at 0 and pi/2 and are the allowed "
              "states; b_{2k}(-q) = b_{2k}(q) for even order, so q may be taken positive. In the sine basis "
              "sin(2mx) the Mathieu eigenproblem is the tridiagonal matrix with diagonal 4m^2 and off-diagonal q, "
              "which is exactly the character-basis matrix above after the shift and scale. The gap is Delta = E_{j=1/2} - E_{j=0}.",
    ))
    E_, a_sym, q_sym = sp.symbols("E a q")
    y = sp.Function("y")
    xx = sp.symbols("x", real=True)
    a_val = 8 / g2 * (E_ + g2 / 2 - 2 / g2)
    q_val = -8 / g2**2
    S = -(g2 / 2) * (sp.diff(y(xx), xx, 2) / 4 + y(xx)) + (2 / g2) * (1 - sp.cos(2 * xx)) * y(xx) - E_ * y(xx)
    Mop = sp.diff(y(xx), xx, 2) + (a_val - 2 * q_val * sp.cos(2 * xx)) * y(xx)
    check(NAME, "with theta = 2x, (H - E) u = -(g^2/8) [y'' + (a - 2 q cos 2x) y] for a = (8/g^2)(E + g^2/2 - 2/g^2), q = -8/g^4: "
                "the identity holds for a generic function y",
          sp.simplify(S + (g2 / 8) * Mop) == 0, eq=eq_math)
    mm = sp.symbols("m", positive=True, integer=True)
    diag_ok = sp.simplify((g2 / 8) * 4 * mm**2 - g2 / 2 - (g2 / 2 * (mm**2 - 1))) == 0     # shift -g^2/2 plus (g^2/8) b with diag 4m^2
    off_ok = sp.simplify((g2 / 8) * q_val + 1 / g2) == 0                                      # (g^2/8) q = -1/g^2
    check(NAME, "the Mathieu matrix (diagonal 4m^2, off-diagonal q = -8/g^4) maps onto the character matrix: "
                "(g^2/8) 4m^2 - g^2/2 = (g^2/2)(m^2 - 1) and (g^2/8) q = -1/g^2", diag_ok and off_ok, eq=eq_math)

    grid = (0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
    err_fd = max(float(np.max(abs(levels_jbasis(g) - levels_fd(g)))) for g in grid + (0.05,))
    err_m = max(float(np.max(abs(levels_jbasis(g) - levels_mathieu(g)))) for g in grid)
    check(NAME, f"lowest four levels: character matrix against finite differences with Richardson extrapolation, "
                f"g^2 in {{0.05, 0.1, ..., 8}}, within 1e-7 (largest difference {err_fd:.1e})", err_fd < 1e-7, eq=eq_math)
    check(NAME, f"lowest four levels: character matrix against scipy's Mathieu b_2k at q = 8/g^4 (q up to 800), "
                f"g^2 in {{0.1, ..., 8}}, within 1e-8 (largest difference {err_m:.1e})", err_m < 1e-8, eq=eq_math)
    conv = max(abs((lambda a, b: (a[1] - a[0]) - (b[1] - b[0]))(levels_jbasis(g, 2, K), levels_jbasis(g, 2, 2 * K)))
               for g in (0.05, 0.25, 1.0, 4.0) for K in (160,))
    check(NAME, f"truncation convergence: the gap at K = 160 and K = 320 characters agree to 1e-9 (largest difference {conv:.1e})",
          conv < 1e-9, eq=eq_math)
    mins = minimize_scalar(lambda lg: gap_num(float(np.exp(lg))), bounds=(np.log(0.2), np.log(5)), method="bounded",
                           options={"xatol": 1e-10})
    g2_min, gap_min = float(np.exp(mins.x)), float(mins.fun)
    gs = np.exp(np.linspace(np.log(1e-2), np.log(1e3), 400))
    gaps = np.array([gap_num(float(g)) for g in gs])
    # ---- strong coupling: exact rational series ---------------------------------------------------------
    M = 14
    coeffs = series_gap(M)
    E0s = strong_series(0, 8)
    higher = "".join((r" + " if coeffs[n] > 0 else r" - ") + sp.latex(abs(sp.Rational(coeffs[n].numerator, coeffs[n].denominator)))
                     + rf"\,\epsilon^{{{n}}}" for n in (4, 6))
    even = [coeffs[n] for n in range(2, M + 1, 2)]
    ratios = [float(abs(even[i + 1] / even[i])) for i in range(len(even) - 1)]
    eq_strong = register(NAME, Equation(
        key="su2_plaq_strong", entry=ENTRY, title="Strong-coupling series of the one-plaquette gap (exact rational coefficients)",
        source_kind="derived", source="derived here from su2_plaq_hamiltonian by Brillouin-Wigner perturbation theory in exact arithmetic",
        latex_str=(r"\Delta=2g^{2}\Bigl[\tfrac34+\tfrac{28}{15}\epsilon^{2}" + higher
                   + r"+\cdots\Bigr],\qquad \epsilon=\frac{1}{2g^{4}},\ \ E_0=2g^{2}\bigl[2\epsilon-\tfrac43\epsilon^{2}+\tfrac{40}{27}\epsilon^{4}-\cdots\bigr]"),
        notes="Units: H/(2 g^2) = j(j+1) + epsilon (2 - chi_{1/2}). Only even orders occur because chi_{1/2} -> -chi_{1/2} "
              "is a unitary change of sign on the odd characters. Second order is 28/15 = 4/3 (ground state, via j = 1/2) "
              "+ 4/3 - 4/5 (j = 1/2, via j = 0 and j = 1). At g^2 = 3 the error is about 1e-8, 2e-12 and 1e-13 after "
              "6, 10 and 14 orders. The ratios of successive even coefficients rise from 1.6 to 4.2 over the 14 orders, "
              "which suggests a finite radius of convergence in epsilon^2 of order 0.2, i.e. g^2 near 1 (an estimate "
              "from seven coefficients, not a proof; the gap minimum sits at g^2 = 0.9875). It is an expansion of one "
              "eigenvalue of one plaquette, not of a field theory.",
    ))
    check(NAME, "gap series starts 3/4 + (28/15) eps^2 and the E_0 series 2 eps - (4/3) eps^2: the hand second-order values "
                "(-1/(3/4); +1/(3/4) - 1/(5/4) for j = 1/2)",
          coeffs[0] == Fr(3, 4) and coeffs[2] == Fr(28, 15) and E0s[1] == 2 and E0s[2] == Fr(-4, 3), eq=eq_strong)
    check(NAME, "all odd orders of the gap series and of E_0 beyond the shift 2 eps vanish (sign symmetry of the chain)",
          all(coeffs[n] == 0 for n in range(1, M + 1, 2)) and all(E0s[n] == 0 for n in (3, 5, 7)), eq=eq_strong)
    errs = {}
    for Mt in (6, 10, 14):
        errs[Mt] = abs(eval_gap_series(series_gap(Mt), 3.0) - gap_num(3.0))
    check(NAME, f"the series converges to the numerical gap at g^2 = 3: errors {errs[6]:.1e}, {errs[10]:.1e}, {errs[14]:.1e} "
                f"for 6, 10, 14 orders, decreasing, last below 1e-10",
          errs[6] > errs[10] > errs[14] and errs[14] < 1e-10, eq=eq_strong)
    check(NAME, "at g^2 = 2 and 4 the 14-order series is within 1e-8 of the numerical gap",
          all(abs(eval_gap_series(coeffs, g) - gap_num(g)) < 1e-8 for g in (2.0, 4.0)), eq=eq_strong)

    # ---- weak coupling ------------------------------------------------------------------------------------
    eq_weak = register(NAME, Equation(
        key="su2_plaq_weak", entry=ENTRY, title="Weak-coupling limit of the one-plaquette gap",
        source_kind="derived", source=f"derived here; the g^4 term from {DLMF}",
        latex_str=(r"\Delta(g)=2\sqrt2-\frac58\,g^{2}-\frac{41\sqrt2}{512}\,g^{4}+O(g^{6}),\qquad "
                   r"\omega=\sqrt2,\ \ \Delta_0=2\omega"),
        notes="Near U = 1 the class angle theta obeys H ~ p^2 g^2/2 + theta^2/g^2: an oscillator of mass 1/g^2 and "
              "frequency sqrt 2 independent of g, restricted to odd states by u(0) = 0, so the first gap is 2 omega. "
              "The quartic term of 2(1 - cos theta)/g^2 gives -5 g^2/8 at first order (computed with exact Hermite "
              "integrals); the g^4 coefficient is read from the large-q expansion of Mathieu's equation and checked "
              "numerically only. Weak coupling here means a small plaquette in a fixed lattice unit, which is not "
              "a continuum limit.",
    ))
    m_, om = sp.symbols("m omega", positive=True)
    om_sol = sp.solve(sp.Eq(sp.Rational(1, 2) * (1 / g2) * om**2, 1 / g2), om)
    check(NAME, "oscillator: kinetic g^2 p^2/2 gives mass 1/g^2, potential theta^2/g^2 = (1/2) m omega^2 theta^2 gives omega = sqrt 2, "
                "and odd states n = 1, 3 are 2 omega apart: 2 sqrt 2", om_sol == [sp.sqrt(2)] and 2 * om_sol[0] == 2 * sp.sqrt(2), eq=eq_weak)
    tt = sp.symbols("tt", real=True)
    alpha = sp.sqrt(2) / g2

    def psi(n):
        return (alpha / sp.pi) ** sp.Rational(1, 4) / sp.sqrt(2**n * sp.factorial(n)) * sp.hermite(n, sp.sqrt(alpha) * tt) * sp.exp(-alpha * tt**2 / 2)
    shift = {}
    for n in (1, 3):
        nrm = sp.integrate(sp.simplify(psi(n) ** 2), (tt, -sp.oo, sp.oo))
        t4 = sp.integrate(sp.simplify(psi(n) ** 2 * tt**4), (tt, -sp.oo, sp.oo))
        shift[n] = sp.simplify(-t4 / (12 * g2))
        want = {1: -5 * g2 / 32, 3: -25 * g2 / 32}[n]
        check(NAME, f"oscillator state n = {n} is normalised and has first-order quartic shift {want} (from -theta^4/(12 g^2), "
                    f"exact Hermite integrals)",
              sp.simplify(nrm - 1) == 0 and sp.simplify(shift[n] - want) == 0, eq=eq_weak)
    check(NAME, "first-order quartic correction to the gap is (-25 + 5) g^2/32 = -5 g^2/8",
          sp.simplify(shift[3] - shift[1] + sp.Rational(5, 8) * g2) == 0, eq=eq_weak)
    check(NAME, "the cosine expansion 2(1 - cos t) = t^2 - t^4/12 + ... used above",
          sp.series(2 * (1 - sp.cos(tt)), tt, 0, 6).removeO() == tt**2 - tt**4 / 12, eq=eq_weak)
    c4 = -41 * np.sqrt(2) / 512
    rel = []
    for g in (0.01, 0.02):
        resid = (gap_num(g) - (2 * np.sqrt(2) - 5 * g / 8)) / g**2
        rel.append(resid / c4)
    check(NAME, f"numerical gap minus 2 sqrt 2 + 5 g^2/8, divided by g^4, is -41 sqrt 2/512 = {c4:.5f} to within 2% at g^2 = 0.01, 0.02 "
                f"(ratios {rel[0]:.4f}, {rel[1]:.4f})", all(abs(r - 1) < 0.02 for r in rel), eq=eq_weak)
    check(NAME, "the numerical gap tends to 2 sqrt 2: at g^2 = 1e-3 it is within 1e-3 of 2 sqrt 2",
          abs(gap_num(1e-3) - 2 * np.sqrt(2)) < 1e-3, eq=eq_weak)
    # ---- the whole curve: no closing, one interior minimum -----------------------------------------------
    eq_min = register(NAME, Equation(
        key="su2_plaq_gap_min", entry=ENTRY, title="The one-plaquette gap never closes and has an interior minimum",
        source_kind="derived", source="derived here from su2_plaq_mathieu, su2_plaq_strong and su2_plaq_weak (numerical)",
        latex_str=(r"\Delta(g)>0\ \ \forall g>0,\qquad \min_{g^{2}}\Delta=" + f"{gap_min:.4f}" + r"\ \text{ at }\ g^{2}\approx"
                   + f"{g2_min:.4f}" + r",\qquad \Delta\to 2\sqrt2\ (g\to0),\ \ \Delta\sim\tfrac32 g^{2}\ (g\to\infty)"),
        notes="Positivity for every g follows from Sturm-Liouville theory (a regular problem on a bounded interval has "
              "simple eigenvalues; quoted, not reproved) and is checked on a grid of 400 values of g^2 from 0.01 to "
              "1000. The minimum is numerical. The gap falls from 2 sqrt 2 at weak coupling, reaches its minimum near "
              "g^2 = 1 where neither expansion converges well, and then rises as 3 g^2/2. A positive gap in a system "
              "with one degree of freedom is not evidence about the mass gap of a field theory.",
    ))
    check(NAME, "the gap is positive on 400 log-spaced values of g^2 in [0.01, 1000]", bool(np.all(gaps > 0)), eq=eq_min)
    check(NAME, "the minimum gap is 2.1962 at g^2 = 0.9875 (to 1e-3), below the weak-coupling value 2 sqrt 2 and "
                "below every gap on the grid",
          abs(gap_min - 2.19625) < 1e-3 and abs(g2_min - 0.98748) < 1e-3 and gap_min < 2 * np.sqrt(2)
          and bool(np.all(gaps >= gap_min - 1e-9)), eq=eq_min)
    check(NAME, "the gap is decreasing on the weak side and increasing on the strong side of the minimum "
                "(g^2 = 0.05 .. 0.9 and 1.1 .. 8, 12 points each side)",
          bool(np.all(np.diff([gap_num(g) for g in np.linspace(0.05, 0.9, 12)]) < 0)
               and np.all(np.diff([gap_num(g) for g in np.linspace(1.1, 8, 12)]) > 0)), eq=eq_min)
    check(NAME, "the strong-coupling gap 3 g^2/2 is within 1e-3 relative of the numerical gap at g^2 = 100 and the weak value "
                "2 sqrt 2 within 1e-3 absolute at g^2 = 1e-3 (the two asymptotes)",
          abs(gap_num(100.0) / 150.0 - 1) < 1e-3 and abs(gap_num(1e-3) - 2 * np.sqrt(2)) < 1e-3, eq=eq_min)

    g2s = (0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
    return {g: gap_num(g) for g in g2s}


if __name__ == "__main__":
    tab = build_and_verify()
    print("g^2     gap")
    for g, v in tab.items():
        print(f"{g:6.3f}  {v:.10f}")
