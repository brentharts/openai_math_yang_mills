"""The BPST instanton and the Hopf fibration: a classical, textbook check next to the Hopf models.

This is the one place where the Hopf material touches gauge theory, so the limits come first.

What is checked: the classical SU(2) Yang-Mills equations on Euclidean R^4 have the explicit
Belavin-Polyakov-Schwartz-Tyupkin solution A^a_mu = 2 eta^a_{mu nu} x_nu / (x^2 + rho^2). Its field strength
is self-dual, its action density is 192 rho^4 / (x^2 + rho^2)^4 in the normalisation below, and its second
Chern number is exactly 1. The pure-gauge map at infinity, x/|x| : S^3 -> SU(2) = S^3, is the identity of S^3
(degree 1), the same S^3 that carries the Hopf fibration over S^2 = SU(2)/U(1).

What is not claimed: nothing here constructs a quantum Yang-Mills theory, proves existence, or bears on a
mass gap. An instanton is a classical solution; its integer is a topological number, not a spectral
quantity. The Hopf insulator's band gap (models/hopf_insulator.py) is a different object again. The
repository states the same limits in YANGMILLS.md.

Source: Belavin, Polyakov, Schwartz and Tyupkin, Phys. Lett. B 59, 85 (1975); 't Hooft symbols as in any
instanton review. Standard textbook material, rederived here; labelled Standard, not new.
"""
import sympy as sp

from eqlib import Equation, check, register

NAME = "hopf_gauge"
ENTRY = "hopf-fibration"
SRC = "Belavin, Polyakov, Schwartz and Tyupkin, Phys. Lett. B 59, 85 (1975)"


def eta(a, mu, nu):
    """'t Hooft symbol eta^a_{mu nu}, indices 0-based: a, mu, nu in 0..2 are space, 3 is the fourth axis."""
    if mu < 3 and nu < 3:
        return sp.LeviCivita(a, mu, nu)
    if nu == 3 and mu < 3:
        return 1 if a == mu else 0
    if mu == 3 and nu < 3:
        return -1 if a == nu else 0
    return 0


def build_and_verify():
    x = sp.symbols("x1 x2 x3 x4", real=True)
    rho = sp.symbols("rho", positive=True)
    r2 = sum(v**2 for v in x)
    den = r2 + rho**2

    A = [[2 * sum(eta(a, mu, nu) * x[nu] for nu in range(4)) / den for mu in range(4)] for a in range(3)]
    F = [[[0] * 4 for _ in range(4)] for _ in range(3)]
    for a in range(3):
        b, c = (a + 1) % 3, (a + 2) % 3
        for mu in range(4):
            for nu in range(4):
                # F^a = d_mu A^a_nu - d_nu A^a_mu + eps^{abc} A^b_mu A^c_nu  (coupling g = 1)
                F[a][mu][nu] = sp.simplify(sp.diff(A[a][nu], x[mu]) - sp.diff(A[a][mu], x[nu])
                                           + sp.LeviCivita(a, b, c) * (A[b][mu] * A[c][nu] - A[c][mu] * A[b][nu]))

    eq_f = register(NAME, Equation(
        key="hg_bpst_field", entry=ENTRY, cites=["bpst-1975"], title="The BPST instanton: A, F and self-duality (classical SU(2))",
        source_kind="textbook", source=SRC,
        latex_str=(r"A^a_\mu=\frac{2\,\eta^a_{\mu\nu}x_\nu}{x^2+\rho^2},\qquad "
                   r"F^a_{\mu\nu}=-\frac{4\,\eta^a_{\mu\nu}\rho^2}{(x^2+\rho^2)^2},\qquad F=\tilde F"),
        notes="Classical, Euclidean, SU(2), coupling g = 1. Standard textbook material: it is not a new result and says "
              "nothing about quantum Yang-Mills existence or a mass gap.",
    ))
    check(NAME, "F^a_{mu nu} = -4 eta^a_{mu nu} rho^2 / (x^2 + rho^2)^2 for all a, mu, nu",
          all(sp.simplify(F[a][m][n] + 4 * eta(a, m, n) * rho**2 / den**2) == 0
              for a in range(3) for m in range(4) for n in range(4)), eq=eq_f)
    dual = lambda a, m, n: sum(sp.LeviCivita(m, n, p, q) * F[a][p][q] for p in range(4) for q in range(4)) / 2
    check(NAME, "F is self-dual: (1/2) eps_{mu nu rho sigma} F^a_{rho sigma} = F^a_{mu nu}",
          all(sp.simplify(dual(a, m, n) - F[a][m][n]) == 0 for a in range(3) for m in range(4) for n in range(4)),
          eq=eq_f)
    # Yang-Mills equation D_mu F_{mu nu} = 0; self-duality plus the Bianchi identity implies it, checked directly
    def DF(a, nu):
        b, c = (a + 1) % 3, (a + 2) % 3
        return sum(sp.diff(F[a][mu][nu], x[mu])
                   + sp.LeviCivita(a, b, c) * (A[b][mu] * F[c][mu][nu] - A[c][mu] * F[b][mu][nu]) for mu in range(4))
    check(NAME, "the Yang-Mills equation D_mu F_{mu nu} = 0 holds",
          all(sp.simplify(DF(a, nu)) == 0 for a in range(3) for nu in range(4)), eq=eq_f)

    dens = sp.simplify(sum(F[a][m][n] * dual(a, m, n) for a in range(3) for m in range(4) for n in range(4)))
    rr = sp.symbols("r", positive=True)
    vol = 2 * sp.pi**2 * rr**3                                   # area of the 3-sphere of radius r times dr
    # dens is radial (checked equal to the closed form below), so the integral is the radial one times the S^3 area
    c2 = sp.integrate((192 * rho**4 / (rr**2 + rho**2) ** 4) * vol, (rr, 0, sp.oo)) / (32 * sp.pi**2)
    eq_c2 = register(NAME, Equation(
        key="hg_bpst_c2", entry=ENTRY, cites=["bpst-1975"], title="The second Chern number of the BPST instanton is 1",
        source_kind="textbook", source=SRC,
        latex_str=(r"F^a_{\mu\nu}\tilde F^a_{\mu\nu}=\frac{192\,\rho^4}{(x^2+\rho^2)^4},\qquad "
                   r"c_2=\frac{1}{32\pi^2}\int_{\mathbb R^4}F^a_{\mu\nu}\tilde F^a_{\mu\nu}\,d^4x=1"),
        notes="Counted: the integral is done exactly, independent of rho. The integer is the winding of the pure-gauge "
              "map at infinity, a topological number of a classical field. It is not a mass, a gap, or a statement about "
              "the quantum theory.",
    ))
    check(NAME, "F^a F~^a = 192 rho^4 / (x^2 + rho^2)^4 pointwise",
          sp.simplify(dens - 192 * rho**4 / den**4) == 0, eq=eq_c2)
    check(NAME, "c_2 = (1/32 pi^2) Int F F~ d^4x = 1 for every rho > 0", sp.simplify(c2) == 1, eq=eq_c2)

    # ---- the boundary map is the identity of S^3 (the Hopf total space) -----------------------
    # degree of g = (x4 + i x.sigma)/|x| from S^3 to SU(2) = S^3: with g in unit-quaternion coordinates (x1..x4)/|x|,
    # g is the identity map of S^3, whose degree is 1; verify via the Jacobian of the map x -> x/|x| on the sphere:
    # for the identity the pull-back of the volume form is the volume form, so the determinant of the 4x4 matrix
    # (n, d_1 n, d_2 n, d_3 n) in any chart equals the sphere volume density.
    th1, th2, ph = sp.symbols("theta1 theta2 phi", real=True)
    n = sp.Matrix([sp.cos(th1), sp.sin(th1) * sp.cos(th2), sp.sin(th1) * sp.sin(th2) * sp.cos(ph),
                   sp.sin(th1) * sp.sin(th2) * sp.sin(ph)])
    J = sp.Matrix.hstack(n, n.diff(th1), n.diff(th2), n.diff(ph))
    vol_density = sp.simplify(J.det())
    eq_deg = register(NAME, Equation(
        key="hg_boundary_degree", entry=ENTRY, title="The pure gauge at infinity is a degree-1 map S^3 -> S^3",
        source_kind="derived", source="derived here; the same S^3 as the Hopf total space (models/hopf_fibration.py)",
        latex_str=(r"A\to g^{-1}dg,\ \ g=\frac{x_4+i\,\vec x\cdot\vec\sigma}{|x|}:\ S^3\to SU(2)\cong S^3,\ \ \deg g=1"),
        notes="SU(2) is the unit quaternions, i.e. the S^3 on which the Hopf map h acts; U(1) acts on its fibres. Degree 1 "
              "because g is the identity of S^3 in the coordinates (x1, x2, x3, x4)/|x|, up to a fixed orientation. This "
              "fixes the sign convention only; it does not connect the instanton number to the Hopf invariant "
              "of any map.",
    ))
    check(NAME, "det(n, d_theta1 n, d_theta2 n, d_phi n) = sin^2(theta1) sin(theta2): the identity map has the sphere's "
                "volume density",
          sp.simplify(vol_density - sp.sin(th1) ** 2 * sp.sin(th2)) == 0, eq=eq_deg)
    total = sp.integrate(sp.sin(th1) ** 2 * sp.sin(th2), (th1, 0, sp.pi), (th2, 0, sp.pi), (ph, 0, 2 * sp.pi))
    check(NAME, "the integral of that density is 2 pi^2, the volume of S^3, so the degree is 2 pi^2 / 2 pi^2 = 1",
          sp.simplify(total) == 2 * sp.pi**2, eq=eq_deg)
