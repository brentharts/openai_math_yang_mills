"""The Hopf fibration, from the companion repository brentharts/HopfFibration.

The repository (hopf_math.py, hopf_figures.py, HOPF_FIGURES.md) builds fibres, their stereographic
projections and linking numbers with numpy. This module does not copy its scripts (the repository is
GPL-2.0, this fork is Apache-2.0). It rederives, in exact symbolic algebra, the closed forms those
figures rely on:

  * the Hopf map h(z1, z2) = (2 Re z1 conj z2, 2 Im z1 conj z2, |z1|^2 - |z2|^2) has |h| = |z|^2,
    is constant on the fibres e^{i t}(z1, z2), and satisfies |<z, w>|^2 = (1 + h(z).h(w)) / 2,
    so the Clifford (great-circle) distance between two fibres is half the angle between their base points;
  * in Hopf coordinates z1 = sin(eta) e^{i xi1}, z2 = cos(eta) e^{i xi2} the pull-back of the area
    form of the unit S^2 is minus twice d(A) (the sign is the orientation convention of this h) with A = Im(conj z dz); then (1/4 pi^2) Int A ^ dA = 1, the
    Hopf invariant of the fibration;
  * the maps f_{p,q} = Hopf(z1^p, z2^q) normalised have Hopf invariant p q (the preimage of a point is the
    (q, p) torus knot and two preimages link p q times);
  * the half-plane {Im z1 = 0, Re z1 > 0} is a section: its stereographic image has closed form
    rho = s / (1 + c sin a), z = c cos a / (1 + c sin a) with s = sin(el/2), c = cos(el/2), and the
    image of each latitude circle is the circle of centre (1/s, 0) and radius c/s.

Not tested: the numerical figures themselves (they are the repository's, with its own test run), the
orientation conventions of the Blender scripts, any claim about physics. The Hopf fibration is a pure
topological fact; nothing here is a statement about Yang-Mills existence or a mass gap.
"""
import sympy as sp

from eqlib import Equation, check, register

NAME = "hopf_fibration"
ENTRY = "hopf-fibration"
REPO = "brentharts/HopfFibration"
DOC = f"{REPO}: HOPF_FIGURES.md"


def build_and_verify():
    x1, x2, x3, x4 = sp.symbols("x1 x2 x3 x4", real=True)
    y1, y2, y3, y4 = sp.symbols("y1 y2 y3 y4", real=True)
    t = sp.symbols("t", real=True)
    z1, z2 = x1 + sp.I * x2, x3 + sp.I * x4
    w1, w2 = y1 + sp.I * y2, y3 + sp.I * y4

    def hopf(a, b):
        ab = sp.expand(a * sp.conjugate(b))
        return (2 * sp.re(ab), 2 * sp.im(ab), sp.expand(a * sp.conjugate(a) - b * sp.conjugate(b)))

    h = hopf(z1, z2)
    n2 = sp.expand(x1**2 + x2**2 + x3**2 + x4**2)

    # ---- the map, its fibres, the Clifford distance -------------------------------------------
    eq_map = register(NAME, Equation(
        key="hf_map", entry=ENTRY, title="The Hopf map and its fibres",
        source_kind="textbook", source="Hopf, Math. Ann. 104 (1931); the convention of hopf_math.py in " + REPO,
        latex_str=(r"h(z_1,z_2)=\bigl(2\,\mathrm{Re}\,z_1\bar z_2,\;2\,\mathrm{Im}\,z_1\bar z_2,\;|z_1|^2-|z_2|^2\bigr),"
                   r"\quad |h|=|z|^2,\quad h(e^{it}z)=h(z)"),
        notes="Maps S^3 in C^2 onto S^2. The fibre over a point is a great circle e^{it}(z1, z2). Quantum-lattice's "
              "d-vector uses conj(z1) z2 instead, which is the reflection y -> -y of this convention.",
    ))
    mag = sp.simplify(sum(c**2 for c in h) - n2**2)
    check(NAME, "|h|^2 = (|z1|^2 + |z2|^2)^2 identically, so h maps the unit S^3 into the unit S^2", mag == 0, eq=eq_map)
    hz = hopf(sp.exp(sp.I * t) * z1, sp.exp(sp.I * t) * z2)
    check(NAME, "h is constant on the circle e^{it}(z1, z2)",
          all(sp.simplify(a - b) == 0 for a, b in zip(hz, h)), eq=eq_map)

    ip = sp.expand(sp.conjugate(z1) * w1 + sp.conjugate(z2) * w2)
    ip2 = sp.expand(ip * sp.conjugate(ip))
    hw = hopf(w1, w2)
    dot = sum(a * b for a, b in zip(h, hw))
    eq_cliff = register(NAME, Equation(
        key="hf_clifford", entry=ENTRY, title="Distance between two fibres is half the base angle",
        source_kind="derived", source=f"derived here; compare {DOC} figure 1 (linking atlas)",
        latex_str=(r"|\langle z,w\rangle|^2=\tfrac12\bigl(|z|^2|w|^2+h(z)\!\cdot\! h(w)\bigr)"
                   r"\;\Rightarrow\;\min_t d_{S^3}(z,e^{it}w)=\tfrac12\,\angle(h(z),h(w))"),
        notes="For unit z, w the angular distance of the fibres is arccos|<z,w>| = Theta/2 with Theta the angle "
              "between the base points, so distinct fibres are everywhere equidistant (Clifford parallel).",
    ))
    check(NAME, "|<z,w>|^2 = (|z|^2|w|^2 + h(z).h(w)) / 2 as polynomials in the eight real coordinates",
          sp.expand(2 * ip2 - (n2 * sp.expand(y1**2 + y2**2 + y3**2 + y4**2) + dot)) == 0, eq=eq_cliff)

    # ---- Hopf coordinates: the pull-back of the area form and the Hopf invariant ---------------
    eta, a1, a2 = sp.symbols("eta xi1 xi2", real=True)
    s, c = sp.sin(eta), sp.cos(eta)
    Z1, Z2 = s * sp.exp(sp.I * a1), c * sp.exp(sp.I * a2)
    H = sp.Matrix(hopf(Z1, Z2)).applyfunc(sp.simplify)
    co = (eta, a1, a2)
    dH = [H.diff(v) for v in co]

    def sigma(i, j):  # pull-back of the unit-sphere area form: h . (d_i h x d_j h)
        return sp.simplify(H.dot(dH[i].cross(dH[j])))

    # A = Im(conj z dz) = sin^2 dxi1 + cos^2 dxi2 ; dA = 2 sin cos d eta ^ (dxi1 - dxi2)
    A = sp.Matrix([0, s**2, c**2])
    dA = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            dA[i, j] = sp.diff(A[j], co[i]) - sp.diff(A[i], co[j])
    eq_inv = register(NAME, Equation(
        key="hf_invariant", entry=ENTRY, title="Hopf invariant 1 from the connection form",
        source_kind="derived", source=f"derived here (Hopf coordinates); compare {DOC} figure 4",
        latex_str=(r"h^*\sigma=-2\,dA,\quad A=\mathrm{Im}(\bar z\,dz),\qquad "
                   r"H=\frac{1}{4\pi^2}\int_{S^3}A\wedge dA=1"),
        notes="The pull-back of the area form sigma of the unit S^2 is -2 dA for this h (the sign is orientation), so "
              "(1/4 pi) h^* sigma = -dA / 2 pi has primitive -A / 2 pi and the Hopf invariant is (1/4 pi^2) Int A ^ dA "
              "(the sign of the primitive cancels in the product). The integral is exact below.",
    ))
    check(NAME, "h^* sigma = -2 dA on S^3, component by component, in Hopf coordinates",
          all(sp.simplify(sigma(i, j) + 2 * dA[i, j]) == 0 for i in range(3) for j in range(3)), eq=eq_inv)
    # A ^ dA = (A_i dA_jk) coefficient of d eta d xi1 d xi2
    top = sp.simplify(A[0] * dA[1, 2] - A[1] * dA[0, 2] + A[2] * dA[0, 1])
    integral = sp.integrate(top, (eta, 0, sp.pi / 2)) * (2 * sp.pi) ** 2
    check(NAME, "(1/4 pi^2) Int A ^ dA over eta in [0, pi/2], xi1, xi2 in [0, 2 pi] equals 1 in absolute value",
          sp.simplify(integral / (4 * sp.pi**2)) in (1, -1), eq=eq_inv)

    # ---- f_{p,q}: Hopf invariant p q ----------------------------------------------------------
    p, q = sp.symbols("p q", positive=True, integer=True)
    S = sp.Function("S")          # S(eta) = |z1'|^2 / |z'|^2 of the normalised (z1^p, z2^q)
    Ap = sp.Matrix([0, p * S(eta), q * (1 - S(eta))])
    dAp = sp.zeros(3, 3)
    for i in range(3):
        for j in range(3):
            dAp[i, j] = sp.diff(Ap[j], co[i]) - sp.diff(Ap[i], co[j])
    topp = sp.simplify(Ap[0] * dAp[1, 2] - Ap[1] * dAp[0, 2] + Ap[2] * dAp[0, 1])
    # the integrand is p q S'(eta) up to sign, whose integral over [0, pi/2] is p q (S(pi/2) - S(0)) = p q
    eq_pq = register(NAME, Equation(
        key="hf_pq", entry=ENTRY, title="The maps f_{p,q} have Hopf invariant p q",
        source_kind="derived", source=f"derived here; compare {DOC} figure 4 (torus knots, linking p q)",
        latex_str=(r"f_{p,q}(z_1,z_2)=h\!\left(\tfrac{z_1^p}{\sqrt{|z_1|^{2p}+|z_2|^{2q}}},"
                   r"\tfrac{z_2^q}{\sqrt{|z_1|^{2p}+|z_2|^{2q}}}\right),\qquad H(f_{p,q})=pq"),
        notes="Here S(eta) is the squared weight of the first coordinate after normalisation; it rises from 0 to 1, "
              "which is all the integral uses. Preimages of a point are (q, p) torus knots and two of them link p q "
              "times (checked numerically, not here).",
    ))
    check(NAME, "the integrand of A ^ dA for f_{p,q} is p q S'(eta) up to sign",
          sp.simplify(sp.Abs(topp) - p * q * sp.Abs(sp.diff(S(eta), eta))) == 0
          or sp.simplify(topp - p * q * sp.diff(S(eta), eta)) == 0
          or sp.simplify(topp + p * q * sp.diff(S(eta), eta)) == 0, eq=eq_pq)
    check(NAME, "for p = q = 1 and S = sin^2 the formula reduces to the Hopf invariant 1 above",
          sp.simplify(sp.integrate((topp.subs({p: 1, q: 1}).subs(S(eta), s**2)).doit(), (eta, 0, sp.pi / 2))) in (1, -1),
          eq=eq_pq)

    # ---- the half-plane section and its closed form -------------------------------------------
    el, az = sp.symbols("el az", real=True)
    sh, ch = sp.sin(el / 2), sp.cos(el / 2)
    # z1 = sh (real), z2 = ch e^{i phi} with phi = -az; stereographic projection from (0,0,0,1)
    phi = -az
    P = sp.Matrix([sh, 0, ch * sp.cos(phi), ch * sp.sin(phi)])
    proj = sp.Matrix(P[:3]) / (1 - P[3])
    rho = sp.sqrt(proj[0]**2 + proj[1]**2)
    zz = proj[2]
    eq_sec = register(NAME, Equation(
        key="hf_section", entry=ENTRY, title="A half-plane section and its stereographic image",
        source_kind="derived", source=f"derived here; compare {DOC} figure 2 (Poincare half-plane section)",
        latex_str=(r"\rho=\frac{s}{1+c\sin a},\quad z=\frac{c\cos a}{1+c\sin a},\quad s=\sin\tfrac{el}{2},\ c=\cos\tfrac{el}{2};"
                   r"\qquad (\rho-1/s)^2+z^2=c^2/s^2"),
        notes="The set {Im z1 = 0, Re z1 > 0} meets every fibre once, so it is a global section over the sphere minus "
              "a point. At fixed elevation the image is a circle of centre (1/s, 0) and radius c/s, a meridian circle "
              "of a torus in the stereographic picture.",
    ))
    # the base point of P: h2 = -2 s c sin(phi) = 2 s c sin(az), so the azimuth is az; el from h3 = s^2 - c^2 = -cos el
    hP = hopf(P[0] + sp.I * P[1], P[2] + sp.I * P[3])
    check(NAME, "the section point has base point (sin el cos az, sin el sin az, -cos el)",
          all(sp.simplify(sp.expand_trig(a - b)) == 0 for a, b in zip(
              hP, (sp.sin(el) * sp.cos(az), sp.sin(el) * sp.sin(az), -sp.cos(el)))), eq=eq_sec)
    rho_closed = sh / (1 + ch * sp.sin(az))
    z_closed = ch * sp.cos(az) / (1 + ch * sp.sin(az))
    check(NAME, "the stereographic image is rho = s/(1 + c sin a), z = c cos a/(1 + c sin a)",
          sp.simplify(proj[0] - rho_closed) == 0 and sp.simplify(proj[1]) == 0 and sp.simplify(zz - z_closed) == 0,
          eq=eq_sec)
    check(NAME, "(rho - 1/s)^2 + z^2 = c^2/s^2: the image of each latitude is a circle",
          sp.simplify(sp.expand_trig((rho_closed - 1 / sh) ** 2 + z_closed**2 - ch**2 / sh**2)) == 0, eq=eq_sec)
