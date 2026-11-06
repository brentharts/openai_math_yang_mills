"""The quantum gravity paper: spin networks on the SU(2) ladder, large generated equations, no tables.

Generates build/qgravity_body.tex, build/qgravity_bibliography.tex and build/qgravity_numbers.tex, which the hand-editable wrapper
qgravity_document.tex pulls in. It was made by copying concise.py and stripping it down: the machinery that builds the large
generated displays (Display, family_row, geometry_row, arrow_row, raw_row, the figures' \resizebox wrapper, esc_title, selected_section)
is kept as it was, the figures that are about tilings, Lie algebras and the Hopf and gauge material are gone, and the quantum gravity
content is new. Everything shown is built from the knowledge graph (knowledgebase.py) and from the registered sympy models, in the
manner of Figures 1 and 2 of RosettaMath's rosettapaper.py:

    \\overbrace{quantity}^{name} = \\underbrace{expression}_{equation} = \\underbrace{...}

Each term of a row is the library equation *solved for the row's quantity*, and the build substitutes that solution back into the
equation and demands an identity (verify_solution). A failed check raises, so nothing unverified is written.

To extend the paper: add equations to the knowledge graph (EXTRA_* in knowledgebase.py), add a row to QG_ROWS (or a new *_ROWS tuple and a
figure function like quantum_geometry_figure), register results in models/qg_spin_network.py, and add a section in write().
"""
from __future__ import annotations

import re
from pathlib import Path

import sympy as sp

import eqlib
import knowledgebase as kb
from document import esc, cite, render_see


class QgError(AssertionError):
    pass


# ---------------------------------------------------------------------------------------
# short names for the underbraces (RosettaMath keeps the same kind of table)
# ---------------------------------------------------------------------------------------
SHORT = {
    "Mass-energy equivalence": "mass-energy",
    "Planck relation": "Planck",
    "Kinetic energy": "kinetic",
    "Thermal energy": "thermal",
    "Mass gap in physical units": "mass gap",
    "Coulomb potential energy": "Coulomb",
    "Newton's second law": "Newton II",
    "Newton's law of gravitation": "Newton",
    "Coulomb's law": "Coulomb",
    "Hooke's law": "Hooke",
    "Schwarzschild radius": "Schwarzschild",
    "Bekenstein-Hawking entropy": "Bekenstein--Hawking",
    "Planck length": "Planck length",
    "Escape velocity": "escape",
    "Boltzmann entropy": "Boltzmann",
    "de Broglie relation": "de Broglie",
    "Compton wavelength": "Compton",
    "Euclidean correlator decay": "correlator",
    "Asymptotic-freedom mass scaling": "asymptotic freedom",
    "Spectre inflation factor": "Spectre",
    "Spectre reciprocal pair": "reciprocal pair",
    "Spectre linear-area relation": "linear--area",
    "Spectre growth rate": "growth rate",
    "Hat inflation factor": "hat",
    "Hat constant as a golden power": "golden power",
    "Spin-model transmutation slope": "O($N$) slope",
    "Spin block coupling": "spin block",
    "Gauge transmutation slope": "gauge slope",
    "Gauge block coupling": "gauge block",
    "Lattice coupling": "lattice coupling",
    "Golden ratio": "golden ratio",
    "Penrose tile ratio": "Penrose tiles",
    "Half-spinor dimension": "half-spinor",
    "E6 fundamental branching": "$27=1+16+10$",
    "Chirality splitting": "chirality splitting",
    "Net chiral class": "net chiral class",
    "Weight vector demand": "weight demand",
    "Standard Model content of the 27": "$27=15+10+2$",
    "Hypercharge": "hypercharge",
    "Gell-Mann-Nishijima formula": "Gell-Mann--Nishijima",
    "Orthogonal algebra dimension": "$\\mathfrak{so}(n)$",
    "Special unitary dimension": "$\\mathfrak{su}(N)$",
    "Product representation dimension": "family module",
    "Hypercharge normalisation": "traces",
    "Weak mixing at unification": "$\\sin^2\\theta_W$",
}

SHORT.update({
    "Spin Casimir": "Casimir",
    "LQG area spectrum": "LQG area",
    "Area of N spin-1/2 punctures": "$N$ punctures",
    "Horizon area": "horizon",
    "Spin-1/2 puncture microstates": "$2^N$ states",
    "Euler characteristic of an Einstein four-manifold": "Gauss--Bonnet",
    "Signature of an Einstein four-manifold": "signature theorem",
    "Minimal sectional curvature": "least plane",
    "Weyl block Casimir": "Weyl block, $j=2$",
    "Penrose mass at equality": "Penrose",
})

# Figure 1: quantities named rather than computed, for the reason RosettaMath names its pivots:
# "every family" does not fit on a page, and some of them (the two inflation factors) only
# look alike. Each pivot lists the equations whose solution for it is shown.
PHYSICS_ROWS = (
    ("energy", ("Mass-energy equivalence", "Planck relation", "Kinetic energy", "Thermal energy",
                "Mass gap in physical units")),
    ("force", ("Newton's second law", "Newton's law of gravitation", "Coulomb's law", "Hooke's law")),
    ("gravitational constant", ("Newton's law of gravitation", "Schwarzschild radius",
                                "Bekenstein-Hawking entropy", "Planck length", "Escape velocity")),
    ("boltzmann constant", ("Boltzmann entropy", "Bekenstein-Hawking entropy", "Thermal energy")),
    ("reduced planck constant", ("Bekenstein-Hawking entropy", "Planck length",
                                 "Mass gap in physical units")),
    ("planck constant", ("Planck relation", "de Broglie relation", "Compton wavelength")),
    ("lattice mass", ("Euclidean correlator decay", "Mass gap in physical units")),
)

# Figure 3: the transmutation rows come from the library; the bridge between the O(N) model
# and Yang-Mills is checked separately in transmutation_bridge().
TRANSMUTATION_ROWS = (
    ("spin transmutation slope", ("Spin-model transmutation slope", "Spin block coupling")),
    ("gauge transmutation slope", ("Gauge transmutation slope", "Gauge block coupling")),
    ("inverse coupling", ("Lattice coupling", "Gauge block coupling")),
    ("lattice mass", ("Euclidean correlator decay", "Mass gap in physical units")),
)

# Pairs that must never share an identity row: the O(N) spin model's scaling form is a statement
# about a different theory from the gauge-theory lattice mass and coupling, and the library
# relates them by a link (a correspondence), not by an equals sign (rosettaphys joins them to
# nothing; see its selftest).
NOT_IDENTICAL = (
    ("Asymptotic-freedom mass scaling", "Euclidean correlator decay"),
    ("Asymptotic-freedom mass scaling", "Mass gap in physical units"),
    ("Asymptotic-freedom mass scaling", "Lattice coupling"),
    ("Spin-model transmutation slope", "Gauge transmutation slope"),
    ("Spin block coupling", "Gauge block coupling"),
)

# Figure 3's correspondences, drawn from the library's links (kind and note come from there).
TRANSMUTATION_ARROWS = (
    ("Asymptotic-freedom mass scaling", "Dimensional transmutation"),
    ("Euclidean correlator decay", "Asymptotic-freedom mass scaling"),
    ("Spin-model transmutation slope", "Gauge transmutation slope"),
    ("Spin block coupling", "Gauge block coupling"),
)


# Figure 2: quantum geometry. Each row is one quantity stated by several library equations: the area of a surface, the Casimir that
# the area is a function of, the number of microstates and the entropy. The first three equations are new in the knowledge graph
# (Spin Casimir, LQG area spectrum, Area of N spin-1/2 punctures, Horizon area, Spin-1/2 puncture microstates); the entropy and the
# Planck length are the library's own.
QG_ROWS = (
    ("area", ("LQG area spectrum", "Area of N spin-1/2 punctures", "Horizon area", "Bekenstein-Hawking entropy")),
    ("casimir eigenvalue", ("Spin Casimir", "LQG area spectrum")),
    ("microstate count", ("Spin-1/2 puncture microstates", "Boltzmann entropy")),
    ("entropy", ("Bekenstein-Hawking entropy", "Boltzmann entropy")),
)
# links of the library, drawn as correspondences and never as equations
QG_ARROWS = (
    ("Spin Casimir", "LQG area spectrum"),
    ("Area of N spin-1/2 punctures", "LQG area spectrum"),
    ("LQG area spectrum", "Bekenstein-Hawking entropy"),
    ("Spin-1/2 puncture microstates", "Boltzmann entropy"),
)

# Figure 3: Einstein four-manifolds (entry 348). The two characteristic identities state the volume and the two Weyl second moments in two
# ways each; the three model spaces are then rows whose two terms are equal numbers once chi, tau and V take the model's values.
EINSTEIN_ROWS = (
    ("normalised volume", ("Euler characteristic of an Einstein four-manifold", "Signature of an Einstein four-manifold")),
    ("self-dual Weyl moment", ("Euler characteristic of an Einstein four-manifold", "Signature of an Einstein four-manifold")),
    ("anti-self-dual Weyl moment", ("Euler characteristic of an Einstein four-manifold", "Signature of an Einstein four-manifold")),
)
EINSTEIN_SPACES = (   # label, assignment of the letters other than the one solved for
    ("round $S^{4}$", {"chi_E": 2, "tau_E": 0, "V": 8 * sp.pi ** 2 / 3, "s_A": 0}),
    (r"Fubini--Study $\mathrm{CP}^{2}$", {"chi_E": 3, "tau_E": 1, "V": 2 * sp.pi ** 2, "s_A": 0}),
    (r"$S^{2}\!\times\! S^{2}$", {"chi_E": 4, "tau_E": 0, "V": 16 * sp.pi ** 2 / 9, "s_A": 1}),
)
EINSTEIN_ARROWS = (
    ("Euler characteristic of an Einstein four-manifold", "Signature of an Einstein four-manifold"),
    ("Spin Casimir", "Weyl block Casimir"),
    ("Weyl block Casimir", "Signature of an Einstein four-manifold"),
)

# registered results shown in the appendix (the ladder's own equations, from su2_lattice_ladder.py); the quantum gravity results
# of models/qg_spin_network.py are displayed in the sections of the paper
SELECTED = (
    ("The ladder: basis, Hamiltonian and plaquette operator", "su2-lattice-gap", ("su2_ladder_basis",),
     "The seven links and two trivalent vertices of the theta graph, the Kogut--Susskind Hamiltonian in lattice units and the matrix "
     "elements of the plaquette operators, from the model \\texttt{su2\\_lattice\\_ladder}. These are the equations on which the "
     "geometry of Section~\\ref{sec:ladder} rests."),
    ("The ladder: weak and strong coupling, and the gap", "su2-lattice-gap",
     ("su2_ladder_weak", "su2_ladder_strong", "su2_ladder_gap"),
     "The weak-coupling limit of the spectrum from two coupled oscillators (the harmonic argument that also gives the weak-coupling "
     "geometry of Section~\\ref{sec:ladder}), the exact strong-coupling splitting, and the gap against the one-plaquette gap. "
     "A spectral gap of a lattice gauge theory on two plaquettes: it says nothing about a gap of any field theory."),
)

# every registered result the paper needs; if one is missing (a partial run of export.py) the paper is skipped, not half-built
QG_MODELS = ("qg_spin_network", "qg_einstein4", "su2_lattice_strip", "su2_peter_weyl", "su2_lattice_stripn", "qg_volume", "qg_penrose", "su2_lattice_lshape", "forest_gauge_4d_theorems", "qg_kerr", "qg_planar_maps")   # the models whose checks the paper counts as its own
REQUIRED = ("su2_ladder_basis", "su2_ladder_weak", "su2_ladder_strong", "su2_ladder_gap",
            "qg_area_spectrum", "qg_casimir_dictionary", "qg_intertwiners", "qg_entropy_immirzi",
            "qg_quantum_tetrahedron", "qg_ladder_geometry",
            "qg_so4_casimirs", "qg_einstein_blocks", "qg_einstein_topology",
            "su2_strip_basis", "su2_stripn_basis", "su2_strip_weak", "su2_strip_strong", "su2_strip_gap",
            "qg_volume_triple_product", "qg_volume_spectra", "qg_volume_semiclassical", "qg_volume_l_shape",
            "qg_penrose_equality", "qg_penrose_quantum",
            "su2_lshape_basis", "su2_lshape_volume", "su2_lshape_gap",
            "fg4t_linefield", "fg4t_partners", "fg4t_table1", "fg4t_theorems",
            "qg_kerr_horizons", "qg_kerr_penrose", "qg_kerr_quantum", "qg_maps_counts", "qg_maps_exponents", "qg_maps_su2k")

# ---------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------
def short(name: str) -> str:
    return SHORT.get(name, name.lower())


def _tex_text(s: str) -> str:
    """Text for \\text{...}: SHORT entries may already carry $...$ math."""
    return s.replace("&", r"\&").replace("%", r"\%").replace("_", r"\_").replace("#", r"\#")


def _compact(tex: str) -> str:
    return tex.replace(r"\frac", r"\tfrac").replace(r"\dfrac", r"\tfrac")


def _positive(expr):
    return expr.xreplace({s: sp.Symbol(s.name, positive=True) for s in expr.free_symbols})


def verify_solution(eq: kb.KBEquation, sol: kb.Solution) -> str:
    """Substituting the solution into the equation must leave an identity."""
    rel = eq.relation()
    residual = sp.simplify(_positive((rel.lhs - rel.rhs).subs(sol.symbol, sol.expr)))
    if residual != 0:
        raise QgError(f"{eq.name}: solving for {sol.quantity} leaves residual {residual}")
    return f"{eq.name}, solved for {sol.quantity}, satisfies the equation identically"


def _term(eq: kb.KBEquation, sol: kb.Solution) -> str:
    return rf"\underbrace{{{_compact(sp.latex(sol.expr))}}}_{{\text{{\tiny {_tex_text(short(eq.name))}}}}}"


def _head(sol: kb.Solution, quantity: str) -> str:
    return rf"\overbrace{{{sp.latex(sol.symbol)}}}^{{\text{{\tiny {_tex_text(quantity)}}}}}"


class Display:
    """Rows and their provenance: the equations and quantities a display brings together."""

    def __init__(self):
        self.rows: list[str] = []
        self.equations: set[str] = set()
        self.quantities: set[str] = set()
        self.checks: list[str] = []
        self.arrows = 0
        self.names: list[str] = []          # in order of first appearance, for the prose line

    def _note(self, eq: kb.KBEquation):
        if eq.name not in self.equations:
            self.names.append(eq.name)
        self.equations.add(eq.name)
        self.quantities.update(eq.quantities())

    def family_row(self, quantity: str, names, tail: str = ""):
        for a, b in NOT_IDENTICAL:
            if a in names and b in names:
                raise QgError(f"{a} and {b} are related by a correspondence, not an identity; "
                                   f"they cannot share the row for {quantity}")
        terms, head = [], None
        for name in names:
            eq = kb.equation(name)
            sol = kb.solve_for(eq, quantity)
            if sol is None:
                raise QgError(f"{name} cannot be solved for {quantity}")
            self.checks.append(verify_solution(eq, sol))
            self._note(eq)
            head = head or _head(sol, quantity)
            terms.append(_term(eq, sol))
        if len({str(t) for t in terms}) != len(terms):
            raise QgError(f"row {quantity} repeats a term")
        self.rows.append(f"{head} &= {' = '.join(terms)}{tail}")

    def geometry_row(self, quantity: str, letter: str, names, assign: dict, label: str):
        """A row whose terms are equal numbers once the other letters take the stated values."""
        terms, head, values = [], None, []
        for name in names:
            eq = kb.equation(name)
            sol = kb.solve_for(eq, quantity)
            if sol is None or sol.symbol.name != letter:
                raise QgError(f"{name}: expected {quantity} to be denoted {letter}")
            self.checks.append(verify_solution(eq, sol))
            self._note(eq)
            sub = {s: assign[s.name] for s in sol.expr.free_symbols if s.name in assign}
            missing = {s.name for s in sol.expr.free_symbols} - {s.name for s in sub}
            if missing:
                raise QgError(f"{name}: no value for {sorted(missing)}")
            values.append(sp.simplify(sp.radsimp(sol.expr.subs(sub))))
            head = head or _head(sol, f"{quantity}, {label}")
            terms.append(_term(eq, sol))
        first = values[0]
        for name, v in zip(names, values):
            if sp.simplify(sp.expand_log(v - first, force=True)) != 0 and abs(sp.N(v - first, 40)) > sp.Float("1e-30"):
                raise QgError(f"{label}: {name} gives {v}, not {first}")
        self.checks.append(f"{label}: all {len(names)} terms equal {first} under "
                           + ", ".join(f"{k} = {v}" for k, v in assign.items()))
        self.rows.append(f"{head} &= {' = '.join(terms)} = {sp.latex(first)}")

    def arrow_row(self, a: str, b: str):
        """A correspondence: two library equations joined by the arrow their link names."""
        link = next(((k, n) for x, k, y, n in kb.library().links if (x, y) == (a, b)), None)
        if link is None:
            raise QgError(f"the library has no link from {a} to {b}")
        kind, _note = link
        arrow, word = {"FIELD": (r"\longleftrightarrow", "same field"),
                       "SPECIALISES": (r"\longrightarrow", "specialises"),
                       "DEFINES": (r"\longrightarrow", "defines")}[kind]
        side = lambda n: (rf"\underbrace{{{_compact(kb.equation(n).latex)}}}"
                          rf"_{{\text{{\tiny {_tex_text(short(n))}}}}}")
        for n in (a, b):
            self._note(kb.equation(n))
        self.arrows += 1
        self.checks.append(f"{a} -- {kind} -> {b}: a link of the library, drawn as a correspondence, not as an equation")
        self.rows.append(rf"{{}} & {side(a)} \;\overset{{\text{{\tiny {word}}}}}{{{arrow}}}\; {side(b)}")

    def raw_row(self, tex: str, checks=()):
        self.rows.append(tex)
        self.checks.extend(checks)

    def tex(self) -> str:
        ids = [r for r in self.rows if not r.startswith("{}")]
        arrows = [r for r in self.rows if r.startswith("{}")]
        body = " \\\\[1.2ex]\n".join(ids)
        if arrows:
            body += " \\\\[2.4ex]\n" + " \\\\[1.6ex]\n".join(arrows)
        tally = (r"\begin{array}{l} \text{identities:} \\ "
                 rf"\text{{{len(self.equations)} equations, {len(ids)} quantities}}")
        if arrows:
            tally += rf" \\[0.6ex] \text{{correspondences:}} \\ \text{{{len(arrows)} links}}"
        tally += r" \end{array}"
        return f"\\left.\n\\begin{{aligned}}\n{body}\n\\end{{aligned}}\n\\;\\right\\}} {tally}"


def _figure(label: str, display: str, caption: str, width: float = 0.98) -> str:
    return (f"\\begin{{figure}}[p]\n\\centering\n"
            f"\\resizebox{{{width}\\textwidth}}{{!}}{{$\\displaystyle\n{display}\n$}}\n"
            f"\\caption{{{caption}}}\n\\label{{{label}}}\n\\end{{figure}}\n")


def _names_line(names) -> str:
    return ", ".join(rf"\href{{{kb.wikipedia_url(kb.equation(n).slug)}}}{{{esc(n)}}}"
                     if kb.wikipedia_url(kb.equation(n).slug) else esc(n) for n in names)

# ---------------------------------------------------------------------------------------
# the three figures
# ---------------------------------------------------------------------------------------
def physics_figure():
    d = Display()
    for quantity, names in PHYSICS_ROWS:
        d.family_row(quantity, names)
    cap = ("The physics of the library, conjoined. Each row is one quantity, overbraced with its name; "
           "each term is a library equation solved for it, underbraced with the equation's name. "
           "The letters are those of the equations, so an equals sign says that each equation "
           "states the same quantity, not that its other letters take one common value.")
    return d, _figure("fig:physics", d.tex(), cap, 0.98)


def quantum_geometry_figure():
    d = Display()
    for quantity, names in QG_ROWS:
        d.family_row(quantity, names)
    # two exact numbers, rederived here from the library's own equations before they are written
    area_lqg = kb.solve_for(kb.equation("LQG area spectrum"), "area").expr
    area_n = kb.solve_for(kb.equation("Area of N spin-1/2 punctures"), "area").expr
    Nn, gam, c2 = sp.Symbol("N"), sp.Symbol("gamma_I"), sp.Symbol("c_2")
    if sp.simplify(area_lqg.subs(c2, sp.Rational(3, 4)) - area_n.subs(Nn, 1)) != 0:
        raise QgError("the area spectrum at c_2 = 3/4 differs from the area of one spin-1/2 puncture")
    d.raw_row(r"\overbrace{A_{\min}}^{\text{\tiny one puncture}} &= "
              r"\underbrace{8\pi\gamma_{I}\ell_P^{2}\sqrt{\tfrac34}}_{\text{\tiny LQG area},\ c_2=3/4} = "
              r"\underbrace{4\sqrt3\,\pi\gamma_{I}\ell_P^{2}}_{\text{\tiny $N=1$ puncture}}",
              ["the area spectrum at c_2 = 3/4 equals the area of N = 1 spin-1/2 puncture, 4 sqrt3 pi gamma l_P^2"])
    bh = kb.equation("Bekenstein-Hawking entropy").relation()
    pl = kb.equation("Planck length").relation()
    s_over_kb = sp.simplify((bh.rhs / sp.Symbol("k_B")).subs(sp.Symbol("A"), area_n).subs(sp.Symbol("ell_P"), pl.rhs))
    if sp.simplify(s_over_kb - sp.sqrt(3) * sp.pi * gam * Nn) != 0:
        raise QgError("S/k_B for N spin-1/2 punctures is not sqrt3 pi gamma N")
    sol = sp.solve(sp.Eq(s_over_kb, Nn * sp.log(2)), gam)
    if len(sol) != 1 or sp.simplify(sol[0] - sp.log(2) / (sp.pi * sp.sqrt(3))) != 0 or abs(float(sol[0]) - 0.127384) > 1e-6:
        raise QgError(f"the Immirzi value from the library's entropy is {sol}, not ln2/(pi sqrt3) = 0.127384")
    d.raw_row(r"\overbrace{\gamma_{I}}^{\text{\tiny spin-1/2 count only}} &= "
              r"\underbrace{\frac{\ln 2}{\pi\sqrt3}}_{\text{\tiny $S/k_B=\sqrt3\pi\gamma_{I}N=N\ln2$}} \approx 0.1274",
              ["gamma = ln2/(pi sqrt3) = 0.127384 from the library's entropy and Planck length with the area of N spin-1/2 punctures "
               "(G, hbar and c cancel)"])
    for a, b in QG_ARROWS:
        d.arrow_row(a, b)
    cap = ("Quantum geometry in the library, conjoined. Each identity row is one quantity (the area of a surface, the Casimir it is a "
           "function of, the number of microstates, the entropy) stated by several library equations, each solved for it. The equals sign "
           "says that each equation states the same quantity, not that the other letters take a common value; in the microstate row the "
           "$2^N$ is the leading exponential of the exact count, which is a Catalan number (Section~\\ref{sec:intertwiners}). The two rows "
           "below are exact numbers, rederived from the library's own entropy and Planck length: $G$, $\\hbar$ and $c$ cancel, leaving "
           "$S/k_B=\\sqrt3\\pi\\gamma_IN$. The arrows are links of the library, not equations: the Casimir defines the area, the "
           "$N$-puncture area specialises it, and the area spectrum and the entropy state the same area only once a counting fixes "
           "$\\gamma_I$, which depends on which spins are counted (Section~\\ref{sec:entropy}).")
    return d, _figure("fig:qgeometry", d.tex(), cap, 0.98)

def einstein_figure():
    d = Display()
    for quantity, names in EINSTEIN_ROWS[:1]:
        d.family_row(quantity, names)
    for label, assign in EINSTEIN_SPACES:
        d.geometry_row("self-dual Weyl moment", "s_S", EINSTEIN_ROWS[1][1], assign, label)
    # the least sectional curvature of each space, from x and y (the depths of the two blocks, computed in qg_einstein4)
    eq = kb.equation("Minimal sectional curvature")
    for label, xy, kmin in ((r"$S^{4}$", (0, 0), 1), (r"$\mathrm{CP}^{2}$", (1, 0), sp.Rational(1, 2)), (r"$S^{2}\!\times\! S^{2}$", (1, 1), 0)):
        sol = kb.solve_for(eq, "minimal sectional curvature")
        d.checks.append(verify_solution(eq, sol))
        d._note(eq)
        val = sp.simplify(sol.expr.subs({sp.Symbol("x"): xy[0], sp.Symbol("y"): xy[1]}))
        if val != kmin:
            raise QgError(f"least sectional curvature of {label}: {val}, expected {kmin}")
        d.checks.append(f"{label}: x = {xy[0]}, y = {xy[1]} give K_min = {val} (the depths x, y are computed in qg_einstein4)")
        d.raw_row(rf"\overbrace{{K_{{m}}}}^{{\text{{\tiny least sectional curvature, {label}}}}} &= "
                  rf"\underbrace{{1-\tfrac{{{xy[0]}+{xy[1]}}}{{2}}}}_{{\text{{\tiny least plane}}}} = {sp.latex(val)}")
    for a, b in EINSTEIN_ARROWS:
        d.arrow_row(a, b)
    cap = ("Einstein four-manifolds in the library, conjoined (entry 348, $\\mathrm{Ric}=3g$). The first row states the volume twice, through the Euler "
           "characteristic and through the signature, in terms of the two Weyl second moments $s_S=\\mathbb E\\,v^2$ and $s_A=\\mathbb E\\,b^2$. "
           "The next three rows state $s_S$ both ways for the three model spaces, with $\\chi$, $\\tau$, $V$ and $s_A$ as computed from their metrics in "
           "Section~\\ref{sec:einstein}: the two terms of each row are the same number, shown at the end of the row ($0$, $1$, $1$), which is the check, since the equations "
           "alone would let any $s_S$ through. The last three rows give the least sectional curvature $1-(x+y)/2$ of the same spaces, which is "
           "$0$ for $S^2\\times S^2$, the planes of zero curvature that strict positivity excludes. The arrows are links of the library: "
           "$\\chi$ and $\\tau$ share the same two moments, and the Weyl block of one chirality is the spin-2 representation, whose Casimir "
           "is the $j(j+1)$ of a spin-network link at $j=2$.")
    return d, _figure("fig:einstein", d.tex(), cap, 0.98)


def penrose_figure():
    d = Display()
    d.family_row("mass", ("Schwarzschild radius", "Penrose mass at equality"))
    d.family_row("area", ("Horizon area", "Penrose mass at equality", "LQG area spectrum"))
    gam = sp.log(2) / (sp.pi * sp.sqrt(3))
    # the smallest quantised area through the Penrose mass, in Planck units: M/M_P = sqrt(A/(16 pi l_P^2))
    area = kb.solve_for(kb.equation("Area of N spin-1/2 punctures"), "area").expr
    val = sp.simplify(sp.sqrt(area.subs({sp.Symbol("N"): 1, sp.Symbol("gamma_I"): gam}) / (16 * sp.pi * sp.Symbol("ell_P") ** 2)))
    if sp.simplify(val - sp.sqrt(sp.log(2) / (4 * sp.pi))) != 0 or abs(float(val) - 0.2349) > 5e-5:
        raise QgError(f"Penrose bound at one spin-1/2 puncture: {val}")
    d.raw_row(r"\overbrace{\tfrac{M_{\min}}{M_P}}^{\text{\tiny one puncture, bound}} &= "
              r"\underbrace{\sqrt{\tfrac{A_{\min}}{16\pi\ell_P^{2}}}}_{\text{\tiny Penrose}} = "
              r"\underbrace{\sqrt{\tfrac{4\sqrt3\pi\gamma_I}{16\pi}}}_{\text{\tiny $N=1$ puncture}} = "
              r"\sqrt{\tfrac{\ln 2}{4\pi}}\approx 0.2349",
              ["the Penrose mass at the area of one spin-1/2 puncture with gamma = ln2/(pi sqrt3) is sqrt(ln2/(4 pi)) M_P = 0.2349 M_P (exact, then numerical)"])
    for a, b in (("Schwarzschild radius", "Penrose mass at equality"), ("LQG area spectrum", "Penrose mass at equality")):
        d.arrow_row(a, b)
    cap = ("The Penrose inequality of entry 260 at equality, in the library. The first row states the mass twice: through the Schwarzschild radius, and as the "
           "mass that saturates the inequality for a given area. The second states the area three ways: the horizon area, the area at which the Penrose mass "
           "equals $M$, and the quantised area of Section~\\ref{sec:geometry}; as in the other figures the equals sign says that each equation states the same "
           "quantity, not that the other letters share a value. The third row evaluates the bound at the area of one spin-$\\tfrac12$ puncture, with "
           "$\\gamma_I=\\ln2/(\\pi\\sqrt3)$, exactly. The arrows are library links: Schwarzschild data saturate the bound, and the area of the bound is the area that has a "
           "spectrum. Putting a quantum area eigenvalue where the theorem has a classical enclosing area is a heuristic (Section~\\ref{sec:penrose}).")
    return d, _figure("fig:penrose", d.tex(), cap, 0.98)


def transmutation_bridge() -> list[str]:
    """Check that the O(N) slope and the Yang-Mills lattice scale correspond.

    With beta = 2N/g^2 (lattice coupling) and b_0 = 11N/(48 pi^2), the Yang-Mills exponent is
    1/(2 b_0 g^2) = r_g beta with r_g = 12 pi^2/(11 N^2), the gauge transmutation slope; the
    O(N) model's exponent is r_s beta with r_s = 2 pi/(N-2). The two are slopes of the same
    kind, which the display joins by a correspondence arrow rather than an equals sign."""
    P = lambda e: _positive(e)
    N, beta, g = sp.symbols("N beta g", positive=True)
    lc = kb.solve_for(kb.equation("Lattice coupling"), "inverse coupling")
    r_g = kb.solve_for(kb.equation("Gauge transmutation slope"), "gauge transmutation slope")
    r_s = kb.solve_for(kb.equation("Spin-model transmutation slope"), "spin transmutation slope")
    b0 = eqlib.REGISTRY["ym_b0"].expr.rhs
    b0 = P(b0.xreplace({x: sp.Symbol("N") for x in b0.free_symbols if x.name == "N_c"}))
    g2 = sp.solve(sp.Eq(beta, P(lc.expr)), g)       # beta = 2N/g^2, so g^2 = 2N/beta
    g2 = [x for x in g2 if x.is_positive]
    if len(g2) != 1:
        raise QgError("cannot solve the lattice coupling for g")
    exponent = (1 / (2 * b0 * g ** 2)).subs(g, g2[0])
    ym = sp.simplify(exponent - P(r_g.expr) * beta)
    if ym != 0:
        raise QgError(f"Yang-Mills exponent differs from r_g*beta by {ym}")
    return ["Yang-Mills exponent 1/(2 b_0 g^2), with beta = 2N/g^2, equals the gauge transmutation slope times beta",
            f"the O(N) exponent uses the spin slope {r_s.expr}, a slope of the same kind with 1/(N-2) in place of 12 pi^2/(11 N^2)"]


def transmutation_figure():
    d = Display()
    for quantity, names in TRANSMUTATION_ROWS:
        d.family_row(quantity, names)
    d.checks.extend(transmutation_bridge())
    for a, b in TRANSMUTATION_ARROWS:
        d.arrow_row(a, b)
    cap = ("Dimensional transmutation across the library. The identity rows are the library's equations "
           "solved for the slopes, the coupling and the gauge-theory lattice mass. The arrows are "
           "correspondences, not identities: the $O(N)$ spin model's scaling form and the Yang--Mills lattice "
           "mass share a shape (a scale that is an exponential of an inverse coupling over a one-loop "
           "coefficient) but not coefficients, couplings or theory, so the library links them and the "
           "display does not equate them. The one exact statement across the arrow is checked: the "
           "Yang--Mills exponent $1/(2b_0g^2)$ equals $r_g\\beta$, with $r_g=12\\pi^2/(11N^2)$, while the "
           "$O(N)$ exponent is $r_s\\beta$ with $r_s=2\\pi/(N-2)$.")
    return d, _figure("fig:transmutation", d.tex(), cap, 0.98)

# ---------------------------------------------------------------------------------------
# the document
# ---------------------------------------------------------------------------------------
_MATH_TOKEN = re.compile(r"(\b\w+(?:_\{[^}]*\}|\^\{[^}]*\})|\b[A-Z]_[A-Z]+\b)")


def esc_title(text: str) -> str:
    """Escape a plain-text equation title for LaTeX, setting sub- and superscripted symbols (f_{p,q}, T^{3/2},
    A_XY) in math mode instead of printing the braces."""
    out = []
    for i, piece in enumerate(_MATH_TOKEN.split(text)):
        if i % 2 == 0:
            out.append(esc(piece))
        elif re.fullmatch(r"[A-Z]_[A-Z]+", piece):
            out.append(f"${piece[0]}_{{\\mathrm{{{piece[2:]}}}}}$")
        else:
            out.append(f"${piece}$")
    return "".join(out)


def selected_section(eqs_by_key) -> tuple[str, set[str]]:
    parts, keys = [], set()
    for title, entry, ekeys, *note in SELECTED:
        parts.append(f"\\subsection*{{{title}}}")
        if note:
            parts.append(note[0] + "\n")
        for k in ekeys:
            eq = eqs_by_key[k]
            cites = kb.citations_for(eq)
            keys.update(cites)
            parts.append(f"\\noindent {esc_title(eq.title)}{cite(cites)}."
                         f"\n\\begin{{equation}}\\label{{ceq:{eq.key}}}\n  {eq.latex(multiline=True)}\n\\end{{equation}}")
            if eq.see:
                parts.append(f"\\noindent\\textit{{See}} {render_see(eq.see)}.")
    return "\n".join(parts) + "\n", keys


def show(eq) -> tuple[str, list[str]]:
    """An inline display of a registered equation: its title and citations, the equation, and how many checks it passed."""
    cites = kb.citations_for(eq)
    tex = (f"\\noindent {esc_title(eq.title)}{cite(cites)}.\n"
           f"\\begin{{equation}}\\label{{qeq:{eq.key}}}\n  {eq.latex(multiline=True)}\n\\end{{equation}}\n"
           f"\\noindent\\textit{{{len(eq.verified)} checks passed (\\texttt{{models/{esc(eq.model)}.py}}).}}\n")
    return tex, cites


def numbers(eqs) -> dict:
    """The counts the abstract and introduction quote, computed at build time."""
    import pkgutil
    import models
    lib = kb.library()
    qg = [e for e in eqs if e.model in QG_MODELS]
    return {
        "nEq": len(lib.equations), "nQty": len(lib.quantities), "nLinks": len(lib.links),
        "nModels": len(list(pkgutil.iter_modules(models.__path__))),
        "nRegEq": len(eqs), "nChecks": len(eqlib.checks()),
        "nQgEq": len(qg), "nQgChecks": sum(len(e.verified) for e in qg),
        "nEinsteinChecks": sum(len(e.verified) for e in eqs if e.model == "qg_einstein4"),
    }


def numbers_tex(nums: dict) -> str:
    return ("% Generated by sympy_models/qgravity.py. Do not edit by hand.\n"
            + "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in nums.items()))


def figref(label: str, pic: bool) -> str:
    return f"Figure~\\ref{{{label}}}" if pic else "the figure (omitted: matplotlib is missing)"


def introduction(nums: dict, figure_checks: int, pic: bool) -> str:
    return (
        "\\section{Introduction}\\label{sec:intro}\n"
        "Loop quantum gravity describes quantum geometry on the Hilbert space of spin networks: graphs whose links carry spins of "
        "$SU(2)$, and whose area and volume operators have discrete spectra~\\cite{rovelli-smolin-1995,ashtekar-lewandowski-1997}. "
        "The Hamiltonian formulation of $SU(2)$ lattice gauge theory~\\cite{kogut-susskind-1975} lives on the same space on a fixed graph, "
        "and the square of the electric field on a link is the Casimir $j(j+1)$ of its spin, the number whose square root the area "
        "operator takes. This repository has already solved two small cases of that lattice theory exactly: one plaquette, and two "
        "plaquettes sharing a link, whose gauge-invariant states are spin networks on the theta graph, with seven links and two "
        "trivalent vertices. This paper reads the second system, the \\emph{ladder}, as the kinematics of loop quantum gravity, writes "
        "the dictionary down, and checks it. Four questions organise it.\n"
        "\\begin{itemize}\\itemsep2pt\n"
        "\\item \\textbf{What does the geometry of the ladder's ground state look like?} The surface that separates the two plaquettes "
        "cuts only the shared link, so its mean-square area is the ground-state value of $c_2(J)$. It is $2/(3g^8)$ at strong coupling "
        "and $3/(\\sqrt{10}\\,g^2)$ at weak coupling, both exact, and both checked against numerical ground states "
        "(Section~\\ref{sec:ladder}).\n"
        "\\item \\textbf{What does counting punctures fix?} For a horizon pierced by $N$ links of spin $\\tfrac12$, the entropy "
        "$S=A/4\\ell_P^2$ is matched by the count of invariant states if the Immirzi parameter is $\\ln2/(\\pi\\sqrt3)$. We derive this from "
        "the library's own entropy and Planck length, and keep it as a property of that count, not of the area spectrum "
        "(Section~\\ref{sec:entropy}).\n"
        "\\item \\textbf{What is the smallest vertex with a volume of its own?} Four spin-$\\tfrac12$ legs have a two-dimensional space of "
        "intertwiners, and the triple-product operator on it has spectrum $\\pm\\sqrt3/4$: one volume, two orientations "
        "(Section~\\ref{sec:intertwiners}).\n"
        "\\item \\textbf{Which gravity entry of the release speaks this paper's language?} Entry 348, on Einstein four-manifolds, is a statement "
        "about the curvature operator of $\\mathrm{Spin}(4)=SU(2)_+\\times SU(2)_-$, whose Weyl blocks are spin-2 objects of the two $SU(2)$s. "
        "We compute those blocks for $S^4$, $\\mathbb{CP}^2$ and $S^2\\times S^2$ from their metrics, check the topology and the numbers the "
        "papers use, and state what corresponds to a spin network and what does not (Section~\\ref{sec:einstein}).\n"
        "\\end{itemize}\n"
        "As in the companion paper, every term of every displayed row is the source equation solved for the row's quantity, substituted "
        "back and required to give an identity, and every number is checked: \\nQgChecks{} checks over the "
        "quantum gravity model, \\nChecks{} over all registered models, and " + str(figure_checks) + " more over the figures pass in the "
        "build, and a failure stops it. The checks test that stated formulas agree with one another and with independent computation "
        "(Weyl's integration formula against numerical Haar quadrature, exact perturbation theory against numerical ground states, explicit "
        "matrices against exact eigenvalues); they do not re-prove any theorem or establish any physics.\n\n"
        "\\paragraph{What this is not.} It is not a model of quantum gravity. The graph is fixed, the Hamiltonian is that of a gauge "
        "theory and not the Hamiltonian constraint of gravity, nothing is invariant under diffeomorphisms, there is no volume operator "
        "beyond the triple product on four-valent intertwiner spaces of Section~\\ref{sec:intertwiners}, and there is no semiclassical or continuum limit. What the "
        "ladder supplies is a Hilbert space of spin networks and an operator that measures the flux through a link, and what is "
        "checked is what that gives.\n")


def conclusion(pic: bool) -> str:
    return (
        "\\section{What is checked, and what is not}\\label{sec:conclusion}\n"
        "The dictionary holds in every case tested. The link Casimir is the squared area in units of $8\\pi\\gamma_I\\ell_P^2$, and the "
        "ladder's electric energy is a sum of seven such squares. The area spectrum is discrete, with a smallest value "
        "$4\\sqrt3\\pi\\gamma_I\\ell_P^2$ and steps that fall to $\\tfrac12$. The invariants of a vertex are counted by Weyl's integration "
        "formula, agree with the triangle rule at trivalent and four-valent vertices, and are Catalan numbers for spin-$\\tfrac12$ legs. "
        "In the ground state of the ladder the surface that separates the plaquettes has a mean-square area that falls as $g^{-8}$ at "
        "strong coupling and grows as $g^{-2}$ at weak coupling, with exact coefficients. And counting spin-$\\tfrac12$ punctures against "
        "the Bekenstein--Hawking entropy gives $\\gamma_I=\\ln2/(\\pi\\sqrt3)$, which is a statement about that count.\n\n"
        "Six limits apply throughout. (i)~The ladder's Hamiltonian is a gauge-theory Hamiltonian, so its ground state is the ground "
        "state of that theory and not a state of quantum gravity; the geometry computed in it is the geometry that Hamiltonian "
        "induces on a fixed graph. (ii)~The value $0.1274$ of the Immirzi parameter comes from counting spin-$\\tfrac12$ punctures only; "
        "counting all spins gives a different number~\\cite{domagala-lewandowski-2004}, which the fork quotes and does not rederive. "
        "(iii)~Of the volume operator only the triple product on four-valent intertwiner spaces is computed, not its prefactor and "
        "not a volume operator on general vertices; the one place a four-valent lattice vertex appears (the L-shaped cluster of three plaquettes) "
        "is computed with the cluster's Hamiltonian, but only for the triple product, at couplings $g^2\\ge1$. (iv)~The strip is checked for "
        "$N=1,2,3$ plaquettes; the strip formula for general $N$ is a prediction, tested at three values. (v)~Nothing here concerns the dynamics of gravity, "
        "diffeomorphism invariance, spin foams or a limit in which the lattice spacing disappears. (vi)~No matter field is coupled to the spin networks; "
        "the group $SU(2)$ is the structure group of the lattice gauge theory (and a factor of $Spin(4)$ in the part on Einstein four-manifolds), not the weak "
        "isospin group of the electroweak theory; the coupling $g$ is a parameter of the Hamiltonian and not a time, and the problem of time of quantum gravity "
        "is not addressed.\n\n"
        "The natural next steps, in increasing distance from what is here: a volume prefactor for the triple product of Section~\\ref{sec:intertwiners}, and the weak-coupling region of the L-shaped cluster; whether $\\langle c_2(J)\\rangle$ at an "
        "interior link of a longer strip approaches a plateau; the plaquette operators $W_a$ and $W_b$, which are traces of holonomies around small loops and so "
        "of the kind a Thiemann-regularised Hamiltonian constraint is built from, composed with the volume, as the first step toward "
        "any dynamics of the geometry; and the entropy count over all spins, to see how far $\\gamma_I$ moves. Each is a function to add to "
        "\\texttt{models/qg\\_spin\\_network.py} and, when it has an exact statement, a row to Figure~\\ref{fig:qgeometry}.\n\n"
        "\\paragraph{Related material.} The checks of this paper are exact symbolic and seeded numerical computations, not machine-checked proofs. "
        "For a Lean~4 formalization of a quantum field theory (the free bosonic field in four Euclidean dimensions, against the Glimm--Jaffe "
        "axioms) see~\\cite{dhmn-2026}. Two works of the author's collaboration are the nearest neighbours in the repository, and neither is used "
        "in a computation here. The Spectre and Chair44 paper~\\cite{hartshorn-romanov} computes the internal space and the boundary growth of "
        "two hierarchical tilings, and its last section tests, with negative or generic outcomes, whether a substitution hierarchy can serve "
        "as an arrow or a clock; it makes no claim about gravity. Romanov's paper on three-dimensional lattice Yang--Mills theory in a forest "
        "gauge~\\cite{zenodo-23783316} fixes the links of a forest to the identity, as the spanning-tree gauge of the ladder does, and proves for "
        "$SU(2)$ at strong coupling and fixed lattice spacing a unique state and a transfer-matrix gap $m\\ge-\\log(7\\beta_W/4)$ for $\\beta_W<4/7$, "
        "and a bound $9.599991$ on the Hessian of the Wilson action in that gauge; by its own account some of its inputs are unrefereed. "
        "Its Euclidean couplings and transfer-matrix gaps are not those of the Hamiltonian $g^2$ and gaps used here, and no number of this paper "
        "is compared with it; the four-dimensional paper of Section~\\ref{sec:forest} extends it, and the checks of that section are those of the model "
        "\\texttt{forest\\_gauge\\_4d\\_theorems}.\n")


def pictures(n_checks: int) -> tuple[str, str]:
    """(figure of the quantum geometry of a link and its count, figure of the cut in the ladder)."""
    geo = (
        "\\begin{figure}[!htbp]\\centering\n  \\includegraphics[width=\\linewidth]{build/fig_qg_geometry.pdf}\n"
        "  \\caption{The area spectrum and the count of punctures. (a)~The single-link area eigenvalues $\\sqrt{j(j+1)}$, in units of "
        "$8\\pi\\gamma_I\\ell_P^2$, against $j$, with the line $j+\\tfrac12$ they approach from below, and (right axis) the step between "
        "neighbouring levels, which is above $\\tfrac12$ and tends to it. (b)~The entropy per puncture $\\tfrac1N\\ln C_{N/2}$ of $N$ "
        "spin-$\\tfrac12$ punctures, from the exact Catalan numbers, with the asymptotic formula of Section~\\ref{sec:entropy}; it tends to "
        "$\\ln2$, which fixes $\\gamma_I=\\ln2/(\\pi\\sqrt3)$ for this count.}\\label{fig:qg-geometry}\n"
        "\\end{figure}\n")
    lad = (
        "\\begin{figure}[!htbp]\\centering\n  \\includegraphics[width=\\linewidth]{build/fig_qg_ladder.pdf}\n"
        "  \\caption{The surface that separates the two plaquettes of the ladder. (a)~The mean-square area $\\langle c_2(J)\\rangle_0$ "
        "in the ground state, against $g^2$ (both axes logarithmic): the numerical ground states, the exact strong-coupling law "
        "$2/(3g^8)$, and the weak-coupling law $3/(\\sqrt{10}\\,g^2)$ from the virial theorem. (b)~$g^2\\langle c_2(J)\\rangle_0$ for the "
        "shared link and $g^2\\langle c_2(j_a)\\rangle_0$ for a link of one plaquette, against $g^2$: both approach their exact "
        "limits (arrows, at $g^2=0$) linearly. The ladder's Hamiltonian is a gauge-theory Hamiltonian, not a gravitational one; "
        "the cut is a surface of its graph, not a horizon.}\\label{fig:qg-ladder}\n"
        "\\end{figure}\n")
    return geo, lad


def write(out_dir: Path, eqs) -> list[Path]:
    """Build the quantum gravity body and its bibliography; raises if any check fails."""
    import kbdoc
    by_key = {e.key: e for e in eqs}
    body, bib, numfile = (out_dir / f"qgravity_{x}.tex" for x in ("body", "bibliography", "numbers"))
    missing = [k for k in REQUIRED if k not in by_key]
    if missing:                                    # a partial run of export.py: do not leave a stale paper behind
        for p in (body, bib, numfile):
            p.unlink(missing_ok=True)
        print(f"qgravity: skipped (needs the models qg_spin_network and su2_lattice_ladder; missing {', '.join(missing)})")
        return []
    phys, ftex = physics_figure()
    qgf, qtex = quantum_geometry_figure()
    tra, ttex = transmutation_figure()
    ein, einx = einstein_figure()
    pen, penx = penrose_figure()
    try:
        import qg_pictures
        pic_checks = qg_pictures.make(out_dir)
    except ImportError as exc:                     # matplotlib missing: the paper builds without Figures 4 and 5
        print(f"qgravity: pictures skipped ({exc})")
        pic_checks = []
    pic = bool(pic_checks)
    fig_checks = sum(len(d.checks) for d in (phys, qgf, tra, ein, pen)) + len(pic_checks)
    nums = numbers(eqs)
    geo_fig, lad_fig = pictures(len(pic_checks)) if pic else ("", "")

    area, c1 = show(by_key["qg_area_spectrum"])
    dic, c2_ = show(by_key["qg_casimir_dictionary"])
    inter, c3 = show(by_key["qg_intertwiners"])
    tet, c4 = show(by_key["qg_quantum_tetrahedron"])
    vq, _ = show(by_key["qg_volume_triple_product"])
    vs, _ = show(by_key["qg_volume_spectra"])
    vsc, _ = show(by_key["qg_volume_semiclassical"])
    vl, _ = show(by_key["qg_volume_l_shape"])
    lb, _ = show(by_key["su2_lshape_basis"])
    lv, _ = show(by_key["su2_lshape_volume"])
    lg, _ = show(by_key["su2_lshape_gap"])
    sb, _ = show(by_key["su2_strip_basis"])
    sw, _ = show(by_key["su2_strip_weak"])
    ss, _ = show(by_key["su2_strip_strong"])
    sg, _ = show(by_key["su2_strip_gap"])
    sn, _ = show(by_key["su2_stripn_basis"])
    basis, c5 = show(by_key["su2_ladder_basis"])
    geom, c6 = show(by_key["qg_ladder_geometry"])
    ent, c7 = show(by_key["qg_entropy_immirzi"])
    so4, c8 = show(by_key["qg_so4_casimirs"])
    blk, c9 = show(by_key["qg_einstein_blocks"])
    top, c10 = show(by_key["qg_einstein_topology"])
    pe, _ = show(by_key["qg_penrose_equality"])
    pq, _ = show(by_key["qg_penrose_quantum"])
    f4l, _ = show(by_key["fg4t_linefield"])
    f4p, _ = show(by_key["fg4t_partners"])
    f4t, _ = show(by_key["fg4t_table1"])
    f4r, _ = show(by_key["fg4t_theorems"])
    kh, _ = show(by_key["qg_kerr_horizons"])
    kp_, _ = show(by_key["qg_kerr_penrose"])
    kq, _ = show(by_key["qg_kerr_quantum"])
    mc, _ = show(by_key["qg_maps_counts"])
    me, _ = show(by_key["qg_maps_exponents"])
    mk, _ = show(by_key["qg_maps_su2k"])

    lines = ["% Generated by sympy_models/qgravity.py. Do not edit by hand.\n"]
    lines.append(introduction(nums, fig_checks, pic))

    lines.append("\\section{The library, conjoined}\\label{sec:physics}\n"
                 "Every statement about quantum gravity ends in the same few constants: $G$, $\\hbar$, $c$, $k_B$, and the Planck length built "
                 "from the first three. Equations that state the same quantity can be set equal once each is solved for it. "
                 f"Figure~\\ref{{fig:physics}} does this for {len(PHYSICS_ROWS)} quantities at once, among them the gravitational constant, which "
                 "appears in Newton's law, the Schwarzschild radius, the Bekenstein--Hawking entropy and the Planck length. "
                 "The equations are " + _names_line(phys.names) + ".\n")
    lines.append(ftex)

    lines.append("\\section{Quantum geometry in the library}\\label{sec:geometry}\n"
                 "The knowledge graph gains five equations for the geometry of a spin network, and four quantities to read them with (spin, "
                 "Casimir eigenvalue, Immirzi parameter, puncture count): the Casimir, the area spectrum, the area of $N$ spin-$\\tfrac12$ "
                 "punctures, the horizon area, and the leading count of puncture microstates. Figure~\\ref{fig:qgeometry} conjoins them with "
                 "the library's entropy and Planck length. The equations are " + _names_line(qgf.names) + ".\n")
    lines.append(qtex)
    lines.append("The area of a surface is a sum over the links that pierce it. A single link of spin $j$ adds an amount that is a function "
                 "of the Casimir alone.\n")
    lines.append(area)
    lines.append("Two properties of the spectrum matter below. It has a smallest nonzero value, and it is not equally spaced: the step from "
                 "one allowed spin to the next is $0.548$ in units of $8\\pi\\gamma_I\\ell_P^2$ from $j=\\tfrac12$ to $1$, $0.522$ from $1$ to "
                 "$\\tfrac32$, and falls to $\\tfrac12$ (" + figref("fig:qg-geometry", pic) + "a). The Immirzi parameter $\\gamma_I$ is a "
                 "free number in the spectrum.\n")
    lines.append("The same Casimir is the electric energy of a link in the lattice gauge theory, which is what lets the ladder serve as a "
                 "model. The next display states that exactly, link by link.\n")
    lines.append(dic)
    lines.append("This is a relation between two operators on one Hilbert space. It does not make the ladder's dynamics gravitational: the "
                 "Hamiltonian constraint of loop quantum gravity is a different operator, and nothing in the ladder's Hamiltonian is built from "
                 "areas in that way.\n")

    lines.append("\\section{Intertwiners and the quantum tetrahedron}\\label{sec:intertwiners}\n"
                 "A vertex of a spin network carries an intertwiner: an invariant vector in the tensor product of the spins that meet there. "
                 "The dimension of that space is the Haar integral of the product of characters, which Weyl's integration formula makes an "
                 "exact integer. At a trivalent vertex it is $1$ or $0$ according to the triangle rule, which is the rule that selects the "
                 "ladder's states and is why they are labelled by link spins alone. At a four-valent vertex it counts the intermediate spins, "
                 "and from there on a vertex has quantum numbers of its own.\n")
    lines.append(inter)
    lines.append("For four spin-$\\tfrac12$ legs the space is two-dimensional, and on it the triple product of three of the spins is a Hermitian "
                 "operator that commutes with the total spin. It is the squared volume of the quantum tetrahedron~\\cite{barbieri-1998} up to a "
                 "constant, and its spectrum is computed here from explicit $16\\times16$ matrices.\n")
    lines.append(tet)
    lines.append("The prefactor that turns $Q$ into a volume is not computed: the statement made is that the smallest vertex with quantum numbers "
                 "of its own has a single nonzero volume squared, $\\sqrt3/4$ in the units of $Q$, with two orientations. The same operator "
                 "is built for any four spins. In the coupled basis it is $i(x_K-x_{K'})$ times an orthogonal matrix of $6j$ symbols (the spectrum "
                 "does not depend on a phase that depends on one intermediate spin alone), and it agrees with the explicit-matrix construction "
                 "to $4\\times10^{-15}$ for eight sets of spins.\n")
    lines.append(vq)
    lines.append(vs)
    lines.append("The spectra are symmetric under $Q\\to-Q$, and vanish when a leg has spin $0$, when the vertex is trivalent. For four spin-$1$ legs "
                 "there is a state of zero volume, as the dimension, three, is odd. The largest eigenvalue, against the triple product of three "
                 "vectors of length $\\sqrt{j(j+1)}$ of the regular tetrahedron, rises from $0.794$ at $j=\\tfrac32$ to $0.921$ at $j=8$ and stays below one: "
                 "a trend, not a proof of a semiclassical limit.\n")
    lines.append(vsc)
    lines.append("\\paragraph{Where the four-valent vertex is.} After gauge fixing, every vertex of the strip of the next section is trivalent, so the strip "
                 "carries no volume. A four-valent lattice vertex first appears for three plaquettes in an L shape ($a$ at the origin, $b$ to its right, $c$ above), "
                 "at the corner that all three share: its four links are the one shared by $a$ and $b$, the top of $b$, the one shared by $a$ and $c$, and the "
                 "right side of $c$. At strong coupling the ground state is the electric vacuum dressed by strings of spin $\\tfrac12$ around plaquettes. "
                 "Only a string around $b$ together with one around $c$ puts spin $\\tfrac12$ on all four legs, with the two loops each "
                 "contracting a pair of legs, and on that vector $Q^2=\\tfrac3{16}$. Its amplitude is the second-order value $4/(9g^8)$, so that "
                 "$\\langle Q^2\\rangle_0=1/(27g^{16})$ at leading order.\n")
    lines.append(vl)
    lines.append("This estimate was first conditional: the amplitude $4/(9g^8)$ is the one two disjoint strings have, and it was not derived for the L shape. "
                 "The Hamiltonian of the L-shaped cluster is now built, and the estimate is tested against it. The basis labels the spins $j_a,j_b,j_c$, the shared "
                 "links $J_{ab},J_{ac}$, and an intermediate spin $K$ at the four-valent vertex, in the channel $(J_{ab}J_{ac})K$, $(j_bj_c)K$. The plaquette operators "
                 "are built vertex by vertex from explicit Clebsch--Gordan coefficients, with the recoupling at the four-valent vertex included, and the signs "
                 "are not guessed: the relation $\\chi_{1/2}(U_p)\\Psi_s=\\sum_{s'}\\langle s'|W_p|s\\rangle\\Psi_{s'}$ is checked pointwise at random group elements "
                 "for $45$ pairs of a state and a plaquette, with the wave functions evaluated by direct tensor contraction.\n")
    lines.append(lb)
    lines.append("The three plaquette operators are real, symmetric and commute in the interior of the truncation (and fail to at its boundary, as they should); "
                 "with a plaquette trivial they reduce in absolute value to the ladder's squared-$6j$ elements; the kinetic matrix of the harmonic limit "
                 "has eigenvalues $4,\\,4\\pm\\sqrt2$, those of the strip. The ground-state energy is $6/g^2-2/g^6+O(g^{-10})$ at strong coupling, and the amplitude of "
                 "the state with strings around $b$ and $c$ is $4/(9g^8)$ (ratio to the numerics tending to $1$ at $g^2=8,12,16$).\n")
    lines.append(lv)
    lines.append("The estimate is confirmed: the ratio $27\\,g^{16}\\langle Q^2\\rangle_0$ is $0.9808$, $0.9939$, $0.99961$ and $0.999976$ at $g^2=3,4,8,16$. "
                 "At intermediate coupling $\\langle Q^2\\rangle_0$ is $1.3\\times10^{-4}$, $1.09\\times10^{-3}$ and $1.34\\times10^{-2}$ at $g^2=2,1.5,1$, growing "
                 "as the coupling falls. The number at $g^2=0.8$ was dropped because it is converged only to about $8\\%$ in the spin cutoff, and those above to about $2\\%$. "
                 "The constant that turns $Q$ into a volume is still not fixed, so this is a statement about the triple product, not about a physical volume.\n")
    lines.append(lg)
    lines.append("The lowest excitation is $2.4398$ at $g^2=1.5$ and $2.0570$ at $g^2=1$, against the harmonic limit $2.2741$; at $g^2=1.5$ it equals the strip's to $10^{-3}$, "
                 "so the corner vertex does not show in the gap there. The weak-coupling limit itself is not reached by these numbers (the gap is below the limit at $g^2=1$), "
                 "and no larger cluster or continuum limit was tried.\n")

    lines.append("\\section{The ladder as a spin network}\\label{sec:ladder}\n"
                 "The two plaquettes sharing a link have seven links. After a spanning-tree gauge fixing, the gauge-invariant states are the "
                 "spin networks of the theta graph: the three outer links of plaquette $a$ carry one spin $j_a$, those of $b$ carry $j_b$, "
                 "and the shared link carries $J$, with the triangle rule at both vertices. The Hamiltonian is the Kogut--Susskind Hamiltonian "
                 "in lattice units, with the plaquette operators given by squared $6j$ symbols.\n")
    lines.append(basis)
    lines.append("A surface that separates the two plaquettes cuts one link, the shared one, so its area is the area spectrum at spin $J$ and "
                 "its mean square in the ground state is $\\langle c_2(J)\\rangle_0$ in units of $(8\\pi\\gamma_I\\ell_P^2)^2$. At strong coupling "
                 "the ground state is the electric vacuum $J=0$, dressed to first order by the two single-plaquette strings; at weak coupling "
                 "it spreads over large spins and the ladder is two coupled oscillators, whose virial theorem fixes the coefficient.\n")
    lines.append(geom)
    lines.append(lad_fig)
    lines.append("The two laws meet the numerical curve at both ends (" + figref("fig:qg-ladder", pic) + "a), and the weak-coupling "
                 "approach is linear in $g^2$ (" + figref("fig:qg-ladder", pic) + "b). The shared link fluctuates more than a "
                 "plaquette's own links by the factor $1.746$, because it carries the difference of the two loop fields. Only the leading "
                 "terms are derived; the first correction in $g^2$ is not.\n")

    lines.append("\\section{A third plaquette: the strip}\\label{sec:strip}\n"
                 "Three plaquettes in a row have ten links. The corner vertices again force equal spins on the outer links of the end plaquettes, "
                 "so the gauge-invariant states are spin networks on a graph of four trivalent vertices and six edges: $j_a$ and $j_c$ for the ends, "
                 "$t$ and $b$ for the top and bottom of the middle plaquette, and $J_1$, $J_2$ for the two shared links. The electric energy is again "
                 "a sum of squared areas, and the plaquette operators shift every link of the boundary of a plaquette by $\\pm\\tfrac12$ with a product "
                 "of one $6j$ symbol per vertex on the loop: two for $a$ and $c$, four for $b$. Unlike the ladder's, these $6j$ symbols are not squared, "
                 "so the signs matter.\n")
    lines.append(sb)
    lines.append("The signs in the display were found by fitting, not derived from a convention: they were fixed against an independent construction "
                 "from Wigner $3j$ symbols that shares no code with the closed forms, and agree with it on all $5047$ transitions checked "
                 "(states up to spin $\\tfrac32$ and a sample up to spin $2$) to $4\\times10^{-16}$. The three plaquette operators are symmetric "
                 "and commute in the interior of the truncated space ($6\\times10^{-16}$), while the commutators reach $0.94$ at the truncation "
                 "boundary, as expected, and with the third plaquette trivial the operators are the ladder's, exactly, on $948$ transitions.\n")
    lines.append("This is the test that the strip formula of the ladder's harmonic argument lacked. The three loop variables are three-dimensional "
                 "oscillators with the tridiagonal kinetic matrix of the display; gauge invariance leaves total angular momentum zero, so the lowest "
                 "excitation is two quanta in the soft mode.\n")
    lines.append(sw)
    lines.append("Against this, the lowest even level of the numerical spectrum, extrapolated to $g\\to0$ from $g^2=0.2$--$0.4$ with a cubic, is "
                 "$2.2749$ ($2.2741$ predicted), and the second is $2.5520$ ($2.5513$ predicted), with the reflection parities of the oscillator modes. "
                 "The next two even levels mix at finite coupling, and only their sum could be compared (agreement to $1.3\\times10^{-3}$). "
                 "At strong coupling the first excited level is threefold degenerate, three strings of spin $\\tfrac12$ at $3g^2/2$, and exact "
                 "second-order perturbation theory splits it.\n")
    lines.append(ss)
    lines.append("The numerical levels at $g^2=6$ and $8$ match these to $5\\times10^{-4}$ with the right parities. Between the limits the three "
                 "gaps are ordered.\n")
    lines.append(sg)
    lines.append("The ordering holds on the grid $g^2=0.25$--$8$; the ratio to the one-plaquette gap is $0.835$ at $g^2=0.25$, against the limit "
                 "$0.8040$, and $0.99995$ at $g^2=8$. The three volumes say nothing about a larger or infinite system, and, like the ladder's, "
                 "the gap here is that of a gauge-theory Hamiltonian on a small lattice.\n")
    lines.append("\\paragraph{Four and five plaquettes.} The same graph, a ladder of squares, gives $3N-3$ labels for $N$ plaquettes, and the display "
                 "gives the Hamiltonian, with a phase $\\sigma_p$ that extends the three-plaquette rule by a factor $(-1)^{(k_t-k_b)/2}$ for each "
                 "interior neighbour. That extension is a conjecture, and it is tested rather than derived.\n")
    lines.append(sn)
    lines.append("The test is an independent construction with no $6j$ symbol and no phase: wavefunctions of the $N$ loop holonomies, with the "
                 "link electric fields from the Peter--Weyl decomposition and gauge invariance enforced by a penalty. It reproduces the ladder "
                 "exactly at $N=2$ and the strip module at $N=3$. At $g^2=2$ the spin-network builder and the independent construction agree on "
                 "the lowest four levels to $2\\times10^{-6}$ at $N=4$ (ground state $3.67192$) and $10^{-3}$ at $N=5$ (ground state $4.5899$, "
                 "limited by the cut), and the plaquette operators commute away from the cut. The weak-coupling gap formula "
                 "$2\\sqrt{2-\\cos(\\pi/(N+1))}$ ($2.1827$ at $N=4$, $2.1297$ at $N=5$) is \\emph{not} tested: at $g^2=0.8$--$0.4$ the $N=4$ "
                 "gap is still falling with the energy cut ($2.84\\to2.49$ at $g^2=0.4$ between cuts $24$ and $32$), and $N=5$ was too slow "
                 "to converge, so no extrapolation is claimed.\n")

    lines.append("\\section{Counting punctures and the Immirzi parameter}\\label{sec:entropy}\n"
                 "If a horizon of area $A$ is pierced by $N$ links of spin $\\tfrac12$, its area is $N$ times the smallest value, and the "
                 "number of invariant states is a Catalan number whose logarithm is $N\\ln2$ at leading order. Setting "
                 "$S=A/4\\ell_P^2$ against that count fixes $\\gamma_I$. The display uses the library's own entropy and Planck length, in "
                 "which $G$, $\\hbar$ and $c$ cancel~\\cite{rovelli-1996}.\n")
    lines.append(ent)
    lines.append(geo_fig)
    lines.append("The count is exact (Section~\\ref{sec:intertwiners}) and the shortfall of the entropy per puncture below $\\ln2$ is "
                 "$\\tfrac{3}{2N}\\ln\\tfrac N2+\\tfrac{\\ln\\pi}{2N}+\\tfrac{9}{4N^2}$ (" + figref("fig:qg-geometry", pic) + "b), still $0.007$ at "
                 "$N=1600$. The value $0.1274$ belongs to this count: counting all spins with the area fixed gives a different number "
                 "(about $0.2375$~\\cite{domagala-lewandowski-2004}, quoted, not rederived here). The area spectrum does not determine "
                 "$\\gamma_I$, and neither does this section derive the entropy of a black hole from dynamics: the horizon area is given.\n")

    lines.append("\\section{Einstein four-manifolds: the geometric \\texorpdfstring{$SU(2)$}{SU(2)}}\\label{sec:einstein}\n"
                 "Entry 348 of the release classifies the closed four-manifolds that are Einstein, $\\mathrm{Ric}=\\lambda g$, and have positive "
                 "sectional curvature: up to scale, the round $S^4$, Fubini--Study $\\mathbb{CP}^2$ and (for the real projective quotient) "
                 "$\\mathbb{RP}^4$~\\cite{oai-einstein4}. Two companion papers carry it to the boundary case, where $S^2\\times S^2$ enters, and to "
                 "an $L^2$ gap: a simply connected nonnegatively curved four-manifold whose trace-free Ricci tensor has small enough "
                 "$L^2$ energy is diffeomorphic to one of the three~\\cite{oai-einstein-gap,oai-zero-plane}. This is the one gravity entry of the "
                 "release whose language is the language of this paper. The frame group of an oriented four-manifold is $\\mathrm{Spin}(4)=SU(2)_+\\times "
                 "SU(2)_-$, the two-forms split as $\\Lambda^2=\\Lambda^+\\oplus\\Lambda^-$ into the adjoints of the two factors, and the "
                 "curvature operator for $\\mathrm{Ric}=3g$ is $\\mathrm{Id}+T$ on $\\Lambda^+$ and $\\mathrm{Id}+U$ on $\\Lambda^-$, with $T$ and $U$ the two "
                 "trace-free Weyl blocks~\\cite{oai-einstein4}. Each block is a symmetric trace-free endomorphism of a three-dimensional space "
                 "on which one $SU(2)$ acts: a spin-2 object. The link of a spin network carries a spin of one $SU(2)$; the self-dual "
                 "(Ashtekar) connection of a Euclidean four-manifold is the $SU(2)_+$ part of its Spin$(4)$ connection. What follows checks the "
                 "algebra of that statement and the numbers the paper and the library use about the three model spaces.\n")
    lines.append(so4)
    lines.append("The Casimirs are computed from explicit $4\\times4$ generators. On the vector they are $(\\tfrac34,\\tfrac34)$, on $\\Lambda^+$ they "
                 "are $(2,0)$, and on a Weyl block, the five-dimensional trace-free part of the symmetric square of $\\Lambda^+$, the Casimir is $6=j(j+1)$ "
                 "at $j=2$: the same function of the spin as in the area spectrum (Section~\\ref{sec:geometry}), evaluated at spin $2$. The last "
                 "line of the display is the algebraic form of a decomposable two-plane having equal self-dual and anti-self-dual parts, which "
                 "is also why the sectional curvature of a plane is a mean of one value from each block. In spin-foam models the corresponding "
                 "statement about a bivector (the simplicity constraint) is what ties the two spins to one spin through the Immirzi "
                 "parameter~\\cite{eprl-2008}; that is context here, and nothing in it is computed.\n")
    lines.append("The curvature blocks of the three spaces are computed from their metrics in a chart, by differentiating twice, and not from a formula for "
                 "the curvature of the model.\n")
    lines.append(blk)
    lines.append("The mixed block between $\\Lambda^+$ and $\\Lambda^-$ vanishes for the three Einstein spaces and not for the unbalanced product "
                 "$S^2(K{=}3)\\times S^2(K{=}1)$, whose Ricci tensor is $\\mathrm{diag}(3,3,1,1)$: that block is the trace-free Ricci tensor, and its "
                 "vanishing is the Einstein equation. The least sectional curvature is $1-(x+y)/2$, where $x$ and $y$ are the depths of the two blocks, "
                 "and the condition $x+y<2$ of the paper is strict sectional positivity: $S^4$ ($0$) and $\\mathbb{CP}^2$ ($1$) satisfy it, $S^2\\times S^2$ "
                 "($2$) has planes of zero curvature. On every nonzero block of the three spaces $|T|^2$ or $|U|^2$ is $6$, so $v^2=b^2=1$: "
                 "the gap $\\mathbb{E}\\,v^2\\ge1$ that the paper takes from~\\cite{gursky-lebrun-1999} holds with equality on each, and $S^2\\times S^2$ "
                 "has $v+b=2$, exactly on the boundary of the condition $v+b<2$. The topology comes out of the same curvature.\n")
    lines.append(top)
    lines.append(einx)
    lines.append("Figure~\\ref{fig:einstein} puts the two characteristic identities in the library's own form. They are the ones that tie the two second "
                 "moments to the volume in the paper's argument, and they are what the three rows of $s_S$ test: the Euler number and the signature, "
                 "read as equations for the same moment, agree on all three spaces. Two numbers of the papers appear as checks of the model: the "
                 "Hitchin--Thorpe quantity $2\\chi-3|\\tau|$ is $4$, $3$, $8$, and the condition $2\\chi-3|\\tau|\\le4$ under which a Gursky--Malchiodi-type "
                 "classification is quoted in the paper holds for the first two only. The energy $D(g)=\\int|\\mathrm{Ric}-\\tfrac R4g|^2\\,d\\mu$ of the "
                 "$L^2$ paper is invariant under rescaling in dimension four, is zero on the three Einstein spaces, and is $64\\pi^2/3$ on the "
                 "unbalanced product. The existence of the universal threshold $\\varepsilon_0$ below which the classification persists is the theorem "
                 "of~\\cite{oai-einstein-gap}; it has no value in the paper, and none is computed here.\n")
    lines.append("\\paragraph{What corresponds to what, and what does not.} (i)~The $SU(2)$ of a link is, in a Euclidean four-manifold, one of the two "
                 "factors of the frame group, and the electric field of a link is a section of the adjoint of that factor, which is $\\Lambda^+$ when "
                 "the factor is $SU(2)_+$. (ii)~The two-dimensional surface whose area the spin network measures has a bivector $\\Sigma$, "
                 "decomposable, with $|\\Sigma_+|=|\\Sigma_-|$; its area is the length of $\\Sigma$ and the area spectrum quantises the length of the "
                 "$\\Sigma_+$ part, the Casimir of Section~\\ref{sec:geometry} (the factor $1/\\sqrt2$ is a normalisation that is not tracked). "
                 "(iii)~The weak-coupling variables of the ladder are $su(2)$-valued, three components per plaquette; on the geometric side "
                 "the curvature $2$-form of $SU(2)_+$ evaluated on a plane is a vector of $\\Lambda^+\\cong\\mathbb R^3$. That the numbers of three "
                 "components agree is the whole of the correspondence at this point, and it is an interpretation. Not claimed: that the ladder, which is a "
                 "two-dimensional lattice problem, or any spin network, carries an Einstein metric; that $S^4$, $\\mathbb{CP}^2$ or $S^2\\times S^2$ is a "
                 "state of any quantum theory here; anything about the Lorentzian theory, in which the group of the Ashtekar--Barbero connection is $SU(2)$ "
                 "for a different reason; or any statement about the classification theorem itself, which is stated in Lean in the release "
                 "(\\texttt{lean/docs/348.md}) and not connected to the checks above. The other gravity entries of the release (the Liouville quantum "
                 "gravity limit of random planar maps, entry 211; the spacetime Penrose inequalities, 260; strong cosmic censorship near Kerr, 264) "
                 "have not been examined: of them the Penrose inequality is the closest to what is computed here, since it bounds a mass "
                 "from below by an area, and the area has a spectrum in Section~\\ref{sec:geometry}, but combining a classical inequality with a quantum "
                 "area eigenvalue is a heuristic that this paper does not make.\n")

    lines.append("\\section{The Penrose inequality and the smallest area}\\label{sec:penrose}\n"
                 "Entry 260 proves that the invariant ADM mass $m$ of asymptotically flat initial data satisfying the dominant energy condition is at least "
                 "$\\sqrt{A_{\\min}/16\\pi}$, where $A_{\\min}$ is the minimum area of a surface enclosing a weakly trapped boundary, in three spatial dimensions "
                 "and with $G=c=1$~\\cite{oai-penrose}. It states charged, four-dimensional and anti-de Sitter versions, and equality for Schwarzschild "
                 "data under horizon hypotheses. The inequality is the one place in the release where a mass is bounded by an area, and the area "
                 "has a spectrum in Section~\\ref{sec:geometry}. The next display restores the units, checks that the library's Schwarzschild radius and horizon "
                 "area saturate the bound, and checks the equality cases of the stated variants as algebraic identities. They are consistency checks "
                 "of formulas, not tests of the theorems.\n")
    lines.append(pe)
    lines.append(penx)
    lines.append("Figure~\\ref{fig:penrose} puts the saturated bound in the library's form. If the area in the bound is taken to be that of $N$ spin-$\\tfrac12$ "
                 "punctures with $\\gamma_I=\\ln2/(\\pi\\sqrt3)$, the bound is $M/M_P\\ge\\sqrt{N\\ln2/4\\pi}$, which is $\\sqrt{S/4\\pi k_B}$ in terms of the "
                 "entropy of Section~\\ref{sec:entropy}: for one puncture, $0.2349$ Planck masses.\n")
    lines.append(pq)
    lines.append("\\paragraph{Why this is a heuristic.} The theorem is about classical initial data, and the left-hand side of the bound is a number computed from "
                 "such data; an area eigenvalue belongs to an operator on a Hilbert space of spin networks. Nothing here derives a statement about a quantum state "
                 "from the theorem, or the converse, and the value of $\\gamma_I$ used belongs to one counting (Section~\\ref{sec:entropy}). The numbers also show where "
                 "the classical reading breaks. For one puncture the area radius is $0.47\\,\\ell_P$, below the Planck length, and the reduced Compton wavelength of "
                 "the saturating mass, $4.3\\,\\ell_P$, is larger than the horizon. The crude criterion that the Schwarzschild radius of the saturating mass exceeds "
                 "its Compton wavelength first holds at $N=10$ punctures. What survives is arithmetic: an exact closed form relating the smallest quantised area, "
                 "the entropy and a mass in Planck units, which a reader can set beside the theorem and decide what it means.\n")

    lines.append("\\section{Strong cosmic censorship near Kerr}\\label{sec:kerr}\n"
                 "Entry 264 proves a local form of strong cosmic censorship near every fixed rotating subextremal Kerr bridge: in a weighted smooth topology, a dense $G_\\delta$ "
                 "set of nearby smooth two-ended asymptotically flat vacuum data has a maximal globally hyperbolic development with no future extension, in the classes "
                 "$C^2\\subset C^1\\subset C^0\\cap W^{1,2}_{\\mathrm{loc}}$ of the three companion papers (the last class, a continuous metric with a locally square-integrable "
                 "weak connection, gives the strongest statement). No symmetry is imposed, and extensions need not satisfy the vacuum equations~\\cite{oai-kerr-w12,oai-kerr-c1,oai-kerr-c2}. "
                 "What the proofs use of Kerr itself is small and explicit: a reduction by scaling to $M=1$, an inner horizon $r_->0$, and the absolute inner surface gravity "
                 "$\\kappa$, which fixes the decay of the lapse as $e^{-\\kappa T}$ and the packet amplitude $A_T=e^{-\\kappa T/4}$, larger than the square root of the lapse. "
                 "The extension is then tested by holonomies of a family of small loops~\\cite{oai-kerr-w12}. The first display computes, in exact algebra, the Kerr quantities "
                 "these statements rest on; the fork read the abstract, introduction and the normalisation of the first paper and not the proofs.\n")
    lines.append(kh)
    lines.append(kp_)
    lines.append("The second display records one exact coincidence. The Penrose inequality of Section~\\ref{sec:penrose} is an equality for Schwarzschild and strict for every rotating Kerr, "
                 "with gap $M^2-A_+/16\\pi=Mr_-/2$: it vanishes exactly when the inner horizon does, which is the condition $r_->0$ that the near-Kerr papers need. The two "
                 "theorems are different statements about different objects, and the fork can say only that they share this hypothesis. The first display also shows why "
                 "the constants of the papers are tied to the fixed centre: $|\\kappa_-|\\sim M/a^2$ grows without bound as $a\\to0$ and both surface gravities vanish as $a\\to M$.\n")
    lines.append(kq)
    lines.append("\\paragraph{What is, and is not, a connection.} Two things join this entry to the paper. The first is arithmetic and heuristic: composed with the area "
                 "spectrum, the classical bound $A_+\\ge8\\pi|J|$ asks for at least five spin-$\\tfrac12$ punctures at the smallest spin $J=\\hbar/2$ and ten at $J=\\hbar$ "
                 "(the same $2\\pi/\\ln2$ as the Compton threshold of Section~\\ref{sec:penrose}), and the Kerr mass bound $M\\ge\\sqrt{|J|}\\,M_P$ then exceeds the Penrose bound for those "
                 "punctures. It inherits every caveat of Section~\\ref{sec:penrose}: a classical bound against an eigenvalue of an operator. It also stops at the outer horizon: the area "
                 "product $A_+A_-=(8\\pi J)^2$ would make $A_-$ a spin-network area too only if $NN'=4\\pi^2J^2/\\ln^22$ were an integer, which is not the case at $J=\\tfrac12$ "
                 "(it is $20.54$), and the inner horizon is not a surface of the states of Section~\\ref{sec:ladder} in any case. The second is structural. The proofs of entry 264 test extensions by "
                 "holonomies of small loops of the Levi-Civita connection, and a spin network is built from holonomies of a connection along the links of a graph. "
                 "The resemblance ends there: the connection of the proof is the Lorentzian one on the spacetime, with only $L^2$ coefficients in a chart; the connection of this paper "
                 "is an $SU(2)$ connection on a spatial slice. Nothing in the displays tests this analogy.\n")

    lines.append("\\section{Liouville quantum gravity and the Catalan counts}\\label{sec:maps}\n"
                 "Entry 211 proves that critical Fortuin--Kasteleyn planar maps converge, for $0<q\\le4$, to Liouville quantum gravity spheres decorated by conformal loop ensembles, "
                 "and for $q>4$ to the Brownian continuum random tree, which is the transition from surface to tree; for FK--Ising ($q=2$, $\\gamma=\\sqrt3$) and for "
                 "spanning-tree-decorated maps the graph metric and the stationary random walk converge as well, and the spectrum "
                 "converges in the Ising case~\\cite{oai-fk-sub,oai-fk-q4,oai-fk-metric,oai-fk-walk,oai-fk-spectral}. The parameters are linked by $q=2+2\\cos(\\pi\\gamma^2/2)$ and "
                 "$\\kappa'=16/\\gamma^2$~\\cite{oai-fk-sub}. This is the entry least close to loop quantum gravity: Liouville gravity is a continuum random surface in two dimensions whose area is a "
                 "multiplicative chaos measure, and the two theories share the initials and little else. The one place where the counting is the same is the one the paper already uses: the Catalan numbers "
                 "that give the dimension of the invariant space of $N$ spin-$\\tfrac12$ legs (Section~\\ref{sec:intertwiners}) and the entropy of Section~\\ref{sec:entropy} "
                 "also count the plane trees of the tree phase and, shuffled in pairs, the tree-decorated maps. The displays check what can be recomputed.\n")
    lines.append(mc)
    lines.append(me)
    lines.append(mk)
    lines.append("\\paragraph{Reading the three displays.} The brute-force counts reproduce Tutte's $2,9,54,378$ rooted maps and Mullin's $2,10,70,588$ tree-rooted maps, "
                 "and Mullin's number is a shuffle of two Catalan structures, one for the tree and one for its dual. The dictionary of the second display turns the entry's "
                 "parameters into central charges and exponents ($q=0,1,2,3,4$ give $c=-2,0,\\tfrac12,\\tfrac45,1$), and the rule that rooted counts grow as $\\mu^nn^{\\gamma_{\\rm str}-2}$ "
                 "reproduces the exponents $-3$ and $-\\tfrac52$ found by brute force. The third display is the connection with this paper. The invariants of $N$ spin-$\\tfrac12$ "
                 "legs are counted by closed walks on a Dynkin diagram $A_{k+1}$ at level $k$ of the quantum group $U_q(\\mathfrak{sl}_2)$, with loop weight $2\\cos(\\pi/(k+2))$; "
                 "setting that weight equal to $\\sqrt q$ gives $q=1,2,3$ at $k=1,2,4$ (percolation, Ising, three-state Potts), and the undeformed $SU(2)$ of this paper, "
                 "loop weight $2=\\dim V_{1/2}$, is the critical end $q=4$ of the range of the entry. The same Temperley--Lieb algebra acts in both, which is why the numbers agree; the fork constructs "
                 "no map from a spin network to a decorated planar map, and the growth rate $\\tfrac12\\ln q$ of the count is the growth rate of a count and not an Immirzi-type constant (the "
                 "$\\ln2$ of Section~\\ref{sec:entropy} is the $k=\\infty$ value, fixed there by a physical argument). A loop weight above $2$, the "
                 "$q>4$ tree phase of the entry, has no unitary $SU(2)_k$ counterpart. Quantum-group deformations of loop quantum gravity, in which $k$ is set by a cosmological constant, are "
                 "a known direction that this paper does not use and for which it cites nothing. The classical formulas used (Tutte, Mullin, the Kardar--Parisi--Zhang relation for $\\gamma_{\\rm str}$) "
                 "are quoted by name and checked by brute force or exact algebra where possible, not taken from a retrieved source.\n")

    lines.append("\\section{Gauge fixing and the curvature of the link: a four-dimensional lattice result}\\label{sec:forest}\n"
                 "A paper dated 30 October 2026 by Bhattacharjee, Bhattacharya, Lockwood and Romanov~\\cite{bbll-forest-4d} studies the Euclidean Wilson action "
                 "on $\\mathbb Z^4$ at fixed lattice spacing. Fixing the variables on a \\emph{line field} (a set $F$ of straight lines, one through every site, "
                 "chosen so that every plaquette contains exactly one link of $F$) to the identity, it proves a unique infinite-volume state and an explicit "
                 "transfer-matrix gap: $m\\ge-\\log(3\\beta_W)$ for $SU(2)$ and $\\beta_W<\\tfrac13$, plaquette blocks that reach $0.345$ for $SU(2)$ and $0.76$ for $SU(3)$, "
                 "and for every compact simple group a range $\\beta_W<2N\\kappa_G/29$ with a decay rate $m_G(\\beta)$ given by the least Ricci curvature $\\kappa_G$. "
                 "It reduces the continuum gap to one cube condition and says that this condition, the slice ratio and the tree bounds are open. The fork read the first "
                 "six pages and the appendix with the line field; it did not read the proofs, and the paper's own account is that these are strong-coupling results "
                 "at a fixed spacing. What the fork could do is rebuild the finite statements: Table~III of the paper (the line field), the partner counts, "
                 "Table~I from root-system data, and the numbers of the theorems. The next displays are those checks. They test the paper's arithmetic and "
                 "combinatorics, not its theorems.\n")
    lines.append(f4l)
    lines.append(f4p)
    lines.append(f4t)
    lines.append(f4r)
    lines.append("Three observations connect this to the rest of the paper. None is a derivation.\n"
                 "\\paragraph{The link space is the Einstein three-sphere, and the smallest spin saturates Lichnerowicz.} The curvature that the "
                 "paper uses for $SU(2)$ is the Ricci curvature of the group manifold $S^3$ with the bi-invariant metric $-\\mathrm{Tr}(XY)$, which is $\\kappa=1$ "
                 "(Einstein with constant $2$ on the unit sphere, here radius $\\sqrt2$). The electric term of the Hamiltonian of Section~\\ref{sec:ladder} is "
                 "the Laplacian of that sphere, with eigenvalue $2j(j+1)$ in this metric. Lichnerowicz's bound $\\lambda_1\\ge\\tfrac{n}{n-1}\\kappa$ is "
                 "$\\tfrac32$ on $S^3$, and the spin-$\\tfrac12$ eigenvalue is $\\tfrac32$: equality, as Obata's theorem requires of the round sphere. "
                 "The check above finds the same comparison strict for every other group in Table~I, among them $SU(3)$ ($8/3$ against $12/7$). "
                 "The spin $\\tfrac12$ is also the one whose area is the smallest, $\\sqrt{j(j+1)}=\\sqrt3/2$ (Section~\\ref{sec:geometry}), "
                 "so the smallest area quantum of a link and the first electric excitation that the curvature argument of the paper controls are the "
                 "same eigenfunction of the same Laplacian. This is a fact about $S^3$, not an explanation of either result.\n"
                 "\\paragraph{Dimension four, twice.} A field with one fixed link at every site and one in every plaquette needs the number of plaquettes through "
                 "a link, $2(d-1)$, to equal the number of plaquettes at a site, $d(d-1)/2$, and the only solution is $d=4$ (checked above). "
                 "The number $d(d-1)/2$ is also the dimension of the space of two-forms, and the splitting $\\Lambda^2=\\Lambda^+\\oplus\\Lambda^-$ into "
                 "two copies of the adjoint of $SU(2)$, on which Section~\\ref{sec:einstein} rests, exists only for $d=4$ as well. The two statements are "
                 "different (one is a count of cells, the other a Hodge star), and no mechanism relates them; the coincidence is recorded because "
                 "the same equation $(d-1)(d-4)=0$ decides both.\n"
                 "\\paragraph{Gauge fixing and the spin-network basis.} The forest gauge solves the Gauss law by setting the holonomies on $F$ to the identity, "
                 "which leaves three free links at every site; the spin-network basis of this paper solves it by coupling spins at vertices to invariants, and the "
                 "ladder of Section~\\ref{sec:ladder} was fixed in a spanning-tree gauge to count its states. These are two ways of removing the same redundancy, "
                 "and nothing is compared beyond that.\n"
                 "\\paragraph{What does not follow.} The gap of Section~\\ref{sec:ladder} is a Hamiltonian gap of a two-plaquette system in units of $g^2$; the gaps "
                 "above are transfer-matrix gaps in lattice units at strong Euclidean coupling. The paper recalls that a gap uniform in lattice units forces trivial "
                 "continuum two-point functions~\\cite{bbll-forest-4d}, so no strong-coupling bound can by itself be a statement about a continuum or about the "
                 "geometry of Section~\\ref{sec:geometry}, and none of the numbers here is compared with a number there. One printed number does not follow from "
                 "the quantity the paper states for it: for $SU(3)$ at $\\beta_W=\\tfrac34$, $-\\tfrac12\\log0.9653=0.017658$ is printed as $0.0177$; the "
                 "model records it. The Zenodo record 10.5281/zenodo.23841532 that was proposed as a further reference could not be retrieved from this "
                 "environment and is not cited.\n")

    lines.append("\\section{Scale and blocking}\\label{sec:transmutation}\n"
                 "The ladder has one coupling, $g$, and refining or coarsening the lattice flows it. Any coarse-graining of spin networks, "
                 "whether block-spin or spin-foam, is a flow of couplings of this kind, so the library's account of the one-loop blocking "
                 "drift is kept as context. The mass gap of a spin model and the lattice scale of a gauge theory both fall as "
                 "$\\exp(-\\mathrm{slope}\\cdot\\beta)$, but they belong to different theories. Figure~\\ref{fig:transmutation} states the "
                 "identities within each theory and draws an arrow, not an equals sign, between the two. Nothing is computed here for a "
                 "gravitational flow. The equations are " + _names_line(tra.names) + ".\n")
    lines.append(ttex)

    lines.append(conclusion(pic))
    sel, cited = selected_section(by_key)
    lines.append("\\appendix\n\\section{Registered results of the ladder}\\label{sec:selected}\n"
                 "The ladder equations on which the paper rests; each is a registered sympy result whose checks passed in the build.\n" + sel)

    text = "\n".join(lines)
    keys = {k for m in re.findall(r"\\cite\w*\{([^}]*)\}", text) for k in m.split(",")}
    unknown = keys - set(kb.references())
    if unknown:
        raise QgError(f"cited but not in the knowledge base: {sorted(unknown)}")
    body.write_text(text, encoding="utf-8")
    numfile.write_text(numbers_tex(nums), encoding="utf-8")
    bib.write_text(kbdoc.bibliography(keys), encoding="utf-8")
    qg = [e for e in eqs if e.model in QG_MODELS]
    total = fig_checks + sum(len(e.verified) for e in qg)
    print(f"qgravity: {len(phys.equations | qgf.equations | tra.equations | ein.equations | pen.equations)} library equations, "
          f"{len(qg)} registered results, {total} checks passed (model {sum(len(e.verified) for e in qg)}, figures {fig_checks})")
    return [body, bib, numfile]
