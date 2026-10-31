"""The concise second document: physics and geometry, large generated equations, no tables.

Generates build/concise_body.tex and build/concise_bibliography.tex, which the hand-editable
wrapper concise_document.tex pulls in. Everything shown is built from the knowledge graph
(knowledgebase.py) and from the registered sympy models, in the manner of Figures 1 and 2 of
RosettaMath's rosettapaper.py:

    \\overbrace{quantity}^{name} = \\underbrace{expression}_{equation} = \\underbrace{...}

Each term of a row is the library equation *solved for the row's quantity*, and the build
substitutes that solution back into the equation and demands an identity (verify_solution).
Rows that also assert the terms are equal numbers (the geometry rows) are checked in sympy
under the stated assignment. A failed check raises, so nothing unverified is written.
"""
from __future__ import annotations

import re
from pathlib import Path

import sympy as sp

import eqlib
import knowledgebase as kb
from document import esc, cite, render_see


class ConciseError(AssertionError):
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

# Figure 2: geometry. A row lists its terms and the values of the letters that are not the row's
# quantity; the build checks that every term then equals the row's first term.
_S = sp.Symbol
GEOMETRY_ROWS = (
    ("inflation factor", "mu",
     ("Spectre inflation factor", "Spectre reciprocal pair", "Spectre linear-area relation",
      "Spectre growth rate"),
     {"g": 4 - sp.sqrt(15), "ell": (sp.sqrt(6) + sp.sqrt(10)) / 2, "gamma": sp.log(4 + sp.sqrt(15))},
     "Spectre"),
    ("inflation factor", "kappa",
     ("Hat inflation factor", "Hat constant as a golden power"),
     {"phi": (1 + sp.sqrt(5)) / 2},
     "hat"),
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

# Figure 4: Lie algebras, chirality and the Standard Model. Each row is one quantity stated by
# several library equations; its letters other than the row's own take the values written here
# (one concrete case: the SO(10) 27, three generations, the left-handed up quark), and every
# term must then equal the first. The arrows are links of the library, not equations.
CHIRALITY_ROWS = (
    ("half-spinor dimension", "S", ("Half-spinor dimension", "E6 fundamental branching", "Chirality splitting"),
     {"n": 10, "D_f": 27, "s": 32}, "$n=10$"),
    ("generation count", "N_gen", ("Net chiral class", "Weight vector demand"),
     {"chi_V": 3, "chi_F": 1, "W": 180, "F": 15}, "three generations"),
    ("Weyl fields per generation", "F", ("Weight vector demand", "Standard Model content of the 27"),
     {"W": 180, "N_gen": 3, "D_f": 27}, "the 27"),
    ("hypercharge", "Y", ("Hypercharge", "Gell-Mann-Nishijima formula"),
     {"q": 1, "Q": sp.Rational(2, 3), "T_3": sp.Rational(1, 2)}, "left-handed up quark"),
)
CHIRALITY_ARROWS = (
    ("Orthogonal algebra dimension", "Special unitary dimension"),
    ("Product representation dimension", "Standard Model content of the 27"),
    ("Weight vector demand", "Net chiral class"),
    ("Hypercharge normalisation", "Weak mixing at unification"),
)

# registered sympy results shown as single large displays, with one line of prose each
SELECTED = (
    ("Bloch's law", "271", ("bloch_theorem", "bloch_correction_corollary")),
    ("The Kosterlitz--Thouless flow", "216", ("kt_limiting_flow", "kt_pole_time", "kt_amplitude")),
    ("Asymptotic scaling of the $O(N)$ mass gap", "215", ("on_mass_af", "on_log_slope")),
    ("The Chair44 and Spectre tilings", "spectre-chair44", ("chair44_contact_formula", "spectre_window_dimension")),
    ("E8, the family grading and the weak mixing angle", "chirality-lisi",
     ("cl_e8_gradings", "cl_z3_z3_table", "cl_weak_mixing")),
    ("Yang--Mills asymptotic freedom", "", ("ym_b0", "ym_lattice_scale", "ym_log_slope")),
    ("The Hopf fibration and the Hopf insulator", "hopf-fibration",
     ("hf_invariant", "hf_pq", "hi_trim_formula", "hg_bpst_c2")),
)


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
        raise ConciseError(f"{eq.name}: solving for {sol.quantity} leaves residual {residual}")
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
                raise ConciseError(f"{a} and {b} are related by a correspondence, not an identity; "
                                   f"they cannot share the row for {quantity}")
        terms, head = [], None
        for name in names:
            eq = kb.equation(name)
            sol = kb.solve_for(eq, quantity)
            if sol is None:
                raise ConciseError(f"{name} cannot be solved for {quantity}")
            self.checks.append(verify_solution(eq, sol))
            self._note(eq)
            head = head or _head(sol, quantity)
            terms.append(_term(eq, sol))
        if len({str(t) for t in terms}) != len(terms):
            raise ConciseError(f"row {quantity} repeats a term")
        self.rows.append(f"{head} &= {' = '.join(terms)}{tail}")

    def geometry_row(self, quantity: str, letter: str, names, assign: dict, label: str):
        """A row whose terms are equal numbers once the other letters take the stated values."""
        terms, head, values = [], None, []
        for name in names:
            eq = kb.equation(name)
            sol = kb.solve_for(eq, quantity)
            if sol is None or sol.symbol.name != letter:
                raise ConciseError(f"{name}: expected {quantity} to be denoted {letter}")
            self.checks.append(verify_solution(eq, sol))
            self._note(eq)
            sub = {s: assign[s.name] for s in sol.expr.free_symbols if s.name in assign}
            missing = {s.name for s in sol.expr.free_symbols} - {s.name for s in sub}
            if missing:
                raise ConciseError(f"{name}: no value for {sorted(missing)}")
            values.append(sp.simplify(sp.radsimp(sol.expr.subs(sub))))
            head = head or _head(sol, f"{quantity}, {label}")
            terms.append(_term(eq, sol))
        first = values[0]
        for name, v in zip(names, values):
            if sp.simplify(sp.expand_log(v - first, force=True)) != 0 and abs(sp.N(v - first, 40)) > sp.Float("1e-30"):
                raise ConciseError(f"{label}: {name} gives {v}, not {first}")
        self.checks.append(f"{label}: all {len(names)} terms equal {first} under "
                           + ", ".join(f"{k} = {v}" for k, v in assign.items()))
        self.rows.append(f"{head} &= {' = '.join(terms)}")

    def arrow_row(self, a: str, b: str):
        """A correspondence: two library equations joined by the arrow their link names."""
        link = next(((k, n) for x, k, y, n in kb.library().links if (x, y) == (a, b)), None)
        if link is None:
            raise ConciseError(f"the library has no link from {a} to {b}")
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


def geometry_figure():
    d = Display()
    for quantity, letter, names, assign, label in GEOMETRY_ROWS:
        d.geometry_row(quantity, letter, names, assign, label)
    # the window dimension of the Spectre tiling, from the registered model, in the two bases
    mu = 4 + sp.sqrt(15)
    ell = (sp.sqrt(6) + sp.sqrt(10)) / 2
    dim = sp.log(5 + 2 * sp.sqrt(6)) / sp.log(mu)
    reg = eqlib.REGISTRY["spectre_window_dimension"].expr
    if sp.N(reg.rhs - dim, 40) != 0 and abs(sp.N(reg.rhs - dim, 40)) > sp.Float("1e-30"):
        raise ConciseError("spectre_window_dimension differs from log(5+2sqrt6)/log(4+sqrt15)")
    alt = 2 * sp.log(sp.sqrt(2) + sp.sqrt(3)) / (2 * sp.log(ell))
    if abs(sp.N(dim - alt, 40)) > sp.Float("1e-30"):
        raise ConciseError("window dimension: the (sqrt2+sqrt3) form differs")
    d.raw_row(r"D &= \underbrace{\frac{\ln(5+2\sqrt{6})}{\ln(4+\sqrt{15})}}_{\text{\tiny window dimension}} "
              r"= \underbrace{\frac{\ln(\sqrt{2}+\sqrt{3})}{\ln \ell}}_{\text{\tiny in the linear factor }"
              r"\ell=(\sqrt6+\sqrt{10})/2}",
              ["window dimension equals ln(5+2sqrt6)/ln(4+sqrt15) and ln(sqrt2+sqrt3)/ln(ell), to 40 digits"])
    d.raw_row(rf"{eqlib.REGISTRY['chair44_contact_formula'].latex().replace('=', '&=', 1)}")
    cap = ("The geometry of the library, conjoined. Top: the inflation factor $\\mu$ of the Spectre "
           "stated four ways, and the hat's $\\kappa=\\varphi^4$; these are numbers, and each term "
           "is checked equal when $g=4-\\sqrt{15}$, $\\ell=(\\sqrt6+\\sqrt{10})/2$, "
           "$\\gamma=\\ln\\mu$ and $\\varphi=(1+\\sqrt5)/2$. The Spectre and the hat are different "
           "tilings, so the two factors are not equated. Bottom: the window dimension, and the "
           "contact count of the Chair44 tiling.")
    return d, _figure("fig:geometry", d.tex(), cap, 0.9)


def chirality_figure():
    d = Display()
    for quantity, letter, names, assign, label in CHIRALITY_ROWS:
        d.geometry_row(quantity, letter, names, assign, label)
    # one more number: E8 counted by the two order-three gradings and the involution
    d.raw_row(r"\overbrace{\dim E_8}^{\text{\tiny three readings}} &= "
              r"\underbrace{86 + 2\cdot 81}_{\text{\tiny }\mathfrak e_6\oplus\mathfrak{su}(3)} = "
              r"\underbrace{120 + 128}_{\text{\tiny }\mathfrak{so}(16)} = "
              r"\underbrace{28 + 28 + 3\cdot 64 }_{\text{\tiny magic square}} = 248",
              ["86 + 2 x 81 = 120 + 128 = 28 + 28 + 3 x 64 = 248 (each decomposition also checked from the roots of E8 in "
               "models/chirality_lisi.py)"])
    if not (86 + 2 * 81 == 120 + 128 == 28 + 28 + 3 * 64 == 248):
        raise ConciseError("the three decompositions of E8 differ")
    for a, b in CHIRALITY_ARROWS:
        d.arrow_row(a, b)
    cap = ("Lie algebras, chirality and the Standard Model. Each identity row is one quantity stated by several "
           "library equations, checked equal in one concrete case written under its name; the last identity "
           "counts $E_8$ three ways, each decomposition rebuilt from the roots in the companion model. The "
           "arrows are links of the library, not equations: the dimensions of $\\mathfrak{so}(n)$ and "
           "$\\mathfrak{su}(N)$ agree at $\\mathfrak{su}(2)=\\mathfrak{so}(3)$ and $\\mathfrak{su}(4)=\\mathfrak{so}(6)$ "
           "and by coincidence elsewhere, and the dimension count for three generations "
           "($W=4\\cdot15\\cdot3=180$ against the $128$ of an eigenspace) bounds the generation count from a "
           "different side than the net chiral class does.")
    return d, _figure("fig:chirality", d.tex(), cap, 0.98)


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
        raise ConciseError("cannot solve the lattice coupling for g")
    exponent = (1 / (2 * b0 * g ** 2)).subs(g, g2[0])
    ym = sp.simplify(exponent - P(r_g.expr) * beta)
    if ym != 0:
        raise ConciseError(f"Yang-Mills exponent differs from r_g*beta by {ym}")
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
def selected_section(eqs_by_key) -> tuple[str, set[str]]:
    parts, keys = [], set()
    for title, entry, ekeys in SELECTED:
        parts.append(f"\\subsection*{{{title}}}")
        for k in ekeys:
            eq = eqs_by_key[k]
            cites = kb.citations_for(eq)
            keys.update(cites)
            parts.append(f"\\noindent {esc(eq.title)}{cite(cites)}."
                         f"\n\\begin{{equation}}\\label{{ceq:{eq.key}}}\n  {eq.latex(multiline=True)}\n\\end{{equation}}")
            if eq.see:
                parts.append(f"\\noindent\\textit{{See}} {render_see(eq.see)}.")
    return "\n".join(parts) + "\n", keys



def numbers(eqs, explored) -> dict:
    """The counts the abstract and introduction quote, computed at build time."""
    import pkgutil
    import models
    lib = kb.library()
    return {
        "nEq": len(lib.equations), "nQty": len(lib.quantities), "nLinks": len(lib.links),
        "nJoins": len(kb.all_joins()), "nFamilies": len(kb.families(2)),
        "nModels": len(list(pkgutil.iter_modules(models.__path__))),
        "nRegEq": len(eqs), "nChecks": len(eqlib.checks()),
        "nExplore": len(explored), "nExploreChecks": sum(len(x.checks) for x in explored),
    }


def numbers_tex(nums: dict) -> str:
    return ("% Generated by sympy_models/concise.py. Do not edit by hand.\n"
            + "".join(f"\\newcommand{{\\{k}}}{{{v}}}\n" for k, v in nums.items()))


def introduction(nums: dict, figure_checks: int) -> str:
    return (
        "\\section{Introduction}\\label{sec:intro}\n"
        "An equation states a relation between quantities, but it is written in letters, and a letter is not a "
        "quantity: $E$ is the energy in $E=mc^2$ and the electric field in another equation. Two equations can be "
        "compared only after each letter has been read as the quantity it denotes. RosettaMath~\\cite{rosettamath} "
        "makes that reading explicit, so that equations stating one quantity can be solved for it and set equal in a "
        "single large display. We apply the method, with sympy doing the algebra, to the physics and geometry of the "
        "OpenAI math release and of the fork's companion repositories, and use the result to ask three questions.\n"
        "\\begin{itemize}\\itemsep2pt\n"
        "\\item \\textbf{Why is the hat's expansion ratio $\\varphi^4$?} Romanov's observation, that the hat "
        "inflation factor equals the fourth power of the golden ratio, and that the same number is the ratio of "
        "unreflected to reflected hats, is traced to its algebra and to a relation between two tables of the "
        "substitution (Section~\\ref{sec:hat}). The Spectre, by contrast, is not golden.\n"
        "\\item \\textbf{Where would a Yang--Mills mass gap sit?} The release contains no gauge-theory result, so we "
        "do not prove one. We write down the chain a gap would have to satisfy, check its consistency, and keep the "
        "$O(N)$ analogy as a correspondence (Section~\\ref{sec:ym}).\n"
        "\\item \\textbf{What does $E_8$ say about the Standard Model's matter?} The chirality repository's "
        "decompositions are rebuilt from the Cartan matrix and tied to the $16$, to $\\sin^2\\theta_W=3/8$ and to a "
        "dimension count for generations (Sections~\\ref{sec:chirality} and~\\ref{sec:explore}).\n"
        "\\end{itemize}\n"
        "Every term of every displayed row is the source equation solved for the row's quantity, substituted back and "
        "required to give an identity, and every number is checked; \\nChecks{} checks over the registered models and "
        f"{figure_checks + nums['nExploreChecks']} more over the figures and explorations pass in the build, and a failure stops it. "
        "The checks test that stated formulas agree with one another and with independent computation; they do not "
        "re-prove any theorem of the release. Identities are kept apart from correspondences, which are drawn as "
        "arrows (Section~\\ref{sec:physics}).\n")


DEDICATED = ("hat_golden", "ym_mass_gap")        # get sections of their own
OMITTED = ("kb_paths",)                          # graph internals: run by the build, not described in the paper


def dedicated_section(x, label: str, title: str, lead: str, tail: str) -> str:
    return (f"\\section{{{title}}}\\label{{{label}}}\n{lead}\n\n{x.prose}\n\n"
            f"\\begin{{equation}}\\label{{eq:explore-{x.name}}}\n  {x.display}\n\\end{{equation}}\n\n{tail}\n"
            f"\\noindent\\textit{{{len(x.checks)} checks passed (\\texttt{{explore/{esc(x.name)}.py}}).}}\n")


def explorations_section(explored) -> str:
    parts = ["\\section{Further explorations}\\label{sec:explore}\n"
             "Each subsection is a script in \\texttt{sympy\\_models/explore/}, run by the build; its checks are "
             "listed in the script's own output.\n"]
    for x in explored:
        if x.name in DEDICATED or x.name in OMITTED:
            continue
        parts.append(f"\\subsection*{{{x.title.replace('&', chr(92) + '&')}}}\n{x.prose}\n")
        if x.display:
            parts.append(f"\\begin{{equation}}\\label{{eq:explore-{x.name}}}\n  {x.display}\n\\end{{equation}}\n")
        parts.append(f"\\noindent\\textit{{{len(x.checks)} checks passed (\\texttt{{explore/{esc(x.name)}.py}}).}}\n")
    return "\n".join(parts)


def conclusion() -> str:
    return (
        "\\section{Conclusion}\\label{sec:conclusion}\n"
        "Romanov's surprise has an answer at the level of algebra. The hat's substitution has spectrum "
        "$\\{1,-1,\\varphi^4,\\varphi^{-4}\\}$, because $x^2-7x+1$ is the minimal polynomial of $\\varphi^4$ and $7$ is a "
        "Lucas number; the linear scale is $\\varphi^2$; and the chirality ratio equals the inflation factor because the "
        "unreflected hats of a metatile are counted by the $H$ row of the substitution while the one reflected hat sits in "
        "each $H$. What remains unexplained is geometric: why the hat's substitution has these entries. The Spectre's "
        "$4+\\sqrt{15}$ is not golden, and the golden $2+\\sqrt5=\\varphi^3$ of its perimeter is a separate fact; sharing "
        "$\\mathbb{Q}(\\sqrt5)$ is not a link between them.\n\n"
        "For Yang--Mills we have a frame and not a result. A coupling-independent physical gap is the statement that "
        "$m_L/(a\\Lambda_L)$ tends to a constant, the gap is that constant times $\\Lambda_L$, and its positivity is the "
        "open problem. The proved two-dimensional gap and the $O(N)$ asymptotics share the linear drift of the "
        "coupling under blocking, and we have kept them as a correspondence. For $E_8$ the checks confirm the "
        "decompositions behind the chirality theorem, show that a weight-vector count does not fit three generations "
        "in $128$ weights, and reproduce $\\sin^2\\theta_W=3/8$ at unification.\n\n"
        "Three limits apply throughout. The checks test internal consistency, not the truth of cited theorems; "
        "several inputs are quoted from repositories, not derived; and the one-loop numbers are illustrations. The "
        "natural next steps are a geometric derivation of the hat's substitution entries from its tile, a check of "
        "whether the Spectre's boundary operator has a closed form, and, if the expert group's Yang--Mills write-up "
        "appears, a comparison of its statements with the chain of Section~\\ref{sec:ym}.\n")


def write(out_dir: Path, eqs) -> list[Path]:
    """Build the concise body and its bibliography; raises if any check fails."""
    import kbdoc
    by_key = {e.key: e for e in eqs}
    phys, ftex = physics_figure()
    geo, gtex = geometry_figure()
    tra, ttex = transmutation_figure()
    chi, ctex = chirality_figure()
    import explore
    explored = explore.run_all()
    nums = numbers(eqs, explored)
    fig_checks = sum(len(d.checks) for d in (phys, geo, tra, chi))

    lines = ["% Generated by sympy_models/concise.py. Do not edit by hand.\n"]
    lines.append(introduction(nums, fig_checks))
    lines.append("\\section{The library, conjoined}\\label{sec:physics}\n"
                 "Every equation of the knowledge graph states some quantity, and equations that state "
                 "the same quantity can be set equal once each is solved for it. "
                 f"Figure~\\ref{{fig:physics}} does this for {len(PHYSICS_ROWS)} quantities at once. "
                 "The equations are "
                 + _names_line(phys.names) + ".\n")
    lines.append(ftex)
    lines.append("\\section{Geometry}\\label{sec:geometry}\n"
                 "The inflation factors of the aperiodic tilings are numbers, so their rows are checked "
                 "equal in sympy, not only solved for. The equations are "
                 + _names_line(geo.names) + ".\n")
    lines.append(gtex)
    byname = {x.name: x for x in explored}
    lines.append(dedicated_section(
        byname["hat_golden"], "sec:hat", "The golden hat",
        "Of the two inflation factors in Figure~\\ref{fig:geometry}, the hat's is the surprising one: it is a power of "
        "the golden ratio, a number the Penrose tilings made familiar, in a tiling that has nothing else of Penrose's. "
        "We take the four-metatile substitution $M$ of the hat (quoted from the CICY repository) and derive what follows "
        "from it exactly.",
        "Two things stay open. First, the entries of $M$ and the hat counts are quoted, and the equality of the "
        "unreflected-hat vector with the $H$ row is a relation we verify but do not explain geometrically. Second, one "
        "published statement is too strong: the Perron frequencies are not all irrational, since $f_H=1/3$ exactly; "
        "the frequencies of $T$, $P$ and $F$ lie in $\\mathbb{Q}(\\sqrt5)\\setminus\\mathbb{Q}$, and one irrational "
        "frequency is already enough to forbid a period."))
    lines.append("\\section{Transmutation}\\label{sec:transmutation}\n"
                 "The mass gap of a spin model and the lattice scale of a gauge theory both fall as "
                 "$\\exp(-\\mathrm{slope}\\cdot\\beta)$, but they belong to different theories. "
                 "Figure~\\ref{fig:transmutation} states the identities within each theory and draws an "
                 "arrow, not an equals sign, between the two. "
                 "The equations are " + _names_line(tra.names) + ".\n")
    lines.append(ttex)
    lines.append(dedicated_section(
        byname["ym_mass_gap"], "sec:ym", "The Yang--Mills mass gap",
        "The Clay problem asks for a quantum $SU(N)$ Yang--Mills theory on $\\mathbb{R}^4$ whose spectrum has a gap. "
        "No manuscript of the release treats a gauge theory, and reports of a four-dimensional existence claim "
        "mention no mass gap, so this section proves nothing about it. It asks what the library can state and check "
        "around it.",
        "What is established, then, is the textbook scaling, the consistency condition, and the proved two-dimensional "
        "analogues: a positive gap for the $O(3)$ lattice model, and exact asymptotics of the mass for $O(4)$. What is not "
        "is the existence of the limit $c=\\lim m_L/(a\\Lambda_L)$ for a gauge theory, its positivity, and any "
        "implication from the $O(N)$ results to $SU(N)$; the shared linear drift is a reason to look, not evidence of a "
        "theorem."))
    lines.append("\\section{Lie algebras and chirality}\\label{sec:chirality}\n"
                 "The companion repository of the chirality theorem puts the Standard Model in the $27$ of "
                 "$E_6$ inside $E_8$. Figure~\\ref{fig:chirality} equates what the library states about the "
                 "half-spinor, the generation count, the fifteen Weyl fields and the hypercharge, each "
                 "in one concrete case, and draws the links between dimension counts that are not equalities. "
                 "The equations are " + _names_line(chi.names) + ".\n")
    lines.append(ctex)
    lines.append(explorations_section(explored))
    lines.append(conclusion())
    sel, cited = selected_section(by_key)
    lines.append("\\appendix\n\\section{Registered results}\\label{sec:selected}\n"
                 "The registered models behind the figures; each is a sympy expression whose checks "
                 "passed in the full equations document.\n" + sel)

    keys = {"rosettamath"} | cited
    for d in (phys, geo, tra, chi):
        for n in d.equations:
            ref = kb.reference_for_source(kb.equation(n).source)
            if ref:
                keys.add(ref)
    body = out_dir / "concise_body.tex"
    bib = out_dir / "concise_bibliography.tex"
    body.write_text("\n".join(lines), encoding="utf-8")
    numfile = out_dir / "concise_numbers.tex"
    numfile.write_text(numbers_tex(nums), encoding="utf-8")
    bib.write_text(kbdoc.bibliography(keys), encoding="utf-8")
    total = fig_checks + sum(len(x.checks) for x in explored)
    print(f"concise: {len(phys.equations | geo.equations | tra.equations | chi.equations)} library equations, {total} checks passed")
    return [body, bib, numfile]
