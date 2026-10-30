"""A knowledge graph of equations, with sympy doing the algebra.

The design is RosettaMath's (https://github.com/brentharts/RosettaMath, rosettaphys.py,
MIT licence, same author), and the data at the bottom of this file is copied from its
knowledge base: its equations, its physical quantities, its per-equation "readings" and
its hand-stated links. What is different here is the algebra. RosettaMath reads LaTeX
with a purpose-built term type; this module reads the same LaTeX into sympy expressions,
so rearranging, eliminating and checking an equation is done by sympy, and the results
are registered as ordinary `eqlib.Equation` records for the equations document.

The model, in four sentences
----------------------------
Nodes are equations, physical quantities (energy, microstate count, ...) and repository
entries. An equation is joined to the quantities its letters denote by a DENOTES edge;
that table (`reads`) is the point of the whole design, because a glyph is not a quantity:
the W of Boltzmann's S = k_B log W counts microstates, the W of W = F d is an energy, and
joining them on the letter gives a false statement that parses and prints beautifully.
Two equations may be joined only on a quantity that both declare, and only if every other
letter they share is read the same way by both. Silence is not agreement: an equation
without a readings table takes part in no join.

What this module does
---------------------
    kb = graph()                      the graph (built once, then cached)
    kb.census()                       node and edge counts by kind
    equation("Ideal gas law")         a KBEquation: .latex, .reads, .relation(), ...
    solve_for(eq, "temperature")      eq rearranged for the symbol that denotes a quantity
    join(a, b)                        two equations eliminated on a shared quantity
    all_joins(), families()           every pairwise join; chains over one quantity
    collisions()                      letters used for two different quantities
    dimension_report()                dimensional homogeneity of every readable equation
    refusals()                        which equations the reader will not rearrange, and why
    selftest()                        structural checks of the tables; raises on failure

What it refuses, on purpose
---------------------------
  * an equation with a derivative, a sum, an integral, a trace or an inequality is kept in
    the graph (with its readings) but is not rearranged;
  * a quantity that occurs twice, or an even power (x**2 = 4 has two roots), is not solved
    for: sympy returns more than one solution and the join is refused;
  * a letter used twice in one equation for the same quantity (m_1, m_2) is not a pivot.

Nothing here proves physics. A join is the statement that two quoted equations, taken
together, imply a third by transitivity; sympy checks the algebra and the dimensions, and
the quotations are only as good as their sources. `Join.notes` says which assumptions
were made (the same letter, read the same way, is taken to be the same variable).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from functools import lru_cache

import sympy as sp

# --------------------------------------------------------------------------------------
# edge vocabulary
# --------------------------------------------------------------------------------------
DENOTES = "denotes"          # equation -> the quantity one of its letters means
SHARES = "shares"            # equation <-> equation, with a quantity in common
CITES = "cites"              # repository entry -> equation (the fork's tools/ scripts)
MODELS = "models"            # eqlib equation -> knowledge-base equation it cross-checks
SPECIALISES = "specialises"  # a special case of a more general law (stated)
LIMIT_OF = "limit-of"        # what a law becomes in some limit (stated)
DEFINES = "defines"          # one equation giving another its meaning (stated)
FIELD = "field"              # equations of the same kind (stated)

EDGE_NOTES = {
    DENOTES: "equation to the quantity one of its letters means",
    SHARES: "two equations that denote a quantity in common",
    CITES: "repository entry to a knowledge-base equation drawn from its script",
    MODELS: "sympy-modelled equation to the knowledge-base equation it cross-checks",
    SPECIALISES: "a special case of a more general law",
    DEFINES: "one equation giving another its meaning",
    LIMIT_OF: "what a law becomes in some limit",
    FIELD: "equations of the same kind",
}

# Quantities too common to be a useful pivot: joining two equations on "mass" says only
# that both have a mass in them, and the two masses are rarely the same body.
COMMON_QUANTITIES = {"pi", "mass", "length", "time", "dimensionless", "amount", "area", "volume"}

# Scripts of this repository that the knowledge base cites, and their CONTENTS.md entries.
SCRIPT_ENTRY = {
    "tools/on_model_mc.py": "215",
    "tools/kt_flow.py": "216",
    "tools/haldane_gap.py": "268",
    "tools/bloch_law.py": "271",
}


class Unreadable(Exception):
    """The LaTeX is outside what the reader will turn into algebra; the message says why."""


# --------------------------------------------------------------------------------------
# reading LaTeX into sympy
# --------------------------------------------------------------------------------------
_REFUSE = [
    (r"\\partial|\\nabla|\\Box|\\square|\\mathcal\{L\}", "derivative or field operator"),
    (r"\\frac\{d[A-Za-z]*\}\{d", "derivative"),
    (r"\\sum|\\prod|\\int|\\oint|\\lim", "sum, product, integral or limit"),
    (r"\\Re|\\Im|\\Tr|\\det|\\mathrm|\\operatorname|\\text|\\dim", "operator or upright text"),
    (r"\\to\b|\\mapsto|\\rightarrow|\\leftarrow", "map or arrow"),
    (r"\\dagger|\\hat|\\vec|\\bar|\\tilde|\\mathbf|\\mathbb|\\cdots|\\ldots|\\dots|\\vdots", "decoration or ellipsis"),
    (r"\\le\b|\\ge\b|\\leq|\\geq|\\approx|\\sim|\\propto|\\equiv|\\neq|\\ne\b|<|>", "not an equation"),
    (r"'", "prime"),
    (r"\\infty", "infinity"),
    (r"\\pm|\\mp", "plus or minus"),
]

_GREEK = {
    "alpha", "beta", "gamma", "delta", "epsilon", "varepsilon", "zeta", "eta", "theta", "vartheta",
    "iota", "kappa", "lambda", "mu", "nu", "xi", "rho", "sigma", "tau", "upsilon", "phi", "varphi",
    "chi", "psi", "omega", "hbar", "ell", "Gamma", "Delta", "Theta", "Lambda", "Xi", "Sigma",
    "Phi", "Psi", "Omega",
}
_FUNCS = {
    "exp": sp.exp, "log": sp.log, "ln": sp.log, "sin": sp.sin, "cos": sp.cos, "tan": sp.tan,
    "sinh": sp.sinh, "cosh": sp.cosh, "tanh": sp.tanh,
}
_TOKEN = re.compile(r"\\[A-Za-z]+|\\.|\d+(?:\.\d+)?|[A-Za-z]|[{}()\[\]^_+\-*/=,|]")
_SPACING = {"\\,", "\\;", "\\!", "\\:", "\\ ", "\\quad", "\\qquad", "\\left", "\\right"}


def _tokens(latex: str) -> list[str]:
    out = []
    for t in _TOKEN.findall(latex):
        if t in _SPACING:
            continue
        out.append(t)
    return out


class _Reader:
    """Recursive-descent reader for the scalar algebra the knowledge base writes."""

    def __init__(self, tokens: list[str]):
        self.t = tokens
        self.i = 0

    # -- plumbing
    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else None

    def take(self, expect: str | None = None):
        tok = self.peek()
        if tok is None or (expect is not None and tok != expect):
            raise Unreadable(f"expected {expect!r}, found {tok!r}")
        self.i += 1
        return tok

    def done(self) -> bool:
        return self.i >= len(self.t)

    # -- grammar
    def relation(self):
        lhs = self.expr()
        self.take("=")
        rhs = self.expr()
        if not self.done():
            raise Unreadable(f"unread tokens from {self.peek()!r}")
        return sp.Eq(lhs, rhs, evaluate=False)

    def expr(self):
        terms = [self.term()]
        while self.peek() in ("+", "-"):
            sign = self.take()
            nxt = self.term()
            terms.append(nxt if sign == "+" else -nxt)
        return sp.Add(*terms)

    _STOP = {None, "+", "-", "=", ")", "}", "]", "^", "_", ",", "/", "*", "|"}

    def term(self):
        value = self.unary()
        while True:
            tok = self.peek()
            if tok in ("*", "\\cdot", "\\times"):
                self.take()
                value = value * self.unary()
            elif tok == "/":
                self.take()
                value = value / self.unary()
            elif tok not in self._STOP and tok not in ("\\cdot", "\\times"):
                value = value * self.power()      # juxtaposition is multiplication
            else:
                return value

    def unary(self):
        if self.peek() == "-":
            self.take()
            return -self.unary()
        if self.peek() == "+":
            self.take()
            return self.unary()
        return self.power()

    def power(self):
        base = self.atom()
        if self.peek() == "^":
            self.take()
            return base ** self.group()
        return base

    def group(self):
        """A braced expression, or one token (the argument of ^)."""
        if self.peek() == "{":
            self.take("{")
            value = self.expr()
            self.take("}")
            return value
        return self.atom()

    def atom(self):
        tok = self.peek()
        if tok is None:
            raise Unreadable("unexpected end")
        if re.fullmatch(r"\d+(?:\.\d+)?", tok):
            self.take()
            return sp.Rational(tok)
        if tok in ("(", "["):
            close = ")" if tok == "(" else "]"
            self.take()
            value = self.expr()
            self.take(close)
            return value
        if tok == "{":
            return self.group()
        if tok == "\\frac":
            self.take()
            num, den = self.group(), self.group()
            return num / den
        if tok == "\\sqrt":
            self.take()
            root = 2
            if self.peek() == "[":
                self.take("[")
                root = self.expr()
                self.take("]")
            return self.group() ** (sp.Integer(1) / root)
        if tok == "\\pi":
            self.take()
            return sp.pi
        if tok.startswith("\\") and tok[1:] in _FUNCS:
            self.take()
            func = _FUNCS[tok[1:]]
            arg = self.group() if self.peek() != "(" else self._paren()
            return func(arg)
        if tok.startswith("\\") and tok[1:] in _GREEK:
            self.take()
            return self._symbol(tok[1:])
        if len(tok) == 1 and tok.isalpha():
            self.take()
            return self._symbol(tok)
        raise Unreadable(f"cannot read {tok!r}")

    def _paren(self):
        self.take("(")
        value = self.expr()
        self.take(")")
        return value

    def _symbol(self, base: str):
        name = base
        if self.peek() == "_":
            self.take()
            if self.peek() == "{":
                self.take("{")
                parts = []
                while self.peek() != "}":
                    parts.append(self.take())
                self.take("}")
            else:
                parts = [self.take()]
            sub = "".join(p.lstrip("\\") for p in parts)
            name = f"{base}_{sub}"
        return sp.Symbol(name)


def read_relation(latex: str) -> sp.Eq:
    """The equation as `sp.Eq(lhs, rhs, evaluate=False)`; raises Unreadable with the reason."""
    for pattern, why in _REFUSE:
        if re.search(pattern, latex):
            raise Unreadable(why)
    if latex.count("=") != 1:
        raise Unreadable("not a single equation")
    return _Reader(_tokens(latex)).relation()


def read_rate(latex: str):
    """Read 'dx/dl = rhs' written as \\frac{dx}{dl} = rhs into (x, l, rhs) with sympy symbols.

    The join machinery refuses derivatives, but the flows of the Kosterlitz-Thouless
    recursion are exactly that, and cross-checking them against a model's explicit solution
    needs the right-hand side as an expression."""
    m = re.fullmatch(r"\s*\\frac\{d([A-Za-z])\}\{d([A-Za-z])\}\s*=\s*(.*)", latex)
    if not m:
        raise Unreadable("not a rate equation dx/dl = ...")
    reader = _Reader(_tokens(m.group(3)))
    rhs = reader.expr()
    if not reader.done():
        raise Unreadable(f"unread tokens from {reader.peek()!r}")
    return sp.Symbol(m.group(1)), sp.Symbol(m.group(2)), rhs


# --------------------------------------------------------------------------------------
# nodes, edges, graph
# --------------------------------------------------------------------------------------
@dataclass(eq=False)
class Node:
    kind: str                 # "equation" | "quantity" | "entry"
    key: str
    label: str
    data: dict = field(default_factory=dict)
    out: list = field(default_factory=list)
    inn: list = field(default_factory=list)

    def __repr__(self):
        return f"<{self.kind} {self.label}>"


@dataclass(eq=False)
class Edge:
    kind: str
    src: Node
    dst: Node
    note: str = ""


class Graph:
    """Nodes and typed edges, with the few walks the rest of the module needs."""

    def __init__(self):
        self.nodes: dict[tuple[str, str], Node] = {}
        self.edges: list[Edge] = []

    def add(self, node: Node) -> Node:
        self.nodes[(node.kind, node.key)] = node
        return node

    def get(self, kind: str, key: str) -> Node | None:
        return self.nodes.get((kind, key))

    def of_kind(self, kind: str) -> list[Node]:
        return [n for n in self.nodes.values() if n.kind == kind]

    def link(self, kind: str, src: Node | None, dst: Node | None, note: str = "") -> Edge | None:
        if src is None or dst is None or src is dst:
            return None
        edge = Edge(kind, src, dst, note)
        self.edges.append(edge)
        src.out.append(edge)
        dst.inn.append(edge)
        return edge

    def neighbours(self, node: Node, kinds=None) -> list[Node]:
        seen, out = set(), []
        for edge in node.out + node.inn:
            if kinds and edge.kind not in kinds:
                continue
            other = edge.dst if edge.src is node else edge.src
            if other is not node and id(other) not in seen:
                seen.add(id(other))
                out.append(other)
        return out

    def path(self, start: Node, goal: Node, kinds=None, limit: int = 6):
        """A shortest chain of nodes from one to another, or None (breadth first)."""
        if start is goal:
            return [start]
        prev = {id(start): None}
        frontier = [start]
        for _ in range(limit):
            nxt = []
            for node in frontier:
                for other in self.neighbours(node, kinds):
                    if id(other) in prev:
                        continue
                    prev[id(other)] = node
                    if other is goal:
                        chain = [other]
                        while prev[id(chain[-1])] is not None:
                            chain.append(prev[id(chain[-1])])
                        return chain[::-1]
                    nxt.append(other)
            frontier = nxt
            if not frontier:
                break
        return None

    def census(self) -> dict:
        kinds: dict[str, int] = {}
        for edge in self.edges:
            kinds[edge.kind] = kinds.get(edge.kind, 0) + 1
        return kinds

    def __repr__(self):
        return f"<Graph {len(self.nodes)} nodes, {len(self.edges)} edges>"


# --------------------------------------------------------------------------------------
# equations and their readings
# --------------------------------------------------------------------------------------
@dataclass(eq=False)
class KBEquation:
    name: str
    field: str
    latex: str
    source: str | None
    slug: str | None
    blurb: str | None
    reads: dict
    node: Node | None = None
    _relation: object = field(default=None, repr=False)
    refusal: str | None = None            # why the reader will not rearrange it, or None

    def __post_init__(self):
        if not self.reads:
            self.refusal = "no readings declared"
            return
        try:
            self._relation = read_relation(self.latex)
        except Unreadable as why:
            self.refusal = str(why)
            return
        letters = {s.name for s in self._relation.free_symbols}
        declared = set(self.reads) - {"pi"}
        if letters != declared:
            self.refusal = ("letters and readings disagree: " + ", ".join(sorted(letters ^ declared)))
            self._relation = None

    @property
    def label(self) -> str:
        return self.name

    def relation(self) -> sp.Eq | None:
        """The equation as sympy, or None if the reader refused it."""
        return self._relation

    @property
    def entry(self) -> str:
        """CONTENTS.md entry of the script this equation was drawn from, or ''."""
        m = re.search(r"(tools/[A-Za-z0-9_]+\.py)", self.source or "")
        return SCRIPT_ENTRY.get(m.group(1), "") if m else ""

    def symbol_for(self, quantity: str) -> sp.Symbol | None:
        """The one symbol that denotes `quantity`, or None if none or more than one does."""
        hits = [s for s, q in self.reads.items() if q == quantity and s != "pi"]
        return sp.Symbol(hits[0]) if len(hits) == 1 else None

    def quantities(self) -> list[str]:
        return sorted({q for s, q in self.reads.items() if s != "pi"})


@dataclass(frozen=True)
class Quantity:
    name: str
    dimension: str
    blurb: str
    slug: str = ""


_BASE = sp.symbols("M L T I K N", positive=True)       # mass, length, time, current, temperature, amount


def _dimension_of(spec: str) -> sp.Expr:
    """'M L^2 T^-2' -> M*L**2/T**2 ; '1' -> 1."""
    out = sp.Integer(1)
    names = {str(b): b for b in _BASE}
    for part in spec.split():
        if part == "1":
            continue
        base, _, power = part.partition("^")
        if base not in names:
            raise ValueError(f"unknown base dimension {base!r} in {spec!r}")
        out *= names[base] ** (sp.Rational(power) if power else 1)
    return out


# --------------------------------------------------------------------------------------
# the library: built once from the data at the bottom of this file
# --------------------------------------------------------------------------------------
class Library:
    def __init__(self):
        self.quantities = {q[0]: Quantity(*q[:4]) for q in QUANTITIES_DATA}
        self.equations: dict[str, KBEquation] = {}
        for name, fld, latex, source, slug, blurb in EQUATIONS_DATA:
            self.equations[name] = KBEquation(name, fld, latex, source, slug, blurb,
                                              dict(READINGS_DATA.get(name, {})))
        self.links = [(a, kind, b, note) for a, kind, b, note in LINKS_DATA]
        self.graph = self._build()
        self._solutions: dict = {}

    def _build(self) -> Graph:
        g = Graph()
        for q in self.quantities.values():
            g.add(Node("quantity", q.name, q.name, {"quantity": q}))
        for eq in self.equations.values():
            eq.node = g.add(Node("equation", eq.name, eq.name, {"equation": eq}))
        for entry in sorted(set(SCRIPT_ENTRY.values())):
            g.add(Node("entry", entry, f"entry {entry}"))

        for eq in self.equations.values():
            for letter, qname in sorted(eq.reads.items()):
                g.link(DENOTES, eq.node, g.get("quantity", qname), letter)
            if eq.entry:
                g.link(CITES, g.get("entry", eq.entry), eq.node, eq.source or "")

        # equations that denote a telling quantity in common
        by_quantity: dict[str, list[KBEquation]] = {}
        for eq in self.equations.values():
            for qname in eq.quantities():
                by_quantity.setdefault(qname, []).append(eq)
        self.rarity = {q: len(v) for q, v in by_quantity.items()}
        seen = set()
        for qname, members in by_quantity.items():
            if qname in COMMON_QUANTITIES:
                continue
            for i, a in enumerate(members):
                for b in members[i + 1:]:
                    pair = (a.name, b.name)
                    if pair in seen:
                        continue
                    seen.add(pair)
                    # note the rarest shared quantity, as RosettaMath does
                    g.link(SHARES, a.node, b.node, qname)

        for a, kind, b, note in self.links:
            g.link(kind.lower().replace("_", "-"), self.equations[a].node, self.equations[b].node, note)
        return g


@lru_cache(maxsize=1)
def library() -> Library:
    return Library()


def graph() -> Graph:
    return library().graph


def equation(name: str) -> KBEquation:
    try:
        return library().equations[name]
    except KeyError:
        raise KeyError(f"no equation named {name!r} in the knowledge base") from None


def equations(field_name: str | None = None) -> list[KBEquation]:
    return [e for e in library().equations.values() if field_name is None or e.field == field_name]


def _as_eq(what) -> KBEquation:
    return what if isinstance(what, KBEquation) else equation(what)


# --------------------------------------------------------------------------------------
# readings: agreement and collisions
# --------------------------------------------------------------------------------------
def agree(a, b, quantity: str) -> bool:
    """Do both equations denote `quantity` with a single, declared letter each?"""
    a, b = _as_eq(a), _as_eq(b)
    return a.symbol_for(quantity) is not None and b.symbol_for(quantity) is not None


def clashes(a, b, quantity: str) -> list[str]:
    """Letters the two equations share and read as different quantities.

    A join states both equations at once, so a shared letter is one variable; two equations
    that agree about the pivot but read some other letter differently (an amount of
    substance and a number density, say) would be made to agree about that too, which
    neither says. Letters one equation does not declare do not count: silence is not a
    reading. (`quantity` is accepted for symmetry with agree().)"""
    a, b = _as_eq(a), _as_eq(b)
    return sorted(s for s in set(a.reads) & set(b.reads) if s != "pi" and a.reads[s] != b.reads[s])


def collisions() -> list[tuple[KBEquation, KBEquation, str, str, str]]:
    """Every letter two equations share but read as different quantities.

    The output is a list of places where physics reuses a character, generated from the
    declared readings rather than remembered. Each distinct (letter, pair of readings) is
    reported once, with the first two equations that show it."""
    eqs = list(library().equations.values())
    seen, out = set(), []
    for i, a in enumerate(eqs):
        for b in eqs[i + 1:]:
            for letter in set(a.reads) & set(b.reads):
                left, right = a.reads[letter], b.reads[letter]
                if letter == "pi" or left == right:
                    continue
                key = (letter, tuple(sorted((left, right))))
                if key not in seen:
                    seen.add(key)
                    out.append((a, b, letter, left, right))
    return sorted(out, key=lambda r: (r[2].lower(), r[3]))


# --------------------------------------------------------------------------------------
# solving, joining, chaining
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Solution:
    quantity: str
    symbol: sp.Symbol
    expr: sp.Expr
    direct: bool            # True when the equation states it alone on one side as quoted


def solve_for(what, quantity: str) -> Solution | None:
    """The equation rearranged for the letter that denotes `quantity`, or None.

    None means the reader refused the equation, the quantity is denoted by no letter or by
    more than one, or sympy finds other than exactly one root. An even power therefore
    refuses: x**2 = 4 has two roots and a chain built on one of them is a false statement."""
    eq = _as_eq(what)
    lib = library()
    key = (eq.name, quantity)
    if key in lib._solutions:
        return lib._solutions[key]
    result = None
    rel = eq.relation()
    sym = eq.symbol_for(quantity)
    if rel is not None and sym is not None:
        try:
            roots = sp.solve(rel, sym, dict=False)
        except (NotImplementedError, ValueError):
            roots = []
        if len(roots) == 1 and not roots[0].has(sym):
            direct = (rel.lhs == sym and not rel.rhs.has(sym)) or (rel.rhs == sym and not rel.lhs.has(sym))
            result = Solution(quantity, sym, sp.simplify(roots[0]) if not direct
                              else (rel.rhs if rel.lhs == sym else rel.lhs), direct)
    lib._solutions[key] = result
    return result


@dataclass
class Join:
    """Two equations eliminated on a shared quantity: `ea = eb`, both sides free of the pivot."""
    a: KBEquation
    b: KBEquation
    pivot: str
    ea: sp.Expr
    eb: sp.Expr
    rearranged: bool
    notes: str = ""

    @property
    def statement(self) -> sp.Eq:
        return sp.Eq(self.ea, self.eb, evaluate=False)

    @property
    def entry(self) -> str:
        """The CONTENTS entry both equations point at (or the only one that does), else ''."""
        entries = {self.a.entry, self.b.entry} - {""}
        return entries.pop() if len(entries) == 1 else ""

    @property
    def key(self) -> str:
        slug = lambda s: re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")
        return f"join_{slug(self.a.name)}__{slug(self.b.name)}__{slug(self.pivot)}"

    @property
    def title(self) -> str:
        return f"{self.a.name} and {self.b.name}, on {self.pivot}"

    def latex(self) -> str:
        """One line per side, each labelled with the equation it came from."""
        tex = lambda s: s.replace("\\", "").replace("&", r"\&").replace("%", r"\%").replace("_", r"\_") \
                         .replace("#", r"\#").replace("$", r"\$").replace("{", r"\{").replace("}", r"\}")
        left = rf"\underbrace{{{sp.latex(self.ea)}}}_{{\text{{{tex(self.a.name)}}}}}"
        right = rf"\underbrace{{{sp.latex(self.eb)}}}_{{\text{{{tex(self.b.name)}}}}}"
        return f"{left}%%BR%%={right}"

    def verify(self) -> list[str]:
        """Machine checks on the join; raises AssertionError if one fails. Returns what passed."""
        passed = []
        for eq, expr in ((self.a, self.ea), (self.b, self.eb)):
            rel = eq.relation()
            sym = eq.symbol_for(self.pivot)
            residual = (rel.lhs - rel.rhs).subs(sym, expr)
            # every quantity is taken positive here, so sqrt(x**2) is x; the solve itself ran
            # without assumptions and refused anything with more than one root
            residual = sp.simplify(residual.xreplace(
                {s: sp.Symbol(s.name, positive=True) for s in residual.free_symbols}))
            assert residual == 0, f"{eq.name}: substituting the pivot leaves {residual}"
            passed.append(f"{eq.name} with {sym} replaced by its solution is an identity")
        assert sp.simplify(self.ea - self.eb) != 0, "the join is x = x"
        passed.append("the two sides are different expressions (the join is not x = x)")
        assert not self.ea.has(self.a.symbol_for(self.pivot)) and not self.eb.has(self.b.symbol_for(self.pivot))
        passed.append(f"the pivot ({self.pivot}) is eliminated from both sides")
        return passed


def shared_pivots(a, b) -> list[str]:
    """Quantities both equations denote, rarest first, without the too-common ones."""
    a, b = _as_eq(a), _as_eq(b)
    common = [q for q in a.quantities() if q in b.quantities() and q not in COMMON_QUANTITIES]
    return sorted(common, key=lambda q: (library().rarity.get(q, 0), q))


def join(a, b, pivot: str | None = None) -> Join | None:
    """Join two equations on a quantity both state. None if they cannot be joined."""
    a, b = _as_eq(a), _as_eq(b)
    for quantity in ([pivot] if pivot else shared_pivots(a, b)):
        if not agree(a, b, quantity) or clashes(a, b, quantity):
            continue
        left, right = solve_for(a, quantity), solve_for(b, quantity)
        if left is None or right is None:
            continue
        if sp.simplify(left.expr - right.expr) == 0:
            continue                                  # x = x
        notes = (f"Pivot: {quantity} (the letter {left.symbol} in the first equation, "
                 f"{right.symbol} in the second). Letters the two equations share and read alike "
                 f"are taken to be the same variable.")
        return Join(a, b, quantity, left.expr, right.expr,
                    rearranged=not (left.direct and right.direct), notes=notes)
    return None


@lru_cache(maxsize=1)
def all_joins() -> tuple[Join, ...]:
    """Every pairwise join the library yields (agreement, no clashes, solvable, not trivial)."""
    eqs = [e for e in library().equations.values() if e.relation() is not None]
    out = []
    for i, a in enumerate(eqs):
        for b in eqs[i + 1:]:
            j = join(a, b)
            if j is not None:
                out.append(j)
    return tuple(out)


def joins_for(what) -> list[Join]:
    eq = _as_eq(what)
    return [j for j in all_joins() if eq in (j.a, j.b)]


@dataclass
class Family:
    quantity: str
    members: list[Solution]
    equations: list[KBEquation]

    def latex(self) -> str:
        parts = [rf"\underbrace{{{sp.latex(s.expr)}}}_{{\text{{{e.name.replace('_', ' ')}}}}}"
                 for s, e in zip(self.members, self.equations)]
        return " = ".join(parts)


def families(minimum: int = 3) -> list[Family]:
    """Quantities that at least `minimum` equations state (directly or after rearranging), with
    distinct expressions."""
    lib = library()
    out = []
    for qname in sorted({q for e in lib.equations.values() for q in e.quantities()}):
        if qname in COMMON_QUANTITIES:
            continue
        members, eqs = [], []
        for e in lib.equations.values():
            sol = solve_for(e, qname)
            if sol is None:
                continue
            if any(sp.simplify(sol.expr - m.expr) == 0 for m in members):
                continue
            members.append(sol)
            eqs.append(e)
        if len(members) >= minimum:
            out.append(Family(qname, members, eqs))
    return sorted(out, key=lambda f: (-len(f.members), f.quantity))


# --------------------------------------------------------------------------------------
# dimensional analysis
# --------------------------------------------------------------------------------------
class DimensionError(Exception):
    pass


def _dim(expr, dims: dict):
    """Dimension of an expression as a monomial in the base symbols; raises DimensionError."""
    if expr.is_Number or expr in (sp.pi,):
        return sp.Integer(1)
    if expr.is_Symbol:
        if expr.name not in dims:
            raise DimensionError(f"no dimension for {expr.name}")
        return dims[expr.name]
    if expr.is_Add:
        ds = [sp.simplify(_dim(a, dims)) for a in expr.args]
        if any(sp.simplify(d / ds[0]) != 1 for d in ds):
            raise DimensionError("terms of a sum have different dimensions")
        return ds[0]
    if expr.is_Mul:
        out = sp.Integer(1)
        for a in expr.args:
            out *= _dim(a, dims)
        return out
    if expr.is_Pow:
        base, ex = expr.args
        db = _dim(base, dims)
        if sp.simplify(_dim(ex, dims)) != 1:
            raise DimensionError("exponent is not dimensionless")
        return db ** ex
    if isinstance(expr, (sp.exp, sp.log, sp.sin, sp.cos, sp.tan, sp.sinh, sp.cosh, sp.tanh)):
        if sp.simplify(_dim(expr.args[0], dims)) != 1:
            raise DimensionError(f"argument of {expr.func.__name__} is not dimensionless")
        return sp.Integer(1)
    raise DimensionError(f"cannot take the dimension of {expr.func.__name__}")


def check_dimensions(what) -> str | None:
    """None if the equation is dimensionally homogeneous; otherwise the reason it is not."""
    eq = _as_eq(what)
    rel = eq.relation()
    if rel is None:
        return "not readable"
    lib = library()
    dims = {}
    for letter, qname in eq.reads.items():
        if letter == "pi":
            continue
        q = lib.quantities.get(qname)
        if q is None:
            return f"quantity {qname!r} is not in the table"
        dims[letter] = _dimension_of(q.dimension)
    try:
        left, right = sp.simplify(_dim(rel.lhs, dims)), sp.simplify(_dim(rel.rhs, dims))
    except DimensionError as why:
        return str(why)
    if sp.simplify(left / right) != 1:
        return f"left side has dimension {left}, right side {right}"
    return None


def dimension_report() -> dict:
    """{'ok': [names], 'bad': [(name, reason)], 'unreadable': [(name, reason)]}."""
    ok, bad, skipped = [], [], []
    for eq in library().equations.values():
        if eq.relation() is None:
            skipped.append((eq.name, eq.refusal or ""))
            continue
        why = check_dimensions(eq)
        (ok.append(eq.name) if why is None else bad.append((eq.name, why)))
    return {"ok": ok, "bad": bad, "unreadable": skipped}


def refusals() -> dict[str, list[str]]:
    """Reason -> equation names the reader will not rearrange."""
    out: dict[str, list[str]] = {}
    for eq in library().equations.values():
        if eq.relation() is None:
            reason = (eq.refusal or "").split(":")[0]
            out.setdefault(reason, []).append(eq.name)
    return dict(sorted(out.items(), key=lambda kv: -len(kv[1])))


# --------------------------------------------------------------------------------------
# bridges to the eqlib registry
# --------------------------------------------------------------------------------------
def join_to_equation(j: Join, entry: str | None = None):
    """The join as an `eqlib.Equation` of kind 'derived', checks recorded on it."""
    import eqlib
    return eqlib.Equation(
        key=j.key, title=j.title, source_kind="derived",
        source=f"knowledgebase.join({j.a.name!r}, {j.b.name!r}); pivot {j.pivot}",
        expr=j.statement, latex_str=j.latex(),
        entry=j.entry if entry is None else entry, notes=j.notes
        + (" At least one side was rearranged from the quoted form." if j.rearranged else ""),
    )


# --------------------------------------------------------------------------------------
# selftest
# --------------------------------------------------------------------------------------
def selftest() -> list[str]:
    """Structural checks of the tables. Raises AssertionError at the first failure; returns
    the names of the checks that passed."""
    lib = library()
    passed = []

    def ok(name, condition, detail=""):
        assert condition, f"{name}: {detail}"
        passed.append(name)

    names = list(lib.equations)
    ok("equation names are unique", len(names) == len(set(names)) == len(EQUATIONS_DATA))
    unknown = sorted({q for r in READINGS_DATA.values() for q in r.values()} - set(lib.quantities) - {"pi"})
    ok("every reading names a known quantity", not unknown, str(unknown[:5]))
    ok("every reading belongs to an equation", set(READINGS_DATA) <= set(names),
       str(sorted(set(READINGS_DATA) - set(names))[:5]))
    ok("every link joins two equations in the library",
       all(a in lib.equations and b in lib.equations for a, _k, b, _n in lib.links))
    bad_dim = [q.name for q in lib.quantities.values() if _try_dim(q.dimension) is None]
    ok("every quantity has a parsable dimension", not bad_dim, str(bad_dim[:5]))
    ok("every join verifies", all(bool(j.verify()) for j in all_joins()))
    ok("no join is a collision", all(not clashes(j.a, j.b, j.pivot) for j in all_joins()))
    return passed


def _try_dim(spec: str):
    try:
        return _dimension_of(spec)
    except Exception:
        return None


# --------------------------------------------------------------------------------------
# the data, copied from RosettaMath's rosettaphys.py (MIT licence, Brent Hartshorn)
# --------------------------------------------------------------------------------------
QUANTITIES_DATA = [('energy',
  'M L^2 T^-2',
  'The capacity to do work, in joules. Kinetic, potential, thermal and rest energy are all this '
  'one quantity, which is why they can be added and why conversion between them is the subject of '
  'most of mechanics.',
  'Energy'),
 ('force', 'M L T^-2', 'A push, in newtons -- the rate at which momentum is delivered.', 'Force'),
 ('mass',
  'M',
  'Resistance to acceleration, and the source of gravity. That those are the same number is the '
  'equivalence principle, not a definition.',
  'Mass'),
 ('length', 'L', 'A distance or a size, in metres.', 'Length'),
 ('area', 'L^2', 'A surface, in square metres.', 'Area'),
 ('volume', 'L^3', 'A region of space, in cubic metres.', 'Volume'),
 ('time', 'T', 'A duration or an instant, in seconds.', 'Time'),
 ('speed',
  'L T^-1',
  'Distance covered per unit time. Velocity is the same quantity with a direction attached.',
  'Speed'),
 ('acceleration', 'L T^-2', 'The rate at which velocity changes.', 'Acceleration'),
 ('momentum',
  'M L T^-1',
  'Mass times velocity: the conserved quantity that makes collisions solvable.',
  'Momentum'),
 ('pressure', 'M L^-1 T^-2', 'Force per unit area, in pascals.', 'Pressure'),
 ('density',
  'M L^-3',
  'Mass per unit volume. One of at least four things rho is used for, and the one it is named '
  'after.',
  'Density'),
 ('temperature', 'K', 'A measure of the energy per degree of freedom, in kelvin.', 'Temperature'),
 ('entropy',
  'M L^2 T^-2 K^-1',
  'A count of the microstates consistent with what is known, in units that make it add to energy '
  'over temperature.',
  'Entropy'),
 ('microstate count',
  '1',
  "A pure number: how many arrangements a system could be in. It is what sits inside Boltzmann's "
  'logarithm, and it is emphatically not an energy, however much its letter looks like one.',
  'Microstate_(statistical_mechanics)'),
 ('frequency', 'T^-1', 'Cycles per second, in hertz.', 'Frequency'),
 ('angular frequency',
  'T^-1',
  'Radians per second. Dimensionally identical to frequency and numerically different by two pi, '
  'which is exactly the sort of distinction dimensional analysis cannot catch.',
  'Angular_frequency'),
 ('wavelength', 'L', 'The distance between repeats of a wave.', 'Wavelength'),
 ('charge', 'I T', 'Electric charge, in coulombs.', 'Electric_charge'),
 ('voltage',
  'M L^2 T^-3 I^-1',
  'Potential difference, in volts: energy per unit charge.',
  'Voltage'),
 ('current', 'I', 'Charge flow per unit time, in amperes.', 'Electric_current'),
 ('resistance', 'M L^2 T^-3 I^-2', 'Opposition to current, in ohms.', 'Electrical_resistance'),
 ('capacitance', 'M^-1 L^-2 T^4 I^2', 'Charge stored per volt, in farads.', 'Capacitance'),
 ('power', 'M L^2 T^-3', 'Energy per unit time, in watts.', 'Power_(physics)'),
 ('amount', 'N', 'A count of particles, in moles.', 'Amount_of_substance'),
 ('number density', 'L^-3', 'Particles per unit volume.', 'Number_density'),
 ('spring constant', 'M T^-2', 'Stiffness: force per unit displacement.', 'Hooke%27s_law'),
 ('gravitational field',
  'L T^-2',
  'The acceleration a gravitational source imposes, independent of what is falling. Dimensionally '
  'an acceleration, and kept apart from one because g and a appear in the same equations meaning '
  'different things.',
  'Gravitational_field'),
 ('flux', 'M T^-3', 'Energy crossing unit area per unit time.', 'Radiant_flux'),
 ('hubble parameter',
  'T^-1',
  'The expansion rate of the universe, an inverse time.',
  'Hubble%27s_law'),
 ('speed of light',
  'L T^-1',
  'Exactly 299792458 m/s, and a conversion factor between space and time rather than a speed '
  'anything travels at.',
  'Speed_of_light'),
 ('planck constant',
  'M L^2 T^-1',
  'The quantum of action. Note that h and hbar differ by two pi and are not interchangeable, which '
  'is the commonest way to be out by a factor of six in quantum mechanics.',
  'Planck_constant'),
 ('reduced planck constant',
  'M L^2 T^-1',
  'h over two pi, which is what appears whenever the angle is in radians -- so, nearly always.',
  'Planck_constant'),
 ('gravitational constant',
  'M^-1 L^3 T^-2',
  "Newton's G, the weakest and least precisely known of the constants.",
  'Gravitational_constant'),
 ('boltzmann constant',
  'M L^2 T^-2 K^-1',
  'The conversion between temperature and energy per degree of freedom.',
  'Boltzmann_constant'),
 ('gas constant',
  'M L^2 T^-2 K^-1 N^-1',
  'The Boltzmann constant per mole rather than per particle.',
  'Gas_constant'),
 ('permittivity',
  'M^-1 L^-3 T^4 I^2',
  'The permittivity of free space, which sets the strength of the electrostatic force.',
  'Vacuum_permittivity'),
 ('elementary charge',
  'I T',
  "The charge on a proton. Shares its letter with Euler's number, which is the single most "
  'dangerous collision in physics notation.',
  'Elementary_charge'),
 ('stefan-boltzmann constant',
  'M T^-3 K^-4',
  'The constant in the fourth-power radiation law.',
  'Stefan%E2%80%93Boltzmann_constant'),
 ('pi', '1', 'The ratio of a circumference to a diameter.', 'Pi'),
 ('dimensionless',
  '1',
  'A pure number: a ratio, a count, an index, or an argument to a logarithm.',
  'Dimensionless_quantity'),
 ('golden ratio',
  '1',
  'The positive root of x^2 = x + 1. It is the inflation factor of the Penrose tilings, and the '
  'reason those tilings cannot repeat: a periodic tiling forces a rational ratio of tile counts, '
  'and this number is not one.',
  'Golden_ratio'),
 ('inflation factor',
  '1',
  'What a substitution rule multiplies a length by each time it is applied. When it is irrational '
  'the tiling it generates has no translational period.',
  'Substitution_tiling'),
 ('patch area',
  'L^2',
  'The area of a finite patch of a tiling. An area, and joinable to nothing else that is one: a '
  'tiling patch can be made as large as you like, so equating it to a horizon area or a soap film '
  'would be a true conditional about nothing. Dimensional agreement is not the test -- the test is '
  'whether the two equations are about the same thing.',
  'Penrose_tiling'),
 ('tile count',
  '1',
  'How many tiles of one kind appear in a patch. Only ratios of these are meaningful, since the '
  'patch can always be made bigger.',
  'Penrose_tiling'),
 ('tile frequency',
  '1',
  'The limiting fraction of a tiling made of one species of tile. Unlike a tile count it '
  'converges, and for a substitution tiling it is a component of the Perron eigenvector of the '
  'substitution matrix -- which is why these come out as algebraic numbers rather than as '
  'measurements.',
  'Substitution_tiling'),
 ('tile area',
  'L^2',
  'The area of a single tile, as opposed to a patch of them. Kept apart from patch area for the '
  'reason given there: a tile has a fixed size and a patch does not, so equating one to the other '
  'would be a statement about nothing.',
  'Aperiodic_tiling'),
 ('linear inflation factor',
  '1',
  'What one round of substitution multiplies a length by. Its square is the inflation factor, '
  'which scales areas and tile counts. Kept apart from it because the two are used for different '
  'things and substituting one for the other gives plausible wrong answers.',
  'Substitution_tiling'),
 ('gap label',
  '1',
  'A number labelling a gap in the spectrum of an operator on a tiling, taken from a ring built '
  'from the inflation factor. Here only its arithmetic is recorded; what the labels count is a '
  'statement about the operator, not about the algebra.',
  'Aperiodic_tiling'),
 ('tile-count growth rate',
  '1',
  'The logarithm of the factor by which the number of tiles grows per round of substitution, in '
  'nats. A growth rate, not an entropy: primitive substitution tilings have zero topological '
  'entropy.',
  'Substitution_tiling'),
 ('genus',
  '1',
  'The number of handles on a closed orientable surface: nought for a sphere, one for a torus. In '
  'string perturbation theory it is also the loop order.',
  'Genus_(mathematics)'),
 ('euler characteristic',
  '1',
  'A whole number attached to a surface that no amount of bending can change. It is what makes '
  'topology quantitative, and it is the exponent the string coupling is raised to.',
  'Euler_characteristic'),
 ('covering degree',
  '1',
  'How many sheets of a covering surface sit over each point of the surface below.',
  'Covering_space'),
 ('mean curvature',
  'L^-1',
  'The average of the two principal curvatures. Setting it to zero is what makes a surface '
  'minimal.',
  'Mean_curvature'),
 ('principal curvature',
  'L^-1',
  'The largest and smallest curvatures of the normal sections at a point on a surface.',
  'Principal_curvature'),
 ('gaussian curvature',
  'L^-2',
  'The product of the principal curvatures. Unlike the mean curvature it survives bending the '
  'surface without stretching it, which is what Gauss called remarkable.',
  'Gaussian_curvature'),
 ('geodesic curvature',
  'L^-1',
  'How far a curve drawn on a surface bends away from being a geodesic of that surface.',
  'Geodesic_curvature'),
 ('slope',
  '1',
  'A partial derivative of a height function: how steeply a surface written as a graph rises in '
  'one direction.',
  'Slope'),
 ('hausdorff dimension',
  '1',
  'A dimension that need not be a whole number, defined by how a covering of the set scales. A '
  'Kakeya set can have measure zero and still have full Hausdorff dimension.',
  'Hausdorff_dimension'),
 ('count',
  '1',
  'A plain whole number of things -- steps in a construction, sheets in a cover.',
  'Natural_number'),
 ('action',
  'M L^2 T^-1',
  'Energy times time: the quantity a physical path makes stationary. Dimensionally identical to '
  'the Planck constant, which is the point of calling h the quantum of action.',
  'Action_(physics)'),
 ('regge slope',
  'L^2',
  'Alpha prime, the one free parameter of the string: the square of the string length, and the '
  'slope of the line that mass-squared and spin lie on.',
  'Regge_theory'),
 ('string tension',
  'L^-2',
  'Energy per unit length of a string, which in natural units is the reciprocal of the Regge slope '
  'up to a factor of two pi.',
  'String_(physics)'),
 ('string coupling',
  '1',
  'The number that controls how readily a worldsheet splits. It is not put in by hand: it is set '
  'by the vacuum value of the dilaton.',
  'String_theory'),
 ('spin',
  '1',
  'Angular momentum in units of hbar, which for the states on a Regge trajectory is a whole '
  'number.',
  'Spin_(physics)'),
 ('partition function',
  '1',
  'A sum over configurations, weighted. Everything thermodynamic, and in string theory everything '
  'perturbative, is a derivative of it.',
  'Partition_function_(statistical_mechanics)'),
 ('spacetime dimension',
  '1',
  'How many coordinates it takes to label an event. A count, but kept apart from other counts on '
  'purpose: pivoting the dimension of spacetime against a number of steps in a construction would '
  'be dimensionally impeccable and entirely meaningless.',
  'Dimension_(vector_space)'),
 ('inverse coupling',
  '1',
  'The lattice parameter beta, which is 2N over the gauge coupling squared: large for weak '
  'coupling. A bare input of a simulation, not a measured quantity.',
  'Lattice_gauge_theory'),
 ('gauge coupling',
  '1',
  'The strength with which a gauge field couples to itself and to matter. Not constant: it runs '
  'with the energy scale, which is why the lattice coupling has to be tuned to reach a continuum '
  'limit.',
  'Coupling_constant'),
 ('color count',
  '1',
  'The N of SU(N): two for the weak-isospin-like toy case the laboratory simulates, three for the '
  'strong interaction.',
  'Special_unitary_group'),
 ('lattice action',
  '1',
  "The Wilson action, a pure number in units of Planck's constant. Kept apart from the continuum "
  'action because one is a sum over squares and the other an integral over spacetime, and the '
  'statement that they agree is a limit.',
  'Wilson_loop'),
 ('plaquette variable',
  '1',
  'The ordered product of the four link variables around an elementary square: the smallest closed '
  'loop on the lattice.',
  'Lattice_gauge_theory'),
 ('link variable',
  '1',
  'The group element that carries a charge along one edge of the lattice, standing in for the '
  'exponential of the gauge field integrated along it.',
  'Lattice_gauge_theory'),
 ('gauge transformation',
  '1',
  'A group element attached to a site, rotating the internal orientation there. Nothing observable '
  'depends on it.',
  'Gauge_theory'),
 ('gauge field strength',
  'L^-2',
  'The curvature of the gauge field, in natural units. For a non-Abelian group it contains a '
  'commutator of the potentials as well as their derivatives.',
  'Field_strength'),
 ('wilson loop',
  '1',
  'The trace of the product of link variables around a closed rectangle, averaged. How it decays '
  'with the size of the rectangle is the lattice test for confinement.',
  'Wilson_loop'),
 ('loop area',
  'L^2',
  'The area enclosed by a Wilson loop, in lattice units on a lattice. Kept apart from area for the '
  'reason patch area and tile area are: a horizon and a Wilson loop are both areas, and equating '
  'them is dimensionally impeccable and meaningless.',
  'Wilson_loop'),
 ('creutz ratio',
  '1',
  'A ratio of four Wilson loops in which perimeter and corner terms cancel, leaving an estimate of '
  'the string tension in lattice units.',
  'Wilson_loop'),
 ('confining string tension',
  'L^-2',
  'The constant in the area law of a Wilson loop: the energy per unit length of the flux tube '
  'between static charges. Not the same quantity as the string tension of string theory, however '
  'similar the name and the dimension.',
  'Confinement'),
 ('lattice spacing',
  'L',
  'The distance between neighbouring sites. Not an input of the simulation but fixed afterwards, '
  'by matching one measured quantity to its physical value.',
  'Lattice_gauge_theory'),
 ('lattice mass',
  '1',
  'A mass in units of the inverse lattice spacing, which is what a correlator on a lattice '
  "measures. Multiplying by the spacing's inverse is what makes it physical.",
  'Lattice_QCD'),
 ('time slice',
  '1',
  'A Euclidean time coordinate counted in lattice spacings: an index, not a duration. Kept apart '
  'from time because a lattice mass times a time slice is dimensionless and a mass times a time is '
  'only so in natural units.',
  'Lattice_gauge_theory'),
 ('correlator',
  '1',
  'The expectation of an operator at one Euclidean time against another, connected so that the '
  'vacuum drops out.',
  'Correlation_function_(quantum_field_theory)'),
 ('amplitude',
  '1',
  'The overlap of an operator with the state that dominates its correlator, squared.',
  'Amplitude'),
 ('autocorrelation',
  '1',
  'The normalised covariance of a Markov chain with a copy of itself shifted by some number of '
  'updates: one at no shift, decaying towards zero.',
  'Autocorrelation'),
 ('autocorrelation time',
  '1',
  'The number of updates it takes the chain to forget where it was, measured in updates rather '
  'than seconds. That distinction is the whole of the sampling-efficiency result.',
  'Autocorrelation'),
 ('sample count',
  '1',
  'How many measurements were taken, correlated or not.',
  'Sample_size_determination'),
 ('effective sample size',
  '1',
  'How many independent measurements a correlated chain is worth.',
  'Effective_sample_size'),
 ('sampling efficiency',
  'T^-1',
  'Effective samples per second of wall-clock time.',
  'Effective_sample_size'),
 ('efficiency ratio',
  '1',
  "One schedule's sampling efficiency over another's, paired from a common starting configuration.",
  'Overrelaxation'),
 ('spin action',
  '1',
  'The energy of a spin configuration in units of the temperature: minus the sum of dot products '
  'of neighbouring unit vectors, times beta.',
  'N-vector_model'),
 ('spin inverse temperature',
  '1',
  'The beta of a spin model: the coupling between neighbouring spins in units of the temperature. '
  'Large means nearly aligned. Not the inverse coupling of a gauge theory, whatever the letter.',
  'N-vector_model'),
 ('spin variable',
  '1',
  'A unit vector with N components sitting on one site of the lattice.',
  'N-vector_model'),
 ('component count',
  '1',
  'The N of the O(N) model: how many components each spin has. Not the number of colours of a '
  'gauge group.',
  'N-vector_model'),
 ('scaling constant',
  '1',
  'The non-universal constant in front of an asymptotic scaling form. The form predicts the '
  'exponent and not this.',
  'Asymptotic_freedom'),
 ('correlation length',
  '1',
  'The distance over which spins stay correlated, in lattice units. The second-moment estimate is '
  'one way to measure it and is not the transfer-matrix gap.',
  'Correlation_length'),
 ('spin susceptibility',
  '1',
  'The sum of the spin-spin correlation over the whole lattice, per site.',
  'Magnetic_susceptibility'),
 ('structure factor',
  '1',
  'The spin correlation at the smallest non-zero momentum on a periodic box.',
  'Structure_factor'),
 ('lattice momentum',
  '1',
  'The smallest non-zero momentum on a periodic box of L sites, 2 sin(pi / L) in inverse lattice '
  'spacings.',
  'Lattice_(group)'),
 ('stiffness deviation',
  '1',
  'How far the stiffness of the XY model is from its universal-jump value, in the variable of the '
  'Kosterlitz flow.',
  'Kosterlitz%E2%80%93Thouless_transition'),
 ('vortex fugacity',
  '1',
  'The weight of a vortex, in the normalisation of the Kosterlitz flow.',
  'Kosterlitz%E2%80%93Thouless_transition'),
 ('rg scale',
  '1',
  'The logarithm of the length scale, the time of a renormalisation-group flow.',
  'Renormalization_group'),
 ('magnon energy', '1', 'The energy of one spin flip, in units of the exchange.', 'Magnon'),
 ('lattice wave number',
  '1',
  'A component of the wave vector, in inverse lattice spacings.',
  'Wave_vector'),
 ('magnon density',
  '1',
  'The number of thermally excited magnons per site in the ideal gas of non-interacting ones.',
  'Magnon'),
 ('scaled inverse temperature',
  '1',
  'The product of inverse temperature and spin, the variable the magnon density is expanded in.',
  'Bloch%27s_law'),
 ('magnon coefficient',
  '1',
  "A coefficient in the large-x expansion of the magnon density; the first is Bloch's.",
  'Bloch%27s_law'),
 ('mixing angle',
  '1',
  'The angle that mixes the bilinear and biquadratic couplings of the spin-one chain.',
  'Heisenberg_model_(quantum)'),
 ('spin operator', '1', 'The spin-one operator on a site, taken as a vector.', 'Spin_(physics)'),
 ('spin hamiltonian',
  'M L^2 T^-2',
  'The energy operator of a spin chain, in units of the exchange.',
  'Hamiltonian_(quantum_mechanics)'),
 ('spin gap',
  'M L^2 T^-2',
  'The energy of the lowest state with one more unit of total spin than the ground state.',
  'Spin_gap'),
 ('spin energy level',
  'M L^2 T^-2',
  'An energy eigenvalue of a spin chain in a given sector of total spin.',
  'Energy_level'),
 ('fine-structure coupling',
  '1',
  'A gauge coupling in the form alpha = g^2 / 4 pi, a function of the energy scale.',
  'Fine-structure_constant'),
 ('beta function coefficient',
  '1',
  'The one-loop slope of the inverse coupling against the logarithm of the scale. Negative is '
  'asymptotic freedom.',
  'Beta_function_(physics)'),
 ('unification scale',
  'M L^2 T^-2',
  'The energy at which the couplings are taken to meet. A modulus in a compactification, not a '
  'computed number.',
  'Grand_Unified_Theory'),
 ('renormalisation scale',
  'M L^2 T^-2',
  'The energy at which a running coupling is evaluated.',
  'Renormalization_group'),
 ('qcd scale',
  'M L^2 T^-2',
  'The energy at which the one-loop strong coupling diverges.',
  'QCD_vacuum'),
 ('flavour count',
  '1',
  'The number of quark flavours, or of chiral flavour pairs in a supersymmetric spectrum.',
  'Flavour_(particle_physics)'),
 ('vacuum polarisation shift',
  '1',
  'The fractional shift between the electromagnetic coupling at the Z mass and at zero momentum. '
  'Its hadronic part is measured.',
  'Fine-structure_constant'),
 ('chiral index',
  '1',
  'The index of the bundle: minus the net number of chiral generations upstairs. An integer from '
  'topology.',
  'Atiyah%E2%80%93Singer_index_theorem'),
 ('generation count',
  '1',
  'The net number of chiral generations.',
  'Generation_(particle_physics)'),
 ('symmetry order',
  '1',
  'The order of the freely acting group of the quotient.',
  'Order_(group_theory)'),
 ('fivebrane class',
  '1',
  'The class of the five-branes needed to cancel the anomaly; it has to be effective.',
  'Green%E2%80%93Schwarz_mechanism'),
 ('chern class',
  '1',
  'A characteristic class of a bundle, here integrated against a Kahler class to give a number.',
  'Chern_class'),
 ('lie algebra dimension', '1', 'The dimension of a Lie algebra as a vector space.', 'Lie_algebra'),
 ('root count', '1', 'The number of roots of a semisimple Lie algebra.', 'Root_system'),
 ('rank', '1', 'The dimension of a Cartan subalgebra.', 'Rank_of_a_Lie_algebra'),
 ('subalgebra dimension',
  '1',
  'The dimension of a subalgebra of a larger Lie algebra.',
  'Lie_algebra'),
 ('representation dimension',
  '1',
  'The complex dimension of a representation.',
  'Group_representation'),
 ('weyl group order', '1', 'The order of the Weyl group of a root system.', 'Weyl_group'),
 ('stabiliser order',
  '1',
  'The order of the stabiliser of a point under a group action.',
  'Group_action'),
 ('orbit size', '1', 'The number of points in one orbit of a group action.', 'Group_action'),
 ('triple count',
  '1',
  'The number of ordered triples drawn from three finite sets.',
  'Cartesian_product'),
 ('triality algebra dimension',
  '1',
  'The dimension of the triality algebra of a division algebra.',
  'Triality'),
 ('division algebra dimension',
  '1',
  'The dimension of R, C, H or O over the reals.',
  'Division_algebra'),
 ('cubic anomaly coefficient',
  '1',
  'The symmetrised cubic trace of generators in a representation: the group factor of the triangle '
  'anomaly.',
  'Gauge_anomaly'),
 ('integer charge',
  '1',
  'The charge of a representation of the Standard Model group in the normalisation where the quark '
  'doublet has charge one.',
  'Standard_Model'),
 ('colour triality', '1', 'The class in Z/3 of how the centre of SU(3) acts.', 'Triality'),
 ('weak centre parity', '1', 'The class in Z/2 of how the centre of SU(2) acts.', 'Weak_isospin'),
 ('hypercharge',
  '1',
  'The U(1) charge of the Standard Model, in the convention Q = T3 + Y.',
  'Weak_hypercharge'),
 ('weak isospin', '1', 'The third component of weak isospin.', 'Weak_isospin'),
 ('hypercharge normalisation',
  '1',
  'The ratio of summed squared hypercharge to summed squared weak isospin over a complete '
  'multiplet.',
  'Weinberg_angle'),
 ('net chiral class',
  '1',
  'The multiplicity of each representation minus that of its dual, as an element of a quotient of '
  'the representation ring.',
  'Chirality_(physics)'),
 ('dual coxeter number',
  '1',
  'The integer h in the dual Coxeter number of a simple Lie algebra: N for SU(N), twelve for E6, '
  'thirty for E8.',
  'Dual_Coxeter_number'),
 ('adjoint index',
  '1',
  'The Dynkin index of the adjoint representation, which is the dual Coxeter number in the usual '
  'normalisation.',
  'Dynkin_index'),
 ('representation index',
  '1',
  'The Dynkin index of the representation the fermions are in: one half for the fundamental of '
  'SU(N).',
  'Dynkin_index'),
 ('spin transmutation slope',
  '1',
  'How fast the log of the mass falls with the inverse temperature of an O(N) model. Kept apart '
  'from the gauge one: they are the same kind of number for different systems, and equating them '
  'would equate their N.',
  'Asymptotic_freedom'),
 ('gauge transmutation slope',
  '1',
  'How fast the log of the mass falls with the inverse coupling of a pure SU(N) lattice gauge '
  'theory, in the convention beta = 2N / g^2.',
  'Asymptotic_freedom'),
 ('blocked inverse coupling',
  '1',
  'The inverse coupling of the lattice after one blocking step, which is a different lattice and '
  'so a different quantity from the one before it.',
  'Renormalization_group'),
 ('blocked spin inverse temperature',
  '1',
  'The inverse temperature of the O(N) model after one blocking step.',
  'Renormalization_group')]

EQUATIONS_DATA = [('Mass-energy equivalence',
  'Relativity',
  'E = m c^2',
  None,
  'Mass%E2%80%93energy_equivalence',
  'Mass and energy are the same thing in different units, and the conversion factor is enormous. '
  'The equation does not say mass turns into energy; it says a system with energy has inertia, and '
  'a body at rest already carries energy m*c^2.'),
 ("Newton's second law",
  'Classical mechanics',
  'F = m a',
  None,
  'Newton%27s_laws_of_motion',
  'Force equals mass times acceleration -- more precisely, force is the rate of change of '
  'momentum. It is the definition that makes mass measurable and turns mechanics into a solvable '
  'differential equation.'),
 ("Newton's law of gravitation",
  'Classical mechanics',
  'F = G \\frac{m_1 m_2}{r^2}',
  None,
  'Newton%27s_law_of_universal_gravitation',
  'An inverse-square attraction between any two masses. The inverse square is not arbitrary: it is '
  'what a flux spreading over the surface of a sphere must do in three dimensions.'),
 ("Coulomb's law",
  'Electromagnetism',
  'F = \\frac{1}{4 \\pi \\epsilon_0} \\frac{q_1 q_2}{r^2}',
  None,
  'Coulomb%27s_law',
  'The electrostatic force between two charges. Structurally identical to Newtonian gravity, but '
  'roughly 10^36 times stronger and able to take either sign -- which is why bulk matter is '
  'electrically neutral and gravity wins at large scales.'),
 ('Schrodinger equation (time-dependent)',
  'Quantum mechanics',
  'i \\hbar \\frac{\\partial \\Psi}{\\partial t} = \\hat{H} \\Psi',
  None,
  'Schr%C3%B6dinger_equation',
  'The equation of motion for a quantum state. The i on the left makes it a wave equation rather '
  'than a diffusion equation: instead of smoothing out, solutions rotate in phase and interfere. '
  'It is first order in time, so the present state determines the entire future.'),
 ('Schrodinger equation (time-independent)',
  'Quantum mechanics',
  '-\\frac{\\hbar^2}{2m} \\nabla^2 \\psi + V \\psi = E \\psi',
  None,
  'Schr%C3%B6dinger_equation',
  'An eigenvalue problem: only certain energies admit well-behaved solutions, and that '
  'discreteness is where quantisation comes from. Atomic energy levels are the eigenvalues of this '
  'equation for a Coulomb potential.'),
 ('Klein-Gordon equation',
  'Quantum field theory',
  '\\left( \\Box + \\frac{m^2 c^2}{\\hbar^2} \\right) \\phi = 0',
  None,
  'Klein%E2%80%93Gordon_equation',
  'The relativistic wave equation for a spinless field. It is second order in time, unlike the '
  'Schrodinger equation, which is what makes it Lorentz invariant -- and also what forced the '
  'reinterpretation of its negative-energy solutions as antiparticles.'),
 ('Dirac equation',
  'Quantum field theory',
  '( i \\gamma^\\mu \\partial_\\mu - m ) \\psi = 0',
  None,
  'Dirac_equation',
  'The relativistic equation for spin-1/2 particles. Dirac insisted on first order in time, which '
  'required the coefficients to be matrices; spin and antimatter both fell out of the algebra '
  'rather than being put in by hand.'),
 ("d'Alembert operator",
  'Relativity',
  '\\Box = \\frac{1}{c^2} \\frac{\\partial^2}{\\partial t^2} - \\nabla^2',
  None,
  "D'Alembert_operator",
  'The definition of the box operator: a second time derivative against a spatial Laplacian, with '
  'the relative minus sign that is the whole difference between Minkowski spacetime and Euclidean '
  'space. Written out this way it is clear why the operator is Lorentz invariant and the Laplacian '
  'alone is not.'),
 ('Wave equation (box form)',
  'Waves',
  '\\Box \\psi = 0',
  None,
  'Wave_equation',
  "The wave equation written with the d'Alembertian. Compressing the time and space derivatives "
  'into one symbol makes the relativistic content visible at a glance: the equation has the same '
  'form in every inertial frame, and its solutions propagate at exactly the speed c buried inside '
  'the operator.'),
 ('Wave equation',
  'Waves',
  '\\frac{\\partial^2 u}{\\partial t^2} = c^2 \\nabla^2 u',
  None,
  'Wave_equation',
  'A second time derivative equal to a Laplacian: disturbances propagate at fixed speed c without '
  'changing shape. Sound, light, and a plucked string all obey it, and the constant c is set by '
  'the medium.'),
 ('Heat / diffusion equation',
  'Thermodynamics',
  '\\frac{\\partial u}{\\partial t} = \\alpha \\nabla^2 u',
  None,
  'Heat_equation',
  'A first time derivative equal to a Laplacian: gradients smooth out irreversibly. Unlike the '
  'wave equation it has a preferred direction of time -- you cannot run it backwards stably, which '
  'is the arrow of time appearing in a differential equation.'),
 ("Laplace's equation",
  'Fields',
  '\\nabla^2 \\phi = 0',
  None,
  'Laplace%27s_equation',
  'The Laplacian vanishes: the field everywhere equals the average of its neighbours. Such '
  'harmonic functions describe static fields in empty space, and they have no local maxima or '
  'minima in the interior -- extremes live on the boundary.'),
 ("Poisson's equation",
  'Fields',
  '\\nabla^2 \\phi = \\rho',
  None,
  'Poisson%27s_equation',
  "Laplace's equation with a source. The density on the right tells the potential how to curve; "
  'solving it recovers the gravitational or electrostatic potential produced by a given '
  'distribution of mass or charge.'),
 ('Helmholtz equation',
  'Waves',
  '\\nabla^2 A + k^2 A = 0',
  None,
  'Helmholtz_equation',
  'What the wave equation becomes after separating out a single frequency. k is the wavenumber; '
  'the equation governs standing waves, waveguides and optical modes.'),
 ("Gauss's law",
  'Electromagnetism',
  '\\nabla \\cdot E = \\frac{\\rho}{\\epsilon_0}',
  None,
  'Gauss%27s_law',
  'The divergence of the electric field is the charge density: field lines begin and end on '
  'charge. Integrated over a closed surface it says the flux out equals the charge enclosed, '
  'regardless of how that charge is arranged.'),
 ("Gauss's law for magnetism",
  'Electromagnetism',
  '\\nabla \\cdot B = 0',
  None,
  'Gauss%27s_law_for_magnetism',
  'The magnetic field has zero divergence: there are no magnetic monopoles. Every field line '
  'closes on itself, which is why cutting a magnet in half gives two magnets rather than a north '
  'and a south pole.'),
 ("Faraday's law of induction",
  'Electromagnetism',
  '\\nabla \\times E = -\\frac{\\partial B}{\\partial t}',
  None,
  'Faraday%27s_law_of_induction',
  'A changing magnetic field creates a circulating electric field. This is the principle behind '
  "every generator and transformer, and the minus sign (Lenz's law) is what keeps the induced "
  'effect opposing the change that caused it.'),
 ('Ampere-Maxwell law',
  'Electromagnetism',
  '\\nabla \\times B = \\mu_0 J + \\mu_0 \\epsilon_0 \\frac{\\partial E}{\\partial t}',
  None,
  'Amp%C3%A8re%27s_circuital_law',
  "Currents and changing electric fields both create circulating magnetic fields. Maxwell's "
  'addition of the second term is what made the equations self-consistent and predicted '
  'electromagnetic waves travelling at the speed of light.'),
 ('Continuity equation',
  'Conservation laws',
  '\\frac{\\partial \\rho}{\\partial t} + \\nabla \\cdot J = 0',
  None,
  'Continuity_equation',
  'The local statement of a conservation law: whatever the density loses, the flux must carry '
  'across the boundary. The same form governs mass, charge, energy and probability.'),
 ('Navier-Stokes equation',
  'Fluid dynamics',
  '\\rho \\left( \\frac{\\partial v}{\\partial t} + v \\cdot \\nabla v \\right) = -\\nabla p + '
  '\\eta \\nabla^2 v',
  None,
  'Navier%E2%80%93Stokes_equations',
  "Newton's second law for a fluid. The nonlinear term in which the velocity advects itself is "
  'what produces turbulence, and why proving that smooth solutions always exist is an open '
  'Millennium Prize problem.'),
 ('Euler-Lagrange equation',
  'Classical mechanics',
  '\\frac{d}{dt} \\frac{\\partial L}{\\partial \\dot{q}} - \\frac{\\partial L}{\\partial q} = 0',
  None,
  'Euler%E2%80%93Lagrange_equation',
  'The condition for the action to be stationary. Feed it a Lagrangian and it hands back the '
  'equations of motion; for a particle in a potential it reduces to F = ma, but it works just as '
  'well in awkward coordinates where forces are painful.'),
 ('Boltzmann entropy',
  'Thermodynamics',
  'S = k_B \\log W',
  None,
  'Boltzmann%27s_entropy_formula',
  "Entropy is the logarithm of the number of microstates. This formula, carved on Boltzmann's "
  'gravestone, is the bridge between the microscopic world of atoms and the macroscopic laws of '
  'thermodynamics.'),
 ('Shannon entropy',
  'Information theory',
  'H = -\\sum p_i \\log p_i',
  None,
  'Entropy_(information_theory)',
  'The average information content of a distribution, measured in bits when the logarithm is base '
  "two. It is Boltzmann's formula again in a different guise, and it sets the hard limit on how "
  'far data can be compressed.'),
 ('Planck relation',
  'Quantum mechanics',
  'E = h \\nu',
  None,
  'Planck_relation',
  "A photon's energy is proportional to its frequency, with Planck's constant as the conversion. "
  'This is the quantum hypothesis in its smallest form: light comes in discrete packets, which is '
  'why the photoelectric effect depends on colour rather than brightness.'),
 ('de Broglie relation',
  'Quantum mechanics',
  '\\lambda = \\frac{h}{p}',
  None,
  'Matter_wave',
  'Every particle has a wavelength inversely proportional to its momentum. For everyday objects it '
  'is unmeasurably small; for electrons it is atom-sized, which is why they diffract and why '
  'electron microscopes work.'),
 ('Heisenberg uncertainty principle',
  'Quantum mechanics',
  '\\Delta x \\Delta p \\geq \\frac{\\hbar}{2}',
  None,
  'Uncertainty_principle',
  'Position and momentum cannot both be sharply defined. This is not a limit on instruments but a '
  'property of Fourier conjugate pairs: a narrow wave packet in space is necessarily broad in '
  'wavenumber.'),
 ('Einstein field equations',
  'General relativity',
  'G_{\\mu\\nu} + \\Lambda g_{\\mu\\nu} = \\frac{8 \\pi G}{c^4} T_{\\mu\\nu}',
  None,
  'Einstein_field_equations',
  'Curvature on the left, matter and energy on the right: spacetime tells matter how to move, '
  'matter tells spacetime how to curve. Ten coupled nonlinear equations, which is why exact '
  'solutions are rare and precious.'),
 ('Schwarzschild radius',
  'General relativity',
  'r_s = \\frac{2 G M}{c^2}',
  None,
  'Schwarzschild_radius',
  'The radius at which the escape velocity reaches the speed of light -- the event horizon of a '
  'non-rotating black hole. For the Sun it is about three kilometres, for the Earth about nine '
  'millimetres.'),
 ('Bekenstein-Hawking entropy',
  'General relativity',
  'S = \\frac{k_B c^3 A}{4 G \\hbar}',
  None,
  'Black_hole_thermodynamics',
  "A black hole's entropy is proportional to the area of its horizon, not its volume. That every "
  'fundamental constant appears at once -- quantum, relativistic, gravitational and thermodynamic '
  '-- is why this formula is treated as a clue to quantum gravity, and it is the origin of the '
  'holographic principle.'),
 ('Lorentz factor',
  'Relativity',
  '\\gamma = \\frac{1}{\\sqrt{1 - \\frac{v^2}{c^2}}}',
  None,
  'Lorentz_factor',
  'The factor by which time dilates and length contracts. It is essentially one at everyday speeds '
  'and diverges as v approaches c, which is why massive objects cannot reach the speed of light.'),
 ('Ideal gas law',
  'Thermodynamics',
  'P V = n R T',
  None,
  'Ideal_gas_law',
  'Pressure times volume equals amount times temperature. An excellent approximation whenever the '
  'molecules are far enough apart to ignore their size and mutual attraction.'),
 ('Stefan-Boltzmann law',
  'Thermodynamics',
  'j = \\sigma T^4',
  None,
  'Stefan%E2%80%93Boltzmann_law',
  'Radiated power scales as the fourth power of temperature. The steepness is why a modest rise in '
  "a star's surface temperature makes it dramatically brighter."),
 ('Plasma frequency',
  'Plasma physics',
  '\\omega_p = \\sqrt{\\frac{n e^2}{\\epsilon_0 m_e}}',
  None,
  'Plasma_oscillation',
  'The natural oscillation frequency of electrons in a plasma. Waves below it cannot propagate and '
  'are reflected -- which is how the ionosphere bounces radio signals around the curve of the '
  'Earth.'),
 ('Debye length',
  'Plasma physics',
  '\\lambda_D = \\sqrt{\\frac{\\epsilon_0 k_B T}{n e^2}}',
  None,
  'Debye_length',
  'The distance over which a charge is screened out by the surrounding plasma. Beyond it the '
  'plasma looks neutral, and a system is only genuinely a plasma if it is much larger than this '
  'length.'),
 ('Vlasov equation',
  'Plasma physics',
  '\\frac{\\partial f}{\\partial t} + v \\cdot \\nabla f + \\frac{F}{m} \\cdot \\nabla_v f = 0',
  None,
  'Vlasov_equation',
  'Tracks the distribution of particles in position-velocity space for a collisionless plasma. It '
  'is a continuity equation in six dimensions, with the electromagnetic force supplied '
  'self-consistently by the particles themselves.'),
 ('Friedmann equation',
  'Cosmology',
  'H^2 = \\frac{8 \\pi G}{3} \\rho - \\frac{k c^2}{a^2} + \\frac{\\Lambda c^2}{3}',
  None,
  'Friedmann_equations',
  'Governs the expansion of a homogeneous universe. The competing terms are matter, spatial '
  'curvature and the cosmological constant, and which one dominates decides whether the universe '
  'expands forever.'),
 ("Euler's identity",
  'Mathematics',
  'e^{i \\pi} + 1 = 0',
  None,
  'Euler%27s_identity',
  'Rotating by half a turn in the complex plane lands you at -1. It links the additive identity, '
  'the multiplicative identity, and the three constants e, i and pi in one line.'),
 ('Pythagorean theorem',
  'Mathematics',
  'a^2 + b^2 = c^2',
  None,
  'Pythagorean_theorem',
  'In a right triangle the squares on the legs sum to the square on the hypotenuse. Generalised, '
  'it is the definition of distance in Euclidean space -- and changing the signs gives the '
  'spacetime interval of relativity.'),
 ('Quadratic formula',
  'Mathematics',
  'x = \\frac{-b \\pm \\sqrt{b^2 - 4ac}}{2a}',
  None,
  'Quadratic_formula',
  'The roots of a quadratic. The discriminant under the root decides everything: positive gives '
  'two real roots, zero a repeated one, negative a complex conjugate pair.'),
 ("Bayes' theorem",
  'Probability',
  'P(A \\mid B) = \\frac{P(B \\mid A) P(A)}{P(B)}',
  None,
  'Bayes%27_theorem',
  'How to update a belief when evidence arrives. The prior is multiplied by how well the '
  'hypothesis predicted the data and renormalised -- the whole of Bayesian inference is repeated '
  'application of this one line.'),
 ('Normal distribution',
  'Probability',
  'f(x) = \\frac{1}{\\sigma \\sqrt{2\\pi}} e^{-\\frac{(x-\\mu)^2}{2\\sigma^2}}',
  None,
  'Normal_distribution',
  'The bell curve. The central limit theorem explains its ubiquity: add up enough independent '
  'small effects and the total is normally distributed almost regardless of what you started '
  'with.'),
 ('Kullback-Leibler divergence',
  'Information theory',
  'D(P \\parallel Q) = \\sum P(x) \\log \\frac{P(x)}{Q(x)}',
  None,
  'Kullback%E2%80%93Leibler_divergence',
  'The relative entropy between two distributions: the extra information cost of using the wrong '
  'model. It is never negative, zero only when the distributions agree, and is not symmetric -- so '
  'it is a divergence, not a distance.'),
 ('Geometric series',
  'Mathematics',
  '\\sum_{n=0}^{\\infty} r^n = \\frac{1}{1-r}',
  None,
  'Geometric_series',
  'A sum with a constant ratio between terms converges whenever that ratio is less than one in '
  'magnitude. It is the workhorse behind perturbation expansions and the resolvent of an '
  'operator.'),
 ('Kinetic energy',
  'Classical mechanics',
  'E_k = \\frac{1}{2} m v^2',
  None,
  'Kinetic_energy',
  'The energy a body has by moving. The half and the square are not decoration: the square is what '
  'makes kinetic energy add as work does, and the half is what falls out of integrating momentum '
  'with respect to velocity.'),
 ('Momentum',
  'Classical mechanics',
  'p = m v',
  None,
  'Momentum',
  'Mass times velocity -- the quantity that is conserved when nothing external pushes, which is '
  'why it, rather than velocity, is what appears in the laws.'),
 ('Work done by a force',
  'Classical mechanics',
  'W = F d',
  None,
  'Work_(physics)',
  'Force times the distance moved along it. This is the bridge between the force laws and the '
  'energy laws: it is how a push becomes a change in energy.'),
 ('Gravitational potential energy',
  'Classical mechanics',
  'U = m g h',
  None,
  'Gravitational_energy',
  'The energy stored by lifting something in a uniform field. The uniform g is an approximation to '
  'the inverse-square law that is excellent near a surface and useless far from one.'),
 ('Surface gravity',
  'Classical mechanics',
  'g = \\frac{G M}{r^2}',
  None,
  'Surface_gravity',
  'What the inverse-square law becomes when one of the two masses is a planet: the acceleration it '
  'imposes, independent of what is falling. That independence is the equivalence principle in its '
  'oldest form.'),
 ("Hooke's law",
  'Classical mechanics',
  'F = -k x',
  None,
  'Hooke%27s_law',
  'A restoring force proportional to displacement. Almost nothing obeys it exactly, and almost '
  'everything obeys it near equilibrium, because it is the first term of any potential expanded '
  'about a minimum.'),
 ('Angular frequency',
  'Waves',
  '\\omega = 2 \\pi f',
  None,
  'Angular_frequency',
  'Frequency in radians per second rather than cycles per second. The factor of two pi is the '
  'single commonest source of a lost constant when moving between a formula and its '
  'implementation.'),
 ('Wave speed',
  'Waves',
  'v = f \\lambda',
  None,
  'Wavelength',
  'Frequency times wavelength. It holds for any wave at all, which is why it links the optical, '
  'acoustic and quantum descriptions of the same travelling disturbance.'),
 ("Ohm's law",
  'Electromagnetism',
  'V = I R',
  None,
  'Ohm%27s_law',
  'Voltage equals current times resistance. It is a property of certain materials rather than a '
  'law of nature, which is why semiconductors ignore it.'),
 ('Electrical power',
  'Electromagnetism',
  'P = V I',
  None,
  'Electric_power',
  "Voltage times current. Combined with Ohm's law it gives the two forms every electronics text "
  'uses next, P = I^2 R and P = V^2 / R.'),
 ('Capacitor charge',
  'Electromagnetism',
  'Q = C V',
  None,
  'Capacitance',
  'Charge stored is capacitance times voltage. The definition of capacitance, really, dressed as a '
  'result.'),
 ('Density',
  'Fluid dynamics',
  '\\rho = \\frac{m}{V}',
  None,
  'Density',
  'Mass per unit volume. Worth stating explicitly because rho is among the most overloaded letters '
  'in physics, and this is the reading the others are named after.'),
 ('Hydrostatic pressure',
  'Fluid dynamics',
  'P = \\rho g h',
  None,
  'Hydrostatic_equilibrium',
  'The pressure at a depth is the weight of the column above it. It does not depend on the shape '
  'of the container, which surprised people for a long time.'),
 ('Thermal energy',
  'Thermodynamics',
  'E = k_B T',
  None,
  'Thermal_energy',
  'The energy scale a temperature corresponds to, per degree of freedom. This is what the '
  'Boltzmann constant is for: it is a unit conversion between kelvin and joules, and in natural '
  'units it is simply one.'),
 ('Heat capacity',
  'Thermodynamics',
  'Q = m c_p \\Delta T',
  None,
  'Heat_capacity',
  'How much heat a temperature change costs. The specific heat c_p is where the material comes in, '
  'and it is one more meaning for an already overloaded c.'),
 ('Compton wavelength',
  'Quantum mechanics',
  '\\lambda_C = \\frac{h}{m c}',
  None,
  'Compton_wavelength',
  "The wavelength of a photon whose energy equals a particle's rest energy -- the scale at which "
  'quantum mechanics and relativity have to be used together, and below which the idea of a single '
  'particle stops working.'),
 ('Planck length',
  'Cosmology',
  '\\ell_P = \\sqrt{\\frac{\\hbar G}{c^3}}',
  None,
  'Planck_length',
  'The only length that can be built from hbar, G and c alone. That it exists at all is the '
  'argument that quantum gravity has a scale; what happens there is not something the formula '
  'knows.'),
 ('Escape velocity',
  'Classical mechanics',
  'v_e = \\sqrt{\\frac{2 G M}{r}}',
  None,
  'Escape_velocity',
  'The speed at which kinetic energy exactly cancels gravitational binding. Set it equal to c and '
  'the Schwarzschild radius drops out, which is a coincidence of the Newtonian derivation rather '
  'than a proof of it.'),
 ("Hubble's law",
  'Cosmology',
  'v = H_0 D',
  None,
  'Hubble%27s_law',
  'Recession speed grows with distance. The constant of proportionality is an inverse time, and '
  'that time is roughly the age of the universe -- which is the whole reason the law mattered.'),
 ('Coulomb potential energy',
  'Electromagnetism',
  'U = \\frac{1}{4 \\pi \\epsilon_0} \\frac{q_1 q_2}{r}',
  None,
  'Electric_potential_energy',
  'The work needed to bring two charges together from far apart. One power of r rather than two: '
  'the force is the derivative of this, and differentiating is where the second power comes from.'),
 ('Golden ratio',
  'Tilings',
  '\\phi = \\frac{1 + \\sqrt{5}}{2}',
  None,
  'Golden_ratio',
  'The positive root of x^2 = x + 1, and the inflation factor of the Penrose tilings. Its '
  'irrationality is not a curiosity here: a tiling with a translational period has only finitely '
  'many tiles per cell, so its tile ratios are rational, and a tiling whose ratio is this number '
  'therefore cannot have one.'),
 ('Inflation quadratic',
  'Tilings',
  '\\lambda^2 = \\lambda + 1',
  None,
  'Substitution_tiling',
  'The characteristic equation of the Fibonacci substitution matrix, and the defining property of '
  'the golden ratio: a length that reproduces itself plus a copy of what came before. The '
  'eigenvalue is what a substitution rule multiplies lengths by, and whether it is rational '
  'decides whether the tiling can repeat. Note the collision: lambda is a wavelength three entries '
  'above this one.'),
 ('Penrose tile ratio',
  'Tilings',
  'N_f = \\phi N_t',
  None,
  'Penrose_tiling',
  'In any Penrose tiling of the plane, fat rhombs outnumber thin ones by the golden ratio in the '
  'limit of a large patch. The counts themselves are unbounded and mean nothing; the ratio is the '
  'invariant, and it is the same for every one of the uncountably many distinct tilings.'),
 ('Inflation of area',
  'Tilings',
  "A' = \\phi^2 A",
  None,
  'Penrose_tiling',
  'One step of inflation scales lengths by the golden ratio and so scales areas by its square. '
  'This is where the tilings meet the rest of this group: an inflation rule is a statement about '
  'how area is budgeted, which is also what a minimal surface and a string worldsheet are about.'),
 ('Spectre inflation quadratic',
  'Tilings',
  '\\mu^2 = 8 \\mu - 1',
  'brentharts/spectre: chirality_e8.substitution_matrix(), nine species',
  'Einstein_problem',
  'The defining relation of the spectre inflation factor, and the exact analogue of x^2 = x + 1 '
  'for the Penrose tilings four entries above. It is the non-trivial factor of the substitution '
  "matrix's characteristic polynomial; the rest of that polynomial is a power of lambda and the "
  'two roots of unity. Same shape of argument as the golden ratio, different quadratic field.'),
 ('Spectre inflation factor',
  'Tilings',
  '\\mu = 4 + \\sqrt{15}',
  'brentharts/spectre: chirality_e8.substitution_matrix(), nine species',
  'Einstein_problem',
  'The larger root, about 7.873, and the factor by which one round of substitution multiplies the '
  'number of tiles. It is irrational, which is the whole obstruction: a periodic tiling has a '
  'rational ratio of tile counts and this number is not one.'),
 ('Spectre reciprocal pair',
  'Tilings',
  'g \\mu = 1',
  'brentharts/spectre: chirality_e8.substitution_matrix(), nine species',
  'Einstein_problem',
  'The two roots of the inflation quadratic multiply to one, because the constant term is one. So '
  'the smaller root is exactly the reciprocal of the larger, and the frequency of the commonest '
  'metatile species is exactly one over the inflation factor. That coincidence is not a '
  'coincidence: it is the constant term of the quadratic, read twice.'),
 ('Spectre triplet frequency',
  'Tilings',
  'g = 4 - \\sqrt{15}',
  'brentharts/spectre: chirality_e8.substitution_matrix(), nine species',
  'Aperiodic_tiling',
  'The limiting fraction of the tiling occupied by each of the three metatile species that every '
  'supertile contains exactly one of. Three species share this value, and they share it for a '
  'reason: the corresponding rows of the substitution matrix are all ones, which is a conservation '
  'law rather than an accident.'),
 ('Spectre frequency ratio',
  'Tilings',
  'v_\\Theta = \\frac{v_\\Gamma}{\\mu}',
  'brentharts/spectre: chirality_e8.substitution_matrix(), nine species',
  'Substitution_tiling',
  "One species appears exactly one inflation step behind another, so its frequency is the other's "
  'divided by the inflation factor. The counts bear this out literally: the sequence for the rarer '
  'species is the sequence for the commoner one shifted by a generation.'),
 ('Spectre count recurrence',
  'Tilings',
  "N'' = 8 N' - N",
  'brentharts/spectre: spectre_tile_counter.py, chirality_e8.substitution_matrix()',
  'Substitution_tiling',
  'The count of a tile species by generation obeys a second-order linear recurrence whose '
  'characteristic polynomial is the inflation quadratic -- for six of the nine species, exactly, '
  'from any starting tile once at most two generations have passed (the factor of lambda to the '
  'fifth in the characteristic polynomial is a transient, and it is gone by then). Started from a '
  'Delta the Gamma counts run 1, 8, 63, 496, 3905, 30744. The other three species, Xi, Pi and Psi, '
  'do not: the matrix also has the eigenvalues plus and minus one, and those leave a constant, or '
  'an alternating, offset in the recurrence for them, which the code handles with auxiliary '
  'sequences. So "tile counts obey this recurrence" is true of most species and false of three, '
  'and the ratio of consecutive terms tends to the inflation factor for all of them.'),
 ('Spectre tile area',
  'Tilings',
  'A = 2 \\sqrt{3} a^2 + 3 a b + \\sqrt{3} b^2',
  'brentharts/spectre: spectre.get_spectre_points()',
  'Einstein_problem',
  'The area of the fourteen-sided tile in terms of its two edge lengths. The expression is not '
  'symmetric in a and b, and that asymmetry is the whole reason the reflected tile is a different '
  'shape rather than the same one turned over.'),
 ('Hat tile area',
  'Tilings',
  'A_H = 8 \\sqrt{3}',
  'brentharts/spectre: spectre.get_spectre_points()',
  'Einstein_problem',
  'The hat is the tile with a = 1 and b = root three. Its sibling the turtle, with the lengths '
  'swapped, has area ten root three, and the spectre with both edges equal has three plus three '
  'root three -- three different areas from one formula, which is what makes the family a family.'),
 ('Spectre linear inflation',
  'Tilings',
  '\\ell = \\frac{\\sqrt{6} + \\sqrt{10}}{2}',
  'brentharts/nariai: nariai_constants.py, derived from the substitution matrix',
  'Substitution_tiling',
  'The factor by which one round of substitution multiplies a length, about 2.806. It is the '
  'square root of the area factor, and it lives in a quartic field rather than the quadratic one: '
  'its minimal polynomial is x^4 - 8x^2 + 1, and the closed form is a sum of two surds that look '
  'unrelated to four plus root fifteen until you square them. Which of the two inflation factors a '
  'calculation wants is not a detail -- lengths scale by this one and areas, tile counts and gap '
  'labels by its square.'),
 ('Spectre linear-area relation',
  'Tilings',
  '\\mu = \\ell^{2}',
  'brentharts/nariai: nariai_constants.py, derived from the substitution matrix',
  'Substitution_tiling',
  'Areas scale as the square of lengths, so the area inflation factor is the square of the linear '
  'one. The content is in what the two are used for: tile counts and the gap labels of the '
  'spectrum are powers of the area factor, and feeding the linear factor into a sum of those '
  'powers gives a wrong answer even when the linear factor itself is right.'),
 ('Spectre linear quartic',
  'Tilings',
  '\\ell^{4} = 8 \\ell^{2} - 1',
  'brentharts/nariai: nariai_constants.py, derived from the substitution matrix',
  'Substitution_tiling',
  'The quartic satisfied by the linear inflation factor. It is the inflation quadratic with the '
  'square of the unknown substituted for the unknown, which is the relation between the two '
  'factors read as a polynomial: the quartic has roots plus and minus the linear factor and plus '
  'and minus its reciprocal.'),
 ('Spectre gap-label recurrence',
  'Tilings',
  "G'' = 8 G' - G",
  'brentharts/nariai: nariai_constants.py, derived from the substitution matrix',
  'Substitution_tiling',
  'The labels of the gaps in the spectrum of an operator on the tiling live in the ring generated '
  'by the reciprocal of the area factor, and successive powers of that reciprocal obey the same '
  'second-order recurrence as the tile counts, which grow like powers of the factor itself. Both '
  'are forced by the one quadratic: it is reciprocal, so a root and its inverse satisfy the same '
  'relation. A fourth-order recurrence appears in an earlier version of the source, and was a '
  'consequence of using a number that is not the inflation factor.'),
 ('Spectre growth rate',
  'Tilings',
  '\\gamma = \\ln \\mu',
  'brentharts/nariai: nariai_constants.py, derived from the substitution matrix',
  'Substitution_tiling',
  'The rate, per round of substitution, at which the number of tiles grows: about 2.063 nats. It '
  'is a growth rate and not an entropy. A primitive substitution tiling is uniquely ergodic and '
  'the translation action on it has zero topological entropy, so quoting this number as one is a '
  'category error. It is also twice the logarithm of the linear factor, and the two are easy to '
  'swap, which doubles any period built from them.'),
 ('Hat inflation quadratic',
  'Tilings',
  '\\kappa^2 = 7 \\kappa - 1',
  'brentharts/CICY: pyCICY.monotile.SUBSTITUTION, four metatiles',
  'Einstein_problem',
  "The defining relation of the hat inflation factor. One coefficient away from the spectre's "
  "quadratic above, and that one coefficient is the whole difference between the two tilings' "
  'arithmetic. The rest of the characteristic polynomial is the same pair of roots of unity in '
  'both.'),
 ('Hat inflation factor',
  'Tilings',
  '\\kappa = \\frac{7 + 3 \\sqrt{5}}{2}',
  'brentharts/CICY: pyCICY.monotile.inflation_factor()',
  'Einstein_problem',
  'The larger root, about 6.854, which is the fourth power of the golden ratio. So the hat tiling '
  'is governed by the same constant as the Penrose tilings after all, four powers up -- while the '
  'spectre, which is the same shape family, is governed by root fifteen instead.'),
 ('Hat constant as a golden power',
  'Tilings',
  '\\kappa = \\phi^4',
  'brentharts/CICY: pyCICY.monotile.hat_chirality(), is_phi4',
  'Golden_ratio',
  'The hat inflation factor is the fourth power of the golden ratio, so the constant governing the '
  'one-tile system is the same one that governs the two-tile system Penrose found fifty years '
  'earlier. This is the entry that ties the two halves of this group together, and being written '
  'down is what lets the library check it: both sides are pinned to numbers elsewhere, so it is '
  'one of the few statements here that can simply be evaluated.'),
 ('Hat chirality ratio',
  'Tilings',
  '\\frac{N_u}{N_r} = \\kappa',
  'brentharts/CICY: pyCICY.monotile.hat_chirality()',
  'Einstein_problem',
  'Unreflected hats outnumber reflected ones by exactly the inflation factor: about one anti-hat '
  'per 6.854 hats. The tiling is therefore almost, but not quite, chiral, and the imbalance is '
  'irrational -- which is the same argument as the Penrose rhombs, since a periodic tiling would '
  'force the ratio to be rational.'),
 ('Euler characteristic of a closed surface',
  'Extremal geometry',
  '\\chi = 2 - 2 g',
  None,
  'Euler_characteristic',
  'Two for a sphere, nought for a torus, and down by two for every handle added. Nothing short of '
  'cutting the surface changes it, which is what makes it worth computing -- and in string '
  'perturbation theory the genus counting handles is also the genus counting loops.'),
 ('Riemann-Hurwitz (unbranched)',
  'Extremal geometry',
  '\\chi_C = n \\chi_B',
  None,
  'Riemann%E2%80%93Hurwitz_formula',
  'An n-sheeted cover of a surface has n times its Euler characteristic, provided nothing is '
  'branched. The full formula subtracts a term for each ramification point; this is the case where '
  'that term is absent, which is the case that says a torus can cover a torus and nothing can '
  'cover a sphere except a sphere.'),
 ('Gauss-Bonnet theorem',
  'Extremal geometry',
  '\\int_M K \\, dA + \\oint_{\\partial M} k_g \\, ds = 2 \\pi \\chi',
  None,
  'Gauss%E2%80%93Bonnet_theorem',
  'Total curvature is topology. The left side is measured with a ruler and depends on the shape; '
  'the right side is a whole number that does not. Bend the surface however you like and the '
  'curvature moves around without the integral changing, which is the sharpest statement in this '
  'group of what a constraint on geometry costs.'),
 ('Mean curvature',
  'Extremal geometry',
  'H = \\frac{\\kappa_1 + \\kappa_2}{2}',
  None,
  'Mean_curvature',
  'The average of the two principal curvatures at a point. Note the collision: H is a Hamiltonian '
  'everywhere else in this library.'),
 ('Minimal surface condition',
  'Extremal geometry',
  'H = 0',
  None,
  'Minimal_surface',
  'A surface is minimal exactly when its mean curvature vanishes everywhere. It does not say the '
  'surface has least area outright -- only that area is stationary under small deformations, which '
  'is the same distinction as between a minimum and a turning point. A soap film finds one of '
  'these because surface tension makes area the energy.'),
 ('Minimal surface equation',
  'Extremal geometry',
  '(1 + u_y^2) u_{xx} - 2 u_x u_y u_{xy} + (1 + u_x^2) u_{yy} = 0',
  None,
  'Minimal_surface',
  'The mean-curvature condition written out for a surface given as a height above the plane. It is '
  'the Euler-Lagrange equation of the area functional, and it is quasilinear rather than linear -- '
  'the coefficients contain the slopes of the very solution being sought, which is why soap films '
  "are hard and Laplace's equation is easy."),
 ('Catenoid',
  'Extremal geometry',
  'r = c \\cosh \\left( \\frac{z}{c} \\right)',
  None,
  'Catenoid',
  'The surface swept by rotating a catenary, and after the plane the first minimal surface anyone '
  'found. Stretch the two rings holding the soap film too far apart and no catenoid joins them: '
  'the solution stops existing and the film snaps to two discs.'),
 ('Area functional',
  'Extremal geometry',
  'A = \\int \\int \\sqrt{1 + |\\nabla u|^2} \\, dx \\, dy',
  None,
  'Plateau%27s_problem',
  "The area of a surface given as a height above the plane. Plateau's problem is to minimise this "
  'with the boundary held fixed, and the minimal surface equation is what comes out of setting its '
  'first variation to zero. The square root is the whole difficulty.'),
 ('Kakeya conjecture',
  'Extremal geometry',
  '\\dim_H K = n',
  None,
  'Kakeya_set',
  'A Kakeya set contains a unit segment in every direction and can have measure zero, which is '
  "Besicovitch's construction. The conjecture is that it cannot be small in the dimensional sense "
  'as well: its Hausdorff dimension must be the full n. Settled in the plane long ago, and in '
  'three dimensions by Wang and Zahl in 2025.'),
 ('Perron tree bound',
  'Extremal geometry',
  'A \\leq \\frac{C}{\\log n}',
  None,
  'Kakeya_set',
  'Bisect a triangle, slide the halves together so they overlap, and repeat: after n steps the '
  'figure still points in all the original directions but its area has fallen like one over log n. '
  'The decay is slow, and unbounded, which is how a needle turns around in arbitrarily little '
  'area.'),
 ('Nambu-Goto action',
  'String theory',
  'S = -T \\int d^2\\sigma \\sqrt{-\\det(h_{ab})}',
  None,
  'Nambu%E2%80%93Goto_action',
  'The action of a string is the area of the sheet it sweeps through spacetime, times its tension. '
  'That is the whole content: a point particle extremises its worldline length, a string '
  'extremises its worldsheet area, and everything else in the subject is consequence. It is also '
  'why minimal surfaces belong in the same group as this.'),
 ('String tension',
  'String theory',
  "T = \\frac{1}{2 \\pi \\alpha'}",
  None,
  'String_(physics)',
  'Tension and the Regge slope are the same parameter inverted. Note the collision: T is a '
  'temperature and a period elsewhere in this library, and alpha prime is not alpha -- the prime '
  'is part of the name, not an operation.'),
 ('Regge trajectory',
  'String theory',
  "J = \\alpha' M^2 + \\alpha_0",
  None,
  'Regge_theory',
  'Spin against mass squared, and the points lie on a line. This was measured in hadron spectra '
  'before anyone proposed a string, and explaining the straightness is what strings were invented '
  'to do: a spinning relativistic string has exactly this relation between its angular momentum '
  'and its energy.'),
 ('String length',
  'String theory',
  "\\ell_s = \\sqrt{\\alpha'}",
  None,
  'String_(physics)',
  'The one length the theory has, in natural units. Every statement about strings being small is a '
  'statement about this number, and nothing in the theory predicts it.'),
 ('Genus expansion',
  'String theory',
  'Z = \\sum_g g_s^{-\\chi} Z_g',
  None,
  'String_theory',
  'The perturbation series of string theory is a sum over surfaces, ordered by how many handles '
  'they have. The exponent is minus the Euler characteristic, so by chi = 2 - 2g each extra handle '
  'costs two powers of the coupling. A loop expansion and a classification of surfaces turn out to '
  'be the same list.'),
 ('Critical dimension',
  'String theory',
  'D = 26',
  None,
  'Critical_dimension',
  'The bosonic string is only consistent in twenty-six spacetime dimensions; the superstring, in '
  'ten. The number is not chosen. It is what makes a quantum anomaly cancel, and getting a '
  'specific integer out of a consistency condition rather than out of observation is the most '
  'unusual thing about the subject.'),
 ('Hagedorn temperature',
  'String theory',
  "T_H = \\frac{1}{4 \\pi \\sqrt{\\alpha'} k_B}",
  None,
  'Hagedorn_temperature',
  'Heat a gas of strings and the number of available states grows exponentially with energy, so '
  'above this temperature the partition function diverges and temperature stops meaning anything. '
  'It is where this group rejoins ordinary physics: the Boltzmann constant appears, and the '
  'quantity on the left is the same temperature as in the gas laws.'),
 ('Lattice coupling',
  'Lattice gauge theory',
  '\\beta = \\frac{2 N}{g^2}',
  None,
  'Lattice_gauge_theory',
  'The one parameter of a pure gauge lattice, written as an inverse coupling: large beta is weak '
  'coupling, and the continuum limit is beta going to infinity. For SU(2) it is 4 / g^2. The '
  'laboratory this group was drawn from treats beta as given and scans it, which is the honest '
  'thing to do on a lattice: the bare coupling is an input, and what it means in physical units is '
  'a result.'),
 ('Wilson plaquette action',
  'Lattice gauge theory',
  'S_W = \\beta \\sum_p \\left( 1 - \\frac{1}{N} \\Re \\mathrm{Tr}\\, U_p \\right)',
  'Elsbernd09/yang-mills-mass-gap-lab: ymlab.wilson_action.wilson_action(), N = 2',
  'Wilson_loop',
  'Sum, over every elementary square of the lattice, of how far the product of the four links '
  'around it is from the identity. A smooth field has plaquettes near one and costs nothing; a '
  'rough one is penalised. Written with the link variables rather than a field strength, which is '
  'the point: gauge invariance is exact, not recovered in a limit, because the trace of a closed '
  'loop does not change under a gauge transformation.'),
 ('Yang-Mills action',
  'Lattice gauge theory',
  'S = \\frac{1}{4 g^2} \\int d^4x \\, F^a_{\\mu\\nu} F^{a\\,\\mu\\nu}',
  None,
  'Yang%E2%80%93Mills_theory',
  'The continuum action the lattice one approximates: the square of the field strength, '
  "integrated, with the coupling out front. It looks like Maxwell's with an index a on F, and the "
  'difference is entirely in what F contains: for a non-Abelian group it has a term quadratic in '
  'the potential, so the field interacts with itself, and that is the whole difficulty of the '
  'mass-gap problem.'),
 ('Link gauge transformation',
  'Lattice gauge theory',
  'U_\\mu(x) \\to G(x) U_\\mu(x) G^\\dagger(x + \\hat{\\mu})',
  'Elsbernd09/yang-mills-mass-gap-lab: ymlab.gauge_transformations.transform_link()',
  'Gauge_theory',
  'A link variable lives between two sites, so a gauge transformation acts at both ends: G where '
  'it starts and the inverse of G where it ends. Every observable worth computing is unchanged by '
  'this, and a program can test that directly -- transform the whole lattice at random and see '
  'whether the action moves. It is the check that catches a wrong sign in a staple.'),
 ('Creutz ratio',
  'Lattice gauge theory',
  '\\chi = -\\ln \\frac{W_{RT} W_{R-1,T-1}}{W_{R,T-1} W_{R-1,T}}',
  'Elsbernd09/yang-mills-mass-gap-lab: ymlab.creutz.creutz_ratio_from_values()',
  'Wilson_loop',
  'A ratio of four Wilson loops built so that the perimeter and corner contributions cancel and '
  'only the area term survives. If the loops fall off with their area it is the string tension; if '
  'they fall off with their perimeter it is zero. It needs both sides of the ratio positive, which '
  'noisy loops often are not, and an implementation has to say what it does then rather than take '
  'a logarithm of noise.'),
 ('Area law',
  'Lattice gauge theory',
  'W = \\exp(-\\sigma A)',
  'Elsbernd09/yang-mills-mass-gap-lab: tests/test_creutz.py '
  'test_creutz_ratio_from_exact_area_law()',
  'Wilson_loop',
  'A Wilson loop that decays with the area it encloses, which is the lattice criterion for '
  'confinement: the energy of a pair of static charges grows linearly with their separation. The '
  'constant in the exponent is the confining string tension. A synthetic loop built to obey this '
  'exactly is what a Creutz-ratio routine has to recover to be trusted on real data.'),
 ('Creutz string tension',
  'Lattice gauge theory',
  '\\chi = \\sigma a^2',
  'Elsbernd09/yang-mills-mass-gap-lab: ymlab.creutz.estimate_string_tension()',
  'Wilson_loop',
  'For an exact area law the Creutz ratio is the string tension in lattice units, and multiplying '
  'by the lattice spacing squared is what puts it back into the units of an inverse area. The '
  'laboratory this comes from averages valid Creutz ratios and calls the result "string-tension '
  'like", which is the right amount of caution: it is an estimate on one finite lattice at one '
  'bare coupling.'),
 ('Lattice effective mass',
  'Lattice gauge theory',
  'm_{eff} = \\ln \\frac{C_{t}}{C_{t+1}}',
  'Elsbernd09/yang-mills-mass-gap-lab: ymlab.mass_gap.effective_mass()',
  'Effective_mass_(solid-state_physics)',
  'The mass a correlator seems to have between one time slice and the next, in units of the '
  'lattice spacing. If one state dominates it is flat; if several do it drifts down towards the '
  'lightest as t grows. Reading the plateau is the whole art of lattice spectroscopy, and a '
  'plateau is a claim about a finite lattice, not about the continuum.'),
 ('Euclidean correlator decay',
  'Lattice gauge theory',
  'C_{t} = A_0 \\exp(-m_{L} t)',
  'Elsbernd09/yang-mills-mass-gap-lab: ymlab.mass_gap.exponential_decay()',
  'Correlation_function_(quantum_field_theory)',
  'In Euclidean time a correlator of a gauge-invariant operator falls off exponentially, at a rate '
  'set by the lightest state the operator couples to. That rate is the mass, and a positive one '
  'for the lightest state is what a mass gap means. Here it is a fit parameter on a periodic '
  'lattice of a few sites, and nothing more.'),
 ('Mass gap in physical units',
  'Lattice gauge theory',
  'E_{gap} = \\frac{\\hbar c \\, m_{L}}{a}',
  None,
  'Mass_gap',
  'A mass measured on the lattice is a pure number, because the lattice has no units but its own '
  'spacing. Putting it in joules needs the spacing in metres, and the spacing is not an input: it '
  'is fixed afterwards by matching one measured quantity. That is why a lattice number without its '
  'scale-setting is not yet a physical statement, and why the continuum mass gap is a limit, not a '
  'value.'),
 ('Integrated autocorrelation time',
  'Lattice gauge theory',
  '\\tau = \\frac{1}{2} + \\sum_{t=1}^{t_c} \\rho_t',
  'Elsbernd09/yang-mills-mass-gap-lab: ymlab.diagnostics.integrated_autocorrelation_time()',
  'Autocorrelation',
  'How many updates of a Markov chain it takes to produce one independent sample, with the sum cut '
  'off where the estimate would start adding noise. The half is the lag zero term. Slow '
  'decorrelation near the continuum limit is the reason lattice simulation is expensive, and the '
  'reason there is a literature on update algorithms.'),
 ('Effective sample size',
  'Lattice gauge theory',
  'N_{eff} = \\frac{N}{2 \\tau}',
  'Elsbernd09/yang-mills-mass-gap-lab: ymlab.diagnostics.effective_sample_size()',
  'Effective_sample_size',
  'A chain of N correlated samples carries as much information about the mean as N over two tau '
  'independent ones. It is what turns "I ran a million sweeps" into a statement about error bars, '
  'and it is the numerator of the efficiency below.'),
 ('Sampling efficiency',
  'Lattice gauge theory',
  '\\eta = \\frac{N_{eff}}{t_{wall}}',
  'Elsbernd09/yang-mills-mass-gap-lab: ymlab.research_campaign.run_schedule_chain()',
  'Effective_sample_size',
  'Independent samples per second of wall-clock time. It divides the same effective sample size by '
  'a measured runtime, and the runtime is the part that makes it machine dependent: a different '
  'computer gives a different number, which is why the autocorrelation time reproduces exactly and '
  'this does not.'),
 ('Efficiency ratio',
  'Lattice gauge theory',
  'R = \\frac{\\eta_{OR}}{\\eta_{M}}',
  'Elsbernd09/yang-mills-mass-gap-lab: ymlab.research_campaign.summarize_campaign()',
  'Overrelaxation',
  'The efficiency of a hybrid schedule against plain Metropolis, taken pair by pair from the same '
  'starting configuration. Above one means the extra overrelaxation sweeps paid for themselves. '
  'The finding this group rests on is that it can be below one while the autocorrelation time is '
  'also lower: fewer updates per independent sample is not the same as less time.'),
 ('O(N) lattice model',
  'Lattice spin models',
  'H = -\\beta \\sum q_{x} \\cdot q_{y}',
  'brentharts/openai_math_yang_mills: tools/on_model_mc.py',
  'N-vector_model',
  'Unit vectors on the sites of a lattice, each pair of neighbours rewarded for pointing the same '
  'way, with the Gibbs weight exp(-H). For three or more components in two dimensions it has no '
  'ordered phase and a finite correlation length at every temperature, which is the property that '
  'makes it a toy for a gauge theory.'),
 ('Asymptotic-freedom mass scaling',
  'Lattice spin models',
  'm = C \\beta^{\\frac{1}{N-2}} \\exp\\left(-\\frac{2 \\pi \\beta}{N-2}\\right)',
  'brentharts/openai_math_yang_mills: tools/on_model_mc.py',
  'Asymptotic_freedom',
  'The mass of the O(N) model in lattice units, as asymptotic freedom says it should behave at '
  'large beta. Only the form is predicted: the constant C is not, and it is corrected by terms in '
  'one over beta that are large at the values simulations reach. A simulated ratio of correlation '
  'lengths approaches the form from below for a long time, which is a statement about the approach '
  'and not a failure of the theory.'),
 ('Second-moment correlation length',
  'Lattice spin models',
  '\\xi^2 = \\frac{\\frac{\\chi}{F} - 1}{p^2}',
  'brentharts/openai_math_yang_mills: tools/on_model_mc.py',
  'Correlation_function_(statistical_mechanics)',
  'A correlation length read off a simulation from two numbers: the susceptibility, and the '
  'structure factor at the smallest non-zero momentum on a periodic box, whose lattice momentum is '
  'p = 2 sin(pi / L) for a box of L sites. It is a mass in the sense of an inverse length, but not '
  'the gap of the transfer matrix; the two differ by a factor of order one, so a comparison of '
  'shapes is honest and a comparison of constants is not.'),
 ('Kosterlitz coupling flow',
  'Lattice spin models',
  '\\frac{dx}{dl} = -y^2',
  'brentharts/openai_math_yang_mills: tools/kt_flow.py',
  'Kosterlitz%E2%80%93Thouless_transition',
  'How the stiffness of the planar XY model drifts from the universal-jump value as the length '
  'scale grows, driven by the density of vortices. With the next entry it is the whole '
  'renormalisation-group picture of the transition, and it conserves x^2 - y^2.'),
 ('Kosterlitz fugacity flow',
  'Lattice spin models',
  '\\frac{dy}{dl} = -x y',
  'brentharts/openai_math_yang_mills: tools/kt_flow.py',
  'Kosterlitz%E2%80%93Thouless_transition',
  'How the vortex fugacity runs: it dies away where the stiffness is high and grows where it is '
  'low. The correlation length diverges with an essential singularity at the transition, not a '
  'power, because the flow takes a time proportional to one over the distance from the critical '
  'line to leave it.'),
 ('Magnon dispersion',
  'Lattice spin models',
  '\\epsilon = 2 \\sum \\left(1 - \\cos k_{a}\\right)',
  'brentharts/openai_math_yang_mills: tools/bloch_law.py',
  'Magnon',
  'The energy of one spin flip spread over a cubic ferromagnet, in units of the exchange. It '
  'vanishes at zero wave number as k squared, and that quadratic edge is why the number of '
  'thermally excited magnons grows like a power of temperature.'),
 ('Bloch magnon density',
  'Lattice spin models',
  'n = A_{0} x^{-\\frac{3}{2}} + A_{1} x^{-\\frac{5}{2}}',
  'brentharts/openai_math_yang_mills: tools/bloch_law.py',
  'Bloch%27s_law',
  'The ideal-magnon density at large inverse temperature: a Bloch term and the first lattice '
  "correction. The leading coefficient is Bloch's, and the next is 3 zeta(5/2) / (128 pi^(3/2)). "
  'It is the density of non-interacting magnons only; that it controls the magnetisation deficit '
  'is what the manuscript proves, and what this entry does not.'),
 ('Bilinear-biquadratic chain',
  'Lattice spin models',
  'H = \\sum \\left[ \\cos\\theta \\, S_{i} \\cdot S_{i+1} + \\sin\\theta \\, \\left(S_{i} \\cdot '
  'S_{i+1}\\right)^2 \\right]',
  'brentharts/openai_math_yang_mills: tools/haldane_gap.py',
  'Heisenberg_model_(quantum)',
  'A ring of spin-one sites with a Heisenberg coupling and a biquadratic one, mixed by an angle. '
  'Theta zero is the Heisenberg chain; the tangent of theta equal to a third is the AKLT point, '
  'whose ground state is known exactly, which is what makes a numerical check of the code '
  'possible.'),
 ('Spin gap',
  'Lattice spin models',
  '\\Delta = E_{1} - E_{0}',
  'brentharts/openai_math_yang_mills: tools/haldane_gap.py',
  'Spin_gap',
  'The energy to flip one unit of total spin above the ground state. For the spin-one chain on '
  'rings of growing length it approaches a positive number from above, which is the Haldane gap; a '
  'list of ring sizes is numerical evidence for that and not a bound, and the extrapolation '
  'assumes a finite-size form.'),
 ('One-loop running',
  'Gauge unification',
  '\\frac{1}{\\alpha_{i}} = \\frac{1}{\\alpha_{G}} + \\frac{b_{i}}{2 \\pi} \\ln\\frac{M_{G}}{\\mu}',
  'brentharts/CICY: pyCICY.theories.running.beta_coefficients()',
  'Renormalization_group',
  'The inverse of a gauge coupling at one loop: the value at the unification scale plus a slope '
  'times the logarithm of the ratio of scales. The slope b is a sum over the particles in the '
  'spectrum and is the only place the spectrum enters. Negative b is asymptotic freedom.'),
 ('Unification prediction',
  'Gauge unification',
  '\\frac{1}{\\alpha_{3}} = \\frac{b_{3} - b_{2}}{b_{1} - b_{2}} \\frac{1}{\\alpha_{1}} + '
  '\\frac{b_{1} - b_{3}}{b_{1} - b_{2}} \\frac{1}{\\alpha_{2}}',
  'brentharts/CICY: pyCICY.theories.couplings.predict_alpha3()',
  'Grand_Unified_Theory',
  'What is left of three running couplings and one meeting point once the meeting scale and value '
  'are eliminated: the strong coupling from the other two. The supersymmetric spectrum lands '
  'within about one per cent at one loop and the non-supersymmetric one is out by forty. It is '
  'blind to the number of generations, because a generation adds the same amount to every b and '
  'only differences enter.'),
 ('Dimensional transmutation',
  'Gauge unification',
  '\\Lambda = M \\exp\\left(\\frac{2 \\pi}{b_{3} \\alpha}\\right)',
  'brentharts/CICY: pyCICY.theories.running.lambda_qcd()',
  'Dimensional_transmutation',
  'The scale at which a running coupling diverges, set by its value at a higher scale and the '
  'one-loop slope. With b negative the exponent is negative and the scale is far below M. Almost '
  "all of a proton's mass is this scale, not a Yukawa coupling. Only the exponent's slope is exact "
  'in a compactification: the coupling and the scale are moduli.'),
 ('One-loop coefficient (QCD)',
  'Gauge unification',
  'b = -\\frac{11}{3} N + \\frac{2}{3} n_{f}',
  'brentharts/CICY: pyCICY.theories.couplings, Standard Model b_3 = -7',
  'Beta_function_(physics)',
  'The slope of the inverse coupling of SU(N) with n_f Dirac fundamentals, in the convention where '
  'negative means asymptotic freedom. Gauge bosons push it down by eleven thirds of N and each '
  'quark flavour pushes it back up by two thirds. With three colours it stays negative up to '
  'sixteen flavours.'),
 ('One-loop coefficient (SUSY)',
  'Gauge unification',
  'b = -3 N + n_{f}',
  'brentharts/CICY: pyCICY.theories.running.beta_coefficients()',
  'Beta_function_(physics)',
  'The same slope when the spectrum is supersymmetric and n_f counts chiral flavour pairs. Three '
  'colours and two flavours per generation give minus nine plus twice the number of generations, '
  'which stays negative for four generations and turns positive at five: no confinement, no '
  'proton.'),
 ('Electromagnetic coupling at the Z',
  'Gauge unification',
  '\\frac{1}{\\alpha_{em}} = \\frac{1}{\\alpha_{2}} + \\frac{5}{3} \\frac{1}{\\alpha_{1}}',
  'brentharts/CICY: pyCICY.theories.couplings.alpha_em_inverse_mz()',
  'Weinberg_angle',
  'Electromagnetism as the combination of the weak and hypercharge couplings that survives '
  'symmetry breaking, with the five thirds that converts the unification normalisation of '
  'hypercharge to the usual one. Computed at the Z mass from the two couplings it comes out within '
  'a tenth of a per cent of the measurement.'),
 ('Thomson limit',
  'Gauge unification',
  '\\frac{1}{\\alpha_{0}} = \\frac{1}{\\alpha_{em}} \\frac{1}{1 - \\epsilon}',
  'brentharts/CICY: pyCICY.theories.couplings.alpha_em_inverse_zero()',
  'Fine-structure_constant',
  'The step from the coupling at the Z mass to the fine-structure constant at zero momentum, '
  'through a shift epsilon, the Delta alpha of the literature. The shift is measured, not '
  'computed: its hadronic part comes from electron-positron data through a dispersion relation '
  'because the momenta are non-perturbative. Any derivation of 137 from pure geometry stops here.'),
 ('Generation count',
  'Gauge unification',
  'N_{gen} = - \\frac{I}{G}',
  'brentharts/CICY: pyCICY.breaking.generations_downstairs()',
  'Atiyah%E2%80%93Singer_index_theorem',
  'The net number of chiral generations downstairs: minus the index of the bundle upstairs divided '
  'by the order of the freely acting group that quotients the manifold. An integer from topology, '
  'with no metric anywhere. It counts only the chiral part; vector-like matter, which the '
  'couplings are sensitive to, does not appear in an index.'),
 ('Heterotic anomaly condition',
  'Gauge unification',
  'W = c_{2}(TX) - c_{2}(V)',
  'brentharts/CICY: pyCICY.bundles.anomaly()',
  'Green%E2%80%93Schwarz_mechanism',
  'What the bundle leaves of the second Chern class of the tangent bundle has to be an effective '
  'class, the class of any five-branes wrapped to make up the difference. Checking it integrated '
  'against the Kahler classes is necessary and not sufficient, and the code says which one it '
  'did.'),
 ('Lie algebra dimension',
  'Lie algebras and chirality',
  'D = R + r',
  'emad-ii/chirality-standard-model: verification/python/check_e6_pair_orderings.py',
  'Root_system',
  'A semisimple Lie algebra is its Cartan subalgebra, one dimension for each unit of rank r, plus '
  'one dimension for every root, so the dimension is the root count plus the rank: E8 has 240 '
  'roots and rank 8, which is 248, and E6 has 72 and 6, which is 78. The chirality repository uses '
  'it to rank the three pair types inside E6. It is for semisimple algebras only; a central torus '
  'adds to the dimension and is not a root.'),
 ('Adjoint branching',
  'Lie algebras and chirality',
  'D = D_{h} + 2 D_{V}',
  'brentharts/chirality-standard-model: lisi/check_e6_su3_bridge.py',
  'Branching_rule',
  'When a compact algebra contains an effective subalgebra whose complement is one irreducible of '
  'complex type, the complement has real dimension twice the complex dimension of V, because its '
  'complexification is V plus its dual. For E8 containing E6 times SU(3) that is 248 = 86 + 2 '
  'times 81, where 86 is 78 + 8 and the 162 is the 27 over the 3 and its conjugate. A dimension '
  'count is a necessary condition on the pair and not a proof that it exists.'),
 ('Product representation dimension',
  'Lie algebras and chirality',
  'D_{V} = D_{a} D_{b}',
  'brentharts/chirality-standard-model: lisi/check_e6_su3_bridge.py',
  'Tensor_product_of_representations',
  'The external tensor product of a representation of one factor with a representation of another '
  'has the product of the two dimensions. The family module is the 27 of E6 times the 3 of SU(3), '
  'which is 81.'),
 ('Orbit size',
  'Lie algebras and chirality',
  'n = \\frac{w}{s}',
  'emad-ii/chirality-standard-model: verification/certificates/e6_extrema_certificate.json',
  'Group_action',
  'The orbit-stabiliser theorem: the orbit of a point has the order of the group divided by the '
  "order of the point's stabiliser. The Weyl group of E6 has order 51840 and the maximal "
  'subsystems D5, A5+A1 and A2 cubed have stabilisers of orders 1920, 1440 and 1296, so their '
  'orbits hold 27, 36 and 40 subsystems.'),
 ('Triple space',
  'Lie algebras and chirality',
  'T = n_{1} n_{2} n_{3}',
  'emad-ii/chirality-standard-model: verification/certificates/e6_extrema_certificate.json',
  'Cartesian_product',
  'One subsystem from each of the three orbits, chosen independently: 27 times 36 times 40 is '
  '38880 ordered triples. The repository exhausts all of them, which is why its claims about '
  'extremal intersections are statements about a finite set and not about a sample.'),
 ('Magic square dimension',
  'Lie algebras and chirality',
  'L = \\tau_{A} + \\tau_{B} + 3 a b',
  'brentharts/chirality-standard-model: lisi/check_magic_square.py',
  'Freudenthal_magic_square',
  'The dimension of the Freudenthal-Tits Lie algebra built from two division algebras: the '
  'dimensions of their triality algebras, plus three times the product of the dimensions of the '
  'algebras. With both the octonions it is 28 + 28 + 3 times 64 = 248, the D4 times D4 count of '
  'E8. The repository computes the triality dimensions 0, 2, 9, 28 as exact nullspaces; the '
  'formula is the standard one and the table it reproduces is the standard table.'),
 ('Cubic anomaly coefficient',
  'Lie algebras and chirality',
  'A_{abc} = \\frac{1}{2} \\mathrm{tr}\\left(T_{a} \\{T_{b}, T_{c}\\}\\right)',
  'emad-ii/chirality-standard-model: verification/python/check_cubic_radical_selection.py',
  'Gauge_anomaly',
  'The group factor in the triangle anomaly of a left-handed Weyl multiplet: the symmetrised trace '
  'of three generators in the representation. It is odd under taking the dual, zero on a self-dual '
  "representation and zero for every representation of SU(2). The repository's second condition is "
  'that the ideal on which it vanishes be nonzero, which E6, with its complex 27 and no cubic '
  'invariant, satisfies.'),
 ('Anomaly of a product representation',
  'Lie algebras and chirality',
  'A_{V} = D_{a} A_{b}',
  'emad-ii/chirality-standard-model: paper/paper.tex, eq. (11)',
  'Gauge_anomaly',
  'Generators acting through the second factor of a tensor product act as the same matrix once for '
  'each basis state of the first, so the trace is multiplied by its dimension: the SU(3) cubic '
  'anomaly of the 27 times 3 is 27 times that of the 3. The repository calls the SU(3) a '
  'background symmetry and says that making it a gauge group needs further anomaly-cancelling '
  'structure.'),
 ('Charge lattice',
  'Lie algebras and chirality',
  'q + 2 t + 3 s \\equiv 0 \\pmod{6}',
  'emad-ii/chirality-standard-model: paper/paper.tex, section on the 27 and the Standard Model',
  'Standard_Model',
  'Which representations descend to the faithful Standard Model group S(U(3) x U(2)): with the '
  'integer charge q normalised so the quark doublet has q = 1, colour triality t in Z/3 and '
  'weak-centre parity s in Z/2, a representation exists on the group exactly when q + 2t + 3s is a '
  'multiple of six. With t = 0 it forces the electric charges of colour singlets to be integers.'),
 ('Hypercharge',
  'Lie algebras and chirality',
  'Y = \\frac{q}{6}',
  'emad-ii/chirality-standard-model: paper/paper.tex, section on the 27 and the Standard Model',
  'Weak_hypercharge',
  'The conventional hypercharge is the integer charge of the lattice above divided by six, which '
  'is what makes the quark doublet 1/6, the right-handed up antiquark -2/3 and the positron 1.'),
 ('Hypercharge normalisation',
  'Lie algebras and chirality',
  'k = \\frac{\\sum Y^{2}}{\\sum T_{3}^{2}}',
  'emad-ii/chirality-standard-model: paper/paper.tex, section on the 27 and the Standard Model',
  'Weinberg_angle',
  'The ratio of the sum of squared hypercharges to the sum of squared weak isospins over a '
  'complete multiplet. For one family of left-handed fermions, which is a 10 and a 5-bar of SU(5), '
  'it is five thirds -- the number that converts the unification normalisation of hypercharge to '
  'the usual one. Here it is a trace over a family and not an assumption, and it is the same with '
  'or without the singlet, which has zero hypercharge.'),
 ('Anomaly-free hypercharge',
  'Lie algebras and chirality',
  '\\sum_{i} Y_{i}^{3} = 0',
  'emad-ii/chirality-standard-model: paper/paper.tex, section on the 27 and the Standard Model',
  'Gauge_anomaly',
  'Over one family of left-handed Weyl fermions the cubes of the hypercharges cancel, and so do '
  'the hypercharges themselves and the mixed sums with SU(3) and SU(2). The paper derives this '
  'from E6 being the radical of the cubic form; the bridge adds the five summands directly and '
  'does not use that argument, so the two are independent.'),
 ('Net chiral class',
  'Lie algebras and chirality',
  '\\chi_{V} = N_{gen} \\chi_{F}',
  'emad-ii/chirality-standard-model: paper/paper.tex, eq. (10)',
  'Chirality_(physics)',
  "The net chirality of the family module as a multiple of one Standard Model family's: for each "
  "representation not equal to its dual, its multiplicity minus its conjugate's. Vector-like pairs "
  'and singlets drop out and the three that remains is the dimension of the SU(3) on the family '
  'index. It counts the chiral part only and fixes no mass or mixing.'),
 ('Special unitary dimension',
  'Lie algebras and chirality',
  'D = N^{2} - 1',
  None,
  'Special_unitary_group',
  'The dimension of SU(N): the traceless N by N Hermitian matrices, N squared less the one trace. '
  'It is three for the weak isospin and eight for colour. Together with the roots-plus-rank entry '
  'it says the roots of SU(N) number N squared less N, which is the count that makes its dual '
  'Coxeter number N.'),
 ('Dual Coxeter number (simply laced)',
  'Lie algebras and chirality',
  'h_{c} = \\frac{R}{r}',
  None,
  'Dual_Coxeter_number',
  'For a simply laced algebra the dual Coxeter number is the number of roots divided by the rank. '
  'E8 has 240 roots and rank 8, so thirty; E6 has 72 and 6, so twelve; SU(N) has N squared less N '
  "roots and rank N less one, so N. It is the integer that sets the gauge bosons' contribution to "
  'the running, and it is why the exceptional groups run faster than SU(3) does: their couplings '
  'change more for the same change of scale. The quotient is exact only when every root has the '
  'same length; for the other algebras it is a different integer.'),
 ('Adjoint index',
  'Gauge unification',
  'T_{a} = h_{c}',
  None,
  'Dynkin_index',
  'The Dynkin index of the adjoint representation is the dual Coxeter number, in the normalisation '
  'where the fundamental of SU(N) has index one half. It is what turns a count of roots into a '
  'coefficient of the beta function.'),
 ('One-loop coefficient (any group)',
  'Gauge unification',
  'b = -\\frac{11}{3} T_{a} + \\frac{4}{3} n_{f} T_{f}',
  None,
  'Beta_function_(physics)',
  'The one-loop slope for any compact gauge group with n_f Dirac fermions in a representation of '
  'index T_f: the gauge bosons through the adjoint index, the fermions through theirs. With the '
  "adjoint index N and the fundamental's one half it is the SU(N) coefficient of this library, "
  'which is the sense in which that one specialises this. The pure-gauge case, n_f = 0, is what a '
  'lattice without quarks sees.'),
 ('Spin-model transmutation slope',
  'Lattice spin models',
  'r_{s} = \\frac{2 \\pi}{N - 2}',
  None,
  'Asymptotic_freedom',
  'The rate at which the logarithm of the mass falls with the inverse temperature in the O(N) '
  'model: the exponent of the asymptotic-freedom scaling form, taken alone. It is the one number '
  'in that form the theory fixes. For three components it is two pi, and it diverges at two, where '
  'there is no asymptotic freedom to speak of.'),
 ('Gauge transmutation slope',
  'Lattice gauge theory',
  'r_{g} = \\frac{12 \\pi^{2}}{11 N^{2}}',
  None,
  'Asymptotic_freedom',
  'The same quantity for a pure SU(N) lattice gauge theory, with the inverse coupling written in '
  'the convention of this library, beta = 2N / g^2. It comes from the one-loop coefficient of the '
  'pure gauge theory, eleven N over forty-eight pi squared, which in these units gives the mass '
  'falling as exp(-r beta). For two colours it is about 2.69 and for three about 1.20, against two '
  'pi for the O(3) model: the same lattice step has to move the coupling about two and a third '
  'times as far for SU(2) as for O(3), and five times as far for SU(3).'),
 ('Spin block coupling',
  'Block renormalisation',
  '\\beta_{1} = \\beta_{0} - \\frac{\\ln \\ell}{r_{s}}',
  None,
  'Renormalization_group',
  'Blocking a lattice by a length factor ell makes the mass ell times larger in the new lattice '
  'units, and the scaling form says which coupling does that: it is lower by the logarithm of ell '
  'over the slope. This is the leading-order statement, in which the tiling enters only through '
  'the logarithm of its inflation factor; the finite-coupling shift a simulation measures is '
  'larger, because the form is corrected at the values reached. Whether an aperiodic inflation is '
  'a valid blocking at all is not something this equation says.'),
 ('Gauge block coupling',
  'Block renormalisation',
  '\\beta_{1} = \\beta_{0} - \\frac{\\ln \\ell}{r_{g}}',
  None,
  'Renormalization_group',
  'The same rule for the SU(N) lattice, with its own slope. The factor ell is the only input from '
  'the tiling: for one Spectre inflation step it is about 2.806, whose logarithm is 1.032.')]

READINGS_DATA = {'Mass-energy equivalence': {'E': 'energy', 'm': 'mass', 'c': 'speed of light'},
 "Newton's second law": {'F': 'force', 'm': 'mass', 'a': 'acceleration'},
 "Newton's law of gravitation": {'F': 'force',
                                 'G': 'gravitational constant',
                                 'm_1': 'mass',
                                 'm_2': 'mass',
                                 'r': 'length'},
 "Coulomb's law": {'F': 'force',
                   'pi': 'pi',
                   'epsilon_0': 'permittivity',
                   'q_1': 'charge',
                   'q_2': 'charge',
                   'r': 'length'},
 'Coulomb potential energy': {'U': 'energy',
                              'pi': 'pi',
                              'epsilon_0': 'permittivity',
                              'q_1': 'charge',
                              'q_2': 'charge',
                              'r': 'length'},
 'Kinetic energy': {'E_k': 'energy', 'm': 'mass', 'v': 'speed'},
 'Momentum': {'p': 'momentum', 'm': 'mass', 'v': 'speed'},
 'Work done by a force': {'W': 'energy', 'F': 'force', 'd': 'length'},
 'Gravitational potential energy': {'U': 'energy',
                                    'm': 'mass',
                                    'g': 'gravitational field',
                                    'h': 'length'},
 'Surface gravity': {'g': 'gravitational field',
                     'G': 'gravitational constant',
                     'M': 'mass',
                     'r': 'length'},
 "Hooke's law": {'F': 'force', 'k': 'spring constant', 'x': 'length'},
 'Escape velocity': {'v_e': 'speed', 'G': 'gravitational constant', 'M': 'mass', 'r': 'length'},
 'Angular frequency': {'omega': 'angular frequency', 'pi': 'pi', 'f': 'frequency'},
 'Wave speed': {'v': 'speed', 'f': 'frequency', 'lambda': 'wavelength'},
 'Planck relation': {'E': 'energy', 'h': 'planck constant', 'nu': 'frequency'},
 'de Broglie relation': {'lambda': 'wavelength', 'h': 'planck constant', 'p': 'momentum'},
 'Compton wavelength': {'lambda_C': 'wavelength',
                        'h': 'planck constant',
                        'm': 'mass',
                        'c': 'speed of light'},
 'Boltzmann entropy': {'S': 'entropy', 'k_B': 'boltzmann constant', 'W': 'microstate count'},
 'Bekenstein-Hawking entropy': {'S': 'entropy',
                                'k_B': 'boltzmann constant',
                                'c': 'speed of light',
                                'A': 'area',
                                'G': 'gravitational constant',
                                'hbar': 'reduced planck constant'},
 'Thermal energy': {'E': 'energy', 'k_B': 'boltzmann constant', 'T': 'temperature'},
 'Ideal gas law': {'P': 'pressure',
                   'V': 'volume',
                   'n': 'amount',
                   'R': 'gas constant',
                   'T': 'temperature'},
 'Heat capacity': {'Q': 'energy', 'm': 'mass', 'c_p': 'dimensionless', 'T': 'temperature'},
 'Stefan-Boltzmann law': {'j': 'flux', 'sigma': 'stefan-boltzmann constant', 'T': 'temperature'},
 'Density': {'rho': 'density', 'm': 'mass', 'V': 'volume'},
 'Hydrostatic pressure': {'P': 'pressure',
                          'rho': 'density',
                          'g': 'gravitational field',
                          'h': 'length'},
 "Ohm's law": {'V': 'voltage', 'I': 'current', 'R': 'resistance'},
 'Electrical power': {'P': 'power', 'V': 'voltage', 'I': 'current'},
 'Capacitor charge': {'Q': 'charge', 'C': 'capacitance', 'V': 'voltage'},
 'Schwarzschild radius': {'r_s': 'length',
                          'G': 'gravitational constant',
                          'M': 'mass',
                          'c': 'speed of light'},
 'Planck length': {'ell_P': 'length',
                   'hbar': 'reduced planck constant',
                   'G': 'gravitational constant',
                   'c': 'speed of light'},
 "Hubble's law": {'v': 'speed', 'H_0': 'hubble parameter', 'D': 'length'},
 'Debye length': {'lambda_D': 'length',
                  'epsilon_0': 'permittivity',
                  'k_B': 'boltzmann constant',
                  'T': 'temperature',
                  'n': 'number density',
                  'e': 'elementary charge'},
 'Plasma frequency': {'omega_p': 'angular frequency',
                      'n': 'number density',
                      'e': 'elementary charge',
                      'epsilon_0': 'permittivity',
                      'm_e': 'mass'},
 'Pythagorean theorem': {'a': 'length', 'b': 'length', 'c': 'length'},
 'Golden ratio': {'phi': 'golden ratio'},
 'Inflation quadratic': {'lambda': 'inflation factor'},
 'Penrose tile ratio': {'N_f': 'tile count', 'phi': 'golden ratio', 'N_t': 'tile count'},
 'Inflation of area': {'A_prime': 'patch area', 'phi': 'golden ratio', 'A': 'patch area'},
 'Spectre inflation quadratic': {'mu': 'inflation factor'},
 'Spectre inflation factor': {'mu': 'inflation factor'},
 'Spectre reciprocal pair': {'g': 'tile frequency', 'mu': 'inflation factor'},
 'Spectre triplet frequency': {'g': 'tile frequency'},
 'Spectre frequency ratio': {'v_Theta': 'tile frequency',
                             'v_Gamma': 'tile frequency',
                             'mu': 'inflation factor'},
 'Spectre count recurrence': {'N_prime_prime': 'tile count',
                              'N_prime': 'tile count',
                              'N': 'tile count'},
 'Spectre tile area': {'A': 'tile area', 'a': 'length', 'b': 'length'},
 'Hat inflation quadratic': {'kappa': 'inflation factor'},
 'Hat inflation factor': {'kappa': 'inflation factor'},
 'Hat constant as a golden power': {'kappa': 'inflation factor', 'phi': 'golden ratio'},
 'Hat chirality ratio': {'N_u': 'tile count', 'N_r': 'tile count', 'kappa': 'inflation factor'},
 'Hat tile area': {'A_H': 'tile area'},
 'Spectre linear inflation': {'ell': 'linear inflation factor'},
 'Spectre linear-area relation': {'mu': 'inflation factor', 'ell': 'linear inflation factor'},
 'Spectre linear quartic': {'ell': 'linear inflation factor'},
 'Spectre gap-label recurrence': {'G_prime_prime': 'gap label',
                                  'G_prime': 'gap label',
                                  'G': 'gap label'},
 'Spectre growth rate': {'gamma': 'tile-count growth rate', 'mu': 'inflation factor'},
 'Euler characteristic of a closed surface': {'chi': 'euler characteristic', 'g': 'genus'},
 'Riemann-Hurwitz (unbranched)': {'chi_C': 'euler characteristic',
                                  'n': 'covering degree',
                                  'chi_B': 'euler characteristic'},
 'Gauss-Bonnet theorem': {'K': 'gaussian curvature',
                          'k_g': 'geodesic curvature',
                          'pi': 'pi',
                          'chi': 'euler characteristic'},
 'Mean curvature': {'H': 'mean curvature',
                    'kappa_1': 'principal curvature',
                    'kappa_2': 'principal curvature'},
 'Minimal surface condition': {'H': 'mean curvature'},
 'Minimal surface equation': {'u_x': 'slope', 'u_y': 'slope'},
 'Catenoid': {'r': 'length', 'c': 'length', 'z': 'length'},
 'Area functional': {'A': 'area'},
 'Kakeya conjecture': {'n': 'hausdorff dimension'},
 'Perron tree bound': {'A': 'area', 'C': 'dimensionless', 'n': 'count'},
 'Nambu-Goto action': {'S': 'action', 'T': 'string tension'},
 'String tension': {'T': 'string tension', 'pi': 'pi', 'alpha_prime': 'regge slope'},
 'Regge trajectory': {'J': 'spin',
                      'alpha_prime': 'regge slope',
                      'M': 'mass',
                      'alpha_0': 'dimensionless'},
 'String length': {'ell_s': 'length', 'alpha_prime': 'regge slope'},
 'Genus expansion': {'Z': 'partition function',
                     'chi': 'euler characteristic',
                     'g_s': 'string coupling'},
 'Critical dimension': {'D': 'spacetime dimension'},
 'Hagedorn temperature': {'T_H': 'temperature',
                          'pi': 'pi',
                          'alpha_prime': 'regge slope',
                          'k_B': 'boltzmann constant'},
 'Lattice coupling': {'beta': 'inverse coupling', 'N': 'color count', 'g': 'gauge coupling'},
 'Wilson plaquette action': {'S_W': 'lattice action',
                             'beta': 'inverse coupling',
                             'N': 'color count',
                             'U_p': 'plaquette variable'},
 'Yang-Mills action': {'S': 'action', 'g': 'gauge coupling', 'F': 'gauge field strength'},
 'Link gauge transformation': {'U': 'link variable', 'G': 'gauge transformation'},
 'Creutz ratio': {'chi': 'creutz ratio',
                  'W_RT': 'wilson loop',
                  'W_R1T1': 'wilson loop',
                  'W_RT1': 'wilson loop',
                  'W_R1T': 'wilson loop'},
 'Area law': {'W': 'wilson loop', 'sigma': 'confining string tension', 'A': 'loop area'},
 'Creutz string tension': {'chi': 'creutz ratio',
                           'sigma': 'confining string tension',
                           'a': 'lattice spacing'},
 'Lattice effective mass': {'m_eff': 'lattice mass', 'C_t': 'correlator', 'C_t1': 'correlator'},
 'Euclidean correlator decay': {'C_t': 'correlator',
                                'A_0': 'amplitude',
                                'm_L': 'lattice mass',
                                't': 'time slice'},
 'Mass gap in physical units': {'E_gap': 'energy',
                                'hbar': 'reduced planck constant',
                                'c': 'speed of light',
                                'm_L': 'lattice mass',
                                'a': 'lattice spacing'},
 'Integrated autocorrelation time': {'tau': 'autocorrelation time', 'rho_t': 'autocorrelation'},
 'Effective sample size': {'N_eff': 'effective sample size',
                           'N': 'sample count',
                           'tau': 'autocorrelation time'},
 'Sampling efficiency': {'eta': 'sampling efficiency',
                         'N_eff': 'effective sample size',
                         't_wall': 'time'},
 'Efficiency ratio': {'R': 'efficiency ratio',
                      'eta_OR': 'sampling efficiency',
                      'eta_M': 'sampling efficiency'},
 'O(N) lattice model': {'H': 'spin action',
                        'beta': 'spin inverse temperature',
                        'q_x': 'spin variable',
                        'q_y': 'spin variable'},
 'Asymptotic-freedom mass scaling': {'m': 'lattice mass',
                                     'C': 'scaling constant',
                                     'beta': 'spin inverse temperature',
                                     'N': 'component count',
                                     'pi': 'pi'},
 'Second-moment correlation length': {'xi': 'correlation length',
                                      'chi': 'spin susceptibility',
                                      'F': 'structure factor',
                                      'p': 'lattice momentum'},
 'Kosterlitz coupling flow': {'x': 'stiffness deviation', 'y': 'vortex fugacity', 'l': 'rg scale'},
 'Kosterlitz fugacity flow': {'x': 'stiffness deviation', 'y': 'vortex fugacity', 'l': 'rg scale'},
 'Magnon dispersion': {'epsilon': 'magnon energy', 'k_a': 'lattice wave number'},
 'Bloch magnon density': {'n': 'magnon density',
                          'A_0': 'magnon coefficient',
                          'A_1': 'magnon coefficient',
                          'x': 'scaled inverse temperature'},
 'Bilinear-biquadratic chain': {'H': 'spin hamiltonian',
                                'theta': 'mixing angle',
                                'S_i': 'spin operator',
                                'S_i1': 'spin operator'},
 'Spin gap': {'Delta': 'spin gap', 'E_1': 'spin energy level', 'E_0': 'spin energy level'},
 'One-loop running': {'alpha_i': 'fine-structure coupling',
                      'alpha_G': 'fine-structure coupling',
                      'b_i': 'beta function coefficient',
                      'pi': 'pi',
                      'M_G': 'unification scale',
                      'mu': 'renormalisation scale'},
 'Unification prediction': {'alpha_3': 'fine-structure coupling',
                            'alpha_2': 'fine-structure coupling',
                            'alpha_1': 'fine-structure coupling',
                            'b_3': 'beta function coefficient',
                            'b_2': 'beta function coefficient',
                            'b_1': 'beta function coefficient'},
 'Dimensional transmutation': {'Lambda': 'qcd scale',
                               'M': 'unification scale',
                               'b_3': 'beta function coefficient',
                               'alpha': 'fine-structure coupling',
                               'pi': 'pi'},
 'One-loop coefficient (QCD)': {'b': 'beta function coefficient',
                                'N': 'color count',
                                'n_f': 'flavour count'},
 'One-loop coefficient (SUSY)': {'b': 'beta function coefficient',
                                 'N': 'color count',
                                 'n_f': 'flavour count'},
 'Electromagnetic coupling at the Z': {'alpha_em': 'fine-structure coupling',
                                       'alpha_2': 'fine-structure coupling',
                                       'alpha_1': 'fine-structure coupling'},
 'Thomson limit': {'alpha_0': 'fine-structure coupling',
                   'alpha_em': 'fine-structure coupling',
                   'epsilon': 'vacuum polarisation shift'},
 'Generation count': {'N_gen': 'generation count', 'I': 'chiral index', 'G': 'symmetry order'},
 'Heterotic anomaly condition': {'W': 'fivebrane class', 'c_2': 'chern class'},
 'Lie algebra dimension': {'D': 'lie algebra dimension', 'R': 'root count', 'r': 'rank'},
 'Adjoint branching': {'D': 'lie algebra dimension',
                       'D_h': 'subalgebra dimension',
                       'D_V': 'representation dimension'},
 'Product representation dimension': {'D_V': 'representation dimension',
                                      'D_a': 'representation dimension',
                                      'D_b': 'representation dimension'},
 'Orbit size': {'n': 'orbit size', 'w': 'weyl group order', 's': 'stabiliser order'},
 'Triple space': {'T': 'triple count',
                  'n_1': 'orbit size',
                  'n_2': 'orbit size',
                  'n_3': 'orbit size'},
 'Magic square dimension': {'L': 'lie algebra dimension',
                            'tau_A': 'triality algebra dimension',
                            'tau_B': 'triality algebra dimension',
                            'a': 'division algebra dimension',
                            'b': 'division algebra dimension'},
 'Cubic anomaly coefficient': {'A_abc': 'cubic anomaly coefficient'},
 'Anomaly of a product representation': {'A_V': 'cubic anomaly coefficient',
                                         'D_a': 'representation dimension',
                                         'A_b': 'cubic anomaly coefficient'},
 'Charge lattice': {'q': 'integer charge', 't': 'colour triality', 's': 'weak centre parity'},
 'Hypercharge': {'Y': 'hypercharge', 'q': 'integer charge'},
 'Hypercharge normalisation': {'k': 'hypercharge normalisation',
                               'Y': 'hypercharge',
                               'T_3': 'weak isospin'},
 'Anomaly-free hypercharge': {'Y_i': 'hypercharge'},
 'Net chiral class': {'chi_V': 'net chiral class',
                      'N_gen': 'generation count',
                      'chi_F': 'net chiral class'},
 'Special unitary dimension': {'D': 'lie algebra dimension', 'N': 'color count'},
 'Dual Coxeter number (simply laced)': {'h_c': 'dual coxeter number',
                                        'R': 'root count',
                                        'r': 'rank'},
 'Adjoint index': {'T_a': 'adjoint index', 'h_c': 'dual coxeter number'},
 'One-loop coefficient (any group)': {'b': 'beta function coefficient',
                                      'T_a': 'adjoint index',
                                      'n_f': 'flavour count',
                                      'T_f': 'representation index'},
 'Spin-model transmutation slope': {'r_s': 'spin transmutation slope',
                                    'N': 'component count',
                                    'pi': 'pi'},
 'Gauge transmutation slope': {'r_g': 'gauge transmutation slope', 'N': 'color count', 'pi': 'pi'},
 'Spin block coupling': {'beta_1': 'blocked spin inverse temperature',
                         'beta_0': 'spin inverse temperature',
                         'ell': 'linear inflation factor',
                         'r_s': 'spin transmutation slope'},
 'Gauge block coupling': {'beta_1': 'blocked inverse coupling',
                          'beta_0': 'inverse coupling',
                          'ell': 'linear inflation factor',
                          'r_g': 'gauge transmutation slope'}}

LINKS_DATA = [('Schrodinger equation (time-independent)',
  'SPECIALISES',
  'Schrodinger equation (time-dependent)',
  'separating the time dependence of a stationary state leaves an eigenvalue problem in space '
  'alone'),
 ('Wave equation (box form)',
  'SPECIALISES',
  'Klein-Gordon equation',
  'the massless case: with m = 0 the Klein-Gordon mass term vanishes and nothing but the '
  "d'Alembertian is left"),
 ("d'Alembert operator",
  'DEFINES',
  'Wave equation (box form)',
  'the box is what makes the one-symbol form of the wave equation mean anything -- it is an '
  'abbreviation, not a new idea'),
 ("Laplace's equation",
  'SPECIALISES',
  "Poisson's equation",
  'the source-free case: Laplace is Poisson with the density set to zero'),
 ('Helmholtz equation',
  'SPECIALISES',
  'Wave equation',
  'separating out a single frequency turns the wave equation into a condition on the spatial part '
  'alone'),
 ('Wave equation',
  'FIELD',
  'Heat / diffusion equation',
  'the same Laplacian, set equal to a second time derivative rather than a first -- which is the '
  'whole difference between propagating and diffusing'),
 ('Boltzmann entropy',
  'DEFINES',
  'Shannon entropy',
  'counting microstates and counting messages are the same sum, differing only in the base of the '
  'logarithm and a constant'),
 ('Bekenstein-Hawking entropy',
  'SPECIALISES',
  'Boltzmann entropy',
  'a black hole obeys the same entropy law as a gas, but the microstates being counted scale with '
  'horizon area rather than with volume'),
 ('Schwarzschild radius',
  'SPECIALISES',
  'Einstein field equations',
  'the radius falls out of the spherically symmetric vacuum solution'),
 ('Planck relation',
  'DEFINES',
  'de Broglie relation',
  'the same Planck constant, read the other way round: if a wave carries energy in quanta then a '
  'particle carries a wavelength'),
 ('Continuity equation',
  'SPECIALISES',
  'Ampere-Maxwell law',
  'taking the divergence of the Ampere-Maxwell law forces charge conservation -- which is why the '
  'displacement current had to be there'),
 ('Lorentz factor',
  'DEFINES',
  'Mass-energy equivalence',
  'gamma is what the rest energy gets multiplied by once the body moves'),
 ('Surface gravity',
  'SPECIALISES',
  "Newton's law of gravitation",
  'dividing out the falling mass leaves an acceleration that does not depend on what is falling, '
  'which is the equivalence principle in its oldest form'),
 ('Gravitational potential energy',
  'LIMIT_OF',
  'Surface gravity',
  'm g h is the work done against a gravity treated as constant, which it is to excellent accuracy '
  'over any height small against the radius'),
 ('Compton wavelength',
  'DEFINES',
  'Planck relation',
  "the wavelength at which a photon carries exactly a particle's rest energy, which is where "
  'relativity and quantum mechanics stop being separable'),
 ('Momentum',
  'SPECIALISES',
  "Newton's second law",
  'force is the rate of change of momentum, and m a is that rate when the mass is what stays '
  'constant'),
 ('Work done by a force',
  'DEFINES',
  'Kinetic energy',
  'integrating force over distance is what produces the half and the square -- the energy laws are '
  'the force laws, integrated'),
 ("Hooke's law",
  'SPECIALISES',
  'Work done by a force',
  'a restoring force linear in displacement stores energy quadratic in it'),
 ('Escape velocity',
  'SPECIALISES',
  'Schwarzschild radius',
  'setting the escape velocity to c gives the Schwarzschild radius exactly -- a coincidence of the '
  'Newtonian derivation rather than a proof of the relativistic one'),
 ('Planck length',
  'DEFINES',
  'Bekenstein-Hawking entropy',
  'the horizon area measured in Planck areas is the microstate count, which is what makes the '
  'entropy a number rather than a scale'),
 ('Thermal energy',
  'DEFINES',
  'Ideal gas law',
  'the gas law is thermal energy counted per mole rather than per particle, which is the whole '
  'difference between R and k_B'),
 ('Electrical power',
  'SPECIALISES',
  "Ohm's law",
  'substituting one into the other gives the two forms every electronics text reaches for next, P '
  '= I^2 R and P = V^2 / R'),
 ('Hydrostatic pressure',
  'SPECIALISES',
  'Density',
  'the pressure at a depth is the weight of the column above it, so the density is doing the same '
  'job it does in its own definition'),
 ('Angular frequency',
  'DEFINES',
  'Planck relation',
  'writing the Planck relation with omega rather than nu is what turns h into hbar; the two pi has '
  'to go somewhere'),
 ('Wave speed',
  'DEFINES',
  'de Broglie relation',
  'a wavelength is only meaningful for something that travels, and this is the relation that says '
  'how fast'),
 ("Hubble's law",
  'FIELD',
  'Friedmann equation',
  "the Friedmann equation is what H obeys; Hubble's law is what H means"),
 ('Golden ratio',
  'SPECIALISES',
  'Inflation quadratic',
  'the surd is the solution of the quadratic; the quadratic is what the number is for'),
 ('Inflation of area',
  'DEFINES',
  'Golden ratio',
  'inflation scales lengths by phi, so it scales areas by phi squared -- the same constant, one '
  'dimension up'),
 ('Spectre inflation factor',
  'SPECIALISES',
  'Spectre inflation quadratic',
  'the surd is the larger root of the quadratic; the quadratic is what the number is for'),
 ('Spectre triplet frequency',
  'SPECIALISES',
  'Spectre inflation quadratic',
  'and this is the smaller root of the same quadratic, which is why the two are reciprocal without '
  'anyone arranging it'),
 ('Spectre reciprocal pair',
  'DEFINES',
  'Spectre inflation quadratic',
  'the constant term of the quadratic is one, so the product of its roots is one -- the reciprocal '
  'relation is read straight off the polynomial'),
 ('Spectre count recurrence',
  'DEFINES',
  'Spectre inflation quadratic',
  'the recurrence has the inflation quadratic as its characteristic polynomial, which is why '
  'counts grow like powers of the larger root while staying whole numbers'),
 ('Spectre inflation quadratic',
  'FIELD',
  'Inflation quadratic',
  'the same statement for the one-tile system that the golden ratio makes for the two-tile one: an '
  'inflation factor is a quadratic irrational either way, and the irrationality is the whole '
  'obstruction to a period'),
 ('Spectre linear-area relation',
  'DEFINES',
  'Spectre linear inflation',
  'areas scale as the square of lengths, so the area factor is the square of this one'),
 ('Spectre linear quartic',
  'DEFINES',
  'Spectre linear-area relation',
  'substitute the square of the unknown into the inflation quadratic and the quartic is what you '
  'have written'),
 ('Spectre gap-label recurrence',
  'DEFINES',
  'Spectre inflation quadratic',
  'the same characteristic polynomial as the tile counts, applied to powers of the reciprocal root '
  'instead of the root'),
 ('Spectre growth rate',
  'SPECIALISES',
  'Spectre inflation factor',
  'the logarithm of the inflation factor, and nothing more than that'),
 ('Asymptotic-freedom mass scaling',
  'FIELD',
  'Dimensional transmutation',
  'both are dimensional transmutation: a scale that is an exponential of two pi times an inverse '
  'coupling, divided by a positive one-loop coefficient -- N minus two for the spin model, the '
  'magnitude of b for the gauge theory. The shape is shared; the coefficients and the '
  'normalisation of the coupling are not, and nothing here says one model is the other'),
 ('Euclidean correlator decay',
  'FIELD',
  'Asymptotic-freedom mass scaling',
  'both state a lattice mass in lattice units -- the laboratory reads it off a correlator, the '
  'spin model predicts it from the coupling -- and they are about different theories, so the graph '
  'relates them without equating them'),
 ('Second-moment correlation length',
  'FIELD',
  'Euclidean correlator decay',
  'the inverse second-moment length is a mass but not in general the transfer-matrix gap that a '
  'correlator decay measures, so the comparison tests the shape of a scaling form and no more'),
 ('Dimensional transmutation',
  'FIELD',
  'Mass gap in physical units',
  'the lattice spacing in the mass-gap formula is what asymptotic scaling fixes from the bare '
  'coupling: the scale Lambda is the reason a fixed physical mass needs a drifting beta as the '
  'spacing shrinks'),
 ('Hat inflation quadratic',
  'FIELD',
  'Spectre inflation quadratic',
  'two monotile substitutions and two quadratics, x^2 - 7x + 1 and x^2 - 8x + 1: they differ in '
  'one coefficient, and that coefficient is the difference between the quadratic fields generated '
  'by root five and root fifteen'),
 ('Hat tile area',
  'SPECIALISES',
  'Spectre tile area',
  'the hat is the member of the family with a = 1 and b = root three; swapping the two lengths '
  'gives the turtle and a different area, which is the asymmetry the formula records'),
 ('Riemann-Hurwitz (unbranched)',
  'SPECIALISES',
  'Euler characteristic of a closed surface',
  'a cover has n times the Euler characteristic of what it covers, so by chi = 2 - 2g it also has '
  'n times the handles, less n - 1 of them'),
 ('Gauss-Bonnet theorem',
  'DEFINES',
  'Euler characteristic of a closed surface',
  'the integral of curvature over a closed surface is 2 pi chi, which is what makes a topological '
  'count something a ruler can measure'),
 ('Mean curvature',
  'DEFINES',
  'Minimal surface condition',
  'minimal means H = 0, and this is what H is'),
 ('Minimal surface equation',
  'SPECIALISES',
  'Minimal surface condition',
  'the same condition written out for a surface given as a height above the plane, where it '
  'becomes a quasilinear partial differential equation'),
 ('Catenoid',
  'SPECIALISES',
  'Minimal surface equation',
  'the rotationally symmetric solution, and the first minimal surface found that was not the '
  'plane'),
 ('Area functional',
  'DEFINES',
  'Minimal surface equation',
  'setting the first variation of the area integral to zero is where the minimal surface equation '
  'comes from'),
 ('Perron tree bound',
  'SPECIALISES',
  'Kakeya conjecture',
  'the Perron construction settles the measure -- it can be made as small as you like -- and '
  'leaves open the dimension, which is what the conjecture is about'),
 ('Nambu-Goto action',
  'FIELD',
  'Area functional',
  'the same variational problem in a different signature: a soap film extremises area in Euclidean '
  'space, a string extremises area in Minkowski spacetime, and the minus sign under the root is '
  'the only difference the mathematics sees'),
 ('String tension',
  'DEFINES',
  'Nambu-Goto action',
  'the T multiplying the worldsheet area is this tension, and writing it as one over two pi alpha '
  'prime is what introduces the only free parameter the string has'),
 ('String length',
  'DEFINES',
  'String tension',
  'alpha prime is a length squared, so the tension is an inverse area -- in natural units, the '
  'only scale in the theory'),
 ('Euler characteristic of a closed surface',
  'DEFINES',
  'Genus expansion',
  'the coupling is raised to minus the Euler characteristic, so chi = 2 - 2g is what turns a sum '
  'over surfaces into an expansion ordered by loops'),
 ('Hagedorn temperature',
  'FIELD',
  'Boltzmann entropy',
  'counting string states is counting microstates, and they grow exponentially with energy -- fast '
  'enough that above this temperature the sum defining the partition function stops converging'),
 ('Yang-Mills action',
  'LIMIT_OF',
  'Wilson plaquette action',
  'as the lattice spacing goes to zero the sum over squares becomes the integral of the squared '
  'field strength, with the coupling out front -- a statement about a limit that no finite lattice '
  'reaches'),
 ('Lattice coupling',
  'DEFINES',
  'Wilson plaquette action',
  'beta is the only parameter the plaquette sum has, and writing it as 2N over g squared is what '
  'ties the lattice to the continuum coupling'),
 ('Link gauge transformation',
  'DEFINES',
  'Wilson plaquette action',
  'the trace of a closed loop of links is unchanged by a gauge transformation, which is why a sum '
  'of such traces is a legitimate action and a sum of single links would not be'),
 ('Wilson plaquette action',
  'DEFINES',
  'Creutz ratio',
  'the plaquette is the smallest Wilson loop; the same trace over larger rectangles is what a '
  'Creutz ratio is built from'),
 ('Creutz string tension',
  'SPECIALISES',
  'Creutz ratio',
  'for loops that obey the area law exactly the ratio is the string tension in lattice units; for '
  'real data it is an estimate of it'),
 ('Area law',
  'DEFINES',
  'Creutz string tension',
  'the exponent of the area law is the quantity the Creutz ratio is trying to isolate, with the '
  'perimeter terms cancelled'),
 ('String tension',
  'FIELD',
  'Area law',
  'a confining flux tube is modelled as a string, and in that model the tension is the confining '
  'string tension -- but whether the model is exact is a research question, so the two stay '
  'distinct quantities and are not joined'),
 ('Lattice effective mass',
  'DEFINES',
  'Euclidean correlator decay',
  'the logarithm of successive time slices recovers the decay rate, and is flat exactly when a '
  'single state dominates the correlator'),
 ('Euclidean correlator decay',
  'DEFINES',
  'Mass gap in physical units',
  'the decay rate is the mass in lattice units; multiplying by hbar c over the spacing is what '
  'turns it into an energy'),
 ('Mass gap in physical units',
  'SPECIALISES',
  'Mass-energy equivalence',
  'the gap is the rest energy of the lightest state, E = m c^2 with the mass being the one the '
  'correlator measures'),
 ('Mass gap in physical units',
  'FIELD',
  'Compton wavelength',
  'the same mass read as a length: the range over which the lightest state can mediate anything, '
  'which is why a gap means short-range'),
 ('Integrated autocorrelation time',
  'DEFINES',
  'Effective sample size',
  'tau is how many correlated updates it takes to make one independent one, so it is what the '
  'sample count is divided by'),
 ('Effective sample size',
  'DEFINES',
  'Sampling efficiency',
  'the efficiency divides this by the measured runtime, which is the one step that makes it '
  'machine dependent'),
 ('Sampling efficiency',
  'DEFINES',
  'Efficiency ratio',
  "the ratio of two schedules' efficiencies from a common starting configuration; above one means "
  'the extra sweeps paid for themselves'),
 ('O(N) lattice model',
  'DEFINES',
  'Asymptotic-freedom mass scaling',
  'the model is what the scaling form is a statement about: the mass of its lightest excitation, '
  'as the temperature is lowered'),
 ('Asymptotic-freedom mass scaling',
  'DEFINES',
  'Second-moment correlation length',
  'a simulation tests the scaling form through a correlation length, and the second-moment '
  'estimate is the one it can afford'),
 ('Kosterlitz coupling flow',
  'SPECIALISES',
  'Kosterlitz fugacity flow',
  'the two flow equations are one system: the stiffness drives the fugacity and the fugacity '
  'drives the stiffness'),
 ('Magnon dispersion',
  'DEFINES',
  'Bloch magnon density',
  'the density is the thermal occupation of the modes this dispersion labels, and its small-k form '
  'fixes the exponent'),
 ('Bilinear-biquadratic chain',
  'DEFINES',
  'Spin gap',
  'the gap is a property of the spectrum of this operator on a ring'),
 ('Spin gap',
  'FIELD',
  'Mass gap in physical units',
  'a spin gap on a lattice and a mass gap in lattice units are both energies above a ground state, '
  'but one is a finite quantum system and the other a field theory, so they stay distinct '
  'quantities'),
 ('Lie algebra dimension',
  'FIELD',
  'Adjoint branching',
  'both count the dimension of one algebra: as its roots plus its rank, and as the pieces of a '
  'decomposition under a subalgebra; for E8 the first gives 240 + 8 and the second 86 + 162, and '
  'they agree'),
 ('Net chiral class',
  'FIELD',
  'Generation count',
  'both count net chiral generations, by different arguments: here the multiplicity of an SU(3) '
  'acting on the family index, in the heterotic compactification an index divided by the order of '
  "a freely acting group; neither argument predicts the other's number, and agreement would be a "
  'fact about the model and not about the library'),
 ('Hypercharge normalisation',
  'DEFINES',
  'Electromagnetic coupling at the Z',
  'the five thirds in the electromagnetic coupling is this ratio for a complete SU(5) multiplet; '
  'the family traced here gives the same number, which fixes the normalisation and says nothing '
  'about the couplings themselves'),
 ('Anomaly of a product representation',
  'SPECIALISES',
  'Cubic anomaly coefficient',
  'the cubic anomaly coefficient for a tensor product, where the factor acting trivially '
  'contributes only its dimension'),
 ('Special unitary dimension',
  'FIELD',
  'Lie algebra dimension',
  'the same count for the one family of algebras: N squared less one is roots plus rank for SU(N), '
  'where the roots number N squared less N and the rank is N less one'),
 ('Dual Coxeter number (simply laced)',
  'DEFINES',
  'Adjoint index',
  'the adjoint index is this integer, in the normalisation where the fundamental of SU(N) has one '
  'half'),
 ('Adjoint index',
  'DEFINES',
  'One-loop coefficient (any group)',
  'the gauge bosons enter the slope through the adjoint index, so it is what carries the size of '
  'the group into the running'),
 ('One-loop coefficient (QCD)',
  'SPECIALISES',
  'One-loop coefficient (any group)',
  "SU(N) with fundamental quarks: the adjoint index is N and the fundamental's is one half, so "
  'four thirds of n_f times one half is the two thirds of the specialised form'),
 ('Asymptotic-freedom mass scaling',
  'DEFINES',
  'Spin-model transmutation slope',
  'the slope is the exponent of the scaling form with the prefactor removed, the one part the '
  'theory fixes'),
 ('One-loop coefficient (any group)',
  'DEFINES',
  'Gauge transmutation slope',
  'the slope is the pure-gauge case of the coefficient, n_f = 0, written against the lattice '
  'coupling beta = 2N / g^2 and not against 1 / g^2'),
 ('Spin-model transmutation slope',
  'FIELD',
  'Gauge transmutation slope',
  'the same role in two systems -- the exponent in which the mass falls with the inverse coupling '
  '-- with different N: components in the one, colours in the other; they are not equated, and the '
  'O(3) value is two pi where the SU(2) value is about 2.69'),
 ('Spin-model transmutation slope',
  'DEFINES',
  'Spin block coupling',
  'blocking by a length factor divides the mass by it, and the slope says what change of coupling '
  'does that'),
 ('Gauge transmutation slope',
  'DEFINES',
  'Gauge block coupling',
  'the same step for the gauge theory, with its own slope'),
 ('Magic square dimension',
  'FIELD',
  'Lie algebra dimension',
  'two ways to the dimension of one algebra: as roots plus rank, and as the two triality algebras '
  "and three times the product of the two division algebras' dimensions; for the octonions on both "
  "sides the second is 28 + 28 + 3 * 64 and the first 240 + 8, and both are E8's 248"),
 ('Special unitary dimension',
  'FIELD',
  'Net chiral class',
  "the same SU(3) arithmetic in two different groups: the laboratory's colour group and the family "
  'symmetry of the chirality theorem both have dimension eight and dual Coxeter number three, '
  'which is where the three of net chirality comes from; nothing here says they are one group, and '
  'the library joins them by no row'),
 ('Spin block coupling',
  'FIELD',
  'Gauge block coupling',
  'the same rule for two systems: both lower the coupling by the logarithm of the length factor '
  'over a slope, and the factor is the only thing they share; the couplings are an inverse '
  'temperature and an inverse gauge coupling, which the library keeps apart and does not equate, '
  'so the two steps are joined by no row'),
 ('Spectre linear inflation',
  'DEFINES',
  'Spin block coupling',
  'the tiling supplies the length factor and nothing else: at one loop the block step depends on '
  'it only through its logarithm, so any scale factor would do, and what the Spectre adds is a '
  'particular irrational one')]
