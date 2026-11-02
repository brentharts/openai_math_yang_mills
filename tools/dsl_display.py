"""Long displays made from the items of the leandsl DSL, in the manner of sympy_models/concise.py.

concise.py builds large equations from the knowledge graph: a row is one quantity, every term is a library
equation *solved for that quantity* and underbraced with its name, and each solution is substituted back into its
equation and must leave an identity. A failed check raises, so nothing unverified is written. This does the same for
the items of a category (each name below is a link to the item's own subsection, through its label):

    d = Energy.display("energy-figure", caption="...")
    d.identity_row("energy", [mass_energy, planck_relation, kinetic_energy, thermal_energy])
    d.value_row("inflation factor", [a, b], {"g": 4 - sqrt(15)}, "Spectre")   # the terms are equal numbers
    d.chain_row(clock_tail)                                                  # a checked chain, as one display
    d.link_row(a, b)                                                         # a connection of the DSL, as an arrow

* ``identity_row(quantity, items)``: the items are equations (``sp.Eq``); ``reads=`` says which of each item's letters
  denotes ``quantity``. The rows are chosen by the author: a reading such as "length" is a kind of quantity, not one
  quantity, so rows are never made by grouping readings.
* A row is refused when two of its items are connected by ``kind="field"`` (the DSL's word for "of the same kind, but a
  correspondence rather than an equality"), when a letter that two items both read is read differently, when an
  equation cannot be solved for the quantity in exactly one way, or when a solution does not satisfy its equation.
* ``chain_row(item)`` draws the item's chain of inequalities as one line, each step underbraced with the lemma it uses
  and whether sympy decided it or it was only sampled; a chain that does not reach the stated conclusion is refused.
Run ``python3 tools/dsl_display.py model.py [--out DIR] [--pdf]`` to make a document of a model's categories and figures.
"""
from __future__ import annotations

import pathlib
import sys

import sympy as sp

import leandsl
from leandsl import Item, tex, STATUS_TAG, EDGE_KINDS, CATEGORIES, lean, Real, Nat, require, chain, le, eq, E, Rational, fact, binom


class DisplayError(AssertionError):
    pass


_TEXREL = {"<=": r"\le", "<": "<", "=": "=", ">=": r"\ge", ">": ">"}
_ARROWS = {"uses": r"\longrightarrow", "specialises": r"\longrightarrow", "limit-of": r"\longrightarrow",
           "defines": r"\longrightarrow", "field": r"\longleftrightarrow", "cites": r"\longrightarrow",
           "shares": r"\longleftrightarrow", "mentions": r"\longrightarrow"}


def _plain(s: str) -> str:
    return s.replace("&", r"\&").replace("%", r"\%").replace("_", r"\_").replace("#", r"\#")


def _compact(t: str) -> str:
    return t.replace(r"\frac", r"\tfrac").replace(r"\dfrac", r"\tfrac")


def _positive(expr):
    return expr.xreplace({s: sp.Symbol(s.name, positive=True) for s in expr.free_symbols})


def _name(it: Item) -> str:
    """The item's name under a term: a link to its subsection, in plain words."""
    return rf"\hyperref[{it.latex_label()}]{{{_plain(it.name.replace('_', ' '))}}}"


def _under(expr_tex: str, label: str, status: str = "") -> str:
    lab = rf"\text{{\tiny {label}}}"
    if status:
        lab = rf"\substack{{{lab}\\ \text{{\tiny [{status}]}}}}"
    return rf"\underbrace{{{expr_tex}}}_{{{lab}}}"


def _symbol_for(item: Item, quantity: str):
    hits = [p for p, _ in item.params if item.reads.get(p.name) == quantity]
    if not hits:
        raise DisplayError(f"{item.name} reads no letter as {quantity!r} (reads: {item.reads or 'nothing'})")
    if len(hits) > 1:
        raise DisplayError(f"{item.name}: {len(hits)} letters read as {quantity!r}: {', '.join(h.name for h in hits)}")
    return hits[0]


def solve_for(item: Item, quantity: str):
    """``(symbol, expression)``: the item's equation rearranged for the letter that reads as ``quantity``. Exactly one
    root is required (an even power gives two, and a display built on one of them would state something false), and the
    root is substituted back: the result must be an identity."""
    if not isinstance(item, Item) or not isinstance(item.concl, sp.Eq):
        raise DisplayError(f"{getattr(item, 'name', item)}: an identity row needs an equation (sp.Eq), not a relation")
    sym = _symbol_for(item, quantity)
    try:
        roots = sp.solve(item.concl, sym, dict=False)
    except (NotImplementedError, ValueError) as ex:
        raise DisplayError(f"{item.name}: cannot be solved for {quantity}: {ex}")
    if len(roots) != 1 or roots[0].has(sym):
        raise DisplayError(f"{item.name}: solving for {quantity} gives {len(roots)} root(s), not exactly one")
    expr = sp.simplify(roots[0])
    residual = sp.simplify(_positive((item.concl.lhs - item.concl.rhs).subs(sym, expr)))
    if residual != 0:
        raise DisplayError(f"{item.name}: solving for {quantity} leaves residual {residual}")
    return sym, expr


class Display:
    """Rows, and what they bring together; ``figure()`` is the LaTeX figure, printed with its category."""

    def __init__(self, label: str, category=None, caption: str = ""):
        self.label, self.caption, self.category = label, caption, category
        self.rows: list[str] = []
        self.items: list[Item] = []
        self.quantities: list[str] = []
        self.checks: list[str] = []
        self.arrows = 0
        self.kinds: dict[str, int] = {}                  # rows by kind: identity, value, chain, link
        if category is not None:
            if "__figures__" not in category.__dict__:
                category.__figures__ = []
            category.__figures__.append(self)

    # ---- references to the figure --------------------------------------------------------------
    def latex_label(self) -> str:
        return "fig:" + leandsl._LABEL_BAD.sub("-", self.label)

    def latex_ref(self, text: str | None = None, page: bool = False) -> str:
        return leandsl._ref(self.latex_label(), text or rf"Figure~\ref{{{self.latex_label()}}}", page)

    # ---- the rows -------------------------------------------------------------------------------
    def _row(self, kind: str, tex_: str):
        self.kinds[kind] = self.kinds.get(kind, 0) + 1
        self.rows.append((kind, tex_))

    def _note(self, it: Item):
        if not any(it is x for x in self.items):
            self.items.append(it)

    def _guard(self, quantity: str, items):
        if len(items) < 2:
            raise DisplayError(f"row {quantity}: a row needs at least two items")
        for i, a in enumerate(items):
            for b in items[i + 1:]:
                if a is b:
                    raise DisplayError(f"row {quantity}: {a.name} appears twice")
                for l in a.links + b.links:
                    if l.kind == "field" and {id(l.src), id(l.dst)} == {id(a), id(b)}:
                        raise DisplayError(f"{a.name} and {b.name} are connected as of the same kind (a correspondence), "
                                           f"not as equal: they cannot share the row for {quantity}")
        seen: dict[str, tuple[str, str]] = {}
        for it in items:
            for letter, reading in it.reads.items():
                if letter in seen and seen[letter][1] != reading:
                    raise DisplayError(f"row {quantity}: {letter} is {seen[letter][1]!r} in {seen[letter][0]} "
                                       f"but {reading!r} in {it.name}")
                seen.setdefault(letter, (it.name, reading))

    def identity_row(self, quantity: str, items, tail: str = ""):
        """One quantity, stated by every item: each term is the item's equation solved for it."""
        items = list(items)
        self._guard(quantity, items)
        terms, head = [], None
        for it in items:
            sym, expr = solve_for(it, quantity)
            self.checks.append(f"{it.name}, solved for {quantity}, satisfies the equation identically")
            self._note(it)
            head = head or rf"\overbrace{{{tex(sym)}}}^{{\text{{\tiny {_plain(quantity)}}}}}"
            terms.append(_under(_compact(tex(expr)), _name(it)))
        if len({t for t in terms}) != len(terms):
            raise DisplayError(f"row {quantity} repeats a term")
        self.quantities.append(quantity)
        self._row("identity", f"{head} &= {' = '.join(terms)}{tail}")
        return self

    def value_row(self, quantity: str, items, assign: dict, label: str):
        """A row whose terms are equal *numbers* once the letters not solved for take the values in ``assign``."""
        items = list(items)
        self._guard(quantity, items)
        terms, head, values = [], None, []
        for it in items:
            sym, expr = solve_for(it, quantity)
            self._note(it)
            sub = {s: assign[s.name] for s in expr.free_symbols if s.name in assign}
            missing = sorted({s.name for s in expr.free_symbols} - {s.name for s in sub})
            if missing:
                raise DisplayError(f"{it.name}: no value for {missing}")
            values.append(sp.simplify(sp.radsimp(expr.subs(sub))))
            head = head or rf"\overbrace{{{tex(sym)}}}^{{\text{{\tiny {_plain(quantity)}, {_plain(label)}}}}}"
            terms.append(_under(_compact(tex(expr)), _name(it)))
        first = values[0]
        for it, v in zip(items, values):
            if abs(sp.N(v - first, 40)) > sp.Float("1e-30"):
                raise DisplayError(f"{label}: {it.name} gives {v}, not {first}")
        self.checks.append(f"{label}: all {len(items)} terms equal {first} under "
                           + ", ".join(f"{k} = {v}" for k, v in assign.items()))
        self.quantities.append(quantity)
        self._row("value", f"{head} &= {' = '.join(terms)}")
        return self

    def chain_row(self, item: Item, which: int = 0):
        """The item's chain of inequalities as one line. Each step names the lemma it uses (else its ``tag=``, else its
        number) and says whether it was decided symbolically or sampled."""
        if not item.chains:
            raise DisplayError(f"{item.name} has no chain")
        ch = item.chains[which]
        if not ch.closes:
            raise DisplayError(f"{item.name}: its chain does not reach the stated conclusion ({ch.problem})")
        parts = [rf"\overbrace{{{tex(ch.start)}}}^{{\text{{\tiny {_name(item)}}}}}"]
        for i, (st, res) in enumerate(zip(ch.steps, ch.results), 1):
            if st.using is not None:
                self._note(st.using)
            lab = _name(st.using) if st.using is not None else (_plain(st.tag) if st.tag else f"step {i}")
            parts.append(rf"{_TEXREL[st.rel]} {_under(_compact(tex(st.expr)), lab, STATUS_TAG[res.status])}")
        self._note(item)
        self.checks.append(f"{item.name}: a chain of {len(ch.steps)} steps, {ch.summary()}; its ends are the two sides of the conclusion")
        self._row("chain", f"{parts[0]} & {' '.join(parts[1:])}")
        return self

    def link_row(self, a: Item, b: Item):
        """A connection of the DSL (``a.connect(b, kind=...)``) drawn as an arrow between two statements."""
        link = next((l for l in a.links if l.dst is b), None)
        if link is None:
            raise DisplayError(f"{a.name} is not connected to {b.name}")
        word = EDGE_KINDS[link.kind][0].replace("is a special case of", "special case of")
        side = lambda it: _under(_compact(tex(it.concl)), _name(it))
        self._note(a)
        self._note(b)
        self.arrows += 1
        self.checks.append(f"{a.name} -- {link.kind} -> {b.name}: a connection of the DSL, drawn as a correspondence, not as an equation")
        arrow = _ARROWS.get(link.kind, r"\longrightarrow")
        self._row("link", rf"{{}} & {side(a)} \;\overset{{\text{{\tiny {_plain(word)}}}}}{{{arrow}}}\; {side(b)}")
        return self

    def path_row(self, items, head: str = "", per_line: int = 2):
        """A path of connections drawn as statements joined by arrows: each item ``uses`` the next (``a.connect(b)``,
        the default kind), and that link must exist. Hypotheses are left out. ``head`` is plain text over the first line."""
        items = list(items)
        if len(items) < 2:
            raise DisplayError("a path needs at least two items")
        for a, b in zip(items, items[1:]):
            if not any(l.dst is b and l.kind == "uses" for l in a.links):
                raise DisplayError(f"{a.name} is not connected to {b.name} with 'uses': not a path")
        if len({id(i) for i in items}) != len(items):
            raise DisplayError("a path may not visit an item twice")
        for it in items:
            self._note(it)
        side = lambda it: _under(_compact(tex(it.concl)), _name(it))
        arrow = r"\;\overset{\text{\tiny uses}}{\longrightarrow}\;"
        lines = []
        for k in range(0, len(items), per_line):
            chunk = " ".join(side(it) + (arrow if (k + j + 1) < len(items) else "") for j, it in enumerate(items[k:k + per_line]))
            lines.append(chunk)
        first = rf"\text{{\small {head}}}" if head else r"{}"
        rows = [f"{first} & {lines[0]}"] + [f"{{}} & {l}" for l in lines[1:]]
        self.checks.append(f"{len(items)} items, each linked to the next with 'uses'")
        self._row("path", " \\\\[1.4ex]\n".join(rows))
        return self

    def category_link_row(self, a, b):
        """A connection between two categories (``A.connect(B, kind=...)``) drawn as an arrow between their titles."""
        link = next((l for l in a.__dict__.get("__links__", []) if l.dst is b and l.kind != "inherits"), None)
        if link is None:
            raise DisplayError(f"{a.__name__} is not connected to {b.__name__}")
        word = EDGE_KINDS[link.kind][0]
        arrow = _ARROWS.get(link.kind, r"\longrightarrow")
        side = lambda c: rf"\text{{{leandsl._cat_ref(c)}}}"
        note = rf"\underset{{\text{{\tiny {leandsl._plain_tex(link.note)}}}}}{{\overset{{\text{{\tiny {_plain(word)}}}}}{{{arrow}}}}}" if link.note \
            else rf"\overset{{\text{{\tiny {_plain(word)}}}}}{{{arrow}}}"
        self.arrows += 1
        self.checks.append(f"{a.__name__} -- {link.kind} -> {b.__name__}: a connection recorded between the categories")
        self._row("category", rf"{{}} & {side(a)} \;{note}\; {side(b)}")
        return self

    # ---- the document ---------------------------------------------------------------------------
    def tex(self) -> str:
        if not self.rows:
            raise DisplayError(f"display {self.label}: no rows")
        ids = [t for k_, t in self.rows if k_ not in ("link", "category")]
        arrows = [t for k_, t in self.rows if k_ in ("link", "category")]
        body = " \\\\[1.2ex]\n".join(ids)
        if arrows:
            body += " \\\\[2.4ex]\n" + " \\\\[1.6ex]\n".join(arrows)
        k = self.kinds
        n_id = k.get("identity", 0) + k.get("value", 0)
        lines = []
        if n_id:
            nit = len({id(i) for i in self.items})
            lines.append(rf"\text{{identities:}} \\ \text{{{n_id} quantit{'ies' if n_id != 1 else 'y'}}}")
        if k.get("chain"):
            lines.append(rf"\text{{chains:}} \\ \text{{{k['chain']} checked}}")
        if k.get("path"):
            lines.append(rf"\text{{paths:}} \\ \text{{{k['path']}}}")
        if self.arrows:
            lines.append(rf"\text{{connections:}} \\ \text{{{self.arrows} link{'s' if self.arrows != 1 else ''}}}")
        lines.append(rf"\text{{items:}} \\ \text{{{len(self.items)} brought together}}")
        tally = r"\begin{array}{l} " + r" \\[0.6ex] ".join(lines) + r" \end{array}"
        return f"\\left.\n\\begin{{aligned}}\n{body}\n\\end{{aligned}}\n\\;\\right\\}} {tally}"

    def figure(self, width: float = 0.98) -> str:
        k, sent = self.kinds, []
        if k.get("identity") or k.get("value"):
            sent.append("Each term of an identity row is an item of the DSL solved for the quantity named over the row and "
                        "underbraced with the item's name (a link); the solution was substituted back into the item's equation "
                        "and left an identity. The letters are those of the items, so an equals sign says that each item states "
                        "the same quantity, not that its other letters take one common value.")
        if k.get("chain"):
            sent.append("A chain is drawn only when its ends are the two sides of the item's statement; each step is "
                        "underbraced with the lemma it uses (a link) or its tag, and [symbolic] says that sympy decided the "
                        "step while [sampled] says it held at random points satisfying the hypotheses, which is evidence, not a proof.")
        if k.get("link"):
            sent.append("An arrow is a connection recorded between two items, not an equation.")
        if k.get("path"):
            sent.append("In a path each statement is followed by one that its proof mentions by name: a textual reference, "
                        "not a checked dependence.")
        if k.get("category"):
            sent.append("An arrow between two categories is a connection recorded between them, not a statement of mathematics.")
        sent.append("Hypotheses are left out of the displays; they are in the items.")
        cap = (self.caption + " " if self.caption else "") + " ".join(sent)
        return (f"\\begin{{figure}}[p]\n\\centering\n\\resizebox{{{width}\\textwidth}}{{!}}{{$\\displaystyle\n{self.tex()}\n$}}\n"
                f"\\caption{{{cap}}}\n\\label{{{self.latex_label()}}}\n\\end{{figure}}\n")


# ---- a document of a model's categories ------------------------------------------------------------
_PREAMBLE = r"""\documentclass[10pt,a4paper]{article}
\usepackage[a4paper,margin=20mm]{geometry}
\usepackage{amsmath,amssymb,graphicx}
\usepackage[colorlinks=true,linkcolor=blue!50!black]{hyperref}
\usepackage{xcolor}
\setlength{\parindent}{0pt}\setlength{\parskip}{4pt}
\setlength{\emergencystretch}{3em}
\providecommand{\lt}{<}\providecommand{\gt}{>}
"""


def document(cats, title: str = "DSL model") -> str:
    body = "\n".join(leandsl.category_tex(c) for c in cats)
    doc = (_PREAMBLE + leandsl.TEX_PREAMBLE + f"\n\\title{{{title}}}\\author{{}}\\date{{}}\n\\begin{{document}}\n\\maketitle\n"
           + body + "\n\\end{document}\n")
    return leandsl.close_refs(doc)[0]


def main(argv=None):
    import argparse, runpy, subprocess
    ap = argparse.ArgumentParser(description="A LaTeX document of the categories (and long displays) of model files")
    ap.add_argument("files", nargs="*")
    ap.add_argument("--out", default="build/display")
    ap.add_argument("--pdf", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    here = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(here))
    sys.modules.setdefault("dsl_display", sys.modules[__name__])
    before = list(CATEGORIES)
    for f in a.files:
        runpy.run_path(f, run_name="dsl_model")
    cats = [c for c in CATEGORIES if c not in before]
    figs = sum(len(c.__dict__.get("__figures__", [])) for c in cats)
    doc = document(cats)
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tf = out / "display.tex"
    tf.write_text(doc)
    ck = leandsl.check_labels(doc)
    print(f"{len(cats)} categories, {figs} figure(s) -> {tf}; labels: {ck['defined']} defined, "
          f"{len(ck['dangling'])} dangling, {len(ck['duplicate'])} duplicate")
    if a.pdf:
        r = None
        for _ in range(3):
            r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error", tf.name], cwd=out, capture_output=True, text=True)
            if r.returncode:
                print("\n".join(l for l in r.stdout.splitlines() if l.startswith("!"))[:1500])
                print("pdflatex failed")
                return 1
        und = [l for l in r.stdout.splitlines() if "undefined" in l and "eference" in l]
        print("pdf ok" + (f", {len(und)} undefined reference warning(s)" if und else ", references resolved"))


def selftest():
    n = [0]

    def ok(c, m=""):
        n[0] += 1
        assert c, m

    saved_cats, saved_reg = list(CATEGORIES), list(leandsl.REGISTRY)
    try:
        E_, m, c, h, nu, v, kB, T = sp.symbols("E m c h nu v k_B T")

        @lean.category(r"\section{Energy}", title="Energy")
        class En:
            pass

        @En.theorems(reads={"E": "energy", "m": "mass", "c": "speed of light"})
        def mass_energy(E: Real, m: Real, c: Real):
            return sp.Eq(E, m * c ** 2)

        @En.theorems(reads={"E": "energy", "h": "planck constant", "nu": "frequency"})
        def planck(E: Real, h: Real, nu: Real):
            return sp.Eq(E, h * nu)

        @En.theorems(reads={"E_k": "energy", "m": "mass", "v": "speed"})
        def kinetic(E_k: Real, m: Real, v: Real):
            return sp.Eq(E_k, m * v ** 2 / 2)

        @En.theorems(reads={"x": "speed", "m": "the other mass"})
        def clashing(x: Real, m: Real):
            return sp.Eq(x, m)

        @En.theorems(reads={"x": "speed"})
        def square(x: Real, y: Real):
            return sp.Eq(y, x ** 2)

        @En.theorems(reads={"y": "speed"})
        def ineq(x: Real, y: Real):
            return x <= y

        d = En.display("energy", caption="Energy.")
        d.identity_row("energy", [mass_energy, planck, kinetic])
        t = d.tex()
        ok(t.count("\\underbrace") == 3 and "\\overbrace{E}^{\\text{\\tiny energy}}" in t and "c^{2} m" in t and "h \\nu" in t
           and "\\tfrac{m v^{2}}{2}" in t, t)
        ok("\\hyperref[cat:En:thm:mass_energy]{mass energy}" in t, t)
        ok(len(d.checks) == 3 and [i.name for i in d.items] == ["mass_energy", "planck", "kinetic"])
        fig = d.figure()
        ok("\\label{fig:energy}" in fig and "\\resizebox" in fig and fig.count("\\caption{Energy. Each term") == 1, fig)
        category_text = leandsl.category_tex(En)
        ok(d.latex_ref() == "\\hyperref[fig:energy]{Figure~\\ref{fig:energy}}")
        ok("\\begin{figure}" in category_text and category_text.count("\\label{fig:energy}") == 1)
        for bad, why in ((lambda: Display("a").identity_row("energy", [mass_energy]), "one item"),
                         (lambda: Display("b").identity_row("energy", [mass_energy, mass_energy]), "twice"),
                         (lambda: Display("c").identity_row("energy", [mass_energy, clashing]), "no reading"),
                         (lambda: Display("d").identity_row("mass", [mass_energy, clashing]), "clash: m is mass and the other mass"),
                         (lambda: Display("e").identity_row("speed", [square, clashing]), "x read twice, same; but square has y not read"),
                         (lambda: Display("f").identity_row("speed", [square, ineq]), "not an equation"),
                         (lambda: Display("g").tex(), "no rows")):
            try:
                bad()
                ok(False, why + " accepted")
            except DisplayError:
                ok(True)
        # an even power has two roots: refused, as in the knowledge base
        try:
            Display("h").identity_row("speed", [square, kinetic])
            ok(False, "two roots accepted")
        except DisplayError as ex:
            ok("root" in str(ex) or "reads no letter" in str(ex), str(ex))
        try:
            solve_for(square, "speed")                       # y = x**2 solved for x: two roots
            ok(False, "two roots accepted")
        except DisplayError as ex:
            ok("2 root(s)" in str(ex), str(ex))
        # a correspondence is not an equality
        mass_energy.connect(planck, kind="field", note="same kind")
        try:
            Display("i").identity_row("energy", [mass_energy, planck])
            ok(False, "field-linked items shared a row")
        except DisplayError as ex:
            ok("correspondence" in str(ex), str(ex))
        d2 = Display("j")
        d2.link_row(mass_energy, planck)
        ok("\\longleftrightarrow" in d2.tex() and "of the same kind as" in d2.tex(), d2.tex())
        try:
            d2.link_row(planck, kinetic)
            ok(False, "a link that is not there was drawn")
        except DisplayError:
            ok(True)
        # value rows: equal numbers under an assignment
        r, g, mu = sp.symbols("r g mu")

        @En.theorems(reads={"mu": "inflation", "g": "gap"})
        def a_form(mu: Real, g: Real):
            return sp.Eq(mu, 4 + sp.sqrt(15) + 0 * g)

        @En.theorems(reads={"mu": "inflation", "g": "gap"})
        def b_form(mu: Real, g: Real):
            return sp.Eq(mu, 1 / (4 - sp.sqrt(15)) + 0 * g)

        @En.theorems(reads={"mu": "inflation", "g": "gap"})
        def wrong(mu: Real, g: Real):
            return sp.Eq(mu, 5 + g)
        Display("k").value_row("inflation", [a_form, b_form], {"g": 1}, "both")
        try:
            Display("l").value_row("inflation", [a_form, wrong], {"g": 1}, "wrong")
            ok(False, "unequal numbers accepted")
        except DisplayError as ex:
            ok("gives" in str(ex), str(ex))
        try:
            Display("m").value_row("inflation", [a_form, wrong], {}, "missing")
            ok(False)
        except DisplayError as ex:
            ok("no value" in str(ex), str(ex))
        ok(True)

        # a chain as one display
        @lean.category(r"\section{Tail}", title="Tail")
        class Tail:
            pass

        @Tail.lemmas(reads={"r": "the exponent"})
        def fact_lower(r: Nat):
            return sp.Pow(sp.Mul(r, sp.Pow(E, -1, evaluate=False), evaluate=False), r, evaluate=False) <= fact(r)

        @Tail.theorems(reads={"r": "the exponent"})
        def tail(N: Nat, r: Nat, B: Real, delta: Real):
            require(B >= 0, delta >= 0, r > 0, 4 * E * B * (N + r) * delta <= r)
            return binom(N + r, r) * (2 * B * delta) ** r <= Rational(1, 2) ** r

        @tail.chain
        def tail_chain(N, r, B, delta):
            nn = N + r
            return chain(binom(nn, r) * (2 * B * delta) ** r,
                         le(nn ** r / fact(r) * (2 * B * delta) ** r, "choose", tag="choose bound"),
                         eq((2 * B * delta * nn) ** r / fact(r), "ring"),
                         le((2 * B * delta * nn) ** r / (r / E) ** r, "stirling", using=fact_lower),
                         eq((2 * B * delta * nn / (r / E)) ** r, "div_pow"),
                         le(Rational(1, 2) ** r, "monotone"))
        d3 = Tail.display("tail")
        d3.chain_row(tail)
        t3 = d3.tex()
        ok(t3.count("\\underbrace") == 5 and "choose bound" in t3 and "\\hyperref[cat:Tail:lem:fact_lower]{fact lower}" in t3
           and "sampled" in t3 and "symbolic" in t3 and "step 2" in t3, t3)
        ok([i.name for i in d3.items] == ["fact_lower", "tail"])
        fc = d3.figure()
        ok("A chain is drawn only when" in fc and "identity" not in fc.split("\\caption")[1] and "An arrow" not in fc, fc)
        ok("\\text{chains:}" in d3.tex() and "identities" not in d3.tex(), d3.tex())
        ok("\\text{identities:}" in d.tex() and "chains" not in d.tex() and "1 quantity" in d.tex(), d.tex())
        # a chain that does not reach the conclusion is not drawn

        @Tail.theorems
        def tail2(x: Real, y: Real):
            require(x >= 0, y >= 0)
            return x ** 2 + y ** 2 >= 3 * x * y

        @tail2.chain
        def tail2_chain(x, y):
            return chain(x ** 2 + y ** 2, le(x ** 2 + y ** 2 + 1, "add one"))
        try:
            Display("n").chain_row(tail2)
            ok(False, "open chain drawn")
        except DisplayError as ex:
            ok("does not reach" in str(ex), str(ex))
        try:
            Display("o").chain_row(mass_energy)
            ok(False)
        except DisplayError:
            ok(True)
        # a link row between items of a category
        fact_lower.connect(tail, kind="specialises")
        d4 = Display("p")
        d4.link_row(fact_lower, tail)
        ok(d4.tex().count("\\underbrace") == 2 and "\\overset{\\text{\\tiny special case of}}" in d4.tex(), d4.tex())
        # paths and category links
        mass_energy.links[:] = []
        planck.links[:] = []
        mass_energy.connect(planck)
        planck.connect(kinetic)
        dp = Display("q")
        dp.path_row([mass_energy, planck, kinetic], head="a path")
        tp = dp.tex()
        ok(tp.count("\\underbrace") == 3 and tp.count("uses}}{\\longrightarrow}") == 2 and "\\\\[1.4ex]" in tp and "\\text{paths:}" in tp, tp)
        ok("a textual reference" in dp.figure() and "category" not in dp.figure().split("\\caption")[1], dp.figure())
        for bad in (lambda: Display("r").path_row([mass_energy]), lambda: Display("s").path_row([mass_energy, kinetic]),
                    lambda: Display("t").path_row([mass_energy, planck, mass_energy])):
            try:
                bad()
                ok(False, "bad path accepted")
            except DisplayError:
                ok(True)
        En.connect(Tail, kind="shares", note="3 shared Lean modules")
        dc = Display("u")
        dc.category_link_row(En, Tail)
        tc = dc.tex()
        ok("\\longleftrightarrow" in tc and "3 shared Lean modules" in tc and "\\hyperref[cat:En]{Energy}" in tc, tc)
        ok("between them" in dc.figure() or "between two categories" in dc.figure())
        try:
            dc.category_link_row(Tail, En)
            ok(False)
        except DisplayError:
            ok(True)
        # the document
        doc = document([En, Tail])
        ck = leandsl.check_labels(doc)
        ok(ck["dangling"] == [] and ck["duplicate"] == [] and "\\begin{figure}" in doc and "\\providecommand{\\dsllabel}" in doc, ck)
    finally:
        CATEGORIES[:] = saved_cats
        leandsl.REGISTRY[:] = saved_reg
    print(f"dsl_display: {n[0]} checks ok")


if __name__ == "__main__":
    sys.modules.setdefault("leandsl", sys.modules.get("leandsl") or __import__("leandsl"))
    main()
