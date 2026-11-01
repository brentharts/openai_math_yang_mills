#!/usr/bin/env python3
r"""sympy -> LaTeX -> Lean 4, with every step checked against the model.

A model in ``sympy_models/`` states an equation as a sympy relation. This tool

  1. prints it as LaTeX (``to_tex``, sympy's printer with multi-letter symbols made unambiguous),
     or takes the equation's hand-written ``latex_str`` when the model supplies one;
  2. reads that LaTeX back with ``texmath.parse`` and turns the tree into sympy again;
  3. compares the result with the model's own expression (``same``: symbolically, else numerically).
     This is the step that makes the rest trustworthy: the Lean is generated from the *LaTeX the
     reader sees*, and only after that LaTeX has been shown to say what the sympy model said;
  4. prints the tree as a Lean 4 theorem (``lean_emit``), with the symbols' sympy assumptions
     as hypotheses and the LaTeX, source and provenance in the doc comment;
  5. optionally gives the file to Lean (``--check``), and reports which theorems Lean accepted
     as statements and which it proved. Without a Lean on the machine nothing is claimed.

Every equation lands in exactly one status, and the report counts them:

    verified   the LaTeX reads back to the model's expression; a Lean theorem was written
    tex-only   hand-written LaTeX with no sympy relation to compare to; written only with --tex-only
    mismatch   the LaTeX does NOT read back to the model's expression (a documentation bug)
    refused    the LaTeX or the expression is outside the subset; the reason is recorded

What a theorem says is a statement, not a result: its body is ``sorry`` (Lean calls that out) unless
``--prove`` found a rational-function identity for ``ring1``/``field_simp`` to close. Most model
equations are laws or definitions, not identities, and stay conjectures. A conjecture that Lean
accepts is a statement that parses, elaborates and is a Prop; nothing more.

    python3 tools/sym2lean.py                         default models, both flavors
    python3 tools/sym2lean.py bloch_law kt_flow       chosen models
    python3 tools/sym2lean.py --list                  census of the registry, nothing written
    python3 tools/sym2lean.py --latex 'E = m c^{2}'   one hand-written statement
    python3 tools/sym2lean.py --check --prove         ask Lean; try ring/field_simp on identities
    python3 tools/sym2lean.py --selftest
"""

from __future__ import annotations

import argparse
import importlib
import pkgutil
import random
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import sympy as sp
from sympy.core.relational import Relational

import common as C
import lean_emit as E
import texmath as T

SYMMODELS = C.ROOT / "sympy_models"
LEAN_PROJECT = C.ROOT / "lean"

if str(SYMMODELS) not in sys.path:          # appended, so tools/ names (common, ...) are never shadowed
    sys.path.append(str(SYMMODELS))
import document as _document  # noqa: E402  (the fork's LaTeX escaper, so text looks as it does there)
DEFAULT_MODELS = ("bloch_law", "kt_flow", "on_scaling", "ym_asymptotic_freedom", "haldane_gap")
STATUSES = ("verified", "tex-only", "mismatch", "refused")


# ----------------------------------------------------------------------------
# sympy -> LaTeX
# ----------------------------------------------------------------------------
@dataclass
class SymInfo:
    r"""One sympy symbol, and how it travels through LaTeX.

    The fork names sympy symbols with LaTeX source (``Symbol('\beta_k')``, ``Symbol('a\Lambda_L')``),
    so a name is not an identifier. ``canon`` is the name the reader and Lean use; ``display`` is the
    override handed to sympy's printer (None: sympy's own printing reads back fine); ``atom`` is set
    for a *composite* name, one quantity that LaTeX shows as juxtaposed letters (``a\Lambda_L``):
    its text is replaced by ``\mathrm{canon}`` before reading, so it is read as the one quantity it
    is and not as a product.
    """
    sym: sp.Symbol
    canon: str
    display: str | None = None
    atom: str | None = None


def _sanitize(name: str) -> str:
    s = re.sub(r"[^A-Za-z0-9]", "", name)
    return s if s and not s[0].isdigit() else "q" + s


def _lone_name(name: str) -> str | None:
    """The canonical name if ``name`` is a single quantity sympy and texmath both print and read the same way."""
    base, _, rest = name.partition("_")
    if rest and not rest.replace(",", "").isalnum():
        return None
    if len(base) == 1 and base.isalpha() and base.isascii():
        return name
    if base in T.GREEK and base != "pi":
        return name
    if base.startswith("\\"):
        try:
            t = T.parse(name)
        except T.Unreadable:
            return None
        return t.name if isinstance(t, T.Sym) else None
    return None


_DISGUISED = re.compile(r"\\(?:lim|sum|prod|int|oint|partial|nabla|infty)(?![A-Za-z])|\\[dt]?frac\{\s*(?:d|\\partial)\b")


def analyse(symbols) -> list[SymInfo]:
    infos: list[SymInfo] = []
    for s in sorted(symbols, key=lambda x: x.name):
        name = s.name
        base, _, rest = name.partition("_")
        if name == "e":
            infos.append(SymInfo(s, "e", r"\mathit{e}"))            # sympy prints Euler's constant as e too
        elif (canon := _lone_name(name)) is not None:
            infos.append(SymInfo(s, canon))
        elif base.isalnum() and base.isascii() and len(base) > 1 and (not rest or rest.replace(",", "").isalnum()):
            sub = ("_{%s}" % rest) if rest else ""
            infos.append(SymInfo(s, name, r"\mathrm{%s}%s" % (base, sub)))   # Nsites, not N s i t e s
        else:
            if _DISGUISED.search(name):
                raise T.Unreadable(
                    f"the symbol named {name!r} is an operator (limit, sum, derivative ...) written as a "
                    "symbol: the model does not define it, so a Lean statement about it would say nothing")
            infos.append(SymInfo(s, _sanitize(name), "{%s}" % name, name))   # composite: one quantity
    seen: dict[str, str] = {}
    for i in infos:
        if i.canon in seen:
            raise T.Unreadable(f"symbols {seen[i.canon]!r} and {i.sym.name!r} both read as {i.canon!r}")
        seen[i.canon] = i.sym.name
    return infos


def symbol_names(infos: list[SymInfo]) -> dict:
    return {i.sym: i.display for i in infos if i.display}


def to_tex(expr, infos: list[SymInfo] | None = None) -> str:
    infos = analyse(expr.free_symbols) if infos is None else infos
    return sp.latex(expr, symbol_names=symbol_names(infos))


def apply_atoms(tex: str, infos: list[SymInfo]) -> str:
    r"""Replace each composite symbol's text by ``\mathrm{canon}``, longest first."""
    repl = []
    for i in infos:
        if i.atom:
            repl += [("{%s}" % i.atom, i.canon), (i.atom, i.canon)]
    for text, canon in sorted(repl, key=lambda r: -len(r[0])):
        tex = tex.replace(text, r"\mathrm{%s}" % canon)
    return tex


_ELEMENTARY = (sp.sin, sp.cos, sp.tan, sp.sinh, sp.cosh, sp.tanh, sp.asin, sp.acos, sp.atan,
               sp.exp, sp.log, sp.Abs)
_NOT_ARITHMETIC = (sp.Derivative, sp.Integral, sp.Sum, sp.Product, sp.Piecewise, sp.Limit, sp.Subs, sp.Order)


def unsupported(expr) -> list[str]:
    """What in a sympy expression the subset cannot state. Known on the sympy side, where a product
    ``zeta (3/2)`` and the function ``zeta(3/2)`` are not confusable, unlike in their LaTeX."""
    bad = set()
    for node in sp.preorder_traversal(expr):
        if isinstance(node, _NOT_ARITHMETIC):
            bad.add(type(node).__name__)
        elif isinstance(node, sp.Function) and not isinstance(node, _ELEMENTARY):
            bad.add(f"{node.func}(...)")
        elif node in (sp.I, sp.oo, -sp.oo, sp.zoo, sp.nan):
            bad.add(str(node))
    return sorted(bad)


# ----------------------------------------------------------------------------
# Is the LaTeX still the model?
# ----------------------------------------------------------------------------
def same(a, b) -> str | None:
    """``'symbolic'`` if a - b simplifies to 0, ``'numeric'`` if it vanishes at random points, else None."""
    d = a - b
    if d == 0:
        return "symbolic"
    if sp.count_ops(d) < 250:
        try:
            if sp.simplify(d) == 0:
                return "symbolic"
        except Exception:
            pass
    rng = random.Random(20261101)
    good = 0
    for _ in range(12):
        sub = {}
        for s in d.free_symbols:
            if s.is_integer:
                sub[s] = rng.randint(1, 6) if (s.is_nonnegative or s.is_positive) else rng.randint(-6, 6)
            elif s.is_positive or s.is_nonnegative:
                sub[s] = rng.uniform(0.3, 2.5)
            elif s.is_negative:
                sub[s] = -rng.uniform(0.3, 2.5)
            else:
                sub[s] = rng.uniform(-2.5, 2.5)
        try:
            v = complex(d.subs(sub).evalf(30))
            scale = 1 + sum(abs(complex(x.subs(sub).evalf(30))) for x in (a, b))
        except Exception:
            continue
        if v != v or scale != scale:        # nan
            continue
        if abs(v) > 1e-9 * scale:
            return None
        good += 1
    return "numeric" if good >= 4 else None


def compare(back, expr) -> tuple[str | None, str]:
    """Is the relation read back from LaTeX the model's relation? (method, note)"""
    if type(back) is not type(expr):
        return None, f"the LaTeX states a {type(back).__name__}, the model a {type(expr).__name__}"
    m1, m2 = same(back.lhs, expr.lhs), same(back.rhs, expr.rhs)
    if m1 and m2:
        return ("numeric" if "numeric" in (m1, m2) else "symbolic"), ""
    if isinstance(expr, sp.Equality):                # a = b read as b = a, or terms moved across
        m = same(back.lhs - back.rhs, expr.lhs - expr.rhs) or same(back.lhs - back.rhs, expr.rhs - expr.lhs)
        if m:
            return m, "sides rearranged relative to the model (lhs - rhs agrees)"
    return None, "lhs/rhs differ from the model's"


# ----------------------------------------------------------------------------
# Binders: what a symbol is, from what sympy was told about it
# ----------------------------------------------------------------------------
def binder_for(info: SymInfo) -> E.Binder:
    sym = info.sym
    if sym.is_real is False or sym.is_imaginary or sym.is_complex is True and sym.is_real is not True:
        raise T.Unreadable(f"symbol {sym.name} is not real (sympy says so); the statements are over ℝ")
    if sym.is_integer:
        kind = "nat" if sym.is_nonnegative else "int"
    else:
        kind = "real"
    assume = []
    if sym.is_positive:
        assume.append("positive")
    elif sym.is_negative:
        assume.append("negative")
    else:
        if sym.is_nonnegative and kind != "nat":
            assume.append("nonnegative")
        if sym.is_nonpositive:
            assume.append("nonpositive")
        if sym.is_nonzero:
            assume.append("nonzero")
    return E.Binder(info.canon, kind, tuple(assume))


def is_field_identity(rel) -> bool:
    """A true identity between rational functions: the one case ring1 / field_simp can close."""
    if not isinstance(rel, sp.Equality):
        return False
    d = rel.lhs - rel.rhs
    for node in sp.preorder_traversal(d):
        if isinstance(node, sp.Pow) and not node.exp.is_Integer:
            return False
        if not isinstance(node, (sp.Add, sp.Mul, sp.Pow, sp.Symbol, sp.Number)):
            return False
    return sp.cancel(sp.together(d)) == 0


# ----------------------------------------------------------------------------
# One equation
# ----------------------------------------------------------------------------
@dataclass
class Row:
    key: str
    title: str = ""
    entry: str = ""
    source_kind: str = ""
    source: str = ""
    notes: str = ""
    latex: str = ""
    status: str = "refused"
    method: str = ""                  # how sympy and the LaTeX were compared
    reason: str = ""                  # why refused / mismatched / any caveat
    binder_tex: list = field(default_factory=list)       # (tex, kind, assume)
    theorems: dict = field(default_factory=dict)         # flavor -> lean_emit.Theorem
    proof_attempted: bool = False
    lean: dict = field(default_factory=dict)             # flavor -> 'proved'|'statement'|'error'|'unchecked'

    def summary(self) -> dict:
        return {"key": self.key, "title": self.title, "entry": self.entry, "status": self.status,
                "method": self.method, "reason": self.reason, "latex": self.latex,
                "theorem": next((t.name for t in self.theorems.values()), None),
                "proof_attempted": self.proof_attempted, "lean": self.lean}


def _squash(text: str) -> str:
    return re.sub(r"\s+", " ", text).strip()


def doc_comment(row: Row, rel_tex: str) -> str:
    head = f"**{row.title or row.key}**" + (f" (entry {row.entry})" if row.entry else "")
    lines = [head]
    if row.source:
        lines.append(f"Source ({row.source_kind}): {_squash(row.source)}")
    lines.append(f"LaTeX: `{_squash(rel_tex)}`")
    lines.append(f"Translated from the model by tools/sym2lean.py; LaTeX vs sympy: {row.method or 'not compared'}.")
    if row.reason:
        lines.append(f"Note: {row.reason}.")
    if row.notes:
        lines.append(f"Caveat from the model: {_squash(row.notes)}")
    return "\n".join(lines)


def process(eq, flavors=("mathlib", "core"), prove: bool = False, tex_only: bool = False) -> Row:
    """One ``eqlib.Equation`` (or anything with the same fields) through the whole pipeline."""
    row = Row(eq.key, getattr(eq, "title", ""), getattr(eq, "entry", ""), getattr(eq, "source_kind", ""),
              getattr(eq, "source", ""), getattr(eq, "notes", ""))
    expr = getattr(eq, "expr", None)
    rel = expr if isinstance(expr, Relational) else None
    latex_str = getattr(eq, "latex_str", None)

    # 0. what the sympy side already knows it cannot state
    infos: list[SymInfo] = []
    if rel is not None:
        bad = unsupported(rel)
        if bad:
            row.latex = _squash(latex_str) if latex_str is not None else ""
            row.reason = f"the model uses {', '.join(bad)}, outside the subset (arithmetic and elementary functions)"
            return row
        try:
            infos = analyse(rel.free_symbols)
        except T.Unreadable as exc:
            row.reason = str(exc)
            return row

    # 1. the LaTeX the reader sees
    try:
        if latex_str is not None:
            tex, euler = eq.latex(), False
        elif rel is not None:
            tex, euler = to_tex(rel, infos), True
        elif expr in (sp.true, sp.false):
            row.reason = ("the model's expression evaluated to the Boolean %s (sympy decided it from the "
                          "symbols' assumptions); build it with evaluate=False or give a latex_str" % expr)
            return row
        else:
            row.reason = "no relation to state: the expression is not an equation or inequality and there is no latex_str"
            return row
    except Exception as exc:
        row.reason = f"could not print the equation: {exc}"
        return row
    row.latex = tex

    # 2. read it back
    by_canon = {i.canon: i for i in infos}
    try:
        tree = T.parse(apply_atoms(tex, infos), euler_e=euler)
        if not isinstance(tree, T.Rel):
            raise T.Unreadable("the LaTeX has no relation (= < > ≤ ≥ ≠), so it states nothing")
    except T.Unreadable as exc:
        row.reason = str(exc)
        return row

    # 3. is it still the model?
    if rel is not None:
        try:
            back = T.to_sympy(tree, {c: i.sym for c, i in by_canon.items()})
        except Exception as exc:
            row.reason = f"could not rebuild sympy from the LaTeX: {exc}"
            return row
        method, note = compare(back, rel)
        if method is None:
            row.status, row.reason = "mismatch", note
            return row
        row.method, row.reason = method, note
        row.status = "verified"
    else:
        row.status = "tex-only"
        row.reason = "no sympy relation to compare the LaTeX with"
        if not tex_only:
            return row

    # 4. binders from sympy's assumptions
    try:
        binders = []
        names = symbol_names(infos)
        for name in T.symbols(tree):
            info = by_canon.get(name)
            b = binder_for(info) if info is not None else E.Binder(name)
            binders.append(b)
            row.binder_tex.append((sp.latex(info.sym, symbol_names=names) if info is not None else name,
                                   b.kind, b.assume))
    except T.Unreadable as exc:
        row.status, row.reason = "refused", str(exc)
        return row

    # 5. the theorem, once per flavor
    tactic = "sorry"
    if prove and rel is not None and is_field_identity(rel):
        row.proof_attempted = True
    for flavor in flavors:
        t = "first | ring1 | (field_simp; ring1) | sorry" if (row.proof_attempted and flavor == "mathlib") else "sorry"
        try:
            row.theorems[flavor] = E.build_theorem(row.key, tree, binders, flavor, doc_comment(row, tex), t)
        except T.Unreadable as exc:
            row.status, row.reason, row.theorems = "refused", str(exc), {}
            return row
    return row


# ----------------------------------------------------------------------------
# The registry
# ----------------------------------------------------------------------------
def load_registry(models: list[str] | None):
    """Run the chosen sympy_models and return their registered equations."""
    if str(SYMMODELS) not in sys.path:
        sys.path.append(str(SYMMODELS))
    import eqlib
    pkg = importlib.import_module("models")
    available = sorted(m.name for m in pkgutil.iter_modules(pkg.__path__))
    chosen = available if models == ["all"] else (models or [m for m in DEFAULT_MODELS if m in available])
    unknown = sorted(set(chosen) - set(available))
    if unknown:
        raise SystemExit(f"unknown model(s): {', '.join(unknown)}; have: {', '.join(available)}")
    for name in chosen:
        eqlib.REGISTRY  # noqa: B018 (import side effects only)
        importlib.import_module(f"models.{name}").build_and_verify()
    return chosen, list(eqlib.REGISTRY.values())


# ----------------------------------------------------------------------------
# Output
# ----------------------------------------------------------------------------
def tex_escape(text: str) -> str:
    """Plain text for LaTeX, by the fork's own rules (``sympy_models/document.esc``), so titles and
    sources look here as they do in the fork's other documents. Characters it cannot typeset become ``?``."""
    try:
        return _document.esc(text)
    except ValueError:
        return _document.esc(text.encode("ascii", "replace").decode())


def tt(text: str) -> str:
    """A key or Lean name in \\texttt, allowed to break after an underscore."""
    return "\\texttt{" + tex_escape(text).replace(r"\_", r"\_\allowbreak{}") + "}"


def longtable(header: list[str], rows: list[list[str]], colspec: str, caption: str) -> str:
    """A table that may break across pages (``common.table`` makes a float, which may not)."""
    head = " & ".join(header) + " \\\\"
    lines = [f"\\begin{{longtable}}{{{colspec}}}", f"\\caption{{{caption}}}\\\\", "\\toprule", head, "\\midrule",
             "\\endfirsthead", "\\toprule", head, "\\midrule", "\\endhead", "\\bottomrule", "\\endfoot"]
    lines += [" & ".join(r) + " \\\\" for r in rows]
    lines.append("\\end{longtable}")
    return "\n".join(lines) + "\n"


_RAGGED = r">{\raggedright\arraybackslash}p{%s}"


def quantities_sentence(row: Row) -> str:
    kinds = {"real": r"\mathbb{R}", "int": r"\mathbb{Z}", "nat": r"\mathbb{N}"}
    items = []
    for tex, kind, assume in row.binder_tex:
        cond = {"positive": ">0", "nonnegative": r"\ge 0", "negative": "<0", "nonpositive": r"\le 0",
                "nonzero": r"\ne 0"}
        extra = ", ".join(f"{tex}{cond[a]}" if a != "nonzero" else f"{tex}{cond[a]}" for a in assume)
        items.append(f"${tex}\\in {kinds[kind]}$" + (f" with ${extra}$" if extra else ""))
    return ", ".join(items)


def tex_fragment(rows: list[Row], flavors, lean_note: str) -> str:
    counts = {s: sum(1 for r in rows if r.status == s) for s in STATUSES}
    out = [
        "% Generated by tools/sym2lean.py. Do not edit by hand.",
        "\\section{Model equations as Lean statements}",
        ("Each equation below is a sympy relation from \\texttt{sympy\\_models/}. It was printed as the LaTeX shown, "
         "read back and compared with the sympy expression, and only then written as a Lean~4 theorem. "
         "A theorem here is a \\emph{statement}: unless marked proved, its body is \\texttt{sorry}, and "
         "Lean accepting it says only that it is a well-formed proposition. "
         "Where an equation has a denominator, Lean's \\texttt{x / 0 = 0} makes the statement easier to "
         "satisfy than the physics is, so the sympy assumptions on each quantity are carried over as hypotheses."),
        "",
        (f"Of {len(rows)} equations: " + ", ".join(f"{counts[s]} {s}" for s in STATUSES) + ". " + tex_escape(lean_note)),
        "",
    ]
    done = [r for r in rows if r.theorems or r.status in ("verified", "mismatch")]
    if done:
        table_rows = []
        for r in done:
            lean = "; ".join(f"{f}: {r.lean.get(f, 'unchecked')}" for f in flavors) if r.theorems else "-"
            table_rows.append([tt(r.key), tex_escape(r.title), r.status + (f" ({r.method})" if r.method else ""),
                               tex_escape(lean)])
        out.append(longtable(["key", "title", "LaTeX against sympy", "Lean"], table_rows,
                             _RAGGED % "3.1cm" + _RAGGED % "5.3cm" + _RAGGED % "2.6cm" + _RAGGED % "3.4cm",
                             "Equations written as Lean theorems, and what Lean made of them."))
    why: dict[str, list[str]] = {}
    for r in rows:
        if r.status in ("refused", "tex-only") and not r.theorems:
            reason = re.sub(r" at position \d+", "", _squash(r.reason))[:140]
            why.setdefault(reason or "(no reason recorded)", []).append(r.key)
    if why:
        reason_rows = [[tex_escape(reason), str(len(keys)), ", ".join(tt(k) for k in keys)]
                       for reason, keys in sorted(why.items(), key=lambda kv: -len(kv[1]))]
        out.append(longtable(["why it was not translated", "n", "equations"], reason_rows,
                             _RAGGED % "5.2cm" + "r" + _RAGGED % "8.2cm",
                             "Equations the pipeline declined, by reason. A refusal names what it could not "
                             "read; none is a guess."))
    out.append("\\clearpage")
    for r in rows:
        if r.status != "verified" and not r.theorems:
            continue
        out.append(f"\\subsection*{{{tex_escape(r.title or r.key)}}}")
        if r.entry:
            out.append(f"Repository entry {tex_escape(r.entry)}. ")
        latex = r.latex.replace("%%BR%%", "\\quad ")
        out.append(C.equation(latex, f"s2l:{r.key}"))
        if r.binder_tex:
            out.append(f"Lean quantifies over {quantities_sentence(r)}.")
        facts = []
        if r.source:
            facts.append(f"Source ({tex_escape(r.source_kind)}): {tex_escape(r.source)}.")
        facts.append(f"LaTeX against sympy: {r.method or 'not compared'}" + (f" ({tex_escape(r.reason)})" if r.reason else "") + ".")
        theorem = next((t.name for t in r.theorems.values()), None)
        if theorem:
            lean = ", ".join(f"{f}: {r.lean.get(f, 'unchecked')}" for f in flavors if f in r.theorems)
            facts.append(f"Lean theorem \\texttt{{{tex_escape(theorem)}}} ({lean}).")
        out.append("\\par\\small " + " ".join(facts) + "\\normalsize\n")
    return "\n".join(out) + "\n"


def write_all(rows: list[Row], flavors, models: list[str], check: bool) -> dict:
    C.ensure_dirs()
    lean_dir = C.BUILD / "lean"
    lean_dir.mkdir(parents=True, exist_ok=True)
    results = {}
    for flavor in flavors:
        theorems = [r.theorems[flavor] for r in rows if flavor in r.theorems]
        title = "Model equations, as Lean statements"
        blurb = ("Generated by `tools/sym2lean.py` from sympy models: sympy -> LaTeX -> Lean 4. Each theorem's "
                 "doc comment carries the LaTeX it was read from and its provenance.\n"
                 f"Models: {', '.join(models) or '(none)'}.")
        src = E.lean_file(theorems, flavor, title, blurb)
        path = lean_dir / f"SymModels_{flavor}.lean"
        path.write_text(src)
        info = {"path": str(path.relative_to(C.ROOT)), "theorems": len(theorems)}
        if check and theorems:
            res = E.run_lean(str(path), flavor, project=str(LEAN_PROJECT))
            if res is None:
                info["checked"] = "no lean on this machine; nothing claimed"
            elif "skipped" in res:
                info["checked"] = "skipped: " + res["skipped"]
            else:
                status = E.theorem_status(src, res["output"])
                for r in rows:
                    if flavor in r.theorems:
                        r.lean[flavor] = status.get(r.theorems[flavor].name, "unchecked")
                info.update(checked="lean " + res["version"], ok=res["ok"], errors=res["errors"][:10],
                            proved=sum(v == "proved" for v in status.values()),
                            statements=sum(v == "statement" for v in status.values()),
                            rejected=sum(v == "error" for v in status.values()))
        results[flavor] = info
    return results


# ----------------------------------------------------------------------------
# Self test
# ----------------------------------------------------------------------------
def selftest() -> int:
    n = 0

    def ok(cond, msg):
        nonlocal n
        if not cond:
            raise SystemExit(f"sym2lean selftest FAILED: {msg}")
        n += 1

    class Eq_:                      # the fields process() reads, without needing eqlib
        def __init__(self, key, expr=None, latex_str=None, **kw):
            self.key, self.expr, self.latex_str = key, expr, latex_str
            self.__dict__.update(kw)

        def latex(self, multiline=False):
            return self.latex_str if self.latex_str is not None else sp.latex(self.expr)

    x, y, T_ = sp.symbols("x y T", positive=True)
    Nsites, e_sym, hbar, kB = sp.symbols("Nsites e hbar k_B", positive=True)
    nn = sp.Symbol("n", integer=True, nonnegative=True)

    cases = [
        ("simple", sp.Eq(x, y**2 / 2)),
        ("exp-and-log", sp.Eq(x, sp.exp(-y) + sp.log(T_))),
        ("multi-letter", sp.Eq(Nsites, x * y)),
        ("symbol e", sp.Eq(e_sym, x**2)),
        ("subscripts", sp.Eq(kB * T_, hbar * x)),
        ("roots and powers", sp.Eq(x, sp.sqrt(y) * T_ ** sp.Rational(3, 2) + y ** -2)),
        ("trig", sp.Eq(x, sp.sin(y) ** 2 + sp.cos(y) ** 2)),
        ("pi and e", sp.Eq(x, 2 * sp.pi * sp.E * y)),
        ("natural exponent", sp.Eq(x, y**nn)),
        ("inequality", sp.Le(x, y + 1)),
        ("floats", sp.Eq(x, sp.Float("0.5") * y)),
        ("abs", sp.Eq(x, sp.Abs(y - T_))),
    ]
    # the fork's convention: a symbol's *name* is LaTeX source
    beta_k = sp.Symbol("\\beta_k", positive=True)
    aL = sp.Symbol("a\\Lambda_L", positive=True)
    smb = sp.Symbol("S-m_S(\\beta)", positive=True)
    cases += [("latex-named symbol", sp.Eq(beta_k, x + y / 2)),
              ("composite name", sp.Eq(aL, sp.exp(-1 / x) * y**2)),
              ("composite on both sides", sp.Eq(smb, x * aL**2 + beta_k)),
              ("greek subscripts", sp.Eq(sp.Symbol("mu_0", positive=True), sp.Symbol("hbar", positive=True) * x))]
    for name, rel in cases:
        row = process(Eq_(name, rel, title=name))
        ok(row.status == "verified", f"{name}: status {row.status} ({row.reason}); latex {row.latex}")
        ok(set(row.theorems) == {"mathlib", "core"}, f"{name}: a theorem per flavor")

    # fuzz: random expressions, printed by sympy and read back by texmath
    rng = random.Random(7)
    syms = [x, y, T_]

    def rand(depth):
        if depth == 0 or rng.random() < 0.25:
            return rng.choice(syms + [sp.Integer(rng.randint(1, 5)), sp.Rational(rng.randint(1, 4), rng.randint(2, 5)), sp.pi])
        a, b = rand(depth - 1), rand(depth - 1)
        k = rng.randrange(8)
        return (a + b, a - b, a * b, a / (b + 1), a ** rng.choice((2, 3, -1, sp.Rational(1, 2))),
                sp.sqrt(a**2 + 1), sp.exp(a / (1 + b**2)), sp.log(a**2 + b**2 + 1))[k]

    bad = []
    for i in range(150):
        rhs = rand(3)
        if rhs.free_symbols == set() and rhs.is_number and rhs.has(sp.zoo):
            continue
        row = process(Eq_(f"fuzz{i}", sp.Eq(sp.Symbol("w", positive=True), rhs, evaluate=False)))
        if row.status != "verified":
            bad.append((sp.latex(rhs), row.status, row.reason))
    ok(not bad, f"fuzz: {len(bad)} of 150 failed, first: {bad[:2]}")

    # a latex_str that does not say what the model says
    row = process(Eq_("liar", sp.Eq(x, y**2), latex_str=r"x = y^{3}"))
    ok(row.status == "mismatch", f"a wrong latex_str must be caught, got {row.status}")
    row = process(Eq_("rearranged", sp.Eq(x, y + T_), latex_str=r"x - y = T"))
    ok(row.status == "verified" and "rearranged" in row.reason, f"rearranged: {row.status} {row.reason}")
    # beyond the subset: refused, with a reason, never guessed
    f = sp.Function("f")
    row = process(Eq_("derivative", sp.Eq(x, sp.Derivative(f(y), y))))
    ok(row.status == "refused" and row.reason, f"derivative: {row.status} {row.reason}")
    row = process(Eq_("applied", sp.Eq(x, f(y) + 1)))
    ok(row.status == "refused" and "function" in row.reason, f"applied function: {row.reason}")
    row = process(Eq_("zeta", sp.Eq(x, sp.zeta(sp.Rational(3, 2)) * y)))
    ok(row.status == "refused" and "zeta" in row.reason, f"an applied special function is named: {row.reason}")
    row = process(Eq_("hand-composite", sp.Eq(aL, sp.exp(-1 / x) * y**2),
                      latex_str=r"a\Lambda_L = \exp\!\Bigl(-\frac{1}{x}\Bigr)\,y^{2}"))
    ok(row.status == "verified", f"hand-written LaTeX naming a composite symbol: {row.status} {row.reason}")
    disguised = sp.Symbol(r"\lim \frac{d\ln(a)}{d(1/g^2)}", positive=True)
    row = process(Eq_("disguised", sp.Eq(disguised, x)))
    ok(row.status == "refused" and "operator" in row.reason, f"an operator named as a symbol: {row.status} {row.reason}")
    row = process(Eq_("little-o", None, latex_str=r"x = y + o\bigl(y^{2}\bigr)"))
    ok(row.status == "refused" and "asymptotic" in row.reason, f"little-o: {row.reason}")
    row = process(Eq_("hand-e", None, latex_str=r"x = e^{y}"))
    ok(row.status == "refused" and "ambiguous" in row.reason, f"hand-written e^: {row.reason}")
    row = process(Eq_("sum", None, latex_str=r"x = \sum_{i} y_i"))
    ok(row.status == "refused" and "binder" in row.reason, f"sum: {row.reason}")
    row = process(Eq_("handwritten", None, latex_str=r"E = m c^{2}"))
    ok(row.status == "tex-only" and not row.theorems, "tex-only is held back by default")
    row = process(Eq_("handwritten", None, latex_str=r"E = m c^{2}"), tex_only=True)
    ok(row.status == "tex-only" and "mathlib" in row.theorems, "tex-only is written with --tex-only")

    # assumptions become hypotheses; identities get a proof attempt
    row = process(Eq_("hyp", sp.Eq(x, hbar / y)), prove=True)
    ok("(h_x : 0 < x)" in row.theorems["mathlib"].text() and "(h_y : 0 < y)" in row.theorems["mathlib"].text(), "hypotheses")
    row = process(Eq_("ident", sp.Eq((x + y) ** 2, x**2 + 2 * x * y + y**2)), prove=True)
    ok(row.proof_attempted and "ring1" in row.theorems["mathlib"].text() and "ring1" not in row.theorems["core"].text(),
       "ring1 only for mathlib identities")
    row = process(Eq_("law", sp.Eq(x, y**2)), prove=True)
    ok(not row.proof_attempted, "a law is not an identity")
    row = process(Eq_("nat", sp.Eq(x, y**nn)))
    ok("(n : ℕ)" in row.theorems["mathlib"].text(), "natural-number symbols are typed ℕ")

    # LaTeX for the document
    ok(tex_escape("a_b & 50%") == r"a\_b \& 50\%", tex_escape("a_b & 50%"))
    ok(tex_escape("T^{3/2}") == r"T\textasciicircum{}\{3/2\}", "braces are not escaped twice")
    ok(tt("on_mass_af") == r"\texttt{on\_\allowbreak{}mass\_\allowbreak{}af}", tt("on_mass_af"))
    return n


# ----------------------------------------------------------------------------
# Command line
# ----------------------------------------------------------------------------
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("models", nargs="*", help="sympy_models names, or 'all' (default: %s)" % ", ".join(DEFAULT_MODELS))
    ap.add_argument("--flavor", choices=("mathlib", "core", "both"), default="both")
    ap.add_argument("--check", action="store_true", help="give the files to Lean if there is one")
    ap.add_argument("--prove", action="store_true", help="try ring1/field_simp on rational-function identities (mathlib)")
    ap.add_argument("--tex-only", action="store_true", help="also write Lean for hand-written LaTeX with no sympy to compare")
    ap.add_argument("--strict", action="store_true", help="exit nonzero if any LaTeX does not match its model")
    ap.add_argument("--list", action="store_true", help="census of the registry; write nothing")
    ap.add_argument("--latex", metavar="STATEMENT", help="translate one hand-written LaTeX statement and print the Lean")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args(argv)

    if args.selftest:
        print(f"sym2lean: {selftest()} checks ok")
        return 0
    flavors = ("mathlib", "core") if args.flavor == "both" else (args.flavor,)

    if args.latex:
        class One:
            key, title, entry, source_kind, source, notes, expr = "latex", "", "", "", "", "", None
            latex_str = args.latex

            @staticmethod
            def latex(multiline=False):
                return args.latex
        row = process(One, flavors, args.prove, tex_only=True)
        if not row.theorems:
            print(f"{row.status}: {row.reason}")
            return 1
        for f in flavors:
            print(f"-- {f}\n{row.theorems[f].text()}\n")
        print("-- hand-written LaTeX has no sympy model to compare with: unverified, all quantities real, no hypotheses")
        return 0

    models, eqs = load_registry(args.models)
    rows = [process(eq, flavors, args.prove, args.tex_only) for eq in eqs]
    counts = {s: sum(1 for r in rows if r.status == s) for s in STATUSES}

    print(f"models: {', '.join(models)}")
    print(f"{len(rows)} equations: " + ", ".join(f"{counts[s]} {s}" for s in STATUSES))
    reasons: dict[str, int] = {}
    for r in rows:
        if r.status in ("refused", "tex-only"):
            reasons[r.reason.split(":")[0][:90]] = reasons.get(r.reason.split(":")[0][:90], 0) + 1
    for why, k in sorted(reasons.items(), key=lambda kv: -kv[1])[:8]:
        print(f"    {k:3d}  {why}")
    for r in rows:
        if r.status == "mismatch":
            print(f"MISMATCH {r.key}: {r.reason}\n    latex: {_squash(r.latex)}")
    if args.list:
        for r in rows:
            print(f"  {r.status:9s} {r.key}")
        return 0

    info = write_all(rows, flavors, models, args.check)
    for f, d in info.items():
        line = f"{f:8s} {d['theorems']} theorems -> {d['path']}"
        if "checked" in d:
            line += f"  [{d['checked']}" + (f"; {d['proved']} proved, {d['statements']} statements, {d['rejected']} rejected" if "proved" in d else "") + "]"
        print(line)
        for err in d.get("errors", []):
            print("    " + err)
    note = ("Lean was not asked." if not args.check else
            "; ".join(f"{f}: {d.get('checked', 'not checked')}" for f, d in info.items()) + ".")
    C.write_tex("sym2lean", tex_fragment(rows, flavors, note))
    C.save_json("sym2lean", {"models": models, "census": counts, "files": info,
                             "rows": [r.summary() for r in rows]})
    print(f"wrote build/tex/sym2lean.tex and build/data/sym2lean.json")
    return 1 if (args.strict and counts["mismatch"]) else 0


if __name__ == "__main__":
    sys.exit(main())
