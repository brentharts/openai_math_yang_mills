"""Arithmetic theorems of the OpenAI release's Lean files, read into sympy and into leandsl items.

A theorem or lemma qualifies when every binder is a number of type ℝ, ℕ, ℤ or ℚ or a hypothesis, and its hypotheses and
conclusion use nothing but variables, numerals, ``+ - * / ^ ⁻¹ |x|``, ``Real.sqrt/log/exp``, ``Nat.choose/factorial``,
``min``/``max`` and the order relations. Anything else (a definition of the release, a quantifier, a sum, a structure, a
variable declared elsewhere in the file) is not read, and the reason is counted. Nothing is typed in: the statement is
parsed from the Lean text, and the Lean text stays beside it.

What is checked, and what is not:

* The reading is a small parser, so each statement is checked against the Lean text two ways: printed back with
  leandsl's Lean printer and parsed again, it must give the same sympy statement; and sympy must find the statement true
  at points that satisfy its hypotheses (the Lean statement is proved in the release, so a counterexample means the
  parser misread it, and the statement is dropped, not shown). ``sampled`` is evidence and no more.
* Natural-number subtraction and division truncate in Lean; a statement that uses either on natural numbers is not read.
* Lean's total functions (x/0 = 0, log 0 = 0, sqrt of a negative = 0) are not modelled: points where sympy would
  divide by zero or leave the reals are skipped, so such a statement is checked on fewer points.
* ``uses`` links between items come from names: the proof of A mentions the name of B, which must be unique among the
  indexed declarations. That is a textual fact, not a checked dependence (a name can occur in a comment or a string).

``python3 tools/oai_lean.py --family 275 [--limit N]`` prints what was read from a family; ``--selftest`` runs the checks.
"""
from __future__ import annotations

import collections
import re
import sys
from dataclasses import dataclass, field

import sympy as sp

import leandsl
from leandsl import Real, Nat, Int, _sym

TYPES = {"ℝ": Real, "ℕ": Nat, "ℤ": Int, "ℚ": Real}
_NORM = {"Real": "ℝ", "Nat": "ℕ", "Int": "ℤ", "Rat": "ℚ"}         # the spelled-out type names
_RANK = {None: 0, "ℕ": 1, "ℤ": 2, "ℚ": 3, "ℝ": 4}


class Unsupported(ValueError):
    """The statement uses something outside the arithmetic subset; the message is the reason."""


# ---- reading the file --------------------------------------------------------------------------------
_COMMENTS = re.compile(r"/-.*?-/|--[^\n]*", re.S)
_DECL = re.compile(r"^(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|noncomputable\s+)*(theorem|lemma)\s+([^\s:({\[⦃]+)", re.M)
_NS = re.compile(r"^(namespace|section|end)\b[ \t]*([^\s]*)", re.M)


@dataclass
class Decl:
    module: str
    kind: str                 # "theorem" | "lemma"
    name: str                 # as written
    qualified: str            # with the enclosing namespaces
    binders: list             # [(open bracket, text)]
    concl: str                # conclusion text, after the colon
    body: str = ""            # the proof text, comments removed
    text: str = ""            # the whole statement as written, comments removed


def _binders(rest: str):
    i, n, out = 0, len(rest), []
    close = {"(": ")", "{": "}", "[": "]", "⦃": "⦄"}
    while i < n:
        c = rest[i]
        if c in close:
            d, j = 1, i + 1
            while j < n and d:
                d += (rest[j] in close) - (rest[j] in close.values())
                j += 1
            out.append((c, rest[i + 1:j - 1]))
            i = j
        elif c == ":":
            return out, rest[i + 1:].strip()
        elif c.isspace():
            i += 1
        else:
            return None
    return None


def declarations(text: str, module: str = "") -> list[Decl]:
    """Every ``theorem``/``lemma`` of a Lean file whose statement is in the form binders ``:`` conclusion."""
    text = _COMMENTS.sub(" ", text)
    spans, stack = [], []                                  # namespace in force at each position
    for m in _NS.finditer(text):
        kw, arg = m.group(1), m.group(2)
        if kw == "namespace":
            stack.append((arg, True))
        elif kw == "section":
            stack.append((arg, False))
        elif stack and (not arg or stack[-1][0] == arg or stack[-1][0].endswith("." + arg)):
            stack.pop()
        spans.append((m.start(), ".".join(a for a, isns in stack if isns and a)))
    out = []
    ds = list(_DECL.finditer(text))
    for k, d in enumerate(ds):
        end = ds[k + 1].start() if k + 1 < len(ds) else len(text)
        blk = text[d.start():end]
        cut = blk.find(":=")
        stmt, body = (blk[:cut], blk[cut:]) if cut >= 0 else (blk, "")
        head = re.match(r"\s*(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|noncomputable\s+)*(?:theorem|lemma)\s+\S+", stmt)
        r = _binders(stmt[head.end():]) if head else None
        if r is None:
            continue
        ns = next((n for pos, n in reversed(spans) if pos < d.start()), "")
        nm = d.group(2)
        qual = nm[len("_root_."):] if nm.startswith("_root_.") else (f"{ns}.{nm}" if ns else nm)
        out.append(Decl(module, d.group(1), nm, qual, r[0], r[1], body, stmt.strip()))
    return out


# ---- the expression parser --------------------------------------------------------------------------
_TOKEN = re.compile(r"""\s*(?:
    (?P<num>\d+(?:\.\d+)?)
  | (?P<id>[A-Za-zα-ωΑ-Ωℓ_][A-Za-z0-9_'α-ωΑ-Ω₀-₉.]*)
  | (?P<sym>≤|≥|≠|∧|⁻¹|[<>=()|^+\-*/:↑,×·ℝℕℤℚ])
)""", re.X)
_FUNCS = {"Real.sqrt": ("ℝ", 1), "Real.log": ("ℝ", 1), "Real.exp": ("ℝ", 1), "Nat.choose": ("ℕ", 2),
          "Nat.factorial": ("ℕ", 1), "Real.rpow": ("ℝ", 2)}
_POLY = {"min": 2, "max": 2, "abs": 1}
_CONST = {"Real.pi": "ℝ", "π": "ℝ"}
_REL = {"≤": sp.Le, "≥": sp.Ge, "<": sp.Lt, ">": sp.Gt, "=": sp.Eq, "≠": sp.Ne}


def _tokens(s: str):
    pos, out = 0, []
    s = s.strip()
    while pos < len(s):
        m = _TOKEN.match(s, pos)
        if not m or m.end() == pos:
            raise Unsupported(f"unread character {s[pos:pos + 1]!r}")
        pos = m.end()
        out.append((m.lastgroup, m.group(m.lastgroup)))
        while pos < len(s) and s[pos].isspace():
            pos += 1
    return out


class _P:
    """A precedence parser for Lean arithmetic: ``+ -`` 65, ``* /`` 70, ``^`` 75 (right), prefix ``-`` 75, ``⁻¹`` tighter."""

    def __init__(self, toks, vars_):
        self.t, self.i, self.vars, self.abs = toks, 0, vars_, 0

    def peek(self):
        return self.t[self.i] if self.i < len(self.t) else (None, None)

    def eat(self, val=None):
        k, v = self.peek()
        if k is None or (val is not None and v != val):
            raise Unsupported(f"expected {val!r}")
        self.i += 1
        return v

    def relation(self):
        lhs = self.expr(0)
        k, v = self.peek()
        if v not in _REL:
            raise Unsupported("not a relation")
        self.i += 1
        return (v, lhs, self.expr(0))

    def expr(self, minp):
        lhs = self.prefix()
        while True:
            k, v = self.peek()
            prec = {"+": 65, "-": 65, "*": 70, "/": 70, "×": 70, "·": 70}.get(v) if k == "sym" else None
            if v == "^" and k == "sym":
                if 75 < minp:
                    return lhs
                self.i += 1
                lhs = ("pow", lhs, self.expr(75))            # right associative
                continue
            if prec is None or prec < minp:
                return lhs
            self.i += 1
            rhs = self.expr(prec + 1)
            lhs = ({"+": "add", "-": "sub", "*": "mul", "/": "div", "×": "mul", "·": "mul"}[v], lhs, rhs)

    def prefix(self):
        k, v = self.peek()
        if k == "sym" and v == "-":
            self.i += 1
            return ("neg", self.expr(75))
        return self.app()

    def app(self):
        k, v = self.peek()
        if k == "id" and (v in _FUNCS or v in _POLY):
            self.i += 1
            arity = _FUNCS[v][1] if v in _FUNCS else _POLY[v]
            return ("fn", v, [self.atom() for _ in range(arity)])
        return self.atom()

    def starts_atom(self):
        k, v = self.peek()
        return k in ("num", "id") or v in ("(", "↑") or (v == "|" and self.abs == 0)

    def atom(self):
        k, v = self.peek()
        if k == "num":
            self.i += 1
            node = ("num", v)
        elif k == "id":
            self.i += 1
            if v in _CONST:
                node = ("const", v)
            elif v in self.vars:
                node = ("var", v)
            else:
                raise Unsupported(f"unknown name {v}")
        elif v == "(":
            self.i += 1
            inner = self.expr(0)
            if self.peek()[1] == ":":
                self.i += 1
                ty = self.eat()
                ty = _NORM.get(ty, ty)
                if ty not in TYPES:
                    raise Unsupported(f"type {ty}")
                inner = ("asc", inner, ty)
            self.eat(")")
            node = inner
        elif v == "|":
            self.i += 1
            self.abs += 1
            inner = self.expr(0)
            self.abs -= 1
            self.eat("|")
            node = ("fn", "abs", [inner])
        elif v == "↑":
            self.i += 1
            return self.atom()
        else:
            raise Unsupported(f"unexpected {v!r}")
        while self.peek()[1] == "⁻¹":
            self.i += 1
            node = ("inv", node)
        return node


def _tree_type(n, vtypes):
    """The type Lean unifies an arithmetic tree to (None while only numerals are seen); function arguments and
    ascribed expressions are elaborated on their own and do not take part."""
    op = n[0]
    if op == "num":
        return None
    if op == "const":
        return "ℝ"
    if op == "var":
        return vtypes[n[1]]
    if op == "asc":
        return n[2]
    if op == "fn":
        if n[1] in _FUNCS:
            return _FUNCS[n[1]][0]
        return _join(_tree_type(a, vtypes) for a in n[2])
    if op == "pow":
        return _tree_type(n[1], vtypes)
    return _join(_tree_type(c, vtypes) for c in n[1:])


def _join(ts):
    best = None
    for t in ts:
        if _RANK[t] > _RANK[best]:
            best = t
    return best


def _check_semantics(n, ty, vtypes):
    """Reject what Lean computes differently from the reals: ℕ/ℤ subtraction, division and inverse truncate."""
    op = n[0]
    if op in ("sub", "neg") and ty in ("ℕ", None) and op == "sub":
        raise Unsupported("subtraction on natural numbers (truncated in Lean)")
    if op == "neg" and ty == "ℕ":
        raise Unsupported("negation of a natural number")
    if op in ("div", "inv") and ty in ("ℕ", "ℤ", None):
        raise Unsupported("division on integers (floor division in Lean)")
    if op in ("num", "const", "var"):
        return
    if op == "asc":
        _check_semantics(n[1], n[2], vtypes)
    elif op == "fn":
        if n[1] in _FUNCS:
            for a in n[2]:
                _check_semantics(a, _FUNCS[n[1]][0], vtypes)
        else:
            for a in n[2]:
                _check_semantics(a, ty, vtypes)
    elif op == "pow":
        _check_semantics(n[1], ty, vtypes)
        e = n[2]
        et = _tree_type(e, vtypes) or "ℕ"                  # an exponent is elaborated on its own: ℕ unless it says otherwise
        _check_semantics(e, et, vtypes)
    else:
        for c in n[1:]:
            _check_semantics(c, ty, vtypes)


def _to_sympy(n, syms, ty):
    op = n[0]
    if op == "num":
        s = n[1]
        return sp.Integer(s) if s.isdigit() else sp.Rational(s)
    if op == "const":
        return sp.pi
    if op == "var":
        return syms[n[1]]
    if op == "asc":
        return _to_sympy(n[1], syms, n[2])
    if op == "neg":
        return -_to_sympy(n[1], syms, ty)
    if op == "inv":
        return sp.Pow(_to_sympy(n[1], syms, ty), -1)
    a = [_to_sympy(c, syms, ty) for c in n[1:] if isinstance(c, tuple)]
    if op == "add":
        return sp.Add(*a)
    if op == "sub":
        return a[0] - a[1]
    if op == "mul":
        return sp.Mul(*a)
    if op == "div":
        return a[0] / a[1]
    if op == "pow":
        return sp.Pow(a[0], a[1])
    if op == "fn":
        f, args = n[1], [_to_sympy(x, syms, ty) for x in n[2]]
        return {"Real.sqrt": lambda x: sp.sqrt(x), "Real.log": sp.log, "Real.exp": sp.exp, "Nat.choose": sp.binomial,
                "Nat.factorial": sp.factorial, "Real.rpow": sp.Pow, "min": sp.Min, "max": sp.Max, "abs": sp.Abs}[f](*args)
    raise Unsupported(f"node {op}")


@dataclass
class Statement:
    params: list              # [(typed symbol, leandsl type)]
    hyps: list                # relations (sympy, unevaluated)
    concl: object
    lean: str                 # the Lean statement text, as written (comments removed)
    hyp_texts: list = field(default_factory=list)


def read_statement(d: Decl) -> Statement:
    """The statement of a declaration as sympy, or Unsupported(reason)."""
    vtypes, order, hyps = {}, [], []
    for br, b in d.binders:
        if br == "[":
            raise Unsupported("instance argument")
        m = re.fullmatch(r"\s*([^:]+?)\s*:\s*(.+?)\s*", b, flags=re.S)
        if not m:
            raise Unsupported("binder without a type")
        names, ty = m.group(1).split(), _NORM.get(m.group(2).strip(), m.group(2).strip())
        if ty in TYPES:
            for nm in names:
                vtypes[nm] = ty
                order.append(nm)
        else:
            hyps.append(ty)
    if not order:
        raise Unsupported("no numeric variable")
    syms = {nm: _sym(nm, TYPES[vtypes[nm]]) for nm in order}
    parsed_h = []
    for h in hyps:
        for part in re.split(r"\s*∧\s*", h):
            parsed_h.append(_relation(part, vtypes, syms))
    c = _relation(d.concl, vtypes, syms)
    return Statement([(syms[nm], TYPES[vtypes[nm]]) for nm in order], parsed_h, c, d.text, hyps)


def _relation(text: str, vtypes, syms):
    if re.search(r"∨|¬|↔|→|∀|∃|fun|∑|∏|∫|:=|\bat\b|%|∣|\bif\b|\bthen\b|\bBool|\bTrue|\bFalse", text):
        raise Unsupported("logic or syntax outside the subset")
    if "∧" in text:
        raise Unsupported("a conjunction as the conclusion")
    toks = _tokens(text)
    p = _P(toks, vtypes)
    rel, lhs, rhs = p.relation()
    if p.i != len(toks):
        raise Unsupported("trailing text")
    ty = _join([_tree_type(lhs, vtypes), _tree_type(rhs, vtypes)])
    if ty is None:
        raise Unsupported("a relation between numerals only")
    _check_semantics(lhs, ty, vtypes)
    _check_semantics(rhs, ty, vtypes)
    return _REL[rel](_to_sympy(lhs, syms, ty), _to_sympy(rhs, syms, ty), evaluate=False)


# ---- checking a statement ----------------------------------------------------------------------------
def lean_text(st: Statement) -> str:
    """The statement printed by leandsl, to be parsed again."""
    return leandsl.lean_expr(st.concl.lhs) + {sp.Le: " ≤ ", sp.Ge: " ≥ ", sp.Lt: " < ", sp.Gt: " > ", sp.Eq: " = ", sp.Ne: " ≠ "}[type(st.concl)] \
        + leandsl.lean_expr(st.concl.rhs)


def check_statement(d: Decl, st: Statement, want: int = 60, tries: int = 12000) -> tuple[str, str]:
    """``(status, detail)``: "sampled" (held at sampled points), "unchecked" (no point satisfying the hypotheses was
    found, or it cannot be evaluated) or "refuted" (false at a point: the parser misread it)."""
    item = leandsl.Item("theorem", d.name, st.params, [h for h in st.hyps], st.concl)
    c = st.concl
    kind = {sp.Le: "<=", sp.Lt: "<", sp.Ge: ">=", sp.Gt: ">", sp.Eq: "=", sp.Ne: "!="}[type(c)]
    if kind == "!=":
        return "unchecked", "a disequality is not sampled"
    triple = {"<=": (c.lhs, c.rhs, "le"), "<": (c.lhs, c.rhs, "<"), ">=": (c.rhs, c.lhs, "le"),
              ">": (c.rhs, c.lhs, "<"), "=": (c.lhs, c.rhs, "=")}[kind]
    res, feasible = leandsl._sample(item, [triple], want, tries)
    n_ok, bad, evaluable = res[0]
    if bad is not None:
        return "refuted", f"false at {', '.join(f'{k}={v}' for k, v in bad.items())}"
    if not evaluable:
        return "unchecked", "cannot be evaluated numerically"
    if n_ok >= 10:
        return "sampled", f"held at {n_ok} points satisfying the hypotheses"
    return "unchecked", f"only {n_ok} point(s) satisfying the hypotheses were found"


@dataclass
class Read:
    decl: Decl
    stmt: Statement | None = None
    status: str = ""          # sampled | unchecked | refuted | unsupported
    detail: str = ""


def read_declaration(d: Decl, sample: bool = True) -> Read:
    try:
        st = read_statement(d)
    except Unsupported as ex:
        return Read(d, None, "unsupported", str(ex))
    except (TypeError, ValueError, ZeroDivisionError, RecursionError, sp.SympifyError) as ex:
        return Read(d, None, "unsupported", f"not read: {type(ex).__name__}")
    # a round trip through leandsl's Lean printer must give the same statement
    try:
        back = _relation(lean_text(st), {leandsl.lean_ident(p.name): ("ℕ" if t is Nat else "ℤ" if t is Int else "ℝ") for p, t in st.params},
                         {leandsl.lean_ident(p.name): p for p, _ in st.params})
        if sp.srepr(back) != sp.srepr(st.concl):
            return Read(d, None, "unsupported", "printing and reading it again gives a different statement")
    except (Unsupported, leandsl.Unprintable) as ex:
        return Read(d, None, "unsupported", f"cannot be printed back: {ex}")
    if not sample:
        return Read(d, st, "unchecked", "not sampled")
    status, detail = check_statement(d, st)
    return Read(d, st if status != "refuted" else None, status, detail)


# ---- into categories ---------------------------------------------------------------------------------
_REF_TOKEN = re.compile(r"[A-Za-z_][A-Za-z0-9_'.]*")


def read_module(text: str, module: str = "", sample: bool = True) -> list[Read]:
    """Every theorem/lemma of one Lean file, read (or with the reason it was not)."""
    return [read_declaration(d, sample) for d in declarations(text, module)]


def evidence_text(r: Read) -> str:
    n = re.search(r"held at (\d+) points", r.detail)
    if r.status == "sampled" and n:
        return (f"Read from the Lean statement below; sympy found it true at {n.group(1)} sampled points that satisfy "
                "the hypotheses (evidence, not a proof; the proof is the release's).")
    return "Read from the Lean statement below; not sampled (" + (r.detail or "no point satisfying the hypotheses was found") + ")."


def ingest(cls, reads: list[Read]) -> dict:
    """Declare the readable statements of ``reads`` as items of the category ``cls``; ``{qualified Lean name: item}``."""
    made = {}
    for r in reads:
        if r.stmt is None:
            continue
        d, st = r.decl, r.stmt
        reg = cls.lemmas if d.kind == "lemma" else cls.theorems
        name = d.qualified
        k = 2
        while name in reg._items:
            name = f"{d.qualified} ({k})"
            k += 1
        made[d.qualified] = reg.ingest(name, st.params, st.hyps, st.concl, lean=d.text, source=d.qualified,
                                       evidence=evidence_text(r))
    return made


def link_uses(items: dict, bodies: dict, index: dict) -> int:
    """``uses`` links from names: the proof of A mentions the name of B. ``items``: qualified name -> Item (what may be
    linked); ``bodies``: qualified name -> proof text; ``index``: short name -> every qualified name that has it
    (a short name that is not unique is not followed). Returns the number of links made."""
    short = {}
    for q in items:
        s = q.split(".")[-1]
        if len(index.get(s, ())) == 1:
            short[s] = q
    n = 0
    for q, it in items.items():
        seen = set()
        for tok in _REF_TOKEN.findall(bodies.get(q, "")):
            s = tok.split(".")[-1]
            tgt = short.get(s)
            if tgt and tgt != q and tgt not in seen and s != q.split(".")[-1]:
                seen.add(tgt)
                it.connect(items[tgt], kind="uses", note="its proof in the Lean file mentions the name")
                n += 1
    return n


# ---- command line and selftest -----------------------------------------------------------------------
def family_reads(cat, num: str, sample: bool = True, modules=None) -> list[Read]:
    """Everything the family's Lean modules declare, read or not (``modules`` restricts which modules)."""
    out = []
    root = cat.root / "lean"
    for m in sorted(modules if modules is not None else cat.family_closure(num)):
        p = root / (m.replace(".", "/") + ".lean")
        if p.exists():
            out += read_module(p.read_text(errors="replace"), m, sample)
    return out


def main(argv=None):
    import argparse
    from oai_catalog import Catalog, ROOT
    ap = argparse.ArgumentParser(description="Read the arithmetic theorems behind a family into sympy")
    ap.add_argument("--family")
    ap.add_argument("--limit", type=int, default=8, help="how many statements to print")
    ap.add_argument("--no-sample", action="store_true")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if not a.family:
        ap.error("--family is required")
    cat = Catalog(ROOT)
    num = f"{int(a.family):03d}"
    rs = family_reads(cat, num, sample=not a.no_sample)
    c = collections.Counter(r.status for r in rs)
    print(f"family {num}: {len(rs)} theorems/lemmas with a plain statement; read: {c['sampled'] + c['unchecked']} "
          f"(sampled {c['sampled']}, not sampled {c['unchecked']}), refuted {c['refuted']}, not read {c['unsupported']}")
    why = collections.Counter(r.detail.split(":")[0][:60] for r in rs if r.status == "unsupported")
    for k, v in why.most_common(6):
        print(f"  not read ({v}): {k}")
    for r in [r for r in rs if r.stmt is not None][:a.limit]:
        print(f"  {r.decl.qualified}: {sp.sstr(r.stmt.concl)}   [{r.status}]")


def selftest():
    n = [0]

    def ok(c, m=""):
        n[0] += 1
        assert c, m

    src = """
import Mathlib
namespace OAI
namespace Foo
/-- a comment with theorem fake : x := y -/
theorem tail_small {E T t : ℝ} {m : ℕ} (hE : 0 ≤ E) (ht : 0 ≤ t) (htT : t ≤ T) (hT : E*T^(m+1) ≤ 1) :
    E*t^(m+1) ≤ 1 := by
  -- uses nothing
  nlinarith [pow_le_pow_left ht htT (m+1)]
lemma abs_four {a b c d E : ℝ} (h1 : |a| ≤ E) (h2 : |b| ≤ E) : |a + b| ≤ 2 * E := by
  exact (abs_add _ _).trans (by linarith)
end Foo
theorem Foo.cast_inv {k : ℕ} (hk : 0 < k) : 0 ≤ 1 - (k:ℝ)⁻¹ := by
  have := tail_small
  sorry
theorem uses_def {x : ℝ} : myDef x ≤ 1 := by simp
theorem nat_sub {a b : ℕ} : a - b ≤ a := by omega
theorem nat_div {a b : ℕ} : a / b ≤ a := by omega
theorem mixed {a b : ℕ} {x : ℝ} (h : (a:ℝ) ≤ x) : (a:ℝ) - b ≤ x := by sorry
theorem quant {x : ℝ} : ∀ y : ℝ, x ≤ y := by sorry
theorem conj {x : ℝ} (h : 1 ≤ x) : 0 ≤ x ∧ 1 ≤ x := by sorry
theorem wrong {x y : ℝ} : x ≤ y := by sorry
theorem neg_pow {x : ℝ} : -x^2 ≤ 0 := by nlinarith
theorem sqrt_sq {x : ℝ} (hx : 0 ≤ x) : Real.sqrt (x^2) = x := by sorry
theorem fn_apps {x y : ℝ} (hx : 0 < x) (hy : 0 < y) : Real.log (x*y) = Real.log x + Real.log y := by sorry
theorem minmax {a b : ℝ} : min a b ≤ max a b := by sorry
theorem choose_ok {n : ℕ} : Nat.choose n 2 ≤ n^2 := by sorry
theorem prec {x : ℝ} (hx : 0 < x) : 2 * x ^ 2 ^ 1 + x⁻¹ * x = 2 * x ^ 2 + 1 := by sorry
end OAI
"""
    ds = declarations(src, "OAI.Demo")
    names = {d.qualified: d for d in ds}
    ok("OAI.Foo.tail_small" in names and "OAI.Foo.abs_four" in names and "OAI.Foo.cast_inv" in names and "OAI.fake" not in names
       and "fake" not in " ".join(names), list(names))
    ok(names["OAI.Foo.tail_small"].binders[0] == ("{", "E T t : ℝ") and names["OAI.Foo.tail_small"].concl.startswith("E*t^(m+1)"))
    ok("nlinarith" in names["OAI.Foo.tail_small"].body and "uses nothing" not in names["OAI.Foo.tail_small"].body)
    rs = {k.split(".")[-1]: read_declaration(d) for k, d in names.items()}
    ts = rs["tail_small"]
    ok(ts.status == "sampled" and len(ts.stmt.hyps) == 4 and len(ts.stmt.params) == 4, (ts.status, ts.detail))
    ok([p.name for p, _ in ts.stmt.params] == ["E", "T", "t", "m"] and ts.stmt.params[3][1] is Nat)
    ok(sp.srepr(ts.stmt.concl).startswith("LessThan") or isinstance(ts.stmt.concl, sp.Le))
    ok(rs["abs_four"].status == "sampled" and isinstance(rs["abs_four"].stmt.concl.lhs, sp.Abs), rs["abs_four"])
    ok(rs["cast_inv"].status in ("sampled", "unchecked") and rs["cast_inv"].stmt is not None, rs["cast_inv"])
    for k, why in (("uses_def", "unknown name"), ("nat_sub", "natural numbers"), ("nat_div", "integers"),
                   ("quant", "logic"), ("conj", "conjunction")):
        ok(rs[k].status == "unsupported" and why in rs[k].detail, (k, rs[k].detail))
    ok(rs["mixed"].status == "sampled", rs["mixed"].detail)            # a - b with a real tree: real subtraction, read
    ok(rs["wrong"].status == "refuted" and rs["wrong"].stmt is None and "false at" in rs["wrong"].detail, rs["wrong"])
    ok(rs["neg_pow"].status == "sampled" and rs["sqrt_sq"].status == "sampled" and rs["fn_apps"].status == "sampled")
    ok(rs["minmax"].status == "sampled" and rs["choose_ok"].status == "sampled", (rs["minmax"].detail, rs["choose_ok"].detail))
    x = sp.Symbol("x", real=True)
    c = rs["neg_pow"].stmt.concl
    ok(sp.simplify(c.lhs + rs["neg_pow"].stmt.params[0][0] ** 2) == 0, c)    # -x^2 is -(x^2)
    ok(rs["prec"].status == "sampled", rs["prec"].detail)                   # x ^ 2 ^ 1 associates to the right
    # the sampler can tell a misreading from a reading: a statement with its inequality reversed is refuted
    bad = Decl("M", "theorem", "rev", "rev", [("{", "x y : ℝ"), ("(", "h : x ≤ y")], "y ≤ x", "", "")
    ok(read_declaration(bad).status == "refuted")
    # into a category, with `uses` from names
    from leandsl import lean, CATEGORIES, REGISTRY
    saved_cats, saved_reg = list(CATEGORIES), list(REGISTRY)
    try:
        @lean.category("m")
        class Mod:
            pass
        made = ingest(Mod, list(rs.values()))
        ok("OAI.Foo.tail_small" in made and "OAI.wrong" not in made and "OAI.uses_def" not in made, list(made))
        it = made["OAI.Foo.tail_small"]
        ok(it.kind == "theorem" and made["OAI.Foo.abs_four"].kind == "lemma")
        ok(it.source == "OAI.Foo.tail_small" and it.lean_stmt.startswith("theorem tail_small") and "sampled points" in it.evidence)
        bodies = {q: d.body for q, d in names.items()}
        index = collections.defaultdict(list)
        for q in names:
            index[q.split(".")[-1]].append(q)
        nl = link_uses(made, bodies, index)
        ok(nl == 1 and [l.dst.name for l in made["OAI.Foo.cast_inv"].links] == ["OAI.Foo.tail_small"], (nl, list(made)))
        t = leandsl.category_tex(Mod)
        ok("sampled points" in t and "\\begin{verbatim}\ntheorem tail_small" in t and "\\hyperref[cat:Mod:thm:OAI.Foo.tail_small]" in t, t[:600])
        # a duplicate name gets a suffix
        again = ingest(Mod, [rs["tail_small"]])
        ok(list(again.values())[0].name == "OAI.Foo.tail_small (2)")
    finally:
        CATEGORIES[:] = saved_cats
        REGISTRY[:] = saved_reg
    print(f"oai_lean: {n[0]} checks ok")


if __name__ == "__main__":
    sys.modules.setdefault("oai_lean", sys.modules["__main__"])
    main()
