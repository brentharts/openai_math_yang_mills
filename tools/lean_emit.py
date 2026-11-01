#!/usr/bin/env python3
r"""Trees from ``texmath.py``, written as Lean 4 source.

Two flavors, because the fork's Lean project and a bare toolchain can answer different
questions:

``mathlib``  ``import Mathlib``; ``ℝ`` is Mathlib's. This is the flavor that belongs next to
             ``lean/OAI`` (the project pins Mathlib), and the only one in which Lean can
             actually *prove* anything (``ring1`` and ``field_simp`` close true field
             identities). It needs a built Mathlib to check.
``core``     No imports. ``Real`` and its operations are declared as axioms, as RosettaMath's
             ``rosettalean.py`` does, so a bare ``lean`` can check that every statement parses,
             elaborates and is a ``Prop``. It proves nothing: the operations have no laws.

What is and is not claimed. Emitting a statement says "this is what the LaTeX says, in Lean".
Whether Lean accepts it is ``run_lean``'s answer, and ``run_lean`` returns None when there is no
Lean to ask, so nothing is then claimed. A theorem whose body is ``sorry`` is a conjecture and
Lean says so with a warning; ``theorem_status`` reads those warnings back so a report can say
which theorems Lean proved and which it merely accepted as statements.

Two facts about Lean's reals that change what a statement means, both written into the output:

* ``x / 0 = 0`` in Lean. A statement with a denominator is true there for reasons that need not
  hold on paper, so a symbol's sympy assumptions (``positive=True``, ``nonzero=True`` ...)
  become hypotheses of the theorem rather than being dropped.
* ``x ^ (y : ℝ)`` is ``Real.rpow``, which is only the familiar power for ``x > 0``. Integer
  exponents are printed as natural or integer powers, which are the familiar ones.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from dataclasses import dataclass, field

from texmath import Call, Const, Num, Op, Rel, Sym, Unreadable

FLAVORS = ("mathlib", "core")

LEAN_RESERVED = {
    "fun", "let", "in", "do", "if", "then", "else", "match", "with", "have", "show", "from", "by", "at",
    "end", "this", "Type", "Prop", "Sort", "at", "open", "def", "theorem", "forall", "exists", "section",
    "namespace", "variable", "universe", "import", "for", "return", "mut", "unless", "calc", "suffices",
    "obtain", "deriving", "inductive", "example", "abbrev", "axiom", "instance", "structure", "class",
    "where", "mutual", "partial", "unsafe", "macro", "syntax", "notation", "local", "private",
    "protected", "noncomputable", "termination_by", "decreasing_by", "try", "catch", "finally", "break",
    "continue", "rfl", "sorry", "only", "ℝ", "ℕ", "ℤ", "Real", "Nat", "Int", "exp", "log", "abs", "max",
    "min", "sin", "cos", "tan", "sqrt", "pi", "id", "e",
}

PREC_ADD, PREC_MUL, PREC_POW, PREC_APP, PREC_ATOM = 65, 70, 75, 90, 100

REAL_FUNCS = {"sin": "Real.sin", "cos": "Real.cos", "tan": "Real.tan", "sinh": "Real.sinh",
              "cosh": "Real.cosh", "tanh": "Real.tanh", "arcsin": "Real.arcsin",
              "arccos": "Real.arccos", "arctan": "Real.arctan", "exp": "Real.exp",
              "log": "Real.log", "sqrt": "Real.sqrt"}

REL_SYMBOLS = {"eq": "=", "ne": "≠", "lt": "<", "gt": ">", "le": "≤", "ge": "≥"}


def ident(name: str) -> str:
    """A Lean-safe identifier for a canonical symbol name. Injective on the names we see."""
    out = re.sub(r"[^A-Za-z0-9_]", "_", name)
    if not out or out[0].isdigit() or out[0] == "_":
        out = "v" + out
    if out in LEAN_RESERVED:
        out += "_"
    return out


@dataclass
class Binder:
    """A quantity a theorem quantifies over: its type and what sympy assumed about it."""
    name: str                       # canonical (texmath) name
    kind: str = "real"              # real | int | nat
    assume: tuple = ()              # subset of: positive nonnegative negative nonpositive nonzero

    @property
    def lean(self) -> str:
        return ident(self.name)

    @property
    def type(self) -> str:
        return {"real": "ℝ", "int": "ℤ", "nat": "ℕ"}[self.kind]

    def hypotheses(self) -> list[tuple[str, str]]:
        v = self.lean
        table = {"positive": f"0 < {v}", "nonnegative": f"0 ≤ {v}", "negative": f"{v} < 0",
                 "nonpositive": f"{v} ≤ 0", "nonzero": f"{v} ≠ 0"}
        return [(f"h_{v}_{a}" if len(self.assume) > 1 else f"h_{v}", table[a]) for a in self.assume]


class Emitter:
    """One statement's worth of printing, remembering which core features it used."""

    def __init__(self, flavor: str, binders: dict[str, Binder]):
        if flavor not in FLAVORS:
            raise ValueError(f"flavor must be one of {FLAVORS}")
        self.flavor = flavor
        self.binders = binders
        self.features: set[str] = set()      # core only: zpow rpow sci natcast intcast abs

    # -- numbers and symbols ---------------------------------------------
    def num(self, n: Num) -> str:
        if "." in n.text:
            self.features.add("sci")
        return n.text

    def sym(self, n: Sym) -> str:
        b = self.binders.get(n.name)
        if b is None:
            raise Unreadable(f"symbol {n.name} has no binder")
        if b.kind == "real":
            return b.lean
        self.features.add("natcast" if b.kind == "nat" else "intcast")
        return f"({b.lean} : ℝ)"

    # -- expressions -----------------------------------------------------
    def go(self, node, need: int = 0) -> str:
        text, prec = self._go(node)
        return f"({text})" if prec < need else text

    def _go(self, node) -> tuple[str, int]:
        if isinstance(node, Num):
            return self.num(node), PREC_ATOM
        if isinstance(node, Sym):
            return self.sym(node), PREC_ATOM
        if isinstance(node, Const):
            return ("Real.pi", PREC_ATOM) if node.name == "pi" else ("Real.exp 1", PREC_APP)
        if isinstance(node, Call):
            return self.call(node)
        if isinstance(node, Op):
            return self.op(node)
        raise Unreadable(f"cannot print {node!r}")

    def call(self, node: Call) -> tuple[str, int]:
        if node.func == "abs":
            if self.flavor == "mathlib":
                return f"|{self.go(node.arg)}|", PREC_ATOM
            self.features.add("abs")
            return f"Real.abs {self.go(node.arg, PREC_ATOM)}", PREC_APP
        return f"{REAL_FUNCS[node.func]} {self.go(node.arg, PREC_ATOM)}", PREC_APP

    def op(self, node: Op) -> tuple[str, int]:
        a = node.args
        if node.op in ("add", "sub"):
            sign = "+" if node.op == "add" else "-"
            return f"{self.go(a[0], PREC_ADD)} {sign} {self.go(a[1], PREC_ADD + 1)}", PREC_ADD
        if node.op in ("mul", "div"):
            sign = "*" if node.op == "mul" else "/"
            return f"{self.go(a[0], PREC_MUL)} {sign} {self.go(a[1], PREC_MUL + 1)}", PREC_MUL
        if node.op == "neg":
            # always bracket the operand, and bracket the negation wherever it is not leftmost:
            # Lean's precedence for prefix minus is easy to get subtly wrong
            return f"-{self.go(a[0], PREC_ATOM)}", PREC_ADD
        if node.op == "pow":
            if isinstance(a[0], Const) and a[0].name == "e":
                # e^x is Real.exp x for every real x (Real.exp_one_rpow); this reads far better
                # than Lean's real power of Real.exp 1, and avoids rpow for the commonest function
                return f"Real.exp {self.go(a[1], PREC_ATOM)}", PREC_APP
            return f"{self.go(a[0], PREC_POW + 1)} ^ {self.exponent(a[1])}", PREC_POW
        raise Unreadable(f"cannot print operator {node.op}")

    def exponent(self, e) -> str:
        """The exponent of a power, typed so Lean picks the familiar power where there is one."""
        if isinstance(e, Num) and "." not in e.text:
            return e.text if self.flavor == "mathlib" else f"({e.text} : Nat)"
        if isinstance(e, Op) and e.op == "neg" and isinstance(e.args[0], Num) and "." not in e.args[0].text:
            self.features.add("zpow")
            return f"(-{e.args[0].text} : ℤ)"
        if isinstance(e, Sym):
            b = self.binders.get(e.name)
            if b is not None and b.kind == "nat":
                return b.lean if self.flavor == "mathlib" else f"({b.lean} : Nat)"
            if b is not None and b.kind == "int":
                self.features.add("zpow")
                return b.lean
        self.features.add("rpow")
        return f"({self.go(e)} : ℝ)"

    # -- statements ------------------------------------------------------
    def statement(self, rel) -> str:
        if not isinstance(rel, Rel):
            raise Unreadable("an expression is not a statement: there is no relation to state")
        # the ascription pins the type of an all-numeral left side to ℝ, not ℕ
        return f"({self.go(rel.left)} : ℝ) {REL_SYMBOLS[rel.op]} {self.go(rel.right)}"


# ----------------------------------------------------------------------------
# Theorems and files
# ----------------------------------------------------------------------------
@dataclass
class Theorem:
    name: str
    statement: str
    binders: list[Binder]
    doc: str = ""
    tactic: str = "sorry"
    features: set = field(default_factory=set)

    def text(self) -> str:
        groups: dict[str, list[str]] = {}
        for b in self.binders:
            groups.setdefault(b.type, []).append(b.lean)
        binder_text = " ".join(f"({' '.join(names)} : {t})" for t, names in groups.items())
        hyps = [f"({h} : {p})" for b in self.binders for h, p in b.hypotheses()]
        head = " ".join(x for x in [f"theorem {self.name}", binder_text, *hyps] if x)
        out = []
        if self.doc:
            # Lean block comments nest: a stray /- inside a doc comment would open one that never closes
            safe = self.doc.replace("/-", "/ -").replace("-/", "- /").strip()
            out.append("/-- " + safe + "\n-/")
        out.append(f"{head} :\n    {self.statement} := by\n  {self.tactic}")
        return "\n".join(out)


def theorem_name(key: str) -> str:
    name = "eq_" + re.sub(r"[^A-Za-z0-9]+", "_", key).strip("_")
    return name


def build_theorem(key: str, rel, binders: list[Binder], flavor: str, doc: str = "",
                  tactic: str = "sorry") -> Theorem:
    """A statement tree as a theorem. Raises Unreadable if it cannot be printed."""
    table = {b.name: b for b in binders}
    em = Emitter(flavor, table)
    stmt = em.statement(rel)
    used = [b for b in binders if _mentions(rel, b.name)]
    return Theorem(theorem_name(key), stmt, used, doc, tactic, em.features)


def _mentions(node, name: str) -> bool:
    from texmath import symbols
    return name in symbols(node)


MATHLIB_HEADER = """import Mathlib

/-!
# {title}

{blurb}
-/

namespace OAI.SymModels

"""

CORE_HEADER = """/-!
# {title}

{blurb}

Flavor `core`: `Real` and its operations are declared as axioms, not defined. This file asks
whether the statements below are well-formed propositions. Proving them needs the field axioms
and the physics, and neither is claimed here.
-/

namespace SymModels

noncomputable section

axiom Real : Type
axiom Real.add : Real -> Real -> Real
axiom Real.sub : Real -> Real -> Real
axiom Real.mul : Real -> Real -> Real
axiom Real.div : Real -> Real -> Real
axiom Real.neg : Real -> Real
axiom Real.sqrt : Real -> Real
axiom Real.log : Real -> Real
axiom Real.exp : Real -> Real
axiom Real.sin : Real -> Real
axiom Real.cos : Real -> Real
axiom Real.tan : Real -> Real
axiom Real.sinh : Real -> Real
axiom Real.cosh : Real -> Real
axiom Real.tanh : Real -> Real
axiom Real.arcsin : Real -> Real
axiom Real.arccos : Real -> Real
axiom Real.arctan : Real -> Real
axiom Real.pi : Real
axiom Real.pow : Real -> Nat -> Real
axiom Real.lit : Nat -> Real
axiom Real.le : Real -> Real -> Prop
axiom Real.lt : Real -> Real -> Prop

instance : Add Real := Add.mk Real.add
instance : Sub Real := Sub.mk Real.sub
instance : Mul Real := Mul.mk Real.mul
instance : Div Real := Div.mk Real.div
instance : Neg Real := Neg.mk Real.neg
instance : LE Real := LE.mk Real.le
instance : LT Real := LT.mk Real.lt
instance : HPow Real Nat Real := HPow.mk Real.pow
instance (n : Nat) : OfNat Real n := OfNat.mk (Real.lit n)
notation "ℝ" => Real
"""

CORE_FEATURES = {
    "zpow": ["axiom Real.zpow : Real -> Int -> Real",
             "instance : HPow Real Int Real := HPow.mk Real.zpow"],
    "rpow": ["axiom Real.rpow : Real -> Real -> Real",
             "instance : HPow Real Real Real := HPow.mk Real.rpow"],
    "sci": ["axiom Real.sci : Nat -> Bool -> Nat -> Real",
            "instance : OfScientific Real := OfScientific.mk Real.sci"],
    "natcast": ["instance : NatCast Real := NatCast.mk Real.lit"],
    "intcast": ["axiom Real.ofInt : Int -> Real",
                "instance : IntCast Real := IntCast.mk Real.ofInt"],
    "abs": ["axiom Real.abs : Real -> Real"],
}
# `ℕ` and `ℤ` in core binders
CORE_NOTATION = ['notation "ℕ" => Nat', 'notation "ℤ" => Int']


def lean_file(theorems: list[Theorem], flavor: str, title: str, blurb: str) -> str:
    """A complete .lean file: header, (core: only the features used), then the theorems."""
    if flavor == "mathlib":
        parts = [MATHLIB_HEADER.format(title=title, blurb=blurb)]
        closing = "\nend OAI.SymModels\n"
    else:
        parts = [CORE_HEADER.format(title=title, blurb=blurb)]
        used = set().union(*(t.features for t in theorems)) if theorems else set()
        parts.append("\n".join(CORE_NOTATION))
        for feature in ("zpow", "rpow", "sci", "natcast", "intcast", "abs"):
            if feature in used:
                parts.append("\n".join(CORE_FEATURES[feature]))
        parts.append("")
        closing = "\nend\n\nend SymModels\n"
    parts.extend(t.text() + "\n" for t in theorems)
    parts.append(closing)
    return "\n".join(parts)


# ----------------------------------------------------------------------------
# Asking Lean
# ----------------------------------------------------------------------------
def run_lean(path: str, flavor: str, project: str | None = None, timeout: int = 1800) -> dict | None:
    """Give the file to Lean. None when there is no Lean to ask (nothing is then claimed).

    A dict otherwise: ``ok`` (exit 0 and no errors), ``version``, ``errors`` (list of lines),
    ``output`` (everything Lean printed), or ``skipped`` with the reason when the flavor
    cannot be checked here (``mathlib`` needs ``lake`` and a built Mathlib in ``project``).
    """
    lean = shutil.which("lean")
    if flavor == "mathlib":
        lake = shutil.which("lake")
        if lake is None:
            return None if lean is None else {"skipped": "lake not found; mathlib flavor needs the lean/ project"}
        if project is None or not os.path.isdir(os.path.join(project, ".lake", "packages", "mathlib")):
            return {"skipped": "Mathlib is not built in the lean/ project "
                               "(run `lake exe cache get` there, then retry)"}
        cmd, cwd = [lake, "env", "lean", os.path.abspath(path)], project
    else:
        if lean is None:
            return None
        cmd, cwd = [lean, os.path.abspath(path)], None
    try:
        run = subprocess.run(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True, timeout=timeout)
        version = subprocess.run([lean or "lean", "--version"], stdout=subprocess.PIPE,
                                 universal_newlines=True).stdout
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"skipped": f"could not run Lean: {exc}"}
    found = re.search(r"version ([0-9][^,\s]*)", version)
    errors = [ln for ln in run.stdout.splitlines() if re.search(r":\s*error\b", ln)]
    return {"ok": run.returncode == 0 and not errors, "version": found.group(1) if found else "?",
            "errors": errors, "output": run.stdout}


def theorem_status(source: str, output: str) -> dict[str, str]:
    """For each theorem: ``proved`` (Lean accepted it with no sorry), ``statement`` (accepted, body is
    sorry), or ``error``. Reads Lean's own messages back against the theorem's line."""
    starts = [(i + 1, m.group(1)) for i, ln in enumerate(source.splitlines())
              if (m := re.match(r"theorem\s+(\S+)", ln))]

    def owner(line: int):
        names = [n for (s, n) in starts if s <= line]
        return names[-1] if names else None

    status = {n: "proved" for _, n in starts}
    for ln in output.splitlines():
        m = re.match(r".*?:(\d+):\d+:\s*(error|warning):\s*(.*)", ln)
        if not m:
            continue
        who = owner(int(m.group(1)))
        if who is None:
            continue
        if m.group(2) == "error":
            status[who] = "error"
        elif "sorry" in m.group(3) and status[who] != "error":
            status[who] = "statement"
    return status


# ----------------------------------------------------------------------------
# Self test
# ----------------------------------------------------------------------------
def selftest() -> int:
    from texmath import parse
    n = 0

    def ok(cond, msg):
        nonlocal n
        if not cond:
            raise SystemExit(f"lean_emit selftest FAILED: {msg}")
        n += 1

    def stmt(src, flavor="mathlib", **kinds):
        tree = parse(src)
        from texmath import symbols
        binders = {s: Binder(s, kinds.get(s, "real")) for s in symbols(tree)}
        return Emitter(flavor, binders).statement(tree)

    ok(stmt(r"y = \frac{1}{2} m v^{2}") == "(y : ℝ) = 1 / 2 * m * v ^ 2", stmt(r"y = \frac{1}{2} m v^{2}"))
    ok(stmt(r"a = b - (c - d)") == "(a : ℝ) = b - (c - d)", "right operand of subtraction is bracketed")
    ok(stmt(r"a = (b + c) d") == "(a : ℝ) = (b + c) * d", "sum inside product")
    ok(stmt(r"a = -b^{2}") == "(a : ℝ) = -(b ^ 2)", stmt(r"a = -b^{2}"))
    ok(stmt(r"a = c \cdot (-b)") == "(a : ℝ) = c * (-b)", stmt(r"a = c \cdot (-b)"))
    ok(stmt(r"a = \sqrt{x + 1}") == "(a : ℝ) = Real.sqrt (x + 1)", "sqrt")
    ok(stmt(r"a = \sin^{2}{\left(x \right)}") == "(a : ℝ) = Real.sin x ^ 2", stmt(r"a = \sin^{2}{\left(x \right)}"))
    ok(stmt(r"a = x^{-2}") == "(a : ℝ) = x ^ (-2 : ℤ)", stmt(r"a = x^{-2}"))
    ok(stmt(r"a = x^{\frac{3}{2}}") == "(a : ℝ) = x ^ (3 / 2 : ℝ)", stmt(r"a = x^{\frac{3}{2}}"))
    ok(stmt(r"a = x^{n}", n="nat") == "(a : ℝ) = x ^ n", "natural exponent stays natural")
    ok(stmt(r"a = n x", n="nat") == "(a : ℝ) = (n : ℝ) * x", "natural symbol cast in a real position")
    ok(stmt(r"a = \mathrm{e}^{x}") == "(a : ℝ) = Real.exp x", stmt(r"a = \mathrm{e}^{x}"))
    ok(stmt(r"a = \mathrm{e}^{-x^{2}}") == "(a : ℝ) = Real.exp (-(x ^ 2))", stmt(r"a = \mathrm{e}^{-x^{2}}"))
    ok(stmt(r"a = \mathrm{e}") == "(a : ℝ) = Real.exp 1", "the constant e")
    ok(stmt(r"2 + 2 = 4") == "(2 + 2 : ℝ) = 4", stmt(r"2 + 2 = 4"))
    ok(stmt(r"x \leq y") == "(x : ℝ) ≤ y", "inequality")
    ok(stmt(r"a = x^{2}", "core") == "(a : ℝ) = x ^ (2 : Nat)", "core flavor types its exponents")
    try:
        Emitter("mathlib", {}).statement(parse(r"x + 1"))
    except Unreadable as exc:
        ok("relation" in str(exc), "an expression is not a statement")
    else:
        ok(False, "an expression must not become a statement")

    # a whole theorem, with hypotheses
    tree = parse(r"\omega = \frac{\hbar}{m}")
    thm = build_theorem("magnon.disp", tree, [Binder("omega"), Binder("hbar", assume=("positive",)),
                                              Binder("m", assume=("positive",))], "mathlib", doc="**Test**")
    text = thm.text()
    ok("theorem eq_magnon_disp (omega hbar m : ℝ) (h_hbar : 0 < hbar) (h_m : 0 < m) :" in text, text)
    ok(text.rstrip().endswith("sorry"), "default body is sorry")
    ok(text.startswith("/-- **Test**\n-/"), "doc comment first")
    nested = build_theorem("n", tree, [Binder("omega"), Binder("hbar"), Binder("m")], "mathlib",
                           doc="a/-b and c-/d").text()
    ok(nested.count("/-") == 1 and nested.count("-/") == 1, "a stray comment delimiter in a doc cannot nest: " + nested)

    # names
    ok(ident("k_B") == "k_B" and ident("exp") == "exp_" and ident("1x") == "v1x", "identifiers")
    ok(theorem_name("bloch.A1") == "eq_bloch_A1", "theorem name")

    # reading Lean's messages back
    src = "theorem a : True := by\n  sorry\n\ntheorem b : True := by\n  trivial\n\ntheorem c : 1 = 2 := by\n  rfl\n"
    out = "f.lean:1:8: warning: declaration uses 'sorry'\nf.lean:7:0: error: bad\n"
    ok(theorem_status(src, out) == {"a": "statement", "b": "proved", "c": "error"}, theorem_status(src, out))

    # a core file declares only the features its theorems used
    t = build_theorem("x", parse(r"a = x^{\frac{3}{2}}"), [Binder("a"), Binder("x")], "core")
    f = lean_file([t], "core", "T", "B")
    ok("axiom Real.rpow" in f and "axiom Real.zpow" not in f and "OfScientific" not in f, "core features on demand")
    ok(f.count("(") == f.count(")"), "balanced brackets in the core file")
    return n


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        print(f"lean_emit: {selftest()} checks ok")
