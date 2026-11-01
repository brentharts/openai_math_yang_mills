#!/usr/bin/env python3
r"""The physics subset of LaTeX, read into a small tree (and back into sympy).

This is the middle of the sympy -> LaTeX -> Lean pipeline (see ``sym2lean.py``). It reads
what a model's equations look like on the page: fractions, powers, roots, subscripted
quantities, Greek letters, the usual functions, and one relation per statement. The tree it
builds is what ``lean_emit.py`` prints as Lean 4, and what ``to_sympy`` turns back into a
sympy expression so the pipeline can check the LaTeX still says what the model said.

The reader is deliberately small, and it refuses rather than guesses. Each of these is
refused with a message that names it, because each has a reading that is silently wrong:

    f(x)         a symbol applied to brackets: a function, or a product?  (sympy's own
                 ``f{\left(x\right)}`` is recognised and refused as a function)
    \sin 2 x     an unbracketed argument followed by more factors
    x 2          juxtaposed numbers
    x^23         TeX reads only the first digit of a bare exponent
    e^{x}        Euler's constant or a symbol named e? write \mathrm{e} or \exp
    \sin^{-1}    reciprocal or inverse?
    \sum \int \partial \nabla \lim ...   binders and operators outside the subset
    \approx \sim \propto                 not equalities, so no exact Lean statement

The notation follows what sympy's printer emits, which ``sym2lean.TexPrinter`` tightens:
``\mathrm{name}`` is one multi-letter symbol, ``\mathrm{e}`` is Euler's constant, ``\exp``
is the exponential. Hand-written LaTeX in that notation reads the same way.

This module is stdlib-only except ``to_sympy``, which imports sympy when called.
"""

from __future__ import annotations

import re
from dataclasses import dataclass


class Unreadable(Exception):
    """This corner of the notation has no reading in the subset. The message names it."""


# ----------------------------------------------------------------------------
# The tree
# ----------------------------------------------------------------------------
@dataclass(frozen=True)
class Num:
    text: str                      # '2', '0.5'


@dataclass(frozen=True)
class Sym:
    name: str                      # canonical: 'k_B', 'hbar', 'mu_0', 'Delta'


@dataclass(frozen=True)
class Const:
    name: str                      # 'pi' or 'e'


@dataclass(frozen=True)
class Call:
    func: str                      # a key of FUNCS (canonical name, ln -> log)
    arg: object


@dataclass(frozen=True)
class Op:
    op: str                        # add sub mul div pow neg
    args: tuple


@dataclass(frozen=True)
class Rel:
    op: str                        # eq ne lt gt le ge
    left: object
    right: object


Node = (Num, Sym, Const, Call, Op, Rel)

# canonical function names the subset knows
FUNCS = ("sin", "cos", "tan", "sinh", "cosh", "tanh", "arcsin", "arccos", "arctan",
         "exp", "log", "sqrt", "abs")
FUNC_ALIASES = {"ln": "log", "asin": "arcsin", "acos": "arccos", "atan": "arctan"}

GREEK = ("alpha beta gamma delta epsilon varepsilon zeta eta theta vartheta iota kappa lambda mu nu xi "
         "rho varrho sigma tau upsilon phi varphi chi psi omega "
         "Gamma Delta Theta Lambda Xi Pi Sigma Upsilon Phi Psi Omega hbar ell").split()
# \pi is a constant, not a quantity; \varpi is not supported

REFUSED_COMMANDS = {
    "sum": "a summation needs a binder", "prod": "a product needs a binder",
    "int": "an integral needs a binder", "oint": "an integral needs a binder",
    "lim": "a limit needs a binder", "partial": "a derivative is not in the subset",
    "nabla": "a differential operator is not in the subset",
    "infty": "infinity is not a real number", "cdots": "an ellipsis hides a pattern",
    "ldots": "an ellipsis hides a pattern", "dots": "an ellipsis hides a pattern",
    "approx": "not an equality, so it has no exact Lean statement",
    "sim": "not an equality, so it has no exact Lean statement",
    "simeq": "not an equality, so it has no exact Lean statement",
    "propto": "not an equality, so it has no exact Lean statement",
    "equiv": "congruence or identity? no single reading",
    "to": "a limit or a map, not arithmetic", "rightarrow": "a limit or a map, not arithmetic",
    "pm": "two values at once", "mp": "two values at once",
    "binom": "binomial coefficients are not in the subset",
    "begin": "environments (matrices, cases, align) are not in the subset",
    "text": None, "mathrm": None, "mathit": None, "operatorname": None,   # handled, listed for clarity
}

SPACING = {"\\,", "\\;", "\\:", "\\!", "\\ ", "\\quad", "\\qquad", "\\displaystyle", "\\textstyle",
           "\\left", "\\right", "\\bigl", "\\bigr", "\\Bigl", "\\Bigr", "\\big", "\\Big", "%%BR%%"}

RELATIONS = {"=": "eq", "<": "lt", ">": "gt", "\\le": "le", "\\leq": "le", "\\ge": "ge",
             "\\geq": "ge", "\\neq": "ne", "\\ne": "ne"}

_TOKEN = re.compile(
    r"""\s+                              # whitespace (dropped)
      | (?P<br>%%BR%%)
      | (?P<cmd>\\[A-Za-z]+|\\[^A-Za-z])
      | (?P<num>\d+(?:\.\d+)?)
      | (?P<ch>[A-Za-z])
      | (?P<p>[{}()\[\]^_+\-*/=<>|,&'!])
    """, re.X)


def tokenize(src: str) -> list[str]:
    out, i = [], 0
    while i < len(src):
        m = _TOKEN.match(src, i)
        if m is None:
            raise Unreadable(f"unreadable character {src[i]!r} at position {i}")
        i = m.end()
        tok = m.group(0)
        if tok.strip() == "" or tok in SPACING:
            continue
        out.append(tok)
    return out


# ----------------------------------------------------------------------------
# The parser
# ----------------------------------------------------------------------------
class _Parser:
    def __init__(self, tokens: list[str], euler_e: bool = False):
        self.t = tokens
        self.i = 0
        self.euler_e = euler_e     # read a bare e as Euler's constant (text sympy printed itself)

    # -- plumbing --------------------------------------------------------
    def peek(self, k: int = 0):
        j = self.i + k
        return self.t[j] if j < len(self.t) else None

    def take(self):
        tok = self.peek()
        if tok is None:
            raise Unreadable("the statement ends where more was expected")
        self.i += 1
        return tok

    def expect(self, tok: str):
        got = self.take()
        if got != tok:
            raise Unreadable(f"expected {tok!r} but found {got!r}")

    def starts_atom(self, tok) -> bool:
        """Can this token begin an operand (so that juxtaposition means a product)?"""
        if tok is None:
            return False
        if tok.isalpha() and len(tok) == 1:
            return True
        if tok in ("(", "{", "["):
            return True
        if tok.startswith("\\"):
            name = tok[1:]
            return (name in GREEK or name in ("pi", "frac", "dfrac", "tfrac", "sqrt", "mathrm", "mathit",
                                              "text", "operatorname") or name in FUNCS
                    or name in FUNC_ALIASES)
        return False

    # -- statements ------------------------------------------------------
    def statement(self):
        left = self.expr()
        tok = self.peek()
        if tok is not None and (tok in RELATIONS or (tok.startswith("\\") and tok[1:] in REFUSED_COMMANDS)):
            if tok in RELATIONS:
                self.take()
                right = self.expr()
                node = Rel(RELATIONS[tok], left, right)
                if self.peek() is not None:
                    raise Unreadable(f"more after the relation: {self.peek()!r} (one relation per statement)")
                return node
            raise Unreadable(f"{tok}: {REFUSED_COMMANDS[tok[1:]]}")
        if tok == "'":
            raise Unreadable("a prime (derivative or a distinct symbol?)")
        if tok is not None:
            raise Unreadable(f"cannot read {tok!r}")
        return left

    # -- sums, products, powers -------------------------------------------
    def expr(self):
        node = self.term()
        while self.peek() in ("+", "-"):
            op = "add" if self.take() == "+" else "sub"
            node = Op(op, (node, self.term()))
        return node

    def term(self):
        # a leading sign applies to the whole product: -\pi\beta is -(pi beta), as it is on the page
        tok = self.peek()
        if tok in ("-", "+"):
            self.take()
            body = self.term_body(self.power())
            return Op("neg", (body,)) if tok == "-" else body
        return self.term_body(self.unary())

    def term_body(self, node):
        while True:
            tok = self.peek()
            if tok in ("*", "\\cdot", "\\times"):
                self.take()
                node = Op("mul", (node, self.unary()))
            elif tok == "/":
                self.take()
                node = Op("div", (node, self.unary()))
            elif tok is not None and tok[0].isdigit():
                raise Unreadable(f"the number {tok} follows a factor with no operator (juxtaposed numbers)")
            elif self.starts_atom(tok):
                node = Op("mul", (node, self.power()))
            else:
                return node

    def unary(self):
        tok = self.peek()
        if tok == "-":
            self.take()
            return Op("neg", (self.unary(),))
        if tok == "+":
            self.take()
            return self.unary()
        return self.power()

    def power(self):
        base = self.atom()
        if self.peek() == "^":
            self.take()
            base = Op("pow", (base, self.script()))
            if self.peek() in ("^", "_"):
                raise Unreadable("a second script on a power (x^a^b or x^a_b)")
        return base

    def script(self):
        """The argument of ^ or _ : a braced expression, or exactly one token."""
        tok = self.peek()
        if tok == "{":
            self.take()
            node = self.expr()
            self.expect("}")
            return node
        if tok is None:
            raise Unreadable("a script ends the statement")
        if tok[0].isdigit() and len(tok) > 1:
            raise Unreadable(f"bare exponent {tok}: TeX reads only the first digit; write ^{{{tok}}}")
        if tok == "-":
            self.take()
            return Op("neg", (self.script(),))
        return self.atom(bare_script=True)

    # -- atoms -----------------------------------------------------------
    def atom(self, bare_script: bool = False):
        tok = self.take()
        if tok[0].isdigit():
            return Num(tok)
        if tok.isalpha() and len(tok) == 1:
            if tok == "e" and self.euler_e:
                return Const("e")
            if tok == "e" and self.peek() == "^":
                raise Unreadable("e^... is ambiguous (Euler's constant or a symbol named e); "
                                 "write \\mathrm{e} or \\exp")
            return self.symbol(tok)
        if tok in ("(", "["):
            node = self.expr()
            self.expect(")" if tok == "(" else "]")
            return node
        if tok == "{":
            node = self.expr()
            self.expect("}")
            return node
        if tok == "|":
            node = self.expr()
            self.expect("|")
            return Call("abs", node)
        if tok == "'":
            raise Unreadable("a prime (derivative or a distinct symbol?)")
        if tok.startswith("\\"):
            return self.command(tok[1:])
        raise Unreadable(f"cannot read {tok!r}")

    def subscript(self) -> str:
        if self.peek() != "_":
            return ""
        self.take()
        tok = self.take()
        parts = []
        if tok == "{":
            while self.peek() != "}":
                p = self.take()
                if p in ("_", "^", "{"):
                    raise Unreadable("nested scripts inside a subscript")
                parts.append(self.subscript_piece(p))
            self.take()
        else:
            parts.append(self.subscript_piece(tok))
        return "_" + "".join(parts)

    def subscript_piece(self, p: str) -> str:
        if p in ("\\mathrm", "\\text", "\\mathit", "\\operatorname"):
            return self.raw_name()
        return self._subscript_piece(p)

    @staticmethod
    def _subscript_piece(p: str) -> str:
        if p.startswith("\\") and p[1:] in GREEK:
            return p[1:]
        if p.isalnum() and p.isascii():
            return p
        if p == ",":
            return ","
        raise Unreadable(f"cannot read {p!r} inside a subscript")

    def symbol(self, name: str):
        name += self.subscript()
        self._not_applied(name)
        if name in ("o", "O") and self.peek() in ("(", "["):
            raise Unreadable(f"{name}(...) reads as the asymptotic notation, which is not arithmetic "
                             "(and as a product it would be silently wrong)")
        return Sym(name)

    def _not_applied(self, name: str):
        # sympy writes an applied function as  f{\left(x \right)}
        if self.peek() == "{" and self.peek(1) in ("(", "\\left", "["):
            raise Unreadable(f"{name} applied to an argument is a function, not a product; "
                             "functions of symbols are not in the subset")

    def group(self) -> object:
        self.expect("{")
        node = self.expr()
        self.expect("}")
        return node

    def raw_name(self) -> str:
        """The text of \\mathrm{...} / \\text{...}."""
        self.expect("{")
        parts = []
        while self.peek() != "}":
            p = self.take()
            if not (p.isalnum() and p.isascii()):
                raise Unreadable(f"cannot read {p!r} inside a name")
            parts.append(p)
        self.take()
        if not parts:
            raise Unreadable("an empty name")
        return "".join(parts)

    def command(self, name: str):
        if name in ("frac", "dfrac", "tfrac"):
            if self.peek() == "{" and self.peek(1) in ("d", "\\partial"):
                raise Unreadable("a derivative written as a fraction is not in the subset")
            num = self.group()
            den = self.group()
            return Op("div", (num, den))
        if name == "sqrt":
            index = None
            if self.peek() == "[":
                self.take()
                index = self.expr()
                self.expect("]")
            arg = self.group()
            if index is None:
                return Call("sqrt", arg)
            return Op("pow", (arg, Op("div", (Num("1"), index))))
        if name == "pi":
            return Const("pi")
        if name in GREEK:
            return self.symbol(name)
        if name in ("mathrm", "mathit", "text", "operatorname"):
            text = self.raw_name()
            canon = FUNC_ALIASES.get(text, text)
            if name == "operatorname" and canon in FUNCS:
                return self.function(canon)
            if name == "mathrm" and text == "e":
                return Const("e")
            if name == "mathrm" and text in FUNCS:
                return self.function(text)
            return self.symbol(text)
        canon = FUNC_ALIASES.get(name, name)
        if canon in FUNCS:
            return self.function(canon)
        if name in REFUSED_COMMANDS and REFUSED_COMMANDS[name]:
            raise Unreadable(f"\\{name}: {REFUSED_COMMANDS[name]}")
        raise Unreadable(f"no reading for \\{name}")

    def function(self, func: str):
        exponent = None
        if self.peek() == "^":
            self.take()
            exponent = self.script()
            if exponent == Op("neg", (Num("1"),)):
                raise Unreadable(f"{func}^{{-1}} is ambiguous (reciprocal or inverse function)")
        tok = self.peek()
        if tok is None:
            raise Unreadable(f"\\{func} has no argument")
        if tok in ("(", "["):                      # \left( ... \right) arrives here with \left dropped
            self.take()
            arg = self.expr()
            self.expect(")" if tok == "(" else "]")
            bare = False
        elif tok == "{":
            arg = self.group()
            bare = False
        else:
            arg = self.atom()
            bare = True
        if bare and (self.starts_atom(self.peek()) or (self.peek() or "")[:1].isdigit()):
            raise Unreadable(f"\\{func} {self._describe(arg)} is followed by more factors; "
                             "bracket the argument")
        node = Call(func, arg)
        return Op("pow", (node, exponent)) if exponent is not None else node

    @staticmethod
    def _describe(node) -> str:
        return node.name if isinstance(node, Sym) else (node.text if isinstance(node, Num) else "...")


def parse(latex: str, euler_e: bool = False):
    """A LaTeX statement or expression as a tree. Raises Unreadable, naming what it met.

    ``euler_e=True`` reads a bare ``e`` as Euler's constant. Use it only for text that sympy's
    printer wrote (which prints a symbol named e as ``\\mathit{e}``); hand-written LaTeX keeps the
    default, where ``e^{x}`` is refused as ambiguous.
    """
    # the fork's line-break marker must go first: its % would otherwise open a comment and
    # silently drop the rest of the equation
    cleaned = latex.replace("%%BR%%", " ")
    cleaned = re.sub(r"(?<!\\)%[^\n]*", "", cleaned)           # LaTeX comments
    cleaned = cleaned.replace("\\\\", " ").strip().rstrip(".,;")
    if not cleaned:
        raise Unreadable("an empty statement")
    p = _Parser(tokenize(cleaned), euler_e)
    return p.statement()


# ----------------------------------------------------------------------------
# Walking the tree
# ----------------------------------------------------------------------------
def symbols(node) -> list[str]:
    """Names of the symbols a tree mentions, in order of first appearance."""
    seen: list[str] = []

    def walk(n):
        if isinstance(n, Sym):
            if n.name not in seen:
                seen.append(n.name)
        elif isinstance(n, Call):
            walk(n.arg)
        elif isinstance(n, Op):
            for a in n.args:
                walk(a)
        elif isinstance(n, Rel):
            walk(n.left)
            walk(n.right)

    walk(node)
    return seen


def show(node) -> str:
    """A compact s-expression, for tests and error messages."""
    if isinstance(node, Num):
        return node.text
    if isinstance(node, Sym):
        return node.name
    if isinstance(node, Const):
        return node.name
    if isinstance(node, Call):
        return f"{node.func}({show(node.arg)})"
    if isinstance(node, (Op, Rel)):
        args = node.args if isinstance(node, Op) else (node.left, node.right)
        return f"({node.op} " + " ".join(show(a) for a in args) + ")"
    raise TypeError(node)


# ----------------------------------------------------------------------------
# Back to sympy
# ----------------------------------------------------------------------------
def to_sympy(node, symmap: dict | None = None):
    """The tree as a sympy expression.

    ``symmap`` maps symbol names to the sympy Symbols to use for them (normally the symbols of
    the original model, so their assumptions carry over). Names not in it become plain Symbols.
    """
    import sympy as sp

    symmap = dict(symmap or {})
    funcs = {"sin": sp.sin, "cos": sp.cos, "tan": sp.tan, "sinh": sp.sinh, "cosh": sp.cosh,
             "tanh": sp.tanh, "arcsin": sp.asin, "arccos": sp.acos, "arctan": sp.atan,
             "exp": sp.exp, "log": sp.log, "sqrt": sp.sqrt, "abs": sp.Abs}
    # not evaluated: a statement about a positive symbol that sympy can decide from its
    # assumptions must still come back as the statement that was written, not as True/False
    rels = {"eq": lambda a, b: sp.Eq(a, b, evaluate=False), "ne": lambda a, b: sp.Ne(a, b, evaluate=False),
            "lt": lambda a, b: sp.Lt(a, b, evaluate=False), "gt": lambda a, b: sp.Gt(a, b, evaluate=False),
            "le": lambda a, b: sp.Le(a, b, evaluate=False), "ge": lambda a, b: sp.Ge(a, b, evaluate=False)}

    def go(n):
        if isinstance(n, Num):
            return sp.Rational(n.text) if "." in n.text else sp.Integer(n.text)
        if isinstance(n, Sym):
            if n.name not in symmap:
                symmap[n.name] = sp.Symbol(n.name)
            return symmap[n.name]
        if isinstance(n, Const):
            return sp.pi if n.name == "pi" else sp.E
        if isinstance(n, Call):
            return funcs[n.func](go(n.arg))
        if isinstance(n, Op):
            a = [go(x) for x in n.args]
            return {"add": lambda: a[0] + a[1], "sub": lambda: a[0] - a[1], "mul": lambda: a[0] * a[1],
                    "div": lambda: a[0] / a[1], "pow": lambda: a[0] ** a[1], "neg": lambda: -a[0]}[n.op]()
        if isinstance(n, Rel):
            return rels[n.op](go(n.left), go(n.right))
        raise TypeError(n)

    return go(node)


# ----------------------------------------------------------------------------
# Self test
# ----------------------------------------------------------------------------
def selftest() -> int:
    n = 0

    def ok(cond, msg):
        nonlocal n
        if not cond:
            raise SystemExit(f"texmath selftest FAILED: {msg}")
        n += 1

    def reads(src, expected):
        got = show(parse(src))
        ok(got == expected, f"{src!r} -> {got}, expected {expected}")

    def refuses(src, needle):
        try:
            parse(src)
        except Unreadable as exc:
            ok(needle in str(exc), f"{src!r} refused with {exc!s}, expected to mention {needle!r}")
        else:
            ok(False, f"{src!r} should have been refused")

    reads(r"a + b c", "(add a (mul b c))")
    reads(r"a - b - c", "(sub (sub a b) c)")
    reads(r"a / b c", "(mul (div a b) c)")
    reads(r"-x^2", "(neg (pow x 2))")
    reads(r"-\pi \beta + c", "(add (neg (mul pi beta)) c)")
    reads(r"a - -b", "(sub a (neg b))")
    reads(r"a \cdot -b c", "(mul (mul a (neg b)) c)")
    reads(r"x^{2}", "(pow x 2)")
    reads(r"\frac{1}{2} m v^{2}", "(mul (mul (div 1 2) m) (pow v 2))")
    reads(r"k_{B} T", "(mul k_B T)")
    reads(r"\hbar \omega_{0}", "(mul hbar omega_0)")
    reads(r"\mathrm{Delta} + \Delta", "(add Delta Delta)")
    reads(r"\sqrt{x + 1}", "sqrt((add x 1))")
    reads(r"\sqrt[3]{x}", "(pow x (div 1 3))")
    reads(r"\sin^{2}{\left(x \right)}", "(pow sin(x) 2)")
    reads(r"\exp{\left(-x \right)}", "exp((neg x))")
    reads(r"\operatorname{asin}{\left(x \right)}", "arcsin(x)")
    reads(r"\left|{x - y}\right|", "abs((sub x y))")
    reads(r"\left(\frac{a}{b}\right)^{2}", "(pow (div a b) 2)")
    reads(r"2 \cdot 3", "(mul 2 3)")
    reads(r"1.0 \cdot 10^{-5} x", "(mul (mul 1.0 (pow 10 (neg 5))) x)")
    reads(r"\mathrm{e}^{x}", "(pow e x)")
    reads(r"\pi r^{2}", "(mul pi (pow r 2))")
    reads(r"y = \frac{1}{x}", "(eq y (div 1 x))")
    reads(r"a \leq b", "(le a b)")
    reads(r"a = b %%BR%% + c", "(eq a (add b c))")
    refuses(r"f{\left(x \right)}", "function")
    refuses(r"\sin 2 x", "bracket the argument")
    refuses(r"x 2", "juxtaposed")
    refuses(r"x^23", "first digit")
    refuses(r"e^{x}", "ambiguous")
    ok(show(parse(r"e^{-x} + \mathit{e}", euler_e=True)) == "(add (pow e (neg x)) e)", "euler_e reads a bare e as the constant")
    refuses(r"\sin^{-1} x", "ambiguous")
    refuses(r"\sum_{i} x_i", "binder")
    refuses(r"\frac{d}{dx} y", "derivative")
    refuses(r"a \approx b", "not an equality")
    refuses(r"a = b = c", "one relation")
    refuses(r"x'", "prime")
    refuses(r"a = b + o\bigl(x^{2}\bigr)", "asymptotic")
    reads(r"x_{\mathrm{max}} + y_{\text{in}}", "(add x_max y_in)")
    refuses(r"", "empty")

    # sympy round trip on a few expressions, using sympy's own printer for the plain cases
    import sympy as sp
    x, y = sp.symbols("x y", positive=True)
    for e in (x**2 / (1 + y), sp.sqrt(x) * y**sp.Rational(3, 2), sp.sin(x) ** 2 + sp.cos(x) ** 2,
              (x - y) / (x + y), 2 * sp.pi * x, sp.Rational(1, 2) * x**3):
        back = to_sympy(parse(sp.latex(e)), {"x": x, "y": y})
        ok(sp.simplify(back - e) == 0, f"round trip of {sp.latex(e)}")
    return n


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        print(f"texmath: {selftest()} checks ok")
    else:
        for line in sys.argv[1:]:
            try:
                print(show(parse(line)))
            except Unreadable as exc:
                print(f"unreadable: {exc}")
