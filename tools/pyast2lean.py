"""Python source -> LaTeX and Lean 4, by parsing it with ``ast`` (no sympy in the path).

A function written in a small, typed subset of Python becomes
  * a Lean ``def`` (``Id.run do`` with ``let mut``, ``while``, ``if``) that has the same text structure,
  * LaTeX pseudo-code with the same structure, and
  * for a function with one top-level ``while``/``for range`` loop, the Hoare-style verification
    conditions read from its ``assert`` statements:

        assert <pre>              # first statements of the function: precondition
        ...
        assert <inv>              # asserts immediately before the loop: its invariant
        decreases(<expr>)         # optional: a variant (also before the loop)
        while <guard>: ...
        assert <post>             # asserts after the loop: postcondition

    -> theorems  f_init (pre -> inv at entry), f_preserve (inv, guard -> inv after one pass),
       f_exit (inv, not guard -> post), f_variant (0 <= v and v decreases).

The loop body is executed *symbolically* over the ast (parallel substitution, ``if`` becomes an
if-expression), so ``a, b = b, a + b`` and branches are handled exactly.

Types: annotations on arguments are required. ``int`` -> ℤ, ``Nat`` -> ℤ plus the hypothesis 0 ≤ n
(never ℕ: truncated subtraction would silently change meaning), ``float`` -> ℝ (an exact real: rounding is
not modelled), ``bool`` -> Prop. Python's ``//`` and ``%`` are floor operations and print as
``Int.fdiv`` / ``Int.fmod``. Anything else (lists, strings, break, nested loops in the VCs, ...) is
*refused with the line number*, not guessed.

What is claimed: the Lean text is printed from the parse tree. It has not been checked by Lean here.
The VCs are stated about the transition built from the same tree; that the ``do`` block and the
transition agree is by construction, not proved. Theorem bodies are ``sorry`` unless ``--prove``.
"""
from __future__ import annotations
import ast, copy, re, sys, textwrap
from dataclasses import dataclass, field

ATOM, APP, POW, MUL, ADD, CMP, NOT, AND, OR, IFE = 100, 90, 75, 70, 65, 50, 40, 35, 30, 0
Z, R, B = "Z", "R", "B"
LT = {Z: "ℤ", R: "ℝ", B: "Prop"}
KEYWORDS = {"end", "fun", "from", "at", "do", "let", "in", "if", "then", "else", "match", "with", "where",
            "open", "def", "theorem", "by", "have", "show", "return", "for", "while", "mut", "type", "set",
            "instance", "class", "structure", "namespace", "section", "variable", "example", "lemma"}
GREEK = {"alpha", "beta", "gamma", "delta", "epsilon", "theta", "lambda", "mu", "nu", "sigma", "phi", "rho",
         "omega", "tau", "eta", "kappa", "xi", "zeta", "chi", "psi"}
LEAN_GREEK = {"alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ε", "theta": "θ",
              "mu": "μ", "nu": "ν", "sigma": "σ", "phi": "φ", "rho": "ρ", "omega": "ω", "tau": "τ",
              "eta": "η", "kappa": "κ", "xi": "ξ", "zeta": "ζ", "chi": "χ", "psi": "ψ"}


class Refuse(Exception):
    def __init__(self, msg, node=None):
        self.line = getattr(node, "lineno", None)
        super().__init__(f"line {self.line}: {msg}" if self.line else msg)


@dataclass
class X:
    lean: str
    tex: str
    prec: int
    t: str


def L(x: X, p: int) -> str:
    return x.lean if x.prec >= p else f"({x.lean})"


def T(x: X, p: int) -> str:
    return x.tex if x.prec >= p else f"\\left({x.tex}\\right)"


def lid(name: str) -> str:
    if name in KEYWORDS:
        return name + "_"
    if name in LEAN_GREEK:
        return LEAN_GREEK[name]
    return name


def tid(name: str) -> str:
    base, _, sub = name.partition("_")
    b = f"\\{base}" if base in GREEK else (base if len(base) == 1 else f"\\mathrm{{{base}}}")
    return b + (f"_{{{sub.replace('_', ',')}}}" if sub else "")


# ---- expressions -------------------------------------------------------------------------------
def coerce(x: X, t: str) -> X:
    if x.t == t or x.t == B:
        return x
    if t == R and x.t == Z:
        return X(f"(({x.lean}) : ℝ)" if x.prec < ATOM and not re.fullmatch(r"\d+", x.lean)
                 else f"({x.lean} : ℝ)", x.tex, x.prec, R)
    raise Refuse(f"cannot use a {x.t} value as {t}")


def num_join(a: X, b: X):
    t = R if R in (a.t, b.t) else Z
    return coerce(a, t), coerce(b, t), t


CMPOPS = {ast.Lt: ("<", "<"), ast.LtE: ("≤", "\\leq"), ast.Gt: (">", ">"), ast.GtE: ("≥", "\\geq"),
          ast.Eq: ("=", "="), ast.NotEq: ("≠", "\\neq")}


def ex(n, env: dict, notes: list) -> X:
    """Translate an expression. ``env`` maps visible names to 'Z'/'R'/'B'."""
    if isinstance(n, ast.Constant):
        v = n.value
        if isinstance(v, bool):
            return X("True" if v else "False", "\\top" if v else "\\bot", ATOM, B)
        if isinstance(v, int):
            return X(str(v), str(v), ATOM, Z) if v >= 0 else X(f"-{-v}", f"-{-v}", ADD + 1, Z)
        if isinstance(v, float):
            s = repr(v)
            if "e" in s or "inf" in s or "nan" in s:
                raise Refuse(f"float literal {s}", n)
            notes.append("float literals and float arithmetic are modelled as exact reals (no rounding)")
            return X(s if v >= 0 else f"-{s[1:]}", s, ATOM if v >= 0 else ADD + 1, R)
        raise Refuse(f"constant {v!r} (only int, float, bool)", n)
    if isinstance(n, ast.Name):
        if n.id in ("pi",):
            return X("Real.pi", "\\pi", ATOM, R)
        if n.id in ("True", "False"):
            return X(n.id, "\\top" if n.id == "True" else "\\bot", ATOM, B)
        if n.id not in env:
            raise Refuse(f"name '{n.id}' is not defined here (assigned in another block, or not an int/float)", n)
        return X(lid(n.id), tid(n.id), ATOM, env[n.id])
    if isinstance(n, ast.Attribute) and isinstance(n.value, ast.Name) and n.value.id == "math" and n.attr == "pi":
        return X("Real.pi", "\\pi", ATOM, R)
    if isinstance(n, ast.UnaryOp):
        o = ex(n.operand, env, notes)
        if isinstance(n.op, ast.USub):
            if o.t == B:
                raise Refuse("negating a bool", n)
            return X(f"-{L(o, APP)}", f"-{T(o, MUL)}", ADD + 1, o.t)
        if isinstance(n.op, ast.UAdd):
            return o
        if isinstance(n.op, ast.Not):
            return X(f"¬{L(o, ATOM)}", f"\\lnot {T(o, ATOM)}", NOT, B)
        raise Refuse("unary operator", n)
    if isinstance(n, ast.BinOp):
        a, b = ex(n.left, env, notes), ex(n.right, env, notes)
        if B in (a.t, b.t):
            raise Refuse("arithmetic on a bool", n)
        op = type(n.op)
        if op in (ast.Add, ast.Sub):
            a, b, t = num_join(a, b)
            s, ts = ("+", "+") if op is ast.Add else ("-", "-")
            return X(f"{L(a, ADD)} {s} {L(b, ADD + 1)}", f"{T(a, ADD)} {ts} {T(b, ADD + 1)}", ADD, t)
        if op is ast.Mult:
            a, b, t = num_join(a, b)
            return X(f"{L(a, MUL)} * {L(b, MUL + 1)}", f"{T(a, MUL)} \\cdot {T(b, MUL + 1)}", MUL, t)
        if op is ast.Div:
            a, b = coerce(a, R), coerce(b, R)
            notes.append("Python `/` is float division: modelled in ℝ; Lean's x / 0 = 0 where Python raises")
            return X(f"{L(a, MUL)} / {L(b, MUL + 1)}", f"\\frac{{{a.tex}}}{{{b.tex}}}", MUL, R)
        if op in (ast.FloorDiv, ast.Mod):
            if a.t != Z or b.t != Z:
                raise Refuse("// and % are only supported on ints", n)
            fn, tx = ("Int.fdiv", "\\left\\lfloor \\frac{%s}{%s} \\right\\rfloor") if op is ast.FloorDiv \
                else ("Int.fmod", None)
            notes.append("// and % are Python floor division/modulo: Int.fdiv / Int.fmod")
            if tx:
                return X(f"{fn} {L(a, ATOM)} {L(b, ATOM)}", tx % (a.tex, b.tex), APP, Z)
            return X(f"Int.fmod {L(a, ATOM)} {L(b, ATOM)}", f"{T(a, MUL)} \\bmod {T(b, MUL + 1)}", APP, Z)
        if op is ast.Pow:
            e = n.right
            if isinstance(e, ast.Constant) and isinstance(e.value, int) and not isinstance(e.value, bool):
                if e.value >= 0:
                    return X(f"{L(a, ATOM)} ^ {e.value}", f"{T(a, ATOM)}^{{{e.value}}}", POW, a.t)
                a = coerce(a, R)
                return X(f"{L(a, ATOM)} ^ (({e.value} : ℤ))", f"{T(a, ATOM)}^{{{e.value}}}", POW, R)
            if b.t == Z:
                notes.append("x ** n with a variable int n prints as x ^ n.toNat: Python returns a float for n < 0")
                return X(f"{L(a, ATOM)} ^ ({b.lean}).toNat", f"{T(a, ATOM)}^{{{b.tex}}}", POW, a.t)
            a = coerce(a, R)
            notes.append("real exponent is Real.rpow: agrees with Python only for a positive base")
            return X(f"{L(a, ATOM)} ^ ({b.lean} : ℝ)", f"{T(a, ATOM)}^{{{b.tex}}}", POW, R)
        raise Refuse(f"operator {op.__name__}", n)
    if isinstance(n, ast.Compare):
        parts, left = [], ex(n.left, env, notes)
        for o, c in zip(n.ops, n.comparators):
            if type(o) not in CMPOPS:
                raise Refuse(f"comparison {type(o).__name__}", n)
            right = ex(c, env, notes)
            a, b, _ = num_join(left, right)
            ls, ts = CMPOPS[type(o)]
            parts.append((f"{L(a, ADD)} {ls} {L(b, ADD)}", f"{T(a, ADD)} {ts} {T(b, ADD)}"))
            left = right
        if len(parts) == 1:
            return X(parts[0][0], parts[0][1], CMP, B)
        return X(" ∧ ".join(p for p, _ in parts), " \\land ".join(t for _, t in parts), AND, B)
    if isinstance(n, ast.BoolOp):
        xs = [ex(v, env, notes) for v in n.values]
        for x in xs:
            if x.t != B:
                raise Refuse("and/or on a non-bool (Python returns an operand, Lean a Prop)", n)
        and_ = isinstance(n.op, ast.And)
        p = AND if and_ else OR
        return X((" ∧ " if and_ else " ∨ ").join(L(x, p + 1) for x in xs),
                 (" \\land " if and_ else " \\lor ").join(T(x, p + 1) for x in xs), p, B)
    if isinstance(n, ast.IfExp):
        c, a, b = ex(n.test, env, notes), ex(n.body, env, notes), ex(n.orelse, env, notes)
        if c.t != B:
            raise Refuse("if-expression test is not a comparison", n)
        a, b, t = (a, b, B) if a.t == B and b.t == B else num_join(a, b)
        return X(f"if {c.lean} then {a.lean} else {b.lean}",
                 f"\\begin{{cases}} {a.tex} & {c.tex} \\\\ {b.tex} & \\text{{otherwise}} \\end{{cases}}", IFE, t)
    if isinstance(n, ast.Call):
        fname = n.func.attr if isinstance(n.func, ast.Attribute) else getattr(n.func, "id", None)
        args = [ex(a, env, notes) for a in n.args]
        if n.keywords:
            raise Refuse("keyword arguments", n)
        if fname == "abs" and len(args) == 1:
            return X(f"|{args[0].lean}|", f"\\left|{args[0].tex}\\right|", ATOM, args[0].t)
        if fname in ("min", "max") and len(args) == 2:
            a, b, t = num_join(*args)
            return X(f"{fname} {L(a, ATOM)} {L(b, ATOM)}", f"\\{fname}\\left({a.tex}, {b.tex}\\right)", APP, t)
        if fname == "sqrt" and len(args) == 1:
            a = coerce(args[0], R)
            return X(f"Real.sqrt {L(a, ATOM)}", f"\\sqrt{{{a.tex}}}", APP, R)
        if fname in ("exp", "log", "sin", "cos", "tan") and len(args) == 1:
            a = coerce(args[0], R)
            return X(f"Real.{fname} {L(a, ATOM)}", f"\\{fname}\\left({a.tex}\\right)", APP, R)
        if fname == "factorial" and len(args) == 1 and args[0].t == Z:
            notes.append("factorial(n) prints as ((n.toNat)! : ℤ): n < 0 is outside Python's domain")
            return X(f"((Int.toNat {L(args[0], ATOM)}).factorial : ℤ)", f"{T(args[0], ATOM)}!", ATOM, Z)
        if fname == "comb" and len(args) == 2 and all(a.t == Z for a in args):
            return X(f"((Nat.choose (Int.toNat {L(args[0], ATOM)}) (Int.toNat {L(args[1], ATOM)})) : ℤ)",
                     f"\\binom{{{args[0].tex}}}{{{args[1].tex}}}", ATOM, Z)
        raise Refuse(f"call to '{fname}' (supported: abs min max sqrt exp log sin cos tan factorial comb)", n)
    raise Refuse(f"{type(n).__name__} expression", n)


# ---- ast helpers -------------------------------------------------------------------------------
class _Subst(ast.NodeTransformer):
    def __init__(self, env):
        self.env = env

    def visit_Name(self, n):
        return copy.deepcopy(self.env[n.id]) if n.id in self.env else n


def subst(node, env):
    return _Subst(env).visit(copy.deepcopy(node))


def names(node) -> set:
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name)}


def desugar(stmts: list) -> list:
    """``for i in range(a, b, c)`` -> ``i = a; while i < b: body; i += c`` (recursively)."""
    out = []
    for s in stmts:
        for field_ in ("body", "orelse"):
            if hasattr(s, field_) and isinstance(getattr(s, field_), list):
                setattr(s, field_, desugar(getattr(s, field_)))
        if isinstance(s, ast.For):
            it = s.iter
            if not (isinstance(s.target, ast.Name) and isinstance(it, ast.Call) and getattr(it.func, "id", "") == "range"
                    and 1 <= len(it.args) <= 3 and not s.orelse):
                raise Refuse("only `for i in range(...)` loops are supported", s)
            a = it.args
            lo, hi = (ast.Constant(0), a[0]) if len(a) == 1 else (a[0], a[1])
            step = a[2] if len(a) == 3 else ast.Constant(1)
            if not (isinstance(step, ast.Constant) and isinstance(step.value, int) and step.value != 0):
                raise Refuse("range step must be a non-zero int literal", s)
            i = s.target.id
            cmp_ = ast.Lt() if step.value > 0 else ast.Gt()
            out.append(ast.copy_location(ast.Assign([ast.Name(i, ast.Store())], lo), s))
            body = s.body + [ast.copy_location(ast.AugAssign(ast.Name(i, ast.Store()), ast.Add(), step), s)]
            test = ast.Compare(ast.Name(i, ast.Load()), [cmp_], [hi])
            w = ast.copy_location(ast.While(test, body, []), s)
            ast.fix_missing_locations(w)
            w._from_for = i
            out.append(w)
        else:
            out.append(s)
    return out


def is_decreases(s) -> bool:
    return isinstance(s, ast.Expr) and isinstance(s.value, ast.Call) and getattr(s.value.func, "id", "") == "decreases"


# ---- statement translation (Lean do-block + LaTeX rows) ----------------------------------------
class Gen:
    def __init__(self, fn: "FnInfo"):
        self.fn, self.notes = fn, fn.notes
        self.lean: list[str] = []
        self.rows: list[tuple[int, str]] = []
        self.snap: dict[int, dict] = {}       # id(While) -> type env at that point
        self.returns: list[X] = []

    def emit(self, ind, lean, tex):
        if lean is not None:
            self.lean.append("  " * ind + lean)
        if tex is not None:
            self.rows.append((ind, tex))

    def block(self, stmts, env, ind, top=False):
        env = dict(env)
        pending: list[ast.stmt] = []          # asserts/decreases before a loop
        for s in stmts:
            if isinstance(s, ast.Assert) or is_decreases(s):
                pending.append(s)
                continue
            self.stmt(s, env, ind, pending)
            pending = [] if not isinstance(s, ast.While) else []
        for s in pending:                      # trailing asserts: postconditions / documentation
            self.stmt(s, env, ind, [])
        return env

    def stmt(self, s, env, ind, pending):
        if isinstance(s, ast.Assert):
            x = ex(s.test, env, self.notes)
            self.emit(ind, f"-- assert: {x.lean}", f"\\{{\\,{x.tex}\\,\\}}")
        elif is_decreases(s):
            x = ex(s.value.args[0], env, self.notes)
            self.emit(ind, f"-- decreases: {x.lean}", f"\\mathrm{{decreases}}\\ {x.tex}")
        elif isinstance(s, ast.Pass):
            pass
        elif isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant):
            pass                                # docstring
        elif isinstance(s, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
            self.assign(s, env, ind)
        elif isinstance(s, ast.If):
            c = ex(s.test, env, self.notes)
            if c.t != B:
                raise Refuse("`if` condition is not a comparison", s)
            self.emit(ind, f"if {c.lean} then", f"\\mathbf{{if}}\\ {c.tex}\\ \\mathbf{{then}}")
            self.block(s.body, env, ind + 1) if s.body else None
            if s.orelse:
                self.emit(ind, "else", "\\mathbf{else}")
                self.block(s.orelse, env, ind + 1)
        elif isinstance(s, ast.While):
            if s.orelse:
                raise Refuse("while-else", s)
            self.snap[id(s)] = dict(env)
            for p in pending:
                self.stmt(p, env, ind, [])
            c = ex(s.test, env, self.notes)
            if c.t != B:
                raise Refuse("`while` condition is not a comparison", s)
            self.emit(ind, f"while {c.lean} do", f"\\mathbf{{while}}\\ {c.tex}\\ \\mathbf{{do}}")
            self.block(s.body, env, ind + 1)
        elif isinstance(s, ast.Return):
            if s.value is None:
                raise Refuse("bare return", s)
            x = ex(s.value, env, self.notes)
            self.returns.append(x)
            self.emit(ind, f"return {x.lean}" if x.t != B else None, None)
            self.emit(ind, None, f"\\mathbf{{return}}\\ {x.tex}")
        elif isinstance(s, (ast.Break, ast.Continue)):
            raise Refuse("break/continue: write the loop guard instead", s)
        else:
            raise Refuse(f"{type(s).__name__} statement", s)

    def assign(self, s, env, ind):
        if isinstance(s, ast.AugAssign):
            if not isinstance(s.target, ast.Name):
                raise Refuse("assignment target must be a name", s)
            tgts, val = [s.target], ast.BinOp(ast.Name(s.target.id, ast.Load()), s.op, s.value)
            ast.copy_location(val, s)
        elif isinstance(s, ast.AnnAssign):
            if s.value is None or not isinstance(s.target, ast.Name):
                raise Refuse("annotated assignment needs a name and a value", s)
            tgts, val = [s.target], s.value
        else:
            if len(s.targets) != 1:
                raise Refuse("chained assignment", s)
            tgts, val = s.targets, s.value
            if isinstance(tgts[0], ast.Tuple):
                tgts = tgts[0].elts
                if not (isinstance(val, ast.Tuple) and len(val.elts) == len(tgts)):
                    raise Refuse("tuple assignment needs a tuple of the same length on the right", s)
                vals = val.elts
            else:
                vals = [val]
        if not isinstance(s, ast.Assign) or not isinstance(s.targets[0], ast.Tuple):
            vals = [val]
        if not all(isinstance(t, ast.Name) for t in tgts):
            raise Refuse("assignment target must be a name", s)
        xs = [ex(v, env, self.notes) for v in vals]
        if any(x.t == B for x in xs):
            raise Refuse("bool-valued variables are not supported", s)
        new = [t.id not in env for t in tgts]
        if any(new) and not all(new):
            raise Refuse("tuple assignment mixing new and existing variables", s)
        final = []
        for t, x in zip(tgts, xs):
            if t.id in env:
                if env[t.id] == Z and x.t == R:
                    raise Refuse(f"'{t.id}' is an int and is assigned a float (Lean variables have one type)", s)
                x = coerce(x, env[t.id])
            final.append(x)
        if len(tgts) == 1:
            t, x = tgts[0], final[0]
            if new[0]:
                env[t.id] = x.t
                self.emit(ind, f"let mut {lid(t.id)} : {LT[x.t]} := {x.lean}", f"{tid(t.id)} \\gets {x.tex}")
            else:
                self.emit(ind, f"{lid(t.id)} := {x.lean}", f"{tid(t.id)} \\gets {x.tex}")
        else:
            lhs = ", ".join(lid(t.id) for t in tgts)
            rhs = ", ".join(f"({x.lean} : {LT[x.t]})" for x in final)
            tl = ", ".join(tid(t.id) for t in tgts)
            tr = ", ".join(x.tex for x in final)
            if all(new):
                for t, x in zip(tgts, final):
                    env[t.id] = x.t
            self.emit(ind, f"{'let mut ' if all(new) else ''}({lhs}) := ({rhs})", f"({tl}) \\gets ({tr})")


# ---- function-level ----------------------------------------------------------------------------
@dataclass
class FnInfo:
    name: str
    node: ast.FunctionDef
    notes: list = field(default_factory=list)
    args: list = field(default_factory=list)        # (name, type, nat)
    pre: list = field(default_factory=list)
    ret: str = Z
    lean_def: str = ""
    rows: list = field(default_factory=list)
    vcs: list = field(default_factory=list)         # (name, lean text, tex text)
    vc_refused: str = ""
    refused: str = ""
    inv_tex: str = ""


def annotation_type(a, node):
    nm = a.id if isinstance(a, ast.Name) else (a.attr if isinstance(a, ast.Attribute) else None)
    if nm == "int":
        return Z, False
    if nm == "Nat":
        return Z, True
    if nm == "float":
        return R, False
    raise Refuse("every argument needs an annotation: int, Nat or float", node)


def group_binders(vars_: list) -> str:
    groups: dict[str, list[str]] = {}
    for n, t in vars_:
        groups.setdefault(LT[t], []).append(lid(n))
    return " ".join(f"({' '.join(v)} : {t})" for t, v in groups.items())


def translate_function(fd: ast.FunctionDef, prove: bool = False) -> FnInfo:
    info = FnInfo(fd.name, fd)
    try:
        _translate(info, prove)
    except Refuse as e:
        info.refused = str(e)
    return info


def _translate(info: FnInfo, prove: bool):
    fd = info.node
    env: dict = {}
    for a in fd.args.args:
        if a.annotation is None:
            raise Refuse(f"argument '{a.arg}' has no annotation", a)
        t, nat = annotation_type(a.annotation, a)
        env[a.arg] = t
        info.args.append((a.arg, t, nat))
    if fd.args.vararg or fd.args.kwarg or fd.args.defaults or fd.args.kwonlyargs:
        raise Refuse("defaults, *args and keyword-only arguments", fd)
    body = desugar(list(fd.body))
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
        body = body[1:]
    pre = []
    while body and isinstance(body[0], ast.Assert):
        pre.append(body.pop(0))
    info.pre = [p.test for p in pre]
    g = Gen(info)
    for p in info.pre:
        g_x = ex(p, env, info.notes)
        g.emit(1, f"-- requires: {g_x.lean}", f"\\{{\\,{g_x.tex}\\,\\}}")
    g.block(body, env, 1, top=True)
    if not g.returns:
        raise Refuse("the function never returns a value", fd)
    ret = fd.returns
    info.ret = annotation_type(ret, fd)[0] if ret is not None else g.returns[0].t
    arg_t = [(n, t) for n, t, _ in info.args]
    noncomp = "noncomputable " if R in [t for _, t in arg_t] + [info.ret] or any("ℝ" in l for l in g.lean) else ""
    info.lean_def = (f"{noncomp}def {lid(info.name)} {group_binders(arg_t)} : {LT[info.ret]} := Id.run do\n"
                     + "\n".join(g.lean) + "\n")
    info.rows = g.rows
    try:
        _vcs(info, g, body, env, prove)
    except Refuse as e:
        info.vc_refused = str(e)


def _vcs(info: FnInfo, g: Gen, body: list, env0: dict, prove: bool):
    loops = [i for i, s in enumerate(body) if isinstance(s, ast.While)]
    if not loops:
        raise Refuse("no top-level loop: no verification conditions")
    if len(loops) > 1 or any(isinstance(n, ast.While) for s in body for n in ast.walk(s) if n is not s):
        raise Refuse("verification conditions are generated for exactly one loop, not nested or several loops")
    k = loops[0]
    w = body[k]
    j = k
    block_ = []
    while j > 0 and (isinstance(body[j - 1], ast.Assert) or is_decreases(body[j - 1])):
        j -= 1
        block_.insert(0, body[j])
    prefix, suffix = body[:j], body[k + 1:]
    invs = [b.test for b in block_ if isinstance(b, ast.Assert)]
    variants = [b.value.args[0] for b in block_ if is_decreases(b)]
    tenv = g.snap[id(w)]
    S = list(tenv.items())                                  # state variables at loop entry
    state = set(tenv)
    # init environment: execute the prefix (straight line only)
    class Ctx:
        asserts: list
    def sx(stmts, e, collect=None):
        e = dict(e)
        for st in stmts:
            if isinstance(st, ast.Assign) or isinstance(st, ast.AnnAssign) or isinstance(st, ast.AugAssign):
                if isinstance(st, ast.AugAssign):
                    tg, vs = [st.target], [ast.BinOp(ast.Name(st.target.id, ast.Load()), st.op, st.value)]
                elif isinstance(st, ast.AnnAssign):
                    tg, vs = [st.target], [st.value]
                elif isinstance(st.targets[0], ast.Tuple):
                    tg, vs = st.targets[0].elts, st.value.elts
                else:
                    tg, vs = st.targets, [st.value]
                new = {t.id: subst(v, e) for t, v in zip(tg, vs)}      # parallel
                e.update(new)
            elif isinstance(st, ast.If):
                c = subst(st.test, e)
                e1, e2 = sx(st.body, e), sx(st.orelse, e)
                for key in set(e1) | set(e2):
                    a, b = e1.get(key, ast.Name(key, ast.Load())), e2.get(key, ast.Name(key, ast.Load()))
                    e[key] = a if ast.dump(a) == ast.dump(b) else ast.IfExp(copy.deepcopy(c), a, b)
            elif isinstance(st, ast.Assert):
                if collect is None:
                    pass
                else:
                    collect.append(subst(st.test, e))
            elif is_decreases(st) or isinstance(st, ast.Pass) or isinstance(st, ast.Expr):
                pass
            elif isinstance(st, ast.Return):
                if collect is not None:
                    break
                raise Refuse("return before the loop")
            else:
                raise Refuse(f"{type(st).__name__} inside the loop/prefix is not handled by the VC generator", st)
        return e

    init = sx(prefix, {})
    step = sx(w.body, {})
    post_list: list = []
    sx(suffix, {}, collect=post_list)
    for x_ in invs + [w.test] + post_list + variants:
        free = names(x_) - state - {"abs", "min", "max", "sqrt", "exp", "log", "sin", "cos", "tan", "factorial", "comb", "pi", "math", "True", "False"}
        if free:
            raise Refuse(f"'{sorted(free)[0]}' is used in the invariant/guard/postcondition but is not defined before the loop")
    notes = info.notes
    tenv_all = dict(tenv)
    pre_x = [ex(p, env0, notes) for p in info.pre]
    inv_x = [ex(i, tenv_all, notes) for i in invs]
    if not inv_x:
        raise Refuse("no invariant: put `assert ...` immediately before the loop")
    name = lid(info.name)
    svars = [(n, t) for n, t in S]
    tac = ("by\n  intros\n  try simp only [%s_inv, %s_post] at *\n  first | omega | nlinarith | sorry" % (name, name)
           if prove else "by\n  sorry")
    conj = lambda xs: " ∧ ".join(L(x, AND + 1) for x in xs)
    tconj = lambda xs: " \\land ".join(T(x, AND + 1) for x in xs)
    inv_def = f"def {name}_inv {group_binders(svars)} : Prop :=\n  {conj(inv_x)}\n"
    pre_names = [n for n, _, nat in info.args if nat]
    arg_vars = [(n, t) for n, t, _ in info.args]

    def hyps(with_pre=True, extra=()):
        hs = [f"(hnat_{lid(n)} : 0 ≤ {lid(n)})" for n in pre_names]
        hs += [f"(hpre{i} : {x.lean})" for i, x in enumerate(pre_x)] if with_pre else []
        return " ".join(hs + list(extra))

    def call_inv(envmap):
        args = []
        for n, t in svars:
            v = envmap.get(n)
            x = ex(v, tenv_all, notes) if v is not None else ex(ast.Name(n, ast.Load()), tenv_all, notes)
            args.append(L(coerce(x, t) if x.t != B else x, ATOM))
        return f"{name}_inv " + " ".join(args)

    tex_inv_sub = lambda e: tconj([ex(subst(i, e), tenv_all, notes) for i in invs])
    vcs = [(f"{name}_inv", inv_def, None)]
    arg_b = group_binders(arg_vars)
    vcs.append((f"{name}_init",
                f"theorem {name}_init {arg_b} {hyps()} :\n    {call_inv(init)} := {tac}\n",
                f"\\text{{pre}} \\Rightarrow {tex_inv_sub(init)}"))
    guard = ex(w.test, tenv_all, notes)
    st_b = group_binders(svars)
    inv_now = f"{name}_inv " + " ".join(lid(n) for n, _ in svars)
    vcs.append((f"{name}_preserve",
                f"theorem {name}_preserve {st_b} {hyps(extra=[f'(hinv : {inv_now})', f'(hguard : {guard.lean})'])} :\n"
                f"    {call_inv(step)} := {tac}\n",
                f"{tconj(inv_x)} \\land {guard.tex} \\Rightarrow {tex_inv_sub(step)}"))
    if post_list:
        post_x = [ex(p, tenv_all, notes) for p in post_list]
        vcs.append((f"{name}_exit",
                    f"theorem {name}_exit {st_b} {hyps(extra=[f'(hinv : {inv_now})', f'(hexit : ¬({guard.lean}))'])} :\n"
                    f"    {conj(post_x)} := {tac}\n",
                    f"{tconj(inv_x)} \\land \\lnot\\left({guard.tex}\\right) \\Rightarrow {tconj(post_x)}"))
    for vi, v in enumerate(variants):
        vx = ex(v, tenv_all, notes)
        vn = ex(subst(v, step), tenv_all, notes)
        if vx.t != Z:
            raise Refuse("decreases(...) must be an int expression")
        vcs.append((f"{name}_variant{vi or ''}",
                    f"theorem {name}_variant{vi or ''} {st_b} {hyps(extra=[f'(hinv : {inv_now})', f'(hguard : {guard.lean})'])} :\n"
                    f"    0 ≤ {vx.lean} ∧ {vn.lean} < {vx.lean} := {tac}\n",
                    f"{tconj(inv_x)} \\land {guard.tex} \\Rightarrow 0 \\leq {vx.tex} \\land {vn.tex} < {vx.tex}"))
    info.vcs = vcs
    info.inv_tex = tconj(inv_x)
    if not post_list:
        info.notes.append("no assert after the loop: no exit condition was generated")


# ---- rendering ---------------------------------------------------------------------------------
def san(s: str) -> str:
    return s.replace("/-", "/ -").replace("-/", "- /")


def lean_text(info: FnInfo, src_name: str = "") -> str:
    if info.refused:
        return f"-- {info.name}: refused ({info.refused})\n"
    doc = ["/--", f"`{info.name}` translated from Python (`{src_name}`) by tools/pyast2lean.py."]
    for n in dict.fromkeys(info.notes):
        doc.append(f"Note: {san(n)}.")
    doc.append("-/")
    out = "\n".join(doc) + "\n" + info.lean_def
    if info.vc_refused:
        out += f"\n-- verification conditions not generated: {info.vc_refused}\n"
    for _, text, _ in info.vcs:
        out += "\n" + text
    return out


def tex_text(info: FnInfo) -> str:
    esc = info.name.replace("_", "\\_")
    if info.refused:
        return f"\\subsection*{{\\texttt{{{esc}}}}}\nRefused: {info.refused.replace('_', chr(92) + '_')}.\n"
    sig = ", ".join(tid(n) for n, _, _ in info.args)
    lines = [f"\\subsection*{{\\texttt{{{esc}}}}}", f"\\[ \\mathrm{{{info.name.replace('_', chr(92) + '_')}}}({sig}) \\]", "\\[\\begin{array}{l}"]
    for ind, t in info.rows:
        lines.append(("\\quad " * max(ind - 1, 0)) + t + " \\\\")
    lines.append("\\end{array}\\]")
    if info.vcs:
        lines.append("Verification conditions (from the \\texttt{assert} statements):")
        for nm, _, tx in info.vcs:
            if tx:
                lines.append(f"\\[ \\text{{{nm[len(lid(info.name)) + 1:]}:}}\\quad {tx} \\]")
    elif info.vc_refused:
        lines.append(f"No verification conditions: {info.vc_refused.replace('_', chr(92) + '_')}.")
    return "\n".join(lines) + "\n"


def translate_source(src: str, prove: bool = False) -> list[FnInfo]:
    tree = ast.parse(textwrap.dedent(src))
    return [translate_function(n, prove) for n in tree.body if isinstance(n, ast.FunctionDef)]


# ---- selftest ----------------------------------------------------------------------------------
SUM_SRC = '''
def sum_to(n: Nat) -> int:
    assert n >= 0
    s = 0
    i = 0
    assert 2 * s == i * (i - 1) and 0 <= i and i <= n + 1 - 1
    decreases(n - i)
    while i < n:
        s = s + i
        i = i + 1
    assert 2 * s == n * (n - 1)
    return s
'''

FIB_SRC = '''
def fib(n: Nat) -> int:
    assert n >= 0
    a, b = 0, 1
    k = 0
    assert a >= 0 and b >= a and k <= n
    while k < n:
        a, b = b, a + b
        k += 1
    assert a >= 0
    return a
'''


def selftest():
    n = [0]

    def ok(c, m=""):
        n[0] += 1
        assert c, m

    s = translate_source(SUM_SRC)[0]
    ok(not s.refused, s.refused)
    t = lean_text(s)
    ok("def sum_to (n : ℤ) : ℤ := Id.run do" in t, t)
    ok("let mut s : ℤ := 0" in t and "while i < n do" in t and "s := s + i" in t, t)
    ok("theorem sum_to_init" in t and "theorem sum_to_preserve" in t and "theorem sum_to_exit" in t
       and "theorem sum_to_variant" in t, t)
    ok("(hnat_n : 0 ≤ n)" in t, t)
    ok("sum_to_inv n (s + i) (i + 1)" in t, t)                   # parallel substitution of the body
    ok("0 ≤ n - i ∧ n - (i + 1) < n - i" in t, t)
    f = translate_source(FIB_SRC)[0]
    ft = lean_text(f)
    ok("(a, b) := ((b : ℤ), (a + b : ℤ))" in ft, ft)
    ok("fib_inv n b (a + b) (k + 1)" in ft or "fib_inv n b (a + b) (k + 1)".replace("n b", "n b") in ft, ft)
    # for-loop desugaring
    fr = translate_source("def tri(n: Nat) -> int:\n  assert n >= 0\n  t = 0\n  assert 2*t == i*(i-1)\n"
                          "  for i in range(0, n):\n    t += i\n  return t\n")[0]
    ok("while i < n do" in lean_text(fr) or fr.refused, lean_text(fr))
    # refusals carry a line number and a reason
    for src, why in (("def f(x: int) -> int:\n    while x > 0:\n        break\n    return x\n", "break/continue"),
                     ("def f(x) -> int:\n    return x\n", "no annotation"),
                     ("def f(x: int) -> int:\n    return [x]\n", "List expression"),
                     ("def f(x: int) -> int:\n    y = 1\n    y = 1.5\n    return x\n", "assigned a float")):
        r = translate_source(src)[0]
        ok(why in r.refused, (src, r.refused))
    ok("line 3" in translate_source("def f(x: int) -> int:\n    while x > 0:\n        break\n    return x\n")[0].refused)
    # python semantics: floor division, no silent reals
    fd = translate_source("def h(a: int, b: int) -> int:\n    return a // b + a % b\n")[0]
    ok("Int.fdiv a b + Int.fmod a b" in lean_text(fd), lean_text(fd))
    fm = translate_source("def m(a: int, b: float) -> float:\n    return a / 2 + b * a\n")[0]
    ok("noncomputable def" in lean_text(fm) and "(a : ℝ)" in lean_text(fm), lean_text(fm))
    # a variable used outside the block that defined it is refused
    r = translate_source("def f(x: int) -> int:\n    if x > 0:\n        y = 1\n    return y\n")[0]
    ok(r.refused and "y" in r.refused, r.refused)
    # tex is produced and contains the invariant and the guard
    tx = tex_text(s)
    ok("\\mathbf{while}" in tx and "\\Rightarrow" in tx, tx)
    print(f"pyast2lean: {n[0]} checks ok")


def main(argv=None):
    import argparse, pathlib
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("files", nargs="*")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--prove", action="store_true", help="try omega/nlinarith before sorry in the VC theorems")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    here = pathlib.Path(__file__).resolve().parent
    files = a.files or sorted(str(p) for p in (here / "py_models").glob("*.py"))
    lean, tex = [], []
    for f in files:
        for info in translate_source(open(f).read(), a.prove):
            lean.append(lean_text(info, pathlib.Path(f).name))
            tex.append(tex_text(info))
            print(f"{info.name}: " + ("refused: " + info.refused if info.refused
                                      else f"{len(info.vcs)} vc part(s)" + (f" ({info.vc_refused})" if info.vc_refused else "")))
    root = here.parent / "build"
    (root / "lean").mkdir(parents=True, exist_ok=True)
    (root / "tex").mkdir(parents=True, exist_ok=True)
    (root / "lean" / "PyAst_mathlib.lean").write_text(
        "import Mathlib\n\nopen scoped Classical\n\nnamespace OAI.PyAst\n\n" + "\n".join(lean) + "\nend OAI.PyAst\n")
    (root / "tex" / "pyast.tex").write_text("\n".join(tex))
    print(f"wrote {root / 'lean' / 'PyAst_mathlib.lean'}")


if __name__ == "__main__":
    main()
