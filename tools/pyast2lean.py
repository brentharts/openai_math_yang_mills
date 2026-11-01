"""Python source -> LaTeX and Lean 4, by parsing it with ``ast`` (no sympy in the path).

A function written in a small, typed subset of Python becomes
  * a Lean ``def`` (``Id.run do`` with ``let mut``, ``while``, ``if``, ``break``, ``continue``) that has the same
    text structure,
  * LaTeX pseudo-code with the same structure, and
  * Hoare-style verification conditions read from its ``assert`` statements:

        assert <pre>              # first statements of the function: precondition
        ...
        assert <inv>              # asserts immediately before a loop: that loop's invariant
        decreases(<expr>)         # optional: a variant (also before the loop)
        while <guard>: ...
        assert <post>             # asserts right after the loop: its exit condition

    -> theorems  f_init (what holds before the loop gives inv), f_preserve (inv, guard -> inv after one pass),
       f_exit (inv, not guard -> post), f_variant (0 <= v and v decreases). Any other assert is a theorem
       f_assertN. A list index a[i] gives f_safeN (0 <= i < len(a)), a call to a function with a precondition
       gives f_safeN (that precondition at the call).

Each loop has its own invariant, so loops may be nested or follow each other (f_L1_init, f_L2_preserve, ...); an
inner loop's conditions assume the outer invariant and guard, except the parts that mention variables the inner
loop changes (state what is needed in the inner invariant). ``break`` makes the exit condition a theorem on
every way out (f_exit_break1, ...); ``continue`` and ``if`` with a jump split the body into paths
(f_preserve_path1, ...). The loop body is executed *symbolically* over the ast (parallel substitution,
``if`` becomes an if-expression), so ``a, b = b, a + b`` and branches are handled exactly.

Types: annotations on arguments are required. ``int`` -> ℤ, ``Nat`` -> ℤ plus the hypothesis 0 ≤ n
(never ℕ: truncated subtraction would silently change meaning), ``float`` -> ℝ (an exact real: rounding is
not modelled), ``bool`` -> Bool (a comparison is ``decide``d when stored; in an assert it is a Prop),
``list[int|float|bool]`` -> Array (a value: a list argument is never changed in place, ``b = list(a)`` copies).
Python's ``//`` and ``%`` are floor operations and print as ``Int.fdiv`` / ``Int.fmod``. Functions of one source
may call each other (callees are printed first; recursion is refused). ``for i in range(...)`` and ``for x in a``
become ``while`` loops, and are refused where that would differ from Python (range bound changed in the body,
loop variable used after the loop). Anything else (strings, dicts, comprehensions, ...) is *refused with the
line number*, not guessed.

What is claimed: the Lean text is printed from the parse tree. It has not been checked by Lean here.
The VCs are stated about the transition built from the same tree; that the ``do`` block and the
transition agree is by construction, not proved. Theorem bodies are ``sorry`` unless ``--prove``.
``--smoke`` evaluates every theorem on a grid of small values: it finds a wrong invariant or a wrong
symbolic execution, it proves nothing.
"""
from __future__ import annotations
import ast, copy, itertools, re, sys, textwrap
from dataclasses import dataclass, field

ATOM, APP, POW, MUL, ADD, CMP, NOT, AND, OR, IFE = 100, 90, 75, 70, 65, 50, 40, 35, 30, 0
Z, R, B = "Z", "R", "B"
AZ, AR, AB = "AZ", "AR", "AB"                       # list[int], list[float], list[bool] -> Array
ARR = {AZ: Z, AR: R, AB: B}                         # array type -> element type
ARR_OF = {Z: AZ, R: AR, B: AB}
LT = {Z: "ℤ", R: "ℝ", B: "Bool", AZ: "Array ℤ", AR: "Array ℝ", AB: "Array Bool"}
KEYWORDS = {"end", "fun", "from", "at", "do", "let", "in", "if", "then", "else", "match", "with", "where",
            "open", "def", "theorem", "by", "have", "show", "return", "for", "while", "mut", "type", "set",
            "instance", "class", "structure", "namespace", "section", "variable", "example", "lemma"}
GREEK = {"alpha", "beta", "gamma", "delta", "epsilon", "theta", "lambda", "mu", "nu", "sigma", "phi", "rho",
         "omega", "tau", "eta", "kappa", "xi", "zeta", "chi", "psi"}
LEAN_GREEK = {"alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ε", "theta": "θ",
              "mu": "μ", "nu": "ν", "sigma": "σ", "phi": "φ", "rho": "ρ", "omega": "ω", "tau": "τ",
              "eta": "η", "kappa": "κ", "xi": "ξ", "zeta": "ζ", "chi": "χ", "psi": "ψ"}


@dataclass
class FuncSig:
    """What a caller needs to know about another translated function."""
    name: str
    params: list                    # (name, type, is_nat)
    ret: str
    pre: list = field(default_factory=list)      # precondition asts over the parameter names


FUNCS: dict[str, FuncSig] = {}       # translated functions of the source being processed
FUNCS_REFUSED: dict[str, str] = {}   # translated functions that were refused, with the reason
INVS: dict[str, tuple] = {}          # "__invK" -> (lean predicate name, [(var, type)], [invariant asts])
BUILTIN_NAMES = {"abs", "min", "max", "sqrt", "exp", "log", "sin", "cos", "tan", "factorial", "comb", "len", "list",
                 "pi", "math", "True", "False", "__set", "__push"}


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
    bl: str = ""        # for t == B: the Lean ``Bool`` term when the value is natively a Bool (variable, call)


def as_bool(x: X) -> str:
    """A Lean ``Bool`` term for a bool-valued X (a Prop is wrapped in ``decide``; Classical is open)."""
    if x.bl:
        return x.bl
    return "true" if x.lean == "True" else "false" if x.lean == "False" else f"(decide ({x.lean}))"


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
            return X("True" if v else "False", "\\top" if v else "\\bot", ATOM, B, "true" if v else "false")
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
            return X(n.id, "\\top" if n.id == "True" else "\\bot", ATOM, B, n.id.lower())
        if n.id not in env:
            raise Refuse(f"name '{n.id}' is not defined here (assigned in another block, or not a supported value)", n)
        if env[n.id] == B:
            return X(f"({lid(n.id)} = true)", tid(n.id), ATOM, B, lid(n.id))
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
        if isinstance(n.op, ast.Mult) and isinstance(n.left, ast.List) and len(n.left.elts) == 1:
            c, cnt = ex(n.left.elts[0], env, notes), ex(n.right, env, notes)      # [c] * n
            if cnt.t != Z or c.t in ARR:
                raise Refuse("[c] * n needs an int count and a scalar element", n)
            cv = as_bool(c) if c.t == B else L(c, ATOM)
            return X(f"(Array.replicate ({cnt.lean}).toNat {cv})", f"\\left[{c.tex}\\right]^{{{cnt.tex}}}",
                     ATOM, ARR_OF[c.t])
        a, b = ex(n.left, env, notes), ex(n.right, env, notes)
        if a.t in ARR or b.t in ARR:
            raise Refuse("arithmetic on a list (only [c] * n, indexing, len, append, list(a) are supported)", n)
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
            if left.t in ARR or right.t in ARR:
                raise Refuse("comparing lists is not supported", n)
            if left.t == B or right.t == B:
                if not (left.t == B and right.t == B and isinstance(o, (ast.Eq, ast.NotEq))):
                    raise Refuse("a bool can only be compared with == / != against a bool", n)
                neg = isinstance(o, ast.NotEq)
                ln, tn = ("¬" if neg else ""), ("\\lnot " if neg else "")
                parts.append((f"{ln}({left.lean} ↔ {right.lean})",
                              f"{tn}\\left({left.tex} \\leftrightarrow {right.tex}\\right)"))
                left = right
                continue
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
        if a.t == b.t and (a.t == B or a.t in ARR):
            t = a.t
        else:
            a, b, t = num_join(a, b)
        return X(f"if {c.lean} then {a.lean} else {b.lean}",
                 f"\\begin{{cases}} {a.tex} & {c.tex} \\\\ {b.tex} & \\text{{otherwise}} \\end{{cases}}", IFE, t)
    if isinstance(n, ast.List):
        if not n.elts:
            raise Refuse("an empty list needs a declared type: write `a: list[int] = []`", n)
        xs = [ex(e, env, notes) for e in n.elts]
        if any(x.t in ARR for x in xs):
            raise Refuse("nested lists are not supported", n)
        if all(x.t == B for x in xs):
            t, items = AB, [as_bool(x) for x in xs]
        elif B in [x.t for x in xs]:
            raise Refuse("a list mixing bools and numbers", n)
        else:
            et = R if R in [x.t for x in xs] else Z
            xs = [coerce(x, et) for x in xs]
            t, items = ARR_OF[et], [x.lean for x in xs]
        return X("#[" + ", ".join(items) + "]", "\\left[" + ", ".join(x.tex for x in xs) + "\\right]", ATOM, t)
    if isinstance(n, ast.Subscript):
        a, i = ex(n.value, env, notes), ex(n.slice, env, notes)
        if a.t not in ARR:
            raise Refuse("indexing something that is not a list", n)
        if i.t != Z:
            raise Refuse("a list index must be an int", n)
        notes.append("a[i] prints as a[i.toNat]!: negative indices (Python counts from the end) and out-of-range "
                     "indices are outside the model; the _safe theorems state 0 ≤ i < len(a)")
        et, tx = ARR[a.t], f"{T(a, ATOM)}\\left[{i.tex}\\right]"
        if et == B:
            return X(f"({L(a, ATOM)}[({i.lean}).toNat]! = true)", tx, ATOM, B, f"{L(a, ATOM)}[({i.lean}).toNat]!")
        return X(f"{L(a, ATOM)}[({i.lean}).toNat]!", tx, ATOM, et)
    if isinstance(n, ast.Call):
        fname = n.func.attr if isinstance(n.func, ast.Attribute) else getattr(n.func, "id", None)
        if fname == "__empty" and isinstance(n.func, ast.Name):
            return X("#[]", "\\left[\\right]", ATOM, n.args[0].value)
        args = [ex(a, env, notes) for a in n.args]
        if n.keywords:
            raise Refuse("keyword arguments", n)
        if isinstance(n.func, ast.Name) and n.func.id in INVS:                    # loop invariant instance
            lean_name, svars, inv_asts = INVS[n.func.id]
            if len(args) != len(svars):
                raise Refuse("internal: invariant arity", n)
            parts, mapping = [], {}
            for (v, vt), x, node in zip(svars, args, n.args):
                x = x if vt in ARR else coerce(x, vt)
                parts.append(as_bool(x) if vt == B else L(x, ATOM))
                mapping[v] = node
            tex_parts = [ex(subst(i, mapping), env, notes) for i in inv_asts]
            return X(f"{lean_name} " + " ".join(parts), " \\land ".join(T(x, AND + 1) for x in tex_parts), APP, B)
        if fname in ("__set", "__push") and isinstance(n.func, ast.Name):
            a = args[0]
            if a.t not in ARR:
                raise Refuse("internal: list update on a non-list", n)
            et = ARR[a.t]
            if fname == "__set":
                i, v = args[1], coerce(args[2], et) if args[2].t != B else args[2]
                vt = as_bool(v) if et == B else L(v, ATOM)
                return X(f"({L(a, ATOM)}.set! ({i.lean}).toNat {vt})",
                         f"{T(a, ATOM)}\\left[{i.tex} \\mapsto {v.tex}\\right]", ATOM, a.t)
            v = coerce(args[1], et) if args[1].t != B else args[1]
            return X(f"({L(a, ATOM)}.push {as_bool(v) if et == B else L(v, ATOM)})",
                     f"{T(a, ATOM)} \\mathbin{{+\\!\\!+}} \\left[{v.tex}\\right]", ATOM, a.t)
        if fname == "len" and len(args) == 1 and args[0].t in ARR:
            return X(f"(({args[0].lean}).size : ℤ)", f"\\left|{args[0].tex}\\right|", ATOM, Z)
        if fname == "list" and len(args) == 1 and args[0].t in ARR:
            return args[0]                                      # a copy: Lean arrays are values
        if isinstance(n.func, ast.Name) and fname in FUNCS:                       # another translated function
            sig = FUNCS[fname]
            if len(args) != len(sig.params):
                raise Refuse(f"'{fname}' takes {len(sig.params)} argument(s), got {len(args)}", n)
            parts = []
            for (pn, pt, _), x in zip(sig.params, args):
                if pt in ARR or x.t in ARR:
                    if x.t != pt:
                        raise Refuse(f"argument '{pn}' of '{fname}' needs a {LT[pt]}", n)
                elif pt == B:
                    if x.t != B:
                        raise Refuse(f"argument '{pn}' of '{fname}' needs a bool", n)
                    x = X(x.lean, x.tex, x.prec, B, x.bl)
                else:
                    x = coerce(x, pt)
                parts.append(as_bool(x) if pt == B else L(x, ATOM))
            call = f"{lid(fname)} " + " ".join(parts) if parts else lid(fname)
            tx = f"\\mathrm{{{fname.replace('_', chr(92) + '_')}}}\\left({', '.join(x.tex for x in args)}\\right)"
            if sig.ret == B:
                return X(f"(({call}) = true)", tx, ATOM, B, f"({call})")
            return X(call, tx, APP, sig.ret)
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
        if isinstance(n.func, ast.Name) and fname in FUNCS_REFUSED:
            raise Refuse(f"calls '{fname}', which was refused: {FUNCS_REFUSED[fname]}", n)
        raise Refuse(f"call to '{fname}' (supported: abs min max sqrt exp log sin cos tan factorial comb len list, "
                     "and the other functions of the same source)", n)
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


def _each(x):
    return x if isinstance(x, list) else [x]


def assigned(x) -> set:
    """Names a statement (list) may change: assigned, element-assigned or appended to."""
    out = set()
    for st in _each(x):
        for n in ast.walk(st):
            if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Store):
                out.add(n.id)
            elif isinstance(n, ast.Subscript) and isinstance(n.ctx, ast.Store):
                out |= {m.id for m in ast.walk(n.value) if isinstance(m, ast.Name)}
            elif isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute) and n.func.attr == "append":
                out |= {m.id for m in ast.walk(n.func.value) if isinstance(m, ast.Name)}
    return out


def loaded(x) -> set:
    return {n.id for st in _each(x) for n in ast.walk(st) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)}


class _ContFix(ast.NodeTransformer):
    """In a desugared ``for``: ``continue`` must still advance the loop variable."""
    def __init__(self, inc):
        self.inc = inc

    def visit_While(self, n):
        return n

    def visit_For(self, n):
        return n

    def visit_Continue(self, n):
        return [copy.deepcopy(self.inc), n]


def _is_assert_group(s) -> bool:
    return isinstance(s, ast.Assert) or is_decreases(s)


def desugar(stmts: list, later: frozenset = frozenset()) -> list:
    """``for i in range(a, b, c)`` -> ``i = a; while i < b: body; i += c`` and ``for x in a`` -> an index loop
    (recursively). The asserts written before the ``for`` become the asserts directly before the ``while``."""
    out = []
    for k, s in enumerate(stmts):
        rest = frozenset(later | loaded(stmts[k + 1:]))
        for field_ in ("body", "orelse"):
            if hasattr(s, field_) and isinstance(getattr(s, field_), list):
                setattr(s, field_, desugar(getattr(s, field_), rest))
        if not isinstance(s, ast.For):
            out.append(s)
            continue
        it = s.iter
        if not isinstance(s.target, ast.Name) or s.orelse:
            raise Refuse("only `for i in range(...)` and `for x in a` loops are supported", s)
        i = s.target.id
        by_range = isinstance(it, ast.Call) and getattr(it.func, "id", "") == "range" and 1 <= len(it.args) <= 3
        if by_range:
            a = it.args
            lo, hi = (ast.Constant(0), a[0]) if len(a) == 1 else (a[0], a[1])
            step = a[2] if len(a) == 3 else ast.Constant(1)
            if not (isinstance(step, ast.Constant) and isinstance(step.value, int) and step.value != 0):
                raise Refuse("range step must be a non-zero int literal", s)
            if (loaded(hi) | loaded(lo)) & (assigned(s.body)):
                raise Refuse("the range bound is changed inside the loop (Python evaluates range(...) once)", s)
            if i in assigned(s.body):
                raise Refuse(f"the loop variable '{i}' is assigned inside the loop", s)
            if i in rest:
                raise Refuse(f"the loop variable '{i}' is used after the loop (Python leaves the last value, "
                             "the while form the bound)", s)
            idx, init, bind = i, lo, []
            test = ast.Compare(ast.Name(i, ast.Load()), [ast.Lt() if step.value > 0 else ast.Gt()], [hi])
        elif isinstance(it, ast.Name):
            if it.id in assigned(s.body):
                raise Refuse(f"the list '{it.id}' is changed inside the loop that iterates over it", s)
            idx, step, init = f"idx_{i}", ast.Constant(1), ast.Constant(0)
            test = ast.Compare(ast.Name(idx, ast.Load()), [ast.Lt()],
                               [ast.Call(ast.Name("len", ast.Load()), [ast.Name(it.id, ast.Load())], [])])
            bind = [ast.copy_location(ast.Assign([ast.Name(i, ast.Store())], ast.Subscript(
                ast.Name(it.id, ast.Load()), ast.Name(idx, ast.Load()), ast.Load())), s)]
        else:
            raise Refuse("only `for i in range(...)` and `for x in a` loops are supported", s)
        group = []
        while out and _is_assert_group(out[-1]):
            group.insert(0, out.pop())
        inc = ast.copy_location(ast.AugAssign(ast.Name(idx, ast.Store()), ast.Add(), step), s)
        fix = _ContFix(inc)
        body = [r for st in s.body for r in _each(fix.visit(st))]
        out.append(ast.copy_location(ast.Assign([ast.Name(idx, ast.Store())], init), s))
        out.extend(group)
        w = ast.copy_location(ast.While(test, bind + body + [inc], []), s)
        ast.fix_missing_locations(w)
        out.append(w)
    return out


def is_decreases(s) -> bool:
    return isinstance(s, ast.Expr) and isinstance(s.value, ast.Call) and getattr(s.value.func, "id", "") == "decreases"


# ---- statement translation (Lean do-block + LaTeX rows) ----------------------------------------
def assign_parts(s):
    """(targets, values, declared type) of an assignment statement; ``[]`` with a declared list type becomes
    ``__empty('AZ')``. Targets are Names or ``name[index]``."""
    decl = None
    if isinstance(s, ast.AugAssign):
        tgts = [s.target]
        load = copy.deepcopy(s.target)
        load.ctx = ast.Load()
        vals = [ast.copy_location(ast.BinOp(load, s.op, s.value), s)]
    elif isinstance(s, ast.AnnAssign):
        if s.value is None or not isinstance(s.target, ast.Name):
            raise Refuse("annotated assignment needs a name and a value", s)
        tgts, vals, decl = [s.target], [s.value], annotation_type(s.annotation, s)[0]
    else:
        if len(s.targets) != 1:
            raise Refuse("chained assignment", s)
        if isinstance(s.targets[0], ast.Tuple):
            tgts = list(s.targets[0].elts)
            if not (isinstance(s.value, ast.Tuple) and len(s.value.elts) == len(tgts)):
                raise Refuse("tuple assignment needs a tuple of the same length on the right", s)
            vals = list(s.value.elts)
        else:
            tgts, vals = [s.targets[0]], [s.value]
    if not all(isinstance(t, ast.Name) or (isinstance(t, ast.Subscript) and isinstance(t.value, ast.Name))
               for t in tgts):
        raise Refuse("assignment target must be a name or a list element", s)
    name_tgts = {t.id for t in tgts if isinstance(t, ast.Name)}
    for t in tgts:
        if isinstance(t, ast.Subscript) and len(tgts) > 1 and (names(t.slice) | {t.value.id}) & name_tgts:
            raise Refuse("a list element target whose index or list is assigned in the same tuple assignment", s)
    vals = [ast.Call(ast.Name("__empty", ast.Load()), [ast.Constant(decl)], []) if isinstance(v, ast.List)
            and not v.elts and decl in ARR else v for v in vals]
    return tgts, vals, decl


class Gen:
    def __init__(self, fn: "FnInfo", ret=None):
        self.fn, self.notes = fn, fn.notes
        self.ret = ret                        # declared return type, or None
        self.lean: list[str] = []
        self.rows: list[tuple[int, str]] = []
        self.snap: dict[int, dict] = {}       # id(While) -> type env at that point
        self.returns: list[X] = []
        self.depth = 0                        # loop nesting, for break/continue

    def emit(self, ind, lean, tex):
        if lean is not None:
            self.lean.append("  " * ind + lean)
        if tex is not None:
            self.rows.append((ind, tex))

    def branch(self, stmts, env, ind):
        """A nested block; a Lean block with only comments in it would not parse, so it gets ``pure ()``."""
        n0 = len(self.lean)
        self.block(stmts, env, ind)
        if all(l.lstrip().startswith("--") for l in self.lean[n0:]):
            self.emit(ind, "pure ()", None)

    def block(self, stmts, env, ind, top=False):
        env = dict(env)
        pending: list[ast.stmt] = []          # asserts/decreases before a loop
        for s in stmts:
            if isinstance(s, ast.Assert) or is_decreases(s):
                pending.append(s)
                continue
            self.stmt(s, env, ind, pending)
            pending = []
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
        elif isinstance(s, ast.Expr) and isinstance(s.value, ast.Call) \
                and isinstance(s.value.func, ast.Attribute) and s.value.func.attr == "append":
            self.append(s, env, ind)
        elif isinstance(s, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
            self.assign(s, env, ind)
        elif isinstance(s, ast.If):
            c = ex(s.test, env, self.notes)
            if c.t != B:
                raise Refuse("`if` condition is not a comparison", s)
            self.emit(ind, f"if {c.lean} then", f"\\mathbf{{if}}\\ {c.tex}\\ \\mathbf{{then}}")
            self.branch(s.body, env, ind + 1)
            if s.orelse:
                self.emit(ind, "else", "\\mathbf{else}")
                self.branch(s.orelse, env, ind + 1)
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
            self.depth += 1
            self.branch(s.body, env, ind + 1)
            self.depth -= 1
        elif isinstance(s, ast.Return):
            if s.value is None:
                raise Refuse("bare return", s)
            x = ex(s.value, env, self.notes)
            d = self.ret
            if d is not None:
                if d in ARR or d == B or x.t in ARR or x.t == B:
                    if x.t != d:
                        raise Refuse(f"returns a {LT[x.t]} but the declared return type is {LT[d]}", s)
                else:
                    x = coerce(x, d)
            self.returns.append(x)
            self.emit(ind, f"return {as_bool(x) if x.t == B else x.lean}", f"\\mathbf{{return}}\\ {x.tex}")
        elif isinstance(s, (ast.Break, ast.Continue)):
            if self.depth == 0:
                raise Refuse("break/continue outside a loop", s)
            w = "break" if isinstance(s, ast.Break) else "continue"
            self.emit(ind, w, f"\\mathbf{{{w}}}")
        else:
            raise Refuse(f"{type(s).__name__} statement", s)

    def check_mutable(self, name, node):
        if name in self.fn.list_params:
            raise Refuse(f"the list argument '{name}' is changed in place: Lean arrays are values and the caller "
                         "would not see it; work on a copy (`b = list(a)`)", node)

    def append(self, s, env, ind):
        c = s.value
        if not (isinstance(c.func.value, ast.Name) and c.func.value.id in env and env[c.func.value.id] in ARR
                and len(c.args) == 1 and not c.keywords):
            raise Refuse("only `a.append(v)` on a list variable is supported", s)
        a = c.func.value.id
        self.check_mutable(a, s)
        x = ex(ast.Call(ast.Name("__push", ast.Load()), [ast.Name(a, ast.Load()), c.args[0]], []), env, self.notes)
        v = ex(c.args[0], env, self.notes)
        self.emit(ind, f"{lid(a)} := {x.lean}", f"{tid(a)} \\gets {tid(a)} \\mathbin{{+\\!\\!+}} \\left[{v.tex}\\right]")

    def assign(self, s, env, ind):
        tgts, vals, decl = assign_parts(s)
        xs = []
        for v in vals:
            x = ex(v, env, self.notes)
            if x.t in ARR and isinstance(v, ast.Name):
                raise Refuse(f"`{v.id}` is a list: assigning it would alias it in Python but copy it in Lean; "
                             "use list(...)", s)
            xs.append(x)
        new = [isinstance(t, ast.Name) and t.id not in env for t in tgts]
        if any(new) and not all(new):
            raise Refuse("tuple assignment mixing new and existing variables", s)
        rhs: dict[str, tuple[str, str]] = {}      # var -> (lean text, type)
        chain: dict[str, str] = {}
        rows = []
        for t, x in zip(tgts, xs):
            if isinstance(t, ast.Name):
                want = env.get(t.id, decl)
                if want in ARR or want == B or x.t in ARR or x.t == B:
                    if want is not None and x.t != want:
                        raise Refuse(f"'{t.id}' is a {LT[want]} and is assigned a {LT[x.t]}", s)
                elif want is not None:
                    if want == Z and x.t == R:
                        raise Refuse(f"'{t.id}' is an int and is assigned a float (Lean variables have one type)", s)
                    x = coerce(x, want)
                rhs[t.id] = (as_bool(x) if x.t == B else x.lean, x.t)
                rows.append((tid(t.id), x.tex))
            else:
                a = t.value.id
                if a not in env or env[a] not in ARR:
                    raise Refuse(f"'{a}' is not a list", s)
                self.check_mutable(a, s)
                i = ex(t.slice, env, self.notes)
                if i.t != Z:
                    raise Refuse("a list index must be an int", s)
                et = ARR[env[a]]
                if et == B or x.t == B:
                    if x.t != et:
                        raise Refuse(f"a {LT[et]} list element is assigned a {LT[x.t]}", s)
                else:
                    if et == Z and x.t == R:
                        raise Refuse("a float is stored in an int list", s)
                    x = coerce(x, et)
                base = chain.get(a, lid(a))
                base = f"({base})" if a in chain else base
                chain[a] = f"{base}.set! ({i.lean}).toNat {as_bool(x) if et == B else L(x, ATOM)}"
                rhs[a] = (chain[a], env[a])
                self.notes.append("a[i] = v prints as a.set! i.toNat v: out-of-range indices are outside the model")
                rows.append((f"{tid(a)}\\left[{i.tex}\\right]", x.tex))
        if len(rhs) == 1:
            (v, (txt, ty)), = rhs.items()
            if all(new):
                env[v] = ty
                self.emit(ind, f"let mut {lid(v)} : {LT[ty]} := {txt}", f"{rows[0][0]} \\gets {rows[0][1]}")
            else:
                self.emit(ind, f"{lid(v)} := {txt}", f"{rows[0][0]} \\gets {rows[0][1]}")
        else:
            lhs = ", ".join(lid(v) for v in rhs)
            rtxt = ", ".join(f"({txt} : {LT[ty]})" for txt, ty in rhs.values())
            if all(new):
                for v, (_, ty) in rhs.items():
                    env[v] = ty
            self.emit(ind, f"{'let mut ' if all(new) else ''}({lhs}) := ({rtxt})",
                      f"({', '.join(r[0] for r in rows)}) \\gets ({', '.join(r[1] for r in rows)})")


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
    vc_data: list = field(default_factory=list)     # (name, symbols, [hypothesis asts], [goal asts]) for --smoke
    list_params: set = field(default_factory=set)   # list arguments (never changed in place)
    calls: set = field(default_factory=set)         # other translated functions this one calls


def annotation_type(a, node):
    nm = a.id if isinstance(a, ast.Name) else (a.attr if isinstance(a, ast.Attribute) else None)
    if nm == "int":
        return Z, False
    if nm == "Nat":
        return Z, True
    if nm == "float":
        return R, False
    if nm == "bool":
        return B, False
    if isinstance(a, ast.Subscript) and getattr(a.value, "id", None) in ("list", "List") \
            and isinstance(a.slice, ast.Name) and a.slice.id in ("int", "float", "bool"):
        return {"int": AZ, "float": AR, "bool": AB}[a.slice.id], False
    raise Refuse("every argument needs an annotation: int, Nat, float, bool, list[int], list[float] or list[bool]", node)


def group_binders(vars_: list) -> str:
    """Lean binders in the given order; only neighbours of the same type share a group, so that the order of the
    variables is the order of the arguments in a call."""
    groups: list[tuple[str, list[str]]] = []
    for n, t in vars_:
        if groups and groups[-1][0] == LT[t]:
            groups[-1][1].append(lid(n))
        else:
            groups.append((LT[t], [lid(n)]))
    return " ".join(f"({' '.join(v)} : {t})" for t, v in groups)


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
        if t in ARR:
            info.list_params.add(a.arg)
    if fd.args.vararg or fd.args.kwarg or fd.args.defaults or fd.args.kwonlyargs:
        raise Refuse("defaults, *args and keyword-only arguments", fd)
    body = desugar(list(fd.body))
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
        body = body[1:]
    pre = []
    while body and isinstance(body[0], ast.Assert):
        pre.append(body.pop(0))
    info.pre = [p.test for p in pre]
    declared = annotation_type(fd.returns, fd)[0] if fd.returns is not None else None
    g = Gen(info, declared)
    for p in info.pre:
        g_x = ex(p, env, info.notes)
        g.emit(1, f"-- requires: {g_x.lean}", f"\\{{\\,{g_x.tex}\\,\\}}")
    g.block(body, env, 1, top=True)
    if not g.returns:
        raise Refuse("the function never returns a value", fd)
    info.ret = declared if declared is not None else g.returns[0].t
    info.calls = {n.func.id for st in fd.body for n in ast.walk(st)
                  if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id in FUNCS}
    arg_t = [(n, t) for n, t, _ in info.args]
    noncomp = "noncomputable " if R in [t for _, t in arg_t] + [info.ret] or any("ℝ" in l for l in g.lean) else ""
    info.lean_def = (f"{noncomp}def {lid(info.name)} {group_binders(arg_t)} : {LT[info.ret]} := Id.run do\n"
                     + "\n".join(g.lean) + "\n")
    info.rows = g.rows
    try:
        _vcs(info, g, body, env, prove)
    except Refuse as e:
        info.vc_refused = str(e)


# ---- verification conditions -------------------------------------------------------------------
@dataclass
class Fact:
    label: str
    e: object                # an ast expression over the symbols of its segment


@dataclass
class Seg:
    """A stretch of symbolic execution. ``E`` maps every variable in scope to an expression over the segment's
    symbols (= the variables in scope when the segment began); ``facts`` are what may be assumed."""
    E: dict
    types: dict
    syms: list
    facts: list

    def copy(self):
        return Seg(dict(self.E), dict(self.types), self.syms, list(self.facts))


def _call(name, *args):
    return ast.Call(ast.Name(name, ast.Load()), list(args), [])


def _cmp(a, op, b):
    return ast.Compare(a, [op], [b])


def _not(x):
    return ast.UnaryOp(ast.Not(), x)


def _and(xs):
    return xs[0] if len(xs) == 1 else ast.BoolOp(ast.And(), xs)


def syms_of(node) -> set:
    """Variable names in an expression (not function names)."""
    skip = {id(n.func) for n in ast.walk(node) if isinstance(n, ast.Call)}
    return {n.id for n in ast.walk(node) if isinstance(n, ast.Name) and id(n) not in skip
            and n.id not in BUILTIN_NAMES and n.id not in INVS}


class _Rename(ast.NodeTransformer):
    def __init__(self, m):
        self.m = m

    def visit_Call(self, n):
        n.args = [self.visit(a) for a in n.args]
        return n

    def visit_Name(self, n):
        return ast.Name(self.m.get(n.id, n.id), n.ctx)


def rename(node, m):
    return _Rename(m).visit(copy.deepcopy(node))


def conjuncts(e) -> list:
    """The conjuncts of a fact; an invariant instance is expanded into the loop's asserts."""
    if isinstance(e, ast.Call) and isinstance(e.func, ast.Name) and e.func.id in INVS:
        _, svars, inv_asts = INVS[e.func.id]
        mp = {v: a for (v, _), a in zip(svars, e.args)}
        return [c for i in inv_asts for c in conjuncts(subst(i, mp))]
    if isinstance(e, ast.BoolOp) and isinstance(e.op, ast.And):
        return [c for v in e.values for c in conjuncts(v)]
    return [e]


def assert_tests(group) -> list:
    return [g.test for g in group if isinstance(g, ast.Assert)]


def group_items(stmts) -> list:
    """Split a block into ('stmt', s), ('loop', while, asserts-before-it) and ('asserts', group) items."""
    items, buf = [], []
    for s in stmts:
        if _is_assert_group(s):
            buf.append(s)
        elif isinstance(s, ast.While):
            items.append(("loop", s, buf))
            buf = []
        else:
            if buf:
                items.append(("asserts", buf))
                buf = []
            items.append(("stmt", s))
    if buf:
        items.append(("asserts", buf))
    return items


class VCGen:
    def __init__(self, info: "FnInfo", g: Gen, prove: bool, total: int):
        self.info, self.g, self.prove, self.total = info, g, prove, total
        self.notes = info.notes
        self.fn = lid(info.name)
        self.vcs: list = []              # (name, lean text, tex text)
        self.data: list = []
        self.seen: set = set()
        self.names: set = set()
        self.counter: dict = {}
        self.k = 0                       # loops seen so far
        self.inv_names: list = []
        self.dropped: list = []

    # -- naming ----------------------------------------------------------------------------------
    def fresh(self, what: str) -> int:
        self.counter[what] = self.counter.get(what, 0) + 1
        return self.counter[what]

    def uniq(self, name: str) -> str:
        n, i = name, 2
        while n in self.names:
            n, i = f"{name}_{i}", i + 1
        self.names.add(n)
        return n

    def lname(self, k: int) -> str:
        return self.fn if self.total == 1 else f"{self.fn}_L{k}"

    def lab(self, what: str, k: int) -> str:
        return what if self.total == 1 else f"{what}{k}"

    def tac(self) -> str:
        if not self.prove:
            return "by\n  sorry"
        unfold = f"try simp only [{', '.join(self.inv_names)}] at *\n  " if self.inv_names else ""
        return f"by\n  intros\n  {unfold}first | omega | nlinarith | sorry"

    # -- theorem text ----------------------------------------------------------------------------
    def emit(self, name: str, S: Seg, goals: list, dedupe: bool = False) -> bool:
        env = dict(S.syms)
        hyps = [(f.label, ex(f.e, env, self.notes)) for f in S.facts if not f.label.startswith("__")]
        gx = [ex(x, env, self.notes) for x in goals]
        parts = [group_binders(S.syms)] + [f"({lb} : {x.lean})" for lb, x in hyps]
        text = (f"theorem {name} {' '.join(p for p in parts if p)} :\n"
                f"    {' ∧ '.join(L(x, AND if len(gx) == 1 else AND + 1) for x in gx)} := {self.tac()}\n")
        key = text.replace(name, "", 1)
        if dedupe and key in self.seen:
            return False
        self.seen.add(key)
        ctx = [x for lb, x in hyps if not lb.startswith(("hnat", "hpre"))]
        lhs = " \\land ".join(T(x, AND + 1) for x in ctx) if ctx else "\\text{pre}"
        self.data.append((name, list(S.syms), [f.e for f in S.facts if not f.label.startswith("__")], list(goals)))
        self.vcs.append((name, text, f"{lhs} \\Rightarrow " + " \\land ".join(T(x, AND if len(gx) == 1 else AND + 1) for x in gx)))
        return True

    # -- segments --------------------------------------------------------------------------------
    def cut(self, S: Seg, M: set, scope: list, extra: list) -> Seg:
        """A new segment whose symbols are the variables of ``scope``. A fact of the old segment is carried over
        (renamed) only if every variable it mentions still has the value it had, i.e. is not in ``M``."""
        back: dict = {}
        for v, _ in scope:
            e = S.E.get(v)
            if isinstance(e, ast.Name) and v not in M and e.id == v:
                back[e.id] = v
        for v, _ in scope:
            e = S.E.get(v)
            if isinstance(e, ast.Name) and v not in M:
                back.setdefault(e.id, v)
        facts = []
        for f in S.facts:
            if all(x in back for x in syms_of(f.e)):
                facts.append(Fact(f.label, rename(f.e, back)))
                continue
            parts = conjuncts(f.e)                              # keep the conjuncts that are still true
            keep = [(i, c) for i, c in enumerate(parts, 1) if all(x in back for x in syms_of(c))]
            facts += [Fact(f"{f.label}_{i}", rename(c, back)) for i, c in keep]
            if len(keep) < len(parts):
                self.dropped.append(f.label)
        return Seg({v: ast.Name(v, ast.Load()) for v, _ in scope}, dict(scope), list(scope), facts + extra)

    # -- side conditions: list indices and calls to other functions ---------------------------------
    def bound(self, a, i):
        return _and([_cmp(ast.Constant(0), ast.LtE(), i), _cmp(i, ast.Lt(), _call("len", a))])

    def safe(self, node, S: Seg):
        conds = []
        for n in ast.walk(node):
            if isinstance(n, ast.Subscript) and isinstance(n.ctx, ast.Load):
                conds.append(self.bound(n.value, n.slice))
            elif isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id in FUNCS:
                sig = FUNCS[n.func.id]
                mp = {p[0]: a for p, a in zip(sig.params, n.args)}
                conds += [_cmp(ast.Constant(0), ast.LtE(), mp[p[0]]) for p in sig.params if p[2] and p[0] in mp]
                conds += [subst(pr, mp) for pr in sig.pre]
        self.safe_conds(conds, S)

    def safe_conds(self, conds, S: Seg):
        seen, out = set(), []
        for c in conds:
            if ast.dump(c) not in seen:
                seen.add(ast.dump(c))
                out.append(subst(c, S.E))
        if out:
            name = self.uniq(f"{self.fn}_safe{self.fresh('safe')}")
            if not self.emit(name, S, out, dedupe=True):        # an identical theorem exists: reuse its number
                self.names.discard(name)
                self.counter["safe"] -= 1

    # -- execution -------------------------------------------------------------------------------
    def block(self, stmts, S: Seg) -> list:
        return self.items(group_items(stmts), S)

    def items(self, items, S: Seg) -> list:
        states, exits, k = [S], [], 0
        while k < len(items):
            it = items[k]
            k += 1
            post = None
            if it[0] == "loop" and k < len(items) and items[k][0] == "asserts":
                post = items[k][1]
                k += 1
            nxt = []
            for st in states:
                if it[0] == "stmt":
                    res = self.stmt(it[1], st)
                elif it[0] == "asserts":
                    res = self.asserts(it[1], st)
                else:
                    res = self.loop(it[1], it[2], post, st)
                for kind, s2 in res:
                    (nxt if kind == "next" else exits).append(s2 if kind == "next" else (kind, s2))
            states = nxt
        return [("next", s) for s in states] + exits

    def asserts(self, group, S: Seg) -> list:
        tests = assert_tests(group)
        if not tests:
            return [("next", S)]
        goals = [subst(t, S.E) for t in tests]
        j = self.fresh("assert")
        self.emit(self.uniq(f"{self.fn}_assert{j}"), S, goals)
        S2 = S.copy()
        S2.facts += [Fact(f"hassert{j}_{i}", gl) for i, gl in enumerate(goals)]
        return [("next", S2)]

    def stmt(self, s, S: Seg) -> list:
        if isinstance(s, ast.Pass) or (isinstance(s, ast.Expr) and isinstance(s.value, ast.Constant)):
            return [("next", S)]
        if isinstance(s, ast.Expr) and isinstance(s.value, ast.Call):          # a.append(v)
            a, arg = s.value.func.value.id, s.value.args[0]
            self.safe(arg, S)
            S2 = S.copy()
            S2.E[a] = _call("__push", S.E[a], subst(arg, S.E))
            return [("next", S2)]
        if isinstance(s, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
            return [("next", self.assign(s, S))]
        if isinstance(s, ast.If):
            return self.if_(s, S)
        if isinstance(s, ast.Return):
            self.safe(s.value, S)
            return [("return", S)]
        if isinstance(s, ast.Break):
            return [("break", S)]
        if isinstance(s, ast.Continue):
            return [("continue", S)]
        raise Refuse(f"{type(s).__name__} statement is not handled by the VC generator", s)

    def assign(self, s, S: Seg) -> Seg:
        tgts, vals, decl = assign_parts(s)
        for v in vals:
            self.safe(v, S)
        for t in tgts:
            if isinstance(t, ast.Subscript):
                self.safe(t.slice, S)
                self.safe_conds([self.bound(t.value, t.slice)], S)
        S2, new = S.copy(), {}
        for t, v in zip(tgts, vals):
            vv = subst(v, S.E)                                  # parallel: every right side sees the old state
            if isinstance(t, ast.Name):
                new[t.id] = vv
                if t.id not in S2.types:
                    S2.types[t.id] = decl or ex(v, S.types, self.notes).t
            else:
                a = t.value.id
                new[a] = _call("__set", new.get(a, S.E[a]), subst(t.slice, S.E), vv)
        S2.E.update(new)
        return S2

    def if_(self, s, S: Seg) -> list:
        self.safe(s.test, S)
        c = subst(s.test, S.E)
        branches = s.body + s.orelse
        if not any(isinstance(n, (ast.Break, ast.Continue, ast.Return, ast.While)) for b in branches for n in ast.walk(b)):
            outs = []
            for cond, stmts in ((c, s.body), (_not(copy.deepcopy(c)), s.orelse)):
                S1 = S.copy()
                S1.facts.append(Fact(f"hpath{self.fresh('hpath')}", cond))
                r = self.block(stmts, S1)
                outs.append(r[0][1])
            S3 = S.copy()                                       # merge: an if-expression for what differs
            for key in S.types:
                a, b = outs[0].E.get(key), outs[1].E.get(key)
                if a is not None and b is not None and ast.dump(a) != ast.dump(b):
                    S3.E[key] = ast.IfExp(copy.deepcopy(c), a, b)
            return [("next", S3)]
        res = []
        for cond, stmts in ((c, s.body), (_not(copy.deepcopy(c)), s.orelse)):
            S1 = S.copy()
            S1.facts.append(Fact(f"hpath{self.fresh('hpath')}", cond))
            res += self.block(stmts, S1)
        return res

    # -- loops -----------------------------------------------------------------------------------
    def loop(self, w, group, post, S: Seg) -> list:
        self.k += 1
        k = self.k
        base = self.lname(k)
        invs = assert_tests(group)
        variants = [g.value.args[0] for g in group if is_decreases(g)]
        if not invs:
            raise Refuse(f"the loop at line {w.lineno} has no invariant: put `assert ...` immediately before it")
        scope = list(self.g.snap[id(w)].items())
        senv = dict(scope)
        M = assigned(w.body)
        inv_id = f"__inv_{self.fn}_{k}"
        INVS[inv_id] = (f"{base}_inv", scope, invs)
        self.inv_names.append(f"{base}_inv")
        inv_x = [ex(i, senv, self.notes) for i in invs]
        self.vcs.append((f"{base}_inv", f"def {base}_inv {group_binders(scope)} : Prop :=\n  "
                         + " ∧ ".join(L(x, AND + 1) for x in inv_x) + "\n", None))
        for v in variants:
            if ex(v, senv, self.notes).t != Z:
                raise Refuse("decreases(...) must be an int expression", w)
        ident = _call(inv_id, *[ast.Name(v, ast.Load()) for v, _ in scope])
        self.emit(self.uniq(f"{base}_init"), S, [_call(inv_id, *[S.E.get(v, ast.Name(v, ast.Load())) for v, _ in scope])])
        hinv, hguard, hexit = self.lab("hinv", k), self.lab("hguard", k), self.lab("hexit", k)
        Sb = self.cut(S, M, scope, [Fact(hinv, ident), Fact(hguard, w.test)]
                      + [Fact(f"__var{k}_{vi}", v) for vi, v in enumerate(variants)])
        Sg = Sb.copy()
        Sg.facts = [f for f in Sg.facts if f.label != hguard]
        self.safe(w.test, Sg)                                   # the guard is evaluated under the invariant only
        res = self.block(w.body, Sb)
        nexts = [st for kind, st in res if kind in ("next", "continue")]
        breaks = [st for kind, st in res if kind == "break"]
        pv, evc, vv = [], [], []
        for j, St in enumerate(nexts):
            pn = f"_path{j + 1}" if len(nexts) > 1 else ""
            pv.append((self.uniq(f"{base}_preserve{pn}"), St,
                       [_call(inv_id, *[St.E.get(v, ast.Name(v, ast.Load())) for v, _ in scope])]))
            for vi, v in enumerate(variants):
                v0 = next((f.e for f in St.facts if f.label == f"__var{k}_{vi}"), None)
                if v0 is None:
                    raise Refuse("the decreases(...) expression mentions a variable that an inner loop changes")
                vv.append((self.uniq(f"{base}_variant{vi or ''}{pn}"), St,
                           [_and([_cmp(ast.Constant(0), ast.LtE(), copy.deepcopy(v0)),
                                  _cmp(subst(v, St.E), ast.Lt(), copy.deepcopy(v0))])]))
        extra = [Fact(hinv, ident), Fact(hexit, _not(copy.deepcopy(w.test)))]
        Sx = self.cut(S, M, scope, extra)                       # the normal exit: invariant and not guard
        Sa = Sx.copy() if not breaks else self.cut(S, M, scope, [])   # what later code may assume
        tests = assert_tests(post) if post else []
        if tests:
            evc.append((self.uniq(f"{base}_exit"), Sx, tests))
            for j, St in enumerate(breaks):
                evc.append((self.uniq(f"{base}_exit_break{j + 1}"), St, [subst(t, St.E) for t in tests]))
            Sa.facts += [Fact(f"hpost{k}_{i}", t) for i, t in enumerate(tests)]
        elif self.total == 1:
            self.notes.append("no assert after the loop: no exit condition was generated")
        for name, St, goals in pv + evc + vv:
            self.emit(name, St, goals)
        return [("next", Sa)]


def _vcs(info: FnInfo, g: Gen, body: list, env0: dict, prove: bool):
    total = sum(isinstance(n, ast.While) for st in body for n in ast.walk(st))
    types = dict(env0)
    facts = [Fact(f"hnat_{lid(n)}", _cmp(ast.Constant(0), ast.LtE(), ast.Name(n, ast.Load())))
             for n, _, nat in info.args if nat]
    facts += [Fact(f"hpre{i}", p) for i, p in enumerate(info.pre)]
    S0 = Seg({n: ast.Name(n, ast.Load()) for n in types}, types, list(types.items()), facts)
    vg = VCGen(info, g, prove, total)
    res = vg.block(body, S0)
    if any(kind in ("break", "continue") for kind, _ in res):
        raise Refuse("break/continue outside a loop")
    if not vg.vcs:
        raise Refuse("nothing to check (no loop, assert after the precondition, list index or call)")
    dropped = sorted({d for d in vg.dropped if not d.startswith("__")})
    if dropped:
        info.notes.append("parts of these context assumptions are not carried into the verification conditions of "
                          f"a later loop because the variables they mention change ({', '.join(dropped)}): state "
                          "what is needed in the invariant")
    info.vcs, info.vc_data = vg.vcs, vg.data
    first = INVS.get(f"__inv_{vg.fn}_1")
    if first:
        info.inv_tex = " \\land ".join(T(ex(i, dict(first[1]), []), AND + 1) for i in first[2])


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
                lines.append(f"\\[ \\text{{{nm[len(lid(info.name)) + 1:].replace('_', chr(92) + '_')}:}}\\quad {tx} \\]")
    elif info.vc_refused:
        lines.append(f"No verification conditions: {info.vc_refused.replace('_', chr(92) + '_')}.")
    return "\n".join(lines) + "\n"


def translate_source(src: str, prove: bool = False) -> list[FnInfo]:
    """Translate every top-level function. Callees come first (Lean needs the definition before its use);
    recursive functions are refused, since a Lean ``def`` of them would need a termination proof."""
    tree = ast.parse(textwrap.dedent(src))
    fds = [n for n in tree.body if isinstance(n, ast.FunctionDef)]
    by_name = {f.name: f for f in fds}
    FUNCS.clear()
    FUNCS_REFUSED.clear()
    callees = {f.name: [n.func.id for n in ast.walk(f) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
                        and n.func.id in by_name] for f in fds}
    state, order, cyclic = {}, [], set()

    def visit(name, stack):
        if state.get(name) == 2:
            return
        if state.get(name) == 1:                         # back edge: everything on the stack from here is a cycle
            cyclic.update(stack[stack.index(name):])
            return
        state[name] = 1
        for c in callees[name]:
            visit(c, stack + [name])
        state[name] = 2
        order.append(name)

    for f in fds:
        visit(f.name, [])
    infos = []
    for name in order:
        fd = by_name[name]
        if name in cyclic:
            info = FnInfo(name, fd)
            info.refused = "recursive or mutually recursive function (the Lean def would need a termination proof)"
        else:
            info = translate_function(fd, prove)
        if info.refused:
            FUNCS_REFUSED[name] = info.refused
        else:
            FUNCS[name] = FuncSig(name, list(info.args), info.ret, list(info.pre))
        infos.append(info)
    return infos


# ---- numeric smoke test of the verification conditions ------------------------------------------
class _ToPy(ast.NodeTransformer):
    """Plain Python for a VC expression: invariant instances are expanded, list updates become helper calls."""
    def visit_Call(self, n):
        n = self.generic_visit(n)
        if isinstance(n.func, ast.Name) and n.func.id in INVS:
            _, svars, inv_asts = INVS[n.func.id]
            body = [subst(i, {v: a for (v, _), a in zip(svars, n.args)}) for i in inv_asts]
            return ast.copy_location(self.visit(ast.BoolOp(ast.And(), body) if len(body) > 1 else body[0]), n)
        if isinstance(n.func, ast.Name) and n.func.id in ("__set", "__push"):
            n.func.id = "_" + n.func.id.strip("_")
        if isinstance(n.func, ast.Name) and n.func.id == "__empty":
            return ast.List([], ast.Load())
        return n


def _smoke_values(t, size):
    if t == Z:
        return list(range(-2, -2 + size))
    if t == R:
        return [-2, -1, -0.5, 0, 0.5, 1, 2, 3]
    if t == B:
        return [False, True]
    et = [v for v in _smoke_values(ARR[t], 16) if v in (-1, 0, 1, 2, False, True, 0.5)]
    out = [[]]
    for ln in (1, 2, 3):
        out += [list(c) for c in itertools.product(et, repeat=ln)]
    return out


def smoke(src: str, infos: list, limit: int = 40000) -> list:
    """Evaluate every generated verification condition on a grid of small values. Returns
    (function, theorem, samples with the hypotheses true, counterexample or None). A sanity check, not a proof."""
    import itertools, math, random
    rng = random.Random(1)
    ns = {"Nat": int, "decreases": lambda *a: None, "math": math, "_set": lambda a, i, v: a[:i] + [v] + a[i + 1:],
          "_push": lambda a, v: a + [v], "sqrt": math.sqrt, "exp": math.exp, "log": math.log, "sin": math.sin,
          "cos": math.cos, "tan": math.tan, "factorial": math.factorial, "comb": math.comb, "pi": math.pi}
    exec(compile(textwrap.dedent(src), "<smoke>", "exec"), ns)
    out = []
    for info in infos:
        for name, syms, hyps, goals in info.vc_data:
            def code(e):
                tree = _ToPy().visit(copy.deepcopy(e))
                for nd in ast.walk(tree):
                    for at in ("lineno", "end_lineno", "col_offset", "end_col_offset"):
                        nd.__dict__.pop(at, None)
                return compile(ast.fix_missing_locations(ast.Expression(tree)), "<vc>", "eval")
            hc, gc = [code(h) for h in hyps], [code(g) for g in goals]
            size = max(6, min(16, int(limit ** (1 / max(1, len(syms))))))       # fewer values when there are many
            doms = [_smoke_values(t, size) for _, t in syms]
            total = 1
            for d in doms:
                total *= len(d)
            combos = (itertools.product(*doms) if total <= limit
                      else (tuple(rng.choice(d) for d in doms) for _ in range(limit)))
            sat, bad = 0, None
            for vals in combos:
                env = dict(ns)
                env.update({n: (list(v) if isinstance(v, list) else v) for (n, _), v in zip(syms, vals)})
                try:
                    if not all(eval(h, env) for h in hc):
                        continue
                    sat += 1
                    if not all(eval(g, env) for g in gc):
                        bad = bad or dict(zip((n for n, _ in syms), vals))
                except (AssertionError, ArithmeticError, IndexError, ValueError, TypeError, RecursionError):
                    continue
            out.append((info.name, name, sat, bad))
    return out


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
    for src, why in (("def f(x: int) -> int:\n    break\n    return x\n", "outside a loop"),
                     ("def f(x) -> int:\n    return x\n", "no annotation"),
                     ("def f(x: int) -> int:\n    return [x]\n", "declared return type"),
                     ("def f(x: int) -> int:\n    return (x, x)\n", "Tuple expression"),
                     ("def f(x: int) -> int:\n    y = 1\n    y = 1.5\n    return x\n", "assigned a float")):
        r = translate_source(src)[0]
        ok(why in r.refused, (src, r.refused))
    ok("line 3" in translate_source("def f(x: int) -> int:\n    y = 1\n    y = 1.5\n    return x\n")[0].refused)
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
    import pathlib
    models = pathlib.Path(__file__).resolve().parent / "py_models"
    by = {}
    for f in sorted(models.glob("*.py")):
        for i in translate_source(f.read_text()):
            by[i.name] = i
    # nested and several loops: one invariant, init, preserve, variant per loop, in the right context
    t = lean_text(by["nested"])
    ok("theorem nested_L1_init" in t and "theorem nested_L2_init" in t and "theorem nested_L2_preserve" in t, t)
    ok("(hinv1_2 : 0 ≤ i)" in t and "(hinv2 : nested_L2_inv n s i j)" in t, t)           # the still-true part of hinv1
    ok("theorem nested_L1_preserve" in t and "nested_L1_inv n s (i + 1)" in t, t)
    t = lean_text(by["two_phase"])
    ok("two_phase_L1_exit" in t and "two_phase_L2_init" in t and "(hpost1_0 : a = n)" in t, t)
    # break / continue: Lean keeps them, the exit must hold on every way out, continue is its own path
    t = lean_text(by["first_neg"])
    ok("      break" in t and "theorem first_neg_exit_break1" in t and "(hexit : ¬(k < " in t, t)
    t = lean_text(by["count_odd"])
    ok("      continue" in t and "count_odd_preserve_path1" in t and "count_odd_preserve_path2" in t, t)
    # a for loop with continue still advances its variable; a branch with only a comment gets pure ()
    fc = translate_source("def f(n: Nat) -> int:\n  t = 0\n  for i in range(n):\n    if i == 2:\n      continue\n"
                          "    t += i\n  return t\n")[0]
    ok("i := i + 1\n      continue" in lean_text(fc), lean_text(fc))
    fp = translate_source("def f(x: int) -> int:\n  if x > 0:\n    assert x > 0\n  return x\n")[0]
    ok("pure ()" in lean_text(fp), lean_text(fp))
    # bool variables are Lean Bool, a Prop is decided when stored; binders keep the argument order
    t = lean_text(by["seen_big"])
    ok("let mut seen : Bool := false" in t and "(n : ℤ) (seen : Bool) (i : ℤ)" in t, t)
    ok("seen_big_inv n (decide (if i > 3 then True else (seen = true))) (i + 1)" in t, t)
    ok("def pick (flag : Bool) (x y : ℤ) : ℤ" in lean_text(by["pick"]), lean_text(by["pick"]))
    # calls: callee before caller, the callee's precondition is a theorem at the call
    t = lean_text(by["sum_sq"])
    ok("sq i" in t and "theorem sum_sq_safe1" in t and "0 ≤ i ∧ i ≥ 0" in t, t)
    names_ = [i.name for i in translate_source("def g(x: int) -> int:\n  return f(x)\ndef f(x: int) -> int:\n  return x\n")]
    ok(names_ == ["f", "g"], names_)
    rec = translate_source("def f(n: Nat) -> int:\n  return f(n)\n")[0]
    ok("recursive" in rec.refused, rec.refused)
    ok("calls 'f', which was refused" in translate_source(
        "def f(n: Nat) -> int:\n  return f(n)\ndef g(n: Nat) -> int:\n  return f(n)\n")[1].refused)
    # lists
    t = lean_text(by["fill_squares"])
    ok("Array.replicate (n).toNat 0" in t and "a.set! (i).toNat (i * i)" in t and "theorem fill_squares_safe1" in t, t)
    t = lean_text(by["swap_first"])
    ok("(b.set! (0).toNat b[(1).toNat]!).set! (1).toNat b[(0).toNat]!" in t, t)
    ok("r := (r.push i)" in lean_text(by["evens"]) and "let mut r : Array ℤ := #[]" in lean_text(by["evens"]))
    ok("let mut x : ℤ := a[(idx_x).toNat]!" in lean_text(by["total"]), lean_text(by["total"]))
    for src, why in (("def f(a: list[int]) -> int:\n  a[0] = 1\n  return 0\n", "changed in place"),
                     ("def f(a: list[int]) -> int:\n  b = a\n  return 0\n", "alias"),
                     ("def f(n: Nat) -> int:\n  for i in range(n):\n    n -= 1\n  return 0\n", "range bound"),
                     ("def f(n: Nat) -> int:\n  for i in range(n):\n    pass\n  return i\n", "used after the loop"),
                     ("def f(a: list[int]) -> int:\n  b = []\n  return 0\n", "needs a declared type"),
                     ("def f(a: list[int]) -> int:\n  return a + a\n", "arithmetic on a list")):
        r = translate_source(src)[0]
        ok(why in r.refused, (src, r.refused))
    # every verification condition of every model holds on a grid of small values, and none is vacuous there
    for f in sorted(models.glob("*.py")):
        src = f.read_text()
        for fn, name, sat, bad in smoke(src, translate_source(src)):
            ok(not bad, (fn, name, bad))
            ok(sat > 0, (fn, name, "hypotheses never satisfied on the grid"))
    print(f"pyast2lean: {n[0]} checks ok")


def main(argv=None):
    import argparse, pathlib
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("files", nargs="*")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--smoke", action="store_true",
                    help="evaluate every verification condition on a grid of small values (a sanity check, not a proof)")
    ap.add_argument("--prove", action="store_true", help="try omega/nlinarith before sorry in the VC theorems")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    here = pathlib.Path(__file__).resolve().parent
    files = a.files or sorted(str(p) for p in (here / "py_models").glob("*.py"))
    lean, tex, failed = [], [], False
    for f in files:
        src = open(f).read()
        infos = translate_source(src, a.prove)
        if a.smoke:
            for fn, name, sat, bad in smoke(src, infos):
                print(f"smoke {fn}.{name}: " + ("FAILS at " + str(bad) if bad else
                      f"{sat} satisfying sample(s)" if sat else "hypotheses never satisfied on the grid (vacuous here)"))
                failed = failed or bool(bad)
        for info in infos:
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
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
