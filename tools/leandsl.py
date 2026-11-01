"""A small Python DSL: write the model as a decorated Python function, get LaTeX and Lean 4.

    from leandsl import lean, Nat, Real, require, binom, fact, E, Rational

    @lean.theorem(source="OAI.Heisenberg.Clocks.clock_tail")
    def clock_tail(N: Nat, r: Nat, B: Real, d: Real):
        require(B >= 0, d >= 0, r > 0, 4*E*B*(N + r)*d <= r)
        return binom(N + r, r) * (2*B*d)**r <= Rational(1, 2)**r

    @lean.define
    def energy(m: Real, v: Real) -> Real:
        return m * v**2 / 2

The function body runs once on sympy symbols (typed by the annotations); it is never executed
numerically. ``require`` collects hypotheses, the returned relation is the conclusion, a returned
expression is a definition body. Use ``lean.raw(...)`` for a Lean term the printer cannot express.

What is and is not claimed: the Lean text is printed straight from the sympy tree. Unlike
sym2lean.py there is no LaTeX read-back, so the status of every item is "generated"; ``--check``
(if a ``lean`` exists) only says whether Lean accepts it. Bodies are ``sorry`` unless ``proof=`` is given.
"""
from __future__ import annotations
import sys, re, inspect
from dataclasses import dataclass, field
import sympy as sp
from sympy.printing.latex import LatexPrinter
from sympy.core.function import AppliedUndef

# ---- types -------------------------------------------------------------------------------------
class _T:
    def __init__(self, name, lean): self.name, self.lean = name, lean
    def __repr__(self): return self.name
Real, Nat, Int = _T("Real", "ℝ"), _T("Nat", "ℕ"), _T("Int", "ℤ")

E, Rational, pi = sp.E, sp.Rational, sp.pi
binom, fact = sp.binomial, sp.factorial

_CTX: list | None = None


def require(*rels):
    """Add hypotheses to the theorem being defined."""
    if _CTX is None:
        raise RuntimeError("require() outside @lean.theorem")
    for r in rels:
        if r in (sp.true, sp.false) or not isinstance(r, sp.core.relational.Relational):
            raise ValueError(f"hypothesis is not an open relation (sympy decided it?): {r}")
        _CTX.append(r)


class Raw(sp.Symbol):
    """A literal Lean term (and LaTeX) the printer cannot derive."""


def raw(lean_text: str, tex: str | None = None, kind: _T = Real):
    s = Raw(lean_text)
    s._tex = tex or lean_text
    s._kind = kind
    return s


# ---- symbol kinds ------------------------------------------------------------------------------
class _TSym(sp.Symbol):
    """Symbol that remembers its Lean type (sympy's Symbol has __slots__, a subclass does not)."""


def _sym(name: str, t: _T) -> sp.Symbol:
    if t is Nat:
        s = _TSym(name, integer=True, nonnegative=True)
    elif t is Int:
        s = _TSym(name, integer=True)
    else:
        s = _TSym(name, real=True)
    s._lean_t = t
    return s


def kind_of(s) -> _T:
    if isinstance(s, Raw):
        return s._kind
    return getattr(s, "_lean_t", Real)


# ---- Lean printer ------------------------------------------------------------------------------
class Unprintable(Exception):
    pass


GREEK = {"alpha": "α", "beta": "β", "gamma": "γ", "delta": "δ", "epsilon": "ε", "theta": "θ",
         "lambda": "lam", "mu": "μ", "sigma": "σ", "phi": "φ", "rho": "ρ", "omega": "ω"}


def lean_ident(name: str) -> str:
    n = GREEK.get(name, name)
    n = re.sub(r"[^A-Za-z0-9_'αβγδεθμσφρω]", "_", n)
    return n if n and not n[0].isdigit() else "x_" + n


def is_natexpr(e) -> bool:
    """True when e can be printed in ℕ with no loss: Nat symbols, integer literals, + * ^ and binom/fact."""
    if isinstance(e, sp.Integer):
        return e >= 0
    if isinstance(e, sp.Symbol):
        return kind_of(e) is Nat
    if isinstance(e, (sp.Add, sp.Mul)):
        return all(is_natexpr(a) for a in e.args)
    if isinstance(e, sp.Pow):
        return is_natexpr(e.base) and isinstance(e.exp, sp.Integer) or (is_natexpr(e.base) and is_natexpr(e.exp))
    if isinstance(e, (sp.binomial, sp.factorial)):
        return all(is_natexpr(a) for a in e.args)
    return False


def _paren(s: str) -> str:
    return s if re.fullmatch(r"[\w.'αβγδεθμσφρω]+|\(.*\)", s) and _balanced(s) else f"({s})"


def _balanced(s: str) -> bool:
    d = 0
    for i, c in enumerate(s):
        d += c == "("
        d -= c == ")"
        if d == 0 and i < len(s) - 1 and not re.fullmatch(r"[\w.'αβγδεθμσφρω]+", s):
            return False
    return True


def lean_expr(e, nat=False) -> str:
    """Print e. ``nat=True`` means the surrounding context is ℕ (e must be natexpr); else ℝ."""
    if nat:
        if not is_natexpr(e):
            raise Unprintable(f"not a natural-number expression: {e}")
        return _nat(e)
    return _real(e)


def _nat(e) -> str:
    if isinstance(e, sp.Integer):
        return str(int(e))
    if isinstance(e, sp.Symbol):
        return lean_ident(e.name)
    if isinstance(e, sp.Add):
        return " + ".join(_paren(_nat(a)) if isinstance(a, sp.Add) else _nat(a) for a in e.args)
    if isinstance(e, sp.Mul):
        return " * ".join(_paren(_nat(a)) if isinstance(a, sp.Add) else _nat(a) for a in e.args)
    if isinstance(e, sp.Pow):
        return f"{_paren(_nat(e.base))} ^ {_paren(_nat(e.exp))}"
    if isinstance(e, sp.binomial):
        return f"Nat.choose {_paren(_nat(e.args[0]))} {_paren(_nat(e.args[1]))}"
    if isinstance(e, sp.factorial):
        return f"Nat.factorial {_paren(_nat(e.args[0]))}"
    raise Unprintable(str(e))


def _real(e) -> str:
    if is_natexpr(e) and not isinstance(e, (sp.Integer, sp.Symbol)):
        return f"(({_nat(e)} : ℕ) : ℝ)"
    if isinstance(e, Raw):
        return e.name
    if isinstance(e, sp.Symbol):
        n = lean_ident(e.name)
        t = kind_of(e)
        return n if t is Real else f"({n} : ℝ)"
    if e is sp.pi:
        return "Real.pi"
    if e is sp.E:
        return "(Real.exp 1)"
    if isinstance(e, sp.Integer):
        return str(int(e)) if e >= 0 else f"(-{-int(e)})"
    if isinstance(e, sp.Rational):
        return f"({e.p} / {e.q} : ℝ)"
    if isinstance(e, sp.Float):
        raise Unprintable("floating-point literal: use Rational")
    if isinstance(e, sp.Add):
        out = ""
        for i, a in enumerate(sp.Add.make_args(e)):
            c, rest = a.as_coeff_Mul()
            if i and c < 0:
                out += " - " + _term(-a)
            else:
                out += (" + " if i else "") + _term(a)
        return out
    if isinstance(e, sp.Mul):
        return _term(e)
    if isinstance(e, sp.Pow):
        b, x = e.args
        if x == sp.Rational(1, 2):
            return f"Real.sqrt {_paren(_real(b))}"
        if x == -1:
            return f"{_atomR(sp.Integer(1))} / {_atomR(b)}"
        if is_natexpr(x):
            return f"{_atomR(b)} ^ {_paren(_nat(x))}"
        if isinstance(x, sp.Integer):               # negative integer power
            return f"{_atomR(b)} ^ (({int(x)} : ℤ))"
        return f"{_atomR(b)} ^ ({_real(x)} : ℝ)"   # Real.rpow: meaning differs for b <= 0
    if isinstance(e, sp.exp):
        return f"Real.exp {_atomR(e.args[0])}"
    if isinstance(e, sp.log):
        return f"Real.log {_atomR(e.args[0])}"
    for cls, nm in ((sp.sin, "sin"), (sp.cos, "cos"), (sp.tan, "tan"), (sp.sinh, "sinh"), (sp.cosh, "cosh"),
                    (sp.atan, "arctan"), (sp.tanh, "tanh")):
        if isinstance(e, cls):
            return f"Real.{nm} {_atomR(e.args[0])}"
    if isinstance(e, sp.Abs):
        return f"|{_real(e.args[0])}|"
    if isinstance(e, sp.binomial):
        return f"(({_nat(e)} : ℕ) : ℝ)"
    if isinstance(e, sp.factorial):
        return f"(({_nat(e)} : ℕ) : ℝ)"
    if isinstance(e, sp.Sum):
        f, (k, lo, hi) = e.function, e.limits[0]
        if lo != 0 or kind_of(k) is not Nat:
            raise Unprintable("only sums  k = 0..n  over a Nat index are supported")
        return f"∑ {lean_ident(k.name)} ∈ Finset.range ({_nat(hi + 1) if is_natexpr(hi + 1) else _nat(hi) + ' + 1'}), {_real(f)}"
    if isinstance(e, AppliedUndef):
        return " ".join([lean_ident(e.func.__name__)] + [_atomR(a) for a in e.args])
    raise Unprintable(f"{type(e).__name__}: {e}")


def _atomR(e) -> str:
    s = _real(e)
    return s if re.fullmatch(r"[\w.'αβγδεθμσφρω]+|\(.*\)", s) and _balanced(s) else f"({s})"


def _term(a) -> str:
    """One term of a sum or a product; a leading unary minus covers the whole product."""
    c, rest = a.as_coeff_Mul()
    if a.is_Mul:
        num, den = sp.fraction(a)
        sign = ""
        if num.could_extract_minus_sign():
            sign, num = "-", -num
        if den != 1:
            return f"{sign}{_prod(num)} / {_atomR(den) if not isinstance(den, sp.Mul) else '(' + _prod(den) + ')'}"
        return sign + _prod(num)
    if a.could_extract_minus_sign():
        return "-" + _atomR(-a)
    return _atomR(a) if isinstance(a, sp.Add) else _real(a)


def _prod(m) -> str:
    if not m.is_Mul:
        return _atomR(m) if isinstance(m, sp.Add) else _real(m)
    return " * ".join(_atomR(a) if isinstance(a, sp.Add) else _real(a) for a in m.args)


REL = {sp.Le: "≤", sp.Lt: "<", sp.Ge: "≥", sp.Gt: ">", sp.Eq: "=", sp.Ne: "≠"}


def lean_rel(r) -> str:
    sym = REL[type(r)]
    l, rr = r.lhs, r.rhs
    nat = is_natexpr(l) and is_natexpr(rr)
    if nat:
        return f"{_nat(l)} {sym} {_nat(rr)}"
    return f"{_real(l)} {sym} {_real(rr)}"


# ---- LaTeX -------------------------------------------------------------------------------------
class TexPrinter(LatexPrinter):
    def _print_Symbol(self, e, style="plain"):
        return e._tex if isinstance(e, Raw) else super()._print_Symbol(e, style)


def tex(e) -> str:
    return TexPrinter({"mul_symbol": "dot" if False else None, "inv_trig_style": "power"}).doprint(e)


# ---- items -------------------------------------------------------------------------------------
@dataclass
class Item:
    kind: str                     # "theorem" | "def"
    name: str
    params: list                  # [(sympy symbol, _T)]
    hyps: list = field(default_factory=list)
    concl: object = None          # relation (theorem) or expression (def)
    ret: _T = Real
    doc: str = ""
    source: str = ""
    proof: str = "sorry"
    func: object = None

    def tex(self) -> str:
        pre = ",\\ ".join(tex(h) for h in self.hyps)
        if self.kind == "def":
            args = ",".join(tex(p) for p, _ in self.params)
            return f"{self.name}({args}) := {tex(self.concl)}"
        return (f"\\text{{if }} {pre}\\text{{, then }} " if pre else "") + tex(self.concl)

    def lean(self) -> str:
        groups: dict[str, list[str]] = {}
        for p, t in self.params:
            groups.setdefault(t.lean, []).append(lean_ident(p.name))
        binders = " ".join(f"({' '.join(n)} : {t})" for t, n in groups.items())
        hs = " ".join(f"(h{i} : {lean_rel(h)})" for i, h in enumerate(self.hyps))
        head = f"/-- {self.doc}".rstrip()
        lines = [head, f"LaTeX: ${self.tex()}$".replace("/-", "/ -").replace("-/", "- /")] \
            if False else None
        comment = ["/--"]
        if self.doc:
            comment.append(self.doc.replace("/-", "/ -").replace("-/", "- /"))
        comment.append("$" + self.tex().replace("/-", "/ -").replace("-/", "- /") + "$")
        if self.source:
            comment.append(f"Modelled on `{self.source}`.")
        comment.append("-/")
        if self.kind == "def":
            body = lean_expr(self.concl) if self.ret is Real else lean_expr(self.concl, nat=True)
            pre = "noncomputable " if self.ret is Real else ""
            return "\n".join(comment) + f"\n{pre}def {lean_ident(self.name)} {binders} : {self.ret.lean} :=\n  {body}\n"
        return "\n".join(comment) + f"\ntheorem {lean_ident(self.name)} {binders} {hs} :\n    {lean_rel(self.concl)} := by\n  {self.proof}\n"


REGISTRY: list[Item] = []


def _annotations(fn):
    sig = inspect.signature(fn, eval_str=True)  # annotations may be strings (from __future__ import annotations)
    out = []
    for n, p in sig.parameters.items():
        t = p.annotation if isinstance(p.annotation, _T) else Real
        out.append((_sym(n, t), t))
    ret = sig.return_annotation if isinstance(sig.return_annotation, _T) else Real
    return out, ret


class _Lean:
    """``lean.theorem`` / ``lean.define`` decorators."""

    raw = staticmethod(raw)

    def theorem(self, fn=None, *, source="", doc=None, proof="sorry", name=None):
        def deco(f):
            global _CTX
            params, _ = _annotations(f)
            _CTX = []
            try:
                concl = f(*[p for p, _ in params])
                hyps = list(_CTX)
            finally:
                _CTX = None
            if not isinstance(concl, sp.core.relational.Relational) or concl in (sp.true, sp.false):
                raise ValueError(f"{f.__name__}: a theorem must return an open relation, got {concl!r}")
            it = Item("theorem", name or f.__name__, params, hyps, concl, doc=doc or (f.__doc__ or "").strip(),
                      source=source, proof=proof, func=f)
            REGISTRY.append(it)
            return it
        return deco(fn) if fn else deco

    def define(self, fn=None, *, doc=None, source="", name=None):
        def deco(f):
            params, ret = _annotations(f)
            nm = name or f.__name__
            body = f(*[p for p, _ in params])
            it = Item("def", nm, params, [], sp.sympify(body), ret, doc=doc or (f.__doc__ or "").strip(),
                      source=source, func=f)
            REGISTRY.append(it)
            fnc = sp.Function(nm)
            it.applied = fnc
            return fnc
        return deco(fn) if fn else deco


lean = _Lean()


def render(items=None, flavor="mathlib") -> str:
    items = REGISTRY if items is None else items
    head = ("import Mathlib\n\nnamespace OAI.Dsl\n\nnoncomputable section\n\n" if flavor == "mathlib" else "")
    body = []
    for it in items:
        try:
            body.append(it.lean())
        except Unprintable as ex:
            body.append(f"-- {it.name}: not printed ({ex})\n")
    tail = "\nend\n\nend OAI.Dsl\n" if flavor == "mathlib" else ""
    return head + "\n".join(body) + tail


def tex_fragment(items=None) -> str:
    out = []
    for it in REGISTRY if items is None else items:
        out.append(f"\\subsection*{{\\texttt{{{it.name.replace('_', chr(92) + '_')}}}}}\n"
                   f"{'Source: ' + chr(92) + 'texttt{' + it.source.replace('_', chr(92) + '_') + '}' + chr(92) + chr(92) if it.source else ''}\n"
                   f"\\[{it.tex()}\\]\n")
    return "\n".join(out)


# ---- selftest ----------------------------------------------------------------------------------
def selftest():
    n = [0]

    def ok(c, m=""):
        n[0] += 1
        assert c, m

    saved = list(REGISTRY)
    REGISTRY.clear()
    try:
        @lean.theorem(source="OAI.Heisenberg.Clocks.clock_tail")
        def clock_tail(N: Nat, r: Nat, B: Real, d: Real):
            require(B >= 0, d >= 0, r > 0, 4*E*B*(N + r)*d <= r)
            return binom(N + r, r) * (2*B*d)**r <= Rational(1, 2)**r

        t = clock_tail.lean()
        ok("(N r : ℕ)" in t and "(B d : ℝ)" in t, t)
        ok("(h2 : r > 0)" in t, t)
        ok("(h3 : 4 * (Real.exp 1) * ((N + r : ℕ) : ℝ) * B * d ≤ (r : ℝ))" in t, t)
        ok("Nat.choose (N + r) r" in t and "≤ (1 / 2 : ℝ) ^ r" in t, t)

        @lean.define
        def energy(m: Real, v: Real) -> Real:
            return m * v**2 / 2

        @lean.theorem
        def energy_nonneg(m: Real, v: Real):
            require(m >= 0)
            return energy(m, v) >= 0
        t2 = energy_nonneg.lean()
        ok("energy m v ≥ 0" in t2, t2)
        ok("def energy (m v : ℝ) : ℝ" in REGISTRY[1].lean(), REGISTRY[1].lean())

        @lean.define
        def gsum(n: Nat) -> Real:
            k = _sym("k", Nat)
            return sp.Sum(1/(k + 1), (k, 0, n))
        ok("Finset.range" in REGISTRY[-1].lean(), REGISTRY[-1].lean())

        try:
            @lean.theorem
            def bad(x: Real):
                return sp.Eq(x, x)
            ok(False, "decided relation accepted")
        except ValueError:
            ok(True)
        try:
            lean_expr(sp.Float(1.5))
            ok(False)
        except Unprintable:
            ok(True)
        ok("$" in clock_tail.lean() and "binom" in clock_tail.tex(), clock_tail.tex())
    finally:
        REGISTRY[:] = saved
    print(f"leandsl: {n[0]} checks ok")


def main(argv=None):
    import argparse, runpy, pathlib
    ap = argparse.ArgumentParser(description="Run DSL model files; write build/lean/Dsl_mathlib.lean and build/tex/dsl.tex")
    ap.add_argument("files", nargs="*")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--check", action="store_true", help="give the Lean file to lean if it is installed")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    here = pathlib.Path(__file__).resolve().parent
    sys.path.insert(0, str(here))
    sys.modules.setdefault("leandsl", sys.modules[__name__])   # models import leandsl; share this registry when run as a script
    for f in a.files or sorted(str(p) for p in (here / "dsl_models").glob("*.py")):
        runpy.run_path(f, run_name="dsl_model")
    root = here.parent / "build"
    (root / "lean").mkdir(parents=True, exist_ok=True)
    (root / "tex").mkdir(parents=True, exist_ok=True)
    (root / "lean" / "Dsl_mathlib.lean").write_text(render())
    (root / "tex" / "dsl.tex").write_text(tex_fragment())
    print(f"{len(REGISTRY)} items -> {root / 'lean' / 'Dsl_mathlib.lean'}")
    if a.check:
        import shutil
        print("lean not installed: nothing checked" if not shutil.which("lean") else "lean found: run it from lean/ with `lake env lean`")


if __name__ == "__main__":
    main()
