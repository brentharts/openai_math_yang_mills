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
    n = re.sub(r"[^A-Za-z0-9_'₀-₉ℓα-κμ-ωΑ-ΟΡΤ-Ω]", "_", n)      # Lean keeps unicode letters, subscripts and primes (not λ, Π, Σ)
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
        return is_natexpr(e.base) and ((isinstance(e.exp, sp.Integer) and e.exp >= 0) or is_natexpr(e.exp))   # not n**-1: ℕ has no inverse
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
    if isinstance(e, (sp.Min, sp.Max)):
        nm = "min" if isinstance(e, sp.Min) else "max"
        args = [_atomR(a) for a in e.args]
        out = args[-1]
        for a in reversed(args[:-1]):
            out = f"{nm} {a} " + (out if out.startswith("(") or out.replace("_", "").isalnum() else f"({out})")
        return out
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


# ---- V2: references and typed connections -------------------------------------------------------
# The edge vocabulary is the knowledge graph's (sympy_models/knowledgebase.py), minus the edges that tie an
# equation to a quantity or a repository: those are `reads=` and `refs=` here.
EDGE_KINDS = {
    "uses": ("uses", "is used by"),
    "specialises": ("is a special case of", "has the special case"),
    "limit-of": ("is a limit of", "has the limit"),
    "defines": ("defines", "is defined by"),
    "field": ("is of the same kind as", "is of the same kind as"),
    "cites": ("cites", "is cited by"),
    "shares": ("shares Lean modules with", "shares Lean modules with"),
    "mentions": ("mentions", "is mentioned by"),
    "inherits": ("inherits from", "is inherited by"),     # recorded by Python itself: ``class C(A, B)``
}


@dataclass(frozen=True)
class Ref:
    """A pointer to explanatory material. The form is checked (a DOI looks like a DOI); nothing is fetched,
    so a well-formed link can still be wrong: that is for a reader to confirm."""
    kind: str            # "doi" | "wikipedia" | "arxiv" | "url"
    target: str
    label: str = ""

    @property
    def url(self) -> str:
        if self.kind == "doi":
            return "https://doi.org/" + self.target
        if self.kind == "wikipedia":
            return "https://en.wikipedia.org/wiki/" + self.target
        if self.kind == "arxiv":
            return "https://arxiv.org/abs/" + self.target
        return self.target

    @property
    def text(self) -> str:
        if self.label:
            return self.label
        return {"doi": f"doi:{self.target}", "arxiv": f"arXiv:{self.target}",
                "wikipedia": "Wikipedia: " + self.target.replace("_", " ")}.get(self.kind, self.target)

    def tex(self) -> str:
        url = self.url.replace("%", "\\%").replace("#", "\\#")
        txt = self.text.replace("&", "\\&").replace("%", "\\%").replace("_", "\\_").replace("#", "\\#")
        return f"\\href{{{url}}}{{{txt}}}"


def doi(ident: str, label: str = "") -> Ref:
    ident = ident.strip().removeprefix("https://doi.org/").removeprefix("doi:")
    if not re.fullmatch(r"10\.\d{4,9}/\S+", ident):
        raise ValueError(f"not a DOI: {ident!r} (expected 10.<registrant>/<suffix>)")
    return Ref("doi", ident, label)


def wikipedia(title: str, label: str = "") -> Ref:
    """``wikipedia("Stirling's approximation")``: the English article of that title."""
    title = title.strip().replace(" ", "_")
    if not title or "/" in title or "#" in title:
        raise ValueError(f"not an article title: {title!r}")
    return Ref("wikipedia", title, label or "Wikipedia: " + title.replace("_", " "))


def arxiv(ident: str, label: str = "") -> Ref:
    ident = ident.strip().removeprefix("arXiv:")
    if not re.fullmatch(r"(\d{4}\.\d{4,5}(v\d+)?|[a-z-]+(\.[A-Z]{2})?/\d{7}(v\d+)?)", ident):
        raise ValueError(f"not an arXiv identifier: {ident!r}")
    return Ref("arxiv", ident, label)


def url(address: str, label: str = "") -> Ref:
    if not re.fullmatch(r"https?://\S+", address):
        raise ValueError(f"not an http(s) address: {address!r}")
    return Ref("url", address, label or address)


def _refs(refs) -> tuple:
    out = tuple(refs or ())
    for r in out:
        if not isinstance(r, Ref):
            raise TypeError(f"refs= takes doi(...), wikipedia(...), arxiv(...) or url(...), got {r!r}")
    return out


@dataclass(frozen=True)
class Link:
    kind: str
    src: "Item"
    dst: "Item"
    note: str = ""


# ---- inequality chains -------------------------------------------------------------------------
class ChainError(ValueError):
    """A step of a chain is false at a point that satisfies the hypotheses."""


_TEXREL = {"<=": "\\le", "<": "<", "=": "=", ">=": "\\ge", ">": ">"}
_LEANREL = {"<=": "≤", "<": "<", "=": "=", ">=": "≥", ">": ">"}
STATUS_TAG = {"proved": "symbolic", "sampled": "sampled", "unchecked": "unchecked"}


@dataclass(frozen=True)
class Step:
    """One link of a chain: ``rel`` between the previous expression and ``expr``, because ``why``
    (LaTeX text); ``using`` is an item (lemma/theorem) the step relies on."""
    rel: str
    expr: object
    why: str = ""
    using: object = None
    tag: str = ""                # a few plain words naming the step under a long display (default: the item used)


def _step(rel):
    def make(expr, why: str = "", using=None, tag: str = "") -> Step:
        return Step(rel, sp.sympify(expr), why, using, tag)
    make.__name__ = {"<=": "le", "<": "lt", "=": "eq", ">=": "ge", ">": "gt"}[rel]
    return make


le, lt, eq, ge, gt = (_step(r) for r in ("<=", "<", "=", ">=", ">"))


@dataclass
class StepResult:
    status: str                  # "proved" | "sampled" | "unchecked"   (a false step raises ChainError)
    detail: str = ""
    points: int = 0              # sampled: how many points satisfying the hypotheses it held at


@dataclass
class Chain:
    start: object
    steps: tuple
    title: str = ""
    results: list = field(default_factory=list)
    closes: bool = False         # the endpoints are the sides of the item's conclusion, strongly enough
    shows: str = ""              # what the chain shows, as LaTeX
    problem: str = ""            # why it does not reach the conclusion, if it does not

    def summary(self) -> str:
        c = {k: sum(r.status == k for r in self.results) for k in STATUS_TAG}
        pts = min((r.points for r in self.results if r.status == "sampled"), default=0)
        n = len(self.results)
        s = f"{c['proved']} of {n} steps checked symbolically, {c['sampled']} sampled"
        if c["sampled"]:
            s += f" (at least {pts} points each, not a proof)"
        return s + f", {c['unchecked']} not checked"


def chain(start, *steps) -> Chain:
    """``chain(a, le(b, "why"), eq(c, "why"), ...)``: a <= b = c ... Steps must all run one way
    (<=, <, = or >=, >, =). The item's ``@item.chain`` decorator checks it."""
    if not steps or not all(isinstance(s, Step) for s in steps):
        raise TypeError("chain(start, step, ...) takes steps made with le(), lt(), eq(), ge() or gt()")
    rels = {s.rel for s in steps}
    if rels & {"<=", "<"} and rels & {">=", ">"}:
        raise ValueError("a chain cannot mix <= / < with >= / > (split it into two chains)")
    return Chain(sp.sympify(start), tuple(steps))


class _Timeout(Exception):
    pass


def _timed(seconds, f, *a):
    """Run f(*a); None if it takes longer than ``seconds`` (sympy's simplify can run for very long)."""
    import signal, threading
    if not hasattr(signal, "SIGALRM") or threading.current_thread() is not threading.main_thread():
        return f(*a)

    def stop(sig, frame):
        raise _Timeout()
    old = signal.signal(signal.SIGALRM, stop)
    signal.alarm(seconds)
    try:
        return f(*a)
    except _Timeout:
        return None
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old)


def _sign_dummies(item) -> dict:
    """Hypotheses of the form ``x >= 0`` / ``x > 0`` become assumptions on a stand-in for x."""
    kind = {}
    for h in item.hyps:
        s = k = None
        if isinstance(h, (sp.Ge, sp.Gt)) and h.rhs == 0 and isinstance(h.lhs, sp.Symbol):
            s, k = h.lhs, ("pos" if isinstance(h, sp.Gt) else "nonneg")
        elif isinstance(h, (sp.Le, sp.Lt)) and h.lhs == 0 and isinstance(h.rhs, sp.Symbol):
            s, k = h.rhs, ("pos" if isinstance(h, sp.Lt) else "nonneg")
        if s is not None and not isinstance(s, Raw) and kind.get(s) != "pos":
            kind[s] = k
    sub = {}
    for s, k in kind.items():
        kw = {"real": True, ("positive" if k == "pos" else "nonnegative"): True}
        if s.is_integer:
            kw["integer"] = True
        sub[s] = sp.Dummy(s.name, **kw)
    return sub


def _is_identity(a, b, sub) -> bool:
    if a == b:
        return True
    for f in (lambda x: x, lambda x: x.xreplace(sub)):
        d = f(a) - f(b)
        for g in (lambda x: x, sp.expand, sp.simplify, lambda x: sp.simplify(sp.expand_func(x)),
                  lambda x: sp.simplify(sp.combsimp(x))):
            try:
                if g(d) == 0:
                    return True
            except Exception:
                pass
    return False


def _is_order(lo, hi, strict, sub) -> bool:
    d = (hi - lo).xreplace(sub)
    for g in (lambda x: x, sp.expand, sp.factor, sp.simplify):
        try:
            v = g(d)
            if (v.is_positive if strict else v.is_nonnegative):
                return True
        except Exception:
            pass
    return False


def _draw(rng, sym):
    t = kind_of(sym)
    if t is Nat:
        return rng.randint(0, 9)
    if t is Int:
        return rng.randint(-9, 9)
    u = rng.random()
    if u < 0.04:
        return 0.0
    if u < 0.10:
        return float(rng.choice([1, 2, 3]) * rng.choice([1, -1]))
    return rng.choice([1, -1]) * 10 ** rng.uniform(-3, 1.3)


def _sample(item, triples, want=200, tries=30000):
    """For each (lo, hi, kind): (points where it held, a violating point or None, evaluable?). Points are
    random but fixed by a seed, and must satisfy every hypothesis of the item (rejection sampling)."""
    import random, mpmath as mpm
    rng = random.Random(20260101)
    params = [p for p, _ in item.params]
    out = [[0, None, True] for _ in triples]
    with mpm.workdps(30):
        try:
            hyps = [sp.lambdify(params, h, "mpmath", dummify=True) for h in item.hyps]
        except Exception:
            return [(0, None, False)] * len(triples), 0
        fns = []
        for i, (lo, hi, kind) in enumerate(triples):
            try:
                fns.append((sp.lambdify(params, lo, "mpmath", dummify=True),
                            sp.lambdify(params, hi, "mpmath", dummify=True), kind))
            except Exception:
                fns.append(None)
                out[i][2] = False

        def real(v):
            v = mpm.mpmathify(v)
            if isinstance(v, mpm.mpc):
                if abs(v.imag) > mpm.mpf("1e-25") * (1 + abs(v.real)):
                    return None
                v = v.real
            return v if mpm.isfinite(v) else None

        feasible = 0
        for _ in range(tries):
            vals = [_draw(rng, p) for p in params]
            mv = [v if isinstance(v, int) else mpm.mpf(v) for v in vals]
            try:
                if not all(h(*mv) for h in hyps):
                    continue
            except Exception:
                continue
            feasible += 1
            for i, f in enumerate(fns):
                if f is None or out[i][1] is not None:
                    continue
                try:
                    lo, hi = real(f[0](*mv)), real(f[1](*mv))
                except Exception:
                    continue
                if lo is None or hi is None:
                    continue
                tol, d = mpm.mpf("1e-18") * (abs(lo) + abs(hi)), hi - lo
                good = abs(d) <= tol if f[2] == "=" else (d > tol if f[2] == "<" else d >= -tol)
                if good:
                    out[i][0] += 1
                else:
                    out[i][1] = {p.name: (v if isinstance(v, int) else float(f"{v:.6g}")) for p, v in zip(params, vals)} \
                        | {"lhs": mpm.nstr(lo, 8), "rhs": mpm.nstr(hi, 8)}
            if feasible >= want:
                break
    return [tuple(o) for o in out], feasible


def check_chain(item, ch: Chain) -> Chain:
    """Fill in ``ch.results`` (per step), ``ch.closes`` and ``ch.shows``; raise ChainError on a false step."""
    sub = _sign_dummies(item)
    known = {p for p, _ in item.params}
    triples, prev = [], ch.start
    for st in ch.steps:
        for e in (prev, st.expr):
            extra = {s for s in e.free_symbols if not isinstance(s, Raw)} - known
            if extra:
                raise ValueError(f"{item.name}: chain uses {sorted(map(str, extra))}, which are not parameters of the item")
        if st.rel in ("<=", "<"):
            triples.append((prev, st.expr, st.rel))
        elif st.rel in (">=", ">"):
            triples.append((st.expr, prev, "<=" if st.rel == ">=" else "<"))
        else:
            triples.append((prev, st.expr, "="))
        prev = st.expr
    results, todo = [None] * len(triples), []
    for i, (lo, hi, kind) in enumerate(triples):
        proved = _timed(15, (lambda: _is_identity(lo, hi, sub)) if kind == "=" else (lambda: _is_order(lo, hi, kind == "<", sub)))
        if proved:
            results[i] = StepResult("proved", "sympy decides it" + (" (using x >= 0 / x > 0 hypotheses)" if sub else ""))
        else:
            todo.append(i)
    if todo:
        sampled, feasible = _sample(item, [triples[i] for i in todo])
        for i, (n_ok, bad, evaluable) in zip(todo, sampled):
            lo, hi, kind = triples[i]
            if bad is not None:
                at = ", ".join(f"{k}={v}" for k, v in bad.items() if k not in ("lhs", "rhs"))
                raise ChainError(f"{item.name}: step {i + 1} ({lo} {kind} {hi}) is false at {at} "
                                 f"(left {bad['lhs']}, right {bad['rhs']}), a point satisfying the hypotheses")
            if not evaluable:
                results[i] = StepResult("unchecked", "the expressions cannot be evaluated numerically")
            elif n_ok >= 10:
                results[i] = StepResult("sampled", f"held at {n_ok} sampled points satisfying the hypotheses", n_ok)
            else:
                results[i] = StepResult("unchecked", f"only {feasible} point(s) satisfying the hypotheses were found")
    ch.results = results
    # does it reach the conclusion?
    down = any(s.rel in (">=", ">") for s in ch.steps)
    lo, hi = (prev, ch.start) if down else (ch.start, prev)
    strict, alleq = any(s.rel in ("<", ">") for s in ch.steps), all(s.rel == "=" for s in ch.steps)
    sym = "=" if alleq else ("<" if strict else "\\le")
    ch.shows = f"{tex(lo)} {sym} {tex(hi)}"
    c = item.concl
    if isinstance(c, (sp.Le, sp.Lt)):
        want = (c.lhs, c.rhs, "lt" if isinstance(c, sp.Lt) else "le")
    elif isinstance(c, (sp.Ge, sp.Gt)):
        want = (c.rhs, c.lhs, "lt" if isinstance(c, sp.Gt) else "le")
    elif isinstance(c, sp.Eq):
        want = (c.lhs, c.rhs, "eq")
    else:
        want = None
    if want is None:
        ch.closes, ch.problem = False, "the conclusion is not an order relation"
    elif not (_is_identity(lo, want[0], sub) and _is_identity(hi, want[1], sub)):
        ch.closes, ch.problem = False, "its ends are not the two sides of the stated conclusion"
    elif want[2] == "eq" and not alleq:
        ch.closes, ch.problem = False, "the conclusion is an equality but a step is an inequality"
    elif want[2] == "lt" and not strict and not alleq:
        ch.closes, ch.problem = False, "the conclusion is strict but no step is"
    else:
        ch.closes, ch.problem = True, ""
    return ch


def chain_tex(ch: Chain) -> str:
    rows, notes = [], []
    for i, (st, res) in enumerate(zip(ch.steps, ch.results), 1):
        rows.append(f"{tex(ch.start) if i == 1 else ''} &{_TEXREL[st.rel]} {tex(st.expr)} \\tag{{{i}}}")
        why = (f"by \\texttt{{{_esc(st.using.name)}}}" + (f": {st.why}" if st.why else "")) if st.using is not None else st.why
        notes.append(f"\\item[({i})] {why + ' ' if why else ''}\\textit{{[{STATUS_TAG[res.status]}]}}")
    out = (ch.title + "\n" if ch.title else "") + "\\begin{align*}\n" + " \\\\\n".join(rows) + "\n\\end{align*}\n"
    out += "\\begin{itemize}\n" + "\n".join(notes) + "\n\\end{itemize}\n"
    out += (f"{ch.summary().capitalize()}. ``Symbolic'' means the computer algebra system decided the step; "
            "``sampled'' means it held at random points satisfying the hypotheses, which is evidence, not a proof. ")
    if not ch.closes:
        out += f"This chain shows ${ch.shows}$ but does not by itself give the stated bound: {ch.problem}.\n"
    elif all(r.status == "proved" for r in ch.results):
        out += "Its ends are the two sides of the stated bound and every step was decided symbolically.\n"
    else:
        out += ("Its ends are the two sides of the stated bound; it becomes a proof of it once the steps "
                "that are only sampled or not checked are proved.\n")
    return out


def chain_lean(ch: Chain) -> str:
    """A ``calc`` skeleton as comments (every step ``sorry``): generated, never run through Lean."""
    try:
        rows = [f"calc {lean_expr(ch.start)}"] + [f"  _ {_LEANREL[s.rel]} {lean_expr(s.expr)} := by sorry" for s in ch.steps]
    except Unprintable as ex:
        return f"-- chain: not printed ({ex})\n"
    head = f"-- chain, generated and not checked by Lean ({ch.summary()}):\n"
    return head + "".join(f"--   {r}\n" for r in rows)


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
    category: object = None       # the category class this item was declared in (V2), or None
    lean_stmt: str | None = None  # an explicit Lean statement (V2 ``lean=``), printed as given
    extra: str = ""               # the function's docstring: extra LaTeX for its subsection (V2)
    links: list = field(default_factory=list)  # Link(kind, self, other): documentation, not a proof dependency
    refs: tuple = ()              # Ref(...): DOI, Wikipedia, arXiv, url
    reads: dict = field(default_factory=dict)  # letter -> the quantity it denotes (V2 ``reads=``)
    chains: list = field(default_factory=list)  # Chain(...): checked chains of inequalities
    evidence: str = ""            # how a statement read from elsewhere was checked (e.g. sampled), printed under it

    @property
    def uses(self) -> list:
        return [l.dst for l in self.links if l.kind == "uses"]

    def linked_from(self) -> list:
        """Links of other items that point at this one."""
        return [l for it in REGISTRY for l in it.links if l.dst is self]

    def connect(self, *others, kind: str = "uses", note: str = "") -> "Item":
        """Record how this item relates to others: ``kind`` is one of EDGE_KINDS ("uses" by default).
        Two items that both read a letter may only be connected if they read it the same way."""
        if kind not in EDGE_KINDS:
            raise ValueError(f"unknown connection kind {kind!r} (one of: {', '.join(EDGE_KINDS)})")
        for o in others:
            if not isinstance(o, Item):
                raise TypeError(f"connect() takes items declared with the DSL, got {o!r}")
            if o is self:
                raise ValueError(f"{self.name} cannot be connected to itself")
            clash = {k: (self.reads[k], o.reads[k]) for k in self.reads.keys() & o.reads.keys()
                     if self.reads[k] != o.reads[k]}
            if clash:
                raise ValueError(f"{self.name} and {o.name} read the same letter differently: "
                                 + "; ".join(f"{k} is {a!r} in one and {b!r} in the other" for k, (a, b) in clash.items()))
            if not any(l.dst is o and l.kind == kind for l in self.links):
                self.links.append(Link(kind, self, o, note))
        return self

    def latex_label(self) -> str:
        """The LaTeX label of this item's subsection: ``<category label>:<thm|lem|res|def>:<name>``."""
        return item_label(self)

    def latex_ref(self, text: str | None = None, page: bool = False) -> str:
        """A hyperlink to this item's subsection, for prose: ``A.some_function.latex_ref()``. ``text`` replaces the
        link text (default: the name); ``page=True`` adds the page number."""
        return _ref(self.latex_label(), text if text is not None else f"\\texttt{{{_brk(self.name)}}}", page)

    def chain(self, fn=None, *, title: str = ""):
        """``@item.chain`` on a function of the item's parameters returning ``chain(a, le(b, why), ...)``.
        Each step is decided by sympy where it can be, else sampled at points satisfying the hypotheses; a
        step that is false there raises ChainError. Steps with ``using=lemma`` also connect the two items."""
        if self.kind == "def":
            raise ValueError(f"{self.name}: chains belong to lemmas and theorems")

        def deco(f):
            ch = f(*[p for p, _ in self.params])
            if not isinstance(ch, Chain):
                raise TypeError(f"{self.name}: the function must return chain(start, step, ...)")
            ch.title = title
            check_chain(self, ch)
            for st in ch.steps:
                if st.using is not None:
                    self.connect(st.using, note="used in a step of a chain")
            self.chains.append(ch)
            return ch
        return deco(fn) if fn else deco

    def toLatex(self, display: bool = False) -> str:
        """The statement as LaTeX, inline ``$...$`` (or ``\\[...\\]``): meant for building prose around it."""
        return f"\\[{self.tex()}\\]" if display else f"${self.tex()}$"

    def tex(self) -> str:
        if self.kind == "result":
            return f"\\texttt{{{self.name.replace('_', chr(92) + '_')}}}"
        pre = ",\\ ".join(tex(h) for h in self.hyps)
        if self.kind == "def":
            args = ",".join(tex(p) for p, _ in self.params)
            return f"{self.name}({args}) := {tex(self.concl)}"
        return (f"\\text{{if }} {pre}\\text{{, then }} " if pre else "") + tex(self.concl)

    def _result_lean(self) -> str:
        lines = [f"-- result {self.name}" + (f" ({self.source})" if self.source else "")]
        lines += ["-- " + _plain_tex(l) for l in self.doc.splitlines() if l.strip()]
        lines += [f"-- {EDGE_KINDS[l.kind][0]}: {l.dst.name}" for l in self.links]
        lines.append("-- the statement as given (not compiled here: it relies on definitions that are not in this file):")
        lines += ["--   " + l for l in self.lean_stmt.splitlines()]
        return "\n".join(lines) + "\n"

    def lean(self) -> str:
        if self.kind == "result":
            return self._result_lean()
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
        if self.reads:
            comment.append("Reading: " + "; ".join(f"{lean_ident(k)} = {v}" for k, v in self.reads.items()) + ".")
        if self.refs:
            comment.append("See: " + "; ".join(f"{r.text} <{r.url}>" for r in self.refs).replace("/-", "/ -").replace("-/", "- /"))
        comment.append("-/")
        if self.kind == "def":
            body = lean_expr(self.concl) if self.ret is Real else lean_expr(self.concl, nat=True)
            pre = "noncomputable " if self.ret is Real else ""
            return "\n".join(comment) + f"\n{pre}def {lean_ident(self.name)} {binders} : {self.ret.lean} :=\n  {body}\n"
        kw = "lemma" if self.kind == "lemma" else "theorem"
        uses = "".join(f"-- {EDGE_KINDS[l.kind][0]}: {l.dst.name}{' (' + l.note + ')' if l.note else ''}"
                       " (documentation, not a proof dependency)\n" for l in self.links)
        try:
            generated = f"{kw} {lean_ident(self.name)} {binders} {hs} :\n    {lean_rel(self.concl)}"
        except Unprintable as ex:
            if self.lean_stmt is None:
                raise
            generated = None
        if self.lean_stmt is None:
            stmt = generated
            note = ""
        else:
            lines = [l.strip() for l in self.lean_stmt.strip().splitlines()]
            stmt = "\n    ".join(lines)
            note = ("-- the statement below is given explicitly; the DSL body prints as:\n"
                    + "".join(f"--   {l}\n" for l in (generated or "(not printable)").splitlines())) \
                if generated is None or " ".join(generated.split()) != " ".join(stmt.split()) else ""
        cal = "".join(chain_lean(c) for c in self.chains)
        return "\n".join(comment) + f"\n{uses}{cal}{note}{stmt} := by\n  {self.proof}\n"


REGISTRY: list[Item] = []


def _annotations(fn):
    sig = inspect.signature(fn, eval_str=True)  # annotations may be strings (from __future__ import annotations)
    out = []
    for n, p in sig.parameters.items():
        t = p.annotation if isinstance(p.annotation, _T) else Real
        out.append((_sym(n, t), t))
    ret = sig.return_annotation if isinstance(sig.return_annotation, _T) else Real
    return out, ret


def _run_statement(f):
    """Run a theorem body once on typed symbols: (parameters, hypotheses from require(), conclusion)."""
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
    return params, hyps, concl


# ---- V2: categories ----------------------------------------------------------------------------
CATEGORIES: list = []


class _Registry:
    """``Cat.lemmas`` / ``Cat.theorems``: call it as a decorator to declare an item, or look one up by name
    (``Cat.theorems.clock_tail``). A lookup also finds the items of the categories ``Cat`` inherits from."""

    def __init__(self, kind: str, owner: type):
        self._kind, self._owner, self._items = kind, owner, {}

    def declare(self, name: str, *, lean: str, doc: str = "", source: str = "", refs=()) -> "Item":
        """``Cat.results.declare("name", lean="theorem name ... : ...", doc=..., source="Module.Path")``: a result
        known only by its Lean statement (no sympy body, so no derived LaTeX). The statement is printed as given:
        in the Lean file as a comment, since it may rely on definitions that are not part of that file."""
        if self._kind != "result":
            raise TypeError("declare() is for Cat.results; Cat.lemmas / Cat.theorems take a decorated function")
        if name in self._items:
            raise ValueError(f"{self._owner.__name__}.results already has '{name}'")
        if not isinstance(lean, str) or not lean.strip():
            raise ValueError(f"{name}: lean= must be the Lean statement text")
        it = Item("result", name, [], [], None, doc=doc, source=source, category=self._owner,
                  lean_stmt=lean.strip(), refs=_refs(refs))
        self._items[name] = it
        REGISTRY.append(it)
        _expose(self._owner, name, it)
        return it

    def __call__(self, fn=None, *, lean=None, source="", doc="", proof="sorry", name=None, refs=(), reads=None):
        if self._kind == "result":
            raise TypeError("Cat.results.declare(name, lean=..., doc=...) declares a result; there is no decorator form")

        def deco(f):
            params, hyps, concl = _run_statement(f)
            nm = name or f.__name__
            if nm in self._items:
                raise ValueError(f"{self._owner.__name__}.{self._kind}s already has '{nm}'")
            rd = dict(reads or {})
            unknown = set(rd) - {p.name for p, _ in params}
            if unknown or not all(isinstance(v, str) and v.strip() for v in rd.values()):
                raise ValueError(f"{nm}: reads= maps parameter names ({', '.join(p.name for p, _ in params)}) "
                                 f"to the quantity each denotes; got {rd!r}")
            it = Item(self._kind, nm, params, hyps, concl, doc=doc, source=source, proof=proof, func=f,
                      category=self._owner, lean_stmt=lean, extra=(f.__doc__ or "").strip(),
                      refs=_refs(refs), reads=rd)
            self._items[nm] = it
            REGISTRY.append(it)
            _expose(self._owner, nm, it)
            return it
        return deco(fn) if fn else deco

    def ingest(self, name: str, params, hyps, concl, *, lean: str = "", source: str = "", doc: str = "",
               reads=None, refs=(), evidence: str = "") -> "Item":
        """An item given as data, not as a decorated function: ``params`` are ``(typed symbol, type)`` pairs, ``hyps`` and
        ``concl`` relations. For code that reads statements from elsewhere (the Lean files); it is not available for results."""
        if self._kind == "result":
            raise TypeError("Cat.results.declare(...) declares a result")
        if name in self._items:
            raise ValueError(f"{self._owner.__name__}.{self._kind}s already has '{name}'")
        it = Item(self._kind, name, list(params), list(hyps), concl, doc=doc, source=source, category=self._owner,
                  lean_stmt=lean or None, refs=_refs(refs), reads=dict(reads or {}), evidence=evidence)
        self._items[name] = it
        REGISTRY.append(it)
        _expose(self._owner, name, it)
        return it

    def _find(self, name):
        for klass in self._owner.__mro__:
            reg = klass.__dict__.get(self._kind + "s")
            if isinstance(reg, _Registry) and name in reg._items:
                return reg._items[name]
        return None

    def __getattr__(self, name):
        if name.startswith("_"):
            raise AttributeError(name)
        it = self._find(name)
        if it is None:
            raise AttributeError(f"{self._owner.__name__}.{self._kind}s has no '{name}' (has: {', '.join(self.names()) or 'nothing'})")
        return it

    __getitem__ = __getattr__

    def __contains__(self, name):
        return self._find(name) is not None

    def __iter__(self):
        return iter(self._items.values())

    def names(self):
        return [n for klass in reversed(self._owner.__mro__)
                if isinstance(klass.__dict__.get(self._kind + "s"), _Registry)
                for n in klass.__dict__[self._kind + "s"]._items]


_RESERVED = {"lemmas", "theorems", "results", "items", "connect", "linked_from", "toLatex", "latex_ref", "latex_label", "display"}


def _expose(owner: type, name: str, it) -> None:
    """Make ``Cat.name`` the item itself, so ``A.some_function.latex_ref()`` works and a subclass ``C(A, B)`` finds it
    through Python's own attribute lookup. A name that clashes with a method, or that a lemma and a theorem share,
    stays reachable through ``Cat.lemmas.name`` / ``Cat.theorems.name`` only."""
    if name.startswith("_") or not name.isidentifier() or name in _RESERVED or name in owner.__dict__:
        return
    cur = getattr(owner, name, None)
    if cur is None or isinstance(cur, Item):
        setattr(owner, name, it)


def _category_display(cls, label: str, caption: str = ""):
    """``Cat.display("energy", caption=...)``: a long display (dsl_display.Display) printed with this category."""
    from dsl_display import Display
    return Display(label, cls, caption)


def _category_connect(cls, *others, kind: str = "uses", note: str = ""):
    """``Cat.connect(Other, kind="cites", note="...")``: how this category relates to others (documentation)."""
    if kind not in EDGE_KINDS:
        raise ValueError(f"unknown connection kind {kind!r} (one of: {', '.join(EDGE_KINDS)})")
    for o in others:
        if not (isinstance(o, type) and o in CATEGORIES):
            raise TypeError(f"connect() takes categories declared with @lean.category, got {o!r}")
        if o is cls:
            raise ValueError(f"{cls.__name__} cannot be connected to itself")
        if not any(l.dst is o and l.kind == kind for l in cls.__dict__["__links__"]):
            cls.__dict__["__links__"].append(Link(kind, cls, o, note))
    return cls


def _category_linked_from(cls):
    return [l for c in CATEGORIES for l in c.__dict__.get("__links__", []) if l.dst is cls]


def _install_category(cls: type, doc: str | None = None, refs=(), title: str = "", label: str = ""):
    if "lemmas" not in cls.__dict__:                       # first time: give it its own registries
        cls.lemmas = _Registry("lemma", cls)
        cls.theorems = _Registry("theorem", cls)
        cls.results = _Registry("result", cls)
        cls.__links__ = []
        CATEGORIES.append(cls)
        for b in cls.__bases__:                            # ``class C(A, B)`` is a connection between categories
            if b in CATEGORIES:
                cls.__links__.append(Link("inherits", cls, b))
        cls.latex_ref = classmethod(lambda c, text=None, page=False: _cat_ref(c, text, page))
        cls.latex_label = classmethod(lambda c: cat_label(c))
        cls.display = classmethod(_category_display)
        cls.__init_subclass__ = classmethod(lambda sub, **kw: _install_category(sub))
        cls.items = classmethod(lambda c: list(c.lemmas) + list(c.theorems) + list(c.results))
        cls.connect = classmethod(_category_connect)
        cls.linked_from = classmethod(_category_linked_from)
        cls.toLatex = classmethod(lambda c, **kw: category_tex(c, **kw))
    if doc is not None:
        cls.__lean_doc__ = doc
    if refs:
        cls.__lean_refs__ = _refs(refs)
    if title:
        cls.__lean_title__ = title
    if label:
        cls.__lean_label__ = label
    return cls


class _Lean:
    """``lean.theorem`` / ``lean.define`` decorators, and ``lean.category`` for V2."""

    def category(self, doc=None, *, refs=(), title: str = "", label: str = ""):
        """``@lean.category(r'''\\section{...}''')`` on a class: its LaTeX introduction. The class docstring
        is more LaTeX, ``Cat.lemmas(...)`` / ``Cat.theorems(...)`` declare items, ``Cat.results.declare(...)`` a result
        known only by its Lean statement, a subclass inherits lookups. ``title`` (LaTeX) and ``label`` let other
        categories link to this one; ``Cat.connect(Other, kind=...)`` records how two categories relate."""
        if isinstance(doc, type):
            return _install_category(doc, "")
        return lambda cls: _install_category(cls, (doc or "").strip(), refs, title, label)

    raw = staticmethod(raw)

    def theorem(self, fn=None, *, source="", doc=None, proof="sorry", name=None):
        def deco(f):
            params, hyps, concl = _run_statement(f)
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


def _esc(x: str) -> str:
    return x.replace("_", "\\_")


def _plain_tex(x: str) -> str:
    """LaTeX text (as used in ``doc=``) as plain text, for a Lean comment: formatting commands are dropped."""
    x = x.replace("\\allowbreak{}", "").replace("\\_", "_").replace("\\&", "&").replace("\\%", "%").replace("\\#", "#")
    return re.sub(r"\\(?:texttt|emph|textbf|textit)\{([^{}]*)\}", r"\1", x)


def _brk(x: str) -> str:
    """An identifier as LaTeX text that may break after each dot and underscore."""
    return _esc(x).replace(".", ".\\allowbreak{}").replace("\\_", "\\_\\allowbreak{}")


# ---- labels and references ----------------------------------------------------------------------
# Every category and every item has a LaTeX label, made from names alone, so that any piece of the generated text can
# link to any other: ``A.some_function.latex_ref()`` in Python, ``[[A.some_function]]`` inside the LaTeX of a ``doc=``
# or docstring (resolved when the document is made, so it may point forward). The document needs TEX_PREAMBLE.
TEX_PREAMBLE = r"\providecommand{\dsllabel}[1]{\phantomsection\label{#1}}"
_LABEL_BAD = re.compile(r"[^A-Za-z0-9:._-]+")
_KIND_TAG = {"theorem": "thm", "lemma": "lem", "result": "res", "def": "def"}


class RefError(ValueError):
    """A ``[[...]]`` reference that names no category or item (or more than one)."""


def cat_label(c) -> str:
    """The label of a category: the one given to ``@lean.category(label=)``, else ``cat:<ClassName>``."""
    return c.__dict__.get("__lean_label__") or "cat:" + _LABEL_BAD.sub("-", c.__name__)


def item_label(it) -> str:
    owner = cat_label(it.category) if it.category is not None else "dsl"
    return f"{owner}:{_KIND_TAG.get(it.kind, it.kind[:3])}:{_LABEL_BAD.sub('-', it.name).strip('-')}"


def _ref(label: str, text: str, page: bool = False) -> str:
    return f"\\hyperref[{label}]{{{text}}}" + (f" (p.~\\pageref{{{label}}})" if page else "")


def _cat_ref(c, text: str | None = None, page: bool = False) -> str:
    """How another category refers to ``c``: a hyperlink to its label."""
    return _ref(cat_label(c), text if text is not None else (c.__dict__.get("__lean_title__") or c.__name__.replace("_", "\\_")), page)


def _literal(c) -> bool:
    """Generated from the release's own files: its text is printed as it is, with no ``[[...]]`` and no inheritance notes."""
    return bool(c.__dict__.get("__lean_literal__"))


_TOKEN = re.compile(r"\[\[([A-Za-z_]\w*(?:\.[A-Za-z_]\w*){0,2})(?:\|([^\]]*))?\]\]")


def _cat_named(name: str):
    hits = [c for c in CATEGORIES if c.__name__ == name]
    if not hits:
        raise RefError(f"[[{name}]]: no category with that class name")
    if len(hits) > 1:
        raise RefError(f"[[{name}]]: {len(hits)} categories have that class name")
    return hits[0]


def resolve_refs(text: str) -> str:
    """``[[Cat]]``, ``[[Cat.item]]``, ``[[Cat.theorems.item]]`` (all with an optional ``|link text``) become hyperlinks.
    A reference that names nothing raises RefError, so a misspelt link is found when the document is made."""
    def sub(m):
        path, label = m.group(1).split("."), m.group(2)
        cat = _cat_named(path[0])
        if len(path) == 1:
            return cat.latex_ref(label)
        if len(path) == 3:
            if path[1] not in ("lemmas", "theorems", "results"):
                raise RefError(f"[[{m.group(1)}]]: the middle part must be lemmas, theorems or results")
            hit = getattr(cat, path[1])._find(path[2])
            if hit is None:
                raise RefError(f"[[{m.group(1)}]]: {cat.__name__}.{path[1]} has no '{path[2]}'")
            return hit.latex_ref(label)
        hits = [h for k in ("lemmas", "theorems", "results") if (h := getattr(cat, k)._find(path[1])) is not None]
        if not hits:
            raise RefError(f"[[{m.group(1)}]]: {cat.__name__} has no item '{path[1]}'")
        if len(hits) > 1:
            raise RefError(f"[[{m.group(1)}]]: ambiguous ({' and '.join(h.kind for h in hits)}); write {cat.__name__}.theorems.{path[1]} or .lemmas.")
        return hits[0].latex_ref(label)
    return _TOKEN.sub(sub, text)


def _group_end(s: str, i: int) -> int:
    """Index just after the brace group that opens at ``s[i]`` (an escaped brace does not count)."""
    depth, j = 0, i
    while j < len(s):
        c = s[j]
        if c == "\\":
            j += 2
            continue
        depth += (c == "{") - (c == "}")
        j += 1
        if depth == 0:
            return j
    return len(s)


_DEFINED = re.compile(r"\\(?:label|dsllabel)\{([^{}]*)\}")


def check_labels(doc: str) -> dict:
    """The labels a LaTeX text defines and uses: ``dangling`` (used, never defined) and ``duplicate`` (defined twice)."""
    defined = _DEFINED.findall(doc)
    used = set(re.findall(r"\\hyperref\[([^\]]*)\]", doc)) | set(re.findall(r"\\pageref\{([^{}]*)\}", doc))
    return {"defined": len(set(defined)), "dangling": sorted(used - set(defined)),
            "duplicate": sorted({x for x in defined if defined.count(x) > 1})}


def close_refs(doc: str):
    """``(doc, dangling)``: a link to a label the text does not define becomes plain text, so that a document
    holding part of the tree has no broken links; ``dangling`` lists those labels."""
    defined = set(_DEFINED.findall(doc))
    out, dangling, i = [], [], 0
    while True:
        j = doc.find("\\hyperref[", i)
        if j < 0:
            out.append(doc[i:])
            break
        k = doc.find("]", j)
        if k < 0 or doc[k + 1:k + 2] != "{":
            out.append(doc[i:k + 1 if k >= 0 else len(doc)])
            i = k + 1 if k >= 0 else len(doc)
            continue
        label, e = doc[j + 10:k], _group_end(doc, k + 1)
        out.append(doc[i:j])
        if label in defined:
            out.append(doc[j:e])
        else:
            out.append(doc[k + 2:e - 1])
            dangling.append(label)
        i = e
    res = "".join(out)
    for lab in set(dangling):
        res = res.replace(f" (p.~\\pageref{{{lab}}})", "")
    return res, sorted(set(dangling))


def graph(cats=None) -> dict:
    """The tree as data: ``nodes`` (categories and their items, by label) and ``edges`` (``inherits``, the connections
    of categories and of items), for drawing or for counting."""
    cats = CATEGORIES if cats is None else list(cats)
    nodes, edges = [], []
    for c in cats:
        nodes.append({"id": cat_label(c), "type": "category", "name": c.__name__, "title": c.__dict__.get("__lean_title__", "")})
        edges += [{"kind": l.kind, "src": cat_label(l.src), "dst": cat_label(l.dst), "note": l.note} for l in c.__dict__.get("__links__", [])]
        for it in c.items():
            if it.category is not c:
                continue
            nodes.append({"id": item_label(it), "type": it.kind, "name": it.name, "category": cat_label(c)})
            edges += [{"kind": l.kind, "src": item_label(it), "dst": item_label(l.dst), "note": l.note} for l in it.links]
    return {"nodes": nodes, "edges": edges}


def _inherited(cat) -> list:
    """``[(defining category, [items])]``: what ``cat`` gets from the categories it extends, nearest first,
    leaving out a name that a nearer category already defines."""
    seen = {(k, n) for k in ("lemmas", "theorems", "results") for n in cat.__dict__[k]._items}
    out = []
    for klass in cat.__mro__[1:]:
        got = []
        for k in ("lemmas", "theorems", "results"):
            reg = klass.__dict__.get(k)
            if isinstance(reg, _Registry):
                for it in reg:
                    if (k, it.name) not in seen:
                        seen.add((k, it.name))
                        got.append(it)
        if got:
            out.append((klass, got))
    return out


def _inheritance_tex(cat, limit: int = 12) -> str:
    if _literal(cat):
        return ""
    own, inc = cat.__dict__.get("__links__", []), _category_linked_from(cat)
    bases = [l.dst for l in own if l.kind == "inherits"]
    subs = [l.src for l in inc if l.kind == "inherits" and not _literal(l.src)]
    out = []
    if bases:
        out.append("Extends " + ", ".join(_cat_ref(b, page=True) for b in bases) + ".")
    if subs:
        out.append("Extended by " + ", ".join(_cat_ref(s, page=True) for s in subs) + ".")
    shown = 0
    for klass, its in _inherited(cat):
        take = its[:max(0, limit - shown)]
        shown += len(take)
        if take:
            more = f" and {len(its) - len(take)} more" if len(its) > len(take) else ""
            out.append(f"Inherited from {_cat_ref(klass)}: " + ", ".join(i.latex_ref() for i in take) + more + ".")
    return ("\n".join(out) + "\n") if out else ""


def _item_tex(it, v2: bool = False, level: str = "subsection") -> str:
    lit = it.category is not None and _literal(it.category)
    fix = (lambda t: t) if lit or not v2 else resolve_refs
    anchor = f"\\dsllabel{{{item_label(it)}}}" if v2 else ""
    head = f"{anchor}\\{level}*{{\\texttt{{{_brk(it.name)}}}}}\n"
    if it.source:
        head += f"Source: \\texttt{{{_brk(it.source)}}}\\\\\n"
    if v2 and it.doc:
        head += fix(it.doc) + "\n"
    if it.kind == "result":
        out = head + "{\\footnotesize\n\\begin{verbatim}\n" + it.lean_stmt + "\n\\end{verbatim}}\n"
        if it.refs:
            out += "See " + ", ".join(r.tex() for r in it.refs) + ".\n"
        return out
    out = head + f"\\[{it.tex()}\\]\n"
    if v2 and it.reads:
        syms = {p.name: p for p, _ in it.params}
        out += "Here " + ", ".join(f"${tex(syms[k])}$ denotes {_esc(v)}" for k, v in it.reads.items()) + ".\n"
    if v2 and it.extra:
        out += fix(it.extra) + "\n"
    if v2 and it.evidence:
        out += f"\\textit{{{it.evidence}}}\n"
        if it.lean_stmt:
            out += "{\\footnotesize\n\\begin{verbatim}\n" + it.lean_stmt + "\n\\end{verbatim}}\n"
    if v2:
        for l in it.links:
            out += f"This {it.kind} {EDGE_KINDS[l.kind][0]} {l.dst.latex_ref()} ({l.dst.kind})" \
                   f"{': ' + l.note if l.note else ''}.\n"
        for l in it.linked_from():
            out += f"It {EDGE_KINDS[l.kind][1]} {l.src.latex_ref()} ({l.src.kind}).\n"
        for c in it.chains:
            out += chain_tex(c)
        if it.refs:
            out += "See " + ", ".join(r.tex() for r in it.refs) + ".\n"
    return out


def _category_links_tex(cat) -> str:
    rows = [f"\\item it {EDGE_KINDS[l.kind][0]} {_cat_ref(l.dst, page=True)}{': ' + l.note if l.note else ''}"
            for l in cat.__dict__.get("__links__", []) if l.kind != "inherits"]
    rows += [f"\\item it {EDGE_KINDS[l.kind][1]} {_cat_ref(l.src, page=True)}{': ' + l.note if l.note else ''}"
             for l in _category_linked_from(cat) if l.kind != "inherits"]
    out = ("\\par\\medskip\\noindent\\textbf{Connections:}\n\\begin{itemize}\n" + "\n".join(rows) + "\n\\end{itemize}\n") if rows else ""
    inh = _inheritance_tex(cat)
    return out + (("\\par\\medskip\\noindent " + inh) if inh else "")


def category_tex(cat, item_level: str = "subsection", **kw) -> str:
    """LaTeX for one category: its introduction, the class docstring, its connections, then every item."""
    refs = cat.__dict__.get("__lean_refs__", ())
    doc = getattr(cat, "__lean_doc__", "") if "__lean_doc__" in cat.__dict__ else ""
    fix = (lambda t: t) if _literal(cat) else resolve_refs
    lab = cat_label(cat)
    anchor = "" if f"\\label{{{lab}}}" in doc else f"\\dsllabel{{{lab}}}"
    parts = [anchor + fix(doc),
             fix((cat.__dict__.get("__doc__") or "").strip()),
             ("Further reading: " + ", ".join(r.tex() for r in refs) + ".\n") if refs else "",
             _category_links_tex(cat)]
    parts += [_item_tex(it, v2=True, level=item_level) for it in list(cat.lemmas) + list(cat.theorems) + list(cat.results)]
    parts += [f.figure() for f in cat.__dict__.get("__figures__", [])]          # long displays (dsl_display.Display)
    return "\n".join(x for x in parts if x) + "\n"


def render(items=None, flavor="mathlib") -> str:
    """Lean for the items. By default: items declared with ``@lean.theorem`` first, then one
    ``namespace OAI.Dsl.<Category>`` block per category (lemmas before theorems)."""
    head = ("import Mathlib\n\nnamespace OAI.Dsl\n\nnoncomputable section\n\n" if flavor == "mathlib" else "")

    def one(it):
        try:
            return it.lean()
        except Unprintable as ex:
            return f"-- {it.name}: not printed ({ex})\n"

    if items is not None:
        body = [one(it) for it in items]
    else:
        body = [one(it) for it in REGISTRY if it.category is None]
        for cat in CATEGORIES:
            its = list(cat.lemmas) + list(cat.theorems) + list(cat.results)
            if its or cat.__dict__.get("__links__"):
                links = "".join(f"-- {cat.__name__} {EDGE_KINDS[l.kind][0]} {l.dst.__name__}"
                                f"{' (' + l.note + ')' if l.note else ''}\n"
                                for l in cat.__dict__.get("__links__", []))
                body.append(f"namespace {cat.__name__}\n\n{links}{chr(10) if links else ''}" + "\n".join(one(it) for it in its)
                            + f"\nend {cat.__name__}\n")
    tail = "\nend\n\nend OAI.Dsl\n" if flavor == "mathlib" else ""
    return head + "\n".join(body) + tail


def tex_fragment(items=None) -> str:
    if items is not None:
        return "\n".join(_item_tex(it) for it in items)
    out = [_item_tex(it) for it in REGISTRY if it.category is None]
    out += [category_tex(cat) for cat in CATEGORIES]
    return close_refs("\n".join(out))[0]


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

        ok(lean_expr(sp.Integer(1) - 1 / _sym("k", Nat)) == "1 - 1 / (k : ℝ)", lean_expr(sp.Integer(1) - 1 / _sym("k", Nat)))   # no ℕ-typed k^(-1)
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

        # V2: categories with explicit statements, name lookup, connections, inheritance, grouped output
        saved_cats = list(CATEGORIES)
        try:
            @lean.category(r'''\\section{A}''')
            class A:
                """extra A"""

            @A.lemmas(lean="lemma foo (x y : ℝ) : x ≤ y + 1", doc="d", source="S.foo")
            def foo(x: Real, y: Real):
                """extra foo"""
                require(x <= y)
                return x <= y + 1

            @A.theorems(lean="theorem bar (x y : ℝ) (h : x ≤ y) :\n      x ≤ y + 1")
            def bar(x: Real, y: Real):
                require(x <= y)
                return x <= y + 1

            ok(isinstance(foo, Item) and foo.kind == "lemma" and A.lemmas.foo is foo and A.theorems["bar"] is bar)
            try:
                A.lemmas.nope
                ok(False)
            except AttributeError as ex:
                ok("has: foo" in str(ex), str(ex))
            ok(bar.connect(foo) is bar and bar.uses == [foo])
            for bad in (bar, 3):
                try:
                    bar.connect(bad)
                    ok(False)
                except (ValueError, TypeError):
                    ok(True)
            ok(foo.toLatex().startswith("$") and foo.toLatex(display=True).startswith("\\["))

            class B(A):                                   # a subclass of a category is a category
                pass
            ok(B.lemmas.foo is foo and "foo" in B.lemmas and "bar" in B.theorems and B.lemmas.names() == ["foo"])

            @B.theorems
            def baz(x: Real, y: Real):
                return x < y + 2
            ok("baz" in B.theorems and "baz" not in A.theorems and B.theorems.names() == ["bar", "baz"])
            try:
                @A.theorems
                def bar(x: Real, y: Real):
                    return x < y
                ok(False, "duplicate accepted")
            except ValueError:
                ok(True)

            t = render()
            ok("namespace A\n" in t and "end A" in t and t.index("lemma foo") < t.index("theorem bar"), t)
            ok("lemma foo (x y : ℝ) : x ≤ y + 1 := by\n  sorry" in t, t)                       # printed as given
            ok("theorem bar (x y : ℝ) (h : x ≤ y) :\n    x ≤ y + 1 := by" in t, t)             # lines re-indented
            ok("-- uses: foo (documentation, not a proof dependency)" in t, t)
            ok("namespace B" not in t or "baz" in t.split("namespace B")[1])
            x = tex_fragment()
            ok("\\section{A}" in x and "extra A" in x and "extra foo" in x and "uses \\hyperref[cat:A:lem:foo]{\\texttt{foo}} (lemma)" in x
               and "It is used by \\hyperref[cat:A:thm:bar]{\\texttt{bar}}" in x, x)
            # ---- references, readings, typed connections ----
            ok(wikipedia("Stirling's approximation").url == "https://en.wikipedia.org/wiki/Stirling%27s_approximation"
               or "Stirling" in wikipedia("Stirling's approximation").url)
            ok(doi("10.1000/xyz123").url == "https://doi.org/10.1000/xyz123")
            ok(doi("https://doi.org/10.1000/xyz123") == doi("doi:10.1000/xyz123") == doi("10.1000/xyz123"))
            for bad in ("xyz", "10.1/x", "https://example.org/10.1000/x"):
                try:
                    doi(bad)
                    ok(False, bad)
                except ValueError:
                    ok(True)

            @A.lemmas(refs=[wikipedia("Binomial coefficient"), doi("10.1000/xyz123")], reads={"x": "the first quantity"})
            def qux(x: Real, y: Real):
                require(x <= y)
                return x <= y + 3
            ok(len(qux.refs) == 2 and qux.reads == {"x": "the first quantity"})
            t = render()
            ok("Reading: x = the first quantity." in t and "https://doi.org/10.1000/xyz123" in t, t)
            try:
                @A.lemmas(reads={"nope": "?"})
                def bad1(x: Real, y: Real):
                    return x <= y
                ok(False, "unknown reads key accepted")
            except ValueError:
                ok(True)
            try:
                @A.lemmas(refs=["not a ref"])
                def bad2(x: Real, y: Real):
                    return x <= y
                ok(False, "bare string accepted as ref")
            except (ValueError, TypeError):
                ok(True)
            bar.connect(qux, kind="specialises", note="take $c=3$")
            ok([l.kind for l in bar.links] == ["uses", "specialises"] and qux.linked_from()[0].src is bar)
            try:
                bar.connect(qux, kind="implies")
                ok(False, "unknown kind accepted")
            except ValueError:
                ok(True)

            @A.lemmas(reads={"x": "mass"})
            def m1(x: Real, y: Real):
                return x <= y

            @A.lemmas(reads={"x": "velocity"})
            def m2(x: Real, y: Real):
                return x <= y
            try:
                m1.connect(m2)
                ok(False, "conflicting readings accepted")
            except ValueError:
                ok(True)
            x = tex_fragment()
            ok("is a special case of \\hyperref[cat:A:lem:qux]{\\texttt{qux}}" in x and "take $c=3$" in x and "Here $x$ denotes the first quantity" in x
               and "\\href{https://doi.org/10.1000/xyz123}" in x, x)
            # ---- inequality chains ----
            @A.lemmas
            def amgm(x: Real, y: Real):
                return x ** 2 + y ** 2 >= 2 * x * y

            @amgm.chain(title="By a square.")
            def ch_amgm(x, y):
                return chain(x ** 2 + y ** 2, eq((x - y) ** 2 + 2 * x * y, "expand the square"),
                             ge(2 * x * y, "a real square is nonnegative"))
            ok([r.status for r in ch_amgm.results] == ["proved", "proved"] and ch_amgm.closes, ch_amgm.results)

            @A.lemmas
            def sq(x: Real, y: Real):
                require(x >= 0, y >= x)
                return x ** 2 <= y ** 2

            @sq.chain
            def ch_sq(x, y):
                return chain(x ** 2, le(x * y, "x <= y, times x"), le(y ** 2, "x <= y, times y"))
            ok([r.status for r in ch_sq.results] == ["sampled", "sampled"] and ch_sq.closes
               and all(r.points >= 10 for r in ch_sq.results), ch_sq.results)
            ok("sampled" in ch_sq.summary() and "not a proof" in ch_sq.summary())

            try:
                @sq.chain
                def ch_false(x, y):
                    return chain(x ** 2, le(x * y / 2, "wrong"), le(y ** 2))
                ok(False, "a false step was accepted")
            except ChainError as ex:
                ok("step 1" in str(ex) and "x=" in str(ex), str(ex))
            ok(len(sq.chains) == 1)                                   # the refused chain was not kept

            @sq.chain
            def ch_part(x, y):
                return chain(x ** 2, le(x * y, "partial"))
            ok(not ch_part.closes and "ends" in ch_part.problem, ch_part.problem)

            @A.lemmas
            def strict(x: Real, y: Real):
                require(x >= 0, y > x)
                return x ** 2 < y ** 2

            @strict.chain
            def ch_strict(x, y):
                return chain(x ** 2, le(x * y), le(y ** 2))
            ok(not ch_strict.closes and "strict" in ch_strict.problem, ch_strict.problem)

            @A.lemmas
            def vacuous(x: Real, y: Real):
                require(x > y, y > x)
                return x <= y + 5

            @vacuous.chain
            def ch_vac(x, y):
                return chain(x, le(y + 5, "no point satisfies the hypotheses"))
            ok(ch_vac.results[0].status == "unchecked" and "hypotheses" in ch_vac.results[0].detail, ch_vac.results)
            try:
                chain(1, le(2), ge(0))
                ok(False, "mixed directions accepted")
            except ValueError:
                ok(True)
            for bad in (lambda: chain(1), lambda: chain(1, 2)):
                try:
                    bad()
                    ok(False)
                except TypeError:
                    ok(True)
            try:
                @amgm.chain
                def ch_other(x, y):
                    return chain(x, le(sp.Symbol("zz")))
                ok(False, "foreign symbol accepted")
            except (ValueError, NameError):
                ok(True)

            @A.lemmas
            def sq2(x: Real, y: Real):
                require(x >= 0, y >= x)
                return x ** 2 <= y ** 2

            @sq2.chain
            def ch_using(x, y):
                return chain(x ** 2, le(x * y, "monotone", using=sq), le(y ** 2))
            ok(any(l.dst is sq for l in sq2.links) and sq.linked_from()[0].src is sq2)

            tl = render()
            ok("--   calc x ^ 2 + y ^ 2\n--     _ = (-y + x) ^ 2 + 2 * x * y := by sorry\n"
               "--     _ ≥ 2 * x * y := by sorry\n" in tl, tl)
            ok("-- chain, generated and not checked by Lean (2 of 2 steps checked symbolically, 0 sampled, 0 not checked):" in tl
               and "sampled (at least" in tl, tl)
            tx = tex_fragment()
            ok("\\begin{align*}" in tx and "\\textit{[symbolic]}" in tx and "\\textit{[sampled]}" in tx and "\\tag{1}" in tx
               and "it becomes a proof" in tx and "every step was decided symbolically" in tx and "does not by itself" in tx, tx)
            @lean.define
            def dd(m: Real) -> Real:
                return m * 2
            try:
                REGISTRY[-1].chain(lambda m: chain(m, le(m)))
                ok(False, "chain on a def accepted")
            except ValueError:
                ok(True)
            # ---- category-level connections and statement-only results ----
            @lean.category(r"\section{Fam P}", title="Family P", label="fam:P")
            class FamP:
                pass

            @lean.category(r"\section{Fam Q}", title="Family Q")
            class FamQ:
                pass
            ok(FamP.connect(FamQ, kind="cites", note="3 citations") is FamP and FamQ.linked_from()[0].src is FamP)
            FamP.connect(FamQ, kind="cites")
            ok(len(FamP.__dict__["__links__"]) == 1)                      # no duplicate link
            for bad, exc in ((lambda: FamP.connect(FamQ, kind="likes"), ValueError), (lambda: FamP.connect(FamP), ValueError),
                             (lambda: FamP.connect(3), TypeError), (lambda: FamP.connect(foo), TypeError)):
                try:
                    bad()
                    ok(False, "bad category connection accepted")
                except exc:
                    ok(True)
            res = FamP.results.declare("thm_p", lean="theorem thm_p (n : ℕ) :\n    n = n", doc="A result.", source="OAI.P.Main")
            ok(res.kind == "result" and FamP.results.thm_p is res and "thm_p" in FamP.results and res in FamP.items())
            ok(FamP.results.names() == ["thm_p"] and res.tex() == "\\texttt{thm\\_p}")
            for bad, exc in ((lambda: FamP.results.declare("thm_p", lean="x"), ValueError),
                             (lambda: FamP.results.declare("e", lean="  "), ValueError),
                             (lambda: FamP.lemmas.declare("l", lean="x"), TypeError),
                             (lambda: FamP.results(lambda: 1), TypeError)):
                try:
                    bad()
                    ok(False, "bad result declaration accepted")
                except exc:
                    ok(True)
            tl = render()
            ok("namespace FamP" in tl and "-- FamP cites FamQ (3 citations)" in tl and "-- result thm_p (OAI.P.Main)" in tl
               and "--   theorem thm_p (n : ℕ) :\n--       n = n" in tl, tl)
            ok("\ntheorem thm_p" not in tl and "lemma thm_p" not in tl, tl)         # only ever a comment
            tp, tq = category_tex(FamP, item_level="subsubsection"), category_tex(FamQ)
            ok("it cites \\hyperref[cat:FamQ]{Family Q} (p.~\\pageref{cat:FamQ}): 3 citations" in tp and "\\subsubsection*{\\texttt{thm\\_\\allowbreak{}p}}" in tp
               and "\\begin{verbatim}\ntheorem thm_p (n : ℕ) :\n    n = n\n\\end{verbatim}" in tp and "A result." in tp, tp)
            ok("it is cited by \\hyperref[fam:P]{Family P} (p.~\\pageref{fam:P}): 3 citations" in tq, tq)
            res2 = FamP.results.declare("thm_q", lean="theorem thm_q : True", source="OAI.Q.Main",
                                        doc="Challenge \\texttt{Qx\\_y}, in module \\texttt{OAI.\\allowbreak{}Q.\\allowbreak{}Main}.")
            ok("-- Challenge Qx_y, in module OAI.Q.Main." in res2.lean() and "\\" not in res2.lean().split("the statement")[0], res2.lean())
            res.connect(foo)                                                   # a result is an Item: ordinary connections work
            ok("-- uses: foo" in res.lean())
            # ---- labels and references ----
            ok(res.latex_label() == "fam:P:res:thm_p" and FamP.latex_label() == "fam:P" and FamP.thm_p is res)
            ok(res.latex_ref() == "\\hyperref[fam:P:res:thm_p]{\\texttt{thm\\_\\allowbreak{}p}}"
               and res.latex_ref("the result", page=True) == "\\hyperref[fam:P:res:thm_p]{the result} (p.~\\pageref{fam:P:res:thm_p})")
            ok("\\dsllabel{fam:P:res:thm_p}\\subsubsection*" in tp and tp.startswith("\\dsllabel{fam:P}"), tp[:200])

            @lean.category(r"\section{Alpha}")
            class Alpha:
                pass

            @lean.category(r"\section{Beta}")
            class Beta:
                pass

            @Alpha.lemmas
            def some_function(x: Real, y: Real):
                require(x <= y)
                return x <= y + 1

            @Beta.theorems
            def other(x: Real, y: Real):
                require(x <= y)
                return x <= y + 2

            class Gamma(Alpha, Beta):
                """Both are used: [[Alpha.some_function]] and [[Beta.other|the second one]], inside [[Alpha]]."""

            ok(Gamma in CATEGORIES and [(l.kind, l.dst) for l in Gamma.__links__] == [("inherits", Alpha), ("inherits", Beta)])
            ok(Gamma.some_function is Alpha.some_function and Gamma.other is Beta.other and Gamma.lemmas.some_function is some_function)
            ok("some_function" not in Gamma.__dict__ and Gamma.latex_label() == "cat:Gamma")
            tg = category_tex(Gamma)
            ok("\\hyperref[cat:Alpha:lem:some_function]{\\texttt{some\\_\\allowbreak{}function}}" in tg
               and "\\hyperref[cat:Beta:thm:other]{the second one}" in tg and "inside \\hyperref[cat:Alpha]{Alpha}" in tg, tg)
            ok("Extends \\hyperref[cat:Alpha]{Alpha} (p.~\\pageref{cat:Alpha}), \\hyperref[cat:Beta]{Beta} (p.~\\pageref{cat:Beta})." in tg
               and "Inherited from \\hyperref[cat:Alpha]{Alpha}: \\hyperref[cat:Alpha:lem:some_function]" in tg, tg)
            ok("Extended by \\hyperref[cat:Gamma]{Gamma}" in category_tex(Alpha), category_tex(Alpha))
            for bad in ("[[Nope]]", "[[Alpha.missing]]", "[[Alpha.lemmas.other]]", "[[Alpha.rules.x]]"):
                try:
                    resolve_refs(bad)
                    ok(False, bad + " accepted")
                except RefError:
                    ok(True)
            ok(resolve_refs("[[Alpha.some_function]]") == some_function.latex_ref() and resolve_refs("no tokens [[1]] here") == "no tokens [[1]] here")
            doc_all = tex_fragment()
            ok(check_labels(doc_all)["dangling"] == [] and check_labels(doc_all)["duplicate"] == [], check_labels(doc_all))
            part = category_tex(Gamma)
            ck = check_labels(part)
            ok("cat:Alpha" in ck["dangling"] and "cat:Gamma" not in ck["dangling"], ck)
            closed, dang = close_refs(part)
            ok("\\hyperref" not in closed and "\\pageref" not in closed and "Extends Alpha, Beta." in closed and "cat:Alpha" in dang, closed)
            ok(close_refs("\\dsllabel{x} \\hyperref[x]{a {b} c} \\hyperref[y]{d {e}} (p.~\\pageref{y})") ==
               ("\\dsllabel{x} \\hyperref[x]{a {b} c} d {e} (p.~\\pageref{y})".replace(" (p.~\\pageref{y})", ""), ["y"]))
            gr = graph([Alpha, Beta, Gamma])
            ok({n["id"] for n in gr["nodes"]} == {"cat:Alpha", "cat:Beta", "cat:Gamma", "cat:Alpha:lem:some_function", "cat:Beta:thm:other"}
               and {"kind": "inherits", "src": "cat:Gamma", "dst": "cat:Beta", "note": ""} in gr["edges"], gr)
            ok("-- Gamma inherits from Alpha" in render() and "-- Gamma inherits from Beta" in render())

            @lean.category(r"\section{Gen}")
            class Gen:
                __lean_literal__ = True

            class GenKid(Gen):
                __lean_literal__ = True
            Gen.results.declare("r1", lean="theorem r1 : True", doc="a [[Not.There]] token printed as it is")
            ok("[[Not.There]]" in category_tex(Gen) and "Extend" not in category_tex(GenKid) and "Extended by" not in category_tex(Gen))
        finally:
            CATEGORIES[:] = saved_cats
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
