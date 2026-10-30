"""Equation registry for sympy_models.

Each model module under models/ builds sympy expressions, wraps the ones worth
keeping in an `Equation`, and verifies them symbolically in `build_and_verify()`.
`export.py` collects every registered equation and writes LaTeX and JSON, so a
later document can pull all the equations from one place.

An Equation has a sympy `expr` (used for computation and, by default, for the
LaTeX), a hand-written `latex_str` (for theorem statements with limits, o(.) terms
or set-theoretic content that sympy cannot express), or both. When both are
given, `latex_str` is what is exported and `expr` is what the checks use.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import sympy as sp

# Source kinds. Keep the provenance honest: the repo's papers, the fork's own
# markdown notes, and textbook facts are different levels of evidence.
SOURCE_KINDS = ("paper", "repo-doc", "textbook", "derived", "knowledge-base")
BR = "%%BR%%"  # optional line-break marker inside latex_str


@dataclass
class Equation:
    key: str                          # unique id, used as the LaTeX label: eq:<key>
    title: str                        # short human title
    source_kind: str                  # one of SOURCE_KINDS
    source: str                       # citation: paper, theorem/equation number, section
    expr: sp.Basic | None = None      # usually sp.Eq(lhs, rhs)
    latex_str: str | None = None      # overrides sp.latex(expr) in the export
    entry: str = ""                   # CONTENTS.md entry number, e.g. "271"
    model: str = ""                   # filled in by register()
    notes: str = ""                   # caveats: what the equation does not say
    verified: list[str] = field(default_factory=list)  # names of passed checks
    cites: list[str] = field(default_factory=list)     # extra bibliography keys (see knowledgebase.references)
    see: list[tuple] = field(default_factory=list)     # (label or None, text, url or None): related objects

    def __post_init__(self):
        if self.expr is None and self.latex_str is None:
            raise ValueError(f"{self.key}: give expr, latex_str or both")

    def latex(self, multiline: bool = False) -> str:
        """LaTeX of the equation. `latex_str` may contain the marker %%BR%% where a line break is
        allowed: it becomes a real break (inside gathered) when multiline=True, a \\quad otherwise."""
        tex = self.latex_str if self.latex_str is not None else sp.latex(self.expr)
        if multiline and BR in tex:
            # the leading {} stops a first '[' being read as gathered's optional alignment argument
            return "\\begin{gathered}\n    {}" + tex.replace(BR, " \\\\\n    {}") + "\n  \\end{gathered}"
        return tex.replace(BR, r"\quad ")

    def to_dict(self) -> dict:
        return {
            "key": self.key,
            "title": self.title,
            "entry": self.entry,
            "latex": self.latex(),
            "source_kind": self.source_kind,
            "source": self.source,
            "model": self.model,
            "notes": self.notes,
            "verified": self.verified,
            "cites": self.cites,
            "see": [list(t) for t in self.see],
        }


REGISTRY: dict[str, Equation] = {}
_CHECKS: list[tuple[str, str, bool]] = []  # (model, check name, passed)


def register(model: str, eq: Equation) -> Equation:
    if eq.source_kind not in SOURCE_KINDS:
        raise ValueError(f"{eq.key}: source_kind must be one of {SOURCE_KINDS}")
    if eq.key in REGISTRY:
        raise ValueError(f"duplicate equation key: {eq.key}")
    eq.model = model
    REGISTRY[eq.key] = eq
    return eq


def check(model: str, name: str, condition: bool, *, eq: Equation | None = None) -> None:
    """Record a verification. A failed check raises, so export never emits an unverified set."""
    _CHECKS.append((model, name, bool(condition)))
    if not condition:
        raise AssertionError(f"[{model}] check failed: {name}")
    if eq is not None:
        eq.verified.append(name)


def is_zero(expr: sp.Basic) -> bool:
    """Symbolic zero test that simplifies first."""
    return sp.simplify(expr) == 0


def checks() -> list[tuple[str, str, bool]]:
    return list(_CHECKS)


def tex_block(eq: Equation) -> str:
    """One numbered equation with provenance as a LaTeX comment and a label."""
    entry = f" (entry {eq.entry})" if eq.entry else ""
    return (
        f"% {eq.title}{entry} [{eq.source_kind}: {eq.source}]\n"
        f"\\begin{{equation}}\\label{{eq:{eq.key}}}\n  {eq.latex(multiline=True)}\n\\end{{equation}}\n"
    )
