"""Exploration scripts: each tests one connection the knowledge graph suggests.

Every module defines run() -> Result. A Result carries a title, a few sentences of prose (LaTeX),
an optional display equation (LaTeX) and the checks that passed. A check that fails raises, so the
concise document is not built from an unverified exploration. Run one on its own:

    python3 sympy_models/explore/su_so.py

or all of them:

    python3 sympy_models/explore/__init__.py     (or: python3 -m explore from sympy_models/)
"""
from __future__ import annotations

import importlib
import sys
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))

ORDER = ("su_so", "spinor_series", "e8_gradings", "generation_bound", "weak_mixing", "hat_golden", "ym_mass_gap", "kb_paths")


@dataclass
class Result:
    name: str
    title: str
    prose: str                       # LaTeX paragraph(s)
    display: str = ""                # LaTeX for an equation environment, or ""
    checks: list[str] = field(default_factory=list)
    data: dict = field(default_factory=dict)

    def check(self, description: str, condition: bool):
        if not condition:
            raise AssertionError(f"{self.name}: {description}")
        self.checks.append(description)


def run_all() -> list[Result]:
    out = []
    for name in ORDER:
        mod = importlib.import_module(f"explore.{name}")
        out.append(mod.run())
    return out


def main(result: Result):
    print(f"{result.title}\n" + "=" * len(result.title))
    for c in result.checks:
        print("  ok ", c)
    for k, v in result.data.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    for r in run_all():
        main(r)
        print()
