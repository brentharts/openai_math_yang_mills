"""Shared helpers for the physics model scripts.

Every script in this folder writes into ``build/`` at the repository root:

    build/figures/*.pdf, *.png   matplotlib figures
    build/tex/*.tex              LaTeX fragments (sympy output, tables, prose)
    build/data/*.json            raw numerical results

``tools/templates/report.tex`` then ``\\input``s the fragments.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # headless; no display needed
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
BUILD = ROOT / "build"
FIG = BUILD / "figures"
TEX = BUILD / "tex"
DATA = BUILD / "data"


def ensure_dirs() -> None:
    for d in (FIG, TEX, DATA):
        d.mkdir(parents=True, exist_ok=True)


def setup_style() -> None:
    plt.rcParams.update(
        {
            "figure.figsize": (6.4, 4.2),
            "axes.grid": True,
            "grid.alpha": 0.3,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "font.size": 10,
            "mathtext.fontset": "cm",
        }
    )


def savefig(fig, name: str) -> None:
    ensure_dirs()
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def _breakable_texttt(text: str) -> str:
    """Allow line breaks after '-' and '/' inside \\texttt{...} (repo folder names are long)."""

    def fix(m: re.Match) -> str:
        body = m.group(1).replace("-", "-\\allowbreak ").replace("/", "/\\allowbreak ")
        return "\\texttt{" + body + "}"

    return re.sub(r"\\texttt\{([^{}]*)\}", fix, text)


def write_tex(name: str, text: str) -> Path:
    ensure_dirs()
    path = TEX / f"{name}.tex"
    path.write_text(_breakable_texttt(text))
    return path


def save_json(name: str, payload: dict) -> Path:
    ensure_dirs()
    path = DATA / f"{name}.json"
    path.write_text(json.dumps(payload, indent=2, default=float))
    return path


def load_json(name: str) -> dict | None:
    path = DATA / f"{name}.json"
    return json.loads(path.read_text()) if path.exists() else None


def equation(latex: str, label: str | None = None) -> str:
    lab = f"\\label{{{label}}}" if label else ""
    return f"\\begin{{equation}}{lab}\n{latex}\n\\end{{equation}}\n"


def table(header: list[str], rows: list[list[str]], colspec: str | None = None,
          caption: str | None = None) -> str:
    colspec = colspec or ("r" * len(header))
    lines = [
        "\\begin{table}[H]\\centering",
        f"\\begin{{tabular}}{{{colspec}}}",
        "\\toprule",
        " & ".join(header) + " \\\\",
        "\\midrule",
    ]
    lines += [" & ".join(r) + " \\\\" for r in rows]
    lines += ["\\bottomrule", "\\end{tabular}"]
    if caption:
        lines.append(f"\\caption{{{caption}}}")
    lines.append("\\end{table}")
    return "\n".join(lines) + "\n"


def figure(name: str, caption: str, width: str = "0.85\\linewidth") -> str:
    return (
        "\\begin{figure}[H]\\centering\n"
        f"\\includegraphics[width={width}]{{{name}.pdf}}\n"
        f"\\caption{{{caption}}}\n"
        "\\end{figure}\n"
    )


def check(condition: bool, message: str) -> None:
    """Fail the build loudly if a numerical sanity check is violated."""
    if not condition:
        print(f"CHECK FAILED: {message}", file=sys.stderr)
        raise SystemExit(1)
    print(f"check ok: {message}")
