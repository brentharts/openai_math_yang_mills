"""LaTeX for the knowledge-graph part of the equations document.

Writes build/knowledge_graph.tex, which equations_document.tex \\input's after the entries. Every
number in it is computed from knowledgebase.py when the file is written; nothing is typed in
(the same rule RosettaMath's rosettapaper.py follows: a paper about a knowledge base is a
second copy of it, and the second copy is the one that goes stale).
"""
from __future__ import annotations

from pathlib import Path

import sympy as sp

import knowledgebase as kb
from document import esc

NUMBERS = ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten"]


def words(n: int, singular: str, plural: str | None = None) -> str:
    """'one equation', 'six equations'; numbers above ten stay digits."""
    plural = plural or singular + "s"
    return f"{NUMBERS[n] if n < len(NUMBERS) else n} {singular if n == 1 else plural}"


def sym(name: str) -> str:
    """A letter of the knowledge base as inline math: 'k_B' -> $k_{B}$, 'hbar' -> $\\hbar$."""
    return "$" + sp.latex(sp.Symbol(name)) + "$"


def _longtable(spec: str, header: str, rows: list[str], caption: str = "") -> str:
    return "\n".join([
        rf"\begin{{longtable}}{{{spec}}}",
        r"\toprule", header + r" \\", r"\midrule", r"\endhead",
        *rows,
        r"\bottomrule",
        rf"\caption{{{caption}}}" if caption else "",
        r"\end{longtable}",
    ])


def census_section() -> str:
    lib = kb.library()
    g = lib.graph
    census = g.census()
    readable = [e for e in lib.equations.values() if e.relation() is not None]
    joins = kb.all_joins()
    rep = kb.dimension_report()
    refused = kb.refusals()
    pivots: dict[str, int] = {}
    for j in joins:
        pivots[j.pivot] = pivots.get(j.pivot, 0) + 1

    edge_rows = "\n".join(
        rf"\texttt{{{esc(k)}}} & {census.get(k, 0)} & {esc(note)} \\"
        for k, note in kb.EDGE_NOTES.items())
    refuse_rows = "\n".join(
        rf"{esc(why)} & {len(names)} & {esc(names[0])} \\" for why, names in refused.items())
    bad_rows = "\n".join(rf"{esc(n)} & {esc(why)} \\" for n, why in rep["bad"])
    pivot_rows = "\n".join(
        rf"{esc(q)} & {n} & {'registered' if q in _registered() else 'counted only'} \\"
        for q, n in sorted(pivots.items(), key=lambda kv: (-kv[1], kv[0])))

    return rf"""
\section{{The knowledge graph}}\label{{sec:kb}}

The equations above are modelled one paper at a time. This part joins them to a library of
{words(len(lib.equations), 'equation')} and {words(len(lib.quantities), 'physical quantity', 'physical quantities')}
copied from RosettaMath's \texttt{{rosettaphys.py}} (MIT licence), where an equation is a node, a
quantity is a node, and an equation is tied to the quantity each of its letters denotes. The algebra is
different: RosettaMath reads the \LaTeX{{}} with its own term type, here it is read into sympy, so
rearranging, eliminating and checking are done by sympy and every composition arrives as an
equation of the registry with its checks. The graph has {len(g.nodes)} nodes and
{len(g.edges)} edges; the tables of this section are computed from it when the document is built.

\begin{{table}}[h]
\centering\small
\begin{{tabular}}{{@{{}}llp{{0.62\linewidth}}@{{}}}}
\toprule
Edge & Count & What it records \\
\midrule
{edge_rows}
\bottomrule
\end{{tabular}}
\caption{{The edge vocabulary. Most edges are derived: \texttt{{denotes}} from the readings tables,
\texttt{{shares}} from \texttt{{denotes}}, \texttt{{cites}} from the \texttt{{source}} of an entry. The
rest are stated by hand in RosettaMath. The edge \texttt{{models}}, which ties a sympy-modelled equation to
the library entries it cross-checks, is added by the model \texttt{{knowledge}}.}}
\end{{table}}

\subsection{{Readings: a glyph is not a quantity}}\label{{sec:kb-readings}}
Each equation declares which quantity each of its letters denotes. Two equations may be joined only on a
quantity both declare, and only if every other letter they share is read the same way; an equation without
readings takes part in no join. {words(len(readable), 'equation')} of {len(lib.equations)} are both readable
and consistent with their readings; the other {len(lib.equations) - len(readable)} stay in the graph and are
not rearranged, for the reasons in Table~\ref{{tab:kb-refused}}.

\begin{{table}}[h]
\centering\footnotesize
\begin{{tabular}}{{@{{}}p{{0.34\linewidth}}rp{{0.44\linewidth}}@{{}}}}
\toprule
Not rearranged because & Count & First example \\
\midrule
{refuse_rows}
\bottomrule
\end{{tabular}}
\caption{{What the reader declines to turn into algebra. A derivative, a sum, a trace or a prime changes what
juxtaposition and equality mean, and treating them as scalar algebra is how a false statement is made to parse.}}
\label{{tab:kb-refused}}
\end{{table}}

\subsection{{Dimensions}}
Each quantity carries a dimension in the basis $M,L,T,I,K,N$, so every readable equation can be tested for
dimensional homogeneity: {words(len(rep['ok']), 'equation')} pass and {words(len(rep['bad']), 'equation')} do not.
{'The failures are listed here as found in the copied library; they are neither repaired nor registered.' if rep['bad'] else ''}

{('\\begin{tabular}{@{}p{0.3\\linewidth}p{0.6\\linewidth}@{}}\\toprule Equation & Reason \\\\\\midrule ' + bad_rows + r' \bottomrule\end{tabular}') if rep['bad'] else ''}

\subsection{{Compositions}}
Eliminating a shared quantity between two equations that both state it gives a new equation. The library yields
{words(len(joins), 'such join')}, and each one has been checked by substituting the pivot back into both
equations (Table~\ref{{tab:kb-pivots}}). Some joins state a relation between the things that both equal the pivot (temperature, entropy, a
transmutation slope); others state only a condition. $mc^2=\tfrac12 mv^2$ is a true consequence of the two
quoted laws on energy and says the rest energy equals the kinetic energy. The model \texttt{{knowledge}} registers
the joins on the pivots chosen for this repository's themes, namely statistical mechanics, lattice scaling and
running couplings, as derived equations in the sections above; the others are counted here and not registered.

\begin{{table}}[h]
\centering\footnotesize
\begin{{tabular}}{{@{{}}p{{0.34\linewidth}}rl@{{}}}}
\toprule
Pivot quantity & Joins & Kind \\
\midrule
{pivot_rows}
\bottomrule
\end{{tabular}}
\caption{{Joins by the quantity eliminated.}}
\label{{tab:kb-pivots}}
\end{{table}}
"""


def _registered() -> set[str]:
    import importlib
    return importlib.import_module("models.knowledge").REGISTERED_PIVOTS


def collisions_section() -> str:
    cols = kb.collisions()
    rows = [rf"{sym(letter)} & {esc(left)} & {esc(right)} & {esc(a.name)} \\"
            for a, _b, letter, left, right in cols]
    table = _longtable(r"@{}lp{0.22\linewidth}p{0.22\linewidth}p{0.3\linewidth}@{}",
                       r"Letter & One reading & The other & First seen in", rows,
                       "Letters the library uses for two different quantities, computed from the declared readings.")
    return rf"""
\subsection{{Collisions}}
Because readings are declared, collisions can be computed: two equations that share a letter and disagree about
it mark a place where physics reuses a character. The library has {words(len(cols), 'distinct collision')}.
A join on the letter alone would equate every pair below; the readings table refuses all of them.

\begingroup\footnotesize
{table}
\endgroup
"""


def families_section() -> str:
    fams = kb.families()
    rows = [rf"{esc(f.quantity)} & {len(f.members)} & {esc('; '.join(e.name for e in f.equations))} \\"
            for f in fams[:12]]
    table = _longtable(r"@{}p{0.2\linewidth}rp{0.62\linewidth}@{}", r"Quantity & Equations & Which", rows,
                       "Quantities stated by three or more equations with different right-hand sides.")
    return rf"""
\subsection{{Families}}
Where three or more equations state the same quantity directly they chain into one statement. The library has
{words(len(fams), 'family', 'families')}; the largest are listed, with the equations that state each.

\begingroup\footnotesize
{table}
\endgroup
"""


def readings_appendix() -> str:
    rows = []
    for eq in kb.library().equations.values():
        rel = eq.relation()
        if rel is None:
            continue
        reads = ", ".join(f"{sym(s)}: {esc(q)}" for s, q in sorted(eq.reads.items()))
        formula = sp.latex(rel)
        rows.append(rf"{esc(eq.name)} & \multicolumn{{1}}{{p{{0.64\linewidth}}}}{{\raggedright ${formula}$ \newline {{\footnotesize {reads}}}}} \\[2pt]")
    table = _longtable(r"@{}p{0.27\linewidth}p{0.66\linewidth}@{}", r"Equation & Statement and readings", rows,
                       "The readable equations of the library, as sympy prints them, with what each letter denotes.")
    return rf"""
\appendix
\section{{The readable equations and their readings}}\label{{app:kb-readings}}
Every equation the reader accepted, as sympy prints it (sympy chooses the layout, not the source), with what
each letter denotes. An equation that is not here takes part in no join.

\begingroup\small
{table}
\endgroup
"""


def write(out_dir: Path) -> Path:
    path = out_dir / "knowledge_graph.tex"
    path.write_text(
        "% Generated by sympy_models/kbdoc.py from knowledgebase.py. Do not edit by hand.\n"
        + census_section() + collisions_section() + families_section() + readings_appendix(),
        encoding="utf-8")
    return path
