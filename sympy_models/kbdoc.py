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
from document import cite, esc

NUMBERS = ["no", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten"]


def words(n: int, singular: str, plural: str | None = None) -> str:
    """'one equation', 'six equations'; numbers above ten stay digits."""
    plural = plural or singular + "s"
    return f"{NUMBERS[n] if n < len(NUMBERS) else n} {singular if n == 1 else plural}"


def sym(name: str) -> str:
    """A letter of the knowledge base as inline math: 'k_B' -> $k_{B}$, 'hbar' -> $\\hbar$."""
    return "$" + sp.latex(sp.Symbol(name)) + "$"


def href(url: str | None, text: str) -> str:
    return rf"\href{{{url}}}{{{text}}}" if url else text


def eq_name(eq: kb.KBEquation) -> str:
    """Name of a library equation, linked to its Wikipedia article."""
    return href(kb.wikipedia_url(eq.slug), esc(eq.name))


def quantity_name(name: str) -> str:
    """Name of a quantity, linked to its row of the glossary."""
    return rf"\hyperref[{kb.quantity_label(name)}]{{{esc(name)}}}"


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
copied from RosettaMath's \texttt{{rosettaphys.py}}~\cite{{rosettamath}} (MIT licence), where an equation is a node, a
quantity is a node, and an equation is tied to the quantity each of its letters denotes. The algebra is
different: RosettaMath reads the \LaTeX{{}} with its own term type, here it is read into sympy, so
rearranging, eliminating and checking are done by sympy and every composition arrives as an
equation of the registry with its checks. The graph has {len(g.nodes)} nodes and
{len(g.edges)} edges: {words(len(g.of_kind('equation')), 'equation')}, {words(len(g.of_kind('quantity')), 'quantity', 'quantities')},
{words(len(g.of_kind('reference')), 'reference')}, {words(len(g.of_kind('entry')), 'repository entry', 'repository entries')} and
{words(len(g.of_kind('wikipedia')), 'Wikipedia article')}. The tables of this section are computed from it when the
document is built. Every equation and quantity of the library has a row in Appendices~\ref{{app:kb-readings}}
and~\ref{{app:kb-quantities}}, where it is a labelled, linked target; names link to Wikipedia, citations to the
bibliography, and entries to their sections (for example Sections~\ref{{sec:entry-216}} and~\ref{{sec:entry-271}}).

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
not rearranged, for the reasons in Table~\ref{{tab:kb-refused}}. All of them are listed, with their readings, in
Appendix~\ref{{app:kb-readings}}.

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


def _formula(eq: kb.KBEquation) -> str:
    """Statement of an equation for the appendix, or a note on why none is shown."""
    rel = eq.relation()
    if rel is not None:
        return f"${sp.latex(rel)}$"
    try:
        num, den, rhs = kb.read_rate(eq.latex)
        return rf"$\frac{{d {num}}}{{d {den}}}={sp.latex(rhs)}$"
    except kb.Unreadable:
        return rf"\textit{{not rearranged: {esc((eq.refusal or '').split(':')[0])}}}"


def readings_appendix() -> str:
    lib = kb.library()
    rows = []
    for eq in lib.equations.values():
        reads = ", ".join(f"{sym(s)}: {quantity_name(q)}" for s, q in sorted(eq.reads.items()) if s != "pi")
        src = kb.reference_for_source(eq.source)
        where = []
        if src:
            where.append(f"from{cite([src])}")
        if eq.entry:
            where.append(rf"entry~\hyperref[sec:entry-{eq.entry}]{{{eq.entry}}}")
        extra = f"\\newline {{\\footnotesize {', '.join(where)}}}" if where else ""
        body = rf"{_formula(eq)}" + (rf" \newline {{\footnotesize {reads}}}" if reads else "") + extra
        rows.append(rf"\refstepcounter{{kbeq}}\label{{{kb.kb_label(eq.name)}}}\thekbeq & {eq_name(eq)} & "
                    rf"{{\raggedright {body}}} \\[2pt]")
    table = _longtable(r"@{}rp{0.22\linewidth}p{0.62\linewidth}@{}", r"No. & Equation & Statement and readings", rows,
                       "The equations of the library, as sympy prints them (sympy chooses the layout, not the source), "
                       "with what each letter denotes. The equation names link to Wikipedia; the quantities to Table~\\ref{app:kb-quantities}.")
    return rf"""
\appendix
\newcounter{{kbeq}}\renewcommand{{\thekbeq}}{{K\arabic{{kbeq}}}}
\newcounter{{kbq}}\renewcommand{{\thekbq}}{{Q\arabic{{kbq}}}}
\section{{The equations of the library and their readings}}\label{{app:kb-readings}}
Every equation of the library, with what each of its letters denotes. An equation without readings, or one the
reader declined, is still listed; it takes part in no join. The entries drawn from a repository cite it, and the
entries drawn from a script of this repository point at the CONTENTS entry of that script.

\begingroup\small
{table}
\endgroup
"""


def quantities_appendix() -> str:
    lib = kb.library()
    users: dict[str, list[str]] = {}
    for eq in lib.equations.values():
        for q in eq.quantities():
            users.setdefault(q, []).append(eq.name)
    rows = []
    for q in sorted(users):
        quantity = lib.quantities[q]
        dim = _dim_tex(quantity.dimension)
        shown = users[q][:8]
        refs = ", ".join(rf"\ref{{{kb.kb_label(n)}}}" for n in shown) + (r", \dots" if len(users[q]) > len(shown) else "")
        rows.append(rf"\refstepcounter{{kbq}}\label{{{kb.quantity_label(q)}}}\thekbq & "
                    rf"{href(kb.wikipedia_url(quantity.slug), esc(q))} & {dim} & {refs} \\")
    table = _longtable(r"@{}rp{0.3\linewidth}p{0.14\linewidth}p{0.36\linewidth}@{}", r"No. & Quantity & Dimension & Denoted in", rows,
                       "The quantities that equations of the library denote, with their dimension in the basis $M,L,T,I,K,N$ "
                       "and the equations (Appendix~\\ref{app:kb-readings}) that denote them.")
    return rf"""
\section{{The physical quantities}}\label{{app:kb-quantities}}
Each quantity links to its Wikipedia article, carries a dimension, and lists the equations that denote it.

\begingroup\small
{table}
\endgroup
"""


def _dim_tex(spec: str) -> str:
    return "$" + sp.latex(kb._dimension_of(spec)) + "$" if spec.strip() != "1" else "$1$"


def bibliography(keys) -> str:
    refs = kb.references()
    items = []
    for key in sorted(set(keys), key=lambda k: (refs[k].kind != "preprint", refs[k].kind != "article", refs[k].title.lower())):
        r = refs[key]
        note = f" {esc(r.note)}" if r.note else ""
        year = "" if str(r.year) in r.note else f" {r.year}."
        doi = f" DOI \\href{{https://doi.org/{r.doi}}}{{{esc(r.doi)}}}." if r.doi else ""
        status = f" \\textit{{{esc(r.status)}}}" if r.status else ""
        items.append(rf"\bibitem{{{key}}} {esc(r.authors)}. \emph{{{esc(r.title)}}}.{year}{note}{doi}{status} \url{{{r.url}}}")
    return ("% Generated by sympy_models/kbdoc.py. Do not edit by hand.\n"
            "{\\small\\raggedright\n\\begin{thebibliography}{99}\n" + "\n".join(items) + "\n\\end{thebibliography}\n}\n")


def cited_keys(eqs) -> set[str]:
    keys = {"rosettamath"}
    for eq in eqs:
        keys.update(kb.citations_for(eq))
        if eq.entry:
            keys.update(r.key for r in kb.manuscripts_for(eq.entry))
    for e in kb.library().equations.values():
        ref = kb.reference_for_source(e.source)
        if ref:
            keys.add(ref)
    return keys


def write(out_dir: Path, eqs=()) -> list[Path]:
    path = out_dir / "knowledge_graph.tex"
    path.write_text(
        "% Generated by sympy_models/kbdoc.py from knowledgebase.py. Do not edit by hand.\n"
        + census_section() + collisions_section() + families_section()
        + readings_appendix() + quantities_appendix(),
        encoding="utf-8")
    bib = out_dir / "bibliography.tex"
    bib.write_text(bibliography(cited_keys(eqs)), encoding="utf-8")
    return [path, bib]
