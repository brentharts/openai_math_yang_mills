"""Generate the LaTeX body of the equations document, grouped by CONTENTS entry.

Writes two fragments under build/ that equations_document.tex (the hand-editable wrapper)
\\input's:

    equations_summary.tex     a table: equations and checks per entry and source kind
    equations_by_entry.tex    one \\section per entry, one \\subsection per equation, with its
                              provenance, caveats and the checks it passed

Entry titles are read from ../CONTENTS.md ("**271. Title.**"). Equations with no entry are
collected in a final "Background" section.
"""
from __future__ import annotations

import re
from collections import Counter, OrderedDict
from pathlib import Path

import eqlib
import knowledgebase as kb

NO_ENTRY = "-"
NO_ENTRY_TITLE = "Background and comparison (no repository entry)"
KIND_LABEL = {
    "paper": "paper",
    "repo-doc": "repo-doc",
    "textbook": "textbook",
    "derived": "derived",
    "knowledge-base": "knowledge-base",
}
KIND_COLUMNS = (("paper", "paper"), ("repo-doc", "repo-doc"), ("textbook", "textbook"),
                ("derived", "derived"), ("knowledge-base", "kb"))

_SPECIAL = {
    "\\": r"\textbackslash{}",
    "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
    "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    # OT1 prints these as inverted marks and a dash if typed directly
    "<": r"\textless{}", ">": r"\textgreater{}", "|": r"\textbar{}",
}
# Non-ASCII characters we know how to typeset. Anything else raises, so it is noticed.
_UNICODE = {
    "\u2013": "--", "\u2014": "---", "\u2212": "-", "\u2018": "`", "\u2019": "'",
    "\u201c": "``", "\u201d": "''", "\u2265": r"$\ge$", "\u2264": r"$\le$", "\u2192": r"$\to$",
    "\u03c0": r"$\pi$", "\u03b6": r"$\zeta$", "\u00d7": r"$\times$", "\u221a": r"$\surd$",
}


def esc(text: str) -> str:
    """Escape plain text for LaTeX (OT1 fonts, UTF-8 input)."""
    out = []
    for ch in text:
        if ch in _SPECIAL:
            out.append(_SPECIAL[ch])
        elif ord(ch) < 128:
            out.append(ch)
        elif ch in _UNICODE:
            out.append(_UNICODE[ch])
        elif ch.isalpha():
            out.append(ch)  # Latin letters are handled by inputenc/T1
        else:
            raise ValueError(f"document.py cannot typeset {ch!r} (U+{ord(ch):04X}); add it to _UNICODE")
    return "".join(out)


def _title_from_markdown(raw: str) -> tuple[str, str]:
    """Convert a CONTENTS.md title to (LaTeX, plain text).

    HTML tags are dropped and $`...`$ becomes inline math in the LaTeX form. The plain form is for
    PDF bookmarks, where math is not allowed (used through \\texorpdfstring)."""
    raw = re.sub(r"</?(i|b|em|sup|sub)>", "", raw)
    pieces = re.split(r"\$`(.+?)`\$", raw)
    tex, plain = [], []
    for i, piece in enumerate(pieces):
        tex.append(f"${piece}$" if i % 2 else esc(piece))
        plain.append(piece.replace("{", "").replace("}", "").replace("\\", ""))
    return "".join(tex).rstrip("."), "".join(plain).rstrip(".")


def entry_titles(contents_path: Path) -> dict[str, tuple[str, str]]:
    titles: dict[str, tuple[str, str]] = {}
    if not contents_path.exists():
        return titles
    pat = re.compile(r"^\*\*(\d{3})\. (.*?)\*\*", re.M)
    for m in pat.finditer(contents_path.read_text(encoding="utf-8")):
        try:
            titles[m.group(1)] = _title_from_markdown(m.group(2))
        except ValueError:
            titles[m.group(1)] = (f"Entry {m.group(1)}", f"Entry {m.group(1)}")
    return titles


def _entry_sort_key(entry: str):
    """CONTENTS entries by number, then companion-note pseudo-entries, then equations with no entry."""
    return (entry == NO_ENTRY, 0 if entry.isdigit() else 1, int(entry) if entry.isdigit() else 0, entry)


def group_by_entry(eqs) -> "OrderedDict[str, list]":
    groups: dict[str, list] = {}
    for e in eqs:
        groups.setdefault(e.entry or NO_ENTRY, []).append(e)
    for items in groups.values():      # transcriptions of the papers first, knowledge-base cross-checks after
        items.sort(key=lambda e: e.source_kind == "knowledge-base")
    return OrderedDict(sorted(groups.items(), key=lambda kv: _entry_sort_key(kv[0])))


def _checks_by_equation(eq) -> list[str]:
    return list(eq.verified)


def render_summary(groups, titles) -> str:
    lines = [
        r"\begin{table}[h]",
        r"\centering\small",
        r"\begin{tabular}{@{}rp{0.36\linewidth}rrrrrrr@{}}",
        r"\toprule",
        r"Entry & Result & Eq. & paper & repo-doc & textbook & derived & kb & Checks \\",
        r"\midrule",
    ]
    totals = Counter()
    for entry, eqs in groups.items():
        kinds = Counter(e.source_kind for e in eqs)
        n_checks = sum(len(e.verified) for e in eqs)
        if entry in kb.PSEUDO_ENTRIES:
            short, long_title = kb.PSEUDO_ENTRIES[entry]
            title, shown = esc(long_title), short
        else:
            title = titles.get(entry, (f"Entry {entry}",) * 2)[0] if entry != NO_ENTRY else NO_ENTRY_TITLE
            shown = entry
        label = r"\textendash" if entry == NO_ENTRY else f"\\hyperref[sec:entry-{entry}]{{{shown}}}"
        lines.append(
            f"{label} & {title} & {len(eqs)} & {kinds['paper']} & {kinds['repo-doc']} & "
            f"{kinds['textbook']} & {kinds['derived']} & {kinds['knowledge-base']} & {n_checks} \\\\"
        )
        totals["eq"] += len(eqs)
        totals["checks"] += n_checks
        for k in KIND_LABEL:
            totals[k] += kinds[k]
    lines += [
        r"\midrule",
        f" & Total & {totals['eq']} & {totals['paper']} & {totals['repo-doc']} & "
        f"{totals['textbook']} & {totals['derived']} & {totals['knowledge-base']} & {totals['checks']} \\\\",
        r"\bottomrule",
        r"\end{tabular}",
        r"\caption{Registered equations and passed checks, by CONTENTS entry and source kind.}",
        r"\end{table}",
    ]
    return "\n".join(lines) + "\n"


def render_body(groups, titles) -> str:
    out = ["% Generated by sympy_models/document.py. Do not edit by hand.\n"]
    for entry, eqs in groups.items():
        if entry == NO_ENTRY:
            heading = esc(NO_ENTRY_TITLE)
            sec_label = "sec:entry-none"
        elif entry in kb.PSEUDO_ENTRIES:
            heading = esc(kb.PSEUDO_ENTRIES[entry][1])
            sec_label = f"sec:entry-{entry}"
        else:
            tex_t, plain_t = titles.get(entry, ("", ""))
            if not tex_t:
                heading = f"Entry {entry}"
            elif tex_t == plain_t:
                heading = f"Entry {entry}: {tex_t}"
            else:
                heading = f"Entry {entry}: \\texorpdfstring{{{tex_t}}}{{{plain_t}}}"
            sec_label = f"sec:entry-{entry}"
        out.append(f"\\section{{{heading}}}\\label{{{sec_label}}}\n")
        papers = kb.manuscripts_for(entry) if entry != NO_ENTRY else []
        if papers:
            items = "; ".join(f"[{esc(r.tag)}]~\\emph{{{esc(r.title)}}}{cite([r.key])}" for r in papers)
            out.append(f"\\noindent\\small\\textit{{Manuscripts.}} {items}.\\normalsize\n")
        for eq in eqs:
            out.append(render_equation(eq))
    return "\n".join(out)


def render_see(items) -> str:
    """Related objects as links: 'text (ref)', with the text a Wikipedia link when there is a URL."""
    out = []
    for label, text, url in items:
        shown = f"\\href{{{url}}}{{{esc(text)}}}" if url else esc(text)
        if label:
            ref = f"\\eqref{{{label}}}" if label.startswith("eq:") else f"\\ref{{{label}}}"
            shown += f"~({ref})"
        out.append(shown)
    return ", ".join(out)


def cite(keys) -> str:
    return f"~\\cite{{{','.join(keys)}}}" if keys else ""


def render_equation(eq) -> str:
    parts = [
        f"% ---- {eq.key} ({eq.model}) ----",
        f"\\subsection{{{esc(eq.title)}}}",
        (f"\\noindent{{\\small\\textsc{{{esc(KIND_LABEL[eq.source_kind])}}}"
         f"\\quad {esc(eq.source)}{cite(kb.citations_for(eq))}\\quad[\\texttt{{{esc(eq.model)}}}]}}"),
        f"\\begin{{equation}}\\label{{eq:{eq.key}}}\n  {eq.latex(multiline=True)}\n\\end{{equation}}",
    ]
    if eq.notes.strip():
        parts.append(f"\\begin{{quote}}\\small\\textit{{Note.}} {esc(eq.notes.strip())}\\end{{quote}}")
    if eq.see:
        parts.append(f"\\begin{{quote}}\\small\\textit{{See.}} {render_see(eq.see)}.\\end{{quote}}")
    checks = _checks_by_equation(eq)
    if checks:
        items = "\n".join(f"  \\item {esc(c)}" for c in checks)
        parts.append(
            "\\ifchecks\n"
            f"\\begin{{quote}}\\footnotesize\\textbf{{Checks passed ({len(checks)}).}}\n"
            "\\begin{itemize}\\itemsep0pt\\parsep0pt\\topsep2pt\\partopsep0pt\n"
            f"{items}\n"
            "\\end{itemize}\\end{quote}\n"
            "\\fi"
        )
    return "\n".join(parts) + "\n"


def write(eqs, out_dir: Path, contents_path: Path) -> list[Path]:
    titles = entry_titles(contents_path)
    groups = group_by_entry(eqs)
    summary = out_dir / "equations_summary.tex"
    body = out_dir / "equations_by_entry.tex"
    summary.write_text(render_summary(groups, titles), encoding="utf-8")
    body.write_text(render_body(groups, titles), encoding="utf-8")
    import kbdoc            # imported here: kbdoc imports esc() from this module
    return [summary, body, *kbdoc.write(out_dir, eqs)]
