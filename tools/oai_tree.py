"""The OpenAI math release as one tree of leandsl categories, generated from the repository's own files.

    root  OpenAICollection
      └ discipline classes        (17, from overview.tex)
          └ family classes        F001 ... F378: summary, manuscripts, Lean scope note, the Comparator
                                  theorems as statement-only results, and edges to other families

Classes are made on demand and cached, so a model file can say

    from oai_tree import family
    @lean.category(...)
    class HeisenbergClockTails(family("271")): ...

and a guide (a fork document such as PHYSICS.md) becomes a class that inherits the families it lists:

    python3 tools/oai_tree.py --guide PHYSICS.md --out DIR [--pdf]

writes DIR/physics.tex (and the PDF, with xelatex) containing the guide, then every family it names, with the
connections that can be read off the release (see oai_catalog.py for what each kind means), then the families
outside the guide that those connections lead to.

The class tree is built from the catalogue only; nothing in it is typed in by hand. What the text says is the
release's own summaries, abstracts and scope notes, converted to LaTeX; the guide's own claims are checked
against the catalogue and every disagreement is listed (``Guide.problems``).
"""
from __future__ import annotations
import re, sys, pathlib, subprocess, collections
from dataclasses import dataclass, field

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import leandsl
from leandsl import lean, CATEGORIES
from oai_catalog import Catalog, ROOT

UPSTREAM = "https://github.com/openai/math/blob/main/"


def _fork_url() -> str:
    try:
        u = subprocess.run(["git", "-C", str(ROOT), "remote", "get-url", "origin"], capture_output=True, text=True,
                           timeout=5).stdout.strip()
        m = re.match(r"(?:https://github\.com/|git@github\.com:)(.+?)(?:\.git)?$", u)
        if m:
            return f"https://github.com/{m.group(1)}/blob/main/"
    except Exception:
        pass
    return UPSTREAM


# ---- markdown (as the release writes it) -> LaTeX -------------------------------------------------
_SPECIAL = {"\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}",
            "~": r"\textasciitilde{}", "^": r"\textasciicircum{}", "$": r"\$"}
_UNI = re.compile("([\u0370-\u03ff\u2070-\u209f\u2100-\u214f\u2190-\u21ff\u2200-\u22ff\u27c0-\u27ef\u2980-\u29ff"
                  "\u2a00-\u2aff\U0001d400-\U0001d7ff]+)")
_LAST_KNOWN: dict = {"families": ()}


def _url_tex(u: str) -> str:
    return u.replace("%", r"\%").replace("#", r"\#")


def resolve_link(url: str, families=None) -> tuple[str, str]:
    """('href'|'ref', target) for a link in the release's markdown: absolute URLs stay, links into the release go
    to the upstream copy, a Lean scope note goes to that family's section when it is in the document."""
    fams = _LAST_KNOWN["families"] if families is None else families
    if re.match(r"https?://", url):
        return "href", url
    path = re.sub(r"^(\.\./)+|^\./", "", url)
    m = re.match(r"(?:lean/)?docs/(\d{3})\.md", path)
    if m and m.group(1) in fams:
        return "ref", f"fam:{m.group(1)}"
    if path.startswith(("preprints/", "lean/ComparatorChallenges/", "ComparatorChallenges/", "lean/docs/", "reasoning_traces/",
                        "overview.pdf", "CONTENTS.md", "figures/")):
        if path.startswith("ComparatorChallenges/"):
            path = "lean/" + path
        return "href", UPSTREAM + path
    return "href", _fork_url() + path        # a page of the fork (YANGMILLS.md, tools/README.md, ...)


def md_to_tex(md: str, families=None) -> str:
    """Inline markdown, with the release's $`math`$ and a little HTML, as LaTeX for xelatex + unicode-math."""
    store: list[str] = []

    def keep(s: str) -> str:
        s = re.sub(r"([_^])\s*(\\math(?:bb|cal|frak|scr|rm|bf|it|sf|tt))\s*([A-Za-z0-9])", r"\1{\2 \3}", s)   # KaTeX: _\mathbb R
        store.append(s)
        return f"\x00{len(store) - 1}\x00"
    md = re.sub(r"\$\$(.+?)\$\$", lambda m: keep("\\[" + m.group(1).strip() + "\\]"), md, flags=re.S)
    md = re.sub(r"\$`(.+?)`\$", lambda m: keep("$" + m.group(1) + "$"), md, flags=re.S)
    md = re.sub(r"(?<!\\)\$([^$\n]+?)(?<!\\)\$", lambda m: keep("$" + m.group(1) + "$"), md)
    md = re.sub(r"`([^`]+)`", lambda m: keep("\\texttt{" + "".join(_SPECIAL.get(c, c) for c in m.group(1)) + "}"), md)

    def plain(t: str) -> str:
        t = "".join(_SPECIAL.get(c, c) for c in t)
        t = t.replace("&emsp;", r"\quad ").replace(r"\&emsp;", r"\quad ")
        for tag, cmd in (("i", r"\emph"), ("em", r"\emph"), ("b", r"\textbf"), ("sub", r"\textsubscript"), ("sup", r"\textsuperscript")):
            t = re.sub(rf"<{tag}>(.*?)</{tag}>", lambda m: cmd + "{" + m.group(1) + "}", t, flags=re.S)
        t = re.sub(r"</?[a-zA-Z][^>]*>", "", t)
        t = re.sub(r"\*\*(\S.*?)\*\*", r"\\textbf{\1}", t)
        t = re.sub(r"(?<![\w*])\*(\S[^*]*?)\*(?![\w*])", r"\\emph{\1}", t)
        return _UNI.sub(lambda m: r"\ensuremath{" + m.group(1) + "}", t)

    out, pos = [], 0
    for m in re.finditer(r"\[((?:[^\[\]]|\[[^\]]*\])*)\]\(([^)\s]+)\)", md):
        out.append(plain(md[pos:m.start()]))
        kind, target = resolve_link(m.group(2), families)
        label = plain(m.group(1))
        out.append(f"\\hyperref[{target}]{{{label}}}" if kind == "ref" else f"\\href{{{_url_tex(target)}}}{{{label}}}")
        pos = m.end()
    out.append(plain(md[pos:]))
    res = "".join(out)
    while "\x00" in res:                                   # placeholders may be nested inside link labels
        res = re.sub(r"\x00(\d+)\x00", lambda m: store[int(m.group(1))], res)
    return res


def plain_title(md: str) -> str:
    """A title for PDF bookmarks: no markup."""
    t = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", md)                 # a markdown link: its text
    t = re.sub(r"\$`?(.+?)`?\$", r"\1", t)
    t = re.sub(r"<[^>]+>", "", t)
    return re.sub(r"[\\{}$`]", "", t)


def blocks_to_tex(md: str, families=None) -> str:
    """Paragraphs, bullet lists and tables of the fork's guide pages as LaTeX (tables become item lists)."""
    out, lines, i = [], md.splitlines(), 0
    while i < len(lines):
        l = lines[i]
        if not l.strip():
            i += 1
        elif l.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not all(re.fullmatch(r":?-{2,}:?", c) for c in cells):
                    rows.append(cells)
                i += 1
            out.append("\\begin{rows}\n" + "\n".join(
                f"\\item \\textbf{{{md_to_tex(r[0], families)}}} \\quad " + " --- ".join(md_to_tex(c, families) for c in r[1:] if c)
                for r in rows[1:]) + "\n\\end{rows}")
        elif re.match(r"\s*[-*] ", l):
            items = []
            while i < len(lines) and re.match(r"\s*[-*] ", lines[i]):
                items.append(re.sub(r"^\s*[-*] ", "", lines[i]))
                i += 1
            out.append("\\begin{itemize}\n" + "\n".join("\\item " + md_to_tex(x, families) for x in items) + "\n\\end{itemize}")
        else:
            para = []
            while i < len(lines) and lines[i].strip() and not lines[i].startswith("|") and not re.match(r"\s*[-*] ", lines[i]):
                para.append(lines[i])
                i += 1
            out.append(md_to_tex(" ".join(para), families) + "\n")
    return "\n".join(out)


# ---- the generated classes ------------------------------------------------------------------------
_CAT: Catalog | None = None
_TREE: dict = {}


def catalog(root=None) -> Catalog:
    global _CAT
    if _CAT is None or (root is not None and pathlib.Path(root) != _CAT.root):
        _CAT = Catalog(root or ROOT)
        _TREE.clear()
    return _CAT


def _ident(s: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "", s.title())


def root_class() -> type:
    if "root" not in _TREE:
        cls = type("OpenAICollection", (), {"__doc__": None, "__oai_generated__": True})
        _TREE["root"] = lean.category("", title="The OpenAI math release", label="oai:root")(cls)
    return _TREE["root"]


def discipline_class(name: str) -> type:
    key = ("d", name)
    if key not in _TREE:
        cls = type("D_" + _ident(name), (root_class(),), {"__oai_generated__": True, "__oai_discipline__": name})
        _TREE[key] = lean.category("", title=md_to_tex(name), label="disc:" + _ident(name))(cls)
    return _TREE[key]


def family_tex(f, include_scope=True) -> str:
    """The body of a family's section: summary, manuscripts with abstracts, the Lean scope note."""
    out = []
    if f.summary:
        out.append(md_to_tex(f.summary))
    if f.manuscripts:
        out.append("\\par\\medskip\\noindent\\textbf{Manuscripts.}\n\\begin{itemize}")
        for m in f.manuscripts:
            out.append(f"\\item \\href{{{_url_tex(UPSTREAM + m.pdf)}}}{{{md_to_tex(m.title)}}}"
                       + (f"\\\\\n{{\\small {md_to_tex(m.abstract)}}}" if m.abstract else ""))
        out.append("\\end{itemize}")
    if include_scope and f.lean_doc:
        out.append("\\par\\medskip\\noindent\\textbf{What the Lean formalization covers.} " + blocks_to_tex(f.lean_doc))
    return "\n".join(out)


def _mod_tex(name: str) -> str:
    """A Lean module path as \\texttt that may break after each dot."""
    return "\\texttt{" + "".join(_SPECIAL.get(c, c) for c in name).replace(".", ".\\allowbreak{}") + "}"


def family(num: str, cat: Catalog | None = None) -> type:
    """The category class of family ``num`` (made on first use): a subclass of its discipline's class."""
    cat = cat or catalog()
    num = f"{int(num):03d}" if str(num).isdigit() else str(num)
    key = ("f", num)
    if key in _TREE:
        return _TREE[key]
    f = cat.families.get(num)
    if f is None:
        raise KeyError(f"no family {num} in the catalogue")
    _LAST_KNOWN["families"] = set(cat.families)
    heading = (f"\\subsection{{\\texorpdfstring{{[{num}] {md_to_tex(f.title)}}}{{[{num}] {plain_title(f.title)}}}}}"
               f"\\label{{fam:{num}}}\n")
    cls = type(f"F{num}", (discipline_class(f.discipline),), {"__oai_generated__": True, "__oai_family__": num})
    _TREE[key] = lean.category(heading + family_tex(f), title=f"[{num}] {md_to_tex(f.title)}", label=f"fam:{num}")(cls)
    for ch in f.challenges:
        for full in ch.theorem_names:
            stmt = ch.statements.get(full)
            if not stmt:
                continue
            same = cls.results._items.get(full)
            if same is not None and same.lean_stmt == stmt:           # the same theorem listed by two challenges
                same.doc += f" Also Comparator challenge \\texttt{{{ch.name}}}."
                continue
            name = full if same is None else f"{full} ({ch.name})"    # a different statement under the same name
            cls.results.declare(name, lean=stmt, source=ch.solution_module or "",
                                doc=f"Comparator challenge \\texttt{{{ch.name}}}, solved in module {_mod_tex(ch.solution_module)}.")
    return _TREE[key]


def connect_families(nums, cat: Catalog | None = None, kinds=("cites", "shares", "mentions")) -> list:
    """Record the catalogue's edges between the given families as category connections; the edges that lead
    outside ``nums`` are returned as (edge, outside family number) so a document can list the neighbours."""
    cat, nums = cat or catalog(), set(nums)
    outside = []
    for e in cat.edges():
        if e.kind not in kinds:
            continue
        a, b = e.src in nums, e.dst in nums
        if a and b:
            note = {"cites": f"{e.weight} bibliography {'entry' if e.weight == 1 else 'entries'} between their manuscripts",
                    "shares": f"{e.weight} shared Lean module{'s' if e.weight != 1 else ''}",
                    "mentions": "the summary says so"}[e.kind]
            family(e.src, cat).connect(family(e.dst, cat), kind=e.kind, note=note)
        elif a or b:
            outside.append((e, e.dst if a else e.src))
    return outside


# ---- guides: the fork's thematic pages ---------------------------------------------------------------
@dataclass
class Row:
    entries: list
    result_md: str
    notes_md: str
    dirs: list


@dataclass
class Section:
    title: str
    intro_md: str
    rows: list


@dataclass
class Guide:
    title: str
    intro_md: str
    sections: list                       # sections with an Entry table
    tail: list                           # (heading, markdown) of the other sections
    problems: list = field(default_factory=list)

    def entries(self) -> list:
        return sorted({e for s in self.sections for r in s.rows for e in r.entries})


def parse_guide(md: str, cat: Catalog) -> Guide:
    """Split a guide page on its ``##`` headings. A section whose table starts ``| Entry |`` is a result section;
    each row's entry numbers are checked against the catalogue (the manuscripts it links must belong to them)."""
    parts = re.split(r"^## ", md, flags=re.M)
    head = parts[0]
    title = re.match(r"# (.*)", head).group(1).strip() if head.startswith("# ") else "Guide"
    g = Guide(title, re.sub(r"^# .*\n", "", head).strip(), [], [])
    for part in parts[1:]:
        heading, _, body = part.partition("\n")
        table = [l for l in body.splitlines() if l.startswith("|")]
        if table and re.match(r"\|\s*Entry\s*\|", table[0]) and re.match(r"\d+\.", heading):
            intro = "\n".join(l for l in body.split("\n\n")[0].splitlines() if not l.startswith("|")).strip()
            sec = Section(re.sub(r"^\d+\.\s*", "", heading).strip(), intro, [])
            for l in table[2:]:
                cells = [c.strip() for c in l.strip().strip("|").split("|")]
                if len(cells) < 3:
                    continue
                row = Row(re.findall(r"\d{3}", cells[0]), cells[1], " | ".join(cells[2:]),
                          list(dict.fromkeys(re.findall(r"preprints/([^/)]+)/", l))))
                if "none" in cells[0].lower() and not row.entries:
                    row.entries = []
                sec.rows.append(row)
                _check_row(g, sec, row, cat)
            g.sections.append(sec)
        else:
            g.tail.append((heading.strip(), body.strip()))
    return g


def _check_row(g: Guide, sec: Section, row: Row, cat: Catalog):
    fams = {cat.manuscripts[d].family for d in row.dirs if d in cat.manuscripts}
    where = f"{sec.title}: {plain_title(row.result_md)[:60]}"
    for e in row.entries:
        if e not in cat.families:
            g.problems.append(f"{where}: entry [{e}] is not a family in the catalogue")
    if row.entries and fams and not fams <= set(row.entries):
        g.problems.append(f"{where}: lists {row.entries} but links manuscripts of {sorted(fams - set(row.entries))}")
    if not row.entries and fams:
        g.problems.append(f"{where}: says no entry, but its manuscript is in family {sorted(fams)}")
    for d in row.dirs:
        if d not in cat.manuscripts:
            g.problems.append(f"{where}: links {d}, "
                              + ("a preprint directory that CONTENTS.md does not list (an older edition?)"
                                 if d in cat.unmapped_dirs else "which is not in preprints/"))


def guide_families(g: Guide, cat: Catalog) -> list:
    nums = set(g.entries())
    for s in g.sections:
        for r in s.rows:
            nums |= {cat.manuscripts[d].family for d in r.dirs if d in cat.manuscripts}
    return sorted(n for n in nums if n in cat.families)


def build_guide(g: Guide, cat: Catalog | None = None, name: str = "Guide", label: str = "guide") -> type:
    """A class inheriting every family the guide names (flat and in ascending order, which always has a consistent
    method resolution order), plus one class per section inheriting that section's families."""
    cat = cat or catalog()
    _LAST_KNOWN["families"] = set(cat.families)
    secs = []
    for k, s in enumerate(g.sections, 1):
        nums = sorted({e for r in s.rows for e in r.entries}
                      | {cat.manuscripts[d].family for r in s.rows for d in r.dirs if d in cat.manuscripts})
        nums = [n for n in nums if n in cat.families]
        body = [f"\\subsection{{{md_to_tex(s.title)}}}\\label{{{label}:{k}}}"]
        if s.intro_md:
            body.append(md_to_tex(s.intro_md))
        body.append("\\begin{rows}")
        for r in s.rows:
            ent = ", ".join(f"\\hyperref[fam:{e}]{{[{e}]}}" if e in cat.families else f"[{e}]" for e in r.entries) or "no entry"
            body.append(f"\\item \\textbf{{{ent}}} \\quad {md_to_tex(r.result_md)} --- {md_to_tex(r.notes_md)}")
        body.append("\\end{rows}")
        bases = tuple(family(n, cat) for n in nums) or (root_class(),)
        cls = type(f"{name}_S{k}", bases, {"__oai_generated__": True})
        secs.append(lean.category("\n".join(body), title=md_to_tex(s.title), label=f"{label}:{k}")(cls))
    allf = tuple(family(n, cat) for n in guide_families(g, cat))
    intro = f"\\section{{{md_to_tex(g.title)}}}\\label{{{label}}}\n" + blocks_to_tex(g.intro_md)
    top = type(name, allf or (root_class(),), {"__oai_generated__": True, "__sections__": secs, "__guide__": g})
    return lean.category(intro, title=md_to_tex(g.title), label=label)(top)


# ---- models that extend a family ----------------------------------------------------------------------
def check_model(cat: Catalog, cls: type) -> list:
    """What a hand-written model claims about the Lean files, checked against them as text: each item's ``source`` must be
    declared in a module behind the family, and a ``uses`` connection is accepted only if the used declaration's name
    occurs in the proof of the other. (Textual checks; Lean has not been run.)"""
    fam = next((c.__oai_family__ for c in cls.__mro__ if c.__dict__.get("__oai_family__")), None)
    if fam is None:
        return [f"{cls.__name__}: does not extend a family class"]
    items = list(cls.lemmas) + list(cls.theorems)
    out, found = [], {}
    for it in items:
        if not it.source:
            continue
        d = cat.declaration(fam, it.source)
        if d is None:
            out.append(f"{cls.__name__}.{it.name}: source {it.source} is not declared in the Lean files behind family {fam}")
        else:
            found[id(it)] = d
    for it in items:
        for l in it.links:
            if l.kind == "uses" and id(it) in found and id(l.dst) in found:
                short = l.dst.source.rsplit(".", 1)[-1]
                proof = found[id(it)][1].split(":=", 1)[-1]
                if not re.search(rf"(?<![\w.'])" + re.escape(short) + r"(?![\w'])", proof):
                    out.append(f"{cls.__name__}.{it.name} is connected as using {l.dst.name}, but {short} does not occur in its proof")
    return out


# ---- the document ---------------------------------------------------------------------------------
PREAMBLE = r"""\documentclass[10pt,a4paper]{article}
\usepackage[a4paper,margin=22mm]{geometry}
\usepackage{fontspec}
\usepackage{amsmath}
\usepackage{unicode-math}
\setmonofont{FreeMono}[Scale=0.95,BoldFont=FreeMonoBold]
\usepackage{xcolor}
\usepackage{fancyvrb}
\newcommand{\fb}[1]{{\fontspec{FreeSerif}#1}}
\usepackage{enumitem}
\setlist{itemsep=2pt,topsep=3pt}
\newlist{rows}{itemize}{1}
\setlist[rows]{label={},leftmargin=0pt,itemsep=4pt}
\usepackage[colorlinks=true,linkcolor=blue!50!black,urlcolor=blue!50!black,citecolor=blue!50!black]{hyperref}
\providecommand{\mathscr}[1]{\mathcal{#1}}
\providecommand{\lt}{<}\providecommand{\gt}{>}
\setlength{\parindent}{0pt}\setlength{\parskip}{4pt}
\setlength{\emergencystretch}{3em}
\setcounter{secnumdepth}{2}
"""

LEGEND = r"""\section*{How to read this document}
Everything here is generated by \texttt{tools/oai\_tree.py} from the files of the repository; nothing is typed in by
hand. The text of each family is the release's own summary, abstracts and Lean scope note. The theorem statements are
the Comparator challenges' statements, exactly as written there, printed as text: they were not compiled in producing
this document. A \emph{connection} between two families is one of

\begin{itemize}
\item \emph{cites}: a manuscript of one cites a manuscript of the other in its bibliography;
\item \emph{shares Lean modules with}: the Lean files behind the two main results import some of the same modules
(modules imported by five or more families are ignored as infrastructure); this is evidence of related machinery,
not of one result depending on the other;
\item \emph{mentions}: the summary names the other result.
\end{itemize}

No connection is inferred from similar subject matter. A family with no Lean scope note has no Lean-based connections.
"""


_NOT_IN_MONO = re.compile("([\U0001d400-\U0001d7ff]+)")       # e.g. the nhds notation 𝓝: FreeMono has no such glyph


def fix_verbatim(tex: str) -> str:
    """Verbatim blocks keep their text exactly; the few characters the monospaced font lacks are set in FreeSerif."""
    def one(m):
        body = m.group(1)
        if not _NOT_IN_MONO.search(body):
            return m.group(0)
        body = _NOT_IN_MONO.sub(lambda k: "\ue000fb\ue001" + k.group(1) + "\ue002", body)
        return "\\begin{Verbatim}[commandchars=\ue000\ue001\ue002]" + body + "\\end{Verbatim}\n"
    return re.sub(r"\\begin\{verbatim\}(.*?)\\end\{verbatim\}", one, tex, flags=re.S)


def render_document(top: type, cat: Catalog, outside=(), title="OpenAI math release: a slice", models=(), problems=()) -> str:
    """The LaTeX document for a guide class: legend, the guide, its families by discipline, the neighbours."""
    g: Guide = top.__guide__
    famcls = [c for c in top.__mro__ if c.__dict__.get("__oai_family__")]
    by_disc = collections.OrderedDict()
    for n in sorted(c.__oai_family__ for c in famcls):
        by_disc.setdefault(cat.families[n].discipline, []).append(n)
    body = [leandsl.category_tex(top)]
    for s in top.__sections__:
        body.append(leandsl.category_tex(s))
    body.append("\\clearpage\\section{The families}\nBy discipline, with their connections and Lean statements.\n")
    order = [d for d in cat.disciplines if d in by_disc]
    for d in order:
        body.append(f"\\subsection*{{\\large {md_to_tex(d)}}}")
        for n in by_disc[d]:
            fc = family(n, cat)
            body.append(leandsl.category_tex(fc, item_level="subsubsection"))
            for m in models:
                if fc in m.__bases__:
                    body.append(f"\\subsubsection*{{Worked model: {m.__dict__.get('__lean_title__') or m.__name__}}}\n"
                                + leandsl.category_tex(m, item_level="paragraph"))
    if outside:
        nb = collections.OrderedDict()
        for e, n in outside:
            nb.setdefault(n, []).append(e)
        body.append("\\clearpage\\section{Neighbouring families outside the guide}\n"
                    "Families the connections above lead to. Only their summaries are given here.\n\\begin{rows}")
        for n, es in sorted(nb.items()):
            f = cat.families[n]
            how = "; ".join(sorted({f"{e.kind} [{e.src if e.dst == n else e.dst}]" + (f" ({e.weight} module{'s' if e.weight != 1 else ''})" if e.kind == 'shares' else "")
                                    for e in es}))
            body.append(f"\\item \\textbf{{[{n}] {md_to_tex(f.title)}}} ({md_to_tex(f.discipline)}). {how}.\\\\\n"
                        f"{{\\small {md_to_tex(f.summary)}}}")
        body.append("\\end{rows}")
    if problems:
        body.append("\\clearpage\\section{Checks of the guide against the catalogue}\n"
                    "Where the guide page and the release's files disagree:\n\\begin{itemize}\n"
                    + "\n".join("\\item " + md_to_tex(p) for p in problems) + "\n\\end{itemize}")
    tail = [f"\\section{{{md_to_tex(h)}}}\n" + blocks_to_tex(b) for h, b in g.tail if h not in ("Related pages",)]
    return fix_verbatim(PREAMBLE + f"\\title{{{md_to_tex(title)}}}\\author{{generated from the repository by \\texttt{{tools/oai\\_tree.py}}}}\n"
            "\\date{}\n\\begin{document}\n\\maketitle\n\\tableofcontents\n" + LEGEND + "\n".join(body) + "\n"
            + "\n".join(tail) + "\n\\end{document}\n")


def build_pdf(tex: pathlib.Path) -> bool:
    ok = True
    for _ in range(2):
        r = subprocess.run(["xelatex", "-interaction=nonstopmode", "-halt-on-error", tex.name], cwd=tex.parent,
                           capture_output=True, text=True)
        ok = ok and r.returncode == 0
        if r.returncode:
            print("\n".join(l for l in r.stdout.splitlines() if l.startswith("!") or "Undefined" in l)[:1500])
            break
    return ok


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--guide", help="a guide page of the fork, e.g. PHYSICS.md")
    ap.add_argument("--out", default="build/oai")
    ap.add_argument("--pdf", action="store_true")
    ap.add_argument("--no-edges", action="store_true", help="skip the Lean import analysis (about 30 s)")
    ap.add_argument("--with-models", action="store_true", help="also load tools/dsl_models/*.py that extend a family")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    if not a.guide:
        ap.error("--guide is required")
    cat = catalog()
    path = pathlib.Path(a.guide)
    path = path if path.exists() else ROOT / a.guide
    g = parse_guide(path.read_text(), cat)
    name = _ident(path.stem) + "Guide"
    top = build_guide(g, cat, name=name, label=path.stem.lower())
    print(f"{path.name}: {len(g.sections)} sections, {sum(len(s.rows) for s in g.sections)} rows, "
          f"{len(top.__mro__) and len(guide_families(g, cat))} families; {len(g.problems)} disagreement(s) with the catalogue")
    outside = []
    if not a.no_edges:
        print("reading the Lean imports of every family (about 30 s) ...")
        outside = connect_families(guide_families(g, cat), cat)
    models = []
    if a.with_models:
        import runpy
        before = set(CATEGORIES)
        for f in sorted((HERE / "dsl_models").glob("*.py")):
            runpy.run_path(str(f), run_name="dsl_model")
        models = [c for c in CATEGORIES if c not in before and not c.__dict__.get("__oai_generated__")]
        for m in models:
            probs = check_model(cat, m)
            g.problems += probs
            print(f"model {m.__name__}: " + ("every source and 'uses' claim found in the Lean text" if not probs else f"{len(probs)} problem(s)"))
    tex = render_document(top, cat, outside, title=g.title, models=models, problems=g.problems)
    out = pathlib.Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    tf = out / (path.stem.lower() + ".tex")
    tf.write_text(tex)
    print(f"wrote {tf}")
    for p in g.problems:
        print("  guide vs catalogue:", p)
    if a.pdf:
        print("pdf ok" if build_pdf(tf) else "xelatex failed")


def selftest():
    n = [0]

    def ok(c, m=""):
        n[0] += 1
        assert c, m
    t = md_to_tex
    ok(t("Proves $`\\overline{\\mathbb Q}`$ & 50% of <i>p</i>-adic ℤ<sub><i>p</i></sub>")
       == "Proves $\\overline{\\mathbb Q}$ \\& 50\\% of \\emph{p}-adic \\ensuremath{ℤ}\\textsubscript{\\emph{p}}", t("Proves $`\\overline{\\mathbb Q}`$ & 50% of <i>p</i>-adic ℤ<sub><i>p</i></sub>"))
    ok(t("a_b and **bold** and `x_y`") == "a\\_b and \\textbf{bold} and \\texttt{x\\_y}", t("a_b and **bold** and `x_y`"))
    ok(t("[Lean](lean/docs/271.md)", families={"271"}) == "\\hyperref[fam:271]{Lean}")
    ok(t("[Lean](lean/docs/999.md)", families={"271"}).startswith("\\href{http") and "docs/999.md" in t("[Lean](lean/docs/999.md)", families=set()))
    ok(t("[p](preprints/A-B-2026/paper.pdf)", families=set()) == "\\href{" + UPSTREAM + "preprints/A-B-2026/paper.pdf}{p}")
    ok(t("[x](https://doi.org/10.1/a_b#c)").startswith("\\href{https://doi.org/10.1/a_b\\#c}"))
    ok(t("[Heisenberg](../ComparatorChallenges/Heisenberg.lean)", families=set()) == "\\href{" + UPSTREAM + "lean/ComparatorChallenges/Heisenberg.lean}{Heisenberg}")
    ok(t("$`\\binom{N}{r}`$ [x](u.md)", families=set()).startswith("$\\binom{N}{r}$ \\href{"))
    ok("\x00" not in t("[a $`x`$ b](https://e.org) $`y`$ `z`"))
    ok(t("$`\\int_\\mathbb R f\\lt 1`$") == "$\\int_{\\mathbb R} f\\lt 1$", t("$`\\int_\\mathbb R f\\lt 1`$"))
    ok(t("$`x^\\mathcal N`$") == "$x^{\\mathcal N}$")
    v = fix_verbatim("a \\begin{verbatim}\nx ∈ 𝓝 y\n\\end{verbatim} b \\begin{verbatim}\nplain ℝ\n\\end{verbatim}")
    ok("\\begin{Verbatim}[commandchars=\ue000\ue001\ue002]\nx ∈ \ue000fb\ue001𝓝\ue002 y\n\\end{Verbatim}\n" in v
       and "\\begin{verbatim}\nplain ℝ\n\\end{verbatim}" in v, v)                      # only blocks that need it change
    ok(plain_title("see [The Title](preprints/dir/) now") == "see The Title now", plain_title("see [The Title](preprints/dir/) now"))
    ok(plain_title("A $`\\overline{\\mathbb Q}`$ <i>b</i>") == "A \\overline\\mathbb Q b".replace("\\", ""), plain_title("A $`\\overline{\\mathbb Q}`$ <i>b</i>"))
    ok(blocks_to_tex("| E | R |\n|---|---|\n| [1] | row one |\n").count("\\item") == 1 and "row one" in blocks_to_tex("| E | R |\n|---|---|\n| [1] | row one |\n"))
    ok("\\item a" in blocks_to_tex("- a\n- b\n") and "\\item b" in blocks_to_tex("- a\n- b\n"))
    # a tiny catalogue: classes, inheritance, edges, guide parsing and the checks against the catalogue
    import tempfile
    sys.path.insert(0, str(HERE))
    import oai_catalog
    with tempfile.TemporaryDirectory() as d:
        r = pathlib.Path(d)
        for p in ("preprints/A-One-1-2026", "preprints/B-Two-2-2026"):
            (r / p).mkdir(parents=True)
        (r / "lean" / "docs").mkdir(parents=True)
        (r / "overview.tex").write_text("\\cataloguesection{Algebra}{2}\n\\resultentry{001}{One}{d}{l}\n\\resultentry{002}{Two}{d}{l}\n")
        (r / "CONTENTS.md").write_text("**001. First.** Sum $`x`$.\n\n&emsp;[Paper A](preprints/A-One-1-2026/p.pdf)\n\nAbs A.\n\n"
                                       "**002. Second.** Plain.\n\n&emsp;[Paper B](preprints/B-Two-2-2026/p.pdf)\n\nAbs B.\n")
        saved = list(CATEGORIES)
        _TREE.clear()
        try:
            c = Catalog(r)
            f1, f2 = family("1", c), family("002", c)
            ok(f1 is family("001", c) and f1.__name__ == "F001" and issubclass(f1, discipline_class("Algebra")) and issubclass(f1, root_class()))
            ok(f1 is not f2 and f1.__oai_family__ == "001")
            try:
                family("999", c)
                ok(False)
            except KeyError:
                ok(True)
            md = ("# Guide\n\nIntro.\n\n## 1. Things\n\nWhat.\n\n| Entry | Result | Notes |\n|---|---|---|\n"
                  "| [001] | [One](preprints/A-One-1-2026/) | fine |\n| [002] | [Two](preprints/A-One-1-2026/) | wrong entry |\n"
                  "| none | [Z](preprints/B-Two-2-2026/) | missing entry |\n| [077] | nothing | not a family |\n\n## Caveats\n\nBe careful.\n")
            g = parse_guide(md, c)
            ok(g.title == "Guide" and len(g.sections) == 1 and len(g.sections[0].rows) == 4 and g.tail[0][0] == "Caveats", g)
            ok(len(g.problems) == 3 and any("002" in p and "001" in p for p in g.problems) and any("says no entry" in p for p in g.problems)
               and any("[077]" in p for p in g.problems), g.problems)
            ok(guide_families(g, c) == ["001", "002"], guide_families(g, c))
            top = build_guide(g, c, name="G1")
            ok(issubclass(top, f1) and issubclass(top, f2) and len(top.__sections__) == 1 and issubclass(top.__sections__[0], f1))
            doc = render_document(top, c, problems=g.problems)
            ok("\\section{Guide}" in doc and "\\label{fam:001}" in doc and "Sum $x$." in doc and "\\hyperref[fam:001]{[001]}" in doc
               and "Checks of the guide" in doc and "\\section{Caveats}" in doc and doc.count("\\begin{document}") == 1, doc)
            (r / "lean" / "docs" / "001.md").write_text("# T\n\n## Scope\n\nS.\n\n| a | [A](../ComparatorChallenges/AA.lean) |\n")
            (r / "lean" / "ComparatorChallenges").mkdir(parents=True)
            (r / "lean" / "OAI" / "X").mkdir(parents=True)
            (r / "lean" / "ComparatorChallenges" / "AA.json").write_text('{"solution_module": "OAI.X.Main1", "theorem_names": ["OAI.AA.t"]}')
            (r / "lean" / "ComparatorChallenges" / "AA.lean").write_text("theorem t : True := trivial\n")
            (r / "lean" / "OAI" / "X" / "Main1.lean").write_text(
                "namespace OAI\nnamespace N\nlemma lem (x : Nat) : x = x := rfl\ntheorem main (x : Nat) : x = x := by\n  exact lem x\n"
                "theorem unrelated : True := trivial\nend N\nend OAI\n")
            c = Catalog(r)
            _TREE.clear()
            f1 = family("001", c)
            ok(list(f1.results)[0].name == "OAI.AA.t" and list(f1.results)[0].source == "OAI.X.Main1")
            M = lean.category("m")(type("M1", (f1,), {}))
            @M.lemmas(source="OAI.N.lem")
            def lem(x: leandsl.Real, y: leandsl.Real):
                return x <= y

            @M.theorems(source="OAI.N.main")
            def main_(x: leandsl.Real, y: leandsl.Real):
                return x <= y + 1

            @M.theorems(source="OAI.N.unrelated")
            def unrel(x: leandsl.Real, y: leandsl.Real):
                return x <= y + 2

            @M.theorems(source="OAI.N.ghost")
            def ghost(x: leandsl.Real, y: leandsl.Real):
                return x <= y + 3
            main_.connect(lem)
            ok(check_model(c, M) == [f"M1.ghost: source OAI.N.ghost is not declared in the Lean files behind family 001"], check_model(c, M))
            unrel.connect(lem)                                               # a claim the Lean text does not bear out
            probs = check_model(c, M)
            ok(len(probs) == 2 and "unrel is connected as using lem, but lem does not occur in its proof" in probs[1], probs)
            ok(check_model(c, discipline_class("Algebra")) == ["D_Algebra: does not extend a family class"])
        finally:
            CATEGORIES[:] = saved
            _TREE.clear()
            _LAST_KNOWN["families"] = ()
    print(f"oai_tree: {n[0]} checks ok")


if __name__ == "__main__":
    sys.modules.setdefault("oai_tree", sys.modules["__main__"])   # model files say `from oai_tree import family`: share the classes
    main()
