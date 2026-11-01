"""Read-only ingester for the OpenAI math release: disciplines, families, manuscripts, Lean scope notes,
Comparator challenges, and the links that can be read off the files themselves.

    python3 tools/oai_catalog.py --stats            # counts, and what could not be linked
    python3 tools/oai_catalog.py --json out.json    # the whole catalogue as JSON
    python3 tools/oai_catalog.py --selftest

Nothing here writes into the repository, and nothing is inferred from topic similarity. Where each fact
comes from:

  discipline, family number, short title   overview.tex   (\\cataloguesection / \\resultentry)
  headline, summary, manuscripts, abstracts CONTENTS.md
  Lean scope note, linked challenges        lean/docs/NNN.md
  solution module, theorem names            lean/ComparatorChallenges/X.json, statement text from X.lean
  "cites" (family -> family)                OAI:<directory> keys in preprints/*/build/*.tex or *.tex
  "shares" (family <-> family)              the OAI.* modules two families' solution modules both import
                                            (transitively); modules imported by 5 or more families are
                                            treated as infrastructure and ignored
  "mentions"                                "result NNN" in a family summary

What an edge means: ``cites`` is a bibliography entry; ``shares`` is shared Lean machinery, which is evidence
of related material but not of one result depending on the other; ``mentions`` is the summary saying so.
Manuscripts whose sources are absent (and the 128 families without a Lean note) simply have no such edges.
"""
from __future__ import annotations
import re, sys, json, pathlib, collections
from dataclasses import dataclass, field

ROOT = pathlib.Path(__file__).resolve().parent.parent


@dataclass
class Manuscript:
    dir: str                  # directory name under preprints/ (also its OAI:<key>)
    title: str
    pdf: str                  # path of the PDF, relative to the repository root
    abstract: str = ""
    family: str = ""


@dataclass
class Challenge:
    name: str                 # stem of the file in lean/ComparatorChallenges
    solution_module: str = ""
    theorem_names: list = field(default_factory=list)
    definition_names: list = field(default_factory=list)
    statements: dict = field(default_factory=dict)   # theorem name -> Lean statement text, when found
    family: str = ""


@dataclass
class Family:
    num: str
    title: str                # CONTENTS headline, e.g. "Milne’s rationality conjecture and algebraic specialization"
    tex_title: str = ""       # overview.tex short title (TeX)
    discipline: str = ""
    summary: str = ""         # CONTENTS description (markdown with $`...`$ math and a little HTML)
    manuscripts: list = field(default_factory=list)  # [Manuscript]
    lean_doc: str = ""        # scope note of lean/docs/NNN.md, markdown ("" when there is none)
    challenges: list = field(default_factory=list)   # [Challenge]


@dataclass(frozen=True)
class Edge:
    kind: str                 # "cites" | "shares" | "mentions"
    src: str                  # family number
    dst: str
    weight: int = 1
    evidence: tuple = ()      # short strings: manuscript directories, module names


class Catalog:
    def __init__(self, root=ROOT):
        self.root = pathlib.Path(root)
        self.disciplines: dict[str, list[str]] = {}      # discipline -> family numbers in catalogue order
        self.families: dict[str, Family] = {}
        self.manuscripts: dict[str, Manuscript] = {}     # by directory
        self.challenges: dict[str, Challenge] = {}
        self.unmapped_dirs: list[str] = []               # preprint directories no CONTENTS entry links to
        self.orphan_challenges: list[str] = []           # challenges no Lean note links to
        self.problems: list[str] = []                    # things the files disagree about
        self._read_overview()
        self._read_contents()
        self._read_docs()
        self._read_challenges()
        self._edges: list[Edge] | None = None

    # ---- the four sources -------------------------------------------------------------------
    def _read_overview(self):
        p = self.root / "overview.tex"
        if not p.exists():
            self.problems.append("overview.tex not found")
            return
        disc = ""
        for line in p.read_text(errors="ignore").splitlines():
            m = re.match(r"\\cataloguesection\{(.*?)\}\{\d+\}", line)
            if m:
                disc = m.group(1)
                self.disciplines.setdefault(disc, [])
                continue
            m = re.match(r"\\resultentry\{(\d+)\}", line)
            if m:
                title = _brace_arg(line, len(m.group(0)) - 1)
                f = Family(m.group(1), tex_title=title, title=title, discipline=disc)
                self.families[f.num] = f
                self.disciplines.setdefault(disc, []).append(f.num)

    def _read_contents(self):
        p = self.root / "CONTENTS.md"
        if not p.exists():
            self.problems.append("CONTENTS.md not found")
            return
        fam = ms = None
        want_abstract = False
        for line in p.read_text(errors="ignore").splitlines():
            m = re.match(r"\*\*(\d{3})\. (.*?)\*\* ?(.*)", line)
            if m:
                fam = self.families.get(m.group(1))
                if fam is None:
                    self.problems.append(f"CONTENTS entry {m.group(1)} is not in overview.tex")
                    fam = self.families[m.group(1)] = Family(m.group(1), title="")
                fam.title = m.group(2).rstrip(".").strip()
                fam.summary = m.group(3).strip()
                ms, want_abstract = None, False
                continue
            m = re.match(r"&emsp;\[(.*?)\]\((preprints/([^/]+)/[^)]+)\)", line)
            if m and fam is not None:
                ms = Manuscript(m.group(3), m.group(1), m.group(2), family=fam.num)
                self.manuscripts[ms.dir] = ms
                fam.manuscripts.append(ms)
                want_abstract = True
                continue
            if want_abstract and ms is not None and line.strip() and not line.startswith("<"):
                ms.abstract = line.strip()
                want_abstract = False
        pre = self.root / "preprints"
        if pre.is_dir():
            self.unmapped_dirs = sorted(d.name for d in pre.iterdir() if d.is_dir() and d.name not in self.manuscripts)

    def _read_docs(self):
        d = self.root / "lean" / "docs"
        if not d.is_dir():
            return
        for p in sorted(d.glob("*.md")):
            fam = self.families.get(p.stem)
            if fam is None:
                self.problems.append(f"lean/docs/{p.name} has no family")
                continue
            text = p.read_text(errors="ignore")
            m = re.search(r"^## Scope\s*\n(.*?)(?=^## |\Z)", text, re.S | re.M)
            fam.lean_doc = (m.group(1) if m else text).strip()
            fam._doc_challenges = re.findall(r"\(\.\./ComparatorChallenges/(\w+)\.lean\)", text)

    def _read_challenges(self):
        c = self.root / "lean" / "ComparatorChallenges"
        if not c.is_dir():
            return
        for p in sorted(c.glob("*.json")):
            try:
                j = json.loads(p.read_text())
            except ValueError:
                self.problems.append(f"{p.name} is not JSON")
                continue
            ch = Challenge(p.stem, j.get("solution_module", ""), list(j.get("theorem_names", [])),
                           list(j.get("definition_names", [])))
            lean = p.with_suffix(".lean")
            if lean.exists():
                src = lean.read_text(errors="ignore")
                for full in ch.theorem_names:
                    st = _statement(src, full)
                    if st:
                        ch.statements[full] = st
            self.challenges[ch.name] = ch
        linked = set()
        for fam in self.families.values():
            for name in getattr(fam, "_doc_challenges", []):
                if name in self.challenges:
                    self.challenges[name].family = fam.num
                    fam.challenges.append(self.challenges[name])
                    linked.add(name)
                else:
                    self.problems.append(f"lean/docs/{fam.num}.md links challenge {name}, which has no JSON")
        self.orphan_challenges = sorted(set(self.challenges) - linked)

    # ---- links read off the files -------------------------------------------------------------
    def edges(self) -> list[Edge]:
        if self._edges is None:
            self._edges = self._citation_edges() + self._mention_edges() + self._sharing_edges()
        return self._edges

    def _citation_edges(self) -> list[Edge]:
        dirs = set(self.manuscripts) | set(self.unmapped_dirs)
        by_dir: dict[tuple, list] = collections.defaultdict(list)
        self.manuscript_cites: dict[str, set] = collections.defaultdict(set)
        for d in sorted(dirs):
            base = self.root / "preprints" / d
            for t in list(base.glob("build/*.tex")) + list(base.glob("*.tex")):
                for key in re.findall(r"OAI:([A-Za-z0-9\-_.]+)", t.read_text(errors="ignore")):
                    key = key.rstrip(".,")
                    if key not in dirs:                     # keys can be glued to a year: accept an unambiguous prefix
                        cands = [x for x in dirs if x.startswith(key[:40])]
                        key = cands[0] if len(cands) == 1 else None
                    if key and key != d:
                        self.manuscript_cites[d].add(key)
        for a, bs in self.manuscript_cites.items():
            for b in bs:
                fa = self.manuscripts[a].family if a in self.manuscripts else ""
                fb = self.manuscripts[b].family if b in self.manuscripts else ""
                if fa and fb and fa != fb:
                    by_dir[(fa, fb)].append(f"{a} -> {b}")
        return [Edge("cites", a, b, len(ev), tuple(ev)) for (a, b), ev in sorted(by_dir.items())]

    def _mention_edges(self) -> list[Edge]:
        out = set()
        for n, f in self.families.items():
            for m in re.finditer(r"results? (\d{3})(?:(?:,| and| or|, and)\s*(?:results? )?(\d{3}))*", f.summary):
                for x in re.findall(r"\d{3}", m.group(0)):
                    if x != n and x in self.families:
                        out.add((n, x))
        return [Edge("mentions", a, b) for a, b in sorted(out)]

    def lean_closure(self, module: str, cache: dict | None = None) -> set:
        cache = {} if cache is None else cache
        seen, stack = set(), [module]
        while stack:
            m = stack.pop()
            if m in seen:
                continue
            seen.add(m)
            if m not in cache:
                p = self.root / "lean" / (m.replace(".", "/") + ".lean")
                cache[m] = re.findall(r"^import\s+(OAI\.[\w.]+)", p.read_text(errors="ignore"), re.M) if p.exists() else []
            stack += cache[m]
        return seen

    def family_closure(self, num: str) -> set:
        """The OAI.* modules the family's solution modules import, transitively (and those modules themselves)."""
        cache = self.__dict__.setdefault("_imports", {})
        mods: set = set()
        for ch in self.families[num].challenges:
            if ch.solution_module:
                mods |= self.lean_closure(ch.solution_module, cache)
        return mods

    def declaration(self, num: str, qualified: str):
        """(module, text) of the ``theorem|lemma|def`` named ``qualified`` (e.g. ``OAI.Heisenberg.Clocks.foo``) in a
        module of the family's import closure, or None. The namespace is read from the ``namespace`` lines of the
        file; the text runs to the next top-level declaration or ``end``. This is a textual check, not Lean's."""
        *prefix, short = qualified.split(".")
        for mod in sorted(self.family_closure(num)):
            p = self.root / "lean" / (mod.replace(".", "/") + ".lean")
            if not p.exists():
                continue
            text = p.read_text(errors="ignore")
            m = re.search(rf"^[ \t]*(?:(?:@\[[^\]]*\]|protected|private|noncomputable)\s+)*(?:theorem|lemma|def)\s+{re.escape(short)}(?![\w'])",
                          text, re.M)
            if not m:
                continue
            ns = re.findall(r"^namespace\s+([\w.]+)", text[:m.start()], re.M)
            if ".".join(prefix) in {".".join(ns[:k]) for k in range(len(ns) + 1)}:
                nxt = re.search(r"^(?:theorem|lemma|def|end|namespace|@\[)", text[m.end():], re.M)
                return mod, text[m.start():m.end() + nxt.start() if nxt else len(text)]
        return None

    def _sharing_edges(self, infrastructure: int = 5) -> list[Edge]:
        closure: dict[str, set] = {}
        for fam in self.families.values():
            mods = self.family_closure(fam.num)
            if mods:
                closure[fam.num] = mods
        self.family_modules = closure
        owners: dict[str, set] = collections.defaultdict(set)
        for f, mods in closure.items():
            for m in mods:
                owners[m].add(f)
        pairs: dict[tuple, list] = collections.defaultdict(list)
        for m, o in owners.items():
            if 2 <= len(o) < infrastructure:
                for a in o:
                    for b in o:
                        if a < b:
                            pairs[(a, b)].append(m)
        return [Edge("shares", a, b, len(ms), tuple(sorted(ms)[:5])) for (a, b), ms in sorted(pairs.items())]

    # ---- views ------------------------------------------------------------------------------
    def related(self, num: str, kinds=("cites", "shares", "mentions")) -> list[Edge]:
        return [e for e in self.edges() if e.kind in kinds and num in (e.src, e.dst)]

    def by_dir(self, d: str) -> Manuscript | None:
        return self.manuscripts.get(d)

    def stats(self) -> str:
        e = self.edges()
        kinds = collections.Counter(x.kind for x in e)
        in_edge = {n for x in e for n in (x.src, x.dst)}
        mc = sum(len(v) for v in self.manuscript_cites.values())
        lines = [
            f"disciplines {len(self.disciplines)}; families {len(self.families)}; manuscripts {len(self.manuscripts)}"
            f" (+{len(self.unmapped_dirs)} preprint directories no CONTENTS entry links to)",
            f"Lean notes {sum(bool(f.lean_doc) for f in self.families.values())}; challenges {len(self.challenges)}"
            f" ({sum(len(c.theorem_names) for c in self.challenges.values())} theorem names,"
            f" {sum(len(c.statements) for c in self.challenges.values())} statements found in the challenge files);"
            f" challenges no note links to: {self.orphan_challenges}",
            f"edges: {dict(kinds)} ({mc} manuscript-level citations); families with any edge: {len(in_edge)}",
        ]
        if self.problems:
            lines.append(f"files disagree: {self.problems[:5]}{' ...' if len(self.problems) > 5 else ''}")
        return "\n".join(lines)

    def to_json(self) -> dict:
        return {
            "disciplines": self.disciplines,
            "families": {n: {"title": f.title, "discipline": f.discipline, "manuscripts": [m.dir for m in f.manuscripts],
                             "challenges": [c.name for c in f.challenges], "has_lean_note": bool(f.lean_doc)}
                         for n, f in self.families.items()},
            "edges": [{"kind": e.kind, "src": e.src, "dst": e.dst, "weight": e.weight, "evidence": list(e.evidence)}
                      for e in self.edges()],
        }


def _brace_arg(s: str, i: int) -> str:
    """The balanced {...} argument of a TeX command starting at s[i] == '{'."""
    depth, out = 0, []
    for ch in s[i:]:
        if ch == "{":
            depth += 1
            if depth == 1:
                continue
        elif ch == "}":
            depth -= 1
            if depth == 0:
                break
        out.append(ch)
    return "".join(out)


def _statement(src: str, full: str) -> str:
    """Text of the ``theorem|lemma`` named ``full`` (e.g. ``OAI.Thomason.main``) up to the ``:=`` that ends its
    statement. A declaration may be written with only part of its name (inside ``namespace OAI``, as
    ``theorem Thomason.main``), so the longest tail of the name is tried first; the last component alone is
    the final resort."""
    parts = full.split(".")
    for k in range(len(parts), 0, -1):
        suffix = r"\.".join(re.escape(x) for x in parts[-k:])
        m = re.search(rf"^[ \t]*(?:(?:@\[[^\]]*\]|protected|private|noncomputable)\s+)*(?:theorem|lemma)\s+(?<![\w.'])(?:[\w']+\.)*{suffix}(?![\w'])",
                      src, re.M)
        if m:
            end = src.find(":=", m.end())
            return src[m.start():end if end != -1 else m.end() + 400].strip()
    return ""


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--root", default=str(ROOT))
    ap.add_argument("--stats", action="store_true")
    ap.add_argument("--json")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args(argv)
    if a.selftest:
        return selftest()
    cat = Catalog(a.root)
    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(cat.to_json(), indent=1))
    if a.stats or not a.json:
        print(cat.stats())


def selftest():
    import tempfile
    n = [0]

    def ok(c, m=""):
        n[0] += 1
        assert c, m

    ok(_brace_arg(r"\resultentry{001}{Zeta {and} L}{x}", len(r"\resultentry{001}")) == "Zeta {and} L")
    ok(_statement("namespace A\ntheorem foo (x : Nat) :\n    x = x := by\n  rfl\n", "A.foo") == "theorem foo (x : Nat) :\n    x = x")
    ok(_statement("theorem foobar : True := trivial", "foo") == "")        # whole-word match only
    two = "theorem Other.main : P := by\n  sorry\ntheorem Thomason.main : Q := by\n  sorry\n"
    ok(_statement(two, "OAI.Thomason.main") == "theorem Thomason.main : Q")  # the qualified tail picks the right one
    ok(_statement("private theorem T.x : R := by\n  sorry", "OAI.T.x") == "private theorem T.x : R")
    with tempfile.TemporaryDirectory() as d:
        r = pathlib.Path(d)
        (r / "preprints" / "A-One-May-1-2026" / "build").mkdir(parents=True)
        (r / "preprints" / "B-Two-May-2-2026" / "build").mkdir(parents=True)
        (r / "preprints" / "C-Old-May-3-2026").mkdir(parents=True)
        (r / "preprints" / "A-One-May-1-2026" / "build" / "p.tex").write_text(r"\bibitem{x} OpenAI, Two, OAI:B-Two-May-2-2026, 2026.")
        (r / "preprints" / "B-Two-May-2-2026" / "build" / "p.tex").write_text("no citations")
        (r / "lean" / "docs").mkdir(parents=True)
        (r / "lean" / "ComparatorChallenges").mkdir(parents=True)
        (r / "lean" / "OAI" / "X").mkdir(parents=True)
        (r / "overview.tex").write_text(
            "\\cataloguesection{Algebra}{2}\n\\resultentry{001}{One}{d}{l}\n\\cataloguesection{Topology}{4}\n\\resultentry{002}{Two}{d}{l}\n")
        (r / "CONTENTS.md").write_text(
            "**001. First family.** Uses result 002 heavily.\n\n&emsp;[Paper A](preprints/A-One-May-1-2026/p.pdf)\n\nAbstract A.\n\n"
            "**002. Second family.** Plain.\n\n&emsp;[Paper B](preprints/B-Two-May-2-2026/p.pdf)\n\nAbstract B.\n")
        (r / "lean" / "docs" / "001.md").write_text("# T\n\n## Scope\n\nScope text.\n\n## Comparator links\n\n| a | [A.lean](../ComparatorChallenges/AA.lean) |\n")
        (r / "lean" / "docs" / "002.md").write_text("# T\n\n## Scope\n\nScope two.\n\n| a | [B.lean](../ComparatorChallenges/BB.lean) |\n")
        for c, mod in (("AA", "OAI.X.Main1"), ("BB", "OAI.X.Main2"), ("ZZ", "OAI.X.Main1")):
            (r / "lean" / "ComparatorChallenges" / f"{c}.json").write_text(json.dumps(
                {"challenge_module": "C." + c, "solution_module": mod, "theorem_names": [f"OAI.{c}.thm"]}))
            (r / "lean" / "ComparatorChallenges" / f"{c}.lean").write_text(f"namespace OAI.{c}\ntheorem thm (n : Nat) : n = n := by sorry\n")
        (r / "lean" / "OAI" / "X" / "Main1.lean").write_text("import OAI.X.Shared\nimport Mathlib\n")
        (r / "lean" / "OAI" / "X" / "Main2.lean").write_text("import OAI.X.Shared\nimport OAI.X.Own\n")
        (r / "lean" / "OAI" / "X" / "Shared.lean").write_text(
            "import Mathlib\nnamespace OAI\nnamespace Sh.Ns\nlemma foo (n : Nat) : n = n := by\n  rfl\ntheorem bar : True := by\n  exact foo 1\nend Sh.Ns\nend OAI\n")
        (r / "lean" / "OAI" / "X" / "Own.lean").write_text("")
        c = Catalog(r)
        ok(list(c.disciplines) == ["Algebra", "Topology"] and c.disciplines["Topology"] == ["002"], c.disciplines)
        f1 = c.families["001"]
        ok(f1.title == "First family" and f1.discipline == "Algebra" and f1.summary.startswith("Uses result 002"), f1)
        ok([m.dir for m in f1.manuscripts] == ["A-One-May-1-2026"] and f1.manuscripts[0].abstract == "Abstract A.")
        ok(f1.lean_doc == "Scope text." and [x.name for x in f1.challenges] == ["AA"], f1.lean_doc)
        ok(c.challenges["AA"].statements == {"OAI.AA.thm": "theorem thm (n : Nat) : n = n"}, c.challenges["AA"])
        ok(c.orphan_challenges == ["ZZ"] and c.unmapped_dirs == ["C-Old-May-3-2026"])
        kinds = {(e.kind, e.src, e.dst) for e in c.edges()}
        ok(("cites", "001", "002") in kinds and ("mentions", "001", "002") in kinds and ("shares", "001", "002") in kinds, kinds)
        sh = [e for e in c.edges() if e.kind == "shares"][0]
        ok(sh.evidence == ("OAI.X.Shared",) and sh.weight == 1, sh)   # the family's own main modules are not shared
        ok(c.related("002") and not c.related("002", kinds=("nope",)))
        ok(json.loads(json.dumps(c.to_json()))["families"]["001"]["challenges"] == ["AA"])
        mod, txt = c.declaration("001", "OAI.Sh.Ns.foo")
        ok(mod == "OAI.X.Shared" and txt.startswith("lemma foo") and "bar" not in txt, txt)      # stops at the next declaration
        ok("foo" in c.declaration("001", "OAI.Sh.Ns.bar")[1])
        ok(c.declaration("001", "OAI.Sh.Other.foo") is None and c.declaration("001", "OAI.Sh.Ns.nope") is None)
        ok(c.declaration("001", "OAI.Sh.Ns.foo") == c.declaration("002", "OAI.Sh.Ns.foo"))      # also in 002's closure
        ok("families 2" in c.stats() and "edges" in c.stats(), c.stats())
    print(f"oai_catalog: {n[0]} checks ok")


if __name__ == "__main__":
    main()
