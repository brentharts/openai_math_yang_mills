"""Connections the knowledge graph makes between islands of the library.

Two questions. First, the equations fall into connected pieces (through shared quantities, links, joins and
references): how many pieces are there, with and without the chirality repository's additions, and which
equations bridge them? Second, between pairs of results from different repositories that nobody linked
directly, what is the shortest chain of graph edges?
"""
import sys
from collections import deque
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import knowledgebase as kb
from explore import Result, main

# edges that carry mathematics: not the ones that only say where an equation was published or explained
MATH_EDGES = ("denotes", "shares", "field", "specialises", "defines", "limit-of")

PAIRS = (
    ("Weak mixing at unification", "Spectre inflation factor"),
    ("Net chiral class", "Chair44 contact count"),
    ("Weight vector demand", "Asymptotic-freedom mass scaling"),
    ("Half-spinor dimension", "Hat inflation factor"),
    ("Triality block roots", "Bloch magnon density"),
)


# a constant such as pi appears in half the library and links nothing in particular
HUBS = frozenset({"pi"})


def shortest(graph, start, goal, limit=8):
    """Breadth-first chain over the mathematical edges, never through a hub constant."""
    prev, frontier = {id(start): None}, [start]
    for _ in range(limit):
        nxt = []
        for node in frontier:
            for other in graph.neighbours(node, MATH_EDGES):
                if id(other) in prev or other.key in HUBS:
                    continue
                prev[id(other)] = node
                if other is goal:
                    chain = [other]
                    while prev[id(chain[-1])] is not None:
                        chain.append(prev[id(chain[-1])])
                    return chain[::-1]
                nxt.append(other)
        frontier = nxt
    return None


def components(graph, kinds=("equation", "quantity"), edge_kinds=None, drop=frozenset()):
    drop = frozenset(drop) | HUBS
    """Connected components over nodes of the given kinds, optionally ignoring some equations (by key)."""
    nodes = [n for n in graph.nodes.values() if n.kind in kinds and n.key not in drop]
    seen, comps = set(), []
    for n in nodes:
        if id(n) in seen:
            continue
        comp, queue = [], deque([n])
        seen.add(id(n))
        while queue:
            x = queue.popleft()
            comp.append(x)
            for y in graph.neighbours(x, edge_kinds):
                if y.kind in kinds and y.key not in drop and id(y) not in seen:
                    seen.add(id(y))
                    queue.append(y)
        comps.append(comp)
    return sorted(comps, key=len, reverse=True)


def run() -> Result:
    r = Result("kb_paths", "Chains between islands of the knowledge graph", prose="")
    g = kb.graph()
    lib = kb.library()
    chi = {e.name for e in lib.equations.values() if e.entry == "chirality-lisi"}
    eq_in = lambda comp: sum(1 for n in comp if n.kind == "equation")
    comps_all = [c for c in components(g, edge_kinds=MATH_EDGES) if eq_in(c)]
    comps_old = [c for c in components(g, edge_kinds=MATH_EDGES, drop=chi) if eq_in(c)]
    r.data["components with the additions (equations per piece, largest 6)"] = [eq_in(c) for c in comps_all[:6]]
    r.data["components without them"] = [eq_in(c) for c in comps_old[:6]]
    r.check("the graph has equation and quantity nodes", len(comps_all) >= 1 and eq_in(comps_all[0]) > 0)
    r.check(f"the chirality additions are {len(chi)} equations, all present in the graph",
            len(chi) >= 10 and all(g.get("equation", n) is not None for n in chi))
    singletons_old = sum(1 for c in comps_old if eq_in(c) == 1 and len(c) == 1)
    r.data["isolated equation nodes before / after"] = (singletons_old, sum(1 for c in comps_all if len(c) == 1))
    # the paths
    found = []
    for a, b in PAIRS:
        na, nb = g.get("equation", a), g.get("equation", b)
        if na is None or nb is None:
            raise AssertionError(f"kb_paths: no node for {a!r} or {b!r}")
        chain = shortest(g, na, nb)
        found.append((a, b, chain))
    reachable = [(a, b, c) for a, b, c in found if c]
    r.check(f"{len(reachable)} of {len(PAIRS)} chosen pairs from different repositories are joined by a chain of at most 8 edges",
            len(reachable) >= 3)
    for a, b, chain in reachable:
        r.data[f"{a} -> {b}"] = " -> ".join(f"{n.label if hasattr(n, 'label') else n.key}" for n in chain)
        r.check(f"the chain from '{a}' to '{b}' has {len(chain) - 1} edges and starts and ends at those equations",
                chain[0].key == a and chain[-1].key == b)
    # which edge kinds carry the chains
    kinds = {}
    for a, b, chain in reachable:
        for x, y in zip(chain, chain[1:]):
            for e in x.out + x.inn:
                if (e.src is y or e.dst is y):
                    kinds[e.kind] = kinds.get(e.kind, 0) + 1
                    break
    r.data["edge kinds along the chains"] = kinds
    r.prose = (r"Over equations and quantities, and only the edges that carry mathematics (quantities shared, fields, specialisations, definitions, limits), the graph with the chirality additions has "
               + str(len(comps_all)) + r" connected pieces that contain an equation (without them, " + str(len(comps_old)) + r"). Chains of at most eight edges join "
               + str(len(reachable)) + r" of the " + str(len(PAIRS)) + r" pairs chosen from different repositories, "
               r"for example the weak mixing angle and the Spectre's inflation factor, which no link states directly: "
               r"each chain is a path through shared quantities and links, never through where an equation was published, and it is a place to look for a connection, not a connection.")
    return r


if __name__ == "__main__":
    main(run())
