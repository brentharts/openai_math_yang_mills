"""N-plaquette SU(2) strip (N = 3, 4, 5) in the spin-network basis: the generalisation of su2_lattice_strip.py.

Graph: a ladder of squares. Labels in k = 2j: j_a, j_c (the three outer links of the end plaquettes), t_i, b_i for the interior
plaquettes i = 2..N-1 and the rungs v_1..v_{N-1}: 3N-3 labels. Each plaquette operator W_p = chi_{1/2}(U_p) shifts every link on its
boundary by +-1/2; an end plaquette has two signed 6j symbols (one per vertex on the loop, third legs t and b of the neighbour),
an interior plaquette four. Phases: the rule of the three-plaquette strip, extended by a factor (-1)^((t_q - b_q)/2) for every
interior neighbour q of an interior plaquette. This is a conjecture fixed by comparison, not a derivation: it is tested here
against an independent Peter-Weyl construction (models/su2_peter_weyl.py, wavefunctions of the N loop holonomies, no 6j symbols,
no phases) for N = 4 and 5, and for N = 3 against su2_lattice_strip.py.
"""
import sys
from pathlib import Path
import itertools
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as sla

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register
from models.su2_lattice_strip import half6j_exact
from models import su2_lattice_strip as strip3
from models import su2_peter_weyl as pw

NAME = "su2_lattice_stripn"
ENTRY = "su2-lattice-gap"


def tri(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) % 2 == 0


def vert(x, xp, y, yp, z):
    sg, sq = half6j_exact(z, yp, xp, x, y)
    return sg * float(sq) ** 0.5


class Net:
    def __init__(self, N):
        self.N = N
    def unpack(self, s):
        """s = (ja, jc, t2..t_{N-1}, b2..b_{N-1}, v1..v_{N-1}) in k -> arrays t[1..N], b[1..N], v[0..N] (index = plaquette/rung)."""
        N = self.N
        ja, jc = s[0], s[1]
        mid = N - 2
        t = [None, ja] + list(s[2:2 + mid]) + [jc]
        b = [None, ja] + list(s[2 + mid:2 + 2 * mid]) + [jc]
        v = [ja] + list(s[2 + 2 * mid:2 + 2 * mid + N - 1]) + [jc]
        return t, b, v
    def pack(self, t, b, v):
        N = self.N
        return (t[1], t[N]) + tuple(t[2:N]) + tuple(b[2:N]) + tuple(v[1:N])
    def valid(self, s):
        t, b, v = self.unpack(s)
        N = self.N
        if min(s) < 0: return False
        for i in range(1, N):
            if not tri(t[i], t[i + 1], v[i]) or not tri(b[i], b[i + 1], v[i]):
                return False
        return True
    def e_el(self, s):
        t, b, v = self.unpack(s)
        N = self.N
        c = lambda k: k * (k + 2) / 4
        e = 3 * c(t[1]) + 3 * c(t[N])
        e += sum(c(t[i]) + c(b[i]) for i in range(2, N))
        e += sum(c(v[i]) for i in range(1, N))
        return e
    def states(self, lam):
        N = self.N
        kmax = int(np.sqrt(4 * lam)) + 2
        out = []
        n = 3 * N - 3
        def rec(prefix):
            # brute force with pruning on energy only; validity checked at end (small spaces)
            if len(prefix) == n:
                if self.valid(prefix) and self.e_el(prefix) <= lam + 1e-12:
                    out.append(prefix)
                return
            for k in range(kmax):
                new = prefix + (k,)
                # crude energy prune: ja,jc weighted 3, others 1
                w = [3, 3] + [1] * (n - 2)
                e = sum(wi * (ki * (ki + 2) / 4) for wi, ki in zip(w, new))
                if e <= lam + 1e-12:
                    rec(new)
        rec(())
        return out
    def amp(self, p, s, s2, rule):
        t, b, v = self.unpack(s); t2, b2, v2 = self.unpack(s2)
        N = self.N
        if p == 1 or p == N:
            if p == 1:
                x, xp, y, yp = t[1], t2[1], v[1], v2[1]
                zt, zb = t[2], b[2]
                ztn, zbn = t[2], b[2]
            else:
                x, xp, y, yp = t[N], t2[N], v[N - 1], v2[N - 1]
                zt, zb = t[N - 1], b[N - 1]
            val = vert(x, xp, y, yp, zt) * vert(x, xp, y, yp, zb) * np.sqrt((x + 1) * (xp + 1) * (y + 1) * (yp + 1))
            return val * rule('end', p, s, s2, self)
        val = (vert(t[p], t2[p], v[p - 1], v2[p - 1], t[p - 1]) * vert(t[p], t2[p], v[p], v2[p], t[p + 1])
               * vert(b[p], b2[p], v[p], v2[p], b[p + 1]) * vert(b[p], b2[p], v[p - 1], v2[p - 1], b[p - 1]))
        dims = (t[p] + 1) * (t2[p] + 1) * (b[p] + 1) * (b2[p] + 1) * (v[p - 1] + 1) * (v2[p - 1] + 1) * (v[p] + 1) * (v2[p] + 1)
        return val * np.sqrt(dims) * rule('mid', p, s, s2, self)
    def transitions(self, p, s):
        t, b, v = self.unpack(s)
        N = self.N
        out = []
        if p == 1 or p == N:
            idx = [0, 2 + 2 * (N - 2)] if p == 1 else [1, 2 + 2 * (N - 2) + N - 2]
            for d in itertools.product((1, -1), repeat=2):
                s2 = list(s)
                for i, dd in zip(idx, d): s2[i] += dd
                s2 = tuple(s2)
                if min(s2) >= 0 and self.valid(s2): out.append(s2)
        else:
            mid = N - 2
            idx = [2 + (p - 2), 2 + mid + (p - 2)]  # t_p, b_p
            idx += [2 + 2 * mid + (p - 2), 2 + 2 * mid + (p - 1)]   # v_{p-1}, v_p
            for d in itertools.product((1, -1), repeat=4):
                s2 = list(s)
                for i, dd in zip(idx, d): s2[i] += dd
                s2 = tuple(s2)
                if min(s2) >= 0 and self.valid(s2): out.append(s2)
        return out
    def hamiltonian(self, lam, g2, rule):
        S = self.states(lam)
        ix = {s: i for i, s in enumerate(S)}
        n = len(S)
        rows, cols, vals = [], [], []
        for s in S:
            for p in range(1, self.N + 1):
                for s2 in self.transitions(p, s):
                    if s2 in ix:
                        a = self.amp(p, s, s2, rule)
                        if a != 0:
                            rows.append(ix[s2]); cols.append(ix[s]); vals.append(a)
        W = sps.csr_matrix((vals, (rows, cols)), shape=(n, n))
        H = sps.diags([0.5 * g2 * self.e_el(s) + 2 * self.N / g2 for s in S]) - W / g2
        return H.tocsr(), W, S



def rule_general(kind, p, s, s2, net):
    """end plaquette: (-1)^((t-b)/2) of its interior neighbour; interior: same/opposite shift sign times the neighbours' factors."""
    t, b, v = net.unpack(s); t2, b2, v2 = net.unpack(s2)
    N = net.N
    if kind == 'end':
        nb = 2 if p == 1 else N - 1
        return -1.0 if ((t[nb] - b[nb]) // 2) % 2 else 1.0
    sg = 1.0 if (t2[p] - t[p]) * (b2[p] - b[p]) > 0 else -1.0
    for q in (p - 1, p + 1):
        if 1 < q < N:
            sg *= -1.0 if ((t[q] - b[q]) // 2) % 2 else 1.0
    return sg


def lowest(N, lam, g2, k=6):
    H, W, S = Net(N).hamiltonian(lam, g2, rule_general)
    if len(S) <= k + 2:
        return np.sort(np.linalg.eigvalsh(H.toarray()))[:k], len(S)
    return np.sort(sla.eigsh(H, k=k, sigma=0, which='LM', return_eigenvectors=False)), len(S)


# values of the lowest levels at g^2 = 2 from the Peter-Weyl engine (computed in build_and_verify, repeated here only for reference)
def build_and_verify():
    eq = register(NAME, Equation(
        key="su2_stripn_basis", entry=ENTRY, title="N-plaquette strip: general-N spin-network Hamiltonian and its independent check",
        source_kind="derived", source="derived here; Kogut-Susskind, Phys. Rev. D 11, 395 (1975); closed 6j forms of Edmonds Table 5",
        latex_str=(r"H_N=\frac{g^{2}}{2}\sum_{\ell}C(j_\ell)+\frac{1}{g^{2}}\sum_{p=1}^{N}\bigl(2-W_p\bigr),\ \ "
                   r"W_p=\sigma_p\sqrt{\textstyle\prod d}\prod_{v\in p}\begin{Bmatrix}x_v&\tfrac12&x_v'\\y_v'&z_v&y_v\end{Bmatrix}"
                   r"%%BR%%\sigma_p=s_p\prod_{q\sim p,\ 1<q<N}(-1)^{(k_{t_q}-k_{b_q})/2}"),
        notes="Generalisation of the three-plaquette strip to N plaquettes with 3N-3 labels. The phase sigma_p is a conjecture fitted "
              "to the three-plaquette rule; it is tested by an independent construction that contains no 6j symbol and no phase.",
    ))
    # N = 3 against the existing strip module
    ref = strip3.levels(2.0, 24, nev=2)[0]
    ref_low = min(ref[1][0], ref[-1][0])
    w3, n3 = lowest(3, 24, 2.0)
    check(NAME, f"the general-N builder at N = 3 reproduces the three-plaquette strip module: lowest level at g^2 = 2 is "
                f"{w3[0]:.6f} against {ref_low:.6f}", abs(w3[0] - ref_low) < 1e-5, eq=eq)
    # N = 4, 5 against Peter-Weyl
    for N, lam, ecut in ((4, 20, 6), (5, 14, 6)):
        w, n = lowest(N, lam, 2.0, k=4)
        s = pw.Strip(N, ecut, ecut)
        wp = s.levels(2.0, 4)
        dev = float(np.max(np.abs(w - wp)))
        check(NAME, f"N = {N}: the lowest four levels at g^2 = 2 from the {n}-state spin-network builder ({w[0]:.5f}, {w[1]:.5f}, "
                    f"{w[2]:.5f}, {w[3]:.5f}) agree with the independent Peter-Weyl construction to {dev:.1e}", dev < 5e-3, eq=eq)
    # plaquette operators commute
    for N in (4, 5):
        net = Net(N)
        H, W, S = net.hamiltonian(16, 1.0, rule_general)
        ix = {s: i for i, s in enumerate(S)}
        n = len(S)
        Ws = []
        for p in range(1, N + 1):
            r, c, v = [], [], []
            for s in S:
                for s2 in net.transitions(p, s):
                    if s2 in ix:
                        a = net.amp(p, s, s2, rule_general)
                        r.append(ix[s2]); c.append(ix[s]); v.append(a)
            Ws.append(sps.csr_matrix((v, (r, c)), shape=(n, n)))
        # test on the states whose neighbours are all inside the cut (energy <= 4)
        inner = [i for i, s in enumerate(S) if net.e_el(s) <= 9]
        P = sps.csr_matrix((np.ones(len(inner)), (inner, range(len(inner)))), shape=(n, len(inner)))
        worst = 0.0
        for p in range(N):
            for q in range(p + 1, N):
                C = (Ws[p] @ Ws[q] - Ws[q] @ Ws[p]) @ P
                worst = max(worst, abs(C).max() if C.nnz else 0.0)
        # Hermiticity
        herm = max(abs(Wp - Wp.T).max() for Wp in Ws)
        check(NAME, f"N = {N}: the {N} plaquette operators commute on the {len(inner)} states away from the cut (largest commutator entry "
                    f"{worst:.1e}) and are symmetric (largest asymmetry {herm:.1e})", worst < 1e-10 and herm < 1e-10, eq=eq)
