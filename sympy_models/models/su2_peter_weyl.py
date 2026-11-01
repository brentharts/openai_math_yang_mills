"""Peter-Weyl engine for the open N-plaquette SU(2) strip: wavefunctions of N loop holonomies, conjugation-invariant (penalty G^2).
Independent of the spin-network closed forms; used only to check them for N = 4, 5."""
import itertools, sys
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as sla
from sympy.physics.quantum.cg import CG
from sympy import Rational, S
from functools import lru_cache

def jmats(k):
    j = k / 2
    ms = np.arange(j, -j - 1, -1)          # m = j .. -j
    d = len(ms)
    Jz = np.diag(ms)
    Jp = np.zeros((d, d))
    for i in range(1, d):
        m = ms[i]
        Jp[i - 1, i] = np.sqrt(j * (j + 1) - m * (m + 1))
    Jx = (Jp + Jp.T) / 2
    Jy = (Jp - Jp.T) / (2j)
    return Jx, Jy, Jz, ms

@lru_cache(None)
def cgc(k1, m1, k2, m2, K, M):
    # <j1 m1 j2 m2 | J M>, arguments doubled
    return float(CG(Rational(k1, 2), Rational(m1, 2), Rational(k2, 2), Rational(m2, 2), Rational(K, 2), Rational(M, 2)).doit())

class Space:
    def __init__(self, N, ecut, kmax):
        self.N = N
        single = []
        for k in range(kmax + 1):
            for m2 in range(-k, k + 1, 2):
                for n2 in range(-k, k + 1, 2):
                    single.append((k, m2, n2))
        # enumerate product states with sum C(j) <= ecut
        C = lambda k: k / 2 * (k / 2 + 1)
        bykc = {}
        for s in single:
            bykc.setdefault(s[0], []).append(s)
        ks = sorted(bykc)
        states = []
        def rec(i, left, cur):
            if i == N:
                states.append(tuple(cur)); return
            for k in ks:
                if C(k) <= left + 1e-12:
                    for s in bykc[k]:
                        rec(i + 1, left - C(k), cur + [s])
        rec(0, ecut, [])
        self.states = states
        self.index = {s: i for i, s in enumerate(states)}
        self.dim = len(states)

def single_ops(kmax):
    """Hermitian L^a (on m) and R^a (on n) as dict k -> list of 3 matrices in the (m,n) product basis of dimension d*d, with index m*d+n; here m index order j..-j."""
    out = {}
    for k in range(kmax + 1):
        Jx, Jy, Jz, ms = jmats(k)
        d = k + 1
        I = np.eye(d)
        L = [np.kron(J, I) for J in (Jx, Jy, Jz)]          # acts on m
        # R: conjugate rep on n:  -J^*  (so that algebra holds)
        R = [np.kron(I, -J.conj()) for J in (Jx, Jy, Jz)]
        out[k] = (L, R)
    return out


# ------------------------------------------------------------------------------------------------
def local_ops(kmax):
    L, R, Cas, chi = {}, {}, {}, {}
    for k in range(kmax + 1):
        Jx, Jy, Jz, ms = jmats(k)
        d = k + 1
        I = np.eye(d)
        L[k] = [sps.csr_matrix(np.kron(J, I)) for J in (Jx, Jy, Jz)]
        R[k] = [sps.csr_matrix(np.kron(I, -J.conj())) for J in (Jx, Jy, Jz)]
        Cas[k] = k / 2 * (k / 2 + 1)
    # chi_{1/2} multiplication: <k',M,N| chi |k,m,n> = sqrt(d_j/d_J) sum_a CG(1/2 a, j m | J M) CG(1/2 a, j n | J N)
    for k in range(kmax + 1):
        for kp in (k - 1, k + 1):
            if kp < 0 or kp > kmax:
                continue
            d, dp = k + 1, kp + 1
            M = np.zeros((dp * dp, d * d))
            ms = [k - 2 * i for i in range(d)]          # doubled m values in the order j..-j
            mps = [kp - 2 * i for i in range(dp)]
            for ia, a in enumerate((1, -1)):             # doubled a = +-1
                for im, m in enumerate(ms):
                    Mv = a + m
                    if abs(Mv) > kp:
                        continue
                    c1 = cgc(1, a, k, m, kp, Mv)
                    if c1 == 0:
                        continue
                    for inn, n in enumerate(ms):
                        Nv = a + n
                        if abs(Nv) > kp:
                            continue
                        c2 = cgc(1, a, k, n, kp, Nv)
                        iM, iN = mps.index(Mv), mps.index(Nv)
                        M[iM * dp + iN, im * d + inn] += np.sqrt((k + 1) / (kp + 1)) * c1 * c2
            chi[(k, kp)] = sps.csr_matrix(M)
    return L, R, Cas, chi


class Strip:
    def __init__(self, N, ecut, kmax):
        self.N, self.kmax = N, kmax
        self.L, self.R, self.Cas, self.chi = local_ops(kmax)
        C = lambda k: k / 2 * (k / 2 + 1)
        tuples = [t for t in itertools.product(range(kmax + 1), repeat=N) if sum(C(k) for k in t) <= ecut + 1e-12]
        self.tuples = tuples
        self.dims = [int(np.prod([(k + 1) ** 2 for k in t])) for t in tuples]
        self.off = np.concatenate([[0], np.cumsum(self.dims)])
        self.dim = int(self.off[-1])
        self.tindex = {t: i for i, t in enumerate(tuples)}

    def embed(self, t, site, ops):
        """kron over sites of identity except ops at 'site' (single matrix)."""
        mats = []
        for s, k in enumerate(t):
            d2 = (k + 1) ** 2
            mats.append(ops if s == site else sps.identity(d2, format="csr"))
        out = mats[0]
        for m in mats[1:]:
            out = sps.kron(out, m, format="csr")
        return out

    def vec_ops(self, t, coeffs):
        """The three components of sum_site (cL * L_site + cR * R_site) for coeffs = {(site,'L'): c, (site,'R'): c}."""
        comps = []
        for a in range(3):
            acc = sps.csr_matrix((int(np.prod([(k + 1) ** 2 for k in t])),) * 2)
            for (site, typ), c in coeffs.items():
                op = (self.L if typ == "L" else self.R)[t[site]][a]
                acc = acc + c * self.embed(t, site, op)
            comps.append(acc)
        return comps

    def link_vectors(self):
        """E_l as coefficient dicts, 0-based sites; the open strip with tree = bottom row and all rungs."""
        N = self.N
        links = []
        for j in range(N):        # b_{j+1}: L_j + sum_{k>j} (L_k + R_k)
            d = {(j, "L"): 1.0}
            for k in range(j + 1, N):
                d[(k, "L")] = 1.0
                d[(k, "R")] = 1.0
            links.append(("b%d" % (j + 1), d))
        for j in range(N):        # t_j = R_j
            links.append(("t%d" % (j + 1), {(j, "R"): 1.0}))
        links.append(("v0", {(0, "R"): 1.0}))
        for j in range(1, N):     # v_j = L_{j} + R_{j+1}  (1-based j) -> 0-based sites j-1, j
            links.append(("v%d" % j, {(j - 1, "L"): 1.0, (j, "R"): 1.0}))
        links.append(("v%d" % N, {(N - 1, "L"): 1.0}))
        return links

    def hamiltonian(self, g2, mu=20.0, sign=+1.0):
        N = self.N
        rows, cols, vals = [], [], []
        diag_blocks = []
        links = self.link_vectors()
        if sign != 1.0:
            # flip the sign of the cross term of v_j by using L - R  (test of sign convention)
            links = [(n, {k: (c if k[1] == "L" else sign * c) for k, c in d.items()}) if n.startswith("v") and len(d) == 2 else (n, d) for n, d in links]
        blocks = []
        for ti, t in enumerate(self.tuples):
            K = None
            for name, d in links:
                comps = self.vec_ops(t, d)
                sq = sum(c @ c for c in comps)
                K = sq if K is None else K + sq
            G = self.vec_ops(t, {**{(s, "L"): 1.0 for s in range(N)}, **{(s, "R"): 1.0 for s in range(N)}})
            G2 = sum(c @ c for c in G)
            blocks.append((K, G2))
        # assemble block diagonal electric + penalty
        Kfull = sps.block_diag([b[0] for b in blocks], format="csr")
        G2full = sps.block_diag([b[1] for b in blocks], format="csr")
        # magnetic: chi_site
        rows, cols, vals = [], [], []
        W = sps.csr_matrix((self.dim, self.dim))
        data = []
        for ti, t in enumerate(self.tuples):
            for site in range(N):
                for kp in (t[site] - 1, t[site] + 1):
                    if kp < 0 or kp > self.kmax:
                        continue
                    t2 = t[:site] + (kp,) + t[site + 1:]
                    if t2 not in self.tindex:
                        continue
                    tj = self.tindex[t2]
                    mats = []
                    for s, k in enumerate(t):
                        d2 = (k + 1) ** 2
                        if s == site:
                            mats.append(self.chi[(k, kp)])
                        else:
                            mats.append(sps.identity(d2, format="csr"))
                    M = mats[0]
                    for m in mats[1:]:
                        M = sps.kron(M, m, format="csr")
                    M = M.tocoo()
                    rows.append(M.row + self.off[tj]); cols.append(M.col + self.off[ti]); vals.append(M.data)
        W = sps.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(self.dim, self.dim))
        self.K, self.G2, self.W = Kfull, G2full, W
        H = 0.5 * g2 * Kfull + (2 * N * sps.identity(self.dim) - W) / g2
        return H, G2full

    def levels(self, g2, nev=6, mu=20.0):
        H, G2 = self.hamiltonian(g2)
        Hp = (H + mu * G2).tocsr()
        w = sla.eigsh(Hp, k=nev, which="SA", ncv=max(40, 3 * nev), tol=1e-10, return_eigenvectors=False)
        return np.sort(w)


def build_and_verify():
    """Helper module with no registered equations; its own check: N = 2 reproduces the ladder's levels exactly (the ladder is a different,
    spin-network construction), which is the validation the N = 4, 5 comparisons in su2_lattice_stripn.py rest on."""
    if __package__ in (None, ""):
        sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))
    from eqlib import check
    from models import su2_lattice_ladder as ladder
    wl, _ = ladder.levels(2.0, 4)
    wp = Strip(2, 6, 6).levels(2.0, 4)
    dev = float(np.max(np.abs(np.sort(wl) - wp)))
    check("su2_peter_weyl", f"the Peter-Weyl engine at N = 2 reproduces the ladder's lowest four levels at g^2 = 2 "
                            f"({wp[0]:.6f}, {wp[1]:.6f}, {wp[2]:.6f}, {wp[3]:.6f}) to {dev:.1e}", dev < 1e-5)
