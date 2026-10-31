"""Three-plaquette SU(2) strip (companion to su2_lattice_ladder.py): the first test of the strip formula.

Three plaquettes a, b, c in a row, open boundary: 8 vertices, 10 links, Kogut-Susskind Hamiltonian in lattice units
H = (g^2/2) sum_l E_l^2 + (1/g^2) sum_p (2 - Tr U_p), the same convention as su2_lattice_ladder.py (constant 6/g^2 here, 4/g^2 there).

Basis. After the gauge reduction the four corner vertices have two links and force equal spins, so the gauge-invariant states are
spin networks on a graph with four trivalent vertices and six edges, labelled (k_a, k_c, k_t, k_b, K1, K2) in k = 2j:
  k_a  the three outer links of plaquette a (all equal),  k_c  the same for c,
  k_t  top link of b,  k_b  bottom link of b,  K1  link shared by a and b,  K2  link shared by b and c,
with the triangle rule at (j_a, t, J1), (t, j_c, J2), (j_a, b, J1), (b, j_c, J2). The ladder had one trivalent vertex pair and one
label per plaquette; here the middle plaquette has four links and so two labels of its own (t and b), which are equal modulo 2
but not equal. The electric energy is (g^2/2)[3C(j_a) + 3C(j_c) + C(t) + C(b) + C(J1) + C(J2)], C(j) = j(j+1).

The plaquette operators. W_p = chi_{1/2}(U_p) is a sum over spin shifts of +-1/2 on the links of the plaquette's boundary with
one factor per vertex on the loop, exactly as in the ladder, but with signs that are not all positive:
  W_a:  (j_a, J1) -> (j_a +- 1/2, J1 +- 1/2),  amplitude  s_ac sqrt(d_a d_a' d_J1 d_J1') {j_a 1/2 j_a'; J1' t J1} {j_a 1/2 j_a'; J1' b J1}
  W_c:  the same with (j_c, J2) and the legs t, b
  W_b:  (t, b, J1, J2) -> each +-1/2, amplitude  s_b sqrt(d_t d_t' d_b d_b' d_J1 d_J1' d_J2 d_J2') times four 6j symbols,
        one for each vertex on the loop: {t 1/2 t'; J1' j_a J1} {t 1/2 t'; J2' j_c J2} {b 1/2 b'; J2' j_c J2} {b 1/2 b'; J1' j_a J1}
with the phases s_ac = (-1)^((k_t - k_b)/2) for W_a and W_c, and s_b = +1 if t and b shift in the same direction, -1 if in opposite
directions, for W_b. The signed 6j symbols with an entry 1/2 are the closed forms of Edmonds Table 5. The ladder's squared 6j is
the special case where both vertices on the loop have the same third leg (then s = +1 and the product is a square).

How the phases were fixed. They are not derived here from a convention; they were found against an independent construction and
then checked on everything below. The independent construction (the "engine" below) writes the gauge-fixed wavefunction in the
three loop holonomies U_a, U_b, U_c with explicit Wigner 3j tensors at the four vertices and the epsilon pairing on the three tree
links, multiplies by chi_{1/2}(U_p) using Clebsch-Gordan series, and projects by Schur orthogonality. It shares no code with the
closed forms. The phases are exactly those that make the closed form equal the engine entry by entry, in the engine's basis.

Checked: the signed closed forms for {a b c; 1/2 e f} against sympy's 6j symbols (exact: sign and rational square), the full
amplitudes against the engine for every transition of every state with spins up to 3/2, and for a sample with spins up to 2 (both
the values, to 1e-12, and the sets of reachable states); the 3j tensors are invariant under random rotations; Hermiticity; that
W_a, W_b, W_c commute on the interior of the truncated space ([W_a,W_b] and [W_b,W_c] are not trivially zero: they act on shared
links); that reflection a <-> c exchanges W_a and W_c; that on the subspace j_c = 0 (where t = b = J2 and the strip is the ladder)
W_a and W_b are exactly the ladder's operators; truncation convergence of the lowest levels in the energy cut; the harmonic
(weak-coupling) limit of the three lowest even and odd levels from three coupled three-dimensional oscillators with kinetic matrix
tridiagonal(4, -1) of eigenvalues 4 - 2cos(k pi/4); the strip formula 2 sqrt(2 - cos(pi/(N+1))) of su2_ladder_weak for N = 3; the exact
strong-coupling splitting of the triply degenerate first excited level by degenerate second-order perturbation theory; the gap
against the ladder's and the one-plaquette gaps.

Not tested: four or more plaquettes; plaquettes in two dimensions (the L-shaped cluster, where a four-valent vertex first appears:
in a strip every vertex stays trivalent after the reduction, so there is no volume operator on the strip); periodic boundaries;
the first correction in g^2 at weak coupling beyond the extrapolation used; higher-order strong coupling; any continuum limit. The
weak-coupling numbers are extrapolations from g^2 = 0.2 .. 0.4, not computations at g -> 0, and the third and fourth even levels
(2 omega_2 and omega_1 + omega_3, nearly degenerate) are checked only through their sum. Nothing here bears on existence or a mass
gap of Yang-Mills on R^4.
"""
import itertools
import sys
from fractions import Fraction as Fr
from functools import lru_cache
from pathlib import Path

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as sla
import sympy as sp
from scipy.linalg import expm
from sympy.physics.wigner import wigner_6j

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register
from models import su2_lattice_ladder as ladder
from models.su2_lattice_gap import gap_num as plaquette_gap

NAME = "su2_lattice_strip"
ENTRY = "su2-lattice-gap"
KS = "Kogut-Susskind, Phys. Rev. D 11, 395 (1975)"
EDM = "Edmonds, Angular Momentum in Quantum Mechanics (1957), Table 5 (6j symbols with an entry 1/2)"

# state = (k_a, k_c, k_t, k_b, K1, K2), spins as k = 2j
M = 128                      # base for the integer key of a state (all spins stay below it)


# ---------------------------------------------------------------------------------------------------------------
# exact signed 6j symbols with an entry 1/2 and exact amplitudes (sign, rational square)
# ---------------------------------------------------------------------------------------------------------------
def _tri_int(a, b, c):
    return abs(a - b) <= c <= a + b and (a + b + c) % 2 == 0


def half6j_exact(a, b, c, e, f):
    """{a b c; 1/2 e f} (k-notation) as (sign, square) with square a Fraction: Edmonds Table 5 and the symmetries."""
    if not (_tri_int(a, b, c) and _tri_int(a, e, f) and abs(e - c) == 1 and abs(f - b) == 1):
        return 0, Fr(0)
    S = (a + b + c) // 2
    if e == c - 1 and f == b + 1:
        sq, sg = Fr(((a + c - b) // 2) * ((a + b - c) // 2 + 1), (b + 1) * (b + 2) * c * (c + 1)), (-1) ** S
    elif e == c - 1 and f == b - 1:
        sq, sg = Fr((S + 1) * (S - a), b * (b + 1) * c * (c + 1)), (-1) ** S
    elif e == c + 1 and f == b - 1:
        sq, sg = Fr(((a + b - c) // 2) * ((a + c - b) // 2 + 1), (c + 1) * (c + 2) * b * (b + 1)), (-1) ** S
    else:
        sq, sg = Fr((S + 2) * (S - a + 1), (b + 1) * (b + 2) * (c + 1) * (c + 2)), -((-1) ** S)
    return (sg if sq else 0), sq


def vert_exact(x, xp, y, yp, z):
    """The 6j at a vertex whose loop legs go x -> x', y -> y' and whose third leg is z: {x 1/2 x'; y' z y}."""
    return half6j_exact(z, yp, xp, x, y)


def amp_exact(which, s, t):
    """(sign, square) of the matrix element <t|W_which|s>: sign * sqrt(square). Zero if the transition is not allowed."""
    ka, kc, kt, kb, K1, K2 = s
    ka2, kc2, kt2, kb2, K12, K22 = t
    if which in "ac":
        x, xp, y, yp = (ka, ka2, K1, K12) if which == "a" else (kc, kc2, K2, K22)
        s1, q1 = vert_exact(x, xp, y, yp, kt)
        s2, q2 = vert_exact(x, xp, y, yp, kb)
        ph = -1 if ((kt - kb) // 2) % 2 else 1
        return ph * s1 * s2, q1 * q2 * (x + 1) * (xp + 1) * (y + 1) * (yp + 1)
    parts = [vert_exact(kt, kt2, K1, K12, ka), vert_exact(kt, kt2, K2, K22, kc),
             vert_exact(kb, kb2, K2, K22, kc), vert_exact(kb, kb2, K1, K12, ka)]
    sg, sq = 1, Fr(1)
    for a_, b_ in parts:
        sg, sq = sg * a_, sq * b_
    ph = 1 if (kt2 - kt) * (kb2 - kb) > 0 else -1
    dims = (kt + 1) * (kt2 + 1) * (kb + 1) * (kb2 + 1) * (K1 + 1) * (K12 + 1) * (K2 + 1) * (K22 + 1)
    return ph * sg, sq * dims


def amp_value(which, s, t):
    sg, sq = amp_exact(which, s, t)
    return sg * sp.sqrt(sp.Rational(sq.numerator, sq.denominator))


def steps(which, s):
    """Candidate targets of W_which from s (spin shifts of +-1 in k, triangle rules at the vertices that carry the loop)."""
    ka, kc, kt, kb, K1, K2 = s
    out = []
    if which == "a":
        for d1, d2 in itertools.product((1, -1), repeat=2):
            ka2, K12 = ka + d1, K1 + d2
            if ka2 >= 0 and K12 >= 0 and ladder.allowed(ka2, kt, K12) and ladder.allowed(ka2, kb, K12):
                out.append((ka2, kc, kt, kb, K12, K2))
    elif which == "c":
        for d1, d2 in itertools.product((1, -1), repeat=2):
            kc2, K22 = kc + d1, K2 + d2
            if kc2 >= 0 and K22 >= 0 and ladder.allowed(kt, kc2, K22) and ladder.allowed(kb, kc2, K22):
                out.append((ka, kc2, kt, kb, K1, K22))
    else:
        for d in itertools.product((1, -1), repeat=4):
            kt2, kb2, K12, K22 = kt + d[0], kb + d[1], K1 + d[2], K2 + d[3]
            if min(kt2, kb2, K12, K22) >= 0 and ladder.allowed(ka, kt2, K12) and ladder.allowed(kt2, kc, K22) \
                    and ladder.allowed(ka, kb2, K12) and ladder.allowed(kb2, kc, K22):
                out.append((ka, kc, kt2, kb2, K12, K22))
    return out


def e_el(s):
    """Electric energy sum in units of g^2/2: 3C(j_a) + 3C(j_c) + C(t) + C(b) + C(J1) + C(J2), as an exact Fraction."""
    ka, kc, kt, kb, K1, K2 = s
    return Fr(3 * ka * (ka + 2) + 3 * kc * (kc + 2) + kt * (kt + 2) + kb * (kb + 2) + K1 * (K1 + 2) + K2 * (K2 + 2), 4)


# ---------------------------------------------------------------------------------------------------------------
# vectorized numerics: states below an electric-energy cut, plaquette operators, parity sectors
# ---------------------------------------------------------------------------------------------------------------
def build_states(lam):
    """All states with e_el <= lam, as an integer array sorted by key. The set is closed under nothing: it is a truncation, and
    since H is a sum of a diagonal part and the W's, its eigenvalues are Rayleigh-Ritz upper bounds that fall as lam grows."""
    c = lambda k: k * (k + 2) / 4
    kmax = int(np.sqrt(4 * lam)) + 2
    out = []
    for ka in range(kmax):
        ea = 3 * c(ka)
        if ea > lam:
            break
        for kc in range(kmax):
            ec = ea + 3 * c(kc)
            if ec > lam:
                break
            for kt in range(kmax):
                et = ec + c(kt)
                if et > lam:
                    break
                for kb in range(kt % 2, kmax, 2):
                    eb = et + c(kb)
                    if eb > lam:
                        break
                    for K1 in range(max(abs(ka - kt), abs(ka - kb)), min(ka + kt, ka + kb) + 1, 2):
                        e1 = eb + c(K1)
                        if e1 > lam:
                            break
                        for K2 in range(max(abs(kc - kt), abs(kc - kb)), min(kc + kt, kc + kb) + 1, 2):
                            if e1 + c(K2) <= lam:
                                out.append((ka, kc, kt, kb, K1, K2))
    S = np.array(out, dtype=np.int64)
    kk = _keys(S)
    o = np.argsort(kk)
    return S[o], kk[o]


def _keys(S):
    k = np.zeros(len(S), dtype=np.int64)
    for i in range(6):
        k = k * M + S[:, i]
    return k


def _tri(a, b, c):
    return (np.abs(a - b) <= c) & (c <= a + b) & (((a + b + c) % 2) == 0)


def half6j(a, b, c, e, f):
    """Signed {a b c; 1/2 e f} for integer arrays (k-notation): the float version of half6j_exact."""
    a, b, c, e, f = (np.asarray(x, dtype=np.int64) for x in (a, b, c, e, f))
    S = (a + b + c) // 2
    ok = _tri(a, b, c) & _tri(a, e, f) & (e >= 0) & (f >= 0)
    sgn = np.where(S % 2 == 0, 1.0, -1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        c1 = np.sqrt(((a + c - b) / 2) * ((a + b - c) / 2 + 1) / ((b + 1) * (b + 2) * c * (c + 1)))
        c2 = np.sqrt((S + 1) * (S - a) / (b * (b + 1) * c * (c + 1)))
        c3 = np.sqrt(((a + b - c) / 2) * ((a + c - b) / 2 + 1) / ((c + 1) * (c + 2) * b * (b + 1)))
        c4 = np.sqrt((S + 2) * (S - a + 1) / ((b + 1) * (b + 2) * (c + 1) * (c + 2)))
    v = np.where((e == c - 1) & (f == b + 1), sgn * c1,
                 np.where((e == c - 1) & (f == b - 1), sgn * c2,
                          np.where((e == c + 1) & (f == b - 1), sgn * c3, -sgn * c4)))
    return np.where(ok, np.nan_to_num(v), 0.0)


def _vert(x, xp, y, yp, z):
    return half6j(z, yp, xp, x, y)


def _amp(which, S, T):
    ka, kc, kt, kb, K1, K2 = (S[:, i] for i in range(6))
    ka2, kc2, kt2, kb2, K12, K22 = (T[:, i] for i in range(6))
    d = lambda k: (k + 1).astype(float)
    if which in "ac":
        x, xp, y, yp = (ka, ka2, K1, K12) if which == "a" else (kc, kc2, K2, K22)
        v = _vert(x, xp, y, yp, kt) * _vert(x, xp, y, yp, kb) * np.sqrt(d(x) * d(xp) * d(y) * d(yp))
        return v * np.where(((kt - kb) // 2) % 2 == 0, 1.0, -1.0)
    v = (_vert(kt, kt2, K1, K12, ka) * _vert(kt, kt2, K2, K22, kc) * _vert(kb, kb2, K2, K22, kc) * _vert(kb, kb2, K1, K12, ka))
    v = v * np.sqrt(d(kt) * d(kt2) * d(kb) * d(kb2) * d(K1) * d(K12) * d(K2) * d(K22))
    return v * np.where((kt2 - kt) * (kb2 - kb) > 0, 1.0, -1.0)


_SHIFTS = {"a": [((0, d1), (4, d2)) for d1 in (1, -1) for d2 in (1, -1)],
           "c": [((1, d1), (5, d2)) for d1 in (1, -1) for d2 in (1, -1)],
           "b": [((2, d0), (3, d1), (4, d2), (5, d3)) for d0 in (1, -1) for d1 in (1, -1) for d2 in (1, -1) for d3 in (1, -1)]}


def plaquette_ops(S, kk):
    """The three plaquette operators as sparse matrices on the states S (sorted by key kk); transitions that leave S are dropped."""
    n = len(S)
    mats = {}
    for w in "abc":
        R, C, V = [], [], []
        for sh in _SHIFTS[w]:
            T = S.copy()
            for col, dlt in sh:
                T[:, col] += dlt
            ok = (T.min(axis=1) >= 0) & (T.max(axis=1) < M)
            tk = np.zeros(n, dtype=np.int64)
            for i in range(6):
                tk = tk * M + np.where(ok, T[:, i], 0)
            pos = np.minimum(np.searchsorted(kk, tk), n - 1)
            ok &= kk[pos] == tk
            idx = np.nonzero(ok)[0]
            v = _amp(w, S[idx], T[idx])
            nz = v != 0
            R.append(idx[nz])
            C.append(pos[idx][nz])
            V.append(v[nz])
        mats[w] = sps.csr_matrix((np.concatenate(V), (np.concatenate(R), np.concatenate(C))), shape=(n, n))
    return mats


def _eel_array(S):
    c = lambda k: k * (k + 2) / 4
    return 3 * c(S[:, 0]) + 3 * c(S[:, 1]) + c(S[:, 2]) + c(S[:, 3]) + c(S[:, 4]) + c(S[:, 5])


def hamiltonian(lam, g2):
    S, kk = build_states(lam)
    W = plaquette_ops(S, kk)
    H = sps.diags(0.5 * g2 * _eel_array(S) + 6 / g2) - (W["a"] + W["b"] + W["c"]) / g2
    return H.tocsr(), S, kk, W


def reflection(S, kk):
    """Index of the reflected state (a <-> c): (k_a, k_c, k_t, k_b, K1, K2) -> (k_c, k_a, k_t, k_b, K2, K1)."""
    rk = _keys(S[:, [1, 0, 2, 3, 5, 4]])
    pos = np.searchsorted(kk, rk)
    assert np.all(kk[pos] == rk)
    return pos


def _sector_basis(perm, par):
    n = len(perm)
    idx = np.arange(n)
    r = np.nonzero(idx < perm)[0]
    m = len(r)
    rows = np.concatenate([r, perm[r]])
    cols = np.concatenate([np.arange(m), np.arange(m)])
    vals = np.concatenate([np.full(m, 1 / np.sqrt(2)), np.full(m, par / np.sqrt(2))])
    if par == 1:
        f = np.nonzero(idx == perm)[0]
        rows, cols, vals = (np.concatenate([rows, f]), np.concatenate([cols, m + np.arange(len(f))]),
                            np.concatenate([vals, np.ones(len(f))]))
        m += len(f)
    return sps.csr_matrix((vals, (rows, cols)), shape=(n, m))


def levels(g2, lam, nev=4, tol=1e-10):
    """Lowest levels of H in the two reflection sectors (a <-> c), as {+1: array, -1: array}, from exact sector bases
    (a Krylov start vector alone does not keep an eigensolver in one sector). Deterministic start vector."""
    H, S, kk, _ = hamiltonian(lam, g2)
    perm = reflection(S, kk)
    out = {}
    for par in (1, -1):
        B = _sector_basis(perm, par)
        Hs = (B.T @ H @ B).tocsr()
        m = Hs.shape[0]
        w = sla.eigsh(Hs, k=nev, which="SA", v0=np.ones(m) / np.sqrt(m), ncv=40, tol=tol, maxiter=20000,
                      return_eigenvectors=False)
        out[par] = np.sort(w)
    return out, len(S)


def excitations(g2, lam, nev=4):
    """The excitation energies above the ground state with their reflection parity, ascending: [(energy, parity), ...]."""
    out, n = levels(g2, lam, nev)
    e0 = out[1][0]
    lv = sorted((float(w - e0), p) for p in (1, -1) for w in out[p])
    return lv[1:], n


# ---------------------------------------------------------------------------------------------------------------
# independent construction from explicit Wigner 3j tensors and Clebsch-Gordan series (used only for checks)
# ---------------------------------------------------------------------------------------------------------------
def _jmats(k):
    j = k / 2
    m = -j + np.arange(k + 1)
    Jp = np.zeros((k + 1, k + 1))
    for i in range(k):
        Jp[i + 1, i] = np.sqrt(j * (j + 1) - m[i] * (m[i] + 1))
    return (Jp + Jp.T) / 2, (Jp - Jp.T) / 2j, np.diag(m)


def wigner_D(k, ang):
    Jx, Jy, Jz = _jmats(k)
    return expm(-1j * (ang[0] * Jx + ang[1] * Jy + ang[2] * Jz))


@lru_cache(None)
def eps_tensor(k):
    """The pairing epsilon^j_{m m'} = (-1)^(j-m) delta_{m,-m'}: D(g) epsilon D(g)^T = epsilon."""
    j = k / 2
    e = np.zeros((k + 1, k + 1))
    for i in range(k + 1):
        e[i, k - i] = (-1) ** int(round(j - (-j + i)))
    return e


@lru_cache(None)
def three_j(k1, k2, k3):
    """Wigner 3j symbol as an array t[m1, m2, m3] (all legs covariant), from the Clebsch-Gordan coefficients of the ladder module."""
    C = ladder.cg(k1, k2, k3)
    t = np.zeros((k1 + 1, k2 + 1, k3 + 1))
    for a in range(k1 + 1):
        for b in range(k2 + 1):
            for Mi in range(k3 + 1):
                if C[a, b, Mi] != 0:
                    t[a, b, k3 - Mi] = (-1) ** ((k1 - k2 - k3) // 2 + Mi) * C[a, b, Mi] / np.sqrt(k3 + 1)
    return t


def three_j_invariance(k1, k2, k3, seed=1):
    ang = np.random.default_rng(seed).normal(size=3)
    t = three_j(k1, k2, k3)
    rot = np.einsum("abc,Aa,Bb,Cc->ABC", t, wigner_D(k1, ang), wigner_D(k2, ang), wigner_D(k3, ang))
    return float(np.abs(rot - t).max()), float(np.linalg.norm(t))


@lru_cache(None)
def _R(k, kp):
    """chi_{1/2}(U) E^k(U)_{mn} = sum_{M,N} R[m,n,M,N] E^{k'}(U)_{MN} with E^k = D^k epsilon^k (Clebsch-Gordan series)."""
    C = ladder.cg(k, 1, kp)
    return np.einsum("maM,pao,pn,oN->mnMN", C, C, eps_tensor(k), eps_tensor(kp), optimize=True)


def check_R(k, seed=5):
    ang = np.random.default_rng(seed).normal(size=3)
    lhs = np.trace(wigner_D(1, ang)) * (wigner_D(k, ang) @ eps_tensor(k))
    rhs = sum(np.einsum("mnMN,MN->mn", _R(k, kp), wigner_D(kp, ang) @ eps_tensor(kp)) for kp in (k - 1, k + 1) if kp >= 0)
    return float(np.abs(lhs - rhs).max())


@lru_cache(None)
def _T(s):
    """The gauge-fixed wavefunction as a tensor T[m,n,r,s,p,q] against E^{j_a}_{mn}(U_a) E^{j_b}_{rs}(U_b) E^{j_c}_{pq}(U_c):
    3j symbols at the four vertices and epsilon on the three tree links (set to the identity by the gauge)."""
    ka, kc, kt, kb, K1, K2 = s
    T = np.einsum("mxy,wpz,nru,sqv,xw,yu,zv->mnrspq", three_j(ka, kt, K1), three_j(kt, kc, K2), three_j(ka, kb, K1),
                  three_j(kb, kc, K2), eps_tensor(kt), eps_tensor(K1), eps_tensor(K2), optimize=True)
    return T, float(np.linalg.norm(T))


_AX = {"a": (0, 1), "b": (2, 3), "c": (4, 5)}
_POS = {"a": 0, "c": 1, "b": 3}


def engine_W(which, s, tmax):
    """<s''|chi_{1/2}(U_which)|s> for all s'' with k_t'' <= tmax, from the explicit construction (Schur orthogonality)."""
    T, nT = _T(s)
    p = _POS[which]
    k = s[p]
    res = {}
    for kp in (k + 1, k - 1):
        if kp < 0:
            continue
        i, j = _AX[which]
        src = list("mnrspq")
        tgt = src.copy()
        tgt[i], tgt[j] = "M", "N"
        TR = np.einsum("".join(src) + "," + src[i] + src[j] + "MN->" + "".join(tgt), T, _R(k, kp), optimize=True)
        ns = list(s)
        ns[p] = kp
        for kt_ in range(0, tmax + 1):
            for K1_ in range(abs(ns[0] - kt_), ns[0] + kt_ + 1, 2):
                if not ladder.allowed(ns[0], ns[3], K1_):
                    continue
                for K2_ in range(abs(kt_ - ns[1]), kt_ + ns[1] + 1, 2):
                    if not ladder.allowed(ns[3], ns[1], K2_):
                        continue
                    s2 = (ns[0], ns[1], kt_, ns[3], K1_, K2_)
                    T2, n2 = _T(s2)
                    v = float(np.sum(T2 * TR)) / (n2 * nT) * np.sqrt((k + 1) / (kp + 1))
                    if abs(v) > 1e-12:
                        res[s2] = v
    return res


def engine_states(kmax):
    out = []
    for ka, kc, kb in itertools.product(range(kmax + 1), repeat=3):
        for kt in range(kmax + 1):
            for K1 in range(abs(ka - kt), ka + kt + 1, 2):
                if not ladder.allowed(ka, kb, K1):
                    continue
                for K2 in range(abs(kt - kc), kt + kc + 1, 2):
                    if ladder.allowed(kb, kc, K2):
                        out.append((ka, kc, kt, kb, K1, K2))
    return out


# ---------------------------------------------------------------------------------------------------------------
# strong coupling: exact degenerate second-order perturbation theory
# ---------------------------------------------------------------------------------------------------------------
STRING_A, STRING_B, STRING_C = (1, 0, 0, 0, 1, 0), (0, 0, 1, 1, 1, 1), (0, 1, 0, 0, 0, 1)


def neighbours_exact(s):
    """Exact matrix elements out of s under all three plaquette operators: {target: sympy number} (a target reached by
    two operators cannot occur: W_a, W_b, W_c change different labels)."""
    out = {}
    for w in "abc":
        for t in steps(w, s):
            v = amp_value(w, s, t)
            if v != 0:
                out[t] = out.get(t, 0) + v
    return out


def strong_matrix():
    """The effective second-order matrix M (units g^-6) on the strings (A, B, C) around plaquette a, b, c and the vacuum shift.

    H = (g^2/2) E + (1/g^2)(6 - W_a - W_b - W_c); the strings have E = 3, so E_D = 3 g^2/2. Then
    M_XY = sum_{m not in D} W_Xm W_mY / ((E_D - E_m)/(g^2/2)) and the vacuum (E = 0) shifts by -sum_m W_0m^2 / (E_m/2) = -2."""
    D = (STRING_A, STRING_B, STRING_C)
    ED = Fr(3)
    assert all(e_el(d) == ED for d in D)
    nb = {d: neighbours_exact(d) for d in D}
    vac = (0, 0, 0, 0, 0, 0)
    nv = neighbours_exact(vac)
    M_ = sp.zeros(3, 3)
    for i, X in enumerate(D):
        for j, Y in enumerate(D):
            tot = 0
            for m, vx in nb[X].items():
                if m in D:
                    raise AssertionError("first-order coupling inside the degenerate level")
                vy = nb[Y].get(m)
                if vy is None:
                    continue
                em = e_el(m)
                assert em != ED
                tot += vx * vy * 2 / sp.Rational((ED - em).numerator, (ED - em).denominator)
            M_[i, j] = sp.simplify(tot)
    vshift = sum(v ** 2 * 2 / sp.Rational(e_el(m).numerator, e_el(m).denominator) for m, v in nv.items())
    return M_, sp.simplify(-vshift)


# ---------------------------------------------------------------------------------------------------------------
def build_and_verify():
    # =====================================================================================================
    # basis and plaquette operators
    # =====================================================================================================
    eq_basis = register(NAME, Equation(
        key="su2_strip_basis", entry=ENTRY, title="Three-plaquette strip: basis, energies and plaquette operators",
        source_kind="derived", source=f"derived here; {KS} for H; {EDM} for the signed 6j closed forms",
        latex_str=(r"|j_a,j_c,t,b,J_1,J_2\rangle,\ \ \triangle(j_a,t,J_1),\,\triangle(t,j_c,J_2),\,\triangle(j_a,b,J_1),\,\triangle(b,j_c,J_2)"
                   r"%%BR%%H=\frac{g^{2}}{2}\bigl[3C(j_a)+3C(j_c)+C(t)+C(b)+C(J_1)+C(J_2)\bigr]"
                   r"+\frac{1}{g^{2}}\bigl(6-W_a-W_b-W_c\bigr)"
                   r"%%BR%%\langle j_a{+}\epsilon,J_1{+}\eta|W_a|j_a,J_1\rangle=(-1)^{(k_t-k_b)/2}\sqrt{d_{j_a}d_{j_a'}d_{J_1}d_{J_1'}}\,"
                   r"\begin{Bmatrix}j_a&\tfrac12&j_a'\\J_1'&t&J_1\end{Bmatrix}\begin{Bmatrix}j_a&\tfrac12&j_a'\\J_1'&b&J_1\end{Bmatrix}"
                   r"%%BR%%W_b:\ \ s_b\sqrt{\textstyle\prod d}\ \begin{Bmatrix}t&\tfrac12&t'\\J_1'&j_a&J_1\end{Bmatrix}"
                   r"\begin{Bmatrix}t&\tfrac12&t'\\J_2'&j_c&J_2\end{Bmatrix}"
                   r"\begin{Bmatrix}b&\tfrac12&b'\\J_2'&j_c&J_2\end{Bmatrix}\begin{Bmatrix}b&\tfrac12&b'\\J_1'&j_a&J_1\end{Bmatrix}"),
        notes="Spin network on a graph of four trivalent vertices and six edges (the corner vertices of the three plaquettes force "
              "equal spins on the outer links of a and of c). W_a, W_b, W_c are chi_{1/2} of the plaquette holonomies and shift every link "
              "on the plaquette's boundary by +-1/2. The signs: s_ac = (-1)^((k_t-k_b)/2) for W_a and W_c; s_b = +1 if t and b shift in "
              "the same direction and -1 if in opposite directions, for W_b. These were fixed against the independent Clebsch-Gordan "
              "construction in the module and hold on every transition checked; they are not derived from a convention. The "
              "6j symbols are signed here (the ladder's are squared because both vertices on its loop have the same third leg).",
    ))
    # --- exact signed 6j against sympy -------------------------------------------------------------------
    n6, bad6 = 0, 0
    for a, b, c in itertools.product(range(0, 6), repeat=3):
        if not _tri_int(a, b, c):
            continue
        for e, f in itertools.product((c - 1, c + 1), (b - 1, b + 1)):
            if e < 0 or f < 0 or not _tri_int(a, e, f):
                continue
            ref = wigner_6j(sp.Rational(a, 2), sp.Rational(b, 2), sp.Rational(c, 2), sp.Rational(1, 2),
                            sp.Rational(e, 2), sp.Rational(f, 2))
            sg, sq = half6j_exact(a, b, c, e, f)
            mine = sg * sp.sqrt(sp.Rational(sq.numerator, sq.denominator))
            bad6 += sp.simplify(mine - ref) != 0
            n6 += 1
    check(NAME, f"the closed forms for the signed 6j symbol {{a b c; 1/2 e f}} (Edmonds Table 5 and symmetries) equal sympy's exact "
                f"symbol, sign and rational square, for all {n6} allowed arguments with spins up to 5/2", bad6 == 0, eq=eq_basis)
    # --- 3j tensors and the Clebsch-Gordan series -------------------------------------------------------
    inv = max(three_j_invariance(*ks)[0] for ks in [(1, 1, 2), (2, 3, 3), (2, 2, 2), (3, 2, 1), (4, 3, 3), (5, 4, 3)])
    nrm = max(abs(three_j_invariance(*ks)[1] - 1) for ks in [(1, 1, 2), (2, 3, 3), (3, 2, 1)])
    rr = max(check_R(k) for k in range(0, 5))
    check(NAME, f"the independent construction is sound: the Wigner 3j tensors are invariant under a random rotation (largest "
                f"deviation {inv:.1e}) and have unit norm ({nrm:.0e}); chi_1/2 E^k = sum R E^k' reproduces the character "
                f"product for k up to 4 (largest deviation {rr:.1e})", inv < 1e-13 and nrm < 1e-13 and rr < 1e-13, eq=eq_basis)
    # --- engine against closed form ----------------------------------------------------------------------
    st3 = engine_states(3)
    sample4 = [s for s in engine_states(4) if max(s[:4]) == 4][::9]
    worst, ntr, setbad = 0.0, 0, 0
    for s in st3 + sample4:
        for w in "abc":
            eng = engine_W(w, s, 8)
            cl = {t: amp_exact(w, s, t) for t in steps(w, s)}
            cs = {t for t, (sg, sq) in cl.items() if sg != 0}
            setbad += cs != set(eng)
            for t, v in eng.items():
                sg, sq = cl[t]
                worst = max(worst, abs(v - sg * float(np.sqrt(float(sq)))))
                ntr += 1
    check(NAME, f"the closed-form amplitudes of W_a, W_b, W_c (signed 6j symbols with the stated phases) equal the independent "
                f"Wigner-3j / Clebsch-Gordan construction for all {ntr} transitions out of {len(st3)} states with spins up to 3/2 "
                f"and {len(sample4)} sampled states with spins up to 2 (largest difference {worst:.1e}); the sets of reachable "
                f"states agree for all {3 * (len(st3) + len(sample4))} (state, plaquette) pairs ({setbad} differ)",
          worst < 1e-12 and setbad == 0 and ntr > 4000, eq=eq_basis)
    # --- float vs exact -----------------------------------------------------------------------------------
    S, kk = build_states(24)
    W = plaquette_ops(S, kk)
    idx = {tuple(int(x) for x in s): i for i, s in enumerate(S)}
    dev, nel = 0.0, 0
    for w in "abc":
        Wc = W[w].tocoo()
        for i, j, v in list(zip(Wc.row, Wc.col, Wc.data))[::7]:
            sg, sq = amp_exact(w, tuple(int(x) for x in S[i]), tuple(int(x) for x in S[j]))
            dev = max(dev, abs(v - sg * float(np.sqrt(float(sq)))))
            nel += 1
    check(NAME, f"the vectorized float operators used for the spectra equal the exact amplitudes on a sample of {nel} elements "
                f"of the {len(S)}-state space e_el <= 24 (largest difference {dev:.1e})", dev < 1e-13, eq=eq_basis)
    # --- structure of the operators ----------------------------------------------------------------------
    symm = max(abs(W[w] - W[w].T).max() for w in "abc")
    perm = reflection(S, kk)
    refl = abs(W["a"][perm][:, perm] - W["c"]).max() + abs(W["b"][perm][:, perm] - W["b"]).max()
    check(NAME, f"W_a, W_b, W_c are symmetric (largest entry of W - W^T: {symm:.1e}) and reflection a <-> c maps W_a to W_c and "
                f"leaves W_b invariant (difference {refl:.1e})", symm == 0 and refl < 1e-14, eq=eq_basis)
    Hh, Sb, kb_, Wb_ = hamiltonian(160, 1.0)
    inner = _eel_array(Sb) <= 25
    cm = {}
    for nm, (X, Y) in {"ab": (Wb_["a"], Wb_["b"]), "bc": (Wb_["b"], Wb_["c"]), "ac": (Wb_["a"], Wb_["c"])}.items():
        Cm = (X @ Y - Y @ X).tocsr()
        cm[nm] = (float(abs(Cm[inner][:, inner]).max()), float(abs(Cm).max()))
    check(NAME, f"the plaquette operators commute: on the {int(inner.sum())} states with e_el <= 25, inside the {len(Sb)}-state space "
                f"e_el <= 160, the largest entries of [W_a,W_b], [W_b,W_c], [W_a,W_c] are {cm['ab'][0]:.1e}, {cm['bc'][0]:.1e}, "
                f"{cm['ac'][0]:.1e}; the same commutators are as large as {cm['ab'][1]:.2f} at the truncation boundary, so the check "
                f"is not vacuous", max(v[0] for v in cm.values()) < 1e-13 and cm["ab"][1] > 0.1, eq=eq_basis)
    # --- reduction to the ladder -----------------------------------------------------------------------
    S2, kk2 = build_states(80)
    W2 = plaquette_ops(S2, kk2)
    sel = np.nonzero((S2[:, 1] == 0) & (S2[:, 2] == S2[:, 3]) & (S2[:, 5] == S2[:, 2]))[0]
    lab = {(int(S2[i, 0]), int(S2[i, 2]), int(S2[i, 4])): j for j, i in enumerate(sel)}
    Wa_r, Wb_r = W2["a"][sel][:, sel], W2["b"][sel][:, sel]
    dl, nl = 0.0, 0
    for (ka, kb, K), j in lab.items():
        for which, t in ladder.steps(ka, kb, K, 400):
            if t in lab:
                mine = (Wa_r if which == "a" else Wb_r)[j, lab[t]]
                dl = max(dl, abs(ladder.element(which, (ka, kb, K), t) - mine))
                nl += 1
    check(NAME, f"j_c = 0 reduces the strip to the ladder: on the {len(sel)} states with k_c = 0 and k_t = k_b = K2 (the strip's three "
                f"outer links of b all carry j_b) W_a and W_b equal the ladder's operators, all {nl} transitions, largest difference "
                f"{dl:.1e}; the electric energies agree by construction (3C(j_b) + C(J1) + 3C(j_a))", dl < 1e-13 and nl > 500, eq=eq_basis)
    ea = levels(0.5, 80, 4)[0]
    eb = levels(0.5, 160, 4)[0]
    dconv = max(float(np.max(abs(ea[p] - eb[p]))) for p in (1, -1))
    check(NAME, f"truncation convergence at g^2 = 0.5: the lowest four levels in each reflection sector at e_el <= 80 and <= 160 agree "
                f"to 5e-5 (largest difference {dconv:.1e}); the larger space is the lower, as Rayleigh-Ritz requires",
          dconv < 5e-5 and all(np.all(eb[p] <= ea[p] + 1e-12) for p in (1, -1)), eq=eq_basis)

    # =====================================================================================================
    # weak coupling: three coupled oscillators
    # =====================================================================================================
    lam_k = [4 - 2 * sp.cos(sp.Rational(k, 4) * sp.pi) for k in (1, 2, 3)]
    om = [sp.sqrt(l / 2) for l in lam_k]
    pred = sorted([(2 * om[0], 1), (om[0] + om[1], -1), (om[0] + om[2], 1), (2 * om[1], 1), (om[1] + om[2], -1), (2 * om[2], 1)],
                  key=lambda x: float(x[0]))
    delta3 = 2 * sp.sqrt(2 - sp.cos(sp.pi / 4))
    eq_weak = register(NAME, Equation(
        key="su2_strip_weak", entry=ENTRY, title="Weak-coupling limit of the three-plaquette strip",
        source_kind="derived", source="derived here from su2_strip_basis (harmonic approximation); numerical check",
        latex_str=(r"\omega_k=\sqrt{2-\cos\tfrac{k\pi}{4}}\ \ (k=1,2,3)=\sqrt{\tfrac{4-\sqrt2}{2}},\ \sqrt2,\ \sqrt{\tfrac{4+\sqrt2}{2}}"
                   r"%%BR%%\Delta_{g\to0}=2\omega_1=2\sqrt{2-\tfrac{\sqrt2}{2}}\approx2.2741,\quad "
                   r"\text{next: }\omega_1{+}\omega_2\ (\text{odd}),\ \ \omega_1{+}\omega_3,\ \ 2\omega_2\ (\text{even})"
                   r"%%BR%%\frac{\Delta_3}{\Delta_{\mathrm{plaquette}}}\to\sqrt{\tfrac12\bigl(2-\tfrac{\sqrt2}{2}\bigr)}\approx0.8040"),
        notes="Near U = 1 each plaquette angle is a three-dimensional oscillator. The electric energy is (g^2/2) E^T K E with the "
              "kinetic matrix K = tridiag(4,-1) in the three loop electric fields (each outer link carries one loop's field, each "
              "shared link the difference of two), eigenvalues 4 - 2cos(k pi/4) = 4 - sqrt2, 4, 4 + sqrt2, and omega = sqrt(lambda/2). "
              "Gauge invariance keeps total angular momentum 0, so the lowest excitation is two quanta in the soft mode, 2 omega_1; "
              "mixed pairs of different modes are singlets as well, with reflection parity (-1)^(k1+k2) for modes k1, k2 (mode k has "
              "parity (-1)^(k+1) under a <-> c). This is the N = 3 case of the strip formula 2 sqrt(2 - cos(pi/(N+1))) stated "
              "in su2_ladder_weak. The 3rd and 4th even levels, omega_1 + omega_3 and 2 omega_2, differ by 0.046 and mix at "
              "finite g^2; only their sum is checked.",
    ))
    K3 = sp.Matrix([[4, -1, 0], [-1, 4, -1], [0, -1, 4]])
    ev = sorted(K3.eigenvals(), key=lambda x: float(x))
    check(NAME, "the kinetic matrix tridiag(4,-1) of order 3 has eigenvalues 4 - sqrt2, 4, 4 + sqrt2, so omega_k = sqrt(2 - cos(k pi/4)) "
                "= 1.1371, 1.4142, 1.6453", all(sp.simplify(a - b) == 0 for a, b in zip(ev, lam_k)) and
          abs(float(om[0]) - 1.13706) < 1e-5 and abs(float(om[1]) - 1.41421) < 1e-5 and abs(float(om[2]) - 1.64533) < 1e-5, eq=eq_weak)
    check(NAME, f"strip formula at N = 3: 2 sqrt(2 - cos(pi/4)) = {float(delta3):.4f} equals 2 omega_1; the sequence N = 1, 2, 3 is "
                f"2.8284, 2.4495, 2.2741, decreasing toward 2 as N grows",
          sp.simplify(delta3 - 2 * om[0]) == 0 and float(delta3) < float(2 * sp.sqrt(2 - sp.cos(sp.pi / 3))) < float(2 * sp.sqrt(2 - sp.cos(sp.pi / 2))),
          eq=eq_weak)
    gs = (0.2, 0.25, 0.3, 0.4)
    data = {g: excitations(g, int(round(40 / g)), 4)[0] for g in gs}
    xs = np.array(gs)
    extrap = []
    for i in range(5):
        extrap.append(float(np.polyfit(xs, np.array([data[g][i][0] for g in gs]), 3)[-1]))
    par_num = [data[0.2][i][1] for i in range(5)]
    target = [float(p[0]) for p in pred[:5]]
    check(NAME, f"reflection parities of the five lowest excitations at g^2 = 0.2 are {par_num}, as the oscillator modes predict "
                f"({[p[1] for p in pred[:5]]}: 2 omega_1, omega_1 + omega_2, omega_1 + omega_3, 2 omega_2, omega_2 + omega_3)",
          par_num == [p[1] for p in pred[:5]], eq=eq_weak)
    check(NAME, f"cubic extrapolation of the lowest two excitations from g^2 = 0.2, 0.25, 0.3, 0.4 to g^2 = 0: {extrap[0]:.4f} and "
                f"{extrap[1]:.4f} against 2 omega_1 = {target[0]:.4f} and omega_1 + omega_2 = {target[1]:.4f} (tolerance 6e-3)",
          abs(extrap[0] - target[0]) < 6e-3 and abs(extrap[1] - target[1]) < 6e-3, eq=eq_weak)
    check(NAME, f"the third and fourth levels (even, nearly degenerate) extrapolate to {extrap[2]:.4f} and {extrap[3]:.4f}, sum "
                f"{extrap[2] + extrap[3]:.4f} against omega_1 + omega_3 + 2 omega_2 = {target[2] + target[3]:.4f} (tolerance 6e-3); "
                f"the fifth, odd, to {extrap[4]:.4f} against omega_2 + omega_3 = {target[4]:.4f} (tolerance 1e-2)",
          abs(extrap[2] + extrap[3] - target[2] - target[3]) < 6e-3 and abs(extrap[4] - target[4]) < 1e-2, eq=eq_weak)

    # =====================================================================================================
    # strong coupling
    # =====================================================================================================
    Mx, vac = strong_matrix()
    eq_strong = register(NAME, Equation(
        key="su2_strip_strong", entry=ENTRY, title="Strong-coupling splitting of the first excited level of the strip",
        source_kind="derived", source="derived here by exact degenerate second-order perturbation theory from su2_strip_basis",
        latex_str=(r"\Delta_i=\frac32\,g^{2}+\frac{\mu_i}{g^{6}}+O(g^{-14}),\ \ \mu_i=\text{eigenvalues of }M+2\cdot\mathbb 1,\ \ "
                   r"M=" + sp.latex(Mx) + r"\ \text{ on }(A,B,C)"),
        notes="The three strings (1/2) around plaquette a, b and c are degenerate at 3 g^2/2 and are coupled only at second order, "
              "through the vacuum and through states with two strings. The odd combination (A - C)/sqrt2 and the two even "
              "combinations separate by a <-> c parity. The vacuum shifts by -2/g^6 (three strings, each |W|^2 / (3/2) = 2/3), "
              "and the gap is the level minus the shifted vacuum. Exact algebraic numbers.",
    ))
    P3 = sp.Matrix([[0, 0, 1], [0, 1, 0], [1, 0, 0]])
    gapm = Mx - vac * sp.eye(3)
    even_b = sp.Matrix([[1 / sp.sqrt(2), 0], [0, 1], [1 / sp.sqrt(2), 0]])
    odd_b = sp.Matrix([1 / sp.sqrt(2), 0, -1 / sp.sqrt(2)])
    Meven = sp.simplify(even_b.T * gapm * even_b)
    modd = sp.simplify((odd_b.T * gapm * odd_b)[0])
    mu_even = sorted([sp.simplify(x) for x in Meven.eigenvals()], key=lambda x: float(x))
    mus = [(float(m_), 1) for m_ in mu_even] + [(float(modd), -1)]
    check(NAME, f"the second-order matrix on (A, B, C) is symmetric, commutes with the reflection (A <-> C), the vacuum shifts by "
                f"{vac} g^-6, and the coupled strings have gap coefficients mu = {[sp.nsimplify(m_) for m_ in mu_even]} (even) and "
                f"{sp.nsimplify(modd)} (odd), all in closed form",
          Mx == Mx.T and P3 * Mx * P3 == Mx and vac == -2 and sp.simplify((P3 * even_b - even_b)).is_zero_matrix, eq=eq_strong)
    obs = {}
    for g in (6.0, 8.0):
        lv, _ = excitations(g, 40, 4)
        obs[g] = lv
    ok, detail = True, []
    for g, lv in obs.items():
        pred_g = sorted([(1.5 * g + m_ / g ** 3, p) for m_, p in mus])
        for (e, p), (pe, pp) in zip(lv[:3], pred_g):
            ok &= abs(e - pe) < 5e-4 and p == pp
        detail.append(f"g^2={g}: " + ", ".join(f"{e:.5f}({p:+d})/{pe:.5f}" for (e, p), (pe, pp) in zip(lv[:3], pred_g)))
    check(NAME, "the three lowest excitations at g^2 = 6 and 8 have the predicted reflection parity and match 3g^2/2 + mu/g^6 to 5e-4: "
          + "; ".join(detail), ok, eq=eq_strong)
    lo = min(mus)
    check(NAME, f"the lowest strong-coupling excitation is the {'even' if lo[1] == 1 else 'odd'} combination, mu = {lo[0]:.5f} g^-6, "
                f"and numerically so at g^2 = 8 (parity {obs[8.0][0][1]:+d}); the one-plaquette string would sit at 3g^2/2 + 14/(15 g^6)",
          obs[8.0][0][1] == lo[1], eq=eq_strong)

    # =====================================================================================================
    # the gap against the ladder and one plaquette
    # =====================================================================================================
    grid = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
    rows = {}
    for g in grid:
        lv, _ = excitations(g, max(int(round(40 / g)), 40), 3)
        wl, _pl = ladder.levels(g, 3)
        rows[g] = (lv[0][0], lv[0][1], float(wl[1] - wl[0]), float(plaquette_gap(g)))
    eq_gap = register(NAME, Equation(
        key="su2_strip_gap", entry=ENTRY, title="The strip gap against the ladder and one-plaquette gaps",
        source_kind="derived", source="derived here from su2_strip_basis (numerical)",
        latex_str=(r"\Delta_{3}(g)\le\Delta_{2}(g)\le\Delta_{1}(g),\qquad "
                   r"\frac{\Delta_3}{\Delta_1}\to\sqrt{1-\tfrac{\sqrt2}{4}}\approx0.8040\ (g\to0),\quad\to1\ (g\to\infty)"),
        notes="Numerical on the grid g^2 = 0.25 .. 8 with truncation checked; Delta_1 is the one-plaquette gap and Delta_2 the "
              "ladder's. At weak coupling the harmonic limit gives 2.8284, 2.4495, 2.2741 for one, two and three plaquettes; "
              "at strong coupling the first excitation is a string around one plaquette whichever system it lives in. "
              "A volume dependence for three volumes only; nothing here is a statement about a larger or infinite system.",
    ))
    check(NAME, "strip gap <= ladder gap <= one-plaquette gap on the grid g^2 = 0.25 .. 8: " +
          ", ".join(f"{g}: {rows[g][0]:.4f} <= {rows[g][2]:.4f} <= {rows[g][3]:.4f}" for g in grid),
          all(rows[g][0] <= rows[g][2] + 1e-9 and rows[g][2] <= rows[g][3] + 1e-9 for g in grid), eq=eq_gap)
    r3 = {g: rows[g][0] / rows[g][3] for g in grid}
    limit3 = float(sp.sqrt(1 - sp.sqrt(2) / 4))
    check(NAME, f"the strip-to-plaquette ratio is {r3[0.25]:.4f} at g^2 = 0.25 (weak-coupling limit {limit3:.4f}) and {r3[8.0]:.5f} at "
                f"g^2 = 8 (strong-coupling limit 1)", abs(r3[0.25] - limit3) < 0.035 and r3[8.0] > 0.99
          and abs(limit3 - float(delta3 / (2 * sp.sqrt(2)))) < 1e-12, eq=eq_gap)
    gs2 = np.array([0.5, 0.7, 0.9, 1.1, 1.3, 1.6])
    gv = np.array([excitations(g, max(int(round(40 / g)), 40), 3)[0][0][0] for g in gs2])
    imin = int(np.argmin(gv))
    check(NAME, f"the strip gap has an interior minimum (scan g^2 = 0.5 .. 1.6: smallest {gv[imin]:.4f} at g^2 = {gs2[imin]}, larger at both "
                f"ends)", 0 < imin < len(gs2) - 1, eq=eq_gap)
    return rows


if __name__ == "__main__":
    tab = build_and_verify()
    print("g^2     strip gap  parity  ladder gap  plaquette gap")
    for g, (a, p, b, c) in tab.items():
        print(f"{g:6.3f}  {a:.6f}   {p:+d}    {b:.6f}    {c:.6f}")
