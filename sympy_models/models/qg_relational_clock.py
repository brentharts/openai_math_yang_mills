"""A relational clock for the two-plaquette ladder: Page-Wootters time in its finite (Feynman-Kitaev) form.

Part II's Hamiltonian H_S (su2_lattice_ladder.py) generates evolution in an external parameter. A constraint theory has no
such parameter. This module adds a clock by hand and asks what "evolution" then means. The clock is a register with N
readings |k>, k = 0..N-1, the system is the ladder truncated at spin kmax, and the constraint is the propagation operator
    C = (1/2) sum_{k=0}^{N-2} [ (|k><k| + |k+1><k+1|) x 1 - |k+1><k| x U - |k><k+1| x U^dagger ],   U = exp(-i H_S tau).
C is Hermitian and positive semidefinite. Its kernel (the physical states) is spanned by the history states
    Psi(psi0) = N^(-1/2) sum_k |k> x U^k psi0 ,
one for every initial state psi0, so dim ker C = dim(system). Conditioning a physical state on the clock reading k returns the
Schroedinger state U^k psi0, and any system observable gives a relational expectation <A>(k) that equals the Heisenberg
value. In this sense time is a correlation inside a timeless state, which is what Page and Wootters proposed. The clock is a
discretisation: the usual continuous clock, with H_C = p on L^2(R), is its limit and is not built here.

Checked, in floating point on the truncated ladder (tolerances stated in the checks): C is Hermitian and positive
semidefinite; every history state is annihilated by C; the kernel of C has exactly the dimension of the system and coincides
with the span of the history states; the conditional state at reading k is U^k psi0; the relational expectation of the plaquette
operator W_a equals the Heisenberg expectation, and the clock readings are uniformly distributed; a periodic clock has,
for generic H_S, no physical state at all, which is why the open register is used; the constraints built from H_S at two
different couplings g have different kernels, so g labels a constraint and is not a time; and the truncation at spin kmax
changes the short-time <W_a> by less than a stated number.

Not tested here: a continuous clock; a clock that interacts with the system (the conditioning is then not a unitary evolution,
Smith and Ahmadi); gravity: the Hamiltonian is the lattice gauge-theory Hamiltonian, and the constraint is linear in it, which is
not the Hamiltonian constraint of loop quantum gravity, whose algebra is quadratic and spatially dependent;
diffeomorphism invariance; the Master Constraint; anything about the spin-network Hilbert space beyond the ladder. This
shows that a gauge-theory Hamiltonian can be given a relational time when a clock is added, and no more.
"""
import sys
from pathlib import Path

import numpy as np
import scipy.linalg as sla

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from eqlib import Equation, check, register
from models.su2_lattice_ladder import build, operators

NAME = "qg_relational_clock"
ENTRY = "spin-network-qg"

KMAX = 4          # spin truncation (k = 2j)
G2 = 1.0          # g^2
TAU = 0.3         # clock step
NCLK = 10         # clock readings


def system(g2=G2, kmax=KMAX):
    H, states, idx = build(kmax, g2)
    Wa, Wb, _, _, _ = operators(kmax)
    return H.toarray(), Wa.toarray(), Wb.toarray(), states


def unitary(H, tau=TAU):
    return sla.expm(-1j * tau * H)


def constraint(U, n=NCLK, periodic=False):
    d = U.shape[0]
    C = np.zeros((n * d, n * d), complex)
    I = np.eye(d)
    steps = n if periodic else n - 1
    for k in range(steps):
        k1 = (k + 1) % n
        C[k * d:(k + 1) * d, k * d:(k + 1) * d] += 0.5 * I
        C[k1 * d:(k1 + 1) * d, k1 * d:(k1 + 1) * d] += 0.5 * I
        C[k1 * d:(k1 + 1) * d, k * d:(k + 1) * d] -= 0.5 * U
        C[k * d:(k + 1) * d, k1 * d:(k1 + 1) * d] -= 0.5 * U.conj().T
    return C


def history(U, psi0, n=NCLK):
    out, v = [], psi0.astype(complex)
    for _ in range(n):
        out.append(v)
        v = U @ v
    return np.concatenate(out) / np.sqrt(n)


def build_and_verify():
    eq = register(NAME, Equation(
        key="qg_relational_clock", entry=ENTRY, title="A relational clock for the ladder: the history states of the propagation constraint",
        source_kind="derived", source="derived here; the discrete clock constraint is the propagation term of the Feynman-Kitaev construction, and conditioning on a clock is the Page-Wootters idea",
        latex_str=(r"C=\tfrac12\sum_{k=0}^{N-2}\Bigl[(|k\rangle\langle k|+|k{+}1\rangle\langle k{+}1|)\otimes\mathbb 1-|k{+}1\rangle\langle k|\otimes U-|k\rangle\langle k{+}1|\otimes U^{\dagger}\Bigr],\quad U=e^{-iH_S\tau}"
                   r"%%BR%%\ker C=\Bigl\{\tfrac1{\sqrt N}\sum_{k}|k\rangle\otimes U^{k}\psi_0\Bigr\},\qquad \dim\ker C=\dim\mathcal H_S"
                   r"%%BR%%\langle k|\otimes\mathbb 1\;\Psi=\tfrac1{\sqrt N}\,U^{k}\psi_0,\qquad \langle A\rangle(k)=\langle\psi_0|U^{\dagger k}AU^{k}|\psi_0\rangle"),
        notes="The system is the two-plaquette ladder of su2_lattice_ladder.py truncated at spin kmax; the clock is an open register of N readings. "
              "The constraint is linear in the lattice Hamiltonian, so this is a model of relational time for a gauge theory with a clock added by hand. "
              "It is not the Hamiltonian constraint of loop quantum gravity, it does not involve diffeomorphisms, and the clock is a discretisation of the "
              "continuous Page-Wootters clock, not that clock.",
    ))
    H, Wa, Wb, states = system()
    d = H.shape[0]
    U = unitary(H)
    C = constraint(U)
    rng = np.random.default_rng(7)

    herm = np.max(np.abs(C - C.conj().T))
    w = np.linalg.eigvalsh((C + C.conj().T) / 2)
    check(NAME, f"the ladder truncated at spin {KMAX / 2:g} at g^2 = {G2:g} has {d} states; its Hamiltonian is Hermitian and the propagation constraint "
                f"C with N = {NCLK} readings is Hermitian (largest asymmetry {herm:.1e}) and positive semidefinite (smallest eigenvalue {w.min():.1e})",
          np.max(np.abs(H - H.T)) < 1e-12 and herm < 1e-12 and w.min() > -1e-12, eq=eq)

    psis = [rng.normal(size=d) + 1j * rng.normal(size=d) for _ in range(5)]
    res = max(np.linalg.norm(C @ history(U, p / np.linalg.norm(p))) for p in psis)
    check(NAME, f"every history state is a physical state: ||C Psi|| <= {res:.1e} for five random initial states", res < 1e-11, eq=eq)

    ker = int(np.sum(w < 1e-9))
    gap = np.sort(w)[ker]
    check(NAME, f"the kernel of C has dimension {ker} = dim of the system {d}; the next eigenvalue is {gap:.3f}", ker == d and gap > 1e-3, eq=eq)

    basis = np.stack([history(U, e) for e in np.eye(d)], axis=1)
    q, _ = np.linalg.qr(basis)
    evals, evecs = np.linalg.eigh((C + C.conj().T) / 2)
    kerspace = evecs[:, :ker]
    dist = np.linalg.norm(q @ q.conj().T - kerspace @ kerspace.conj().T)
    check(NAME, f"the span of the history states is exactly the kernel of C (projector difference {dist:.1e})", dist < 1e-9, eq=eq)

    p0 = psis[0] / np.linalg.norm(psis[0])
    Psi = history(U, p0)
    fid = []
    for k in range(NCLK):
        cond = Psi[k * d:(k + 1) * d] * np.sqrt(NCLK)
        ref = np.linalg.matrix_power(U, k) @ p0
        fid.append(abs(np.vdot(ref, cond)))
    check(NAME, f"conditioning on the clock reading k gives the Schroedinger state U^k psi0 (smallest fidelity over k = 0..{NCLK - 1}: {min(fid):.12f})",
          min(fid) > 1 - 1e-12, eq=eq)

    # relational expectation of a plaquette operator against the Heisenberg one, from an initial state with weight on several levels
    p1 = np.zeros(d, complex)
    p1[0] = 1.0                   # the strong-coupling vacuum |0,0,0>
    p1[min(3, d - 1)] += 0.5
    p1 /= np.linalg.norm(p1)
    Ps = history(U, p1)
    rel, heis, prob = [], [], []
    for k in range(NCLK):
        blk = Ps[k * d:(k + 1) * d]
        pk = np.vdot(blk, blk).real
        prob.append(pk)
        rel.append((np.vdot(blk, Wa @ blk).real) / pk)
        v = np.linalg.matrix_power(U, k) @ p1
        heis.append(np.vdot(v, Wa @ v).real)
    check(NAME, f"the relational expectation of the plaquette operator W_a given the clock reading equals the Heisenberg expectation for all {NCLK} readings "
                f"(largest difference {max(abs(a - b) for a, b in zip(rel, heis)):.1e}); the clock readings are uniform (largest deviation from 1/N: "
                f"{max(abs(p - 1 / NCLK) for p in prob):.1e}); <W_a> moves from {heis[0]:.4f} to {max(heis, key=lambda x: abs(x - heis[0])):.4f}, so the test is not vacuous",
          max(abs(a - b) for a, b in zip(rel, heis)) < 1e-11 and max(abs(p - 1 / NCLK) for p in prob) < 1e-12
          and max(abs(h - heis[0]) for h in heis) > 1e-3, eq=eq)

    # a periodic clock
    Cp = constraint(U, periodic=True)
    wp = np.linalg.eigvalsh((Cp + Cp.conj().T) / 2)
    kerp = int(np.sum(wp < 1e-9))
    Eh = np.linalg.eigvalsh(H)
    phases = np.exp(-1j * NCLK * TAU * Eh)
    expected = int(np.sum(np.abs(phases - 1) < 1e-9))
    check(NAME, f"a periodic register has physical states only for the eigenvalues of H_S with exp(-i E N tau) = 1: its kernel has dimension {kerp}, "
                f"the number of such eigenvalues is {expected}, and the open register has {d}", kerp == expected and kerp < d, eq=eq)

    # g is not a time
    H2, _, _, _ = system(g2=1.1)
    C2 = constraint(unitary(H2))
    w2, v2 = np.linalg.eigh((C2 + C2.conj().T) / 2)
    k2 = v2[:, :int(np.sum(w2 < 1e-9))]
    P1 = kerspace @ kerspace.conj().T
    P2 = k2 @ k2.conj().T
    ov = np.linalg.norm(P1 @ P2) ** 2 / d
    check(NAME, f"the coupling labels a constraint, not a time: the kernels of C at g^2 = 1.0 and g^2 = 1.1 are different subspaces (normalised overlap {ov:.3f} < 1)",
          ov < 0.999, eq=eq)

    # truncation
    Hb, Wab, _, _ = system(kmax=KMAX + 2)
    db = Hb.shape[0]
    v0 = np.zeros(db, complex)
    v0[0] = 1.0
    v0[min(3, db - 1)] += 0.5
    v0 /= np.linalg.norm(v0)
    Ub = unitary(Hb)
    big = []
    for k in range(4):
        v = np.linalg.matrix_power(Ub, k) @ v0
        big.append(np.vdot(v, Wab @ v).real)
    small = [heis[k] for k in range(4)]
    diff = max(abs(a - b) for a, b in zip(big, small))
    check(NAME, f"truncation: raising the cut from spin {KMAX / 2:g} to {(KMAX + 2) / 2:g} changes <W_a> at the first four readings by at most {diff:.1e}",
          diff < 5e-3, eq=eq)
    return eq


if __name__ == "__main__":
    build_and_verify()
    print("ok")
