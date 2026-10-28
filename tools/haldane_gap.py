#!/usr/bin/env python3
"""Exact diagonalisation of spin chains: the spin-1 Haldane gap (CONTENTS entry 268).

What the repo states (entry 268): for the pure antiferromagnetic Heisenberg spin-1 chain on even
periodic rings the spectral gap stays uniformly positive as the chain grows (the "Haldane gap
conjecture"); a companion paper treats odd open chains with an endpoint field.

Model computed here
-------------------
    H(theta) = sum_i [ cos(theta) S_i . S_{i+1} + sin(theta) (S_i . S_{i+1})^2 ],  periodic ring of L sites.

theta = 0 is the Heisenberg chain of the repo paper. theta = atan(1/3) is the AKLT point, whose ground
state is known exactly (E_0 = -2L/3), which gives a strong test of the code. The Haldane phase
extends over -pi/4 < theta < pi/4.

The gap is  Delta(L) = E_0(S^z = 1) - E_0(S^z = 0),  from Lanczos in S^z sectors. The Hamiltonian is
built from an arbitrary two-site operator, so one code path serves spin 1/2 (used for validation),
the bilinear term and the biquadratic term.

sympy is used to diagonalise the two-site bond operator, and to show symbolically that the AKLT
Hamiltonian is a sum of projectors onto total bond spin 2.

Finite rings give a sequence of numbers; this is numerical evidence about Delta(L), not a proof of a
uniform positive lower bound, which is what the repo paper establishes.
"""

from __future__ import annotations

import sys

import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla
import sympy as sp

import common as C

LITERATURE_GAP = 0.4105  # infinite-chain Haldane gap (DMRG: White & Huse 1993; Golinelli et al.), approximate


# ----------------------------------------------------------------------------
# Local operators and two-site bond operators
# ----------------------------------------------------------------------------
def spin_ops(spin2: int):
    """Spin matrices for spin S = spin2/2, basis ordered m = -S, ..., +S. Returns Sx, Sy, Sz."""
    S = spin2 / 2.0
    d = spin2 + 1
    m = np.arange(-S, S + 1)
    Sz = np.diag(m)
    Sp = np.zeros((d, d))
    for k in range(d - 1):
        Sp[k + 1, k] = np.sqrt(S * (S + 1) - m[k] * (m[k] + 1))
    Sm = Sp.T
    return (Sp + Sm) / 2, (Sp - Sm) / 2j, Sz


def bond_ops(spin2: int):
    Sx, Sy, Sz = spin_ops(spin2)
    SS = (np.kron(Sx, Sx) + np.kron(Sy, Sy) + np.kron(Sz, Sz)).real
    return SS, SS @ SS


# ----------------------------------------------------------------------------
# Sector Hamiltonians
# ----------------------------------------------------------------------------
class Sector:
    """Basis of S^z = sz sector of an L-site ring of local dimension d = spin2 + 1."""

    def __init__(self, L: int, spin2: int, sz2: int):
        self.L, self.d = L, spin2 + 1
        codes = np.arange(self.d**L, dtype=np.int64)
        digits = np.empty((self.d**L, L), dtype=np.int8)
        tmp = codes.copy()
        for i in range(L):
            digits[:, i] = tmp % self.d
            tmp //= self.d
        # S^z (times 2) = sum_i (2 m_i), with 2m = 2*digit - spin2
        total2 = (2 * digits.astype(np.int64) - spin2).sum(axis=1)
        keep = total2 == sz2
        self.codes = codes[keep]
        self.digits = digits[keep]
        self.pow = self.d ** np.arange(L, dtype=np.int64)

    def bond_matrix(self, op: np.ndarray) -> sps.csr_matrix:
        """Sum over ring bonds of a two-site operator (indices a*d+b: a on site i, b on site i+1)."""
        L, d, n = self.L, self.d, len(self.codes)
        rows, cols, vals = [], [], []
        for i in range(L):
            j = (i + 1) % L
            for a in range(d):
                for b in range(d):
                    idx_in = np.nonzero((self.digits[:, i] == a) & (self.digits[:, j] == b))[0]
                    if idx_in.size == 0:
                        continue
                    for a2 in range(d):
                        for b2 in range(d):
                            amp = op[a2 * d + b2, a * d + b]
                            if abs(amp) < 1e-14:
                                continue
                            new = self.codes[idx_in] + (a2 - a) * self.pow[i] + (b2 - b) * self.pow[j]
                            idx_out = np.searchsorted(self.codes, new)
                            assert np.all(self.codes[idx_out] == new), "bond operator leaves the sector"
                            rows.append(idx_out)
                            cols.append(idx_in)
                            vals.append(np.full(idx_in.size, amp))
        return sps.csr_matrix(
            (np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))


def lowest(H: sps.csr_matrix, k: int = 1) -> np.ndarray:
    n = H.shape[0]
    if n <= 1500:
        return np.linalg.eigvalsh(H.toarray())[:k]
    vals = spla.eigsh(H, k=k, which="SA", tol=1e-10, return_eigenvectors=False)
    return np.sort(vals)


def ring_energies(L: int, spin2: int, thetas: list[float]):
    """For each theta return (E0 in Sz=0, E0 in Sz=1) of H(theta); spin2 must be even for Sz sectors 0, 1."""
    SS, SS2 = bond_ops(spin2)
    out = []
    mats = {}
    for sz2 in (0, 2):
        sec = Sector(L, spin2, sz2)
        mats[sz2] = (sec.bond_matrix(SS), sec.bond_matrix(SS2))
    for th in thetas:
        e = []
        for sz2 in (0, 2):
            H1, H2 = mats[sz2]
            e.append(float(lowest(np.cos(th) * H1 + np.sin(th) * H2)[0]))
        out.append(tuple(e))
    return out


def spin_half_energy(L: int) -> float:
    """Heisenberg spin-1/2 ring ground-state energy (validation of the engine)."""
    SS, _ = bond_ops(1)
    sec = Sector(L, 1, 0)
    return float(lowest(sec.bond_matrix(SS))[0])


# ----------------------------------------------------------------------------
# Symbolic part
# ----------------------------------------------------------------------------
def symbolic():
    s2 = sp.sqrt(2)
    Sz = sp.diag(-1, 0, 1)
    Sp = sp.Matrix([[0, 0, 0], [s2, 0, 0], [0, s2, 0]])
    Sm = Sp.T
    Sx = (Sp + Sm) / 2
    Sy = (Sp - Sm) / (2 * sp.I)
    kp = sp.kronecker_product
    SS = sp.simplify(kp(Sx, Sx) + kp(Sy, Sy) + kp(Sz, Sz))
    eig = SS.eigenvals()
    C.check(eig == {-2: 1, -1: 3, 1: 5}, "S.S on two spin-1 sites has eigenvalues -2 (x1), -1 (x3), +1 (x5)")
    I9 = sp.eye(9)
    P2 = (SS + 2 * I9) * (SS + I9) / 6  # projector onto the eigenvalue +1 space (total spin J = 2)
    C.check(sp.simplify(P2 * P2 - P2) == sp.zeros(9, 9) and P2.trace() == 5,
            "P_2 = (S.S+2)(S.S+1)/6 is a rank-5 projector (total bond spin 2)")
    aklt = SS + SS * SS / 3
    C.check(sp.simplify(aklt - (2 * P2 - sp.Rational(2, 3) * I9)) == sp.zeros(9, 9),
            "S.S + (S.S)^2/3 = 2 P_2 - 2/3")
    return eig


# ----------------------------------------------------------------------------
def main(full: bool = False) -> None:
    import matplotlib.pyplot as plt

    C.ensure_dirs()
    C.setup_style()
    eig = symbolic()

    # --- validation of the engine ---------------------------------------------------------------
    C.check(abs(spin_half_energy(4) + 2.0) < 1e-9, "spin-1/2 ring L=4: E0 = -2 (exact)")
    C.check(abs(spin_half_energy(6) + 2.8027756377) < 1e-6, "spin-1/2 ring L=6: E0 = -2.80278 (known)")
    th_aklt = float(np.arctan(1 / 3))
    # H(theta) = cos(theta) * [S.S + tan(theta) (S.S)^2], so at tan(theta) = 1/3 it is cos(theta) times the
    # AKLT Hamiltonian, whose exact ring ground energy is -2L/3. cos(atan(1/3)) = 3/sqrt(10).
    for L in (6, 8):
        e0 = ring_energies(L, 2, [th_aklt])[0][0]
        exact = -2.0 * L / np.sqrt(10.0)
        C.check(abs(e0 - exact) < 1e-8, f"AKLT point, ring L={L}: E0 = -2L/sqrt(10) = {exact:.6f} (exact)")

    # --- Heisenberg gap vs L ----------------------------------------------------------------------
    Ls = [4, 6, 8, 10, 12] + ([14] if full else [])
    gaps, e0s = [], []
    for L in Ls:
        e_s0, e_s1 = ring_energies(L, 2, [0.0])[0]
        gaps.append(e_s1 - e_s0)
        e0s.append(e_s0 / L)
        print(f"L={L:2d}: E0/L={e_s0 / L:.8f}  gap={e_s1 - e_s0:.6f}", flush=True)
    gaps = np.array(gaps)
    C.check(abs(e0s[Ls.index(12)] + 1.401484) < 0.01, "ground-state energy per site at L=12 within 1% of infinite chain (-1.40148)")
    C.check(bool(np.all(np.diff(gaps) < 0)), "periodic-ring gap decreases monotonically with L")
    C.check(bool(np.all(gaps > LITERATURE_GAP)), "ring gaps stay above the infinite-chain gap (approach from above)")
    dd = np.abs(np.diff(gaps))
    C.check(dd[-1] < dd[0], "finite-size corrections to the gap shrink with L")
    # Finite-size form for a massive excitation on a ring: Delta(L) = Delta_inf + A exp(-L/xi)/sqrt(L).
    # Three parameters, fitted to the three largest sizes (an exact solve, so it carries no error bar).
    from scipy.optimize import curve_fit

    def model(L_, D, A, xi):
        return D + A * np.exp(-L_ / xi) / np.sqrt(L_)

    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # "covariance could not be estimated": expected for an exact fit
        pfit, _ = curve_fit(model, np.array(Ls[-3:], float), gaps[-3:], p0=[0.41, 2.0, 6.0], maxfev=20000)
    extrap, A_fit, xi_fit = (float(v) for v in pfit)
    C.check(abs(extrap / LITERATURE_GAP - 1) < 0.03,
            f"finite-size fit gives Delta_inf = {extrap:.4f}, within 3% of the literature value ~{LITERATURE_GAP} "
            f"(correlation length from fit: {xi_fit:.2f})")
    aitken = extrap  # keep the downstream variable name used in figures and tables

    # --- gap as function of theta -----------------------------------------------------------------------
    thetas = [-0.8, -0.6, -0.4, -0.2, 0.0, 0.1, 0.2, th_aklt, 0.45, 0.6, 0.7, 0.78]
    scan = {}
    for L in (8, 10):
        res = ring_energies(L, 2, thetas)
        scan[L] = [e1 - e0 for e0, e1 in res]
        print(f"theta scan L={L} done", flush=True)
    aklt_gap = scan[10][thetas.index(th_aklt)]

    # --- figure ---------------------------------------------------------------------------------------------
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.2))
    ax[0].plot(Ls, gaps, "o-", label=r"$\Delta(L)$, periodic ring")
    ax[0].axhline(LITERATURE_GAP, color="k", ls="--", lw=1, label=rf"literature $\approx{LITERATURE_GAP}$")
    ax[0].axhline(aitken, color="C3", ls=":", lw=1.5, label=rf"fit $\Delta_\infty+Ae^{{-L/\xi}}/\sqrt{{L}}$: {aitken:.4f}")
    ax[0].set_xlabel("$L$")
    ax[0].set_ylabel(r"$\Delta=E_0(S^z{=}1)-E_0(S^z{=}0)$")
    ax[0].set_title(r"Spin-1 Heisenberg ring, $\theta=0$")
    ax[0].set_ylim(0.3, 0.9)
    ax[0].legend(fontsize=8)

    for L, st in ((8, "o--"), (10, "s-")):
        ax[1].plot(thetas, scan[L], st, ms=4, label=f"$L={L}$")
    ax[1].axvspan(-np.pi / 4, np.pi / 4, color="C2", alpha=0.12, label=r"Haldane phase $|\theta|<\pi/4$")
    ax[1].axvline(th_aklt, color="C1", ls=":", lw=1)
    ax[1].text(th_aklt + 0.02, 0.05, "AKLT", fontsize=8, color="C1")
    ax[1].axvline(0, color="0.5", ls=":", lw=1)
    ax[1].text(0.02, 0.02, "Heisenberg", fontsize=8, color="0.4")
    ax[1].set_xlabel(r"$\theta$")
    ax[1].set_ylabel(r"$\Delta$")
    ax[1].set_title("Gap across the bilinear-biquadratic family")
    ax[1].text(0.2, 0.92, "$\\Delta=0$ at right edge: lowest level is a\nspin-1 multiplet (finite-size crossing)",
               fontsize=8, color="0.3")
    ax[1].legend(fontsize=8, loc="lower left")
    C.savefig(fig, "haldane_gap")

    # --- data and LaTeX ---------------------------------------------------------------------------------------
    C.save_json("haldane_gap", {"L": Ls, "gap": gaps.tolist(), "e0_per_site": e0s, "aitken": float(aitken),
                                "aklt_gap_L10": float(aklt_gap), "fit_A": A_fit, "fit_xi": xi_fit,
                                "theta_scan": {str(L): {"theta": thetas, "gap": scan[L]} for L in scan}})
    rows = [[f"{L}", f"{e:.6f}", f"{g:.5f}"] for L, e, g in zip(Ls, e0s, gaps)]
    eigtex = ",\\ ".join(f"{sp.latex(k)}\\ (\\times{v})" for k, v in sorted(eig.items()))
    stats = "including $L=14$" if full else "up to $L=12$"

    tex = r"""
\section{The spin-1 Haldane gap}
\label{sec:haldane}

\paragraph{Source in the repository.}
CONTENTS entry 268; preprints \texttt{The-periodic-spin-one-Haldane-gap-September-24-2026} and
\texttt{A-boundary-field-gap-for-the-spin-one-Heisenberg-chain-September-24-2026}. For the pure antiferromagnetic Heisenberg chain on even periodic
rings the first paper proves that the spectral gap remains uniformly positive as the chain grows.

\paragraph{Model computed here.} A periodic ring of $L$ spin-1 sites with
""" + C.equation(
        r"H(\theta)=\sum_i\Big[\cos\theta\,\mathbf S_i\!\cdot\!\mathbf S_{i+1}+\sin\theta\,(\mathbf S_i\!\cdot\!\mathbf S_{i+1})^2\Big],", "eq:H"
    ) + r"""
where $\theta=0$ is the Heisenberg chain and $\theta=\arctan\frac13$ is the AKLT point. The gap is
$\Delta(L)=E_0(S^z{=}1)-E_0(S^z{=}0)$, computed by Lanczos in $S^z$ sectors (""" + stats + r""").

\paragraph{Symbolic bond algebra (sympy).} On two spin-1 sites, $\mathbf S_1\!\cdot\!\mathbf S_2=\frac12[J(J+1)-4]$ has eigenvalues
""" + C.equation(eigtex + r"\quad\text{for }J=0,1,2.", "eq:bond") + r"""
Hence $P_2=\frac16(\mathbf S\!\cdot\!\mathbf S+2)(\mathbf S\!\cdot\!\mathbf S+1)$ is the rank-5 projector onto $J=2$, and the identity
""" + C.equation(r"\mathbf S\!\cdot\!\mathbf S+\tfrac13(\mathbf S\!\cdot\!\mathbf S)^2=2P_2-\tfrac23", "eq:aklt") + r"""
(checked symbolically) shows that $H_{\rm AKLT}=\sum_i[\mathbf S_i\!\cdot\!\mathbf S_{i+1}+\frac13(\mathbf S_i\!\cdot\!\mathbf S_{i+1})^2]$
is a sum of projectors up to a constant, so $E_0=-\frac{2L}{3}$ exactly on a ring. Since $H(\theta)=\cos\theta\,H_{\rm AKLT}$ at $\tan\theta=\frac13$,
the code must give $E_0=-\frac{2L}{3}\cos\theta=-\frac{2L}{\sqrt{10}}$ there, and it does to $10^{-8}$ for $L=6,8$. It also reproduces $E_0=-2$ and $-2.80278$
for the spin-$\frac12$ rings with $L=4,6$. Together these validate the bond-operator machinery used for the Heisenberg term.

""" + C.figure(
        "haldane_gap",
        r"Left: gap of the periodic spin-1 Heisenberg ring against $L$, with the infinite-chain literature value and a finite-size fit. "
        r"Right: $E_0(S^z{=}1)-E_0(S^z{=}0)$ across the bilinear--biquadratic family for two ring sizes, with the Haldane phase shaded. "
        r"For $\theta\gtrsim0.65$ these small rings reach $\Delta=0$ before $\pi/4$: the lowest level becomes a spin-1 multiplet, a finite-size "
        r"level crossing, and the plotted quantity no longer measures a singlet--triplet gap.",
        "\\linewidth",
    ) + C.table(
        [r"$L$", r"$E_0/L$", r"$\Delta(L)$"], rows,
        caption=(r"Heisenberg ring ($\theta=0$). The infinite-chain DMRG gap is approximately %.4f (White and Huse, 1993); "
                 r"a finite-size fit to the three largest rings gives %.4f." % (LITERATURE_GAP, aitken)),
        colspec="rrr",
    ) + r"""
\paragraph{Caveat.} Exact diagonalisation gives $\Delta(L)$ for a handful of sizes; the repository paper proves a uniform positive lower bound for
all even $L$. The extrapolation fits an assumed form $\Delta_\infty+Ae^{-L/\xi}/\sqrt L$ exactly through three points, and the quoted literature value is from DMRG, not from this repository.
In the $\theta$ scan, the gap vanishes for $\theta\gtrsim0.65$ on the $L=8,10$ rings because the lowest level there is a spin-1 multiplet
(the minima of the $S^z=0$ and $S^z=1$ sectors coincide); this is a finite-size level crossing near the critical region, so the scan is only
meaningful as a gap in the interior of the Haldane window. Outside $|\theta|<\pi/4$ the quantity $E_0(S^z{=}1)-E_0(S^z{=}0)$ is not the physical gap
of the other phases.
"""
    C.write_tex("haldane", tex)
    print("wrote build/tex/haldane.tex")


if __name__ == "__main__":
    main(full="--full" in sys.argv)
