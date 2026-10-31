"""sin^2(theta_W) = 3/8 at unification, run down to the Z mass with the library's one-loop law.

The knowledge graph links 'Hypercharge normalisation' (k, the ratio of traces, 5/3 over one generation) to
'Weak mixing at unification' (s_W = 1/(1 + k)) and that to 'One-loop running'. This script closes the loop at
one loop: given alpha_em^-1 and alpha_s at M_Z, require the three couplings to meet, with the library's relation
1/alpha_i(mu) = 1/alpha_G + (b_i / 2 pi) ln(M_G / mu), and solve for the weak mixing parameter at M_Z, the
unification scale and the unified coupling. The three equations are linear in (s_W, ln(M_G/M_Z), 1/alpha_G).
It does this for the Standard Model and for the minimal supersymmetric one. One loop, no thresholds, no
two-loop terms: it reproduces the classic contrast and nothing finer.
"""
import sys
from pathlib import Path

import sympy as sp

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import knowledgebase as kb
from explore import Result, main

ALPHA_EM_INV = sp.Rational(127952, 1000)      # alpha_em^-1 at M_Z (MSbar)
ALPHA_S = sp.Rational(1179, 10000)
S2_MEASURED = sp.Rational(23122, 100000)      # MSbar value at M_Z
MZ = sp.Rational(911876, 10000)               # GeV
B = {"SM": (sp.Rational(41, 10), sp.Rational(-19, 6), sp.Integer(-7)),
     "MSSM": (sp.Rational(33, 5), sp.Integer(1), sp.Integer(-3))}


def texexp(x):
    m, e = f"{x:.1e}".split("e")
    return rf"{m}\times10^{{{int(e)}}}"


def unify(b, k=sp.Rational(5, 3)):
    s2, t, aG = sp.symbols("s2 t aG")
    two_pi = 2 * sp.pi
    # alpha_1^-1 = (1/k) cos^2 alpha_em^-1 (k = 5/3 gives 3/5); alpha_2^-1 = s2 alpha_em^-1; alpha_3^-1 = 1/alpha_s
    eqs = [sp.Eq((1 / k) * (1 - s2) * ALPHA_EM_INV, aG + b[0] * t / two_pi),
           sp.Eq(s2 * ALPHA_EM_INV, aG + b[1] * t / two_pi),
           sp.Eq(1 / ALPHA_S, aG + b[2] * t / two_pi)]
    sol = sp.solve(eqs, [s2, t, aG], dict=True)[0]
    return sol[s2], sol[t], sol[aG]


def run() -> Result:
    r = Result("weak_mixing", "The weak mixing angle: from 3/8 at unification to the Z mass", prose="")
    # the library's laws
    mix = kb.equation("Weak mixing at unification")
    sW = kb.solve_for(mix, "weak mixing parameter")
    r.check("the library's s_W = 1/(1 + k) gives 3/8 at k = 5/3", sW.expr.subs(sp.Symbol("k"), sp.Rational(5, 3)) == sp.Rational(3, 8))
    run_eq = kb.equation("One-loop running").relation()
    b, aG_, ai, MG_, mu_ = (sp.Symbol(x) for x in ("b_i", "alpha_G", "alpha_i", "M_G", "mu"))
    aGv, biv, Lv = sp.symbols("aGv biv Lv", positive=True)
    ai_run = 1 / (1 / aGv + biv * Lv / (2 * sp.pi))                     # the library's relation solved for alpha_i
    resid = (run_eq.lhs - run_eq.rhs).subs({b: biv, aG_: aGv, ai: ai_run, MG_: sp.exp(Lv), mu_: 1})
    r.check("the form 1/alpha_i(mu) = 1/alpha_G + (b_i/2 pi) ln(M_G/mu) satisfies the library's 'One-loop running' identically",
            sp.simplify(resid) == 0)
    k = sp.Rational(5, 3)
    results = {}
    for name, bb in B.items():
        s2, t, aG = unify(bb, k)
        MG = MZ * sp.exp(t)
        results[name] = dict(s2=float(s2), MG=float(MG), aG_inv=float(aG), t=float(t))
    r.data.update({f"{n}: sin^2 theta_W(M_Z), M_G [GeV], 1/alpha_G": (round(v["s2"], 4), f"{v['MG']:.2e}", round(v["aG_inv"], 2))
                   for n, v in results.items()})
    r.data["measured sin^2 theta_W(M_Z)"] = float(S2_MEASURED)
    sm, ss = results["SM"], results["MSSM"]
    r.check("the three couplings meet in both cases by construction (three linear equations, three unknowns); the "
            "solution has positive M_G > M_Z and positive 1/alpha_G",
            all(v["t"] > 0 and v["aG_inv"] > 0 for v in results.values()))
    r.check("MSSM running predicts sin^2 theta_W(M_Z) within 0.01 of the measured value",
            abs(ss["s2"] - float(S2_MEASURED)) < 0.01)
    r.check("the Standard Model's prediction is further from the measured value than the MSSM's",
            abs(sm["s2"] - float(S2_MEASURED)) > abs(ss["s2"] - float(S2_MEASURED)))
    r.check("both predictions lie below the unification value 3/8, as running from 3/8 requires",
            sm["s2"] < 0.375 and ss["s2"] < 0.375)
    # consistency: run each solution back up to M_G and read sin^2 there from the three inverse couplings
    for name, bb in B.items():
        s2, t, aG = unify(bb, k)
        inv = [(1 / k) * (1 - s2) * ALPHA_EM_INV - bb[0] * t / (2 * sp.pi),
               s2 * ALPHA_EM_INV - bb[1] * t / (2 * sp.pi), 1 / ALPHA_S - bb[2] * t / (2 * sp.pi)]
        em_inv = inv[1] + k * inv[0]                                  # 1/alpha_em = 1/alpha_2 + k/alpha_1
        r.check(f"{name}: run to M_G, the three inverse couplings all equal 1/alpha_G and sin^2 theta_W = alpha_2^-1/alpha_em^-1 = 3/8",
                all(sp.simplify(x - aG) == 0 for x in inv) and sp.simplify(inv[1] / em_inv - sp.Rational(3, 8)) == 0)
    r.prose = (r"One-loop running with the library's law, from the measured $\alpha_{\mathrm{em}}^{-1}$ and $\alpha_s$ at $M_Z$ and the "
               r"requirement that the three couplings meet, predicts $\sin^2\theta_W(M_Z)="
               + f"{sm['s2']:.3f}" + r"$ in the Standard Model and $" + f"{ss['s2']:.3f}" + r"$ with the supersymmetric "
               r"coefficients, against the measured $0.2312$; the unification scales are $" + texexp(sm['MG']) + r"$ and $"
               + texexp(ss['MG']) + r"$ GeV. It is the classic one-loop contrast and nothing finer: thresholds, two-loop "
               r"terms and the supersymmetric scale are ignored, and the $3/8$ is the value at the meeting point, not an input.")
    r.display = r"s_W=\frac{1}{1+k}=\frac38\quad\xrightarrow{\ \text{one loop}\ }\quad s_W(M_Z)\approx %.3f\ (\text{SM}),\ %.3f\ (\text{MSSM})" % (sm["s2"], ss["s2"])
    return r


if __name__ == "__main__":
    main(run())
