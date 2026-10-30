"""The Hilbert-Smith theorem (entry 304), recorded as a statement.

  [E304]  The Hilbert-Smith conjecture in every finite dimension
          preprints/The-Hilbert-Smith-conjecture-in-every-finite-dimension-September-23-2026/

SPECTRE_CLOCK_AND_CHAIR44.md Section 3 concludes that this theorem has no interaction with the Chair44
tiling hull: Theorem 1.2 restricts actions on connected manifolds, and a tiling hull is not a manifold. The
statements are kept here so the document and the graph can cite the manuscript and link the discussion. No
computation is attached, so no check is claimed.
"""
from eqlib import Equation, register

NAME = "hilbert_smith"
ENTRY = "304"
E304 = "[E304] Hilbert-Smith paper"


def build_and_verify():
    register(NAME, Equation(
        key="e304_p_adic_exclusion", entry=ENTRY,
        title="No faithful p-adic action on a manifold",
        source_kind="paper", source=f"{E304} Theorem 1.2",
        latex_str=(r"\mathbb Z_p=\varprojlim \mathbb Z/p^{k}\mathbb Z\ \curvearrowright\ M^{n}\ \text{(connected, Hausdorff, "
                   r"second countable, no boundary)}\ \Longrightarrow\ \text{the action has nonzero kernel}"),
        notes="Every jointly continuous action factors through a finite quotient. Theorem 1.1 then gives: a locally "
              "compact second-countable group acting faithfully on such a manifold is a Lie group. The fork's notes "
              "record that a tiling hull is not a manifold (it is locally a Euclidean factor times a totally "
              "disconnected one), so a 2-adic odometer on the Chair44 hull is consistent with the theorem.",
    ))
