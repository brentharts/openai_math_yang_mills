# A pose-level model from the Chair44 atlas

This page describes [tools/chair44_model.py](tools/chair44_model.py), a first computational step on the idea in [HIERARCHIES.md](HIERARCHIES.md#7-an-actionable-idea-flickers-24-state-model): a lattice-style model with the 24 proper cubic rotations as its states. Citation keys are those of [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md#citation-key): [T] Tsiokos (arXiv:2609.19214), [HR] Hartshorn–Romanov (Research Square, unrefereed), [F] Flicker.

Run it with `make chair44` (about five seconds) or `make` for the full report. The section it adds to the PDF is `build/tex/chair44.tex`.

## What the model is, and what it is not

[F] proposes a 24-state cubic-lattice model with nearest-neighbour penalties for illegal contacts. **That penalty table is not reproduced here**, because it was not available to me. The script instead builds the model that [T]'s contact atlas supports on its own:

- a state is a rotation G, one of the 24, together with an integer translation t, acting on the unit chair (seven unit cubes a + [0,1]³, a ∈ {0,1}³ \ {(1,1,1)}) by x ↦ Gx + t;
- the energy is E = (overlapping cells) + (face-contacting chair pairs whose relative pose is not one of the 44 atlas poses);
- ground states are the configurations with E = 0.

Because chairs touch at offsets such as (2,2,2) and (4,0,0), this is a model on poses, not a nearest-neighbour model on unit steps of ℤ³. Whether it is equivalent to [F]'s model is an open question here.

## Data and the one inferred convention

The 44 poses (Figure 7 of [T]) and the eight child poses (Table 1) were **transcribed from a text extraction of the paper**. The extraction gives frames and offsets but not the action on the chair. The script uses x ↦ Gx + t about the origin on unit cubes, which I chose because it was the only convention I tried (unit or doubled cubes, rotation about the origin or the chair centre, G or Gᵀ) under which all 44 poses are overlap-free and in face contact with the root chair. That is an inference from fit, not a statement from the paper.

The legality test of [T] (pyramid profiles on a shared panel) and the panel labels of its Figure 3 were not in the extracted text, so the atlas is used as given and is not re-derived.

## What the script checks (each stops the build on failure)

| Check | Result |
|---|---|
| 44 distinct poses, all det G = +1 | passes |
| each of the 44 is overlap-free and face-touching the root | passes under the convention above |
| 24 proper cubic rotations; sympy: the group on the 4 body diagonals has order 24 | passes |
| the 8 child poses of Table 1 fill the doubled chair exactly | 56 cells, no overlap |
| iterated refinement (G,t) ↦ {(GH, 2t+Gu)}, levels 1–3 | 8, 64, 512 chairs; no overlap; patch fills the 2ⁿ-scaled chair |
| every face contact in those patches is an atlas contact | E = 0: 16, 183 and 1663 contacts, none illegal |

The exact fill of the doubled chair and E = 0 at three levels are strong tests of the transcription: a wrong digit in a child pose or an atlas offset would be very likely to break them.

## Findings that are reported, not asserted

- The atlas has **19 distinct rotations** among the 24 and is **closed under inversion** (all 44 poses have their inverse in the atlas).
- The substitution patches realise **30 of the 44 poses** (21 at level 1, 30 at levels 2 and 3), counting both orientations of each contact. [HR] state in their Corollary 1 that only 30 of the 44 contacts occur; the counts agree. I have only the secondary summary of [HR] for this.
- A level-2 patch already uses **all 24 rotation states**.
- **Zero-temperature random growth** (add a chair at an atlas pose of a placed chair if it overlaps nothing and creates no illegal contact; 24 runs, cap 300 chairs) **never got stuck**, and the first run used all 44 poses and all 24 states. A rejected candidate stays rejected as chairs are added, so an empty candidate pool would have been an exact dead end; none occurred.

The last point says the pairwise atlas rule, as I have used it, does not confine local growth to the 30 poses of the hierarchical tilings, at least up to 300 chairs. It does not show that the grown patches extend to infinite tilings, and it does not contradict [T]: the legality in [T] may involve more than the pairwise atlas, and the pairwise rule may be missing that. This is the point the model could be used to probe.

## Caveats

- No temperature, no thermodynamic limit and no Monte Carlo in the Metropolis sense; the "growth" is a zero-temperature greedy procedure.
- The atlas is assumed to be the complete set of legal contacts.
- Nothing here tests [F]'s face-colour or arrow rules, or says whether [F]'s 24-state model has Chair tilings as ground states.
- A growth run that reaches the cap shows only that random local growth was not blocked, not that a global tiling exists.

## Next steps

1. Obtain the Figure 3 panel labels (or the Lean/Python sources of [T]) and re-derive the 44 contacts from the pyramid geometry, as an independent check of the transcription.
2. Add a finite-temperature Metropolis on a finite patch, with a penalty per illegal contact, and look at whether the hierarchical structure appears as the temperature falls.
3. Check whether the 14 poses outside the 30 are excluded by contacts with two or more neighbours (a three-body condition), which would explain why growth does not stay hierarchical.

## Related pages

- [SPECTRE_CHAIR44.md](SPECTRE_CHAIR44.md), [MONOTILES.md](MONOTILES.md), [HIERARCHIES.md](HIERARCHIES.md), [PERIODIC_TILING.md](PERIODIC_TILING.md)
- [tools/README.md](tools/README.md), for how the scripts and the PDF are built.
