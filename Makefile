# Physics models for the OpenAI math release (fork tooling).
#
# Variables: PY=python3   FULL=1 (same as `make full`)

PY       ?= python3
FULL     ?=
ARGS     := $(if $(FULL),--full,)
BUILD    := build
STAMPS   := $(BUILD)/stamps
PDF      := $(BUILD)/physics_models.pdf
TEMPLATE := tools/templates/physics_models.tex

SCRIPTS     := bloch_law on_model_mc kt_flow haldane_gap chair44_model
STAMP_FILES := $(addprefix $(STAMPS)/,$(addsuffix .stamp,$(SCRIPTS)))
COMMON      := tools/common.py

.DEFAULT_GOAL := all
.PHONY: all help setup models pdf full bloch on kt haldane chair44 figures combinatorics zeta sympy sympy-pdf hopf su2 sym2lean sym2lean-pdf sym2lean-test dsl pyast physics-lean clean

all: pdf

help:
	@echo "make            run the models, then build $(PDF)"
	@echo "make setup      install Python dependencies (tools/requirements.txt)"
	@echo "make models     run all model scripts (figures, LaTeX fragments, data in $(BUILD)/)"
	@echo "make pdf        build the PDF (runs any models that are out of date)"
	@echo "make full       longer Monte Carlo runs and one more ring size (several minutes)"
	@echo "make bloch | on | kt | haldane | chair44    run a single model"
	@echo "make sympy      run the sympy models in sympy_models/ and export their equations (build in sympy_models/build/)"
	@echo "make hopf       check only the Hopf models (fibration, insulator, BPST); see HOPF_GAUGE.md"
	@echo "make su2        check only the SU(2) lattice-gap models (2D string tension, one-plaquette gap, two-plaquette ladder); see SU2_LATTICE_GAP.md"
	@echo "make sympy-pdf  also build sympy_models/build/equations_document.pdf, grouped by CONTENTS entry"
	@echo "make figures     regenerate the committed figures/*.png used by SPECTRE_CLOCK_AND_CHAIR44.md (about 90 s)"
	@echo "make combinatorics  exact checks and figures for COMBINATORICS_REPRESENTATIONS.md (about 1 min; needs networkx)"
	@echo "make zeta           zero-counting check and strip figures for ZETA_ZEROS.md (about 30 s; needs mpmath)"
	@echo "make sym2lean       sympy -> LaTeX -> Lean 4: write $(BUILD)/lean/*.lean from the sympy models (SYM2LEAN_ARGS='bloch_law --check --prove')"
	@echo "make sym2lean-pdf   also typeset the report $(BUILD)/sym2lean_document.pdf"
	@echo "make sym2lean-test  selftests of the three sym2lean modules (no models run)"
	@echo "make dsl         @lean decorator DSL: tools/dsl_models/*.py -> $(BUILD)/lean/Dsl_mathlib.lean"
	@echo "make pyast       Python (ast) -> LaTeX + Lean with loop-invariant VCs: tools/py_models/*.py -> $(BUILD)/lean/PyAst_mathlib.lean"
	@echo "make physics-lean  PHYSICS.md families, their papers and Lean arithmetic statements, one xelatex document: $(BUILD)/oai/physics.pdf (several minutes; GUIDE=PHYSICS.md)"
	@echo "make clean      remove $(BUILD)/"
	@echo "Variables: PY=python3  FULL=1  SYM2LEAN_ARGS="

setup:
	$(PY) -m pip install -r tools/requirements.txt

# --- models ---------------------------------------------------------------------------------
# One stamp per script, so editing a script (or common.py) re-runs only what depends on it.
$(STAMPS)/%.stamp: tools/%.py $(COMMON)
	@mkdir -p $(STAMPS)
	$(PY) tools/$*.py $(ARGS)
	@touch $@

models: $(STAMP_FILES)

bloch:   $(STAMPS)/bloch_law.stamp
on:      $(STAMPS)/on_model_mc.stamp
kt:      $(STAMPS)/kt_flow.stamp
haldane: $(STAMPS)/haldane_gap.stamp
chair44: $(STAMPS)/chair44_model.stamp

# Figures for the markdown page SPECTRE_CLOCK_AND_CHAIR44.md; written to figures/ (committed), not part of the PDF.
figures:
	$(PY) tools/hierarchy_figures.py $(ARGS)

# Exact checks and figures for COMBINATORICS_REPRESENTATIONS.md; writes figures/comb_*.png, not part of the PDF.
combinatorics:
	$(PY) tools/combinatorics_checks.py all --figures --nmax=22

# Zero-counting check and figures for ZETA_ZEROS.md; writes figures/zeta_*.png, not part of the PDF.
zeta:
	$(PY) tools/zeta_strip_figure.py --nzeros=200 --figures

# sympy models: equations checked symbolically and exported as LaTeX (see sympy_models/README.md).
sympy:
	$(PY) sympy_models/export.py

sympy-pdf:
	$(PY) sympy_models/export.py --pdf

# sympy -> LaTeX -> Lean 4 (tools/sym2lean.py): every equation is printed, read back, compared with its sympy
# expression, and only then written as a Lean theorem. `--check` gives the files to Lean if there is one.
SYM2LEAN_ARGS ?=
sym2lean:
	$(PY) tools/sym2lean.py $(SYM2LEAN_ARGS)

sym2lean-pdf: sym2lean
	@command -v pdflatex >/dev/null || { echo "pdflatex not found (install TeX Live)"; exit 1; }
	cp tools/templates/sym2lean_document.tex $(BUILD)/
	cd $(BUILD) && pdflatex -interaction=nonstopmode -halt-on-error sym2lean_document.tex > sym2lean.log 2>&1 \
	  || { tail -40 sym2lean.log; exit 1; }
	cd $(BUILD) && pdflatex -interaction=nonstopmode -halt-on-error sym2lean_document.tex >> sym2lean.log 2>&1 \
	  || { tail -40 sym2lean.log; exit 1; }
	@echo "built $(BUILD)/sym2lean_document.pdf"

sym2lean-test:
	$(PY) tools/texmath.py --selftest
	$(PY) tools/lean_emit.py --selftest
	$(PY) tools/sym2lean.py --selftest
	$(PY) tools/leandsl.py --selftest
	$(PY) tools/pyast2lean.py --selftest

# @lean DSL (tools/leandsl.py): models in tools/dsl_models/*.py -> build/lean/Dsl_mathlib.lean
dsl:
	$(PY) tools/leandsl.py

# Python source -> LaTeX + Lean via ast, loop invariants from asserts (tools/pyast2lean.py); PYAST_ARGS=--prove
PYAST_ARGS ?=
pyast:
	$(PY) tools/pyast2lean.py $(PYAST_ARGS)

# One document for a guide page: its families, the Lean arithmetic statements read from their files, the chains between them
# (tools/oai_tree.py --lean-items, tools/oai_lean.py). Needs xelatex; the first run reads every family's Lean imports (about 30 s).
GUIDE ?= PHYSICS.md
OAI_ARGS ?=
physics-lean:
	@command -v xelatex >/dev/null || { echo "xelatex not found (install TeX Live)"; exit 1; }
	$(PY) tools/oai_tree.py --guide $(GUIDE) --out $(BUILD)/oai --pdf --with-models --lean-items --graph $(BUILD)/oai/graph.json $(OAI_ARGS)

# Run only the checks of the Hopf models (no export; a partial export currently stops in concise.py). See HOPF_GAUGE.md.
hopf:
	cd sympy_models && $(PY) -c "import importlib, eqlib; [importlib.import_module('models.' + m).build_and_verify() for m in ('hopf_fibration', 'hopf_insulator', 'hopf_gauge')]; print(len(eqlib.checks()), 'checks passed in 3 Hopf models')"

# Build only the quantum gravity paper's inputs (full export; see QUANTUM_GRAVITY.md).
qg:
	$(PY) sympy_models/export.py --pdf

su2:
	cd sympy_models && $(PY) -c "import importlib, eqlib; [importlib.import_module('models.' + m).build_and_verify() for m in ('su2_lattice_gap', 'su2_lattice_ladder')]; print(len(eqlib.checks()), 'checks passed in the 2 SU(2) lattice-gap models')"

# --- summary and PDF --------------------------------------------------------------------------
$(BUILD)/tex/summary.tex: tools/make_summary.py $(STAMP_FILES)
	$(PY) tools/make_summary.py

# pdflatex runs twice (cross-references) inside build/, so the template's relative paths resolve there.
$(PDF): $(TEMPLATE) $(BUILD)/tex/summary.tex $(STAMP_FILES)
	@command -v pdflatex >/dev/null || { echo "pdflatex not found (install TeX Live)"; exit 1; }
	cp $(TEMPLATE) $(BUILD)/physics_models.tex
	cd $(BUILD) && pdflatex -interaction=nonstopmode -halt-on-error physics_models.tex > pdflatex.log 2>&1 \
	  || { tail -40 pdflatex.log; exit 1; }
	cd $(BUILD) && pdflatex -interaction=nonstopmode -halt-on-error physics_models.tex >> pdflatex.log 2>&1 \
	  || { tail -40 pdflatex.log; exit 1; }
	@echo "built $(PDF)"

pdf: $(PDF)

# -B forces every target to rebuild, and FULL=1 is passed to the scripts in this same invocation.
full:
	$(MAKE) -B FULL=1 pdf

clean:
	rm -rf $(BUILD) sympy_models/build
	find sympy_models -name __pycache__ -type d -prune -exec rm -rf {} +
	find tools -name __pycache__ -type d -prune -exec rm -rf {} +
