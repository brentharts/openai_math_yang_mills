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
.PHONY: all help setup models pdf full bloch on kt haldane chair44 figures combinatorics zeta sympy sympy-pdf hopf clean

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
	@echo "make sympy-pdf  also build sympy_models/build/equations_document.pdf, grouped by CONTENTS entry"
	@echo "make figures     regenerate the committed figures/*.png used by SPECTRE_CLOCK_AND_CHAIR44.md (about 90 s)"
	@echo "make combinatorics  exact checks and figures for COMBINATORICS_REPRESENTATIONS.md (about 1 min; needs networkx)"
	@echo "make zeta           zero-counting check and strip figures for ZETA_ZEROS.md (about 30 s; needs mpmath)"
	@echo "make clean      remove $(BUILD)/"
	@echo "Variables: PY=python3  FULL=1"

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

# Run only the checks of the Hopf models (no export; a partial export currently stops in concise.py). See HOPF_GAUGE.md.
hopf:
	cd sympy_models && $(PY) -c "import importlib, eqlib; [importlib.import_module('models.' + m).build_and_verify() for m in ('hopf_fibration', 'hopf_insulator', 'hopf_gauge')]; print(len(eqlib.checks()), 'checks passed in 3 Hopf models')"

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
