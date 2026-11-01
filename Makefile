# Tomato drought transcriptome meta-analysis - end-to-end pipeline.
#
# Every rule logs the exact command, the git commit, and tool versions to
# logs/<target>.log, because an output is only reproducible if the parameters that
# produced it are recorded with it.
#
#   make reference      download and index ITAG4.1
#   make quant          quantify every library in config/samplesheet.csv  (needs network)
#   make de             per-stratum differential expression, design chosen per stratum
#   make meta           random-effects meta-analysis, leaf, concurrent-control strata
#   make moderators     stress-method subgroup + duration meta-regression
#   make replication    Moon et al. 2026 candidate replication
#   make inference      panel inference + variance-model sensitivity
#   make geneset        competitive gene-set test vs an expression-matched null
#   make loo            leave-one-stratum-out influence analysis
#   make figures        every manuscript figure, from the tables above
#   make all

SHELL      := /bin/bash
.SHELLFLAGS := -o pipefail -c   # a failure inside a `| tee` pipe must fail the rule
PY         := python
QUANT_ENV  ?= tomato-drought-quant
ANALYSIS_ENV ?= tomato-drought-analysis
REF        := ref
MATRICES   := matrices
DE         := de
META       := meta
FIGS       := figures
MIN_REPS   ?= 3
EXPORT_DIR ?= export
MIN_STRATA ?= 8
THREADS    ?= 4

# A literal comma inside $(call ...) is read by make as an argument separator and
# silently truncates the command, so any comma-separated flag value must use this.
COMMA := ,

# The rotation/CAMERA step is the only R step. Point RSCRIPT at an interpreter whose
# library has limma and edgeR; `make all RSCRIPT=/path/to/Rscript` is the usual form.
RSCRIPT    ?= Rscript

define log
@mkdir -p logs
@{ echo "### $@"; date -u +'%Y-%m-%dT%H:%M:%SZ'; \
   echo "commit: $$(git rev-parse --short HEAD 2>/dev/null || echo untracked)"; \
   echo "command: $(1)"; \
   $(PY) -c "import sys;print('python:',sys.version.split()[0])"; \
   } >> logs/$@.log
@set +e; $(1) 2>&1 | tee -a logs/$@.log; statuses=("$${PIPESTATUS[@]}"); rc=$${statuses[0]}; \
  if [[ $$rc -eq 0 && $${statuses[1]} -ne 0 ]]; then rc=$${statuses[1]}; fi; \
  echo "ended_utc: $$(date -u +'%Y-%m-%dT%H:%M:%SZ') | exit_code: $$rc" | tee -a logs/$@.log; \
  exit $$rc
endef

.NOTPARALLEL:
.PHONY: robustness prose sensitivity complete out_of_sample background_sensitivity currency seal all reference quant de meta moderators genotype replication tissue coordination \
        inference geneset loo figures manifest clean

# `all` runs everything downstream of quantification. `quant` is separate because it
# needs network access and downloads ~1.3 TB; run it once, then `make all`.
# `figures` is last and depends on the rest: every panel is drawn from a table
# produced earlier in the same run, so a figure cannot disagree with the text.
all: manifest meta moderators genotype replication inference geneset tissue coordination contrast rotation loo figures
	@$(MAKE) seal

# Complete workflow additionally runs the separate optional analyses before final sealing.
# Explicitly serialised: failure of any phase prevents the final seal.
complete:
	@$(MAKE) all
	@$(MAKE) sensitivity
	@$(MAKE) robustness
	@$(MAKE) out_of_sample
	@$(MAKE) background_sensitivity
	@$(MAKE) calibration
	@test -s rep/out/PUBLIC_rotation_combined.csv
	@test -s meta_rotation/background_sensitivity/PUBLIC_rotation_background_sensitivity.csv
	@test -s meta_calibration/PUBLIC_calibration_summary.csv
	@test -s meta_sensitivity/PUBLIC_cohort_conflict_set_level.csv
	@test -s meta_sensitivity/PUBLIC_lfcse_ceiling_sensitivity.csv
	@test -s meta_sensitivity/PUBLIC_vif_bootstrap.csv
	@test -s meta_calibration/k13/PUBLIC_calibration_summary.csv
	@$(MAKE) seal

# Written first, so the hashes describe the code that the rest of this run used.
# It lands in logs/ rather than $(META)/, which the meta rule wipes on every rerun.
manifest:
	@mkdir -p logs docs
	$(call log,$(PY) scripts/run_manifest.py --out logs/RUN_MANIFEST.json \
	    --inputs $(MATRICES) --extra-inputs rep/matrices rep/stratum_design.csv \
	    --rscript "$(RSCRIPT)" --note "$(RUN_NOTE)")
	@# docs/ is what the distributed archive carries; copying it here rather than by
	@# hand is what stops the shipped manifest describing an older run than the results.
	@cp logs/RUN_MANIFEST.json docs/RUN_MANIFEST.json

reference:
	$(call log,bash scripts/00_prepare_reference.sh $(REF))

quant: reference
	$(call log,$(PY) scripts/01_quantify.py --sheet config/samplesheet.csv \
	    --index $(REF)/salmon_itag41_k31 --out quant --threads $(THREADS))
	$(call log,$(PY) scripts/02_build_matrices.py --quant quant \
	    --sheet config/samplesheet.csv --tx2gene $(REF)/tx2gene.tsv --out $(MATRICES))

de:
	rm -rf $(DE) $(DE)_leaf_primary $(DE)_leaf_baseline $(DE)_leaf_unresolved
	$(call log,$(PY) scripts/03_de_design_aware.py --matrices $(MATRICES) \
	    --samples config/samplesheet.csv --out $(DE) --min-reps $(MIN_REPS))

# The primary cohort is leaf strata whose control arm is concurrent with the drought
# arm. Baseline-referenced strata answer a different question and are pooled separately.
meta: de
	rm -rf $(META) $(META)_baseline $(META)_unresolved
	$(call log,$(PY) scripts/split_by_estimand.py --de $(DE) --tissue leafshoot \
	    --out-primary $(DE)_leaf_primary --out-sensitivity $(DE)_leaf_baseline \
	    --out-unresolved $(DE)_leaf_unresolved)
	$(call log,$(PY) scripts/04_meta_analysis.py --de $(DE)_leaf_primary \
	    --out $(META) --min-studies $(MIN_STRATA))
	$(call log,$(PY) scripts/04_meta_analysis.py --de $(DE)_leaf_baseline \
	    --out $(META)_baseline --min-studies 3)
	$(call log,$(PY) scripts/04_meta_analysis.py --de $(DE)_leaf_unresolved \
	    --out $(META)_unresolved --min-studies 3)

moderators: meta
	$(call log,$(PY) scripts/05_subgroup.py --de $(DE)_leaf_primary \
	    --cov config/stratum_covariates.csv --moderator stress_method \
	    --out $(META)_subgroup_method)
	$(call log,$(PY) scripts/07_metaregression_duration.py --long $(META)/per_study_effects_long.csv.gz \
	    --durations config/duration_curation.csv --out $(META)/metareg_duration.csv)
	$(call log,$(PY) scripts/07_metaregression_duration.py --long $(META)/per_study_effects_long.csv.gz \
	    --durations config/duration_curation.csv --stress-method water_withholding \
	    --out $(META)/metareg_duration_withholding_only.csv)

genotype: de
	$(call log,$(PY) scripts/06_genotype_moderator.py --arm-effects config/genotype_arm_effects.csv.gz \
	    --wildtype-de $(DE) --panel config/rfo_panel_full.csv --out genotype_results)


# Moon et al. candidate-gene comparison (not the out-of-sample cohort; see out_of_sample).
replication: meta
	$(call log,$(PY) scripts/09_moon_replication.py \
	    --meta $(META)/PUBLIC_meta_drought_random_effects.csv \
	    --candidates config/moon_candidates.csv --metadegs config/moon_metaDEGs.csv \
	    --out $(META)/moon_replication.csv --min-strata $(MIN_STRATA))

# Per-gene inference for the pre-specified panel under three variance models. The
# panel file must be the full pathway set: BH within a panel chosen for significance
# is circular, and the whole point of this rule is an honest within-panel threshold.
inference: meta
	$(call log,$(PY) scripts/10_panel_inference.py --long $(META)/per_study_effects_long.csv.gz \
	    --panel config/rfo_panel_full.csv --out $(META) --min-strata $(MIN_STRATA))

geneset: meta
	$(call log,$(PY) scripts/11_geneset_test.py \
	    --meta $(META)/PUBLIC_meta_drought_random_effects.csv \
	    --long $(META)/per_study_effects_long.csv.gz --min-overlap $(MIN_STRATA) \
	    --sets config/geneset_sets.csv --panel config/rfo_panel_full.csv --out $(META))

loo: meta
	$(call log,$(PY) scripts/12_leave_one_out.py --long $(META)/per_study_effects_long.csv.gz \
	    --panel config/rfo_panel_full.csv --out $(META) --min-strata $(MIN_STRATA))

tissue: geneset
	$(call log,$(PY) scripts/15_tissue_comparison.py --de $(DE) \
	    --leaf-geneset $(META)/PUBLIC_geneset_test.csv \
	    --sets config/geneset_sets.csv --panel config/rfo_panel_full.csv \
	    --out $(META)_tissue)

coordination: geneset
	$(call log,$(PY) scripts/16_coordination.py \
	    --long $(META)/per_study_effects_long.csv.gz \
	    --meta $(META)/PUBLIC_meta_drought_random_effects.csv \
	    --sets config/geneset_sets.csv --pairs config/coordination_pairs.csv \
	    --panel config/rfo_panel_full.csv --out $(META))

# Stratum-level version of the competitive test: immune to inter-gene correlation by
# construction, but assumes the strata are independent. PRJNA472702 and PRJNA644865
# re-sequence eight shared BioSamples, so they are declared as one dependent group.
# Sample-level competitive (CAMERA) and self-contained (fry) tests. Reads the design
# table the DE step wrote, so it fits exactly the libraries that were modelled there.
rotation: de
	$(call log,$(RSCRIPT) scripts/18_competitive_rotation.R \
	    --design de/PUBLIC_stratum_design.csv --matrix-dir matrices \
	    --sets config/geneset_sets.csv --estimand concurrent_control --tissue leafshoot \
	    --dependent-group PRJNA472702$(COMMA)PRJNA644865 --out $(META)_rotation)

# Pre-registered out-of-sample cohort: two studies outside the corpus, quantified by the
# same reference and pipeline, scored with the same competitive test. Separate from `all`
# because it needs the two extra count matrices under rep/matrices.
out_of_sample:
	$(call log,$(RSCRIPT) scripts/18_competitive_rotation.R \
	    --design rep/stratum_design.csv --matrix-dir rep/matrices \
	    --sets config/geneset_sets.csv --estimand concurrent_control --tissue leafshoot \
	    --out rep/out)

# Sensitivity of the sample-level test to the background definition: default camera
# background against one with every non-target panel gene removed.
background_sensitivity:
	$(call log,$(RSCRIPT) scripts/18_competitive_rotation.R \
	    --design de/PUBLIC_stratum_design.csv --matrix-dir matrices \
	    --sets config/geneset_sets.csv --estimand concurrent_control --tissue leafshoot \
	    --dependent-group PRJNA472702$(COMMA)PRJNA644865 --drop-nontarget-panel \
	    --out $(META)_rotation/background_sensitivity)
	$(call log,$(PY) scripts/21_merge_background_sensitivity.py \
	    --default $(META)_rotation/PUBLIC_rotation_combined.csv \
	    --dropped $(META)_rotation/background_sensitivity/PUBLIC_rotation_combined.csv \
	    --out $(META)_rotation/background_sensitivity/PUBLIC_rotation_background_sensitivity.csv)

# Closes the provenance chain at the far end: content hashes of every result file and
# of the count matrices the run consumed, written into the same manifest the run opened.
# `all` depends on it, so a complete run cannot finish unsealed.
seal:
	$(call log,$(PY) scripts/run_manifest.py --out logs/RUN_MANIFEST.json --seal \
	    --inputs $(MATRICES) --extra-inputs rep/matrices rep/stratum_design.csv \
	    --required-output de/PUBLIC_stratum_design.csv \
	    --required-output meta/PUBLIC_meta_drought_random_effects.csv)
	@cp logs/RUN_MANIFEST.json docs/RUN_MANIFEST.json

# Does every exported table still match the pipeline output that produces it?
currency:
	$(call log,$(PY) scripts/20_check_export_currency.py --export-dir $(EXPORT_DIR) --repo .)

# Calibration of the across-stratum Stouffer combination under a competitive null.
# Deliberately OUTSIDE `all`: it validates the method rather than producing a result the
# manuscript's tables or figures read, and it takes about 40 minutes on its own.
calibration:
	$(call log,$(RSCRIPT) scripts/19_calibration.R --reps 400 --strata 8 \
	    --genes 6000 --set-size 7 --per-arm 6 --set-correlation 0.55 \
	    --seed 1 --out $(META)_calibration)
# The simulation above runs at 8 strata, which is the floor a gene must clear, not the
# number the real combination uses. Repeating it at the cohort's actual 13 strata removes
# that mismatch; both are reported because the 8-stratum run is the one earlier text cites.
	$(call log,$(RSCRIPT) scripts/19_calibration.R --reps 400 --strata 13 \
	    --genes 6000 --set-size 7 --per-arm 6 --set-correlation 0.55 \
	    --seed 1 --out meta_calibration/k13)

contrast: meta
	$(call log,$(PY) scripts/17_study_contrast.py \
	    --long $(META)/per_study_effects_long.csv.gz \
	    --sets config/geneset_sets.csv --panel config/rfo_panel_full.csv \
	    --samplesheet config/samplesheet.csv --matrix-dir matrices \
	    --dependent-group PRJNA472702$(COMMA)PRJNA644865 \
	    --restrict-timepoints "PRJNA657834:T1" \
	    --out $(META))

figures: inference geneset loo moderators genotype tissue coordination
	$(call log,$(PY) scripts/13_figures.py --meta $(META) --de $(DE) \
	    --subgroup $(META)_subgroup_method --genotype genotype_results \
	    --tissue $(META)_tissue --sets config/geneset_sets.csv \
	    --panel config/rfo_panel_full.csv --out $(FIGS))


clean:
	rm -rf $(DE) $(DE)_leaf_primary $(DE)_leaf_baseline $(DE)_leaf_unresolved \
	    $(META) $(META)_baseline $(META)_unresolved $(META)_subgroup_method \
	    genotype_results \
	    $(FIGS) logs

# Cohort-definition sensitivity: record-conflict and shared-BioSample exclusions
# applied to BOTH the gene-level refit and the competitive CAMERA recombination.
# No prerequisites on purpose. Every target here is .PHONY, and the `meta` rule starts by
# deleting $(META), so naming it as a prerequisite would re-run it and wipe the inference,
# gene-set, leave-one-out, contrast and coordination tables that other rules wrote into the
# same directory. `complete` already sequences this after `all`; the inputs are checked.
sensitivity:
	@test -s $(META)/per_study_effects_long.csv.gz
	@test -s $(META)_rotation/PUBLIC_rotation_per_stratum.csv
	$(call log,$(PY) scripts/22_cohort_conflict_sensitivity.py \
	    --long $(META)/per_study_effects_long.csv.gz \
	    --rotation $(META)_rotation/PUBLIC_rotation_per_stratum.csv \
	    --min-studies $(MIN_STRATA) --out meta_sensitivity)

# Third link of the provenance chain. The SHA-256 inventory proves a copy is faithful and
# 20_check_export_currency.py proves an exported table matches the pipeline output; neither
# checks that the numbers written in the prose are the numbers in those tables, which is
# where the 1 November audit found 32 stale values. Non-zero exit stops the build.
prose:
	$(call log,$(PY) scripts/23_check_prose_against_tables.py \
	    --claims config/prose_claims.csv --tables $(EXPORT_DIR)/04_donmus_tablolar \
	    --docs $(EXPORT_DIR)/01_makale_taslagi/BMCPB_results_v1.md \
	          $(EXPORT_DIR)/02_bulgular_ve_dergi/BULGULAR.md \
	    --out logs/prose_check.csv)

# Two gaps the audit recorded as open: the standard-error ceiling was documented but its
# effect never shown, and the variance-inflation correction was reported as a point
# estimate although it moves the p-value by three orders of magnitude.
robustness:
	@test -d $(DE)_leaf_primary
	@test -s $(META)/per_study_effects_long.csv.gz
	@test -s $(META)/PUBLIC_meta_drought_random_effects.csv
	$(call log,$(PY) scripts/24_lfcse_sensitivity.py \
	    --de $(DE)_leaf_primary --sets config/geneset_sets.csv \
	    --min-studies $(MIN_STRATA) --out meta_sensitivity)
	$(call log,$(PY) scripts/25_vif_bootstrap.py \
	    --meta $(META)/PUBLIC_meta_drought_random_effects.csv \
	    --long $(META)/per_study_effects_long.csv.gz --sets config/geneset_sets.csv \
	    --panel config/rfo_panel_full.csv --min-overlap $(MIN_STRATA) \
	    --out meta_sensitivity)
