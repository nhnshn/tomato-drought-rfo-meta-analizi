#!/usr/bin/env python3
"""Check that every exported FROZEN_* table is byte-identical to the pipeline
output it claims to be, and that the exported copy carries the cohort the
manuscript reports.

A SHA-256 inventory proves the copy is faithful to the artifact it came from.
It does NOT prove the artifact came from the current run -- that is the gap
this script closes.  Exit status is non-zero if any table is stale.
"""
import argparse, hashlib, json, os, sys
import pandas as pd

# FROZEN_* export name -> path of the pipeline output that produces it.
# None means the table is produced outside the Makefile (recorded, not regenerable).
MAP = {
    "FROZEN_meta_leaf_primary.csv":           "meta/PUBLIC_meta_drought_random_effects.csv",
    "FROZEN_per_study_effects_long.csv.gz":   "meta/per_study_effects_long.csv.gz",
    "FROZEN_meta_root_cohort.csv":            "meta_tissue/root/PUBLIC_meta_drought_random_effects.csv",
    "FROZEN_cohort_conflict_gene_level.csv":  "meta_sensitivity/PUBLIC_cohort_conflict_gene_level.csv",
    "FROZEN_cohort_conflict_set_level.csv":   "meta_sensitivity/PUBLIC_cohort_conflict_set_level.csv",
    "FROZEN_lfcse_ceiling_sensitivity.csv":   "meta_sensitivity/PUBLIC_lfcse_ceiling_sensitivity.csv",
    "FROZEN_vif_bootstrap.csv":               "meta_sensitivity/PUBLIC_vif_bootstrap.csv",
    "FROZEN_calibration_summary_k13.csv":     "meta_calibration/k13/PUBLIC_calibration_summary.csv",
    "FROZEN_calibration_pvalues_k13.csv":     "meta_calibration/k13/calibration_pvalues.csv",
    "FROZEN_inference_summary.csv":           "meta/PUBLIC_inference_summary.csv",
    "FROZEN_panel_inference.csv":             "meta/PUBLIC_panel_inference.csv",
    "FROZEN_panel_leave_one_out.csv":         "meta/PUBLIC_panel_leave_one_out.csv",
    "FROZEN_panel_direction_consistency.csv": "meta/PUBLIC_panel_direction_consistency.csv",
    "FROZEN_leave_one_stratum_out.csv":       "meta/PUBLIC_leave_one_stratum_out_genomewide.csv",
    "FROZEN_geneset_test.csv":                "meta/PUBLIC_geneset_test.csv",
    "FROZEN_coordination_test.csv":           "meta/PUBLIC_coordination_test.csv",
    "FROZEN_coordination_between_sets.csv":   "meta/PUBLIC_coordination_between_sets.csv",
    "FROZEN_panel_correlation_matrix.csv":    "meta/PUBLIC_panel_correlation_matrix.csv",
    "FROZEN_study_contrasts.csv":             "meta/PUBLIC_study_contrasts.csv",
    "FROZEN_study_contrast_summary.csv":      "meta/PUBLIC_study_contrast_summary.csv",
    "FROZEN_study_contrast_loo.csv":          "meta/PUBLIC_study_contrast_loo.csv",
    "FROZEN_target_measurability.csv":        "meta/PUBLIC_target_measurability.csv",
    "FROZEN_de_summary.csv":                  "de/PUBLIC_de_summary.csv",
    "FROZEN_stratum_design.csv":              "de/PUBLIC_stratum_design.csv",
    "FROZEN_subgroup_stress_method.csv":      "meta_subgroup_method/PUBLIC_subgroup_by_stress_method.csv",
    "FROZEN_metaregression_duration.csv":     "meta/metareg_duration.csv",
    "FROZEN_genotype_response.csv":           "genotype_results/PUBLIC_genotype_response_preservation.csv",
    "FROZEN_moon_replication.csv":            "meta/moon_replication.csv",
    "FROZEN_moon_metaDEG_summary.csv":        "meta/moon_replication_metaDEG_summary.csv",
    "FROZEN_rotation_combined.csv":           "meta_rotation/PUBLIC_rotation_combined.csv",
    "FROZEN_rotation_per_stratum.csv":        "meta_rotation/PUBLIC_rotation_per_stratum.csv",
    "FROZEN_calibration_summary.csv":         "meta_calibration/PUBLIC_calibration_summary.csv",
    "FROZEN_calibration_pvalues.csv":         "meta_calibration/calibration_pvalues.csv",
    "FROZEN_replication_combined.csv":        "rep/out/PUBLIC_rotation_combined.csv",
    "FROZEN_replication_per_stratum.csv":     "rep/out/PUBLIC_rotation_per_stratum.csv",
    "FROZEN_tissue_geneset_test.csv":         "meta_tissue/PUBLIC_tissue_geneset_test.csv",
    "FROZEN_tissue_leave_one_out.csv":        "meta_tissue/PUBLIC_tissue_leave_one_out.csv",
    "FROZEN_tissue_per_stratum.csv":          "meta_tissue/PUBLIC_tissue_per_stratum.csv",
    "FROZEN_replication_set_coverage.csv":        None,
    "FROZEN_replication_quant_log.tsv":           None,
    "FROZEN_study_contrast_record_sensitivity.csv": None,
    "FROZEN_rotation_background_sensitivity.csv": "meta_rotation/background_sensitivity/PUBLIC_rotation_background_sensitivity.csv",
    "FROZEN_export_currency.csv":            None,
    "FROZEN_RUN_MANIFEST.json":               "logs/RUN_MANIFEST.json",
}

# Cohort markers the manuscript depends on: column -> expected maximum.
COHORT = {
    "FROZEN_meta_leaf_primary.csv":  ("n_strata", 13),
    "FROZEN_panel_inference.csv":    ("k", 13),
    "FROZEN_panel_leave_one_out.csv":("k", 13),
}

def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--export-dir", required=True, help="exported package root")
    ap.add_argument("--repo", required=True, help="pipeline working tree")
    ap.add_argument("--expect-strata", type=int, default=13)
    ap.add_argument("--out")
    a = ap.parse_args()

    found = {}
    for root, _, files in os.walk(a.export_dir):
        for f in files:
            # Figures are results too. Checking only FROZEN_* tables is how a figure
            # drawn from a superseded cohort survives an export that reports no problems.
            if f.startswith("FROZEN_") or f.endswith(".png"):
                found[f] = os.path.join(root, f)

    rows = []
    for name, path in sorted(found.items()):
        if name.endswith(".png"):
            full = os.path.join(a.repo, "figures", name)
            if not os.path.exists(full):
                rows.append(dict(table=name, verdict="NO-OUTPUT", detail="no figures/" + name))
            else:
                rows.append(dict(table=name,
                                 verdict="CURRENT" if sha(path) == sha(full) else "STALE",
                                 detail="figures/" + name))
            continue
        if name not in MAP:
            rows.append(dict(table=name, verdict="UNMAPPED", detail="not in MAP"))
            continue
        src = MAP[name]
        if src is None:
            rows.append(dict(table=name, verdict="AD-HOC",
                             detail="produced outside the Makefile; not regenerable by make"))
            continue
        full = os.path.join(a.repo, src)
        if not os.path.exists(full):
            rows.append(dict(table=name, verdict="NO-OUTPUT", detail=f"missing {src}"))
            continue
        same = sha(path) == sha(full)
        detail = src
        if name in COHORT and same:
            col, _ = COHORT[name]
            t = pd.read_csv(path)
            if col in t.columns:
                got = int(t[col].max())
                if got != a.expect_strata:
                    rows.append(dict(table=name, verdict="WRONG-COHORT",
                                     detail=f"{col} max {got}, expected {a.expect_strata}"))
                    continue
                detail += f" ({col} max {got})"
        rows.append(dict(table=name, verdict="CURRENT" if same else "STALE", detail=detail))

    df = pd.DataFrame(rows)
    if a.out:
        df.to_csv(a.out, index=False)
    print(df.to_string(index=False))
    bad = df[df.verdict.isin(["STALE", "WRONG-COHORT", "UNMAPPED", "NO-OUTPUT"])]
    print(f"\n{len(df)} tables | CURRENT {int((df.verdict=='CURRENT').sum())} "
          f"| AD-HOC {int((df.verdict=='AD-HOC').sum())} | problems {len(bad)}")
    return 1 if len(bad) else 0

if __name__ == "__main__":
    sys.exit(main())
