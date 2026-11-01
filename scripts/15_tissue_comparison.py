#!/usr/bin/env python
"""Is the RFO response leaf-specific? Tissue-resolved pathway test.

The primary cohort is concurrent-control leaf/shoot only, so it cannot by itself
distinguish a whole-plant RFO response from a leaf-specific one. This step repeats
the competitive gene-set test in the other organs that reach the same evidence bar.

Only root reaches three concurrent-control strata, the minimum this project accepts
for a random-effects pooled estimate; fruit/seed (2), flower/anther (1) and guard
cell/other (1) are reported per stratum, descriptively, and are NOT pooled.

Because the root cohort rests on three strata, every root claim is accompanied by a
leave-one-stratum-out refit: with k=3 a single study can carry a pooled effect, and
a set-level p-value does not reveal that on its own.

Usage:
  python 15_tissue_comparison.py --de de --leaf-meta meta/PUBLIC_meta_drought_random_effects.csv \
      --leaf-geneset meta/PUBLIC_geneset_test.csv --out meta_tissue
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

PUBLIC_PREFIX = "PUBLIC_"
HERE = os.path.dirname(os.path.abspath(__file__))
# leaf/shoot is listed so its strata appear in the per-stratum description; its
# pooled test is the primary cohort and is read in, not recomputed.
TISSUES = ["leafshoot", "root", "fruitseed", "floweranther", "other"]
ALREADY_POOLED = "leafshoot"
MIN_STRATA = 3


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"failed: {' '.join(cmd)}\n{r.stderr}")
    return r.stdout


def split(de, tissue, root):
    d = {k: os.path.join(root, tissue, k) for k in ("p", "s", "u")}
    for v in d.values():
        # cleared, not reused: a stratum left by an earlier cohort would otherwise be
        # picked up by the glob below and silently pooled into this run
        shutil.rmtree(v, ignore_errors=True)
        os.makedirs(v)
    run([sys.executable, os.path.join(HERE, "split_by_estimand.py"), "--de", de,
         "--tissue", tissue, "--out-primary", d["p"], "--out-sensitivity", d["s"],
         "--out-unresolved", d["u"]])
    return sorted(glob.glob(os.path.join(d["p"], "*.csv")))


def meta_and_test(de_dir, out_dir, sets, panel, min_studies):
    os.makedirs(out_dir, exist_ok=True)
    run([sys.executable, os.path.join(HERE, "04_meta_analysis.py"), "--de", de_dir,
         "--out", out_dir, "--min-studies", str(min_studies)])
    # min-overlap follows the cohort's own coverage rule: in a three-stratum cohort a
    # pair cannot share eight strata, and leaving the default would silently set every
    # inter-gene correlation to NaN and switch the correlation correction off.
    run([sys.executable, os.path.join(HERE, "11_geneset_test.py"),
         "--meta", os.path.join(out_dir, f"{PUBLIC_PREFIX}meta_drought_random_effects.csv"),
         "--long", os.path.join(out_dir, "per_study_effects_long.csv.gz"),
         "--min-overlap", str(min_studies),
         "--sets", sets, "--panel", panel, "--out", out_dir])
    return pd.read_csv(os.path.join(out_dir, f"{PUBLIC_PREFIX}geneset_test.csv"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--de", default="de")
    ap.add_argument("--leaf-geneset", default="meta/PUBLIC_geneset_test.csv")
    ap.add_argument("--sets", default="config/geneset_sets.csv")
    ap.add_argument("--panel", default="config/rfo_panel_full.csv")
    ap.add_argument("--out", default="meta_tissue")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    work = os.path.join(args.out, "_strata")

    # list, not set: a set's iteration order varies with the process hash seed, which
    # reorders the summation in the stratum means and changes their last bits.
    bios = list(dict.fromkeys(pd.read_csv(args.sets).query("set == 'rfo_biosynthetic'").gene))

    # ---- per-stratum description, every concurrent-control non-leaf stratum ----
    per_stratum, pooled_dirs = [], {}
    for t in TISSUES:
        files = split(args.de, t, work)
        pooled_dirs[t] = (os.path.join(work, t, "p"), len(files))
        for f in files:
            st = os.path.basename(f).replace(PUBLIC_PREFIX, "").replace("_droughtVScontrol.csv", "")
            d = pd.read_csv(f)
            gcol = "gene" if "gene" in d.columns else d.columns[0]
            d = d.dropna(subset=["log2FoldChange"]).set_index(gcol)
            sub = d.reindex([g for g in bios if g in d.index])
            per_stratum.append(dict(
                stratum=st, tissue=t, n_genes_measured=len(d),
                n_biosynthetic_measured=int(sub.log2FoldChange.notna().sum()),
                mean_biosynthetic_log2FC=float(sub.log2FoldChange.mean()),
                genome_median_log2FC=float(d.log2FoldChange.median()),
                pooled=(len(files) >= MIN_STRATA and t != ALREADY_POOLED)))
    PS = pd.DataFrame(per_stratum)
    PS.to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}tissue_per_stratum.csv"), index=False)

    # ---- pooled set-level test wherever MIN_STRATA is reached ----
    frames = [pd.read_csv(args.leaf_geneset).assign(
        tissue="leafshoot", n_strata_in_cohort=np.nan)]
    for t, (de_dir, n) in pooled_dirs.items():
        if n < MIN_STRATA or t == ALREADY_POOLED:
            continue
        R = meta_and_test(de_dir, os.path.join(args.out, t), args.sets, args.panel, n)
        frames.append(R.assign(tissue=t, n_strata_in_cohort=n))
    S = pd.concat(frames, ignore_index=True)
    keep = ["tissue", "gene_set", "family", "n_strata_in_cohort", "n_members_declared",
            "n_members_tested", "mean_pooled_log2FC", "null_mean", "null_sd",
            "sd_above_null", "mean_intergene_r", "n_intergene_pairs", "vif",
            "sd_above_null_vif", "alternative", "n_draws", "seed", "p_empirical",
            "p_empirical_vif", "p_bh_within_family", "genome_median_pooled"]
    S[keep].to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}tissue_geneset_test.csv"), index=False)

    # ---- leave-one-stratum-out on every pooled non-leaf cohort ----
    loo = []
    for t, (de_dir, n) in pooled_dirs.items():
        if n < MIN_STRATA or t == ALREADY_POOLED:
            continue
        files = sorted(glob.glob(os.path.join(de_dir, "*.csv")))
        for drop in files:
            st = os.path.basename(drop).replace(PUBLIC_PREFIX, "").replace("_droughtVScontrol.csv", "")
            with tempfile.TemporaryDirectory() as td:
                sub_de = os.path.join(td, "de")
                os.makedirs(sub_de)
                for f in files:
                    if f != drop:
                        os.symlink(os.path.abspath(f), os.path.join(sub_de, os.path.basename(f)))
                # Coverage rule for a small cohort: a gene must be measured in
                # EVERY stratum of that cohort, so the bar follows the refit size.
                R = meta_and_test(sub_de, os.path.join(td, "m"), args.sets, args.panel, n - 1)
                for _, r in R[R.gene_set.isin(["rfo_biosynthetic", "gt8_control"])].iterrows():
                    loo.append(dict(tissue=t, dropped_stratum=st, n_strata_remaining=n - 1,
                                    gene_set=r.gene_set, n_members_tested=r.n_members_tested,
                                    mean_pooled_log2FC=r.mean_pooled_log2FC,
                                    sd_above_null=r.sd_above_null,
                                    mean_intergene_r=r.mean_intergene_r, vif=r.vif,
                                    sd_above_null_vif=r.sd_above_null_vif,
                                    p_empirical=r.p_empirical,
                                    p_empirical_vif=r.p_empirical_vif))
    L = pd.DataFrame(loo)
    L.to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}tissue_leave_one_out.csv"), index=False)

    print("== concurrent-control strata outside leaf/shoot ==")
    print(PS[["stratum", "tissue", "n_biosynthetic_measured", "mean_biosynthetic_log2FC",
              "genome_median_log2FC", "pooled"]].round(3).to_string(index=False))
    print("\n== pooled biosynthetic-arm test by tissue ==")
    b = S[S.gene_set == "rfo_biosynthetic"]
    print(b[["tissue", "n_strata_in_cohort", "n_members_tested", "mean_pooled_log2FC",
             "mean_intergene_r", "vif", "sd_above_null_vif",
             "p_empirical_vif"]].round(3).to_string(index=False))
    if len(L):
        print("\n== leave-one-stratum-out, non-leaf pooled cohorts ==")
        print(L[L.gene_set == "rfo_biosynthetic"][
            ["tissue", "dropped_stratum", "n_members_tested", "mean_pooled_log2FC",
             "vif", "sd_above_null_vif", "p_empirical_vif"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
