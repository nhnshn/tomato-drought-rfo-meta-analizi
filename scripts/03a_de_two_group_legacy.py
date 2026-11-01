#!/usr/bin/env python3
"""Per-study differential expression: drought vs well-watered control.

House convention for public bulk RNA-seq: ONE multi-level group factor plus an
explicit numeric contrast vector (never an interaction term whose reference level
the encoder picks), BH-controlled FDR, PUBLIC_-prefixed outputs, and a per-gene
expression floor so an unexpressed gene is distinguishable from an unchanged one.

  python 03_de_per_study.py --matrices matrices --out de --min-total-count 10
"""
import argparse, glob, os, re
import numpy as np, pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

PUBLIC_PREFIX = "PUBLIC_"

def group_contrast_vector(dds, plus, minus):
    cols = list(dds.obsm["design_matrix"].columns)
    v = np.zeros(len(cols))
    for lvl, sign in [(p, +1.0) for p in plus] + [(m, -1.0) for m in minus]:
        hits = [i for i, c in enumerate(cols) if c.endswith(f"group[T.{lvl}]") or c == f"group_{lvl}"]
        if hits:
            v[hits[0]] += sign / (len(plus) if sign > 0 else len(minus))
        elif not any(c.endswith(f"[T.{lvl}]") for c in cols):
            pass  # reference level: absorbed in the intercept, weight stays on the others
    if abs(v.sum()) > 1e-9:
        raise ValueError(f"contrast weights sum to {v.sum():.3f}, not 0 - that tests a mean, not a contrast")
    return v, dict(zip(cols, v))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrices", default="matrices")
    ap.add_argument("--out", default="de")
    ap.add_argument("--min-total-count", type=int, default=10)
    ap.add_argument("--min-reps", type=int, default=3)
    ap.add_argument("--split-by", default="tissue_h",
                    help="comma-separated columns defining a stratum within a study. Default tissue_h: "
                         "one contrast per study x tissue (33 strata, 495 reps). Adding 'timepoint' "
                         "splits further but only 25 strata then keep >=3v3.")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    summary = []
    for cfile in sorted(glob.glob(os.path.join(args.matrices, "*_counts.csv.gz"))):
        study = os.path.basename(cfile).replace("_counts.csv.gz", "")
        meta = pd.read_csv(os.path.join(args.matrices, f"{study}_samples.csv"))
        counts = pd.read_csv(cfile, index_col=0)
        meta = meta[meta.sample_id.isin(counts.columns)].copy()
        meta = meta[meta.treatment_class.isin(["control", "drought"])]
        if "include_in_core" in meta.columns:
            meta = meta[meta.include_in_core.astype(str).str.lower().isin(["true", "1"])]
        # one multi-level group factor; a timepoint/tissue split stays inside the level name
        split_cols = [c for c in args.split_by.split(",") if c and c in meta.columns]
        parts = split_cols + ["treatment_class"]
        meta["group"] = meta[parts].astype(str).agg("_".join, axis=1).str.replace(r"[^0-9A-Za-z]+", "", regex=True)
        for strat, sm in meta.groupby(split_cols or [lambda _: "all"]):
            tag = study + "__" + "_".join(re.sub(r"\W+", "", str(x)) for x in (strat if isinstance(strat, tuple) else (strat,)))
            n_c = (sm.treatment_class == "control").sum(); n_d = (sm.treatment_class == "drought").sum()
            if n_c < args.min_reps or n_d < args.min_reps:
                summary.append(dict(stratum=tag, n_control=n_c, n_drought=n_d, status="skipped: <min reps"))
                continue
            X = counts[sm.sample_id].T.round().astype(int)
            keep = X.sum(axis=0) >= args.min_total_count
            # A stratum can arrive with nothing mappable - e.g. PRJNA395638 is ABI SOLiD
            # colorspace data, which salmon cannot map against a nucleotide index, so every
            # count is zero. Record it as a QC exclusion instead of crashing PyDESeq2.
            if int(keep.sum()) == 0:
                summary.append(dict(stratum=tag, n_control=n_c, n_drought=n_d, genes_tested=0,
                                    status="excluded: no gene reaches the count floor "
                                           "(check platform/library type)"))
                print(f"{tag}: EXCLUDED - zero genes above the count floor", flush=True)
                continue
            dds = DeseqDataSet(counts=X.loc[:, keep], metadata=sm.set_index("sample_id"),
                               design="~group", refit_cooks=True, n_cpus=1, quiet=True)
            dds.deseq2()
            lv_d = sm.loc[sm.treatment_class == "drought", "group"].iloc[0]
            lv_c = sm.loc[sm.treatment_class == "control", "group"].iloc[0]
            st = DeseqStats(dds, contrast=["group", lv_d, lv_c], n_cpus=1, quiet=True)
            st.summary()
            res = st.results_df.copy()
            res["max_count"] = X.loc[:, keep].max(axis=0)      # unexpressed != unchanged
            res["n_control"], res["n_drought"], res["stratum"], res["study"] = n_c, n_d, tag, study
            res.index.name = "gene"
            res.to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}{tag}_droughtVScontrol.csv"))
            sig = int((res.padj < 0.05).sum())
            summary.append(dict(stratum=tag, n_control=n_c, n_drought=n_d, genes_tested=int(keep.sum()),
                                n_padj_lt_05=sig, status="ok"))
            print(f"{tag}: {n_c}v{n_d}, {keep.sum()} genes tested, {sig} at FDR<0.05", flush=True)
    pd.DataFrame(summary).to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}de_summary.csv"), index=False)

if __name__ == "__main__":
    main()
