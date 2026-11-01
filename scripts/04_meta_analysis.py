#!/usr/bin/env python3
"""Random-effects meta-analysis of per-study drought log2 fold-changes.

Inverse-variance weighting; DerSimonian-Laird tau^2; Cochran's Q and I^2 reported
per gene. Fixed-effect estimates are reported alongside so a reader can see where
the two diverge. Genes are kept when they were testable in >= --min-studies strata.

  python 04_meta_analysis.py --de de --out meta --min-studies 5
"""
import argparse, glob, os
import numpy as np, pandas as pd
from scipy import stats

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--de", default="de")
    ap.add_argument("--out", default="meta")
    ap.add_argument("--min-studies", type=int, default=5)
    ap.add_argument("--max-lfcse", type=float, default=5.0, help="drop absurd SEs from near-zero counts")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    frames = []
    for f in sorted(glob.glob(os.path.join(args.de, "PUBLIC_*_droughtVScontrol.csv"))):
        d = pd.read_csv(f)
        d = d.dropna(subset=["log2FoldChange", "lfcSE"])
        d = d[(d.lfcSE > 0) & (d.lfcSE < args.max_lfcse)]
        frames.append(d[["gene", "stratum", "study", "log2FoldChange", "lfcSE", "max_count"]])
    long = pd.concat(frames, ignore_index=True)
    # mtime=0: gzip stores the write time in its header, so without this two
    # identical runs produce different bytes and a checksum cannot verify a rerun.
    long.to_csv(os.path.join(args.out, "per_study_effects_long.csv.gz"), index=False,
                compression={"method": "gzip", "mtime": 0})

    rows = []
    for gene, g in long.groupby("gene"):
        k = len(g)
        if k < args.min_studies:
            continue
        y, se = g.log2FoldChange.to_numpy(), g.lfcSE.to_numpy()
        w = 1.0 / se**2
        beta_f = float((w * y).sum() / w.sum()); se_f = float(np.sqrt(1 / w.sum()))
        Q = float((w * (y - beta_f) ** 2).sum())
        C = w.sum() - (w**2).sum() / w.sum()
        tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
        I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
        wr = 1.0 / (se**2 + tau2)
        beta_r = float((wr * y).sum() / wr.sum()); se_r = float(np.sqrt(1 / wr.sum()))
        z = beta_r / se_r
        rows.append(dict(gene=gene, n_strata=k, beta_fixed=beta_f, se_fixed=se_f,
                         beta_random=beta_r, se_random=se_r, z=z,
                         pval=float(2 * stats.norm.sf(abs(z))),
                         Q=Q, df=k - 1, pval_Q=float(stats.chi2.sf(Q, k - 1)),
                         tau2=tau2, I2=I2,
                         prop_same_direction=float((np.sign(y) == np.sign(beta_r)).mean()),
                         median_max_count=float(g.max_count.median())))
    if not rows:
        raise SystemExit(f"no gene was testable in >= {args.min_studies} strata "
                         f"({long.stratum.nunique()} strata available) - lower --min-studies")
    res = pd.DataFrame(rows).sort_values("pval")
    m = res.pval.notna()
    res.loc[m, "padj"] = stats.false_discovery_control(res.loc[m, "pval"], method="bh")
    res.to_csv(os.path.join(args.out, "PUBLIC_meta_drought_random_effects.csv"), index=False)
    print(f"{len(res)} genes meta-analysed (>= {args.min_studies} strata); "
          f"{int((res.padj < 0.05).sum())} at random-effects FDR < 0.05; "
          f"median I2 = {res.I2.median():.0f}%")

if __name__ == "__main__":
    main()
