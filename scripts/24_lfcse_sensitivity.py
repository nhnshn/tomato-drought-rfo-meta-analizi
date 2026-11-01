#!/usr/bin/env python
"""Sensitivity of the primary meta-analysis to the standard-error ceiling.

04_meta_analysis.py drops a gene x stratum observation when its lfcSE is >= --max-lfcse
(default 5). The filter is documented in Methods but its effect was never shown, and it
is informative censoring rather than a neutral clean-up: the observations it removes are
the low-count ones, so the ceiling changes which strata contribute to which gene, and
therefore the per-gene k that the >= --min-studies floor is applied to.

This re-pools the primary cohort at several ceilings and reports what moves: the number
of genes that clear the floor, the FDR<0.05 count under DerSimonian-Laird, and the
biosynthetic arm's set mean, which is the quantity the pathway claim rests on.
"""
import argparse, os
import numpy as np, pandas as pd
from scipy import stats


def pool(df, min_studies):
    rows = []
    for gene, g in df.groupby("gene", sort=False):
        k = len(g)
        if k < min_studies:
            continue
        y, se = g.log2FoldChange.to_numpy(), g.lfcSE.to_numpy()
        w = 1.0 / se ** 2
        bf = (w * y).sum() / w.sum()
        Q = (w * (y - bf) ** 2).sum()
        C = w.sum() - (w ** 2).sum() / w.sum()
        tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
        wr = 1.0 / (se ** 2 + tau2)
        br = (wr * y).sum() / wr.sum()
        sr = np.sqrt(1 / wr.sum())
        rows.append((gene, k, br, 2 * stats.norm.sf(abs(br / sr))))
    r = pd.DataFrame(rows, columns=["gene", "k", "mu", "p"])
    r["padj"] = stats.false_discovery_control(r.p, method="bh")
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--long", required=True)
    ap.add_argument("--sets", required=True)
    ap.add_argument("--min-studies", type=int, default=8)
    ap.add_argument("--ceilings", default="2,3,5,10,inf")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    long = pd.read_csv(a.long)
    sets = pd.read_csv(a.sets)
    arm = [g for g in sets[sets["set"] == "rfo_biosynthetic"].gene]

    rows = []
    for c in a.ceilings.split(","):
        cap = float(c)
        sub = long[(long.lfcSE > 0) & (long.lfcSE < cap)]
        r = pool(sub, a.min_studies)
        armr = r[r.gene.isin(arm)]
        rows.append(dict(max_lfcse=c, n_obs=len(sub), n_genes=len(r),
                         n_bh05=int((r.padj < 0.05).sum()),
                         arm_members_tested=len(armr),
                         arm_mean_log2FC=float(armr.mu.mean()) if len(armr) else np.nan,
                         genome_median_log2FC=float(r.mu.median())))
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(a.out, "PUBLIC_lfcse_ceiling_sensitivity.csv"), index=False)
    print(R.to_string(index=False))


if __name__ == "__main__":
    main()
