#!/usr/bin/env python
"""Sensitivity of the primary meta-analysis to the standard-error ceiling.

04_meta_analysis.py drops a gene x stratum observation whose lfcSE is at or above
--max-lfcse (default 5). The filter is documented in Methods but its effect was never
shown, and it is informative censoring rather than neutral clean-up: the observations it
removes are the low-count ones, so the ceiling changes which strata contribute to which
gene and therefore the per-gene k that the >= --min-studies floor is applied to.

This reads the PER-STRATUM DIFFERENTIAL-EXPRESSION TABLES, not the pooled long file the
meta step writes. That distinction is the whole point: the long file has already had the
ceiling applied, so testing ceilings against it can only ever reproduce the default and
would make a looser ceiling look inert by construction. Starting from the DE tables keeps
the observations above the default ceiling available, so "looser" is a real condition.

Reported per ceiling: observations retained and removed, genes clearing the measurability
floor, FDR<0.05 count under DerSimonian-Laird, and the biosynthetic arm's set mean, which
is the quantity the pathway claim rests on.
"""
import argparse, glob, os
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
    ap.add_argument("--de", required=True,
                    help="per-stratum DE directory, BEFORE the ceiling is applied")
    ap.add_argument("--sets", required=True)
    ap.add_argument("--min-studies", type=int, default=8)
    ap.add_argument("--ceilings", default="2,3,5,10,inf")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    files = sorted(glob.glob(os.path.join(a.de, "PUBLIC_*_droughtVScontrol.csv")))
    if not files:
        raise SystemExit(f"no DE tables under {a.de}")
    raw = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    raw = raw.dropna(subset=["log2FoldChange", "lfcSE"])
    raw = raw[raw.lfcSE > 0]
    arm = list(pd.read_csv(a.sets).query("set == 'rfo_biosynthetic'").gene)

    rows = []
    for c in a.ceilings.split(","):
        cap = float(c)
        sub = raw[raw.lfcSE < cap]
        r = pool(sub, a.min_studies)
        armr = r[r.gene.isin(arm)]
        rows.append(dict(max_lfcse=c, n_obs_kept=len(sub), n_obs_removed=len(raw) - len(sub),
                         n_genes=len(r), n_bh05=int((r.padj < 0.05).sum()),
                         arm_members_tested=len(armr),
                         arm_mean_log2FC=float(armr.mu.mean()) if len(armr) else np.nan,
                         genome_median_log2FC=float(r.mu.median())))
    R = pd.DataFrame(rows)
    R.insert(0, "source", "per-stratum DE tables (no ceiling pre-applied)")
    R.to_csv(os.path.join(a.out, "PUBLIC_lfcse_ceiling_sensitivity.csv"), index=False)
    print(f"{len(raw)} unfiltered observations, max lfcSE {raw.lfcSE.max():.1f}")
    print(R.drop(columns="source").to_string(index=False))


if __name__ == "__main__":
    main()
