#!/usr/bin/env python3
"""Subgroup meta-analysis: does a moderator explain the between-study heterogeneity?

For every gene, the strata are split by a moderator (stress method, platform, ...).
Within each level a DerSimonian-Laird random-effects estimate is computed; the levels are
then compared with Cochran's Q_between on the level estimates. A moderator that explains
heterogeneity shows lower within-level I^2 than the pooled I^2 and a significant Q_between.

  python 05_subgroup.py --de de_leafshoot --cov subgroup/stratum_covariates.csv \
      --moderator stress_method --min-strata-per-level 3 --out subgroup/leaf_stress_method
"""
import argparse, glob, os
import numpy as np, pandas as pd
from scipy import stats

def dl_pool(y, se):
    """DerSimonian-Laird random-effects pooling. Returns beta, se, tau2, I2, Q."""
    k = len(y)
    w = 1.0 / se**2
    beta_f = (w * y).sum() / w.sum()
    Q = float((w * (y - beta_f) ** 2).sum())
    if k > 1:
        C = w.sum() - (w**2).sum() / w.sum()
        tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
        I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    else:
        tau2, I2 = 0.0, 0.0
    wr = 1.0 / (se**2 + tau2)
    return float((wr * y).sum() / wr.sum()), float(np.sqrt(1 / wr.sum())), tau2, I2, Q

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--de", required=True)
    ap.add_argument("--cov", required=True)
    ap.add_argument("--moderator", required=True)
    ap.add_argument("--min-strata-per-level", type=int, default=3)
    ap.add_argument("--min-total-strata", type=int, default=8)
    ap.add_argument("--max-lfcse", type=float, default=5.0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    cov = pd.read_csv(args.cov).set_index("stratum")
    frames = []
    for f in sorted(glob.glob(os.path.join(args.de, "PUBLIC_*_droughtVScontrol.csv"))):
        d = pd.read_csv(f).dropna(subset=["log2FoldChange", "lfcSE"])
        d = d[(d.lfcSE > 0) & (d.lfcSE < args.max_lfcse)]
        frames.append(d[["gene", "stratum", "log2FoldChange", "lfcSE", "max_count"]])
    long = pd.concat(frames, ignore_index=True)
    long[args.moderator] = long.stratum.map(cov[args.moderator])
    long = long.dropna(subset=[args.moderator])

    counts = cov.loc[cov.index.isin(long.stratum.unique()), args.moderator].value_counts()
    levels = [l for l, n in counts.items() if n >= args.min_strata_per_level and str(l) != "unknown"]
    print(f"moderator '{args.moderator}' levels kept: "
          + ", ".join(f"{l} ({counts[l]} strata)" for l in levels))
    if len(levels) < 2:
        raise SystemExit("fewer than two levels clear --min-strata-per-level; nothing to compare")
    long = long[long[args.moderator].isin(levels)]

    rows = []
    for gene, g in long.groupby("gene"):
        if g.stratum.nunique() < args.min_total_strata:
            continue
        per, ok = {}, True
        for lv in levels:
            s = g[g[args.moderator] == lv]
            if len(s) < args.min_strata_per_level:
                ok = False; break
            b, se, tau2, I2, _ = dl_pool(s.log2FoldChange.to_numpy(), s.lfcSE.to_numpy())
            per[lv] = (b, se, I2, len(s))
        if not ok:
            continue
        b_all, se_all, tau2_all, I2_all, _ = dl_pool(g.log2FoldChange.to_numpy(), g.lfcSE.to_numpy())
        bs = np.array([per[l][0] for l in levels]); ses = np.array([per[l][1] for l in levels])
        wb = 1.0 / ses**2
        mu = (wb * bs).sum() / wb.sum()
        Qb = float((wb * (bs - mu) ** 2).sum())
        rec = dict(gene=gene, n_strata=int(g.stratum.nunique()),
                   beta_overall=b_all, I2_overall=I2_all,
                   Q_between=Qb, df_between=len(levels) - 1,
                   p_between=float(stats.chi2.sf(Qb, len(levels) - 1)),
                   I2_within_median=float(np.median([per[l][2] for l in levels])))
        for l in levels:
            key = str(l).replace(" ", "_")
            rec[f"beta_{key}"], rec[f"se_{key}"], rec[f"I2_{key}"], rec[f"k_{key}"] = per[l]
        rows.append(rec)

    res = pd.DataFrame(rows)
    if res.empty:
        raise SystemExit("no gene had enough strata in every level")
    res["padj_between"] = stats.false_discovery_control(res.p_between, method="bh")
    res = res.sort_values("p_between")
    res.to_csv(os.path.join(args.out, "PUBLIC_subgroup_by_%s.csv" % args.moderator), index=False)
    print(f"{len(res)} genes tested; {int((res.padj_between < 0.05).sum())} show a significant "
          f"subgroup difference (FDR<0.05)")
    print(f"median I2 overall {res.I2_overall.median():.0f}% -> median within-level "
          f"{res.I2_within_median.median():.0f}%")

if __name__ == "__main__":
    main()
