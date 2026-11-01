#!/usr/bin/env python
"""Cohort-definition sensitivity: how much of the pathway result survives when the
strata whose deposited records conflict, or which share BioSamples, are removed.

The shipped package tested this only for the stratum-level matched contrast
(study_contrast_record_sensitivity). The competitive CAMERA test, which carries the
smallest reported p-value, was never re-run under the same exclusions. This script
closes that gap and runs it from the current pipeline outputs rather than from a
previous reanalysis.

Two quantities per scenario:
  * gene level  - DerSimonian-Laird and HKSJ refits from per_study_effects_long,
                  keeping the >= --min-studies floor used by the primary analysis
  * set level   - Stouffer recombination of the per-stratum one-sided CAMERA
                  p-values already computed by 18_competitive_rotation.R

A scenario whose stratum count falls below the gene-level floor is reported as
not computable, never as "no significant genes".
"""
import argparse, json, os
import numpy as np, pandas as pd
from scipy import stats

SCENARIOS = {
    "full": [],
    "drop_record_conflicts": ["PRJNA1164923", "PRJNA1307411"],
    "drop_conflicts_plus_shared_biosample": ["PRJNA1164923", "PRJNA1307411", "PRJNA644865"],
    "drop_conflicts_shared_and_timing": ["PRJNA1164923", "PRJNA1307411",
                                         "PRJNA644865", "PRJNA1381306"],
}
BASIS = {
    "full": "all primary strata",
    "drop_record_conflicts": "BioSample/title record conflicts removed",
    "drop_conflicts_plus_shared_biosample": "record conflicts plus the stratum sharing 8 BioSamples",
    "drop_conflicts_shared_and_timing": "the above plus the weakest control-timing evidence",
}


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
        p_norm = 2 * stats.norm.sf(abs(br / sr))
        # Hartung-Knapp-Sidik-Jonkman variance on the same weights
        q = (wr * (y - br) ** 2).sum() / (k - 1)
        se_hk = np.sqrt(q / wr.sum())
        p_hk = 2 * stats.t.sf(abs(br / se_hk), df=k - 1) if k > 1 else np.nan
        rows.append((gene, k, br, sr, p_norm, p_hk))
    r = pd.DataFrame(rows, columns=["gene", "k", "mu", "se", "p_normal", "p_hksj"])
    if len(r):
        r["bh_normal"] = stats.false_discovery_control(r.p_normal, method="bh")
        r["bh_hksj"] = stats.false_discovery_control(r.p_hksj, method="bh")
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--long", required=True)
    ap.add_argument("--rotation", required=True, help="PUBLIC_rotation_per_stratum.csv")
    ap.add_argument("--min-studies", type=int, default=8)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    long = pd.read_csv(a.long)
    long = long[(long.lfcSE > 0) & (long.lfcSE < 5)]
    rot = pd.read_csv(a.rotation)
    rot["study"] = rot.stratum.str.split("__").str[0]

    gene_rows, set_rows = [], []
    for name, drop in SCENARIOS.items():
        sub = long[~long.study.isin(drop)]
        k_strata = sub.stratum.nunique()
        if k_strata < a.min_studies:
            gene_rows.append(dict(scenario=name, n_strata=k_strata, n_genes=np.nan,
                                  n_bh_normal=np.nan, n_bh_hksj=np.nan,
                                  status=f"not computable: {k_strata} strata is below the "
                                         f"{a.min_studies}-stratum floor",
                                  excluded=";".join(drop), basis=BASIS[name]))
        else:
            r = pool(sub, a.min_studies)
            gene_rows.append(dict(scenario=name, n_strata=k_strata, n_genes=len(r),
                                  n_bh_normal=int((r.bh_normal < 0.05).sum()),
                                  n_bh_hksj=int((r.bh_hksj < 0.05).sum()),
                                  status="ok", excluded=";".join(drop), basis=BASIS[name]))
        rs = rot[~rot.study.isin(drop)]
        for gs, g in rs.groupby("gene_set"):
            p = g.camera_p_up.to_numpy()
            z = stats.norm.isf(p)
            Z = z.sum() / np.sqrt(len(z))
            set_rows.append(dict(scenario=name, gene_set=gs, k=len(p),
                                 n_up=int((g.direction == "Up").sum()),
                                 stouffer_z=Z, combined_p_up=float(stats.norm.sf(Z)),
                                 excluded=";".join(drop)))
    G = pd.DataFrame(gene_rows)
    S = pd.DataFrame(set_rows)
    S["bh_within_scenario"] = np.nan
    for sc, idx in S.groupby("scenario").groups.items():
        S.loc[idx, "bh_within_scenario"] = stats.false_discovery_control(
            S.loc[idx, "combined_p_up"], method="bh")
    G.to_csv(os.path.join(a.out, "PUBLIC_cohort_conflict_gene_level.csv"), index=False)
    S.to_csv(os.path.join(a.out, "PUBLIC_cohort_conflict_set_level.csv"), index=False)
    print(G.to_string(index=False))
    print()
    print(S[S.gene_set == "rfo_biosynthetic"].to_string(index=False))


if __name__ == "__main__":
    main()
