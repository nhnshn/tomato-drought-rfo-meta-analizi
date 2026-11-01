#!/usr/bin/env python
"""Per-gene inference for a pre-specified pathway panel, plus the inference-model
sensitivity that qualifies it.

Three variance models are reported side by side, because on this corpus the number
of genes that clear FDR depends on the uncertainty model far more than on the
tau-squared estimator:

  normal-z       DerSimonian-Laird tau^2, Wald test against the normal
  HKSJ           Hartung-Knapp-Sidik-Jonkman variance inflation, t with k-1 df
  modified HKSJ  the same with the inflation floored at 1 (Roever et al. 2015,
                 doi:10.1186/s12874-015-0091-1), which stops an unusually
                 homogeneous gene from being handed a *smaller* interval

Paule-Mandel tau^2 is estimated as well, so that "the estimator moved it" can be
told apart from "the uncertainty model moved it". The two questions are separate
and the distinction decides how the genome-wide count may be described.

A sign test across strata is reported alongside: it reads only the direction of
each stratum's effect, so it is untouched by all of the above.

Multiplicity is corrected twice, for two different questions. Genome-wide BH
answers "is this gene remarkable among all tested genes". Panel-restricted BH
answers "is this gene remarkable within the pre-specified pathway". The panel must
be fixed in advance for the second to mean anything -- pass the full pathway panel,
never the subset that happened to reach significance in an earlier run.

Usage:
  python 10_panel_inference.py --long meta/per_study_effects_long.csv.gz \
      --panel config/rfo_panel_full.csv --out meta
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

PUBLIC_PREFIX = "PUBLIC_"


def dl_tau2(y, v):
    """DerSimonian-Laird tau^2 (method of moments), truncated at zero."""
    k = len(y)
    w = 1.0 / v
    fe = (w * y).sum() / w.sum()
    Q = (w * (y - fe) ** 2).sum()
    C = w.sum() - (w ** 2).sum() / w.sum()
    return max(0.0, (Q - k + 1) / C) if C > 0 else 0.0


def pm_tau2(y, v, tol=1e-6, it=100):
    """Paule-Mandel tau^2 by Newton iteration on the weighted Q statistic.

    Solves sum w_i (y_i - mu)^2 = k - 1 for tau^2, where w depends on tau^2.
    Less biased than DerSimonian-Laird at the small k typical of a meta-analysis
    over a dozen strata, which is why Moon et al. (2026) use it.
    """
    t2 = 0.0
    for _ in range(it):
        w = 1 / (v + t2)
        mu = (w * y).sum() / w.sum()
        F = (w * (y - mu) ** 2).sum() - (len(y) - 1)
        if abs(F) < tol:
            break
        d = -((w ** 2) * (y - mu) ** 2).sum()
        if d == 0:
            break
        t2 = max(0.0, t2 - F / d)
    return t2


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--long", required=True,
                    help="per_study_effects_long.csv.gz from 04_meta_analysis.py")
    ap.add_argument("--panel", required=True,
                    help="CSV with gene[,annotation] - the PRE-SPECIFIED pathway panel")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-strata", type=int, default=8,
                    help="a gene is pooled only where it was testable in this many strata")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    eff = pd.read_csv(args.long)
    panel_df = pd.read_csv(args.panel)
    panel = list(panel_df.gene)
    ann = dict(zip(panel_df.gene, panel_df.get("annotation", pd.Series(dtype=str))))

    # ---- genome-wide, three variance models + Paule-Mandel ------------------
    rows = []
    for gene, g in eff.groupby("gene"):
        k = len(g)
        if k < args.min_strata:
            continue
        y = g.log2FoldChange.to_numpy()
        v = g.lfcSE.to_numpy() ** 2
        if not (np.isfinite(y).all() and np.isfinite(v).all() and (v > 0).all()):
            continue
        t2 = dl_tau2(y, v)
        wr = 1 / (v + t2)
        mu = (wr * y).sum() / wr.sum()
        se = np.sqrt(1 / wr.sum())
        # HKSJ: inflate the variance by the weighted residual scatter, test on t
        q = (wr * (y - mu) ** 2).sum() / (k - 1)
        t2_pm = pm_tau2(y, v)
        w_pm = 1 / (v + t2_pm)
        mu_pm = (w_pm * y).sum() / w_pm.sum()
        se_pm = np.sqrt(1 / w_pm.sum())
        rows.append((gene, k, mu, se, t2, q,
                     2 * stats.norm.sf(abs(mu / se)),
                     2 * stats.t.sf(abs(mu / (se * np.sqrt(q))), k - 1),
                     2 * stats.t.sf(abs(mu / (se * np.sqrt(max(1.0, q)))), k - 1),
                     mu_pm, t2_pm, 2 * stats.norm.sf(abs(mu_pm / se_pm))))
    H = pd.DataFrame(rows, columns=["gene", "k", "mu", "se", "tau2_dl", "q",
                                    "p_normal", "p_hksj", "p_hksj_mod",
                                    "mu_pm", "tau2_pm", "p_pm"])
    for c in ["normal", "hksj", "hksj_mod", "pm"]:
        H["bh_" + c] = multipletests(H["p_" + c], method="fdr_bh")[1]
    H.to_csv(os.path.join(args.out, "inference_sensitivity_hksj.csv"), index=False)

    # ---- the same genes, BH restricted to the pre-specified panel -----------
    P = H[H.gene.isin(panel)].copy()
    for c in ["normal", "hksj", "hksj_mod", "pm"]:
        P["panelBH_" + c] = multipletests(P["p_" + c], method="fdr_bh")[1]
    P["annotation"] = P.gene.map(lambda g: ann.get(g, ""))
    P = P.sort_values("mu", ascending=False)
    P.to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}panel_inference.csv"), index=False)

    # ---- sign test over strata: no variance model at all --------------------
    rows = []
    for g in panel:
        s = eff[eff.gene == g]
        if len(s) < args.min_strata:
            continue
        pos = int((s.log2FoldChange > 0).sum())
        k = len(s)
        rows.append(dict(gene=g, annotation=str(ann.get(g, ""))[:34], k=k, n_positive=pos,
                         median_lfc=s.log2FoldChange.median(),
                         p_sign=stats.binomtest(pos, k, 0.5).pvalue))
    D = pd.DataFrame(rows)
    if len(D):
        D["padj_sign"] = multipletests(D.p_sign, method="fdr_bh")[1]
        D = D.sort_values("p_sign")
    D.to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}panel_direction_consistency.csv"),
             index=False)

    # ---- one-line summary of what each model concludes ----------------------
    summary = pd.DataFrame([
        dict(scope="genome-wide", model="normal-z (DerSimonian-Laird tau2)",
             n_tested=len(H), n_fdr05=int((H.bh_normal < 0.05).sum()),
             smallest_bh=float(H.bh_normal.min())),
        dict(scope="genome-wide", model="normal-z (Paule-Mandel tau2)",
             n_tested=len(H), n_fdr05=int((H.bh_pm < 0.05).sum()),
             smallest_bh=float(H.bh_pm.min())),
        dict(scope="genome-wide", model="HKSJ",
             n_tested=len(H), n_fdr05=int((H.bh_hksj < 0.05).sum()),
             smallest_bh=float(H.bh_hksj.min())),
        dict(scope="genome-wide", model="modified HKSJ",
             n_tested=len(H), n_fdr05=int((H.bh_hksj_mod < 0.05).sum()),
             smallest_bh=float(H.bh_hksj_mod.min())),
        dict(scope="panel", model="normal-z (DerSimonian-Laird tau2)",
             n_tested=len(P), n_fdr05=int((P.panelBH_normal < 0.05).sum()),
             smallest_bh=float(P.panelBH_normal.min()) if len(P) else np.nan),
        dict(scope="panel", model="HKSJ",
             n_tested=len(P), n_fdr05=int((P.panelBH_hksj < 0.05).sum()),
             smallest_bh=float(P.panelBH_hksj.min()) if len(P) else np.nan),
        dict(scope="panel", model="sign test over strata",
             n_tested=len(D), n_fdr05=int((D.padj_sign < 0.05).sum()) if len(D) else 0,
             smallest_bh=float(D.padj_sign.min()) if len(D) else np.nan),
    ])
    summary.to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}inference_summary.csv"),
                   index=False)

    print(f"genes pooled in >= {args.min_strata} strata: {len(H)}")
    print(summary.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
    print(f"\nmedian q (between-study scale factor): {H.q.median():.2f} | "
          f"q > 1 in {100 * (H.q > 1).mean():.0f}% of genes")
    print(f"panel genes tested: {len(P)} of {len(panel)}")


if __name__ == "__main__":
    main()
