#!/usr/bin/env python
"""Leave-one-stratum-out influence analysis.

A corpus assembled from public studies is not evenly weighted: a stratum with deep
coverage and tight standard errors can dominate the pooled estimate of thousands of
genes. This refits the whole meta-analysis once per stratum, each time with that
stratum removed, and reports both

  - the genome-wide significant-gene count, which on this corpus is NOT stable and
    is therefore reportable only as a screen, and
  - per-gene stability for the pathway panel, as the number of refits in which each
    gene stays significant.

The second is the one that supports a claim about a gene. A gene significant in
every refit does not depend on any single study; a gene significant in a handful
does, and the refit count says so plainly.

The pooling here is DerSimonian-Laird with normal-z, matching 04_meta_analysis.py,
and is vectorised over genes so that k+1 full refits stay cheap.

Usage:
  python 12_leave_one_out.py --long meta/per_study_effects_long.csv.gz \
      --panel config/rfo_panel_full.csv --out meta
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

PUBLIC_PREFIX = "PUBLIC_"


def dl_pool(y, se, mask, min_k):
    """Vectorised DerSimonian-Laird pooling over a gene x stratum matrix.

    `mask` marks usable cells. Returns pooled effect, BH-adjusted p (computed only
    among genes meeting `min_k`, since BH over a mixed coverage set is not a single
    family), and the coverage flag.
    """
    safe_se = np.where(mask, se, 1.0)
    w = np.where(mask, 1 / safe_se ** 2, 0.0)
    k = mask.sum(1)
    sw = w.sum(1)
    y0 = np.nan_to_num(y)
    fe = np.divide((w * y0).sum(1), sw, out=np.full(len(y), np.nan), where=sw > 0)
    Q = (w * (y0 - fe[:, None]) ** 2 * mask).sum(1)
    C = sw - np.divide((w ** 2).sum(1), sw, out=np.zeros_like(sw), where=sw > 0)
    tau2 = np.clip(np.divide(Q - (k - 1), C, out=np.zeros_like(Q), where=C > 0), 0, None)
    ws = np.where(mask, 1 / (safe_se ** 2 + tau2[:, None]), 0.0)
    sws = ws.sum(1)
    mu = np.divide((ws * y0).sum(1), sws, out=np.full(len(y), np.nan), where=sws > 0)
    semu = np.divide(1, np.sqrt(sws), out=np.full(len(y), np.nan), where=sws > 0)
    p = 2 * stats.norm.sf(np.abs(mu / semu))
    keep = k >= min_k
    padj = np.full(len(y), np.nan)
    if keep.any():
        padj[keep] = multipletests(p[keep], method="fdr_bh")[1]
    return mu, padj, keep


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--long", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-strata", type=int, default=8)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    eff = pd.read_csv(args.long)
    strata = sorted(eff.stratum.unique())
    Y = eff.pivot_table(index="gene", columns="stratum", values="log2FoldChange") \
           .reindex(columns=strata)
    S = eff.pivot_table(index="gene", columns="stratum", values="lfcSE") \
           .reindex(columns=strata)
    Yv, Sv = Y.to_numpy(), S.to_numpy()
    ok = np.isfinite(Yv) & np.isfinite(Sv) & (Sv > 0)

    mu0, padj0, keep0 = dl_pool(Yv, Sv, ok, args.min_strata)
    base_sig = set(Y.index[(padj0 < 0.05) & keep0])
    print(f"full cohort: {len(strata)} strata, {int(keep0.sum())} genes pooled, "
          f"{len(base_sig)} significant")

    panel = [g for g in pd.read_csv(args.panel).gene if g in Y.index]
    gi = {g: i for i, g in enumerate(Y.index)}

    rows, tracks = [], {g: [] for g in panel}
    for j, st in enumerate(strata):
        m = ok.copy()
        m[:, j] = False
        mu, padj, keep = dl_pool(Yv, Sv, m, args.min_strata)
        sig = set(Y.index[(padj < 0.05) & keep])
        rows.append(dict(left_out=st, n_genes_pooled=int(keep.sum()), n_sig=len(sig),
                         retained_from_full=len(sig & base_sig),
                         pct_of_full=100 * len(sig & base_sig) / len(base_sig)
                         if base_sig else np.nan))
        for g in panel:
            i = gi[g]
            tracks[g].append((padj[i], mu[i], keep[i]))

    L = pd.DataFrame(rows).sort_values("n_sig")
    L.to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}leave_one_stratum_out_genomewide.csv"),
             index=False)

    prows = []
    for g in panel:
        t = tracks[g]
        padjs = np.array([x[0] for x in t], dtype=float)
        mus = np.array([x[1] for x in t], dtype=float)
        testable = np.array([bool(x[2]) for x in t])
        prows.append(dict(
            gene=g,
            padj_full=float(padj0[gi[g]]), mu_full=float(mu0[gi[g]]),
            n_refits=len(t), n_testable=int(testable.sum()),
            n_sig_refits=int(np.nansum(padjs[testable] < 0.05)) if testable.any() else 0,
            min_padj=float(np.nanmin(padjs)) if testable.any() else np.nan,
            max_padj=float(np.nanmax(padjs)) if testable.any() else np.nan,
            mu_min=float(np.nanmin(mus)) if testable.any() else np.nan,
            mu_max=float(np.nanmax(mus)) if testable.any() else np.nan,
            most_influential_stratum=strata[int(np.nanargmax(padjs))]
            if testable.any() else ""))
    P = pd.DataFrame(prows).sort_values("n_sig_refits", ascending=False)
    P.to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}panel_leave_one_out.csv"), index=False)

    print("\ngenome-wide significant-gene count by dropped stratum:")
    print(L.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    print(f"\nrange across refits: {L.n_sig.min()} to {L.n_sig.max()} "
          f"against {len(base_sig)} for the full cohort")
    print("\npanel stability (refits significant / refits testable):")
    print(P[["gene", "mu_full", "padj_full", "n_sig_refits", "n_testable",
             "most_influential_stratum"]]
          .to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
