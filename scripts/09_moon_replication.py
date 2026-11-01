#!/usr/bin/env python3
"""09_moon_replication.py - replication test against Moon et al. 2026
(Appl Biol Chem, doi:10.1186/s13765-026-01137-7).

Two comparisons at different resolutions, because they answer different questions.

1. The ten RT-qPCR candidates of their Table S6, each with a stated pooled
   direction. Ten named genes is a weak test -- it can only ever give a sign test
   over ten trials -- but it is the comparison they themselves foreground.

2. All 5,681 of their heterogeneity-filtered meta-DEGs (Table S3). Correlating two
   independent pooled effect estimates over thousands of genes is a far stronger
   check of whether two meta-analyses of the same literature agree, and it does not
   depend on anyone's significance threshold.

Note on independence: their corpus is four BioProjects, of which PRJEB83804 is also
in ours. The overlap is one stratum of fourteen, so agreement is close to but not
fully independent, and that is stated with the result rather than glossed.

Replication is declared per gene on criteria fixed before looking at the result:
the pooled random-effects estimate has the reported direction, it is significant at
FDR < 0.05 here, and it was testable in at least --min-strata strata. A gene that
agrees in direction without reaching significance is reported as exactly that --
agreement in direction is not replication, and failure to reach significance in a
cohort with different power is not refutation.

  python 09_moon_replication.py --meta meta/PUBLIC_meta_drought_random_effects.csv \
      --candidates config/moon_candidates.csv --metadegs config/moon_metaDEGs.csv \
      --out meta/moon_replication.csv
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy import stats


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--candidates", required=True,
                    help="their Table S6: gene_base, meta_direction, annotation")
    ap.add_argument("--metadegs", default=None,
                    help="their Table S3: gene_base, random_log2FC_PM, random_fdr_PM, "
                         "high_confidence_core")
    ap.add_argument("--out", default="moon_replication.csv")
    ap.add_argument("--min-strata", type=int, default=8)
    args = ap.parse_args()

    m = pd.read_csv(args.meta).set_index("gene")
    cand = pd.read_csv(args.candidates)
    name_col = ("common_name_or_annotation" if "common_name_or_annotation" in cand.columns
                else cand.columns[1])

    rows = []
    for _, c in cand.iterrows():
        g = c.gene_base
        direction = +1 if str(c.meta_direction).strip().lower() == "up" else -1
        if g not in m.index:
            rows.append(dict(gene=g, name=c[name_col], reported_direction=direction,
                             status="not testable here",
                             note="below the expression floor in too many strata"))
            continue
        r = m.loc[g]
        same = int(np.sign(r.beta_random)) == direction
        rows.append(dict(
            gene=g, name=c[name_col], reported_direction=direction,
            reported_log2FC=c.get("random_meta_log2FC_PM", np.nan),
            beta_random=r.beta_random, padj=r.padj, n_strata=r.n_strata, I2=r.I2,
            prop_same_direction=r.prop_same_direction, direction_matches=same,
            status=("replicated" if same and r.padj < 0.05 and r.n_strata >= args.min_strata
                    else "same direction, not significant" if same
                    else "opposite direction, significant" if r.padj < 0.05
                    else "opposite direction, not significant")))
    out = pd.DataFrame(rows)
    out.to_csv(args.out, index=False)

    t = out[out.direction_matches.notna()]
    if len(t):
        k, n = int(t.direction_matches.sum()), len(t)
        print(f"Table S6 candidates: {n} of {len(cand)} testable | direction agreement "
              f"{k}/{n} (sign test p = {stats.binomtest(k, n, 0.5).pvalue:.3f}) | "
              f"replicated at FDR<0.05: {int((out.status == 'replicated').sum())}")
        print(out[["gene", "name", "reported_direction", "beta_random", "padj", "status"]]
              .to_string(index=False, float_format=lambda x: f"{x:.3f}", max_colwidth=34))

    # ---- the stronger comparison: their whole meta-DEG set -------------------
    if not args.metadegs:
        return
    s3 = pd.read_csv(args.metadegs)
    keep = [c for c in ["gene_base", "random_log2FC_PM", "random_fdr_PM",
                        "high_confidence_core", "I2"] if c in s3.columns]
    s3 = s3[keep].dropna(subset=["gene_base"])
    mine = m.reset_index()[["gene", "beta_random", "padj", "n_strata"]]
    J = s3.merge(mine, left_on="gene_base", right_on="gene", how="inner")
    J = J.dropna(subset=["random_log2FC_PM", "beta_random"])
    if not len(J):
        print("\nno overlap with the meta-DEG table")
        return

    r = stats.pearsonr(J.random_log2FC_PM, J.beta_random)
    rs = stats.spearmanr(J.random_log2FC_PM, J.beta_random)
    agree_n = int((np.sign(J.random_log2FC_PM) == np.sign(J.beta_random)).sum())
    print(f"\ntheir heterogeneity-filtered meta-DEGs: {len(s3)} | testable here: {len(J)}")
    print(f"effect-size correlation: Pearson r = {r[0]:+.3f} (p = {r[1]:.1e}) | "
          f"Spearman rho = {rs[0]:+.3f}")
    print(f"direction agreement over all {len(J)}: {100 * agree_n / len(J):.1f}% "
          f"(binomial p = {stats.binomtest(agree_n, len(J), 0.5).pvalue:.1e})")

    summary = [dict(subset="all testable meta-DEGs", n=len(J), pearson_r=r[0],
                    spearman_rho=rs[0], pct_direction_agreement=100 * agree_n / len(J),
                    our_median_abs_effect=float(J.beta_random.abs().median()),
                    our_n_fdr05=int((J.padj < 0.05).sum()))]
    sig = J[J.random_fdr_PM < 0.05] if "random_fdr_PM" in J else J.iloc[:0]
    if len(sig):
        a = int((np.sign(sig.random_log2FC_PM) == np.sign(sig.beta_random)).sum())
        print(f"their FDR<0.05 subset ({len(sig)}): direction agreement "
              f"{100 * a / len(sig):.1f}%, our median |effect| "
              f"{sig.beta_random.abs().median():.2f}")
        summary.append(dict(subset="their FDR<0.05 subset", n=len(sig),
                            pearson_r=float(stats.pearsonr(sig.random_log2FC_PM,
                                                           sig.beta_random)[0]),
                            spearman_rho=float(stats.spearmanr(sig.random_log2FC_PM,
                                                               sig.beta_random)[0]),
                            pct_direction_agreement=100 * a / len(sig),
                            our_median_abs_effect=float(sig.beta_random.abs().median()),
                            our_n_fdr05=int((sig.padj < 0.05).sum())))
    core = J[J.high_confidence_core == True] if "high_confidence_core" in J else J.iloc[:0]
    if len(core):
        a = int((np.sign(core.random_log2FC_PM) == np.sign(core.beta_random)).sum())
        print(f'their "high-confidence core" ({len(core)}): direction agreement '
              f"{100 * a / len(core):.1f}%, our FDR<0.05 {int((core.padj < 0.05).sum())}")
        summary.append(dict(subset="their high-confidence core", n=len(core),
                            pearson_r=float(stats.pearsonr(core.random_log2FC_PM,
                                                           core.beta_random)[0]),
                            spearman_rho=float(stats.spearmanr(core.random_log2FC_PM,
                                                               core.beta_random)[0]),
                            pct_direction_agreement=100 * a / len(core),
                            our_median_abs_effect=float(core.beta_random.abs().median()),
                            our_n_fdr05=int((core.padj < 0.05).sum())))

    base = os.path.splitext(args.out)[0]
    J.to_csv(f"{base}_metaDEG_genes.csv", index=False)
    pd.DataFrame(summary).to_csv(f"{base}_metaDEG_summary.csv", index=False)


if __name__ == "__main__":
    main()
