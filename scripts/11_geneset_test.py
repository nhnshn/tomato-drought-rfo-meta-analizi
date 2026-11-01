#!/usr/bin/env python
"""Competitive gene-set test against an expression-matched background.

Why this test carries the pathway-level claim. Per-gene inference on this corpus is
fragile: under the HKSJ variance no single gene in the panel clears FDR, and the
genome-wide significant-gene count swings by an order of magnitude when one stratum
is dropped. A set-level test asks a different and more robust question -- is the
pathway as a whole displaced relative to comparable genes -- and it answers it with
ONE test, so there is no multiplicity to spend.

"Competitive" means the comparison is against other genes, not against zero. The
null is built by resampling genes matched on expression level, because pooled effect
sizes are not expression-independent: lowly expressed genes carry larger standard
errors and wider pooled effects, so an unmatched null would flatter any set that
happens to be weakly expressed. Each target gene is matched to its
`--neighbours` nearest genes on log10(median max count + 1) and one is drawn at
random; repeating that `n_draws` times gives the distribution of the set mean under
"genes like these, but unrelated".

Panel members are excluded from the sampling pool, so a pathway cannot be tested
against itself.

The test is one-sided for a set with a directional hypothesis and two-sided for a
control set, as declared per set in the sets file -- not chosen after seeing the
result.

Caveat recorded with the result: genes in a pathway are co-expressed, so the draws
are more independent than the target set is. That inflates significance. The
reported effect is therefore also expressed in null standard deviations, where a
threefold inflation of the null spread can be applied by eye.

Usage:
  python 11_geneset_test.py --meta meta/PUBLIC_meta_drought_random_effects.csv \
      --sets config/geneset_sets.csv --panel config/rfo_panel_full.csv --out meta
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from intergene import effect_matrix, mc_pvalue, pairwise_r, vif  # noqa: E402

PUBLIC_PREFIX = "PUBLIC_"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meta", required=True,
                    help="PUBLIC_meta_drought_random_effects.csv from 04_meta_analysis.py")
    ap.add_argument("--sets", required=True,
                    help="CSV: set,gene,alternative,n_draws,seed[,description]")
    ap.add_argument("--panel", required=True,
                    help="pathway panel to exclude from the sampling pool")
    ap.add_argument("--long", required=True,
                    help="per-stratum effects table, for the inter-gene correlation estimate")
    ap.add_argument("--min-overlap", type=int, default=8,
                    help="strata a gene pair must share before its correlation is used")
    ap.add_argument("--out", required=True)
    ap.add_argument("--neighbours", type=int, default=300,
                    help="size of the expression-matched candidate pool per target gene")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    long = pd.read_csv(args.long)
    M = pd.read_csv(args.meta)
    M = M.dropna(subset=["beta_random", "median_max_count"])
    M["logexp"] = np.log10(M.median_max_count + 1)
    panel = set(pd.read_csv(args.panel).gene)
    sets = pd.read_csv(args.sets)

    # the effect matrix the inter-gene correlation is estimated from
    Z, Mask = effect_matrix(long, M.gene.to_numpy())
    rowpos = {gn: i for i, gn in enumerate(M.gene)}

    pool_df = M[~M.gene.isin(panel)]
    pool_exp = pool_df.logexp.to_numpy()
    pool_beta = pool_df.beta_random.to_numpy()

    results, members = [], []
    for set_name, sg in sets.groupby("set", sort=False):
        genes = [g for g in sg.gene if g in set(M.gene)]
        if not genes:
            print(f"{set_name}: no member testable, skipped")
            continue
        alternative = str(sg.alternative.iloc[0])
        n_draws = int(sg.n_draws.iloc[0])
        seed = int(sg.seed.iloc[0])

        tgt = M[M.gene.isin(genes)]
        obs = float(tgt.beta_random.mean())

        # Matching is deterministic, so the candidate pools are built once; the draw
        # order below is identical to drawing inside the loop.
        cand = []
        for le in tgt.logexp:
            idx = np.abs(pool_exp - le).argsort()[:args.neighbours]
            cand.append(pool_beta[idx])

        rng = np.random.default_rng(seed)
        null = np.empty(n_draws)
        for d in range(n_draws):
            null[d] = np.mean([c[rng.integers(len(c))] for c in cand])

        centre = null.mean()
        if alternative == "greater":
            n_extreme = int((null >= obs).sum())
        elif alternative == "less":
            n_extreme = int((null <= obs).sum())
        else:
            n_extreme = int((np.abs(null - centre) >= abs(obs - centre)).sum())
        sd = float(null.std())

        # Correlation-aware null: the draws are independent genes, the set is not.
        # Rescaling the draws about their centre by sqrt(VIF) widens the null to the
        # variance an equicorrelated set of this size would actually have.
        idx = np.array([rowpos[gn] for gn in genes])
        rbar, n_pairs = pairwise_r(Z, Mask, idx, args.min_overlap)
        v = vif(rbar, len(genes))
        null_c = centre + (null - centre) * np.sqrt(v)
        if alternative == "greater":
            n_extreme_c = int((null_c >= obs).sum())
        elif alternative == "less":
            n_extreme_c = int((null_c <= obs).sum())
        else:
            n_extreme_c = int((np.abs(null_c - centre) >= abs(obs - centre)).sum())
        p_emp = mc_pvalue(n_extreme, n_draws)

        results.append(dict(
            gene_set=set_name, n_members_declared=len(sg), n_members_tested=len(genes),
            mean_pooled_log2FC=obs, null_mean=centre, null_sd=sd,
            null_pct2_5=float(np.percentile(null, 2.5)),
            null_pct97_5=float(np.percentile(null, 97.5)),
            sd_above_null=(obs - centre) / sd if sd > 0 else np.nan,
            sd_above_null_3x_inflated=(obs - centre) / (3 * sd) if sd > 0 else np.nan,
            alternative=alternative, n_draws=n_draws, seed=seed,
            n_draws_as_extreme=n_extreme, p_empirical=p_emp,
            mean_intergene_r=rbar, n_intergene_pairs=n_pairs, vif=v,
            sd_above_null_vif=(obs - centre) / (sd * np.sqrt(v)) if sd > 0 else np.nan,
            n_draws_as_extreme_vif=n_extreme_c,
            p_empirical_vif=mc_pvalue(n_extreme_c, n_draws),
            genome_median_pooled=float(M.beta_random.median()),
            neighbours=args.neighbours,
            family=str(sg.family.iloc[0]) if "family" in sg else set_name,
            description=str(sg.description.iloc[0]) if "description" in sg else ""))
        for g in genes:
            r = M[M.gene == g].iloc[0]
            members.append(dict(gene_set=set_name, gene=g,
                                beta_random=float(r.beta_random),
                                median_max_count=float(r.median_max_count),
                                n_strata=int(r.n_strata) if "n_strata" in r else np.nan))

    R = pd.DataFrame(results)
    # Multiplicity is spent only within a declared family of sets. The pre-specified
    # pathway test and the negative control are each their own family (one test, no
    # correction); the sub-arm decomposition is one family of five, so those five are
    # BH-corrected against each other. A set whose family holds a single test keeps
    # its empirical p. An empirical p of exactly 0 means no draw was as extreme and is
    # bounded below by 1/n_draws; the bound, not 0, is what is carried into BH.
    # the correlation-aware p is the one corrected and reported as primary
    R["p_for_bh"] = R.p_empirical_vif
    R["p_bh_within_family"] = np.nan
    for fam, idx in R.groupby("family").groups.items():
        p = R.loc[idx, "p_for_bh"].to_numpy()
        o = np.argsort(p)
        n = len(p)
        q = np.empty(n)
        q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
        R.loc[idx, "p_bh_within_family"] = np.minimum(q, 1.0)
    R.to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}geneset_test.csv"), index=False)
    pd.DataFrame(members).to_csv(
        os.path.join(args.out, f"{PUBLIC_PREFIX}geneset_members.csv"), index=False)

    for _, r in R.iterrows():
        print(f"{r.gene_set} ({r.n_members_tested} of {r.n_members_declared} testable)")
        print(f"  mean pooled log2FC = {r.mean_pooled_log2FC:+.3f}")
        print(f"  expression-matched null: mean {r.null_mean:+.3f}, sd {r.null_sd:.3f}, "
              f"97.5th pct {r.null_pct97_5:+.3f}")
        print(f"  independent-gene null: {r.alternative} p = {r.p_empirical:.4f} "
              f"({r.n_draws_as_extreme} of {r.n_draws} draws as extreme), "
              f"{r.sd_above_null:+.1f} null sd")
        print(f"  inter-gene r = {r.mean_intergene_r:+.3f} over {r.n_intergene_pairs} pairs "
              f"-> VIF {r.vif:.2f}, null sd widened to {r.null_sd * np.sqrt(r.vif):.3f}")
        print(f"  CORRELATION-AWARE p = {r.p_empirical_vif:.4f} "
              f"({r.n_draws_as_extreme_vif} of {r.n_draws} widened draws as extreme), "
              f"{r.sd_above_null_vif:+.1f} null sd  [BH {r.p_bh_within_family:.4f}]")
    if len(R):
        print(f"\ngenome-wide median pooled effect: {R.genome_median_pooled.iloc[0]:+.3f}")


if __name__ == "__main__":
    main()
