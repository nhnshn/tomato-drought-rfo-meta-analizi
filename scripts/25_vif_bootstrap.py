#!/usr/bin/env python
"""Partial sensitivity of the variance-inflation correction to the strata it is estimated from.

The pathway p-value is produced by widening an expression-matched null by sqrt(VIF), with
VIF = 1 + (m-1)*rbar. rbar is a point estimate from at most 13 effect values per gene pair,
and the correction it drives moves the p-value by roughly three orders of magnitude, so
reporting the corrected p as a single number states a precision the data do not support.

WHAT THIS RESAMPLES, AND WHAT IT HOLDS FIXED. Strata are resampled with replacement; for
each resample rbar and VIF are recomputed from the resampled effect matrix and the
expression-matched null is redrawn. The observed set mean and the pooled per-gene effects
of the background genes are NOT recomputed -- they are read once from the meta table and
held fixed. The output is therefore a sensitivity distribution over the correlation and
the inflation factor, conditional on the fitted meta-analysis, and NOT a confidence
interval for the pathway p-value. A full bootstrap would refit the per-gene effects,
re-derive set membership and the background, and recompute the test statistic inside every
resample, and would additionally have to treat the strata that share BioSamples as one
cluster rather than as independent draws; that is a different and much larger analysis.

Resampling is over strata rather than genes because strata are the unit the meta-analysis
treats as exchangeable and a gene bootstrap would break the set membership under test.
Note that the two strata sharing eight BioSamples are drawn independently here, and a
stratum drawn twice contributes twice to the minimum-overlap count, both of which make the
spread somewhat optimistic.
"""
import argparse, os, sys
import numpy as np, pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from intergene import mc_pvalue, pairwise_r, vif  # noqa: E402


def effect_matrix_from(long, genes):
    W = long.pivot_table(index="gene", columns="stratum", values="log2FoldChange")
    W = W.reindex(genes)
    X = W.to_numpy(dtype=float)
    return np.where(np.isfinite(X), X, 0.0), np.isfinite(X).astype(float), list(W.columns)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", required=True)
    ap.add_argument("--long", required=True)
    ap.add_argument("--sets", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--set-name", default="rfo_biosynthetic")
    ap.add_argument("--min-overlap", type=int, default=8)
    ap.add_argument("--neighbours", type=int, default=300)
    ap.add_argument("--boot", type=int, default=400)
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    M = pd.read_csv(a.meta).dropna(subset=["beta_random", "median_max_count"])
    M["logexp"] = np.log10(M.median_max_count + 1)
    long = pd.read_csv(a.long)
    panel = set(pd.read_csv(a.panel).gene)
    sets = pd.read_csv(a.sets)
    members = [g for g in sets[sets["set"] == a.set_name].gene if g in set(M.gene)]

    pool_df = M[~M.gene.isin(panel)]
    pool_exp, pool_beta = pool_df.logexp.to_numpy(), pool_df.beta_random.to_numpy()
    tgt = M[M.gene.isin(members)]
    cand = [pool_beta[np.abs(pool_exp - le).argsort()[:a.neighbours]] for le in tgt.logexp]

    Z, Mask, strata = effect_matrix_from(long[long.gene.isin(members)], members)
    idx = np.arange(len(members))
    rng = np.random.default_rng(a.seed)

    rbars, vifs, ps = [], [], []
    for _ in range(a.boot):
        take = rng.integers(0, len(strata), len(strata))
        rb, npair = pairwise_r(Z[:, take], Mask[:, take], idx, a.min_overlap)
        if not np.isfinite(rb):
            continue
        v = vif(rb, len(members))
        null = np.array([np.mean([c[rng.integers(len(c))] for c in cand])
                         for _ in range(a.draws)])
        centre = null.mean()
        obs = float(tgt.beta_random.mean())
        widened = centre + (null - centre) * np.sqrt(v)
        ps.append(mc_pvalue(int((widened >= obs).sum()), a.draws))
        rbars.append(rb); vifs.append(v)

    q = lambda x: (float(np.percentile(x, 2.5)), float(np.percentile(x, 50)),
                   float(np.percentile(x, 97.5)))
    rl, rm, rh = q(rbars); vl, vm, vh = q(vifs); pl, pm, ph = q(ps)
    R = pd.DataFrame([dict(gene_set=a.set_name,
                           quantity="partial sensitivity: correlation and VIF resampled "
                                    "over strata; meta-effects and set mean held fixed; "
                                    "not a confidence interval",
                           n_members=len(members),
                           n_bootstrap=len(ps), n_draws=a.draws, seed=a.seed,
                           rbar_median=rm, rbar_lo=rl, rbar_hi=rh,
                           vif_median=vm, vif_lo=vl, vif_hi=vh,
                           p_median=pm, p_lo=pl, p_hi=ph)])
    R.to_csv(os.path.join(a.out, "PUBLIC_vif_bootstrap.csv"), index=False)
    print(f"{a.set_name}: {len(members)} members, {len(ps)} stratum resamples")
    print("  partial sensitivity (meta-effects held fixed), not a confidence interval")
    print(f"  rbar  median {rm:+.3f}  95% [{rl:+.3f}, {rh:+.3f}]")
    print(f"  VIF   median {vm:.2f}   95% [{vl:.2f}, {vh:.2f}]")
    print(f"  p     median {pm:.4f}  95% [{pl:.4f}, {ph:.4f}]")


if __name__ == "__main__":
    main()
