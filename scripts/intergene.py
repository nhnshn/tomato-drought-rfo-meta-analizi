#!/usr/bin/env python
"""Inter-gene correlation shared by the gene-set and coordination steps.

A competitive gene-set test that draws its null genes independently assumes the
target set is independent too. Pathway genes are not: on this corpus the
biosynthetic arm's members co-vary across strata at mean pairwise r = 0.55, while
a matched random draw sits at ~0. The null is therefore too narrow and the p-value
too small -- the inter-gene correlation problem Wu & Smyth (2012, Nucleic Acids
Research, doi:10.1093/nar/gks461) identify for competitive tests.

The correction used here is theirs: inflate the variance of the set mean by

    VIF = 1 + (m - 1) * rbar

for a set of m tested members with mean pairwise correlation rbar, which is the
exact variance ratio for an equicorrelated set. Applied to an empirical null it
rescales the draws about their centre by sqrt(VIF), widening the null without
assuming it is normal.

rbar is estimated from the same per-stratum effect matrix the meta-analysis pools,
with pairwise-complete deletion and a minimum overlap; strata play the role
samples play in CAMERA. It is an estimate from at most a handful of points per
pair, so VIF is itself uncertain and the uncorrected statistic is reported beside
the corrected one rather than replaced by it.
"""
import numpy as np


def effect_matrix(long_df, genes_index, value="log2FoldChange"):
    """Gene x stratum matrix of per-stratum effects, as (Z filled with 0, presence mask)."""
    W = long_df.pivot_table(index="gene", columns="stratum", values=value)
    W = W.reindex(genes_index)
    X = W.to_numpy(dtype=float)
    mask = np.isfinite(X).astype(float)
    return np.where(np.isfinite(X), X, 0.0), mask


def pairwise_r(Z, M, idx, min_overlap):
    """Mean pairwise Pearson r among rows `idx`, pairwise-complete.

    Filling missing entries with 0 makes each cross-product sum over exactly the
    strata both genes share, so the moments below are the pairwise-complete ones.
    Returns (mean r, number of pairs meeting `min_overlap`).
    """
    z, m = Z[idx], M[idx]
    n = m @ m.T
    sx = z @ m.T
    sxy = z @ z.T
    sxx = (z * z) @ m.T
    with np.errstate(invalid="ignore", divide="ignore"):
        num = n * sxy - sx * sx.T
        den = np.sqrt((n * sxx - sx ** 2) * (n * sxx.T - sx.T ** 2))
        r = num / den
    iu = np.triu_indices(len(idx), 1)
    rv, nv = r[iu], n[iu]
    ok = (nv >= min_overlap) & np.isfinite(rv)
    if not ok.any():
        return np.nan, 0
    return float(rv[ok].mean()), int(ok.sum())


def vif(rbar, m):
    """Variance inflation factor for the mean of m equicorrelated members.

    Floored at 1: a negative rbar would otherwise shrink the null and manufacture
    significance, which is not a correction we are willing to apply.
    """
    if not np.isfinite(rbar) or m < 2:
        return 1.0
    return float(max(1.0, 1.0 + (m - 1) * rbar))


def mc_pvalue(n_extreme, n_draws):
    """Monte Carlo p with the plus-one correction.

    An empirical p of exactly 0 claims a resolution the simulation does not have;
    (b + 1) / (B + 1) is the standard unbiased form (Phipson & Smyth 2010) and is
    what gets reported and corrected.
    """
    return (n_extreme + 1) / (n_draws + 1)
