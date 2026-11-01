#!/usr/bin/env python
"""Do the pathway genes move TOGETHER across strata, or merely in the same direction?

The competitive test in 11_geneset_test.py asks whether a set's MEAN pooled effect is
displaced. A set can pass that test without being coordinated: each member could be
induced in a different subset of studies. Coordination is a different claim -- that the
members rise and fall together from stratum to stratum -- and it is the claim implied by
calling a pathway "co-regulated".

Statistic: the mean pairwise Pearson correlation of per-stratum log2 fold changes among
the members of a set. Pairs are computed with pairwise-complete deletion and require
`--min-overlap` shared strata, which keeps genes that the meta-analysis could test but
that are not measured in every stratum.

The null is built exactly as in 11_geneset_test.py: each target gene is matched to its
`--neighbours` nearest genes on log10(median max count + 1), one is drawn at random, and
the same statistic is computed on the drawn set. This matters more here than for the mean
test, because stratum-level technical structure correlates ALL genes to some degree; the
matched null carries that same structure, so an excess over it is coordination specific to
the set. Panel members are excluded from the sampling pool.

Coverage is not matched explicitly -- the matching variable is expression, as in the mean
test -- so the output carries the mean number of strata per member for the set and for the
null, as a diagnostic the reader can check.

Multiplicity follows the same family convention as the mean test.

Usage:
  python 16_coordination.py --long meta/per_study_effects_long.csv.gz \
      --meta meta/PUBLIC_meta_drought_random_effects.csv \
      --sets config/geneset_sets.csv --panel config/rfo_panel_full.csv --out meta
"""
import argparse
import os

import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from intergene import mc_pvalue, pairwise_r  # noqa: E402

PUBLIC_PREFIX = "PUBLIC_"


def cross_r(Z, M, ia, ib, min_overlap):
    """Mean Pearson r over every pair with one member in `ia` and one in `ib`."""
    za, ma = Z[ia], M[ia]
    zb, mb = Z[ib], M[ib]
    n = ma @ mb.T
    sx = za @ mb.T
    sy = (zb @ ma.T).T
    sxy = za @ zb.T
    sxx = (za * za) @ mb.T
    syy = ((zb * zb) @ ma.T).T
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (n * sxy - sx * sy) / np.sqrt((n * sxx - sx ** 2) * (n * syy - sy ** 2))
    ok = (n >= min_overlap) & np.isfinite(r)
    # a gene shared by both sets would contribute a spurious r = 1
    shared = np.isin(ia[:, None], ib[None, :])
    ok &= ~shared
    if not ok.any():
        return np.nan, 0
    return float(r[ok].mean()), int(ok.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--long", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--sets", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--neighbours", type=int, default=300)
    ap.add_argument("--min-overlap", type=int, default=8)
    ap.add_argument("--pairs", default=None,
                    help="declared set pairs for the between-set coordination test")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    long = pd.read_csv(args.long)
    W = long.pivot_table(index="gene", columns="stratum", values="log2FoldChange")
    M_meta = pd.read_csv(args.meta).dropna(subset=["beta_random", "median_max_count"])
    # the gene universe is exactly the one the meta-analysis could test
    W = W.reindex(M_meta.gene.to_numpy()).dropna(how="all")
    M_meta = M_meta.set_index("gene").loc[W.index]
    genes = W.index.to_numpy()
    X = W.to_numpy(dtype=float)
    Mask = np.isfinite(X).astype(float)
    Z = np.where(np.isfinite(X), X, 0.0)
    logexp = np.log10(M_meta.median_max_count.to_numpy() + 1)
    n_strata = Mask.sum(1)
    pos = {g: i for i, g in enumerate(genes)}

    panel = set(pd.read_csv(args.panel).gene)
    pool = np.array([i for i, g in enumerate(genes) if g not in panel])
    pool_exp = logexp[pool]

    sets = pd.read_csv(args.sets)
    rows = []
    for set_name, sg in sets.groupby("set", sort=False):
        idx = np.array([pos[g] for g in sg.gene if g in pos])
        if len(idx) < 2:
            print(f"{set_name}: fewer than two members testable, skipped")
            continue
        alternative = str(sg.alternative.iloc[0])
        n_draws = int(sg.n_draws.iloc[0])
        seed = int(sg.seed.iloc[0])

        obs, n_pairs = pairwise_r(Z, Mask, idx, args.min_overlap)

        cand = [pool[np.abs(pool_exp - logexp[i]).argsort()[:args.neighbours]] for i in idx]
        rng = np.random.default_rng(seed)
        null = np.empty(n_draws)
        null_cov = np.empty(n_draws)
        for d in range(n_draws):
            pick = np.array([c[rng.integers(len(c))] for c in cand])
            null[d], _ = pairwise_r(Z, Mask, pick, args.min_overlap)
            null_cov[d] = n_strata[pick].mean()
        good = np.isfinite(null)
        null = null[good]
        centre = float(null.mean())
        sd = float(null.std())
        if alternative == "greater":
            n_extreme = int((null >= obs).sum())
        elif alternative == "less":
            n_extreme = int((null <= obs).sum())
        else:
            n_extreme = int((np.abs(null - centre) >= abs(obs - centre)).sum())

        rows.append(dict(
            gene_set=set_name, n_members_declared=len(sg), n_members_tested=len(idx),
            n_pairs=n_pairs, mean_pairwise_r=obs, null_mean=centre, null_sd=sd,
            null_pct2_5=float(np.percentile(null, 2.5)),
            null_pct97_5=float(np.percentile(null, 97.5)),
            sd_above_null=(obs - centre) / sd if sd > 0 else np.nan,
            sd_above_null_3x_inflated=(obs - centre) / (3 * sd) if sd > 0 else np.nan,
            mean_n_strata_set=float(n_strata[idx].mean()),
            mean_n_strata_null=float(np.nanmean(null_cov)),
            alternative=alternative, n_draws=int(good.sum()), seed=seed,
            n_draws_as_extreme=n_extreme, p_empirical=mc_pvalue(n_extreme, int(good.sum())),
            min_overlap=args.min_overlap, neighbours=args.neighbours,
            family=str(sg.family.iloc[0]) if "family" in sg else set_name,
            description=str(sg.description.iloc[0]) if "description" in sg else ""))

    R = pd.DataFrame(rows)
    R["p_for_bh"] = R.p_empirical
    R["p_bh_within_family"] = np.nan
    for fam, ix in R.groupby("family").groups.items():
        p = R.loc[ix, "p_for_bh"].to_numpy()
        o = np.argsort(p)
        n = len(p)
        q = np.empty(n)
        q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
        R.loc[ix, "p_bh_within_family"] = np.minimum(q, 1.0)
    R.to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}coordination_test.csv"), index=False)

    # the panel correlation matrix itself, for the figure and for inspection
    pidx = np.array([pos[g] for g in pd.read_csv(args.panel).gene if g in pos])
    pn = genes[pidx]
    C = np.full((len(pidx), len(pidx)), np.nan)
    for a in range(len(pidx)):
        for b in range(len(pidx)):
            if a == b:
                C[a, b] = 1.0
            elif a < b:
                r, k = pairwise_r(Z, Mask, np.array([pidx[a], pidx[b]]), args.min_overlap)
                C[a, b] = C[b, a] = r
    pd.DataFrame(C, index=pn, columns=pn).to_csv(
        os.path.join(args.out, f"{PUBLIC_PREFIX}panel_correlation_matrix.csv"))


    # ---- between-set coordination, for declared pairs only ----
    if args.pairs:
        members = {k: np.array([pos[g] for g in sg.gene if g in pos])
                   for k, sg in sets.groupby("set", sort=False)}
        crows = []
        for _, pr in pd.read_csv(args.pairs).iterrows():
            ia, ib = members[pr.set_a], members[pr.set_b]
            obs, n_pairs = cross_r(Z, Mask, ia, ib, args.min_overlap)
            ca = [pool[np.abs(pool_exp - logexp[i]).argsort()[:args.neighbours]] for i in ia]
            cb = [pool[np.abs(pool_exp - logexp[i]).argsort()[:args.neighbours]] for i in ib]
            rng = np.random.default_rng(int(pr.seed))
            nd = int(pr.n_draws)
            null = np.empty(nd)
            for d in range(nd):
                pa = np.array([c[rng.integers(len(c))] for c in ca])
                pb = np.array([c[rng.integers(len(c))] for c in cb])
                null[d], _ = cross_r(Z, Mask, pa, pb, args.min_overlap)
            null = null[np.isfinite(null)]
            centre, sd = float(null.mean()), float(null.std())
            if pr.alternative == "greater":
                ne = int((null >= obs).sum())
            elif pr.alternative == "less":
                ne = int((null <= obs).sum())
            else:
                ne = int((np.abs(null - centre) >= abs(obs - centre)).sum())
            crows.append(dict(
                set_a=pr.set_a, set_b=pr.set_b, n_pairs=n_pairs, mean_pairwise_r=obs,
                null_mean=centre, null_sd=sd,
                sd_above_null=(obs - centre) / sd if sd > 0 else np.nan,
                sd_above_null_3x_inflated=(obs - centre) / (3 * sd) if sd > 0 else np.nan,
                alternative=pr.alternative, n_draws=len(null), seed=int(pr.seed),
                n_draws_as_extreme=ne, p_empirical=mc_pvalue(ne, len(null)),
                family=pr.family, description=pr.description))
        X2 = pd.DataFrame(crows)
        X2["p_for_bh"] = X2.p_empirical
        X2["p_bh_within_family"] = np.nan
        for fam, ix in X2.groupby("family").groups.items():
            p = X2.loc[ix, "p_for_bh"].to_numpy()
            o = np.argsort(p)
            n = len(p)
            q = np.empty(n)
            q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
            X2.loc[ix, "p_bh_within_family"] = np.minimum(q, 1.0)
        X2.to_csv(os.path.join(args.out, f"{PUBLIC_PREFIX}coordination_between_sets.csv"),
                  index=False)
        print()
        for _, r in X2.iterrows():
            print(f"{r.set_a} x {r.set_b} ({r.n_pairs} pairs): r = {r.mean_pairwise_r:+.3f}, "
                  f"null {r.null_mean:+.3f} sd {r.null_sd:.3f}, "
                  f"p = {r.p_empirical:.4f}, BH {r.p_bh_within_family:.4f}, "
                  f"{r.sd_above_null:+.1f} sd")

    for _, r in R.iterrows():
        print(f"{r.gene_set} ({r.n_members_tested} members, {r.n_pairs} pairs)")
        print(f"  mean pairwise r = {r.mean_pairwise_r:+.3f}")
        print(f"  expression-matched null: mean {r.null_mean:+.3f}, sd {r.null_sd:.3f}")
        print(f"  {r.alternative} empirical p = {r.p_empirical:.4f}, "
              f"BH within family {r.p_bh_within_family:.4f}")
        print(f"  {r.sd_above_null:+.1f} null sd ({r.sd_above_null_3x_inflated:+.1f} if "
              f"the null spread is inflated threefold); coverage "
              f"{r.mean_n_strata_set:.1f} vs null {r.mean_n_strata_null:.1f} strata")


if __name__ == "__main__":
    main()
