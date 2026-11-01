#!/usr/bin/env python3
"""07_metaregression_duration.py - random-effects meta-regression of per-stratum
drought log2 fold-changes on log10(stress duration in days).

tau^2 by DerSimonian-Laird with the Knapp-Hartung small-sample correction, which
matters here: with ~10 strata the uncorrected Wald standard errors are
anti-conservative (IntHout et al. 2014, doi:10.1186/1471-2288-14-25). Genes are
fitted only where every included stratum has an estimate, so the design matrix is
identical across genes.

Durations come from config/duration_curation.csv, which records per stratum where
the number came from and, where it is missing, why it could not be established.
Strata without a duration are dropped rather than imputed.

  python 07_metaregression_duration.py --long meta/per_study_effects_long.csv.gz \
      --durations config/duration_curation.csv --out metareg_duration.csv
"""
import argparse
import numpy as np, pandas as pd
from pymare.estimators import DerSimonianLaird
from statsmodels.stats.multitest import multipletests


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--long", required=True)
    ap.add_argument("--durations", required=True)
    ap.add_argument("--out", default="metareg_duration.csv")
    ap.add_argument("--stress-method", default=None,
                    help="restrict to one stress_method level, e.g. water_withholding, to keep "
                         "duration from standing in for assay type")
    ap.add_argument("--correction", default="knapp-hartung",
                    choices=["knapp-hartung", "knapp-hartung-conservative", "wald"])
    args = ap.parse_args()

    eff = pd.read_csv(args.long)
    dur = pd.read_csv(args.durations).dropna(subset=["duration_days"])
    dur = dur[dur.stratum.isin(set(eff.stratum))]
    if args.stress_method:
        dur = dur[dur.stress_method == args.stress_method]
    strata = list(dur.stratum)
    if len(strata) < 5:
        raise SystemExit(f"only {len(strata)} strata carry a duration; meta-regression not attempted")

    s = eff[eff.stratum.isin(strata)]
    y = s.pivot_table(index="gene", columns="stratum", values="log2FoldChange")
    se = s.pivot_table(index="gene", columns="stratum", values="lfcSE")
    keep = y.notna().all(axis=1) & se.notna().all(axis=1) & (se > 0).all(axis=1)
    y, se = y[keep][strata], se[keep][strata]
    d = dur.set_index("stratum").loc[strata, "duration_days"].to_numpy()
    X = np.c_[np.ones(len(strata)), np.log10(d)]

    r = DerSimonianLaird(small_sample_correction=args.correction).fit(
        y=y.to_numpy().T, v=(se.to_numpy().T) ** 2, X=X).summary()
    fe = r.get_fe_stats()
    out = pd.DataFrame({
        "gene": y.index,
        "intercept_at_1day": np.asarray(fe["est"])[0],
        "slope_per_10x_duration": np.asarray(fe["est"])[1],
        "slope_se": np.asarray(fe["se"])[1],
        "p": np.asarray(fe["p"])[1],
        "tau2": np.asarray(r.tau2).ravel()})
    out["padj"] = multipletests(out.p, method="fdr_bh")[1]
    out.to_csv(args.out, index=False)
    print(f"{len(strata)} strata ({d.min():.3g}-{d.max():.3g} d), {len(out)} genes, "
          f"{(out.padj < 0.05).sum()} with FDR<0.05 duration slope "
          f"[correction: {args.correction}]")


if __name__ == "__main__":
    main()
