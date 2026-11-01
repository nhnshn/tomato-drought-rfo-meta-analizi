#!/usr/bin/env python
"""06_genotype_moderator.py — DE within mutant/transgenic arms, then a WT-vs-altered
genotype moderator test restricted to studies contributing BOTH arms."""
import sys, numpy as np, pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import stats
from statsmodels.stats.multitest import multipletests

def de_one(counts, meta, label):
    """counts: genes x samples; meta: index=samples with column 'group' (control/drought)."""
    counts = counts.loc[:, meta.index]
    keep = counts.sum(axis=1) > 0
    maxc = counts.max(axis=1)
    dds = DeseqDataSet(counts=counts[keep].T.astype(int), metadata=meta,
                       design="~group", refit_cooks=True, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["group", "drought", "control"], quiet=True)
    st.summary()
    r = st.results_df.reset_index().rename(columns={"index": "gene"})
    r["stratum"] = label
    r["max_count"] = r.gene.map(maxc)
    return r[["gene", "stratum", "log2FoldChange", "lfcSE", "pvalue", "padj", "max_count"]]

def dl_pool(y, se):
    """DerSimonian-Laird random-effects pooling of one gene."""
    w = 1.0 / se**2
    fe = (w * y).sum() / w.sum()
    Q = (w * (y - fe) ** 2).sum()
    k = len(y)
    C = w.sum() - (w**2).sum() / w.sum()
    tau2 = max(0.0, (Q - (k - 1)) / C) if C > 0 else 0.0
    ws = 1.0 / (se**2 + tau2)
    mu = (ws * y).sum() / ws.sum()
    semu = np.sqrt(1.0 / ws.sum())
    I2 = max(0.0, (Q - (k - 1)) / Q) * 100 if Q > 0 else 0.0
    return mu, semu, tau2, I2, k

def moderator(df, grp_col):
    """Q_between test across two genotype classes for each gene."""
    out = []
    for g, sub in df.groupby("gene"):
        parts = {}
        for lvl, s in sub.groupby(grp_col):
            if len(s) < 2:
                continue
            parts[lvl] = dl_pool(s.log2FoldChange.values, s.lfcSE.values)
        if len(parts) < 2:
            continue
        mus = np.array([p[0] for p in parts.values()])
        ses = np.array([p[1] for p in parts.values()])
        w = 1 / ses**2
        Qb = (w * (mus - (w * mus).sum() / w.sum()) ** 2).sum()
        p = stats.chi2.sf(Qb, len(mus) - 1)
        rec = {"gene": g, "Q_between": Qb, "p_between": p}
        for lvl, pr in parts.items():
            rec[f"pooled_{lvl}"] = pr[0]
            rec[f"se_{lvl}"] = pr[1]
            rec[f"n_strata_{lvl}"] = pr[4]
        out.append(rec)
    o = pd.DataFrame(out)
    o["padj_between"] = multipletests(o.p_between, method="fdr_bh")[1]
    return o.sort_values("p_between")


# ---------------------------------------------------------------- command line

def reliability_corrected_slope(d, min_reliability=0.30):
    """Fraction of the wild-type drought response retained in one altered-genotype arm.

    Regressing the altered-arm effects on the wild-type effects is biased downwards
    because the wild-type effects carry measurement error, so the slope is divided by
    the reliability ratio (true effect variance / observed effect variance) estimated
    from the wild-type standard errors. Below `min_reliability` the wild-type effects
    are mostly noise and the corrected slope is division by almost nothing, so it is
    returned as NaN rather than as a number.
    """
    import statsmodels.api as sm
    x, y, sx = d.lfc_wt.to_numpy(), d.log2FoldChange.to_numpy(), d.se_wt.to_numpy()
    slope_raw = sm.OLS(y, sm.add_constant(x)).fit().params[1]
    var_obs = x.var(ddof=1)
    var_err = float(np.mean(sx ** 2))
    reliability = (var_obs - var_err) / var_obs if var_obs > 0 else float("nan")
    ok = reliability >= min_reliability
    return dict(n_genes=len(d), slope_raw=float(slope_raw), reliability=float(reliability),
                slope_corrected=float(slope_raw / reliability) if ok else np.nan,
                status="estimable" if ok else
                       f"not estimable: wild-type effect variance is dominated by measurement "
                       f"error (reliability {reliability:.2f} < {min_reliability})")


def pathway_offset(d, panel_genes):
    """Does the panel sit off the arm's own attenuation line, beyond overall damping?"""
    import statsmodels.api as sm
    X = sm.add_constant(d.lfc_wt.to_numpy())
    fit = sm.WLS(d.delta.to_numpy(), X, weights=1 / d.se_delta.to_numpy() ** 2).fit()
    resid = d.delta - fit.predict(X)
    z = resid / resid.std()
    inp = z[d.gene.isin(panel_genes)]
    rest = z[~d.gene.isin(panel_genes)]
    if len(inp) < 5:
        return dict(n_panel=len(inp), median_z_panel=np.nan, p_panel_vs_genome=np.nan)
    u = stats.mannwhitneyu(inp, rest, alternative="two-sided")
    return dict(n_panel=int(len(inp)), median_z_panel=float(inp.median()),
                p_panel_vs_genome=float(u.pvalue))


def main():
    import argparse, glob, os
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm-effects", required=True,
                    help="per-arm DE table for the altered genotypes (gene, stratum, study, "
                         "tissue, genotype, log2FoldChange, lfcSE, max_count)")
    ap.add_argument("--wildtype-de", required=True,
                    help="directory of per-stratum wild-type DE tables from 03_de_design_aware.py")
    ap.add_argument("--panel", required=True, help="pathway gene panel CSV with a 'gene' column")
    ap.add_argument("--out", default="genotype_results")
    ap.add_argument("--min-reliability", type=float, default=0.30)
    ap.add_argument("--min-wt-count", type=int, default=20,
                    help="wild-type expression floor; genes below it carry no usable effect")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    arms = pd.read_csv(args.arm_effects)
    arms["wt_stratum"] = arms.study + "__" + arms.tissue.str.replace(r"[^0-9A-Za-z]+", "", regex=True)
    wt = []
    for st in sorted(arms.wt_stratum.unique()):
        f = os.path.join(args.wildtype_de, f"PUBLIC_{st}_droughtVScontrol.csv")
        if not os.path.exists(f):
            raise SystemExit(f"no wild-type table for {st}: {f}")
        t = pd.read_csv(f)[["gene", "log2FoldChange", "lfcSE", "max_count"]]
        t.columns = ["gene", "lfc_wt", "se_wt", "maxc_wt"]
        t["wt_stratum"] = st
        wt.append(t)
    d = arms.merge(pd.concat(wt, ignore_index=True), on=["gene", "wt_stratum"])
    d["delta"] = d.log2FoldChange - d.lfc_wt
    d["se_delta"] = np.sqrt(d.lfcSE ** 2 + d.se_wt ** 2)
    d["p_delta"] = 2 * stats.norm.sf(np.abs(d.delta / d.se_delta))
    panel = set(pd.read_csv(args.panel).gene)

    rows = []
    for arm, s in d.groupby("stratum"):
        s = s[(s.maxc_wt >= args.min_wt_count) & s.lfc_wt.notna()
              & s.log2FoldChange.notna() & (s.se_wt > 0)].copy()
        s["padj_delta"] = multipletests(s.p_delta, method="fdr_bh")[1]
        s.to_csv(os.path.join(args.out, f"PUBLIC_delta_{arm}.csv"), index=False)
        rec = dict(arm=arm)
        rec.update(reliability_corrected_slope(s, args.min_reliability))
        if rec["status"] == "estimable":
            rec.update(pathway_offset(s, panel))
        rows.append(rec)
        print(f"{arm}: slope_corrected="
              f"{rec['slope_corrected'] if rec['slope_corrected'] == rec['slope_corrected'] else 'NA'}"
              f" ({rec['status'].split(':')[0]})", flush=True)
    out = pd.DataFrame(rows).sort_values("slope_corrected", na_position="last")
    out.to_csv(os.path.join(args.out, "PUBLIC_genotype_response_preservation.csv"), index=False)
    print(f"\n{(out.status == 'estimable').sum()} of {len(out)} arms estimable")


if __name__ == "__main__":
    main()
