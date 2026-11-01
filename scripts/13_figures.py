#!/usr/bin/env python
"""Every figure in the manuscript, regenerated from the frozen run's own outputs.

Figures are built here and nowhere else. A figure drawn in an interactive session
drifts from the tables as soon as the pipeline is rerun, and a reader cannot tell
which run a stale figure came from -- so each panel reads the same CSV the text
cites, and counts that appear in a title are computed from the data at render time
rather than typed in. A claim-title that cannot be computed is written as a
description instead.

Usage:
  python 13_figures.py --meta meta --de de --subgroup meta_subgroup_method \
      --genotype genotype_results --panel config/rfo_panel_full.csv \
      --out figures
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GREY = "#8c8c8c"
FOCAL = "#b2182b"        # the pathway under test
BG = "#6a9ec0"           # everything else
POS = "#2166ac"
NEG = "#b2182b"
BASE, ANNOT, TICK = 8, 7, 6


def apply_style():
    plt.rcParams.update({
        "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.size": BASE, "axes.titlesize": BASE, "axes.labelsize": BASE,
        "legend.fontsize": ANNOT, "xtick.labelsize": TICK, "ytick.labelsize": TICK,
        "axes.titlelocation": "left", "axes.titleweight": "regular",
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "legend.frameon": False,
        "font.family": "sans-serif",
    })


def letter(ax, c):
    ax.text(-0.08, 1.06, c, transform=ax.transAxes, fontsize=BASE + 2,
            fontweight="bold", va="bottom", ha="right")


def short(s, n=30):
    s = str(s)
    return s if len(s) <= n else s[:n - 1] + "\u2026"


# --------------------------------------------------------------------------- 1
def fig_cohort(de_dir, out):
    d = pd.read_csv(os.path.join(de_dir, "PUBLIC_de_summary.csv"))
    d = d[d.status == "ok"].copy()
    d["tissue"] = d.stratum.str.split("__").str[-1]
    leaf = d[d.tissue == "leafshoot"].copy()
    order = {"concurrent_control": 0, "baseline_referenced": 1, "unresolved_timing": 2}
    leaf["o"] = leaf.estimand.map(order)
    leaf = leaf.sort_values(["o", "n_drought"], ascending=[True, True])
    colours = {"concurrent_control": FOCAL, "baseline_referenced": BG,
               "unresolved_timing": GREY}
    nice = {"concurrent_control": "concurrent control",
            "baseline_referenced": "baseline-referenced",
            "unresolved_timing": "control timing unresolved"}

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6),
                             gridspec_kw={"width_ratios": [2.1, 1]})
    ax = axes[0]
    y = np.arange(len(leaf))
    ax.barh(y - 0.2, leaf.n_control, height=0.38, color="white",
            edgecolor=[colours[e] for e in leaf.estimand], linewidth=0.9)
    ax.barh(y + 0.2, leaf.n_drought, height=0.38,
            color=[colours[e] for e in leaf.estimand], linewidth=0)
    ax.set_yticks(y)
    ax.set_yticklabels([s.replace("__leafshoot", "") for s in leaf.stratum])
    ax.set_xlabel("libraries (open = control, filled = drought)")
    n_pri = int((leaf.estimand == "concurrent_control").sum())
    ax.set_title(f"{n_pri} of {len(leaf)} leaf strata support a\nconcurrent-control contrast")
    ax.margins(y=0.02)
    letter(ax, "a")
    for e, c in colours.items():
        if (leaf.estimand == e).any():
            ax.scatter([], [], marker="s", s=18, color=c, label=nice[e])
    ax.legend(loc="lower right", handletextpad=0.4)

    ax = axes[1]
    g = d.groupby(["tissue", "estimand"]).size().unstack(fill_value=0)
    g = g.reindex(columns=[k for k in order if k in g.columns], fill_value=0)
    bottom = np.zeros(len(g))
    for e in g.columns:
        ax.bar(np.arange(len(g)), g[e], bottom=bottom, color=colours[e], width=0.62)
        bottom += g[e].to_numpy()
    ax.set_xticks(np.arange(len(g)))
    ax.set_xticklabels([t.replace("leafshoot", "leaf/shoot").replace("floweranther", "flower")
                        for t in g.index], rotation=30, ha="right")
    ax.set_ylabel("strata")
    ax.set_title("Strata by tissue")
    letter(ax, "b")

    fig.tight_layout()
    p = os.path.join(out, "tomato_drought_data_landscape.png")
    fig.savefig(p)
    plt.close(fig)
    return p


# --------------------------------------------------------------------------- 2
def fig_rfo(meta_dir, panel_file, out):
    P = pd.read_csv(os.path.join(meta_dir, "PUBLIC_panel_inference.csv"))
    M = pd.read_csv(os.path.join(meta_dir, "PUBLIC_meta_drought_random_effects.csv"))
    G = pd.read_csv(os.path.join(meta_dir, "PUBLIC_geneset_test.csv"))
    ann = dict(zip(pd.read_csv(panel_file).gene,
                   pd.read_csv(panel_file).get("annotation", "")))
    gw = dict(zip(M.gene, M.padj))
    P = P.sort_values("mu")
    P["padj_gw"] = P.gene.map(gw)
    n_sig = int((P.padj_gw < 0.05).sum())

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 4.0),
                             gridspec_kw={"width_ratios": [1.45, 1]})
    ax = axes[0]
    y = np.arange(len(P))
    lo = P.mu - 1.96 * P.se
    hi = P.mu + 1.96 * P.se
    ax.hlines(y, lo, hi, color=GREY, linewidth=1.0, zorder=1)
    sig = (P.padj_gw < 0.05).to_numpy()
    ax.scatter(P.mu[sig], y[sig], s=26, color=FOCAL, zorder=3, label="FDR < 0.05 genome-wide")
    ax.scatter(P.mu[~sig], y[~sig], s=22, facecolors="white", edgecolors=FOCAL,
               linewidths=0.9, zorder=3, label="not significant")
    ax.axvline(0, color="black", linewidth=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([short(ann.get(g, g), 32) for g in P.gene])
    ax.set_xlabel("pooled log$_2$ fold change (drought vs control)")
    ax.set_title(f"{n_sig} of {len(P)} pathway genes clear\ngenome-wide FDR")
    ax.margins(y=0.02)
    ax.legend(loc="lower right", handletextpad=0.3)
    letter(ax, "a")

    ax = axes[1]
    row = G[G.gene_set == "rfo_biosynthetic"]
    if len(row):
        r = row.iloc[0]
        rng = np.random.default_rng(0)
        # redraw the null's shape from its reported moments only for display
        # The drawn null must be the one the reported p came from. The panel title quotes
        # p_empirical_vif, so the curve has to use the correlation-widened spread, not the
        # raw independent-gene null_sd that the uncorrected test used.
        sd_eff = r.null_sd * np.sqrt(r.vif)
        xs = np.linspace(r.null_mean - 4 * sd_eff,
                         max(r.mean_pooled_log2FC * 1.08, r.null_mean + 4 * sd_eff), 400)
        dens = np.exp(-0.5 * ((xs - r.null_mean) / sd_eff) ** 2)
        ax.fill_between(xs, dens, color=BG, alpha=0.5, linewidth=0)
        ax.plot(xs, dens, color=BG, linewidth=1.0)
        ax.axvline(r.mean_pooled_log2FC, color=FOCAL, linewidth=1.6)
        ax.annotate(f"observed\n{r.mean_pooled_log2FC:+.2f}",
                    xy=(r.mean_pooled_log2FC, 0.78), xytext=(-6, 0),
                    textcoords="offset points", ha="right", va="center",
                    color=FOCAL, fontsize=ANNOT)
        ax.annotate("expression-matched\nrandom sets",
                    xy=(r.null_mean, 0.45), xytext=(10, 14), textcoords="offset points",
                    ha="left", va="center", color=BG, fontsize=ANNOT,
                    arrowprops=dict(arrowstyle="-", color=BG, linewidth=0.7))
        ax.set_yticks([])
        ax.set_xlabel("mean pooled log$_2$ fold change of the set")
        # the correlation-aware p is the reported one; at the simulation's resolution
        # floor it is quoted as a bound rather than rounded to zero
        floor = 1.0 / (r.n_draws + 1)
        p_txt = (f"p $\\leq$ {floor:.0e}" if r.p_empirical_vif <= floor
                 else f"p = {r.p_empirical_vif:.4f}")
        ax.set_title(f"The biosynthetic arm as a set:\n{p_txt} against matched genes, "
                     f"VIF {r.vif:.1f}")
        ax.spines["left"].set_visible(False)
    letter(ax, "b")

    fig.tight_layout()
    p = os.path.join(out, "rfo_pathway_drought.png")
    fig.savefig(p)
    plt.close(fig)
    return p


# --------------------------------------------------------------------------- 3
def fig_gols(meta_dir, out):
    eff = pd.read_csv(os.path.join(meta_dir, "per_study_effects_long.csv.gz"))
    M = pd.read_csv(os.path.join(meta_dir, "PUBLIC_meta_drought_random_effects.csv"))
    gols = {"Solyc01g100830": "SlGolS1", "Solyc02g084980": "SlGolS (AtGolS1 ortholog)",
            "Solyc02g062590": "SlGolS (AtGolS2 ortholog)", "Solyc01g079170": "SlGolS (weakest)"}
    Mi = M.set_index("gene")
    present = [g for g in gols if g in Mi.index]
    n_sig = sum(1 for g in present if float(Mi.loc[g, "padj"]) < 0.05)

    fig, ax = plt.subplots(figsize=(5.4, 3.5))
    for i, g in enumerate(present):
        s = eff[eff.gene == g]
        ax.scatter(s.log2FoldChange, np.full(len(s), i) + np.random.default_rng(i)
                   .uniform(-0.14, 0.14, len(s)), s=14, color=BG, alpha=0.7,
                   linewidths=0, zorder=2)
        mu = float(Mi.loc[g, "beta_random"])
        se = float(Mi.loc[g, "se_random"])
        padj = float(Mi.loc[g, "padj"])
        ax.hlines(i, mu - 1.96 * se, mu + 1.96 * se, color=FOCAL, linewidth=1.2, zorder=3)
        ax.scatter([mu], [i], s=34, color=FOCAL if padj < 0.05 else "white",
                   edgecolors=FOCAL, linewidths=1.0, zorder=4)
        ax.text(mu, i + 0.30, f"FDR {padj:.2f}", color=FOCAL, fontsize=ANNOT,
                ha="center", va="bottom")
    ax.axvline(0, color="black", linewidth=0.7)
    ax.set_yticks(range(len(present)))
    ax.set_yticklabels([f"$\\it{{{gols[g].split(' ')[0]}}}$\n{g}" for g in present])
    ax.set_xlabel("log$_2$ fold change (points = strata, bar = pooled 95% CI)")
    ax.set_title(f"Galactinol synthases are induced in most strata;\n"
                 f"{n_sig} of {len(present)} clears genome-wide FDR")
    # The topmost gene's FDR label sits 0.30 above its row, so a symmetric margin
    # leaves it overlapping the two-line title. Give the top band explicit room.
    ax.set_ylim(-0.55, len(present) - 1 + 0.85)
    fig.tight_layout()
    p = os.path.join(out, "gols_family_drought.png")
    fig.savefig(p)
    plt.close(fig)
    return p


# --------------------------------------------------------------------------- 4
def fig_subgroup(subgroup_dir, meta_dir, panel_file, out):
    f = [x for x in os.listdir(subgroup_dir) if x.startswith("PUBLIC_subgroup_by_")]
    S = pd.read_csv(os.path.join(subgroup_dir, f[0]))
    # list, not set: iteration order must not depend on the process hash seed
    panel = list(dict.fromkeys(pd.read_csv(panel_file).gene))
    lv = [c[len("beta_"):] for c in S.columns if c.startswith("beta_")
          and c not in ("beta_overall",)]
    a, b = lv[0], lv[1]

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5))
    ax = axes[0]
    ax.hist(S.I2_overall.dropna(), bins=40, color=GREY, alpha=0.65, linewidth=0,
            label="ignoring the moderator")
    ax.hist(S.I2_within_median.dropna(), bins=40, color=BG, alpha=0.75, linewidth=0,
            label="within moderator levels")
    ax.set_xlabel("I$^2$ (%)")
    ax.set_ylabel("genes")
    ax.set_title(f"Stress method absorbs heterogeneity:\nmedian I$^2$ "
                 f"{S.I2_overall.median():.0f}% \u2192 {S.I2_within_median.median():.0f}%")
    ax.legend(loc="upper left")
    letter(ax, "a")

    ax = axes[1]
    sub = S[S.gene.isin(panel)].dropna(subset=[f"beta_{a}", f"beta_{b}"])
    ax.scatter(S[f"beta_{a}"], S[f"beta_{b}"], s=3, color=GREY, alpha=0.25,
               linewidths=0, rasterized=True)
    ax.scatter(sub[f"beta_{a}"], sub[f"beta_{b}"], s=26, color=FOCAL, zorder=3,
               label="RFO pathway")
    lim = np.nanpercentile(np.abs(np.r_[S[f"beta_{a}"], S[f"beta_{b}"]]), 99.5)
    ax.plot([-lim, lim], [-lim, lim], color="black", linewidth=0.6, linestyle="--")
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel(f"pooled log$_2$FC, {a.replace('_', ' ')}")
    ax.set_ylabel(f"pooled log$_2$FC, {b.replace('_', ' ')}")
    n_bet = int((S.padj_between < 0.05).sum())
    ax.set_title(f"{n_bet} of {len(S):,} genes differ\nbetween methods (FDR < 0.05)")
    ax.legend(loc="upper left")
    letter(ax, "b")

    fig.tight_layout()
    p = os.path.join(out, "subgroup_heterogeneity.png")
    fig.savefig(p)
    plt.close(fig)
    return p


# --------------------------------------------------------------------------- 5
def fig_duration(meta_dir, panel_file, out):
    D = pd.read_csv(os.path.join(meta_dir, "metareg_duration.csv"))
    W = None
    wf = os.path.join(meta_dir, "metareg_duration_withholding_only.csv")
    if os.path.exists(wf):
        W = pd.read_csv(wf)
        # When every duration-carrying stratum already uses water withholding the
        # restriction removes nothing, and two identical panels would imply a
        # sensitivity analysis that was never separable. Show one, and say why.
        if len(W) == len(D) and np.allclose(W.p.to_numpy(), D.p.to_numpy(),
                                            equal_nan=True):
            W = "identical"
    # list, not set: iteration order must not depend on the process hash seed
    panel = list(dict.fromkeys(pd.read_csv(panel_file).gene))

    if isinstance(W, str):
        fig, ax0 = plt.subplots(figsize=(4.6, 3.4))
        axes = [ax0]
        panels = [(D, "every stratum carrying a curated duration uses water\n"
                      "withholding, so restricting by method changes nothing")]
    else:
        fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
        panels = [(D, "all strata with a curated duration"),
                  (W, "water withholding only")]
    for ax, (tab, lab) in zip(axes, panels):
        if tab is None:
            ax.axis("off")
            continue
        ax.scatter(tab.slope_per_10x_duration, -np.log10(tab.p), s=3, color=GREY,
                   alpha=0.25, linewidths=0, rasterized=True)
        sub = tab[tab.gene.isin(panel)]
        ax.scatter(sub.slope_per_10x_duration, -np.log10(sub.p), s=26, color=FOCAL,
                   zorder=3, label="RFO pathway")
        n_sig = int((tab.padj < 0.05).sum())
        ax.axvline(0, color="black", linewidth=0.6)
        ax.set_xlabel("slope per 10\u00d7 duration (log$_2$FC)")
        ax.set_title(f"{lab}\n{n_sig} of {len(tab):,} genes with an FDR-significant slope")
    axes[0].set_ylabel("\u2212log$_{10}$ p")
    axes[0].legend(loc="upper left")
    if len(axes) > 1:
        letter(axes[0], "a")
        letter(axes[1], "b")
    fig.tight_layout()
    p = os.path.join(out, "duration_metaregression.png")
    fig.savefig(p)
    plt.close(fig)
    return p


# --------------------------------------------------------------------------- 6
def fig_genotype(geno_dir, out):
    G = pd.read_csv(os.path.join(geno_dir, "PUBLIC_genotype_response_preservation.csv"))
    est = G[G.status.astype(str).str.contains("estimable", case=False, na=False)
            | G.slope_corrected.notna()].copy()
    est = est.dropna(subset=["slope_corrected"]).sort_values("slope_corrected")
    ne = G[~G.index.isin(est.index)]

    fig, ax = plt.subplots(figsize=(6.0, 3.4))
    y = np.arange(len(est))
    ax.hlines(y, 0, est.slope_corrected, color=GREY, linewidth=0.9)
    ax.scatter(est.slope_corrected, y, s=30, color=FOCAL, zorder=3)
    ax.axvline(1.0, color="black", linewidth=0.7, linestyle="--")
    ax.text(1.0, len(est) - 0.4, "wild-type response", fontsize=ANNOT, ha="center",
            va="bottom", color="black")
    for i, (_, r) in enumerate(est.iterrows()):
        ax.text(r.slope_corrected + 0.02, i, f"{r.slope_corrected:.2f}", fontsize=ANNOT,
                va="center", ha="left", color=FOCAL)
    ax.set_yticks(y)
    ax.set_yticklabels([short(a, 26) for a in est.arm])
    ax.set_xlabel("response retained, relative to the matched wild-type arm")
    note = f"; {len(ne)} not estimable" if len(ne) else ""
    ax.set_title(f"{int((est.slope_corrected > 0.7).sum())} of {len(est)} altered-genotype arms "
                 f"retain >70%\nof the wild-type drought response{note}")
    ax.set_xlim(0, max(1.12, est.slope_corrected.max() * 1.12))
    ax.margins(y=0.06)
    fig.tight_layout()
    p = os.path.join(out, "genotype_moderator.png")
    fig.savefig(p)
    plt.close(fig)
    return p
# Row order, row labels and tissue labels for Figure 4 (tissue_subarm_rfo.png).
# These three names were referenced by fig_tissue but never defined in the shipped
# tarball, so the figure target aborted with NameError and the figure could not be
# regenerated. Restored from the axis text of the exported figure; the five steps are
# the pathway partition in catalytic order, followed by the non-RFO control.
ARM_ORDER = ["rfo_galactinol_synthase", "rfo_raffinose_stachyose_synthase",
             "rfo_alpha_galactosidase_catabolic", "inositol_1P_synthase",
             "inositol_monophosphatase", "gt8_control"]
ARM_LABEL = {"rfo_galactinol_synthase": "galactinol synthase\n(UDP-Gal + myo-inositol)",
             "rfo_raffinose_stachyose_synthase": "raffinose / stachyose\nsynthase",
             "rfo_alpha_galactosidase_catabolic": "$\\alpha$-galactosidase\n(RFO catabolism)",
             "inositol_1P_synthase": "myo-inositol-1-P\nsynthase",
             "inositol_monophosphatase": "myo-inositol\nmonophosphatase",
             "gt8_control": "GT8 control\n(non-RFO)"}
TISSUE_LABEL = {"leafshoot": "leaf / shoot", "root": "root",
                "fruitseed": "fruit / seed", "floweranther": "flower / anther",
                "other": "guard cell / other"}


def fig_tissue(tissue_dir, out):
    S = pd.read_csv(os.path.join(tissue_dir, "PUBLIC_tissue_geneset_test.csv"))
    P = pd.read_csv(os.path.join(tissue_dir, "PUBLIC_tissue_per_stratum.csv"))
    tissues = [t for t in ["leafshoot", "root"] if t in set(S.tissue)]
    cols = {"leafshoot": POS, "root": "#1b7837"}

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5),
                             gridspec_kw={"width_ratios": [1.25, 1.0]})
    ax = axes[0]
    arms = [a for a in ARM_ORDER if a in set(S.gene_set)]
    h = 0.34
    for i, arm in enumerate(arms):
        y0 = len(arms) - 1 - i
        for j, t in enumerate(tissues):
            r = S[(S.gene_set == arm) & (S.tissue == t)]
            if not len(r):
                continue
            r = r.iloc[0]
            y = y0 + (0.5 - j) * h
            # same rule as Fig 2b: the band must be the widened null the BH-adjusted
            # p_empirical_vif was computed against, or a filled marker can sit outside
            # a band that the test it reports did not actually reject.
            sd_eff = r.null_sd * np.sqrt(r.vif)
            ax.plot([r.null_mean - 1.96 * sd_eff, r.null_mean + 1.96 * sd_eff],
                    [y, y], color="0.75", lw=4, solid_capstyle="butt", zorder=1)
            sig = r.p_bh_within_family < 0.05
            ax.plot([r.mean_pooled_log2FC], [y], "o", ms=5.5, color=cols[t],
                    mfc=cols[t] if sig else "white", mew=1.2, zorder=3)
            ax.annotate(f"{r.mean_pooled_log2FC:+.2f}", (r.mean_pooled_log2FC, y),
                        xytext=(5, 0), textcoords="offset points", fontsize=TICK,
                        va="center", color=cols[t])
    ax.axvline(0, color="0.4", lw=0.8, zorder=0)
    ax.set_yticks(range(len(arms)))
    ax.set_yticklabels([ARM_LABEL[a] for a in reversed(arms)], fontsize=TICK)
    ax.set_xlabel("mean pooled log$_2$ fold change of the set")
    ax.set_ylim(-0.6, len(arms) - 0.4)
    ax.set_title("Which step of the pathway moves", fontsize=BASE)
    hs = [plt.Line2D([], [], marker="o", ls="", color=cols[t], label=TISSUE_LABEL[t])
          for t in tissues]
    hs += [plt.Line2D([], [], marker="o", ls="", mfc="white", mec="0.3", color="0.3",
                      label="not significant"),
           plt.Line2D([], [], color="0.75", lw=4, label="95% of the\nexpression-matched null")]
    ax.legend(handles=hs, loc="lower right", fontsize=TICK - 0.5, handletextpad=0.4,
              borderpad=0.3, labelspacing=0.35)
    letter(ax, "a")

    ax = axes[1]
    order = [t for t in ["leafshoot", "root", "fruitseed", "floweranther", "other"]
             if t in set(P.tissue)]
    rng = np.random.default_rng(0)
    for i, t in enumerate(order):
        y0 = len(order) - 1 - i
        v = P[P.tissue == t].mean_biosynthetic_log2FC.to_numpy()
        jit = rng.uniform(-0.17, 0.17, len(v)) if len(v) > 1 else np.zeros(len(v))
        pooled = bool(P[P.tissue == t].pooled.iloc[0]) or t == "leafshoot"
        ax.scatter(v, y0 + jit, s=22, facecolor=cols.get(t, "0.55"),
                   edgecolor="white", linewidth=0.5, zorder=3,
                   alpha=0.95 if pooled else 0.6)
        ax.plot([np.mean(v)], [y0], "|", ms=16, mew=1.8, color="0.25", zorder=4)
    ax.axvline(0, color="0.4", lw=0.8, zorder=0)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(
        [f"{TISSUE_LABEL[t]}\n({(P.tissue == t).sum()} "
         f"{'stratum' if (P.tissue == t).sum() == 1 else 'strata'})"
         for t in reversed(order)], fontsize=TICK)
    ax.set_xlabel("stratum mean log$_2$ fold change,\nRFO biosynthetic arm")
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_title("Every concurrent-control stratum", fontsize=BASE)
    letter(ax, "b")

    fig.tight_layout()
    p = os.path.join(out, "tissue_subarm_rfo.png")
    fig.savefig(p)
    plt.close(fig)
    return p


# -------------------------------------------------------------------------- 10
STEP_SHORT = {"rfo_galactinol_synthase": "GolS",
              "rfo_raffinose_stachyose_synthase": "RS/STS",
              "rfo_alpha_galactosidase_catabolic": "$\\alpha$-Gal",
              "inositol_1P_synthase": "MIPS",
              "inositol_monophosphatase": "IMP",
              "rfo_biosynthetic": "biosynthetic arm",
              "gt8_control": "GT8 control"}


def fig_coordination(meta_dir, sets_file, out):
    C = pd.read_csv(os.path.join(meta_dir, "PUBLIC_panel_correlation_matrix.csv"), index_col=0)
    R = pd.read_csv(os.path.join(meta_dir, "PUBLIC_coordination_test.csv"))
    X = pd.read_csv(os.path.join(meta_dir, "PUBLIC_coordination_between_sets.csv"))
    sets = pd.read_csv(sets_file)

    steps = ["rfo_galactinol_synthase", "rfo_raffinose_stachyose_synthase",
             "rfo_alpha_galactosidase_catabolic", "inositol_1P_synthase",
             "inositol_monophosphatase"]
    order, bounds = [], []
    for s in steps:
        mem = [g for g in sets[sets["set"] == s].gene if g in C.index]
        order += mem
        bounds.append(len(order))
    Cm = C.loc[order, order].to_numpy(dtype=float)

    fig = plt.figure(figsize=(7.2, 3.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.05], wspace=0.55)

    ax = fig.add_subplot(gs[0, 0])
    im = ax.imshow(Cm, cmap="RdBu_r", vmin=-1, vmax=1)
    for b in bounds[:-1]:
        ax.axhline(b - 0.5, color="k", lw=0.9)
        ax.axvline(b - 0.5, color="k", lw=0.9)
    lab = [g.replace("Solyc", "") for g in order]
    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    ax.set_xticklabels(lab, rotation=90, fontsize=4.2)
    ax.set_yticklabels(lab, fontsize=4.2)
    ax.tick_params(length=1.5, pad=1)
    prev = 0
    for s, b in zip(steps, bounds):
        ax.text(len(order) - 0.2, (prev + b - 1) / 2, STEP_SHORT[s], fontsize=TICK - 0.5,
                va="center", ha="left")
        prev = b
    cb = fig.colorbar(im, ax=ax, fraction=0.042, pad=0.22, location="left")
    cb.set_label("Pearson $r$ across strata", fontsize=TICK)
    cb.ax.tick_params(labelsize=TICK - 1)
    cb.ax.yaxis.set_label_position("left")
    cb.ax.yaxis.set_ticks_position("left")
    ax.set_title("Pairwise co-variation of the 21 panel genes", fontsize=BASE)
    letter(ax, "a")

    ax = fig.add_subplot(gs[0, 1])
    items = []
    for s in ["rfo_biosynthetic"] + steps + ["gt8_control"]:
        r = R[R.gene_set == s]
        if len(r):
            items.append((STEP_SHORT[s], r.iloc[0], False))
    # the label must name BOTH sides: not every declared pair starts from the arm
    cross_short = dict(STEP_SHORT, rfo_biosynthetic="arm")
    for _, r in X.iterrows():
        items.append((f"{cross_short[r.set_a]} \u00d7 {cross_short[r.set_b]}", r, True))
    n = len(items)
    for i, (lbl, r, cross) in enumerate(items):
        y = n - 1 - i
        lo = r.null_mean - 1.96 * r.null_sd
        hi = r.null_mean + 1.96 * r.null_sd
        ax.plot([lo, hi], [y, y], color="0.78", lw=4.5, solid_capstyle="butt", zorder=1)
        sig = r.p_bh_within_family < 0.05
        col = NEG if r.mean_pairwise_r < 0 else POS
        ax.plot([r.mean_pairwise_r], [y], "o", ms=5.5, color=col,
                mfc=col if sig else "white", mew=1.2, zorder=3)
        neg = r.mean_pairwise_r < 0
        ax.annotate(f"{r.mean_pairwise_r:+.2f}", (r.mean_pairwise_r, y),
                    xytext=(-6 if neg else 6, 0), textcoords="offset points",
                    fontsize=TICK, va="center", ha="right" if neg else "left", color=col)
    ax.axvline(0, color="0.4", lw=0.8, zorder=0)
    ax.axhline(len(X) - 0.5, color="0.6", lw=0.7, ls=":")
    ax.set_yticks(range(n))
    ax.set_yticklabels([l for l, _, _ in items][::-1], fontsize=TICK)
    ax.set_xlabel("mean pairwise $r$ of per-stratum log$_2$ fold change")
    ax.set_xlim(-0.75, 1.12)
    # one empty slot below the last row so the legend cannot sit on a data point
    ax.set_ylim(-1.9, n - 0.3)
    ax.set_title("Observed against an expression-matched null", fontsize=BASE)
    hs = [plt.Line2D([], [], color="0.78", lw=4.5, label="95% of the null"),
          plt.Line2D([], [], marker="o", ls="", color="0.3", label="BH < 0.05 in its family"),
          plt.Line2D([], [], marker="o", ls="", mfc="white", mec="0.3", color="0.3",
                     label="not significant")]
    ax.legend(handles=hs, loc="lower right", fontsize=TICK - 0.5, handletextpad=0.4,
              borderpad=0.3, labelspacing=0.35)
    letter(ax, "b")

    p = os.path.join(out, "rfo_coordination.png")
    fig.savefig(p)
    plt.close(fig)
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--de", required=True)
    ap.add_argument("--subgroup", required=True)
    ap.add_argument("--genotype", required=True)
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--sets", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)
    apply_style()

    made = [
        fig_cohort(args.de, args.out),
        fig_rfo(args.meta, args.panel, args.out),
        fig_gols(args.meta, args.out),
        fig_subgroup(args.subgroup, args.meta, args.panel, args.out),
        fig_duration(args.meta, args.panel, args.out),
        fig_genotype(args.genotype, args.out),
        fig_tissue(args.tissue, args.out),
        fig_coordination(args.meta, args.sets, args.out),
    ]
    for p in made:
        print("wrote", p)


if __name__ == "__main__":
    main()
